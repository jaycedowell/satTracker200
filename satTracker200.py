#!/usr/bin/env python2.7

import os
import sys
import copy
import time
import math
import ephem
import curses
import getopt
import serial
import urllib
import threading
from datetime import datetime, timedelta
from ConfigParser import SafeConfigParser


def usage(exitCode=None):
	print """satTracker200.py - Satellite predictor/tracker for the LX200 classic 
telescope.

Usage: satTracker200.py [OPTIONS] tle_file|download [tle_file [...]]

Options:
-h, --help           Display this help message
-c, --config-file    Site configuration file to use 
                     (default = satTracker200.config)
-p, --predict        Run in predictor mode (default = no, tracking mode)
-d, --date           In predictor mode, the start date as YYYY/MM/DD for the
                     predictions (default = now)
-t, --time           In predictor mode, the start time as HH:MM:SS for the 
                     predictions (default = now)
-i, --interval       In predictor mode, the duration in minutes to create 
                     predictions for (default = 180)
-m, --mag-limit      In predictor mode, filters the event list of passes
                     brighter thant the provided limit (default = no limit)
                     
Note:  If the first argument to the command is "download", the script will
       connect to the CelesTrak website and download the latest versions of
       the 'visual.txt' and 'science.txt' TLEs to the current directory and
       parse those.
"""
	
	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	config['configFile'] = 'satTracker200.config'
	config['predictorMode'] = False
	config['predictorDate'] = None
	config['predictorTime'] = None
	config['predictorInterval'] = 180.0
	config['predictorMagLimit'] = None
	
	try:
		opts, args = getopt.getopt(args, "hc:pd:t:i:m:", ["help", "config-file=", "predict", "date=", "time=", "interval=", "mag-limit="])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
	
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		elif opt in ('-c', '--config-file'):
			config['configFile'] = str(value)
		elif opt in ('-p', '--predict'):
			config['predictorMode'] = True
		elif opt in ('-d', '--date'):
			config['predictorDate'] = value
		elif opt in ('-t', '--time'):
			config['predictorTime'] = value
		elif opt in ('-i', '--interval'):
			config['predictorInterval'] = float(value)
		elif opt in ('-m', '--mag-limit'):
			config['predictorMagLimit'] = float(value)
		else:
			assert False
			
	# Add in arguments
	config['args'] = args
	
	# Parse the configuration file
	cFile = SafeConfigParser()
	cFile.read(config['configFile'])
	## Site parameters
	config['lat'] = cFile.getfloat('Site', 'latitude')
	config['lon'] = cFile.getfloat('Site', 'longitude')
	config['elev'] = cFile.getfloat('Site', 'elevation')
	config['utcOffset'] = cFile.getfloat('Site', 'utcOffset')
	## Serial parameters
	config['port'] = cFile.get('LX200', 'port')
	config['baud'] = cFile.get('LX200', 'baudrate')
	## Tracking parameters
	config['updateInterval'] = cFile.getfloat('Tracking', 'updateInterval')
	config['trackOffsetStep'] = cFile.getfloat('Tracking', 'trackOffsetStep')
	config['prepOffsetStep'] = cFile.getfloat('Tracking', 'perpOffsetStep')
	## Convert the trackOffsetStep from a float to a timedelta instance
	s = int(config['trackOffsetStep'])
	m = int((s-config['trackOffsetStep'])*1e6)
	config['trackOffsetStep'] = timedelta(seconds=s, microseconds=m)
	
	# Validate
	if len(config['args']) == 0:
		raise RuntimeError("Must provide at least one TLE file to use")
	if config['predictorDate'] is None and config['predictorTime'] is not None:
		raise RuntimeError("Must enter both a time and a date for the predictor mode")
	if config['predictorDate'] is not None and config['predictorTime'] is None:
		raise RuntimeError("Must enter both a date and a time for the predictor mode")
		
	# Return configuration
	return config


def getBearingFromPoint(ra1, dec1, ra2, dec2, ):
	"""
	Given ephem.Angle instances for two sets of RA/dec values, compute the 
	initial bearing and distance of the second point from the first.
	
	From:
	http://www.movable-type.co.uk/scripts/latlong.html
	"""
	
	# Distance (via PyEphem)
	distance = ephem.separation((ra1,dec1), (ra2,dec2))
	
	# Bearing
	bearing = math.atan2(math.sin(ra2-ra1)*math.cos(dec2), math.cos(dec1)*math.sin(dec2)-math.sin(dec1)*math.cos(dec2)*math.cos(ra2-ra1))
	
	return ephem.degrees(bearing), ephem.degrees(distance)


def getPointFromBearing(ra, dec, bearing, distance):
	"""
	Given ephem.Angle instances for RA, dec, calculate the RA/Dec pair that 
	lies at bearing 'bearing' a distance 'distance' away.
	
	From:
	http://www.movable-type.co.uk/scripts/latlong.html
	"""
	
	newDec = math.asin( math.sin(dec)*math.cos(distance) + math.cos(dec)*math.sin(distance)*math.cos(bearing) )
	newRA  = ra + math.atan2(math.sin(bearing)*math.sin(distance)*math.cos(dec), math.cos(distance) - math.sin(dec)*math.sin(newDec))
	if abs(math.cos(dec)) < 1e-12:
			newRA = bearing
			
	return ephem.hours(newRA), ephem.degrees(newDec)


class LX200(object):
	"""
	Minimal RS-232 interface to the LX200 classic telescope.
	"""
	
	def __init__(self, device, baud=9600, timeout=1.0):
		"""
		Create the connection to the telescope given a device name.  
		Optionally, set the baud rate and timeout for the connection.
		"""
		
		# Fixed parameters
		self.device = device
		self.timeout = timeout
		
		# Open the port
		self._open(baud=9600)
		
		# Change the baud rate, if necessary
		if baud != 9600:
			self.setBaudRate(baud)
			
		# The default state is for low-precision (HH:MM.M and sDD*MM) 
		# coordinates
		self.highPrecision = False
		
		# Semaphore for backgrounded moveToPosition() calls
		self.lock = threading.Semaphore()
		
	def _open(self, baud=9600):
		"""
		Internal function to open the serial port and make sure the
		telescope is there.
		"""
		
		# Save the baud rate in case a setBaudRate() is called
		self.baud = baud
		
		# Open the port
		self.port = serial.Serial(self.device, baud, timeout=self.timeout)
		
		# Make sure there is a telescope on the other end
		self.port.write(chr(0x06))
		status = self.port.read(1)
		if status not in ('L', 'A', 'P'):
			raise RuntimeError("Error communicating with the LX200 on '%s'" % self.device)
			
	def _ra2str(self, ra, deg=False):
		"""
		Internal function to take an RA value in HH.HHHHH or DDD.DDDDD 
		format and convert it to a HH:MM.M or HH:MM:SS string.
		"""
		
		if deg:
			ra /= 15.0
		ra %= 24.0
		
		if self.highPrecision:
			raH, raM, raS = int(ra), int(ra*60)%60, int(ra*3600)%60
			value = "%02i:%02i:%02i" % (raH, raM, raS)
		else:
			raH, raM = int(ra), (ra*60)%60
			value = "%02i:%04.1f" % (raH, raM)
			
		return value
		
	def _dec2str(self, dec):
		"""
		Internal function to take an dec value in sDD.DDDDD format and
		# convert it to a sDD*MM or sDD*MM:SS string.
		"""
		
		decG = '-' if dec < 0 else '+'
		if self.highPrecision:
			decD, decM, decS = int(abs(dec)), int(abs(dec)*60)%60, int(abs(dec)*3600)%60
			value = "%s%02i%s%02i:%02i" % (decG, decD, chr(223), decM, decS)
		else:
			decD, decM = int(abs(dec)), (abs(dec)*60)%60
			value = "%s%02i%s%02i" % (decG, decD, chr(223), decM)
			
		return value
		
	def _str2coord(self, data, forceLow=False):
		"""
		Internal function to take a coordinate string (either low or high
		precision) and convert it to a HH.HHHHHH or sDD.DDDDDD value.
		"""
		
		data = data.replace(chr(223), ':')
		if self.highPrecision and not forceLow:
			d, m, s = data.split(':', 2)
			d, m, s = int(d), int(m), int(s)
			value = abs(d) + m/60. + s/3600.
		else:
			d, m = data.split(':', 1)
			d, m = int(d), float(m)
			value = abs(d) + m/60.0
		if d < 0:
			value *= -1
		return value
		
	def _readNumber(self):
		"""
		Internal function to read a single byte from the serial port.
		"""
		
		c = self.port.read(1)
		try:
			c = int(c, 10)
		except ValueError:
			pass
		return c
		
	def _readString(self):
		"""
		Read a set of bytes from the serial port until a '#' is found.
		The data up to but not including the '#' are returned.
		"""
		
		c = self.port.read(1)
		try:
			while c[-1] != '#':
				c += self.port.read(1)
		except IndexError:
			pass
		return c[:-1]
		
	def setBaudRate(self, baud):
		"""
		Change the baud rate for the RS-232 connection to the telescope.
		Returns True if the change was successful, False otherwise.
		"""
		
		rates = {56700: 1, 38400: 2, 28800: 3, 19200: 4, 
			    14400: 5,  9600: 6,  4800: 7,  2400: 8, 
			    1200: 9}
			    
		if baud == self.baud:
			return True
			
		try:
			self.port.write('#:SB%i#' % rates[baud])
			status = bool( self._readNumber() )
			if status:
				self.port.close()
				self._open(baud=baud)
				return True
			else:
				return False
				
		except KeyError:
			return False
			
	def getDateTime(self):
		"""
		Get the date/time according to the telescope as a datetime 
		instance in UTC.
		"""
		
		self.port.write('#:GC#')
		tDate = self._readString()
		self.port.write('#:GL#')
		tTime = self._readString()
		self.port.write('#:GG#')
		tOff = self._readString()
		
		dt = datetime.strptime("%s %s" % (tDate, tTime), "%m/%d/%y %H:%M:%S")
		
		tOff = float(tOff)
		gOff = -1 if tOff < 0 else 1
		hOff = int(abs(tOff))
		mOff = int(abs(tOff)*60) % 60
		sOff = int(abs(tOff)*3600) % 60
		td = timedelta(hours=hOff, minutes=mOff, seconds=sOff)
		dt += gOff*td
		
		return dt
		
	def setDateTime(self, dt):
		"""
		Set the date/time on the telescope using the provided UTC 
		datetime instance.
		
		Note:  This function overrides any existing UTC offset and
		forces it to 0.
		"""
		
		tDate = dt.strftime("%m/%d/%y")
		tTime = dt.strftime("%H:%M:%S")
		
		self.port.write('#:SC%s#' % tDate)
		statusD = bool( self._readNumber() )
		if statusD:
			msg = self._readString()
		self.port.write('#:SL%s#' % tTime)
		statusT = bool( self._readNumber() )
		self.port.write('#:SG+00.0#')
		statusO = bool( self._readNumber() )
		
		return statusD and statusT and statusO
		
	def getObserver(self, elevation=1706.9):
		"""
		Query the telescope for its current location and return it as an
		ephem.Observer instance.  The 'elevation' keyword is used to set
		the elevation of the location, in meters, above sea level.
		"""
		
		self.port.write('#:Gt#')
		tLat = self._readString()
		self.port.write('#:Gg#')
		tLng = self._readString()
		
		tLat = self._str2coord(tLat, forceLow=True)
		tLng = self._str2coord(tLng, forceLow=True)*-1.0		# For the LX200, W is positive
		
		obs = ephem.Observer()
		obs.lat = tLat*math.pi/180
		obs.lon = tLng*math.pi/180
		obs.elevation = elevation
		
		return obs
		
	def setHighPrecision(self):
		if not self.highPrecision:
			self.port.write('#:U#')
			self.highPrecision = True
			
	def setLowPrecision(self):
		if self.highPrecision:
			self.port.write('#:U#')
			self.highPrecision = False
			
	def setSlewRateFast(self):
		self.port.write('#:RS#')
		
	def setSlewRateFind(self):
		self.port.write('#:RM#')
		
	def setSlewRateCenter(self):
		self.port.write('#:RC#')
		
	def setSlewRateGuide(self):
		self.port.write('#:RG#')
		
	def getCurrentPointing(self):
		"""
		Query the telescope for the current poitning location and return
		it as a two-element tuple of RA, in decimal hours, and dec, in
		decimal degrees.
		"""
		
		self.port.write('#:GR#')
		ra = self._readString()
		self.port.write('#:GD#')
		dec = self._readString()
		
		ra = self._str2coord(ra)
		dec = self._str2coord(dec)
		return ra, dec
		
	def _moveToPosition(self, raStr, decStr, fast=False):
		"""
		Internal function to help with running moveToPosition() in the 
		background.
		"""
		
		if self.lock.acquire(False):
			status = True
			self.port.write('#:Sr%s#' % raStr)
			if not fast:
				status &= bool( self._readNumber() )
			self.port.write('#:Sd%s#' % decStr)
			if not fast:
				status &= bool( self._readNumber() )
				
			if status:
				self.port.write('#:MS#')
				if not fast:
					status &= not bool( self._readNumber() )
				if not status:
					msg = self._readString()
				else:
					if not fast:
						self.port.write('#:D#')
						dist = self._readString()
						dist = sum([1 for c in dist if c != ' '])
						while dist > 2:
							time.sleep(0.001)
							self.port.write('#:D#')
							dist = self._readString()
							dist = sum([1 for c in dist if c != ' '])
			self.lock.release()
			
		else:
			status = False
			
		return status
		
	def moveToPosition(self, ra, dec, fast=False, background=False):
		"""
		Given a RA value in decimal hours and a declination value in 
		decimal degrees, slew the telescope there provided that it is
		observable.  Returns True if the target was set the and the
		move was initiated, False otherwise.
		
		There are two keywords that control how this functions behaves:
		'fast' and 'background'.  Setting 'fast' to True disables the 
		on-line error detection in order to increase the command rate.
		'background' causes the function to spawn a background thread
		for the movement so that this function is non-blocking.
		"""
		
		ra = self._ra2str(ra)
		dec = self._dec2str(dec)
		
		if background:
			t = threading.Thread(target=self._moveToPosition, args=(ra,dec), kwargs={'fast':fast})
			t.start()
			status = True
			
		else:
			status = self._moveToPosition(ra, dec, fast=fast)
			
		return status
		
	def haltCurrentSlew(self):
		"""
		Stop the current slew.
		"""
		
		self.port.write('#:Q#')


class EarthSatellitePlus(ephem.EarthSatellite):
	"""
	Sub-class of the ephem.EarthSatellite class that adds in additional
	information about the magnitude of the satellite.
	
	.. note::  The 'bearing' attribute is only computed relative to the 
	the previous call to compute().  This field is only useful during 
	tracking operations.
	"""
	
	# Rising/setting
	rising = False
	
	# Bearing information
	bearing = 0.0
	_lastRA = None
	_lastDec = None
	_lastAlt = None
	
	# Magnitude information
	stdMag = 99.0
	_sol = ephem.Sun()
	magnitude = 99.0
	
	def compute(self, observer):
		"""
		Wrapper around ephem.EarthSatellite.compute() function that 
		adds in visibility, magntiude, and satellite bearing 
		computations.
		"""
		
		## Save the 
		try:
			self._lastRA, self._lastDec, self._lastAlt = self.ra, self.dec, self.alt
		except (RuntimeError, AttributeError):
			self._lastRA, self._lastDec, self._lastAlt = None, None, None
			
		## Locate the satellite
		ephem.EarthSatellite.compute(self, observer)
		
		## Compute the location of the Sun so that we can estimate the 
		## brightness of the satellite.
		self._sol.compute(observer)
		B = math.pi - ephem.separation(self._sol, self)
		self.magnitude = self.stdMag - 15.0 + 5.0*math.log10(self.range/1e3)
		self.magnitude -= 2.5*math.log10( math.sin(B) + (math.pi-B)*math.cos(B) )
		if self.magnitude >= 15.0 or self.eclipsed:
			self.magnitude = 15.0
			
		## Calculate the bearing if this is not the first time
		## that we have tried
		try:
			self.bearing = getBearingFromPoint(self._lastRA, self._lastDec, self.ra, self.dec)[0]
			self.rising = True if self.alt > self._lastAlt else False
		except (ValueError, TypeError):
			self.bearing = 0.0
			
	def setStdMag(self, stdMag):
		"""
		Set the "standard magntiude" for this satellite to use for 
		brightness calculations.
		"""
		
		self.stdMag = stdMag
		
	def fillFromPyEphem(self, sat):
		"""
		Initialize this instance using an existing 
		ephem.EarthSatellite instance.
		"""
		
		for attr in ('name', 'catalog_number', '_epoch', '_n', '_orbit', '_drag', '_decay', '_e'):
			setattr(self, attr, getattr(sat, attr, None))
		for attr in ('_inc', '_raan', '_ap', '_M'):
			setattr(self, attr, getattr(sat, attr, None)*180/math.pi)


def passPredictor(observer, satellites, date=None, time=None, utcOffset=0.0, duration=180.0, magnitudeCut=None):
	"""
	Function to predict passes of satellites and return them as a list of 
	tuples.  The entry in each tuple is:
	  * satellite name
	  * maximimum brightness
	  * rise local date/time
	  * rise azimuth
	  * maximum altitude local date/time
	  * maximum altitude
	  * set local date/time
	  * set azimuth
	"""
	
	# Deal with the magntiude magnitude
	if magnitudeCut is None:
		magnitudeCut = 15.0
		
	# Set the date/time
	if date is None or time is None:
		## Current date/time in local (UTC + proivded UTC offset)
		dt = datetime.utcnow()
		oG = -1 if utcOffset < 0 else 1
		oS = abs(utcOffset*3600.0)
		oM = int((oS-int(oS))*1e6)
		oS = int(oS)
		dt += oG*timedelta(seconds=oS, microseconds=oM)
	else:
		## Input date/time
		dt = datetime.strptime("%s %s" % (date, time), "%Y/%m/%d %H:%M:%S")
	## Save the original observer date so that we can restore it later
	origDate = 1.0*observer.date
	## Set the date to local
	observer.date = dt.strftime("%Y/%m/%d %H:%M:%S")
	## Correct for the provided UTC offset
	observer.date -= utcOffset/24.0
	
	# Save the start time so that we can keep up with it
	tStart = observer.date*1.0
	
	# Loop through satellites to come up with a sequence of events
	events = []
	for sat in satellites:
		## Reset the start time
		observer.date = tStart
		
		## Run until we are outside of the window to search
		while (observer.date-tStart) <= (duration/60.0/24.0):
			### Compute the next pass
			try:
				rTime, rAz, mTime, mEl, sTime, sAz = observer.next_pass(sat)
			except ValueError:
				#### If this satellite is circumpolar, set some dummy values
				rTime, rAz = observer.date, ephem.degrees('0:00:00')
				mTime, mEl = observer.date+duration/2.0/60.0/24.0, ephem.degrees('1:00:00')
				sTime, sAz = observer.date+duration/60.0/24.0, ephem.degrees('0:00:00')
				
			### Does the satellite actually rise?
			if rTime is None or mTime is None or sTime is None:
				observer.date += 2*duration/60.0/24.0
				continue
				
			### Make sure that this pass is valid
			visible = False
			maxBright = 15.0
			for d in (rTime, mTime, sTime):
				#### Compute the satellite's particulars for this time
				observer.date = d
				sat.compute(observer)
				
				#### Is it visible?
				visible |= not sat.eclipsed
				
				#### Is it brighter than our initial guess?
				if sat.magnitude < maxBright:
					maxBright = sat.magnitude
					
			if visible and maxBright <= magnitudeCut and (rTime-tStart) <= (duration/60.0/24.0):
				### Save this one and update the time to local
				rTimeLocal = ephem.date( rTime + utcOffset/24.0 )
				mTimeLocal = ephem.date( mTime + utcOffset/24.0 )
				sTimeLocal = ephem.date( sTime + utcOffset/24.0 )
				events.append( (sat.name, maxBright, rTimeLocal, rAz, mTimeLocal, mEl, sTimeLocal, sAz) )
				
			### Step forward in time to look for the next pass
			observer.date = sTime + 5/60.0/24.0
			
	# Sort the passes by rise time
	events.sort(key=lambda x: x[2])
	
	# Restore the original date
	observer.date = origDate
	
	# Return
	return events


class SatellitePositionTracker(object):
	"""
	Class for tracking an ensemble of satellites and using a telescope to see them.
	"""
	
	def __init__(self, observer, satellites, interval=0.1, lx200=None):
		"""
		Initialize the SatellitePositions instace using an ephem.Observer
		instance, a list of ephem.EarthSatellite instances, and, optionally:
		  * An update interval in seconds using the 'interval' keyword and
		  * A LX200 instance for tracking using the 'lx200' keyword.
		"""
		
		# Copy the observer and convert the list of ephem.EarthSatellite 
		# instances to EarthSatellitePlus instances.  Also, save the
		# interval and telescope information
		self.observer = observer
		self.satellites = []
		for sat in satellites:
			if type(sat) is not EarthSatellitePlus:
				newSat = EarthSatellitePlus()
				newSat.fillFromPyEphem(sat)
				sat = newSat
			self.satellites.append( sat )
		self.interval = float(interval)
		self.lx200 = lx200
		
		# State variables for fine-tuning the telescope tracking
		self.tracking = None
		self.timeOffset = timedelta()
		self.perpOffset = 0.0
		
		# Threading state variales
		self.thread = None
		self.alive = threading.Event()
		
	def updateObserver(self, observer):
		"""
		Change the location of the observer using an ephem.Observer
		instance.  If the background tracking thread is running is
		also restarts the thread.
		"""
		
		self.observer = copy.copy(observer)
		
		if self.thread is not None:
			self.stop()
			self.start()
			
	def setTimeOffset(self, seconds):
		"""
		Set the tracking time offset to using a timedelta instance.
		"""
		
		self.timeOffset = seconds
		
	def adjustTimeOffset(self, seconds):
		"""
		Apply an additional adjustment to the tracking time offset 
		using a timedelta instance.
		"""
		
		self.timeOffset += seconds
		
	def getTimeOffset(self):
		"""
		Return the current tracking time offset as a timedelta instance.
		"""
		
		return self.timeOffset
		
	def setPerpOffset(self, degrees):
		"""
		Set the perpendicular track offset to a value in decimal degrees.
		"""
		
		self.perpOffset = degrees*math.pi/180
		
	def adjustPerpOffset(self, degrees):
		"""
		Apply an additional adjustment to the perpendicular track offset
		used a value in decimal degrees.
		"""
		
		self.perpOffset += degrees*math.pi/180
		
	def getPerpOffset(self):
		"""
		Return the perpendicular track offset in decimal degrees.
		"""
		
		return self.perpOffset*180/math.pi
		
	def startTracking(self, catalog_number):
		"""
		Begin tracking the satellite with the specified NORAD ID number.
		"""
		
		self.tracking = catalog_number
		
	def stopTracking(self):
		"""
		Halt satellite tracking and stop the current telescope slew.
		"""
		
		self.tracking = None
		if self.lx200 is not None:
			self.lx200.haltCurrentSlew()
			
	def getTracking(self):
		"""
		Return the NORAD ID number of the satellite being currently 
		tracked.  If nothing is being tracked return None.
		"""
		
		return self.tracking
		
	def start(self):
		"""
		Start the tracking background thread.
		"""
		
		if self.thread is not None:
			self.stop()
			
		self.thread = threading.Thread(target=self.trackSats, name='trackSats')
		self.thread.setDaemon(1)
		self.alive.set()
		self.thread.start()
		time.sleep(1)
		
	def stop(self):
		"""
		Stop the tracking background thread.
		"""
		
		if self.thread is not None:
			self.alive.clear()
			
			self.thread.join()
			self.thread = None
			
	def trackSats(self):
		"""
		Background thread function that computes the locations of all of 
		the satellites provided and provides poitning information to the
		telescope for the satellite being tracked.
		"""
		
		while self.alive.isSet():
			# Start and current times.  The current time is further adjusted
			# by the 'timeOffset' attribute to help with tracking.
			tStart = time.time()
			tNow = datetime.utcnow() + self.timeOffset
			
			# Update the observer with the current time
			self.observer.date = tNow.strftime("%Y/%m/%d %H:%M:%S.%f")
			
			# Loop through the satellites and update them
			for sat in self.satellites:
				sat.compute(self.observer)
				
				if sat.alt > 0:
					## If the satellie is up, check and see if it is 
					## the one that we should be tracking.
					if sat.catalog_number == self.tracking:
						### Is there a telescope to use?
						if self.lx200 is not None:
							#### Apply a perpendicular correction to the
							#### track to help with tracking
							ra, dec = getPointFromBearing(sat.ra, sat.dec, sat.bearing+math.pi/2, self.perpOffset)
							
							#### Radians -> hours/degrees
							ra = ra*12/math.pi
							dec = dec*180/math.pi
							
							#### Command the telescope
							self.lx200.moveToPosition(ra, dec, fast=True, background=True)
							
				else:
					## If it is no longer visible, check and see if it is
					## the satellite that we were tracking so that we can
					## stop.
					if sat.catalog_number == self.tracking:
						self.tracking = None
						
			# Final time to figure out how much time we spent calculating 
			# positions and how long we need to sleep to get ready for the
			# next iteration
			tStop = time.time()
			sleepCount = 0.0
			sleepTime = self.interval - (tStop - tStart)
			
			# Sleep in 10 ms increments
			while self.alive.isSet() and sleepCount < sleepTime:
				time.sleep(0.01)
				sleepCount += 0.01


def main(args):
	# Parse the command line
	config = parseOptions(args)
	
	# Grab the list of filenames to look for TLEs in
	filenames = config['args']
	
	# Special case:  If the first filename in 'download', go to CelesTrak and
	# download the latest
	if filenames[0] == 'download':
		filenames = []
		
		## Top ~100 brightest satellites
		try:
			urllib.urlretrieve("http://celestrak.com/NORAD/elements/visual.txt", "visual.txt")
			filenames.append( 'visual.txt' )
		except Exception as e:
			print "ERROR: Download of 'visual.txt' failed: %s" % str(e)
			
		## Earth & space science satellites
		try:
			urllib.urlretrieve("http://celestrak.com/NORAD/elements/science.txt", "science.txt")
			filenames.append( 'science.txt' )
		except Exception as e:
			print "ERROR: Download of 'science.txt' failed: %s" % str(e)
			
	# Read in the TLEs
	satellites = []
	for filename in filenames:
		fh = open(filename, 'r')
		data = fh.readlines()
		fh.close()
		
		## Convert to EarthSatellitePlus instances
		for i in xrange(len(data)/3):
			sat = EarthSatellitePlus()
			sat.fillFromPyEphem( ephem.readtle(data[3*i+0], data[3*i+1], data[3*i+2]) )
			satellites.append( sat )
			
	# Load in the standard magnitude if we happen to have a qs.mag file to read
	## NORAD ID -> list index dictionary to help with mapping the
	## QuickSat qs.mag file information to the right EarthSatellitePlus
	## instance.
	lookup = {}
	for i in xrange(len(satellites)):
		lookup[satellites[i].catalog_number] = i
		
	## Try to added the QuickSat qs.mag information to the various
	## instances.
	try:
		fh = open('qs.mag', 'r')
		for line in fh:
			cNum = int(line[:5], 10)
			if cNum == 1 or cNum == 99999:
				continue
			try:
				mag = float(line[33:37])
			except ValueError:
				pass
			try:
				satellites[ lookup[cNum] ].setStdMag( mag )
			except KeyError:
				pass
		fh.close()
	except IOError:
		pass
		
	# Report on what we found
	nSats = len(satellites)
	if nSats == 0:
		print "ERROR: Found no valid TLEs, exiting"
		sys.exit(1)
		
	print "Loaded %i TLEs from %i files" % (nSats, len(filenames))
	
	# Try and setup the telescope.  This not only opens up the port but also
	# makes sure that the time is right.
	try:
		tel = LX200(config['port'], baud=config['baud'])
		
		## Make sure that the time is "reasonable"
		tTime = tel.getDateTime()
		cTime = datetime.utcnow()
		tcOffset = cTime - tTime
		
		if abs(tcOffset) > timedelta(seconds=10):
			tel.setDateTime(cTime)
			
			tTime = tel.getDateTime()
			cTime = datetime.utcnow()
			tcOffset = cTime - tTime
			
		print "LX200 set to %s, computer at %s" % (tTime, cTime)
		print "-> Difference is %s" % tcOffset
		
		# Set the slew rate
		tel.setSlewRateFast()
		
		# Set high precision coordinates
		tel.setHighPrecision()
		
	except Exception as e:
		print "ERROR: %s" % str(e)
		tel = None
		tcOffset = timedelta()
		
	# Set the observer according to the telescope, if found.  Otherwise, use
	# the default values.
	try:
		obs = tel.getObserver(elevation=config['elev'])
	except Exception as e:
		print "ERROR: %s" % str(e)
		obs = ephem.Observer()
		obs.lat = config['lat']*math.pi/180
		obs.lon = config['lon']*math.pi/180
		obs.elevation = config['elev']
	print "Observer set to %.4f %s, %.4f %s @ %.1f m" % (abs(obs.lat)*180/math.pi, 'N' if obs.lat >= 0 else 'S', abs(obs.lon)*180/math.pi, 'E' if obs.lon >= 0 else 'W', obs.elevation)
	
	# Select how to run
	if config['predictorMode']:
		# Predictor mode
		
		# Work with the date and time provided so that we can display 
		# something helpful with the pass predictions
		d, t = config['predictorDate'], config['predictorTime']
		if d is None and t is None:
			## Use the current UTC time and convert it to local using 
			## the provided UTC offset
			dt = datetime.utcnow()
			oG = -1 if config['utcOffset'] < 0 else 1
			oS = abs(config['utcOffset']*3600.0)
			oM = int((oS-int(oS))*1e6)
			oS = int(oS)
			dt += oG*timedelta(seconds=oS, microseconds=oM)
			dt = dt.strftime("%Y/%m/%d %H:%M:%S")
		else:
			## Use the provided date and time
			dt = "%s %s" % (d, t)
			
		# Compute the passes - this should return values in local time
		events = passPredictor(obs, satellites, date=d, time=t, utcOffset=config['utcOffset'],
					duration=config['predictorInterval'], magnitudeCut=config['predictorMagLimit'])
					
		print " "
		print "Satellite Pass Predictions for %s:" % dt
		print "%24s  %5s  %13s  %13s  %13s" % ('Name', 'Mag.', 'Rise Time', 'Max. El. Time', 'Set Time')
		for event in events:
			mag = event[1]
			if mag < 15:
				mag = "%5.1f" % mag
			else:
				mag = " --- "
			print "%24s  %5s  %13s  %13s  %13s" % (event[0], mag, str(event[2])[5:], str(event[4])[5:], str(event[6])[5:])
			
	else:
		# Tracking mode
		
		# Convert the provided UT offset to a timedelta instance that we 
		# can use to update the clock
		oG = -1 if config['utcOffset'] < 0 else 1
		oS = abs(config['utcOffset']*3600.0)
		oM = int((oS-int(oS))*1e6)
		oS = int(oS)
		utcOffset = oG*timedelta(seconds=oS, microseconds=oM)
		
		# Start the tracking thread
		trkr = SatellitePositionTracker(obs, satellites, interval=config['updateInterval'], lx200=tel)
		trkr.start()
		
		# Setup the TUI
		stdscr = curses.initscr()
		curses.noecho()
		curses.cbreak()
		stdscr.keypad(1)
		stdscr.nodelay(1)
		
		# Main TUI loop
		## State control variables and such
		act = 1
		trk = -1
		msg = ''
		oldMsg = ''
		msgCount = 0
		empty = ''
		while len(empty) < (5+2+24+2+12+2+12+2+12+2+5):
			empty += ' '
		empty += '\n'
		magLimit = 6.0
		
		while True:
			## Current time
			tNow = datetime.utcnow()
			tNowLocal = tNow + utcOffset
			
			## Figure out how many satellites are currently visible and brighter
			## than the current magntiude limit
			trkChange = False
			nVis = sum( [1 for sat in trkr.satellites if sat.alt > 0 and sat.magnitude <= magLimit] )
			
			## Interact with the user's keypresses
			c = stdscr.getch()
			if c == ord('q'):
				### 'q' to exit the main loop
				break
			elif c == curses.KEY_UP:
				### Up arrow to change what is selected
				act -= 1
				if act < 1:
					act = 1
			elif c == curses.KEY_DOWN:
				### Down array to change what is selected
				act += 1
				if act > min([nVis, 15]):
					act = min([nVis, 15])
			elif c == ord('a'):
				### 'a' to adjust the tracking time offset - negative
				trkr.adjustTimeOffset( -config['trackOffsetStep'] )
				off = trkr.getTimeOffset()
				msg = 'Time offset now %+.1f s' % (off.days*86400+off.seconds+off.microseconds/1e6)
			elif c == ord('s'):
				### 's' to adjust the tracking time offset - positive
				trkr.adjustTimeOffset( config['trackOffsetStep'] )
				off = trkr.getTimeOffset()
				msg = 'Time offset now %+.1f s' % (off.days*86400+off.seconds+off.microseconds/1e6)
			elif c == ord('z'):
				### 'z' to adjust the perpendicular track offset - negative
				trkr.adjustPerpOffset( -config['perpOffsetStep'] )
				off = trkr.getPerpOffset()
				msg = 'Perpendicular offset now %+.1f degrees' % off
			elif c == ord('w'):
				### 'w' to adjust the perpendicular track offset - positive
				trkr.adjustPerpOffset( config['perpOffsetStep'] )
				off = trkr.getPerpOffset()
				msg = 'Perpendicular offset now %+.1f degrees' % off
			elif c == ord('p'):
				off1 = trkr.getTimeOffset()
				off2 = trkr.getPerpOffset()
				msg = 'Time offset is %+.1f s; Perpendicular offset is %+.1f degrees' % (off1.days*86400+off1.seconds+off1.microseconds/1e6, off2)
			elif c == ord('k'):
				### 'k' to adjust the magnitude limit - negative
				magLimit -= 0.5
				magLimit = max([magLimit, -2.0])
				msg = 'Magntiude limit now <= %.1f mag' % magLimit
			elif c == ord('l'):
				### 'k' to adjust the magnitude limit - negative
				magLimit += 0.5
				magLimit = min([magLimit, 15.0])
				if magLimit == 15.0:
					msg = 'Magntiude limit now disabled'
				else:
					msg = 'Magntiude limit now <= %.1f mag' % magLimit
			elif c == ord('t'):
				### 't' to toggle telesope tracking on/off
				if trk == -1:
					trk = act
					trkChange = True
				else:
					trk = -1
					trkChange = True
					
			## Additional state check to make sure that we aren't still trying
			## to track a satellite that has set
			if not trkChange:
				if trk != -1 and trkr.getTracking() is None:
					trk = -1
					trkChange = True
					
			## Header for the satellite data
			output = "%5s  %24s  %12s  %12s  %12s  %5s\n" % ('ID', 'Name', 'Azimuth', 'Elevation', 'Status', 'Mag.')
			stdscr.addstr(0, 0, output, curses.A_UNDERLINE)
			
			## Satellite information
			k = 1
			for j in xrange(nSats):
				### Make the satellite easy-to-access
				sat = trkr.satellites[j]
				
				### Is it visible and bright enought?
				if sat.alt > 0 and sat.magnitude <= magLimit:
					#### Create the output line
					output = "%5i  %24s  %12s  %12s  %8s-%3s  %5s\n" % (sat.catalog_number, sat.name, sat.az, sat.alt, 'Eclipsed' if sat.eclipsed else 'In Sun', 'Asc' if sat.rising else 'Dsc', '%5.1f' % sat.magnitude if sat.magnitude < 15 else ' --- ')
					
					#### See if we need to enable/disable tracking
					if trkChange:
						if k == trk:
							if sat.catalog_number != trkr.getTracking():
								trkr.setTimeOffset( timedelta() )
								trkr.setPerpOffset( 0.0 )
								trkr.startTracking(sat.catalog_number)
								msg = 'Now tracking NORAD ID #%i' % sat.catalog_number
								
						elif trk == -1 and k == 1:
							trkr.stopTracking()
							msg = 'Tracking stopped'
							
					#### Set the text scheme to use for displaying
					if trkr.satellites[j].catalog_number == trkr.getTracking():
						## Tracking
						options = (curses.A_UNDERLINE, curses.A_REVERSE|curses.A_UNDERLINE)
					else:
						## Non-tracking
						options = (curses.A_NORMAL, curses.A_REVERSE)
						
					#### Add the line to the screen, provided it will fit
					if k < 14:
						if k == act:
							stdscr.addstr(k, 0, output, options[1])
						else:
							stdscr.addstr(k, 0, output, options[0])
						k += 1
					else:
						break
						
			## Pad out the message stored in 'msg' generated from interaction 
			## with the user/telescope
			while len(msg) < (5+2+24+2+12+2+12+2+12+2+5):
				msg += ' '
			if msg[-1] != '\n':
				msg += '\n'
				
			## See if it is time to clear off the message
			if msg == oldMsg:
				msgCount += 1
				if msgCount == 101:
					msg = empty
					msgCount = 0
			oldMsg = msg
			
			## Pad out the current UTC time
			output = "  UTC: "
			output += tNow.strftime("%Y/%m/%d %H:%M:%S")
			output += ".%01i" % (float(tNow.strftime("%f"))/1e5,)
			output += "                         "
			output += "LT: "
			output += tNowLocal.strftime("%Y/%m/%d %H:%M:%S")
			output += ".%01i" % (float(tNowLocal.strftime("%f"))/1e5,)
			while len(output) < (5+2+24+2+12+2+12+2+12+2+5):
				output += ' '
			output += '\n'
			
			## Final time/message/help information
			stdscr.addstr(k+0, 0, output)
			stdscr.addstr(k+1, 0, msg)
			stdscr.addstr(k+2, 0, 'Keys:\n')
			stdscr.addstr(k+3, 0, '  t   - Track currently selected satellite\n')
			stdscr.addstr(k+4, 0, '  p   - Print current tracking offsets\n')
			stdscr.addstr(k+5, 0, '  a/s - Decrease/Increase track lead by 1 second\n')
			stdscr.addstr(k+6, 0, '  z/w - Decrease/Increase track offset by 0.2 degrees\n')
			stdscr.addstr(k+7, 0, '  k/l - Decrease/Increase the magntiude limit by 0.5 mag\n')
			stdscr.addstr(k+8, 0, '  q   - Exit\n')
			
			## Pad out the window to deal with satellite setting
			while k+9 < 22:
				stdscr.addstr(k+9, 0, empty)
				k += 1
				
			## Display and sleep until the next iteration
			stdscr.refresh()
			time.sleep(0.1)
			
		# Close out the curses session
		curses.nocbreak()
		stdscr.keypad(0)
		curses.echo()
		curses.endwin()
		
		# Stop the tracking thread
		trkr.stop()


if __name__ == "__main__":
	main(sys.argv[1:])
	