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
from datetime import datetime, timedelta
from ConfigParser import SafeConfigParser
import traceback
try:
	import cStringIO as StringIO
except ImportError:
	import StringIO


_deg2rad = math.pi/180.0
_rad2deg = 180.0/math.pi
_rad2hr  = 12.0/math.pi


def usage(exitCode=None):
	print """satTracker200.py - Satellite predictor/tracker for the LX200 classic 
telescope.

Usage: satTracker200.py [OPTIONS] tle_file [tle_file [...]]

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
                     predictions for (default = 180 = 3 hours)
-m, --mag-limit      In predictor mode, filters the event list of passes
                     brighter than the provided limit (default = no limit)
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
	config['trackOffsetStep'] = abs( cFile.getfloat('Tracking', 'trackOffsetStep') )
	config['crossTrackOffsetStep'] = abs( cFile.getfloat('Tracking', 'crossTrackOffsetStep') )
	## Convert the trackOffsetStep from a float to a timedelta instance
	s = int(config['trackOffsetStep'])
	u = int((config['trackOffsetStep']-s)*1e6)
	config['trackOffsetStep'] = timedelta(seconds=s, microseconds=u)
	
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
	
	def __init__(self, device, baud=9600, timeout=None):
		"""
		Create the connection to the telescope given a device name.  
		Optionally, set the baud rate and timeout for the connection.
		"""
		
		# Fixed parameters
		self._device = device
		self._timeout = timeout
		
		# Open the port
		self._open(baud=9600)
		
		# Change the baud rate, if necessary
		if baud != 9600:
			self.setBaudRate(baud)
			
		# The default state is for low-precision (HH:MM.M and sDD*MM) 
		# coordinates
		self._highPrecision = False
		
	def _open(self, baud=9600):
		"""
		Internal function to open the serial port and make sure the
		telescope is there.
		"""
		
		# Save the baud rate in case a setBaudRate() is called
		self._baud = baud
		
		# Open the port
		self.port = serial.Serial(self._device, baudrate=baud, bytesize=serial.EIGHTBITS, parity=serial.PARITY_NONE,
								stopbits=serial.STOPBITS_ONE, xonxoff=False, rtscts=False, dsrdtr=None,
								timeout=self._timeout)
								
		# Clear it out
		self.port.flushInput()
		self.port.flushOutput()
		
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
		
		if self._highPrecision:
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
		if self._highPrecision:
			decD, decM, decS = int(abs(dec)), int(abs(dec)*60)%60, int(abs(dec)*3600)%60
			value = "%s%02i%s%02i:%02i" % (decG, decD, chr(223), decM, decS)
		else:
			decD, decM = int(abs(dec)), (abs(dec)*60)%60
			value = "%s%02i%s%02i" % (decG, decD, chr(223), decM)
			
		return value
		
	def _str2coord(self, data, forceLow=False, forceHigh=False):
		"""
		Internal function to take a coordinate string (either low or high
		precision) and convert it to a HH.HHHHHH or sDD.DDDDDD value.
		"""
		
		data = data.replace(chr(223), ':')
		if (self._highPrecision and not forceLow) or forceHigh:
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
			    
		if baud == self._baud:
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
		obs.lat = tLat*_deg2rad
		obs.lon = tLng*_deg2rad
		obs.elevation = elevation
		
		return obs
		
	def setHighPrecision(self):
		"""
		Set the coordinate format to high precision mode (HH:MM:SS or 
		sDD:MM:SS).
		"""
		
		# Are we already there?
		if not self._highPrecision:
			self.port.write('#:U#')
			self._highPrecision = True
			
	def setLowPrecision(self):
		"""
		Set the coordinate format to low precision mode (default; HH:MM.M
		or sDD:MM).
		"""
		
		# Are the already there?
		if self._highPrecision:
			self.port.write('#:U#')
			self._highPrecision = False
			
	def setMaximumSlewRate(self, degPerSec):
		"""
		Set the maximum slew rate in degrees per second.  This is the speed
		used in the "Slew" (fastest) mode.  Valid values are: 
		  * 2, 
		  * 3, 
		  * 4, 
		  * 5, 
		  * 6, 
		  * 7, and
		  * 8.
		"""
		
		# Validate the rate
		if int(degPerSec) not in (2, 3, 4, 5, 6, 7, 8):
			raise ValueError("Invalid slew rate: %i deg/s" % int(degPerSec))
			
		self.port.write('#:Sw%i#' % int(degPerSec))
		return bool( self._readNumber() )
		
	def setSlewRateMax(self):
		"""
		Set the slew rate to "Slew" (fastest).
		"""
		
		self.port.write('#:RS#')
		
	def setSlewRateFind(self):
		"""
		Set the slew rate to "Find" (second fastest).
		"""
		
		self.port.write('#:RM#')
		
	def setSlewRateCenter(self):
		"""
		Set the slew rate to "Center" (second slowest).
		"""
		
		self.port.write('#:RC#')
		
	def setSlewRateGuide(self):
		"""
		Set the slew rate to "Guide" (slowest).
		"""
		
		self.port.write('#:RG#')
		
	def getCurrentPointing(self):
		"""
		Query the telescope for the current pointing location and return
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
		
	def getSiderealTime(self):
		"""
		Query the telescope for the current sidereal time and return it in
		decimal hours.
		"""
		
		self.port.write('#:GS#')
		lst = self._readString()
		
		lst = self._str2coord(lst, forceHigh=True)
		return lst
		
	def getElevationLimits(self):
		"""
		Get the current slew elevation limits, in degrees, and return them
		as a two-element tuple of lower limit, upper limit.
		"""
		
		self.port.write('#:Gh#')
		ll = self._readString()
		self.port.write('#:Go#')
		ul = self._readString()
		
		ll = ll.replace(chr(223), '')
		ll = int(ll)
		ul = ul.replace(chr(223), '')
		ul = int(ul)
		
		return ll, ul
		
	def setElevationLimits(self, lowerLimit, upperLimit):
		"""
		Set the telescope slew elevation limits.
		"""
		
		# Validate
		if lowerLimit < 0 or lowerLimit > 90:
			raise ValueError("Invalid lower limit: %i" % int(lowerLimit))
		if upperLimit < 0 or upperLimit > 90:
			raise ValueError("Invalid upper limit: %i" % int(upperLimit))
		if upperLimit <= lowerLimit:
			raise ValueError("Invalid limits: %i to %i" % (int(lowerLimit), int(upperLimit)))
			
		# Set
		status = True
		self.port.write("#:Sh%02i%s#" % (int(lowerLimit), chr(223)))
		status &= bool( self._readNumber() )
		self.port.write("#:So%02i%s#" % (int(upperLimit), chr(223)))
		status &= bool( self._readNumber() )
		
		return status
		
	def isSlewing(self):
		"""
		Return if the telescope is currently slewing to a target or not.
		Returns True if it is, False otherwise.
		"""
		
		self.port.write('#:D#')
		dist = self._readString()
		dist = sum([1 for c in dist if c == chr(255)])
		if dist > 2:
			return True
		else:
			return False
			
	def _moveToPosition(self, raStr, decStr, fast=False):
		"""
		Internal function to help with running moveToPosition() in the 
		background.
		"""
		
		# Set the RA/dec
		status = True
		self.port.write('#:Sr%s#' % raStr)
		if not fast:
			status &= bool( self._readNumber() )
		self.port.write('#:Sd%s#' % decStr)
		if not fast:
			status &= bool( self._readNumber() )
			
		# If that has worked, start the slew
		if status:
			self.port.write('#:MS#')
			if not fast:
				status &= not bool( self._readNumber() )
			if not status:
				msg = self._readString()
				
		# Flush the serial connection if we are in fast mode
		if fast:
			self.port.flushInput()
			self.port.flushOutput()
			
		return status
		
	def moveToPosition(self, ra, dec, fast=False, blocking=True):
		"""
		Given a RA value in decimal hours and a declination value in 
		decimal degrees, slew the telescope there provided that it is
		observable.  Returns True if the target was set the and the
		move was initiated, False otherwise.
		
		There are two keywords that control how this functions behaves:
		'fast' and 'blocking'.  Setting 'fast' to True disables the 
		on-line error detection in order to increase the command rate.
		'blocking' causes the function to poll the telescope for motion
		and waits until the current slew has finished before moving.
		"""
		
		# Convert to a string
		ra = self._ra2str(ra)
		dec = self._dec2str(dec)
		
		if blocking:
			# Wait until we are ready to move again
			while self.isSlewing():
				time.sleep(0.05)
				
			# Move
			self.haltCurrentSlew()
			status = self._moveToPosition(ra, dec, fast=fast)
		else:
			if self.isSlewing():
				status = False
			else:
				self.haltCurrentSlew()
				status = self._moveToPosition(ra, dec, fast=fast)
				
		return status
		
	def haltCurrentSlew(self):
		"""
		Stop the current slew.
		"""
		
		self.port.write('#:Q#')
		for d in ('n', 's'):
			self.port.write('#:Q%s#' % d)
			
	def resetTarget(self):
		"""
		Reset the target to clear a previous GoTo command.
		"""
		
		ra, dec = self.getCurrentPointing()
		ra, dec = self._ra2str(ra), self._dec2str(dec)
		self._moveToPosition(ra, dec, fast=False)


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
	
	# Signal loss information
	_nextSet = None
	los = 0.0
	
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
			
		## Calculate how long before we expect to lose the satellite
		if self.alt > 0:
			if self._nextSet is None:
				try:
					event = observer.next_pass(self)
					self._nextSet = event[4]
					if self._nextSet is None:
						self._nextSet = observer.date - 1.0
						
				except ValueError:
					self._nextSet = observer.date + 1.0
					
			self.los = (self._nextSet-observer.date)*86400.0
			
		else:
			self._nextSet = None
			self.los = 0.0
		
	def setStdMag(self, stdMag):
		"""
		Set the "standard magnitude" for this satellite to use for 
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
			setattr(self, attr, getattr(sat, attr, None)*_rad2deg)


def passPredictor(observer, satellites, date=None, time=None, utcOffset=0.0, duration=180.0, magnitudeCut=None):
	"""
	Function to predict passes of satellites and return them as a list of 
	tuples.  The entry in each tuple is:
	  * satellite name
	  * maximum brightness
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
		oU = int((oS-int(oS))*1e6)
		oS = int(oS)
		dt += oG*timedelta(seconds=oS, microseconds=oU)
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
		
		## Run until we are outside of the window to search (or until it looks
		## like we got stuck)
		j = 0
		while (observer.date-tStart) <= (duration/60.0/24.0) and j < 100:
			### Compute the next pass
			try:
				rTime, rAz, mTime, mEl, sTime, sAz = observer.next_pass(sat)
			except ValueError:
				#### If this satellite is geosynchronous, set some dummy values 
				#### using its current position
				sat.compute(observer)
				rTime, rAz = observer.date, sat.az
				mTime, mEl = observer.date+duration/2.0/60.0/24.0, sat.alt
				sTime, sAz = observer.date+duration/60.0/24.0, sat.az
				
			### Does the satellite actually rise?
			if rTime is None or mTime is None or sTime is None:
				observer.date += 2*duration/60.0/24.0
				continue
				
			### Does the rise time make sense?  This cases is needed to deal 
			### with a problem in PyEphem where the rise time is strange if 
			### satellite is already above the horizon at the current time
			if rTime > mTime and observer.date == tStart:
				observer.date -= 15/60.0/24.0
				rTime, rAz, junk1, junk2, junk3, junk4 = observer.next_pass(sat)
				
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
					
			if visible and maxBright <= magnitudeCut and (rTime-tStart) <= (duration/60.0/24.0) and mEl > 0:
				### Save this one and update the time to local
				rTimeLocal = ephem.date( rTime + utcOffset/24.0 )
				mTimeLocal = ephem.date( mTime + utcOffset/24.0 )
				sTimeLocal = ephem.date( sTime + utcOffset/24.0 )
				events.append( (sat.name, maxBright, rTimeLocal, rAz, mTimeLocal, mEl, sTimeLocal, sAz) )
				
			### Step forward in time to look for the next pass
			observer.date = sTime + 5/60.0/24.0
			j += 1
			
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
	
	def __init__(self, observer, satellites, telescope=None):
		"""
		Initialize the SatellitePositions instance using an ephem.Observer
		instance, a list of ephem.EarthSatellite instances, and, optionally:
		  * An update interval in seconds using the 'interval' keyword and
		  * An LX200-like telescope instance for tracking using the 
		    'telescope' keyword.
		"""
		
		# Copy the observer and convert the list of ephem.EarthSatellite 
		# instances to EarthSatellitePlus instances.  Also, save the
		# telescope information
		self.observer = observer
		self.satellites = []
		for sat in satellites:
			if type(sat) is not EarthSatellitePlus:
				newSat = EarthSatellitePlus()
				newSat.fillFromPyEphem(sat)
				sat = newSat
			self.satellites.append( sat )
		self.telescope = telescope
		
		# State variables for progressive updates to keep down the execution
		# time of update()
		self.nTiers = len(self.satellites)/50+1
		self.currentTier = 0
		self.tiers = [i % self.nTiers for i in xrange(len(self.satellites))]
		
		# State variables for fine-tuning the telescope tracking
		self.tracking = None
		self.timeOffset = timedelta()
		self.crossTrackOffset = 0.0
		
	def getObserver(self):
		"""
		Return a copy of the ephem.Observer instance being used for 
		computation.
		"""
		
		return copy.copy(self.observer)
		
	def setObserver(self, observer):
		"""
		Change the location of the observer using an ephem.Observer
		instance.  If the background tracking thread is running is
		also restarts the thread.
		"""
		
		self.observer = copy.copy(observer)
		
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
		
	def setCrossTrackOffset(self, degrees):
		"""
		Set the perpendicular track offset to a value in decimal degrees.
		"""
		
		self.crossTrackOffset = degrees*_deg2rad
		
	def adjustCrossTrackOffset(self, degrees):
		"""
		Apply an additional adjustment to the cross track offset using a 
		value in decimal degrees.
		"""
		
		self.crossTrackOffset += degrees*_deg2rad
		
	def getCrossTrackOffset(self):
		"""
		Return the perpendicular track offset in decimal degrees.
		"""
		
		return self.crossTrackOffset*_rad2deg
		
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
		if self.telescope is not None:
			self.telescope.haltCurrentSlew()
			
	def getTracking(self):
		"""
		Return the NORAD ID number of the satellite being currently 
		tracked.  If nothing is being tracked return None.
		"""
		
		return self.tracking
		
	def resetTracking(self):
		"""
		Reset the telescope target information to clear tracking problems.
		"""
		
		if self.telescope is not None:
			self.telescope.resetTarget()
			
	def getNumberVisible(self, magnitudeCut=None):
		"""
		Return the number of satellites currently visible.  This means 
		above the horizon and brighter than the specified magnitude cut.
		"""
		
		# Deal with the magnitude magnitude
		if magnitudeCut is None:
			magnitudeCut = 15.0
			
		nVis = sum( [1 for sat in self.satellites if sat.alt > 0 and sat.magnitude <= magnitudeCut] )
		return nVis
		
	def update(self, full=False):
		"""
		Computes the locations of all satellites in the current computation
		tier and provides pointing information to the telescope for the 
		satellite being tracked.  Returns the number of seconds used.
		
		.. note:: If a full upate is required, set the 'full' keyword to True.
		"""
		
		# Marker
		tStart = time.time()
		
		# The current time is further adjusted by the 'timeOffset' attribute
		# to help with tracking.
		tNow = datetime.utcnow() + self.timeOffset
		
		# Update the observer with the current time
		self.observer.date = tNow.strftime("%Y/%m/%d %H:%M:%S.%f")
		
		# Loop through the satellites and update them
		for tier,sat in zip(self.tiers, self.satellites):
			if not full:
				if tier != self.currentTier and sat.catalog_number != self.tracking:
					continue
					
			sat.compute(self.observer)
			
			if sat.alt > 0:
				## If the satellite is up, check and see if it is 
				## the one that we should be tracking.
				if sat.catalog_number == self.tracking:
					### Is there a telescope to use?
					if self.telescope is not None:
						#### Apply a perpendicular correction to the
						#### track to help with tracking
						ra, dec = getPointFromBearing(sat.ra, sat.dec, sat.bearing+math.pi/2, self.crossTrackOffset)
						
						#### Radians -> hours/degrees
						ra = ra*_rad2hr
						dec = dec*_rad2deg
						
						#### Command the telescope
						self.telescope.moveToPosition(ra, dec, blocking=False)
						
			else:
				## If it is no longer visible, check and see if it is
				## the satellite that we were tracking so that we can
				## stop.
				if sat.catalog_number == self.tracking:
					self.stopTracking()
					
		# Final time to figure out how much time we spent calculating 
		# positions.
		tStop = time.time()
		
		# Update the tier being computed
		self.currentTier = (self.currentTier + 1) % self.nTiers
		
		# Done
		return tStop-tStart


def main(args):
	# Parse the command line
	config = parseOptions(args)
	
	# Grab the list of filenames to look for TLEs in
	filenames = config['args']
	
	# Read in the TLEs and create a dictionary that helps map the NORAD ID
	# number to a particular entry.  This helps with assigned intrinsic 
	# magnitudes to satellites from the QuickSat qs.mag file.
	satellites = []
	lookup = {}
	for filename in filenames:
		fh = open(filename, 'r')
		data = fh.readlines()
		fh.close()
		
		## Is this geo.txt?  If so, filter out everything but the "bright" ones.  If we
		## don't do this we end up tracking far too many things at once and the telescope
		## command rate drops.
		if os.path.basename(filename) == 'geo.txt':
			toKeep = [26038, 26608, 26724, 28626, 28644, 28903, 29520, 32018, 39616]
		else:
			toKeep = None
			
		## Convert to EarthSatellitePlus instances
		for i in xrange(len(data)/3):
			sat = EarthSatellitePlus()
			try:
				sat.fillFromPyEphem( ephem.readtle(data[3*i+0], data[3*i+1], data[3*i+2]) )
			except ValueError:
				print "ERROR: Cannot parse TLE:"
				print "  %s" % data[3*i+0].rstrip()
				print "  %s" % data[3*i+1].rstrip()
				print "  %s" % data[3*i+2].rstrip()
				print "-> skipping"
				continue
				
			### Is this in the geo.txt "bright" list?
			try:
				if sat.catalog_number not in toKeep:
					continue
			except TypeError:
				pass
				
			### Have we already loaded this satellite from a different file?
			if sat.catalog_number not in lookup.keys():
				satellites.append( sat )
				lookup[sat.catalog_number] = len(satellites)-1
				
	# Load in the standard magnitude if we happen to have a qs.mag file to read
	# and add the information to the various
	# instances.
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
	except IOError as e:
		print "ERROR: Error reading 'qs.mag': %s" % str(e)
		pass
		
	# Report on what we found
	nSats = len(satellites)
	if nSats == 0:
		print "ERROR: Found no valid TLEs, exiting"
		sys.exit(1)
		
	print "Loaded %i TLEs from %i file%s" % (nSats, len(filenames), 's' if len(filenames) > 1 else '')
	
	# Try and setup the telescope.  This not only opens up the port but also
	# makes sure that the time is right.
	try:
		lx200 = LX200(config['port'], baud=config['baud'])
		
		## Make sure that the time is "reasonable"
		tTime = lx200.getDateTime()
		cTime = datetime.utcnow()
		tcOffset = cTime - tTime
		
		if abs(tcOffset) > timedelta(seconds=10):
			lx200.setDateTime(cTime)
			
			tTime = lx200.getDateTime()
			cTime = datetime.utcnow()
			tcOffset = cTime - tTime
			
		print "LX200 set to %s, computer at %s" % (tTime, cTime)
		print "-> Difference is %s" % tcOffset
		lx200StatusString = "LX200 set to %s, computer at %s" % (tTime, cTime)
		
		# Set the slew rate to maximum
		lx200.setMaximumSlewRate(8)
		lx200.setSlewRateMax()
		
		# Set high precision coordinates
		lx200.setHighPrecision()
		
	except Exception as e:
		print "ERROR: %s" % str(e)
		lx200 = None
		lx200StatusString = 'Telescope not connected'
		
	# Set the observer according to the telescope, if found.  Otherwise, use
	# the default values.
	try:
		obs = lx200.getObserver(elevation=config['elev'])
	except Exception as e:
		print "ERROR: %s" % str(e)
		obs = ephem.Observer()
		obs.lat = config['lat']*_deg2rad
		obs.lon = config['lon']*_deg2rad
		obs.elevation = config['elev']
	print "Observer set to %.4f %s, %.4f %s @ %.1f m" % (abs(obs.lat)*_rad2deg, 'N' if obs.lat >= 0 else 'S', abs(obs.lon)*_rad2deg, 'E' if obs.lon >= 0 else 'W', obs.elevation)
	
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
			oU = int((oS-int(oS))*1e6)
			oS = int(oS)
			dt += oG*timedelta(seconds=oS, microseconds=oU)
			dt = dt.strftime("%Y/%m/%d %H:%M:%S")
		else:
			## Use the provided date and time
			dt = "%s %s" % (d, t)
			
		# Compute the passes - this should return values in local time
		events = passPredictor(obs, satellites, date=d, time=t, utcOffset=config['utcOffset'],
					duration=config['predictorInterval'], magnitudeCut=config['predictorMagLimit'])
					
		print " "
		print "Satellite Pass Predictions for %s:" % dt
		print "%24s  %5s  %13s  %15s  %13s" % ('Name', 'Mag.', 'Rise Time', 'Max. Time/El.', 'Set Time')
		print "-" * (24+2+5+2+13+2+15+2+13)
		for event in events:
			mag = event[1]
			if mag < 15:
				mag = "%5.1f" % mag
			else:
				mag = " --- "
			rTime = str(event[2])[5:]
			mTime = str(event[4]).split(None, 1)[1]
			mEl = event[5]*_rad2deg
			sTime = str(event[6])[5:]
			print "%24s  %5s  %13s  %8s @ %4.1f  %13s" % (event[0], mag, rTime, mTime, mEl, sTime)
			
	else:
		# Tracking mode
		
		# Convert the provided UT offset to a timedelta instance that we 
		# can use to update the clock
		oG = -1 if config['utcOffset'] < 0 else 1
		oS = abs(config['utcOffset']*3600.0)
		oU = int((oS-int(oS))*1e6)
		oS = int(oS)
		utcOffset = oG*timedelta(seconds=oS, microseconds=oU)
		
		# Initialize the tracker and run the first full update
		trkr = SatellitePositionTracker(obs, satellites, telescope=lx200)
		trkr.update(full=True)
		
		# Setup the TUI
		stdscr = curses.initscr()
		curses.start_color()
		curses.noecho()
		curses.cbreak()
		stdscr.keypad(1)
		stdscr.nodelay(1)
		size = stdscr.getmaxyx()
		if size[0] < 24 or size[1] < 80:
			curses.nocbreak()
			stdscr.keypad(0)
			curses.echo()
			curses.endwin()
			raise RuntimeError("Terminal size needs to be at least 24 by 80")
		
		if curses.has_colors():
			curses.init_pair(1, curses.COLOR_WHITE,  curses.COLOR_BLACK)
			curses.init_pair(2, curses.COLOR_RED,    curses.COLOR_BLACK)
			curses.init_pair(3, curses.COLOR_BLACK,  curses.COLOR_GREEN)
			curses.init_pair(4, curses.COLOR_CYAN,   curses.COLOR_BLACK)
			
			std = curses.color_pair(1)
			sel = curses.color_pair(2) | curses.A_BOLD
			inf = curses.color_pair(3)
			ntf = curses.color_pair(4)
			
		else:
			std = curses.A_NORMAL
			sel = curses.A_NORMAL
			inf = curses.A_NORMAL
			ntf = curses.A_NORMAL
			
		# Main TUI loop
		## State control variables and such
		act = 0					# Line number of the satellite the selector is on
		trk = -1					# Line number of the satellite being tracked
		info = -1					# Whether or not to get info about a satellite
		magLimit = 6.0				# Magnitude limit for what to display
		msg = lx200StatusString		# Message string
		oldMsg = ''				# Message string - previous value (to help with clearing old messages)
		msgCount = 0				# Message string display counter (to help with clearing old messages)
		empty = ''				# Empty string (to help with clearing old messages/the window)
		while len(empty) < (5+2+24+2+12+2+12+2+12+2+5):
			empty += ' '
			
		try:
			while True:
				## Current time
				tNow = datetime.utcnow()
				tNowLocal = tNow + utcOffset
			
				## Update the tracker
				tElapsed = trkr.update()
			
				## Figure out how many satellites are currently visible and brighter
				## than the current magnitude limit
				trkChange = False
				nVis = trkr.getNumberVisible(magnitudeCut=magLimit)
			
				## Interact with the user's key presses - one key at a time
				c = stdscr.getch()
				curses.flushinp()
				if c == ord('Q'):
					### 'q' to exit the main loop after stopping the movement
					if trkr.getTracking() is not None:
						trkr.stopTracking()
					if lx200 is not None:
						lx200.haltCurrentSlew()
					break
				elif c == curses.KEY_UP:
					### Up arrow to change what is selected
					act -= 1
					act = max([act, 0])
				elif c == curses.KEY_DOWN:
					### Down array to change what is selected
					act += 1
					act = min([act, nVis-1, 12])
				elif c == ord('a'):
					### 'a' to adjust the tracking time offset - negative
					trkr.adjustTimeOffset( -config['trackOffsetStep'] )
					off = trkr.getTimeOffset()
					msg = 'Time offset now %+.1f s' % (off.days*86400+off.seconds+off.microseconds/1e6)
				elif c == ord('A'):
					### 'A' to adjust the tracking time offset - negative x10
					trkr.adjustTimeOffset( -10*config['trackOffsetStep'] )
					off = trkr.getTimeOffset()
					msg = 'Time offset now %+.1f s' % (off.days*86400+off.seconds+off.microseconds/1e6)
				elif c == ord('s'):
					### 's' to adjust the tracking time offset - positive
					trkr.adjustTimeOffset( config['trackOffsetStep'] )
					off = trkr.getTimeOffset()
					msg = 'Time offset now %+.1f s' % (off.days*86400+off.seconds+off.microseconds/1e6)
				elif c == ord('S'):
					### 's' to adjust the tracking time offset - positive x10
					trkr.adjustTimeOffset( 10*config['trackOffsetStep'] )
					off = trkr.getTimeOffset()
					msg = 'Time offset now %+.1f s' % (off.days*86400+off.seconds+off.microseconds/1e6)
				elif c == ord('z'):
					### 'z' to adjust the perpendicular track offset - negative
					trkr.adjustCrossTrackOffset( -config['crossTrackOffsetStep'] )
					off = trkr.getCrossTrackOffset()
					msg = 'Cross track offset now %+.1f degrees' % off
				elif c == ord('Z'):
					### 'Z' to adjust the perpendicular track offset - negative x10
					trkr.adjustCrossTrackOffset( -10*config['crossTrackOffsetStep'] )
					off = trkr.getCrossTrackOffset()
					msg = 'Cross track offset now %+.1f degrees' % off
				elif c == ord('w'):
					### 'w' to adjust the perpendicular track offset - positive
					trkr.adjustCrossTrackOffset( config['crossTrackOffsetStep'] )
					off = trkr.getCrossTrackOffset()
					msg = 'Cross track offset now %+.1f degrees' % off
				elif c == ord('W'):
					### 'W' to adjust the perpendicular track offset - positive x10
					trkr.adjustCrossTrackOffset( 10*config['crossTrackOffsetStep'] )
					off = trkr.getCrossTrackOffset()
					msg = 'Cross track offset now %+.1f degrees' % off
				elif c == ord('p'):
					off1 = trkr.getTimeOffset()
					off2 = trkr.getCrossTrackOffset()
					msg = 'Time offset is %+.1f s; Cross track offset is %+.1f degrees' % (off1.days*86400+off1.seconds+off1.microseconds/1e6, off2)
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
				elif c == ord('r'):
					### 'r' to reset the telescope target
					trkr.resetTracking()
					msg = 'Resetting telescope target'
				elif c == ord('t'):
					### 't' to toggle telesope tracking on/off
					if trk == -1:
						trk = act
						trkChange = True
					else:
						trk = -1
						trkChange = True
				elif c == ord('i'):
					### 'i' to get information about the currently selected satellite
					info = act
					
				## Additional state check to make sure that we aren't still trying
				## to track a satellite that has set/been eclipsed
				if not trkChange:
					if trk != -1 and trkr.getTracking() is None:
						trk = -1
						trkChange = True
					
				## Header for the satellite data
				output = "%5s  %24s  %12s  %12s  %12s  %5s\n" % ('ID', 'Name', 'Azimuth', 'Elevation', 'Status', 'Mag.')
				stdscr.addstr(0, 0, output, inf)
			
				## Have we stopped a track?  If so, let the user know.
				if trkChange:
					if trk == -1:
						trkr.stopTracking()
						msg = 'Tracking stopped'
						trkChange = False
					
				## Satellite information
				k = 0
				for j in xrange(nSats):
					### Make the satellite easy-to-access
					sat = trkr.satellites[j]
				
					### Is it visible and bright enough?
					if sat.alt > 0 and sat.magnitude <= magLimit:
						#### Create the output line
						output = "%5i  %24s  %12s  %12s  %8s-%3s  %5s" % (sat.catalog_number, sat.name, sat.az, sat.alt, 'Eclipsed' if sat.eclipsed else 'In Sun', 'Asc' if sat.rising else 'Dsc', '%5.1f' % sat.magnitude if sat.magnitude < 15 else ' --- ')
					
						#### See if we need to enable/disable tracking
						if trkChange:
							if k == trk:
								if sat.catalog_number != trkr.getTracking():
									trkr.setTimeOffset( timedelta() )
									trkr.setCrossTrackOffset( 0.0 )
									trkr.startTracking(sat.catalog_number)
									msg = 'Now tracking \'%s\' (NORAD ID #%i)' % (sat.name, sat.catalog_number)
								
						#### See if we need to poll and print info
						if k == info:
							## Is this one active?
							isTracked = ''
							if sat.catalog_number == trkr.getTracking():
								isTracked = ' (tracking)'
							
							## Loss of signal information
							los = sat.los
							los = '%i:%02i:%04.1f' % (int(los/3600.0), int(los/60.0)%60, los%60)
						
							msg = '%s%s: LoS %s, range %.1f km, velocity %.1f km/s' % (sat.name, isTracked, los, sat.range/1e3, sat.range_velocity/1e3)
							info = -1
						
						#### Add the line to the screen, provided it will fit
						if k <= 12:
							## Standard flag for normal satellites
							displayFlags = std
							## Make the satellite currently be traced more obvious
							if trkr.satellites[j].catalog_number == trkr.getTracking():
								displayFlags = sel | curses.A_BOLD
							## Make the satellite currently selected reversed
							if k == act:
								displayFlags |= curses.A_REVERSE
							
							stdscr.addstr(k+1,  0, output, displayFlags)
							k += 1
						else:
							break
						
				## Pad out the message stored in 'msg' generated from interaction 
				## with the user/telescope
				while len(msg) < (5+2+24+2+12+2+12+2+12+2+5):
					msg += ' '
				
				## See if it is time to clear off the message
				if msg == oldMsg:
					msgCount += 1
					if msgCount == 201:
						msg = empty
						msgCount = 0
				oldMsg = msg
			
				## Pad out the current UTC time
				output = "UTC: "
				output += tNow.strftime("%Y/%m/%d %H:%M:%S")
				output += ".%01i" % (float(tNow.strftime("%f"))/1e5,)
				output += "           "
				output += "%5.3f s" % tElapsed
				output += "           "
				output += "LT: "
				output += tNowLocal.strftime("%Y/%m/%d %H:%M:%S")
				output += ".%01i" % (float(tNowLocal.strftime("%f"))/1e5,)
				while len(output) < (5+2+24+2+12+2+12+2+12+2+5):
					output += ' '
				
				## Final time/message/help information
				stdscr.addstr(k+1,  0, output, inf)
				stdscr.addstr(k+2,  0, msg, ntf)
				stdscr.addstr(k+3,  0, 'Keys')
				stdscr.clrtoeol()
				stdscr.addstr(k+4,  0, '  t   - Start/stop tracking of the currently selected satellite')
				stdscr.clrtoeol()
				stdscr.addstr(k+5,  0, '  r   - Reset failed telescope slew')
				stdscr.clrtoeol()
				stdscr.addstr(k+6,  0, '  a/s - Decrease/Increase track offset by %.1f second (x10 with shift)' % ( config['trackOffsetStep'].seconds+config['trackOffsetStep'].microseconds/1e6))
				stdscr.clrtoeol()
				stdscr.addstr(k+7,  0, '  z/w - Decrease/Increase cross track offset by %.1f degrees (x10 with shift)' % config['crossTrackOffsetStep'])
				stdscr.clrtoeol()
				stdscr.addstr(k+8,  0, '  k/l - Decrease/Increase the magnitude limit by 0.5 mag')
				stdscr.clrtoeol()
				stdscr.addstr(k+9,  0, '  p   - Print current tracking offsets')
				stdscr.clrtoeol()
				stdscr.addstr(k+10, 0, '  Q   - Exit')
				stdscr.clrtoeol()
			
				## Pad out the window to deal with satellite setting
				stdscr.clrtobot()
			
				## Refresh the display
				stdscr.refresh()
			
				## Sleep until the next iteration
				tSleep = config['updateInterval'] - tElapsed
				time.sleep( max([tSleep, 0.001]) )
			
		except KeyboardInterrupt:
			pass
			
		except Exception as error:
			exc_type, exc_value, exc_traceback = sys.exc_info()
			fileObject = StringIO.StringIO()
			traceback.print_tb(exc_traceback, file=fileObject)
			tbString = fileObject.getvalue()
			fileObject.close()
			
		# Close out the curses session
		curses.nocbreak()
		stdscr.keypad(0)
		curses.echo()
		curses.endwin()
		
		# Final reporting
		try:
			## Error
			print "%s: failed with %s at line %i" % (os.path.basename(__file__), str(error), traceback.tb_lineno(exc_traceback))
			for line in tbString.split('\n'):
				print line
		except NameError:
			pass


if __name__ == "__main__":
	main(sys.argv[1:])
	