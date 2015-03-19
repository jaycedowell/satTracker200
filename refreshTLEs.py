#!/usr/bin/env python

import os
import sys
import time
import getopt
import urllib


def usage(exitCode=None):
	print """refreshTLEs.py - Look at the CelesTrak TLEs in the current directory
and refresh old ones.

Usage: refreshTLEs.py [OPTIONS]

Options:
-h, --help           Display this help message
-f, --force          Force re-downloading the TLEs
-a, --age            Age limit in days for refreshing a file 
                     (default = 4 days)
"""
	
	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	config['force'] = False
	config['age'] = 4*86400
	
	try:
		opts, args = getopt.getopt(args, "hfa:", ["help", "force", "age="])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
	
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		elif opt in ('-f', '--force'):
			config['force'] = True
		elif opt in ('-a', '--age'):
			config['age'] = float(value)*86400
		else:
			assert False
			
	# Add in arguments
	config['args'] = args
	
	# Return configuration
	return config
	


def main(args):
	# Parse the command line
	config = parseOptions(args)
	
	# Figure out what do to
	toRefresh = []
	for filename in ('visual.txt', 'science.txt'):
		if os.path.exists(filename):
			## Get the age of the file
			mtime = os.stat(filename)[7]
			age = time.time() - mtime
			
			print "File '%s' last modified %.1f days ago" % (filename, age/86400.)
			
			## Does this file need to be refreshed?
			if age > config['age'] or config['force']:
				toRefresh.append( filename )
				print "-> adding to update list"
			else:
				print "-> skipping"
				
		else:
			toRefresh.append( filename )
			print "File '%s' does not exist" % filename
			print "-> adding to update list"
	print " "
	
	# Do it
	for filename in toRefresh:
		## Download the file and extract its contents
		print "Downloading '%s' from 'http://celestrak.com'" % filename
		
		url = "http://celestrak.com/NORAD/elements/%s" % filename
		
		t0 = time.time()
		urllib.urlretrieve(url, filename)
		t1 = time.time()
		sz = os.path.getsize(filename)
		print "-> downloaded %i bytes in %.3f s (%.1f kB/s)" % (sz, t1-t0, sz/1024./(t1-t0))


if __name__ == "__main__":
	main(sys.argv[1:])
	