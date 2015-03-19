#!/usr/bin/env python

import os
import sys
import time
import getopt
import urllib
import zipfile


def usage(exitCode=None):
	print """downloadQSMag.py - Download and extract the QuickSat intrinsic
magntiude file

Usage: downloadQSMag.py [OPTIONS]

Options:
-h, --help           Display this help message
-f, --force          Force re-downloading and extracting the data
"""
	
	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	config['force'] = False
	
	try:
		opts, args = getopt.getopt(args, "hf", ["help", "force"])
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
	if os.path.exists('qs.mag') and not config['force']:
		## Nothing to do
		print "The 'qs.mag' file already exists, skipping."
		
	else:
		## Download the file and extract its contents
		print "Downloading 'qsmag.zip'"
		t0 = time.time()
		urllib.urlretrieve("https://www.prismnet.com/~mmccants/programs/qsmag.zip", "qsmag.zip")
		t1 = time.time()
		sz = os.path.getsize("qsmag.zip")
		print "-> downloaded %i bytes in %.3f s (%.1f kB/s)" % (sz, t1-t0, sz/1024./(t1-t0))
		
		print "Extracting to 'qs.mag'"
		qsmag = zipfile.ZipFile('qsmag.zip')
		
		fh = open('qs.mag', 'w')
		fh.write( qsmag.read('qs.mag') )
		fh.close()
		sz = os.path.getsize('qs.mag')
		print "-> extracted %s bytes to 'qs.mag'" % sz
		
		## Cleanup
		os.unlink('qsmag.zip')


if __name__ == "__main__":
	main(sys.argv[1:])
	