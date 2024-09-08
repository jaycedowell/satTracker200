#!/usr/bin/env python3

import os
import sys
import time
import zipfile
import argparse
import urllib.request import urlretrieve


def main(args):
    # Figure out what do to
    if os.path.exists('qs.mag') and not args.force:
        ## Nothing to do
        print("The 'qs.mag' file already exists, skipping.")
        
    else:
        ## Download the file and extract its contents
        print("Downloading 'qsmag.zip'")
        t0 = time.time()
        urlretrieve("https://www.prismnet.com/~mmccants/programs/qsmag.zip", "qsmag.zip")
        t1 = time.time()
        sz = os.path.getsize("qsmag.zip")
        print("-> downloaded %i bytes in %.3f s (%.1f kB/s)" % (sz, t1-t0, sz/1024./(t1-t0)))
        
        print("Extracting to 'qs.mag'")
        qsmag = zipfile.ZipFile('qsmag.zip')
        
        fh = open('qs.mag', 'w')
        fh.write( qsmag.read('qs.mag') )
        fh.close()
        sz = os.path.getsize('qs.mag')
        print("-> extracted %s bytes to 'qs.mag'" % sz)
        
        ## Cleanup
        os.unlink('qsmag.zip')


if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description='download and extract the QuickSat intrinsic magntiude file',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
	parser.add_argument('-f', '--force', action='store_true',
						help='force re-downloading and extracting the data')
	args = parser.parse_args()
	main(args)
