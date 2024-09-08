#!/usr/bin/env python3

import os
import sys
import time
import argparse
from urllib.request import urlretrieve


def main(args):
    # Figure out what do to
    toRefresh = []
    for filename in ('visual.txt', 'science.txt', 'resource.txt', 'geo.txt', 'starlink.txt'):
        if os.path.exists(filename):
            ## Get the age of the file
            mtime = os.stat(filename)[8]
            age = time.time() - mtime
            
            print("File '%s' last modified %.1f days ago" % (filename, age/86400.))
            
            ## Does this file need to be refreshed?
            if age > args.age or args.force:
                toRefresh.append( filename )
                print("-> adding to update list")
            else:
                print("-> skipping")
                
        else:
            toRefresh.append( filename )
            print("File '%s' does not exist" % filename)
            print("-> adding to update list")
    print(" ")
    
    # Do it
    for filename in toRefresh:
        ## Download the file and extract its contents
        print("Downloading '%s' from 'http://celestrak.com'" % filename)
        
        url = "http://celestrak.com/NORAD/elements/%s" % filename
        
        t0 = time.time()
        urlretrieve(url, filename)
        t1 = time.time()
        sz = os.path.getsize(filename)
        print("-> downloaded %i bytes in %.3f s (%.1f kB/s)" % (sz, t1-t0, sz/1024./(t1-t0)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='look at the CelesTrak TLEs in the current directory and refresh old ones', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('-f', '--force', action='store_true',
                        help='force re-downloading the TLEs')
    parser.add_argument('-a', '--age', type=int, default=2,
                        help='age limit in days for refreshing a file')
    args = parser.parse_args()
    args.age *= 86400    # days -> seconds
    main(args)
