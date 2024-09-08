#!/usr/bin/env python3

import os
import sys
import time
import zipfile
import argparse
from urllib.request import urlretrieve


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
        print(f"-> downloaded {sz} bytes in {t1-t0:.3f} s ({sz/1024/(t1-t0):.1f} kB/s)")
        
        print("Extracting to 'qs.mag'")
        qsmag = zipfile.ZipFile('qsmag.zip')
        
        with open('qs.mag', 'wb') as fh:
            fh.write( qsmag.read('qs.mag') )
        sz = os.path.getsize('qs.mag')
        print(f"-> extracted {sz} bytes to 'qs.mag'")
        
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
