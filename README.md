# satTracker200
===============

## Description
Python-based satellite pass prediction/tracking software for the LX200 classic telescope.

## Requirements
 * Python >=2.7 and <3.0
 * PyEphem
 * PySerial
 
## Included Files
 * satTracker200.py - Main script that deals with satellite pass prediction and tracking
 * satTracker200.config.example - Example satTracker200.py configuration file
 * downloadQSMag.py - Utility to download and extract the QuickSat intrinsic magnitude file
 * refreshTLEs.py - Utility to help keep two-line element set (TLE) files up to date

## Setup
 1) Copy the 'satTracker200.config.example' file to 'satTracker200.config' and update it for your site.  You will want to change the latitude, longitude, elevation, utcOffset, and port parameters.
 
 2) Run 'downloadQSMag.py' to grab the satellite magntiude list.
 
 3) Run 'refreshTLEs.py' to download the latest TLEs from www.celestrak.com.

## Usage
### Tracking Mode (Default)

	python2.7 satTracker200.py visual.txt science.txt

### Prediction Mode

	python2.7 satTracker200.py -p visual.txt science.txt

For a complete list of options for the included scripts, see the integrated help.
