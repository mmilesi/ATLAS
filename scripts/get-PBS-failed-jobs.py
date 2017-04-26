#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='Search PBS logfiles for failed jobs')
parser.add_argument('PBSdir', metavar='PBSdir',type=str, help='Path to directory containing PBS log files.')
args = parser.parse_args()

match = "Errno"
cmd = 'grep -nrHI "{0}" {1}'.format(match,args.PBSdir)

import subprocess

# This simply prints out lines where grep succeeded
subprocess.call(cmd,shell=True)

# This saves the grep match to a string for further usage
# child = subprocess.Popen('grep -nrHI "{0}" {1}'.format(match,args.PBSdir),stdout=subprocess.PIPE,shell=True)
# output = child.communicate()[0]
# print output: