#!/usr/bin/env python

import glob, os, sys, subprocess, shutil, argparse

samplescsv = os.path.abspath(os.path.curdir) + "/HTopMultilepAnalysis/PlotUtils/Files/samples2015_HTopMultilep_25ns.csv"

sys.path.append(os.path.abspath(os.path.curdir)+"/HTopMultilepAnalysis/PlotUtils/")
from Core import NTupleTools, DatasetManager, listifyInputFiles

datasets = DatasetManager.DatasetManager()
sampledict = datasets.getListSamples(samplescsv,genericPath=True)

import normaliseTrees

parser = argparse.ArgumentParser(description='Weight ROOT trees created by PBS job w/ Xsec weight')

parser.add_argument('dest', metavar='dest',type=str,
                   help='Base directory where to store output for the process in question')
parser.add_argument('DSID', metavar='DSID',type=str,
                   help='Dataset ID of the process in question')

args = parser.parse_args()

if __name__ == '__main__':

    if not os.path.exists(args.dest):
        os.makedirs(args.dest)
    for s in sampledict:
        groupdir = args.dest + "/" + s["group"]
        if not os.path.exists(groupdir):
            os.makedirs(groupdir)

    for s in sampledict:

        if args.DSID == s["ID"] or ( args.DSID == "Data" and not s["ID"] ):

            knownDSID = True

            separator = "."
            if not s["ID"]:
                separator = ""

            INTREE  = args.DSID + "/data-output" + "/" + args.DSID + ".root"
            INHIST  = args.DSID + "/hist-" + args.DSID + ".root"
            if ( args.DSID == "Data" ):
                INTREE  = args.DSID +"/data-output" + "/list-local-HTopGroupNTup.root"
                INHIST  = args.DSID + "/hist-list-local-HTopGroupNTup.root"

            OUTTREE = args.dest + "/" + s["group"] + "/" + s["ID"] + separator + s["name"] + ".root"
            OUTHIST = args.dest + "/" + s["group"] + "/hist-" + s["ID"] + separator + s["name"] + ".root"

            print("Moving :\n{0}\nto:\n{1}".format(INTREE,OUTTREE))
            print("Moving :\n{0}\nto:\n{1}".format(INHIST,OUTHIST))

            shutil.move(INTREE,OUTTREE)
            shutil.move(INHIST,OUTHIST)

            normaliseTrees.applyWeight(OUTTREE,s,isdata=bool(args.DSID == "Data" and not s["ID"]))

            break

    if not knownDSID:
        print("Simply removing {0} b/c corresponding DSID is unknown...".format(args.DSID))
    shutil.rmtree(os.path.abspath(os.path.curdir) + "/" + args.DSID)
