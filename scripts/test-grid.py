#!/usr/bin/env python

import glob, os, sys, subprocess
import datetime

current_time = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

username   = "mmilesi"
prodtag    = "v7-MiniNTup-02"
infilepath = "user.dhohn.410000.PowhegPythiaEvtGen.DAOD_HIGG8D1.e3698_s2608_s2183_r7267_r6282_p2559.29.03.16.Nominal_output.root/"
configpath = "$ROOTCOREBIN/user_scripts/HTopMultilepAnalysis/jobOptions_HTopMultilepMiniNtupMaker.py"
outdir     = "output_grid_DxAOD-2015-13TeV_" + str(current_time)
destSE     = "AUSTRALIA-ATLAS_LOCALGROUPDISK"
exclSE     = "ANALY_IHEP,ANALY_JINR,ANALY_IHEP_GLEXEC,ANALY_RRC-KI-HPC,ANALY_RRC-KI-T1,IHEP_MCORE,IHEP_PROD,RRC-KI-HPC2,RRC-KI-T1,RRC-KI-T1_MCORE,RRC-KI-T1_TEST"

gridDSname = "user." + username + ".HTopMultilep." + prodtag + ".%in:name[3]%.%in:name[10]%"

xAH_run = "xAH_run.py -vv --files {0} --config {1} --inputDQ2 --submitDir {2} prun --optGridMergeOutput=1 --optGridNFilesPerJob=1.0 --optGridDestSE={3} --optGridOutputSampleName={4} --optGridExcludedSite={5}".format(infilepath, configpath, outdir, destSE, gridDSname, exclSE)

print("Executing command:\n{0}".format(xAH_run))
subprocess.call(xAH_run,shell=True)
