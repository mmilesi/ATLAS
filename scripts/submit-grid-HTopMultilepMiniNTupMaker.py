#!/usr/bin/env python

import glob, os, sys, subprocess
import datetime

current_time = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

username   = "mmilesi"
prodtag    = "v14-03-MiniNTup"
treename   = "nominal"
configpath = "$ROOTCOREBIN/user_scripts/HTopMultilepAnalysis/jobOptions_HTopMultilepMiniNTupMaker.py"
outdir     = "output_grid_DxAOD-13TeV_" + str(current_time)
destSE     = "AUSTRALIA-ATLAS_LOCALGROUPDISK"
exclSE     = "ANALY_IHEP,ANALY_JINR,ANALY_IHEP_GLEXEC,ANALY_RRC-KI-HPC,ANALY_RRC-KI-T1,IHEP_MCORE,IHEP_PROD,RRC-KI-HPC2,RRC-KI-T1,RRC-KI-T1_MCORE,RRC-KI-T1_TEST"

xAH_run = None

# ------------------------------------------------------------

# For MC
#
infilepath = "HTopMultilepAnalysis/doc/list-grid-HTopGroupNTup-MC.txt"
gridDSname = "user." + username + ".HTopMultilep." + prodtag + ".%in:name[3]%.%in:name[10]%"

xAH_run = "xAH_run.py -vv --files {0} --config {1} --treeName {2} --inputList --inputDQ2 --submitDir {3} prun --optGridMergeOutput=1 --optGridNFilesPerJob=1.0 --optGridDestSE={4} --optGridOutputSampleName={5} --optGridExcludedSite={6}".format(infilepath, configpath, treename, outdir, destSE, gridDSname, exclSE)

print("Executing command:\n{0}".format(xAH_run))
subprocess.call(xAH_run,shell=True)

# ------------------------------------------------------------

# For DATA
#
infilepath = "HTopMultilepAnalysis/doc/list-grid-HTopGroupNTup-DATA.txt"
gridDSname = "user." + username + ".HTopMultilep." + prodtag + ".%in:name[3]%.%in:name[4]%"

xAH_run = "xAH_run.py -vv --files {0} --config {1} --treeName {2} --inputList --inputDQ2 --submitDir {3} prun --optGridMergeOutput=1 --optGridNFilesPerJob=1.0 --optGridDestSE={4} --optGridOutputSampleName={5} --optGridExcludedSite={6}".format(infilepath, configpath, treename, outdir, destSE, gridDSname, exclSE)

print("Executing command:\n{0}".format(xAH_run))
subprocess.call(xAH_run,shell=True)

