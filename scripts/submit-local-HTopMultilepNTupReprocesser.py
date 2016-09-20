#!/usr/bin/env python

import glob, os, sys, subprocess, shutil, argparse

parser = argparse.ArgumentParser(description='Run HTopMultilepNTupReprocesser algorithm interactively')

parser.add_argument('source', metavar='source', type=str,
                    help='The source (base) directory where the input is stored.')
parser.add_argument('destination', metavar='destination', type=str,
                    help='The base directory where the output will be stored. Subdirectories for different sample groups will be created automatically by the job.')
parser.add_argument('--nevents', dest='nevents', action='store', default=0, type=int,
                    help='The number of events to be processed. Default is 0 (i.e., ALL events)')
parser.add_argument('--treename', dest='treename', action='store', default="physics", type=str,
                    help='The name of the input TTree. Default is \"physics\"')
parser.add_argument('--closure', dest='closure', action='store_true', default=False,
                    help='Run on ttbar to perform MM closure test. Default is False, i.e., the code will run on data.')

args = parser.parse_args()

samplescsv = os.path.abspath(os.path.curdir) + "/HTopMultilepAnalysis/PlotUtils/Files/samples2015_HTopMultilep_25ns.csv"

sys.path.append(os.path.abspath(os.path.curdir)+"/HTopMultilepAnalysis/PlotUtils/")
from Core import NTupleTools, DatasetManager, listifyInputFiles

datasets = DatasetManager.DatasetManager()
sampledict = datasets.getListSamples(samplescsv,genericPath=True)

motherdir = args.source
destdir   = args.destination

infilelist = []

if not args.closure:
    infilelist.append(motherdir + "/Data/physics_Main.root")
else:
    infilelist.append(motherdir + "/tops/410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.root")
    #infilelist.append(motherdir + "/tops/410187.Sherpa_NNPDF30NNLO_ttbar_SingleLeptonP_MEPS_NLO.root")
    #infilelist.append(motherdir + "/tops/410188.Sherpa_NNPDF30NNLO_ttbar_SingleLeptonM_MEPS_NLO.root")

# -------------------------------------------------------------------------------------------------------

configpath = "$ROOTCOREBIN/user_scripts/HTopMultilepAnalysis/jobOptions_HTopMultilepNTupReprocesser.py"
treename   = args.treename
nevents    = args.nevents

for infile in infilelist:

  xAH_run = None

  outdir = ( infile.split("/") )[-1]
  outdir = ( outdir.split(".") )[-2]
  xAH_run = "xAH_run.py -vv --files {0} --config {1} --treeName {2} --submitDir {3} --nevents {4} --force direct".format(infile,configpath,treename,outdir,nevents)

  print("Executing command:\n{0}".format(xAH_run))
  subprocess.call(xAH_run,shell=True)

  # Move output file(s) from job directory to the proper one,
  # but first, change the file name to be readable in KG's FW!
  #
  for s in sampledict:

     if outdir == s["name"]:

	outputfilepath = outdir +"/data-output"

        separator = "."
        if not s["ID"]:
          separator = ""

	INTREE  = outputfilepath + "/" + s["group"] + ".root"

        if not os.path.exists(destdir + "/" + s["group"]):
            os.makedirs(destdir + "/" + s["group"])

	OUTTREE = destdir + "/" + s["group"] + "/" + s["ID"] + separator + s["name"] + ".root"

	print("Moving :\n{0}\nto:\n{1}".format(INTREE,OUTTREE))

	shutil.move(INTREE,OUTTREE)

        break

  shutil.rmtree(os.path.abspath(os.path.curdir) + "/" + outdir)
