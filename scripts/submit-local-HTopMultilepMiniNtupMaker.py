#!/usr/bin/env python

import glob, os, sys, subprocess

samplescsv = os.path.abspath(os.path.curdir) + "/HTopMultilepAnalysis/PlotUtils/Files/samples2015_HTopMultilep_25ns.csv"

sys.path.append(os.path.abspath(os.path.curdir)+"/HTopMultilepAnalysis/PlotUtils/")
from Core import NTupleTools, DatasetManager, listifyInputFiles

datasets = DatasetManager.DatasetManager()
sampledict = datasets.getListSamples(samplescsv,genericPath=True)

# -----------------------------
# 25ns_v7
# -----------------------------

infilelist = [
"HTopMultilepAnalysis/doc/list-local-HTopGroupNTup.txt", # for data
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361063/361063.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361064/361064.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361065/361065.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361066/361066.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361067/361067.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361068/361068.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361069/361069.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361070/361070.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361071/361071.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361072/361072.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361073/361073.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361074/361074.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361077/361077.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361078/361078.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361079/361079.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361080/361080.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361081/361081.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361082/361082.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361083/361083.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361084/361084.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361085/361085.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361086/361086.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361087/361087.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/410000/410000.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/410066/410066.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/410067/410067.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/410068/410068.root"
]

# -------------------------------------------------------------------------------------------------------

configpath = "$ROOTCOREBIN/user_scripts/HTopMultilepAnalysis/jobOptions_HTopMultilepMiniNtupMaker.py"
treename   = "nominal"
nevents    = 0

motherdir = os.path.abspath(os.path.curdir) + "/test"
if not os.path.exists(motherdir):
    os.makedirs(motherdir)
for s in sampledict:
    groupdir = motherdir + "/" + s["group"]
    if not os.path.exists(groupdir):
        os.makedirs(groupdir)

for infile in infilelist:

  xAH_run = None

  # In case of DATA, read the list of infiles from a txt file and execute one single job
  #
  if infile[-4:] == ".txt":
     outdir = "Data"
     xAH_run = "xAH_run.py -vv --files {0} --config {1} --inputList --treeName {2} --submitDir {3} --nevents {4} --force direct".format(infile,configpath,treename,outdir,nevents)
  else :
     outdir = ( infile.split("/") )[-1]
     outdir = outdir.replace(".root","")
     xAH_run = "xAH_run.py -vv --files {0} --config {1} --treeName {2} --submitDir {3} --nevents {4} --force direct".format(infile,configpath,treename,outdir,nevents)

  print("Executing command:\n{0}".format(xAH_run))
  subprocess.call(xAH_run,shell=True)

  # Move output file(s) from job directory to the proper one,
  # but first, change the file name to be readable in KG's FW!
  #
  for s in sampledict:

     if outdir == s["ID"] or ( outdir == "Data" and not s["ID"] ):
	outputfilepath = outdir +"/data-output"

        separator = "."
        if not s["ID"]:
          separator = ""

	INTREE  = outputfilepath + "/" + outdir + ".root"
	INHIST  =  outdir + "/hist-" + outdir + ".root"
	if ( outdir == "Data" ):
	   INTREE  = outputfilepath + "/list-local-HTopGroupNTup.root"
	   INHIST  =  outdir + "/hist-list-local-HTopGroupNTup.root"

	OUTTREE = motherdir + "/" + s["group"] + "/" + s["ID"] + separator + s["name"] + ".root"
	OUTHIST = motherdir + "/" + s["group"] + "/hist-" + s["ID"] + separator + s["name"] + ".root"

	print("Moving :\n{0}\nto:\n{1}".format(INTREE,OUTTREE))
	print("Moving :\n{0}\nto:\n{1}".format(INHIST,OUTHIST))

	os.rename(INTREE,OUTTREE)
	os.rename(INHIST,OUTHIST)
