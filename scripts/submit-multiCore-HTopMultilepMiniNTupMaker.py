#!/usr/bin/env python

import glob, os, sys, subprocess, shutil, argparse

import multiprocessing

parser = argparse.ArgumentParser(description='Run HTopMultilepMiniNTupMaker algorithm interactively on multiple cores')

parser.add_argument('destination', metavar='destination', type=str,
                    help='The base directory where the output will be stored. Subdirectories for different sample groups will be created automatically by the job.')
parser.add_argument("--prod_ID", dest="prod_ID", action="store", default="25ns_v29/01", type=str,
                    help="The NTup production tag, e.g. 25ns_v19, 25ns_v24/02, ...  (default: prod_ID=25ns_v29/01)")
parser.add_argument('--nparallel', dest='nparallel', action='store', default=1, type=int,
                    help='The maximum number of parallel processes to be executed. The number to choose depends on the features of the machine you are using. Default is 1 (i.e., no parallel jobs)')
parser.add_argument('--nevents', dest='nevents', action='store', default=0, type=int,
                    help='The number of events to be processed. Default is 0 (i.e., ALL events)')
parser.add_argument('--skip', dest='skip', action='store', default=0, type=int,
                    help='Skip all entries before the one selected. Use in combination w/ --nevents=1 if you wish to debug a single event. Default is 0 (i.e., do NOT skip any entry)')
parser.add_argument('--treename', dest='treename', action='store', default="nominal", type=str,
                    help='The name of the input TTree. Default is \"nominal\"')

args = parser.parse_args()

samplescsv = os.path.abspath(os.path.curdir) + "/HTopMultilepAnalysis/PlotUtils/Files/samples_HTopMultilep_Priority1.csv"

sys.path.append(os.path.abspath(os.path.curdir)+"/HTopMultilepAnalysis/PlotUtils/")
from Core import NTupleTools, DatasetManager, listifyInputFiles

datasets = DatasetManager.DatasetManager()
sampledict = datasets.getListSamples(samplescsv,genericPath=True)

import normaliseTrees

def generateCmdList(samples):

    cmdlist = []
    for sample in samples:

        xAH_run = None
        infile  = None

        # In case of DATA, read the list of infiles from a txt file and execute one single job

        if sample[-4:] == ".txt":
            outdir = "Data"
	    infile = sample
    	    xAH_run = "xAH_run.py -vv --files {0} --config {1} --inputList --treeName {2} --submitDir {3} --nevents {4} --skip {5} --force direct".format(infile,configpath,treename,outdir,nevents,skip)
        else :
	    outdir = sample
	    infile = sample_path + "/" + sample + "/" + sample + ".root"
    	    xAH_run = "xAH_run.py -vv --files {0} --config {1} --treeName {2} --submitDir {3} --nevents {4} --skip {5} --force direct".format(infile,configpath,treename,outdir,nevents,skip)

	cmdlist.append(xAH_run)

    return cmdlist

def listchunks(l, n):
    n = max(1, n)
    return [l[i:i + n] for i in range(0, len(l), n)]


def miniNTuplise(sample):

    print("Executing command:\n{0}".format(sample))

    subprocess.call(sample,shell=True)

    # Move output file(s) from job directory to the proper one,
    # but first, change the file name to be readable in KG's FW!

    knownDSID = False

    # Get output directory name from the executed command string

    outdir = sample[sample.find("--submitDir")+len("--submitDir")+1:sample.find("--nevents")-1]

    for s in sampledict:

     if outdir == s["ID"] or ( outdir == "Data" and not s["ID"] ):

        knownDSID = True
	outputfilepath = outdir +"/data-output"

        separator = "."
        if not s["ID"]:
          separator = ""

	INTREE  = outputfilepath + "/" + outdir + ".root"
	INHIST  = outdir + "/hist-" + outdir + ".root"
	if ( outdir == "Data" ):
	   INTREE  = outputfilepath + "/list-local-HTopGroupNTup.root"
	   INHIST  = outdir + "/hist-list-local-HTopGroupNTup.root"

	OUTTREE = motherdir + "/" + s["group"] + "/" + s["ID"] + separator + s["name"] + ".root"
	OUTHIST = motherdir + "/" + s["group"] + "/hist-" + s["ID"] + separator + s["name"] + ".root"

	print("Moving :\n{0}\nto:\n{1}".format(INTREE,OUTTREE))
	print("Moving :\n{0}\nto:\n{1}".format(INHIST,OUTHIST))

	shutil.move(INTREE,OUTTREE)
	shutil.move(INHIST,OUTHIST)

	normaliseTrees.applyWeight(OUTTREE,s,isdata=bool(outdir == "Data" and not s["ID"]))

        break

    if not knownDSID:
        print("Simply removing {0} b/c corresponding DSID is unknown...".format(outdir))
    shutil.rmtree(os.path.abspath(os.path.curdir) + "/" + outdir)


if __name__ == '__main__':

    #sample_path = "/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/" + args.prod_ID + "/Nominal"
    #sample_path = "/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/" + args.prod_ID + "/ttbar_ttgamma/Nominal"
    sample_path = "/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/" + args.prod_ID + "/Nominal_PLICFT"
    #sample_path =  "/afs/cern.ch/user/m/mmilesi/work/public/ttH/GroupNTup/25ns_" + args.prod_ID + "/Nominal"

    infilelist = [
## Modify this file to include ROOT files w/ data runs
#
#"HTopMultilepAnalysis/doc/list-local-HTopGroupNTup.txt",
#
## MC: just list DSIDs
#
#"410000",
#"410155",
"410501",
    ]

    # -------------------------------------------------------------------------------------------------------

    configpath = "$ROOTCOREBIN/user_scripts/HTopMultilepAnalysis/jobOptions_HTopMultilepMiniNTupMaker.py"
    treename   = args.treename
    nevents    = args.nevents
    skip       = args.skip
    motherdir  = args.destination

    if not os.path.exists(motherdir):
	os.makedirs(motherdir)
    for s in sampledict:
	groupdir = motherdir + "/" + s["group"]
	if not os.path.exists(groupdir):
	    os.makedirs(groupdir)

    list_commands = generateCmdList(infilelist)

    MAX_PARALLEL = args.nparallel

    print listchunks(list_commands,MAX_PARALLEL)

    for chunk in listchunks(list_commands,MAX_PARALLEL):

	print("Processing samples: ")
	print("\n".join("{0} - {1}".format(idx,elem[elem.find("--submitDir")+len("--submitDir")+1:elem.find("--nevents")-1]) for idx, elem in enumerate(chunk)))
	p = multiprocessing.Pool(MAX_PARALLEL)
	p.map(miniNTuplise,chunk)
	p.close()
	p.join()
