#!/usr/bin/env python

import glob, os, sys, subprocess, shutil, argparse

import multiprocessing

parser = argparse.ArgumentParser(description='Run HTopMultilepMiniNTupMaker algorithm interactively on multiple cores')

parser.add_argument('destination', metavar='destination', type=str,
                    help='The base directory where the output will be stored. Subdirectories for different sample groups will be created automatically by the job.')
parser.add_argument('version', metavar='version', type=str,
                    help='The NTuple version to be used (e.g., v17, v18 ...)')
parser.add_argument('--nparallel', dest='nparallel', action='store', default=1, type=int,
                    help='The maximum number of parallel processes to be executed. The number to choose depends on the features of the machine you are using. Default is 1 (i.e., no parallel jobs)')
parser.add_argument('--nevents', dest='nevents', action='store', default=0, type=int,
                    help='The number of events to be processed. Default is 0 (i.e., ALL events)')
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
        #
        if sample[-4:] == ".txt":
            outdir = "Data"
	    infile = sample
    	    xAH_run = "xAH_run.py -vv --files {0} --config {1} --inputList --treeName {2} --submitDir {3} --nevents {4} --force direct".format(infile,configpath,treename,outdir,nevents)
        else :
	    outdir = sample
	    infile = sample_path + "/" + sample + "/" + sample + ".root"
    	    xAH_run = "xAH_run.py -vv --files {0} --config {1} --treeName {2} --submitDir {3} --nevents {4} --force direct".format(infile,configpath,treename,outdir,nevents)

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
    #
    knownDSID = False

    # Get output directory name from the executed command string
    #
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

    sample_path = "/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_" + args.version + "/Nominal"

    infilelist = [
#"HTopMultilepAnalysis/doc/list-local-HTopGroupNTup.txt",
#"341989",
#"341992",
#"341995",
#"341998",
#"342001",
#"342004",
#"342284",
#"342285",
#"343365",
#"343366",
#"343367",
#"361063",
#"361064",
#"361065",
#"361066",
#"361067",
#"361068",
#"361069",
#"361070",
#"361071",
#"361072",
#"361073",
#"361077",
#"361091",
#"361092",
#"361093",
#"361094",
#"361095",
#"361096",
#"361097",
#"361372",
#"361373",
#"361374",
#"361375",
#"361376",
#"361377",
#"361378",
#"361379",
#"361380",
#"361381",
#"361382",
#"361383",
#"361384",
#"361385",
#"361386",
#"361387",
#"361388",
#"361389",
#"361390",
#"361391",
#"361392",
#"361393",
#"361394",
#"361395",
#"361396",
#"361397",
#"361398",
#"361399",
#"361400",
#"361401",
#"361402",
#"361403",
#"361404",
#"361405",
#"361406",
#"361407",
#"361408",
#"361409",
#"361410",
#"361411",
#"361412",
#"361413",
#"361414",
#"361415",
#"361416",
#"361417",
#"361418",
#"361419",
#"361420",
#"361421",
#"361422",
#"361423",
#"361424",
#"361425",
#"361426",
#"361427",
#"361428",
#"361429",
#"361430",
#"361431",
#"361432",
#"361433",
#"361434",
#"361435",
#"361436",
#"361437",
#"361438",
#"361439",
#"361440",
#"361441",
#"361442",
#"361443",
#"361468",
#"361469",
#"361470",
#"361471",
#"361472",
#"361473",
#"361474",
#"361475",
#"361476",
#"361477",
#"361478",
#"361479",
#"361480",
#"361481",
#"361482",
#"361483",
#"361484",
#"361485",
#"361486",
#"361487",
#"361488",
#"361489",
#"361490",
#"361491",
#"361620",
#"361621",
#"361622",
#"361623",
#"361624",
#"361625",
#"361626",
#"361627",
"410000",
#"410011",
#"410012",
#"410013",
#"410014",
#"410025",
#"410026",
#"410049",
#"410073",
#"410074",
#"410075",
#"410155",
#"410156",
#"410157",
#"410215",
#"410218",
#"410219",
#"410220",
    ]

    # -------------------------------------------------------------------------------------------------------

    configpath = "$ROOTCOREBIN/user_scripts/HTopMultilepAnalysis/jobOptions_HTopMultilepMiniNTupMaker.py"
    treename   = args.treename
    nevents    = args.nevents
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
