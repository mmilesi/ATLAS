#!/usr/bin/env python

import glob, os, sys, subprocess, shutil

import multiprocessing

samplescsv = os.path.abspath(os.path.curdir) + "/HTopMultilepAnalysis/PlotUtils/Files/samples2015_HTopMultilep_25ns.csv"

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

    #sample_path = "/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v18/Nominal"
    sample_path = "/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v19/Nominal"
    #sample_path = "/afs/cern.ch/user/m/mmilesi/work/private/HTopMultileptonsTestSamples/25ns_v13/Nominal",

    infilelist = [
## data
#"HTopMultilepAnalysis/doc/list-local-HTopGroupNTup.txt",
## ttbar noallhad
"410000",
## ttH
"343365",
"343366",
"343367",
## ttW
"410155",
## ttZ
"410218",
"410219",
"410220",
## tZ noallhad
"410050",
## 4 tops
"410080",
## ttWW
"410081",
## VV
"361063",
"361064",
"361065",
"361066",
"361067",
"361068",
"361069",
"361070",
"361071",
"361072",
"361073",
"361077",
"361081",
"361082",
"361083",
"361084",
"361085",
"361086",
"361087",
## tHbj
"341989",
"341992",
"341995",
## WtH
"341998",
"342001",
"342004",
## VVV
"361620",
"361621",
"361622",
"361623",
"361624",
"361625",
"361626",
"361627",
## WtZ
"410215",
## single top
"410011",
"410012",
## tW
"410013",
"410014",
## Z+jets
"361372",
"361373",
"361374",
"361375",
"361376",
"361377",
"361378",
"361379",
"361380",
"361381",
"361382",
"361383",
"361384",
"361385",
"361386",
"361387",
"361388",
"361389",
"361390",
"361391",
"361392",
"361393",
"361394",
"361395",
"361396",
"361397",
"361398",
"361399",
"361400",
"361401",
"361402",
"361403",
"361404",
"361405",
"361406",
"361407",
"361408",
"361409",
"361410",
"361411",
"361412",
"361413",
"361414",
"361415",
"361416",
"361417",
"361418",
"361419",
"361420",
"361421",
"361423",
"361424",
"361425",
"361426",
"361427",
"361428",
"361429",
"361430",
"361431",
"361432",
"361433",
"361434",
"361435",
"361436",
"361437",
"361438",
"361439",
"361440",
"361441",
"361442",
"361443",
"361468",
"361469",
"361470",
"361471",
"361472",
"361473",
"361474",
"361475",
"361476",
"361477",
"361478",
"361479",
"361480",
"361481",
"361482",
"361483",
"361484",
"361485",
"361486",
"361487",
"361488",
"361489",
"361490",
"361491",
## W+jets
"361300",
"361301",
"361302",
"361303",
"361304",
"361305",
"361306",
"361307",
"361308",
"361309",
"361310",
"361311",
"361312",
"361313",
"361314",
"361315",
"361316",
"361317",
"361318",
"361319",
"361320",
"361321",
"361322",
"361323",
"361324",
"361325",
"361326",
"361327",
"361328",
"361329",
"361330",
"361331",
"361332",
"361333",
"361334",
"361335",
"361336",
"361337",
"361338",
"361339",
"361340",
"361341",
"361342",
"361343",
"361344",
"361345",
"361346",
"361347",
"361348",
"361349",
"361350",
"361351",
"361352",
"361353",
"361354",
"361355",
"361356",
"361357",
"361358",
"361359",
"361360",
"361361",
"361362",
"361363",
"361364",
"361365",
"361366",
"361367",
"361368",
"361369",
"361370",
"361371",
#
## Others
#
#"301535",
#"301536",
#"301890",
#"301891",
#"301892",
#"301893",
#"301894",
#"301895",
#"301896",
#"301897",
#"301898",
#"301899",
#"301900",
#"301901",
#"301902",
#"301903",
#"301904",
#"301905",
#"301906",
#"301907",
#"301908",
#"301909",
#"301910",
#"341177",
#"341270",
#"341271",
#"341988",
#"341990",
#"341991",
#"341994",
#"341996",
#"341997",
#"341999",
#"342000",
#"342005",
#"342170",
#"342171",
#"342172",
#"342284",
#"342285",
#"343266",
#"343267",
#"343268",
#"343269",
#"343270",
#"343271",
#"343272",
#"343273",
#"343274",
#"361074",
#"361075",
#"361076",
#"361078",
#"361079",
#"361080",
#"361088",
#"361089",
#"361090",
#"361091",
#"361092",
#"361093",
#"361094",
#"361095",
#"361096",
#"361097",
#"361100",
#"361101",
#"361102",
#"361103",
#"361104",
#"361105",
#"361106",
#"361107",
#"361108",
#"361500",
#"361501",
#"361502",
#"361503",
#"361504",
#"361505",
#"361506",
#"361507",
#"361508",
#"361509",
#"361510",
#"361511",
#"361512",
#"361513",
#"361514",
#"361520",
#"361521",
#"361523",
#"361525",
#"361526",
#"361527",
#"361528",
#"361529",
#"361531",
#"361533",
#"361534",
#"361628",
#"361629",
#"361630",
#"361631",
#"361632",
#"361633",
#"361634",
#"361636",
#"361637",
#"361638",
#"361639",
#"361640",
#"361641",
#"361642",
#"363102",
#"363103",
#"363104",
#"363105",
#"363106",
#"363107",
#"363108",
#"363109",
#"363110",
#"363111",
#"363112",
#"363113",
#"363114",
#"363115",
#"363116",
#"363117",
#"363118",
#"363119",
#"363120",
#"363121",
#"363122",
#"363331",
#"363332",
#"363333",
#"363334",
#"363335",
#"363336",
#"363337",
#"363338",
#"363339",
#"363340",
#"363341",
#"363342",
#"363343",
#"363344",
#"363345",
#"363346",
#"363347",
#"363348",
#"363349",
#"363350",
#"363351",
#"363352",
#"363353",
#"363354",
#"363361",
#"363362",
#"363363",
#"363364",
#"363365",
#"363366",
#"363367",
#"363368",
#"363369",
#"363370",
#"363371",
#"363372",
#"363373",
#"363374",
#"363375",
#"363376",
#"363377",
#"363378",
#"363379",
#"363380",
#"363381",
#"363382",
#"363383",
#"363384",
#"363385",
#"363386",
#"363387",
#"363388",
#"363389",
#"363390",
#"363391",
#"363392",
#"363393",
#"363394",
#"363395",
#"363396",
#"363397",
#"363398",
#"363399",
#"363400",
#"363401",
#"363402",
#"363403",
#"363404",
#"363405",
#"363406",
#"363407",
#"363408",
#"363409",
#"363410",
#"363411",
#"363436",
#"363437",
#"363438",
#"363439",
#"363440",
#"363441",
#"363442",
#"363443",
#"363444",
#"363445",
#"363446",
#"363447",
#"363448",
#"363449",
#"363450",
#"363451",
#"363452",
#"363453",
#"363454",
#"363455",
#"363456",
#"363457",
#"363458",
#"363459",
#"363460",
#"363461",
#"363462",
#"363463",
#"363464",
#"363465",
#"363466",
#"363467",
#"363468",
#"363469",
#"363470",
#"363471",
#"363472",
#"363473",
#"363474",
#"363475",
#"363476",
#"363477",
#"363478",
#"363479",
#"363480",
#"363481",
#"363482",
#"363483",
#"410001",
#"410002",
#"410007",
#"410009",
#"410015",
#"410016",
#"410025",
#"410026",
#"410049",
#"410066",
#"410067",
#"410068",
#"410073",
#"410074",
#"410075",
#"410111",
#"410112",
#"410113",
#"410114",
#"410115",
#"410116",
#"410121",
#"410142",
#"410143",
#"410144",
#"410156",
#"410157",
#"410159",
#"410187",
#"410188",
#"410189",
#"410500",
    ]

    # -------------------------------------------------------------------------------------------------------

    configpath = "$ROOTCOREBIN/user_scripts/HTopMultilepAnalysis/jobOptions_HTopMultilepMiniNTupMaker.py"
    treename   = "nominal"
    nevents    = 0

    #motherdir = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v14_testDilepTrig"
    #motherdir = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v17_test"
    #motherdir = "/coepp/cephfs/mel/mmilesi/ttH/MiniNTup/25ns_v14/25ns_v14_Direct"
    #motherdir = "/coepp/cephfs/mel/mmilesi/ttH/MiniNTup/25ns_v14/25ns_v14_Direct_DLT"
    #motherdir = "/coepp/cephfs/mel/mmilesi/ttH/MiniNTup/25ns_v17/25ns_v17_Direct"
    #motherdir = "/coepp/cephfs/mel/mmilesi/ttH/MiniNTup/25ns_v18/25ns_v18_Direct"
    #motherdir = "/coepp/cephfs/mel/mmilesi/ttH/MiniNTup/25ns_v18/25ns_v18_Skim"
    motherdir = "/coepp/cephfs/mel/mmilesi/ttH/MiniNTup/25ns_v19/25ns_v19_Skim"

    if not os.path.exists(motherdir):
        os.makedirs(motherdir)
    for s in sampledict:
        groupdir = motherdir + "/" + s["group"]
        if not os.path.exists(groupdir):
            os.makedirs(groupdir)

    list_commands = generateCmdList(infilelist)

    MAX_PARALLEL = 6

    print listchunks(list_commands,MAX_PARALLEL)

    for chunk in listchunks(list_commands,MAX_PARALLEL):

        print("Processing samples: ")
        print("\n".join("{0} - {1}".format(idx,elem[elem.find("--submitDir")+len("--submitDir")+1:elem.find("--nevents")-1]) for idx, elem in enumerate(chunk)))
        p = multiprocessing.Pool(MAX_PARALLEL)
        p.map(miniNTuplise,chunk)
	p.close()
        p.join()

