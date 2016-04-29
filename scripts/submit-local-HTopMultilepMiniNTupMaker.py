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
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/410007/410007.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/410066/410066.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/410067/410067.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/410068/410068.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/410011/410011.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/410012/410012.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/410013/410013.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/410014/410014.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/410016/410016.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/410049/410049.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/410050/410050.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/410073/410073.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/410074/410074.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/410075/410075.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/410080/410080.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/410081/410081.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/410111/410111.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/410112/410112.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/410113/410113.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/410114/410114.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/410115/410115.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/410116/410116.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/341177/341177.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/341270/341270.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/341271/341271.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361374/361374.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361377/361377.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361380/361380.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361383/361383.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361386/361386.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361389/361389.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361392/361392.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361395/361395.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361398/361398.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361401/361401.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361404/361404.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361407/361407.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361410/361410.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361413/361413.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361416/361416.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361419/361419.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361422/361422.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361425/361425.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361428/361428.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361431/361431.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361434/361434.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361437/361437.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361440/361440.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361443/361443.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361373/361373.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361376/361376.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361379/361379.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361382/361382.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361385/361385.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361388/361388.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361391/361391.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361394/361394.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361400/361400.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361403/361403.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361406/361406.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361409/361409.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361412/361412.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361415/361415.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361418/361418.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361421/361421.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361424/361424.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361427/361427.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361430/361430.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361433/361433.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361436/361436.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361439/361439.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361442/361442.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361372/361372.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361375/361375.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361378/361378.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361381/361381.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361384/361384.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361387/361387.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361390/361390.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361393/361393.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361396/361396.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361399/361399.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361402/361402.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361405/361405.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361408/361408.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361411/361411.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361414/361414.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361417/361417.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361420/361420.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361423/361423.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361426/361426.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361429/361429.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361432/361432.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361435/361435.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361438/361438.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361441/361441.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361468/361468.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361469/361469.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361470/361470.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361471/361471.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361472/361472.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361473/361473.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361474/361474.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361475/361475.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361476/361476.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361477/361477.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361478/361478.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361479/361479.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361480/361480.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361481/361481.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361482/361482.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361483/361483.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361484/361484.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361485/361485.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361486/361486.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361487/361487.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361488/361488.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361489/361489.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361490/361490.root",
"/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v7/Nominal/361491/361491.root",
]

# -------------------------------------------------------------------------------------------------------

configpath = "$ROOTCOREBIN/user_scripts/HTopMultilepAnalysis/jobOptions_HTopMultilepMiniNTupMaker.py"
treename   = "nominal"
nevents    = 0

motherdir = os.path.abspath(os.path.curdir) + "/25ns_v7"
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
