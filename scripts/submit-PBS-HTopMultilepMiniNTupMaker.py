#!/usr/bin/env python

""" submit-PBS-HTopMultilepMiniNTupMaker.py: PBS NTup skimming submission script for HTopMultilepAnalysis """

__author__     = "Marco Milesi"
__email__      = "marco.milesi@cern.ch"
__maintainer__ = "Marco Milesi"

import glob, os, sys, subprocess, shutil, string, argparse

parser = argparse.ArgumentParser(description='PBS NTup skimming submission script for HTopMultilepAnalysis')

parser.add_argument('--prod_ID', dest='prod_ID', action='store', default='25ns_v19', type=str,
                    help='The NTup production tag (default: prod_ID=\'25ns_v19\')')
parser.add_argument('--nevents', dest='nevents', action='store', default=0, type=int,
                    help='The number of events to run on (default: nevents=0 [ALL])')
parser.add_argument('--release', dest='release', action='store', default='2.3.51', type=str,
                    help='The ASG release to be used (default: release=\'2.3.51\')')
parser.add_argument('--treename', dest='treename', action='store', default='nominal', type=str,
                    help='The input TTree name (default: treename=\'nominal\')')
parser.add_argument('--queue', dest='queue', action='store', default='short', type=str,
                    help='The PBS batch queue type to be used (\'short\',\'long\', default: queue=\'short\')')

args = parser.parse_args()

def copy_source(sub_dir = string.Template("$TMPDIR").substitute(os.environ), force = False):

    if os.path.exists(sub_dir) and force:
        shutil.rmtree(sub_dir)

    print("Creating submission directory:\n{0},\nand copying source code in it...\n".format(sub_dir))
    if not os.path.exists(sub_dir):
        os.makedirs(sub_dir)
    else:
        print("Good! Submission directory already exists...")

    tarballname = "xAH.tar.gz"
    print("Making a tarball of code and moving it to submission directory...")
    if not os.path.isfile(sub_dir+"/"+tarballname):
        subprocess.call(["tar","-zcf",tarballname,"xAODAnaHelpers/","HTopMultilepAnalysis/","RootCoreBin/","--exclude-vcs"])
        shutil.move("./"+tarballname,sub_dir+"/"+tarballname)
    else:
        print("Good! Tarball of code already exists in submission directory...")


def create_jobs(params):

    print("Creating job submission scripts...")

    job_list = []
    for sample in params["sample_list"]:

	inputpath = params["inputpath"] + "/Nominal"

	# Special treatment for DATA: need an input .txt list of datasets
	#
        if ( sample[-4:] == ".txt" ):
	    infile = sample
	    params.update({"inlistcmd" : "--inputList"})
            sample = "Data"
        else:
	    infile    = inputpath + "/{0}/{1}.root".format(sample,sample)
	    params.update({"inlistcmd" : ""})

	outdir = sample

        # Update items to parameters dictionary for *this* sample
        #
        params.update({"sample" : sample, "infile" : infile, "outdir" : outdir})

        steer_job_name = "job_xAH_{0}.pbs".format(sample)
        steer_job = open(steer_job_name, "w")
        steer_job_script = params["steer_job"]
        steer_job_script = steer_job_script.format(**params)

        steer_job.write(steer_job_script)
        steer_job.close()

        job_list.append( steer_job_name )

    return job_list

def submit_jobs(jobs):

    for job in jobs:
        print("Submitting steering job script {0} to PBS node...".format(job))
        subprocess.call(["qsub",job])

if __name__ == '__main__':

    steer_job_script="""#!/bin/bash

# --- Start PBS Directives ---

# Submit to the long/short queue
#PBS -q {queue}

# Job Name
#PBS -N skim_xAH_{sample}.pbs

# Combine the standard output and error output files
#PBS -j oe

# Email on abort and exit
#PBS -m ae
#PBS -M m.milesi@student.unimelb.edu.au

# --- End PBS Directives ---

# Run job

# Run the job from current working directory
cd $PBS_O_WORKDIR

echo "Running on host" $(hostname)
echo "Time is" $(date)
echo "Current directory is" $(pwd)
echo "This jobs runs on the following processors:"
NODES=$(cat $PBS_NODEFILE)
echo $NODES
NPROCS=$(wc -l < $PBS_NODEFILE)
echo "This job has allocated $NPROCS nodes"
echo "Running batch job for sample {sample}"
echo ""
TMP=$PWD/tmp.{sample}
echo "Creating temporary directory for this job: $TMP"
if [ -d "$TMP" ]; then
    echo "Directory already found! Removing it first..."
    rm -rf $TMP
fi
mkdir $TMP
echo "Copying tarballed code and cd'ing into job directory..."
rsync -arvxSH xAH.tar.gz $TMP/
cd $TMP/
echo "Opening tarball..."
tar -zxf xAH.tar.gz
rm -rf xAH.tar.gz
echo "Setting up ATLAS software..."
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh
echo "Setting up RootCore and ASG..."
source $ATLAS_LOCAL_RCSETUP_PATH/rcSetup.sh Base,{release}
rc find_packages
rc compile
echo "Printing env variables:"
env
echo "Running xAH job:"
python $PWD/RootCoreBin/bin/x86_64-slc6-gcc49-opt/xAH_run.py -vv --files {infile} {inlistcmd} --config {configpath} --treeName {treename} --submitDir {outdir} --nevents {nevents} --force direct
python $PWD/HTopMultilepAnalysis/scripts/weightTrees-PBS.py {dest} {sample}
echo "Removing temporary directory..."
rm -rf $TMP
exit 0
    """

    samplelist = [
## data
"HTopMultilepAnalysis/doc/list-local-HTopGroupNTup.txt",
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
"361091",
"361092",
"361093",
"361094",
"361095",
"361096",
"361097",
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
#
# v20-02, v20-03
#
#"341990",
#"341271",
#"341989",
#"341177",
#"341988",
#"341270",
#"342001",
#"341992",
#"341995",
#"341999",
#"341997",
#"341998",
#"342284",
#"342170",
#"342172",
#"342285",
#"342171",
#"342004",
#"343367",
#"343266",
#"343268",
#"343267",
#"343365",
#"343366",
#"361065",
#"361068",
#"361063",
#"361067",
#"361066",
#"361064",
#"361071",
#"361074",
#"361073",
#"361070",
#"361072",
#"361069",
#"361078",
#"361077",
#"361079",
#"361076",
#"361075",
#"361080",
#"361081",
#"361086",
#"361082",
#"361084",
#"361083",
#"361085",
#"361087",
#"361089",
#"361090",
#"361088",
#"361091",
#"361092",
#"361095",
#"410000",
#"361096",
#"361094",
#"361093",
#"361097",
#"410002",
#"410001",
#"410007",
#"410009",
#"410011",
#"410012",
#"410016",
#"410013",
#"410014",
#"410026",
#"410025",
#"410015",
#"410073",
#"410068",
#"410049",
#"410050",
#"410067",
#"410066",
#"410080",
#"410081",
#"410074",
#"410111",
#"410075",
#"410112",
#"410113",
#"410121",
#"410120",
#"410115",
#"410114",
#"410116",
#"410143",
#"410144",
#"410156",
#"410142",
#"410157",
#"410155",
#"410188",
#"410159",
#"410187",
#"410218",
#"410215",
#"410189",
#"410500",
#"410220",
#"410219",
    ]

    prod_ID  = args.prod_ID
    nevents  = args.nevents
    release  = args.release
    treename = args.treename
    queue    = args.queue

    job_params = {
	"release"     :  release,
	"sample_list" :  samplelist,
	"inlistcmd"   :  "", # An xAH_run commmand to specify to read an input .txt file with a list of samples --> will be set by the script to non-dummy value only for DATA
	"sub_dir"     :  "/coepp/cephfs/mel/mmilesi/ttH/PBS/" + prod_ID + "_PBS",  # The path to the submission node (NB: must NOT be under /home : PBS cannot read/write into it!!)
	"inputpath"   :  "/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/" + prod_ID,
	"configpath"  :  "$ROOTCOREBIN/user_scripts/HTopMultilepAnalysis/jobOptions_HTopMultilepMiniNTupMaker.py",
	"treename"    :  treename,
	"nevents"     :  int(nevents),
	"steer_job"   :  steer_job_script,
	"dest"        :  "/coepp/cephfs/mel/mmilesi/ttH/MiniNTup/" + prod_ID + "/" + prod_ID + "_Skim_PBS",
	"queue"       :  queue,
    }

    copy_source(sub_dir = job_params["sub_dir"], force = False)

    print("cd'ing to submission directory:\n{0}".format(job_params["sub_dir"]))
    os.chdir(os.path.abspath(job_params["sub_dir"]))
    subprocess.call(["ls","-l"])

    jobs = create_jobs(job_params)

    submit_jobs(jobs)



