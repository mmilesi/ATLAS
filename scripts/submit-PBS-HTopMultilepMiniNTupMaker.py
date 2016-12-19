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
parser.add_argument('--tag', dest='tag', action='store', default=None, type=str,
                    help='Add a tag identifier to the submission directory and the destination directory (default: None)')

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

        # Append an increasing index to destination directory if the input one already exists
        # (if retry index "0" already exists, will try with idx "1", and so on...)
        #
        if os.path.exists(params["dest"]):
            append_idx = 0
            while True:
                new_dest = params["dest"] + "_Retry_" + str(append_idx)
                if not os.path.exists(new_dest):
                    params.update({"dest": new_dest})
                    break
                append_idx = append_idx + 1

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
"HTopMultilepAnalysis/doc/list-local-HTopGroupNTup.txt",
"341989",
"341992",
"341995",
"341998",
"342001",
"342004",
"342284",
"342285",
"343365",
"343366",
"343367",
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
"361091",
"361092",
"361093",
"361094",
"361095",
"361096",
"361097",
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
"361422",
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
"361620",
"361621",
"361622",
"361623",
"361624",
"361625",
"361626",
"361627",
"410000",
"410011",
"410012",
"410013",
"410014",
"410025",
"410026",
"410049",
"410073",
"410074",
"410075",
"410155",
"410156",
"410157",
"410215",
"410218",
"410219",
"410220",
    ]

    prod_ID  = args.prod_ID.split('/')
    print "prod_ID: ", prod_ID
    nevents  = args.nevents
    release  = args.release
    treename = args.treename
    queue    = args.queue
    tag      = args.tag
    if tag:
        tag = "_" + tag

    job_params = {
	"release"     :  release,
	"sample_list" :  samplelist,
	"inlistcmd"   :  "", # An xAH_run commmand to specify to read an input .txt file with a list of samples --> will be set by the script to non-dummy value only for DATA
	"sub_dir"     :  "/coepp/cephfs/mel/mmilesi/ttH/PBS/" + prod_ID[0] + "_PBS" + tag,  # The path to the submission node (NB: must NOT be under /home : PBS cannot read/write into it!!)
	"inputpath"   :  "/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/" + args.prod_ID,
	"configpath"  :  "$ROOTCOREBIN/user_scripts/HTopMultilepAnalysis/jobOptions_HTopMultilepMiniNTupMaker.py",
	"treename"    :  treename,
	"nevents"     :  int(nevents),
	"steer_job"   :  steer_job_script,
	"dest"        :  "/coepp/cephfs/mel/mmilesi/ttH/MiniNTup/" + prod_ID[0] + "/" + prod_ID[0] + "_Skim_PBS" + tag,
	"queue"       :  queue,
    }

    copy_source(sub_dir = job_params["sub_dir"], force = False)

    print("cd'ing to submission directory:\n{0}".format(job_params["sub_dir"]))
    os.chdir(os.path.abspath(job_params["sub_dir"]))
    subprocess.call(["ls","-l"])

    jobs = create_jobs(job_params)

    submit_jobs(jobs)



