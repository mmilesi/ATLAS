#!/usr/bin/env python

""" submit-PBS-ARRAY-HTopMultilepMiniNTupMaker.py: PBS NTup skimming submission script for HTopMultilepAnalysis """

__author__     = "Marco Milesi"
__email__      = "marco.milesi@cern.ch"
__maintainer__ = "Marco Milesi"

import glob, os, sys, subprocess, shutil, string, argparse

parser = argparse.ArgumentParser(description="PBS NTup skimming submission script for HTopMultilepAnalysis. Will submit a single main PBS job made up of an array of subjobs, one for each sample ID.")

parser.add_argument("--prod_ID", dest="prod_ID", action="store", default="25ns_v19", type=str,
                    help="The NTup production tag, e.g. 25ns_v19, 25ns_v24/02, ...  (default: prod_ID=25ns_v19)")
parser.add_argument("--nevents", dest="nevents", action="store", default=0, type=int,
                    help="The number of events to run on (default: nevents=0 [ALL])")
parser.add_argument("--release", dest="release", action="store", default="2.4.22", type=str,
                    help="The ASG release to be used (default: release=\"2.4.22\")")
parser.add_argument("--treename", dest="treename", action="store", default="nominal", type=str,
                    help="The input TTree name (default: treename=\"nominal\")")
parser.add_argument("--queue", dest="queue", action="store", default="short", type=str,
                    help="The PBS batch queue type to be used (\"short\",\"long\", default: queue=\"short\")")
parser.add_argument("--tag", dest="tag", action="store", default=None, type=str,
                    help="Add a tag identifier to the submission directory and the destination directory (default: None)")
parser.add_argument("--showsamples", dest="showsamples", action="store", const="-1", default=None, type=str, nargs='?',
                    help="Show on screen the list of sample IDs, with their corresponding indexes in the job array, then exit. If no argument is given to the option ( or using --showsamples -1 ), will print the entire list. If the user gives one ore more comma-separated indexes as option, only the corresponding samples will be printed.)")

args = parser.parse_args()

def copy_source(subdir = string.Template("$TMPDIR").substitute(os.environ), force = False):

    if os.path.exists(subdir) and force:
        shutil.rmtree(subdir)

    print("Creating submission directory:\n{0},\nand copying source code in it...".format(subdir))
    if not os.path.exists(subdir):
        os.makedirs(subdir)
    else:
        print("Good! Submission directory already exists...")

    tarballname = "xAH.tar.gz"
    print("Making a tarball of code and moving it to submission directory...")
    if not os.path.isfile(subdir+"/"+tarballname):
        subprocess.call(["tar","-zcf",tarballname,"xAODAnaHelpers/","HTopMultilepAnalysis/","RootCoreBin/","--exclude-vcs"])
        shutil.move("./"+tarballname,subdir+"/"+tarballname)
    else:
        print("Good! Tarball of code already exists in submission directory...")


def create_jobs(params):

    job_list = []

    params.update( {"njobs" : len(params["samplelist"]) })
    params.update( {"upperjobidx" : len(params["samplelist"])-1 })

    # Transform the sample list into a single string of space-separated sample IDs
    # to pass as an argument to the PBS job script.

    samples_str = " ".join( params["samplelist"] )
    params.update( {"samples" : samples_str} )

    # Here would go an hypotethical loop over systematics...

    steer_job_name = "job_xAH.pbs"
    steer_job = open(steer_job_name, "w")
    steer_job_script = params["steer_job"]
    steer_job_script = steer_job_script.format(**params)

    steer_job.write(steer_job_script)
    steer_job.close()

    print("Creating job submission scripts w/ parameters:\n")
    for key,value in params.iteritems():
        if key == "steer_job" : continue
        print("{0} : {1}".format(key,value))
    print("")

    job_list.append( steer_job_name )

    return job_list

def submit_jobs(jobs):

    for job in jobs:
        print("\nSubmitting steering job script {0} to PBS node...".format(job))
        subprocess.call(["qsub",job])

if __name__ == '__main__':

    # This is the PBS script that will be actually executed on the batch node

    steer_job_script="""#!/bin/bash

# --- Start PBS Directives ---

# Submit to the long/short queue
#PBS -q {queue}

# Job Name
#PBS -N skim_xAH

# Set the number of subjobs in the array
#PBS -t 0-{upperjobidx}

# Combine the standard output and error output files
#PBS -j oe

# Email on abort and exit
#PBS -m ae
#PBS -M m.milesi@student.unimelb.edu.au

# --- End PBS Directives ---

# Run job

# Run the job from current working directory
cd $PBS_O_WORKDIR

echo "Running array of batch xAH jobs - subjob index: PBS_ARRAYID="$PBS_ARRAYID
echo ""
echo "Running on host:" $(hostname)
echo ""
echo "Time is:" $(date)
echo ""
echo "Current directory is:" $(pwd)
echo ""
echo "This job runs on the following processors:"
echo ""
NODES=$(cat $PBS_NODEFILE)
echo $NODES
NPROCS=$(wc -l < $PBS_NODEFILE)
echo ""
echo "This job has allocated $NPROCS nodes"
echo ""
TMP=`mktemp -d $TMPDIR/mmilesi.$PBS_ARRAYID.XXXX`
echo "Creating temporary directory for this subjob: $TMP"
if [ -d "$TMP" ]; then
    echo "Directory already found! Removing it first..."
    rm -rf $TMP
fi
mkdir $TMP
echo "Copying tarballed code and cd'ing into temp subjob directory..."
echo ""
rsync -arvxSH xAH.tar.gz $TMP/
cd $TMP/
echo "Opening tarball..."
echo ""
tar -zxf xAH.tar.gz
rm -rf xAH.tar.gz
echo ""
echo "Setting up ATLAS software..."
echo ""
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh
echo ""
echo "Setting up RootCore and ASG..."
echo ""
source $ATLAS_LOCAL_RCSETUP_PATH/rcSetup.sh Base,{release}
rc find_packages
rc compile
echo ""
echo "Printing env variables:"
echo ""
env
echo ""

# The following script wraps around the xAH_run.py script, picking the sample for the PBS subjob from the input sample list via the env variable os.getenv('PBS_ARRAYID'). It also executes the script to reweight the final trees.

python $PWD/HTopMultilepAnalysis/scripts/wrapper-xAH-PBS.py --inputpath {inputpath} --outputpath {outputpath} --samplelist {samples} --config {configpath} --treename {treename} --nevents {nevents}

echo "Removing temporary directory..."
rm -rf $TMP
exit 0
    """

    # List of samples to be processed (one sample will correspond to one subjob in the array)

    samplelist = [
# DATA
#
"00276262",
"00276329",
# MC
#
"410000",
"410011",
"410012",
    ]

    if args.showsamples:
        showlist = [ int(s) for s in args.showsamples.split(',') ]
        print("List of sample IDs and their indexes in job list:\n")
        for idx, s in enumerate(samplelist):
            if not ( len(showlist) == 1 and showlist[0] == -1 ):
                if idx not in showlist: continue
            print("job[{0}]={1}".format(idx,s))
        sys.exit()

    prod_ID  = args.prod_ID.split('/')
    print "prod_ID: ", prod_ID
    nevents  = args.nevents
    release  = args.release
    treename = args.treename
    queue    = args.queue
    tag      = args.tag
    if tag:
	tag = "_" + tag

    # Set the job parameters

    job_params = {
	"release"     :  release,
	"samplelist"  :  samplelist,
	"subdir"      :  "/coepp/cephfs/mel/mmilesi/ttH/PBS/" + prod_ID[0] + "_PBS" + tag,  # The path to the submission node (NB: must NOT be under /home : PBS cannot read/write into it!!)
	"inputpath"   :  "/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/" + args.prod_ID,
	"outputpath"  :  "/coepp/cephfs/mel/mmilesi/ttH/MiniNTup/" + prod_ID[0] + "/" + prod_ID[0] + "_Skim_PBS" + tag,
	"configpath"  :  "$ROOTCOREBIN/user_scripts/HTopMultilepAnalysis/jobOptions_HTopMultilepMiniNTupMaker.py",
	"treename"    :  treename,
	"nevents"     :  int(nevents),
	"steer_job"   :  steer_job_script,
	"queue"       :  queue,
    }

    # Copying source code onto the submission node, and cd in there

    copy_source(subdir = job_params["subdir"], force = False)

    print("cd'ing to submission directory:\n{0}".format(job_params["subdir"]))
    os.chdir(os.path.abspath(job_params["subdir"]))
    subprocess.call(["ls","-l"])

    # Create a list of PBS job scripts.
    # Producing a list is useful if for example one wants to submit several main jobs for each systematic,
    # where each job will consist of an array of subjobs (one for each sample ID).

    jobs = create_jobs(job_params)

    # Finally, execute the PBS script

    submit_jobs(jobs)















































