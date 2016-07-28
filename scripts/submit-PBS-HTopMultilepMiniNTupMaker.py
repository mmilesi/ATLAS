#!/usr/bin/env python

import glob, os, sys, subprocess, shutil, string

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

        infile = params["inputpath"] + "/{0}/{1}.root".format(sample,sample)
        outdir = sample

        # Add items to parameters dictionary for *this* sample
        #
        params.update({"sample" : int(sample), "infile" : infile, "outdir" : outdir})

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
#PBS -q short

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
python $PWD/RootCoreBin/bin/x86_64-slc6-gcc49-opt/xAH_run.py -vv --files {infile} --config {configpath} --treeName {treename} --submitDir {outdir} --nevents {nevents} --force direct
exit 0
    """

    samplelist = [
        "341270",
        "343365",
        #"343366",
        #"343367",
    ]

    job_params = {
        "release"     : "2.3.51",
        "sample_list" :  samplelist,
        "sub_dir"     :  "/coepp/cephfs/mel/mmilesi/test_batch",  # The path to the submission node (NB: must NOT be under /home : PBS cannot read/write into it!!)
        #"inputpath"  :  "/afs/cern.ch/user/m/mmilesi/work/private/HTopMultileptonsTestSamples/25ns_v15/Nominal",
        "inputpath"   :  "/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v19/Nominal",
        "configpath"  :  "$ROOTCOREBIN/user_scripts/HTopMultilepAnalysis/jobOptions_HTopMultilepMiniNTupMaker.py",
        "treename"    :  "nominal",
        "nevents"     :  10000,
        "steer_job"   :  steer_job_script,
    }

    copy_source(sub_dir = job_params["sub_dir"], force = False)

    print("cd'ing to submission directory:\n{0}".format(job_params["sub_dir"]))
    os.chdir(os.path.abspath(job_params["sub_dir"]))
    subprocess.call(["ls","-l"])

    jobs = create_jobs(job_params)

    submit_jobs(jobs)
