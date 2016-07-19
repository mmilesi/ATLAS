#!/usr/bin/env python

import glob, os, sys, subprocess, shutil, string

def copy_rootcore(sample_list, params, home_dir):

    for sample in sample_list:
        tempdir = params["temp_RC_dir"]+"/xAH.{0}".format(sample)
        print("Creating temp directory:\n{0},\nand copying RootCore source in it...\n".format(tempdir))

        if os.path.exists(tempdir):
            shutil.rmtree(tempdir)

	copy_cmd = "rc grid_submit_nobuild {0}".format(tempdir)
        subprocess.call(copy_cmd,shell=True)

def create_jobs(sample_list, params, exec_job_script, steer_job_script):

    print("Creating job submission scripts...")

    job_list = []
    for sample in sample_list:

        this_temp_dir = params["temp_RC_dir"] + "/xAH.{0}".format(sample)
        this_file     = params["inputpath"] + "/{0}/{1}.root".format(sample,sample)
	outdir        = sample

        # Add items to parameters dictionary for *this* sample
        #
        params.update({"sample" : int(sample), "temp_RC_dir" : this_temp_dir, "infile" : this_file, "outdir" : outdir})

	exec_job_name = "run_xAH_job_{0}.py".format(sample)
        exec_job = open(exec_job_name, "w")
        exec_job_script = (exec_job_script).format(**params)
        exec_job.write(exec_job_script)
        exec_job.close()
        subprocess.call(["chmod", "755", exec_job_name], shell=False)

	steer_job_name = "job_xAH_{0}.pbs".format(sample)
        steer_job = open(steer_job_name, "w")
        steer_job_script = (steer_job_script).format(sample=sample,name=exec_job_name)
        steer_job.write(steer_job_script)
        steer_job.close()

	job_list.append( (steer_job_name, exec_job_name) )

    return job_list

def submit_jobs(job_list, home_dir, sub_dir):

    print("Moving submission scripts to submission node...")
    for job in job_list:
        shutil.move(os.path.abspath(home_dir)+"/"+job[0],os.path.abspath(sub_dir)+"/"+job[0])
        shutil.move(os.path.abspath(home_dir)+"/"+job[1],os.path.abspath(sub_dir)+"/"+job[1])
   
    print("Now cd'ing to submission node:\n{0}\n".format(sub_dir))
    os.chdir(os.path.abspath(sub_dir))
    
    for job in job_list:
        print("Submitting steering job script to PBS node: {0} (==> will be executing script: {1})".format(job[0],job[1]))
        subprocess.call("qsub "+job[0],shell=True)

if __name__ == '__main__':

    exec_job_script = """#!/bin/python

import os, subprocess

def runJob(temp_RC_dir,infile,configpath,treename,outdir,nevents):

    print("Setting up RootCore in temp directory: %s" % (temp_RC_dir))

    setupATLAS = "source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh"

    setup_cmd = "%s; source %s/RootCore/scripts/grid_run_nobuild.sh %s" % (setupATLAS,temp_RC_dir,temp_RC_dir)

    cmd = "python %s/RootCoreBin/bin/x86_64-slc6-gcc49-opt/xAH_run.py" % (temp_RC_dir)
    xAH_run = "%s; %s -vv --files %s --config %s --treeName %s --submitDir %s --nevents %s --force direct" % (setup_cmd,cmd,infile,configpath,treename,outdir,nevents)

    print("Running xAH job: %s" % (xAH_run))
    
    subprocess.call(xAH_run,shell=True)

if __name__ == '__main__':

    runJob( temp_RC_dir="{temp_RC_dir}", infile="{infile}", configpath="{configpath}", treename="{treename}", outdir="{outdir}", nevents={nevents})

    """

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

# Run the job from current working directory
cd $PBS_O_WORKDIR

# Run job

echo "Now into directory:"
pwd
echo "Running batch job for sample {sample}:"
echo "Printing env variables:"
env
python {name}
exit 0
    """

    # The path to the directory where RootCore is originally set, and where the packages reside
    #
    #home_dir = "/afs/cern.ch/user/m/mmilesi/ttH/RUN2/HTopMultilepAnalysisCode/trunk"
    home_dir = "/imports/home/mmilesi/PhD/ttH_MultiLeptons/RUN2/HTopMultilepAnalysisCode/trunk"
    
    # The path to the submission node (NB: must NOT be under /home : PBS cannot read/write into it!!)
    #
    sub_dir = "/coepp/cephfs/mel/mmilesi/test_batch"
    #sub_dir = string.Template("$TMPDIR/mmilesi").substitute(os.environ)
    
    samplelist = [
        #"341270",
        "343365",
        #"343366",
        #"343367",
    ]

    job_params = {
        "temp_RC_dir" : sub_dir,
	#"inputpath"   : "/afs/cern.ch/user/m/mmilesi/work/private/HTopMultileptonsTestSamples/25ns_v15/Nominal",
	"inputpath"   : "/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v17/Nominal",
        "configpath"  :  "$ROOTCOREBIN/user_scripts/HTopMultilepAnalysis/jobOptions_HTopMultilepMiniNTupMaker.py",
        "treename"    :  "nominal",
        "nevents"     :  10000,
    }

    copy_rootcore(samplelist,job_params,home_dir)

    jobs = create_jobs(samplelist,job_params,exec_job_script,steer_job_script)
    
    submit_jobs(jobs,home_dir,sub_dir)
