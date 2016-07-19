#!/usr/bin/env python

import glob, os, sys, subprocess, shutil, string

def copyRootCore(sample_list, params, home_dir):

    for sample in sample_list:
        tempdir = params["temp_RC_dir"]+"/xAH.{0}".format(sample)
        print("Creating temp directory:\n{0}...".format(tempdir))

        if not os.path.exists(tempdir):
            os.makedirs(tempdir)
        else:
            shutil.rmtree(tempdir)
            os.makedirs(tempdir)

        copy_cmd  = "rsync -arzSH --exclude='.git/' {0}/RootCoreBin {1}/xAODAnaHelpers {2}/HTopMultilepAnalysis {3}/".format(home_dir,home_dir,home_dir,tempdir)
        subprocess.call(copy_cmd,shell=True)

def create_jobs(sample_list, params, exec_job_script, steer_job_script):

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

def submit_jobs(job_list):

    for job_par in job_list:
       print("Submitting steering job script: {0} ==> executing script: {1}".format(job_par[0],job_par[1]))
       subprocess.call("qsub "+job_par[0],shell=True)

if __name__ == '__main__':

    exec_job_script = """#!/bin/python

import glob, os, sys, subprocess, shutil, string

def runJob(temp_RC_dir,infile,configpath,treename,outdir,nevents):

    print("cd'ing to temp directory: %s" % (temp_RC_dir))
    os.chdir(os.path.abspath(temp_RC_dir))

    setupATLAS = "source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh"
    rcSetup    = "source $ATLAS_LOCAL_RCSETUP_PATH/rcSetup.sh"

    setup_cmd = "%s; %s Base,2.3.51; rc make_par" % (setupATLAS,rcSetup)

    cmd = "python %s/RootCoreBin/bin/x86_64-slc6-gcc49-opt/xAH_run.py" % (temp_RC_dir)
    xAH_run = "%s; %s -vv --files %s --config %s --treeName %s --submitDir %s --nevents %s --force direct" % (setup_cmd,cmd,infile,configpath,treename,outdir,nevents)

    print("Running job: %s" % (xAH_run))
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

echo "Running batch job for sample {sample}:"
env
pwd
python {name}
exit 0
    """

    #home_dir = "/afs/cern.ch/user/m/mmilesi/ttH/RUN2/HTopMultilepAnalysisCode/trunk"
    home_dir = "/imports/home/mmilesi/PhD/ttH_MultiLeptons/RUN2/HTopMultilepAnalysisCode/trunk"

    samplelist = [
        #"341270",
        "343365",
        #"343366",
        #"343367",
    ]

    job_params = {
        #"temp_RC_dir" : string.Template("$TMPDIR/mmilesi").substitute(os.environ),
        "temp_RC_dir" : "/coepp/cephfs/mel/mmilesi/test_batch",
	#"inputpath"   : "/afs/cern.ch/user/m/mmilesi/work/private/HTopMultileptonsTestSamples/25ns_v15/Nominal",
	"inputpath"   : "/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v17/Nominal",
        "configpath"  :  "$ROOTCOREBIN/user_scripts/HTopMultilepAnalysis/jobOptions_HTopMultilepMiniNTupMaker.py",
        "treename"    :  "nominal",
        "nevents"     :  10000,
    }

    copyRootCore(samplelist,job_params,home_dir)

    jobs = create_jobs(samplelist,job_params,exec_job_script, steer_job_script)
    submit_jobs(jobs)
