import subprocess

def create_jobs(exec_job_script, steer_job_script):

    job_list = []
    for idx in range(5):
        exec_job_name = "exec_job_{0}.py".format(idx)
        exec_job = open(exec_job_name, "w")
        THIS_exec_job_script = exec_job_script % (idx)
        exec_job.write(THIS_exec_job_script)
        exec_job.close()

	steer_job_name = "job_{0}.pbs".format(idx)
        steer_job = open(steer_job_name, "w")
        THIS_steer_job_script = steer_job_script % (idx,idx,exec_job_name)
        steer_job.write(THIS_steer_job_script)
        steer_job.close()

	job_list.append( (steer_job_name,exec_job_name) )
    return job_list


def submit_jobs(job_list):

    for job_par in job_list:
       print("Submitting steering job script: {0} ==> executing script: {1}".format(job_par[0],job_par[1]))
       subprocess.call("qsub "+job_par[0],shell=True)

if __name__ == '__main__':

#    exec_job_script="""#!/bin/bash
#echo $PATH
#echo "Starting job %s:"
#COUNTER=0
#while [ $COUNTER -lt 5 ]; do
#  sleep 60
#  let COUNTER=COUNTER+1
#  echo "${COUNTER} mins have passed..."
#done
#    """

    exec_job_script="""#!/bin/python
import time
print("Starting job %s")
counter = 0
while counter < 1:
    time.sleep(60)
    counter += 1
    print("{0} mins have passed...".format(counter))
    """

    steer_job_script="""#!/bin/bash

# --- Start PBS Directives ---
# Inherit current user environment
#PBS -V

# Submit to the short queue
#PBS -q short

# Job Name
#PBS -N test_sleep_%s.pbs

# Combine the standard output and error output files
#PBS -j oe

# Email on abort and exit
#PBS -m ae
#PBS -M m.milesi@student.unimelb.edu.au
# --- End PBS Directives ---

# Run the job from current working directory
cd $PBS_O_WORKDIR

# Run job
echo "Current PBS working directory:"
echo " "
echo ${PBS_O_WORKDIR}
echo " "
cd /coepp/cephfs/mel/mmilesi/test_batch
echo "Executing batch job %s from directory:"
pwd
python %s
exit 0
    """

    jobs = create_jobs(exec_job_script, steer_job_script)
    submit_jobs(jobs)
