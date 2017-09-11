#!/usr/bin/env python

""" submit-PBS-ARRAY-ttHMultiGFW2.py: PBS NTup skimming submission script for GFW2 """

__author__     = "Marco Milesi"
__email__      = "marco.milesi@cern.ch"
__maintainer__ = "Marco Milesi"

import glob, os, sys, subprocess, shutil, string, argparse

parser = argparse.ArgumentParser(description="PBS NTup skimming submission script for GFW2. Will submit a single main PBS job made up of an array of subjobs, one for each sample ID.")

parser.add_argument("--prod_ID", dest="prod_ID", action="store", default="25ns_v19", type=str,
                    help="The NTup production tag, e.g. 25ns_v19, 25ns_v24/02, ...  (default: prod_ID=25ns_v19)")
parser.add_argument("--closure", dest="closure", action="store_true", default=False,
                    help="Run jobs on MC ttbar,ttgamma samples (for closure test)")
parser.add_argument("--nevents", dest="nevents", action="store", default=0, type=int,
                    help="The number of events to run on (default: nevents=0 [ALL])")
parser.add_argument("--release", dest="release", action="store", default="2.4.30", type=str,
                    help="The ASG release to be used (default: release=\"2.4.30\")")
parser.add_argument("--queue", dest="queue", action="store", default="long", type=str,
                    help="The PBS batch queue type to be used (\"short\",\"long\", default: queue=\"long\")")
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

    tarballname = "ttHMultiGFW2.tar.gz"
    print("Making a tarball of code and moving it to submission directory...")
    if not os.path.isfile(subdir+"/"+tarballname):
        subprocess.call(["tar","-zcf",tarballname,"ttHMultilepton/","RootCoreBin/","Ext_XGBoost/","Ext_yaml/","--exclude-vcs"])
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

    steer_job_name = "job_ttHMultiGFW2.pbs"
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
#PBS -N skim_ttHMultiGFW2

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

echo "Running array of batch ttHMultiGFW2 jobs - subjob index: PBS_ARRAYID="$PBS_ARRAYID
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
echo "Creating temporary directory for this subjob"
TMP=`mktemp -d $TMPDIR/mmilesi.$PBS_ARRAYID.XXXX`
echo "Temporary directory is $TMP"
if [ ! -d "$TMP" ]; then
    echo "Directory not created! Bad! Exiting"
    exit 1
fi
echo "Copying tarballed code and cd'ing into temp subjob directory..."
echo ""
rsync -arvxSH ttHMultiGFW2.tar.gz $TMP/
cd $TMP/
echo "Opening tarball..."
echo ""
tar -zxf ttHMultiGFW2.tar.gz
rm -rf ttHMultiGFW2.tar.gz
echo ""
echo "Setting up ATLAS software..."
echo ""
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh
echo ""
echo "Setting up RootCore and ASG..."
echo ""
source $ATLAS_LOCAL_RCSETUP_PATH/rcSetup.sh Top,{release}
# rc clean
rc find_packages
rc compile
echo ""
echo "Printing env variables:"
echo ""
env
echo ""

# The following script wraps around the testGFW2Tool script, picking the sample for the PBS subjob from the input sample list via the env variable os.getenv('PBS_ARRAYID'). It also executes the script to reweight the final trees.

ls -l /coepp/cephfs/mel/mmilesi/ttH/
python $PWD/ttHMultilepton/scripts/wrapper-ttHMultiGFW2-PBS.py --inputpath {inputpath} --outputpath {outputpath} --samplelist {samples} --config {configpath} --nevents {nevents}

echo "Removing temporary directory..."
rm -rf $TMP
exit 0
    """

    # List of samples to be processed (one sample will correspond to one subjob in the array)

    samplelist = [
"00276262",
"00276329",
"00276336",
"00276416",
"00276511",
"00276689",
"00276778",
"00276790",
"00276952",
"00276954",
"00278880",
"00278912",
"00278968",
"00279169",
"00279259",
"00279279",
"00279284",
"00279345",
"00279515",
"00279598",
"00279685",
"00279813",
"00279867",
"00279928",
"00279932",
"00279984",
"00280231",
"00280273",
"00280319",
"00280368",
"00280423",
"00280464",
"00280500",
"00280520",
"00280614",
"00280673",
"00280753",
"00280853",
"00280862",
"00280950",
"00280977",
"00281070",
"00281074",
"00281075",
"00281317",
"00281385",
"00281411",
"00282625",
"00282631",
"00282712",
"00282784",
"00282992",
"00283074",
"00283155",
"00283270",
"00283429",
"00283608",
"00283780",
"00284006",
"00284154",
"00284213",
"00284285",
"00284420",
"00284427",
"00284484",
"00297730",
"00298595",
"00298609",
"00298633",
"00298687",
"00298690",
"00298771",
"00298773",
"00298862",
"00298967",
"00299055",
"00299144",
"00299147",
"00299184",
"00299243",
"00299584",
"00300279",
"00300345",
"00300415",
"00300418",
"00300487",
"00300540",
"00300571",
"00300600",
"00300655",
"00300687",
"00300784",
"00300800",
"00300863",
"00300908",
"00301912",
"00301918",
"00301932",
"00301973",
"00302053",
"00302137",
"00302265",
"00302269",
"00302300",
"00302347",
"00302380",
"00302391",
"00302393",
"00302737",
"00302831",
"00302872",
"00302919",
"00302925",
"00302956",
"00303007",
"00303079",
"00303201",
"00303208",
"00303264",
"00303266",
"00303291",
"00303304",
"00303338",
"00303421",
"00303499",
"00303560",
"00303638",
"00303832",
"00303846",
"00303892",
"00303943",
"00304006",
"00304008",
"00304128",
"00304178",
"00304198",
"00304211",
"00304243",
"00304308",
"00304337",
"00304409",
"00304431",
"00304494",
"00305380",
"00305543",
"00305571",
"00305618",
"00305671",
"00305674",
"00305723",
"00305727",
"00305735",
"00305777",
"00305811",
"00305920",
"00306269",
"00306278",
"00306310",
"00306384",
"00306419",
"00306442",
"00306448",
"00306451",
"00307126",
"00307195",
"00307259",
"00307306",
"00307354",
"00307358",
"00307394",
"00307454",
"00307514",
"00307539",
"00307569",
"00307601",
"00307619",
"00307656",
"00307710",
"00307716",
"00307732",
"00307861",
"00307935",
"00308047",
"00308084",
"00309375",
"00309390",
"00309440",
"00309516",
"00309640",
"00309674",
"00309759",
"00310015",
"00310247",
"00310249",
"00310341",
"00310370",
"00310405",
"00310468",
"00310473",
"00310634",
"00310691",
"00310738",
"00310809",
"00310863",
"00310872",
"00310969",
"00311071",
"00311170",
"00311244",
"00311287",
"00311321",
"00311365",
"00311402",
"00311473",
"00311481",
    ]

    if args.closure:
        del samplelist[:]
        samplelist.extend(["410082","410501","410503"])
        #samplelist.extend(["410082","410000"])
        #samplelist.extend(["410082","410250","410251","410252"])

    if args.showsamples:
        showlist = [ int(s) for s in args.showsamples.split(',') ]
        print("List of sample IDs and their indexes in job list:\n")
        for idx, s in enumerate(samplelist):
            if not ( len(showlist) == 1 and showlist[0] == -1 ):
                if idx not in showlist: continue
            print("job[{0}]={1}".format(idx,s))
        os.sys.exit()

    prod_ID  = args.prod_ID.split('/')
    print "prod_ID: ", prod_ID
    nevents  = args.nevents
    release  = args.release
    queue    = args.queue
    tag      = args.tag
    if tag:
	tag = "_" + tag
    else:
        tag = ""

    conf_id = "data"
    if args.closure:
        conf_id = "MC"

    # Set the job parameters

    job_params = {
	"release"     :  release,
	"samplelist"  :  samplelist,
	"subdir"      :  "/coepp/cephfs/mel/mmilesi/ttH/PBS/GFW2MiniNTup/" + prod_ID[0] + "_PBS" + tag,  # The path to the submission node (NB: must NOT be under /home : PBS cannot read/write into it!!)
	"inputpath"   :  "/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/" + args.prod_ID,
	"outputpath"  :  "/coepp/cephfs/mel/mmilesi/ttH/GFW2MiniNTup/" + prod_ID[0] + "/Fakes/" + prod_ID[0] + "_Skim_PBS" + tag,
	"configpath"  :  "$ROOTCOREBIN/../ttHMultilepton/data/generic_" + conf_id + "_Melbourne.conf",
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

