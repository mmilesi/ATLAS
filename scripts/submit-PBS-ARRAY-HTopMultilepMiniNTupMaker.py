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
## DATA
##
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
## MC
##
"341989",
"341992",
"341995",
"341998",
"342001",
"342004",
"342284",
"342285",
##"343365",
##"343366",
"343367",
"361063",
"361064",
"361065",
"361066",
"361067",
"361068",
"361069",
##"361070",
"361071",
"361072",
"361073",
"361077",
##"361091",
##"361092",
##"361093",
##"361094",
"361095",
##"361096",
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
##"361382",
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
##"361399",
"361400",
"361401",
##"361402",
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
##"361422",
"361423",
"361424",
##"361425",
"361426",
##"361427",
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
##"361479",
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
##"410000",
"410501",
"410011",
"410012",
"410013",
"410014",
"410025",
##"410026",
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
##"361300",
##"361301",
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
##"361321",
"361322",
"361323",
##"361324",
"361325",
"361326",
"361327",
"361328",
"361329",
"361330",
"361331",
"361332",
##"361333",
"361334",
"361335",
##"361336",
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
##"361350",
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

