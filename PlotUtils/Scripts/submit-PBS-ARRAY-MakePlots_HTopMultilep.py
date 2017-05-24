#!/usr/bin/env python

""" submit-PBS-ARRAY-MakePlots_HTopMultilep.py: PBS submission script for making plots """

__author__     = "Marco Milesi"
__email__      = "marco.milesi@cern.ch"
__maintainer__ = "Marco Milesi"

import glob, os, sys, subprocess, shutil, string, argparse

parser = argparse.ArgumentParser(description="PBS plot making submission script for HTopMultilepAnalysis. Will submit a single main PBS job made up of an array of subjobs, one for each variable we want to plot.")

parser.add_argument("--outputdir", dest="outputdir", action="store", default="PLOTS_TEST", type=str,
                    help="The base directory where output of each job will be stored,, e.g. PLOTS_25ns_v27,...  (default: PLOTS_TEST)")
parser.add_argument("-f","--forcetarball", dest="forcetarball", action="store_true", default=False,
                    help="Force recreation of package tarball")
parser.add_argument("--optstr", dest="optsr", action="store", type=str,
                    help='A string representing the command options for Plotter/MakePlots_HTopMultilep.py, including the input directory, the channel for which plots should be made, etc. To see all the options, just type \'python Plotter/MakePlots_HTopMultilep.py --help\'')
parser.add_argument("--tag", dest="tag", action="store", default=None, type=str,
                    help="Add a tag identifier to the submission directory (default: None)")
parser.add_argument("--release", dest="release", action="store", default="2.4.22", type=str,
                    help="The ASG release to be used (default: release=\"2.4.22\")")
parser.add_argument("--queue", dest="queue", action="store", default="long", type=str,
                    help="The PBS batch queue type to be used (\"short\",\"long\", default: queue=\"long\")")
parser.add_argument("--showvars", dest="showvars", action="store", const="-1", default=None, type=str, nargs='?',
                    help="Show on screen the list of variable names, with their corresponding indexes in the job array, then exit. If no argument is given to the option ( or using --showvars -1 ), will print the entire list. If the user gives one ore more comma-separated indexes as option, only the corresponding variable names will be printed.)")
parser.add_argument("--dry", dest="dry", action="store_true", default=False,
                    help="Dry-run")

args = parser.parse_args()

def copy_source(subdir = string.Template("$TMPDIR").substitute(os.environ), force = False):

    if os.path.exists(subdir) and force:
        shutil.rmtree(subdir)

    print("Creating submission directory:\n{0},\nand copying source code in it...".format(subdir))
    if not os.path.exists(subdir):
        os.makedirs(subdir)
    else:
        print("Good! Submission directory already exists...")

    tarballname = "HTop.tar.gz"

    # The directory where the inut package lives

    # WARNING: the following works only if RootCore has been set up previously on the submission node...
    # rcdir = string.Template("$ROOTCOREBIN").substitute(os.environ)
    # tokens = rcdir.split('/')
    # pckgdir = "/".join( "{0}".format(tk) for tk in tokens[:-1] )

    pckgdir = "/imports/home/mmilesi/PhD/ttH_MultiLeptons/RUN2/HTopMultilepAnalysisCode/trunk"

    # Remove some ROOT object and shared libraries files before tarballing and copying

    for filetype in ["*.so","*.d","*.pcm"]:
        map(os.remove, glob.glob(pckgdir+"/RootCoreBin/user_scripts/HTopMultilepAnalysis/ROOT_TTreeFormulas/"+filetype))
        map(os.remove, glob.glob(pckgdir+"/HTopMultilepAnalysis/scripts/ROOT_TTreeFormulas/"+filetype))

    os.chdir(os.path.abspath(pckgdir))

    print("Making a tarball of input package and moving it to submission directory...")
    if not os.path.isfile(subdir+"/"+tarballname) or args.forcetarball:
        subprocess.call(["tar","-zcf",tarballname,"HTopMultilepAnalysis/","RootCoreBin/","--exclude-vcs","--exclude=\'Output*\'"])
        shutil.move("./"+tarballname,subdir+"/"+tarballname)
    else:
        print("Good! Tarball of code already exists in submission directory...")


def create_jobs(params):

    job_list = []

    params.update( {"njobs" : len(params["varlist"]) })
    params.update( {"upperjobidx" : len(params["varlist"])-1 })

    # Transform the plotting variable list into a single string of space-separated variable names
    # to pass as an argument to the PBS job script.

    vars_str = " ".join( params["varlist"] )
    params.update( {"vars" : vars_str} )

    # Here would go an hypotethical loop over systematics...

    steer_job_name = "job_MakePlots_HTopMultilep.pbs"
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
#PBS -N makePlots_HTopMultlep

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

echo "Running array of batch plotting jobs - subjob index: PBS_ARRAYID="$PBS_ARRAYID
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
    echo "Directory w/ this name already found! Removing it first, and recreate it..."
    rm -rf $TMP
fi
mkdir $TMP
echo ""
echo "Temp dir: " $TMP
echo ""
echo "Copying tarballed package and cd'ing into temp subjob directory..."
echo ""
rsync -arvxSH HTop.tar.gz $TMP/
cd $TMP/
echo "Opening tarball..."
echo ""
tar -zxf HTop.tar.gz
rm -rf HTop.tar.gz
echo ""
echo "Setting up ATLAS software..."
echo ""
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh
echo ""
echo "Setting up RootCore and ASG..."
echo ""
source $ATLAS_LOCAL_RCSETUP_PATH/rcSetup.sh Base,{release}
# rc find_packages
# rc compile
echo ""
echo "Printing env variables:"
echo ""
env
echo ""

# The following script wraps around the MakePlots_HTopMultilep.py script, picking the plotting variable for the PBS subjob from the input plotting variable list via the env variable os.getenv('PBS_ARRAYID'). It also copies the output to a more suitable directory.

# Clear the cache first
ls -l /coepp/cephfs/mel/mmilesi/ttH/

python $PWD/HTopMultilepAnalysis/PlotUtils/Scripts/wrapper-MakePlots_HTopMultilep-PBS.py --optstr="{optstr}" --varlist {vars} --outputpath {outputpath}

echo "Removing temporary directory..."
# rm -rf $TMP
exit 0
    """

    # List of samples to be processed (one sample will correspond to one subjob in the array)

    varlist = [
        "Integral",
        # "Integral_LOGY",
        "ElProbePt",
        "ElProbePt_LOGY",
        "ElProbePt_RealEffBinning",
        "ElProbePt_FakeEffBinning",
        "ElProbePt_RealEffBinning_LOGY",
        "ElProbePt_FakeEffBinning_LOGY",
        "ElProbeEta",
        "ElProbeDistanceClosestJet",
        # "ElProbeDistanceClosestBJet",
        # "ElProbeDistanceOtherLep",
        # "ElProbeNJets",
        # "ElProbeNBJets",
        "ElProbeNBJets_VS_ElProbePt",
        "ElProbeDistanceClosestJet_VS_ElProbePt",
        # "ElProbeEta_VS_ElProbePt",
        "MuProbePt",
        "MuProbePt_LOGY",
        "MuProbePt_RealEffBinning",
        "MuProbePt_FakeEffBinning",
        "MuProbePt_RealEffBinning_LOGY",
        "MuProbePt_FakeEffBinning_LOGY",
        "MuProbeEta",
        "MuProbeDistanceClosestJet",
        # "MuProbeDistanceClosestBJet",
        # "MuProbeDistanceOtherLep",
        # "MuProbeNJets",
        # "MuProbeNBJets",
        # "MuProbeNBJets_VS_MuProbePt",
        "MuProbeDistanceClosestJet_VS_MuProbePt",
        #
        # "LepFlavours",
        #
        # "Integral",
        # "NJets2j3j",
        # "NJets4j",
        # "El0Pt",
        # "El1Pt",
        # "Mu0Pt",
        # "Mu1Pt",
        # "El0Eta",
        # "El1Eta",
        # "Mu0Eta",
        # "Mu1Eta",
        # "El0DeltaRClosestJet",
        # "El1DeltaRClosestJet",
        # "Mu0DeltaRClosestJet",
        # "Mu1DeltaRClosestJet",
        # "NBJets",
        # "Mll01_inc",
        # "MET_FinalTrk",
        # "deltaRLep0Lep1",
        # "deltaPhiLep0Lep1",
        # "TotLepCharge",
        #
        # "NN_Rebinned",
        # "RNN_Rebinned",
        # "NN_ttV",
        # "NN_top",
        # "RNN_ttV",
        # "RNN_top",
        #
        # "Lep1Type",
        # "Lep1Origin",
        # 'Lep1Type_VS_Lep1Origin',
        # "Lep1Origin_VS_NJets",
        # "Lep1Origin_VS_Lep1Pt",
        # "Lep1Origin_VS_Lep1Eta",
        # "Lep1Origin_VS_Lep1EtaBE2",
        # "Lep1Origin_VS_Lep1DistanceClosestJet",
        # "Lep2Type",
        # "Lep2Origin",
        # 'Lep2Type_VS_Lep2Origin',
        # "Lep2Origin_VS_NJets",
        # "Lep2Origin_VS_Lep2Pt",
        # "Lep2Origin_VS_Lep2Eta",
        # "Lep2Origin_VS_Lep2EtaBE2",
        # "Lep2Origin_VS_Lep2DistanceClosestJet",
        #
        # "MassMuMuEl_LOGY",
        # "MassMuMu_LOGY",
        # "mTElMET_LOGY",
        # "ElProbePt",
        # "ElProbePt_LOGY",
        # "ElProbePt_FakeEffBinning",
        # "ElProbePt_RealEffBinning_LOGY",
        # "ElProbeEta",
        # "ElProbeEta_LOGY",
        ]

    if args.showvars:
        showvars = [ int(s) for s in args.showvars.split(',') ]
        print("List of plotting variables and their indexes in job list:\n")
        for idx, s in enumerate(varlist):
            if not ( len(showvars) == 1 and showvars[0] == -1 ):
                if idx not in showvars: continue
            print("job[{0}]={1}".format(idx,s))
        sys.exit()

    outputdir  = args.outputdir
    print("outputdir: {0}".format(outputdir))
    release  = args.release
    queue    = args.queue
    optstr   = args.optsr
    print("optstr: {0}".format(optstr))
    tag      = args.tag
    if tag:
        tag = "_" + tag
    else:
        tag = ""

    # Set the job parameters

    job_params = {
	"release"     :  release,
	"varlist"     :  varlist,
	"subdir"      :  "/coepp/cephfs/mel/mmilesi/ttH/PBS/Plotting/" + outputdir + "_PBS" + tag,  # The path to the submission node (NB: must NOT be under /home : PBS cannot read/write into it!!)
	"outputpath"  :  "/coepp/cephfs/mel/mmilesi/ttH/PlotVault/" + outputdir,
        "optstr"      :  optstr,
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
    # where each job will consist of an array of subjobs (one for each plotting variable).

    jobs = create_jobs(job_params)

    # Finally, execute the PBS script

    if not args.dry:
        submit_jobs(jobs)
