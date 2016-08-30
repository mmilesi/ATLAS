#!/bin/bash

# -------------------------------
# 20.7.6.7 - 25 ns - MC (skimmed)
# -------------------------------

inDS="mc15_13TeV.410155.aMcAtNloPythia8EvtGen_MEN30NLO_A14N23LO_ttW.merge.DAOD_HIGG8D1.e5070_s2726_r7772_r7676_p2719"
#infilepath="/data/mmilesi/HTopMultileptonsTestSamples/HIGG8D1_20.7.6.7/${inDS}/*root*"
infilepath="/afs/cern.ch/user/m/mmilesi/work/private/HTopMultileptonsTestSamples/HIGG8D1_20.7.6.7/${inDS}/*root*"

# -----------------------
# 20.7.6.5 - 25 ns - data
# -----------------------

#inDS="data16_13TeV.00303726.physics_Main.merge.DAOD_HIGG8D1.f716_m1620_p2689"
#infilepath="/data/mmilesi/HTopMultileptonsTestSamples/HIGG8D1_20.7.6.5/${inDS}/*root*"
#infilepath="/afs/cern.ch/user/m/mmilesi/work/private/HTopMultileptonsTestSamples/HIGG8D1_20.7.6.5/${inDS}/*root*"

# ------------------------------------------------------------------------------------

# tokenize inDS using '.' as separator
#
tokens=(${inDS//./ })

configpath="$ROOTCOREBIN/user_scripts/HTopMultilepAnalysis/jobOptions_HTopMultilep.py"

current_time="$(date +'%d-%m-%Y-%T')"
outdir=output_local_DxAOD-2016-13TeV_${tokens[2]}_${current_time}
nevents=1000

echo ""
echo "Input file path :"
echo ""
echo ${infilepath}
echo ""
echo "Configuring job with :"
echo ""
echo ${configpath}
echo ""
echo "Output will be stored in :"
echo ""
echo ${outdir}
echo ""

xAH_run.py -vv --files ${infilepath} --config ${configpath} --submitDir ${outdir} --nevents ${nevents} --force direct
