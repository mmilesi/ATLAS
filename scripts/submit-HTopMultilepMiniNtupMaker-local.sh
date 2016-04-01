#!/bin/bash

# -----------------------------
# 25ns_v6
# -----------------------------

inDS="group.ntup.25ns_v6_410000"
infilepath="/coepp/cephfs/mel/mmilesi/ttH/multileptons_ntuple_run2/25ns_v6/Nominal/410000/410000.root"

# ------------------------------------------------------------------------------------

# tokenize inDS using '.' as separator
#
tokens=(${inDS//./ })

configpath="$ROOTCOREBIN/user_scripts/HTopMultilepAnalysis/jobOptions_HTopMultilepMiniNtupMaker.py"
current_time="$(date +'%d-%m-%Y-%T')"
outdir=output_local_DxAOD-2015-13TeV_${tokens[2]}_${current_time}
nevents=500

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

xAH_run.py -vv --files ${infilepath} --config ${configpath} --treeName "nominal" --submitDir ${outdir} --nevents ${nevents} direct
