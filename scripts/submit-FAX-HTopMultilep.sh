#!/bin/bash 

# ----------------
# 20.1.8.2 - 25 ns
# ----------------

inDS=mc15_13TeV.341270.aMcAtNloHerwigppEvtGen_UEEE5_CTEQ6L1_CT10ME_ttH125_semilep.merge.DAOD_HIGG8D1.e4277_s2608_s2183_r6869_r6282_p2434/

# ------------------------------------------------------------------------------------

configpath=HTopMultilepAnalysis/scripts/jobOptions_HTopMultilep.py
current_time="$(date +'%d-%m-%Y-%T')"
outdir=output_FAX_DxAOD-2015-13TeV_${current_time}
nevents=0

echo ""
echo "Input DS :"
echo ""
echo ${inDS}
echo ""
echo "Configuring job with :"
echo ""
echo ${configpath}
echo ""
echo "Output will be stored in :"
echo ""
echo ${outdir}
echo ""

xAH_run.py --files ${inDS} --config ${configpath} --submitDir ${outdir} --inputDQ2 --nevents ${nevents} direct
