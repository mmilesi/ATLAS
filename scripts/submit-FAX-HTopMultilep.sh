#!/bin/bash

# when using DQ2 CONTAINERS, put a trailing '/'
#
#inDS="mc15_13TeV.341270.aMcAtNloHerwigppEvtGen_UEEE5_CTEQ6L1_CT10ME_ttH125_semilep.merge.DAOD_HIGG8D1.e4277_s2608_s2183_r6869_r6282_p2434/"

# when using DQ2 DATASETs, do *not* put a trailing '/'
#
inDS="mc15_13TeV.361505.MadGraphPythia8EvtGen_A14NNPDF23LO_Zmumu_Np0.merge.DAOD_HIGG8D1.e3898_s2608_s2183_r6869_r6282_p2438_tid06808706_00"

# ------------------------------------------------------------------------------------

configpath="$ROOTCOREBIN/user_scripts/HTopMultilepAnalysis/jobOptions_HTopMultilep.py"
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

#xAH_run.py -vv --files ${inDS} --config ${configpath} --submitDir ${outdir} --inputDQ2 --nevents ${nevents} direct
#
# specify --inputTag if you wish to run on specific file(s) in the input dataset
#
xAH_run.py -vv --files ${inDS} --config ${configpath} --submitDir ${outdir} --inputDQ2 --inputTag 'DAOD_HIGG8D1.06808706._000002.pool.root.1' --nevents ${nevents} direct
