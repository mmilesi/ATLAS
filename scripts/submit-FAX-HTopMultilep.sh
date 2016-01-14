#!/bin/bash

# when using DQ2 CONTAINERS, put a trailing '/'
#
#inDS="mc15_13TeV.341270.aMcAtNloHerwigppEvtGen_UEEE5_CTEQ6L1_CT10ME_ttH125_semilep.merge.DAOD_HIGG8D1.e4277_s2608_s2183_r6869_r6282_p2434/"
#inDS="data15_13TeV.00280464.physics_Main.merge.DAOD_HIGG8D1.f629_m1504_p2432/"
inDS="mc15_13TeV.341270.aMcAtNloHerwigppEvtGen_UEEE5_CTEQ6L1_CT10ME_ttH125_semilep.merge.DAOD_HIGG8D1.e4277_s2608_s2183_r6869_r6282_p2501/"

# when using DQ2 DATASETS, do *not* put a trailing '/'
#
#inDS="mc15_13TeV.410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.merge.DAOD_HIGG8D1.e3698_s2608_s2183_r6765_r6282_p2434_tid06670680_00"

# ------------------------------------------------------------------------------------

# tokenize inDS using '.' as separator
#
tokens=(${inDS//./ })

configpath="$ROOTCOREBIN/user_scripts/HTopMultilepAnalysis/jobOptions_HTopMultilep.py"
current_time="$(date +'%d-%m-%Y-%T')"
outdir=output_FAX_DxAOD-2015-13TeV_${tokens[2]}_${current_time}
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

xAH_run.py -vv --files ${inDS} --config ${configpath} --submitDir ${outdir} --inputDQ2 --nevents ${nevents} direct
#
# specify --inputTag if you wish to run on specific file(s) in the input dataset
#
#xAH_run.py -vv --files ${inDS} --config ${configpath} --submitDir ${outdir} --inputDQ2 --inputTag 'DAOD_HIGG8D1.06808706._000002.pool.root.1' --nevents ${nevents} direct
