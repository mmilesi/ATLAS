#!/bin/bash

# -----------------------------
# full xAOD - mc15 13TeV - 25ns
# -----------------------------

#inDS="mc15_13TeV.410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.merge.AOD.e3698_s2608_s2183_r6765_r6282"
#infilepath="/data/mmilesi/HTopMultileptonsTestSamples/MC15/r6282/${inDS}/*root*"


# ---------------------
# 20.1.9.3 - 25 ns - MC
# ---------------------

#inDS="mc15_13TeV.410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.merge.DAOD_HIGG8D1.e3698_s2608_s2183_r7267_r6282_p2501"
#inDS="mc15_13TeV.341270.aMcAtNloHerwigppEvtGen_UEEE5_CTEQ6L1_CT10ME_ttH125_semilep.merge.DAOD_HIGG8D1.e4277_s2608_s2183_r6869_r6282_p2501"
#infilepath="/data/mmilesi/HTopMultileptonsTestSamples/HIGG8D1_20.1.9.3/${inDS}/*root*"

# -------------------------------
# 20.1.9.8 - 25 ns - MC (skimmed)
# -------------------------------

#inDS="mc15_13TeV.410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.merge.DAOD_HIGG8D1.e3698_s2608_s2183_r7267_r6282_p2559"
#inDS="mc15_13TeV.341270.aMcAtNloHerwigppEvtGen_UEEE5_CTEQ6L1_CT10ME_ttH125_semilep.merge.DAOD_HIGG8D1.e4277_s2608_s2183_r6869_r6282_p2559"
#infilepath="/data/mmilesi/HTopMultileptonsTestSamples/HIGG8D1_20.1.9.8/${inDS}/*root*"

# -----------------------
# 20.1.9.8 - 25 ns - data
# -----------------------

inDS="data15_13TeV.00281070.physics_Main.merge.DAOD_HIGG8D1.f629_m1504_p2559"
infilepath="/data/mmilesi/HTopMultileptonsTestSamples/HIGG8D1_20.1.9.8/${inDS}/*root*"

# -----------------------
# 20.1.6.3 - 25 ns - data
# -----------------------
#
#inDS="data15_13TeV.00276329.physics_Main.merge.DAOD_HIGG8D1.f620_m1480_p2411/DAOD_HIGG8D1.06323611._000001.pool.root.1"
#infilepath="/data/mmilesi/HTopMultileptonsTestSamples/HIGG8D1_20.1.6.3/${inDS}"

# ------------------------------------------------------------------------------------

# tokenize inDS using '.' as separator
#
tokens=(${inDS//./ })

configpath="$ROOTCOREBIN/user_scripts/HTopMultilepAnalysis/jobOptions_HTopMultilep.py"
current_time="$(date +'%d-%m-%Y-%T')"
outdir=output_local_DxAOD-2015-13TeV_${tokens[2]}_${current_time}
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

xAH_run.py -vv --files ${infilepath} --config ${configpath} --submitDir ${outdir} --nevents ${nevents} direct
