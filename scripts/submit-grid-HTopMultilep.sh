#!/bin/bash

username=mmilesi

prodtag=024a_DxAOD
#CFChallenge_v04_DxAOD
#024a_DxAOD
#CFChallenge_v03_DxAOD

infilepath="HTopMultilepAnalysis/doc/list-grid-DxAOD-2015-13TeV.txt"

configpath="$ROOTCOREBIN/user_scripts/HTopMultilepAnalysis/jobOptions_HTopMultilep.py"

current_time="$(date +'%d-%m-%Y-%T')"
outdir=output_grid_DxAOD-2015-13TeV_${current_time}

#inSE=INFN-T1_DATADISK,MWT2_DATADISK

destSE=AUSTRALIA-ATLAS_LOCALGROUPDISK
exclSE=ANALY_IHEP,ANALY_DESY-HH,ANALY_FZK,ANALY_FZK_SHORT,ANALY_GRIF-IRFU,ANALY_MWT2_SL6,ANALY_INFN-MILANO-ATLASC,ANALY_GLASGOW_SL6,ANALY_INFN-T1,ANALY_LPSC,ANALY_SLAC_SHORT_1HR,ANALY_MCGILL,ANALY_TRIUMF,ANALY_AGLT2_SL6,ANALY_BNL_SHORT,ANALY_SWT2_CPB,ANALY_MANC_SL6_SHORT,ANALY_SiGNET,ANALY_INFN-FRASCATI,ANALY_IFIC,ANALY_VICTORIA,ANALY_INFN-NAPOLI-RECAS

gridDSname="user.${username}.HTopMultilep.${prodtag}.%in:name[2]%.%in:name[3]%"

xAH_run.py -vv --files ${infilepath} --config ${configpath} --inputList --inputDQ2 --submitDir ${outdir} prun --optGridMergeOutput=1 --optGridNFilesPerJob=1.0 --optGridDestSE=${destSE} --optGridOutputSampleName=${gridDSname} --optGridExcludedSite=${exclSE} #--optGridSite=${inSE}

