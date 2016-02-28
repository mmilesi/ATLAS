#!/bin/bash

username=mmilesi

prodtag=027_TruthTP_Retry2_DxAOD
#028_TruthTP_DxAOD
#025a_TruthTP_DxAOD
#CFChallenge_v04_DxAOD
#024a_DxAOD
#CFChallenge_v03_DxAOD

infilepath="HTopMultilepAnalysis/doc/list-grid-DxAOD-2015-13TeV.txt"

configpath="$ROOTCOREBIN/user_scripts/HTopMultilepAnalysis/jobOptions_HTopMultilep.py"

current_time="$(date +'%d-%m-%Y-%T')"
outdir=output_grid_DxAOD-2015-13TeV_${current_time}

#inSE=INFN-T1_DATADISK,MWT2_DATADISK

destSE=AUSTRALIA-ATLAS_LOCALGROUPDISK
exclSE=ANALY_IHEP,ANALY_JINR,ANALY_IHEP_GLEXEC,ANALY_RRC-KI-HPC,ANALY_RRC-KI-T1,IHEP_MCORE,IHEP_PROD,RRC-KI-HPC2,RRC-KI-T1,RRC-KI-T1_MCORE,RRC-KI-T1_TEST

gridDSname="user.${username}.HTopMultilep.${prodtag}.%in:name[2]%.%in:name[3]%"

xAH_run.py -vv --files ${infilepath} --config ${configpath} --inputList --inputDQ2 --submitDir ${outdir} prun --optGridMergeOutput=1 --optGridNFilesPerJob=1.0 --optGridDestSE=${destSE} --optGridOutputSampleName=${gridDSname} --optGridExcludedSite=${exclSE} #--optGridSite=${inSE}

