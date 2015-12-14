#!/bin/bash

username=mmilesi

prodtag=023b_DxAOD

infilepath="HTopMultilepAnalysis/doc/list-grid-DxAOD-2015-13TeV.txt"

configpath="$ROOTCOREBIN/user_scripts/HTopMultilepAnalysis/jobOptions_HTopMultilep.py"

current_time="$(date +'%d-%m-%Y-%T')"
outdir=output_grid_DxAOD-2015-13TeV_${current_time}

destSE=AUSTRALIA-ATLAS_LOCALGROUPDISK
exclSE=ANALY_IHEP

gridDSname="user.${username}.HTopMultilep.${prodtag}.%in:name[2]%.%in:name[3]%"

xAH_run.py -vv --files ${infilepath} --config ${configpath} --inputList --inputDQ2 --submitDir ${outdir} prun --optGridMergeOutput=1 --optGridNFilesPerJob=1.0 --optGridDestSE=${destSE} --optGridOutputSampleName=${gridDSname} --optGridExcludedSite=${exclSE}

