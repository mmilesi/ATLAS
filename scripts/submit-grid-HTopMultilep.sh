#!/bin/bash 

username=mmilesi

prodtag=021_testLeak1_DxAOD

infilepath=HTopMultilepAnalysis/doc/list-grid-DxAOD-2015-13TeV.txt

configpath=HTopMultilepAnalysis/scripts/jobOptions_HTopMultilep.py

current_time="$(date +'%d-%m-%Y-%T')"
outdir=output_grid_DxAOD-2015-13TeV_${current_time}

destSE=AUSTRALIA-ATLAS_LOCALGROUPDISK

#exclSE=

gridDSname=user.${username}.HTopMultilep.${prodtag}.%in:name[2]%.%in:name[3]%

xAH_run.py --files ${infilepath} --config ${configpath} --inputList --submitDir ${outdir} prun --optGridMergeOutput=1 --optGridNFilesPerJob=1.0 --optGridDestSE=${destSE} --optGridOutputSampleName=${gridDSname} #--optGridExcludedSite=${exclSE}
