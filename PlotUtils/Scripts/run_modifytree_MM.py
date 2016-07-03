#!/usr/bin/env python

""" run_modifytree_MM.py: simple script to execute a ROOT macro by calling the CINT interpreter with TROOT::ProcessLine() """

__author__     = "Marco Milesi, Francesco Nuti"
__email__      = "marco.milesi@cern.ch, francesco.nuti@cern.ch"
__maintainer__ = "Marco Milesi"

import os, subprocess, sys, time, shlex, copy

sys.path.append(os.path.abspath(os.path.curdir))

from ROOT import gROOT

gROOT.SetBatch(True)

#oldpath  = '/data/mmilesi/ttH/MergedDatasets/Merged_v029/ICHEP_NOLEPISO/Merged_Melb15_ttH_029_NoLepIso_DxAOD_DATA/'
#newpath  = '/data/mmilesi/ttH/MergedDatasets/Merged_v029/ICHEP_NOLEPISO/Merged_Melb15_ttH_029_NoLepIso_DxAOD_DATA_MM_WEIGHTED/'

# -------------------------------------------------------------------------------------------------------------

#oldpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v7/25ns_v7_Data_Original/"
#newpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v7/25ns_v7_Data_MM_WEIGHTED/"
#rr_dir = "/afs/cern.ch/user/m/mmilesi/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMRates_25ns_v7_FinalSelection_NominalBinning/Rates_YesSub_LHInput"

#oldpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v7/25ns_v7_410000_Original/"
#newpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v7/25ns_v7_410000_MM_WEIGHTED/"
#rr_dir = "/afs/cern.ch/user/m/mmilesi/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMClosureRates_25ns_v7_FinalSelection_NominalBinning/Rates_YesSub_LHInput"

# -------------------------------------------------------------------------------------------------------------

#oldpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v7/25ns_v7_Data_Original/"
#newpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v7/25ns_v7_Data_MM_WEIGHTED_AVGMUFAKE/"
#rr_dir = "/afs/cern.ch/user/m/mmilesi/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMRates_25ns_v7_FinalSelection_NominalBinning/Rates_YesSub_AvgMuFake"

#oldpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v7/25ns_v7_Data_Original_DDQMisID/"
#newpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v7/25ns_v7_Data_DDQMisID_MM_WEIGHTED_AVGMUFAKE/"
#rr_dir  = "/afs/cern.ch/user/m/mmilesi/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMRates_25ns_v7_FinalSelection_DDQmisID"

#oldpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v7/25ns_v7_410000_Original/"
#newpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v7/25ns_v7_410000_MM_WEIGHTED_AVGMUFAKE/"
#rr_dir = "/afs/cern.ch/user/m/mmilesi/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMClosureRates_25ns_v7_FinalSelection_NominalBinning/Rates_YesSub_AvgMuFake"

# -------------------------------------------------------------------------------------------------------------

#oldpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v14/25ns_v14_Direct_410000_Original/"
#oldpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v14/25ns_v14_Direct_DLT_410000_Original/"

#newpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v14/25ns_v14_Direct_410000_MM_WEIGHTED_EtaPt/"
#newpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v14/25ns_v14_Direct_410000_MM_WEIGHTED_Pt/"
#rr_dir  = "/afs/cern.ch/user/m/mmilesi/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMClosureRates_25ns_v14"

#newpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v14/25ns_v14_Direct_410000_MM_WEIGHTED_Pt_2015/"
#newpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v14/25ns_v14_Direct_410000_MM_WEIGHTED_Pt_2015_ForceFinalRealBin/"
#newpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v14/25ns_v14_Direct_410000_MM_WEIGHTED_Pt_2015_FinerBinningMuFakeEff/"

#newpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v14/25ns_v14_Direct_DLT_410000_MM_WEIGHTED_Pt_2015/"

#rr_dir  = "/afs/cern.ch/user/m/mmilesi/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMClosureRates_25ns_v14_2015"
#rr_dir  = "/afs/cern.ch/user/m/mmilesi/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMClosureRates_25ns_v14_DLT_2015"

# -------------------------------------------------------------------------------------------------------------

#oldpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v15/25ns_v15_Direct_Data_Original/"
#newpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v15/25ns_v15_Direct_Data_MM_WEIGHTED_Pt_TEST/"
#rr_dir  = "/afs/cern.ch/user/m/mmilesi/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMRates_25ns_v15"

#oldpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v15/25ns_v15_Direct_410000_Original/"
#newpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v15/25ns_v15_Direct_410000_MM_WEIGHTED_Pt/"
#rr_dir  = "/afs/cern.ch/user/m/mmilesi/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMClosureRates_25ns_v15/Rates_NominalBinning"
#newpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v15/25ns_v15_Direct_410000_MM_WEIGHTED_Pt/"

#oldpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v15/25ns_v15_Direct_410000_Original/"
#newpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v15/25ns_v15_Direct_410000_MM_WEIGHTED_Pt_FinerBinningMuFakeEff/"
#rr_dir  = "/afs/cern.ch/user/m/mmilesi/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMClosureRates_25ns_v15/Rates_FineMuonFakeEffBinning"

#oldpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v15/25ns_v15_Direct_Data_Original/"
#newpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v15/25ns_v15_Direct_Data_MM_WEIGHTED_Pt_LH/"
#rr_dir  = "/afs/cern.ch/user/m/mmilesi/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMRates_LHFit_25ns_v15"

oldpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v15/25ns_v15_Direct_Data_Original/"
newpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v15/25ns_v15_Direct_Data_MM_WEIGHTED_Pt_ScaledFakeEl/"
rr_dir  = "/afs/cern.ch/user/m/mmilesi/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMRates_25ns_v15"

# -------------------------------------------------------------------------------------------------------------

doMM_PtOnly = True
addMM    = 'YES' # Set to 'YES' if MM weight branch does not exist yet
nentries = 'ALL' #ALL

if not os.path.exists(newpath):
    os.makedirs(newpath)

if doMM_PtOnly:
  gROOT.LoadMacro("modifyttree_MM_PtOnly.cxx+g")
else:
  gROOT.LoadMacro("modifyttree_MM.cxx+g")

group_list = os.listdir(oldpath)
group_list = group_list[:]

do_closure = None

for group in group_list:
    if not os.path.isdir(oldpath+group):
        continue
    if not os.path.exists(newpath+group):
        os.makedirs(newpath+group)
    sample_list = os.listdir(oldpath+group+'/')
    for sample in sample_list:
        if "hist-" in sample:
            continue
        print group+'/'+sample
        infile=oldpath+group+'/'+sample
        outfile=newpath+group+'/'+sample

        if "physics_Main" in sample: do_closure = "NO"
        else:                        do_closure = "YES"

        command_line = None
        if doMM_PtOnly:
	  command_line = 'modifyttree_MM_PtOnly(\"'+infile+'\",\"'+outfile+'\",\"'+addMM+'\",\"'+rr_dir+'\",\"'+do_closure+'\",\"'+nentries+'\")'
        else:
          command_line = 'modifyttree_MM(\"'+infile+'\",\"'+outfile+'\",\"'+addMM+'\",\"'+rr_dir+'\",\"'+do_closure+'\",\"'+nentries+'\")'

        print command_line
        gROOT.ProcessLine(command_line);

