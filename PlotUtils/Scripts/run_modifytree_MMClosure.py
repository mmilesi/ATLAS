#!/usr/bin/env python

""" run_modifytree_MMClosure.py: simple script to execute a ROOT macro by calling the CINT interpreter with TROOT::ProcessLine() """

__author__     = "Marco Milesi, Francesco Nuti"
__email__      = "marco.milesi@cern.ch, francesco.nuti@cern.ch"
__maintainer__ = "Marco Milesi"

import os, subprocess, sys, time, shlex, copy

sys.path.append(os.path.abspath(os.path.curdir))

from ROOT import gROOT

gROOT.SetBatch(True)

#oldpath  = '/data/mmilesi/ttH/MergedDatasets/Merged_Melb15_ttH_027_DxAOD_TTbar_MMClosure_TruthTP/'
#newpath  = '/data/mmilesi/ttH/MergedDatasets/Merged_Melb15_ttH_027_DxAOD_TTbar_MMClosure_TruthTP_WEIGHTED/'

#oldpath  = '/data/mmilesi/ttH/MergedDatasets/Merged_Melb15_ttH_028_DxAOD_TTbar_MMClosure_TruthTP/'
#newpath  = '/data/mmilesi/ttH/MergedDatasets/Merged_Melb15_ttH_028_DxAOD_TTbar_MMClosure_TruthTP_WEIGHTED/'
#newpath  = '/data/mmilesi/ttH/MergedDatasets/Merged_Melb15_ttH_028_DxAOD_TTbar_MMClosure_TruthTP_WEIGHTED_RatesHighNJet/'
#newpath  = '/data/mmilesi/ttH/MergedDatasets/Merged_Melb15_ttH_028_DxAOD_TTbar_MMClosure_TruthTP_WEIGHTED_new/'

oldpath  = '/data/mmilesi/ttH/MergedDatasets/Merged_v029/ICHEP_BASELINE/Merged_Melb15_ttH_029_Baseline_DxAOD_TTbar_MMClosure_TruthTP/'
#newpath  = '/data/mmilesi/ttH/MergedDatasets/Merged_v029/ICHEP_BASELINE/Merged_Melb15_ttH_029_Baseline_DxAOD_TTbar_MMClosure_TruthTP_WEIGHTED/'
newpath  = '/data/mmilesi/ttH/MergedDatasets/Merged_v029/ICHEP_BASELINE/Merged_Melb15_ttH_029_Baseline_DxAOD_TTbar_MMClosure_TruthTP_WEIGHTED_EffNoLepIso/'

#oldpath  = '/data/mmilesi/ttH/MergedDatasets/Merged_v029/ICHEP_NOLEPISO/Merged_Melb15_ttH_029_NoLepIso_DxAOD_TTbar_MMClosure_TruthTP/'
#newpath  = '/data/mmilesi/ttH/MergedDatasets/Merged_v029/ICHEP_NOLEPISO/Merged_Melb15_ttH_029_NoLepIso_DxAOD_TTbar_MMClosure_TruthTP_WEIGHTED/'

#oldpath  = '/data/mmilesi/ttH/MergedDatasets/Merged_v029/ICHEP_NOLEPIP/Merged_Melb15_ttH_029_NoLepIP_DxAOD_TTbar_MMClosure_TruthTP/'
#newpath  = '/data/mmilesi/ttH/MergedDatasets/Merged_v029/ICHEP_NOLEPIP/Merged_Melb15_ttH_029_NoLepIP_DxAOD_TTbar_MMClosure_TruthTP_WEIGHTED/'


treename = 'physics'
nentries = 'ALL' #ALL

if not os.path.exists(newpath):
    os.makedirs(newpath)

gROOT.LoadMacro("modifyttree_MMClosure.cxx+g")
group_list = os.listdir(oldpath)
group_list = group_list[:]
for group in group_list:
    if not os.path.isdir(oldpath+group):
        continue
    if not os.path.exists(newpath+group):
        os.makedirs(newpath+group)
    sample_list = os.listdir(oldpath+group+'/')
    for sample in sample_list:
        print group+'/'+sample
        infile=oldpath+group+'/'+sample
        outfile=newpath+group+'/'+sample

        command_line = 'modifyttree_MMClosure(\"'+infile+'\",\"'+nentries+'\",\"'+treename+'\",\"'+outfile+'\")'
        print command_line
        gROOT.ProcessLine(command_line);

        #command_line = 'root -l -q \"modifyttree_MMClosure.cxx+g(\"'+infile+'\",\"'+treename+'\",\"'+outfile+'\")\"'
        #print command_line
        #args = shlex.split(command_line)
        #print args
        #subprocess.call(args)
        #subprocess.call(['root', '-l', '-q', '"modifyttree_MMClosure.cxx+g(\"'+infile+'\",\"'+treename+'\",\"'+outfile+'\")"'])
        #f = TFile.Open(path+group+'/'+sample)
        #h_tot = f.Get('TotalEvents')
        #t = f.Get('physics')
