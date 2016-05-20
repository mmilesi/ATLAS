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

oldpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v7/25ns_v7_Data_Original/"
newpath = "/afs/cern.ch/user/m/mmilesi/work/private/ttH/MiniNTup/25ns_v7/25ns_v7_Data_MM_WEIGHTED/"

# Set only if MM weight branch does not exist yet
#
addMM = True

treename = 'physics'
nentries = '100' #ALL

if not os.path.exists(newpath):
    os.makedirs(newpath)

gROOT.LoadMacro("modifyttree_MM.cxx+g")
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

        command_line = None
        if addMM:
            command_line = 'modifyttree_AddMM(\"'+infile+'\",\"'+nentries+'\",\"'+treename+'\",\"'+outfile+'\")'
        else:
            command_line = 'modifyttree_MM(\"'+infile+'\",\"'+nentries+'\",\"'+treename+'\",\"'+outfile+'\")'

        print command_line
        gROOT.ProcessLine(command_line);

