#!/usr/bin/python
#
# ***********************************************************************
# Steering macro for modifytree_MMClosure.C:
# loop over the merged samples to change the content of branches in TTree
#
# Authors:
#  Francesco Nuti ( francesco.nuti@cern.ch )
#
# ***********************************************************************

from math import sqrt, pow
from ROOT import TCanvas, TFile, TGraph, TGraphErrors, TColor, TAttFill, TStyle, TLegend, TH1D, gROOT, TF1, TTree
import array
import copy
import os, subprocess, sys, time, shlex
sys.path.append(os.path.abspath(os.path.curdir))

gROOT.SetBatch(True)

oldpath  = '/data/mmilesi/ttH/MergedDatasets/Merged_Melb15_ttH_028_DxAOD_TTbar_MMClosure_TruthTP/'
newpath  = '/data/mmilesi/ttH/MergedDatasets/Merged_Melb15_ttH_028_DxAOD_TTbar_MMClosure_TruthTP_WEIGHTED/'

treename = 'physics'
nentries = 'ALL' #ALL

if not os.path.exists(newpath):
    os.makedirs(newpath)

gROOT.LoadMacro("modifyttree_MMClosure.C+g")
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

        #command_line = 'root -l -q \"modifyttree_MMClosure.C+g(\"'+infile+'\",\"'+treename+'\",\"'+outfile+'\")\"'
        #print command_line
        #args = shlex.split(command_line)
        #print args
        #subprocess.call(args)
        #subprocess.call(['root', '-l', '-q', '"modifyttree_MMClosure.C+g(\"'+infile+'\",\"'+treename+'\",\"'+outfile+'\")"'])
        #f = TFile.Open(path+group+'/'+sample)
        #h_tot = f.Get('TotalEvents')
        #t = f.Get('physics')
