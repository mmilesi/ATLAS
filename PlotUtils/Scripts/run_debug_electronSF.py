from math import sqrt, pow
from ROOT import TCanvas, TFile, TGraph, TGraphErrors, TColor, TAttFill, TStyle, TLegend, TH1D, gROOT, TF1, TTree
import array
import copy
import os, subprocess, sys, time, shlex
sys.path.append(os.path.abspath(os.path.curdir))

gROOT.SetBatch(True)

oldpath  = '/data/mmilesi/ttH/MergedDatasets/test_diboson/'
treename = 'physics'
nentries = 'ALL'

gROOT.LoadMacro("debug_electronSF.cxx+g")
group_list = os.listdir(oldpath)
group_list = group_list[:]
for group in group_list:
    if not os.path.isdir(oldpath+group): 
        continue
    sample_list = os.listdir(oldpath+group+'/')
    for sample in sample_list:
        print group+'/'+sample
        infile=oldpath+group+'/'+sample
        command_line = 'debug_electronSF(\"'+infile+'\",\"'+nentries+'\",\"'+treename+'\")'
        print command_line
        gROOT.ProcessLine(command_line);
