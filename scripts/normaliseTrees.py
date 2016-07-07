#!/usr/bin/env python

import os, sys

from ROOT import TFile, TH1, TH1D, TTree

def applyWeight(filepath, s, treename='physics', isdata=False):

  if isdata:
    print("Sample is DATA, no need to apply normalisation weight!")
    return

  f = TFile.Open(filepath,"UPDATE")
  t = f.Get(treename)

  if not t:
    print ("WARNING: tree {0} cannot be found during tree registration. Sample will not be weighted...".format(treename))
    return 

  print("Weighting tree {0} w/ Xsec weight - sample ID[{1}]".format(filepath,s["ID"]))
  weight = float(s["xsection"]) * float(s["efficiency"]) * float(s["kfactor"]) * 1e3 # to get the weight in fb (the Xsec is in pb)
  h = f.Get("TotalEventsW")
  if not h:
     print ("WARNING: histogram named TotalEventsW in file {1} couldn't be found!".format(filepath))
     return False
  weight /= h.GetBinContent(2)
  t.SetWeight(weight)
  t.Write(t.GetName(),t.kOverwrite)
  print("===> w = {0}".format(weight))
      
  f.Close()
