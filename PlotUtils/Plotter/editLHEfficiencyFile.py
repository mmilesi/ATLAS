#!/usr/bin/env python

from ROOT import ROOT, TFile

def copyMuon():

  path_LH = "../OutputPlots_MMRates_LHFit_25ns_v15/LHRates_25ns_v15/LH_mumu/"
  file_r_LH = TFile(path_LH + "LH_efficiencies_real_mu_mumu.root")
  file_f_LH = TFile(path_LH + "LH_efficiencies_fake_mu_mumu.root")
  
  hist_r_LH = file_r_LH.Get("r_hist") 
  hist_f_LH = file_f_LH.Get("f_hist") 
  
  hist_r_LH.SetName("Mu_ProbePt_Real_Efficiency_observed")
  hist_f_LH.SetName("Mu_ProbePt_Fake_Efficiency_observed")

  hist_r_LH.SetDirectory(0)
  hist_f_LH.SetDirectory(0)
  
  return hist_r_LH, hist_f_LH

def copyElectron():

  path_LH = "../OutputPlots_MMRates_LHFit_25ns_v15/LHRates_25ns_v15/LH_inclusive/"
  file_r_LH = TFile(path_LH + "LH_efficiencies_real_el_incl.root")
  file_f_LH = TFile(path_LH + "LH_efficiencies_fake_el_incl.root")
  
  hist_r_LH = file_r_LH.Get("r_hist") 
  hist_f_LH = file_f_LH.Get("f_hist") 
  
  hist_r_LH.SetName("El_ProbePt_Real_Efficiency_observed")
  hist_f_LH.SetName("El_ProbePt_Fake_Efficiency_observed")

  hist_r_LH.SetDirectory(0)
  hist_f_LH.SetDirectory(0)
  
  return hist_r_LH, hist_f_LH

if __name__ == "__main__":

  f = TFile("../OutputPlots_MMRates_LHFit_25ns_v15/Rates.root","RECREATE")

  h_r_mu, h_f_mu = copyMuon()
  h_r_el, h_f_el = copyElectron()

  f.cd()
  h_r_mu.Write()
  h_f_mu.Write() 
  h_r_el.Write()
  h_f_el.Write() 
  
  f.Close()
