#!/usr/bin/env python

from ROOT import ROOT, TFile

def copyMuon():

  path_LH = "../MMClosure_v21_RightDLTTrigMatching_LikelihoodFit/OutputPlots_MMClosureRates_LHFit_NoCorr_SFmuSFel_DLT_25ns_v21/LeptonEfficiencies_LH/LH_mumu/"
  #path_LH = "../MMClosure_v21_RightDLTTrigMatching_LikelihoodFit/OutputPlots_MMClosureRates_LHFit_NoCorr_SFmuSFel_SLT_25ns_v21/LeptonEfficiencies_LH/LH_mumu/"
  
  file_r_LH = TFile(path_LH + "LH_efficiencies_real_mu_mumu.root")
  file_f_LH = TFile(path_LH + "LH_efficiencies_fake_mu_mumu.root")

  hist_r_LH = file_r_LH.Get("r_hist")
  hist_f_LH = file_f_LH.Get("f_hist")

  hist_r_LH.SetName("Real_Mu_Pt_Efficiency_expectedbkg")
  hist_f_LH.SetName("Fake_Mu_Pt_Efficiency_expectedbkg")

  hist_r_LH.SetDirectory(0)
  hist_f_LH.SetDirectory(0)

  return hist_r_LH, hist_f_LH

def copyElectron():

  #path_LH = "../MMClosure_v21_RightDLTTrigMatching_LikelihoodFit/OutputPlots_MMClosureRates_LHFit_NoCorr_SFmuSFel_DLT_25ns_v21/LeptonEfficiencies_LH/LH_elel/"
  path_LH = "../MMClosure_v21_RightDLTTrigMatching_LikelihoodFit/OutputPlots_MMClosureRates_LHFit_NoCorr_SFmuINCLel_DLT_25ns_v21/LeptonEfficiencies_LH/LH_incl/"
  
  #path_LH = "../MMClosure_v21_RightDLTTrigMatching_LikelihoodFit/OutputPlots_MMClosureRates_LHFit_NoCorr_SFmuSFel_SLT_25ns_v21/LeptonEfficiencies_LH/LH_elel/"  
  
  #file_r_LH = TFile(path_LH + "LH_efficiencies_real_el_elel.root")
  #file_f_LH = TFile(path_LH + "LH_efficiencies_fake_el_elel.root")
  
  file_r_LH = TFile(path_LH + "LH_efficiencies_real_el_incl.root")
  file_f_LH = TFile(path_LH + "LH_efficiencies_fake_el_incl.root")

  hist_r_LH = file_r_LH.Get("r_hist")
  hist_f_LH = file_f_LH.Get("f_hist")

  hist_r_LH.SetName("Real_El_Pt_Efficiency_expectedbkg")
  hist_f_LH.SetName("Fake_El_Pt_Efficiency_expectedbkg")

  hist_r_LH.SetDirectory(0)
  hist_f_LH.SetDirectory(0)

  return hist_r_LH, hist_f_LH

if __name__ == "__main__":

  #f = TFile("../MMClosure_v21_RightDLTTrigMatching_LikelihoodFit/OutputPlots_MMClosureRates_LHFit_NoCorr_SFmuSFel_DLT_25ns_v21/LeptonEfficiencies_LH/LeptonEfficiencies_Files/SF_mu_SF_el/LeptonEfficiencies_LH.root","RECREATE")
  
  f = TFile("../MMClosure_v21_RightDLTTrigMatching_LikelihoodFit/OutputPlots_MMClosureRates_LHFit_NoCorr_SFmuINCLel_DLT_25ns_v21/LeptonEfficiencies_LH/LeptonEfficiencies_Files/SF_mu_INCL_el/LeptonEfficiencies_LH.root","RECREATE")
  
  #f = TFile("../MMClosure_v21_RightDLTTrigMatching_LikelihoodFit/OutputPlots_MMClosureRates_LHFit_NoCorr_SFmuSFel_SLT_25ns_v21/LeptonEfficiencies_LH/LeptonEfficiencies_Files/SF_mu_SF_el/LeptonEfficiencies_LH.root","RECREATE")

  h_r_mu, h_f_mu = copyMuon()
  h_r_el, h_f_el = copyElectron()

  f.cd()
  h_r_mu.Write()
  h_f_mu.Write()
  h_r_el.Write()
  h_f_el.Write()

  f.Close()
