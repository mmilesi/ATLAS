#!/usr/bin/env python

from ROOT import ROOT, TFile

baseline   = False
tm_eff     = False
not_tm_eff = True

def copyMuon():

  # Real muon eff from T&P
  #
  if baseline:
    path_r_mu = "./PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIsoTagPt30_ForceProbeToBeFake/"
  elif tm_eff:
    path_r_mu = "./PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIsoTagPt30_ForceProbeToBeFake_TRIGMATCH_EFF/"
  elif not_tm_eff:
    path_r_mu = "./PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIsoTagPt30_ForceProbeToBeFake_NOT_TRIGMATCH_EFF/"

  file_r_mu = TFile(path_r_mu + "LeptonEfficiencies.root")

  hist_r_mu = file_r_mu.Get("Real_Mu_Pt_Efficiency_expectedbkg")

  hist_r_mu.SetDirectory(0)

  # Fake muon eff from likelihood
  #
  if baseline:
    path_f_mu = "./PLOTS_25ns_v24/MMClosure_v24_LikelihoodFit/OutputPlots_MMClosureRates_LHFit_NoCorr_INCLUSIVE_FLAV_DLT_25ns_v24/"
  elif tm_eff:
    path_f_mu = "./PLOTS_25ns_v24/MMClosure_v24_LikelihoodFit/OutputPlots_MMClosureRates_LHFit_NoCorr_INCLUSIVE_FLAV_DLT_25ns_v24_TRIGMATCH_EFF/"
  elif not_tm_eff:
    path_f_mu = "./PLOTS_25ns_v24/MMClosure_v24_LikelihoodFit/OutputPlots_MMClosureRates_LHFit_NoCorr_INCLUSIVE_FLAV_DLT_25ns_v24_NOT_TRIGMATCH_EFF/"

  file_f_mu = TFile(path_f_mu + "LeptonEfficiencies_LH/LH_mumu/LH_efficiencies_fake_mu_mumu.root")

  hist_f_mu = file_f_mu.Get("f_hist")

  hist_f_mu.SetName("Fake_Mu_Pt_Efficiency_expectedbkg")

  hist_f_mu.SetDirectory(0)

  return hist_r_mu, hist_f_mu

def copyElectron():

  # Real electron eff from T&P
  #
  if baseline:
    path_r_el = "./PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIsoTagPt30_ForceProbeToBeFake/"
  elif tm_eff:
    path_r_el = "./PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIsoTagPt30_ForceProbeToBeFake_TRIGMATCH_EFF/"
  elif not_tm_eff:
    path_r_el = "./PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIsoTagPt30_ForceProbeToBeFake_NOT_TRIGMATCH_EFF/"

  file_r_el = TFile(path_r_el + "LeptonEfficiencies.root")

  hist_r_el = file_r_el.Get("Real_El_Pt_Efficiency_expectedbkg")

  hist_r_el.SetDirectory(0)

  # Fake electron eff from T&P
  #
  hist_f_el = file_r_el.Get("Fake_El_Pt_Efficiency_expectedbkg")

  hist_f_el.SetDirectory(0)

  return hist_r_el, hist_f_el


if __name__ == "__main__":

  append = ""
  if tm_eff:
    append = "_TRIGMATCH_EFF"
  elif not_tm_eff:
    append = "_NOT_TRIGMATCH_EFF"

  f = TFile("./PLOTS_25ns_v24/CombinedEfficiencies" + append + "/LeptonEfficiencies.root","RECREATE")

  h_r_mu, h_f_mu = copyMuon()
  h_r_el, h_f_el = copyElectron()

  f.cd()
  h_r_mu.Write()
  h_f_mu.Write()
  h_r_el.Write()
  h_f_el.Write()

  f.Close()
