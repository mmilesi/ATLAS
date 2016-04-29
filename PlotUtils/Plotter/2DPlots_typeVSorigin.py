 #!/usr/bin/python

import array
import os
import sys

from array import array

from ROOT import gROOT, gDirectory, gStyle, TH1D, TH2D, TFile, TCanvas, TColor, TLegend, TLatex

gROOT.Reset()
gROOT.LoadMacro("AtlasStyle.C")
from ROOT import SetAtlasStyle
SetAtlasStyle()

gROOT.SetBatch(True)

# Are we using group ntuples?
#
useGroupNTup = True

# Normalise histograms to unity
#
doNorm = True

# Are we using truth cuts?
#
useTruth = False

# for fancy 2-dim histograms!
#
def set_fancy_2D_style():

  icol = 0
  gStyle.SetFrameBorderMode(icol);
  gStyle.SetFrameFillColor(icol);
  gStyle.SetCanvasBorderMode(icol);
  gStyle.SetCanvasColor(icol);
  gStyle.SetPadBorderMode(icol);
  gStyle.SetPadColor(icol);
  gStyle.SetStatColor(icol);
  gStyle.SetOptTitle(0);
  gStyle.SetOptStat(0);
  gStyle.SetOptFit(0);

  ncontours=999

  s = array('d', [0.00, 0.34, 0.61, 0.84, 1.00])
  r = array('d', [0.00, 0.00, 0.87, 1.00, 0.51])
  g = array('d', [0.00, 0.81, 1.00, 0.20, 0.00])
  b = array('d', [0.51, 1.00, 0.12, 0.00, 0.00])

  npoints = len(s)
  TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
  gStyle.SetNumberContours(ncontours)

def plot_2D_FakeCR():

  myfilename = ("/coepp/cephfs/mel/mmilesi/ttH/MergedDatasets/Merged_v030/Merged_Melb15_ttH_030_DxAOD_p2559/tops/410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.root","/coepp/cephfs/mel/mmilesi/ttH/MiniNTup/25ns_v7/25ns_v7/tops/410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.root")[bool(useGroupNTup)]
  #myfilename = ("/coepp/cephfs/mel/mmilesi/ttH/MergedDatasets/Merged_v031/Merged_Melb15_ttH_031-01_DxAOD_p2559/tops/410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.root","TO_BE_ADDED")[bool(useGroupNTup)]
  
  myfile     = TFile(myfilename)
  mytree     = gDirectory.Get("physics")

  gROOT.LoadMacro("$ROOTCOREBIN/user_scripts/HTopMultilepAnalysis/ROOT_TTreeFormulas/largeEtaEvent.cxx+")
  from ROOT import largeEtaEvent

  mytree_weight     = "3.209 * mcEventWeight * weight_pileup * weight_lepton_trig_HTop[0] * weight_lepton_reco_HTop[0] * weight_lepton_iso_HTop[0] * weight_lepton_ID_HTop[0] * weight_lepton_TTVA_HTop[0] * weight_jet_JVT_HTop[0] * weight_jet_MV2c20_SFFix77[0]"
  trig_dec	    = "( passHLT == 1 )"
  trigmatching      = "( ( lep_isTrigMatched[0] == 1 && ( ( lep_flavour[0] == 11 && lep_pt[0] > 25e3 ) || ( lep_flavour[0] == 13 && lep_pt[0] > 21e3 ) ) ) || ( lep_isTrigMatched[1] == 1 && ( ( lep_flavour[1] == 11 && lep_pt[1] > 25e3 ) || ( lep_flavour[1] == 13 && lep_pt[1] > 21e3 ) ) ) )"
  lep_tag_trigmatch = "( lep_tag_isTightSelected[0] == 1 && lep_tag_isTrigMatched[0] == 1 && ( ( lep_tag_flavour[0] == 11 && lep_tag_pt[0] > 25e3 ) || ( lep_tag_flavour[0] == 13 && lep_tag_pt[0] > 21e3 ) ) )"
  el_eta	    = "( nel == 0 || ( nel > 0 && Max$( TMath::Abs(el_caloCluster_eta) ) < 1.37 ) )"
  nbjets	    = "( njets_mv2c20_Fix77 > 0 )"
  njets 	    = "( njets > 1 && njets <= 4 )"
  nleptons	    = "( nlep == 2 && Min$( lep_pt ) > 10e3 )"
  flavour_cut	    = "( 1 )"
  same_sign	    = "( isSS01 == 1 )"
  tau_veto	    = "( ntau == 0 )"
  zpeak_veto	    = "( ( nmuon == 1 && nel == 1 ) || ( TMath::Abs( mll01 - 91.187e3 ) > 7.5e3 ) )"
  truth_cut         = "( 1 )"
  if useTruth:
    # Select non-prompt or QMisID probe leptons
    #
    truth_cut	    = "( isMC==0 || ( isMC==1 && ( ( lep_probe_truthType[0] != 6 && lep_probe_truthType[0] != 2 ) || ( lep_probe_isChFlip[0] == 1 ) ) ) )"

  if useGroupNTup:
    mytree_weight     = "3.209 * mcWeightOrg * pileupEventWeight_090 * weight_tag * weight_probe * JVT_EventWeight * MV2c20_77_EventWeight"
    trig_dec	      = "( passEventCleaning == 1 && ( isMC == 1 && HLT_e24_lhmedium_L1EM18VH == 1 ) || ( isMC == 0 && HLT_e24_lhmedium_L1EM20VH == 1 ) || ( HLT_e60_lhmedium == 1 ) || ( HLT_e120_lhloose == 1 ) || ( HLT_mu20_iloose_L1MU15 == 1 ) || ( HLT_mu50 == 1 ) )"
    trigmatching      = "( ( lep_isTrigMatch_0 == 1 && ( ( TMath::Abs( lep_ID_0 ) == 11 && lep_Pt_0 > 25e3 ) || ( TMath::Abs( lep_ID_0 ) == 13 && lep_Pt_0 > 21e3 ) ) ) || ( lep_isTrigMatch_1 == 1 && ( ( TMath::Abs( lep_ID_1 ) == 11 && lep_Pt_1 > 25e3 ) || ( TMath::Abs( lep_ID_1 ) == 13 && lep_Pt_1 > 21e3 ) ) ) )"
    lep_tag_trigmatch = "( lep_Tag_isTightSelected == 1 && lep_Tag_isTrigMatch == 1 && ( ( TMath::Abs( lep_Tag_ID ) == 11 && lep_Tag_Pt > 25e3 ) || ( TMath::Abs( lep_Tag_ID ) == 13 && lep_Tag_Pt > 21e3 ) ) )"
    el_eta            = "( largeEtaEvent( nelectrons,lep_ID_0,lep_ID_1,lep_EtaBE2_0,lep_EtaBE2_1 ) == 0 )"
    nbjets	      = "( nJets_OR_MV2c20_77 > 0 )"
    njets	      = "( nJets_OR > 1 && nJets_OR <= 4 )"
    nleptons	      = "( nleptons == 2 && ( lep_Pt_0 > 10e3 && lep_Pt_1 > 10e3 ) )"
    flavour_cut       = "( 1 )"
    same_sign	      = "( isSS01 == 1 )"
    tau_veto	      = "( nTaus_OR_Pt25 == 0 )"
    zpeak_veto	      = "( ( nmuons == 1 && nelectrons == 1 ) || ( TMath::Abs( Mll01 - 91.187e3 ) > 7.5e3 ) )"
    truth_cut         = "( 1 )"
    if useTruth:
      # Select non-prompt or QMisID probe leptons
      #
      truth_cut	      = "( isMC==0 || ( isMC==1 && lep_Probe_isPrompt != 1 ) )"

  hist_list = []

  basecut =  "(" + trig_dec + " && " + trigmatching + " && " + lep_tag_trigmatch + " && " + el_eta + " && " + njets + " && " + nbjets + " && " + nleptons + " && " + flavour_cut + " && " + same_sign + " && " + tau_veto + " && "  + zpeak_veto + " && " + truth_cut + " && "

  varProbeType   = ('lep_probe_truthType','TO_BE_ADDED')[bool(useGroupNTup)]
  varProbeOrigin = ('lep_probe_truthOrigin','TO_BE_ADDED')[bool(useGroupNTup)]
  varProbeChFlip = ('lep_probe_isChFlip','lep_Probe_isBremsElec')[bool(useGroupNTup)]

  varTagType     = ('lep_tag_truthType','TO_BE_ADDED')[bool(useGroupNTup)]
  varTagOrigin   = ('lep_tag_truthOrigin','TO_BE_ADDED')[bool(useGroupNTup)]
  varTagChFlip   = ('lep_tag_isChFlip','lep_Tag_isBremsElec')[bool(useGroupNTup)]

  # ------------
  #    Muons
  # ------------

  sel_FAKE_MU_PROBE_T = basecut + "( ( isProbeMuEvent == 1 ) && ( muon_probe_isTightSelected[0] == 1 ) ) )"
  sel_FAKE_MU_PROBE_L = basecut + "( ( isProbeMuEvent == 1 ) && ( muon_probe_isTightSelected[0] == 0 ) ) )"

  if useGroupNTup:
    sel_FAKE_MU_PROBE_T = basecut + "( ( TMath::Abs( lep_Probe_ID ) == 13 ) && ( lep_Probe_isTightSelected == 1 ) ) )"
    sel_FAKE_MU_PROBE_L = basecut + "( ( TMath::Abs( lep_Probe_ID ) == 13 ) && ( lep_Probe_isTightSelected == 0 ) ) )"


  print "sel_FAKE_MU_PROBE_T: \n", sel_FAKE_MU_PROBE_T
  print "sel_FAKE_MU_PROBE_L: \n", sel_FAKE_MU_PROBE_L

  h_FAKE_MU_PROBE_T_type_VS_origin = TH2D("FAKE_MU_PROBE_T_type_VS_origin", "FAKE_MU_PROBE_T_type_VS_origin", 18, 0.0, 18.0, 40, 0.0, 40.0)
  hist_list.append(h_FAKE_MU_PROBE_T_type_VS_origin)
  h_FAKE_MU_PROBE_T_type_VS_origin.GetXaxis().SetTitle("Probe lep (T) truthType")
  h_FAKE_MU_PROBE_T_type_VS_origin.GetYaxis().SetTitle("Probe lep (T) truthOrigin")
  h_FAKE_MU_PROBE_L_type_VS_origin = TH2D("FAKE_MU_PROBE_L_type_VS_origin", "FAKE_MU_PROBE_L_type_VS_origin", 18, 0.0, 18.0, 40, 0.0, 40.0)
  hist_list.append(h_FAKE_MU_PROBE_L_type_VS_origin)
  h_FAKE_MU_PROBE_L_type_VS_origin.GetXaxis().SetTitle("Probe lep (L!T) truthType")
  h_FAKE_MU_PROBE_L_type_VS_origin.GetYaxis().SetTitle("Probe lep (L!T) truthOrigin")

  mytree.Project("FAKE_MU_PROBE_T_type_VS_origin", varProbeOrigin + ":" + varProbeType, "%s * (%s)" %(mytree_weight, sel_FAKE_MU_PROBE_T), "text" )
  mytree.Project("FAKE_MU_PROBE_L_type_VS_origin", varProbeOrigin + ":" + varProbeType, "%s * (%s)" %(mytree_weight, sel_FAKE_MU_PROBE_L), "text" )

  # Normalise to unity
  #
  if doNorm:
    if h_FAKE_MU_PROBE_T_type_VS_origin.Integral() > 0:
      h_FAKE_MU_PROBE_T_type_VS_origin.Scale(1.0/h_FAKE_MU_PROBE_T_type_VS_origin.Integral())
    if h_FAKE_MU_PROBE_L_type_VS_origin.Integral() > 0:
      h_FAKE_MU_PROBE_L_type_VS_origin.Scale(1.0/h_FAKE_MU_PROBE_L_type_VS_origin.Integral())

  # -----------------------------------------

  h_FAKE_MU_TAG_T_type_VS_origin = TH2D("FAKE_MU_TAG_T_type_VS_origin", "FAKE_MU_TAG_T_type_VS_origin", 18, 0.0, 18.0, 40, 0.0, 40.0)
  hist_list.append(h_FAKE_MU_TAG_T_type_VS_origin)
  h_FAKE_MU_TAG_T_type_VS_origin.GetXaxis().SetTitle("Tag lep ( Probe (T) ) truthType")
  h_FAKE_MU_TAG_T_type_VS_origin.GetYaxis().SetTitle("Tag lep ( Probe (T) ) truthOrigin")
  h_FAKE_MU_TAG_L_type_VS_origin = TH2D("FAKE_MU_TAG_L_type_VS_origin", "FAKE_MU_TAG_L_type_VS_origin", 18, 0.0, 18.0, 40, 0.0, 40.0)
  hist_list.append(h_FAKE_MU_TAG_L_type_VS_origin)
  h_FAKE_MU_TAG_L_type_VS_origin.GetXaxis().SetTitle("Tag lep ( Probe (L!T) ) truthType")
  h_FAKE_MU_TAG_L_type_VS_origin.GetYaxis().SetTitle("Tag lep ( Probe (L!T) ) truthOrigin")

  mytree.Project("FAKE_MU_TAG_T_type_VS_origin", varTagOrigin + ":" + varTagType, "%s * (%s)" %(mytree_weight, sel_FAKE_MU_PROBE_T), "text" )
  mytree.Project("FAKE_MU_TAG_L_type_VS_origin", varTagOrigin + ":" + varTagType, "%s * (%s)" %(mytree_weight, sel_FAKE_MU_PROBE_L), "text" )

  # Normalise to unity
  #
  if doNorm:
    if h_FAKE_MU_TAG_T_type_VS_origin.Integral() > 0:
      h_FAKE_MU_TAG_T_type_VS_origin.Scale(1.0/h_FAKE_MU_TAG_T_type_VS_origin.Integral())
    if h_FAKE_MU_TAG_L_type_VS_origin.Integral() > 0:
      h_FAKE_MU_TAG_L_type_VS_origin.Scale(1.0/h_FAKE_MU_TAG_L_type_VS_origin.Integral())

  #--------------------------------------------------------------------------------
  # probe truthOrigin VS tag truthOrigin
  """
  h_FAKE_MU_PROBE_T_ProbeOrigin_VS_TagOrigin = TH2D("FAKE_MU_PROBE_T_ProbeOrigin_VS_TagOrigin", "FAKE_MU_PROBE_T_ProbeOrigin_VS_TagOrigin", 40, 0.0, 40.0, 40, 0.0, 40.0)
  hist_list.append(h_FAKE_MU_PROBE_T_ProbeOrigin_VS_TagOrigin)
  h_FAKE_MU_PROBE_T_ProbeOrigin_VS_TagOrigin.GetXaxis().SetTitle("Tag lep (T) truthOrigin")
  h_FAKE_MU_PROBE_T_ProbeOrigin_VS_TagOrigin.GetYaxis().SetTitle("Probe lep (T) truthOrigin")
  h_FAKE_MU_PROBE_L_ProbeOrigin_VS_TagOrigin = TH2D("FAKE_MU_PROBE_L_ProbeOrigin_VS_TagOrigin", "FAKE_MU_PROBE_L_ProbeOrigin_VS_TagOrigin", 40, 0.0, 40.0, 40, 0.0, 40.0)
  hist_list.append(h_FAKE_MU_PROBE_L_ProbeOrigin_VS_TagOrigin)
  h_FAKE_MU_PROBE_L_ProbeOrigin_VS_TagOrigin.GetXaxis().SetTitle("Tag lep (L!T) truthOrigin")
  h_FAKE_MU_PROBE_L_ProbeOrigin_VS_TagOrigin.GetYaxis().SetTitle("Probe lep (L!T) truthOrigin")

  mytree.Project("FAKE_MU_PROBE_T_ProbeOrigin_VS_TagOrigin", varProbeOrigin + ":" + varTagOrigin, "%s * (%s)" %(mytree_weight, sel_FAKE_MU_PROBE_T), "text" )
  mytree.Project("FAKE_MU_PROBE_L_ProbeOrigin_VS_TagOrigin", varProbeOrigin + ":" + varTagOrigin, "%s * (%s)" %(mytree_weight, sel_FAKE_MU_PROBE_L), "text" )

  h_FAKE_MU_PROBE_T_ProbeOrigin_VS_TagOrigin.Scale(1.0/h_FAKE_MU_PROBE_T_ProbeOrigin_VS_TagOrigin.Integral())
  h_FAKE_MU_PROBE_L_ProbeOrigin_VS_TagOrigin.Scale(1.0/h_FAKE_MU_PROBE_L_ProbeOrigin_VS_TagOrigin.Integral())
  """

  # ------------
  #  Electrons
  # ------------

  sel_FAKE_EL_PROBE_T = basecut + "( ( isProbeElEvent == 1 ) && ( el_probe_isTightSelected[0] == 1 ) ) )"
  sel_FAKE_EL_PROBE_L = basecut + "( ( isProbeElEvent == 1 ) && ( el_probe_isTightSelected[0] == 0 ) ) )"

  if useGroupNTup:
    sel_FAKE_EL_PROBE_T = basecut + "( ( TMath::Abs( lep_Probe_ID ) == 11 ) && ( lep_Probe_isTightSelected == 1 ) ) )"
    sel_FAKE_EL_PROBE_L = basecut + "( ( TMath::Abs( lep_Probe_ID ) == 11 ) && ( lep_Probe_isTightSelected == 0 ) ) )"

  print "sel_FAKE_EL_PROBE_T: \n", sel_FAKE_EL_PROBE_T
  print "sel_FAKE_EL_PROBE_L: \n", sel_FAKE_EL_PROBE_L

  h_FAKE_EL_PROBE_T_type_VS_origin = TH2D("FAKE_EL_PROBE_T_type_VS_origin", "FAKE_EL_PROBE_T_type_VS_origin", 18, 0.0, 18.0, 40, 0.0, 40.0)
  hist_list.append(h_FAKE_EL_PROBE_T_type_VS_origin)
  h_FAKE_EL_PROBE_T_type_VS_origin.GetXaxis().SetTitle("Probe lep (T) truthType")
  h_FAKE_EL_PROBE_T_type_VS_origin.GetYaxis().SetTitle("Probe lep (T) truthOrigin")
  h_FAKE_EL_PROBE_L_type_VS_origin = TH2D("FAKE_EL_PROBE_L_type_VS_origin", "FAKE_EL_PROBE_L_type_VS_origin", 18, 0.0, 18.0, 40, 0.0, 40.0)
  hist_list.append(h_FAKE_EL_PROBE_L_type_VS_origin)
  h_FAKE_EL_PROBE_L_type_VS_origin.GetXaxis().SetTitle("Probe lep (L!T) truthType")
  h_FAKE_EL_PROBE_L_type_VS_origin.GetYaxis().SetTitle("Probe lep (L!T) truthOrigin")

  h_FAKE_EL_PROBE_T_isChFlip = TH1D("FAKE_EL_PROBE_T_isChFlip", "FAKE_EL_PROBE_T_isChFlip", 2, -0.5, 1.5)
  hist_list.append(h_FAKE_EL_PROBE_T_isChFlip)
  h_FAKE_EL_PROBE_T_isChFlip.GetXaxis().SetTitle("Probe lep (T) isChFlip")
  h_FAKE_EL_PROBE_L_isChFlip = TH1D("FAKE_EL_PROBE_L_isChFlip", "FAKE_EL_PROBE_L_isChFlip", 2, -0.5, 1.5)
  hist_list.append(h_FAKE_EL_PROBE_L_isChFlip)
  h_FAKE_EL_PROBE_L_isChFlip.GetXaxis().SetTitle("Probe lep (L!T) isChFlip")

  mytree.Project("FAKE_EL_PROBE_T_type_VS_origin", varProbeOrigin + ":" + varProbeType, "%s * (%s)" %(mytree_weight, sel_FAKE_EL_PROBE_T), "text" )
  mytree.Project("FAKE_EL_PROBE_L_type_VS_origin", varProbeOrigin + ":" + varProbeType, "%s * (%s)" %(mytree_weight, sel_FAKE_EL_PROBE_L), "text" )
  #mytree.Project("FAKE_EL_PROBE_T_isChFlip", varProbeChFlip, "%s * (%s)" %(mytree_weight, sel_FAKE_EL_PROBE_T), "text" )
  #mytree.Project("FAKE_EL_PROBE_L_isChFlip", varProbeChFlip, "%s * (%s)" %(mytree_weight, sel_FAKE_EL_PROBE_L), "text" )

  # Normalise to unity
  #
  if doNorm:
    if h_FAKE_EL_PROBE_T_type_VS_origin.Integral() > 0 :
      h_FAKE_EL_PROBE_T_type_VS_origin.Scale(1.0/h_FAKE_EL_PROBE_T_type_VS_origin.Integral())
    if h_FAKE_EL_PROBE_L_type_VS_origin.Integral() > 0 :
      h_FAKE_EL_PROBE_L_type_VS_origin.Scale(1.0/h_FAKE_EL_PROBE_L_type_VS_origin.Integral())
    if h_FAKE_EL_PROBE_T_isChFlip.Integral() > 0 :
      h_FAKE_EL_PROBE_T_isChFlip.Scale(1.0/h_FAKE_EL_PROBE_T_isChFlip.Integral())
    if h_FAKE_EL_PROBE_L_isChFlip.Integral() > 0 :
      h_FAKE_EL_PROBE_L_isChFlip.Scale(1.0/h_FAKE_EL_PROBE_L_isChFlip.Integral())

  # -----------------------------------------

  h_FAKE_EL_TAG_T_type_VS_origin = TH2D("FAKE_EL_TAG_T_type_VS_origin", "FAKE_EL_TAG_T_type_VS_origin", 18, 0.0, 18.0, 40, 0.0, 40.0)
  hist_list.append(h_FAKE_EL_TAG_T_type_VS_origin)
  h_FAKE_EL_TAG_T_type_VS_origin.GetXaxis().SetTitle("Tag lep ( Probe (T) ) truthType")
  h_FAKE_EL_TAG_T_type_VS_origin.GetYaxis().SetTitle("Tag lep ( Probe (T) ) truthOrigin")
  h_FAKE_EL_TAG_L_type_VS_origin = TH2D("FAKE_EL_TAG_L_type_VS_origin", "FAKE_EL_TAG_L_type_VS_origin", 18, 0.0, 18.0, 40, 0.0, 40.0)
  hist_list.append(h_FAKE_EL_TAG_L_type_VS_origin)
  h_FAKE_EL_TAG_L_type_VS_origin.GetXaxis().SetTitle("Tag lep ( Probe (L!T) ) truthType")
  h_FAKE_EL_TAG_L_type_VS_origin.GetYaxis().SetTitle("Tag lep ( Probe (L!T) ) truthOrigin")

  mytree.Project("FAKE_EL_TAG_T_type_VS_origin", varTagOrigin + ":" + varTagType, "%s * (%s)" %(mytree_weight, sel_FAKE_EL_PROBE_T), "text" )
  mytree.Project("FAKE_EL_TAG_L_type_VS_origin", varTagOrigin + ":" + varTagType, "%s * (%s)" %(mytree_weight, sel_FAKE_EL_PROBE_L), "text" )

  # Normalise to unity
  #
  if doNorm:
    if h_FAKE_EL_TAG_T_type_VS_origin.Integral() >0 :
      h_FAKE_EL_TAG_T_type_VS_origin.Scale(1/h_FAKE_EL_TAG_T_type_VS_origin.Integral())
    if h_FAKE_EL_TAG_L_type_VS_origin.Integral() >0 :
      h_FAKE_EL_TAG_L_type_VS_origin.Scale(1/h_FAKE_EL_TAG_L_type_VS_origin.Integral())

  #--------------------------------------------------------------------------------
  # probe truthOrigin VS tag truthOrigin
  """
  h_FAKE_EL_PROBE_T_ProbeOrigin_VS_TagOrigin = TH2D("FAKE_EL_PROBE_T_ProbeOrigin_VS_TagOrigin", "FAKE_EL_PROBE_T_ProbeOrigin_VS_TagOrigin", 40, 0.0, 40.0, 40, 0.0, 40.0)
  hist_list.append(h_FAKE_EL_PROBE_T_ProbeOrigin_VS_TagOrigin)
  h_FAKE_EL_PROBE_T_ProbeOrigin_VS_TagOrigin.GetXaxis().SetTitle("Tag lep (T) truthOrigin")
  h_FAKE_EL_PROBE_T_ProbeOrigin_VS_TagOrigin.GetYaxis().SetTitle("Probe lep (T) truthOrigin")
  h_FAKE_EL_PROBE_L_ProbeOrigin_VS_TagOrigin = TH2D("FAKE_EL_PROBE_L_ProbeOrigin_VS_TagOrigin", "FAKE_EL_PROBE_L_ProbeOrigin_VS_TagOrigin", 40, 0.0, 40.0, 40, 0.0, 40.0)
  hist_list.append(h_FAKE_EL_PROBE_L_ProbeOrigin_VS_TagOrigin)
  h_FAKE_EL_PROBE_L_ProbeOrigin_VS_TagOrigin.GetXaxis().SetTitle("Tag lep (L!T) truthOrigin")
  h_FAKE_EL_PROBE_L_ProbeOrigin_VS_TagOrigin.GetYaxis().SetTitle("Probe lep (L!T) truthOrigin")

  mytree.Project("FAKE_EL_PROBE_T_ProbeOrigin_VS_TagOrigin",  varProbeOrigin + ":" + varTagOrigin, "%s * (%s)" %(mytree_weight, sel_FAKE_EL_PROBE_T), "text" )
  mytree.Project("FAKE_EL_PROBE_L_ProbeOrigin_VS_TagOrigin",  varProbeOrigin + ":" + varTagOrigin, "%s * (%s)" %(mytree_weight, sel_FAKE_EL_PROBE_L), "text" )

  h_FAKE_EL_PROBE_T_ProbeOrigin_VS_TagOrigin.Scale(1.0/h_FAKE_EL_PROBE_T_ProbeOrigin_VS_TagOrigin.Integral())
  h_FAKE_EL_PROBE_L_ProbeOrigin_VS_TagOrigin.Scale(1.0/h_FAKE_EL_PROBE_L_ProbeOrigin_VS_TagOrigin.Integral())
  """

  # ------------------------------------------------------------------

  set_fancy_2D_style()

  gStyle.SetPaintTextFormat("2.1f")

  for hist in hist_list:

     #c = TCanvas("c1","Temp",50,50,700,900)
     c = TCanvas("c1","Temp",50,50,800,900)

     legend = TLegend(0.2,0.8,0.7,0.85); # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
     legend.SetHeader("2 Lep SS Fake CR")
     legend.SetBorderSize(0)  # no border
     legend.SetFillStyle(0); # Legend transparent background
     legend.SetTextSize(0.04) # Increase entry font size!
     legend.SetTextFont(42)   # Helvetica
     leg_ATLAS = TLatex()
     leg_lumi  = TLatex()
     leg_ATLAS.SetTextSize(0.03)
     leg_ATLAS.SetNDC()
     leg_lumi.SetTextSize(0.03)
     leg_lumi.SetNDC()

     hist.Draw("colz text")

     legend.Draw()
     leg_ATLAS.DrawLatex(0.2,0.75,"#bf{#it{ATLAS}} Work In Progress")
     leg_lumi.DrawLatex(0.2,0.7,"#sqrt{s} = 13 TeV, #int L dt = 3.2 fb^{-1}")

     canvasname = ( hist.GetName(), hist.GetName() + "_TruthCut" )[bool(useTruth)]
     c.SaveAs(canvasname + ".png")

# ------------------------------------------------------------------

def plot_2D_RealCR():

  myfilename = ("/coepp/cephfs/mel/mmilesi/ttH/MergedDatasets/Merged_v030/Merged_Melb15_ttH_030_DxAOD_p2559/tops/410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.root","/coepp/cephfs/mel/mmilesi/ttH/MiniNTup/25ns_v7/25ns_v7/tops/410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.root")[bool(useGroupNTup)]
  #myfilename = ("/coepp/cephfs/mel/mmilesi/ttH/MergedDatasets/Merged_v031/Merged_Melb15_ttH_031-01_DxAOD_p2559/tops/410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.root","TO_BE_ADDED")[bool(useGroupNTup)]

  myfile     = TFile(myfilename)
  mytree     = gDirectory.Get("physics")

  gROOT.LoadMacro("$ROOTCOREBIN/user_scripts/HTopMultilepAnalysis/ROOT_TTreeFormulas/largeEtaEvent.cxx+")
  from ROOT import largeEtaEvent

  mytree_weight     = "3.209 * mcEventWeight * weight_pileup * weight_lepton_trig_HTop[0] * weight_lepton_reco_HTop[0] * weight_lepton_iso_HTop[0] * weight_lepton_ID_HTop[0] * weight_lepton_TTVA_HTop[0] * weight_jet_JVT_HTop[0] * weight_jet_MV2c20_SFFix77[0]"
  trig_dec	    = "( passHLT == 1 )"
  trigmatching      = "( ( lep_isTrigMatched[0] == 1 && ( ( lep_flavour[0] == 11 && lep_pt[0] > 25e3 ) || ( lep_flavour[0] == 13 && lep_pt[0] > 21e3 ) ) ) || ( lep_isTrigMatched[1] == 1 && ( ( lep_flavour[1] == 11 && lep_pt[1] > 25e3 ) || ( lep_flavour[1] == 13 && lep_pt[1] > 21e3 ) ) ) )"
  lep_tag_trigmatch = "( lep_tag_isTightSelected[0] == 1 && lep_tag_isTrigMatched[0] == 1 && ( ( lep_tag_flavour[0] == 11 && lep_tag_pt[0] > 25e3 ) || ( lep_tag_flavour[0] == 13 && lep_tag_pt[0] > 21e3 ) ) )"
  el_eta	    = "( nel == 0 || ( nel > 0 && Max$( TMath::Abs(el_caloCluster_eta) ) < 1.37 ) )"
  nbjets	    = "( njets_mv2c20_Fix77 > 0 )"
  njets 	    = "( njets > 1 && njets <= 4 )"
  nleptons	    = "( nlep == 2 && Min$( lep_pt ) > 10e3 )"
  flavour_cut	    = "( 1 )"
  opposite_sign	    = "( isSS01 != 1 )"
  tau_veto	    = "( ntau == 0 )"
  zpeak_veto	    = "( ( nmuon == 1 && nel == 1 ) || ( TMath::Abs( mll01 - 91.187e3 ) > 7.5e3 ) )"
  truth_cut         = "( 1 )"
  if useTruth:
    # Select non-prompt or QMisID probe leptons
    #
    truth_cut	    = "( isMC==0 || ( isMC==1 && ( ( lep_probe_truthType[0] != 6 && lep_probe_truthType[0] != 2 ) || ( lep_probe_isChFlip[0] == 1 ) ) ) )"

  if useGroupNTup:
    mytree_weight     = "3.209 * mcWeightOrg * pileupEventWeight_090 * weight_tag * weight_probe * JVT_EventWeight * MV2c20_77_EventWeight"
    trig_dec	      = "( passEventCleaning == 1 && ( isMC == 1 && HLT_e24_lhmedium_L1EM18VH == 1 ) || ( isMC == 0 && HLT_e24_lhmedium_L1EM20VH == 1 ) || ( HLT_e60_lhmedium == 1 ) || ( HLT_e120_lhloose == 1 ) || ( HLT_mu20_iloose_L1MU15 == 1 ) || ( HLT_mu50 == 1 ) )"
    trigmatching      = "( ( lep_isTrigMatch_0 == 1 && ( ( TMath::Abs( lep_ID_0 ) == 11 && lep_Pt_0 > 25e3 ) || ( TMath::Abs( lep_ID_0 ) == 13 && lep_Pt_0 > 21e3 ) ) ) || ( lep_isTrigMatch_1 == 1 && ( ( TMath::Abs( lep_ID_1 ) == 11 && lep_Pt_1 > 25e3 ) || ( TMath::Abs( lep_ID_1 ) == 13 && lep_Pt_1 > 21e3 ) ) ) )"
    lep_tag_trigmatch = "( lep_Tag_isTightSelected == 1 && lep_Tag_isTrigMatch == 1 && ( ( TMath::Abs( lep_Tag_ID ) == 11 && lep_Tag_Pt > 25e3 ) || ( TMath::Abs( lep_Tag_ID ) == 13 && lep_Tag_Pt > 21e3 ) ) )"
    el_eta            = "( largeEtaEvent( nelectrons,lep_ID_0,lep_ID_1,lep_EtaBE2_0,lep_EtaBE2_1 ) == 0 )"
    nbjets	      = "( nJets_OR_MV2c20_77 > 0 )"
    njets	      = "( nJets_OR > 1 && nJets_OR <= 4 )"
    nleptons	      = "( nleptons == 2 && ( lep_Pt_0 > 10e3 && lep_Pt_1 > 10e3 ) )"
    flavour_cut       = "( 1 )"
    opposite_sign     = "( isSS01 != 1 )"
    tau_veto	      = "( nTaus_OR_Pt25 == 0 )"
    zpeak_veto	      = "( ( nmuons == 1 && nelectrons == 1 ) || ( TMath::Abs( Mll01 - 91.187e3 ) > 7.5e3 ) )"
    truth_cut         = "( 1 )"
    if useTruth:
      # Select non-prompt or QMisID probe leptons
      #
      truth_cut	      = "( isMC==0 || ( isMC==1 && lep_Probe_isPrompt != 1 ) )"

  hist_list = []

  basecut =  "(" + trig_dec + " && " + trigmatching + " && " + lep_tag_trigmatch + " && " + el_eta + " && " + njets + " && " + nbjets + " && " + nleptons + " && " + flavour_cut + " && " + opposite_sign + " && " + tau_veto + " && "  + zpeak_veto + " && " + truth_cut + " && "

  varProbeType   = ('lep_probe_truthType','TO_BE_ADDED')[bool(useGroupNTup)]
  varProbeOrigin = ('lep_probe_truthOrigin','TO_BE_ADDED')[bool(useGroupNTup)]
  varProbeChFlip = ('lep_probe_isChFlip','lep_Probe_isBremsElec')[bool(useGroupNTup)]

  varTagType     = ('lep_tag_truthType','TO_BE_ADDED')[bool(useGroupNTup)]
  varTagOrigin   = ('lep_tag_truthOrigin','TO_BE_ADDED')[bool(useGroupNTup)]
  varTagChFlip   = ('lep_tag_isChFlip','lep_Tag_isBremsElec')[bool(useGroupNTup)]

  # ------------
  #    Muons
  # ------------

  sel_REAL_MU_PROBE_T = basecut + "( ( isProbeMuEvent == 1 ) && ( muon_probe_isTightSelected[0] == 1 ) ) )"
  sel_REAL_MU_PROBE_L = basecut + "( ( isProbeMuEvent == 1 ) && ( muon_probe_isTightSelected[0] == 0 ) ) )"

  if useGroupNTup:
    sel_REAL_MU_PROBE_T = basecut + "( ( TMath::Abs( lep_Probe_ID ) == 13 ) && ( lep_Probe_isTightSelected == 1 ) ) )"
    sel_REAL_MU_PROBE_L = basecut + "( ( TMath::Abs( lep_Probe_ID ) == 13 ) && ( lep_Probe_isTightSelected == 0 ) ) )"

  print "sel_REAL_MU_PROBE_T: \n", sel_REAL_MU_PROBE_T
  print "sel_REAL_MU_PROBE_L: \n", sel_REAL_MU_PROBE_L

  h_REAL_MU_PROBE_T_type_VS_origin = TH2D("REAL_MU_PROBE_T_type_VS_origin", "REAL_MU_PROBE_T_type_VS_origin", 18, 0.0, 18.0, 40, 0.0, 40.0)
  hist_list.append(h_REAL_MU_PROBE_T_type_VS_origin)
  h_REAL_MU_PROBE_T_type_VS_origin.GetXaxis().SetTitle("Probe lep (T) truthType")
  h_REAL_MU_PROBE_T_type_VS_origin.GetYaxis().SetTitle("Probe lep (T) truthOrigin")
  h_REAL_MU_PROBE_L_type_VS_origin = TH2D("REAL_MU_PROBE_L_type_VS_origin", "REAL_MU_PROBE_L_type_VS_origin", 18, 0.0, 18.0, 40, 0.0, 40.0)
  hist_list.append(h_REAL_MU_PROBE_L_type_VS_origin)
  h_REAL_MU_PROBE_L_type_VS_origin.GetXaxis().SetTitle("Probe lep (L!T) truthType")
  h_REAL_MU_PROBE_L_type_VS_origin.GetYaxis().SetTitle("Probe lep (L!T) truthOrigin")

  mytree.Project("REAL_MU_PROBE_T_type_VS_origin", varProbeOrigin + ":" + varProbeType, "%s * (%s)" %(mytree_weight, sel_REAL_MU_PROBE_T), "text" )
  mytree.Project("REAL_MU_PROBE_L_type_VS_origin", varProbeOrigin + ":" + varProbeType, "%s * (%s)" %(mytree_weight, sel_REAL_MU_PROBE_L), "text" )

  # Normalise to unity
  #
  if doNorm:
    if h_REAL_MU_PROBE_T_type_VS_origin.Integral() > 0:
      h_REAL_MU_PROBE_T_type_VS_origin.Scale(1.0/h_REAL_MU_PROBE_T_type_VS_origin.Integral())
    if h_REAL_MU_PROBE_L_type_VS_origin.Integral() > 0:
      h_REAL_MU_PROBE_L_type_VS_origin.Scale(1.0/h_REAL_MU_PROBE_L_type_VS_origin.Integral())

  # -----------------------------------------

  h_REAL_MU_TAG_T_type_VS_origin = TH2D("REAL_MU_TAG_T_type_VS_origin", "REAL_MU_TAG_T_type_VS_origin", 18, 0.0, 18.0, 40, 0.0, 40.0)
  hist_list.append(h_REAL_MU_TAG_T_type_VS_origin)
  h_REAL_MU_TAG_T_type_VS_origin.GetXaxis().SetTitle("Tag lep ( Probe (T) ) truthType")
  h_REAL_MU_TAG_T_type_VS_origin.GetYaxis().SetTitle("Tag lep ( Probe (T) ) truthOrigin")
  h_REAL_MU_TAG_L_type_VS_origin = TH2D("REAL_MU_TAG_L_type_VS_origin", "REAL_MU_TAG_L_type_VS_origin", 18, 0.0, 18.0, 40, 0.0, 40.0)
  hist_list.append(h_REAL_MU_TAG_L_type_VS_origin)
  h_REAL_MU_TAG_L_type_VS_origin.GetXaxis().SetTitle("Tag lep ( Probe (L!T) ) truthType")
  h_REAL_MU_TAG_L_type_VS_origin.GetYaxis().SetTitle("Tag lep ( Probe (L!T) ) truthOrigin")

  mytree.Project("REAL_MU_TAG_T_type_VS_origin", varTagOrigin + ":" + varTagType, "%s * (%s)" %(mytree_weight, sel_REAL_MU_PROBE_T), "text" )
  mytree.Project("REAL_MU_TAG_L_type_VS_origin", varTagOrigin + ":" + varTagType, "%s * (%s)" %(mytree_weight, sel_REAL_MU_PROBE_L), "text" )

  # Normalise to unity
  if h_REAL_MU_TAG_T_type_VS_origin.Integral() > 0:
    h_REAL_MU_TAG_T_type_VS_origin.Scale(1.0/h_REAL_MU_TAG_T_type_VS_origin.Integral())
  if h_REAL_MU_TAG_L_type_VS_origin.Integral() > 0:
    h_REAL_MU_TAG_L_type_VS_origin.Scale(1.0/h_REAL_MU_TAG_L_type_VS_origin.Integral())

  # ------------
  #  Electrons
  # ------------

  sel_REAL_EL_PROBE_T = basecut + "( ( isProbeElEvent == 1 ) && ( el_probe_isTightSelected[0] == 1 ) ) )"
  sel_REAL_EL_PROBE_L = basecut + "( ( isProbeElEvent == 1 ) && ( el_probe_isTightSelected[0] == 0 ) ) )"

  if useGroupNTup:
    sel_REAL_EL_PROBE_T = basecut + "( ( TMath::Abs( lep_Probe_ID ) == 11 ) && ( lep_Probe_isTightSelected == 1 ) ) )"
    sel_REAL_EL_PROBE_L = basecut + "( ( TMath::Abs( lep_Probe_ID ) == 11 ) && ( lep_Probe_isTightSelected == 0 ) ) )"

  print "sel_REAL_EL_PROBE_T: \n", sel_REAL_EL_PROBE_T
  print "sel_REAL_EL_PROBE_L: \n", sel_REAL_EL_PROBE_L

  h_REAL_EL_PROBE_T_type_VS_origin = TH2D("REAL_EL_PROBE_T_type_VS_origin", "REAL_EL_PROBE_T_type_VS_origin", 18, 0.0, 18.0, 40, 0.0, 40.0)
  hist_list.append(h_REAL_EL_PROBE_T_type_VS_origin)
  h_REAL_EL_PROBE_T_type_VS_origin.GetXaxis().SetTitle("Probe lep (T) truthType")
  h_REAL_EL_PROBE_T_type_VS_origin.GetYaxis().SetTitle("Probe lep (T) truthOrigin")
  h_REAL_EL_PROBE_L_type_VS_origin = TH2D("REAL_EL_PROBE_L_type_VS_origin", "REAL_EL_PROBE_L_type_VS_origin", 18, 0.0, 18.0, 40, 0.0, 40.0)
  hist_list.append(h_REAL_EL_PROBE_L_type_VS_origin)
  h_REAL_EL_PROBE_L_type_VS_origin.GetXaxis().SetTitle("Probe lep (L!T) truthType")
  h_REAL_EL_PROBE_L_type_VS_origin.GetYaxis().SetTitle("Probe lep (L!T) truthOrigin")

  h_REAL_EL_PROBE_T_isChFlip = TH1D("REAL_EL_PROBE_T_isChFlip", "REAL_EL_PROBE_T_isChFlip", 2, -0.5, 1.5)
  hist_list.append(h_REAL_EL_PROBE_T_isChFlip)
  h_REAL_EL_PROBE_T_isChFlip.GetXaxis().SetTitle("Probe lep (T) isChFlip")
  h_REAL_EL_PROBE_L_isChFlip = TH1D("REAL_EL_PROBE_L_isChFlip", "REAL_EL_PROBE_L_isChFlip", 2, -0.5, 1.5)
  hist_list.append(h_REAL_EL_PROBE_L_isChFlip)
  h_REAL_EL_PROBE_L_isChFlip.GetXaxis().SetTitle("Probe lep (L!T) isChFlip")

  mytree.Project("REAL_EL_PROBE_T_type_VS_origin",  varProbeOrigin + ":" + varProbeType, "%s * (%s)" %(mytree_weight, sel_REAL_EL_PROBE_T), "text" )
  mytree.Project("REAL_EL_PROBE_L_type_VS_origin",  varProbeOrigin + ":" + varProbeType, "%s * (%s)" %(mytree_weight, sel_REAL_EL_PROBE_L), "text" )
  #mytree.Project("REAL_EL_PROBE_T_isChFlip", varProbeChFlip, "%s * (%s)" %(mytree_weight, sel_REAL_EL_PROBE_T), "text" )
  #mytree.Project("REAL_EL_PROBE_L_isChFlip", varProbeChFlip, "%s * (%s)" %(mytree_weight, sel_REAL_EL_PROBE_L), "text" )

  # Normalise to unity
  #
  if doNorm:
    if h_REAL_EL_PROBE_T_type_VS_origin.Integral() > 0 :
      h_REAL_EL_PROBE_T_type_VS_origin.Scale(1.0/h_REAL_EL_PROBE_T_type_VS_origin.Integral())
    if h_REAL_EL_PROBE_L_type_VS_origin.Integral() > 0 :
      h_REAL_EL_PROBE_L_type_VS_origin.Scale(1.0/h_REAL_EL_PROBE_L_type_VS_origin.Integral())
    if h_REAL_EL_PROBE_T_isChFlip.Integral() > 0 :
      h_REAL_EL_PROBE_T_isChFlip.Scale(1.0/h_REAL_EL_PROBE_T_isChFlip.Integral())
    if h_REAL_EL_PROBE_L_isChFlip.Integral() > 0 :
      h_REAL_EL_PROBE_L_isChFlip.Scale(1.0/h_REAL_EL_PROBE_L_isChFlip.Integral())

  # -----------------------------------------

  h_REAL_EL_TAG_T_type_VS_origin = TH2D("REAL_EL_TAG_T_type_VS_origin", "REAL_EL_TAG_T_type_VS_origin", 10, 0.0, 10.0, 40, 0.0, 40.0)
  hist_list.append(h_REAL_EL_TAG_T_type_VS_origin)
  h_REAL_EL_TAG_T_type_VS_origin.GetXaxis().SetTitle("Tag lep ( Probe (T) ) truthType")
  h_REAL_EL_TAG_T_type_VS_origin.GetYaxis().SetTitle("Tag lep ( Probe (T) ) truthOrigin")
  h_REAL_EL_TAG_L_type_VS_origin = TH2D("REAL_EL_TAG_L_type_VS_origin", "REAL_EL_TAG_L_type_VS_origin", 10, 0.0, 10.0, 40, 0.0, 40.0)
  hist_list.append(h_REAL_EL_TAG_L_type_VS_origin)
  h_REAL_EL_TAG_L_type_VS_origin.GetXaxis().SetTitle("Tag lep ( Probe (L!T) ) truthType")
  h_REAL_EL_TAG_L_type_VS_origin.GetYaxis().SetTitle("Tag lep ( Probe (L!T) ) truthOrigin")

  mytree.Project("REAL_EL_TAG_T_type_VS_origin", "lep_tag_truthOrigin:lep_tag_truthType", "%s * (%s)" %(mytree_weight, sel_REAL_EL_PROBE_T), "text" )
  mytree.Project("REAL_EL_TAG_L_type_VS_origin", "lep_tag_truthOrigin:lep_tag_truthType", "%s * (%s)" %(mytree_weight, sel_REAL_EL_PROBE_L), "text" )

  # Normalise to unity
  #
  if doNorm:
    if h_REAL_EL_TAG_T_type_VS_origin.Integral() >0 :
      h_REAL_EL_TAG_T_type_VS_origin.Scale(1.0/h_REAL_EL_TAG_T_type_VS_origin.Integral())
    if h_REAL_EL_TAG_L_type_VS_origin.Integral() >0 :
      h_REAL_EL_TAG_L_type_VS_origin.Scale(1.0/h_REAL_EL_TAG_L_type_VS_origin.Integral())

  # ------------------------------------------------------------------

  set_fancy_2D_style()

  gStyle.SetPaintTextFormat("2.1f")

  for hist in hist_list:

     #c = TCanvas("c1","Temp",50,50,700,900)
     c = TCanvas("c1","Temp",50,50,800,900)

     legend = TLegend(0.2,0.8,0.7,0.85); # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
     legend.SetHeader("2 Lep OS Real CR")
     legend.SetBorderSize(0)  # no border
     legend.SetFillStyle(0); # Legend transparent background
     legend.SetTextSize(0.04) # Increase entry font size!
     legend.SetTextFont(42)   # Helvetica
     leg_ATLAS = TLatex()
     leg_lumi  = TLatex()
     leg_ATLAS.SetTextSize(0.03)
     leg_ATLAS.SetNDC()
     leg_lumi.SetTextSize(0.03)
     leg_lumi.SetNDC()

     hist.Draw("colz text")

     legend.Draw()
     leg_ATLAS.DrawLatex(0.2,0.75,"#bf{#it{ATLAS}} Work In Progress")
     leg_lumi.DrawLatex(0.2,0.7,"#sqrt{s} = 13 TeV, #int L dt = 3.2 fb^{-1}")

     canvasname = ( hist.GetName(), hist.GetName() + "_TruthCut" )[bool(useTruth)]
     c.SaveAs(canvasname + ".png")

# -----------------------------------------------------------------------------------

plot_2D_RealCR()
plot_2D_FakeCR()






