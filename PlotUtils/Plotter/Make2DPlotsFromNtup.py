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
  
  myfile = TFile("/data/mmilesi/ttH/MergedDatasets/Merged_Melb15_ttH_021_DxAOD/tops/PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.root")
  mytree = gDirectory.Get("physics")
  mytree_weight = mytree.GetWeight()/1e3
  
  trig_dec  = "( passHLT == 1 )" #"( passHLT == 1 && ( ( isMC == 1 && passedTriggers == \"HLT_e24_lhmedium_L1EM18VH\" ) || ( isMC == 0 && passedTriggers == \"HLT_e24_lhmedium_L1EM20VH\" ) || passedTriggers == \"HLT_e60_lhmedium\" || passedTriggers == \"HLT_e120_lhloose\" || passedTriggers == \"HLT_mu20_iloose_L1MU15\" || passedTriggers == \"HLT_mu50\" ) )"
  lep_tag_trigmatch = "( lep_tag_isTrigMatched[0] == 1 && ( ( lep_tag_flavour[0] == 11 && lep_tag_pt[0] > 28e3 ) || ( lep_tag_flavour[0] == 13 && lep_tag_pt[0] > 24e3 ) ) )"
  el_tag_eta   = "( TMath::Abs(el_tag_eta[0]) < 1.37 )"
  nbjets       = "( njets_mv2c20_Fix77 > 0 )"
  njets        = "( njets > 0 && njets < 4 )"
  nleptons     = "( nlep == 2 && ( lep_pt[0] > 20e3 && lep_pt[1] > 20e3 ) )"
  same_sign    = "( isSS01 == 1 )"
  # veto charge flips
  #
  #ch_flip_veto = "( lep_isChFlip[0] == 0 && lep_isChFlip[1] == 0 )" 
  ch_flip_veto = "( 1 )"  
  # require at least 1 !prompt lepton
  #
  non_prompt   = "( ( lep_truthType[0] != 6 && lep_truthType[0] != 2 ) || ( lep_truthType[1] != 6 && lep_truthType[1] != 2 ) )"   
  # require at least one T lepton in the event
  #  
  is_tight_event = "1"#"( isNonTightEvent == 0 )"  
  # require tag lepton to be T
  #
  tight_tag   = "( lep_tag_isTightSelected[0] == 1 )"
  
  hist_list = []
  
  # ------------
  #    Muons
  # ------------
   
  # look only at OF region  
  sel_FAKE_MU_PROBE_T = "(" + trig_dec + " && " + lep_tag_trigmatch + " && " + nleptons + " && " + same_sign + " && " +  is_tight_event + " && " + tight_tag + " && " + "( (nmuon == 1 && nel == 1) && ( isProbeMuEvent == 1 ) && ( muon_probe_isTightSelected[0] == 1 ) )" + " && " + el_tag_eta + " && " + nbjets + " && " + njets  + " && " + ch_flip_veto + " && " + non_prompt + ")"
  sel_FAKE_MU_PROBE_L = "(" + trig_dec + " && " + lep_tag_trigmatch + " && " + nleptons + " && " + same_sign + " && " +  is_tight_event + " && " + tight_tag + " && " + "( (nmuon == 1 && nel == 1) && ( isProbeMuEvent == 1 ) && ( muon_probe_isTightSelected[0] == 0 ) )" + " && " + el_tag_eta + " && " + nbjets + " && " + njets + " && " + ch_flip_veto + " && " + non_prompt + ")"

  # look only at SF region
  #sel_FAKE_MU_PROBE_T = "(" + trig_dec + " && " + lep_tag_trigmatch + " && " + nleptons + " && " + same_sign + " && " +  is_tight_event + " && " + tight_tag + " && " + "( ( nmuon == 2 ) && ( isProbeMuEvent == 1 ) && ( muon_probe_isTightSelected[0] == 1 ) )" + " && " + nbjets + " && " + njets + " && " + ch_flip_veto + " && " + non_prompt + ")"
  #sel_FAKE_MU_PROBE_L = "(" + trig_dec + " && " + lep_tag_trigmatch + " && " + nleptons + " && " + same_sign + " && " +  is_tight_event + " && " + tight_tag + " && " + "( ( nmuon == 2 ) && ( isProbeMuEvent == 1 ) && ( muon_probe_isTightSelected[0] == 0 ) )" + " && " + nbjets + " && " + njets + " && "  + ch_flip_veto + " && " + non_prompt + ")"
  
  print "sel_FAKE_MU_PROBE_T: \n", sel_FAKE_MU_PROBE_T
  print "sel_FAKE_MU_PROBE_L: \n", sel_FAKE_MU_PROBE_L
  
  h_FAKE_MU_PROBE_T_type_VS_origin = TH2D("FAKE_MU_PROBE_T_type_VS_origin", "FAKE_MU_PROBE_T_type_VS_origin", 10, 0.0, 10.0, 40, 0.0, 40.0)
  hist_list.append(h_FAKE_MU_PROBE_T_type_VS_origin)
  h_FAKE_MU_PROBE_T_type_VS_origin.GetXaxis().SetTitle("Probe lep (T) truthType")
  h_FAKE_MU_PROBE_T_type_VS_origin.GetYaxis().SetTitle("Probe lep (T) truthOrigin")
  h_FAKE_MU_PROBE_L_type_VS_origin = TH2D("FAKE_MU_PROBE_L_type_VS_origin", "FAKE_MU_PROBE_L_type_VS_origin", 10, 0.0, 10.0, 40, 0.0, 40.0)
  hist_list.append(h_FAKE_MU_PROBE_L_type_VS_origin)
  h_FAKE_MU_PROBE_L_type_VS_origin.GetXaxis().SetTitle("Probe lep (L!T) truthType")
  h_FAKE_MU_PROBE_L_type_VS_origin.GetYaxis().SetTitle("Probe lep (L!T) truthOrigin")
  
  mytree.Project("FAKE_MU_PROBE_T_type_VS_origin", "lep_probe_truthOrigin:lep_probe_truthType", "%s * (%s)" %(mytree_weight, sel_FAKE_MU_PROBE_T), "text" )
  mytree.Project("FAKE_MU_PROBE_L_type_VS_origin", "lep_probe_truthOrigin:lep_probe_truthType", "%s * (%s)" %(mytree_weight, sel_FAKE_MU_PROBE_L), "text" )

  # Normalise to unit
  h_FAKE_MU_PROBE_T_type_VS_origin.Scale(1/h_FAKE_MU_PROBE_T_type_VS_origin.Integral())
  h_FAKE_MU_PROBE_L_type_VS_origin.Scale(1/h_FAKE_MU_PROBE_L_type_VS_origin.Integral())

  # -----------------------------------------

  h_FAKE_MU_TAG_T_type_VS_origin = TH2D("FAKE_MU_TAG_T_type_VS_origin", "FAKE_MU_TAG_T_type_VS_origin", 10, 0.0, 10.0, 40, 0.0, 40.0)
  hist_list.append(h_FAKE_MU_TAG_T_type_VS_origin)
  h_FAKE_MU_TAG_T_type_VS_origin.GetXaxis().SetTitle("Tag lep ( Probe (T) ) truthType")
  h_FAKE_MU_TAG_T_type_VS_origin.GetYaxis().SetTitle("Tag lep ( Probe (T) ) truthOrigin")
  h_FAKE_MU_TAG_L_type_VS_origin = TH2D("FAKE_MU_TAG_L_type_VS_origin", "FAKE_MU_TAG_L_type_VS_origin", 10, 0.0, 10.0, 40, 0.0, 40.0)
  hist_list.append(h_FAKE_MU_TAG_L_type_VS_origin)
  h_FAKE_MU_TAG_L_type_VS_origin.GetXaxis().SetTitle("Tag lep ( Probe (L!T) ) truthType")
  h_FAKE_MU_TAG_L_type_VS_origin.GetYaxis().SetTitle("Tag lep ( Probe (L!T) ) truthOrigin")
  
  mytree.Project("FAKE_MU_TAG_T_type_VS_origin", "lep_tag_truthOrigin:lep_tag_truthType", "%s * (%s)" %(mytree_weight, sel_FAKE_MU_PROBE_T), "text" )
  mytree.Project("FAKE_MU_TAG_L_type_VS_origin", "lep_tag_truthOrigin:lep_tag_truthType", "%s * (%s)" %(mytree_weight, sel_FAKE_MU_PROBE_L), "text" )

  # Normalise to unit
  if h_FAKE_MU_TAG_T_type_VS_origin.Integral() > 0:
    h_FAKE_MU_TAG_T_type_VS_origin.Scale(1/h_FAKE_MU_TAG_T_type_VS_origin.Integral())
  if h_FAKE_MU_TAG_L_type_VS_origin.Integral() > 0:
    h_FAKE_MU_TAG_L_type_VS_origin.Scale(1/h_FAKE_MU_TAG_L_type_VS_origin.Integral())

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
  
  mytree.Project("FAKE_MU_PROBE_T_ProbeOrigin_VS_TagOrigin", "lep_probe_truthOrigin:lep_tag_truthOrigin", "%s * (%s)" %(mytree_weight, sel_FAKE_MU_PROBE_T), "text" )
  mytree.Project("FAKE_MU_PROBE_L_ProbeOrigin_VS_TagOrigin", "lep_probe_truthOrigin:lep_tag_truthOrigin", "%s * (%s)" %(mytree_weight, sel_FAKE_MU_PROBE_L), "text" )

  h_FAKE_MU_PROBE_T_ProbeOrigin_VS_TagOrigin.Scale(1/h_FAKE_MU_PROBE_T_ProbeOrigin_VS_TagOrigin.Integral())
  h_FAKE_MU_PROBE_L_ProbeOrigin_VS_TagOrigin.Scale(1/h_FAKE_MU_PROBE_L_ProbeOrigin_VS_TagOrigin.Integral())
  """

  # ------------
  #  Electrons
  # ------------
  
  # look only at OF region
  sel_FAKE_EL_PROBE_T = "(" + trig_dec + " && " + lep_tag_trigmatch + " && " + nleptons + " && " + same_sign + " && " +  is_tight_event + " && " + tight_tag + " && " + "( (nmuon == 1 && nel == 1) && ( isProbeElEvent == 1 ) && ( el_probe_isTightSelected[0] == 1 ) )" + " && " + nbjets + " && " + njets + " && "  + ch_flip_veto + " && " + non_prompt + ")"
  sel_FAKE_EL_PROBE_L = "(" + trig_dec + " && " + lep_tag_trigmatch + " && " + nleptons + " && " + same_sign + " && " +  is_tight_event + " && " + tight_tag + " && " + "( (nmuon == 1 && nel == 1) && ( isProbeElEvent == 1 ) && ( el_probe_isTightSelected[0] == 0 ) )" + " && " + nbjets + " && " + njets + " && " + ch_flip_veto + " && " + non_prompt + ")"

  # look only at SF region
  #sel_FAKE_EL_PROBE_T = "(" + trig_dec + " && " + lep_tag_trigmatch + " && " + nleptons + " && " + same_sign + " && " +  is_tight_event + " && " + tight_tag + " && " + "( ( nel == 2 ) && ( isProbeElEvent == 1 ) && ( el_probe_isTightSelected[0] == 1 ) )" + " && " + el_tag_eta + " && " + nbjets + " && " + njets + " && "  + ch_flip_veto + " && " + non_prompt + ")"
  #sel_FAKE_EL_PROBE_L = "(" + trig_dec + " && " + lep_tag_trigmatch + " && " + nleptons + " && " + same_sign + " && " +  is_tight_event + " && " + tight_tag + " && " + "( ( nel == 2 ) && ( isProbeElEvent == 1 ) && ( el_probe_isTightSelected[0] == 0 ) )" + " && " + el_tag_eta + " && " + nbjets + " && " + njets + " && " + ch_flip_veto + " && " + non_prompt + ")"
  
  print "sel_FAKE_EL_PROBE_T: \n", sel_FAKE_EL_PROBE_T
  print "sel_FAKE_EL_PROBE_L: \n", sel_FAKE_EL_PROBE_L
  
  h_FAKE_EL_PROBE_T_type_VS_origin = TH2D("FAKE_EL_PROBE_T_type_VS_origin", "FAKE_EL_PROBE_T_type_VS_origin", 10, 0.0, 10.0, 40, 0.0, 40.0)
  hist_list.append(h_FAKE_EL_PROBE_T_type_VS_origin)
  h_FAKE_EL_PROBE_T_type_VS_origin.GetXaxis().SetTitle("Probe lep (T) truthType")
  h_FAKE_EL_PROBE_T_type_VS_origin.GetYaxis().SetTitle("Probe lep (T) truthOrigin")  
  h_FAKE_EL_PROBE_L_type_VS_origin = TH2D("FAKE_EL_PROBE_L_type_VS_origin", "FAKE_EL_PROBE_L_type_VS_origin", 10, 0.0, 10.0, 40, 0.0, 40.0)
  hist_list.append(h_FAKE_EL_PROBE_L_type_VS_origin)
  h_FAKE_EL_PROBE_L_type_VS_origin.GetXaxis().SetTitle("Probe lep (L!T) truthType")
  h_FAKE_EL_PROBE_L_type_VS_origin.GetYaxis().SetTitle("Probe lep (L!T) truthOrigin")  
  
  h_FAKE_EL_PROBE_T_isChFlip = TH1D("FAKE_EL_PROBE_T_isChFlip", "FAKE_EL_PROBE_T_isChFlip", 2, -0.5, 1.5)
  hist_list.append(h_FAKE_EL_PROBE_T_isChFlip)
  h_FAKE_EL_PROBE_T_isChFlip.GetXaxis().SetTitle("Probe lep (T) isChFlip")
  h_FAKE_EL_PROBE_L_isChFlip = TH1D("FAKE_EL_PROBE_L_isChFlip", "FAKE_EL_PROBE_L_isChFlip", 2, -0.5, 1.5)
  hist_list.append(h_FAKE_EL_PROBE_L_isChFlip)
  h_FAKE_EL_PROBE_L_isChFlip.GetXaxis().SetTitle("Probe lep (L!T) isChFlip")
  
  mytree.Project("FAKE_EL_PROBE_T_type_VS_origin", "lep_probe_truthOrigin:lep_probe_truthType", "%s * (%s)" %(mytree_weight, sel_FAKE_EL_PROBE_T), "text" )
  mytree.Project("FAKE_EL_PROBE_L_type_VS_origin", "lep_probe_truthOrigin:lep_probe_truthType", "%s * (%s)" %(mytree_weight, sel_FAKE_EL_PROBE_L), "text" )
  mytree.Project("FAKE_EL_PROBE_T_isChFlip", "lep_probe_isChFlip", "%s * (%s)" %(mytree_weight, sel_FAKE_EL_PROBE_T), "text" )
  mytree.Project("FAKE_EL_PROBE_L_isChFlip", "lep_probe_isChFlip", "%s * (%s)" %(mytree_weight, sel_FAKE_EL_PROBE_L), "text" )

  # Normalise to unit
  h_FAKE_EL_PROBE_T_type_VS_origin.Scale(1/h_FAKE_EL_PROBE_T_type_VS_origin.Integral())
  h_FAKE_EL_PROBE_L_type_VS_origin.Scale(1/h_FAKE_EL_PROBE_L_type_VS_origin.Integral())

  if h_FAKE_EL_PROBE_T_isChFlip.Integral() > 0 :
    h_FAKE_EL_PROBE_T_isChFlip.Scale(1/h_FAKE_EL_PROBE_T_isChFlip.Integral())
  if h_FAKE_EL_PROBE_L_isChFlip.Integral() > 0 :
    h_FAKE_EL_PROBE_L_isChFlip.Scale(1/h_FAKE_EL_PROBE_L_isChFlip.Integral())

  # -----------------------------------------

  h_FAKE_EL_TAG_T_type_VS_origin = TH2D("FAKE_EL_TAG_T_type_VS_origin", "FAKE_EL_TAG_T_type_VS_origin", 10, 0.0, 10.0, 40, 0.0, 40.0)
  hist_list.append(h_FAKE_EL_TAG_T_type_VS_origin)
  h_FAKE_EL_TAG_T_type_VS_origin.GetXaxis().SetTitle("Tag lep ( Probe (T) ) truthType")
  h_FAKE_EL_TAG_T_type_VS_origin.GetYaxis().SetTitle("Tag lep ( Probe (T) ) truthOrigin")
  h_FAKE_EL_TAG_L_type_VS_origin = TH2D("FAKE_EL_TAG_L_type_VS_origin", "FAKE_EL_TAG_L_type_VS_origin", 10, 0.0, 10.0, 40, 0.0, 40.0)
  hist_list.append(h_FAKE_EL_TAG_L_type_VS_origin)
  h_FAKE_EL_TAG_L_type_VS_origin.GetXaxis().SetTitle("Tag lep ( Probe (L!T) ) truthType")
  h_FAKE_EL_TAG_L_type_VS_origin.GetYaxis().SetTitle("Tag lep ( Probe (L!T) ) truthOrigin")
  
  mytree.Project("FAKE_EL_TAG_T_type_VS_origin", "lep_tag_truthOrigin:lep_tag_truthType", "%s * (%s)" %(mytree_weight, sel_FAKE_EL_PROBE_T), "text" )
  mytree.Project("FAKE_EL_TAG_L_type_VS_origin", "lep_tag_truthOrigin:lep_tag_truthType", "%s * (%s)" %(mytree_weight, sel_FAKE_EL_PROBE_L), "text" )

  # Normalise to unit
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
  
  mytree.Project("FAKE_EL_PROBE_T_ProbeOrigin_VS_TagOrigin", "lep_probe_truthOrigin:lep_tag_truthOrigin", "%s * (%s)" %(mytree_weight, sel_FAKE_EL_PROBE_T), "text" )
  mytree.Project("FAKE_EL_PROBE_L_ProbeOrigin_VS_TagOrigin", "lep_probe_truthOrigin:lep_tag_truthOrigin", "%s * (%s)" %(mytree_weight, sel_FAKE_EL_PROBE_L), "text" )

  h_FAKE_EL_PROBE_T_ProbeOrigin_VS_TagOrigin.Scale(1/h_FAKE_EL_PROBE_T_ProbeOrigin_VS_TagOrigin.Integral())
  h_FAKE_EL_PROBE_L_ProbeOrigin_VS_TagOrigin.Scale(1/h_FAKE_EL_PROBE_L_ProbeOrigin_VS_TagOrigin.Integral())
  """

  # ------------------------------------------------------------------
  
  set_fancy_2D_style()
  
  gStyle.SetPaintTextFormat("2.1f")
  
  for hist in hist_list:
     
     c = TCanvas("c1","Temp",50,50,700,900)

     legend = TLegend(0.2,0.8,0.7,0.85); # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
     legend.SetHeader("2 Lep SS Fake CR")
     legend.SetBorderSize(0)  # no border
     legend.SetFillColor(0)   # Legend background should be white
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
     leg_lumi.DrawLatex(0.2,0.7,"#sqrt{s} = 13 TeV, #int L dt = 3.3 fb^{-1}")

     c.SaveAs(hist.GetName() + ".png")

# ------------------------------------------------------------------

def plot_2D_RealCR():
  
  myfile = TFile("/coepp/cephfs/mel/mmilesi/ttH/MergedDatasets/Merged_v030/Merged_Melb15_ttH_030_DxAOD_p2559/tops/410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.root")
  mytree = gDirectory.Get("physics")
  mytree_weight = "3.209 * mcEventWeight * weight_pileup * weight_lepton_trig_HTop[0] * weight_lepton_reco_HTop[0] * weight_lepton_iso_HTop[0] * weight_lepton_ID_HTop[0] * weight_lepton_TTVA_HTop[0] * weight_jet_JVT_HTop[0] * weight_jet_MV2c20_SFFix77[0]"  
  
  trig_dec          = "( passHLT == 1 )"
  trigmatching      = "( ( lep_isTrigMatched[0] == 1 && ( ( lep_flavour[0] == 11 && lep_pt[0] > 25e3 ) || ( lep_flavour[0] == 13 && lep_pt[0] > 21e3 ) ) ) || ( lep_isTrigMatched[1] == 1 && ( ( lep_flavour[1] == 11 && lep_pt[1] > 25e3 ) || ( lep_flavour[1] == 13 && lep_pt[1] > 21e3 ) ) ) )"
  lep_tag_trigmatch = "( lep_tag_isTightSelected[0] == 1 && lep_tag_isTrigMatched[0] == 1 && ( ( lep_tag_flavour[0] == 11 && lep_tag_pt[0] > 25e3 ) || ( lep_tag_flavour[0] == 13 && lep_tag_pt[0] > 21e3 ) ) )"
  el_tag_eta        = "( nel == 0 || ( nel > 0 && Max$( TMath::Abs(el_caloCluster_eta) ) < 1.37 ) )"
  nbjets            = "( njets_mv2c20_Fix77 > 0 )"
  njets             = "( njets > 1 && njets <= 4 )"
  nleptons          = "( nlep == 2 && Min$( lep_pt ) > 10e3 )"
  flavour_cut       = "( nmuon == 1 && nel == 1 )"
  opposite_sign     = "( isSS01 != 1 )"
  tau_veto          = "( ntau == 0 )"
  zpeak_cut         = "( ( nmuon == 1 && nel == 1 ) || ( TMath::Abs( mll01 - 91.187e3 ) > 7.5e3 ) )"
  truth_cut         = "( isMC==0 || ( isMC==1 && ( ( lep_probe_truthType[0] != 6 && lep_probe_truthType[0] != 2 ) || ( lep_probe_isChFlip[0] != 0 ) ) ) )"
  
  hist_list = []

  basecut =  "(" + trig_dec + " && " + trigmatching + " && " + lep_tag_trigmatch + " && " + el_tag_eta + " && " + njets + " && " + nbjets + " && " + nleptons + " && " + flavour_cut + " && " + opposite_sign + " && " + tau_veto + " && "  + zpeak_cut + " && " + truth_cut + " && " 
  
  # ------------
  #    Muons
  # ------------
   
  sel_REAL_MU_PROBE_T = basecut + "( ( isProbeMuEvent == 1 ) && ( muon_probe_isTightSelected[0] == 1 ) ) )"
  sel_REAL_MU_PROBE_L = basecut + "( ( isProbeMuEvent == 1 ) && ( muon_probe_isTightSelected[0] == 0 ) ) )"
  
  print "sel_REAL_MU_PROBE_T: \n", sel_REAL_MU_PROBE_T
  print "sel_REAL_MU_PROBE_L: \n", sel_REAL_MU_PROBE_L
  
  h_REAL_MU_PROBE_T_type_VS_origin = TH2D("REAL_MU_PROBE_T_type_VS_origin", "REAL_MU_PROBE_T_type_VS_origin", 10, 0.0, 10.0, 40, 0.0, 40.0)
  hist_list.append(h_REAL_MU_PROBE_T_type_VS_origin)
  h_REAL_MU_PROBE_T_type_VS_origin.GetXaxis().SetTitle("Probe lep (T) truthType")
  h_REAL_MU_PROBE_T_type_VS_origin.GetYaxis().SetTitle("Probe lep (T) truthOrigin")
  h_REAL_MU_PROBE_L_type_VS_origin = TH2D("REAL_MU_PROBE_L_type_VS_origin", "REAL_MU_PROBE_L_type_VS_origin", 10, 0.0, 10.0, 40, 0.0, 40.0)
  hist_list.append(h_REAL_MU_PROBE_L_type_VS_origin)
  h_REAL_MU_PROBE_L_type_VS_origin.GetXaxis().SetTitle("Probe lep (L!T) truthType")
  h_REAL_MU_PROBE_L_type_VS_origin.GetYaxis().SetTitle("Probe lep (L!T) truthOrigin")
  
  mytree.Project("REAL_MU_PROBE_T_type_VS_origin", "lep_probe_truthOrigin:lep_probe_truthType", "%s * (%s)" %(mytree_weight, sel_REAL_MU_PROBE_T), "text" )
  mytree.Project("REAL_MU_PROBE_L_type_VS_origin", "lep_probe_truthOrigin:lep_probe_truthType", "%s * (%s)" %(mytree_weight, sel_REAL_MU_PROBE_L), "text" )

  # Normalise to unit
  h_REAL_MU_PROBE_T_type_VS_origin.Scale(1.0/h_REAL_MU_PROBE_T_type_VS_origin.Integral())
  h_REAL_MU_PROBE_L_type_VS_origin.Scale(1.0/h_REAL_MU_PROBE_L_type_VS_origin.Integral())

  # -----------------------------------------

  h_REAL_MU_TAG_T_type_VS_origin = TH2D("REAL_MU_TAG_T_type_VS_origin", "REAL_MU_TAG_T_type_VS_origin", 10, 0.0, 10.0, 40, 0.0, 40.0)
  hist_list.append(h_REAL_MU_TAG_T_type_VS_origin)
  h_REAL_MU_TAG_T_type_VS_origin.GetXaxis().SetTitle("Tag lep ( Probe (T) ) truthType")
  h_REAL_MU_TAG_T_type_VS_origin.GetYaxis().SetTitle("Tag lep ( Probe (T) ) truthOrigin")
  h_REAL_MU_TAG_L_type_VS_origin = TH2D("REAL_MU_TAG_L_type_VS_origin", "REAL_MU_TAG_L_type_VS_origin", 10, 0.0, 10.0, 40, 0.0, 40.0)
  hist_list.append(h_REAL_MU_TAG_L_type_VS_origin)
  h_REAL_MU_TAG_L_type_VS_origin.GetXaxis().SetTitle("Tag lep ( Probe (L!T) ) truthType")
  h_REAL_MU_TAG_L_type_VS_origin.GetYaxis().SetTitle("Tag lep ( Probe (L!T) ) truthOrigin")
  
  mytree.Project("REAL_MU_TAG_T_type_VS_origin", "lep_tag_truthOrigin:lep_tag_truthType", "%s * (%s)" %(mytree_weight, sel_REAL_MU_PROBE_T), "text" )
  mytree.Project("REAL_MU_TAG_L_type_VS_origin", "lep_tag_truthOrigin:lep_tag_truthType", "%s * (%s)" %(mytree_weight, sel_REAL_MU_PROBE_L), "text" )

  # Normalise to unit
  if h_REAL_MU_TAG_T_type_VS_origin.Integral() > 0:
    h_REAL_MU_TAG_T_type_VS_origin.Scale(1.0/h_REAL_MU_TAG_T_type_VS_origin.Integral())
  if h_REAL_MU_TAG_L_type_VS_origin.Integral() > 0:
    h_REAL_MU_TAG_L_type_VS_origin.Scale(1.0/h_REAL_MU_TAG_L_type_VS_origin.Integral())

  # ------------
  #  Electrons
  # ------------
  
  sel_REAL_EL_PROBE_T = basecut + "( ( isProbeElEvent == 1 ) && ( el_probe_isTightSelected[0] == 1 ) ) )"
  sel_REAL_EL_PROBE_L = basecut + "( ( isProbeElEvent == 1 ) && ( el_probe_isTightSelected[0] == 0 ) ) )"
  
  print "sel_REAL_EL_PROBE_T: \n", sel_REAL_EL_PROBE_T
  print "sel_REAL_EL_PROBE_L: \n", sel_REAL_EL_PROBE_L
  
  h_REAL_EL_PROBE_T_type_VS_origin = TH2D("REAL_EL_PROBE_T_type_VS_origin", "REAL_EL_PROBE_T_type_VS_origin", 10, 0.0, 10.0, 40, 0.0, 40.0)
  hist_list.append(h_REAL_EL_PROBE_T_type_VS_origin)
  h_REAL_EL_PROBE_T_type_VS_origin.GetXaxis().SetTitle("Probe lep (T) truthType")
  h_REAL_EL_PROBE_T_type_VS_origin.GetYaxis().SetTitle("Probe lep (T) truthOrigin")  
  h_REAL_EL_PROBE_L_type_VS_origin = TH2D("REAL_EL_PROBE_L_type_VS_origin", "REAL_EL_PROBE_L_type_VS_origin", 10, 0.0, 10.0, 40, 0.0, 40.0)
  hist_list.append(h_REAL_EL_PROBE_L_type_VS_origin)
  h_REAL_EL_PROBE_L_type_VS_origin.GetXaxis().SetTitle("Probe lep (L!T) truthType")
  h_REAL_EL_PROBE_L_type_VS_origin.GetYaxis().SetTitle("Probe lep (L!T) truthOrigin")  
  
  h_REAL_EL_PROBE_T_isChFlip = TH1D("REAL_EL_PROBE_T_isChFlip", "REAL_EL_PROBE_T_isChFlip", 2, -0.5, 1.5)
  hist_list.append(h_REAL_EL_PROBE_T_isChFlip)
  h_REAL_EL_PROBE_T_isChFlip.GetXaxis().SetTitle("Probe lep (T) isChFlip")
  h_REAL_EL_PROBE_L_isChFlip = TH1D("REAL_EL_PROBE_L_isChFlip", "REAL_EL_PROBE_L_isChFlip", 2, -0.5, 1.5)
  hist_list.append(h_REAL_EL_PROBE_L_isChFlip)
  h_REAL_EL_PROBE_L_isChFlip.GetXaxis().SetTitle("Probe lep (L!T) isChFlip")
  
  mytree.Project("REAL_EL_PROBE_T_type_VS_origin", "lep_probe_truthOrigin:lep_probe_truthType", "%s * (%s)" %(mytree_weight, sel_REAL_EL_PROBE_T), "text" )
  mytree.Project("REAL_EL_PROBE_L_type_VS_origin", "lep_probe_truthOrigin:lep_probe_truthType", "%s * (%s)" %(mytree_weight, sel_REAL_EL_PROBE_L), "text" )
  mytree.Project("REAL_EL_PROBE_T_isChFlip", "lep_probe_isChFlip", "%s * (%s)" %(mytree_weight, sel_REAL_EL_PROBE_T), "text" )
  mytree.Project("REAL_EL_PROBE_L_isChFlip", "lep_probe_isChFlip", "%s * (%s)" %(mytree_weight, sel_REAL_EL_PROBE_L), "text" )

  # Normalise to unit
  h_REAL_EL_PROBE_T_type_VS_origin.Scale(1/h_REAL_EL_PROBE_T_type_VS_origin.Integral())
  h_REAL_EL_PROBE_L_type_VS_origin.Scale(1/h_REAL_EL_PROBE_L_type_VS_origin.Integral())

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

  # Normalise to unit
  if h_REAL_EL_TAG_T_type_VS_origin.Integral() >0 :
    h_REAL_EL_TAG_T_type_VS_origin.Scale(1.0/h_REAL_EL_TAG_T_type_VS_origin.Integral())
  if h_REAL_EL_TAG_L_type_VS_origin.Integral() >0 :
    h_REAL_EL_TAG_L_type_VS_origin.Scale(1.0/h_REAL_EL_TAG_L_type_VS_origin.Integral())

  # ------------------------------------------------------------------
  
  set_fancy_2D_style()
  
  gStyle.SetPaintTextFormat("2.1f")
  
  for hist in hist_list:
     
     c = TCanvas("c1","Temp",50,50,700,900)

     legend = TLegend(0.2,0.8,0.7,0.85); # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
     legend.SetHeader("2 Lep OS Real CR")
     legend.SetBorderSize(0)  # no border
     legend.SetFillColor(0)   # Legend background should be white
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
     leg_lumi.DrawLatex(0.2,0.7,"#sqrt{s} = 13 TeV, #int L dt = 3.3 fb^{-1}")

     c.SaveAs(hist.GetName() + ".png")

# -----------------------------------------------------------------------------------

plot_2D_RealCR()  
#plot_2D_FakeCR()  
  
  
  
  
  
  
