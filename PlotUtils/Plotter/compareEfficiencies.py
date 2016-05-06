#!/usr/bin/python

import array
import os
import sys

from array import array

from ROOT import gROOT, gDirectory, gStyle, TH1D, TH2D, TFile, TCanvas, TColor, TLegend, TLatex, kRed, kBlue

gROOT.Reset()
gROOT.LoadMacro("AtlasStyle.C")
from ROOT import SetAtlasStyle
SetAtlasStyle()

gROOT.SetBatch(True)


import argparse
parser = argparse.ArgumentParser(description="Compare efficiencies measured with T&P and Likelihood fit")
parser.add_argument("--flavour", metavar="FLAVOUR", dest="flavour", default="mu", type=str,
		     help="Lepton flavour (\"mu\",\"el\") - default is \"mu\"")
args = parser.parse_args()

def plotter():
  
  file_TP = TFile("../OutputPlots_MMRates_25ns_v7_FinalSelection_LHComparison/Rates_NoSub_LHInput/Rates.root")
  
  hist_TP = ( file_TP.Get("El_ProbePt_Real_Efficiency_observed"), file_TP.Get("Mu_ProbePt_Real_Efficiency_observed") )[bool(args.flavour is "mu")]

  hist_TP.SetLineWidth(2)
  hist_TP.SetMarkerSize(1.0)
 
  hist_TP.GetYaxis().SetTitle("Real efficiency")

  hist_TP.GetXaxis().SetTitleOffset(1.0)
  hist_TP.GetYaxis().SetTitleOffset(1.0)

  hist_TP.GetXaxis().SetRangeUser(0.4,1.0)
  hist_TP.GetYaxis().SetRangeUser(0.4,1.0)

  # -------------------------------------------

  file_LH = ( TFile("../OutputPlots_MMRates_LHFit_2DPt/LH_efficiencies_el.root"), TFile("../OutputPlots_MMRates_LHFit_2DPt/LH_efficiencies_mu.root") )[bool(args.flavour is "mu")]

  hist_LH = file_LH.Get("r_hist")
  
  hist_LH.SetLineWidth(2)
  hist_LH.SetMarkerSize(1.0)

  hist_LH.GetYaxis().SetTitle("Real efficiency")

  hist_LH.GetXaxis().SetTitleOffset(1.0)
  hist_LH.GetYaxis().SetTitleOffset(1.0)

  hist_LH.GetXaxis().SetRangeUser(0.4,1.0)
  hist_LH.GetYaxis().SetRangeUser(0.4,1.0)

  # -------------------------------------------

  hist_TP.SetLineColor(kRed)
  hist_TP.SetMarkerColor(kRed)
  
  hist_LH.SetLineColor(kBlue)
  hist_LH.SetMarkerColor(kBlue)

  # -------------------------------------------

  c = TCanvas("c1","Temp",50,50,1000,900)
     
  legend = TLegend(0.45,0.5,0.925,0.8) # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
  legend.SetHeader("Real effficiency")
  legend.SetBorderSize(0)  # no border
  legend.SetFillStyle(0) # Legend transparent background
  legend.SetTextSize(0.04) # Increase entry font size!
  legend.SetTextFont(42)   # Helvetica
  leg_ATLAS = TLatex()
  leg_lumi  = TLatex()
  leg_ATLAS.SetTextSize(0.03)
  leg_ATLAS.SetNDC()
  leg_lumi.SetTextSize(0.03)
  leg_lumi.SetNDC()

  hist_TP.Draw("E0")
  hist_LH.Draw("E0,SAME")
	  
  legend.AddEntry(hist_TP, "Tag & Probe", "P")
  legend.AddEntry(hist_LH, "Likelihood", "P")
     
  legend.Draw()
  leg_ATLAS.DrawLatex(0.6,0.35,"#bf{#it{ATLAS}} Work In Progress")
  leg_lumi.DrawLatex(0.6,0.27,"#sqrt{s} = 13 TeV, #int L dt = 3.2 fb^{-1}")

  canvasname = ( "real_eff_el_TP_LH", "real_eff_mu_TP_LH" )[bool(args.flavour is "mu")]
  c.SaveAs( canvasname + ".png" )

# ----------------

plotter()
