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
		     help="Lepton flavour (mu,el) - default is mu")
args = parser.parse_args()


# ---------------------------
# for fancy 2-dim histograms!
# ---------------------------

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


def plot2Dhist():

  observed_sel = ["TT","TL","LT","LL"]
  charges      = ["OS","SS"]
  flavours     = ["ElEl","MuMu"]

  #input_path = "../OutputPlots_MMRates_LHFit_2DPt/"
  input_path = "../OutputPlots_MMRates_LHFit_25ns_v7_FinalSelection_NominalBinning/"

  set_fancy_2D_style()

  for ch in charges:
    for flav in flavours:
      for sl in observed_sel:

        path = input_path + ch + "_" + flav + "_" + sl + "/" + ch + "_" + flav + "_" + sl + "_Lep0Pt_VS_Lep1Pt.root"

        myfile = TFile(path)
        myhist = myfile.Get("observed")

        c = TCanvas("c1","Temp",50,50,1000,800)

        header = None
	if ch == "OS": header = "2 Lep SS Real CR"
	elif ch == "SS": header = "2 Lep SS Fake CR"

     	legend = TLegend(0.2,0.8,0.4,0.85); # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
     	legend.SetHeader(header)
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

     	myhist.Draw("colz text")

        legend.Draw()
        leg_ATLAS.DrawLatex(0.2,0.75,"#bf{#it{ATLAS}} Work In Progress")
        leg_lumi.DrawLatex(0.2,0.7,"#sqrt{s} = 13 TeV, #int L dt = 3.2 fb^{-1}")

        canvasname = ch + "_" + flav + "_" + sl + "_observed"
        c.SaveAs(canvasname + ".png")

# ---------------
# Real efficiency
# ---------------

def plotter_r():

  #file_TP_path = "../OutputPlots_MMRates_25ns_v7_FinalSelection_LHComparison/Rates_NoSub_LHInput_6x6bins/Rates.root"

  #file_TP_path = "../OutputPlots_MMRates_25ns_v7_FinalSelection_NominalBinning/Rates_NoSub_LHInput/Rates.root"
  file_TP_path = "../OutputPlots_MMRates_25ns_v7_FinalSelection_NominalBinning/Rates_YesSub_LHInput/Rates.root"

  file_TP = TFile(file_TP_path)

  hist_TP = ( file_TP.Get("El_ProbePt_Real_Efficiency_observed"), file_TP.Get("Mu_ProbePt_Real_Efficiency_observed") )[bool(args.flavour == "mu")]

  print("Reading histogram {0} from file {1}".format(hist_TP.GetName(), file_TP_path))

  hist_TP.SetLineWidth(2)
  hist_TP.SetMarkerSize(1.0)

  hist_TP.GetYaxis().SetTitle("efficiency")

  hist_TP.GetXaxis().SetTitleOffset(1.0)
  hist_TP.GetYaxis().SetTitleOffset(1.0)

  hist_TP.GetXaxis().SetRangeUser(0.4,1.0)
  hist_TP.GetYaxis().SetRangeUser(0.4,1.0)

  # -------------------------------------------

  #LH_init_path = "../OutputPlots_MMRates_LHFit_2DPt/6x6bins/"

  #LH_init_path = "../OutputPlots_MMRates_LHFit_25ns_v7_FinalSelection_NominalBinning/FittedEfficiencies_NoSub/"
  LH_init_path = "../OutputPlots_MMRates_LHFit_25ns_v7_FinalSelection_NominalBinning/FittedEfficiencies_YesSub/"

  file_LH_path = ( LH_init_path + "LH_efficiencies_real_el.root", LH_init_path + "LH_efficiencies_real_mu.root" )[bool(args.flavour == "mu")]

  file_LH = TFile(file_LH_path)

  hist_LH = file_LH.Get("r_hist")

  print("Reading histogram {0} from file {1}".format(hist_LH.GetName(), file_LH_path))

  hist_LH.SetLineWidth(2)
  hist_LH.SetMarkerSize(1.0)

  hist_LH.GetYaxis().SetTitle("efficiency")

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
  legend.SetHeader("#epsilon_{real}")
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

  canvasname = ( "real_eff_el_TP_LH", "real_eff_mu_TP_LH" )[bool(args.flavour == "mu")]
  c.SaveAs( canvasname + ".png" )

# ---------------
# Fake efficiency
# ---------------

def plotter_f():

  #file_TP_path = "../OutputPlots_MMRates_25ns_v7_FinalSelection_LHComparison/Rates_NoSub_LHInput_3x3bins/Rates.root"

  #file_TP_path = "../OutputPlots_MMRates_25ns_v7_FinalSelection_NominalBinning/Rates_NoSub_LHInput/Rates.root"
  file_TP_path = "../OutputPlots_MMRates_25ns_v7_FinalSelection_NominalBinning/Rates_YesSub_LHInput/Rates.root"

  file_TP = TFile(file_TP_path)

  hist_TP = ( file_TP.Get("El_ProbePt_Fake_Efficiency_observed"), file_TP.Get("Mu_ProbePt_Fake_Efficiency_observed") )[bool(args.flavour == "mu")]

  print("Reading histogram {0} from file {1}".format(hist_TP.GetName(), file_TP_path))

  hist_TP.SetLineWidth(2)
  hist_TP.SetMarkerSize(1.0)

  hist_TP.GetYaxis().SetTitle("efficiency")

  hist_TP.GetXaxis().SetTitleOffset(1.0)
  hist_TP.GetYaxis().SetTitleOffset(1.0)

  hist_TP.GetXaxis().SetRangeUser(0.0,1.0)
  hist_TP.GetYaxis().SetRangeUser(0.0,1.0)

  # -------------------------------------------

  #LH_init_path = "../OutputPlots_MMRates_LHFit_2DPt/3x3bins/"

  #LH_init_path = "../OutputPlots_MMRates_LHFit_25ns_v7_FinalSelection_NominalBinning/FittedEfficiencies_NoSub/"
  LH_init_path = "../OutputPlots_MMRates_LHFit_25ns_v7_FinalSelection_NominalBinning/FittedEfficiencies_YesSub/"

  file_LH_path = ( LH_init_path + "LH_efficiencies_fake_el.root", LH_init_path + "LH_efficiencies_fake_mu.root" )[bool(args.flavour == "mu")]

  file_LH = TFile(file_LH_path)

  hist_LH = file_LH.Get("f_hist")

  print("Reading histogram {0} from file {1}".format(hist_LH.GetName(), file_LH_path))


  hist_LH.SetLineWidth(2)
  hist_LH.SetMarkerSize(1.0)

  hist_LH.GetYaxis().SetTitle("efficiency")

  hist_LH.GetXaxis().SetTitleOffset(1.0)
  hist_LH.GetYaxis().SetTitleOffset(1.0)

  hist_LH.GetXaxis().SetRangeUser(0.0,1.0)
  hist_LH.GetYaxis().SetRangeUser(0.0,1.0)

  # -------------------------------------------

  hist_TP.SetLineColor(kRed)
  hist_TP.SetMarkerColor(kRed)

  hist_LH.SetLineColor(kBlue)
  hist_LH.SetMarkerColor(kBlue)

  # -------------------------------------------

  c = TCanvas("c1","Temp",50,50,1000,900)

  legend = TLegend(0.45,0.5,0.925,0.8) # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
  legend.SetHeader("#epsilon_{fake}")
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

  canvasname = ( "fake_eff_el_TP_LH", "fake_eff_mu_TP_LH" )[bool(args.flavour == "mu")]
  c.SaveAs( canvasname + ".png" )


# ----------------

#plotter_r()
#plotter_f()
plot2Dhist()
