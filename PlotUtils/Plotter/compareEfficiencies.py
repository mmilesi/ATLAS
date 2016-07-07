#!/usr/bin/python

import array
import os
import sys

#from array import array

from ROOT import gROOT, gDirectory, gStyle, gPad, TH1D, TH2D, TFile, TCanvas, TColor, TLegend, TLatex, kRed, kBlue, kAzure, kCyan

gROOT.Reset()
gROOT.LoadMacro("$HOME/RootUtils/AtlasStyle.C")
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

  #file_TP_path = "$HOME/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMRates_25ns_v7_FinalSelection_NominalBinning/Rates_NoSub_LHInput/Rates.root"
  file_TP_path = "$HOME/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMRates_25ns_v7_FinalSelection_NominalBinning/Rates_YesSub_LHInput/Rates.root"
  # CLOSURE
  #file_TP_path = "$HOME/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMClosureRates_25ns_v7_FinalSelection_NominalBinning/Rates_YesSub_LHInput/Rates.root"

  file_TP = TFile(file_TP_path)

  hist_TP = ( file_TP.Get("El_ProbePt_Real_Efficiency_observed"), file_TP.Get("Mu_ProbePt_Real_Efficiency_observed") )[bool(args.flavour == "mu")]
  # CLOSURE
  #hist_TP = ( file_TP.Get("El_ProbePt_Real_Efficiency_expected"), file_TP.Get("Mu_ProbePt_Real_Efficiency_expected") )[bool(args.flavour == "mu")]

  print("Reading histogram {0} from file {1}".format(hist_TP.GetName(), file_TP_path))

  hist_TP.SetLineWidth(2)
  hist_TP.SetMarkerSize(1.0)

  hist_TP.GetYaxis().SetTitle("efficiency")

  hist_TP.GetXaxis().SetTitleOffset(1.0)
  hist_TP.GetYaxis().SetTitleOffset(1.0)

  hist_TP.GetXaxis().SetRangeUser(0.4,1.0)
  hist_TP.GetYaxis().SetRangeUser(0.4,1.0)

  # -------------------------------------------

  #LH_init_path = "$HOME/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMRates_LHFit_25ns_v7_FinalSelection_NominalBinning/FittedEfficiencies_NoSub/"
  LH_init_path = "$HOME/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMRates_LHFit_25ns_v7_FinalSelection_NominalBinning/FittedEfficiencies_YesSub/"
  # CLOSURE
  #LH_init_path = "$HOME/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMClosureRates_LHFit_25ns_v7_FinalSelection_NominalBinning/FittedEfficiencies_YesSub/"

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

  #file_TP_path = "$HOME/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMRates_25ns_v7_FinalSelection_NominalBinning/Rates_NoSub_LHInput/Rates.root"
  file_TP_path = "$HOME/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMRates_25ns_v7_FinalSelection_NominalBinning/Rates_YesSub_LHInput/Rates.root"
  # CLOSURE
  #file_TP_path = "$HOME/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMClosureRates_25ns_v7_FinalSelection_NominalBinning/Rates_YesSub_LHInput/Rates.root"

  file_TP = TFile(file_TP_path)

  hist_TP = ( file_TP.Get("El_ProbePt_Fake_Efficiency_observed"), file_TP.Get("Mu_ProbePt_Fake_Efficiency_observed") )[bool(args.flavour == "mu")]
  # CLOSURE
  #hist_TP = ( file_TP.Get("El_ProbePt_Fake_Efficiency_expected"), file_TP.Get("Mu_ProbePt_Fake_Efficiency_expected") )[bool(args.flavour == "mu")]

  print("Reading histogram {0} from file {1}".format(hist_TP.GetName(), file_TP_path))

  hist_TP.SetLineWidth(2)
  hist_TP.SetMarkerSize(1.0)

  hist_TP.GetYaxis().SetTitle("efficiency")

  hist_TP.GetXaxis().SetTitleOffset(1.0)
  hist_TP.GetYaxis().SetTitleOffset(1.0)

  hist_TP.GetXaxis().SetRangeUser(0.0,1.0)
  hist_TP.GetYaxis().SetRangeUser(0.0,1.0)

  # -------------------------------------------

  #LH_init_path = "$HOME/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMRates_LHFit_25ns_v7_FinalSelection_NominalBinning/FittedEfficiencies_NoSub/"
  LH_init_path = "$HOME/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMRates_LHFit_25ns_v7_FinalSelection_NominalBinning/FittedEfficiencies_YesSub/"
  # CLOSURE
  #LH_init_path = "$HOME/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMClosureRates_LHFit_25ns_v7_FinalSelection_NominalBinning/FittedEfficiencies_YesSub/"

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

def plotter_r_flavours():

  #file_TP_path = "$HOME/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMRates_25ns_v7_FinalSelection_NominalBinning/Rates_NoSub_LHInput/Rates.root"
  #file_TP_path = "./OutputPlots_MMRates_25ns_v7_FinalSelection_NominalBinning/Rates_YesSub_LHInput/Rates.root"
  
  # CLOSURE
  #file_TP_path = "$HOME/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMClosureRates_25ns_v7_FinalSelection_NominalBinning/Rates_YesSub_LHInput/Rates.root"
  file_TP_path = "$HOME/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMClosureRates_25ns_v14_DLT_2015/Rates.root"

  file_TP = TFile(file_TP_path)

  #hist_TP = ( file_TP.Get("El_ProbePt_Real_Efficiency_observed"), file_TP.Get("Mu_ProbePt_Real_Efficiency_observed") )[bool(args.flavour == "mu")]
  # CLOSURE
  hist_TP = ( file_TP.Get("El_ProbePt_Real_Efficiency_expected"), file_TP.Get("Mu_ProbePt_Real_Efficiency_expected") )[bool(args.flavour == "mu")]

  print("Reading histogram {0} from file {1}".format(hist_TP.GetName(), file_TP_path))

  hist_TP.SetLineWidth(2)
  hist_TP.SetMarkerSize(1.0)

  hist_TP.GetXaxis().SetTitleOffset(1.0)
  hist_TP.GetYaxis().SetTitleOffset(1.0)

  lepton = ("e","#mu")[bool(args.flavour == "mu")]

  hist_TP.GetXaxis().SetTitle("p_{T}^{"+ lepton +"} [GeV]")
  hist_TP.GetYaxis().SetTitle("efficiency")

  if args.flavour == "el":
    hist_TP.GetYaxis().SetRangeUser(0.4,1.0)
  elif args.flavour == "mu":
    hist_TP.GetYaxis().SetRangeUser(0.65,1.0)

  hist_TP.SetLineColor(kRed)
  hist_TP.SetMarkerColor(kRed)

  #c = TCanvas("c1","Temp",50,50,1000,900)
  c = TCanvas("c1","Temp",50,50,1300,800)

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
  legend.AddEntry(hist_TP, "Tag & Probe", "P")

  # -------------------------------------------

  #LH_init_path = "$HOME/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMRates_LHFit_25ns_v7_FinalSelection_NominalBinning/FittedEfficiencies_YesSub_DiffFlavours/"
  LH_init_path = "$HOME/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMClosureRates_LHFit_25ns_v14_DLT_2015/LHClosureRates_25ns_v14_DLT_2015/"
  
  hist_LH_list = []
  #flav_comp_list = ["mumu","elel","of","incl"]
  #flav_comp_list = ["mumu"]
  flav_comp_list = ["incl"]

  if ( args.flavour == "mu" and "elel" in flav_comp_list ): flav_comp_list.remove("elel")
  if ( args.flavour == "el" and "mumu" in flav_comp_list ): flav_comp_list.remove("mumu")

  for idx, flavcomp in enumerate(flav_comp_list,start=0):

    file_LH_path = LH_init_path + "Fit_" + flavcomp + "/LH_efficiencies_real_" + args.flavour + "_" + flavcomp + ".root"

    file_LH = TFile(file_LH_path)

    hist_LH = file_LH.Get("r_hist")

    print("Reading histogram {0} from file {1} - index : {2}".format(hist_LH.GetName(), file_LH_path, idx))

    # A dirty hack to make multiple plots superimposed with a small offset

    offsets = []
    for ibin in range(0,hist_LH.GetNbinsX()+1):
      offset = ( idx + 1 ) * 0.15 * hist_LH.GetBinWidth(ibin)
      offsets.append(offset)

    print("offsets per bin: {0}".format(offsets))

    binlist = []
    for ibin in range(0,hist_LH.GetNbinsX()+1):
      binlist.append( hist_LH.GetXaxis().GetBinUpEdge(ibin) + offsets[ibin] )

    print("shifted bins: {0}".format(binlist))

    binlist_arr = array.array("d",binlist)

    shiftedhist = TH1D( hist_LH.GetName() + "_shift", hist_LH.GetTitle(), len(binlist)-1, binlist_arr )

    for ibin in range(1, hist_LH.GetNbinsX()+1):
      shiftedhist.SetBinContent( ibin, hist_LH.GetBinContent(ibin) )
      shiftedhist.SetBinError( ibin, hist_LH.GetBinError(ibin) )

    shiftedhist.SetLineWidth(2)
    shiftedhist.SetMarkerSize(1.0)

    shiftedhist.GetXaxis().SetTitle("p_{T}^{"+ lepton +"} [GeV]")
    shiftedhist.GetYaxis().SetTitle("efficiency")

    shiftedhist.GetXaxis().SetTitleOffset(1.0)
    shiftedhist.GetYaxis().SetTitleOffset(1.0)

    if args.flavour == "el":
      shiftedhist.GetYaxis().SetRangeUser(0.4,1.0)
    elif args.flavour == "mu":
      shiftedhist.GetYaxis().SetRangeUser(0.65,1.0)

    if idx == 0:
      shiftedhist.SetLineColor(kBlue)
      shiftedhist.SetMarkerColor(kBlue)
    elif idx == 1:
      shiftedhist.SetLineColor(kBlue+3)
      shiftedhist.SetMarkerColor(kBlue+3)
    elif idx == 2:
      shiftedhist.SetLineColor(kAzure+1)
      shiftedhist.SetMarkerColor(kAzure+1)
    elif idx == 3:
      shiftedhist.SetLineColor(kCyan)
      shiftedhist.SetMarkerColor(kCyan)

    shiftedhist.SetDirectory(0)

    pair = (flavcomp,shiftedhist)
    hist_LH_list.append(pair)

  for idx, histpair in enumerate(hist_LH_list,start=0):
    #if idx == 0:
    #  histpair[1].Draw("E0")
    #else:
    #  histpair[1].Draw("E0,SAME")
    histpair[1].Draw("E0,SAME")

    legend.AddEntry(histpair[1], "Likelihood - " +  histpair[0], "P")
    legend.Draw()

  leg_ATLAS.DrawLatex(0.6,0.35,"#bf{#it{ATLAS}} Work In Progress")
  leg_lumi.DrawLatex(0.6,0.27,"#sqrt{s} = 13 TeV, #int L dt = 3.2 fb^{-1}")

  canvasname = ( "real_eff_el_TP_LH", "real_eff_mu_TP_LH" )[bool(args.flavour == "mu")]
  c.SaveAs( canvasname + ".png" )

# ----------------

def plotter_f_flavours():

  #file_TP_path = "$HOME/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMRates_25ns_v7_FinalSelection_NominalBinning/Rates_NoSub_LHInput/Rates.root"
  #file_TP_path = "./OutputPlots_MMRates_25ns_v7_FinalSelection_NominalBinning/Rates_YesSub_LHInput/Rates.root"
  
  #file_TP_path = "./OutputPlots_MMRates_25ns_v7_FinalSelection_DDQmisID/Rates.root"
 
  # CLOSURE
  #file_TP_path = "$HOME/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMClosureRates_25ns_v7_FinalSelection_NominalBinning/Rates_YesSub_LHInput/Rates.root"
  file_TP_path = "$HOME/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMClosureRates_25ns_v14_DLT_2015/Rates.root"

  file_TP = TFile(file_TP_path)

  #hist_TP = ( file_TP.Get("El_ProbePt_Fake_Efficiency_observed"), file_TP.Get("Mu_ProbePt_Fake_Efficiency_observed") )[bool(args.flavour == "mu")]
  # CLOSURE
  hist_TP = ( file_TP.Get("El_ProbePt_Fake_Efficiency_expected"), file_TP.Get("Mu_ProbePt_Fake_Efficiency_expected") )[bool(args.flavour == "mu")]

  print("Reading histogram {0} from file {1}".format(hist_TP.GetName(), file_TP_path))

  hist_TP.SetLineWidth(2)
  hist_TP.SetMarkerSize(1.0)

  hist_TP.GetXaxis().SetTitleOffset(1.0)
  hist_TP.GetYaxis().SetTitleOffset(1.0)

  lepton = ("e","#mu")[bool(args.flavour == "mu")]

  hist_TP.GetXaxis().SetTitle("p_{T}^{"+ lepton +"} [GeV]")
  hist_TP.GetYaxis().SetTitle("efficiency")

  hist_TP.GetYaxis().SetRangeUser(0.0,1.0)

  hist_TP.SetLineColor(kRed)
  hist_TP.SetMarkerColor(kRed)

  #c = TCanvas("c1","Temp",50,50,1000,900)
  c = TCanvas("c1","Temp",50,50,1300,800)

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
  legend.AddEntry(hist_TP, "Tag & Probe", "P")

  # -------------------------------------------

  #LH_init_path = "$HOME/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMRates_LHFit_25ns_v7_FinalSelection_NominalBinning/FittedEfficiencies_YesSub_DiffFlavours/"
  LH_init_path = "$HOME/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMClosureRates_LHFit_25ns_v14_DLT_2015/LHClosureRates_25ns_v14_DLT_2015/"

  hist_LH_list = []
  #flav_comp_list = ["mumu","elel","of","incl"]
  #flav_comp_list = ["mumu"]
  flav_comp_list = ["incl"]

  if ( args.flavour == "mu" and "elel" in flav_comp_list ): flav_comp_list.remove("elel")
  if ( args.flavour == "el" and "mumu" in flav_comp_list ): flav_comp_list.remove("mumu")

  for idx, flavcomp in enumerate(flav_comp_list,start=0):

    file_LH_path = LH_init_path + "Fit_" + flavcomp + "/LH_efficiencies_fake_" + args.flavour + "_" + flavcomp + ".root"

    file_LH = TFile(file_LH_path)

    hist_LH = file_LH.Get("f_hist")

    print("Reading histogram {0} from file {1} - index : {2}".format(hist_LH.GetName(), file_LH_path, idx))

    # A dirty hack to make multiple plots superimposed with a small offset

    offsets = []
    for ibin in range(0,hist_LH.GetNbinsX()+1):
      offset = ( idx + 1 ) * 0.15 * hist_LH.GetBinWidth(ibin)
      # hardcoded...so ugly, any better ideas?
      if ibin == 4:
        offset = ( idx + 1 ) * 0.05 * hist_LH.GetBinWidth(ibin)
      if ibin == 5:
        offset = ( idx + 1 ) * 0.01 * hist_LH.GetBinWidth(ibin)
      offsets.append(offset)

    print("offsets per bin: {0}".format(offsets))

    binlist = []
    for ibin in range(0,hist_LH.GetNbinsX()+1):
      binlist.append( hist_LH.GetXaxis().GetBinUpEdge(ibin) + offsets[ibin] )

    print("shifted bins: {0}".format(binlist))

    binlist_arr = array.array("d",binlist)

    shiftedhist = TH1D( hist_LH.GetName() + "_shift", hist_LH.GetTitle(), len(binlist)-1, binlist_arr )

    for ibin in range(1, hist_LH.GetNbinsX()+1):
      shiftedhist.SetBinContent( ibin, hist_LH.GetBinContent(ibin) )
      shiftedhist.SetBinError( ibin, hist_LH.GetBinError(ibin) )

    shiftedhist.SetLineWidth(2)
    shiftedhist.SetMarkerSize(1.0)

    lepton = ("e","#mu")[bool(args.flavour == "mu")]
    shiftedhist.GetXaxis().SetTitle("p_{T}^{"+ lepton +"} [GeV]")
    shiftedhist.GetYaxis().SetTitle("efficiency")

    shiftedhist.GetXaxis().SetTitleOffset(1.0)
    shiftedhist.GetYaxis().SetTitleOffset(1.0)

    shiftedhist.GetYaxis().SetRangeUser(0.0,1.0)

    if idx == 0:
      shiftedhist.SetLineColor(kBlue)
      shiftedhist.SetMarkerColor(kBlue)
    elif idx == 1:
      shiftedhist.SetLineColor(kBlue+3)
      shiftedhist.SetMarkerColor(kBlue+3)
    elif idx == 2:
      shiftedhist.SetLineColor(kAzure+1)
      shiftedhist.SetMarkerColor(kAzure+1)
    elif idx == 3:
      shiftedhist.SetLineColor(kCyan)
      shiftedhist.SetMarkerColor(kCyan)

    shiftedhist.SetDirectory(0)

    pair = (flavcomp,shiftedhist)
    hist_LH_list.append(pair)

  gPad.SetLogx()
  hist_TP.GetXaxis().SetMoreLogLabels()

  for idx, histpair in enumerate(hist_LH_list,start=0):
    histpair[1].GetXaxis().SetMoreLogLabels()
    #if idx == 0:
    #  histpair[1].Draw("E0")
    #else:
    #  histpair[1].Draw("E0,SAME")
    histpair[1].Draw("E0,SAME")

    legend.AddEntry(histpair[1], "Likelihood - " +  histpair[0], "P")
    legend.Draw()

  leg_ATLAS.DrawLatex(0.6,0.35,"#bf{#it{ATLAS}} Work In Progress")
  leg_lumi.DrawLatex(0.6,0.27,"#sqrt{s} = 13 TeV, #int L dt = 3.2 fb^{-1}")

  canvasname = ( "fake_eff_el_TP_LH", "fake_eff_mu_TP_LH" )[bool(args.flavour == "mu")]
  c.SaveAs( canvasname + ".png" )

# ----------------

#plotter_r()
#plotter_f()
plotter_r_flavours()
plotter_f_flavours()
#plot2Dhist()
