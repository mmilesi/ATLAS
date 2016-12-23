#!/usr/bin/python

import array
import os
import sys

import argparse

list_flavCRs = ["mumu","elel","of","incl"]

luminosities = { "GRL v73 - Moriond 2016 GRL":3.209,  # March 2016
                 "ICHEP 2015+2016 DS":13.20768,       # August 2016
                 "POST-ICHEP 2015+2016 DS":22.07036,  # October 2016
                 "FULL 2015+2016 DS":36.4702          # December 2016
               }

triggers = ["SLT","DLT"]

parser = argparse.ArgumentParser(description="Compare efficiencies measured with T&P and Likelihood fit")
parser.add_argument("--flavour", metavar="FLAVOUR", dest="flavour", default="mu", type=str,
		     help="Lepton flavour for the efficiency (mu,el) - default is mu")
parser.add_argument("--flavRealCR", dest="flavRealCR", action="store", default="incl", type=str, nargs="+",
                    help="The flavour composition of the CR where REAL efficiencies were measured. Use space-separated args. Full list of available options:\n{0}".format(list_flavCRs))
parser.add_argument("--flavFakeCR", dest="flavFakeCR", action="store", default="incl", type=str, nargs="+",
                    help="The flavour composition of the CR where FAKE efficiencies were measured. Use space-separated args. Full list of available options:\n{0}".format(list_flavCRs))
parser.add_argument("--closure", dest="closure", action="store_true", default=False,
                    help="Run on ttbar to perform MM closure test. Default is False, i.e., the code will run on data.")
parser.add_argument("--lumi", dest="lumi", action="store", type=float, default=luminosities["FULL 2015+2016 DS"],
                    help="The luminosity of the dataset. Pick one of these values: ==> " + ",".join( "{0} ({1})".format( lumi, tag ) for tag, lumi in luminosities.iteritems() ) + ". Default is {0}".format(luminosities["FULL 2015+2016 DS"] ) )
parser.add_argument('--trigger', dest='trigger', action='store', default=triggers[0], type=str, nargs='+',
                    help='The trigger strategy to be used. Choose one of:\n{0}.\nIf this option is not specified, default is {1}'.format(triggers,triggers[0]))
parser.add_argument("--triggerEff", dest="triggerEff", action="store_true", default=False,
                    help="Make comparison plots for trigger efficiency.")

args = parser.parse_args()

from ROOT import gROOT, gDirectory, gStyle, gPad, TPad, TH1, TH1D, TH2D, TFile, TCanvas, TColor, TLegend, TLatex, TLine, kRed, kBlue, kAzure, kCyan, kBlack, kMagenta

gROOT.Reset()
gROOT.LoadMacro("$HOME/RootUtils/AtlasStyle.C")
from ROOT import SetAtlasStyle
SetAtlasStyle()

gROOT.SetBatch(True)

TH1.SetDefaultSumw2()

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


def plotter_flavours( eff_type ):

  if eff_type == "real":
      hist_eff_type = "Real"
  elif eff_type == "fake":
      hist_eff_type = "Fake"

  basepath = "$HOME/PhD/ttH_MultiLeptons/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/"

  # ---------------------------
  # HTop Tag & Probe efficiency
  # ---------------------------

  file_TP_path = basepath + "blahblah"

  if args.closure:
      if "SLT" in args.trigger:
          file_TP_path = basepath + "PLOTS_25ns_v23/MMClosure_v23_HTopTP/OutputPlots_MMClosureRates_HTopTagProbe_NoCorr_SLT_SFmuSFel_25ns_v23/LeptonEfficiencies.root"
      elif "DLT" in args.trigger:
          #file_TP_path = basepath + "PLOTS_25ns_v23/MMClosure_v23_HTopTP/OutputPlots_MMClosureRates_HTopTagProbe_NoCorr_DLT_SFmuSFel_25ns_v23/LeptonEfficiencies.root"
          file_TP_path = basepath + "PLOTS_25ns_v23/MMClosure_v23_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_DLT_SFmuSFel_Pt_AMBISOLVING_25ns_v23/LeptonEfficiencies.root"

  print("------------------------------------------------------------------------------------------")
  print("HTop Tag & Probe efficiency - Opening file:\n{0}".format(file_TP_path))
  print("------------------------------------------------------------------------------------------")

  file_TP = TFile(file_TP_path)

  if args.closure:
      hist_TP = ( file_TP.Get(hist_eff_type+"_El_Pt_Efficiency_expectedbkg"), file_TP.Get(hist_eff_type+"_Mu_Pt_Efficiency_expectedbkg") )[bool(args.flavour == "mu")]
  else:
      hist_TP = ( file_TP.Get(hist_eff_type+"_El_Pt_Efficiency_observed_sub"), file_TP.Get(hist_eff_type+"_Mu_Pt_Efficiency_observed_sub") )[bool(args.flavour == "mu")]

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
  legend.SetHeader("#epsilon_{" + eff_type + "}")
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

  #hist_TP.Draw("E0")
  #legend.AddEntry(hist_TP, "HTop Tag & Probe (ICHEP)", "P")

  # ---------------------------
  # New Tag & Probe efficiency
  # ---------------------------

  file_SUSY_TP_path = basepath + "blahblah"

  if args.closure:
      if "SLT" in args.trigger:
          #file_SUSY_TP_path = basepath + "PLOTS_25ns_v23/MMClosure_v23_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_SFmuSFel_Pt_AMBISOLVING_25ns_v23/LeptonEfficiencies.root"
          #file_SUSY_TP_path = basepath + "PLOTS_25ns_v23/MMClosure_v23_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_OFmuOFel_Pt_AMBISOLVING_25ns_v23/LeptonEfficiencies.root"
          #file_SUSY_TP_path = basepath + "PLOTS_25ns_v23/MMClosure_v23_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_SFmuSFel_massClosestBJet_AMBISOLVING_25ns_v23/LeptonEfficiencies.root"
          #file_SUSY_TP_path = basepath + "PLOTS_25ns_v23/MMClosure_v23_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_SFmuSFel_deltaRClosestBJet_AMBISOLVING_25ns_v23/LeptonEfficiencies.root"
          #file_SUSY_TP_path = basepath + "PLOTS_25ns_v23/MMClosure_v23_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v23/LeptonEfficiencies.root"
          #file_SUSY_TP_path = basepath + "PLOTS_25ns_v23/MMClosure_v23_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v23_OLD_EL_BINNING/LeptonEfficiencies.root"
          #
          #file_SUSY_TP_path = basepath + "PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24/LeptonEfficiencies.root"
          #file_SUSY_TP_path = basepath + "PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_ForceProbeToBeFake/LeptonEfficiencies.root"
          #file_SUSY_TP_path = basepath + "PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIso/LeptonEfficiencies.root"
          #file_SUSY_TP_path = basepath + "PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIsoTagPt30/LeptonEfficiencies.root"
          file_SUSY_TP_path = basepath + "PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIsoTagPt30_ForceProbeToBeFake/LeptonEfficiencies.root"
      elif "DLT" in args.trigger:
          #file_SUSY_TP_path = basepath + "PLOTS_25ns_v23/MMClosure_v23_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_DLT_SFmuSFel_Pt_AMBISOLVING_25ns_v23/LeptonEfficiencies.root"
          #file_SUSY_TP_path = basepath + "PLOTS_25ns_v23/MMClosure_v23_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_DLT_SFmuSFel_massClosestBJet_AMBISOLVING_25ns_v23/LeptonEfficiencies.root"
          file_SUSY_TP_path = basepath + "PLOTS_25ns_v23/MMClosure_v23_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_DLT_SFmuSFel_deltaRClosestBJet_AMBISOLVING_25ns_v23/LeptonEfficiencies.root"


  print("------------------------------------------------------------------------------------------")
  print("New Tag & Probe efficiency - Opening file:\n{0}".format(file_SUSY_TP_path))
  print("------------------------------------------------------------------------------------------")

  file_SUSY_TP = TFile(file_SUSY_TP_path)

  if args.closure:
      hist_SUSY_TP = ( file_SUSY_TP.Get(hist_eff_type+"_El_Pt_Efficiency_expectedbkg"), file_SUSY_TP.Get(hist_eff_type+"_Mu_Pt_Efficiency_expectedbkg") )[bool(args.flavour == "mu")]
  else:
      hist_SUSY_TP = ( file_SUSY_TP.Get(hist_eff_type+"_El_Pt_Efficiency_observed_sub"), file_SUSY_TP.Get(hist_eff_type+"_Mu_Pt_Efficiency_observed_sub") )[bool(args.flavour == "mu")]

  print("Reading histogram {0} from file {1}".format(hist_SUSY_TP.GetName(), file_SUSY_TP_path))

  hist_SUSY_TP.SetLineColor(kMagenta)
  hist_SUSY_TP.SetMarkerColor(kMagenta)

  hist_SUSY_TP.GetXaxis().SetTitle("p_{T}^{"+ lepton +"} [GeV]")
  hist_SUSY_TP.GetYaxis().SetTitle("#varepsilon")

  if eff_type == "real":
    hist_SUSY_TP.GetYaxis().SetRangeUser(0.5,1.0)
  elif eff_type == "fake":
    if args.flavour == "el":
      hist_SUSY_TP.GetYaxis().SetRangeUser(0.0,0.7)
    if args.flavour == "mu":
      hist_SUSY_TP.GetYaxis().SetRangeUser(0.0,0.5)

  hist_SUSY_TP.Draw("E0")
  legend.AddEntry(hist_SUSY_TP, "Tag & Probe", "P")

  # ----------------------------
  # TRUTH Tag & Probe efficiency
  # ----------------------------

  file_TRUTH_TP_path = basepath + "blahblah"

  if args.closure:
      if "SLT" in args.trigger:
          #file_TRUTH_TP_path = basepath + "PLOTS_25ns_v23/MMClosure_v23_TruthTP/OutputPlots_MMClosureRates_TruthTagProbe_NoCorr_SLT_SFmuSFel_25ns_v23/LeptonEfficiencies.root"
          #file_TRUTH_TP_path = basepath + "PLOTS_25ns_v23/MMClosure_v23_TruthTP/OutputPlots_MMClosureRates_TruthTagProbe_NoCorr_SLT_OFmuOFel_25ns_v23/LeptonEfficiencies.root"
          #file_TRUTH_TP_path = basepath + "PLOTS_25ns_v23/MMClosure_v23_TruthTP/OutputPlots_MMClosureRates_TruthTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_25ns_v23/LeptonEfficiencies.root"
          #file_TRUTH_TP_path = basepath + "PLOTS_25ns_v23/MMClosure_v23_TruthTP/OutputPlots_MMClosureRates_TruthTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_25ns_v23_OLD_EL_BINNING/LeptonEfficiencies.root"
          #
          file_TRUTH_TP_path = basepath + "PLOTS_25ns_v24/MMClosure_v24_TruthTP/OutputPlots_MMClosureRates_TruthTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_25ns_v24/LeptonEfficiencies.root"
      elif "DLT" in args.trigger:
          file_TRUTH_TP_path = basepath + "PLOTS_25ns_v23/MMClosure_v23_TruthTP/OutputPlots_MMClosureRates_TruthTagProbe_NoCorr_DLT_SFmuSFel_25ns_v23/LeptonEfficiencies.root"

  print("------------------------------------------------------------------------------------------")
  print("TRUTH Tag & Probe efficiency - Opening file:\n{0}".format(file_TRUTH_TP_path))
  print("------------------------------------------------------------------------------------------")

  file_TRUTH_TP = TFile(file_TRUTH_TP_path)

  if args.closure:
      hist_TRUTH_TP = ( file_TRUTH_TP.Get(hist_eff_type+"_El_Pt_Efficiency_expectedbkg"), file_TRUTH_TP.Get(hist_eff_type+"_Mu_Pt_Efficiency_expectedbkg") )[bool(args.flavour == "mu")]
  else:
      hist_TRUTH_TP = ( file_TRUTH_TP.Get(hist_eff_type+"_El_Pt_Efficiency_observed_sub"), file_TRUTH_TP.Get(hist_eff_type+"_Mu_Pt_Efficiency_observed_sub") )[bool(args.flavour == "mu")]

  print("Reading histogram {0} from file {1}".format(hist_TRUTH_TP.GetName(), file_TRUTH_TP_path))

  hist_TRUTH_TP.SetLineWidth(2)
  hist_TRUTH_TP.SetLineColor(kBlack)
  hist_TRUTH_TP.SetLineStyle(2)
  hist_TRUTH_TP.SetMarkerSize(1.0)
  hist_TRUTH_TP.SetMarkerColor(kBlack)
  hist_TRUTH_TP.SetMarkerStyle(24)

  hist_TRUTH_TP.Draw("E0 SAME")
  legend.AddEntry(hist_TRUTH_TP, "Truth", "P")

  # ---------------------
  # Likelihood efficiency
  # ---------------------

  #"""
  LH_init_path = basepath + "blahblah"

  if args.closure:
      if "SLT" in args.trigger:
          #LH_init_path = basepath + "PLOTS_25ns_v23/MMClosure_v23_LikelihoodFit/OutputPlots_MMClosureRates_LHFit_NoCorr_SFmuSFel_SLT_25ns_v21/LeptonEfficiencies_LH/"
          #LH_init_path = basepath + "PLOTS_25ns_v24/MMClosure_v24_LikelihoodFit/OutputPlots_MMClosureRates_LHFit_NoCorr_INCLUSIVE_FLAV_SLT_25ns_v24/LeptonEfficiencies_LH/"
          LH_init_path = basepath + "PLOTS_25ns_v24/MMClosure_v24_LikelihoodFit/OutputPlots_MMClosureRates_LHFit_NoCorr_INCLUSIVE_FLAV_DLT_25ns_v24/LeptonEfficiencies_LH/"
      elif "DLT" in args.trigger:
          LH_init_path = basepath + "PLOTS_25ns_v23/MMClosure_v23_LikelihoodFit/OutputPlots_MMClosureRates_LHFit_NoCorr_SFmuSFel_DLT_25ns_v21/LeptonEfficiencies_LH/"

  hist_LH_list = []

  if eff_type == "real":
    flav_comp_list = args.flavRealCR
  elif eff_type == "fake":
    flav_comp_list = args.flavFakeCR

  if ( args.flavour == "mu" and "elel" in flav_comp_list ): flav_comp_list.remove("elel")
  if ( args.flavour == "el" and "mumu" in flav_comp_list ): flav_comp_list.remove("mumu")

  for idx, flavcomp in enumerate(flav_comp_list,start=0):

    file_LH_path = LH_init_path + "LH_" + flavcomp + "/LH_efficiencies_" + eff_type + "_" + args.flavour + "_" + flavcomp + ".root"

    print("\t------------------------------------------------------------------------------------------")
    print("\tLikelihood efficiency - Opening file:\n{0}".format(file_LH_path))
    print("\t------------------------------------------------------------------------------------------")

    file_LH = TFile(file_LH_path)

    if eff_type == "real":
         hist_LH_name = "r_hist"
    elif eff_type == "fake":
         hist_LH_name = "f_hist"

    hist_LH = file_LH.Get(hist_LH_name)

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

  #gPad.SetLogx()
  #hist_TP.GetXaxis().SetMoreLogLabels()

  for idx, histpair in enumerate(hist_LH_list,start=0):
    #histpair[1].GetXaxis().SetMoreLogLabels()
    #if idx == 0:
    #  histpair[1].Draw("E0")
    #else:
    #  histpair[1].Draw("E0,SAME")

    # TEMP: plot LH efficiency only for fake muon
    if not ( eff_type == "fake" and args.flavour == "mu"): continue

    histpair[1].Draw("E0,SAME")

    if histpair[0] == "mumu": flag = "#mu#mu"
    else: flag = histpair[0]

    legend.AddEntry(histpair[1], "Likelihood - " +  flag, "P")

  #"""

  legend.Draw()

  leg_ATLAS.DrawLatex(0.6,0.35,"#bf{#it{ATLAS}} Work In Progress")
  leg_lumi.DrawLatex(0.6,0.27,"#sqrt{{s}} = 13 TeV, #int L dt = {0:.1f} fb^{{-1}}".format(args.lumi))

  # Add a vertical line to highlight relevant bins

  if eff_type == "real":
    low_y_edge = 0.5
    up_y_edge  = 1.0
  elif eff_type == "fake":
    low_y_edge = 0.0
    if args.flavour == "el":
      up_y_edge  = 0.7
    elif args.flavour == "mu":
      up_y_edge  = 0.5

  if "SLT" in args.trigger:

      if args.flavour == "el":

        refl_vert0 = TLine(26.0,0,26.0,up_y_edge)
        refl_vert0.SetLineStyle(2)
        refl_vert0.SetLineWidth(2)
        refl_vert0.Draw("SAME")

        refl_vert1 = TLine(60.0,low_y_edge,60.0,up_y_edge)
        refl_vert1.SetLineStyle(2)
        refl_vert1.SetLineWidth(2)
        refl_vert1.Draw("SAME")

        refl_vert2 = TLine(140.0,low_y_edge,140.0,up_y_edge)
        refl_vert2.SetLineStyle(2)
        refl_vert2.SetLineWidth(2)
        refl_vert2.Draw("SAME")

      elif args.flavour == "mu":

        refl_vert0 = TLine(26.0,low_y_edge,26.0,up_y_edge)
        refl_vert0.SetLineStyle(2)
        refl_vert0.SetLineWidth(2)
        refl_vert0.Draw("SAME")

        refl_vert1 = TLine(50.0,low_y_edge,50.0,up_y_edge)
        refl_vert1.SetLineStyle(2)
        refl_vert1.SetLineWidth(2)
        refl_vert1.Draw("SAME")

  elif "DLT" in args.trigger:

      if args.flavour == "el":

        refl_vert0 = TLine(17.0,low_y_edge,17.0,up_y_edge)
        refl_vert0.SetLineStyle(2)
        refl_vert0.SetLineWidth(2)
        refl_vert0.Draw("SAME")

      elif args.flavour == "mu":

        refl_vert0 = TLine(22.0,low_y_edge,22.0,up_y_edge)
        refl_vert0.SetLineStyle(2)
        refl_vert0.SetLineWidth(2)
        refl_vert0.Draw("SAME")

  canvasname = ( eff_type + "_eff_el_TP_LH", eff_type + "_eff_mu_TP_LH" )[bool(args.flavour == "mu")]
  c.SaveAs( canvasname + "_" + args.trigger[0] + ".png" )


def plotter_triggereff( eff_type ):

    if "SLT" in args.trigger:
        basepath = "$HOME/PhD/ttH_MultiLeptons/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_TTbar_RealFakeLep_TriggerEff_SLT/EfficiencyPlots_TriggerEff/BasicPlots/"
    elif "DLT" in args.trigger:
        basepath = "$HOME/PhD/ttH_MultiLeptons/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_TTbar_RealFakeLep_TriggerEff_DLT/EfficiencyPlots_TriggerEff/BasicPlots/"

    filename_L     = "RealFake_L_TriggerEfficiency.root"
    filename_T     = "RealFake_T_TriggerEfficiency.root"
    filename_AntiT = "RealFake_AntiT_TriggerEfficiency.root"

    file_L     = TFile(basepath + filename_L)
    file_T     = TFile(basepath + filename_T)
    file_AntiT = TFile(basepath + filename_AntiT)

    if   args.flavour == "mu": flavour = "Mu"
    elif args.flavour == "el": flavour = "El"

    if   eff_type == "real": efficiency = "Real"
    elif eff_type == "fake": efficiency = "Fake"

    variable = "Pt"

    append = ("_observed_sub","_expectedbkg")[bool(args.closure)]

    histname_L     = efficiency + "_" + flavour + "_" + variable + "_L_TriggerEfficiency" + append
    histname_T     = efficiency + "_" + flavour + "_" + variable + "_T_TriggerEfficiency" + append
    histname_AntiT = efficiency + "_" + flavour + "_" + variable + "_AntiT_TriggerEfficiency" + append

    print("Try to get histogram:\n{0}\nfrom file:\n{1}\n".format(histname_L,basepath + filename_L))
    print("Try to get histogram:\n{0}\nfrom file:\n{1}\n".format(histname_T,basepath + filename_T))
    print("Try to get histogram:\n{0}\nfrom file:\n{1}\n".format(histname_AntiT,basepath + filename_AntiT))

    file_L.cd()
    hist_L = file_L.Get(histname_L)
    hist_L.SetDirectory(0)

    file_T.cd()
    hist_T = file_T.Get(histname_T)
    hist_T.SetDirectory(0)

    file_AntiT.cd()
    hist_AntiT = file_AntiT.Get(histname_AntiT)
    hist_AntiT.SetDirectory(0)

    hist_L.SetLineStyle(1)
    hist_L.SetLineColor(kBlack)

    hist_T.SetLineStyle(2)
    hist_T.SetLineColor(kBlack)
    hist_T.SetMarkerStyle(24)

    delta_eff = hist_L.Clone("DeltaEff")
    delta_eff.SetXTitle(hist_L.GetXaxis().GetTitle())
    delta_eff.SetYTitle("#Delta#varepsilon/#varepsilon [%]")
    delta_eff.GetXaxis().SetTitleSize(0.15)
    delta_eff.GetYaxis().SetTitleSize(0.15)
    delta_eff.GetXaxis().SetTitleOffset(0.90)
    delta_eff.GetYaxis().SetTitleOffset(0.35)
    delta_eff.GetXaxis().SetLabelSize(0.15)
    delta_eff.GetYaxis().SetLabelSize(0.12)
    delta_eff.GetYaxis().SetRangeUser(-50.0,20.0)
    delta_eff.GetYaxis().SetNdivisions(505) #(5)
    delta_eff.SetLineColor(kRed)
    delta_eff.SetMarkerColor(kRed)
    delta_eff.SetMarkerSize(1)

    delta_eff.Add(hist_T, -1)
    delta_eff.Divide(hist_T)
    delta_eff.Scale(100.0)

    legend = TLegend(0.6,0.4,0.8,0.55) # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
    legend.SetHeader(flavour + " - " + efficiency)
    legend.SetBorderSize(0)	# no border
    legend.SetFillStyle(0)	# Legend transparent background
    legend.SetTextSize(0.04)	# Increase entry font size!
    legend.SetTextFont(42)	# Helvetica

    legend.AddEntry(hist_T, "Trigger efficiency - T", "P")
    legend.AddEntry(hist_L, "Trigger efficiency - L", "P")

    leg_ATLAS = TLatex()
    leg_lumi  = TLatex()
    leg_ATLAS.SetTextSize(0.04)
    leg_ATLAS.SetNDC()
    leg_lumi.SetTextSize(0.04)
    leg_lumi.SetNDC()

    # ---------------------------------------------------------------

    c = TCanvas("c1","Temp",50,50,800,600)

    pad1 = TPad("pad1", "", 0, 0.25, 1, 1)
    pad2 = TPad("pad2", "", 0, 0,   1, 0.25)
    pad1.SetBottomMargin(0.02)
    pad2.SetBottomMargin(0.4)
    pad1.Draw()
    pad2.Draw()

    pad1.cd()

    hist_L.GetXaxis().SetLabelSize(0)
    hist_L.GetXaxis().SetLabelOffset(999)

    hist_L.Draw("E0")
    hist_T.Draw("E0 SAME")

    if "SLT" in args.trigger:

        if args.flavour == "el":

	    refl_vert0 = TLine(26.0,hist_L.GetMinimum(),26.0,hist_L.GetMaximum())
            refl_vert0.SetLineStyle(2)
            refl_vert0.SetLineWidth(2)
            refl_vert0.Draw("SAME")

            refl_vert1 = TLine(60.0,hist_L.GetMinimum(),60.0,hist_L.GetMaximum())
            refl_vert1.SetLineStyle(2)
            refl_vert1.SetLineWidth(2)
            refl_vert1.Draw("SAME")

            refl_vert2 = TLine(140.0,hist_L.GetMinimum(),140.0,hist_L.GetMaximum())
            refl_vert2.SetLineStyle(2)
            refl_vert2.SetLineWidth(2)
            refl_vert2.Draw("SAME")

	elif args.flavour == "mu":

	    refl_vert0 = TLine(26.0,hist_L.GetMinimum(),26.0,hist_L.GetMaximum())
            refl_vert0.SetLineStyle(2)
            refl_vert0.SetLineWidth(2)
            refl_vert0.Draw("SAME")

	    refl_vert1 = TLine(50.0,hist_L.GetMinimum(),50.0,hist_L.GetMaximum())
            refl_vert1.SetLineStyle(2)
            refl_vert1.SetLineWidth(2)
            refl_vert1.Draw("SAME")

    elif "DLT" in args.trigger:

	if args.flavour == "el":

	    refl_vert0 = TLine(17.0,hist_L.GetMinimum(),17.0,hist_L.GetMaximum())
            refl_vert0.SetLineStyle(2)
            refl_vert0.SetLineWidth(2)
            refl_vert0.Draw("SAME")

	elif args.flavour == "mu":

	    refl_vert0 = TLine(22.0,hist_L.GetMinimum(),22.0,hist_L.GetMaximum())
            refl_vert0.SetLineStyle(2)
            refl_vert0.SetLineWidth(2)
            refl_vert0.Draw("SAME")

    legend.Draw()
    leg_ATLAS.DrawLatex(0.6,0.27,"#bf{#it{ATLAS}} Work In Progress")
    leg_lumi.DrawLatex(0.6,0.2,"#sqrt{{s}} = 13 TeV, #int L dt = {0:.1f} fb^{{-1}}".format(args.lumi))

    pad2.cd()

    delta_eff.Draw("E0")

    refl = TLine(delta_eff.GetBinLowEdge(1), 0.0, delta_eff.GetBinLowEdge(delta_eff.GetNbinsX()+1), 0.0)
    refl.SetLineStyle(2)
    refl.SetLineWidth(2)
    refl.Draw("SAME")

    outpath = basepath
    if outpath[-1] == '/':
      outpath = outpath[:-1]

    c.SaveAs( outpath + "/TriggerEfficiencyRatio_" + efficiency + "_" + flavour + "_" + variable + ".png" )

# ----------------

if __name__ == "__main__":

    if args.triggerEff:
        plotter_triggereff("real")
        plotter_triggereff("fake")
    else:
        plotter_flavours("real")
        plotter_flavours("fake")

    #plot2Dhist()
