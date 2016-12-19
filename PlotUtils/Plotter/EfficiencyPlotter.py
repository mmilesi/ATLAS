#!/usr/bin/env python

""" EfficiencyPlotter.py: plot efficiencies"""

__author__     = "Marco Milesi"
__email__      = "marco.milesi@cern.ch"
__maintainer__ = "Marco Milesi"

import os, sys, array, math

sys.path.append(os.path.abspath(os.path.curdir))

from ROOT import ROOT, gROOT, Double, TPad, TLine, TH1, TH1D, TFile, TCanvas, TLegend, TLatex, TGraphAsymmErrors, TEfficiency
from ROOT import kBlue, kOrange, kPink, kGreen, kRed, kYellow, kTeal, kMagenta, kViolet, kAzure, kCyan, kSpring, kGray, kBlack, kWhite
from ROOT import kFullCircle, kCircle, kOpenTriangleUp, kDot

gROOT.Reset()
gROOT.LoadMacro("$HOME/RootUtils/AtlasStyle.C")
from ROOT import SetAtlasStyle
SetAtlasStyle()

TH1.SetDefaultSumw2()

gROOT.SetBatch(True)

class Plot:

    luminosity = 100

    legend = TLegend(0.45,0.5,0.925,0.8) # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
    legend.SetBorderSize(0)  # no border
    legend.SetFillStyle(0) # Legend transparent background
    legend.SetTextSize(0.03) # Increase entry font size!
    legend.SetTextFont(42)   # Helvetica

    legendATLAS = TLatex()
    legendLumi  = TLatex()
    legendATLAS.SetTextSize(0.03)
    legendATLAS.SetNDC()
    legendLumi.SetTextSize(0.03)
    legendLumi.SetNDC()

    reflines = []

    def __init__( self, sourceName, sourcePath, properties={} ):

        f = TFile(sourcePath)

        if not f:
	   sys.exit("ERROR: file\n{0}\ncannot be found!".format(sourcePath))

        self.__hist = f.Get(sourceName)

	if not self.__hist:
	   sys.exit("ERROR: histogram:\n{0}\ncannot be found in file:\n{1}".format(sourceName,sourcePath))

	self.__hist.SetDirectory(0)

	self.__props  = properties


    def setPropery( self, propID, propValue ):

        self.__props[propID] = propValue

    def makePlot( self, canvas, drawOpt ):

        canvas.cd()

	if self.__props.get("xAxisTitle") : self.__hist.GetXaxis().SetTitle( self.__props["xAxisTitle"] )
	if self.__props.get("yAxisTitle") : self.__hist.GetYaxis().SetTitle( self.__props["yAxisTitle"] )

	if self.__props.get("xAxisRange") : self.__hist.GetXaxis().SetRangeUser( self.__props["xAxisRange"][0], self.__props["xAxisRange"][1] )
	if self.__props.get("yAxisRange") : self.__hist.GetYaxis().SetRangeUser( self.__props["yAxisRange"][0], self.__props["yAxisRange"][1] )

	if self.__props.get("colour") :
	    self.__hist.SetLineColor(self.__props["colour"])
	    self.__hist.SetMarkerColor(self.__props["colour"])

	if self.__props.get("lineStyle")   : self.__hist.SetLineStyle(self.__props["lineStyle"])
	if self.__props.get("lineWidth")   : self.__hist.SetLineWidth(self.__props["lineWidth"])

	if self.__props.get("markerStyle") : self.__hist.SetMarkerStyle(self.__props["markerStyle"])
	if self.__props.get("markerSize")  : self.__hist.SetMarkerSize(self.__props["markerSize"])

        if self.__props.get("legend")      : Plot.legend.AddEntry(self.__hist, self.__props["legend"], "P")

	self.__hist.Draw( drawOpt )


class MultiPlot:

    def __init__( self, plots=[] ):

	self.__plotlist  = plots

    def makeMultiPlot( self, savePath, saveName ):

        tokens = saveName.split('_')

        c = TCanvas("c_"+tokens[0]+"_"+tokens[1],"Efficiency",50,50,1300,800)

        for idx, plot in enumerate(self.__plotlist):

	    if idx == 0:
                plot.makePlot( c, "E0")
	    else:
                plot.makePlot( c, "E0 SAME")

        for refl in Plot.reflines:
            refl.SetLineStyle(2)
	    refl.Draw("SAME")

        Plot.legend.Draw()

	Plot.legendATLAS.DrawLatex(0.6,0.35,"#bf{#it{ATLAS}} Work In Progress")
        Plot.legendLumi.DrawLatex(0.6,0.27,"#sqrt{{s}} = 13 TeV, #int L dt = {0:.1f} fb^{{-1}}".format(Plot.luminosity))

        for ext in ["png","eps"]:
	    c.SaveAs( savePath + "/" + saveName + "." + ext )

# -------------------------------------------------------------------------

def plotFakeElectron():

    plotlist = []

    Plot.legend.SetHeader("#varepsilon_{fake} - Electrons")

    p0_props = {
                "legend"      : "Truth",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.0,0.9),
                "colour"      : kBlack,
                "lineStyle"   : 2,
                "lineWidth"   : 2,
                "markerStyle" : 24,
                "markerSize"  : 1,
               }

    p0 = Plot("Fake_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_TruthTP/OutputPlots_MMClosureRates_TruthTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_25ns_v24/LeptonEfficiencies.root", p0_props)

    plotlist.append(p0)

    p1_props = {
                "legend"      : "Tag & Probe (OF, #mu tag)",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.0,0.7),
                "colour"      : kOrange+10,
                "lineStyle"   : 1,
                "lineWidth"   : 2,
                "markerStyle" : 20,
                "markerSize"  : 1,
               }

    p1 = Plot("Fake_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIsoTagPt30_ForceProbeToBeFake/LeptonEfficiencies.root", p1_props)

    plotlist.append(p1)

    p2_props = {
                "legend"      : "Tag & Probe (OF, #mu tag) - probe e TM",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.0,0.7),
                "colour"      : kOrange+10,
                "lineStyle"   : 3,
                "lineWidth"   : 2,
                "markerStyle" : 24,
                "markerSize"  : 1,
               }

    p2 = Plot("Fake_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIsoTagPt30_ForceProbeToBeFake_TRIGMATCH_EFF/LeptonEfficiencies.root", p2_props)

    plotlist.append(p2)

    p3_props = {
                "legend"      : "Tag & Probe (OF, #mu tag)- probe e NOT TM",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.0,0.7),
                "colour"      : kOrange+10,
                "lineStyle"   : 3,
                "lineWidth"   : 2,
                "markerStyle" : 30,
                "markerSize"  : 3,
               }

    p3 = Plot("Fake_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIsoTagPt30_ForceProbeToBeFake_NOT_TRIGMATCH_EFF/LeptonEfficiencies.root", p3_props)

    plotlist.append(p3)

    # Add vertical reference lines for trigger thresholds.
    # Use the first histogram in the list to set the range

    refl_0 = TLine(26.0,p0_props["yAxisRange"][0],26.0,p0_props["yAxisRange"][1])
    Plot.reflines.append(refl_0)
    refl_1 = TLine(60.0,p0_props["yAxisRange"][0],60.0,p0_props["yAxisRange"][1])
    Plot.reflines.append(refl_1)
    refl_2 = TLine(140.0,p0_props["yAxisRange"][0],140.0,p0_props["yAxisRange"][1])
    Plot.reflines.append(refl_2)

    multiP = MultiPlot( plots=plotlist )

    multiP.makeMultiPlot( "./PLOTS_25ns_v24/CombinedEfficiencies_Plots", "Fake_El_Compare" )

    # Clear before returning

    del plotlist[:]

    Plot.legend.Clear()

    del Plot.reflines[:]

def plotFakeElectron_anyProbe():

    plotlist = []

    Plot.legend.SetHeader("#varepsilon_{fake} - Electrons")

    p0_props = {
                "legend"      : "Truth",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.0,0.9),
                "colour"      : kBlack,
                "lineStyle"   : 2,
                "lineWidth"   : 2,
                "markerStyle" : 24,
                "markerSize"  : 1,
               }

    p0 = Plot("Fake_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_TruthTP/OutputPlots_MMClosureRates_TruthTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_25ns_v24/LeptonEfficiencies.root", p0_props)

    plotlist.append(p0)

    p1_props = {
                "legend"      : "Tag & Probe (OF, #mu tag) - probe e fake",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.0,0.7),
                "colour"      : kOrange+10,
                "lineStyle"   : 1,
                "lineWidth"   : 2,
                "markerStyle" : 20,
                "markerSize"  : 1,
               }

    p1 = Plot("Fake_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIsoTagPt30_ForceProbeToBeFake/LeptonEfficiencies.root", p1_props)

    plotlist.append(p1)

    p2_props = {
                "legend"      : "Tag & Probe (OF, #mu tag)",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.0,0.7),
                "colour"      : kMagenta,
                "lineStyle"   : 1,
                "lineWidth"   : 2,
                "markerStyle" : 29,
                "markerSize"  : 2,
               }

    p2 = Plot("Fake_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIsoTagPt30/LeptonEfficiencies.root", p2_props)

    plotlist.append(p2)

    # Add vertical reference lines for trigger thresholds.
    # Use the first histogram in the list to set the range

    refl_0 = TLine(26.0,p0_props["yAxisRange"][0],26.0,p0_props["yAxisRange"][1])
    Plot.reflines.append(refl_0)
    refl_1 = TLine(60.0,p0_props["yAxisRange"][0],60.0,p0_props["yAxisRange"][1])
    Plot.reflines.append(refl_1)
    refl_2 = TLine(140.0,p0_props["yAxisRange"][0],140.0,p0_props["yAxisRange"][1])
    Plot.reflines.append(refl_2)

    multiP = MultiPlot( plots=plotlist )

    multiP.makeMultiPlot( "./PLOTS_25ns_v24/CombinedEfficiencies_Plots", "Fake_El_Compare_anyProbe" )

    # Clear before returning

    del plotlist[:]

    Plot.legend.Clear()

    del Plot.reflines[:]


def plotFakeMuon():

    plotlist = []

    Plot.legend.SetHeader("#varepsilon_{fake} - Muons")

    p0_props = {
                "legend"      : "Truth",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.0,0.3),
                "colour"      : kBlack,
                "lineStyle"   : 2,
                "lineWidth"   : 2,
                "markerStyle" : 24,
                "markerSize"  : 1,
               }

    p0 = Plot("Fake_Mu_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_TruthTP/OutputPlots_MMClosureRates_TruthTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_25ns_v24/LeptonEfficiencies.root", p0_props)

    plotlist.append(p0)

    p1_props = {
                "legend"      : "Likelihood (#mu#mu)",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.0,0.3),
                "colour"      : kOrange+10,
                "lineStyle"   : 1,
                "lineWidth"   : 2,
                "markerStyle" : 20,
                "markerSize"  : 1,
               }

    p1 = Plot("f_hist","./PLOTS_25ns_v24/MMClosure_v24_LikelihoodFit/OutputPlots_MMClosureRates_LHFit_NoCorr_INCLUSIVE_FLAV_DLT_25ns_v24/LeptonEfficiencies_LH/LH_mumu/LH_efficiencies_fake_mu_mumu.root", p1_props)

    plotlist.append(p1)

    # Add vertical reference lines for trigger thresholds.
    # Use the first histogram in the list to set the range

    refl_0 = TLine(26.0,p0_props["yAxisRange"][0],26.0,p0_props["yAxisRange"][1])
    Plot.reflines.append(refl_0)
    refl_1 = TLine(50.0,p0_props["yAxisRange"][0],50.0,p0_props["yAxisRange"][1])
    Plot.reflines.append(refl_1)

    multiP = MultiPlot( plots=plotlist )

    multiP.makeMultiPlot( "./PLOTS_25ns_v24/CombinedEfficiencies_Plots", "Fake_Mu_Compare" )

    # Clear before returning

    del plotlist[:]

    Plot.legend.Clear()

    del Plot.reflines[:]


def plotRealElectron():

    plotlist = []

    Plot.legend.SetHeader("#varepsilon_{real} - Electrons")

    p0_props = {
                "legend"      : "Truth",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.3,1.0),
                "colour"      : kBlack,
                "lineStyle"   : 2,
                "lineWidth"   : 2,
                "markerStyle" : 24,
                "markerSize"  : 1,
               }

    p0 = Plot("Real_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_TruthTP/OutputPlots_MMClosureRates_TruthTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_25ns_v24/LeptonEfficiencies.root", p0_props)

    plotlist.append(p0)

    p1_props = {
                "legend"      : "Tag & Probe (OF)",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.5,1.0),
                "colour"      : kBlue,
                "lineStyle"   : 1,
                "lineWidth"   : 2,
                "markerStyle" : 20,
                "markerSize"  : 1,
               }

    p1 = Plot("Real_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIsoTagPt30_ForceProbeToBeFake/LeptonEfficiencies.root", p1_props)

    plotlist.append(p1)

    p2_props = {
                "legend"      : "Tag & Probe (OF) - probe TM",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.0,0.7),
                "colour"      : kBlue,
                "lineStyle"   : 3,
                "lineWidth"   : 2,
                "markerStyle" : 24,
                "markerSize"  : 1,
               }

    p2 = Plot("Real_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIsoTagPt30_ForceProbeToBeFake_TRIGMATCH_EFF/LeptonEfficiencies.root", p2_props)

    plotlist.append(p2)

    p3_props = {
                "legend"      : "Tag & Probe (OF) - probe NOT TM",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.0,0.7),
                "colour"      : kBlue,
                "lineStyle"   : 3,
                "lineWidth"   : 2,
                "markerStyle" : 30,
                "markerSize"  : 3,
               }

    p3 = Plot("Real_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIsoTagPt30_ForceProbeToBeFake_NOT_TRIGMATCH_EFF/LeptonEfficiencies.root", p3_props)

    plotlist.append(p3)

    # Add vertical reference lines for trigger thresholds.
    # Use the first histogram in the list to set the range

    refl_0 = TLine(26.0,p0_props["yAxisRange"][0],26.0,p0_props["yAxisRange"][1])
    Plot.reflines.append(refl_0)
    refl_1 = TLine(60.0,p0_props["yAxisRange"][0],60.0,p0_props["yAxisRange"][1])
    Plot.reflines.append(refl_1)
    refl_2 = TLine(140.0,p0_props["yAxisRange"][0],140.0,p0_props["yAxisRange"][1])
    Plot.reflines.append(refl_2)

    multiP = MultiPlot( plots=plotlist )

    multiP.makeMultiPlot( "./PLOTS_25ns_v24/CombinedEfficiencies_Plots", "Real_El_Compare" )

    # Clear before returning

    del plotlist[:]

    Plot.legend.Clear()

    del Plot.reflines[:]


def plotRealMuon():

    plotlist = []

    Plot.legend.SetHeader("#varepsilon_{real} - Muons")

    p0_props = {
                "legend"      : "Truth",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.5,1.0),
                "colour"      : kBlack,
                "lineStyle"   : 2,
                "lineWidth"   : 2,
                "markerStyle" : 24,
                "markerSize"  : 1,
               }

    p0 = Plot("Real_Mu_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_TruthTP/OutputPlots_MMClosureRates_TruthTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_25ns_v24/LeptonEfficiencies.root", p0_props)

    plotlist.append(p0)

    p1_props = {
                "legend"      : "Tag & Probe (OF)",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.5,1.0),
                "colour"      : kBlue,
                "lineStyle"   : 1,
                "lineWidth"   : 2,
                "markerStyle" : 20,
                "markerSize"  : 1,
               }

    p1 = Plot("Real_Mu_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIsoTagPt30_ForceProbeToBeFake/LeptonEfficiencies.root", p1_props)

    plotlist.append(p1)

    p2_props = {
                "legend"      : "Tag & Probe (OF) - probe TM",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.0,0.7),
                "colour"      : kBlue,
                "lineStyle"   : 3,
                "lineWidth"   : 2,
                "markerStyle" : 24,
                "markerSize"  : 1,
               }

    p2 = Plot("Real_Mu_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIsoTagPt30_ForceProbeToBeFake_TRIGMATCH_EFF/LeptonEfficiencies.root", p2_props)

    plotlist.append(p2)

    p3_props = {
                "legend"      : "Tag & Probe (OF) - probe NOT TM",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.0,0.7),
                "colour"      : kBlue,
                "lineStyle"   : 3,
                "lineWidth"   : 2,
                "markerStyle" : 30,
                "markerSize"  : 3,
               }

    p3 = Plot("Real_Mu_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIsoTagPt30_ForceProbeToBeFake_NOT_TRIGMATCH_EFF/LeptonEfficiencies.root", p3_props)

    plotlist.append(p3)

    # Add vertical reference lines for trigger thresholds.
    # Use the first histogram in the list to set the range

    refl_0 = TLine(26.0,p0_props["yAxisRange"][0],26.0,p0_props["yAxisRange"][1])
    Plot.reflines.append(refl_0)
    refl_1 = TLine(50.0,p0_props["yAxisRange"][0],50.0,p0_props["yAxisRange"][1])
    Plot.reflines.append(refl_1)

    multiP = MultiPlot( plots=plotlist )

    multiP.makeMultiPlot( "./PLOTS_25ns_v24/CombinedEfficiencies_Plots", "Real_Mu_Compare" )

    # Clear before returning

    del plotlist[:]

    Plot.legend.Clear()

    del Plot.reflines[:]

# -------------------------------------------------------------------------

if __name__ == "__main__":

    Plot.luminosity = 36.4702

    plotFakeElectron()
    plotFakeElectron_anyProbe()
    plotFakeMuon()
    plotRealElectron()
    plotRealMuon()
