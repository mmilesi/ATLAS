#!/usr/bin/env python

from ROOT import kBlue, kOrange, kPink, kGreen, kRed, kYellow, kTeal, kMagenta, kViolet, kAzure, kCyan, kSpring, kGray, kBlack, kWhite
from ROOT import TLine

from EfficiencyPlotterClasses import Plot, MultiPlot

Plot.luminosity = 36.4702

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

    p0 = Plot("Fake_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24_ElNoIso/MMClosure_v24_TruthTP/OutputPlots_MMClosureRates_TruthTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_25ns_v24/LeptonEfficiencies.root", p0_props)

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

    p1 = Plot("Fake_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24_ElNoIso/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_ForceProbeToBeFake/LeptonEfficiencies.root", p1_props)

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

    p2 = Plot("Fake_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24_ElNoIso/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_ForceProbeToBeFake_TRIGMATCH_EFF/LeptonEfficiencies.root", p2_props)

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

    p3 = Plot("Fake_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24_ElNoIso/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_ForceProbeToBeFake_NOT_TRIGMATCH_EFF/LeptonEfficiencies.root", p3_props)

    plotlist.append(p3)

    # Add vertical reference lines for trigger thresholds.
    # Use the first histogram in the list to set the range

    refl_0 = TLine(26.0,p0_props["yAxisRange"][0],26.0,p0_props["yAxisRange"][1])
    Plot.reflines.append(refl_0)
    refl_1 = TLine(60.0,p0_props["yAxisRange"][0],60.0,p0_props["yAxisRange"][1])
    Plot.reflines.append(refl_1)
    #refl_2 = TLine(140.0,p0_props["yAxisRange"][0],140.0,p0_props["yAxisRange"][1])
    #Plot.reflines.append(refl_2)

    multiP = MultiPlot( plots=plotlist )

    multiP.makeMultiPlot( "./PLOTS_25ns_v24_ElNoIso/CombinedEfficiencies_Plots", "Fake_El_Compare" )

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

    p0 = Plot("Fake_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24_ElNoIso/MMClosure_v24_TruthTP/OutputPlots_MMClosureRates_TruthTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_25ns_v24/LeptonEfficiencies.root", p0_props)

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

    p1 = Plot("Fake_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24_ElNoIso/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_ForceProbeToBeFake/LeptonEfficiencies.root", p1_props)

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

    p2 = Plot("Fake_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24_ElNoIso/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24/LeptonEfficiencies.root", p2_props)

    plotlist.append(p2)

    # Add vertical reference lines for trigger thresholds.
    # Use the first histogram in the list to set the range

    refl_0 = TLine(26.0,p0_props["yAxisRange"][0],26.0,p0_props["yAxisRange"][1])
    Plot.reflines.append(refl_0)
    refl_1 = TLine(60.0,p0_props["yAxisRange"][0],60.0,p0_props["yAxisRange"][1])
    Plot.reflines.append(refl_1)
    #refl_2 = TLine(140.0,p0_props["yAxisRange"][0],140.0,p0_props["yAxisRange"][1])
    #Plot.reflines.append(refl_2)

    multiP = MultiPlot( plots=plotlist )

    multiP.makeMultiPlot( "./PLOTS_25ns_v24_ElNoIso/CombinedEfficiencies_Plots", "Fake_El_Compare_anyProbe" )

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

    p0 = Plot("Real_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24_ElNoIso/MMClosure_v24_TruthTP/OutputPlots_MMClosureRates_TruthTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_25ns_v24/LeptonEfficiencies.root", p0_props)

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

    p1 = Plot("Real_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24_ElNoIso/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_ForceProbeToBeFake/LeptonEfficiencies.root", p1_props)

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

    p2 = Plot("Real_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24_ElNoIso/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_ForceProbeToBeFake_TRIGMATCH_EFF/LeptonEfficiencies.root", p2_props)

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

    p3 = Plot("Real_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24_ElNoIso/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_ForceProbeToBeFake_NOT_TRIGMATCH_EFF/LeptonEfficiencies.root", p3_props)

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

    multiP.makeMultiPlot( "./PLOTS_25ns_v24_ElNoIso/CombinedEfficiencies_Plots", "Real_El_Compare" )

    # Clear before returning

    del plotlist[:]

    Plot.legend.Clear()

    del Plot.reflines[:]


def plotFakeElectron_BaselineElIso():

    plotlist = []

    Plot.legend.SetHeader("#varepsilon_{fake} - Electrons")

    p0_props = {
                "legend"      : "Tag & Probe (OF, #mu tag) - w/o baseline e iso",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.0,0.7),
                "colour"      : kOrange+10,
                "lineStyle"   : 1,
                "lineWidth"   : 2,
                "markerStyle" : 20,
                "markerSize"  : 1,
               }

    p0 = Plot("Fake_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24_ElNoIso/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_ForceProbeToBeFake/LeptonEfficiencies.root", p0_props)

    plotlist.append(p0)

    p1_props = {
                "legend"      : "Tag & Probe (OF, #mu tag)",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.0,0.7),
                "colour"      : kMagenta,
                "lineStyle"   : 1,
                "lineWidth"   : 2,
                "markerStyle" : 29,
                "markerSize"  : 2,
               }

    p1 = Plot("Fake_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIsoTagd0sig15_ForceProbeToBeFake_REBINNED/LeptonEfficiencies.root", p1_props)

    plotlist.append(p1)

    # Add vertical reference lines for trigger thresholds.
    # Use the first histogram in the list to set the range

    refl_0 = TLine(26.0,p0_props["yAxisRange"][0],26.0,p0_props["yAxisRange"][1])
    Plot.reflines.append(refl_0)
    refl_1 = TLine(60.0,p0_props["yAxisRange"][0],60.0,p0_props["yAxisRange"][1])
    Plot.reflines.append(refl_1)
    #refl_2 = TLine(140.0,p0_props["yAxisRange"][0],140.0,p0_props["yAxisRange"][1])
    #Plot.reflines.append(refl_2)

    multiP = MultiPlot( plots=plotlist )

    multiP.makeMultiPlot( "./PLOTS_25ns_v24_ElNoIso/CombinedEfficiencies_Plots", "Fake_El_Compare_ElIso" )

    # Clear before returning

    del plotlist[:]

    Plot.legend.Clear()

    del Plot.reflines[:]


def plotRealElectron_BaselineElIso():

    plotlist = []

    Plot.legend.SetHeader("#varepsilon_{real} - Electrons")

    p0_props = {
                "legend"      : "Tag & Probe (OF) - w/o baseline e iso",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.5,1.0),
                "colour"      : kBlue,
                "lineStyle"   : 1,
                "lineWidth"   : 2,
                "markerStyle" : 20,
                "markerSize"  : 1,
               }

    p0 = Plot("Real_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24_ElNoIso/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_ForceProbeToBeFake/LeptonEfficiencies.root", p0_props)

    plotlist.append(p0)

    p1_props = {
                "legend"      : "Tag & Probe (OF)",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.5,1.0),
                "colour"      : kMagenta,
                "lineStyle"   : 1,
                "lineWidth"   : 2,
                "markerStyle" : 20,
                "markerSize"  : 1,
	       }

    p1 = Plot("Real_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIsoTagd0sig15_ForceProbeToBeFake_REBINNED/LeptonEfficiencies.root", p1_props)

    plotlist.append(p1)

    # Add vertical reference lines for trigger thresholds.
    # Use the first histogram in the list to set the range

    refl_0 = TLine(26.0,p0_props["yAxisRange"][0],26.0,p0_props["yAxisRange"][1])
    Plot.reflines.append(refl_0)
    refl_1 = TLine(60.0,p0_props["yAxisRange"][0],60.0,p0_props["yAxisRange"][1])
    Plot.reflines.append(refl_1)
    #refl_2 = TLine(140.0,p0_props["yAxisRange"][0],140.0,p0_props["yAxisRange"][1])
    #Plot.reflines.append(refl_2)

    multiP = MultiPlot( plots=plotlist )

    multiP.makeMultiPlot( "./PLOTS_25ns_v24_ElNoIso/CombinedEfficiencies_Plots", "Real_El_Compare_ElIso" )

    # Clear before returning

    del plotlist[:]

    Plot.legend.Clear()

    del Plot.reflines[:]
