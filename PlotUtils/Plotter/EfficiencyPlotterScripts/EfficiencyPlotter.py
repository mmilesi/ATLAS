#!/usr/bin/env python

import os

from ROOT import kBlue, kOrange, kPink, kGreen, kRed, kYellow, kTeal, kMagenta, kViolet, kAzure, kCyan, kSpring, kGray, kBlack, kWhite
from ROOT import TLine

from EfficiencyPlotterClasses import Plot, MultiPlot

luminosity = 36.1

def plotFakeElectron():

    plotlist = []

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

    p0 = Plot("Fake_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_TruthTP/OutputPlots_MMClosureRates_TruthTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_25ns_v24_REBINNED/LeptonEfficiencies.root", p0_props)

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

    p1 = Plot("Fake_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIsoTagd0sig15_ForceProbeToBeFake_REBINNED/LeptonEfficiencies.root", p1_props)

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

    p2 = Plot("Fake_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIsoTagd0sig15_ForceProbeToBeFake_TRIGMATCH_EFF_REBINNED/LeptonEfficiencies.root", p2_props)

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

    p3 = Plot("Fake_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIsoTagd0sig15_ForceProbeToBeFake_NOT_TRIGMATCH_EFF_REBINNED/LeptonEfficiencies.root", p3_props)

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
    multiP.luminosity = luminosity
    multiP.buildLegend( header="#varepsilon_{fake} - Electrons", lcoords=[0.5,0.7,0.93,0.91] )
    multiP.setCanvasCoords([50,50,800,800])
    # multiP.setPlotTitle("My title",(0.6,0.88))

    multiP.makeMultiPlot( "./PLOTS_25ns_v24/CombinedEfficiencies_Plots", "Fake_El_Compare" )

    # Clear before returning

    del plotlist[:]


def plotFakeMuon():

    plotlist = []

    Plot.legend.SetHeader("#varepsilon_{fake} - Muons")

    p0_props = {
                "legend"      : "Truth",
                "yAxisTitle"  : "#varepsilon",
                #"yAxisRange"  : (0.0,0.3),
                "yAxisRange"  : (0.0,0.6),
                "colour"      : kBlack,
                "lineStyle"   : 2,
                "lineWidth"   : 2,
                "markerStyle" : 24,
                "markerSize"  : 1,
               }

    p0 = Plot("Fake_Mu_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_TruthTP/OutputPlots_MMClosureRates_TruthTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_25ns_v24_REBINNED/LeptonEfficiencies.root", p0_props)

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

    p2_props = {
                "legend"      : "Likelihood (#mu#mu) - (HLT_mu26_ivarmedium||HLT_mu50) TM",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.0,0.7),
                "colour"      : kOrange+10,
                "lineStyle"   : 3,
                "lineWidth"   : 2,
                "markerStyle" : 24,
                "markerSize"  : 1,
               }

    p2 = Plot("f_hist","./PLOTS_25ns_v24/MMClosure_v24_LikelihoodFit/OutputPlots_MMClosureRates_LHFit_NoCorr_INCLUSIVE_FLAV_DLT_25ns_v24_TRIGMATCH_EFF/LeptonEfficiencies_LH/LH_mumu/LH_efficiencies_fake_mu_mumu.root", p2_props)

    plotlist.append(p2)

    p3_props = {
                "legend"      : "Likelihood (#mu#mu) - (HLT_mu26_ivarmedium||HLT_mu50) NOT TM",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.0,0.7),
                "colour"      : kOrange+10,
                "lineStyle"   : 3,
                "lineWidth"   : 2,
                "markerStyle" : 30,
                "markerSize"  : 3,
               }

    p3 = Plot("f_hist","./PLOTS_25ns_v24/MMClosure_v24_LikelihoodFit/OutputPlots_MMClosureRates_LHFit_NoCorr_INCLUSIVE_FLAV_DLT_25ns_v24_NOT_TRIGMATCH_EFF/LeptonEfficiencies_LH/LH_mumu/LH_efficiencies_fake_mu_mumu.root", p3_props)

    plotlist.append(p3)

    # Add vertical reference lines for trigger thresholds.
    # Use the first histogram in the list to set the range

    refl_0 = TLine(26.0,p0_props["yAxisRange"][0],26.0,p0_props["yAxisRange"][1])
    Plot.reflines.append(refl_0)
    refl_1 = TLine(50.0,p0_props["yAxisRange"][0],50.0,p0_props["yAxisRange"][1])
    Plot.reflines.append(refl_1)

    multiP = MultiPlot( plots=plotlist )
    multiP.luminosity = luminosity
    multiP.buildLegend( header="#varepsilon_{fake} - Muons", lcoords=[0.5,0.7,0.93,0.91] )
    multiP.setCanvasCoords([50,50,800,800])
    # multiP.setPlotTitle("My title",(0.6,0.88))

    multiP.makeMultiPlot( "./PLOTS_25ns_v24/CombinedEfficiencies_Plots", "Fake_Mu_Compare" )

    # Clear before returning

    del plotlist[:]


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

    #p0 = Plot("Real_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_TruthTP/OutputPlots_MMClosureRates_TruthTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_25ns_v24_REBINNED/LeptonEfficiencies.root", p0_props)

    #plotlist.append(p0)

    p1_props = {
                "legend"      : "Tag & Probe (OF)",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.0,1.0),
                "colour"      : kBlue,
                "lineStyle"   : 1,
                "lineWidth"   : 2,
                "markerStyle" : 20,
                "markerSize"  : 1,
               }

    #p1 = Plot("Real_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIsoTagd0sig15_ForceProbeToBeFake_REBINNED/LeptonEfficiencies.root", p1_props)
    p1 = Plot("Real_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v27/OutputPlots_MMClosureRates_25ns_v27/LeptonEfficiencies.root", p1_props)

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

    #p2 = Plot("Real_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIsoTagd0sig15_ForceProbeToBeFake_TRIGMATCH_EFF_REBINNED/LeptonEfficiencies.root", p2_props)

    #plotlist.append(p2)

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

    #p3 = Plot("Real_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIsoTagd0sig15_ForceProbeToBeFake_NOT_TRIGMATCH_EFF_REBINNED/LeptonEfficiencies.root", p3_props)

    #plotlist.append(p3)

    # Add vertical reference lines for trigger thresholds.
    # Use the first histogram in the list to set the range

    # refl_0 = TLine(26.0,p0_props["yAxisRange"][0],26.0,p0_props["yAxisRange"][1])
    # Plot.reflines.append(refl_0)
    # refl_1 = TLine(60.0,p0_props["yAxisRange"][0],60.0,p0_props["yAxisRange"][1])
    # Plot.reflines.append(refl_1)
    # refl_2 = TLine(140.0,p0_props["yAxisRange"][0],140.0,p0_props["yAxisRange"][1])
    # Plot.reflines.append(refl_2)

    multiP = MultiPlot( plots=plotlist )
    multiP.luminosity = luminosity
    multiP.buildLegend( header="#varepsilon_{real} - Electrons", lcoords=[0.5,0.7,0.93,0.91] )
    multiP.setCanvasCoords([50,50,800,800])
    # multiP.setPlotTitle("My title",(0.6,0.88))

    outdir = "./PLOTS_25ns_v27/NonPrompt_VS_PhotonConv"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    multiP.makeMultiPlot( outdir, "Real_El_Compare_NonPrompt_PhotonConv" )

    #multiP.makeMultiPlot( "./PLOTS_25ns_v24/CombinedEfficiencies_Plots", "Real_El_Compare" )

    # Clear before returning

    del plotlist[:]


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

    p0 = Plot("Real_Mu_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_TruthTP/OutputPlots_MMClosureRates_TruthTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_25ns_v24_REBINNED/LeptonEfficiencies.root", p0_props)

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

    p1 = Plot("Real_Mu_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIsoTagd0sig15_ForceProbeToBeFake_REBINNED/LeptonEfficiencies.root", p1_props)

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

    p2 = Plot("Real_Mu_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIsoTagd0sig15_ForceProbeToBeFake_TRIGMATCH_EFF_REBINNED/LeptonEfficiencies.root", p2_props)

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

    p3 = Plot("Real_Mu_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v24/MMClosure_v24_SUSYTP/OutputPlots_MMClosureRates_SUSYTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_OF_AMBISOLVING_25ns_v24_TightTagIsoTagd0sig15_ForceProbeToBeFake_NOT_TRIGMATCH_EFF_REBINNED/LeptonEfficiencies.root", p3_props)

    plotlist.append(p3)

    # Add vertical reference lines for trigger thresholds.
    # Use the first histogram in the list to set the range

    refl_0 = TLine(26.0,p0_props["yAxisRange"][0],26.0,p0_props["yAxisRange"][1])
    Plot.reflines.append(refl_0)
    refl_1 = TLine(50.0,p0_props["yAxisRange"][0],50.0,p0_props["yAxisRange"][1])
    Plot.reflines.append(refl_1)

    multiP = MultiPlot( plots=plotlist )
    multiP.luminosity = luminosity
    multiP.buildLegend( header="#varepsilon_{real} - Muons", lcoords=[0.5,0.7,0.93,0.91] )
    multiP.setCanvasCoords([50,50,800,800])
    # multiP.setPlotTitle("My title",(0.6,0.88))

    multiP.makeMultiPlot( "./PLOTS_25ns_v24/CombinedEfficiencies_Plots", "Real_Mu_Compare" )

    # Clear before returning

    del plotlist[:]


def plotFakeElectron_NonPromptAndPhotonConvVSPhotonConv():

    plotlist = []

    p0_props = {
                "legend"      : "All fakes",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.0,0.7),
                "colour"      : kOrange+10,
                "lineStyle"   : 1,
                "lineWidth"   : 2,
                "markerStyle" : 20,
                "markerSize"  : 1,
               }

    p0 = Plot("3L_Fake_El_Pt_Efficiency_observed_sub","./PLOTS_25ns_v27/OutputPlots_MMRates_25ns_v27_3LTP_v2/LeptonEfficiencies.root", p0_props)

    plotlist.append(p0)

    p1_props = {
                "legend"      : "#gamma conversion",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.0,0.7),
                #"colour"      : kOrange+7,
                "colour"      : kGreen+3,
                "lineStyle"   : 1,
                "lineWidth"   : 2,
                "markerStyle" : 20,
                "markerSize"  : 1,
               }

    p1 = Plot("Fake_El_Pt_Efficiency_observed_sub","./PLOTS_25ns_v27/OutputPlots_PhotonConvElecRates_ZMuMuElTP_25ns_v27/LeptonEfficiencies.root", p1_props)

    plotlist.append(p1)

    # Add vertical reference lines for trigger thresholds.
    # Use the first histogram in the list to set the range

    refl_0 = TLine(26.0,p0_props["yAxisRange"][0],26.0,p0_props["yAxisRange"][1])
    #Plot.reflines.append(refl_0)
    refl_1 = TLine(60.0,p0_props["yAxisRange"][0],60.0,p0_props["yAxisRange"][1])
    #Plot.reflines.append(refl_1)
    #refl_2 = TLine(140.0,p0_props["yAxisRange"][0],140.0,p0_props["yAxisRange"][1])
    #Plot.reflines.append(refl_2)

    multiP = MultiPlot( plots=plotlist )

    outdir = "./PLOTS_25ns_v27/NonPromptAndPhotonconv_VS_PhotonConv_DD"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    multiP.makeMultiPlot( outdir, "Fake_El_Compare_NonPromptAndPhotonConv_PhotonConv" )
    multiP.luminosity = luminosity
    multiP.buildLegend( header="Data", lcoords=[0.5,0.7,0.93,0.91] )
    multiP.setCanvasCoords([50,50,800,800])
    # multiP.setPlotTitle("My title",(0.6,0.88))

    # Clear before returning

    del plotlist[:]


def plotFakeElectron_NonPromptVSPhotonConv():

    plotlist = []

    p0_props = {
                "legend"      : "Non-prompt (#mu^{#pm}e^{#pm}, CR)",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.0,0.1),
                "colour"      : kOrange+10,
                "lineStyle"   : 2,
                "lineWidth"   : 2,
                "markerStyle" : 20,
                "markerSize"  : 1.4,
                #"yaxisLog"    : True
               }

    p0 = Plot("Fake_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v29/OutputPlots_MMClosureRates_25ns_v29_LEP_FROM_JET/LeptonEfficiencies.root", p0_props)

    plotlist.append(p0)

    p1_props = {
                "legend"      : "#gamma conv. (#mu^{#pm}e^{#pm}, CR)",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.0,0.1),
                "colour"      : kOrange-2,
                "lineStyle"   : 2,
                "lineWidth"   : 2,
                "markerStyle" : 20,
                "markerSize"  : 1.4,
               }

    p1 = Plot("Fake_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v29/OutputPlots_MMClosureRates_25ns_v29_PHOTON_CONV/LeptonEfficiencies.root", p1_props)

    plotlist.append(p1)

    p2_props = {
                "legend"      : "<#varepsilon> (#mu^{#pm}e^{#pm}, CR)",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.0,0.1),
                "colour"      : kBlack,
                "lineStyle"   : 2,
                "lineWidth"   : 2,
                "markerStyle" : 20,
                "markerSize"  : 1.4,
               }

    p2 = Plot("Fake_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v29/OutputPlots_MMClosureRates_25ns_v29/LeptonEfficiencies_410501_410082_FixedRescaling_26_07_17.root", p2_props)

    plotlist.append(p2)

    p3_props = {
                "legend"      : "<#varepsilon> (e^{#pm}e^{#pm}, SR)",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.0,0.1),
                "colour"      : kBlue,
                "lineStyle"   : 1,
                "lineWidth"   : 3,
                "markerStyle" : 20,
                "markerSize"  : 1.4,
               }

    p3 = Plot("RESCALED_2L_ee_Fake_El_Pt_Efficiency_expectedbkg","./PLOTS_25ns_v29/OutputPlots_MMClosureRates_25ns_v29/LeptonEfficiencies_410501_410082_FixedRescaling_26_07_17.root", p3_props)

    plotlist.append(p3)

    multiP = MultiPlot( plots=plotlist )
    multiP.luminosity = luminosity
    multiP.buildLegend( header="Simulation (#it{t#bar{t}})", lcoords=[0.5,0.7,0.93,0.91] )
    multiP.setCanvasCoords([50,50,800,800])
    # multiP.setPlotTitle("My title",(0.6,0.88))

    outdir = "./PLOTS_25ns_v29/NonPrompt_VS_PhotonConv_TTBarTTGamma"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    multiP.makeMultiPlot( outdir, "Fake_El_Compare_NonPrompt_PhotonConv" )

    # Clear before returning

    del plotlist[:]

