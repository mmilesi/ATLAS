#!/usr/bin/env python

import os

from ROOT import kBlue, kOrange, kPink, kGreen, kRed, kYellow, kTeal, kMagenta, kViolet, kAzure, kCyan, kSpring, kGray, kBlack, kWhite
from ROOT import TLine

from EfficiencyPlotterClasses import Plot, MultiPlot

luminosity = 36.1
ATLASlabel = "Internal"
#ATLASlabel = ""

# eff_filepath = "./PLOTS_25ns_v29/OutputPlots_MMRates_25ns_v29_TopAsymmConv/LeptonEfficiencies.root"
eff_filepath = "./PLOTS_25ns_v29/OutputPlots_MMRates_25ns_v29_TopAsymmConv/LeptonEfficiencies_15GeV.root"
# eff_filepath = "./PLOTS_25ns_v29/OutputPlots_MMClosureRates_25ns_v29/LeptonEfficiencies.root"

# outdirbase =  "./PLOTS_25ns_v29/OutputPlots_MMRates_25ns_v29_TopAsymmConv/EfficiencyPlots/QualityPlots"
outdirbase =  "./PLOTS_25ns_v29/OutputPlots_MMRates_25ns_v29_TopAsymmConv/EfficiencyPlots/QualityPlots_15Gev"
# outdirbase =  "./PLOTS_25ns_v29/OutputPlots_MMClosureRates_25ns_v29/EfficiencyPlots/QualityPlots"

def plotRealEfficiency():

    plotlist = []

    p0_props = {
                "legend"      : "#bf{e}",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.45,1.1),
                "xAxisTitle"  : "p_{T} [GeV]",
                "xAxisLog"    : True,
                "colour"      : kBlue,
                "lineStyle"   : 2,
                "lineWidth"   : 3,
                "markerStyle" : 24,
                "markerSize"  : 1.4,
               }

    p0 = Plot("Real_El_Pt_Efficiency_observed_sub",eff_filepath, p0_props)
    p0.addSystematics(syslist=["ND_FakesOS"])

    plotlist.append(p0)

    p1_props = {
                "legend"      : "#bf{#mu}",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.45,1.1),
                "colour"      : kRed,
                "lineStyle"   : 1,
                "lineWidth"   : 3,
                "markerStyle" : 24,
                "markerSize"  : 1.4,
               }

    p1 = Plot("Real_Mu_Pt_Efficiency_observed_sub",eff_filepath, p1_props)
    p1.addSystematics(syslist=["ND_FakesOS"])

    plotlist.append(p1)

    multiP = MultiPlot( plots=plotlist )
    multiP.luminosity = luminosity
    #multiP.buildLegend( header="#it{Data} - e^{#pm}#mu^{#mp}, #mu^{#pm}e^{#mp}" )
    multiP.buildLegend( header="#bf{Data - e^{#pm}#mu^{#mp}, #mu^{#pm}e^{#mp}}", lcoords=[0.62,0.3,0.9,0.5] ) #Unbolded
    multiP.ATLASlabel = ATLASlabel
    multiP.setCanvasCoords([50,50,800,800])
    # multiP.setPlotTitle("Real efficiency",(0.6,0.88))

    outdir = outdirbase + ATLASlabel
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    multiP.makeMultiPlot( outdir, "RealEfficiency_ElMu" )

    # Clear before returning

    del plotlist[:]


def plotFakeEfficiency_El():

    plotlist = []

    p0_props = {
                "legend"      : "#bf{N_{b-tags} = 1}",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.0,0.15),
                "xAxisTitle"  : "p_{T} [GeV]",
                #"xAxisLog"    : True,
                "colour"      : kBlue+2,
                "lineStyle"   : 2,
                "lineWidth"   : 3,
                "markerStyle" : 24,
                "markerSize"  : 1.4,
               }

    p0 = Plot("RESCALED_2L_ee_Fake_El_NBJets_VS_Pt_Efficiency_observed_sub_projPt_1",eff_filepath, p0_props)
    p0.addSystematics(syslist=["ND_VV","ND_TTV","ND_OtherPromptSS","N_QMisID","D_QMisID","ND_ALPHA"])

    plotlist.append(p0)

    p1_props = {
                "legend"      : "#bf{N_{b-tags} #geq 2}",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.0,0.2),
                "colour"      : kBlue-7,
                "lineStyle"   : 1,
                "lineWidth"   : 3,
                "markerStyle" : 24,
                "markerSize"  : 1.4,
               }

    p1 = Plot("RESCALED_2L_ee_Fake_El_NBJets_VS_Pt_Efficiency_observed_sub_projPt_2",eff_filepath, p1_props)
    p1.addSystematics(syslist=["ND_VV","ND_TTV","ND_OtherPromptSS","N_QMisID","D_QMisID","ND_ALPHA"])

    plotlist.append(p1)

    multiP = MultiPlot( plots=plotlist )
    multiP.luminosity = luminosity
    # multiP.buildLegend( header="#bf{Data - #mu^{#pm}e^{#pm}}",lcoords=[0.65,0.7,0.93,0.9] ) # Unbolded
    multiP.buildLegend( header="#bf{Data - #mu^{#pm}e^{#pm}}",lcoords=[0.64,0.72,0.92,0.928] ) # Unbolded

    multiP.ATLASlabel = ATLASlabel
    multiP.setCanvasCoords([50,50,800,800])
    # multiP.setPlotTitle("Electron fake rate",(0.2,0.7))

    outdir = outdirbase + ATLASlabel
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    multiP.makeMultiPlot( outdir, "FakeRate_El" )

    # Clear before returning

    del plotlist[:]


def plotFakeEfficiency_Mu():

    plotlist = []

    p0_props = {
                "legend"      : "#bf{0 #leq min(#Delta R_{#mu,j}) #leq 1}",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.0,1.1),
                "xAxisTitle"  : "p_{T} [GeV]",
                #"xAxisLog"    : True,
                "colour"      : kRed+1,
                "lineStyle"   : 2,
                "lineWidth"   : 3,
                "markerStyle" : 24,
                "markerSize"  : 1.4,
               }

    p0 = Plot("Fake_Mu_DistanceClosestJet_VS_Pt_Efficiency_observed_sub_projPt_1",eff_filepath, p0_props)
    p0.addSystematics(syslist=["ND_VV","ND_TTV","ND_OtherPromptSS"])

    plotlist.append(p0)

    p1_props = {
                "legend"      : "#bf{1 < min(#Delta R_{#mu,j}) #leq 5}",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.0,1.0),
                "colour"      : kRed-9,
                "lineStyle"   : 1,
                "lineWidth"   : 3,
                "markerStyle" : 24,
                "markerSize"  : 1.4,
               }

    p1 = Plot("Fake_Mu_DistanceClosestJet_VS_Pt_Efficiency_observed_sub_projPt_2",eff_filepath, p1_props)
    p1.addSystematics(syslist=["ND_VV","ND_TTV","ND_OtherPromptSS"])

    plotlist.append(p1)

    multiP = MultiPlot( plots=plotlist )
    multiP.luminosity = luminosity
    #multiP.buildLegend( header="#bf{Data - #mu^{#pm}#mu^{#pm}}",lcoords=[0.53,0.7,0.95,0.9] ) # Unbolded
    multiP.buildLegend( header="#bf{Data - #mu^{#pm}#mu^{#pm}}",lcoords=[0.53,0.728,0.95,0.928] ) # Unbolded
    multiP.ATLASlabel = ATLASlabel
    multiP.setCanvasCoords([50,50,800,800])
    # multiP.setPlotTitle("Muon fake rate",(0.2,0.7))

    outdir = outdirbase + ATLASlabel
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    multiP.makeMultiPlot( outdir, "FakeRate_Mu" )

    # Clear before returning

    del plotlist[:]

# ----------------------------------------------------------------

def plotRealEfficiency_Closure():

    plotlist = []

    p0_props = {
                "legend"      : "#bf{e}",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.5,1.1),
                "xAxisTitle"  : "p_{T} [GeV]",
                "xAxisLog"    : True,
                "colour"      : kBlue,
                "lineStyle"   : 2,
                "lineWidth"   : 3,
                "markerStyle" : 24,
                "markerSize"  : 1.4,
               }

    p0 = Plot("Real_El_Pt_Efficiency_expectedbkg",eff_filepath, p0_props)

    plotlist.append(p0)

    p1_props = {
                "legend"      : "#bf{#mu}",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.5,1.1),
                "colour"      : kRed,
                "lineStyle"   : 1,
                "lineWidth"   : 3,
                "markerStyle" : 24,
                "markerSize"  : 1.4,
               }

    p1 = Plot("Real_Mu_Pt_Efficiency_expectedbkg",eff_filepath, p1_props)

    plotlist.append(p1)

    multiP = MultiPlot( plots=plotlist )
    multiP.luminosity = luminosity
    multiP.buildLegend( header="#bf{#it{t#bar{t}} - e^{#pm}#mu^{#mp}, #mu^{#pm}e^{#mp}}", lcoords=[0.62,0.3,0.9,0.5] )
    multiP.ATLASlabel = ATLASlabel
    multiP.setCanvasCoords([50,50,800,800])
    # multiP.setPlotTitle("Real efficiency",(0.6,0.88))

    outdir = outdirbase + ATLASlabel
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    multiP.makeMultiPlot( outdir, "RealEfficiency_ElMu" )

    # Clear before returning

    del plotlist[:]


def plotFakeEfficiency_Closure():

    plotlist = []

    p0_props = {
                "legend"      : "#bf{e}",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.0,0.065),
                "xAxisTitle"  : "p_{T} [GeV]",
                #"xAxisLog"    : True,
                "colour"      : kBlue,
                "lineStyle"   : 2,
                "lineWidth"   : 3,
                "markerStyle" : 24,
                "markerSize"  : 1.4,
               }

    p0 = Plot("RESCALED_2L_ee_Fake_El_NBJets_VS_Pt_Efficiency_expectedbkg_projPt_inclusive",eff_filepath, p0_props)

    plotlist.append(p0)

    p1_props = {
                "legend"      : "#bf{#mu}",
                "yAxisTitle"  : "#varepsilon",
                "yAxisRange"  : (0.0,0.065),
                "colour"      : kRed,
                "lineStyle"   : 1,
                "lineWidth"   : 3,
                "markerStyle" : 24,
                "markerSize"  : 1.4,
               }

    p1 = Plot("Fake_Mu_DistanceClosestJet_VS_Pt_Efficiency_expectedbkg_projPt_inclusive",eff_filepath, p1_props)

    plotlist.append(p1)

    multiP = MultiPlot( plots=plotlist )
    multiP.luminosity = luminosity
    multiP.buildLegend( header="#bf{#it{t#bar{t}} - Fake CRs}",lcoords=[0.53,0.728,0.95,0.928] ) # Unbolded
    multiP.ATLASlabel = ATLASlabel
    multiP.setCanvasCoords([50,50,800,800])
    # multiP.setPlotTitle("Muon fake rate",(0.2,0.7))

    outdir = outdirbase + ATLASlabel
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    multiP.makeMultiPlot( outdir, "FakeRate_ElMu" )

    # Clear before returning

    del plotlist[:]


