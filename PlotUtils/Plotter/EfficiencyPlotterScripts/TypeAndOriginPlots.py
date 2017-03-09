#!/usr/bin/env python

import os

from ROOT import kBlue, kOrange, kPink, kGreen, kRed, kYellow, kTeal, kMagenta, kViolet, kAzure, kCyan, kSpring, kGray, kBlack, kWhite
from ROOT import TLine, TCanvas

from EfficiencyPlotterClasses import Plot, MultiPlot

Plot.luminosity = 36.4702

def plotTypeVSOrigin( normFactor=0, **kwargs ):

    plotlist = []

    if kwargs["flavour"] == "Mu":
        header   = "muons"
    elif kwargs["flavour"] == "El":
        header   = "electrons"

    Plot.legend.SetHeader("Fake {0}".format(header))

    variable = "{0}ProbeType_VS_{1}ProbeOrigin".format(kwargs["flavour"],kwargs["flavour"])

    #basePath = os.path.abspath(os.curdir) + "/" + "PLOTS_" + kwargs["prodID"] + "/MMRates_DATA/" + "OutputPlots_MMRates_DATAMC_" + kwargs["prodID"] + "_LeptonMVA_FakeOrig_" + kwargs["jetSelection"]
    basePath = os.path.abspath(os.curdir) + "/" + "OutputPlots_MMRates_DATAMC_" + kwargs["prodID"] + "_LeptonMVA_FakeOrig_" + kwargs["jetSelection"]

    inputPath = basePath + "/" + kwargs["lepSelection"] + "/"
    inputName = kwargs["lepSelection"] + "_" + variable + ".root"

    p0_props = {
                #"legend"      : "t#bar{t} - Loose sel.",
                "xAxisTitle"  : "truthType",
                "yAxisTitle"  : "truthOrigin",
                "xAxisLabels" : [("Unknown",0),("UnknownElectron",1),("IsoElectron",2),("NonIsoElectron",3),("BkgElectron",4),("UnknownMuon",5),("IsoMuon",6),("NonIsoMuon",7),("BkgMuon",8),("UnknownTau",9),("IsoTau",10),("NonIsoTau",11),("BkgTau",12),("UnknownPhoton",13),("IsoPhoton",14),("NonIsoPhoton",15),("BkgPhoton",16),("Hadron",17),("Neutrino",18),("NuclFrag",19),("NonPrimary",20)], #("GenParticle",21),("SUSYParticle",22),("BBbarMesonPart",23),("BottomMesonPart",24),("CCbarMesonPart",25),("CharmedMesonPart",26),("BottomBaryonPart",27),("CharmedBaryonPart",28),("StrangeBaryonPart",29),("LightBaryonPart",30),("StrangeMesonPart",31),("LightMesonPart",32),("BJet",33),("CJet",34),("LJet",35),("GJet",36),("TauJet",37),("UnknownJet",38)],
                "yAxisLabels" : [("NonDefined",0),("SingleElec",1),("SingleMuon",2),("SinglePhot",3),("SingleTau",4),("PhotonConv",5),("DalitzDec",6),("ElMagProc",7),("Mu",8),("TauLep",9),("top",10),("QuarkWeakDec",11),("WBoson",12),("ZBoson",13),("Higgs",14),("HiggsMSSM",15),("HeavyBoson",16),("WBosonLRSM",17),("NuREle",18),("NuRMu",19),("NuRTau",20),("LQ",21),("SUSY",22),("LightMeson",23),("StrangeMeson",24),("CharmedMeson",25),("BottomMeson",26),("CCbarMeson",27),("JPsi",28),("BBbarMeson",29),("LightBaryon",30),("StrangeBaryon",31),("CharmedBaryon",32),("BottomBaryon",33),("PionDecay",34),("KaonDecay",35),("BremPhot",36),("PromptPhot",37),("UndrPhot",38),("ISRPhot",39),("FSRPhot",40)], #,("NucReact",41),("PiZero",42),("DiBoson",43),("ZorHeavyBoson",44),("QCD",45)],
                "normFactor"   : normFactor,
                "drawGrid"     : True,
               }

    p0 = Plot(kwargs["sample"], inputPath  + inputName, p0_props)

    extension = ""
    if "normFactor" in p0_props:
        if p0_props["normFactor"]:
            extension = "_NORM"

    canvasID = kwargs["lepSelection"] + "_" + variable + "_" + kwargs["sample"] + "_" + kwargs["prodID"] + "_" + kwargs["jetSelection"] + extension

    c = TCanvas(canvasID,canvasID,50,50,800,900)
    c.cd()
    p0.makePlot()

    #Plot.legend.Draw()
    #Plot.legendATLAS.DrawLatex(0.6,0.35,"#bf{#it{ATLAS}} Work In Progress")
    #Plot.legendLumi.DrawLatex(0.6,0.27,"#sqrt{{s}} = 13 TeV, #int L dt = {0:.1f} fb^{{-1}}".format(Plot.luminosity))

    savePath = basePath + "/TypeOriginPlots_" + kwargs["prodID"] + "/"
    saveName = canvasID

    for ext in ["png","eps","root"]:
        c.SaveAs( savePath + saveName + "." + ext )

    Plot.legend.Clear()


def plotTypeVSNjets( normFactor=0, **kwargs ):

    plotlist = []

    if kwargs["flavour"] == "Mu":
        header   = "muons"
    elif kwargs["flavour"] == "El":
        header   = "electrons"

    Plot.legend.SetHeader("Fake {0}".format(header))

    variable = "{0}ProbeType_VS_NJets".format(kwargs["flavour"])

    #basePath = os.path.abspath(os.curdir) + "/" + "PLOTS_" + kwargs["prodID"] + "/MMRates_DATA/" + "OutputPlots_MMRates_DATAMC_" + kwargs["prodID"] + "_LeptonMVA_FakeOrig_" + kwargs["jetSelection"]
    basePath = os.path.abspath(os.curdir) + "/" + "OutputPlots_MMRates_DATAMC_" + kwargs["prodID"] + "_LeptonMVA_FakeOrig_" + kwargs["jetSelection"]

    inputPath = basePath + "/" + kwargs["lepSelection"] + "/"
    inputName = kwargs["lepSelection"] + "_" + variable + ".root"

    p0_props = {
                #"legend"      : "t#bar{t} - Loose sel.",
                "xAxisTitle"  : "truthType",
                "yAxisTitle"  : "Jet Multiplicity",
                "xAxisLabels" : [("Unknown",0),("UnknownElectron",1),("IsoElectron",2),("NonIsoElectron",3),("BkgElectron",4),("UnknownMuon",5),("IsoMuon",6),("NonIsoMuon",7),("BkgMuon",8),("UnknownTau",9),("IsoTau",10),("NonIsoTau",11),("BkgTau",12),("UnknownPhoton",13),("IsoPhoton",14),("NonIsoPhoton",15),("BkgPhoton",16),("Hadron",17),("Neutrino",18),("NuclFrag",19),("NonPrimary",20)], #("GenParticle",21),("SUSYParticle",22),("BBbarMesonPart",23),("BottomMesonPart",24),("CCbarMesonPart",25),("CharmedMesonPart",26),("BottomBaryonPart",27),("CharmedBaryonPart",28),("StrangeBaryonPart",29),("LightBaryonPart",30),("StrangeMesonPart",31),("LightMesonPart",32),("BJet",33),("CJet",34),("LJet",35),("GJet",36),("TauJet",37),("UnknownJet",38)],
                "normFactor"   : normFactor,
                "drawGrid"     : True,
               }

    p0 = Plot(kwargs["sample"], inputPath  + inputName, p0_props)

    extension = ""
    if "normFactor" in p0_props:
        if p0_props["normFactor"]:
            extension = "_NORM"

    canvasID = kwargs["lepSelection"] + "_" + variable + "_" + kwargs["sample"] + "_" + kwargs["prodID"] + "_" + kwargs["jetSelection"] + extension

    c = TCanvas(canvasID,canvasID,50,50,900,700)
    c.cd()
    p0.makePlot()

    #Plot.legend.Draw()
    #Plot.legendATLAS.DrawLatex(0.6,0.35,"#bf{#it{ATLAS}} Work In Progress")
    #Plot.legendLumi.DrawLatex(0.6,0.27,"#sqrt{{s}} = 13 TeV, #int L dt = {0:.1f} fb^{{-1}}".format(Plot.luminosity))

    savePath = basePath + "/TypeOriginPlots_" + kwargs["prodID"] + "/"
    saveName = canvasID

    for ext in ["png","eps","root"]:
        c.SaveAs( savePath + saveName + "." + ext )

    Plot.legend.Clear()


def plotOriginVSNjets( normFactor=0, **kwargs ):

    plotlist = []

    if kwargs["flavour"] == "Mu":
        header   = "muons"
    elif kwargs["flavour"] == "El":
        header   = "electrons"

    Plot.legend.SetHeader("Fake {0}".format(header))

    variable = "{0}ProbeOrigin_VS_NJets".format(kwargs["flavour"])

    #basePath = os.path.abspath(os.curdir) + "/" + "PLOTS_" + kwargs["prodID"] + "/MMRates_DATA/" + "OutputPlots_MMRates_DATAMC_" + kwargs["prodID"] + "_LeptonMVA_FakeOrig_" + kwargs["jetSelection"]
    basePath = os.path.abspath(os.curdir) + "/" + "OutputPlots_MMRates_DATAMC_" + kwargs["prodID"] + "_LeptonMVA_FakeOrig_" + kwargs["jetSelection"]

    inputPath = basePath + "/" + kwargs["lepSelection"] + "/"
    inputName = kwargs["lepSelection"] + "_" + variable + ".root"

    p0_props = {
                #"legend"      : "t#bar{t} - Loose sel.",
                "xAxisTitle"  : "truthOrigin",
                "yAxisTitle"  : "Jet multiplicity",
                "xAxisLabels" : [("NonDefined",0),("SingleElec",1),("SingleMuon",2),("SinglePhot",3),("SingleTau",4),("PhotonConv",5),("DalitzDec",6),("ElMagProc",7),("Mu",8),("TauLep",9),("top",10),("QuarkWeakDec",11),("WBoson",12),("ZBoson",13),("Higgs",14),("HiggsMSSM",15),("HeavyBoson",16),("WBosonLRSM",17),("NuREle",18),("NuRMu",19),("NuRTau",20),("LQ",21),("SUSY",22),("LightMeson",23),("StrangeMeson",24),("CharmedMeson",25),("BottomMeson",26),("CCbarMeson",27),("JPsi",28),("BBbarMeson",29),("LightBaryon",30),("StrangeBaryon",31),("CharmedBaryon",32),("BottomBaryon",33),("PionDecay",34),("KaonDecay",35),("BremPhot",36),("PromptPhot",37),("UndrPhot",38),("ISRPhot",39),("FSRPhot",40)], #,("NucReact",41),("PiZero",42),("DiBoson",43),("ZorHeavyBoson",44),("QCD",45)],
                "normFactor"   : normFactor,
                "drawGrid"     : True,
               }

    p0 = Plot(kwargs["sample"], inputPath  + inputName, p0_props)

    extension = ""
    if "normFactor" in p0_props:
        if p0_props["normFactor"]:
            extension = "_NORM"

    canvasID = kwargs["lepSelection"] + "_" + variable + "_" + kwargs["sample"] + "_" + kwargs["prodID"] + "_" + kwargs["jetSelection"] + extension
    c = TCanvas(canvasID,canvasID,50,50,1000,700)
    c.cd()
    p0.makePlot()

    doLeptonOriginFracPlots = ( not "normFactor" in p0_props or not p0_props["normFactor"] )
    if doLeptonOriginFracPlots:
        canvasID_stack = kwargs["lepSelection"] + "_LepOriginFrac_VS_NJets_" + kwargs["sample"] + "_" + kwargs["prodID"] + "_" + kwargs["jetSelection"] + extension
        cstack = TCanvas(canvasID_stack,canvasID_stack,50,50,1000,700)
        cstack.cd()
        mystack, mystackleg = p0.makeLeptonOriginFracPlots(canvasID_stack)
        mystack.Draw()
        mystack.GetXaxis().SetTitle("Jet multiplicity")
        mystack.GetYaxis().SetTitle("Fake lepton origin fraction")
        mystack.GetXaxis().SetRangeUser(2,mystack.GetHistogram().GetNbinsX()-1)
        mystackleg.Draw()
        saveName_stack = canvasID_stack

    #Plot.legend.Draw()
    #Plot.legendATLAS.DrawLatex(0.6,0.35,"#bf{#it{ATLAS}} Work In Progress")
    #Plot.legendLumi.DrawLatex(0.6,0.27,"#sqrt{{s}} = 13 TeV, #int L dt = {0:.1f} fb^{{-1}}".format(Plot.luminosity))

    savePath = basePath + "/TypeOriginPlots_" + kwargs["prodID"] + "/"
    saveName = canvasID

    for ext in ["png","eps","root"]:
        c.SaveAs( savePath + saveName + "." + ext )
        if doLeptonOriginFracPlots:
            cstack.SaveAs( savePath + saveName_stack + "." + ext )

    Plot.legend.Clear()


