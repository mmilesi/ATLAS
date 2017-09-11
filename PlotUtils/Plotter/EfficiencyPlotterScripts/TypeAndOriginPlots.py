#!/usr/bin/env python

import os

from ROOT import kBlue, kOrange, kPink, kGreen, kRed, kYellow, kTeal, kMagenta, kViolet, kAzure, kCyan, kSpring, kGray, kBlack, kWhite
from ROOT import TLine, TCanvas, TFile, TH1D, THStack, TLegend

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

    variable_2L = "Lep1Origin_VS_NJets"
    variable_3L = "Lep2Origin_VS_NJets"

    basePath = os.path.abspath(os.curdir) + "/" + "PLOTS_" + kwargs["prodID"] + "_v2/" + "OutputPlots_FakeOriginFrac_TTBarTTBarGamma_" + kwargs["prodID"]

    inputPath_2L   = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_ALLNJ_VR_TT" + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_ALLNJ_VR_TT_"
    inputPath_3LSR = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_3L_SR_TT" + "/" + kwargs["flavour"] + "FakeOriginFrac_3L_SR_TT_"

    inputName_2L = variable_2L + ".root"
    inputName_3L = variable_3L + ".root"

    p0_props = {
                #"legend"      : "t#bar{t} - Loose sel.",
                "xAxisTitle"  : "truthOrigin",
                "yAxisTitle"  : "N_{jets}",
                "xAxisLabels" : [("NonDefined",0),("SingleElec",1),("SingleMuon",2),("SinglePhot",3),("SingleTau",4),("PhotonConv",5),("DalitzDec",6),("ElMagProc",7),("Mu",8),("TauLep",9),("top",10),("QuarkWeakDec",11),("WBoson",12),("ZBoson",13),("Higgs",14),("HiggsMSSM",15),("HeavyBoson",16),("WBosonLRSM",17),("NuREle",18),("NuRMu",19),("NuRTau",20),("LQ",21),("SUSY",22),("LightMeson",23),("StrangeMeson",24),("CharmedMeson",25),("BottomMeson",26),("CCbarMeson",27),("JPsi",28),("BBbarMeson",29),("LightBaryon",30),("StrangeBaryon",31),("CharmedBaryon",32),("BottomBaryon",33),("PionDecay",34),("KaonDecay",35),("BremPhot",36),("PromptPhot",37),("UndrPhot",38),("ISRPhot",39),("FSRPhot",40)], #,("NucReact",41),("PiZero",42),("DiBoson",43),("ZorHeavyBoson",44),("QCD",45)],
                "normFactor"   : normFactor,
                "drawGrid"     : True,
               }

    extension = ""
    if "normFactor" in p0_props:
        if p0_props["normFactor"]:
            extension = "_NORM"
    doLeptonOriginFracPlots = ( not "normFactor" in p0_props or not p0_props["normFactor"] )

    pathlist   = [inputPath_2L+inputName_2L,inputPath_3LSR+inputName_3L]
    varlist    = [variable_2L,variable_3L]
    regionlist = ["2LSS_ALLNJ_VR_TT","3L_SR_TT"]

    for idx, path in enumerate(pathlist):

        p0 = Plot(kwargs["sample"], path, p0_props)

        canvasID = varlist[idx] + "_" + kwargs["flavour"] + "_" + kwargs["prodID"] + "_" + regionlist[idx] + extension
        c = TCanvas(canvasID,canvasID,50,50,1000,700)
        c.cd()
        p0.makePlot()

        if doLeptonOriginFracPlots:
            canvasID_stack = "FakeLepOriginFrac" + "_" + kwargs["flavour"] + "_" + varlist[idx] + "_" + kwargs["prodID"] + "_" + regionlist[idx]
            cstack = TCanvas(canvasID_stack,canvasID_stack,50,50,1000,700)
            cstack.cd()
            mystack, mystackleg = p0.makeLeptonOriginFracPlots(canvasID_stack)
            mystack.Draw()
            mystack.GetXaxis().SetTitle("N_{jets}")
            mystack.GetYaxis().SetTitle("Fake {0} origin fraction".format(header))
            mystackleg.Draw()
            saveName_stack = canvasID_stack

        # Plot.legend.Draw()
        # Plot.legendATLAS.DrawLatex(0.6,0.35,"#bf{#it{ATLAS}} Work In Progress")
        # Plot.legendLumi.DrawLatex(0.6,0.27,"#sqrt{{s}} = 13 TeV, #int L dt = {0:.1f} fb^{{-1}}".format(Plot.luminosity))

        savePath = basePath + "/" + kwargs["flavour"] + "_" + "FakeLepOriginFrac_VS_NJets_" + kwargs["prodID"] + "/"

        if not os.path.exists(savePath):
            os.makedirs(savePath)

        saveName = canvasID
        for ext in ["png","pdf","root"]:
            c.SaveAs( savePath + saveName + "." + ext )
            if doLeptonOriginFracPlots:
                cstack.SaveAs( savePath + saveName_stack + "." + ext )

    Plot.legend.Clear()

def plotOriginVSNBjets( normFactor=0, **kwargs ):

    plotlist = []

    if kwargs["flavour"] == "Mu":
        header   = "muons"
    elif kwargs["flavour"] == "El":
        header   = "electrons"

    Plot.legend.SetHeader("Fake {0}".format(header))

    variable_2L = "Lep1Origin_VS_NBJets"
    variable_3L = "Lep2Origin_VS_NBJets"

    basePath = os.path.abspath(os.curdir) + "/" + "PLOTS_" + kwargs["prodID"] + "_v2/" + "OutputPlots_FakeOriginFrac_TTBarTTBarGamma_" + kwargs["prodID"]

    pathlist = []

    inputName_2L = variable_2L + ".root"
    inputName_3L = variable_3L + ".root"

    inputPath_2LVR = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_LOWNJ_VR_TT" + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_"
    pathlist.append(inputPath_2LVR+inputName_2L)
    inputPath_2LSR = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_SR_TT" + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_SR_TT_"
    pathlist.append(inputPath_2LSR+inputName_2L)
    inputPath_3LSR = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_3L_SR_TT" + "/" + kwargs["flavour"] + "FakeOriginFrac_3L_SR_TT_"
    pathlist.append(inputPath_3LSR+inputName_3L)
    if kwargs["flavour"] == "El":
        inputPath_2LVR_ElEl = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_ElEl" + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_ElEl_"
        pathlist.append(inputPath_2LVR_ElEl+inputName_2L)
        inputPath_2LVR_OF   = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_OF" + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_OF_"
        pathlist.append(inputPath_2LVR_OF+inputName_2L)
        inputPath_2LSR_ElEl = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_SR_TT_ElEl" + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_SR_TT_ElEl_"
        pathlist.append(inputPath_2LSR_ElEl+inputName_2L)
        inputPath_2LSR_OF   = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_SR_TT_OF" + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_SR_TT_OF_"
        pathlist.append(inputPath_2LSR_OF+inputName_2L)

    p0_props = {
                #"legend"      : "t#bar{t} - Loose sel.",
                "xAxisTitle"  : "truthOrigin",
                "yAxisTitle"  : "N_{b-tags}",
                "xAxisLabels" : [("NonDefined",0),("SingleElec",1),("SingleMuon",2),("SinglePhot",3),("SingleTau",4),("PhotonConv",5),("DalitzDec",6),("ElMagProc",7),("Mu",8),("TauLep",9),("top",10),("QuarkWeakDec",11),("WBoson",12),("ZBoson",13),("Higgs",14),("HiggsMSSM",15),("HeavyBoson",16),("WBosonLRSM",17),("NuREle",18),("NuRMu",19),("NuRTau",20),("LQ",21),("SUSY",22),("LightMeson",23),("StrangeMeson",24),("CharmedMeson",25),("BottomMeson",26),("CCbarMeson",27),("JPsi",28),("BBbarMeson",29),("LightBaryon",30),("StrangeBaryon",31),("CharmedBaryon",32),("BottomBaryon",33),("PionDecay",34),("KaonDecay",35),("BremPhot",36),("PromptPhot",37),("UndrPhot",38),("ISRPhot",39),("FSRPhot",40)], #,("NucReact",41),("PiZero",42),("DiBoson",43),("ZorHeavyBoson",44),("QCD",45)],
                "normFactor"   : normFactor,
                "drawGrid"     : True,
               }

    extension = ""
    if "normFactor" in p0_props:
        if p0_props["normFactor"]:
            extension = "_NORM"
    doLeptonOriginFracPlots = ( not "normFactor" in p0_props or not p0_props["normFactor"] )

    varlist    = [variable_2L,variable_2L,variable_3L]
    regionlist = ["2LSS_LOWNJ_VR_TT","2LSS_SR_TT","3L_SR_TT"]
    if kwargs["flavour"] == "El":
        regionlist.extend(["2LSS_LOWNJ_VR_TT_ElEl","2LSS_LOWNJ_VR_TT_OF","2LSS_SR_TT_ElEl","2LSS_SR_TT_OF"])
        varlist.extend([variable_2L,variable_2L,variable_2L,variable_2L])

    for idx, path in enumerate(pathlist):

        p0 = Plot(kwargs["sample"], path, p0_props)

        canvasID = varlist[idx] + "_" + kwargs["flavour"] + "_" + kwargs["prodID"] + "_" + regionlist[idx] + extension
        c = TCanvas(canvasID,canvasID,50,50,1000,700)
        c.cd()
        p0.makePlot()

        if doLeptonOriginFracPlots:
            canvasID_stack = "FakeLepOriginFrac" + "_" + kwargs["flavour"] + "_" + varlist[idx] + "_" + kwargs["prodID"] + "_" + regionlist[idx]
            cstack = TCanvas(canvasID_stack,canvasID_stack,50,50,1000,700)
            cstack.cd()
            mystack, mystackleg = p0.makeLeptonOriginFracPlots(canvasID_stack)
            mystack.Draw()
            mystack.GetXaxis().SetTitle("N_{b-tags}")
            mystack.GetYaxis().SetTitle("Fake {0} origin fraction".format(header))
            mystackleg.Draw()
            saveName_stack = canvasID_stack

        # Plot.legend.Draw()
        # Plot.legendATLAS.DrawLatex(0.6,0.35,"#bf{#it{ATLAS}} Work In Progress")
        # Plot.legendLumi.DrawLatex(0.6,0.27,"#sqrt{{s}} = 13 TeV, #int L dt = {0:.1f} fb^{{-1}}".format(Plot.luminosity))

        savePath = basePath + "/" + kwargs["flavour"] + "_" + "FakeLepOriginFrac_VS_NBJets_" + kwargs["prodID"] + "/"

        if not os.path.exists(savePath):
            os.makedirs(savePath)

        saveName = canvasID
        for ext in ["png","pdf","root"]:
            c.SaveAs( savePath + saveName + "." + ext )
            if doLeptonOriginFracPlots:
                cstack.SaveAs( savePath + saveName_stack + "." + ext )

    Plot.legend.Clear()


def plotOriginVSPt( normFactor=0, **kwargs ):

    plotlist = []

    header = "leptons"

    Plot.legend.SetHeader("Fake {0}".format(header))

    variable_2L = "LepFakeOrigin_VS_LepFakePt"

    basePath = os.path.abspath(os.curdir) + "/" + "PLOTS_" + kwargs["prodID"] + "/" + "OutputPlots_FakeOriginFrac_TTBarTTBarGamma_" + kwargs["prodID"] + "_SplitFlavours_TruthTP"

    pathlist    = []
    varlist     = []
    regionlist  = []
    flaglist    = []

    inputName_2L = variable_2L + ".root"

    inputPath_2LVR_of = basePath + "/" + "LepFake" + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_MuEl" + "/" + "LepFake" + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_MuEl_"
    pathlist.append(inputPath_2LVR_of+inputName_2L)
    varlist.append(variable_2L)
    regionlist.append("2LSS_LOWNJ_VR_TT")
    flaglist.append(("me (OF), CR","MuEl_CR"))

    inputPath_2LSR_ee = basePath + "/" + "ElFake" + "FakeOriginFrac_2LSS_SR_TT_ElEl" + "/" + "ElFake" + "FakeOriginFrac_2LSS_SR_TT_ElEl_"
    pathlist.append(inputPath_2LSR_ee+inputName_2L)
    varlist.append(variable_2L)
    regionlist.append("2LSS_SR_TT")
    flaglist.append(("ee, SR","ElEl_SR"))

    inputPath_2LSR_of= basePath + "/" + "LepFake" + "FakeOriginFrac_2LSS_SR_TT_OF" + "/" + "LepFake" + "FakeOriginFrac_2LSS_SR_TT_OF_"
    pathlist.append(inputPath_2LSR_of+inputName_2L)
    varlist.append(variable_2L)
    regionlist.append("2LSS_SR_TT")
    flaglist.append(("OF, SR","OF_SR"))

    inputPath_2LVR_ee = basePath + "/" + "ElFake" + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_ElEl" + "/" + "ElFake" + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_ElEl_"
    pathlist.append(inputPath_2LVR_ee+inputName_2L)
    varlist.append(variable_2L)
    regionlist.append("2LSS_LOWNJ_VR_TT")
    flaglist.append(("ee, CR","ElEl_CR"))

    inputPath_2LVR_mm = basePath + "/" + "MuFake" + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_MuMu" + "/" + "MuFake" + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_MuMu_"
    pathlist.append(inputPath_2LVR_mm+inputName_2L)
    varlist.append(variable_2L)
    regionlist.append("2LSS_LOWNJ_VR_TT")
    flaglist.append(("mm, CR","MuMu_CR"))

    inputPath_2LSR_mm = basePath + "/" + "MuFake" + "FakeOriginFrac_2LSS_SR_TT_MuMu" + "/" + "MuFake" + "FakeOriginFrac_2LSS_SR_TT_MuMu_"
    pathlist.append(inputPath_2LSR_mm+inputName_2L)
    varlist.append(variable_2L)
    regionlist.append("2LSS_LOWNJ_VR_TT")
    flaglist.append(("mm, SR","MuMu_SR"))

    p0_props = {
                #"legend"      : "t#bar{t} - Loose sel.",
                "xAxisTitle"  : "truthOrigin",
                "yAxisTitle"  : "p_{T}^{l}",
                "xAxisLabels" : [("NonDefined",0),("SingleElec",1),("SingleMuon",2),("SinglePhot",3),("SingleTau",4),("PhotonConv",5),("DalitzDec",6),("ElMagProc",7),("Mu",8),("TauLep",9),("top",10),("QuarkWeakDec",11),("WBoson",12),("ZBoson",13),("Higgs",14),("HiggsMSSM",15),("HeavyBoson",16),("WBosonLRSM",17),("NuREle",18),("NuRMu",19),("NuRTau",20),("LQ",21),("SUSY",22),("LightMeson",23),("StrangeMeson",24),("CharmedMeson",25),("BottomMeson",26),("CCbarMeson",27),("JPsi",28),("BBbarMeson",29),("LightBaryon",30),("StrangeBaryon",31),("CharmedBaryon",32),("BottomBaryon",33),("PionDecay",34),("KaonDecay",35),("BremPhot",36),("PromptPhot",37),("UndrPhot",38),("ISRPhot",39),("FSRPhot",40)], #,("NucReact",41),("PiZero",42),("DiBoson",43),("ZorHeavyBoson",44),("QCD",45)],
                "normFactor"   : normFactor,
                "drawGrid"     : True,
               }

    extension = ""
    if "normFactor" in p0_props:
        if p0_props["normFactor"]:
            extension = "_NORM"
    doLeptonOriginFracPlots = ( not "normFactor" in p0_props or not p0_props["normFactor"] )

    makeStackOriginFractions = False

    for idx, path in enumerate(pathlist):

        p0 = Plot(kwargs["sample"], path, p0_props)

        canvasID = varlist[idx] + "_" + "Lep" + "_" + kwargs["prodID"] + "_" + regionlist[idx] + extension
        c = TCanvas(canvasID,canvasID,50,50,1000,700)
        c.cd()
        p0.makePlot()

        if doLeptonOriginFracPlots:
            if makeStackOriginFractions:
                # THIS IS TO MAKE STACK PLOTS
                #
                canvasID_stack = "FakeLepOriginFrac" + "_" + "Lep" + "_" + varlist[idx] + "_" + kwargs["prodID"] + "_" + regionlist[idx] + "_" + flaglist[idx][1]
                cstack = TCanvas(canvasID_stack,canvasID_stack,50,50,1000,700)
                cstack.cd()
                mystack, mystackleg = p0.makeLeptonOriginFracPlots(canvasID_stack)
                mystack.Draw()
                mystack.GetXaxis().SetTitle("p_{T}^{l}")
                mystack.GetYaxis().SetTitle("Fake {0} origin fraction".format(header))
                mystackleg.Draw()
                saveName_stack = canvasID_stack
                print("\nRegion: {0}\nPhoton conversion fraction VS. Pt: {1}\n".format(flaglist[idx][0],p0.conversion_frac_VS_Y))
            else:
                # THIS IS TO MAKE CONVERSION FRACTION ONLY PLOT, WITH UNCERTAINTIES
                #
                canvasIDfrac = "FakeLepOriginFrac" + "_" + "Lep" + "_" + varlist[idx] + "_" + kwargs["prodID"] + "_" + regionlist[idx] + "_" + flaglist[idx][1]
                cfrac = TCanvas(canvasIDfrac,canvasIDfrac,50,50,1000,700)
                cfrac.cd()
                print("\nRegion: {0}\n".format(flaglist[idx][0]))
                hist, leg = p0.makeConversionFracHist(canvasIDfrac)
                hist.GetYaxis().SetRangeUser(0,1)
                hist.SetLineColor(kOrange)
                hist.SetMarkerColor(kOrange)
                hist.GetXaxis().SetTitle("p_{T}(l)")
                hist.GetXaxis().SetTitleOffset(1.0)
                hist.GetYaxis().SetTitle("#gamma conversion fraction")
                hist.GetYaxis().SetTitleOffset(1.1)
                leg.Draw()
                saveName_frac = canvasIDfrac

        # Plot.legend.Draw()
        # Plot.legendATLAS.DrawLatex(0.6,0.35,"#bf{#it{ATLAS}} Work In Progress")
        # Plot.legendLumi.DrawLatex(0.6,0.27,"#sqrt{{s}} = 13 TeV, #int L dt = {0:.1f} fb^{{-1}}".format(Plot.luminosity))

        savePath = basePath + "/" + "Lep" + "_" + "FakeLepOriginFrac_VS_Pt_" + kwargs["prodID"] + "/"

        if not os.path.exists(savePath):
            os.makedirs(savePath)

        saveName = canvasID
        if makeStackOriginFractions:
            for ext in ["png","pdf","root"]:
                c.SaveAs( savePath + saveName + "_" + flaglist[idx][1] + "." + ext )
                if doLeptonOriginFracPlots:
                    cstack.SaveAs( savePath + saveName_stack + "_" + flaglist[idx][1] + "." + ext )
        else:
            hist_outputfile = TFile(savePath + saveName_frac + "_" + flaglist[idx][1] + "_HIST.root","RECREATE")
            hist.Write()


    Plot.legend.Clear()

def plotOriginVSEta( normFactor=0, **kwargs ):

    plotlist = []

    if kwargs["flavour"] == "Mu":
        header   = "muons"
    elif kwargs["flavour"] == "El":
        header   = "electrons"

    Plot.legend.SetHeader("Fake {0}".format(header))

    variable_2L = "Lep1Origin_VS_Lep1Eta"
    variable_3L = "Lep2Origin_VS_Lep2Eta"

    basePath = os.path.abspath(os.curdir) + "/" + "PLOTS_" + kwargs["prodID"] + "/" + "OutputPlots_FakeOriginFrac_TTBarTTBarGamma_" + kwargs["prodID"]

    inputPath_2LVR = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_LOWNJ_VR_TT" + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_"
    inputPath_2LSR = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_SR_TT" + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_SR_TT_"
    inputPath_3LSR = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_3L_SR_TT" + "/" + kwargs["flavour"] + "FakeOriginFrac_3L_SR_TT_"

    inputName_2L = variable_2L + ".root"
    inputName_3L = variable_3L + ".root"

    p0_props = {
                #"legend"      : "t#bar{t} - Loose sel.",
                "xAxisTitle"  : "truthOrigin",
                "yAxisTitle"  : "#eta^{lep}",
                "xAxisLabels" : [("NonDefined",0),("SingleElec",1),("SingleMuon",2),("SinglePhot",3),("SingleTau",4),("PhotonConv",5),("DalitzDec",6),("ElMagProc",7),("Mu",8),("TauLep",9),("top",10),("QuarkWeakDec",11),("WBoson",12),("ZBoson",13),("Higgs",14),("HiggsMSSM",15),("HeavyBoson",16),("WBosonLRSM",17),("NuREle",18),("NuRMu",19),("NuRTau",20),("LQ",21),("SUSY",22),("LightMeson",23),("StrangeMeson",24),("CharmedMeson",25),("BottomMeson",26),("CCbarMeson",27),("JPsi",28),("BBbarMeson",29),("LightBaryon",30),("StrangeBaryon",31),("CharmedBaryon",32),("BottomBaryon",33),("PionDecay",34),("KaonDecay",35),("BremPhot",36),("PromptPhot",37),("UndrPhot",38),("ISRPhot",39),("FSRPhot",40)], #,("NucReact",41),("PiZero",42),("DiBoson",43),("ZorHeavyBoson",44),("QCD",45)],
                "normFactor"   : normFactor,
                "drawGrid"     : True,
               }

    extension = ""
    if "normFactor" in p0_props:
        if p0_props["normFactor"]:
            extension = "_NORM"
    doLeptonOriginFracPlots = ( not "normFactor" in p0_props or not p0_props["normFactor"] )

    pathlist   = [inputPath_2LVR+inputName_2L,inputPath_2LSR+inputName_2L,inputPath_3LSR+inputName_3L]
    varlist    = [variable_2L,variable_2L,variable_3L]
    regionlist = ["2LSS_LOWNJ_VR_TT","2LSS_SR_TT","3L_SR_TT"]

    for idx, path in enumerate(pathlist):

        p0 = Plot(kwargs["sample"], path, p0_props)

        canvasID = varlist[idx] + "_" + kwargs["flavour"] + "_" + kwargs["prodID"] + "_" + regionlist[idx] + extension
        c = TCanvas(canvasID,canvasID,50,50,1000,700)
        c.cd()
        p0.makePlot()

        if doLeptonOriginFracPlots:
            canvasID_stack = "FakeLepOriginFrac" + "_" + kwargs["flavour"] + "_" + varlist[idx] + "_" + kwargs["prodID"] + "_" + regionlist[idx]
            cstack = TCanvas(canvasID_stack,canvasID_stack,50,50,1000,700)
            cstack.cd()
            mystack, mystackleg = p0.makeLeptonOriginFracPlots(canvasID_stack)
            mystack.Draw()
            mystack.GetXaxis().SetTitle("#eta^{lep}")
            mystack.GetYaxis().SetTitle("Fake {0} origin fraction".format(header))
            mystackleg.Draw()
            saveName_stack = canvasID_stack

        # Plot.legend.Draw()
        # Plot.legendATLAS.DrawLatex(0.6,0.35,"#bf{#it{ATLAS}} Work In Progress")
        # Plot.legendLumi.DrawLatex(0.6,0.27,"#sqrt{{s}} = 13 TeV, #int L dt = {0:.1f} fb^{{-1}}".format(Plot.luminosity))

        savePath = basePath + "/" + kwargs["flavour"] + "_" + "FakeLepOriginFrac_VS_Eta_" + kwargs["prodID"] + "/"

        if not os.path.exists(savePath):
            os.makedirs(savePath)

        saveName = canvasID
        for ext in ["png","pdf","root"]:
            c.SaveAs( savePath + saveName + "." + ext )
            if doLeptonOriginFracPlots:
                cstack.SaveAs( savePath + saveName_stack + "." + ext )

    Plot.legend.Clear()

def plotOriginVSEtaBE2( normFactor=0, **kwargs ):

    plotlist = []

    if kwargs["flavour"] == "Mu":
        header   = "muons"
    elif kwargs["flavour"] == "El":
        header   = "electrons"

    Plot.legend.SetHeader("Fake {0}".format(header))

    variable_2L = "Lep1Origin_VS_Lep1EtaBE2"
    variable_3L = "Lep2Origin_VS_Lep2EtaBE2"

    basePath = os.path.abspath(os.curdir) + "/" + "PLOTS_" + kwargs["prodID"] + "_v2/" + "OutputPlots_FakeOriginFrac_TTBarTTBarGamma_" + kwargs["prodID"]

    inputPath_2LVR = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_LOWNJ_VR_TT" + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_"
    inputPath_2LSR = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_SR_TT" + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_SR_TT_"
    inputPath_3LSR = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_3L_SR_TT" + "/" + kwargs["flavour"] + "FakeOriginFrac_3L_SR_TT_"

    inputName_2L = variable_2L + ".root"
    inputName_3L = variable_3L + ".root"

    p0_props = {
                #"legend"      : "t#bar{t} - Loose sel.",
                "xAxisTitle"  : "truthOrigin",
                "yAxisTitle"  : "#eta^{lep}",
                "xAxisLabels" : [("NonDefined",0),("SingleElec",1),("SingleMuon",2),("SinglePhot",3),("SingleTau",4),("PhotonConv",5),("DalitzDec",6),("ElMagProc",7),("Mu",8),("TauLep",9),("top",10),("QuarkWeakDec",11),("WBoson",12),("ZBoson",13),("Higgs",14),("HiggsMSSM",15),("HeavyBoson",16),("WBosonLRSM",17),("NuREle",18),("NuRMu",19),("NuRTau",20),("LQ",21),("SUSY",22),("LightMeson",23),("StrangeMeson",24),("CharmedMeson",25),("BottomMeson",26),("CCbarMeson",27),("JPsi",28),("BBbarMeson",29),("LightBaryon",30),("StrangeBaryon",31),("CharmedBaryon",32),("BottomBaryon",33),("PionDecay",34),("KaonDecay",35),("BremPhot",36),("PromptPhot",37),("UndrPhot",38),("ISRPhot",39),("FSRPhot",40)], #,("NucReact",41),("PiZero",42),("DiBoson",43),("ZorHeavyBoson",44),("QCD",45)],
                "normFactor"   : normFactor,
                "drawGrid"     : True,
               }

    extension = ""
    if "normFactor" in p0_props:
        if p0_props["normFactor"]:
            extension = "_NORM"
    doLeptonOriginFracPlots = ( not "normFactor" in p0_props or not p0_props["normFactor"] )

    pathlist   = [inputPath_2LVR+inputName_2L,inputPath_2LSR+inputName_2L,inputPath_3LSR+inputName_3L]
    varlist    = [variable_2L,variable_2L,variable_3L]
    regionlist = ["2LSS_LOWNJ_VR_TT","2LSS_SR_TT","3L_SR_TT"]

    for idx, path in enumerate(pathlist):

        p0 = Plot(kwargs["sample"], path, p0_props)

        canvasID = varlist[idx] + "_" + kwargs["flavour"] + "_" + kwargs["prodID"] + "_" + regionlist[idx] + extension
        c = TCanvas(canvasID,canvasID,50,50,1000,700)
        c.cd()
        p0.makePlot()

        if doLeptonOriginFracPlots:
            canvasID_stack = "FakeLepOriginFrac" + "_" + kwargs["flavour"] + "_" + varlist[idx] + "_" + kwargs["prodID"] + "_" + regionlist[idx]
            cstack = TCanvas(canvasID_stack,canvasID_stack,50,50,1000,700)
            cstack.cd()
            mystack, mystackleg = p0.makeLeptonOriginFracPlots(canvasID_stack)
            mystack.Draw()
            mystack.GetXaxis().SetTitle("#eta^{lep}")
            mystack.GetYaxis().SetTitle("Fake {0} origin fraction".format(header))
            mystackleg.Draw()
            saveName_stack = canvasID_stack

        # Plot.legend.Draw()
        # Plot.legendATLAS.DrawLatex(0.6,0.35,"#bf{#it{ATLAS}} Work In Progress")
        # Plot.legendLumi.DrawLatex(0.6,0.27,"#sqrt{{s}} = 13 TeV, #int L dt = {0:.1f} fb^{{-1}}".format(Plot.luminosity))

        savePath = basePath + "/" + kwargs["flavour"] + "_" + "FakeLepOriginFrac_VS_EtaBE2_" + kwargs["prodID"] + "/"

        if not os.path.exists(savePath):
            os.makedirs(savePath)

        saveName = canvasID
        for ext in ["png","pdf","root"]:
            c.SaveAs( savePath + saveName + "." + ext )
            if doLeptonOriginFracPlots:
                cstack.SaveAs( savePath + saveName_stack + "." + ext )

    Plot.legend.Clear()


def plotOriginVSDistanceLepClosestJet( normFactor=0, **kwargs ):

    plotlist = []

    if kwargs["flavour"] == "Mu":
        header   = "muons"
    elif kwargs["flavour"] == "El":
        header   = "electrons"

    Plot.legend.SetHeader("Fake {0}".format(header))

    variable_2L = "Lep1Origin_VS_Lep1DistanceClosestJet"
    variable_3L = "Lep2Origin_VS_Lep2DistanceClosestJet"

    basePath = os.path.abspath(os.curdir) + "/" + "PLOTS_" + kwargs["prodID"] + "_v2/" + "OutputPlots_FakeOriginFrac_TTBarTTBarGamma_" + kwargs["prodID"]

    inputPath_2LVR = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_LOWNJ_VR_TT" + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_"
    inputPath_2LSR = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_SR_TT" + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_SR_TT_"
    inputPath_3LSR = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_3L_SR_TT" + "/" + kwargs["flavour"] + "FakeOriginFrac_3L_SR_TT_"

    inputName_2L = variable_2L + ".root"
    inputName_3L = variable_3L + ".root"

    p0_props = {
                #"legend"      : "t#bar{t} - Loose sel.",
                "xAxisTitle"  : "truthOrigin",
                "yAxisTitle"  : "#DeltaR(lep, closest jet)",
                "xAxisLabels" : [("NonDefined",0),("SingleElec",1),("SingleMuon",2),("SinglePhot",3),("SingleTau",4),("PhotonConv",5),("DalitzDec",6),("ElMagProc",7),("Mu",8),("TauLep",9),("top",10),("QuarkWeakDec",11),("WBoson",12),("ZBoson",13),("Higgs",14),("HiggsMSSM",15),("HeavyBoson",16),("WBosonLRSM",17),("NuREle",18),("NuRMu",19),("NuRTau",20),("LQ",21),("SUSY",22),("LightMeson",23),("StrangeMeson",24),("CharmedMeson",25),("BottomMeson",26),("CCbarMeson",27),("JPsi",28),("BBbarMeson",29),("LightBaryon",30),("StrangeBaryon",31),("CharmedBaryon",32),("BottomBaryon",33),("PionDecay",34),("KaonDecay",35),("BremPhot",36),("PromptPhot",37),("UndrPhot",38),("ISRPhot",39),("FSRPhot",40)], #,("NucReact",41),("PiZero",42),("DiBoson",43),("ZorHeavyBoson",44),("QCD",45)],
                "normFactor"   : normFactor,
                "drawGrid"     : True,
               }

    extension = ""
    if "normFactor" in p0_props:
        if p0_props["normFactor"]:
            extension = "_NORM"
    doLeptonOriginFracPlots = ( not "normFactor" in p0_props or not p0_props["normFactor"] )

    pathlist   = [inputPath_2LVR+inputName_2L,inputPath_2LSR+inputName_2L,inputPath_3LSR+inputName_3L]
    varlist    = [variable_2L,variable_2L,variable_3L]
    regionlist = ["2LSS_LOWNJ_VR_TT","2LSS_SR_TT","3L_SR_TT"]

    for idx, path in enumerate(pathlist):

        p0 = Plot(kwargs["sample"], path, p0_props)

        canvasID = varlist[idx] + "_" + kwargs["flavour"] + "_" + kwargs["prodID"] + "_" + regionlist[idx] + extension
        c = TCanvas(canvasID,canvasID,50,50,1000,700)
        c.cd()
        p0.makePlot()

        if doLeptonOriginFracPlots:
            canvasID_stack = "FakeLepOriginFrac" + "_" + kwargs["flavour"] + "_" + varlist[idx] + "_" + kwargs["prodID"] + "_" + regionlist[idx]
            cstack = TCanvas(canvasID_stack,canvasID_stack,50,50,1000,700)
            cstack.cd()
            mystack, mystackleg = p0.makeLeptonOriginFracPlots(canvasID_stack)
            mystack.Draw()
            mystack.GetXaxis().SetTitle("#DeltaR(lep, closest jet)")
            mystack.GetYaxis().SetTitle("Fake {0} origin fraction".format(header))
            mystackleg.Draw()
            saveName_stack = canvasID_stack

        # Plot.legend.Draw()
        # Plot.legendATLAS.DrawLatex(0.6,0.35,"#bf{#it{ATLAS}} Work In Progress")
        # Plot.legendLumi.DrawLatex(0.6,0.27,"#sqrt{{s}} = 13 TeV, #int L dt = {0:.1f} fb^{{-1}}".format(Plot.luminosity))

        savePath = basePath + "/" + kwargs["flavour"] + "_" + "FakeLepOriginFrac_VS_DistanceLepClosestJet_" + kwargs["prodID"] + "/"

        if not os.path.exists(savePath):
            os.makedirs(savePath)

        saveName = canvasID
        for ext in ["png","pdf","root"]:
            c.SaveAs( savePath + saveName + "." + ext )
            if doLeptonOriginFracPlots:
                cstack.SaveAs( savePath + saveName_stack + "." + ext )

    Plot.legend.Clear()


def plotOriginVSDistanceOtherLep( normFactor=0, **kwargs ):

    plotlist = []

    if kwargs["flavour"] == "Mu":
        header   = "muons"
    elif kwargs["flavour"] == "El":
        header   = "electrons"

    Plot.legend.SetHeader("Fake {0}".format(header))

    variable_2L = "Lep1Origin_VS_DistanceOtherLep"
    variable_3L = "Lep2Origin_VS_DistanceOtherLep"

    basePath = os.path.abspath(os.curdir) + "/" + "PLOTS_" + kwargs["prodID"] + "_v2/" + "OutputPlots_FakeOriginFrac_TTBarTTBarGamma_" + kwargs["prodID"]

    pathlist = []

    inputName_2L = variable_2L + ".root"
    inputName_3L = variable_3L + ".root"

    inputPath_2LVR = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_LOWNJ_VR_TT" + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_"
    pathlist.append(inputPath_2LVR+inputName_2L)
    inputPath_2LSR = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_SR_TT" + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_SR_TT_"
    pathlist.append(inputPath_2LSR+inputName_2L)
    inputPath_3LSR = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_3L_SR_TT" + "/" + kwargs["flavour"] + "FakeOriginFrac_3L_SR_TT_"
    pathlist.append(inputPath_3LSR+inputName_3L)
    if kwargs["flavour"] == "El":
        inputPath_2LVR_ElEl = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_ElEl" + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_ElEl_"
        pathlist.append(inputPath_2LVR_ElEl+inputName_2L)
        inputPath_2LVR_OF   = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_OF" + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_OF_"
        pathlist.append(inputPath_2LVR_OF+inputName_2L)
        inputPath_2LSR_ElEl = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_SR_TT_ElEl" + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_SR_TT_ElEl_"
        pathlist.append(inputPath_2LSR_ElEl+inputName_2L)
        inputPath_2LSR_OF   = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_SR_TT_OF" + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_SR_TT_OF_"
        pathlist.append(inputPath_2LSR_OF+inputName_2L)
        # Split 1bjet, >= 2 bjet
        inputPath_2LVR_1BJET_ElEl = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_1BJET_LOWNJ_VR_TT_ElEl" + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_1BJET_LOWNJ_VR_TT_ElEl_"
        pathlist.append(inputPath_2LVR_1BJET_ElEl+inputName_2L)
        inputPath_2LVR_1BJET_OF = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_1BJET_LOWNJ_VR_TT_OF" + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_1BJET_LOWNJ_VR_TT_OF_"
        pathlist.append(inputPath_2LVR_1BJET_OF+inputName_2L)
        inputPath_2LVR_2BJET_ElEl = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_2BJET_LOWNJ_VR_TT_ElEl" + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_2BJET_LOWNJ_VR_TT_ElEl_"
        pathlist.append(inputPath_2LVR_2BJET_ElEl+inputName_2L)
        inputPath_2LVR_2BJET_OF = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_2BJET_LOWNJ_VR_TT_OF" + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_2BJET_LOWNJ_VR_TT_OF_"
        pathlist.append(inputPath_2LVR_2BJET_OF+inputName_2L)
        inputPath_2LSR_1BJET_ElEl = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_1BJET_SR_TT_ElEl" + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_1BJET_SR_TT_ElEl_"
        pathlist.append(inputPath_2LSR_1BJET_ElEl+inputName_2L)
        inputPath_2LSR_1BJET_OF = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_1BJET_SR_TT_OF" + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_1BJET_SR_TT_OF_"
        pathlist.append(inputPath_2LSR_1BJET_OF+inputName_2L)
        inputPath_2LSR_2BJET_ElEl = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_2BJET_SR_TT_ElEl" + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_2BJET_SR_TT_ElEl_"
        pathlist.append(inputPath_2LSR_2BJET_ElEl+inputName_2L)
        inputPath_2LSR_2BJET_OF = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_2BJET_SR_TT_OF" + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_2BJET_SR_TT_OF_"
        pathlist.append(inputPath_2LSR_2BJET_OF+inputName_2L)

    p0_props = {
                #"legend"      : "t#bar{t} - Loose sel.",
                "xAxisTitle"  : "truthOrigin",
                "yAxisTitle"  : "#DeltaR(l_{0},l_{1})",
                "xAxisLabels" : [("NonDefined",0),("SingleElec",1),("SingleMuon",2),("SinglePhot",3),("SingleTau",4),("PhotonConv",5),("DalitzDec",6),("ElMagProc",7),("Mu",8),("TauLep",9),("top",10),("QuarkWeakDec",11),("WBoson",12),("ZBoson",13),("Higgs",14),("HiggsMSSM",15),("HeavyBoson",16),("WBosonLRSM",17),("NuREle",18),("NuRMu",19),("NuRTau",20),("LQ",21),("SUSY",22),("LightMeson",23),("StrangeMeson",24),("CharmedMeson",25),("BottomMeson",26),("CCbarMeson",27),("JPsi",28),("BBbarMeson",29),("LightBaryon",30),("StrangeBaryon",31),("CharmedBaryon",32),("BottomBaryon",33),("PionDecay",34),("KaonDecay",35),("BremPhot",36),("PromptPhot",37),("UndrPhot",38),("ISRPhot",39),("FSRPhot",40)], #,("NucReact",41),("PiZero",42),("DiBoson",43),("ZorHeavyBoson",44),("QCD",45)],
                "normFactor"   : normFactor,
                "drawGrid"     : True,
               }

    extension = ""
    if "normFactor" in p0_props:
        if p0_props["normFactor"]:
            extension = "_NORM"
    doLeptonOriginFracPlots = ( not "normFactor" in p0_props or not p0_props["normFactor"] )

    varlist    = [variable_2L,variable_2L,variable_3L]
    regionlist = ["2LSS_LOWNJ_VR_TT","2LSS_SR_TT","3L_SR_TT"]
    if kwargs["flavour"] == "El":
        regionlist.extend(["2LSS_LOWNJ_VR_TT_ElEl","2LSS_LOWNJ_VR_TT_OF","2LSS_SR_TT_ElEl","2LSS_SR_TT_OF"])
        varlist.extend([variable_2L,variable_2L,variable_2L,variable_2L])
        # Split 1bjet, >= 2 bjet
        regionlist.extend(["2LSS_1BJET_LOWNJ_VR_TT_ElEl","2LSS_1BJET_LOWNJ_VR_TT_OF","2LSS_2BJET_LOWNJ_VR_TT_ElEl","2LSS_2BJET_LOWNJ_VR_TT_OF","2LSS_1BJET_SR_TT_ElEl","2LSS_1BJET_SR_TT_OF","2LSS_2BJET_SR_TT_ElEl","2LSS_2BJET_SR_TT_OF"])
        varlist.extend([variable_2L,variable_2L,variable_2L,variable_2L,variable_2L,variable_2L,variable_2L,variable_2L])

    for idx, path in enumerate(pathlist):

        print("Reading input from : {0}".format(path))

        p0 = Plot(kwargs["sample"], path, p0_props)

        canvasID = varlist[idx] + "_" + kwargs["flavour"] + "_" + kwargs["prodID"] + "_" + regionlist[idx] + extension
        c = TCanvas(canvasID,canvasID,50,50,1000,700)
        c.cd()
        p0.makePlot()

        if doLeptonOriginFracPlots:
            canvasID_stack = "FakeLepOriginFrac" + "_" + kwargs["flavour"] + "_" + varlist[idx] + "_" + kwargs["prodID"] + "_" + regionlist[idx]
            cstack = TCanvas(canvasID_stack,canvasID_stack,50,50,1000,700)
            cstack.cd()
            mystack, mystackleg = p0.makeLeptonOriginFracPlots(canvasID_stack)
            mystack.Draw()
            mystack.GetXaxis().SetTitle("#DeltaR(l_{0},l_{1})")
            mystack.GetYaxis().SetTitle("Fake {0} origin fraction".format(header))
            mystackleg.Draw()
            saveName_stack = canvasID_stack

        # Plot.legend.Draw()
        # Plot.legendATLAS.DrawLatex(0.6,0.35,"#bf{#it{ATLAS}} Work In Progress")
        # Plot.legendLumi.DrawLatex(0.6,0.27,"#sqrt{{s}} = 13 TeV, #int L dt = {0:.1f} fb^{{-1}}".format(Plot.luminosity))

        savePath = basePath + "/" + kwargs["flavour"] + "_" + "FakeLepOriginFrac_VS_DistanceLepOtherLep_" + kwargs["prodID"] + "/"

        if not os.path.exists(savePath):
            os.makedirs(savePath)

        saveName = canvasID
        for ext in ["png","pdf","root"]:
            c.SaveAs( savePath + saveName + "." + ext )
            if doLeptonOriginFracPlots:
                cstack.SaveAs( savePath + saveName_stack + "." + ext )

    Plot.legend.Clear()


def plotOriginVSBDTG( normFactor=0, **kwargs ):

    plotlist = []

    if kwargs["flavour"] == "Mu":
        header   = "muons"
    elif kwargs["flavour"] == "El":
        header   = "electrons"

    Plot.legend.SetHeader("Fake {0}".format(header))

    variable_2L = "El1Origin_VS_BDTGScore"

    # flav = "OF"
    flav = "ElEl"

    basePath = os.path.abspath(os.curdir) + "/" + "PLOTS_" + kwargs["prodID"] + "/" + "OutputPlots_MM_MMClosureTest_HIGHNJ_" + kwargs["prodID"] + "_25_07_17_MCClosure_2L3L_RESCALED"

    pathlist = []

    inputName_2L = variable_2L + ".root"

    inputPath_2L = basePath + "/" + flav + "SS_SR_HighJet_DataDriven_Closure" + "/" + flav + "SS_SR_HighJet_DataDriven_Closure_"
    pathlist.append(inputPath_2L+inputName_2L)

    p0_props = {
                #"legend"      : "t#bar{t} - Loose sel.",
                "xAxisTitle"  : "truthOrigin",
                "yAxisTitle"  : "BDTG",
                "xAxisLabels" : [("NonDefined",0),("SingleElec",1),("SingleMuon",2),("SinglePhot",3),("SingleTau",4),("PhotonConv",5),("DalitzDec",6),("ElMagProc",7),("Mu",8),("TauLep",9),("top",10),("QuarkWeakDec",11),("WBoson",12),("ZBoson",13),("Higgs",14),("HiggsMSSM",15),("HeavyBoson",16),("WBosonLRSM",17),("NuREle",18),("NuRMu",19),("NuRTau",20),("LQ",21),("SUSY",22),("LightMeson",23),("StrangeMeson",24),("CharmedMeson",25),("BottomMeson",26),("CCbarMeson",27),("JPsi",28),("BBbarMeson",29),("LightBaryon",30),("StrangeBaryon",31),("CharmedBaryon",32),("BottomBaryon",33),("PionDecay",34),("KaonDecay",35),("BremPhot",36),("PromptPhot",37),("UndrPhot",38),("ISRPhot",39),("FSRPhot",40)], #,("NucReact",41),("PiZero",42),("DiBoson",43),("ZorHeavyBoson",44),("QCD",45)],
                "normFactor"   : normFactor,
                "drawGrid"     : True,
               }

    extension = ""
    if "normFactor" in p0_props:
        if p0_props["normFactor"]:
            extension = "_NORM"
    doLeptonOriginFracPlots = ( not "normFactor" in p0_props or not p0_props["normFactor"] )

    varlist    = [variable_2L]
    regionlist = ["2LSS_HIGHNJ"]

    for idx, path in enumerate(pathlist):

        p0 = Plot(kwargs["sample"], path, p0_props)

        canvasID = varlist[idx] + "_" + flav + "_" + kwargs["prodID"] + "_" + regionlist[idx] + extension
        c = TCanvas(canvasID,canvasID,50,50,1000,700)
        c.cd()
        #p0.makePlot()

        if doLeptonOriginFracPlots:
            # THIS IS TO MAKE STACK PLOTS
            #
            # canvasID_stack = "FakeLepOriginFrac" + "_" + flav + "_" + varlist[idx] + "_" + kwargs["prodID"] + "_" + regionlist[idx]
            # cstack = TCanvas(canvasID_stack,canvasID_stack,50,50,1000,700)
            # cstack.cd()
            # mystack, mystackleg = p0.makeLeptonOriginFracPlots(canvasID_stack)
            # mystack.Draw()
            # mystack.GetXaxis().SetTitle("BDTG")
            # mystack.GetYaxis().SetTitle("Fake {0} origin fraction".format(header))
            # mystackleg.Draw()
            # saveName_stack = canvasID_stack
            # print("\nRegion: {0}\nPhoton conversion fraction VS. BDTG: {1}\n".format(regionlist[idx],p0.conversion_frac_VS_Y))
            #
            # THIS IS TO MAKE CONVERSION FRACTION ONLY PLOT, WITH UNCERTAINTIES
            #
            canvasIDfrac = "FakeLepOriginFrac" + "_" + flav + "_" + varlist[idx] + "_" + kwargs["prodID"] + "_" + regionlist[idx]
            cfrac = TCanvas(canvasIDfrac,canvasIDfrac,50,50,1000,700)
            cfrac.cd()
            hist, leg = p0.makeConversionFracHist(canvasIDfrac)
            hist.GetXaxis().SetTitle("BDTG")
            hist.GetYaxis().SetTitle("Fake {0} origin fraction".format(header))
            leg.Draw()
            saveName_frac = canvasIDfrac

        # Plot.legend.Draw()
        # Plot.legendATLAS.DrawLatex(0.6,0.35,"#bf{#it{ATLAS}} Work In Progress")
        # Plot.legendLumi.DrawLatex(0.6,0.27,"#sqrt{{s}} = 13 TeV, #int L dt = {0:.1f} fb^{{-1}}".format(Plot.luminosity))

        savePath = basePath + "/" + flav + "_" + "FakeLepOriginFrac_VS_BDTG_" + kwargs["prodID"] + "/"

        if not os.path.exists(savePath):
            os.makedirs(savePath)

        saveName = canvasID
        for ext in ["png","pdf","root"]:
            #c.SaveAs( savePath + saveName + "." + ext )
            if doLeptonOriginFracPlots:
                #cstack.SaveAs( savePath + saveName_stack + "." + ext )
                cfrac.SaveAs( savePath + saveName_frac + "." + ext )
        hist_outputfile = TFile(savePath + saveName_frac + "_HIST.root","RECREATE")
        hist.Write()


    Plot.legend.Clear()


def plotFakeOriginFrac2L3L( **kwargs ):

    plotlist = []

    if kwargs["flavour"] == "Mu":
        header   = "muon"
    elif kwargs["flavour"] == "El":
        header   = "electron"

    Plot.legend.SetHeader("Fake {0}".format(header))

    variable_2L = "Lep1Origin"
    variable_3L = "Lep2Origin"

    basePath = os.path.abspath(os.curdir) + "/" + "PLOTS_" + kwargs["prodID"] + "OutputPlots_FakeOriginFrac_TTBarTTBarGamma_" + kwargs["prodID"] + "_SplitFlavours"

    inputPath_2LVR = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_LOWNJ_VR_TT" + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_"
    inputPath_2LSR = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_SR_TT" + "/" + kwargs["flavour"] + "FakeOriginFrac_2LSS_SR_TT_"
    inputPath_3LSR = basePath + "/" + kwargs["flavour"] + "FakeOriginFrac_3L_SR_TT" + "/" + kwargs["flavour"] + "FakeOriginFrac_3L_SR_TT_"

    inputName_2L = variable_2L + ".root"
    inputName_3L = variable_3L + ".root"

    f_2LVR = TFile(inputPath_2LVR+inputName_2L)
    h_2LVR = f_2LVR.Get(kwargs["sample"])
    h_2LVR.SetDirectory(0)

    f_2LSR = TFile(inputPath_2LSR+inputName_2L)
    h_2LSR = f_2LSR.Get(kwargs["sample"])
    h_2LSR.SetDirectory(0)

    f_3LSR = TFile(inputPath_3LSR+inputName_3L)
    h_3LSR = f_3LSR.Get(kwargs["sample"])
    h_3LSR.SetDirectory(0)

    h_list = [ h_2LVR, h_2LSR, h_3LSR ]

    # Fake lepton origin fraction wrt. 2LVR, 2LSR, 3LSR

    histfakes_BF        = TH1D("histfakes_BF","histfakes_BF",3,-0.5,2.5) # B-hadrons in jets (mesons/baryons)
    histfakes_CF        = TH1D("histfakes_CF","histfakes_CF",3,-0.5,2.5) # C-hadrons in jets (mesons/baryons)
    histfakes_HFRes     = TH1D("histfakes_HFRes","histfakes_HFRes",3,-0.5,2.5) # B,C resonances (J/psi, Upsilon...)
    histfakes_LF        = TH1D("histfakes_LF","histfakes_LF",3,-0.5,2.5) # Light hadrons in jets
    histfakes_PhConv    = TH1D("histfakes_PhConv","histfakes_PhConv",3,-0.5,2.5) # Photon conversions
    histfakes_Other     = TH1D("histfakes_Other","histfakes_Other",3,-0.5,2.5) # Other fakes (mis-id jets, leptons from generic pi/K...)
    histfakes_Unknown   = TH1D("histfakes_Unknown","histfakes_Unknown",3,-0.5,2.5) # Unknown fakes (failure of MCTruthClassifier

    stacklegend = TLegend(0.23,0.25,0.43,0.55)
    stacklegend.SetBorderSize(1)
    stacklegend.SetFillColor(kWhite)
    stacklegend.SetTextSize(0.03)
    stacklegend.SetTextFont(42)

    histfakes_list = [ (histfakes_BF,kRed), (histfakes_CF,kRed-9), (histfakes_HFRes,kPink-2), (histfakes_LF,kOrange+1), (histfakes_PhConv,kYellow), (histfakes_Other,kPink+1), (histfakes_Unknown,kAzure+1) ]

    for h in histfakes_list:
        h[0].SetLineWidth(2)
        h[0].SetLineStyle(1)
        h[0].SetLineColor(1)
        h[0].SetFillColor(h[1])

    # p0_props = {
    #             #"legend"      : " t#bar{t} + t#bar{t}#gamma",
    #             "xAxisTitle"  : "truthOrigin",
    #             "xAxisLabels" : [("NonDefined",0),("SingleElec",1),("SingleMuon",2),("SinglePhot",3),("SingleTau",4),("PhotonConv",5),("DalitzDec",6),("ElMagProc",7),("Mu",8),("TauLep",9),("top",10),("QuarkWeakDec",11),("WBoson",12),("ZBoson",13),("Higgs",14),("HiggsMSSM",15),("HeavyBoson",16),("WBosonLRSM",17),("NuREle",18),("NuRMu",19),("NuRTau",20),("LQ",21),("SUSY",22),("LightMeson",23),("StrangeMeson",24),("CharmedMeson",25),("BottomMeson",26),("CCbarMeson",27),("JPsi",28),("BBbarMeson",29),("LightBaryon",30),("StrangeBaryon",31),("CharmedBaryon",32),("BottomBaryon",33),("PionDecay",34),("KaonDecay",35),("BremPhot",36),("PromptPhot",37),("UndrPhot",38),("ISRPhot",39),("FSRPhot",40)], #,("NucReact",41),("PiZero",42),("DiBoson",43),("ZorHeavyBoson",44),("QCD",45)],
    #             "normFactor"   : normFactor,
    #            }

    # Loop over truthOrigin hists: 2LVR, 2LSR, 3LSR

    for idx, h in enumerate(h_list):

        offset = 1 # (to account for underflow bin, which has idx=0)

        # Get the tot. fakes for *this* hist

        fakes_TOT_h = h.Integral( 0, h.GetNbinsX()+1 )

        # Get the HF fakes for *this* hist

        fakes_BF_h    = h.Integral( 26+offset,26+offset ) + h.Integral(33+offset,33+offset)
        fakes_CF_h    = h.Integral( 25+offset,25+offset ) + h.Integral(32+offset,32+offset)
        fakes_HFRes_h = h.Integral( 27+offset,29+offset )

        # Get the LF fakes for *this* hist

        fakes_LF_h = h.Integral( 23+offset,24+offset ) + h.Integral( 30+offset,31+offset )

        # Get the photon conversion fakes for *this* hist

        fakes_PhConv_h = h.Integral( 5+offset,5+offset )

        # Get the "Unknown" fakes for *this* Y

        fakes_Unknown_h = h.Integral( 0+offset,0+offset )

        # Get the other fakes for *this* hist

        fakes_Other_h = fakes_TOT_h - ( fakes_BF_h + fakes_CF_h + fakes_HFRes_h + fakes_LF_h + fakes_PhConv_h + fakes_Unknown_h )

        # Set the bin content for the fake lepton origin fraction hists for *this* hist

        if fakes_TOT_h:
            fakes_BF_frac_h      = fakes_BF_h/fakes_TOT_h
            fakes_CF_frac_h      = fakes_CF_h/fakes_TOT_h
            fakes_HFRes_frac_h   = fakes_HFRes_h/fakes_TOT_h
            fakes_LF_frac_h      = fakes_LF_h/fakes_TOT_h
            fakes_PhConv_frac_h  = fakes_PhConv_h/fakes_TOT_h
            fakes_Unknown_frac_h = fakes_Unknown_h/fakes_TOT_h
            fakes_Other_frac_h   = fakes_Other_h/fakes_TOT_h
        else:
            fakes_BF_frac_h = fakes_CF_frac_h = fakes_HFRes_frac_h = fakes_LF_frac_h = fakes_PhConv_frac_h = fakes_Unknown_frac_h = fakes_Other_frac_h = 0

        print("bin[{0}] ".format(idx))
        print("\ttot fakes = {0}".format(fakes_TOT_h))
        print("\t-) BF fakes = {0} ({1:.2f})".format(fakes_BF_h,fakes_BF_frac_h))
        print("\t-) CF fakes = {0} ({1:.2f})".format(fakes_CF_h,fakes_CF_frac_h))
        print("\t-) HFRes fakes = {0} ({1:.2f})".format(fakes_HFRes_h,fakes_HFRes_frac_h))
        print("\t-) LF fakes = {0} ({1:.2f})".format(fakes_LF_h,fakes_LF_frac_h))
        print("\t-) PhConv fakes = {0} ({1:.2f})".format(fakes_PhConv_h,fakes_PhConv_frac_h))
        print("\t-) Unknown fakes = {0} ({1:.2f})".format(fakes_Unknown_h,fakes_Unknown_frac_h))
        print("\t-) Other fakes = {0} ({1:.2f})".format(fakes_Other_h,fakes_Other_frac_h))

        histfakes_BF.SetBinContent( idx+offset, fakes_BF_frac_h )
        histfakes_CF.SetBinContent( idx+offset, fakes_CF_frac_h )
        histfakes_HFRes.SetBinContent( idx+offset, fakes_HFRes_frac_h )
        histfakes_LF.SetBinContent( idx+offset, fakes_LF_frac_h )
        histfakes_PhConv.SetBinContent( idx+offset, fakes_PhConv_frac_h )
        histfakes_Unknown.SetBinContent( idx+offset, fakes_Unknown_frac_h )
        histfakes_Other.SetBinContent( idx+offset, fakes_Other_frac_h )

    # Add histograms w/ fake origin fractions into a stack plot

    stacklegend.AddEntry(histfakes_BF, "B-Had Fakes", "F")
    stacklegend.AddEntry(histfakes_CF, "C-Had Fakes", "F")
    stacklegend.AddEntry(histfakes_HFRes, "J/#psi,#Upsilon Fakes", "F")
    stacklegend.AddEntry(histfakes_LF, "L-Had Fakes", "F")
    stacklegend.AddEntry(histfakes_PhConv, "#gamma conversion" , "F")
    stacklegend.AddEntry(histfakes_Other, "Other Fakes", "F")
    stacklegend.AddEntry(histfakes_Unknown, "Unknown" , "F")

    stack = THStack("FakeLepOriginFrac_VS_2L3LCat_STACK","FakeLepOriginFrac_VS_2L3LCat_STACK")
    stack.Add(histfakes_BF)
    stack.Add(histfakes_CF)
    stack.Add(histfakes_HFRes)
    stack.Add(histfakes_LF)
    stack.Add(histfakes_PhConv)
    stack.Add(histfakes_Other)
    stack.Add(histfakes_Unknown)

    canvasID_stack = kwargs["flavour"] + "_" + "FakeLepOriginFrac_VS_2L3LCat_" + kwargs["sample"] + "_" + kwargs["prodID"]
    cstack = TCanvas(canvasID_stack,canvasID_stack,50,50,1000,700)

    cstack.cd()
    stack.Draw()
    stack.GetXaxis().SetTitle("Category")
    stack.GetYaxis().SetTitle("Fake {0} origin fraction".format(header))
    stack.GetXaxis().SetRangeUser(-0.5,stack.GetHistogram().GetNbinsX()-0.5)

    stack.GetXaxis().SetBinLabel(1,"2lSS, CR")
    stack.GetXaxis().SetBinLabel(2,"2lSS, SR")
    stack.GetXaxis().SetBinLabel(3,"3l, SR")
    stack.GetYaxis().SetLabelSize(stack.GetXaxis().GetLabelSize())

    stacklegend.Draw()
    saveName_stack = canvasID_stack

    savePath = basePath + "/" + kwargs["flavour"] + "_" + "FakeLepOriginFrac_VS_2L3LCat_" + kwargs["prodID"] + "/"

    if not os.path.exists(savePath):
        os.makedirs(savePath)

    for ext in ["png","pdf","root"]:
        cstack.SaveAs( savePath + saveName_stack + "." + ext )


 # def plotFakeOriginFrac2LSplitFlavours( **kwargs ):

 #    plotlist = []

 #    Plot.legend.SetHeader("Fake {0}".format("leptons"))

 #    variable    = "Lep1Origin"
 #    variable_em = "Lep0Origin"

 #    basePath = os.path.abspath(os.curdir) + "/" + "PLOTS_" + kwargs["prodID"] + "/" + "OutputPlots_FakeOriginFrac_TTBarTTBarGamma_" + kwargs["prodID"] + "_SplitFlavours"

 #    inputPath_ee_VR = basePath + "/" + "El" + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_ElEl" + "/" + "El" + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_ElEl_"
 #    inputPath_ee_SR = basePath + "/" + "El" + "FakeOriginFrac_2LSS_SR_TT_ElEl" + "/" + "El" + "FakeOriginFrac_2LSS_SR_TT_ElEl_"
 #    inputPath_em_VR = basePath + "/" + "El" + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_ElMu" + "/" + "El" + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_ElMu_"
 #    inputPath_em_SR = basePath + "/" + "El" + "FakeOriginFrac_2LSS_SR_TT_ElMu" + "/" + "El" + "FakeOriginFrac_2LSS_SR_TT_ElMu_"
 #    inputPath_me_VR = basePath + "/" + "El" + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_MuEl" + "/" + "El" + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_MuEl_"
 #    inputPath_me_SR = basePath + "/" + "El" + "FakeOriginFrac_2LSS_SR_TT_MuEl" + "/" + "El" + "FakeOriginFrac_2LSS_SR_TT_MuEl_"
 #    inputPath_mm_VR = basePath + "/" + "Mu" + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_MuMu" + "/" + "Mu" + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_MuMu_"
 #    inputPath_mm_SR = basePath + "/" + "Mu" + "FakeOriginFrac_2LSS_SR_TT_MuMu" + "/" + "Mu" + "FakeOriginFrac_2LSS_SR_TT_MuMu_"

 #    inputName = variable + ".root"
 #    inputName_em = variable_em + ".root"

 #    f_ee_VR = TFile(inputPath_ee_VR+inputName)
 #    h_ee_VR = f_ee_VR.Get(kwargs["sample"])
 #    h_ee_VR.SetDirectory(0)
 #    f_ee_SR = TFile(inputPath_ee_SR+inputName)
 #    h_ee_SR = f_ee_SR.Get(kwargs["sample"])
 #    h_ee_SR.SetDirectory(0)

 #    f_em_VR = TFile(inputPath_em_VR+inputName_em)
 #    h_em_VR = f_em_VR.Get(kwargs["sample"])
 #    h_em_VR.SetDirectory(0)
 #    f_em_SR = TFile(inputPath_em_SR+inputName_em)
 #    h_em_SR = f_em_SR.Get(kwargs["sample"])
 #    h_em_SR.SetDirectory(0)

 #    f_me_VR = TFile(inputPath_me_VR+inputName)
 #    h_me_VR = f_me_VR.Get(kwargs["sample"])
 #    h_me_VR.SetDirectory(0)
 #    f_me_SR = TFile(inputPath_me_SR+inputName)
 #    h_me_SR = f_me_SR.Get(kwargs["sample"])
 #    h_me_SR.SetDirectory(0)

 #    f_mm_VR = TFile(inputPath_mm_VR+inputName)
 #    h_mm_VR = f_mm_VR.Get(kwargs["sample"])
 #    h_mm_VR.SetDirectory(0)
 #    f_mm_SR = TFile(inputPath_mm_SR+inputName)
 #    h_mm_SR = f_mm_SR.Get(kwargs["sample"])
 #    h_mm_SR.SetDirectory(0)

 #    h_list = [ h_ee_VR, h_ee_SR, h_em_VR, h_em_SR, h_me_VR, h_me_SR, h_mm_VR, h_mm_SR ]

 #    # Fake lepton origin fraction wrt. 2LVR, 2LSR, 3LSR

 #    histfakes_BF        = TH1D("histfakes_BF","histfakes_BF",8,-0.5,7.5) # B-hadrons in jets (mesons/baryons)
 #    histfakes_CF        = TH1D("histfakes_CF","histfakes_CF",8,-0.5,7.5) # C-hadrons in jets (mesons/baryons)
 #    histfakes_HFRes     = TH1D("histfakes_HFRes","histfakes_HFRes",8,-0.5,7.5) # B,C resonances (J/psi, Upsilon...)
 #    histfakes_LF        = TH1D("histfakes_LF","histfakes_LF",8,-0.5,7.5) # Light hadrons in jets
 #    histfakes_PhConv    = TH1D("histfakes_PhConv","histfakes_PhConv",8,-0.5,7.5) # Photon conversions
 #    histfakes_Other     = TH1D("histfakes_Other","histfakes_Other",8,-0.5,7.5) # Other fakes (mis-id jets, leptons from generic pi/K...)
 #    histfakes_Unknown   = TH1D("histfakes_Unknown","histfakes_Unknown",8,-0.5,7.5) # Unknown fakes (failure of MCTruthClassifier)

 #    stacklegend = TLegend(0.80,0.25,0.95,0.55)
 #    stacklegend.SetBorderSize(1)
 #    stacklegend.SetFillColor(kWhite)
 #    stacklegend.SetTextSize(0.03)
 #    stacklegend.SetTextFont(42)

 #    histfakes_list = [ (histfakes_BF,kRed), (histfakes_CF,kRed-9), (histfakes_HFRes,kPink-2), (histfakes_LF,kOrange+1), (histfakes_PhConv,kYellow), (histfakes_Other,kPink+1), (histfakes_Unknown,kAzure+1) ]

 #    for h in histfakes_list:
 #        h[0].SetLineWidth(2)
 #        h[0].SetLineStyle(1)
 #        h[0].SetLineColor(1)
 #        h[0].SetFillColor(h[1])

 #    for idx, h in enumerate(h_list):

 #        offset = 1 # (to account for underflow bin, which has idx=0)

 #        # Get the tot. fakes for *this* hist

 #        fakes_TOT_h = h.Integral( 0, h.GetNbinsX()+1 )

 #        # Get the HF fakes for *this* hist

 #        fakes_BF_h    = h.Integral( 26+offset,26+offset ) + h.Integral(33+offset,33+offset)
 #        fakes_CF_h    = h.Integral( 25+offset,25+offset ) + h.Integral(32+offset,32+offset)
 #        fakes_HFRes_h = h.Integral( 27+offset,29+offset )

 #        # Get the LF fakes for *this* hist

 #        fakes_LF_h = h.Integral( 23+offset,24+offset ) + h.Integral( 30+offset,31+offset )

 #        # Get the photon conversion fakes for *this* hist

 #        fakes_PhConv_h = h.Integral( 5+offset,5+offset )

 #        # Get the "Unknown" fakes for *this* Y

 #        fakes_Unknown_h = h.Integral( 0+offset,0+offset )

 #        # Get the other fakes for *this* hist

 #        fakes_Other_h = fakes_TOT_h - ( fakes_BF_h + fakes_CF_h + fakes_HFRes_h + fakes_LF_h + fakes_PhConv_h + fakes_Unknown_h )

 #        # Set the bin content for the fake lepton origin fraction hists for *this* hist

 #        if fakes_TOT_h:
 #            fakes_BF_frac_h      = fakes_BF_h/fakes_TOT_h
 #            fakes_CF_frac_h      = fakes_CF_h/fakes_TOT_h
 #            fakes_HFRes_frac_h   = fakes_HFRes_h/fakes_TOT_h
 #            fakes_LF_frac_h      = fakes_LF_h/fakes_TOT_h
 #            fakes_PhConv_frac_h  = fakes_PhConv_h/fakes_TOT_h
 #            fakes_Unknown_frac_h = fakes_Unknown_h/fakes_TOT_h
 #            fakes_Other_frac_h   = fakes_Other_h/fakes_TOT_h
 #        else:
 #            fakes_BF_frac_h = fakes_CF_frac_h = fakes_HFRes_frac_h = fakes_LF_frac_h = fakes_PhConv_frac_h = fakes_Unknown_frac_h = fakes_Other_frac_h = 0

 #        print("bin[{0}] ".format(idx))
 #        print("\ttot fakes = {0}".format(fakes_TOT_h))
 #        print("\t-) BF fakes = {0} ({1:.2f})".format(fakes_BF_h,fakes_BF_frac_h))
 #        print("\t-) CF fakes = {0} ({1:.2f})".format(fakes_CF_h,fakes_CF_frac_h))
 #        print("\t-) HFRes fakes = {0} ({1:.2f})".format(fakes_HFRes_h,fakes_HFRes_frac_h))
 #        print("\t-) LF fakes = {0} ({1:.2f})".format(fakes_LF_h,fakes_LF_frac_h))
 #        print("\t-) PhConv fakes = {0} ({1:.2f})".format(fakes_PhConv_h,fakes_PhConv_frac_h))
 #        print("\t-) Unknown fakes = {0} ({1:.2f})".format(fakes_Unknown_h,fakes_Unknown_frac_h))
 #        print("\t-) Other fakes = {0} ({1:.2f})".format(fakes_Other_h,fakes_Other_frac_h))

 #        histfakes_BF.SetBinContent( idx+offset, fakes_BF_frac_h )
 #        histfakes_CF.SetBinContent( idx+offset, fakes_CF_frac_h )
 #        histfakes_HFRes.SetBinContent( idx+offset, fakes_HFRes_frac_h )
 #        histfakes_LF.SetBinContent( idx+offset, fakes_LF_frac_h )
 #        histfakes_PhConv.SetBinContent( idx+offset, fakes_PhConv_frac_h )
 #        histfakes_Unknown.SetBinContent( idx+offset, fakes_Unknown_frac_h )
 #        histfakes_Other.SetBinContent( idx+offset, fakes_Other_frac_h )

 #    # Add histograms w/ fake origin fractions into a stack plot

 #    stacklegend.AddEntry(histfakes_PhConv, "#gamma conversion" , "F")
 #    stacklegend.AddEntry(histfakes_BF, "B-Had Fakes", "F")
 #    stacklegend.AddEntry(histfakes_CF, "C-Had Fakes", "F")
 #    stacklegend.AddEntry(histfakes_HFRes, "J/#psi,#Upsilon Fakes", "F")
 #    stacklegend.AddEntry(histfakes_LF, "L-Had Fakes", "F")
 #    stacklegend.AddEntry(histfakes_Other, "Other Fakes", "F")
 #    stacklegend.AddEntry(histfakes_Unknown, "Unknown" , "F")

 #    stack = THStack("FakeLepOriginFrac_VS_2LSplitFlavours_STACK","FakeLepOriginFrac_VS_2LSplitFlavours_STACK")
 #    stack.Add(histfakes_PhConv)
 #    stack.Add(histfakes_BF)
 #    stack.Add(histfakes_CF)
 #    stack.Add(histfakes_HFRes)
 #    stack.Add(histfakes_LF)
 #    stack.Add(histfakes_Other)
 #    stack.Add(histfakes_Unknown)

 #    canvasID_stack = "FakeLepOriginFrac_VS_2LSplitFlavours_" + kwargs["sample"] + "_" + kwargs["prodID"]
 #    cstack = TCanvas(canvasID_stack,canvasID_stack,50,50,1000,700)

 #    cstack.cd()
 #    stack.Draw()
 #    stack.GetXaxis().SetTitle("Category")
 #    stack.GetYaxis().SetTitle("Fake {0} origin fraction".format("lepton"))
 #    stack.GetXaxis().SetRangeUser(-0.5,stack.GetHistogram().GetNbinsX()-0.5)
 #    stack.GetXaxis().SetBinLabel(1,"ee, CR")
 #    stack.GetXaxis().SetBinLabel(2,"ee, SR")
 #    stack.GetXaxis().SetBinLabel(3,"e#mu, CR")
 #    stack.GetXaxis().SetBinLabel(4,"e#mu, SR")
 #    stack.GetXaxis().SetBinLabel(5,"#mue, CR")
 #    stack.GetXaxis().SetBinLabel(6,"#mue, SR")
 #    stack.GetXaxis().SetBinLabel(7,"#mu#mu, CR")
 #    stack.GetXaxis().SetBinLabel(8,"#mu#mu, SR")
 #    stack.GetYaxis().SetLabelSize(stack.GetXaxis().GetLabelSize())

 #    stacklegend.Draw()
 #    saveName_stack = canvasID_stack

 #    savePath = basePath + "/" + "FakeLepOriginFrac_VS_2LSplitFlavours_" + kwargs["prodID"] + "/"

 #    if not os.path.exists(savePath):
 #        os.makedirs(savePath)

 #    for ext in ["png","pdf","root"]:
 #        cstack.SaveAs( savePath + saveName_stack + "." + ext )

def plotFakeOriginFrac2LSplitFlavours( splitOF=False, **kwargs ):

    plotlist = []

    Plot.legend.SetHeader("Fake {0}".format("leptons"))

    variable = "LepFakeOrigin"

    # v28
    # basePath = os.path.abspath(os.curdir) + "/" + "PLOTS_" + kwargs["prodID"] + "/" + "OutputPlots_FakeOriginFrac_TTBarTTBarGamma_" + kwargs["prodID"] + "_SplitFlavours_TruthTP_HiggsApproval"
    # basePath = os.path.abspath(os.curdir) + "/" + "PLOTS_" + kwargs["prodID"] + "/" + "OutputPlots_FakeOriginFrac_TTBarTTBarGamma_" + kwargs["prodID"] + "_SplitFlavours_TruthTP_FINAL" # <--- use this!

    # v29
    # PP8
    basePath = os.path.abspath(os.curdir) + "/" + "PLOTS_" + kwargs["prodID"] + "/" + "OutputPlots_FakeOriginFrac_TTBarTTBarGamma_" + kwargs["prodID"] + "_SplitFlavours_TruthTP"
    # basePath = os.path.abspath(os.curdir) + "/" + "PLOTS_" + kwargs["prodID"] + "/" + "OutputPlots_FakeOriginFrac_TTBarTTBarGamma_" + kwargs["prodID"] + "_SplitFlavours_TruthTP_20GeV"
    # basePath = os.path.abspath(os.curdir) + "/" + "PLOTS_" + kwargs["prodID"] + "/" + "OutputPlots_FakeOriginFrac_TTBarTTBarGamma_" + kwargs["prodID"] + "_SplitFlavours_TruthTP_Corrections"
    # PP6
    # basePath = os.path.abspath(os.curdir) + "/" + "PLOTS_" + kwargs["prodID"] + "/" + "OutputPlots_FakeOriginFrac_TTBarTTBarGamma_" + kwargs["prodID"] + "_SplitFlavours_TruthTP_PP6"
    # basePath = os.path.abspath(os.curdir) + "/" + "PLOTS_" + kwargs["prodID"] + "/" + "OutputPlots_FakeOriginFrac_TTBarTTBarGamma_" + kwargs["prodID"] + "_SplitFlavours_TruthTP_PP6_Corrections"

    # ee
    inputPath_ee_VR = basePath + "/" + "ElFake" + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_ElEl" + "/" + "ElFake" + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_ElEl_"
    inputPath_ee_SR = basePath + "/" + "ElFake" + "FakeOriginFrac_2LSS_SR_TT_ElEl" + "/" + "ElFake" + "FakeOriginFrac_2LSS_SR_TT_ElEl_"
    # TEMP...
    # basePath_ee = os.path.abspath(os.curdir) + "/" + "PLOTS_" + kwargs["prodID"] + "/" + "OutputPlots_FakeOriginFrac_TTBarTTBarGamma_" + kwargs["prodID"] + "_SplitFlavours"
    # inputPath_ee_VR = basePath_ee + "/" + "El" + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_ElEl" + "/" + "El" + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_ElEl_"
    # inputPath_ee_SR = basePath_ee + "/" + "El" + "FakeOriginFrac_2LSS_SR_TT_ElEl" + "/" + "El" + "FakeOriginFrac_2LSS_SR_TT_ElEl_"
    # em
    inputPath_em_VR = basePath + "/" + "LepFake" + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_ElMu" + "/" + "LepFake" + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_ElMu_"
    inputPath_em_SR = basePath + "/" + "LepFake" + "FakeOriginFrac_2LSS_SR_TT_ElMu" + "/" + "LepFake" + "FakeOriginFrac_2LSS_SR_TT_ElMu_"
    # me
    inputPath_me_VR = basePath + "/" + "LepFake" + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_MuEl" + "/" + "LepFake" + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_MuEl_"
    inputPath_me_SR = basePath + "/" + "LepFake" + "FakeOriginFrac_2LSS_SR_TT_MuEl" + "/" + "LepFake" + "FakeOriginFrac_2LSS_SR_TT_MuEl_"
    # OF (em,me)
    if not splitOF:
        # --> use me for VR!
        inputPath_of_VR = basePath + "/" + "LepFake" + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_MuEl" + "/" + "LepFake" + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_MuEl_"
    else:
        inputPath_of_VR = basePath + "/" + "LepFake" + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_OF" + "/" + "LepFake" + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_OF_"
    inputPath_of_SR = basePath + "/" + "LepFake" + "FakeOriginFrac_2LSS_SR_TT_OF" + "/" + "LepFake" + "FakeOriginFrac_2LSS_SR_TT_OF_"
    # mm
    inputPath_mm_VR = basePath + "/" + "MuFake" + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_MuMu" + "/" + "MuFake" + "FakeOriginFrac_2LSS_LOWNJ_VR_TT_MuMu_"
    inputPath_mm_SR = basePath + "/" + "MuFake" + "FakeOriginFrac_2LSS_SR_TT_MuMu" + "/" + "MuFake" + "FakeOriginFrac_2LSS_SR_TT_MuMu_"

    inputName = variable + ".root"

    f_ee_VR = TFile(inputPath_ee_VR+inputName)
    h_ee_VR = f_ee_VR.Get(kwargs["sample"])
    h_ee_VR.SetDirectory(0)
    f_ee_SR = TFile(inputPath_ee_SR+inputName)
    h_ee_SR = f_ee_SR.Get(kwargs["sample"])
    h_ee_SR.SetDirectory(0)
    # TEMP
    # inputName_ee = "Lep1Origin.root"
    # f_ee_VR = TFile(inputPath_ee_VR+inputName_ee)
    # h_ee_VR = f_ee_VR.Get(kwargs["sample"])
    # h_ee_VR.SetDirectory(0)
    # f_ee_SR = TFile(inputPath_ee_SR+inputName_ee)
    # h_ee_SR = f_ee_SR.Get(kwargs["sample"])
    # h_ee_SR.SetDirectory(0)

    if splitOF:

        f_em_VR = TFile(inputPath_em_VR+inputName)
        h_em_VR = f_em_VR.Get(kwargs["sample"])
        h_em_VR.SetDirectory(0)
        f_em_SR = TFile(inputPath_em_SR+inputName)
        h_em_SR = f_em_SR.Get(kwargs["sample"])
        h_em_SR.SetDirectory(0)

        f_me_VR = TFile(inputPath_me_VR+inputName)
        h_me_VR = f_me_VR.Get(kwargs["sample"])
        h_me_VR.SetDirectory(0)
        f_me_SR = TFile(inputPath_me_SR+inputName)
        h_me_SR = f_me_SR.Get(kwargs["sample"])
        h_me_SR.SetDirectory(0)

    f_of_VR = TFile(inputPath_of_VR+inputName)
    h_of_VR = f_of_VR.Get(kwargs["sample"])
    h_of_VR.SetDirectory(0)
    f_of_SR = TFile(inputPath_of_SR+inputName)
    h_of_SR = f_of_SR.Get(kwargs["sample"])
    h_of_SR.SetDirectory(0)

    f_mm_VR = TFile(inputPath_mm_VR+inputName)
    h_mm_VR = f_mm_VR.Get(kwargs["sample"])
    h_mm_VR.SetDirectory(0)
    f_mm_SR = TFile(inputPath_mm_SR+inputName)
    h_mm_SR = f_mm_SR.Get(kwargs["sample"])
    h_mm_SR.SetDirectory(0)

    h_list = []
    h_list.extend([(h_ee_VR,"ee, CR"), (h_ee_SR,"ee, Pre-MVA")])
    if splitOF:
        h_list.extend([(h_em_VR,"em, CR"), (h_em_SR,"em, Pre-MVA")])
        h_list.extend([(h_me_VR,"me, CR"), (h_me_SR,"me, Pre-MVA")])
    h_list.extend([(h_of_VR,"OF, CR"), (h_of_SR,"OF, Pre-MVA")])
    h_list.extend([(h_mm_VR,"mm, CR"), (h_mm_SR,"mm, Pre-MVA")])

    # Fake lepton origin fraction wrt. 2LVR, 2LSR, 3LSR

    nbins = 6
    binmin = -0.5
    binmax = 5.5

    if splitOF:
        nbins = 10
        binmin = -0.5
        binmax = 9.5

    histfakes_BF        = TH1D("histfakes_BF","histfakes_BF",nbins,binmin,binmax) # B-hadrons in jets (mesons/baryons)
    histfakes_CF        = TH1D("histfakes_CF","histfakes_CF",nbins,binmin,binmax) # C-hadrons in jets (mesons/baryons)
    histfakes_HFRes     = TH1D("histfakes_HFRes","histfakes_HFRes",nbins,binmin,binmax) # B,C resonances (J/psi, Upsilon...)
    histfakes_LF        = TH1D("histfakes_LF","histfakes_LF",nbins,binmin,binmax) # Light hadrons in jets
    histfakes_PhConv    = TH1D("histfakes_PhConv","histfakes_PhConv",nbins,binmin,binmax) # Photon conversions
    histfakes_Other     = TH1D("histfakes_Other","histfakes_Other",nbins,binmin,binmax) # Other fakes (mis-id jets, leptons from generic pi/K...)
    histfakes_Unknown   = TH1D("histfakes_Unknown","histfakes_Unknown",nbins,binmin,binmax) # Unknown fakes (failure of MCTruthClassifier)

    stacklegend = TLegend(0.80,0.25,0.95,0.55)
    stacklegend.SetBorderSize(1)
    stacklegend.SetFillColor(kWhite)
    stacklegend.SetTextSize(0.03)
    stacklegend.SetTextFont(42)

    histfakes_list = [ (histfakes_BF,kRed), (histfakes_CF,kRed-9), (histfakes_HFRes,kPink-2), (histfakes_LF,kOrange+1), (histfakes_PhConv,kYellow), (histfakes_Other,kPink+1), (histfakes_Unknown,kAzure+1) ]

    for h in histfakes_list:
        h[0].SetLineWidth(2)
        h[0].SetLineStyle(1)
        h[0].SetLineColor(1)
        h[0].SetFillColor(h[1])

    for idx, H in enumerate(h_list):

        h    = H[0]
        name = H[1]

        offset = 1 # (to account for underflow bin, which has idx=0)

        # Get the tot. fakes for *this* hist

        fakes_TOT_h = h.Integral( 0, h.GetNbinsX()+1 )

        # Get the HF fakes for *this* hist

        fakes_BF_h    = h.Integral( 26+offset,26+offset ) + h.Integral(33+offset,33+offset)
        fakes_CF_h    = h.Integral( 25+offset,25+offset ) + h.Integral(32+offset,32+offset)
        fakes_HFRes_h = h.Integral( 27+offset,29+offset )

        # Get the LF fakes for *this* hist

        fakes_LF_h = h.Integral( 23+offset,24+offset ) + h.Integral( 30+offset,31+offset )

        # Get the photon conversion fakes for *this* hist

        fakes_PhConv_h = h.Integral( 5+offset,5+offset )

        # Get the "Unknown" fakes for *this* Y

        fakes_Unknown_h = h.Integral( 0+offset,0+offset )

        # Get the other fakes for *this* hist

        fakes_Other_h = fakes_TOT_h - ( fakes_BF_h + fakes_CF_h + fakes_HFRes_h + fakes_LF_h + fakes_PhConv_h + fakes_Unknown_h )

        # Set the bin content for the fake lepton origin fraction hists for *this* hist

        if fakes_TOT_h:
            fakes_BF_frac_h      = fakes_BF_h/fakes_TOT_h
            fakes_CF_frac_h      = fakes_CF_h/fakes_TOT_h
            fakes_HFRes_frac_h   = fakes_HFRes_h/fakes_TOT_h
            fakes_LF_frac_h      = fakes_LF_h/fakes_TOT_h
            fakes_PhConv_frac_h  = fakes_PhConv_h/fakes_TOT_h
            fakes_Unknown_frac_h = fakes_Unknown_h/fakes_TOT_h
            fakes_Other_frac_h   = fakes_Other_h/fakes_TOT_h
        else:
            fakes_BF_frac_h = fakes_CF_frac_h = fakes_HFRes_frac_h = fakes_LF_frac_h = fakes_PhConv_frac_h = fakes_Unknown_frac_h = fakes_Other_frac_h = 0

        print("bin[{0}] : {1}".format(idx,name))
        print("\ttot fakes = {0:.2f}".format(fakes_TOT_h))
        print("\t-) B-Had fakes = {0:.2f} ({1:.2f}%)".format(fakes_BF_h,fakes_BF_frac_h*1e2))
        print("\t-) C-Had fakes = {0:.2f} ({1:.2f}%)".format(fakes_CF_h,fakes_CF_frac_h*1e2))
        print("\t-) Jpsi/Upsilon fakes = {0:.2f} ({1:.2f}%)".format(fakes_HFRes_h,fakes_HFRes_frac_h*1e2))
        print("\t-) L-Had fakes = {0:.2f} ({1:.2f}%)".format(fakes_LF_h,fakes_LF_frac_h*1e2))
        print("\t-) PhConv fakes = {0:.2f} ({1:.2f}%)".format(fakes_PhConv_h,fakes_PhConv_frac_h*1e2))
        print("\t-) Other fakes = {0:.2f} ({1:.2f}%)".format(fakes_Other_h,fakes_Other_frac_h*1e2))
        print("\t-) Unknown fakes = {0:.2f} ({1:.2f}%)".format(fakes_Unknown_h,fakes_Unknown_frac_h*1e2))

        histfakes_BF.SetBinContent( idx+offset, fakes_BF_frac_h )
        histfakes_CF.SetBinContent( idx+offset, fakes_CF_frac_h )
        histfakes_HFRes.SetBinContent( idx+offset, fakes_HFRes_frac_h )
        histfakes_LF.SetBinContent( idx+offset, fakes_LF_frac_h )
        histfakes_PhConv.SetBinContent( idx+offset, fakes_PhConv_frac_h )
        histfakes_Unknown.SetBinContent( idx+offset, fakes_Unknown_frac_h )
        histfakes_Other.SetBinContent( idx+offset, fakes_Other_frac_h )

    # Add histograms w/ fake origin fractions into a stack plot

    stacklegend.AddEntry(histfakes_PhConv, "#gamma conversion" , "F")
    stacklegend.AddEntry(histfakes_BF, "B-Had Fakes", "F")
    stacklegend.AddEntry(histfakes_CF, "C-Had Fakes", "F")
    stacklegend.AddEntry(histfakes_HFRes, "J/#psi,#Upsilon Fakes", "F")
    stacklegend.AddEntry(histfakes_LF, "L-Had Fakes", "F")
    stacklegend.AddEntry(histfakes_Other, "Other Fakes", "F")
    stacklegend.AddEntry(histfakes_Unknown, "Unknown" , "F")

    stack = THStack("FakeLepOriginFrac_VS_2LSplitFlavours_STACK","FakeLepOriginFrac_VS_2LSplitFlavours_STACK")
    stack.Add(histfakes_PhConv)
    stack.Add(histfakes_BF)
    stack.Add(histfakes_CF)
    stack.Add(histfakes_HFRes)
    stack.Add(histfakes_LF)
    stack.Add(histfakes_Other)
    stack.Add(histfakes_Unknown)

    canvasID_stack = "FakeLepOriginFrac_VS_2LSplitFlavours_" + kwargs["sample"] + "_" + kwargs["prodID"]

    xcanv = 1200 if splitOF else 1000
    cstack = TCanvas(canvasID_stack,canvasID_stack,50,50,xcanv,700)

    cstack.cd()
    stack.Draw()
    stack.GetXaxis().SetTitle("Category")
    stack.GetYaxis().SetTitle("Fake {0} origin fraction".format("lepton"))
    stack.GetXaxis().SetRangeUser(-0.5,stack.GetHistogram().GetNbinsX()-0.5)
    if splitOF:
        stack.GetXaxis().SetBinLabel(1,"ee, CR")
        stack.GetXaxis().SetBinLabel(2,"ee, SR")
        stack.GetXaxis().SetBinLabel(3,"e#mu, CR")
        stack.GetXaxis().SetBinLabel(4,"e#mu, SR")
        stack.GetXaxis().SetBinLabel(5,"#mue, CR")
        stack.GetXaxis().SetBinLabel(6,"#mue, SR")
        stack.GetXaxis().SetBinLabel(7,"OF, CR")
        stack.GetXaxis().SetBinLabel(8,"OF, SR")
        stack.GetXaxis().SetBinLabel(9,"#mu#mu, CR")
        stack.GetXaxis().SetBinLabel(10,"#mu#mu, SR")
    else:
        stack.GetXaxis().SetBinLabel(1,"ee, CR")
        stack.GetXaxis().SetBinLabel(2,"ee, SR")
        stack.GetXaxis().SetBinLabel(3,"OF, CR")
        stack.GetXaxis().SetBinLabel(4,"OF, SR")
        stack.GetXaxis().SetBinLabel(5,"#mu#mu, CR")
        stack.GetXaxis().SetBinLabel(6,"#mu#mu, SR")
    stack.GetYaxis().SetLabelSize(stack.GetXaxis().GetLabelSize())

    stacklegend.Draw()
    saveName_stack = canvasID_stack

    if splitOF:
        savePath = basePath + "/" + "FakeLepOriginFrac_VS_2LSplitFlavours_ee_em_me_mm_" + kwargs["prodID"] + "/"
    else:
        savePath = basePath + "/" + "FakeLepOriginFrac_VS_2LSplitFlavours_" + kwargs["prodID"] + "/"

    if not os.path.exists(savePath):
        os.makedirs(savePath)

    for ext in ["png","pdf","root"]:
        cstack.SaveAs( savePath + saveName_stack + "." + ext )
