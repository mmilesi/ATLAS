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


def plotFakeOriginFrac2L3L( **kwargs ):

    plotlist = []

    if kwargs["flavour"] == "Mu":
        header   = "muon"
    elif kwargs["flavour"] == "El":
        header   = "electron"

    Plot.legend.SetHeader("Fake {0}".format(header))

    variable_2L = "Lep1Origin"
    variable_3L = "Lep2Origin"

    basePath = os.path.abspath(os.curdir) + "/" + "PLOTS_" + kwargs["prodID"] + "/" + "OutputPlots_FakeOriginFrac_TTBarTTBarGamma_" + kwargs["prodID"]

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

    histfakes_BF        = TH1D("histfakes_BF","histfakes_BF",3,-0.5,2.5)
    histfakes_CF        = TH1D("histfakes_CF","histfakes_CF",3,-0.5,2.5)
    histfakes_LF        = TH1D("histfakes_LF","histfakes_LF",3,-0.5,2.5)
    histfakes_PhConv    = TH1D("histfakes_PhConv","histfakes_PhConv",3,-0.5,2.5)
    histfakes_Other     = TH1D("histfakes_Other","histfakes_Other",3,-0.5,2.5)

    stacklegend = TLegend(0.5,0.3,0.75,0.6) # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
    stacklegend.SetBorderSize(1)
    stacklegend.SetFillColor(kWhite)
    stacklegend.SetTextSize(0.03)
    stacklegend.SetTextFont(42)

    histfakes_list = [ (histfakes_BF,kRed), (histfakes_CF,kRed-9), (histfakes_LF,kOrange+1), (histfakes_PhConv,kYellow), (histfakes_Other,kAzure+1) ]

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

        fakes_BF_h = h.Integral( 26+offset,26+offset ) + h.Integral( 29+offset,29+offset ) + h.Integral(33+offset,33+offset)
        fakes_CF_h = h.Integral( 25+offset,25+offset ) + h.Integral( 27+offset,28+offset ) + h.Integral(32+offset,32+offset)

        # Get the LF fakes for *this* hist

        fakes_LF_h = h.Integral( 23+offset,24+offset ) + h.Integral( 30+offset,31+offset )

        # Get the photon conversion fakes for *this* hist

        fakes_PhConv_h = h.Integral( 5+offset,5+offset )

        # Get the other fakes for *this* hist

        fakes_Other_h = fakes_TOT_h - ( fakes_BF_h + fakes_CF_h + fakes_LF_h + fakes_PhConv_h )

        # Set the bin content for the fake lepton origin fraction hists for *this* hist

        if fakes_TOT_h:
            fakes_BF_frac_h     = fakes_BF_h/fakes_TOT_h
            fakes_CF_frac_h     = fakes_CF_h/fakes_TOT_h
            fakes_LF_frac_h     = fakes_LF_h/fakes_TOT_h
            fakes_PhConv_frac_h = fakes_PhConv_h/fakes_TOT_h
            fakes_Other_frac_h  = fakes_Other_h/fakes_TOT_h
        else:
            fakes_BF_frac_h = fakes_CF_frac_h = fakes_LF_frac_h = fakes_PhConv_frac_h = fakes_Other_frac_h = 0

        print("bin[{0}] ".format(idx))
        print("\ttot fakes = {0}".format(fakes_TOT_h))
        print("\t-) BF fakes = {0} ({1:.2f})".format(fakes_BF_h,fakes_BF_frac_h))
        print("\t-) CF fakes = {0} ({1:.2f})".format(fakes_CF_h,fakes_CF_frac_h))
        print("\t-) LF fakes = {0} ({1:.2f})".format(fakes_LF_h,fakes_LF_frac_h))
        print("\t-) PhConv fakes = {0} ({1:.2f})".format(fakes_PhConv_h,fakes_PhConv_frac_h))
        print("\t-) Other fakes = {0} ({1:.2f})".format(fakes_Other_h,fakes_Other_frac_h))

        histfakes_BF.SetBinContent( idx+offset, fakes_BF_frac_h )
        histfakes_CF.SetBinContent( idx+offset, fakes_CF_frac_h )
        histfakes_LF.SetBinContent( idx+offset, fakes_LF_frac_h )
        histfakes_PhConv.SetBinContent( idx+offset, fakes_PhConv_frac_h )
        histfakes_Other.SetBinContent( idx+offset, fakes_Other_frac_h )

    # Add histograms w/ fake origin fractions into a stack plot

    stacklegend.AddEntry(histfakes_BF, "BF Fakes", "F")
    stacklegend.AddEntry(histfakes_CF, "CF Fakes", "F")
    stacklegend.AddEntry(histfakes_LF, "LF Fakes", "F")
    stacklegend.AddEntry(histfakes_PhConv, "#gamma conversion" , "F")
    stacklegend.AddEntry(histfakes_Other, "Other Fakes", "F")

    stack = THStack("FakeLepOriginFrac_VS_2L3LCat_STACK","FakeLepOriginFrac_VS_2L3LCat_STACK")
    stack.Add(histfakes_BF)
    stack.Add(histfakes_CF)
    stack.Add(histfakes_LF)
    stack.Add(histfakes_PhConv)
    stack.Add(histfakes_Other)

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

