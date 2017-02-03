#!/usr/bin/env python

""" EfficiencyPlotterClasses.py: plot efficiencies"""

__author__     = "Marco Milesi"
__email__      = "marco.milesi@cern.ch"
__maintainer__ = "Marco Milesi"

import os, sys, array, math

sys.path.append(os.path.abspath(os.path.curdir))

from ROOT import ROOT, gROOT, gStyle, gPad, Double, TPad, TLine, TH1, TH1D, TH2, THStack, TFile, TCanvas, TLegend, TLatex, TGraphAsymmErrors, TEfficiency
from ROOT import TPaletteAxis, TColor, kBlue, kOrange, kPink, kGreen, kRed, kYellow, kTeal, kMagenta, kViolet, kAzure, kCyan, kSpring, kGray, kBlack, kWhite
from ROOT import kFullCircle, kCircle, kOpenTriangleUp, kDot

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

        self.name = self.__hist.GetName()
        self.is2D = isinstance(self.__hist,TH2)

	self.__props  = properties

        # Set special properties for 2D histograms

        if self.is2D:
            self.set2DStyle(opt="BasicRainBow")
            from ROOT import kTRUE, kFALSE
            self.__hist.GetXaxis().SetNdivisions(self.__hist.GetNbinsX(),0,0,kTRUE)
            self.__hist.GetYaxis().SetNdivisions(self.__hist.GetNbinsY(),0,0,kTRUE)
            if self.__props["xAxisTitle"] in ["truthType","truthOrigin"]:
                self.__hist.GetXaxis().SetLabelSize(0.02)
            if self.__props["yAxisTitle"] in ["truthType","truthOrigin"]:
                self.__hist.GetYaxis().SetLabelSize(0.02)

        # Set here default properties (if they were not set by the user)

        if not self.__props.get("drawOpt"):
            if self.is2D:
                self.__props["drawOpt"] = "COLZ1"
            else:
                self.__props["drawOpt"] = "E0"

    def setProperty( self, propID, propValue ):

        self.__props[propID] = propValue

    def getProperty( self, propID ):

        if self.__props.get(propID):
            return self.__props[propID]

    def set2DStyle( self, opt = "BasicRainBow" ):

        gStyle.SetPadRightMargin(0.2) # Leave more space to the right side of the current Pad to show the histogram scale

        if opt == "FancyRainBow":
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

            s = array.array('d', [0.00, 0.34, 0.61, 0.84, 1.00])
            r = array.array('d', [0.00, 0.00, 0.87, 1.00, 0.51])
            g = array.array('d', [0.00, 0.81, 1.00, 0.20, 0.00])
            b = array.array('d', [0.51, 1.00, 0.12, 0.00, 0.00])

            npoints = len(s)
            TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
            gStyle.SetNumberContours(ncontours)

        if opt == "BasicRainBow":
            gStyle.SetPalette(1) # This resets the color palette to a simple Rainbow Color Map w/ 50 colors. See https://root.cern.ch/doc/master/classTColor.html

    def setAlphanumLabels( self, axis="X", labels=[], labelsopt="h" ):

        if axis == "X":
            n = self.__hist.GetNbinsX()
        elif axis == "Y":
            n = self.__hist.GetNbinsY()

        if n != len(labels):
            sys.exit("ERROR: {0}-axis - n = {1}, nlabels = {2}".format(axis,n,labels))

        for idx in range(1,n+1):
            if axis == "X":
                self.__hist.GetXaxis().SetBinLabel(idx,labels[idx-1][0])
            elif axis == "Y":
                self.__hist.GetYaxis().SetBinLabel(idx,labels[idx-1][0])

        if axis == "X":
            self.__hist.GetXaxis().LabelsOption(labelsopt)
        elif axis == "Y":
            self.__hist.GetYaxis().LabelsOption(labelsopt)


    def makeLeptonOriginFracPlots( self, histID=None ):

        # Fake lepton origin fraction wrt. njets

        histfakes_HF        = TH1D("histfakes_HF"+ "_"+histID,"histfakes_HF", self.__hist.GetNbinsY(),-0.5,self.__hist.GetNbinsY()-0.5)
        histfakes_LF        = TH1D("histfakes_LF"+"_"+histID,"histfakes_LF", self.__hist.GetNbinsY(),-0.5,self.__hist.GetNbinsY()-0.5)
        histfakes_PhConv    = TH1D("histfakes_PhConv"+"_"+histID,"histfakes_PhConv", self.__hist.GetNbinsY(),-0.5,self.__hist.GetNbinsY()-0.5)
        histfakes_Other     = TH1D("histfakes_Other"+"_"+histID,"histfakes_Other", self.__hist.GetNbinsY(),-0.5,self.__hist.GetNbinsY()-0.5)

        stacklegend = TLegend(0.5,0.3,0.75,0.6) # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
        stacklegend.SetBorderSize(1)
        stacklegend.SetFillColor(kWhite)
        stacklegend.SetTextSize(0.03)
        stacklegend.SetTextFont(42)

        histfakes_HF.SetLineWidth(3)
        histfakes_LF.SetLineWidth(3)
        histfakes_PhConv.SetLineWidth(3)
        histfakes_Other.SetLineWidth(3)

        histfakes_HF.SetLineStyle(1)
        histfakes_LF.SetLineStyle(1)
        histfakes_PhConv.SetLineStyle(1)
        histfakes_Other.SetLineStyle(1)

        histfakes_HF.SetLineColor(1)
        histfakes_LF.SetLineColor(1)
        histfakes_PhConv.SetLineColor(1)
        histfakes_Other.SetLineColor(1)

        histfakes_HF.SetFillColor(kRed)
        histfakes_LF.SetFillColor(kOrange+1)
        histfakes_PhConv.SetFillColor(kYellow)
        histfakes_Other.SetFillColor(kAzure+1)

        # Loop over jet multiplicity bins of the 2D hist

        for biny in range( 1, self.__hist.GetNbinsY()+1 ):

            offset = 1 # (to account for underflow bin, which has idx=0)

            # Get the tot. fakes for *this* nr. of jets

            fakes_TOT_biny = self.__hist.Integral( 0,  self.__hist.GetNbinsX()+1, biny, biny )

            # Get the HF fakes for *this* nr. of jets

            fakes_HF_biny = self.__hist.Integral( 25+offset,29+offset, biny, biny ) + self.__hist.Integral( 32+offset,33+offset, biny, biny )

            # Get the LF fakes for *this* nr. of jets

            fakes_LF_biny = self.__hist.Integral( 23+offset,24+offset, biny, biny ) + self.__hist.Integral( 30+offset,31+offset, biny, biny )

            # Get the photon conversion fakes for *this* nr. of jets

            fakes_PhConv_biny = self.__hist.Integral( 5+offset,5+offset, biny, biny )

            # Get the other fakes for *this* nr. of jets

            fakes_Other_biny = fakes_TOT_biny - ( fakes_HF_biny + fakes_LF_biny + fakes_PhConv_biny )

            # Set the bin content for the fake lepton origin fraction hists for *this* nr. jet bin

            if fakes_TOT_biny:
                fakes_HF_frac_biny     = fakes_HF_biny/fakes_TOT_biny
                fakes_LF_frac_biny     = fakes_LF_biny/fakes_TOT_biny
                fakes_PhConv_frac_biny = fakes_PhConv_biny/fakes_TOT_biny
                fakes_Other_frac_biny  = fakes_Other_biny/fakes_TOT_biny
            else:
                fakes_HF_frac_biny = fakes_LF_frac_biny = fakes_PhConv_frac_biny = fakes_Other_frac_biny = 0

            # print("bin[{0}] - njets = {1}".format(biny, biny-1))
            # print("\ttot fakes = {0}".format(fakes_TOT_biny))
            # print("\t-) HF fakes = {0} ({1:.2f})".format(fakes_HF_biny,fakes_HF_frac_biny))
            # print("\t-) LF fakes = {0} ({1:.2f})".format(fakes_LF_biny,fakes_LF_frac_biny))
            # print("\t-) PhConv fakes = {0} ({1:.2f})".format(fakes_PhConv_biny,fakes_PhConv_frac_biny))
            # print("\t-) Other fakes = {0} ({1:.2f})".format(fakes_Other_biny,fakes_Other_frac_biny))

            histfakes_HF.SetBinContent( biny, fakes_HF_frac_biny )
            histfakes_LF.SetBinContent( biny, fakes_LF_frac_biny )
            histfakes_PhConv.SetBinContent( biny, fakes_PhConv_frac_biny )
            histfakes_Other.SetBinContent( biny, fakes_Other_frac_biny )

        # Add histograms w/ fake origin fractions into a stack plot

        stacklegend.AddEntry(histfakes_HF, "HF Fakes", "F")
        stacklegend.AddEntry(histfakes_LF, "LF Fakes", "F")
        stacklegend.AddEntry(histfakes_PhConv, "#gamma conversion" , "F")
        stacklegend.AddEntry(histfakes_Other, "Other Fakes", "F")

        stack = THStack("LepOriginFrac_VS_NJets_STACK","LepOriginFrac_VS_NJets_STACK")
        stack.Add(histfakes_HF)
        stack.Add(histfakes_LF)
        stack.Add(histfakes_PhConv)
        stack.Add(histfakes_Other)

        return stack, stacklegend

    def makePlot( self ):

	if self.__props.get("xAxisTitle") : self.__hist.GetXaxis().SetTitle( self.__props["xAxisTitle"] )
	if self.__props.get("yAxisTitle") : self.__hist.GetYaxis().SetTitle( self.__props["yAxisTitle"] )

	if self.__props.get("xAxisRange") : self.__hist.GetXaxis().SetRangeUser( self.__props["xAxisRange"][0], self.__props["xAxisRange"][1] )
	if self.__props.get("yAxisRange") : self.__hist.GetYaxis().SetRangeUser( self.__props["yAxisRange"][0], self.__props["yAxisRange"][1] )

        if self.__props.get("normFactor"):
            normfactor = self.__props["normFactor"]
            if normfactor: self.__hist.Scale( normfactor / self.__hist.Integral() )

        if self.__props.get("legend")      : Plot.legend.AddEntry(self.__hist, self.__props["legend"], "P")

        if self.__props.get("setBinVals")  :
            for tup in self.__props["setBinVals"]:
                self.__hist.SetBinContent(tup[0],tup[1])
        if self.__props.get("setBinErrs")  :
            for tup in self.__props["setBinErrs"]:
                self.__hist.SetBinError(tup[0],tup[1])

	if not self.is2D:

            if self.__props.get("colour") :
                self.__hist.SetLineColor(self.__props["colour"])
                self.__hist.SetMarkerColor(self.__props["colour"])
            if self.__props.get("lineStyle")   : self.__hist.SetLineStyle(self.__props["lineStyle"])
            if self.__props.get("lineWidth")   : self.__hist.SetLineWidth(self.__props["lineWidth"])
            if self.__props.get("markerStyle") : self.__hist.SetMarkerStyle(self.__props["markerStyle"])
            if self.__props.get("markerSize")  : self.__hist.SetMarkerSize(self.__props["markerSize"])

        # Draw the histogram on the Pad!

	self.__hist.Draw( self.__props["drawOpt"] )

        # Deal with alphanumeric axis labels

        if self.__props.get("xAxisLabels"):
            self.setAlphanumLabels( axis="X", labels=self.__props["xAxisLabels"], labelsopt="v" )
            self.__hist.GetXaxis().SetTitleOffset(2.1) # increase a bit the axis title offset
            self.__hist.GetXaxis().SetLabelSize(0.03)
            self.__hist.GetXaxis().SetTitleSize(0.04)
        if self.__props.get("yAxisLabels"):
            self.setAlphanumLabels( axis="Y", labels=self.__props["yAxisLabels"], labelsopt="h" )
            self.__hist.GetYaxis().SetTitleOffset(2.255) # increase a bit the axis title offset
            self.__hist.GetYaxis().SetLabelSize(0.03)
            self.__hist.GetYaxis().SetTitleSize(0.04)

        # For 2D histograms, this will change the size of the labels on the COLZ axis
        # NB: must be done after calling Draw()!

        if self.is2D:
            gPad.Update()
            palette = self.__hist.FindObject("palette")
            palette.GetAxis().SetLabelSize(0.05)

        if self.__props.get("drawGrid"):
            if self.__props["drawGrid"]:
                gPad.Update()
                gPad.SetGrid()

class MultiPlot:

    def __init__( self, plots=[] ):

	self.__plotlist  = plots

    def makeMultiPlot( self, savePath, saveName ):

        tokens = saveName.split('_')

        c = TCanvas("c_"+tokens[0]+"_"+tokens[1],"Efficiency",50,50,1300,800)

        for idx, plot in enumerate(self.__plotlist):

            if plot.is2D:
                print("WARNING! Cannot plot a TH2 histogram together with other histograms. Skipping {0}".format(plot.name))
                continue

	    if idx == 0:
                #plot.makePlot(c)
                plot.makePlot()
	    else:
                plot.setProperty("drawOpt", str(plot.getProperty("drawOpt")) + " SAME" )
                #plot.makePlot(c)
                plot.makePlot()

        for refl in Plot.reflines:
            refl.SetLineStyle(2)
	    refl.Draw("SAME")

        Plot.legend.Draw()

	Plot.legendATLAS.DrawLatex(0.6,0.35,"#bf{#it{ATLAS}} Work In Progress")
        Plot.legendLumi.DrawLatex(0.6,0.27,"#sqrt{{s}} = 13 TeV, #int L dt = {0:.1f} fb^{{-1}}".format(Plot.luminosity))

        for ext in ["png","eps","root"]:
	    c.SaveAs( savePath + "/" + saveName + "." + ext )

