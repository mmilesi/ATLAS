#!/usr/bin/env python

""" EfficiencyPlotterClasses.py: plot efficiencies"""

__author__     = "Marco Milesi"
__email__      = "marco.milesi@cern.ch"
__maintainer__ = "Marco Milesi"

import os, sys, array, math

sys.path.append(os.path.abspath(os.path.curdir))

from ROOT import ROOT, gROOT, gStyle, gPad, Double, TPad, TLine, TH1, TH1D, TH2, TFile, TCanvas, TLegend, TLatex, TGraphAsymmErrors, TEfficiency
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


    def makePlot( self, canvas ):

        canvas.cd()

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
                plot.makePlot(c)
	    else:
                plot.setProperty("drawOpt", str(plot.getProperty("drawOpt")) + " SAME" )
                plot.makePlot(c)

        for refl in Plot.reflines:
            refl.SetLineStyle(2)
	    refl.Draw("SAME")

        Plot.legend.Draw()

	Plot.legendATLAS.DrawLatex(0.6,0.35,"#bf{#it{ATLAS}} Work In Progress")
        Plot.legendLumi.DrawLatex(0.6,0.27,"#sqrt{{s}} = 13 TeV, #int L dt = {0:.1f} fb^{{-1}}".format(Plot.luminosity))

        for ext in ["png","eps","root"]:
	    c.SaveAs( savePath + "/" + saveName + "." + ext )

