#!/usr/bin/env python

""" EfficiencyPlotterClasses.py: plot efficiencies"""

__author__     = "Marco Milesi"
__email__      = "marco.milesi@cern.ch"
__maintainer__ = "Marco Milesi"

import os, sys, array, math

sys.path.append(os.path.abspath(os.path.curdir))

from ROOT import ROOT, gROOT, Double, TPad, TLine, TH1, TH1D, TFile, TCanvas, TLegend, TLatex, TGraphAsymmErrors, TEfficiency
from ROOT import kBlue, kOrange, kPink, kGreen, kRed, kYellow, kTeal, kMagenta, kViolet, kAzure, kCyan, kSpring, kGray, kBlack, kWhite
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

        for ext in ["png","eps","root"]:
	    c.SaveAs( savePath + "/" + saveName + "." + ext )

