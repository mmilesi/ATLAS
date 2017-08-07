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
    # legend.SetTextFont(42)   # Helvetica

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

        self.conversion_frac_VS_Y = []

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

        # Fake lepton origin fraction wrt. X

        # Trick to make sure also the input histograms w/ variable bin size are handled correctly
        binsY = []
        for biny in range( 1, self.__hist.GetNbinsY()+2 ):
            lowedge = self.__hist.GetYaxis().GetBinLowEdge(biny)
            binsY.append(lowedge)
        arr_binsY = array.array("d", binsY)

        print "hist N binsY: ", self.__hist.GetNbinsY()
        print "Y axis bin lims: ", binsY, " - N bins: ", len(binsY)-1

        histfakes_BF        = TH1D("histfakes_BF"+ "_"+histID,"histfakes_BF", self.__hist.GetNbinsY(),arr_binsY) # B-hadrons in jets (mesons/baryons)
        histfakes_CF        = TH1D("histfakes_CF"+ "_"+histID,"histfakes_CF", self.__hist.GetNbinsY(),arr_binsY) # C-hadrons in jets (mesons/baryons)
        histfakes_HFRes     = TH1D("histfakes_HFRes"+ "_"+histID,"histfakes_HFRes", self.__hist.GetNbinsY(),arr_binsY) # B,C resonances (J/psi, Upsilon...)
        histfakes_LF        = TH1D("histfakes_LF"+"_"+histID,"histfakes_LF", self.__hist.GetNbinsY(),arr_binsY) # Light hadrons in jets
        histfakes_PhConv    = TH1D("histfakes_PhConv"+"_"+histID,"histfakes_PhConv", self.__hist.GetNbinsY(),arr_binsY) # Photon conversions
        histfakes_Other     = TH1D("histfakes_Other"+"_"+histID,"histfakes_Other", self.__hist.GetNbinsY(),arr_binsY) # Other fakes (mis-id jets, leptons from generic pi/K...)
        histfakes_Unknown   = TH1D("histfakes_Unknown"+"_"+histID,"histfakes_Unknown", self.__hist.GetNbinsY(),arr_binsY) # Unknown fakes (failure of MCTruthClassifier)

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

        # Loop over var Y bins of the 2D hist

        for biny in range( 1, self.__hist.GetNbinsY()+2 ):

            offset = 1 # (to account for underflow bin, which has idx=0)

            # Get the tot. fakes for *this* Y

            fakes_TOT_biny = self.__hist.Integral( 0,  self.__hist.GetNbinsX()+1, biny, biny )

            # Get the HF fakes for *this* Y

            fakes_BF_biny    = self.__hist.Integral( 26+offset,26+offset, biny, biny ) + self.__hist.Integral(33+offset,33+offset, biny, biny )
            fakes_CF_biny    = self.__hist.Integral( 25+offset,25+offset, biny, biny ) + self.__hist.Integral(32+offset,32+offset, biny, biny )
            fakes_HFRes_biny = self.__hist.Integral( 27+offset,29+offset, biny, biny )

            # Get the LF fakes for *this* Y

            fakes_LF_biny = self.__hist.Integral( 23+offset,24+offset, biny, biny ) + self.__hist.Integral( 30+offset,31+offset, biny, biny )

            # Get the photon conversion fakes for *this* Y

            fakes_PhConv_biny = self.__hist.Integral( 5+offset,5+offset, biny, biny )

            # Get the "Unknown" fakes for *this* Y

            fakes_Unknown_biny = self.__hist.Integral( 0+offset,0+offset, biny, biny )

            # Get the other fakes for *this* Y

            fakes_Other_biny = fakes_TOT_biny - ( fakes_BF_biny + fakes_CF_biny + fakes_HFRes_biny + fakes_LF_biny + fakes_PhConv_biny + fakes_Unknown_biny )

            # Set the bin content for the fake lepton origin fraction hists for *this* Y bin

            if fakes_TOT_biny:
                fakes_BF_frac_biny      = fakes_BF_biny/fakes_TOT_biny
                fakes_CF_frac_biny      = fakes_CF_biny/fakes_TOT_biny
                fakes_HFRes_frac_biny   = fakes_HFRes_biny/fakes_TOT_biny
                fakes_LF_frac_biny      = fakes_LF_biny/fakes_TOT_biny
                fakes_PhConv_frac_biny  = fakes_PhConv_biny/fakes_TOT_biny
                fakes_Unknown_frac_biny = fakes_Unknown_biny/fakes_TOT_biny
                fakes_Other_frac_biny   = fakes_Other_biny/fakes_TOT_biny
            else:
                fakes_BF_frac_biny = fakes_CF_frac_biny = fakes_HFRes_frac_biny = fakes_LF_frac_biny = fakes_PhConv_frac_biny = fakes_Unknown_frac_biny = fakes_Other_frac_biny = 0

            if False:
                print("varY - bin[{0}]".format(biny))
                print("\ttot fakes = {0}".format(fakes_TOT_biny))
                print("\t-) BF fakes = {0} ({1:.2f})".format(fakes_BF_biny,fakes_BF_frac_biny))
                print("\t-) CF fakes = {0} ({1:.2f})".format(fakes_CF_biny,fakes_CF_frac_biny))
                print("\t-) HFRes fakes = {0} ({1:.2f})".format(fakes_HFRes_biny,fakes_HFRes_frac_biny))
                print("\t-) LF fakes = {0} ({1:.2f})".format(fakes_LF_biny,fakes_LF_frac_biny))
                print("\t-) PhConv fakes = {0} ({1:.2f})".format(fakes_PhConv_biny,fakes_PhConv_frac_biny))
                print("\t-) Unknown fakes = {0} ({1:.2f})".format(fakes_Unknown_biny,fakes_Unknown_frac_biny))
                print("\t-) Other fakes = {0} ({1:.2f})".format(fakes_Other_biny,fakes_Other_frac_biny))

            histfakes_BF.SetBinContent( biny, fakes_BF_frac_biny )
            histfakes_CF.SetBinContent( biny, fakes_CF_frac_biny )
            histfakes_HFRes.SetBinContent( biny, fakes_HFRes_frac_biny )
            histfakes_LF.SetBinContent( biny, fakes_LF_frac_biny )
            histfakes_PhConv.SetBinContent( biny, fakes_PhConv_frac_biny )
            histfakes_Unknown.SetBinContent( biny, fakes_Unknown_frac_biny )
            histfakes_Other.SetBinContent( biny, fakes_Other_frac_biny )

            self.conversion_frac_VS_Y.append((biny,round(fakes_PhConv_frac_biny,3)))

        # Add histograms w/ fake origin fractions into a stack plot

        stacklegend.AddEntry(histfakes_BF, "B-Had Fakes", "F")
        stacklegend.AddEntry(histfakes_CF, "C-Had Fakes", "F")
        stacklegend.AddEntry(histfakes_HFRes, "J/#psi,#Upsilon Fakes", "F")
        stacklegend.AddEntry(histfakes_LF, "L-Had Fakes", "F")
        stacklegend.AddEntry(histfakes_PhConv, "#gamma conversion" , "F")
        stacklegend.AddEntry(histfakes_Other, "Other Fakes", "F")
        stacklegend.AddEntry(histfakes_Unknown, "Unknown" , "F")

        stack = THStack("LepOriginFrac_VS_Y_STACK","LepOriginFrac_VS_Y_STACK")
        stack.Add(histfakes_BF)
        stack.Add(histfakes_CF)
        stack.Add(histfakes_HFRes)
        stack.Add(histfakes_LF)
        stack.Add(histfakes_PhConv)
        stack.Add(histfakes_Other)
        stack.Add(histfakes_Unknown)

        return stack, stacklegend

    def makeConversionFracHist( self, histID=None ):

        # Photon conversion fraction wrt. X
        # Uncertainties from binmial distribution

        legend = TLegend(0.23,0.25,0.43,0.55)
        legend.SetBorderSize(1)
        legend.SetFillColor(kWhite)
        legend.SetTextSize(0.03)
        legend.SetTextFont(42)

        offset = 1 # (to account for underflow bin, which has idx=0)

        basename = self.__hist.GetName()

        # D: all fakes
        hist_fakes_TOT_projY = self.__hist.ProjectionY( basename+"_py_ALL_FAKES" )
        # N: conversion fakes
        hist_fakes_CONV_projY = self.__hist.ProjectionY( basename+"_py_CONV_FAKES", 5+offset, 5+offset )

        hist_fakes_CONV_FRAC = hist_fakes_CONV_projY.Clone("CONV_FRAC")
        hist_fakes_CONV_FRAC.Divide(hist_fakes_CONV_projY,hist_fakes_TOT_projY,1.0,1.0,"B")

        legend.AddEntry(hist_fakes_CONV_FRAC, "#gamma conversion", "F")

        print("Conversion fraction VS. Y:")
        for ibin in range(1,hist_fakes_CONV_FRAC.GetSize()-1):

            ibin_lowedge = hist_fakes_CONV_FRAC.GetXaxis().GetBinLowEdge(ibin)
            ibin_upedge  = hist_fakes_CONV_FRAC.GetXaxis().GetBinUpEdge(ibin)
            f_gamma = hist_fakes_CONV_FRAC.GetBinContent(ibin)
            f_gamma_err = hist_fakes_CONV_FRAC.GetBinError(ibin)

            print("\tbin: {0} [{1:.3f},{2:.3f}] - f_gamma = {3:.2f} +- {4:.2f}".format(ibin,ibin_lowedge,ibin_upedge,f_gamma,f_gamma_err))

        return hist_fakes_CONV_FRAC, legend


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
            if self.__props.get("lineStyle")    : self.__hist.SetLineStyle(self.__props["lineStyle"])
            if self.__props.get("lineWidth")    : self.__hist.SetLineWidth(self.__props["lineWidth"])
            if self.__props.get("markerStyle")  : self.__hist.SetMarkerStyle(self.__props["markerStyle"])
            if self.__props.get("markerSize")   : self.__hist.SetMarkerSize(self.__props["markerSize"])
            if self.__props.get("markerColour") : self.__hist.SetMarkerColor(self.__props["markerColour"])

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
            if palette:
                palette.GetAxis().SetLabelSize(0.05)

        if self.__props.get("drawGrid"):
            if self.__props["drawGrid"]:
                gPad.Update()
                gPad.SetGrid()

class MultiPlot:

    def __init__( self, plots=[] ):

	self.__plotlist  = plots

    def makeMultiPlot( self, savePath, saveName ):

        if not os.path.exists(savePath):
            print("Creating directory: {0}".format(savePath))
            os.makedirs(savePath)

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

        for ext in ["png","pdf","root"]:
	    c.SaveAs( savePath + "/" + saveName + "." + ext )

