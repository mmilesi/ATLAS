#!/usr/bin/python

import os
import sys

from array import array

from ROOT import gROOT, gDirectory, gStyle, TH1D, TFile, TCanvas, TColor, TLegend, TLatex

gROOT.Reset()
gROOT.LoadMacro("/home/mmilesi/RootUtils/AtlasStyle.C")
from ROOT import SetAtlasStyle
SetAtlasStyle()

gROOT.SetBatch(True)

def plot():

    region = "ElEl"
    flavour = "El"
    CR = "RealCR"
    tightness = "T"
    variable = "Pt"

    subdir = "".join((region,CR,flavour,tightness))
    path = "../OutputPlots_MMRates_v028/" + subdir + "/"

    list_plots = ["".join((flavour,"Probe",variable,".root")),"".join((flavour,"0",variable,".root")),"".join((flavour,"1",variable,".root"))]

    list_files = []
    for item in list_plots:
        item = "".join((path,subdir,"_",item))
        list_files.append(item)

    hist_list = []
    for idx,item in enumerate(list_files):
        print("Now reading: {0}, index {1}".format(item,idx))
        myfile = TFile(item)
        myhist = myfile.Get("observed")
        hist_list.append(myhist)
        # This will keep the histogram alive in memory even after the TFile gets destroyed
        myhist.SetDirectory(0)

    c = TCanvas("c","Temp",50,50,700,900)

    legend = TLegend(0.6,0.8,0.8,0.85); # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
    header ="Lepton " + variable + " - " + subdir
    legend.SetHeader(header)
    legend.AddEntry(None, "", "")
    legend.SetBorderSize(0)  # no border
    legend.SetFillColor(0) # legend_L background should be white
    legend.SetTextSize(15) # Increase entry font size!
    legend.SetTextFont(43) # Helvetica

    outfile = open( path + "Yields_" + subdir + "_" + variable + ".txt", "w")

    for idx,myhist in enumerate(hist_list):

        print("hist name {0}".format(myhist.GetName()))

        outfile.write("%s - %s \n" %( list_plots[idx],myhist.GetName() ) )
        outfile.write("{ Nr. events: %s }; \n" %( myhist.Integral() ) )
        # get the bin number for pT = 60 GeV
        bmin = myhist.FindBin(60.0)
        # get the overflow bin
        bmax = myhist.GetNbinsX()+1
        integral_range = myhist.Integral(bmin,bmax)
        outfile.write("{ Nr. events (pT > 60 GeV): %s }; \n" %( integral_range ) )

        myhist.SetLineColor(2*idx+2)
        myhist.SetMarkerColor(2*idx+2)

        legend.AddEntry(myhist, list_plots[idx], "L")
        legend.AddEntry(None, "", "")

        if idx == 0:
            myhist.Draw("E0")
        else:
            myhist.Draw("E0 SAME")

    legend.Draw()

    c.SaveAs( path + "CompareData_" + subdir + "_" + variable + ".png")

    outfile.close()

plot()
