#!/usr/bin/python

import array, math

from ROOT import gROOT, gDirectory, gStyle, gPad, TH1D, TFile, TCanvas, TColor, TLegend, TLatex, kBlack, kBlue, kOrange

gROOT.Reset()
gROOT.LoadMacro("$HOME/RootUtils/AtlasStyle.C")
from ROOT import SetAtlasStyle
SetAtlasStyle()

gROOT.SetBatch(True)

file_N = TFile("$ROOTCOREBIN/data/HTopMultilepAnalysis/External/QMisID_Pt_rates_Tight_v15.root")
file_D = TFile("$ROOTCOREBIN/data/HTopMultilepAnalysis/External/QMisID_Pt_rates_Loose_v15.root")

hist_N_orig = file_N.Get("LikelihoodPt")
hist_D_orig = file_D.Get("LikelihoodPt")

xbins_N = []
for ibinx in range(1,hist_N_orig.GetNbinsX()+1):
  xbins_N.append(hist_N_orig.GetBinLowEdge(ibinx))

xbins_D = []
for ibinx in range(1,hist_D_orig.GetNbinsX()+1):
  xbins_D.append(hist_D_orig.GetBinLowEdge(ibinx))

print "xbins_N = ", xbins_N
print "xbins_D = ", xbins_D

vxbins_N = array.array("d", xbins_N)
vxbins_D = array.array("d", xbins_D)

hist_N = TH1D("N","N",len(vxbins_N)-1,vxbins_N)
print("QMisID rate T - Numerator:")
for ibinx in range(1,hist_N_orig.GetNbinsX()):

  content = hist_N_orig.GetBinContent(ibinx)
  error   = hist_N_orig.GetBinError(ibinx)
  perc_error = abs(error) * 100 / content

  print("bin {0} - {1} +- {2} ({3} %)".format(ibinx,content,error,perc_error))
  hist_N.SetBinContent(ibinx,content)
  hist_N.SetBinError(ibinx,error)

hist_D = TH1D("D","D",len(vxbins_D)-1,vxbins_D)
print("QMisID rate L - Denominator:")
for ibinx in range(1,hist_D_orig.GetNbinsX()):

  content = hist_D_orig.GetBinContent(ibinx)
  error   = hist_D_orig.GetBinError(ibinx)
  perc_error = abs(error) * 100 / content

  print("bin {0} - {1} +- {2} ({3} %)".format(ibinx,content,error,perc_error))
  hist_D.SetBinContent(ibinx,content)
  hist_D.SetBinError(ibinx,error)

ratio = hist_N.Clone()
ratio.Divide(hist_N,hist_D,1.0,1.0,"B")
#ratio.Divide(hist_N,hist_D)

ratio.GetYaxis().SetTitle("(QMisID r_{T})/(QMisID r_{L})")
ratio.GetXaxis().SetTitle("p_{T} [GeV]")

ratio.GetYaxis().SetRangeUser(0.0,1.5)

c = TCanvas("c1","Temp",50,50,1000,800)

legend = TLegend(0.2,0.8,0.4,0.85); # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
legend.SetHeader("QMisID rate T/L ratio")
legend.SetBorderSize(0)  # no border
legend.SetFillStyle(0); # Legend transparent background
legend.SetTextSize(0.04) # Increase entry font size!
legend.SetTextFont(42)   # Helvetica
leg_ATLAS = TLatex()
leg_lumi  = TLatex()
leg_ATLAS.SetTextSize(0.03)
leg_ATLAS.SetNDC()
leg_lumi.SetTextSize(0.03)
leg_lumi.SetNDC()

ratio.Draw()

legend.Draw()
leg_ATLAS.DrawLatex(0.2,0.75,"#bf{#it{ATLAS}} Work In Progress")
leg_lumi.DrawLatex(0.2,0.7,"#sqrt{s} = 13 TeV, #int L dt = 6.7 fb^{-1}")

c.SaveAs("QMisID_Pt_rates_ratio_T_L.png")

# -------------------------------------

#file_eff = TFile("../OutputPlots_MMRates_25ns_v15/MMRates_25ns_v15_ScaledElFake/Rates.root")
file_eff = TFile("../OutputPlots_MMRates_25ns_v15/MMRates_25ns_v15_ScaledElFake_WEIGHTED_AVG/Rates.root")
#file_eff = TFile("../OutputPlots_MMRates_25ns_v15/MMRates_25ns_v15_ScaledElFake_ARITHMETIC_AVG/Rates.root")

hist_F_orig = file_eff.Get("El_ProbePt_Fake_Efficiency_observed")
hist_R_orig = file_eff.Get("El_ProbePt_RealQMisIDBinning_Efficiency_observed")

hist_QMisID = hist_R_orig.Clone()

print("ratio - GetNbinsX() = {0}".format(ratio.GetNbinsX()))
print("hist_R_orig - GetNbinsX() = {0}".format(hist_R_orig.GetNbinsX()))

for ibinx in range(1,hist_QMisID.GetNbinsX()+1):

  eff_r = hist_QMisID.GetBinContent(ibinx)
  err_eff_r = hist_QMisID.GetBinError(ibinx)

  scale = ratio.GetBinContent(ibinx)
  err_scale = ratio.GetBinError(ibinx)
  perc_err_scale = abs(err_scale) * 100 / scale

  eff_QMisID = scale * eff_r
  err_eff_QMisID = math.sqrt( (eff_r*eff_r) * (err_scale*err_scale) + (scale*scale) * (err_eff_r*err_eff_r) )

  print("bin {0}\n\teff_r: {1} +- {2}\n\tscale : {3} +- {4} ({5} %)\n\teff_QMisID: {6} +- {7}".format(ibinx,eff_r,err_eff_r,scale,err_scale,perc_err_scale,eff_QMisID,err_eff_QMisID))

  hist_QMisID.SetBinContent(ibinx, eff_QMisID)
  hist_QMisID.SetBinError(ibinx, err_eff_QMisID)

c2 = TCanvas("c2","Temp",50,50,1000,800)

legend2 = TLegend(0.4,0.3,0.5,0.35); # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
legend2.SetBorderSize(0)  # no border
legend2.SetFillStyle(0); # Legend transparent background
legend2.SetTextSize(0.04) # Increase entry font size!
legend2.SetTextFont(42)   # Helvetica

hist_F_orig.GetYaxis().SetRangeUser(0.0,0.7)
hist_F_orig.GetYaxis().SetTitle("Efficiency")

hist_F_orig.SetLineColor(kBlue)
hist_F_orig.SetMarkerColor(kBlue)
hist_F_orig.Draw()

hist_QMisID.SetLineColor(kOrange+7)
hist_QMisID.SetMarkerColor(kOrange+7)
hist_QMisID.Draw("SAME")

legend2.AddEntry(hist_F_orig, "#epsilon_{fake}", "P")
legend2.AddEntry(hist_QMisID, "#epsilon_{QMisID}", "P")

legend2.Draw()
leg_ATLAS.DrawLatex(0.6,0.35,"#bf{#it{ATLAS}} Work In Progress")
leg_lumi.DrawLatex(0.6,0.3,"#sqrt{s} = 13 TeV, #int L dt = 6.7 fb^{-1}")

c2.SaveAs("FakeEff_VS_QMisIDEff.png")

