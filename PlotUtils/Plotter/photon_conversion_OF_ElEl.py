# Photon conversion fraction - ee vs. OF (inclusive Nbjets)

frac_ee = [(1, 0.0), (2, 0.176), (3, 0.32), (4, 0.379), (5, 0.657), (6, 0.748), (7, 1.0)]
frac_of = [(1, 0.0), (2, 0.22), (3, 0.228), (4, 0.279), (5, 0.47), (6, 0.562), (7, 0.0)]

# Photon conversion fraction - ee vs. OF (Nbjet = 1)

# frac_ee = [(1, 0.0), (2, 0.068), (3, 0.201), (4, 0.331), (5, 0.559), (6, 0.748), (7, 0.0)]
# frac_of = [(1, 0.0), (2, 0.182), (3, 0.15), (4, 0.222), (5, 0.372), (6, 0.494), (7, 0.0)]

# Photon conversion fraction - ee vs. OF (Nbjet = 2)

# frac_ee = [(1, 0.0), (2, 0.883), (3, 0.795), (4, 0.658), (5, 0.829), (6, 0.0), (7, 1.0)]
# frac_of = [(1, 0.0), (2, 0.495), (3, 0.662), (4, 0.567), (5, 1.0), (6, 1.0), (7, 0.0)]

# alpha = ( frac_ee - frac_of ) / frac_of

alpha = []
for idx, elem in enumerate(frac_ee):
    f_ee  = elem[1]
    f_of  = frac_of[idx][1]
    r = 0
    n = f_ee - f_of
    d = f_of
    if n and d:
        r = n / d
    ratio = (idx+1,round(r,3))
    alpha.append(ratio)

print("alpha = {0}".format(alpha))

import array
from ROOT import gROOT, TFile, TCanvas, TH1F, TLatex, TLine, kRed

gROOT.LoadMacro("$HOME/RootUtils/AtlasStyle.C")
from ROOT import SetAtlasStyle
SetAtlasStyle()

gROOT.SetBatch(True)

c = TCanvas("alpha","alpha",50,50,800,600)
c.SetFrameFillColor(0)
c.SetFrameFillStyle(0)
c.SetFrameBorderMode(0)

n = len(alpha)
binlowedges = [10,15,20,26,35,60,210]
binlowedges_arr = array.array("f", binlowedges )

alpha_h = TH1F("alpha","alpha",n-1,binlowedges_arr)
alpha_h.GetXaxis().SetTitle("p_{T}^{e} [GeV]")
alpha_h.GetYaxis().SetTitle("1+#alpha")
alpha_h.SetMaximum(2.0)
alpha_h.SetMinimum(0.0)

for bin in range(1,alpha_h.GetSize()-1):
    content = 1 + alpha[bin-1][1]
    error = 0.15 * content
    print bin, content
    alpha_h.SetBinContent(bin,content)
    alpha_h.SetBinError(bin,error)

alpha_h.SetLineColor(kRed)
alpha_h.SetLineWidth(2)
alpha_h.SetMarkerColor(kRed)
alpha_h.SetMarkerSize(1)
alpha_h.SetMarkerStyle(20)
alpha_h.Draw("EP")

refl = TLine(alpha_h.GetXaxis().GetBinLowEdge(1),1.0,alpha_h.GetXaxis().GetBinLowEdge(n),1.0)
refl.SetLineStyle(2)
refl.Draw("SAME")

leg_alpha = TLatex()
leg_alpha.SetTextSize(0.06)
leg_alpha.SetNDC()
leg_alpha.DrawLatex(0.7,0.3,"#alpha=#frac{f_{#gamma}^{ee}-f_{#gamma}^{OF}}{f_{#gamma}^{OF}}")

outputpath = "../PLOTS_25ns_v27_v2/OutputPlots_FakeOriginFrac_TTBarTTBarGamma_25ns_v27"

for ext in ["pdf","png"]:
    c.SaveAs(outputpath+"/alpha."+ext)

outputfile = TFile(outputpath+"/"+"alpha.root","RECREATE")
outputfile.cd()
alpha_h.Write()
