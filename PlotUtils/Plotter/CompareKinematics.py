 #!/usr/bin/python

import array
import os
import sys
import argparse

# -------------------------------
# Parser for command line options
# -------------------------------

parser = argparse.ArgumentParser(description='Simple macro to compare kinematics/generic variables')

#***********************************
# positional arguments (compulsory!)
#***********************************
parser.add_argument('inputFile', metavar='inputFile',type=str,
                   help='full path to the input file')
#*******************
# optional arguments
#*******************
parser.add_argument('--channel', dest='channel', action='store', default='SF', type=str,
		    help='the channel chosen (OF, SF - default is SF)')
parser.add_argument('--doChFlipVeto', dest='doChFlipVeto', action='store_true',default=False,
                    help='apply veto on events w/ charge flip leptons')
parser.add_argument('--outdirname', dest='outdirname', action='store', default='', type=str,
		    help='specify a name to append to the output directory')
parser.add_argument('--doNorm', dest='doNorm', action='store_true',default=False,
                    help='normalise histograms to unity')

args = parser.parse_args()

from array import array

from ROOT import gROOT, gDirectory, gStyle, TH1D, TH2D, TFile, TCanvas, TColor, TLegend, TLatex

gROOT.Reset()
gROOT.LoadMacro("$HOME/RootUtils/AtlasStyle.C")
from ROOT import SetAtlasStyle
SetAtlasStyle()

gROOT.SetBatch(True)


# for fancy 2-dim histograms!
#
def set_fancy_2D_style():

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

  s = array('d', [0.00, 0.34, 0.61, 0.84, 1.00])
  r = array('d', [0.00, 0.00, 0.87, 1.00, 0.51])
  g = array('d', [0.00, 0.81, 1.00, 0.20, 0.00])
  b = array('d', [0.51, 1.00, 0.12, 0.00, 0.00])

  npoints = len(s)
  TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
  gStyle.SetNumberContours(ncontours)

def plot():

  myfile = TFile(args.inputFile)

  mytree = gDirectory.Get("physics")
  mytree_weight = mytree.GetWeight()/1e3

  trig_dec  = "( passHLT == 1 )"
  lep_tag_trigmatch = "( lep_tag_isTrigMatched[0] == 1 && ( ( lep_tag_flavour[0] == 11 && lep_tag_pt[0] > 25e3 ) || ( lep_tag_flavour[0] == 13 && lep_tag_pt[0] > 22e3 ) ) )"
  lep_probe_trigmatch = "( lep_probe_isTrigMatched[0] == 1 && ( ( lep_probe_flavour[0] == 11 && lep_probe_pt[0] > 25e3 ) || ( lep_probe_flavour[0] == 13 && lep_probe_pt[0] > 22e3 ) ) )"
  el_tag_eta   = "( TMath::Abs(el_tag_eta[0]) < 1.37 )"
  nbjets       = "( njets_mv2c20_Fix77 > 0 )"
  njets        = "( njets > 0 && njets < 4 )"
  nleptons     = "( nlep == 2 && ( lep_pt[0] > 20e3 && lep_pt[1] > 20e3 ) )"
  same_sign    = "( isSS01 == 1 )"
  # veto charge flips
  #
  ch_flip_veto = "( 1 )"
  if args.doChFlipVeto:
     ch_flip_veto = "( lep_isChFlip[0] == 0 && lep_isChFlip[1] == 0 )"
  # require at least 1 !prompt lepton
  #
  non_prompt   = "( ( lep_truthType[0] != 6 && lep_truthType[0] != 2 ) || ( lep_truthType[1] != 6 && lep_truthType[1] != 2 ) )"
  # require at least one T lepton in the event
  #
  is_tight_event = "1" #"( isNonTightEvent == 0 )"
  # require tag lepton to be T
  #
  tight_tag   = "( lep_tag_isTightSelected[0] == 1 )"

  hist_list_probe_MUPROBEEVT = {}
  hist_list_tag_MUPROBEEVT = {}
  hist_list_probe_ELPROBEEVT = {}
  hist_list_tag_ELPROBEEVT = {}

  # -----------------------------------------------------

  dirname = "OutputPlots_TRUTH"
  if args.outdirname:
     dirname += ( "_" + args.outdirname )

  try:
      os.makedirs(dirname)
  except:
      pass

  # -----------------------------------------------------

  # ----------------
  #    Muon probe
  # ----------------

  sel_FAKE_MUPROBEEVT_T = None
  sel_FAKE_MUPROBEEVT_L = None

  # look only at OF region
  if args.channel == "OF":
    sel_FAKE_MUPROBEEVT_T = "(" + trig_dec + " && " + lep_tag_trigmatch + " && " + lep_probe_trigmatch + " && " + nleptons + " && " + same_sign + " && " +  is_tight_event + " && " + tight_tag + " && " + "( (nmuon == 1 && nel == 1) && ( isProbeMuEvent == 1 ) && ( muon_probe_isTightSelected[0] == 1 ) )" + " && " + el_tag_eta + " && " + nbjets + " && " + njets  + " && " + ch_flip_veto + " && " + non_prompt + ")"
    sel_FAKE_MUPROBEEVT_L = "(" + trig_dec + " && " + lep_tag_trigmatch + " && " + lep_probe_trigmatch + " && " + nleptons + " && " + same_sign + " && " +  is_tight_event + " && " + tight_tag + " && " + "( (nmuon == 1 && nel == 1) && ( isProbeMuEvent == 1 ) && ( muon_probe_isTightSelected[0] == 0 ) )" + " && " + el_tag_eta + " && " + nbjets + " && " + njets + " && " + ch_flip_veto + " && " + non_prompt + ")"
  # look only at SF region
  elif args.channel == "SF":
    sel_FAKE_MUPROBEEVT_T = "(" + trig_dec + " && " + lep_tag_trigmatch + " && " + lep_probe_trigmatch + " && " + nleptons + " && " + same_sign + " && " +  is_tight_event + " && " + tight_tag + " && " + "( ( nmuon == 2 ) && ( isProbeMuEvent == 1 ) && ( muon_probe_isTightSelected[0] == 1 ) )" + " && " + nbjets + " && " + njets + " && " + ch_flip_veto + " && " + non_prompt + ")"
    sel_FAKE_MUPROBEEVT_L = "(" + trig_dec + " && " + lep_tag_trigmatch + " && " + lep_probe_trigmatch + " && " + nleptons + " && " + same_sign + " && " +  is_tight_event + " && " + tight_tag + " && " + "( ( nmuon == 2 ) && ( isProbeMuEvent == 1 ) && ( muon_probe_isTightSelected[0] == 0 ) )" + " && " + nbjets + " && " + njets + " && "  + ch_flip_veto + " && " + non_prompt + ")"

  print "sel_FAKE_MUPROBEEVT_T: \n", sel_FAKE_MUPROBEEVT_T
  print "sel_FAKE_MUPROBEEVT_L: \n", sel_FAKE_MUPROBEEVT_L

  # histograms for probe passing T
  #
  h_probelep_FAKE_MUPROBEEVT_T_pt = TH1D("probelep_FAKE_MUPROBEEVT_T_pt", "probelep_FAKE_MUPROBEEVT_T_pt", 30, 20.0, 200.0)
  hist_list_probe_MUPROBEEVT[h_probelep_FAKE_MUPROBEEVT_T_pt] = ["T","pt"]
  h_probelep_FAKE_MUPROBEEVT_T_ptvarcone30_over_pt = TH1D("probelep_FAKE_MUPROBEEVT_T_ptvarcone30_over_pt", "probelep_FAKE_MUPROBEEVT_T_ptvarcone30_over_pt", 20, 0.0, 0.1)
  hist_list_probe_MUPROBEEVT[h_probelep_FAKE_MUPROBEEVT_T_ptvarcone30_over_pt] = ["T","ptvarcone30_over_pt"]
#  h_probelep_FAKE_MUPROBEEVT_T_ptvarcone20_over_pt = TH1D("probelep_FAKE_MUPROBEEVT_T_ptvarcone20_over_pt", "probelep_FAKE_MUPROBEEVT_T_ptvarcone20_over_pt", 20, 0.0, 0.1)
#  hist_list_probe_MUPROBEEVT[h_probelep_FAKE_MUPROBEEVT_T_ptvarcone20_over_pt] = ["T","ptvarcone20_over_pt"]
  h_probelep_FAKE_MUPROBEEVT_T_topoetcone20_over_pt = TH1D("probelep_FAKE_MUPROBEEVT_T_topoetcone20_over_pt", "probelep_FAKE_MUPROBEEVT_T_topoetcone20_over_pt", 70, -0.2, 0.5)
  hist_list_probe_MUPROBEEVT[h_probelep_FAKE_MUPROBEEVT_T_topoetcone20_over_pt] = ["T","topoetcone20_over_pt"]

  h_probelep_FAKE_MUPROBEEVT_T_pt.GetXaxis().SetTitle("pT [GeV]")
  h_probelep_FAKE_MUPROBEEVT_T_ptvarcone30_over_pt.GetXaxis().SetTitle("ptvarcone30/pT")
#  h_probelep_FAKE_MUPROBEEVT_T_ptvarcone20_over_pt.GetXaxis().SetTitle("ptvarcone20/pT")
  h_probelep_FAKE_MUPROBEEVT_T_topoetcone20_over_pt.GetXaxis().SetTitle("topoetcone20/pT")

  # histograms for probe passing L!T
  #
  h_probelep_FAKE_MUPROBEEVT_L_pt = TH1D("probelep_FAKE_MUPROBEEVT_L_pt", "probelep_FAKE_MUPROBEEVT_L_pt", 30, 20.0, 200.0)
  hist_list_probe_MUPROBEEVT[h_probelep_FAKE_MUPROBEEVT_L_pt]  = ["L","pt"]
  h_probelep_FAKE_MUPROBEEVT_L_ptvarcone30_over_pt = TH1D("probelep_FAKE_MUPROBEEVT_L_ptvarcone30_over_pt", "probelep_FAKE_MUPROBEEVT_L_ptvarcone30_over_pt", 20, 0.0, 0.1)
  hist_list_probe_MUPROBEEVT[h_probelep_FAKE_MUPROBEEVT_L_ptvarcone30_over_pt] = ["L","ptvarcone30_over_pt"]
#  h_probelep_FAKE_MUPROBEEVT_L_ptvarcone20_over_pt = TH1D("probelep_FAKE_MUPROBEEVT_L_ptvarcone20_over_pt", "probelep_FAKE_MUPROBEEVT_L_ptvarcone20_over_pt", 20, 0.0, 0.1)
#  hist_list_probe_MUPROBEEVT[h_probelep_FAKE_MUPROBEEVT_L_ptvarcone20_over_pt] = ["L","ptvarcone20_over_pt"]
  h_probelep_FAKE_MUPROBEEVT_L_topoetcone20_over_pt = TH1D("probelep_FAKE_MUPROBEEVT_L_topoetcone20_over_pt", "probelep_FAKE_MUPROBEEVT_L_topoetcone20_over_pt", 70, -0.2, 0.5)
  hist_list_probe_MUPROBEEVT[h_probelep_FAKE_MUPROBEEVT_L_topoetcone20_over_pt] = ["L","topoetcone20_over_pt"]

  h_probelep_FAKE_MUPROBEEVT_L_pt.GetXaxis().SetTitle("pT [GeV]")
  h_probelep_FAKE_MUPROBEEVT_L_ptvarcone30_over_pt.GetXaxis().SetTitle("ptvarcone30/pT")
#  h_probelep_FAKE_MUPROBEEVT_L_ptvarcone20_over_pt.GetXaxis().SetTitle("ptvarcone20/pT")
  h_probelep_FAKE_MUPROBEEVT_L_topoetcone20_over_pt.GetXaxis().SetTitle("topoetcone20/pT")

  # histograms for tag (passing T: default)
  #
  h_taglep_FAKE_MUPROBEEVT_T_pt = TH1D("taglep_FAKE_MUPROBEEVT_T_pt", "taglep_FAKE_MUPROBEEVT_T_pt", 30, 20.0, 200.0)
  hist_list_tag_MUPROBEEVT[h_taglep_FAKE_MUPROBEEVT_T_pt] = ["T","pt"]
  h_taglep_FAKE_MUPROBEEVT_T_ptvarcone30_over_pt = TH1D("taglep_FAKE_MUPROBEEVT_T_ptvarcone30_over_pt", "taglep_FAKE_MUPROBEEVT_T_ptvarcone30_over_pt", 20, 0.0, 0.1)
  hist_list_tag_MUPROBEEVT[h_taglep_FAKE_MUPROBEEVT_T_ptvarcone30_over_pt] = ["T","ptvarcone30_over_pt"]
  h_taglep_FAKE_MUPROBEEVT_T_ptvarcone20_over_pt = TH1D("taglep_FAKE_MUPROBEEVT_T_ptvarcone20_over_pt", "taglep_FAKE_MUPROBEEVT_T_ptvarcone20_over_pt", 20, 0.0, 0.1)
  hist_list_tag_MUPROBEEVT[h_taglep_FAKE_MUPROBEEVT_T_ptvarcone20_over_pt] = ["T","ptvarcone20_over_pt"]
  h_taglep_FAKE_MUPROBEEVT_T_topoetcone20_over_pt = TH1D("taglep_FAKE_MUPROBEEVT_T_topoetcone20_over_pt", "taglep_FAKE_MUPROBEEVT_T_topoetcone20_over_pt", 70, -0.2, 0.5)
  hist_list_tag_MUPROBEEVT[h_taglep_FAKE_MUPROBEEVT_T_topoetcone20_over_pt] = ["T","topoetcone20_over_pt"]

  h_taglep_FAKE_MUPROBEEVT_T_pt.GetXaxis().SetTitle("pT [GeV]")
  h_taglep_FAKE_MUPROBEEVT_T_ptvarcone30_over_pt.GetXaxis().SetTitle("ptvarcone30/pT")
  h_taglep_FAKE_MUPROBEEVT_T_topoetcone20_over_pt.GetXaxis().SetTitle("topoetcone20/pT")

  mytree.Project("probelep_FAKE_MUPROBEEVT_T_pt", "muon_probe_pt[0]/1e3", "%s * (%s)" %(mytree_weight, sel_FAKE_MUPROBEEVT_T), "" )
  mytree.Project("probelep_FAKE_MUPROBEEVT_L_pt", "muon_probe_pt[0]/1e3", "%s * (%s)" %(mytree_weight, sel_FAKE_MUPROBEEVT_L), "" )
  mytree.Project("probelep_FAKE_MUPROBEEVT_T_ptvarcone30_over_pt", "muon_probe_ptvarcone30[0]/muon_probe_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_MUPROBEEVT_T), "" )
  mytree.Project("probelep_FAKE_MUPROBEEVT_L_ptvarcone30_over_pt", "muon_probe_ptvarcone30[0]/muon_probe_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_MUPROBEEVT_L), "" )
#  mytree.Project("probelep_FAKE_MUPROBEEVT_T_ptvarcone20_over_pt", "muon_probe_ptvarcone20[0]/muon_probe_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_MUPROBEEVT_T), "" )
#  mytree.Project("probelep_FAKE_MUPROBEEVT_L_ptvarcone20_over_pt", "muon_probe_ptvarcone20[0]/muon_probe_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_MUPROBEEVT_L), "" )
  mytree.Project("probelep_FAKE_MUPROBEEVT_T_topoetcone20_over_pt", "muon_probe_topoetcone20[0]/muon_probe_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_MUPROBEEVT_T), "" )
  mytree.Project("probelep_FAKE_MUPROBEEVT_L_topoetcone20_over_pt", "muon_probe_topoetcone20[0]/muon_probe_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_MUPROBEEVT_L), "" )

  if args.channel == "OF":
     # probe is muon, tag is electron
     #
     mytree.Project("taglep_FAKE_MUPROBEEVT_T_pt", "el_tag_pt[0]/1e3", "%s * (%s)" %(mytree_weight, sel_FAKE_MUPROBEEVT_T), "" )
     mytree.Project("taglep_FAKE_MUPROBEEVT_T_ptvarcone30_over_pt", "el_tag_ptvarcone30[0]/el_tag_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_MUPROBEEVT_T), "" )
     mytree.Project("taglep_FAKE_MUPROBEEVT_T_ptvarcone20_over_pt", "el_tag_ptvarcone20[0]/el_tag_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_MUPROBEEVT_T), "" )
     mytree.Project("taglep_FAKE_MUPROBEEVT_T_topoetcone20_over_pt", "el_tag_topoetcone20[0]/el_tag_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_MUPROBEEVT_T), "" )
  elif args.channel == "SF":
     # probe is muon, tag is muon
     #
     mytree.Project("taglep_FAKE_MUPROBEEVT_T_pt", "muon_tag_pt[0]/1e3", "%s * (%s)" %(mytree_weight, sel_FAKE_MUPROBEEVT_T), "" )
     mytree.Project("taglep_FAKE_MUPROBEEVT_T_ptvarcone30_over_pt", "muon_tag_ptvarcone30[0]/muon_tag_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_MUPROBEEVT_T), "" )
#     mytree.Project("taglep_FAKE_MUPROBEEVT_T_ptvarcone20_over_pt", "muon_tag_ptvarcone20[0]/muon_tag_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_MUPROBEEVT_T), "" )
     mytree.Project("taglep_FAKE_MUPROBEEVT_T_topoetcone20_over_pt", "muon_tag_topoetcone20[0]/muon_tag_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_MUPROBEEVT_T), "" )


  # Normalise to unity
  if args.doNorm:
     for hist in hist_list_probe_MUPROBEEVT.keys():
	if hist.Integral() == 0:
	  continue
	hist.Scale(1/hist.Integral())
	hist.GetYaxis().SetTitle("Norm. events")
	#hist.GetYaxis().SetRangeUser(0,1)
     for hist in hist_list_tag_MUPROBEEVT.keys():
	if hist.Integral() == 0:
	  continue
        hist.Scale(1/hist.Integral())
	hist.GetYaxis().SetTitle("Norm. events")
        #hist.GetYaxis().SetRangeUser(0,1)

  # now plot
  #
  print "\tLooking at events where the PROBE is: {0} and the channel is {1}\n".format("MUON", args.channel)
  for  hist_tag in hist_list_tag_MUPROBEEVT.keys():

     c = TCanvas("c","Temp",50,50,700,900)
     c.Divide(1,2)

     c.cd(1)
     legend_T = TLegend(0.6,0.8,0.8,0.85); # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
     legend_T.SetHeader("2 Lep SS Fake CR - T probe")
     legend_T.SetBorderSize(0)  # no border
     legend_T.SetFillColor(0) # Legend background should be white
     legend_T.SetTextSize(0.04) # Increase entry font size!
     legend_T.SetTextFont(42) # Helvetica

     c.cd(2)
     legend_L = TLegend(0.6,0.8,0.8,0.85); # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
     legend_L.SetHeader("2 Lep SS Fake CR - L!T probe")
     legend_L.SetBorderSize(0)  # no border
     legend_L.SetFillColor(0) # legend_L background should be white
     legend_L.SetTextSize(0.04) # Increase entry font size!
     legend_L.SetTextFont(42) # Helvetica

     var_tag    = hist_list_tag_MUPROBEEVT[hist_tag][1]
     type_tag   = hist_list_tag_MUPROBEEVT[hist_tag][0]

     hist_tag.SetLineColor(1)
     hist_tag.SetLineWidth(2)

     print "\tPlotting variable: {0}\n".format(var_tag)

     print "\ttag histogram name: {0}\n".format(hist_tag.GetName())

     draw_pad1 = False
     draw_pad2 = False

     for hist_probe in hist_list_probe_MUPROBEEVT.keys():

       print "\t\tprobe histogram name: {0}".format(hist_probe.GetName())

       hist_probe.SetLineColor(2)
       hist_probe.SetLineWidth(2)

       var_probe  = hist_list_probe_MUPROBEEVT[hist_probe][1]
       type_probe = hist_list_probe_MUPROBEEVT[hist_probe][0]

       # must be looking at the same variable!
       #
       if  var_probe != var_tag:
         continue

       print "\t\tfound a match in variable! {0} = {1}\n".format(var_tag,var_probe)

       if type_probe == "T":
          print "\t\t  probe is of type: {0}".format(type_probe)
          print "\n\t\t  --> drawing: \t {0} \t {1} \t on pad(1) of the same canvas\n".format(hist_tag.GetName(),hist_probe.GetName())
          c.cd(1)
          legend_T.AddEntry(None, "", "")
          legend_T.AddEntry(hist_tag, "REAL lepton", "L")
	  legend_T.AddEntry(None, "", "")
          legend_T.AddEntry(hist_probe, "FAKE lepton", "L")
          hist_tag.Draw("HIST")
          hist_probe.Draw("HIST SAME")
          draw_pad1 = True

       elif type_probe == "L":
          print "\t\t  probe is of type: {0}".format(type_probe)
          print "\n\t\t  --> drawing: \t {0} \t {1} \t on pad(2) of the same canvas\n".format(hist_tag.GetName(),hist_probe.GetName())
          c.cd(2)
          legend_L.AddEntry(None, "", "")
	  legend_L.AddEntry(hist_tag, "REAL lepton", "L")
	  legend_L.AddEntry(None, "", "")
          legend_L.AddEntry(hist_probe, "FAKE lepton", "L")
          hist_tag.Draw("HIST")
          hist_probe.Draw("HIST SAME")
          draw_pad2 = True

       if ( draw_pad1 and draw_pad2 ):
	 c.cd(1)
         legend_T.Draw()
	 c.cd(2)
	 legend_L.Draw()
         plotname = dirname + "/FAKE_MUPROBEEVT_" + var_probe
         c.SaveAs(plotname+".png")

     del c

  # -----------------------------------------------------

  # ---------------------
  #    Electron probe
  # ---------------------

  sel_FAKE_ELPROBEEVT_T = None
  sel_FAKE_ELPROBEEVT_L = None

  # look only at OF region
  if args.channel == "OF":
    sel_FAKE_ELPROBEEVT_T = "(" + trig_dec + " && " + lep_tag_trigmatch + " && " + lep_probe_trigmatch + " && " + nleptons + " && " + same_sign + " && " +  is_tight_event + " && " + tight_tag + " && " + "( (nmuon == 1 && nel == 1) && ( isProbeElEvent == 1 ) && ( el_probe_isTightSelected[0] == 1 ) )" + " && " + nbjets + " && " + njets + " && "  + ch_flip_veto + " && " + non_prompt + ")"
    sel_FAKE_ELPROBEEVT_L = "(" + trig_dec + " && " + lep_tag_trigmatch + " && " + lep_probe_trigmatch + " && " + nleptons + " && " + same_sign + " && " +  is_tight_event + " && " + tight_tag + " && " + "( (nmuon == 1 && nel == 1) && ( isProbeElEvent == 1 ) && ( el_probe_isTightSelected[0] == 0 ) )" + " && " + nbjets + " && " + njets + " && " + ch_flip_veto + " && " + non_prompt + ")"



  # look only at SF region
  elif args.channel == "SF":
    sel_FAKE_ELPROBEEVT_T = "(" + trig_dec + " && " + lep_tag_trigmatch + " && " + lep_probe_trigmatch + " && " + nleptons + " && " + same_sign + " && " +  is_tight_event + " && " + tight_tag + " && " + "( ( nel == 2 ) && ( isProbeElEvent == 1 ) && ( el_probe_isTightSelected[0] == 1 ) )" + " && " + el_tag_eta + " && " + nbjets + " && " + njets + " && "  + ch_flip_veto + " && " + non_prompt + ")"
    sel_FAKE_ELPROBEEVT_L = "(" + trig_dec + " && " + lep_tag_trigmatch + " && " + lep_probe_trigmatch + " && " + nleptons + " && " + same_sign + " && " +  is_tight_event + " && " + tight_tag + " && " + "( ( nel == 2 ) && ( isProbeElEvent == 1 ) && ( el_probe_isTightSelected[0] == 0 ) )" + " && " + el_tag_eta + " && " + nbjets + " && " + njets + " && " + ch_flip_veto + " && " + non_prompt + ")"

  print "sel_FAKE_ELPROBEEVT_T: \n", sel_FAKE_ELPROBEEVT_T
  print "sel_FAKE_ELPROBEEVT_L: \n", sel_FAKE_ELPROBEEVT_L

  # histograms for probe passing T
  #
  h_probelep_FAKE_ELPROBEEVT_T_pt = TH1D("probelep_FAKE_ELPROBEEVT_T_pt", "probelep_FAKE_ELPROBEEVT_T_pt", 30, 20.0, 200.0)
  hist_list_probe_ELPROBEEVT[h_probelep_FAKE_ELPROBEEVT_T_pt] = ["T","pt"]
#  h_probelep_FAKE_ELPROBEEVT_T_ptvarcone30_over_pt = TH1D("probelep_FAKE_ELPROBEEVT_T_ptvarcone30_over_pt", "probelep_FAKE_ELPROBEEVT_T_ptvarcone30_over_pt", 20, 0.0, 0.1)
#  hist_list_probe_ELPROBEEVT[h_probelep_FAKE_ELPROBEEVT_T_ptvarcone30_over_pt] = ["T","ptvarcone30_over_pt"]
  h_probelep_FAKE_ELPROBEEVT_T_ptvarcone20_over_pt = TH1D("probelep_FAKE_ELPROBEEVT_T_ptvarcone20_over_pt", "probelep_FAKE_ELPROBEEVT_T_ptvarcone20_over_pt", 20, 0.0, 0.1)
  hist_list_probe_ELPROBEEVT[h_probelep_FAKE_ELPROBEEVT_T_ptvarcone20_over_pt] = ["T","ptvarcone20_over_pt"]
  h_probelep_FAKE_ELPROBEEVT_T_topoetcone20_over_pt = TH1D("probelep_FAKE_ELPROBEEVT_T_topoetcone20_over_pt", "probelep_FAKE_ELPROBEEVT_T_topoetcone20_over_pt", 70, -0.2, 0.5)
  hist_list_probe_ELPROBEEVT[h_probelep_FAKE_ELPROBEEVT_T_topoetcone20_over_pt] = ["T","topoetcone20_over_pt"]

  h_probelep_FAKE_ELPROBEEVT_T_pt.GetXaxis().SetTitle("pT [GeV]")
#  h_probelep_FAKE_ELPROBEEVT_T_ptvarcone30_over_pt.GetXaxis().SetTitle("ptvarcone30/pT")
  h_probelep_FAKE_ELPROBEEVT_T_ptvarcone20_over_pt.GetXaxis().SetTitle("ptvarcone20/pT")
  h_probelep_FAKE_ELPROBEEVT_T_topoetcone20_over_pt.GetXaxis().SetTitle("topoetcone20/pT")

  # histograms for probe passing L!T
  #
  h_probelep_FAKE_ELPROBEEVT_L_pt = TH1D("probelep_FAKE_ELPROBEEVT_L_pt", "probelep_FAKE_ELPROBEEVT_L_pt", 30, 20.0, 200.0)
  hist_list_probe_ELPROBEEVT[h_probelep_FAKE_ELPROBEEVT_L_pt]  = ["L","pt"]
#  h_probelep_FAKE_ELPROBEEVT_L_ptvarcone30_over_pt = TH1D("probelep_FAKE_ELPROBEEVT_L_ptvarcone30_over_pt", "probelep_FAKE_ELPROBEEVT_L_ptvarcone30_over_pt", 20, 0.0, 0.1)
#  hist_list_probe_ELPROBEEVT[h_probelep_FAKE_ELPROBEEVT_L_ptvarcone30_over_pt] = ["L","ptvarcone30_over_pt"]
  h_probelep_FAKE_ELPROBEEVT_L_ptvarcone20_over_pt = TH1D("probelep_FAKE_ELPROBEEVT_L_ptvarcone20_over_pt", "probelep_FAKE_ELPROBEEVT_L_ptvarcone20_over_pt", 20, 0.0, 0.1)
  hist_list_probe_ELPROBEEVT[h_probelep_FAKE_ELPROBEEVT_L_ptvarcone20_over_pt] = ["L","ptvarcone20_over_pt"]
  h_probelep_FAKE_ELPROBEEVT_L_topoetcone20_over_pt = TH1D("probelep_FAKE_ELPROBEEVT_L_topoetcone20_over_pt", "probelep_FAKE_ELPROBEEVT_L_topoetcone20_over_pt", 70, -0.2, 0.5)
  hist_list_probe_ELPROBEEVT[h_probelep_FAKE_ELPROBEEVT_L_topoetcone20_over_pt] = ["L","topoetcone20_over_pt"]

  h_probelep_FAKE_ELPROBEEVT_L_pt.GetXaxis().SetTitle("pT [GeV]")
#  h_probelep_FAKE_ELPROBEEVT_L_ptvarcone30_over_pt.GetXaxis().SetTitle("ptvarcone30/pT")
  h_probelep_FAKE_ELPROBEEVT_L_ptvarcone20_over_pt.GetXaxis().SetTitle("ptvarcone20/pT")
  h_probelep_FAKE_ELPROBEEVT_L_topoetcone20_over_pt.GetXaxis().SetTitle("topoetcone20/pT")

  # histograms for tag (passing T: default)
  #
  h_taglep_FAKE_ELPROBEEVT_T_pt = TH1D("taglep_FAKE_ELPROBEEVT_T_pt", "taglep_FAKE_ELPROBEEVT_T_pt", 30, 20.0, 200.0)
  hist_list_tag_ELPROBEEVT[h_taglep_FAKE_ELPROBEEVT_T_pt] = ["T","pt"]
  h_taglep_FAKE_ELPROBEEVT_T_ptvarcone30_over_pt = TH1D("taglep_FAKE_ELPROBEEVT_T_ptvarcone30_over_pt", "taglep_FAKE_ELPROBEEVT_T_ptvarcone30_over_pt", 20, 0.0, 0.1)
  hist_list_tag_ELPROBEEVT[h_taglep_FAKE_ELPROBEEVT_T_ptvarcone30_over_pt] = ["T","ptvarcone30_over_pt"]
  h_taglep_FAKE_ELPROBEEVT_T_ptvarcone20_over_pt = TH1D("taglep_FAKE_ELPROBEEVT_T_ptvarcone20_over_pt", "taglep_FAKE_ELPROBEEVT_T_ptvarcone20_over_pt", 20, 0.0, 0.1)
  hist_list_tag_ELPROBEEVT[h_taglep_FAKE_ELPROBEEVT_T_ptvarcone20_over_pt] = ["T","ptvarcone20_over_pt"]
  h_taglep_FAKE_ELPROBEEVT_T_topoetcone20_over_pt = TH1D("taglep_FAKE_ELPROBEEVT_T_topoetcone20_over_pt", "taglep_FAKE_ELPROBEEVT_T_topoetcone20_over_pt", 70, -0.2, 0.5)
  hist_list_tag_ELPROBEEVT[h_taglep_FAKE_ELPROBEEVT_T_topoetcone20_over_pt] = ["T","topoetcone20_over_pt"]

  h_taglep_FAKE_ELPROBEEVT_T_pt.GetXaxis().SetTitle("pT [GeV]")
  h_taglep_FAKE_ELPROBEEVT_T_ptvarcone30_over_pt.GetXaxis().SetTitle("ptvarcone30/pT")
  h_taglep_FAKE_ELPROBEEVT_T_ptvarcone20_over_pt.GetXaxis().SetTitle("ptvarcone20/pT")
  h_taglep_FAKE_ELPROBEEVT_T_topoetcone20_over_pt.GetXaxis().SetTitle("topoetcone20/pT")

  mytree.Project("probelep_FAKE_ELPROBEEVT_T_pt", "el_probe_pt[0]/1e3", "%s * (%s)" %(mytree_weight, sel_FAKE_ELPROBEEVT_T), "" )
  mytree.Project("probelep_FAKE_ELPROBEEVT_L_pt", "el_probe_pt[0]/1e3", "%s * (%s)" %(mytree_weight, sel_FAKE_ELPROBEEVT_L), "" )
#  mytree.Project("probelep_FAKE_ELPROBEEVT_T_ptvarcone30_over_pt", "el_probe_ptvarcone30[0]/el_probe_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_ELPROBEEVT_T), "" )
#  mytree.Project("probelep_FAKE_ELPROBEEVT_L_ptvarcone30_over_pt", "el_probe_ptvarcone30[0]/el_probe_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_ELPROBEEVT_L), "" )
  mytree.Project("probelep_FAKE_ELPROBEEVT_T_ptvarcone20_over_pt", "el_probe_ptvarcone20[0]/el_probe_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_ELPROBEEVT_T), "" )
  mytree.Project("probelep_FAKE_ELPROBEEVT_L_ptvarcone20_over_pt", "el_probe_ptvarcone20[0]/el_probe_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_ELPROBEEVT_L), "" )
  mytree.Project("probelep_FAKE_ELPROBEEVT_T_topoetcone20_over_pt", "el_probe_topoetcone20[0]/el_probe_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_ELPROBEEVT_T), "" )
  mytree.Project("probelep_FAKE_ELPROBEEVT_L_topoetcone20_over_pt", "el_probe_topoetcone20[0]/el_probe_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_ELPROBEEVT_L), "" )

  if args.channel == "OF":
     # probe is electron, tag is muon
     #
     mytree.Project("taglep_FAKE_ELPROBEEVT_T_pt", "muon_tag_pt[0]/1e3", "%s * (%s)" %(mytree_weight, sel_FAKE_ELPROBEEVT_T), "" )
     mytree.Project("taglep_FAKE_ELPROBEEVT_T_ptvarcone30_over_pt", "muon_tag_ptvarcone30[0]/muon_tag_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_ELPROBEEVT_T), "" )
     mytree.Project("taglep_FAKE_ELPROBEEVT_T_ptvarcone20_over_pt", "muon_tag_ptvarcone20[0]/muon_tag_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_ELPROBEEVT_T), "" )
     mytree.Project("taglep_FAKE_ELPROBEEVT_T_topoetcone20_over_pt", "muon_tag_topoetcone20[0]/muon_tag_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_ELPROBEEVT_T), "" )
  elif args.channel == "SF":
     # probe is electron, tag is electron
     #
     mytree.Project("taglep_FAKE_ELPROBEEVT_T_pt", "el_tag_pt[0]/1e3", "%s * (%s)" %(mytree_weight, sel_FAKE_ELPROBEEVT_T), "" )
#     mytree.Project("taglep_FAKE_ELPROBEEVT_T_ptvarcone30_over_pt", "el_tag_ptvarcone30[0]/el_tag_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_ELPROBEEVT_T), "" )
     mytree.Project("taglep_FAKE_ELPROBEEVT_T_ptvarcone20_over_pt", "el_tag_ptvarcone20[0]/el_tag_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_ELPROBEEVT_T), "" )
     mytree.Project("taglep_FAKE_ELPROBEEVT_T_topoetcone20_over_pt", "el_tag_topoetcone20[0]/el_tag_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_ELPROBEEVT_T), "" )


  # Normalise to unity
  if args.doNorm:
     for hist in hist_list_probe_ELPROBEEVT.keys():
	if hist.Integral() == 0:
	  continue
        hist.Scale(1/hist.Integral())
	hist.GetYaxis().SetTitle("Norm. events")
	#hist.GetYaxis().SetRangeUser(0,1)
     for hist in hist_list_tag_ELPROBEEVT.keys():
	if hist.Integral() == 0:
	  continue
        hist.Scale(1/hist.Integral())
	hist.GetYaxis().SetTitle("Norm. events")
	#hist.GetYaxis().SetRangeUser(0,1)

  # now plot
  #
  print "\tLooking at events where the PROBE is: {0} and the channel is {1}\n".format("ELECTRON", args.channel)
  for  hist_tag in hist_list_tag_ELPROBEEVT.keys():

     c = TCanvas("c","Temp",50,50,700,900)
     c.Divide(1,2)

     c.cd(1)
     legend_T = TLegend(0.6,0.8,0.8,0.85); # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
     legend_T.SetHeader("2 Lep SS Fake CR - T probe")
     legend_T.SetBorderSize(0)  # no border
     legend_T.SetFillColor(0) # Legend background should be white
     legend_T.SetTextSize(0.04) # Increase entry font size!
     legend_T.SetTextFont(42) # Helvetica

     c.cd(2)
     legend_L = TLegend(0.6,0.8,0.8,0.85); # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
     legend_L.SetHeader("2 Lep SS Fake CR - L!T probe")
     legend_L.SetBorderSize(0)  # no border
     legend_L.SetFillColor(0) # legend_L background should be white
     legend_L.SetTextSize(0.04) # Increase entry font size!
     legend_L.SetTextFont(42) # Helvetica

     var_tag    = hist_list_tag_ELPROBEEVT[hist_tag][1]
     type_tag   = hist_list_tag_ELPROBEEVT[hist_tag][0]

     hist_tag.SetLineColor(1)
     hist_tag.SetLineWidth(2)

     print "\tPlotting variable: {0}\n".format(var_tag)

     print "\ttag histogram name: {0}\n".format(hist_tag.GetName())

     draw_pad1 = False
     draw_pad2 = False

     for hist_probe in hist_list_probe_ELPROBEEVT.keys():

       print "\t\tprobe histogram name: {0}".format(hist_probe.GetName())

       hist_probe.SetLineColor(2)
       hist_probe.SetLineWidth(2)

       var_probe  = hist_list_probe_ELPROBEEVT[hist_probe][1]
       type_probe = hist_list_probe_ELPROBEEVT[hist_probe][0]

       # must be looking at the same variable!
       #
       if  var_probe != var_tag:
         continue

       print "\t\tfound a match in variable! {0} = {1}\n".format(var_tag,var_probe)

       if type_probe == "T":
          print "\t\t  probe is of type: {0}".format(type_probe)
          print "\n\t\t  --> drawing: \t {0} \t {1} \t on pad(1) of the same canvas\n".format(hist_tag.GetName(),hist_probe.GetName())
          c.cd(1)
	  legend_T.AddEntry(None, "", "")
          legend_T.AddEntry(hist_tag, "REAL lepton", "L")
	  legend_T.AddEntry(None, "", "")
          legend_T.AddEntry(hist_probe, "FAKE lepton", "L")
          hist_tag.Draw("HIST")
          hist_probe.Draw("HIST SAME")
          draw_pad1 = True

       elif type_probe == "L":
          print "\t\t  probe is of type: {0}".format(type_probe)
          print "\n\t\t  --> drawing: \t {0} \t {1} \t on pad(2) of the same canvas\n".format(hist_tag.GetName(),hist_probe.GetName())
          c.cd(2)
	  legend_L.AddEntry(None, "", "")
	  legend_L.AddEntry(hist_tag, "REAL lepton", "L")
	  legend_L.AddEntry(None, "", "")
          legend_L.AddEntry(hist_probe, "FAKE lepton", "L")
          hist_tag.Draw("HIST")
          hist_probe.Draw("HIST SAME")
          draw_pad2 = True

       if ( draw_pad1 and draw_pad2 ):
	 c.cd(1)
         legend_T.Draw()
	 c.cd(2)
	 legend_L.Draw()
         plotname = dirname + "/FAKE_ELPROBEEVT_" + var_probe
         c.SaveAs(plotname+".png")

     del c

# -----------------------------------------------------------------------------------

def plot_2D():

  myfile = TFile(args.inputFile)

  mytree = gDirectory.Get("physics")
  mytree_weight = mytree.GetWeight()/1e3

  trig_dec  = "( passHLT == 1 )"
  lep_tag_trigmatch = "( lep_tag_isTrigMatched[0] == 1 && ( ( lep_tag_flavour[0] == 11 && lep_tag_pt[0] > 25e3 ) || ( lep_tag_flavour[0] == 13 && lep_tag_pt[0] > 22e3 ) ) )"
  lep_probe_trigmatch = "( lep_probe_isTrigMatched[0] == 1 && ( ( lep_probe_flavour[0] == 11 && lep_probe_pt[0] > 25e3 ) || ( lep_probe_flavour[0] == 13 && lep_probe_pt[0] > 22e3 ) ) )"
  el_tag_eta   = "( TMath::Abs(el_tag_eta[0]) < 1.37 )"
  nbjets       = "( njets_mv2c20_Fix77 > 0 )"
  njets        = "( njets > 0 && njets < 4 )"
  nleptons     = "( nlep == 2 && ( lep_pt[0] > 20e3 && lep_pt[1] > 20e3 ) )"
  same_sign    = "( isSS01 == 1 )"
  # veto charge flips
  #
  ch_flip_veto = "( 1 )"
  if args.doChFlipVeto:
     ch_flip_veto = "( lep_isChFlip[0] == 0 && lep_isChFlip[1] == 0 )"
  # require at least 1 !prompt lepton
  #
  non_prompt   = "( ( lep_truthType[0] != 6 && lep_truthType[0] != 2 ) || ( lep_truthType[1] != 6 && lep_truthType[1] != 2 ) )"
  # require at least one T lepton in the event
  #
  is_tight_event = "1" #"( isNonTightEvent == 0 )"
  # require tag lepton to be T
  #
  tight_tag   = "( lep_tag_isTightSelected[0] == 1 )"

  hist_list_tag_MUPROBEEVT = {}
  hist_list_tag_ELPROBEEVT = {}

  set_fancy_2D_style()

  gStyle.SetPaintTextFormat("2.1f")

  # -----------------------------------------------------

  dirname = "OutputPlots_TRUTH"
  if args.outdirname:
     dirname += ( "_" + args.outdirname )

  try:
      os.makedirs(dirname)
  except:
      pass

  # -----------------------------------------------------

  # ----------------
  #    Muon probe
  # ----------------

  sel_FAKE_MUPROBEEVT_T = None

  # look only at OF region
  if args.channel == "OF":
    sel_FAKE_MUPROBEEVT_T = "(" + trig_dec + " && " + lep_tag_trigmatch + " && " + lep_probe_trigmatch + " && " + nleptons + " && " + same_sign + " && " +  is_tight_event + " && " + tight_tag + " && " + "( (nmuon == 1 && nel == 1) && ( isProbeMuEvent == 1 ) && ( muon_probe_isTightSelected[0] == 1 ) )" + " && " + el_tag_eta + " && " + nbjets + " && " + njets  + " && " + ch_flip_veto + " && " + non_prompt + ")"
  # look only at SF region
  elif args.channel == "SF":
    sel_FAKE_MUPROBEEVT_T = "(" + trig_dec + " && " + lep_tag_trigmatch + " && " + lep_probe_trigmatch + " && " + nleptons + " && " + same_sign + " && " +  is_tight_event + " && " + tight_tag + " && " + "( ( nmuon == 2 ) && ( isProbeMuEvent == 1 ) && ( muon_probe_isTightSelected[0] == 1 ) )" + " && " + nbjets + " && " + njets + " && " + ch_flip_veto + " && " + non_prompt + ")"

  print "sel_FAKE_MUPROBEEVT_T: \n", sel_FAKE_MUPROBEEVT_T

  # histograms for tag (passing T: default)
  #
  h_taglep_FAKE_MUPROBEEVT_T_pt = TH2D("taglep_FAKE_MUPROBEEVT_T_pt", "taglep_FAKE_MUPROBEEVT_T_pt", 30, 20.0, 200.0, 30, 20.0, 200.0)
  hist_list_tag_MUPROBEEVT[h_taglep_FAKE_MUPROBEEVT_T_pt] = ["T","pt"]
  h_taglep_FAKE_MUPROBEEVT_T_ptvarcone30_over_pt = TH2D("taglep_FAKE_MUPROBEEVT_T_ptvarcone30_over_pt", "taglep_FAKE_MUPROBEEVT_T_ptvarcone30_over_pt", 20, 0.0, 0.1, 20, 0.0, 0.1)
  hist_list_tag_MUPROBEEVT[h_taglep_FAKE_MUPROBEEVT_T_ptvarcone30_over_pt] = ["T","ptvarcone30_over_pt"]
  h_taglep_FAKE_MUPROBEEVT_T_ptvarcone20_over_pt = TH2D("taglep_FAKE_MUPROBEEVT_T_ptvarcone20_over_pt", "taglep_FAKE_MUPROBEEVT_T_ptvarcone20_over_pt", 20, 0.0, 0.1, 20, 0.0, 0.1)
  hist_list_tag_MUPROBEEVT[h_taglep_FAKE_MUPROBEEVT_T_ptvarcone20_over_pt] = ["T","ptvarcone20_over_pt"]
  h_taglep_FAKE_MUPROBEEVT_T_topoetcone20_over_pt = TH2D("taglep_FAKE_MUPROBEEVT_T_topoetcone20_over_pt", "taglep_FAKE_MUPROBEEVT_T_topoetcone20_over_pt", 20, -0.1, 0.1, 20, -0.1, 0.1)
  hist_list_tag_MUPROBEEVT[h_taglep_FAKE_MUPROBEEVT_T_topoetcone20_over_pt] = ["T","topoetcone20_over_pt"]

  h_taglep_FAKE_MUPROBEEVT_T_pt.GetXaxis().SetTitle("REAL pT [GeV]")
  h_taglep_FAKE_MUPROBEEVT_T_ptvarcone30_over_pt.GetXaxis().SetTitle("REAL ptvarcone30/pT")
  h_taglep_FAKE_MUPROBEEVT_T_topoetcone20_over_pt.GetXaxis().SetTitle("REAL topoetcone20/pT")
  h_taglep_FAKE_MUPROBEEVT_T_pt.GetYaxis().SetTitle("FAKE pT [GeV]")
  h_taglep_FAKE_MUPROBEEVT_T_ptvarcone30_over_pt.GetYaxis().SetTitle("FAKE ptvarcone30/pT")
  h_taglep_FAKE_MUPROBEEVT_T_topoetcone20_over_pt.GetYaxis().SetTitle("FAKE topoetcone20/pT")

  if args.channel == "OF":
     # probe is muon, tag is electron
     #
     mytree.Project("taglep_FAKE_MUPROBEEVT_T_pt", "muon_probe_pt[0]/1e3:el_tag_pt[0]/1e3", "%s * (%s)" %(mytree_weight, sel_FAKE_MUPROBEEVT_T), "" )
     mytree.Project("taglep_FAKE_MUPROBEEVT_T_ptvarcone30_over_pt", "muon_probe_ptvarcone30[0]/muon_probe_pt[0]:el_tag_ptvarcone30[0]/el_tag_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_MUPROBEEVT_T), "" )
     mytree.Project("taglep_FAKE_MUPROBEEVT_T_ptvarcone20_over_pt", "muon_probe_ptvarcone20[0]/muon_probe_pt[0]:el_tag_ptvarcone20[0]/el_tag_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_MUPROBEEVT_T), "" )
     mytree.Project("taglep_FAKE_MUPROBEEVT_T_topoetcone20_over_pt", "muon_probe_topoetcone20[0]/muon_probe_pt[0]:el_tag_topoetcone20[0]/el_tag_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_MUPROBEEVT_T), "" )
  elif args.channel == "SF":
     # probe is muon, tag is muon
     #
     mytree.Project("taglep_FAKE_MUPROBEEVT_T_pt", "muon_probe_pt[0]/1e3:muon_tag_pt[0]/1e3", "%s * (%s)" %(mytree_weight, sel_FAKE_MUPROBEEVT_T), "" )
     mytree.Project("taglep_FAKE_MUPROBEEVT_T_ptvarcone30_over_pt", "muon_probe_ptvarcone30[0]/muon_probe_pt[0]:muon_tag_ptvarcone30[0]/muon_tag_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_MUPROBEEVT_T), "" )
#     mytree.Project("taglep_FAKE_MUPROBEEVT_T_ptvarcone20_over_pt", "muon_probe_ptvarcone20[0]/muon_probe_pt[0]:muon_tag_ptvarcone20[0]/muon_tag_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_MUPROBEEVT_T), "" )
     mytree.Project("taglep_FAKE_MUPROBEEVT_T_topoetcone20_over_pt", "muon_probe_topoetcone20[0]/muon_probe_pt[0]:muon_tag_topoetcone20[0]/muon_tag_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_MUPROBEEVT_T), "" )


  # Normalise to unity
  if args.doNorm:
     for hist in hist_list_tag_MUPROBEEVT.keys():
	if hist.Integral() == 0:
	  continue
        hist.Scale(1/hist.Integral())

  # now plot
  #
  print "\tLooking at events where the PROBE is: {0} and the channel is {1}\n".format("MUON", args.channel)
  for  hist_tag in hist_list_tag_MUPROBEEVT.keys():

     c = TCanvas("c","Temp",50,50,900,900)

     legend_T = TLegend(0.4,0.8,0.6,0.85); # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
     legend_T.SetHeader("2 Lep SS Fake CR - T probe")
     legend_T.SetBorderSize(0)  # no border
     legend_T.SetFillColor(0) # Legend background should be white
     legend_T.SetTextSize(0.04) # Increase entry font size!
     legend_T.SetTextFont(42) # Helvetica

     var_tag    = hist_list_tag_MUPROBEEVT[hist_tag][1]
     type_tag   = hist_list_tag_MUPROBEEVT[hist_tag][0]

     print "\tPlotting variable: {0}\n".format(var_tag)

     print "\ttag histogram name: {0}\n".format(hist_tag.GetName())

     hist_tag.Draw("colz")

     legend_T.Draw()
     plotname = dirname + "/FAKE_MUPROBEEVT_" + var_tag
     c.SaveAs(plotname+".png")

     del c

  # -----------------------------------------------------

  # ---------------------
  #    Electron probe
  # ---------------------

  sel_FAKE_ELPROBEEVT_T = None

  # look only at OF region
  if args.channel == "OF":
    sel_FAKE_ELPROBEEVT_T = "(" + trig_dec + " && " + lep_tag_trigmatch + " && " + lep_probe_trigmatch + " && " + nleptons + " && " + same_sign + " && " +  is_tight_event + " && " + tight_tag + " && " + "( (nmuon == 1 && nel == 1) && ( isProbeElEvent == 1 ) && ( el_probe_isTightSelected[0] == 1 ) )" + " && " + nbjets + " && " + njets + " && "  + ch_flip_veto + " && " + non_prompt + ")"
  # look only at SF region
  elif args.channel == "SF":
    sel_FAKE_ELPROBEEVT_T = "(" + trig_dec + " && " + lep_tag_trigmatch + " && " + lep_probe_trigmatch + " && " + nleptons + " && " + same_sign + " && " +  is_tight_event + " && " + tight_tag + " && " + "( ( nel == 2 ) && ( isProbeElEvent == 1 ) && ( el_probe_isTightSelected[0] == 1 ) )" + " && " + el_tag_eta + " && " + nbjets + " && " + njets + " && "  + ch_flip_veto + " && " + non_prompt + ")"

  print "sel_FAKE_ELPROBEEVT_T: \n", sel_FAKE_ELPROBEEVT_T

  # histograms for tag (passing T: default)
  #
  h_taglep_FAKE_ELPROBEEVT_T_pt = TH2D("taglep_FAKE_ELPROBEEVT_T_pt", "taglep_FAKE_ELPROBEEVT_T_pt", 30, 20.0, 200.0, 30, 20.0, 200.0)
  hist_list_tag_ELPROBEEVT[h_taglep_FAKE_ELPROBEEVT_T_pt] = ["T","pt"]
  h_taglep_FAKE_ELPROBEEVT_T_ptvarcone30_over_pt = TH2D("taglep_FAKE_ELPROBEEVT_T_ptvarcone30_over_pt", "taglep_FAKE_ELPROBEEVT_T_ptvarcone30_over_pt", 20, 0.0, 0.1, 20, 0.0, 0.1)
  hist_list_tag_ELPROBEEVT[h_taglep_FAKE_ELPROBEEVT_T_ptvarcone30_over_pt] = ["T","ptvarcone30_over_pt"]
  h_taglep_FAKE_ELPROBEEVT_T_ptvarcone20_over_pt = TH2D("taglep_FAKE_ELPROBEEVT_T_ptvarcone20_over_pt", "taglep_FAKE_ELPROBEEVT_T_ptvarcone20_over_pt", 20, 0.0, 0.1, 20, 0.0, 0.1)
  hist_list_tag_ELPROBEEVT[h_taglep_FAKE_ELPROBEEVT_T_ptvarcone20_over_pt] = ["T","ptvarcone20_over_pt"]
  h_taglep_FAKE_ELPROBEEVT_T_topoetcone20_over_pt = TH2D("taglep_FAKE_ELPROBEEVT_T_topoetcone20_over_pt", "taglep_FAKE_ELPROBEEVT_T_topoetcone20_over_pt", 20, -0.1, 0.1, 20, -0.1, 0.1)
  hist_list_tag_ELPROBEEVT[h_taglep_FAKE_ELPROBEEVT_T_topoetcone20_over_pt] = ["T","topoetcone20_over_pt"]

  h_taglep_FAKE_ELPROBEEVT_T_pt.GetXaxis().SetTitle("REAL pT [GeV]")
  h_taglep_FAKE_ELPROBEEVT_T_ptvarcone20_over_pt.GetXaxis().SetTitle("REAL ptvarcone20/pT")
  h_taglep_FAKE_ELPROBEEVT_T_topoetcone20_over_pt.GetXaxis().SetTitle("REAL topoetcone20/pT")
  h_taglep_FAKE_ELPROBEEVT_T_pt.GetYaxis().SetTitle("FAKE pT [GeV]")
  h_taglep_FAKE_ELPROBEEVT_T_ptvarcone20_over_pt.GetYaxis().SetTitle("FAKE ptvarcone20/pT")
  h_taglep_FAKE_ELPROBEEVT_T_topoetcone20_over_pt.GetYaxis().SetTitle("FAKE topoetcone20/pT")

  if args.channel == "OF":
     # probe is el, tag is muon
     #
     mytree.Project("taglep_FAKE_ELPROBEEVT_T_pt", "el_probe_pt[0]/1e3:muon_tag_pt[0]/1e3", "%s * (%s)" %(mytree_weight, sel_FAKE_ELPROBEEVT_T), "" )
     mytree.Project("taglep_FAKE_ELPROBEEVT_T_ptvarcone30_over_pt", "el_probe_ptvarcone30[0]/el_probe_pt[0]:muon_tag_ptvarcone30[0]/muon_tag_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_ELPROBEEVT_T), "" )
     mytree.Project("taglep_FAKE_ELPROBEEVT_T_ptvarcone20_over_pt", "el_probe_ptvarcone20[0]/el_probe_pt[0]:muon_tag_ptvarcone20[0]/muon_tag_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_ELPROBEEVT_T), "" )
     mytree.Project("taglep_FAKE_ELPROBEEVT_T_topoetcone20_over_pt", "el_probe_topoetcone20[0]/el_probe_pt[0]:muon_tag_topoetcone20[0]/muon_tag_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_ELPROBEEVT_T), "" )
  elif args.channel == "SF":
     # probe is el, tag is el
     #
     mytree.Project("taglep_FAKE_ELPROBEEVT_T_pt", "el_probe_pt[0]/1e3:el_tag_pt[0]/1e3", "%s * (%s)" %(mytree_weight, sel_FAKE_ELPROBEEVT_T), "" )
     #mytree.Project("taglep_FAKE_ELPROBEEVT_T_ptvarcone30_over_pt", "el_probe_ptvarcone30[0]/el_probe_pt[0]:el_tag_ptvarcone30[0]/el_tag_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_ELPROBEEVT_T), "" )
     mytree.Project("taglep_FAKE_ELPROBEEVT_T_ptvarcone20_over_pt", "el_probe_ptvarcone20[0]/el_probe_pt[0]:el_tag_ptvarcone20[0]/el_tag_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_ELPROBEEVT_T), "" )
     mytree.Project("taglep_FAKE_ELPROBEEVT_T_topoetcone20_over_pt", "el_probe_topoetcone20[0]/el_probe_pt[0]:el_tag_topoetcone20[0]/el_tag_pt[0]", "%s * (%s)" %(mytree_weight, sel_FAKE_ELPROBEEVT_T), "" )


  # Normalise to unity
  if args.doNorm:
     for hist in hist_list_tag_ELPROBEEVT.keys():
	if hist.Integral() == 0:
	  continue
        hist.Scale(1/hist.Integral())

  # now plot
  #
  print "\tLooking at events where the PROBE is: {0} and the channel is {1}\n".format("ELECTRON", args.channel)
  for  hist_tag in hist_list_tag_ELPROBEEVT.keys():

     c = TCanvas("c","Temp",50,50,900,900)

     legend_T = TLegend(0.4,0.8,0.6,0.85); # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
     legend_T.SetHeader("2 Lep SS Fake CR - T probe")
     legend_T.SetBorderSize(0)  # no border
     legend_T.SetFillColor(0) # Legend background should be white
     legend_T.SetTextSize(0.04) # Increase entry font size!
     legend_T.SetTextFont(42) # Helvetica

     var_tag    = hist_list_tag_ELPROBEEVT[hist_tag][1]
     type_tag   = hist_list_tag_ELPROBEEVT[hist_tag][0]

     print "\tPlotting variable: {0}\n".format(var_tag)

     print "\ttag histogram name: {0}\n".format(hist_tag.GetName())

     hist_tag.Draw("colz")

     legend_T.Draw()
     plotname = dirname + "/FAKE_ELPROBEEVT_" + var_tag
     c.SaveAs(plotname+".png")

     del c

# -----------------------------------------------------------------------------------

#plot()

plot_2D()





