#!/usr/bin/python

import array
import os
import sys

sys.path.append(os.path.abspath(os.path.curdir))

# -------------------------------
# Parser for command line options
# -------------------------------
import argparse

parser = argparse.ArgumentParser(description="Plotting python macro for deriving real/fake lepton efficiencies/rates for MM.")

#***********************************
# positional arguments (compulsory!)
#***********************************
parser.add_argument("inputDir", metavar="inputDir",type=str,
		  help="path to the directory containing subdirs w/ input files")
#*******************
# optional arguments
#*******************
parser.add_argument("--flavourComp", metavar="FLAVOUR_COMP", dest="flavourComp", default="", type=str,
		  help="Flavour composition of the two leptons in CR (*empty_string*, ElEl, MuMu, OF) - default is *empty_string*")
parser.add_argument("--usePrediction", metavar="DATA_TYPE", dest="usePrediction", action="store", default="DATA", type=str,
		    help="use Monte-Carlo (MC) or data (DATA) to derive efficiencies/rates - default is DATA")
parser.add_argument("--debug", dest="debug", action="store_true",
		    help="run in debug mode")
parser.add_argument("--doBkgSub", dest="doBkgSub", action="store_true",
		    help="subtract prompt MC in fake region (only if using DATA)")
parser.add_argument("--rebinEta", dest="rebinEta", action="store_true",default=False,
		    help="do rebinning in eta")
parser.add_argument("--rebinPt", dest="rebinPt", action="store_true",default=False,
		    help="do rebinning in pT")
parser.add_argument("--doAvg", dest="doAvg", action="store_true",default=False,
		    help="get average efficiencies/rates (i.e, make 1 bin)")
parser.add_argument("--saveOnlyRates", dest="saveOnlyRates", action="store_true",
		    help="save only rates")
args = parser.parse_args()

if not args.inputDir.endswith("/"):
   args.inputDir += "/"

# ---------------------------------------------------------------------------------
# NB: the following string definitions are to be chosen according to the
#     output name of the files (produced by the MakePlots_HTopMultilep.py scripts)
#     used for comuting the T/!T(=L) ratio to derive rates.
#
#     An example name could be:
#
#        MuElFakeCRMuL_ElProbePt
#
#        channels[h] + lep_types[i] + CR + leptons[j] + selections[k] + leptons[l] + variables[m]
# ---------------------------------------------------------------------------------

# -------------------------------------------------------
# for each channel, store which leptons can be associated
# -------------------------------------------------------
dict_channels_lep = {
 	             "ElEl": ["El"],
                     "OF"  : ["El","Mu"],
                     ""    : ["El","Mu"],
                     "MuMu": ["Mu"]
		    }

list_lep         = dict_channels_lep[args.flavourComp]
list_types       = ["Fake","Real"]
list_variables   = ["ProbeEta","ProbePt"] #,"ProbeNJets"]
list_selections  = ["T","L"]
list_prediction  = ["expected", "observed"]   # expected --> use MC distribution for probe lepton to derive the rate (to be used only as a cross check, and in closure test)
                                              # observed --> use DATA distribution for probe lepton to derive the rate - need to subtract the prompt/ch-flips here!
list_out_samples = ["factor","factorbkgsub","rate","ratebkgsub"]

# ----------------------------------------------------------------
# if you want to subtract backgrounds to estimate the fake rates
#
# ( prompt subtraction, charge flip subtraction )
# ----------------------------------------------------------------

from ROOT import gROOT, TH1, TFile, TGraphAsymmErrors, Double

gROOT.SetBatch(True)

TH1.SetDefaultSumw2()

def RateToEfficiency( v ):
    v = float(v)
    if v < 0:
        v = 0.
    return v/(v+1)

hists  = {}
graphs = {}
yields = {}
fin    = []

channel_str = ""
if ( args.flavourComp is not "" ):
   #channel_str = "_" + args.flavourComp
   channel_str = args.flavourComp
   print " channel string: ", channel_str

for iLep in list_lep:

   print "looking at lepton of flavour: " , iLep

   for iVar in list_variables:

      print "\t looking at variable: " , iVar

      for iType in list_types:

          print "\t\t looking at rate of type: ", iType

          for iSel in list_selections:

              print "\t\t\t object  is: ", iSel

              fname = args.inputDir + channel_str + iType + "CR" + iLep + iSel + "/" + channel_str + iType + "CR" + iLep + iSel + "_" + iLep + iVar + ".root"

	      print "\t\t\t input filename: ", fname

              fin.append( TFile(fname) )

	      for iPred in list_prediction:

                  print "\t\t\t\t checking prediction from: ", iPred

                  # L[-1] can be used to access the last item in a list
          	  htmp = fin[-1].Get( iPred )

		  # let the macro know the name of the histogram in the input ROOT file
		  #
		  # do not even bother for "observed" hist if usePrediction = MC
		  #
		  if ( args.usePrediction == "MC" ) and ( iPred == "observed" ) :
		     continue

	  	  if not htmp:
		     sys.exit("\t\t\t\t\t error: histogram htmp does not exist")

                  histname = iLep + "_" + iVar + "_" + iType + "_" + iSel + "_" + iPred

                  print "\t\t\t\t\t output histname (Num,Denom): ", histname

		  # make a clone of this histogram by default
		  #
	          hists[histname]  = htmp.Clone( histname )

                  if iVar == "ProbeEta":

		     # rebin in eta ONLY for rates in data
		     if ( args.usePrediction == "DATA" ) and args.rebinEta:

		       if iLep == "Mu":
                  	 nBIN  = 5
                         xbins = [ 0.0 , 0.1 , 0.7, 1.3 , 1.9, 2.5 ]
                  	 # the rebinning method automatically creates a new histogram
                  	 #
		     	 vxbins = array.array("d", xbins)
		     	 print "\t\t\t\t\t vxbins: ",vxbins
                  	 hists[histname]  = htmp.Rebin( nBIN, histname, vxbins )

		       elif iLep == "El":

                  	 nBIN  = 6
                         xbins = [ 0.0 , 0.5 , 0.8 , 1.37 , 1.52 , 2.0 , 2.6]
                  	 # the rebinning method automatically creates a new histogram
                  	 #
		     	 vxbins = array.array("d", xbins)
		     	 print "\t\t\t\t\t vxbins: ",vxbins
                  	 hists[histname]  = htmp.Rebin( nBIN, histname, vxbins )

		  elif iVar == "ProbePt":

		     if args.rebinPt:

		       if iType == "Fake":

                          if iLep == "Mu":
                             nBIN  = 3
                             xbins = [25,35,50,200]
			  elif iLep == "El":
                             nBIN  = 3
                             xbins = [25,40,60,200]

                       elif iType == "Real":

                  	  #nBIN  = 12
                  	  #xbins = [25,30,35,40,45,50,60,70,80,90,100,150,200]
			  nBIN  = 4
                  	  xbins = [25,30,40,60,200]

                       vxbins = array.array("d", xbins)
		       print "\t\t\t\t\t vxbins: ",vxbins
                       # the rebinning method automatically creates a new histogram
                       #
                       hists[histname]  = htmp.Rebin( nBIN, histname, vxbins )

                  if args.doAvg:

                     nBIN = 1
                     if iVar == "ProbePt":
                        xbins = [25,200]
                     elif iVar == "ProbeEta":
                        if iLep == "El":
                           xbins = [0.0,2.6]
                        elif iLep == "Mu":
                           xbins == [0.0,2.5]
                     elif iVar == "ProbeNJets":
                        xbins = [2,10]
                     vxbins = array.array("d", xbins)
                     print "\t\t\t\t\t vxbins: ",vxbins
                     # the rebinning method automatically creates a new histogram
                     #
                     hists[histname]  = htmp.Rebin( nBIN, histname, vxbins )

		  yields[histname] = hists[histname].Integral()

	      # -------------------------------------------------
              # subtract prompt bkg from data in fakes histograms
	      # -------------------------------------------------

	      if args.doBkgSub :

		 if ( args.usePrediction == "MC" ):
		     sys.exit("trying to subtract MC when using --usePrediction=", args.usePrediction ," option. Please use DATA instead, or switch off prompt background subtraction")

		 name = iLep + "_"+ iVar +"_"+ iType + "_" +  iSel + "_" + list_prediction[1]

		 if ( iType == "Fake" ):

		     print "\t\t\t\t subtracting prompt/ch-flip MC to data in Fake CR..."

		     hists[name].Add( hists[ iLep + "_" + iVar + "_" + iType + "_" + iSel + "_" + list_prediction[0] ], -1 )

		 # ------------------------------------------------------------
		 # set bin content to zero if subtraction gives negative yields
		 #
		 # compute the yield for this histogram
		 # ------------------------------------------------------------

		 for iBin in range( 0, hists[name].GetNbinsX()+2):

		     if hists[name].GetBinContent( iBin ) < 0:
                        hists[name].SetBinContent( iBin, 0 )

                 yields[name] = hists[name].Integral()

	  # -----------------
	  # compute the rate:
	  # T / !T
	  #
	  # and the efficiency:
	  # T / ( !T + T )
	  # -----------------

	  histname =  iLep + "_" + iVar +"_"+ iType

          print "\t\t histogram to be used in the ratio: ", histname

          # ---------------------------------------------------
	  # use MC based estimate, but only if selected by user
          # ---------------------------------------------------

	  append_str = "observed"
          if ( args.usePrediction == "MC" ):
	      append_str = "expected"

          hists[histname + "_Rate_" + append_str]  = hists[histname + "_T_" + append_str].Clone(histname + "_Rate_" + append_str)
          hists[histname + "_Rate_" + append_str].Divide(hists[histname+"_L_" + append_str])
          yields[histname + "_Rate_" + append_str] = yields[histname + "_T_" + append_str] / yields[histname + "_L_" + append_str]

          # For efficiency, make sure the errors are computed correctly
          # (NB: in this case, numerator and denominator are not independent sets of events!)
          #
          hist_pass = hists[histname + "_T_" + append_str]
          hist_tot  = hists[histname + "_T_" + append_str] + hists[histname+"_L_" + append_str]
          hist_eff  = hist_pass.Clone(histname + "_Efficiency_" + append_str)
          hist_eff.Divide(hist_pass, hist_tot,1.0,1.0,"B")
          hists[histname + "_Efficiency_" + append_str] = hist_eff
          #
          # even better, use TGraphAsymmErrors
          # (handles the errors in case eff is 100%)
          #
          g_efficiency = TGraphAsymmErrors(hist_eff)
          g_efficiency.Divide(hist_pass,hist_tot,"cl=0.683 b(1,1) mode")
          graphs[histname + "_Efficiency_" + append_str + "_graph"] = g_efficiency

          print "\t\t --> ratio hist name: ", histname + "_Rate_" + append_str

if args.doAvg:
   channel_str += "Avg"

outfile = open( args.inputDir + channel_str + "Rates.txt", "w")
outfile.write( "Efficiencies/Rates for FF amd Matrix Method \n")

foutname = args.inputDir + channel_str + "Rates.root"
fout = TFile( foutname,"RECREATE" )
fout.cd()

for g in sorted( graphs.keys() ):

   print "saving graph: ", graphs[g].GetName()

   graphs[g].Write()

   Eff=[]

   if (args.debug ) :
       print "saving efficiencies (and errors) from graph: ", graphs[g].GetName()

   for ipoint in range( 0, graphs[g].GetN() ):

       x = Double(0)
       y = Double(0)
       graphs[g].GetPoint(ipoint,x,y)
       set = [ ipoint, y, graphs[g].GetErrorYhigh(ipoint), graphs[g].GetErrorYlow(ipoint) ]
       Eff.append( set )

   outfile.write("%s \n" %(g) )
   for set in Eff:
       outfile.write("{ %s }; \n" %( "Bin nr: " + str(set[0]) + ", efficiency = " + str(round(set[1],3)) + " + " + str(round(set[2],3)) + " - " + str(round(set[3],3)) ) )


for h in sorted( hists.keys() ):

   print "saving histogram: ", hists[h].GetName()

   hists[h].Write()
   Rate=[]

   if "_Rate_" in h:

      if (args.debug ) :
   	 print "saving rates (and errors) for histogram: ", hists[h].GetName()

      Rtot=yields[h]

      for ibin in range( 1, hists[h].GetNbinsX()+1 ):

         set = [ ibin, hists[h].GetBinContent(ibin), hists[h].GetBinError(ibin)]
  	 Rate.append( set )

      outfile.write("%s \n" %(h) )
      for set in Rate:
          outfile.write("{ %s }; \n" %( "Bin nr: " + str(set[0]) + ", rate = " + str(round(set[1],3)) + " +- " + str(round(set[2],3)) ) )

outfile.close()
fout.Write()
fout.Close()

for f in range( 0,len(fin) ):
   fin[f].Close()
