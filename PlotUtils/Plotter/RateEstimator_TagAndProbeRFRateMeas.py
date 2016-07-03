#!/usr/bin/env python

""" RateEstimator_TagAndProbeRFRateMeas.py: measure r/f efficiencies/rates for Matrix Method and Fake Factor Method """

__author__     = "Marco Milesi, Francesco Nuti"
__email__      = "marco.milesi@cern.ch, francesco.nuti@cern.ch"
__maintainer__ = "Marco Milesi"

import os, sys, array, math

sys.path.append(os.path.abspath(os.path.curdir))

from ROOT import ROOT, gROOT, TH1, TH1D, TFile, TGraphAsymmErrors, Double

gROOT.SetBatch(True)
TH1.SetDefaultSumw2()

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
		  help="subtract MC to data (only if using DATA)")
parser.add_argument("--rebinEta", dest="rebinEta", action="store_true",default=False,
		  help="do rebinning in eta")
parser.add_argument("--rebinPt", dest="rebinPt", action="store_true",default=False,
		  help="do rebinning in pT")
parser.add_argument("--doAvgMuFake", dest="doAvgMuFake", action="store_true",default=False,
		  help="get average efficiency for muon fakes only")
parser.add_argument("--doAvg", dest="doAvg", action="store_true",default=False,
		  help="get average efficiencies/rates (i.e, make 1 bin)")
parser.add_argument("--saveOnlyRates", dest="saveOnlyRates", action="store_true",
		  help="save only rates")
parser.add_argument("--useLogPlots", dest="useLogPlots", action="store_true",default=False,
		  help="Read plots with logarithmic Y scale")

args = parser.parse_args()

def getQMisIDRatesRatio( N, D ):

  xbins_N = []
  for ibinx in range(1,N.GetNbinsX()+1):
    xbins_N.append(N.GetBinLowEdge(ibinx))

  xbins_D = []
  for ibinx in range(1,D.GetNbinsX()+1):
    xbins_D.append(D.GetBinLowEdge(ibinx))

  if args.debug:
    print "xbins_N = ", xbins_N
    print "xbins_D = ", xbins_D

  vxbins_N = array.array("d", xbins_N)
  vxbins_D = array.array("d", xbins_D)

  new_N = TH1D("N","N",len(vxbins_N)-1,vxbins_N)
  if args.debug: print("QMisID rate T - Numerator:")
  for ibinx in range(1,N.GetNbinsX()):
    content = N.GetBinContent(ibinx)
    error   = N.GetBinError(ibinx)
    perc_error = abs(error) * 100 / content
    if args.debug: print("bin {0} - {1} +- {2} ({3} %)".format(ibinx,content,error,perc_error))
    new_N.SetBinContent(ibinx,content)
    new_N.SetBinError(ibinx,error)

  new_D = TH1D("D","D",len(vxbins_D)-1,vxbins_D)
  if args.debug: print("QMisID rate L - Denominator:")
  for ibinx in range(1,D.GetNbinsX()):
    content = D.GetBinContent(ibinx)
    error   = D.GetBinError(ibinx)
    perc_error = abs(error) * 100 / content
    if args.debug: print("bin {0} - {1} +- {2} ({3} %)".format(ibinx,content,error,perc_error))
    new_D.SetBinContent(ibinx,content)
    new_D.SetBinError(ibinx,error)

  ratio = new_N.Clone()
  ratio.Divide(new_N,new_D,1.0,1.0,"B")
  return ratio


def scaleEff( hist_r, hist_f, hist_data, hist_QMisID ):

  file_N = TFile("$ROOTCOREBIN/data/HTopMultilepAnalysis/External/QMisID_Pt_rates_Tight_v15.root")
  file_D = TFile("$ROOTCOREBIN/data/HTopMultilepAnalysis/External/QMisID_Pt_rates_Loose_v15.root")

  hist_N = file_N.Get("LikelihoodPt")
  hist_D = file_D.Get("LikelihoodPt")

  hist_scales = getQMisIDRatesRatio( hist_N, hist_D)

  hist_eff_QMisID   = hist_r.Clone()
  hist_eff_f_w      = hist_r.Clone()

  if args.debug:
    print("Histogram w/ scale factors - GetNbinsX() = {0}".format(hist_scales.GetNbinsX()))
    print("Input REAL efficiency hist - GetNbinsX() = {0}".format(hist_r.GetNbinsX()))

  for ibinx in range(1,hist_r.GetNbinsX()+1):

    eff_r     = hist_r.GetBinContent(ibinx)
    err_eff_r = hist_r.GetBinError(ibinx)

    scale	   = hist_scales.GetBinContent(ibinx)
    err_scale	   = hist_scales.GetBinError(ibinx)
    perc_err_scale = abs(err_scale) * 100 / scale

    eff_QMisID = scale * eff_r
    err_eff_QMisID = math.sqrt( (eff_r*eff_r) * (err_scale*err_scale) + (scale*scale) * (err_eff_r*err_eff_r) )

    hist_eff_QMisID.SetBinContent(ibinx, eff_QMisID)
    hist_eff_QMisID.SetBinError(ibinx, err_eff_QMisID)

    # Now get the updated FAKE efficiency, by taking the weighted average <input eff fake, QMisID eff>
    # The weight is the QMisID fraction wrt the total data for each bin
    # (The sum of weight is correctly normalised to 1)
    #
    w = hist_QMisID.GetBinContent(ibinx) / hist_data.GetBinContent(ibinx)

    # Deal with the (ill) case where the QMisID yield is higher than data
    #
    if w > 1.0: w = 1.0

    eff_f     = hist_f.GetBinContent(ibinx)
    err_eff_f = hist_f.GetBinError(ibinx)

    eff_w   = ( 1.0 - w ) * eff_f + w * eff_QMisID
    err_eff_w = ( 1.0 - w + ( 1.0 / scale ) * w ) * err_eff_f

    # Add in quadrature the difference w/ nominal fake efficiency as a systematic
    #
    sys_err_eff_w = abs(eff_f - eff_w)
    tot_err_eff_w = math.sqrt( (err_eff_w)*(err_eff_w) + (sys_err_eff_w)*(sys_err_eff_w) )

    print("bin {0} - [{1},{2}] GeV \n\teff_r: {3} +- {4}\n\tscale : {5} +- {6} ({7} %)\n\teff_QMisID: {8} +- {9}\n\teff_f: {10} +- {11}\n\tweighted eff_f: {12} +- {13} ({14})".format(ibinx,hist_r.GetBinLowEdge(ibinx),hist_r.GetBinLowEdge(ibinx+1),eff_r,err_eff_r,scale,err_scale,perc_err_scale,eff_QMisID,err_eff_QMisID,eff_f,err_eff_f,eff_w,err_eff_w,tot_err_eff_w))

    # If efficiency unphysically > 1, force it to 1
    #
    if eff_w > 1.0 : eff_w = 1.0

    hist_eff_f_w.SetBinContent(ibinx, eff_w)
    hist_eff_f_w.SetBinError(ibinx, tot_err_eff_w)

  hist_eff_QMisID.SetDirectory(0)
  hist_eff_f_w.SetDirectory(0)

  return hist_eff_QMisID, hist_eff_f_w

def RateToEfficiency( v ):
    v = float(v)
    if v < 0:
	v = 0.
    return v/(v+1)

if __name__ == "__main__":

  if not args.inputDir.endswith("/"):
     args.inputDir += "/"

  # ---------------------------------------------------------------------------------
  # NB: the following string definitions are to be chosen according to the
  #     output name of the files (produced by the MakePlots_HTopMultilep.py scripts)
  #     used for comuting the T/!T(=L) ratio to derive rates.
  #
  #     An example name could be:
  #
  #	 MuElFakeCRMuL_ElProbePt
  #
  #	 channels[h] + lep_types[i] + CR + leptons[j] + selections[k] + leptons[l] + variables[m]
  # ---------------------------------------------------------------------------------

  if args.rebinEta or args.rebinPt or args.doAvg:
     print("\n******************************************************************************************************\n\nWill perform histogram rebinning...")
     print("\nWARNING!\n(From TH1 docs) If ngroup is not an exact divider of the number of bins, the top limit of the rebinned histogram is reduced to the upper edge of the last bin that can make a complete group.")
     print("The remaining bins are added to the overflow bin. Statistics will be recomputed from the new bin contents.\n")
     print("If rebinning to one single bin, this might lead to an \"empty\" histogram (as everything will end up in the overflow bin)")
     print("\n******************************************************************************************************\n")

  # -------------------------------------------------------
  # for each channel, store which leptons can be associated
  # -------------------------------------------------------
  dict_channels_lep = {
  		     "ElEl": ["El"],
  		     "OF"  : ["El","Mu"],
  		     ""    : ["El","Mu"],
  		     "MuMu": ["Mu"]
  		    }

  list_lep	   = dict_channels_lep[args.flavourComp]
  list_types	   = ["Real","RealQMisIDBinning","Fake"]
  list_variables   = ["ProbePt"] #,"ProbeEta"] #"ProbeNJets"]
  list_selections  = ["T","AntiT"]
  list_prediction  = ["expected", "observed"]   # expected --> use MC distribution for probe lepton to derive the rate (to be used only as a cross check, and in closure test)
  					      # observed --> use DATA distribution for probe lepton to derive the rate - need to subtract the prompt/ch-flips here!
  list_out_samples = ["factor","factorbkgsub","rate","ratebkgsub"]

  hists  = {}
  graphs = {}
  yields = {}
  fin    = []

  channel_str = ""
  if ( args.flavourComp is not "" ):
     #channel_str = "_" + args.flavourComp
     channel_str = args.flavourComp
     print " channel string: ", channel_str

  log_plot = ("","_LOGY")[bool(args.useLogPlots)]

  # For later use
  #
  hist_eff_electron_real = None

  for iLep in list_lep:

     print "looking at lepton of flavour: " , iLep

     for iVar in list_variables:

        print "\t looking at variable: " , iVar

        for iType in list_types:

          print "\t\t looking at rate of type: ", iType

          for iSel in list_selections:

              print "\t\t\t object  is: ", iSel

              iActualType = iType
              if iType == "RealQMisIDBinning": iActualType = "Real"

              fname = args.inputDir + channel_str + iActualType + "CR" + iLep + iSel + log_plot + "/" + channel_str + iActualType + "CR" + iLep + iSel + "_" + iLep + iVar + ".root"

  	      print "\t\t\t input filename: ", fname

              fin.append( TFile(fname) )

              htmp_QMisID = None
  	      if  ( args.usePrediction == "DATA" ) and ( iType == "Fake" ) and ( iLep == "El" ) and ( iVar == "ProbePt" ):
                htmp_QMisID = fin[-1].Get("chargeflipbkg")
  	        if not htmp_QMisID:
  	          sys.exit("ERROR: histogram w/ name chargeflipbkg does not exist in input file")

  	      for iPred in list_prediction:

          	  standard_rebin = False

          	  print "\t\t\t\t checking prediction from: ", iPred

          	  # L[-1] can be used to access the last item in a list
          	  htmp = fin[-1].Get( iPred ) # will be either "expected", or "observed"

  	          # let the macro know the name of the histogram in the input ROOT file
  	          #
  	          # do not even bother for "observed" hist if usePrediction = MC
  	          #
  	          if ( args.usePrediction == "MC" ) and ( iPred == "observed" ) :
  	             continue

  	          if not htmp:
  	             sys.exit("ERROR: histogram w/ name %s does not exist in input file", iPred)

          	  histname = iLep + "_" + iVar + "_" + iType + "_" + iSel + "_" + iPred

          	  print "\t\t\t\t\t output histname (Num,Denom): ", histname

  	          # make a clone of this histogram (expected/observed) by default
  	          #
  	  	  hists[histname]  = htmp.Clone( histname )

  	          # make a clone also of the QMisID histogram (do it only once)
  	          #
	          histname_QMisID = histname + "_QMisID"
  	          hists[histname_QMisID] = None
  	          if htmp_QMisID and histname_QMisID not in hists:
  	            hists[histname_QMisID]  = htmp_QMisID.Clone( histname_QMisID )

          	  if iVar == "ProbeEta":

  	            if args.rebinEta:

  	              if iLep == "Mu":
          		nBIN  = 5
          		xbins = [ 0.0 , 0.1 , 0.7, 1.3 , 1.9, 2.5 ]
          		# the rebinning method automatically creates a new histogram
          		#
  	        	vxbins = array.array("d", xbins)
  	        	print "\t\t\t\t\t vxbins: ",vxbins
          		hists[histname]  = htmp.Rebin( nBIN, histname, vxbins )
          		if htmp_QMisID : hists[histname_QMisID]  = htmp_QMisID.Rebin( nBIN, histname_QMisID, vxbins )

  	              elif iLep == "El":

          		nBIN  = 6
          		xbins = [ 0.0 , 0.5 , 0.8 , 1.37 , 1.52 , 2.0 , 2.6]
          		# the rebinning method automatically creates a new histogram
          		#
  	        	vxbins = array.array("d", xbins)
  	        	print "\t\t\t\t\t vxbins: ",vxbins
          		hists[histname]  = htmp.Rebin( nBIN, histname, vxbins )
          		if htmp_QMisID : hists[histname_QMisID]  = htmp_QMisID.Rebin( nBIN, histname_QMisID, vxbins )

  	          elif iVar == "ProbePt":

  	            if args.rebinPt:

  	              if iType == "Fake":

          		if iLep == "Mu":

          		   if args.doAvgMuFake:

  	        	     # Make one single bin
          		     #
          		     #nBIN = htmp.GetNbinsX()
          		     #xbins  = [htmp.GetBinLowEdge(1),htmp.GetBinLowEdge(htmp.GetNbinsX()+1)]
          		     #standard_rebin = True

          		     nBIN  = 4
          		     xbins = [10,15,20,25,200]

          		   else:
          		     # nominal binning
          		     #
          		     nBIN  = 5
          		     xbins = [10,15,20,25,35,200] # merged bin [50,200] - this binning works a bit better for ttbar MC at medium-high pT range

  	        	elif iLep == "El":

  	        	  # standard binning
  	        	  #
  	        	  nBIN  = 5
          		  xbins = [10,15,20,25,40,200] # merged bin [60,200]

          	      elif iType == "Real":

  	        	# standard binning
  	        	#
  	        	nBIN  = 7
          		xbins = [10,15,20,25,30,40,60,200]

                      elif iType == "RealQMisIDBinning":

  	        	# Need this binning to get the weighted fake eff from QMisID eff.
  	        	#
  	        	nBIN  = 5
          		xbins = [10,15,20,25,40,200]

          	      vxbins = array.array("d", xbins)
          	      print "\t\t\t\t\t vxbins: ",vxbins
          	      # the rebinning method automatically creates a new histogram
          	      #
          	      if not standard_rebin:
  	                hists[histname]	= htmp.Rebin( nBIN, histname, vxbins )
          		if htmp_QMisID: hists[histname_QMisID]  = htmp_QMisID.Rebin( nBIN, histname_QMisID, vxbins )
          	      else :
  	                hists[histname] = htmp.Rebin( nBIN, histname )
          		if htmp_QMisID : hists[histname_QMisID]  = htmp_QMisID.Rebin( nBIN, histname_QMisID )

          	  if args.doAvg:

  	            nBIN = htmp.GetNbinsX()
  	            xbins  = [htmp.GetBinLowEdge(1),htmp.GetBinLowEdge(htmp.GetNbinsX()+1)]
          	    vxbins = array.array("d", xbins)
          	    print "\t\t\t\t\t vxbins: ",vxbins
  	            hists[histname]  = htmp.Rebin( nBIN, histname )
          	    if htmp_QMisID : hists[histname_QMisID]  = htmp_QMisID.Rebin( nBIN, histname_QMisID )

  	          if ( args.rebinEta or args.rebinPt or args.doAvg) and args.debug:
  	             print("\t\t\t\t\t Integral BEFORE rebinning: {0}".format(htmp.Integral(0,htmp.GetNbinsX()+1)))
  	             print("\t\t\t\t\t Integral AFTER rebinning: {0}".format(hists[histname].Integral(0,hists[histname].GetNbinsX()+1)))

  	          # -------------------------------------------------------------------
  	          # compute the yield for this histogram, considering also the overflow
  	          # -------------------------------------------------------------------

  	          yields[histname] = hists[histname].Integral(0,hists[histname].GetNbinsX()+1)

  	      # -------------------------------------------------
              # subtract prompt bkg from data in fakes histograms
  	      # -------------------------------------------------

              # Before doing this, save a copy of the original histogram (expected/observed)
  	      #
	      hist_original = hists[histname].Clone( histname )

  	      if args.doBkgSub :

  	         if ( args.usePrediction == "MC" ):
  	           sys.exit("trying to subtract MC when using --usePrediction=", args.usePrediction ," option. Please use DATA instead, or switch off prompt background subtraction")

  	         name = None

  	         if ( iType == "Fake" ):

  	           name = iLep + "_"+ iVar +"_"+ iType + "_" +  iSel + "_" + list_prediction[1] # --> "observed"

  	           print "\t\t\t\t subtracting events w/ prompt/ch-flip probe lepton to data in Fake CR..."

          	   hist_sub = hists[ iLep + "_" + iVar + "_" + iType + "_" + iSel + "_" + list_prediction[0] ] # --> "expected"

  	           if args.debug:
  	              print "\t\t\t\t Integral before sub: ", hists[name].Integral(0,hists[histname].GetNbinsX()+1), " - hist name: ", hists[name].GetName()

  	           hists[name].Add( hist_sub, -1 )

  	           if args.debug:
  	              print "\t\t\t\t Integral after sub: ", hists[name].Integral(0,hists[histname].GetNbinsX()+1)

  	         elif ( iType == "Real" ) or ( iType == "RealQMisIDBinning" ):

                   iActualType = iType

  	           name = iLep + "_"+ iVar +"_"+ iActualType + "_" +  iSel + "_" + list_prediction[1] # --> "observed"

  	           print "\t\t\t\t subtracting events w/ !prompt/ch-flip probe lepton to data in Real CR..."

          	   hist_sub = hists[ iLep + "_" + iVar + "_" + iActualType + "_" + iSel + "_" + list_prediction[0] ] # --> "expected"

  	           if args.debug:
  	              print "\t\t\t\t Integral before sub: ", hists[name].Integral(0,hists[histname].GetNbinsX()+1), " - hist name: ", hists[name].GetName()

  	           hists[name].Add( hist_sub, -1 )

  	           if args.debug:
  	              print "\t\t\t\t Integral after sub: ", hists[name].Integral(0,hists[histname].GetNbinsX()+1)


  	         # ----------------------------------------------------------------------
  	         # set bin content to zero if subtraction gives negative yields
  	         #
  	         # re-compute the yield for this histogram, considering also the overflow
  	         # ----------------------------------------------------------------------

  	         for iBin in range( 0, hists[name].GetNbinsX()+2):

  	           if hists[name].GetBinContent( iBin ) < 0:
                       hists[name].SetBinContent( iBin, 0 )

                 yields[name] = hists[name].Integral(0,hists[histname].GetNbinsX()+1)

  	  # ------------------
  	  # compute the rate:
  	  # T / !T
  	  #
  	  # and the efficiency:
  	  # T / ( !T + T ) = T / L
  	  # ------------------

  	  histname =  iLep + "_" + iVar +"_"+ iType

          print "\t\t histogram to be used in the ratio: ", histname

          # ---------------------------------------------------
  	  # use MC based estimate, but only if selected by user
          # ---------------------------------------------------

  	  append_str = "observed"
          if ( args.usePrediction == "MC" ):
  	    append_str = "expected"

  	  print "Numerator T: tot. yield = ",  yields[histname + "_T_" + append_str]
          for bin in range(1,hists[histname + "_T_" + append_str].GetNbinsX()+1):
            print("\t Bin nr: {0}, [{1},{2}] - yield = {3}".format(bin,hists[histname + "_T_" + append_str].GetBinLowEdge(bin),hists[histname + "_T_" + append_str].GetBinLowEdge(bin+1),hists[histname + "_T_" + append_str].GetBinContent(bin)))
  	  print "Denominator AntiT: tot. yield = ",  yields[histname + "_AntiT_" + append_str]
          for bin in range(1,hists[histname + "_AntiT_" + append_str].GetNbinsX()+1):
            print("\t Bin nr: {0}, [{1},{2}] - yield = {3}".format(bin,hists[histname + "_AntiT_" + append_str].GetBinLowEdge(bin),hists[histname + "_AntiT_" + append_str].GetBinLowEdge(bin+1),hists[histname + "_AntiT_" + append_str].GetBinContent(bin)))

          hists[histname + "_Rate_" + append_str]  = hists[histname + "_T_" + append_str].Clone(histname + "_Rate_" + append_str)
          hists[histname + "_Rate_" + append_str].Divide(hists[histname + "_AntiT_" + append_str])

          yields[histname + "_Rate_" + append_str] = yields[histname + "_T_" + append_str] / yields[histname + "_AntiT_" + append_str]

          # For efficiency, make sure the errors are computed correctly
          # (NB: in this case, numerator and denominator are not independent sets of events!)
          #
          hist_pass = hists[histname + "_T_" + append_str]
          hist_tot  = hists[histname + "_T_" + append_str] + hists[histname + "_AntiT_" + append_str]
          hist_eff  = hist_pass.Clone(histname + "_Efficiency_" + append_str)
          hist_eff.Divide(hist_pass, hist_tot,1.0,1.0,"B")

	  # Cache the electron REAL efficiency (only in DATA, and for pT dependence)
	  #
  	  if ( args.usePrediction == "DATA" ) and ( iType == "RealQMisIDBinning" ) and ( iLep == "El" ) and ( iVar == "ProbePt" ):
	    hist_eff_electron_real = hist_eff

  	  # Get the actual electron FAKE efficiency by scaling the REAL efficiency appropriately by the QMisID rate ratio
  	  #
  	  hist_eff_QMisID = hist_eff_fake_scaled = None

  	  if ( args.usePrediction == "DATA" ) and ( iType == "Fake" ) and ( iLep == "El" ) and ( iVar == "ProbePt" ):

	     if not hist_eff_electron_real:
  	       sys.exit("ERROR: histogram hist_eff_electron_real does not exist")

  	     hist_eff_QMisID, hist_eff_fake_scaled = scaleEff( hist_eff_electron_real, hist_eff, hist_original, hists[histname_QMisID] )

  	     eff_QMisID_name         = (histname + "_Efficiency_" + append_str).replace("Fake","QMisID")
  	     eff_fake_scaled_name    = (histname + "_Efficiency_" + append_str).replace("Fake","ScaledFake")

  	     hist_eff_QMisID.SetName(eff_QMisID_name)
  	     hist_eff_fake_scaled.SetName(eff_fake_scaled_name)

             hists[eff_QMisID_name]	 = hist_eff_QMisID
             hists[eff_fake_scaled_name] = hist_eff_fake_scaled

          hists[histname + "_Efficiency_" + append_str] = hist_eff
          #
          # even better, use TGraphAsymmErrors
          # (handles the errors in case eff is 100%)
          #
          g_efficiency = TGraphAsymmErrors(hist_eff)
          g_efficiency.Divide(hist_pass,hist_tot,"cl=0.683 b(1,1) mode")
          graphs[histname + "_Efficiency_" + append_str + "_graph"] = g_efficiency

  	  print "Denominator AntiT+T: tot. yield = ", hist_tot.Integral(0,hist_tot.GetNbinsX()+1)
          for bin in range(1,hist_tot.GetNbinsX()+1):
            print("\t Bin nr: {0}, [{1},{2}] - yield = {3}".format(bin,hist_tot.GetBinLowEdge(bin),hist_tot.GetBinLowEdge(bin+1),hist_tot.GetBinContent(bin)))

          print "\t\t --> RATE hist name: ", histname + "_Rate_" + append_str
          print "\t\t --> EFFICIENCY hist name: ", histname + "_Efficiency_" + append_str

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

     if not hists[h]:
       continue

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
