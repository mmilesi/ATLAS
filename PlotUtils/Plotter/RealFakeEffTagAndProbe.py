#!/usr/bin/env python

""" RealFakeEffTagAndProbe.py: measure r/f efficiencies/factors for Matrix Method and Fake Factor Method with T&P method"""

__author__     = "Marco Milesi"
__email__      = "marco.milesi@cern.ch"
__maintainer__ = "Marco Milesi"

import os, sys, array, math

sys.path.append(os.path.abspath(os.path.curdir))

# -------------------------------
# Parser for command line options
# -------------------------------
import argparse

parser = argparse.ArgumentParser(description="Module for deriving real/fake lepton efficiencies/factors for MM/FF")

#***********************************
# positional arguments (compulsory!)
#***********************************

parser.add_argument("inputpath", metavar="inputpath",type=str,
                  help="path to the directory containing subdirs w/ input files")

#*******************
# optional arguments
#*******************

parser.add_argument('--variables', dest='variables', action='store', type=str, nargs='*',
                  help='List of variables to be considered. Use a space-separated list. If unspecified, will consider pT only.')
parser.add_argument("--channel", metavar="channel", default="", type=str,
                  help="Flavour composition of the two leptons in CR to be considered (\"ElEl\", \"MuMu\", \"OF\"). If unspecified, will consider all combinations.")
parser.add_argument("--closure", dest="closure", action="store_true",default=False,
                  help="Estimate efficiencies using MonteCarlo (for closure test)")
parser.add_argument("--debug", dest="debug", action="store_true",default=False,
                  help="Run in debug mode")
parser.add_argument("--verbose", dest="verbose", action="store_true",default=False,
                  help="Run in verbose mode")
parser.add_argument("--nosub", dest="nosub", action="store_true",default=False,
                  help="Do not subtract backgrounds to data (NB: subtraction is disabled by default when running w/ option --closure)")
parser.add_argument('--systematics', dest='systematics', action='store', default="", type=str, nargs='*',
                  help='Option to pass a list of systematic variations to be considered. Use a space-separated list.')
parser.add_argument("--log", dest="log", action="store_true",default=False,
                  help="Read plots with logarithmic Y scale.")
parser.add_argument('--rebin', dest='rebin', action='store', type=str, nargs='+',
                  help='Option to pass a bin range for rebinning. Use space-separated sets of options and numbers (comma-separated) to define new bins, specifying whether it should apply to real/fake muons/electrons for a given variable in the following way (eg.):\n --rebin Real,El,Pt,10,20,200 Fake,Mu,Eta,0.0,1.3,2.5')
parser.add_argument('--averagehist', dest='averagehist', action='store', type=str, nargs='+',
                  help='Option to get average histograms (i.e., 1 single bin over the full range). Use space-separated sets of options (comma-separated) to tell which histograms should be averaged out, e.g.:\n --averagehist Real,El,Pt Fake,Mu,Eta.\nIf option ALL is specified, all histograms get averaged out.')
parser.add_argument("--factors", dest="factors", action="store_true",default=False,
                  help="Calculate factors (pass/!pass) in addition to efficiencies.")
parser.add_argument("--outfilename", metavar="outfilename", default=None, type=str,
                  help="Name of the output file(s). If unspecified, default is \"LeptonEfficiencies\"")
parser.add_argument("--outpath", metavar="outpath", default=None, type=str,
                  help="Name of directory where to store outputs. If unspecified, default is the input directory.")
parser.add_argument("--plots", dest="plots", action="store_true", default=False,
                  help="Produce efficiency plots.")

args = parser.parse_args()

from ROOT import ROOT, gROOT, Double, TH1, TH1D, TFile, TCanvas, TLegend, TLatex, TGraphAsymmErrors, TEfficiency, kFullCircle, kCircle, kOpenTriangleUp, kDot, kBlue, kOrange

gROOT.Reset()
gROOT.LoadMacro("$HOME/RootUtils/AtlasStyle.C")
from ROOT import SetAtlasStyle
SetAtlasStyle()

TH1.SetDefaultSumw2()

class RealFakeEffTagAndProbe:

    def __init__( self, closure=False, factors=False, variables=[], systematics=[], efficiency=None, nosub=False ):

	self.closure = closure
	self.factors = factors

        self.tp_lep = "Probe"

    	#self.__channels     = {"" : ["El","Mu"], "ElEl": ["El"], "MuMu": ["Mu"], "OF" : ["El","Mu"]}
    	self.__channels     = {"" : ["El"], "ElEl": ["El"], "MuMu": ["Mu"], "OF" : ["El","Mu"]}
        self.__leptons      = []
    	self.__efficiencies = ["Fake"] # ["Real","Fake"]
    	self.__variables    = ["Pt"]
    	self.__selections   = ["L","T","AntiT"]
    	self.__processes      = []
    	self.__processes_sub  = []
    	self.__systematics    = [""]
    	self.__systematicsdirections = ["","up","dn"]


	self.leptons_greek = {"El":"e","Mu":"#mu"}
	self.leptons_full  = {"El":"Electrons","Mu":"Muons"}

        self.__outputpath = os.path.abspath(os.path.curdir)
        self.__outputfile = None
        self.__outputfile_yields = None

	if not self.closure:
	    self.__processes.append("observed")
	    #self.__processes.extend(["expectedbkg","dibosonbkg","ttbarzbkg","raretopbkg","qmisidbkg","wjetsbkg","ttbarwbkg","ttbarbkg","zjetsbkg","singletopbkg","allsimbkg"])
            if not nosub:
                self.__processes_sub.extend(["qmisidbkg","allsimbkg"])
	else:
	    self.__processes.append("expectedbkg")

        if efficiency:
            self.__efficinecies.append(efficiency)

        if variables:
            for var in variables:
                if var not in self.__variables:
                    self.__variables.append(var)

        if systematics:
                self.__systematics.extend(systematics)

        # The following dictionary associates a known process to a list of affecting systematics

        self.__process_syst_dict = {"qmisidbkg":["QMisID"], "allsimbkg":["AllSimStat"]}

    	# -----------------------------------------
    	# these dictionaries will store the inputs
    	# -----------------------------------------

        self.histkeys = []

        self.loose_hists      = {}
        self.tight_hists      = {}
        self.antitight_hists  = {}
    	self.loose_yields     = {}
    	self.tight_yields     = {}
    	self.antitight_yields = {}

    	# -----------------------------------------
    	# these dictionaries will store the outputs
    	# -----------------------------------------

        self.histfactors       = {}
    	self.histefficiencies  = {}
    	self.graphefficiencies = {}
    	self.tefficiencies     = {}

	# ---------------
    	# General options
    	# ---------------

	self.debug = False
        self.verbose = False
	self.log   = False

        self.lumi = 13.2
	self.extensionlist = ["eps","png"]

    def addProcess( self, processlist=None ):
	self.__processes.extend(processlist)

    def __fillNDHistDicts__( self, key, sel, hist ):

        if sel == "T":
            self.tight_hists[key] = hist
        elif sel == "L":
            self.loose_hists[key] = hist
        elif sel == "AntiT":
            self.antitight_hists[key] = hist


    def __getStatsVariedHist__( self, hist, var_direction ):

        var_hist = hist.Clone(hist.GetName())

        for ibin in range(0, hist.GetNbinsX()+2):
            value    = hist.GetBinContent(ibin)
            stat_err = hist.GetBinError(ibin)
            if var_direction == "up":
                var_hist.SetBinContent(ibin, value + stat_err )
                var_hist.SetBinError(ibin, stat_err)
            if var_direction == "dn":
                var_hist.SetBinContent(ibin, value - stat_err )
                var_hist.SetBinError(ibin, stat_err)

        return var_hist


    def readInputs( self, inputpath=None, channel=None, log=False ):

        self.__outputpath = inputpath

        if inputpath and inputpath.endswith("/"):
            inputpath = inputpath[:-1]

        log_suffix = ("","_LOGY")[bool(self.log)]

        self.__leptons = self.__channels[channel]

        print("Leptons to be considered:")
        print("\n".join("{0}".format(lep) for lep in self.__leptons))
        print("********************************************")
        print("Efficiencies to be considered:")
        print("\n".join("{0}".format(eff) for eff in self.__efficiencies))
        print("********************************************")
        print("Variables to be considered:")
        print("\n".join("{0}".format(var) for var in self.__variables))
        print("********************************************")

        filename = None

        for lep in self.__leptons:

            for eff in self.__efficiencies:

                actual_eff = eff

                for var in self.__variables:

                    for sel in self.__selections:

                        filename = ( inputpath + "/" + channel + actual_eff + "CR" + lep + sel + log_suffix + "/" + channel + actual_eff + "CR" + lep + sel + "_" + lep + self.tp_lep + var + ".root" )

                        thisfile = TFile(filename)
                        if not thisfile:
                            sys.exit("ERROR: file:\n{0}\ndoes not exist!".format(filename))

                        for proc in self.__processes:

			    thishist = thisfile.Get(proc)

                            key = "_".join( (actual_eff,lep,var,proc) )

			    if not key in self.histkeys:
			        self.histkeys.append(key)

                            thishist.SetName(key)
                            thishist.SetDirectory(0)
			    self.__fillNDHistDicts__(key, sel, thishist)

                        for subproc in self.__processes_sub:

			    for sys in self.__systematics:

			    	for sysdir in self.__systematicsdirections:

			    	    append = ("_" +sys + "sys" + "_" + sysdir,"")[bool(not sys and not sysdir)]

				    sysprocname = subproc + append

			    	    if thisfile.GetListOfKeys().Contains(sysprocname) or ( sys == "AllSimStat" and sysdir ):

                                        thissyshist = None

                                        if sys == "AllSimStat":
                                            thissyshist = self.__getStatsVariedHist__( thisfile.Get("allsimbkg"), sysdir )
                                        else:
                                            thissyshist = thisfile.Get(sysprocname)

					keyappend = ("_" + sys + "_" + sysdir,"")[bool(not sys and not sysdir)]
			    		syskey = "_".join( (actual_eff,lep,var,subproc) ) + keyappend

			    		if not syskey in self.histkeys:
			    		    self.histkeys.append(syskey)

			    		thissyshist.SetName(syskey)
			    		thissyshist.SetDirectory(0)
			    		self.__fillNDHistDicts__(syskey, sel, thissyshist)


    def __subtract__ ( self, sel, proc_key, proc_sub_key, bin_idx=None, proc_sub_base_key=None ):

	hist = None
	sub_hist = None
	sub_basehist = None

        # Case 1):
	# Subtract the entire histogram (identified by "proc_sub_key") from the "proc_key" histogram

        if ( bin_idx == None and not proc_sub_base_key ) or ( bin_idx != None and  proc_sub_key == proc_sub_base_key ):

           if sel == "T":
	       hist    = self.tight_hists.get(proc_key)
	       sub_hist = self.tight_hists.get(proc_sub_key)
	   elif sel == "L":
	       hist    = self.loose_hists.get(proc_key)
	       sub_hist = self.loose_hists.get(proc_sub_key)
	   elif sel == "AntiT":
	       hist    = self.antitight_hists.get(proc_key)
	       sub_hist = self.antitight_hists.get(proc_sub_key)
	   if sub_hist:
	       if self.verbose:
	           print("\t{0:.3f} ({1}) - {2:.3f} ({3})".format(hist.Integral(0,hist.GetNbinsX()+1),hist.GetName(),sub_hist.Integral(0,sub_hist.GetNbinsX()+1),sub_hist.GetName()))

               hist.Add( sub_hist, -1 )

	       if self.verbose:
	           print("\t ==> = {0:.3f}".format(hist.Integral(0,hist.GetNbinsX()+1)))

        # Case 2):
	# Make a specific subtraction for the bin in question

	elif ( bin_idx != None and proc_sub_key != proc_sub_base_key ):

           if sel == "T":
	       hist        = self.tight_hists.get(proc_key)
	       sub_hist     = self.tight_hists.get(proc_sub_key)
	       sub_basehist = self.tight_hists.get(proc_sub_base_key)
	   elif sel == "L":
	       hist        = self.loose_hists.get(proc_key)
	       sub_hist     = self.loose_hists.get(proc_sub_key)
	       sub_basehist = self.loose_hists.get(proc_sub_base_key)
	   elif sel == "AntiT":
	       hist        = self.antitight_hists.get(proc_key)
	       sub_hist     = self.antitight_hists.get(proc_sub_key)
	       sub_basehist = self.antitight_hists.get(proc_sub_base_key)

	   if sub_hist and sub_basehist:
	       if self.verbose:
	           print("\t{0:.3f} ({1}) - {2:.3f} ({3})".format(hist.Integral(0,hist.GetNbinsX()+1),hist.GetName(),sub_basehist.Integral(0,sub_basehist.GetNbinsX()+1),sub_basehist.GetName()))
	           print("\tbin {0} - {1:.3f} ({2}) - {3:.3f} ({4})".format(bin_idx,hist.GetBinContent(bin_idx),hist.GetName(),sub_hist.GetBinContent(bin_idx),sub_hist.GetName()))

               bin_idx_sub     = hist.GetBinContent(bin_idx) - sub_hist.GetBinContent(bin_idx)
	       bin_idx_sub_err = math.sqrt( pow( hist.GetBinError(bin_idx),2.0) - pow( sub_hist.GetBinError(bin_idx),2.0 ) )

	       # Firstly, subtract the base histogram...

               hist.Add( sub_basehist, -1 )

	       # ...then change the bin in question!

	       hist.SetBinContent( bin_idx, bin_idx_sub )
	       hist.SetBinError( bin_idx, bin_idx_sub_err )

	       if self.verbose:
	           print("\t ==> = {0:.3f}".format(hist.Integral(0,hist.GetNbinsX()+1)))

        return hist


    def subtractHistograms ( self ):

        for key in self.histkeys:

	    tokens = key.split("_")

	    # If the histogram is not data, do not subtract anything

	    if not ( len(tokens) == 4 and "observed" in tokens[3] ) : continue

	    if self.verbose:
	    	print("Subtracting to {0}...".format(key))

	    base = "_".join( (tokens[0],tokens[1],tokens[2]) )

	    for sys in self.__systematics:

	        keyappend_sys = ("_" + sys,"")[bool(not sys)]

	    	for sysdir in self.__systematicsdirections:

		    # Make sure you get either the nominal, or the systematics

		    if not sys and sysdir: continue
		    if sys and not sysdir: continue

	    	    keyappend_sys_sysdir = keyappend_sys

		    if ( sys and sysdir ):
	    	        keyappend_sys_sysdir += "_" + sysdir

	    	    subkeysys = key + "_sub" + keyappend_sys_sysdir

		    if ( not sys and not sysdir):

                        if not subkeysys in self.histkeys:
                            self.histkeys.append(subkeysys)

	    	        self.tight_hists[subkeysys]	= self.tight_hists[key].Clone(subkeysys)
	    	        self.loose_hists[subkeysys]	= self.loose_hists[key].Clone(subkeysys)
	    	        self.antitight_hists[subkeysys] = self.antitight_hists[key].Clone(subkeysys)

		    else:
		    	for ibin in range(1,self.tight_hists[key].GetNbinsX()+2):

		     	    subkeysys_ibin = subkeysys + "_" + str(ibin)

                            if not subkeysys_ibin in self.histkeys:
                                self.histkeys.append(subkeysys_ibin)

	    	    	    self.tight_hists[subkeysys_ibin]	 = self.tight_hists[key].Clone(subkeysys_ibin)
	    	    	    self.loose_hists[subkeysys_ibin]	 = self.loose_hists[key].Clone(subkeysys_ibin)
	    	    	    self.antitight_hists[subkeysys_ibin] = self.antitight_hists[key].Clone(subkeysys_ibin)

	            for subproc in self.__processes_sub:

			# Key for the nominal processes to be subtracted

			subprockey_base = "_".join( (base,subproc) )

			# Key for the systematically varied processes to be subtracted (NB: systematic match is enforced)

			subprockey_sys_sysdir = subprockey_base
			if ( self.__process_syst_dict.get(subproc) ) and sys in self.__process_syst_dict.get(subproc):
			    subprockey_sys_sysdir += keyappend_sys_sysdir

			# Nominal subtraction

			if ( not sys and not sysdir):

                            if self.verbose:
                                print "A -B ==>"
                                print "subkeysys (A): ", subkeysys
                                print "subprockey (B): ", subprockey_sys_sysdir

		    	    self.tight_hists[subkeysys]     = self.__subtract__( "T",     subkeysys, subprockey_sys_sysdir )
            	    	    self.loose_hists[subkeysys]     = self.__subtract__( "L",     subkeysys, subprockey_sys_sysdir )
            	    	    self.antitight_hists[subkeysys] = self.__subtract__( "AntiT", subkeysys, subprockey_sys_sysdir )

			else:

			    for ibin in range(1,self.tight_hists[key].GetNbinsX()+2):

				subkeysys_ibin = subkeysys + "_" + str(ibin)

                                if self.verbose:
                                    print "\n\tA -B ==>"
                                    print "\tsubkeysys (A): ", subkeysys_ibin
                                    print "\tsubprockey (B): ", subprockey_sys_sysdir

		    	   	self.tight_hists[subkeysys_ibin]     = self.__subtract__( "T",     subkeysys_ibin, subprockey_sys_sysdir, ibin, subprockey_base )
            	    	   	self.loose_hists[subkeysys_ibin]     = self.__subtract__( "L",     subkeysys_ibin, subprockey_sys_sysdir, ibin, subprockey_base )
            	    	   	self.antitight_hists[subkeysys_ibin] = self.__subtract__( "AntiT", subkeysys_ibin, subprockey_sys_sysdir, ibin, subprockey_base )


	    if self.verbose:
            	print("**********************************************************")


    def rebinHistograms ( self, rebinlist=None, averagehistlist=None ):

        if not rebinlist and not averagehistlist:
            print("Will not do any rebinning...")
            return

        if rebinlist or averagehistlist:
            print("\nWARNING!\n")
	    print("(From TH1 docs) If ngroup is not an exact divider of the number of bins, the top limit of the rebinned histogram is reduced to the upper edge of the last bin that can make a complete group.")
            print("The remaining bins are added to the overflow bin. Statistics will be recomputed from the new bin contents.\n")
            print("If rebinning to one single bin, this might lead to an \"empty\" histogram (as everything will end up in the overflow bin)\n")

        if rebinlist:
            for key in self.histkeys:

	        # By construction, the following will be a list w/ N items, where:
	        #
	        # tokens[0] = efficiency ("Real","Fake"...)
	        # tokens[1] = lepton ("El","Mu"...)
	        # tokens[2] = variable ("Pt","Eta"...)
	        # tokens[3] = process ("observed","expectedbkg"...)

	        tokens = key.split("_")

	        print("Current tokens:")
	        print tokens

                for rebinitem in rebinlist:

	            rebinitem = rebinitem.split(",")

	            if ( tokens[0] in rebinitem ) and ( tokens[1] in rebinitem ) and ( tokens[2] in rebinitem ):

                        nbins = len(rebinitem[3:])-1

                        print("\t===> Rebinning matching histograms with the following values:")
                        print "\tbin edges: ", rebinitem[3:], ", number of bins = ", nbins

                        bins = [ float(binedge) for binedge in rebinitem[3:] ]
	                arr_bins = array.array("d", bins)

                        self.tight_hists[key]     = self.tight_hists[key].Rebin( nbins, key, arr_bins )
                        self.loose_hists[key]     = self.loose_hists[key].Rebin( nbins, key, arr_bins )
                        self.antitight_hists[key] = self.antitight_hists[key].Rebin( nbins, key, arr_bins )

        if averagehistlist:
            for key in self.histkeys:

		tokens = key.split("_")

		print("Current tokens:")
	        print tokens

		for averageitem in averagehistlist:

		    averageitem = averageitem.split(",")

		    if ( averageitem[0] == "ALL" ) or ( ( tokens[0] in averageitem ) and ( tokens[1] in averageitem ) and ( tokens[2] in averageitem ) ):

		        print("\t===> Taking average on whole bin range of matching histograms")
                        nbins_tight      = self.tight_hists[key].GetNbinsX()
                        nbins_loose      = self.loose_hists[key].GetNbinsX()
                        nbins_antitight  = self.antitight_hists[key].GetNbinsX()
                        self.tight_hists[key]     = self.tight_hists[key].Rebin( nbins_tight, key )
                        self.loose_hists[key]     = self.loose_hists[key].Rebin( nbins_loose, key )
                        self.antitight_hists[key] = self.antitight_hists[key].Rebin( nbins_antitight, key )


    def __yields_and_integral__( self, histogram ):

        yields = [ histogram.GetBinContent(ibin) for ibin in range( 1, histogram.GetNbinsX()+2 ) ]
        yields.append( histogram.Integral(0,histogram.GetNbinsX()+1) )
        return yields


    def storeYields ( self ):

        # Make sure "tight" (numerator) bins have *never* more events
	# than "loose" (denominator) bins. That is, set tight=loose in such cases.
	# Small fluctuations can happen b/c of subtraction and numerical precision (?)

        for key, tight_hist in self.tight_hists.iteritems():
             for ibin in range( 1, tight_hist.GetNbinsX()+2 ):
	         tight_yield = tight_hist.GetBinContent(ibin)
		 loose_yield = self.loose_hists[key].GetBinContent(ibin)
		 if tight_yield > loose_yield:
		     self.loose_hists[key].SetBinContent(ibin,tight_yield)

        # Reset any bin which became negative after subtraction to zero,
        # then get the yields for every bin and the integarl, and store it in a list

        for key, hist in self.tight_hists.iteritems():
            for ibin in range( 1, hist.GetNbinsX()+2 ):
                if hist.GetBinContent(ibin) < 0:
                    hist.SetBinContent(ibin,0.0)
            self.tight_yields[key] = self.__yields_and_integral__(hist)

        for key, hist in self.loose_hists.iteritems():
            for ibin in range( 1, hist.GetNbinsX()+2 ):
                if hist.GetBinContent(ibin) < 0:
                    hist.SetBinContent(ibin,0.0)
            self.loose_yields[key] = self.__yields_and_integral__(hist)

        for key, hist in self.antitight_hists.iteritems():
            for ibin in range( 1, hist.GetNbinsX()+2 ):
                if hist.GetBinContent(ibin) < 0:
                    hist.SetBinContent(ibin,0.0)
            self.antitight_yields[key] = self.__yields_and_integral__(hist)


    def sanityCheck( self ):
    	print("\n\n")
	for key, num_list in self.tight_yields.iteritems():
    	     numerator   = num_list[-1]
    	     denominator = self.loose_yields[key][-1]
    	     if denominator < numerator:
    		 print("WARNING! Histogram:{0} ==> L = {1:.3f}, T = {2:.3f}, !T = {3:.3f} ==> L < T ??".format(key, denominator, numerator, self.antitight_yields[key][-1] ))
                 print("Numerator T: integral = {0}".format(self.tight_yields[key][-1]))
                 print("\t" + " ; ".join( "({0},{1:.3f})".format(bin,value) for bin,value in enumerate(self.tight_yields[key][:-1])))
                 print("Denominator L: integral = {0}".format(self.loose_yields[key][-1]))
                 print("\t" + " ; ".join( "({0},{1:.3f})".format(bin,value) for bin,value in enumerate(self.loose_yields[key][:-1])))
                 print("Anti T: integral = {0}".format(self.antitight_yields[key][-1]))
                 print("\t" + " ; ".join( "({0},{1:.3f})".format(bin,value) for bin,value in enumerate(self.antitight_yields[key][:-1])))


    # NB: When computing an efficiency, need to make sure the errors are computed correctly!
    # In this case, numerator and denominator are not independent sets of events! The efficiency is described by a binomial PDF.

    def computeEfficiencies( self, variation ):

        print("\nCalculating EFFICIENCIES...\n")

        nominal_key = None

        # NB: here is crucial to loop over the alphabetically-sorted keys!

        for key in sorted(self.histkeys):

            tokens = key.split("_")

            if ( len(tokens) < 5 or tokens[3] not in self.__processes ): continue

            if ( variation == "nominal"):
                if ( len(tokens) > 5 ): continue
                nominal_key = key

            # If checking systematics, need to store the nominal key first

            if ( variation != "nominal" and len(tokens) == 5 ):
                nominal_key = key
                continue

            if self.debug:
                print("\nnominal key: {0}\nvar key: {1}\n".format(nominal_key,key))

            # Define numerator (pass) and denominator (total)

            h_pass = h_tot = None

            append = ""
            if variation == "nominal":
                h_pass = self.tight_hists[nominal_key]
                h_tot  = self.tight_hists[nominal_key] + self.antitight_hists[nominal_key]
            elif variation == "numerator":
                append = "numerator"
                h_pass = self.tight_hists[key]
                h_tot  = self.tight_hists[key] + self.antitight_hists[nominal_key]
            elif variation == "denominator":
                append = "denominator"
                h_pass = self.tight_hists[nominal_key]
                h_tot  = self.tight_hists[nominal_key] + self.antitight_hists[key]

	    ratiolist = []
	    for idx, elem in enumerate(self.tight_yields[key]):
	        n = elem
		d = self.loose_yields[key][idx]
		r = -1 # just an unphysical value
		if d:
		    r = n/d
		ratiolist.append(r)

            if self.verbose:
                print("*****************************************************")
                print("Histogram: {0}\n".format(key))
                print("Numerator T: integral = {0}".format(self.tight_yields[key][-1]))
                print("\t" + " ; ".join( "({0},{1:.3f})".format(bin,value) for bin,value in enumerate(self.tight_yields[key][:-1])))
                print("Denominator L (AntiT+T): integral = {0}".format(self.loose_yields[key][-1]))
                print("\t" + " ; ".join( "({0},{1:.3f})".format(bin,value) for bin,value in enumerate(self.loose_yields[key][:-1])))
                print("AntiT: integral = {0}".format(self.antitight_yields[key][-1]))
                print("\t" + " ; ".join( "({0},{1:.3f})".format(bin,value) for bin,value in enumerate(self.antitight_yields[key][:-1])))
                print("Ratio T/L: integral = {0}".format(ratiolist[-1]))
                print("\t" + " ; ".join( "({0},{1:.3f})".format(bin,value) for bin,value in enumerate(ratiolist[:-1])))
                print("")

            # 1.
            #
	    # The TH1::Divide method with the option "B" calculates binomial errors using the "normal" approximation
            # (NB: the approximation fails when eff = 0 or 1. In such cases, TEfficiency or TGraphAsymmErrors should be used, since they know how to handle such cases)
	    #
            key_heff = "_".join( (tokens[0],tokens[1],tokens[2],"Efficiency",tokens[3],tokens[4],append) )
            if key_heff.endswith("_"):
                key_heff = key_heff[:-1]
            if len(tokens) > 6:
                key_heff = key_heff + "_" + "_".join( ("{0}".format(other_tokens) for other_tokens in tokens[5:]) )

            print "key for efficiency: ", key_heff
            #continue

            h_efficiency  = h_pass.Clone(key_heff)
            h_efficiency.Divide(h_pass,h_tot,1.0,1.0,"B")

            # 2.
            # The TEfficiency class handles the special cases not covered by TH1::Divide
            #
	    t_efficiency = None
            key_teff = key_heff.replace( "Efficiency", "TEfficiency" )


	    if TEfficiency.CheckConsistency(h_pass, h_tot,"w"):
	        t_efficiency = TEfficiency(h_pass, h_tot)
                t_efficiency.SetName(key_teff)

                t_efficiency.SetConfidenceLevel(0.683)

                # 2.a
                #
                # Use TEfficiency, with the frequentist Clopper-Pearson confidence interval at XX% CL (set before)
                # (this handles the eff = 0, 1 case)
                #
                # DOES NOT SEEM TO WORK FOR WEIGHTED HISTOGRAMS --> IT REDUCES TO THE NORMAL APPROX
                #
	        #t_efficiency.SetStatisticOption(TEfficiency.kFCP)

                # 2.b
                #
                # Use TEfficiency, with the Bayesian uniform prior at XX% CL (set before)
                # (This is the same as the TGraphAsymmErrors below)
                #
                # In order to get the same value for efficiency as in a frequentist approach (as from TH1::Divide("B")),
                # the MODE should be used as an estimator. This works as long as a uniform prior is chosen.
                #
                # Please refer to the TEfficiency class docs for details:
                # https://root.cern.ch/doc/master/classTEfficiency.html
                #
	        t_efficiency.SetStatisticOption(TEfficiency.kBUniform)
                t_efficiency.SetPosteriorMode()

            # 3.
            #
            # Calculate efficiency using TGraphAsymmErrors
            # (uses a Bayesian uniform prior beta(1,1) for the efficiency with 68% CL. It handles the errors in case eff is 0 or 1)
            #
            key_geff = key_heff
            g_efficiency = TGraphAsymmErrors(h_efficiency)
            g_efficiency.Divide(h_pass,h_tot,"cl=0.683 b(1,1) mode")

            # Save the efficiencies in the proper dictionaries
            #
            self.histefficiencies[key_heff]  = h_efficiency
            self.graphefficiencies[key_geff] = g_efficiency
       	    self.tefficiencies[key_teff]     = t_efficiency

    def computeFactors( self, variation ):

	print("\nCalculating FACTORS...\n")

        nominal_key = None

        for key in sorted(self.histkeys):

	    tokens = key.split("_")

            if ( len(tokens) < 5 or tokens[3] not in self.__processes ): continue

            if ( variation == "nominal"):
                if ( len(tokens) > 5 ): continue
                nominal_key = key

            # If checking systematics, need to store the nominal key

            if ( variation != "nominal" and len(tokens) == 5 ):
                nominal_key = key
                continue

            if self.debug:
                print("\nnominal key: {0}\nvar key: {1}\n".format(nominal_key,key))

	    # Define numerator (pass) and denominator (not-pass)
	    # These are two set of *independent* events, so just doing the hist ratio will be ok.

            h_pass = h_notpass = None

            append = ""
            if variation == "nominal":
                h_pass     = self.tight_hists[nominal_key]
                h_notpass  = self.antitight_hists[nominal_key]
            elif variation == "numerator":
                append = "numerator"
                h_pass     = self.tight_hists[key]
                h_notpass  = self.antitight_hists[nominal_key]
            elif variation == "denominator":
                append = "denominator"
                h_pass     = self.tight_hists[nominal_key]
                h_notpass  = self.antitight_hists[key]

	    ratiolist = []
	    for idx, elem in enumerate(self.tight_yields[key]):
	        n = elem
		d = self.antitight_yields[key][idx]
		r = -1 # just an unphysical value
		if d:
		    r = n/d
		ratiolist.append(r)

	    if self.verbose:
                print("*****************************************************")
                print("Histogram: {0}\n".format(key))
                print("Numerator T: integral = {0}".format(self.tight_yields[key][-1]))
                print("\t" + " ; ".join( "({0},{1:.3f})".format(bin,value) for bin,value in enumerate(self.tight_yields[key][:-1])))
                print("Denominator AntiT: integral = {0}".format(self.antitight_yields[key][-1]))
                print("\t" + " ; ".join( "({0},{1:.3f})".format(bin,value) for bin,value in enumerate(self.antitight_yields[key][:-1])))
                print("Ratio T/AntiT: integral = {0}".format(ratiolist[-1]))
                print("\t" + " ; ".join( "({0},{1:.3f})".format(bin,value) for bin,value in enumerate(ratiolist[:-1])))
                print("")

            key_hfactor = "_".join( (tokens[0],tokens[1],tokens[2],"Factor",tokens[3],tokens[4],append) )
            if key_hfactor.endswith("_"):
                key_hfactor = key_hfactor[:-1]

            if len(tokens) > 6:
                key_hfactor = key_hfactor + "_" + "_".join( ("{0}".format(other_tokens) for other_tokens in tokens[5:]) )

            h_factor  = h_pass.Clone(key_hfactor)
            h_factor.Divide(h_pass,h_notpass)

            # Save the factors in the proper dictionary
            #
            self.histfactors[key_hfactor] = h_factor

    def saveOutputs( self, filename="LeptonEfficiencies", outputpath=None ):

	if not outputpath:
	    outputpath = self.__outputpath

	self.__outputpath = outputpath

        print("\nStoring output files:\n{0}\n{1}\nin directory: {2}".format(filename+".root",filename+".txt",os.path.abspath(outputpath)))

        self.__outputfile_yields = open(self.__outputpath+"/"+filename+".txt","w")
        self.__outputfile_yields.write( "Efficiencies/Factors for Fake Factor amd Matrix Method\n")

        self.__outputfile = TFile(self.__outputpath+"/"+filename+".root","RECREATE")
        self.__outputfile.cd()

    	for key, h in sorted(self.histefficiencies.iteritems()):
    	    if not h: continue
    	    if self.debug: print("\nSaving histogram: {0}".format(key))
    	    h.Write()
    	    eff=[]
    	    for ibin in range( 1, h.GetNbinsX()+2 ):
    	    	myset = [ ibin, h.GetBinLowEdge(ibin), h.GetBinLowEdge(ibin+1), h.GetBinContent(ibin), h.GetBinError(ibin)]
    	    	eff.append( myset )
    	    self.__outputfile_yields.write("%s:\n" %(key) )
    	    for myset in eff:
    	    	self.__outputfile_yields.write("{ %s }; \n" %( "Bin nr: " + str(myset[0]) + " [" + str(round(myset[1],3)) + "," + str(round(myset[2],3)) + "], efficiency (from TH1::Divide(\"B\")) = " + str(round(myset[3],3)) + " +- " + str(round(myset[4],3)) ) )

    	for key, t in sorted(self.tefficiencies.iteritems()):
    	    if not t: continue
    	    if self.debug: print("\nSaving TEfficiency: {0}".format(key))
    	    t.Write()
    	    teff=[]
    	    for ibin in range( 1, t.GetTotalHistogram().GetNbinsX()+2 ):
    	    	myset = [ ibin, t.GetTotalHistogram().GetBinLowEdge(ibin), t.GetTotalHistogram().GetBinLowEdge(ibin+1), t.GetEfficiency(ibin), t.GetEfficiencyErrorUp(ibin), t.GetEfficiencyErrorLow(ibin)]
    	    	teff.append( myset )
    	    self.__outputfile_yields.write("%s:\n" %(key) )
    	    for myset in teff:
    	    	self.__outputfile_yields.write("{ %s }; \n" %( "Bin nr: " + str(myset[0]) + " [" + str(round(myset[1],3)) + "," + str(round(myset[2],3)) + "], efficiency (from TEfficiency) = " + str(round(myset[3],3)) + " + " + str(round(myset[4],3)) + " - " + str(round(myset[5],3)) ) )

    	for key, g in sorted(self.graphefficiencies.iteritems()):
    	    if not g: continue
    	    if self.debug: print("\nSaving graph: {0}".format(key))
    	    g.Write()
    	    geff=[]
    	    for ipoint in range( 0, g.GetN()+1 ):
    	    	x = Double(0)
    	    	y = Double(0)
    	    	g.GetPoint(ipoint,x,y)
    	    	myset = [ ipoint+1, y, g.GetErrorYhigh(ipoint), g.GetErrorYlow(ipoint) ]
    	    	geff.append( myset )
    	    self.__outputfile_yields.write("%s:\n" %(key) )
    	    for myset in geff:
    	    	self.__outputfile_yields.write("{ %s }; \n" %( "Bin nr: " + str(myset[0]) + ", efficiency (from TGraphAsymmErrors) = " + str(round(myset[1],3)) + " + " + str(round(myset[2],3)) + " - " + str(round(myset[3],3)) ) )

	if self.factors:
    	    for key, h in sorted(self.histfactors.iteritems()):
    	    	if not h: continue
    	    	if self.debug: print("\nSaving histogram: {0}".format(key))
    	    	h.Write()
    	    	eff=[]
    	    	for ibin in range( 1, h.GetNbinsX()+2 ):
    	    	    myset = [ ibin, h.GetBinLowEdge(ibin), h.GetBinLowEdge(ibin+1), h.GetBinContent(ibin), h.GetBinError(ibin)]
    	    	    eff.append( myset )
    	    	self.__outputfile_yields.write("%s:\n" %(key) )
    	    	for myset in eff:
    	    	    self.__outputfile_yields.write("{ %s }; \n" %( "Bin nr: " + str(myset[0]) + " [" + str(round(myset[1],3)) + "," + str(round(myset[2],3)) + "], factor (from TH1::Divide()) = " + str(round(myset[3],3)) + " +- " + str(round(myset[4],3)) ) )

	self.__outputfile_yields.close()
        self.__outputfile.Write()
        self.__outputfile.Close()


    def plotMaker( self, run_batch=True ):

        if run_batch:
	    gROOT.SetBatch(True)

        proc = ("observed","expectedbkg")[bool(self.closure)]

        for var in self.__variables:

	    for lep in self.__leptons:

	        c = TCanvas("c1","Efficiencies")
                c.SetFrameFillColor(0)
                c.SetFrameFillStyle(0)
                c.SetFrameBorderMode(0)

                legend = TLegend(0.45,0.5,0.925,0.8) # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
                legend.SetHeader(self.leptons_full[lep])
                legend.SetBorderSize(0)     # no border
                legend.SetFillStyle(0)      # Legend transparent background
                legend.SetTextSize(0.035)   # Increase entry font size!
                legend.SetTextFont(42)      # Helvetica

                leg_ATLAS  = TLatex()
                leg_lumi   = TLatex()
                leg_ATLAS.SetTextSize(0.04)
                leg_ATLAS.SetNDC()
                leg_lumi.SetTextSize(0.04)
                leg_lumi.SetNDC()

	        for idx_eff, eff in enumerate(self.__efficiencies):

		    for idx_proc, proc in enumerate(self.__processes):

                    	key  = "_".join((eff,lep,var,"Efficiency",proc))

			print "\tplotting histogram: ", key

		    	hist = self.histefficiencies[key]

		    	hist.GetYaxis().SetRangeUser(0,1)

		    	hist.GetYaxis().SetTitle("#varepsilon")
	   	    	hist.GetXaxis().SetTitleOffset(1.0)
	   	    	hist.GetYaxis().SetTitleOffset(1.0)

	            	hist.SetLineStyle(1)
                    	hist.SetMarkerStyle(kFullCircle)

			if not ( proc == "observed" ):
			    hist.SetLineStyle(idx_proc+3)
                    	    hist.SetMarkerStyle(kCircle)

			if eff == "Real":
		    	   hist.SetLineColor(kBlue)
		    	   hist.SetMarkerColor(kBlue)
		    	else:
		    	   hist.SetLineColor(kOrange)
		    	   hist.SetMarkerColor(kOrange)

                    	legend.AddEntry(hist,eff+" - "+proc, "P")

		    	if not idx_eff and not idx_proc:
			   hist.Draw("E0")
		    	else:
		    	   hist.Draw("E0,SAME")

                legend.Draw()
                leg_ATLAS.DrawLatex(0.6,0.35,"#bf{#it{ATLAS}} Work In Progress");
                leg_lumi.DrawLatex(0.6,0.27,"#sqrt{{s}} = 13 TeV, #int L dt = {0} fb^{{-1}}".format(str(self.lumi)));

                canvas_filename = "_".join(("RealFake",lep,var,"Efficiency",proc))

		for extension in self.extensionlist:
		    c.SaveAs(self.__outputpath+"/"+canvas_filename+"."+extension)


    def __set_fancy_2D_style( self ):

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


    def __factorToEfficiency__(self, f):
        if f < 0:
            f = 0.0
	e = f/(f+1)
        return e

    def __efficiencyToFactor__(self, e):
	f = e/(1-e)
        return f

    def checkRebin(self):

    	if self.debug:
    	    print("\n\nTIGHT histograms dictionary:\n")
    	    print("\tkey\t\thistname\t\tnbins\n")
    	    for key, value in sorted( self.tight_hists.iteritems() ):
    		print("\t{0}\t{1}\t{2}".format(key, value.GetName(), value.GetNbinsX()))

    	    print("\n\nLOOSE histograms dictionary:\n")
    	    print("\tkey\t\thistname\t\tnbins\n")
    	    for key, value in sorted( self.loose_hists.iteritems() ):
    		print("\t{0}\t{1}\t{2}".format(key, value.GetName(), value.GetNbinsX()))

    	    print("\n\nANTI-TIGHT histograms dictionary:\n")
    	    print("\tkey\t\thistname\t\tnbins\n")
    	    for key, value in sorted( self.antitight_hists.iteritems() ):
    		print("\t{0}\t{1}\t{2}".format(key, value.GetName(), value.GetNbinsX()))


    def checkYields(self, debug_msg=None):

	if self.debug:

	    if debug_msg:
	        print("\n{0}".format(debug_msg))

    	    print("\n\nTIGHT yields dictionary:\n")
    	    print("\tkey\t\tyields (per bin)\t\tintegral\n")
    	    for key, value in sorted( self.tight_yields.iteritems() ):
    		print("\t{0}".format(key) + "\t[" + ",".join( "{0:.2f}".format(i) for i in value[:-1] ) + "]" + "\t{0:.2f}".format(value[-1]) )

    	    print("\n\nLOOSE yields dictionary:\n")
    	    print("\tkey\t\tyields (per bin)\t\tintegral\n")
    	    for key, value in sorted( self.loose_yields.iteritems() ):
    		print("\t{0}".format(key) + "\t[" + ",".join( "{0:.2f}".format(i) for i in value[:-1] ) + "]" + "\t{0:.2f}".format(value[-1]) )

    	    print("\n\nANTI-TIGHT yields dictionary:\n")
    	    print("\tkey\t\tyields (per bin)\t\tintegral\n")
    	    for key, value in sorted( self.antitight_yields.iteritems() ):
    		print("\t{0}".format(key) + "\t[" + ",".join( "{0:.2f}".format(i) for i in value[:-1] ) + "]" + "\t{0:.2f}".format(value[-1]) )


# --------------------------------------------------------------------------------------------------------------

if __name__ == "__main__":

    eff = RealFakeEffTagAndProbe( closure=args.closure, factors=args.factors, variables=args.variables, systematics=args.systematics, nosub=args.nosub )

    eff.debug = args.debug
    eff.log   = args.log

    #eff.addProcess(processlist=["expectedbkg"])

    eff.readInputs( inputpath=args.inputpath, channel=args.channel )

    eff.rebinHistograms( rebinlist=args.rebin, averagehistlist=args.averagehist )
    eff.checkRebin()

    eff.storeYields()
    eff.checkYields("events BEFORE subtraction")

    eff.subtractHistograms()

    eff.storeYields()
    eff.checkYields("events AFTER subtraction")

    eff.computeEfficiencies(variation="nominal")
    print("\n")
    eff.computeEfficiencies(variation="numerator")
    print("\n")
    eff.computeEfficiencies(variation="denominator")

    if eff.factors:
        eff.computeFactors("nominal")
        print("\n")
        eff.computeFactors("numerator")
        print("\n")
        eff.computeFactors("denominator")

    eff.saveOutputs( filename=args.outfilename )

    if args.plots:
        eff.plotMaker()

