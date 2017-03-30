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

g_avaialble_systematics = ["QMisID","TTV","VV","OtherPromptSS","FakesOS"]

g_luminosities = { "GRL v73 - Moriond 2016 GRL":3.209,  # March 2016
                   "ICHEP 2015+2016 DS":13.20768,       # August 2016
                   "POST-ICHEP 2015+2016 DS":22.07036,  # October 2016
                   "FULL 2015+2016 DS":36.0746          # December 2016
                 }

g_selections   = ["L","T","AntiT"]

g_efficiencies = ["ALL","Real","Fake"]

g_leptons      = ["ALL","El","Mu"]

parser.add_argument("--lumi", dest="lumi", action="store", type=float, default=g_luminosities["FULL 2015+2016 DS"],
                  help="The luminosity of the dataset. Pick one of these values: ==> " + ",".join( "{0} ({1})".format( lumi, tag ) for tag, lumi in g_luminosities.iteritems() ) + ". Default is {0}".format(g_luminosities["FULL 2015+2016 DS"] ) )
parser.add_argument('--variables', dest='variables', action='store', type=str, nargs='+', default=["Pt"],
                  help='List of variables to be considered. Use a space-separated list. If using a 2D input distribution, pass a string of format \"varxname&&varyname\" (remember to escape the \"&\" symbol in the shell!). If unspecified, will consider \"Pt\" only.')
parser.add_argument('--efficiency', dest='efficiency', action='store', default=[g_efficiencies[0]], type=str, nargs='+', choices=g_efficiencies,
                  help='The efficiency type to be measured. Can pass multiple space-separated arguments to this command-line option (picking among the above list). If this option is not specified, default will be \'{0}\''.format(g_efficiencies[0]))
parser.add_argument("--do3L", dest="do3L", action="store_true",default=False,
                  help="Save efficiencies to be used in 3L SR.")
parser.add_argument("--channel", metavar="channel", default="", type=str,
                  help="Flavour composition of the two leptons in CR to be considered (\"ElEl\", \"MuMu\", \"OF\"). If unspecified, will consider all combinations.")
parser.add_argument('--lepton', dest='lepton', action='store', default=[g_leptons[0]], type=str, nargs='+', choices=g_leptons,
                  help='The lepton flavour chosen. Can pass multiple space-separated arguments to this command-line option (picking among the above list). If this option is not specified, default will be \'{0}\''.format(g_leptons[0]))
parser.add_argument("--closure", dest="closure", action="store_true",default=False,
                  help="Estimate efficiencies using MonteCarlo (for closure test)")
parser.add_argument("--debug", dest="debug", action="store_true",default=False,
                  help="Run in debug mode")
parser.add_argument("--verbose", dest="verbose", action="store_true",default=False,
                  help="Run in verbose mode")
parser.add_argument("--nosub", dest="nosub", action="store_true",default=False,
                  help="Do not subtract backgrounds to data (NB: subtraction is disabled by default when running w/ option --closure)")
parser.add_argument('--systematics', dest='systematics', action='store', default="", type=str, nargs='*',
                  help='Option to pass a list of systematic variations to be considered. Use a space-separated list. Currently available systematics: {0}'.format(g_avaialble_systematics))
parser.add_argument("--log", dest="log", action="store_true",default=False,
                  help="Read plots with logarithmic Y scale.")
parser.add_argument('--rebin', dest='rebin', action='store', type=str, nargs='+',
                  help='Option to pass a bin range for rebinning. Use space-separated sets of options and numbers (comma-separated) to define new bins, specifying whether it should apply to real/fake muons/electrons for a given variable in the following way (eg.):\n --rebin Real,El,Pt,10,20,200 Fake,Mu,Eta,0.0,1.3,2.5')
parser.add_argument('--averagehist', dest='averagehist', action='store', type=str, nargs='+',
                  help='Option to get average histograms (i.e., 1 single bin over the full range). Use space-separated sets of options (comma-separated) to tell which histograms should be averaged out, e.g.:\n --averagehist Real,El,Pt Fake,Mu,Eta.\nIf option ALL is specified, all histograms get averaged out.')
parser.add_argument("--factors", dest="factors", action="store_true",default=False,
                  help="Calculate factors (pass/!pass) in addition to efficiencies.")
parser.add_argument("--outfilename", metavar="outfilename", default="LeptonEfficiencies", type=str,
                  help="Name of the output file(s). If unspecified, default is \"LeptonEfficiencies\"")
parser.add_argument("--outpath", metavar="outpath", default=None, type=str,
                  help="Name of directory where to store outputs. If unspecified, default is the input directory.")
parser.add_argument("--plots", dest="plots", action="store_true", default=False,
                  help="Produce efficiency plots.")
parser.add_argument("--triggerEff", dest="triggerEff", action="store", default=None, const=g_selections[0], type=str, nargs="?",
                  help="Measure trigger efficiency for a given lepton selection. The lepton selection can be specified as an optional command line argument to this option (Choose between [" + ",".join( "{0}".format( s ) for s in g_selections ) + "]). If no option is specified, default selection will be {0}".format(g_selections[0]))
parser.add_argument("--probeAssignEff", dest="probeAssignEff", action="store_true", default=False,
                  help="Measure probe assignment efficiency.")
parser.add_argument("--update", dest="update", action="store_true", default=False,
                  help="Update existing ROOT output file (Gets overwritten by default).")

args = parser.parse_args()

from ROOT import ROOT, gROOT, gStyle, Double, gPad, TPad, TLine, TH1, TH1D, TH2, TH2D, TFile, TCanvas, TLegend, TLatex, TGraphAsymmErrors, TEfficiency, kFullCircle, kCircle, kOpenTriangleUp, kDot, kBlue, kOrange, kPink, kGreen, kRed, kYellow, kTeal, kMagenta, kViolet, kAzure, kCyan, kSpring, kGray, kBlack, kWhite

from Plotter.BackgroundTools import set_fancy_2D_style

gROOT.Reset()
gROOT.LoadMacro("$HOME/RootUtils/AtlasStyle.C")
from ROOT import SetAtlasStyle
SetAtlasStyle()

TH1.SetDefaultSumw2()

class RealFakeEffTagAndProbe:

    def __init__( self, closure=False, factors=False, variables=[], systematics=[], efficiencies=None, leptons=None, nosub=False ):

	self.closure    = closure
	self.nosub      = nosub
	self.factors    = factors
	self.triggerEff     = None
	self.probeAssignEff = None

        self.tp_lep = "Probe"

    	self.selections       = {"D":"L","N":"T","AntiN":"AntiT"}
        self.__efficiencies   = []

    	self.__channels       = {"" : ["El","Mu"], "ElEl": ["El"], "MuMu": ["Mu"], "OF" : ["El","Mu"]}
        self.__leptons        = []
    	self.__variables      = []
    	self.__processes      = []
    	self.__processes_sub  = []
    	self.__systematics    = [""]
    	self.__systematicsdirections = ["","up","dn"]

	self.leptons_greek = {"El":"e","Mu":"#mu"}
	self.leptons_full  = {"El":"Electrons","Mu":"Muons"}

        self.__outputpath = os.path.abspath(os.path.curdir)
        self.__outputfile = None
        self.__outputfile_yields = None

	self.__averagehists = False

	if not self.closure:
	    self.__processes.append("observed")
            if not self.nosub:
                self.__processes_sub.extend(["qmisidbkg","ttbarwbkg","ttbarzbkg","dibosonbkg","raretopbkg","ttbarbkg","fakesbkg"])
                #self.__processes_sub.extend(["qmisidbkg","dibosonbkg","ttbarzbkg","raretopbkg","wjetsbkg","ttbarwbkg","ttbarbkg","zjetsbkg","singletopbkg"])
	else:
	    self.__processes.append("expectedbkg")

        if leptons:
            for lep in leptons:
                print("lep: {0}".format(lep))
                if lep == "ALL":
                    self.__leptons.extend(["El","Mu"])
                else:
                    self.__leptons.append(lep)

        if efficiencies:
            for eff in efficiencies:
                print("eff: {0}".format(eff))
                if eff == "ALL":
                    self.__efficiencies.extend(["Real","Fake"])
                else:
                    self.__efficiencies.append(eff)

        if variables:
            for var in variables:
                if var not in self.__variables:
                    self.__variables.append(var)

        if systematics:
                self.__systematics.extend(systematics)

        # The following dictionary associates a known process to a list of affecting systematics

        self.__process_syst_dict = {"qmisidbkg":["QMisID"],"ttbarwbkg":["TTV"],"ttbarzbkg":["TTV"],"dibosonbkg":["VV"],"raretopbkg":["OtherPromptSS"],"ttbarbkg":["OtherPromptSS"],"fakesbkg":["FakesOS"]}
        self.__syst_color_dict   = {"QMisID_N":kGreen+3,
                                    "QMisID_D":kGreen-7,
                                    "TTV_ND":kYellow,
                                    "VV_ND":kGreen-7,
                                    "OtherPromptSS_ND":kGray,
                                    "FakesOS_ND":kViolet-4,
                                    }

    	# -----------------------------------------
    	# these dictionaries will store the inputs
    	# -----------------------------------------

        self.histkeys = []

        self.denominator_hists    = {}
        self.numerator_hists      = {}
        self.antinumerator_hists  = {}
    	self.denominator_yields   = {}
    	self.numerator_yields     = {}
    	self.antinumerator_yields = {}

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

	self.debug   = False
        self.verbose = False
	self.log     = False

        self.lumi = 36.4
	self.extensionlist = ["eps","png","root"]

    def addProcess( self, processlist=None ):
	self.__processes.extend(processlist)

    def __fillNDHistDicts__( self, key, eventset, hist ):

        if not key in self.histkeys:
            self.histkeys.append(key)

        if self.verbose:
            print("Saving histogram for N,D,AntiN w/ key: {0}".format(key))

        if eventset == "N":
            self.numerator_hists[key] = hist
        elif eventset == "D":
            self.denominator_hists[key] = hist
        elif eventset == "AntiN":
            self.antinumerator_hists[key] = hist


    def __getStatsVariedHist__( self, hist, var_direction ):

        var_hist = hist.Clone(hist.GetName())

        for ibin in range(0, hist.GetSize()):
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

        if self.__leptons == "ALL":
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

                    for sel_key, sel_value in self.selections.iteritems():

                        if not "&&" in var:
                            filename = ( inputpath + "/" + channel + actual_eff + "CR" + lep + sel_value + log_suffix + "/" + channel + actual_eff + "CR" + lep + sel_value + "_" + lep + self.tp_lep + var + ".root" )
                        else:
                            vars2D = var.split('&&')
                            varX   = vars2D[0]
                            varY   = vars2D[1]
                            filename = ( inputpath + "/" + channel + actual_eff + "CR" + lep + sel_value + log_suffix + "/" + channel + actual_eff + "CR" + lep + sel_value + "_" + lep + self.tp_lep + varX + "_VS_" + lep + self.tp_lep + varY + ".root" )

                        thisfile = TFile(filename)
                        if not thisfile:
                            sys.exit("ERROR: file:\n{0}\ndoes not exist!".format(filename))

                        for proc in self.__processes:

			    thishist = thisfile.Get(proc)

                            key = "_".join( (actual_eff,lep,var,proc) )

                            thishist.SetName(key)
                            thishist.SetDirectory(0)
			    self.__fillNDHistDicts__(key, sel_key, thishist)

                        for subproc in self.__processes_sub:

			    for sys in self.__systematics:

                                if sys == "QMisID" and not ( eff == "Fake" and lep == "El" ):
                                    print("Skipping {0} systematics for {1},{2}".format(sys,eff,lep))
                                    continue
                                if sys in ["TTV","VV","OtherPromptSS"] and not ( eff == "Fake" ):
                                    print("Skipping {0} systematics for {1},{2}".format(sys,eff,lep))
                                    continue
                                if sys == "FakesOS" and not ( eff == "Real" ):
                                    print("Skipping {0} systematics for {1},{2}".format(sys,eff,lep))
                                    continue

			    	for sysdir in self.__systematicsdirections:

			    	    append = ("_" +sys + "sys" + "_" + sysdir,"")[bool(not sys and not sysdir)]

				    sysprocname = subproc + append

			    	    if thisfile.GetListOfKeys().Contains(sysprocname):

                                        thissyshist = thisfile.Get(sysprocname)

					keyappend = ("_" + sys + "_" + sysdir,"")[bool(not sys and not sysdir)]
			    		syskey = "_".join( (actual_eff,lep,var,subproc) ) + keyappend

			    		thissyshist.SetName(syskey)
			    		thissyshist.SetDirectory(0)
			    		self.__fillNDHistDicts__(syskey, sel_key, thissyshist)


        self.storeYields()


    def __subtract__ ( self, sel, proc_key, proc_sub_key, bin_idx=None, proc_sub_base_key=None ):

	hist = None
	sub_hist = None
	sub_basehist = None

        # Case 1):
	# Subtract the entire histogram (identified by "proc_sub_key") from the "proc_key" histogram

        if ( bin_idx == None and not proc_sub_base_key ) or ( bin_idx != None and proc_sub_key == proc_sub_base_key ):

           if sel == "N":
	       hist     = self.numerator_hists.get(proc_key)
	       sub_hist = self.numerator_hists.get(proc_sub_key)
	   elif sel == "D":
	       hist     = self.denominator_hists.get(proc_key)
	       sub_hist = self.denominator_hists.get(proc_sub_key)
	   elif sel == "AntiN":
	       hist     = self.antinumerator_hists.get(proc_key)
	       sub_hist = self.antinumerator_hists.get(proc_sub_key)
	   if sub_hist:
               integral_pre_sub = hist.Integral(0,hist.GetNbinsX()+1)
               hist.Add( sub_hist, -1 )
               integral_post_sub = hist.Integral(0,hist.GetNbinsX()+1)
	       if self.verbose:
	           print("\t{0} : {1:.3f} ({2} [A]) - {3:.3f} ({4} [B]) = {5:.3f}".format(sel,integral_pre_sub,hist.GetName(),sub_hist.Integral(0,sub_hist.GetNbinsX()+1),sub_hist.GetName(),integral_post_sub))
	   else:
               print("\tCouldn't find histogram for: proc_sub_key = {0}".format(proc_sub_key))


        # Case 2):
	# Make a specific subtraction for the bin in question

	elif ( bin_idx != None and proc_sub_key != proc_sub_base_key ):

           if sel == "N":
	       hist         = self.numerator_hists.get(proc_key)
	       sub_hist     = self.numerator_hists.get(proc_sub_key)
	       sub_basehist = self.numerator_hists.get(proc_sub_base_key)
	   elif sel == "D":
	       hist         = self.denominator_hists.get(proc_key)
	       sub_hist     = self.denominator_hists.get(proc_sub_key)
	       sub_basehist = self.denominator_hists.get(proc_sub_base_key)
	   elif sel == "AntiN":
	       hist         = self.antinumerator_hists.get(proc_key)
	       sub_hist     = self.antinumerator_hists.get(proc_sub_key)
	       sub_basehist = self.antinumerator_hists.get(proc_sub_base_key)

	   if sub_hist and sub_basehist:

               integral_pre_sub = hist.Integral(0,hist.GetNbinsX()+1)
               thisbin_pre_sub  = hist.GetBinContent(bin_idx)

               bin_sub     = hist.GetBinContent(bin_idx) - sub_hist.GetBinContent(bin_idx)
	       bin_sub_err = math.sqrt( pow( hist.GetBinError(bin_idx),2.0) + pow( sub_hist.GetBinError(bin_idx),2.0 ) )

	       # Firstly, subtract the base histogram...

               hist.Add( sub_basehist, -1 )

	       # ...then change the bin in question!

	       hist.SetBinContent( bin_idx, bin_sub )
	       hist.SetBinError( bin_idx, bin_sub_err )

               integral_post_sub = hist.Integral(0,hist.GetNbinsX()+1)
               thisbin_post_sub  = hist.GetBinContent(bin_idx)

	       if self.verbose:
	           print("\t{0} : {1:.3f} ({2} [A]) - {3:.3f} ({4} [B]) = {5:.3f}".format(sel,integral_pre_sub,hist.GetName(),sub_basehist.Integral(0,sub_basehist.GetNbinsX()+1),sub_basehist.GetName(),integral_post_sub))
	           print("\t\tbin {0} - {1:.3f} ({2} [A]) - {3:.3f} ({4} [B]) = {5:.3f}".format(bin_idx,thisbin_pre_sub,hist.GetName(),sub_hist.GetBinContent(bin_idx),sub_hist.GetName(),thisbin_post_sub))

	   else:
               if not sub_hist:
                   print("\tCouldn't find histogram for: proc_sub_key = {0}".format(proc_sub_key))
               if not sub_basehist:
                   print("\tCouldn't find histogram for: proc_sub_base_key = {0}".format(proc_sub_base_key))

        return hist


    def subtractHistograms ( self ):

        if self.nosub: return

        for key in self.histkeys:

	    tokens = key.split("_")

	    # If the histogram is not data, do not subtract anything

	    if not ( len(tokens) == 4 and "observed" in tokens[3] ) : continue

	    base = "_".join( (tokens[0],tokens[1],tokens[2]) )

	    for sys in self.__systematics:

                if sys == "QMisID" and not all( s in key for s in ["Fake","El"] ):
                    print("Skipping {0} systematics for {1}".format(sys,key))
                    continue
                if sys in ["TTV","VV","OtherPromptSS"] and not ( "Fake" in key ):
                    print("Skipping {0} systematics for {1}".format(sys,key))
                    continue
                if sys == "FakesOS" and not ( "Real" in key ):
                    print("Skipping {0} systematics for {1}".format(sys,key))
                    continue

	        keyappend_sys = ("_" + sys,"")[bool(not sys)]

	    	for sysdir in self.__systematicsdirections:

		    # Make sure you get either the nominal, or the systematics

		    if not sys and sysdir: continue
		    if sys and not sysdir: continue

	    	    keyappend_sys_sysdir = keyappend_sys

		    if ( sys and sysdir ):
	    	        keyappend_sys_sysdir += "_" + sysdir

	    	    subkeysys = key + "_sub" + keyappend_sys_sysdir

                    if self.verbose:
                        print("\nSubtracting to {0}...".format(key+keyappend_sys_sysdir))

                    if not subkeysys in self.histkeys:
                        self.histkeys.append(subkeysys)

	    	    self.numerator_hists[subkeysys]	= self.numerator_hists[key].Clone(subkeysys)
	    	    self.denominator_hists[subkeysys]   = self.denominator_hists[key].Clone(subkeysys)
	    	    self.antinumerator_hists[subkeysys] = self.antinumerator_hists[key].Clone(subkeysys)

                    if sys or sysdir:
                        for ibin in range(1,self.numerator_hists[key].GetSize()):

                            subkeysys_ibin = subkeysys + "_" + str(ibin)

                            if not subkeysys_ibin in self.histkeys:
                                self.histkeys.append(subkeysys_ibin)

                                self.numerator_hists[subkeysys_ibin]     = self.numerator_hists[key].Clone(subkeysys_ibin)
                                self.denominator_hists[subkeysys_ibin]   = self.denominator_hists[key].Clone(subkeysys_ibin)
                                self.antinumerator_hists[subkeysys_ibin] = self.antinumerator_hists[key].Clone(subkeysys_ibin)

	            for subproc in self.__processes_sub:

			# Key for the nominal processes to be subtracted

			subprockey_base = "_".join( (base,subproc) )

			# Key for the systematically varied processes to be subtracted (NB: systematic match is enforced)

			subprockey_sys_sysdir = subprockey_base
			if ( self.__process_syst_dict.get(subproc) ) and sys in self.__process_syst_dict.get(subproc):
			    subprockey_sys_sysdir += keyappend_sys_sysdir

                        if self.verbose:
                            print("\n\tProcess to be subtracted: {0} (Base process: {1})".format(subprockey_sys_sysdir,subprockey_base))
                            print("")

		    	self.numerator_hists[subkeysys]     = self.__subtract__( "N",     subkeysys, subprockey_sys_sysdir )
            	    	self.denominator_hists[subkeysys]   = self.__subtract__( "D",     subkeysys, subprockey_sys_sysdir )
            	    	self.antinumerator_hists[subkeysys] = self.__subtract__( "AntiN", subkeysys, subprockey_sys_sysdir )

                        if sys or sysdir:

                            # Subtract sys bin by bin

                            for ibin in range(1,self.numerator_hists[key].GetSize()):

                                subkeysys_ibin = subkeysys + "_" + str(ibin)

                                if self.verbose: print("")
                                self.numerator_hists[subkeysys_ibin]     = self.__subtract__( "N",     subkeysys_ibin, subprockey_sys_sysdir, ibin, subprockey_base )
                                self.denominator_hists[subkeysys_ibin]   = self.__subtract__( "D",     subkeysys_ibin, subprockey_sys_sysdir, ibin, subprockey_base )
                                self.antinumerator_hists[subkeysys_ibin] = self.__subtract__( "AntiN", subkeysys_ibin, subprockey_sys_sysdir, ibin, subprockey_base )

            self.storeYields()

	    if self.verbose:
            	print("**********************************************************\n")


    def rebinHistograms ( self, rebinlist=None, averagehistlist=None ):

        if not rebinlist and not averagehistlist:
            print("Will not do any rebinning...")
            return

        if rebinlist or averagehistlist:
            print("\nWARNING!\n")
            print("Will be calling:\n\nTH1::Rebin(Int_t ngroup = 2, const char * newname = \"\", const Double_t * xbins = 0 )\n")
	    print("(From TH1 docs):\n")
            print("1) If ngroup is not an exact divider of the number of bins, the top limit of the rebinned histogram is reduced to the upper edge of the last bin that can make a complete group.")
            print("   The remaining bins are added to the overflow bin. Statistics will be recomputed from the new bin contents.")
            print("   If rebinning to one single bin, this might lead to an \"empty\" histogram (as everything will end up in the overflow bin)\n")
            print("2) If rebinning to variable bin size (i.e., xbins != 0), the bin edges specified in xbins should correspond to bin edges in the original histogram.")
            print("   If a bin edge in the new histogram is in the middle of a bin in the original histogram, all entries in the split bin in the original histogram will be transfered to the lower of the two possible bins in the new histogram. This is probably not what you want.")
            print("   Therefore, make sure the new binning is compatible w/ the binning of the input histogram!\n")

        if rebinlist:
            for key in self.histkeys:

	        # By construction, the following will be a list w/ N items, where:
	        #
	        # tokens[0] = efficiency ("Real","Fake"...)
	        # tokens[1] = lepton ("El","Mu"...)
	        # tokens[2] = variable ("Pt","Eta","Eta&&Pt"...)
	        # tokens[3] = process ("observed","expectedbkg"...)

	        tokens = key.split("_")

                if self.verbose:
	            print("Current tokens:")
	            print tokens

                if "&&" in tokens[2]:
                    print("Variable: {0}. Cannot rebin w/ variable bin size a TH2. Consider passing an input TH2 already made w/ appropriate binning...".format(tokens[2]))
                    continue

                for rebinitem in rebinlist:

	            rebinitem = rebinitem.split(",")

	            if ( tokens[0] in rebinitem ) and ( tokens[1] in rebinitem ) and ( tokens[2] in rebinitem ):

                        nbins = len(rebinitem[3:])-1

                        if self.debug:
                            print("\t===> Rebinning matching histograms with the following values:")
                            print "\tbin edges: ", rebinitem[3:], ", number of bins = ", nbins

                        bins = [ float(binedge) for binedge in rebinitem[3:] ]
	                arr_bins = array.array("d", bins)

                        self.numerator_hists[key]     = self.numerator_hists[key].Rebin( nbins, key, arr_bins )
                        self.denominator_hists[key]   = self.denominator_hists[key].Rebin( nbins, key, arr_bins )
                        self.antinumerator_hists[key] = self.antinumerator_hists[key].Rebin( nbins, key, arr_bins )

        if averagehistlist:

	    self.__averagehists = True

	    for key in self.histkeys:

		tokens = key.split("_")

		print("Current tokens:")
	        print tokens

		for averageitem in averagehistlist:

		    averageitem = averageitem.split(",")

		    if ( averageitem[0] == "ALL" ) or ( ( tokens[0] in averageitem ) and ( tokens[1] in averageitem ) and ( tokens[2] in averageitem ) ):

		        print("\t===> Taking average on whole bin range of matching histograms")

                        if not isinstance(self.numerator_hists[key],TH2):
                            nbins_numerator      = self.numerator_hists[key].GetNbinsX()
                            nbins_denominator    = self.denominator_hists[key].GetNbinsX()
                            nbins_antinumerator  = self.antinumerator_hists[key].GetNbinsX()
                            self.numerator_hists[key]     = self.numerator_hists[key].Rebin( nbins_numerator, key )
                            self.denominator_hists[key]   = self.denominator_hists[key].Rebin( nbins_denominator, key )
                            self.antinumerator_hists[key] = self.antinumerator_hists[key].Rebin( nbins_antinumerator, key )
                        else:
                            nbins_numerator_x      = self.numerator_hists[key].GetXaxis().GetNbins()
                            nbins_denominator_x    = self.denominator_hists[key].GetXaxis().GetNbins()
                            nbins_antinumerator_x  = self.antinumerator_hists[key].GetXaxis().GetNbins()
                            nbins_numerator_y      = self.numerator_hists[key].GetYaxis().GetNbins()
                            nbins_denominator_y    = self.denominator_hists[key].GetYaxis().GetNbins()
                            nbins_antinumerator_y  = self.antinumerator_hists[key].GetYaxis().GetNbins()
                            self.numerator_hists[key]     = self.numerator_hists[key].Rebin2D( nbins_numerator_x, nbins_numerator_y, key )
                            self.denominator_hists[key]   = self.denominator_hists[key].Rebin2D( nbins_denominator_x, nbins_denominator_y, key )
                            self.antinumerator_hists[key] = self.antinumerator_hists[key].Rebin2D( nbins_antinumerator_x, nbins_antinumerator_y, key )

        self.storeYields()


    def __getAverageHist__( self, histogram ):

        if not isinstance(histogram,TH2):
            nbins = histogram.GetNbinsX()
            h_avg = histogram.Rebin( nbins, histogram.GetName()+"_AVG" )
        else:
            nbinsx = histogram.GetXaxis().GetNbins()
            nbinsy = histogram.GetYaxis().GetNbins()
            h_avg = histogram.Rebin2D( nbinsx, nbinsy, histogram.GetName()+"_AVG" )

        return h_avg


    def __yields_and_integral__( self, histogram ):

        yields = [ histogram.GetBinContent(ibin) for ibin in range( 1, histogram.GetSize() ) ]

        integral = 0.0
        if isinstance(histogram,TH1) and not isinstance(histogram,TH2):
            integral = histogram.Integral(0,histogram.GetNbinsX()+1)
        else:
            integral = histogram.Integral(0,histogram.GetXaxis().GetNbins()+1,0,histogram.GetYaxis().GetNbins()+1)

        yields.append(integral)

        return yields


    def storeYields ( self ):

        # Make sure "numerator" (eg. TIGHT) bins have *never* more events
	# than "denominator" (eg LOOSE) bins. That is, set numerator=denominator in such cases.
	# Small fluctuations can happen b/c of subtraction and numerical precision (?)

        for key, numerator_hist in self.numerator_hists.iteritems():
             for ibin in range( 1, numerator_hist.GetSize() ):
	         numerator_yield   = numerator_hist.GetBinContent(ibin)
		 denominator_yield = self.denominator_hists[key].GetBinContent(ibin)
		 if numerator_yield > denominator_yield:
		     self.denominator_hists[key].SetBinContent(ibin,numerator_yield)

        # Reset to zero any bin which became negative after subtraction,
        # then get the yields for every bin and the integral, and store it in a list

        for key, hist in self.numerator_hists.iteritems():

            # print("Storing yields for N - key: {0} - integral: {1:.3f}".format(key,hist.Integral()))

            for ibin in range( 1, hist.GetSize() ):
                if hist.GetBinContent(ibin) < 0:
                    hist.SetBinContent(ibin,0.0)
            self.numerator_yields[key] = self.__yields_and_integral__(hist)

        for key, hist in self.denominator_hists.iteritems():

            # print("Storing yields for D - key: {0} - integral: {1:.3f}".format(key,hist.Integral()))

            for ibin in range( 1, hist.GetSize() ):
                if hist.GetBinContent(ibin) < 0:
                    hist.SetBinContent(ibin,0.0)
            self.denominator_yields[key] = self.__yields_and_integral__(hist)

        for key, hist in self.antinumerator_hists.iteritems():

            # print("Storing yields for AntiN - key: {0} - integral: {1:.3f}".format(key,hist.Integral()))

            for ibin in range( 1, hist.GetSize() ):
                if hist.GetBinContent(ibin) < 0:
                    hist.SetBinContent(ibin,0.0)
            self.antinumerator_yields[key] = self.__yields_and_integral__(hist)


    def sanityCheck( self ):
    	print("\n\n")
	for key, num_list in self.numerator_yields.iteritems():
    	     numerator   = num_list[-1]
    	     denominator = self.denominator_yields[key][-1]
    	     if denominator < numerator:
    		 print("WARNING! Histogram:{0} ==> D = {1:.3f}, N = {2:.3f}, !N = {3:.3f} ==> D < N ??".format(key, denominator, numerator, self.antinumerator_yields[key][-1] ))
                 print("Numerator : integral = {0}".format(self.numerator_yields[key][-1]))
                 print("\t" + " ; ".join( "({0},{1:.3f})".format(bin,value) for bin,value in enumerate(self.numerator_yields[key][:-1])))
                 print("Denominator : integral = {0}".format(self.denominator_yields[key][-1]))
                 print("\t" + " ; ".join( "({0},{1:.3f})".format(bin,value) for bin,value in enumerate(self.denominator_yields[key][:-1])))
                 print("AntiNumerator : integral = {0}".format(self.antinumerator_yields[key][-1]))
                 print("\t" + " ; ".join( "({0},{1:.3f})".format(bin,value) for bin,value in enumerate(self.antinumerator_yields[key][:-1])))


    def __manageBoundaries__( self, histogram ):

        if isinstance(histogram,TH1) and not isinstance(histogram,TH2):

            # 1. Merge OFlow bin in last visible bin

	    h_last_idx  = histogram.GetNbinsX()
	    h_oflow_idx = histogram.GetNbinsX()+1
	    h_last      = histogram.GetBinContent( h_last_idx )
	    h_oflow     = histogram.GetBinContent( h_oflow_idx )
	    h_last_err  = histogram.GetBinError( h_last_idx )
	    h_oflow_err = histogram.GetBinError( h_oflow_idx )

	    h_merged     = h_last + h_oflow
	    h_merged_err = math.sqrt( pow(h_last_err,2.0) + pow(h_oflow_err,2.0) )

	    histogram.SetBinContent( h_last_idx, h_merged )
	    histogram.SetBinError( h_last_idx, h_merged_err )

	    # 2. Set OFlow bin content and error equal to the just modified last visible bin content and error

	    histogram.SetBinContent( h_oflow_idx, histogram.GetBinContent( h_last_idx ) )
	    histogram.SetBinError( h_oflow_idx, histogram.GetBinError( h_last_idx )  )

        else:

            # NB: for a TH2, the overflow is the top-right corner of the rectangular grid defining the histogram itself.
            # Therefore, we need to apply the same strategy of TH1 for each x and y 1D slice.

            # 1.A Loop over x-axis slices, and fix the y-projetcion overflows

            ybin_last_idx  = histogram.GetYaxis().GetNbins()
            ybin_oflow_idx = histogram.GetYaxis().GetNbins()+1

            for xbin in range(1,histogram.GetXaxis().GetNbins()+2):

                y_last      = histogram.GetBinContent( xbin, ybin_last_idx )
                y_oflow     = histogram.GetBinContent( xbin, ybin_oflow_idx )
                y_last_err  = histogram.GetBinError( xbin, ybin_last_idx )
                y_oflow_err = histogram.GetBinError( xbin, ybin_oflow_idx )

                y_merged     = y_last + y_oflow
                y_merged_err = math.sqrt( pow(y_last_err,2.0) + pow(y_oflow_err,2.0) )

                histogram.SetBinContent( xbin, ybin_last_idx, y_merged )
                histogram.SetBinError( xbin, ybin_last_idx, y_merged_err )

                histogram.SetBinContent( xbin, ybin_oflow_idx, histogram.GetBinContent( ybin_last_idx ) )
                histogram.SetBinError( xbin, ybin_oflow_idx, histogram.GetBinError( ybin_last_idx ) )

            # 1.B Loop over y-axis slices, and fix the x-projetcion overflows

            xbin_last_idx  = histogram.GetXaxis().GetNbins()
            xbin_oflow_idx = histogram.GetXaxis().GetNbins()+1

            for ybin in range(1,histogram.GetYaxis().GetNbins()+2):

                x_last      = histogram.GetBinContent( xbin_last_idx, ybin )
                x_oflow     = histogram.GetBinContent( xbin_oflow_idx, ybin )
                x_last_err  = histogram.GetBinError( xbin_last_idx, ybin )
                x_oflow_err = histogram.GetBinError( xbin_oflow_idx, ybin )

                x_merged     = x_last + x_oflow
                x_merged_err = math.sqrt( pow(x_last_err,2.0) + pow(x_oflow_err,2.0) )

                histogram.SetBinContent( xbin_last_idx, ybin, x_merged )
                histogram.SetBinError( xbin_last_idx, ybin, x_merged_err )

                histogram.SetBinContent( xbin_oflow_idx, ybin, histogram.GetBinContent( xbin_last_idx ) )
                histogram.SetBinError( xbin_oflow_idx, ybin, histogram.GetBinError( xbin_last_idx ) )

        return histogram


    def __makeProjectionHists__( self, histogram, eventset ):

        for key in sorted(self.histkeys):

            if not "&&" in key or "_proj" in key: continue

            tokens = key.split('_')

            vars2D = tokens[2].split('&&')
            varX   = vars2D[0]
            varY   = vars2D[1]

            tot_integral_proj_x = tot_integral_proj_y = 0.0

            if self.verbose:
                print("\tkey: {0}".format(key))
                print("\n\t\tvarX: {0} - varY: {1}\n".format(varX, varY))
                print("\t\tprojecting TH2 onto: {0}".format(varX))

            for i,ybin in enumerate( range(1,histogram.GetYaxis().GetNbins()+2) ):

                ybin_upedge  = histogram.GetYaxis().GetBinUpEdge(ybin)
                ybin_lowedge = histogram.GetYaxis().GetBinLowEdge(ybin)

                if self.verbose: print("\t\t\tybin: {0} - edges: [{1},{2}]".format(ybin,ybin_lowedge,ybin_upedge))

                proj_x_key = "{0}_proj{1}_{2}".format(key,varX,ybin)

                proj_x = histogram.ProjectionX(proj_x_key+"_"+eventset,ybin,ybin)
                proj_x.SetLineColor(i+1)
                proj_x.SetMarkerColor(i+1)

                integral = proj_x.Integral(0,proj_x.GetSize()-1)
                tot_integral_proj_x += integral

                if self.verbose: print("\t\t\tProjection key: {0} - Integral(): {1:.3f}".format(proj_x_key,integral))

                proj_x.SetDirectory(0)
                self.__fillNDHistDicts__(proj_x_key, eventset, proj_x )

            if self.verbose: print("\t\tprojecting TH2 onto: {0}".format(varY))

            for i,xbin in enumerate( range(1,histogram.GetXaxis().GetNbins()+2) ):

                xbin_upedge  = histogram.GetXaxis().GetBinUpEdge(xbin)
                xbin_lowedge = histogram.GetXaxis().GetBinLowEdge(xbin)

                if self.verbose: print("\t\t\txbin: {0} - edges: [{1},{2}]".format(xbin,xbin_lowedge,xbin_upedge))

                proj_y_key = "{0}_proj{1}_{2}".format(key,varY,xbin)

                proj_y = histogram.ProjectionY(proj_y_key+"_"+eventset,xbin,xbin)
                proj_y.SetLineColor(i+1)
                proj_y.SetMarkerColor(i+1)

                integral = proj_y.Integral(0,proj_y.GetSize()-1)
                tot_integral_proj_y += integral

                if self.verbose: print("\t\t\tProjection key: {0} - Integral(): {1:.3f}".format(proj_y_key,integral))

                proj_y.SetDirectory(0)
                self.__fillNDHistDicts__(proj_y_key, eventset, proj_y )

            if self.verbose: print("\n\t\tTot. integral projX: {0:.3f} - projY: {1:.3f}".format(tot_integral_proj_x,tot_integral_proj_y))


    # NB: When computing an efficiency, need to make sure the errors are computed correctly!
    # In this case, numerator and denominator are not independent sets of events! The efficiency is described by a binomial PDF.

    def computeEfficiencies( self, variation ):

        if self.closure and variation != "nominal":
            print("\n\tCLOSURE: do nominal only...")
            return
        if self.nosub and variation != "nominal":
            print("\n\No subtraction activated: do nominal only...")
            return

        for key in sorted(self.histkeys):
            if variation != "nominal": continue
            self.__makeProjectionHists__( self.numerator_hists[key], "N" )
            self.__makeProjectionHists__( self.denominator_hists[key], "D" )
            self.__makeProjectionHists__( self.antinumerator_hists[key], "AntiN" )

        #self.storeYields()

        print("\nCalculating EFFICIENCIES - {0}...".format(variation))

        nominal_key = None

        # NB: here is crucial to loop over the alphabetically-sorted keys!

        for key in sorted(self.histkeys):

            tokens = key.split("_")

            if any( sys in tokens for sys in ["TTV","VV","OtherPromptSS","FakesOS"] ) :
                if not variation == "ND": continue
                if tokens[-1].isdigit(): continue
            if any( sys in tokens for sys in ["QMisID"] ):
                if not variation in ["N","D"]: continue

	    min_tokens = True
	    if self.closure or self.nosub:
	        min_tokens = ( len(tokens) >= 4 )
	    else:
	        min_tokens = ( len(tokens) >= 5 )

            if ( not min_tokens or tokens[3] not in self.__processes ): continue

            if variation == "nominal" :
                if ( "_proj" not in key ) and len(tokens) > 5 : continue
                nominal_key = key

            # If checking systematics, need to store the nominal key first

            if ( variation != "nominal" and len(tokens) == 5 ):
                nominal_key = key
                continue

            # Extension for projection histograms

            proj_extension = ""
            for i, token in enumerate(tokens):
                if "proj" in token:
                    proj_extension = "_" + "_".join(tokens[i:])
                    break

            if self.debug:
                print("\n\t*****************************************************")
                print("\n\tCalculating efficiency for:\n\tnominal key: {0}\n\tvar key: {1}\n".format(nominal_key,key))

            # Define pass (N) and total (D=N+!N)

            h_pass = h_tot = None

            key_N = key_AntiN = nominal_key
            append = ""

            if variation == "N":
                append = "N"
                key_N     = key
                key_AntiN = nominal_key
            elif variation == "D":
                append = "D"
                key_N     = nominal_key
                key_AntiN = key
            elif variation == "ND":
                append = "ND"
                key_N     = key
                key_AntiN = key

            # Crucial to clone histograms here!

            h_pass = self.numerator_hists[key_N].Clone()
            h_tot  = h_pass.Clone()
            h_tot.Add(self.antinumerator_hists[key_AntiN].Clone())

	    # ----------------------------------------------------------------------------------

	    # Wrap around boundaries to take overflows into account

            h_pass = self.__manageBoundaries__(h_pass)
            h_tot  = self.__manageBoundaries__(h_tot)

            # ----------------------------------------------------------------------------------

            h_pass_AVG = self.__getAverageHist__(h_pass)
            h_tot_AVG  = self.__getAverageHist__(h_tot)

            h_pass_AVG = self.__manageBoundaries__(h_pass_AVG)
            h_tot_AVG  = self.__manageBoundaries__(h_tot_AVG)

            # ----------------------------------------------------------------------------------

	    ratiolist = []
	    for idx, elem in enumerate(self.numerator_yields[key_N]):
	        n = elem
		d = n + self.antinumerator_yields[key_AntiN][idx]
		r = -1 # just an unphysical value
		if d:
		    r = n/d
		ratiolist.append(r)

            if self.verbose:
                print("\t-----------------------------------------------------")
                print("\tHistogram: {0}\n".format(key))
                print("\tNumerator : integral = {0:.3f}".format(self.numerator_yields[key_N][-1]))
                print("\t\t" + " ; ".join( "({0},{1:.3f})".format(bin,value) for bin,value in enumerate(self.numerator_yields[key_N][:-1])))
                print("\tDenominator  (N+!N): integral = {0:.3f}".format( self.numerator_yields[key_N][-1] + self.antinumerator_yields[key_AntiN][-1]))
                print("\t\t" + " ; ".join( "({0},{1:.3f})".format(bin,value) for bin,value in enumerate( [ sum(x) for x in zip(  self.numerator_yields[key_N][:-1], self.antinumerator_yields[key_AntiN][:-1] ) ] ) ) )
                print("\tAntiNumerator: integral = {0:.3f}".format(self.antinumerator_yields[key_AntiN][-1]))
                print("\t\t" + " ; ".join( "({0},{1:.3f})".format(bin,value) for bin,value in enumerate(self.antinumerator_yields[key_AntiN][:-1])))
                print("\tRatio N/D (<eff>): integral = {0:.3f}".format(ratiolist[-1]))
                print("\t\t" + " ; ".join( "({0},{1:.3f})".format(bin,value) for bin,value in enumerate(ratiolist[:-1])))
                print("\t-----------------------------------------------------")
                print("\t")

            # 1.
            #
	    # The TH1::Divide method with the option "B" calculates binomial errors using the "normal" approximation
            # (NB: the approximation fails when eff = 0 or 1. In such cases, TEfficiency or TGraphAsymmErrors should be used, since they know how to handle such cases)

	    if self.closure or self.nosub:
                key_heff = "_".join( (tokens[0],tokens[1],tokens[2],"Efficiency",tokens[3]+proj_extension,append) )
            else:
	        key_heff = "_".join( (tokens[0],tokens[1],tokens[2],"Efficiency",tokens[3]+proj_extension,tokens[4],append) )

            if key_heff.endswith("_"):
                key_heff = key_heff[:-1]
            if len(tokens) > 6:
                key_heff = key_heff + "_" + "_".join( ("{0}".format(other_tokens) for other_tokens in tokens[5:]) )

            h_efficiency  = h_pass.Clone(key_heff)
            h_efficiency.Divide(h_pass,h_tot,1.0,1.0,"B")

            h_efficiency_AVG  = h_pass_AVG.Clone(key_heff+"_AVG")
            h_efficiency_AVG.Divide(h_pass_AVG,h_tot_AVG,1.0,1.0,"B")

            # 2.
            # The TEfficiency class handles the special cases not covered by TH1::Divide

	    t_efficiency = None
            key_teff = key_heff.replace( "Efficiency", "TEfficiency" )

	    if isinstance(h_pass,TH1D) and TEfficiency.CheckConsistency(h_pass, h_tot,"w"):
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

	        t_efficiency.SetStatisticOption(TEfficiency.kBUniform)
                t_efficiency.SetPosteriorMode()

            # 3.
            #
            # Calculate efficiency using TGraphAsymmErrors
            # (uses a Bayesian uniform prior beta(1,1) for the efficiency with 68% CL. It handles the errors in case eff is 0 or 1)

            key_geff = key_heff
            g_efficiency = None
            if isinstance(h_efficiency,TH1D):
                g_efficiency = TGraphAsymmErrors(h_efficiency)
                g_efficiency.Divide(h_pass,h_tot,"cl=0.683 b(1,1) mode")

            # Save the efficiencies in the proper dictionaries

            self.histefficiencies[key_heff]         = h_efficiency
            self.histefficiencies[key_heff+"_AVG"]  = h_efficiency_AVG
            self.graphefficiencies[key_geff]        = g_efficiency
       	    self.tefficiencies[key_teff]            = t_efficiency

            print ("\tkey for efficiency: {0}".format(key_heff))
            print ("\tkey for efficiency (AVG): {0}".format(key_heff+"_AVG"))

    def computeFactors( self, variation ):

	print("\nCalculating FACTORS...")

        nominal_key = None

        for key in sorted(self.histkeys):

	    if self.closure and variation != "nominal":
	        print("\n\tCLOSURE: do nominal only...")
		return
	    if self.nosub and variation != "nominal":
	        print("\n\No subtraction activated: do nominal only...")
		return

	    tokens = key.split("_")

	    min_tokens = True
	    if self.closure or self.nosub:
	        min_tokens = ( len(tokens) >= 4 )
	    else:
	        min_tokens = ( len(tokens) >= 5 )

            if ( not min_tokens or tokens[3] not in self.__processes ): continue

            if ( variation == "nominal"):
                if ( len(tokens) > 5 ): continue
                nominal_key = key

            # If checking systematics, need to store the nominal key

            if ( variation != "nominal" and len(tokens) == 5 ):
                nominal_key = key
                continue

            if self.debug:
                print("\n\t*****************************************************")
                print("\n\tCalculating fake factor for:\n\tnominal key: {0}\n\tvar key: {1}\n".format(nominal_key,key))

	    # Define numerator (pass) and denominator (not-pass)
	    # These are two set of *independent* events, so just doing the hist ratio will be ok.

            h_pass = h_notpass = None

            append = ""
            if variation == "nominal":
                h_pass     = self.numerator_hists[nominal_key]
                h_notpass  = self.antinumerator_hists[nominal_key]
            elif variation == "N":
                append = "N"
                h_pass     = self.numerator_hists[key]
                h_notpass  = self.antinumerator_hists[nominal_key]
            elif variation == "D":
                append = "D"
                h_pass     = self.numerator_hists[nominal_key]
                h_notpass  = self.antinumerator_hists[key]

	    # ----------------------------------------------------------------------------------

	    # Wrap around boundaries to take overflows into account

            h_pass    = self.__manageBoundaries__(h_pass)
            h_notpass = self.__manageBoundaries__(h_notpass)

            # ----------------------------------------------------------------------------------

	    ratiolist = []
	    for idx, elem in enumerate(self.numerator_yields[key]):
	        n = elem
		d = self.antinumerator_yields[key][idx]
		r = -1 # just an unphysical value
		if d:
		    r = n/d
		ratiolist.append(r)

	    if self.verbose:
                print("\t-----------------------------------------------------")
                print("\tHistogram: {0}\n".format(key))
                print("\tNumerator : integral = {0}".format(self.numerator_yields[key][-1]))
                print("\t\t" + " ; ".join( "({0},{1:.3f})".format(bin,value) for bin,value in enumerate(self.numerator_yields[key][:-1])))
                print("\tAntiNumerator: integral = {0}".format(self.antinumerator_yields[key][-1]))
                print("\t\t" + " ; ".join( "({0},{1:.3f})".format(bin,value) for bin,value in enumerate(self.antinumerator_yields[key][:-1])))
                print("\tRatio N/!N: integral = {0}".format(ratiolist[-1]))
                print("\t\t" + " ; ".join( "({0},{1:.3f})".format(bin,value) for bin,value in enumerate(ratiolist[:-1])))
                print("\t-----------------------------------------------------")
                print("\t")

	    if self.closure or self.nosub:
                key_hfactor = "_".join( (tokens[0],tokens[1],tokens[2],"Factor",tokens[3],append) )
            else:
	        key_hfactor = "_".join( (tokens[0],tokens[1],tokens[2],"Factor",tokens[3],tokens[4],append) )

            if key_hfactor.endswith("_"):
                key_hfactor = key_hfactor[:-1]

            if len(tokens) > 6:
                key_hfactor = key_hfactor + "_" + "_".join( ("{0}".format(other_tokens) for other_tokens in tokens[5:]) )

            h_factor  = h_pass.Clone(key_hfactor)
            h_factor.Divide(h_pass,h_notpass)

            # Save the factors in the proper dictionary

            self.histfactors[key_hfactor] = h_factor

            print ("\tkey for factor: {0}".format(key_hfactor))


    def saveOutputs( self, filename="LeptonEfficiencies", outputpath=None ):

	if not outputpath:
	    outputpath = self.__outputpath

	self.__outputpath = outputpath

        if self.__averagehists:
	    filename += "Avg"

        if self.nosub:
            filename += "_NO_SUB"

        print("\nStoring output files:\n{0}\n{1}\nin directory: {2}\n".format(filename+".root",filename+".txt",os.path.abspath(outputpath)))

        self.__outputfile_yields = open(self.__outputpath+"/"+filename+".txt","w")
        self.__outputfile_yields.write( "Efficiencies/Factors for Fake Factor amd Matrix Method\n")

        fileopt = "RECREATE"
        if args.update:
            fileopt = "UPDATE"

        self.__outputfile = TFile(self.__outputpath+"/"+filename+".root",fileopt)
        self.__outputfile.cd()

        if self.debug: print("\nSaving histograms:")
    	for key, h in sorted(self.histefficiencies.iteritems()):
    	    if not h: continue
    	    if self.debug: print("\t{0}".format(key))

            # Make sure histogram name contains symbols that ROOT can parse ...

            hname = h.GetName()
            if "&&" in hname:
                hname = hname.replace("&&","_VS_")
                h.SetName(hname)
            if args.do3L:
                h.SetName("3L_"+hname)
    	    h.Write()

    	    eff=[]
            is2DHist = ( isinstance(h,TH2) )
    	    for ibin in range( 1, h.GetSize() ):
                if not is2DHist:
                    myset = [ ibin, h.GetBinLowEdge(ibin), h.GetBinLowEdge(ibin+1), h.GetBinContent(ibin), h.GetBinError(ibin)]
                else:
                    myset = [ ibin, h.GetBinContent(ibin), h.GetBinError(ibin)]
    	    	eff.append( myset )
    	    self.__outputfile_yields.write("%s:\n" %(key) )
    	    for myset in eff:
                if not is2DHist:
                    self.__outputfile_yields.write("{ %s }; \n" %( "Bin nr: " + str(myset[0]) + " [" + str(round(myset[1],3)) + "," + str(round(myset[2],3)) + "], efficiency (from TH1::Divide(\"B\")) = " + str(round(myset[3],3)) + " +- " + str(round(myset[4],3)) ) )
                else:
                    self.__outputfile_yields.write("{ %s }; \n" %( "Bin nr: " + str(myset[0]) + ", efficiency (from TH1::Divide(\"B\")) = " + str(round(myset[1],3)) + " +- " + str(round(myset[2],3)) ) )

    	# for key, t in sorted(self.tefficiencies.iteritems()):
    	#     if not t: continue
    	#     if self.debug: print("\nSaving TEfficiency: {0}".format(key))
    	#     t.Write()
    	#     teff=[]
    	#     for ibin in range( 1, t.GetTotalHistogram().GetSize() ):
    	#     	myset = [ ibin, t.GetTotalHistogram().GetBinLowEdge(ibin), t.GetTotalHistogram().GetBinLowEdge(ibin+1), t.GetEfficiency(ibin), t.GetEfficiencyErrorUp(ibin), t.GetEfficiencyErrorLow(ibin)]
    	#     	teff.append( myset )
    	#     self.__outputfile_yields.write("%s:\n" %(key) )
    	#     for myset in teff:
    	#     	self.__outputfile_yields.write("{ %s }; \n" %( "Bin nr: " + str(myset[0]) + " [" + str(round(myset[1],3)) + "," + str(round(myset[2],3)) + "], efficiency (from TEfficiency) = " + str(round(myset[3],3)) + " + " + str(round(myset[4],3)) + " - " + str(round(myset[5],3)) ) )

    	# for key, g in sorted(self.graphefficiencies.iteritems()):
    	#     if not g: continue
    	#     if self.debug: print("\nSaving graph: {0}".format(key))
    	#     g.Write()
    	#     geff=[]
    	#     for ipoint in range( 0, g.GetN()+1 ):
    	#     	x = Double(0)
    	#     	y = Double(0)
    	#     	g.GetPoint(ipoint,x,y)
    	#     	myset = [ ipoint+1, y, g.GetErrorYhigh(ipoint), g.GetErrorYlow(ipoint) ]
    	#     	geff.append( myset )
    	#     self.__outputfile_yields.write("%s:\n" %(key) )
    	#     for myset in geff:
    	#     	self.__outputfile_yields.write("{ %s }; \n" %( "Bin nr: " + str(myset[0]) + ", efficiency (from TGraphAsymmErrors) = " + str(round(myset[1],3)) + " + " + str(round(myset[2],3)) + " - " + str(round(myset[3],3)) ) )

	if self.factors:
    	    for key, h in sorted(self.histfactors.iteritems()):
    	    	if not h: continue
    	    	# if self.debug: print("\nSaving histogram: {0}".format(key))
    	    	h.Write()
    	    	eff=[]
    	    	for ibin in range( 1, h.GetSize() ):
    	    	    myset = [ ibin, h.GetBinLowEdge(ibin), h.GetBinLowEdge(ibin+1), h.GetBinContent(ibin), h.GetBinError(ibin)]
    	    	    eff.append( myset )
    	    	self.__outputfile_yields.write("%s:\n" %(key) )
    	    	for myset in eff:
    	    	    self.__outputfile_yields.write("{ %s }; \n" %( "Bin nr: " + str(myset[0]) + " [" + str(round(myset[1],3)) + "," + str(round(myset[2],3)) + "], factor (from TH1::Divide()) = " + str(round(myset[3],3)) + " +- " + str(round(myset[4],3)) ) )

	self.__outputfile_yields.close()
        self.__outputfile.Write()
        self.__outputfile.Close()

    def __getLimits__( self, histlist, scale_up=1.0, scale_dn=1.0, shift_up=0.0, shift_dn=0.0, ratio=False ):
        ymin = 1e9
	ymax = 0.0
	for h in histlist:
	    this_ymin = h.GetBinContent( h.GetMinimumBin() )
	    this_ymax = h.GetBinContent( h.GetMaximumBin() )
	    if ratio:
		shifted_dn = [ ( h.GetBinContent(bin) - h.GetBinError(bin)/2.0 ) if ( h.GetBinError(bin) ) else 1 for bin in range(1,h.GetSize()) ]
	        shifted_up = [ ( h.GetBinContent(bin) + h.GetBinError(bin)/2.0 ) if ( h.GetBinError(bin) ) else 1 for bin in range(1,h.GetSize()) ]
		this_ymin = min(shifted_dn)
		this_ymax = max(shifted_up)

	    if this_ymin < ymin:
	        ymin = this_ymin
	    if this_ymax >= ymax:
	        ymax = this_ymax

	if not ratio:
	    if ymin < 0: ymin = 0.0
	    if ymax > 1: ymax = 1.0
	else:
	   if  abs( ymax - 1 ) <= 0.05:
	       ymax = 1.05
	       shift_up = 0
	   if  abs( ymin - 1 ) <= 0.05:
	       ymin = 0.95
	       shift_dn = 0

	return scale_dn * ( ymin - shift_dn ), scale_up * ( ymax + shift_up )


    def __save1DProjections__( self, var, hist2Dkey, canvasname, savepath ):

        if "AVG" in hist2Dkey: return

        c = TCanvas(canvasname+"_Projections","Projections",50,50,1000,600)
        c.SetFrameFillColor(0)
        c.SetFrameFillStyle(0)
        c.SetFrameBorderMode(0)

        c.Divide(2,1)

        legendx = TLegend(0.45,0.25,0.925,0.55) # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
        legendx.SetBorderSize(0)     # no border
        legendx.SetFillStyle(0)      # Legend transparent background
        legendx.SetTextSize(0.035)   # Increase entry font size!
        legendx.SetTextFont(42)      # Helvetica

        legendy = TLegend(0.45,0.25,0.925,0.55) # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
        legendy.SetBorderSize(0)     # no border
        legendy.SetFillStyle(0)      # Legend transparent background
        legendy.SetTextSize(0.035)   # Increase entry font size!
        legendy.SetTextFont(42)      # Helvetica

        normfactor = self.histefficiencies[hist2Dkey+"_AVG"].GetBinContent(1,1)

    	for key, h in sorted(self.histefficiencies.iteritems()):

            vars2D = var.split("_VS_")

            if "AVG" in key: continue
            if not all( v in key for v in vars2D ): continue
            if not "proj" in key: continue

            tokens = key.split("_")

            slicefactor = self.histefficiencies[key+"_AVG"].GetBinContent(1)

            scalefactor = 0.0
            if slicefactor: scalefactor = normfactor/slicefactor

            if self.verbose: print("\tkey: {0} - <TOT eff> = {1:.3f}, <SLICE eff> = {2:.3f}, scale factor = {3:.3f}".format(key,normfactor,slicefactor,scalefactor))

            h.Scale(scalefactor)

            h.SetMaximum(1.0)

            if self.verbose: print("\t\t==> new integral = {0:.3f}".format(h.Integral(0,h.GetSize()-1)))

            varX = "proj" + vars2D[0]
            varY = "proj" + vars2D[1]

            if varX in key:

                c.cd(1)

                slice_lowedge =  self.histefficiencies[hist2Dkey].GetYaxis().GetBinLowEdge(int(tokens[-1]))
                slice_upedge =  self.histefficiencies[hist2Dkey].GetYaxis().GetBinUpEdge(int(tokens[-1]))

                legendx.AddEntry(h,"{0} - [{1},{2}]".format(vars2D[1],slice_lowedge,slice_upedge), "P")

                if not gPad.GetListOfPrimitives().GetSize():
                    h.Draw("E0")
                else:
                    h.Draw("E0 SAME")

                legendx.Draw()


            if varY in key:

                c.cd(2)

                slice_lowedge =  self.histefficiencies[hist2Dkey].GetXaxis().GetBinLowEdge(int(tokens[-1]))
                slice_upedge =  self.histefficiencies[hist2Dkey].GetXaxis().GetBinUpEdge(int(tokens[-1]))

                legendy.AddEntry(h,"{0} - [{1},{2}]".format(vars2D[0],slice_lowedge,slice_upedge), "P")

                if not gPad.GetListOfPrimitives().GetSize():
                    h.Draw("E0")
                else:
                    h.Draw("E0 SAME")

                legendy.Draw()

        for extension in self.extensionlist:
            c.SaveAs(savepath+"/BasicPlots/"+canvasname+"_Projections."+extension)


    def plotMaker( self ):

	gROOT.SetBatch(True)

        proc = ("observed","expectedbkg")[bool(self.closure)]

        sub = ("w/","w/o")[bool(self.nosub)]
        proc_dict = {"observed":"Data ({0} subtraction)".format(sub), "expectedbkg":"simulation"}

        savepath = self.__outputpath+"/EfficiencyPlots"+("","_Avg")[self.__averagehists]

        if self.nosub:
            savepath += "_NO_SUB"

        if self.triggerEff:
	    savepath += "_TriggerEff"

	if not os.path.exists(savepath+"/BasicPlots"):
	    os.makedirs(savepath+"/BasicPlots")

	if self.triggerEff:
            trigeff_file = TFile(savepath+"/BasicPlots/RealFake_"+self.triggerEff+"_TriggerEfficiency.root","RECREATE")

	if self.probeAssignEff:
            probeassigneff_file = TFile(savepath+"/BasicPlots/RealFake_ProbeAssignEfficiency.root","RECREATE")

        leg_ATLAS  = TLatex()
        leg_lumi   = TLatex()
        leg_ATLAS.SetTextSize(0.03)
        leg_ATLAS.SetNDC()
        leg_lumi.SetTextSize(0.03)
        leg_lumi.SetNDC()

        for var in self.__variables:

            is2DHist = ( "&&" in var )

            if is2DHist:

                for lep in self.__leptons:

                    legend = TLegend(0.45,0.5,0.925,0.8) # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
                    legend.SetBorderSize(0)     # no border
                    legend.SetFillStyle(0)      # Legend transparent background
                    legend.SetTextSize(0.035)   # Increase entry font size!
                    #legend.SetTextFont(42)      # Helvetica

                    legend.SetHeader(self.leptons_full[lep])

                    for idx_eff, eff in enumerate(self.__efficiencies):

                        c = TCanvas("c_"+lep+"_"+eff,"Efficiencies",50,50,800,600)
                        c.SetFrameFillColor(0)
                        c.SetFrameFillStyle(0)
                        c.SetFrameBorderMode(0)

                        key  = "_".join((eff,lep,var,"Efficiency",proc,"sub"))
                        if self.closure or self.nosub:
                            key  = "_".join((eff,lep,var,"Efficiency",proc))

                        print("\tplotting histogram: {0}".format(key))

                        hist = self.histefficiencies[key]

                        set_fancy_2D_style()
                        gPad.SetRightMargin(0.2)
                        hist.Draw("COLZ1 text")

                        # legend.AddEntry(hist,eff+" - "+proc_dict[proc], "P")
                        # legend.Draw()
                        leg_ATLAS.DrawLatex(0.2,0.82,"#bf{#it{ATLAS}} Work In Progress");
                        leg_lumi.DrawLatex(0.2,0.77,"#sqrt{{s}} = 13 TeV, #int L dt = {0:.1f} fb^{{-1}}".format(self.lumi));

                        var = var.replace("&&","_VS_")
                        canvas_filename = "_".join((eff,lep,var,"Efficiency",proc))

                        c.Update()

                        for extension in self.extensionlist:
                            c.SaveAs(savepath+"/BasicPlots/"+canvas_filename+"."+extension)

                        self.__save1DProjections__(var, key, canvas_filename, savepath)

                        c.Clear()
                        c_avg = c.Clone(c.GetName()+"_AVG")
                        hist_avg = self.histefficiencies[key+"_AVG"]

                        set_fancy_2D_style()
                        gPad.SetRightMargin(0.2)
                        hist_avg.Draw("COLZ1 text")

                        c_avg.Update()

                        for extension in self.extensionlist:
                            c_avg.SaveAs(savepath+"/BasicPlots/"+canvas_filename+"_AVG."+extension)


            else:

                for lep in self.__leptons:

                    legend = TLegend(0.45,0.5,0.925,0.8) # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
                    legend.SetBorderSize(0)     # no border
                    legend.SetFillStyle(0)      # Legend transparent background
                    legend.SetTextSize(0.035)   # Increase entry font size!
                    #legend.SetTextFont(42)      # Helvetica

                    legend.SetHeader(self.leptons_full[lep])

                    c = TCanvas("c_"+lep,"Efficiencies")
                    c.SetFrameFillColor(0)
                    c.SetFrameFillStyle(0)
                    c.SetFrameBorderMode(0)

                    for idx_eff, eff in enumerate(self.__efficiencies):

                        key  = "_".join((eff,lep,var,"Efficiency",proc,"sub"))
                        if self.closure or self.nosub:
                            key  = "_".join((eff,lep,var,"Efficiency",proc))

                        print("\tplotting histogram: {0}".format(key))

                        hist     = self.histefficiencies[key]
                        hist_AVG = self.histefficiencies[key+"_AVG"]

                        hist.GetYaxis().SetRangeUser(0,1)
                        hist.GetYaxis().SetTitle("#varepsilon")
                        hist.GetXaxis().SetTitleOffset(1.0)
                        hist.GetYaxis().SetTitleOffset(1.0)
                        hist.SetLineStyle(1)
                        hist.SetMarkerStyle(kCircle)

                        hist_AVG.SetLineStyle(2)

                        if not idx_eff:
                            hist.SetLineColor(kBlue)
                            hist_AVG.SetLineColor(kBlue)
                            hist.SetMarkerColor(kBlue)
                        else:
                            hist.SetLineColor(kOrange+7)
                            hist_AVG.SetLineColor(kOrange+7)
                            hist.SetMarkerColor(kOrange+7)

                        if not idx_eff:
                            hist.Draw("E0")
                            hist_AVG.Draw("HIST SAME")
                        else:
                            hist.Draw("E0,SAME")
                            hist_AVG.Draw("HIST SAME")

                        legend.AddEntry(hist,eff+" - "+proc_dict[proc], "P")
                        legend.AddEntry(hist_AVG,eff+" - <#varepsilon> - "+proc_dict[proc], "L")

                        if self.triggerEff:
                            copy_hist_name = hist.GetName()
                            copy_hist_name = copy_hist_name.replace("Efficiency",self.triggerEff+"_TriggerEfficiency")
                            copy_hist = hist.Clone(copy_hist_name)
                            trigeff_file.cd()
                            copy_hist.Write()

                        if self.probeAssignEff:
                            copy_hist_name = hist.GetName()
                            copy_hist_name = copy_hist_name.replace("Efficiency","ProbeAssignEfficiency")
                            copy_hist = hist.Clone(copy_hist_name)
                            probeassigneff_file.cd()
                            copy_hist.Write()

                    legend.Draw()
                    leg_ATLAS.DrawLatex(0.6,0.35,"#bf{#it{ATLAS}} Work In Progress");
                    leg_lumi.DrawLatex(0.6,0.27,"#sqrt{{s}} = 13 TeV, #int L dt = {0:.1f} fb^{{-1}}".format(self.lumi));

                    canvas_filename = "_".join(("RealFake",lep,var,"Efficiency",proc))

                    if self.triggerEff:
                        canvas_filename = canvas_filename.replace("Efficiency",self.triggerEff+"_TriggerEfficiency")

                    if self.probeAssignEff:
                        canvas_filename = canvas_filename.replace("Efficiency","ProbeAssignEfficiency")

                    for extension in self.extensionlist:
                        c.SaveAs(savepath+"/BasicPlots/"+canvas_filename+"."+extension)


    def plotMakerSys( self ):

	proc = ("observed","expectedbkg")[bool(self.closure)]

        proc_dict = {"observed":"Data (w/ subtraction)", "expectedbkg":"simulation"}

	savepath = self.__outputpath+"/EfficiencyPlots"+("","_Avg")[self.__averagehists]

        if not os.path.exists(savepath+"/SplitSys"):
	    os.makedirs(savepath+"/SplitSys")
        if not os.path.exists(savepath+"/CombinedSys"):
	    os.makedirs(savepath+"/CombinedSys")

        for var in self.__variables:

            is2DHist = ( "&&" in var )
            if is2DHist: continue

  	    for lep in self.__leptons:

	        for idx_eff, eff in enumerate(self.__efficiencies):

	            key_nominal = "_".join( (eff,lep,var,"Efficiency",proc) )
	            if not self.nosub:
	                key_nominal += "_sub"

	            hist_nominal = self.histefficiencies[key_nominal]
          	    hist_nominal.SetLineStyle(2)
		    hist_nominal.SetLineWidth(2)
          	    hist_nominal.SetLineColor(kBlack)
	  	    hist_nominal.SetMarkerStyle(kFullCircle)
          	    hist_nominal.SetMarkerColor(kBlack)

		    # Store histograms in a list for future convenience

		    histlist = [hist_nominal]

		    if self.verbose:
    		        print("plotting nominal histogram: {0}".format(key_nominal))

	            hists_sys_numerator   = []
	            hists_sys_denominator = []
	            hists_sys_numden      = []

  	     	    for sys in self.__systematics[1:]:

                        if sys == "QMisID" and not all( s in key_nominal for s in ["Fake","El"] ):
                            print("Skipping {0} systematics for {1}".format(sys,key_nominal))
                            continue
                        if sys in ["TTV","VV","OtherPromptSS"] and not ( "Fake" in key_nominal ):
                            print("Skipping {0} systematics for {1}".format(sys,key_nominal))
                            continue
                        if sys == "FakesOS" and not ( "Real" in key_nominal ):
                            print("Skipping {0} systematics for {1}".format(sys,key_nominal))
                            continue

	     		for sysdir in self.__systematicsdirections[1:]:

                            if sys in ["QMisID"]: # Systematics which are uncorrelated among bins, and between N, D subsets
                                for bin in range( 1, hist_nominal.GetSize()):

                                    keyappend_num   = "_".join(("N",sys,sysdir,str(bin)))
                                    keyappend_denom = "_".join(("D",sys,sysdir,str(bin)))

                                    if self.verbose:
                                        print "\tstoring efficiency: ", "_".join((key_nominal,keyappend_num))

                                    hists_sys_numerator.append( self.histefficiencies["_".join((key_nominal,keyappend_num))] )
                                    hists_sys_denominator.append( self.histefficiencies["_".join((key_nominal,keyappend_denom))] )
                                    histlist.extend([ self.histefficiencies["_".join((key_nominal,keyappend_num))], self.histefficiencies["_".join((key_nominal,keyappend_denom))] ])
                            else:
                                    keyappend_numden = "_".join(("ND",sys,sysdir))

                                    if self.verbose:
                                        print "\tstoring efficiency: ", "_".join((key_nominal,keyappend_numden))

                                    hists_sys_numden.append( self.histefficiencies["_".join((key_nominal,keyappend_numden))] )
                                    histlist.extend([ self.histefficiencies["_".join((key_nominal,keyappend_numden))], self.histefficiencies["_".join((key_nominal,keyappend_numden))] ])


          	    legend = TLegend(0.45,0.5,0.925,0.8) # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
                    legend.SetHeader(eff + " - " + self.leptons_full[lep])
          	    legend.SetBorderSize(0)	# no border
          	    legend.SetFillStyle(0)	# Legend transparent background
          	    legend.SetTextSize(0.035)	# Increase entry font size!
          	    #legend.SetTextFont(42)	# Helvetica

                    leg_ATLAS  = TLatex()
                    leg_lumi   = TLatex()
                    leg_ATLAS.SetTextSize(0.04)
                    leg_ATLAS.SetNDC()
                    leg_lumi.SetTextSize(0.04)
                    leg_lumi.SetNDC()

	  	    legend.AddEntry(hist_nominal,"#varepsilon_{{{0}}} - nominal (stat. unc.)".format(eff), "P")

	  	    # Store the sum in quadrature of all the UP and DN variations wrt nominal to eventually make a single sys error band
		    # Assuming that all uncertainties in each bin are not correlated (which is the case!)

		    all_sys = {} # for each bin nr. store the sum in quadrature of sys for that bin

		    count_sys_num = [0 for i in range(len(self.__systematics[1:]))]

	  	    for h in hists_sys_numerator:

			h.SetLineStyle(1)
			h.SetLineWidth(2)
	  		h.SetMarkerStyle(kCircle)

			for idx, sys in enumerate(self.__systematics[1:]):
			    if sys in h.GetName():

				# print "\tsys hist name: ", h.GetName()
                                tokens = h.GetName().split("_")
				bin = tokens[-1] # get the bin number
				sys_var = abs( hist_nominal.GetBinContent(int(bin)) - h.GetBinContent(int(bin)) )
				# print("\t\tvariation[{0}] = {1}".format(bin,sys_var))

				if not all_sys.get(bin):
				    all_sys[bin] = pow(sys_var,2.0)
				else:
				    all_sys[bin] += pow(sys_var,2.0)
				# print("\t\tsum var[{0}] = {1}".format(bin,all_sys[bin]))

				h.SetLineColor(self.__syst_color_dict[sys+"_N"])
			        h.SetMarkerColor(self.__syst_color_dict[sys+"_N"])
			        if not count_sys_num[idx]:
			            legend.AddEntry(h,"{0} sub. sys. (TT)".format(sys), "P")
				    count_sys_num[idx] += 1

		    count_sys_denom = [0 for i in range(len(self.__systematics[1:]))]

	  	    for h in hists_sys_denominator:

			h.SetLineStyle(1)
			h.SetLineWidth(2)
	  		h.SetMarkerStyle(kOpenTriangleUp)

			for idx, sys in enumerate(self.__systematics[1:]):
			    if sys in h.GetName():

				# print "\tsys hist name: ", h.GetName()
                                tokens = h.GetName().split("_")
				bin = tokens[-1] # get the bin number
				sys_var = abs( hist_nominal.GetBinContent(int(bin)) - h.GetBinContent(int(bin)) )
				# print("\t\tvariation[{0}] = {1}".format(bin,sys_var))

				if not all_sys.get(bin):
				    all_sys[bin] = pow(sys_var,2.0)
				else:
				    all_sys[bin] += pow(sys_var,2.0)
				# print("\t\tsum var[{0}] = {1}".format(bin,all_sys[bin]))

			        h.SetLineColor(self.__syst_color_dict[sys+"_D"])
			        h.SetMarkerColor(self.__syst_color_dict[sys+"_D"])
			        if not count_sys_denom[idx]:
			            legend.AddEntry(h,"{0} sub. sys. (T#slash{{T}})".format(sys), "P")
				    count_sys_denom[idx] += 1

		    count_sys_numden = [0 for i in range(len(self.__systematics[1:]))]

	  	    for h in hists_sys_numden:

			h.SetLineStyle(1)
			h.SetLineWidth(2)
	  		h.SetMarkerStyle(kOpenTriangleUp)

			for idx, sys in enumerate(self.__systematics[1:]):
			    if sys in h.GetName():

				#print "\tsys hist name: ", h.GetName()
                                tokens = h.GetName().split("_")

                                for bin in range(1,hist_nominal.GetNbinsX()+2):
                                    sys_var = h.GetBinContent(int(bin)) - hist_nominal.GetBinContent(int(bin))
                                    bin = str(bin)
                                    #print("\t\tvariation[{0}] = {1}".format(bin,sys_var))
                                    if not all_sys.get(bin):
                                        all_sys[bin] = pow(sys_var,2.0)
                                    else:
                                        all_sys[bin] += pow(sys_var,2.0)
				    #print("\t\tsum var[{0}] = {1}".format(bin,all_sys[bin]))

                                    h.SetLineColor(self.__syst_color_dict[sys+"_ND"])
                                    h.SetMarkerColor(self.__syst_color_dict[sys+"_ND"])
                                    if not count_sys_numden[idx]:
                                        legend.AddEntry(h,"{0} sub. sys.".format(sys), "P")
                                        count_sys_numden[idx] += 1

	            # Now get sqrt of sum in quadrature of systematics

		    for key, sys_var in sorted(all_sys.iteritems()):
		        all_sys[key] = math.sqrt(all_sys[key])
			#print("\tallsys bin[{0}] = {1}".format(key,all_sys[key]))

	  	    c = TCanvas("c_"+ var + "_" + lep + "_" + eff,"Efficiencies")
          	    c.SetFrameFillColor(0)
          	    c.SetFrameFillStyle(0)
          	    c.SetFrameBorderMode(0)

          	    pad1 = TPad("pad1", "", 0, 0.25, 1, 1)
          	    pad2 = TPad("pad2", "", 0, 0,   1, 0.25)
          	    pad1.SetBottomMargin(0.02)
          	    pad2.SetBottomMargin(0.4)
          	    pad2.SetGridy(1)
          	    pad1.Draw()
          	    pad2.Draw()

		    # Be smarter when setting y range

		    scale_dn = scale_up = 1.0
		    if eff == "Fake":
		        scale_dn = 0.80
		        scale_up = 1.20
		    elif eff == "Real":
		        scale_dn = 0.95
		        scale_up = 1.05
		    ymin, ymax = self.__getLimits__(histlist, scale_up, scale_dn)
		    #hist_nominal.GetYaxis().SetRangeUser(ymin,ymax)
                    hist_nominal.SetMaximum(ymax*1.2)
                    hist_nominal.SetMinimum(0)

                    # Make a clone of nominal, and divide by itself

		    rationom = hist_nominal.Clone(hist_nominal.GetName()+"_Ratio")
                    rationom.Add(rationom,-1)
		    rationom.Divide(hist_nominal)
		    for bin in 	range(1,rationom.GetSize()):
		        error = 0.0
			if hist_nominal.GetBinContent(bin) > 0:
		            error = hist_nominal.GetBinError(bin) / hist_nominal.GetBinContent(bin)
		        rationom.SetBinError( bin, error )

		    rationom.SetLineStyle(1)
		    rationom.SetLineWidth(2)
		    rationom.SetMarkerSize(0)
             	    rationom.SetYTitle("#Delta#varepsilon/#varepsilon")
             	    rationom.GetXaxis().SetTitleSize(0.15)
             	    rationom.GetYaxis().SetTitleSize(0.15)
             	    rationom.GetXaxis().SetTitleOffset(1.0)
             	    rationom.GetYaxis().SetTitleOffset(0.35)
             	    rationom.GetXaxis().SetLabelSize(0.15)
             	    rationom.GetYaxis().SetLabelSize(0.12)
             	    rationom.GetYaxis().SetNdivisions(505)#(503)#(5)
	   	    rationom.SetFillColor(kGray+3)
           	    rationom.SetLineColor(kGray+3)
                    rationom.SetFillStyle(3356)
                    gStyle.SetHatchesLineWidth(1)
                    gStyle.SetHatchesSpacing(0.6)
           	    #rationom.SetFillStyle(3004)

		    ratiolist = []
	  	    for h in hists_sys_numerator:
			if self.verbose:
			    print "hist sys: ", h.GetName()
			    print("\t num   = [" + ",".join( "{0:.3f}".format(x) for x in [ h.GetBinContent(i) for i in range(1,h.GetSize()) ] ) + "]" )
			    print("\t denom = [" + ",".join( "{0:.3f}".format(x) for x in [ hist_nominal.GetBinContent(i) for i in range(1,hist_nominal.GetSize()) ] ) + "]" )
			ratio = h.Clone(h.GetName())
                        ratio.Add(hist_nominal,-1)
                        ratio.Divide(hist_nominal)
			if self.verbose:
			    print("\t ratio = [" + ",".join( "{0:.3f}".format(x) for x in [ ratio.GetBinContent(i) for i in range(1,ratio.GetSize()) ] ) + "]" )
			ratiolist.append(ratio)

	  	    for h in hists_sys_denominator:
			if self.verbose:
			    print "hist sys: ", h.GetName()
			    print("\t num   = [" + ",".join( "{0:.3f}".format(x) for x in [ h.GetBinContent(i) for i in range(1,h.GetSize()) ] ) + "]" )
			    print("\t denom = [" + ",".join( "{0:.3f}".format(x) for x in [ hist_nominal.GetBinContent(i) for i in range(1,hist_nominal.GetSize()) ] ) + "]" )
			ratio = h.Clone(h.GetName())
                        ratio.Add(hist_nominal,-1)
                        ratio.Divide(hist_nominal)
			if self.verbose:
			    print("\t ratio = [" + ",".join( "{0:.3f}".format(x) for x in [ ratio.GetBinContent(i) for i in range(1,ratio.GetSize()) ] ) + "]" )
			ratiolist.append(ratio)

	  	    for h in hists_sys_numden:
			if self.verbose:
			    print "hist sys: ", h.GetName()
			    print("\t num   = [" + ",".join( "{0:.3f}".format(x) for x in [ h.GetBinContent(i) for i in range(1,h.GetSize()) ] ) + "]" )
			    print("\t denom = [" + ",".join( "{0:.3f}".format(x) for x in [ hist_nominal.GetBinContent(i) for i in range(1,hist_nominal.GetSize()) ] ) + "]" )
			ratio = h.Clone(h.GetName())
                        ratio.Add(hist_nominal,-1)
                        ratio.Divide(hist_nominal)
			if self.verbose:
			    print("\t ratio = [" + ",".join( "{0:.3f}".format(x) for x in [ ratio.GetBinContent(i) for i in range(1,ratio.GetSize()) ] ) + "]" )
			ratiolist.append(ratio)

		    ratio_ymin, ratio_ymax = self.__getLimits__(ratiolist, ratio=True)
		    #rationom.GetYaxis().SetRangeUser(ratio_ymin, ratio_ymax)
		    ###rationom.SetMaximum(ratio_ymax*1.4)
		    #rationom.SetMinimum(ratio_ymin*(1.0/1.2))
		    ##rationom.SetMinimum(ratio_ymin*1.4)

                    rationom.SetMaximum(0.6)
                    rationom.SetMinimum(-0.6)

		    pad1.cd()
		    # Remove X axis labels from top pad
		    hist_nominal.GetXaxis().SetLabelSize(0)
                    hist_nominal.GetXaxis().SetLabelOffset(999)
	  	    hist_nominal.Draw("E0")
	  	    for h in histlist[1:]:
	  	        h.Draw("HIST, SAME")
	            legend.Draw()
                    leg_ATLAS.DrawLatex(0.6,0.35,"#bf{#it{ATLAS}} Work In Progress");
                    leg_lumi.DrawLatex(0.6,0.27,"#sqrt{{s}} = 13 TeV, #int L dt = {0:.1f} fb^{{-1}}".format(self.lumi));

		    pad2.cd()
		    rationom.Draw("E2")
		    for r in ratiolist:
			r.Draw("HIST SAME")

                    refl = TLine(rationom.GetBinLowEdge(1), 0.0, rationom.GetBinLowEdge(rationom.GetNbinsX()+1), 0.0)
                    refl.SetLineStyle(2)
		    refl.SetLineColor(kRed)
                    refl.SetLineWidth(2)
		    refl.Draw("SAME")

		    canvas_filename = "_".join((eff,lep,var,"Efficiency",proc,"Systematics"))

		    for extension in self.extensionlist:
		        c.SaveAs(savepath+"/SplitSys/"+canvas_filename+"."+extension)

                    # Reset axis labels to default (otherwise next plots won't have labels)

		    hist_nominal.GetXaxis().SetLabelSize()
                    hist_nominal.GetXaxis().SetLabelOffset()

	            # ---------------------------------------------------------------------------------------------------------------

	            # Make a set of plots for nominal + *symmetrised* combination of systematic uncertainties

	  	    c_allsys = TCanvas("c_"+ var + "_" + lep + "_" + eff + "_allsys","Efficiencies_CombinedSystematics")
          	    c_allsys.SetFrameFillColor(0)
          	    c_allsys.SetFrameFillStyle(0)
          	    c_allsys.SetFrameBorderMode(0)

          	    legend_allsys = TLegend(0.45,0.5,0.925,0.8) # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
                    legend_allsys.SetHeader(eff + " - " + self.leptons_full[lep])
          	    legend_allsys.SetBorderSize(0)	# no border
          	    legend_allsys.SetFillStyle(0)	# Legend transparent background
          	    legend_allsys.SetTextSize(0.035)	# Increase entry font size!
          	    #legend_allsys.SetTextFont(42)	# Helvetica

		    hist_allsys = hist_nominal.Clone(hist_nominal.GetName()+"_AllSys")
	   	    hist_allsys.SetFillColor(kOrange-2)
           	    hist_allsys.SetLineColor(kOrange-2)
           	    hist_allsys.SetFillStyle(1001)

	  	    legend_allsys.AddEntry(hist_nominal,"#varepsilon_{{{0}}} - nominal (stat. unc.)".format(eff), "P")
		    legend_allsys.AddEntry(hist_allsys,"Combined sys.", "F")

		    for bin in 	range(1,hist_allsys.GetSize()):
			# A dirty hack: if the error is zero, drawing w/ option E2 seems
			# to be ingnored and HIST gets used instead.
			# Thus, just set a tiny error on the hist
			#
                        # print "bin: ", bin
			error = 0.0001
			if all_sys[str(bin)]:
			    error = all_sys[str(bin)]
		        hist_allsys.SetBinError(bin,error)

		    ymin_allsys, ymax_allsys = self.__getLimits__([hist_allsys,hist_nominal], scale_up, scale_dn)
		    #hist_allsys.GetYaxis().SetRangeUser(ymin_allsys,ymax_allsys)
                    hist_allsys.SetMaximum(ymax_allsys*1.2)
                    hist_allsys.SetMinimum(0)

		    ratio_allsys = hist_nominal.Clone(hist_nominal.GetName()+"_Ratio_AllSys")
                    ratio_allsys.Add(hist_nominal,-1)
		    ratio_allsys.Divide(hist_nominal)

                    #print("Setting error on ratio_allsys:")
		    for bin in 	range(1,ratio_allsys.GetSize()):
			# A dirty hack: if the error is zero, drawing w/ option E2 seems
			# to be ingnored and HIST gets used instead.
			# Thus, just set a tiny error on the hist
			#
			error = 0.0001
			if hist_nominal.GetBinContent(bin) > 0 and all_sys[str(bin)]:
			    #error = 2.0 * ( all_sys[str(bin)] / hist_nominal.GetBinContent(bin) )
                            #print("all_sys[{0}] = {1:.3f}, nominal[{0}] = {2:.3f}".format(bin,all_sys[str(bin)],hist_nominal.GetBinContent(bin)))
			    error = all_sys[str(bin)] / hist_nominal.GetBinContent(bin)
                            #print("error = ( all_sys[{0}] ) / nominal[{0}] = {1:.3f}".format(bin,error))
		        ratio_allsys.SetBinError( bin, error )

		    ratio_allsys.SetLineStyle(1)
		    ratio_allsys.SetLineColor(kBlack)
		    ratio_allsys.SetLineWidth(2)
		    ratio_allsys.SetMarkerSize(0)
             	    ratio_allsys.SetYTitle("#Delta#varepsilon/#varepsilon")
             	    ratio_allsys.GetXaxis().SetTitleSize(0.15)
             	    ratio_allsys.GetYaxis().SetTitleSize(0.15)
             	    ratio_allsys.GetXaxis().SetTitleOffset(1.0)
             	    ratio_allsys.GetYaxis().SetTitleOffset(0.35)
             	    ratio_allsys.GetXaxis().SetLabelSize(0.15)
             	    ratio_allsys.GetYaxis().SetLabelSize(0.12)
             	    ratio_allsys.GetYaxis().SetNdivisions(5)#(503)#(5)
	   	    ratio_allsys.SetFillColor(kOrange-2)
           	    ratio_allsys.SetFillStyle(1001)

                    # Need to re-get the nominal ratio hist with stats errors

		    ratio_ymin_allsys, ratio_ymax_allsys = self.__getLimits__([ratio_allsys,rationom], shift_up=0.2, shift_dn=0.2, ratio=True)

		    print ("ratio_ymin_allsys = {0:.1f}, ratio_ymax_allsys = {1:.1f}".format(ratio_ymin_allsys, ratio_ymax_allsys))
		    #ratio_allsys.GetYaxis().SetRangeUser(round(ratio_ymin_allsys,1), round(ratio_ymax_allsys,1))
                    ratio_allsys.SetMaximum(ratio_ymax_allsys*1.2)
                    ratio_allsys.SetMinimum(ratio_ymin_allsys*(1.0/1.2))

          	    pad1_allsys = TPad("pad1", "", 0, 0.25, 1, 1)
          	    pad2_allsys = TPad("pad2", "", 0, 0,   1, 0.25)
          	    pad1_allsys.SetBottomMargin(0.02)
          	    pad2_allsys.SetBottomMargin(0.4)
          	    pad2_allsys.SetGridy(1)
          	    pad1_allsys.Draw()
          	    pad2_allsys.Draw()

		    pad1_allsys.cd()
		    # Remove X axis labels from top pad
		    hist_allsys.GetXaxis().SetLabelSize(0)
                    hist_allsys.GetXaxis().SetLabelOffset(999)
		    hist_allsys.Draw("E2")
		    hist_nominal.Draw("E0 SAME")

                    legend_allsys.Draw()
		    leg_ATLAS.DrawLatex(0.6,0.35,"#bf{#it{ATLAS}} Work In Progress")
                    leg_lumi.DrawLatex(0.6,0.27,"#sqrt{{s}} = 13 TeV, #int L dt = {0:.1f} fb^{{-1}}".format(self.lumi))

		    print("NOMINAL: bincontent    = [" + ",".join( "{0:.3f}".format(x) for x in [ hist_nominal.GetBinContent(ibin) for ibin in range(1,hist_nominal.GetSize()) ] ) + "]" )
		    print("NOMINAL: binerror (+-) = [" + ",".join( "{0:.3f}".format(x) for x in [ hist_nominal.GetBinError(ibin) for ibin in range(1,hist_nominal.GetSize()) ] ) + "]" )
		    print("ALLSYS:  bincontent    = [" + ",".join( "{0:.3f}".format(x) for x in [ hist_allsys.GetBinContent(ibin) for ibin in range(1,hist_allsys.GetSize()) ] ) + "]" )
		    print("ALLSYS:  binerror (+-) = [" + ",".join( "{0:.3f}".format(x) for x in [ hist_allsys.GetBinError(ibin) for ibin in range(1,hist_allsys.GetSize()) ] ) + "]" )

		    pad2_allsys.cd()
		    rationom.Draw("E2")
		    ratio_allsys.Draw("E2 SAME")
		    refl.Draw("SAME")

		    print("RATIO NOMINAL: bincontent    = [" + ",".join( "{0:.3f}".format(x) for x in [ rationom.GetBinContent(ibin) for ibin in range(1,rationom.GetSize()) ] ) + "]" )
		    print("RATIO NOMINAL: binerror (+-) = [" + ",".join( "{0:.3f}".format(x) for x in [ rationom.GetBinError(ibin) for ibin in range(1,rationom.GetSize()) ] ) + "]" )
		    print("RATIO ALLSYS:  bincontent    = [" + ",".join( "{0:.3f}".format(x) for x in [ ratio_allsys.GetBinContent(ibin) for ibin in range(1,ratio_allsys.GetSize()) ] ) + "]" )
		    print("RATIO ALLSYS:  binerror (+-) = [" + ",".join( "{0:.3f}".format(x) for x in [ ratio_allsys.GetBinError(ibin) for ibin in range(1,ratio_allsys.GetSize()) ] ) + "]" )

		    canvas_allsys_filename = "_".join((eff,lep,var,"Efficiency",proc,"CombinedSystematics"))

		    for extension in self.extensionlist:
		        c_allsys.SaveAs(savepath+"/CombinedSys/"+canvas_allsys_filename+"."+extension)


    def __factorToEfficiency__(self, f):
        if f < 0:
            f = 0.0
	e = f/(f+1)
        return e

    def __efficiencyToFactor__(self, e):
	f = e/(1-e)
        return f

    def checkRebin(self, debug_msg=None):

    	if self.debug:

	    if debug_msg:
	        print("\n{0}".format(debug_msg))

    	    print("\n\tNUMERATOR histograms dictionary:\n")
    	    print("\tkey\t\thistname\t\tnbins\n")
    	    for key, value in sorted( self.numerator_hists.iteritems() ):
    		print("\t{0}\t{1}\t{2}".format(key, value.GetName(), value.GetSize()-2))

    	    print("\n\tDENOMINATOR histograms dictionary:\n")
    	    print("\tkey\t\thistname\t\tnbins\n")
    	    for key, value in sorted( self.denominator_hists.iteritems() ):
    		print("\t{0}\t{1}\t{2}".format(key, value.GetName(), value.GetSize()-2))

    	    print("\n\tANTI-NUMERATOR histograms dictionary:\n")
    	    print("\tkey\t\thistname\t\tnbins\n")
    	    for key, value in sorted( self.antinumerator_hists.iteritems() ):
    		print("\t{0}\t{1}\t{2}".format(key, value.GetName(), value.GetSize()-2))


    def checkYields(self, debug_msg=None):

	if self.debug:

	    if debug_msg:
	        print("\n{0}".format(debug_msg))

    	    print("\n\tNUMERATOR yields dictionary:\n")
    	    print("\tkey\t\tyields (per bin)\t\tintegral\n")
    	    for key, value in sorted( self.numerator_yields.iteritems() ):
    		print("\t{0}".format(key) + "\t[" + ",".join( "{0:.3f}".format(i) for i in value[:-1] ) + "]" + "\t{0:.3f}".format(value[-1]) )

    	    print("\n\tDENOMINATOR yields dictionary:\n")
    	    print("\tkey\t\tyields (per bin)\t\tintegral\n")
    	    for key, value in sorted( self.denominator_yields.iteritems() ):
    		print("\t{0}".format(key) + "\t[" + ",".join( "{0:.3f}".format(i) for i in value[:-1] ) + "]" + "\t{0:.3f}".format(value[-1]) )

    	    print("\n\tANTI-NUMERATOR yields dictionary:\n")
    	    print("\tkey\t\tyields (per bin)\t\tintegral\n")
    	    for key, value in sorted( self.antinumerator_yields.iteritems() ):
    		print("\t{0}".format(key) + "\t[" + ",".join( "{0:.3f}".format(i) for i in value[:-1] ) + "]" + "\t{0:.3f}".format(value[-1]) )


# --------------------------------------------------------------------------------------------------------------

if __name__ == "__main__":

    print args

    for sys in args.systematics:
        if not sys in g_avaialble_systematics:
	    print("\nWARNING!\nSystematic {0} is not supported!!".format(sys))

    eff = RealFakeEffTagAndProbe( closure=args.closure, factors=args.factors, variables=args.variables, efficiencies=args.efficiency, leptons=args.lepton, systematics=args.systematics, nosub=args.nosub )

    if args.triggerEff:
        print("Measuring trigger efficiency for selection {0} ...\n".format(args.triggerEff))
	eff.triggerEff = args.triggerEff
    	eff.selections = {"D":eff.triggerEff+"AnyTM","N":eff.triggerEff+"TM","AntiN":eff.triggerEff+"AntiTM"}

    if args.probeAssignEff:
        print("Measuring probe lepton assignment efficiency...\n")
	eff.probeAssignEff = args.probeAssignEff
        if "Real" in eff.efficiencies: eff.efficiencies.remove("Real")
        eff.selections = {"D":"ProbeTTMMatchedToAny","N":"ProbeTTMMatchedToFake","AntiN":"ProbeTTMMatchedToReal"}
    	#eff.selections = {"D":"ProbeMatchedToAny","N":"ProbeMatchedToFake","AntiN":"ProbeMatchedToReal"}

    eff.lumi  = args.lumi
    eff.debug = args.debug
    eff.verbose = args.verbose
    eff.log   = args.log

    #eff.addProcess(processlist=["expectedbkg"])

    eff.readInputs( inputpath=args.inputpath, channel=args.channel )

    eff.rebinHistograms( rebinlist=args.rebin, averagehistlist=args.averagehist )
    eff.checkRebin("Checking rebinning...")

    eff.checkYields("Checking event yields BEFORE subtraction")
    eff.subtractHistograms()
    if not args.nosub:
        eff.checkYields("Checking event yields AFTER subtraction")

    eff.computeEfficiencies(variation="nominal")
    eff.computeEfficiencies(variation="N")
    eff.computeEfficiencies(variation="D")
    eff.computeEfficiencies(variation="ND")

    if eff.factors:
        eff.computeFactors("nominal")
        eff.computeFactors("N")
        eff.computeFactors("D")
        eff.computeFactors("ND")

    eff.saveOutputs( filename=args.outfilename, outputpath=args.outpath )

    if args.plots:
        eff.plotMaker()
        if args.systematics and not args.closure and not args.nosub:
	    eff.plotMakerSys()
