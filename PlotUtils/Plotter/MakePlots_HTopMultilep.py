#!/usr/bin/env python

""" MakePlots_HTopMultilep.py: plotting script for the HTopMultilep RunII analysis """

__author__     = "Marco Milesi, Francesco Nuti"
__email__      = "marco.milesi@cern.ch, francesco.nuti@cern.ch"
__maintainer__ = "Marco Milesi"

import os
import sys
import math

sys.path.append(os.path.abspath(os.path.curdir))

# -------------------------------
# Parser for command line options
# -------------------------------
import argparse

parser = argparse.ArgumentParser(description='Plotting script for the HTopMultilep RunII analysis')

#***********************************
# positional arguments (compulsory!)
#***********************************
parser.add_argument('inputDir', metavar='inputDir',type=str,
                   help='Path to the directory containing input files')
parser.add_argument('samplesCSV', metavar='samplesCSV',type=str,
                   help='Path to the csv file containing the processes of interest with their cross sections and other metadata')
#*******************
# optional arguments
#*******************

list_available_channel    = ['TwoLepSR','ThreeLepSR','FourLepSR','MMRates','MMRates_LHFit',
			     'TwoLepLowNJetCR', 'ThreeLepLowNJetCR',
			     'WZonCR', 'WZoffCR', 'WZHFonCR', 'WZHFoffCR',
			     'ttWCR', 'ttZCR','ZSSpeakCR', 'DataMC', 'MMClosureTest',
			     'MMClosureRates','CutFlowChallenge']
list_available_fakemethod = ['MC','MM','FF','ABCD']
list_available_flavcomp   = ['OF','SF','INCLUSIVE']

parser.add_argument('--debug', dest='debug', action='store_true', default=False,
                    help='Run in debug mode')
parser.add_argument('--channel', dest='channel', action='store', default='TwoLepSR', type=str,
		    help='The channel chosen ({0})'.format(list_available_channel))
parser.add_argument('--ratesFromMC', dest='ratesFromMC', action='store_true', default=False,
                    help='Extract rates from pure simulation. Use w/ option --channel=MMRates')
parser.add_argument('--useMCChFlip', dest='useMCChFlip', action='store_true', default=False,
                    help='Use Monte-Carlo based estimate of QMisID')
parser.add_argument('--MCCompRF', dest='MCCompRF', action='store_true', default=False,
                    help='Plot simulation estimate in real/fake efficiency CRs. Use w/ option --channel=MMRates')
parser.add_argument('--outdirname', dest='outdirname', action='store', default='', type=str,
		    help='Specify a name to append to the output directory')
parser.add_argument('--fakeMethod', dest='fakeMethod', action='store', default='MC', type=str,
		    help='The fake estimation method chosen ({0})'.format(list_available_fakemethod))
parser.add_argument('--lepFlavComp', dest='lepFlavComp', action='store', default='INCLUSIVE', type=str,
                    help='Flavour composition of the dilepton pair used for efficiency measurement. Use w/ option --channel=MMRates,MMClosureRates. Default is \'INCLUSIVE\' ({0})'.format(list_available_flavcomp))
parser.add_argument('--doLogScaleX', dest='doLogScaleX', action='store_true',
                    help='Use log scale on the X axis')
parser.add_argument('--doLogScaleY', dest='doLogScaleY', action='store_true',
                    help='Use log scale on the Y axis')
parser.add_argument('--doSyst', dest='doSyst', action='store_true',
                    help='Run systematics')
parser.add_argument('--noSignal', action='store_true', dest='noSignal',
                    help='Exclude signal')
parser.add_argument('--noWeights', action='store_true', dest='noWeights', default=False,
                    help='Do not apply weights, except for mcEventWeight and Xsec*lumi...')
parser.add_argument('--noStandardPlots', action='store_true', dest='noStandardPlots',
                    help='exclude all standard plots')
parser.add_argument('--doEPS', action='store_true', dest='doEPS',
                    help='Make .eps output files')
parser.add_argument('--doChFlipRate', dest='doChFlipRate', action='store_true',
                    help='Measure charge flip rate in MC (to be used with --channel=MMClosureRates)')
parser.add_argument('--doUnblinding', dest='doUnblinding', action='store_true', default=False,
                    help='Unblind data in SRs')
parser.add_argument('--printEventYields', dest='printEventYields', action='store_true', default=False,
                    help='Prints out event yields in tabular form (NB: can be slow)')

args = parser.parse_args()

# -------------------------------
# Important to run without popups
# -------------------------------
from ROOT import gROOT

gROOT.SetBatch(True)

# -----------------
# Some ROOT imports
# -----------------
from ROOT import TH1I, TMath, TFile, TAttFill, TColor, kBlack, kWhite, kGray, kBlue, kRed, kYellow, kAzure, kTeal, kSpring, kOrange, kGreen, kCyan, kViolet, kMagenta, kPink, Double

# ---------------------------------------------------------------------
# Importing all the tools and the definitions used to produce the plots
# ---------------------------------------------------------------------
from Plotter.BackgroundTools_ttH import loadSamples, Category, Background, Process, VariableDB, Variable, Cut, Systematics, Category
# ---------------------------------------------------------------------------
# Importing the classes for the different processes.
# They contains many info on the normalization and treatment of the processes
# ---------------------------------------------------------------------------
from Plotter.ttH2015_Background import MyCategory, TTHBackgrounds2015

try:
    args.channel in list_available_channel
except ValueError:
    sys.exit('the channel specified (', args.channel ,') is incorrect. Must be one of: ', list_available_channel)

doTwoLepSR      	= bool( args.channel == 'TwoLepSR' )
doThreeLepSR    	= bool( args.channel == 'ThreeLepSR' )
doFourLepSR     	= bool( args.channel == 'FourLepSR' )
doMMRates          	= bool( args.channel == 'MMRates' )
doMMRatesLHFit          = bool( args.channel == 'MMRates_LHFit')
doTwoLepLowNJetCR       = bool( args.channel == 'TwoLepLowNJetCR' )
doThreeLepLowNJetCR     = bool( args.channel == 'ThreeLepLowNJetCR' )
doWZonCR                = bool( args.channel == 'WZonCR' )
doWZoffCR               = bool( args.channel == 'WZoffCR' )
doWZHFonCR              = bool( args.channel == 'WZHFonCR' )
doWZHFoffCR             = bool( args.channel == 'WZHFoffCR' )
dottWCR                 = bool( args.channel == 'ttWCR' )
dottZCR                 = bool( args.channel == 'ttZCR' )
doZSSpeakCR             = bool( args.channel == 'ZSSpeakCR' )
doDataMCCR              = bool( args.channel == 'DataMC' )
doMMClosureTest         = bool( args.channel == 'MMClosureTest' )
doMMClosureRates        = bool( args.channel == 'MMClosureRates' )
doCFChallenge           = bool( args.channel == 'CutFlowChallenge' )

# -----------------------------------------
# a comprehensive flag for all possible SRs
# -----------------------------------------
doSR = (doTwoLepSR or doThreeLepSR or doFourLepSR)

# -----------------------------------------
# a comprehensive flag for the low-Njet CR
# -----------------------------------------
doLowNJetCR = (doTwoLepLowNJetCR or doThreeLepLowNJetCR)

# ------------------------------------------
# a comprehensive flag for all the other CRs
# ------------------------------------------
doOtherCR = (doWZonCR or doWZoffCR or doWZHFonCR or doWZHFoffCR or dottWCR or dottZCR or doZSSpeakCR or doMMRates or doMMRatesLHFit or doDataMCCR or doMMClosureTest or doMMClosureRates or doCFChallenge)

# ------------------------------------------------
# make standard plots unless differently specified
# ------------------------------------------------
doStandardPlots = False if (args.noStandardPlots) else (doSR or doLowNJetCR or doOtherCR)

# ----------------------------
# Check fake estimation method
# ----------------------------

try:
    args.fakeMethod in list_available_fakemethod
except ValueError:
    sys.exit('the Fake Method specified (', args.fakeMethod ,') is incorrect. Must be one of ', list_available_fakemethod)

doMM   = bool( args.fakeMethod == 'MM' )
doFF   = bool( args.fakeMethod == 'FF' )
doABCD = bool( args.fakeMethod == 'ABCD' )

# -----------------------------------------------------
# Check lepton flavour composition of the dilepton pair
# for efficiency measurement
# -----------------------------------------------------

try:
    args.lepFlavComp in list_available_flavcomp
except ValueError:
    sys.exit('ERROR: the lepton flavour composition of the 2Lep pair is incorrect. Must be one of ', list_available_flavcomp)

# ----------------------------------------------------
# When in debug mode, print out all the input commands
# ----------------------------------------------------
if ( args.debug ):
    print args

# --------------------------
# Retrieve the input samples
# --------------------------
inputs = loadSamples(
                        # path of the data to be processed
                        inputdir    = args.inputDir,
                        samplescsv  = args.samplesCSV,
                        nomtree     = 'physics',
                        # name of the trees that contains values for shifted systematics
                        systrees =  [
                                        ##'METSys',
                                        #'ElEnResSys',
		                        #'ElES_LowPt',
					#'ElES_Zee',
					#'ElES_R12',
					#'ElES_PS',
                                        ##'EESSys',
					#'MuSys',
                                        ##'METResSys',
                                        ##'METScaleSys',
                                        #'JES_Total',
                                        #'JER',
                                    ],
                    )

# ------------------------------------------------------
# Here you include all names of variables to be plotted,
# with min, max, number of bins and ntuple name.
# ------------------------------------------------------
vardb = VariableDB()

doRelaxedBJetCut = False # if true, be inclusive in bjet multiplicity, otherwise, will require at least 1 BTagged jet

# -----------------------------------------------------
# The list of event-level TCuts
#
# WARNING:
# To avoid unexpected behaviour,
# ALWAYS enclose the cut string in '()'!!!
#
# -----------------------------------------------------

# ---------------------
# General cuts
# ---------------------

vardb.registerCut( Cut('DummyCut',    '( 1 )') )
vardb.registerCut( Cut('IsMC',        '( isMC == 1 )') )
# To ask for an event be passing an OR of triggers
#
#gROOT.LoadMacro("$ROOTCOREBIN/user_scripts/HTopMultilepAnalysis/ROOT_TTreeFormulas/containsAny.cxx+")
#from ROOT import containsAny

#vardb.registerCut( Cut('TrigDec',     '( passHLT == 1 && ( Sum$( ( isMC == 1 && passedTriggers == \"HLT_e24_lhmedium_L1EM18VH\" ) + ( isMC == 0 && passedTriggers == \"HLT_e24_lhmedium_L1EM20VH\" ) + ( passedTriggers == \"HLT_e60_lhmedium\" ) + ( passedTriggers == \"HLT_e120_lhloose\" ) + ( passedTriggers == \"HLT_mu20_iloose_L1MU15\" ) + ( passedTriggers == \"HLT_mu50\" ) ) > 0 ) )') )
# the following does not work b/c TTreeFormula does not accept strings as arguments
#vardb.registerCut( Cut('TrigDec',    '( passHLT == 1 && ( containsAny(passedTriggers,\"HLT_e24_lhmedium_L1EM18VH,HLT_e24_lhmedium_L1EM20VH,HLT_e60_lhmedium,HLT_e120_lhloose,HLT_mu20_iloose_L1MU15,HLT_mu50\") ) )') )
#
# To ask for an event be passing any of the saved triggers
#
vardb.registerCut( Cut('TrigDec',     '( passHLT == 1 )') )
vardb.registerCut( Cut('LargeNBJet',  '( njets_mv2c20_Fix77 > 1 )') )
vardb.registerCut( Cut('VetoLargeNBJet',  '( njets_mv2c20_Fix77 < 4 )') )
vardb.registerCut( Cut('BJetVeto',    '( njets_mv2c20_Fix77 == 0 )') )
vardb.registerCut( Cut('OneBJet',     '( njets_mv2c20_Fix77 == 1 )') )
vardb.registerCut( Cut('TauVeto',     '( ntau == 0 )') )
vardb.registerCut( Cut('OneTau',      '( ntau == 1 )') )

# ---------------------
# 3lep cuts
# ---------------------

vardb.registerCut( Cut('3Lep_NLep',         '( nlep == 3 && ( lep_3lepClosestSS_pt[0] > 20e3 && lep_3lepOtherSS_pt[0] > 20e3 && lep_3lepOS_pt[0] > 10e3 ) )') )
vardb.registerCut( Cut('3Lep_Charge',       '( TMath::Abs( Sum$(lep_charge) ) == 1 )') )
vardb.registerCut( Cut('3Lep_TightLeptons', '( lep_3lepClosestSS_isTightSelected[0] == 1 && lep_3lepOtherSS_isTightSelected[0] == 1 && TMath::Abs(lep_3lepOS_trkz0sintheta[0]) < 0.5 && ( ( lep_3lepOS_flavour[0] == 13 && TMath::Abs(lep_3lepOS_trkd0sig[0]) < 3.0 ) || ( lep_3lepOS_flavour[0] == 11 && TMath::Abs(lep_3lepOS_trkd0sig[0]) < 5.0 ) ) )') )
vardb.registerCut( Cut('3Lep_TrigMatch',   '( ( lep_isTrigMatched[0] == 1 && ( ( lep_flavour[0] == 11 && lep_pt[0] > 25e3 ) || ( lep_flavour[0] == 13 && lep_pt[0] > 21e3 ) ) ) || ( lep_isTrigMatched[1] == 1 && ( ( lep_flavour[1] == 11 && lep_pt[1] > 25e3 ) || ( lep_flavour[1] == 13 && lep_pt[1] > 21e3 ) ) ) || ( lep_isTrigMatched[1] == 1 && ( ( lep_flavour[2] == 11 && lep_pt[2] > 25e3 ) || ( lep_flavour[2] == 13 && lep_pt[2] > 21e3 ) ) ) )') )
vardb.registerCut( Cut('3Lep_ZVeto',        '( ( isOSPairSF01 == 0 || ( isOSPairSF01 == 1 &&  TMath::Abs( mOSPair01 - 91e3 ) < 10e3 ) ) || ( isOSPairSF02 == 0 || ( isOSPairSF02 == 1 &&  TMath::Abs( mOSPair02 - 91e3 ) < 10e3 ) ) )') )
vardb.registerCut( Cut('3Lep_MinZCut',      '( ( isOSPairSF01 == 0 || ( isOSPairSF01 == 1 &&  mOSPair01 > 12e3 ) ) || ( isOSPairSF02 == 0 || ( isOSPairSF02 == 1 &&  mOSPair02 > 12e3 ) ) )') )
vardb.registerCut( Cut('3Lep_NJets',        '( ( njets_mv2c20_Fix77 > 0 && njets > 3 ) || ( njets_mv2c20_Fix77 > 1 && njets == 3) )') )

# ---------------------
# 4lep cuts
# ---------------------

vardb.registerCut( Cut('4Lep_NJets',      '( njets >= 2 )') )
vardb.registerCut( Cut('4Lep',        '( nlep == 4 && lep_pt[0] > 25e3 && lep_pt[1] > 15e3 && lep_isTightSelected[0] == 1 && lep_isTightSelected[1] == 1 && lep_isTightSelected[2] == 1 && lep_isTightSelected[3] == 1 && Sum$( lep_charge ) == 0 && ( ( mJPsiCand_ee > 10e3 || mJPsiCand_ee < 0.0 ) && ( mJPsiCand_mm > 10e3 || mJPsiCand_mm < 0.0 ) ) )') )

# ---------------------
# 2Lep SS + 1 tau cuts
# ---------------------

vardb.registerCut( Cut('2Lep1Tau_NLep',         '( nlep == 2 )') )
vardb.registerCut( Cut('2Lep1Tau_TightLeptons', '( Sum$( lep_isTightSelected ) == 2 )') )
vardb.registerCut( Cut('2Lep1Tau_pT',           '( Min$( lep_pt ) > 15e3 )') )
vardb.registerCut( Cut('2Lep1Tau_TrigMatch',    '( ( lep_isTrigMatched[0] == 1 && lep_pt[0] > 25e3 ) || ( lep_isTrigMatched[1] == 1 && lep_pt[1] > 25e3 ) )') )
vardb.registerCut( Cut('2Lep1Tau_SS',           '( isSS01 == 1 )') )
vardb.registerCut( Cut('2Lep1Tau_1Tau',         '( ntau == 1 && ( lep_charge[0] * tau_charge[0] ) < 0 )') )
vardb.registerCut( Cut('2Lep1Tau_Zsidescut',    '( nel <= 1 || ( nel == 2 && TMath::Abs( mll01 - 91.187e3 ) > 10e3 ) )' )  )
vardb.registerCut( Cut('2Lep1Tau_NJet_SR',      '( njets >= 4 )') )
vardb.registerCut( Cut('2Lep1Tau_NJet_CR',      '( njets > 1 && njets < 4 )') )
vardb.registerCut( Cut('2Lep1Tau_NBJet',        '( njets_mv2c20_Fix77 > 0 )') )

# ---------------------
# 2Lep SS + 0 tau cuts
# ---------------------
vardb.registerCut( Cut('2Lep_TrigMatch',   '( ( lep_isTrigMatched[0] == 1 && ( ( lep_flavour[0] == 11 && lep_pt[0] > 25e3 ) || ( lep_flavour[0] == 13 && lep_pt[0] > 21e3 ) ) ) || ( lep_isTrigMatched[1] == 1 && ( ( lep_flavour[1] == 11 && lep_pt[1] > 25e3 ) || ( lep_flavour[1] == 13 && lep_pt[1] > 21e3 ) ) ) )') )
vardb.registerCut( Cut('2Lep_TrigMatchDataMC', '( ( lep_isTrigMatched[0] == 1 || lep_isTrigMatched[1] == 1 ) )') )
if doRelaxedBJetCut:
    print("\nUsing relaxed nr. bjet cut: INCLUSIVE bjet multiplicity...\n")
    vardb.registerCut( Cut('2Lep_NBJet',      '( njets_mv2c20_Fix77 >= 0 )') )
else:
  vardb.registerCut( Cut('2Lep_NBJet',       	      '( njets_mv2c20_Fix77 > 0 )') )
vardb.registerCut( Cut('2Lep_NBJet_SR',      	      '( njets_mv2c20_Fix77 > 0 )') )
vardb.registerCut( Cut('2Lep_NJet_SR',       	      '( njets >= 4 )') )
vardb.registerCut( Cut('2Lep_NJet_CR',       	      '( njets > 1 && njets < 4 )') ) # use this for Moriond
vardb.registerCut( Cut('2Lep_NJet_CR_ttW',   	      '( njets > 1 && njets < 4 )') )
vardb.registerCut( Cut('2Lep_NJet_CR_SStt',  	      '( njets < 4 )') )
vardb.registerCut( Cut('2Lep_SS',            	      '( isSS01 == 1 )') )
vardb.registerCut( Cut('2Lep_OS',            	      '( isSS01 != 1 )') )
vardb.registerCut( Cut('2Lep_NLep_MMRates',  	      '( nlep == 2 && Min$( lep_pt ) > 10e3 )') )
vardb.registerCut( Cut('2Lep_NLep',          	      '( nlep == 2 && Min$( lep_pt ) > 25e3 )') )
vardb.registerCut( Cut('2Lep_NLep_Relaxed',  	      '( nlep == 2 && lep_pt[0] > 25e3 && lep_pt[1] > 10e3 )') )
vardb.registerCut( Cut('2Lep_SF_Event',      	      '( nmuon == 2 || nel == 2 )') )
vardb.registerCut( Cut('2Lep_MuMu_Event',    	      '( nmuon == 2 && nel == 0 )') )
vardb.registerCut( Cut('2Lep_ElEl_Event',    	      '( nel == 2 && nmuon == 0 )') )
vardb.registerCut( Cut('2Lep_OF_Event',      	      '( nmuon == 1 && nel == 1 )') )
vardb.registerCut( Cut('2Lep_MuEl_Event',    	      '( nmuon == 1 && nel == 1 && lep_flavour[0] == 13 )') )
vardb.registerCut( Cut('2Lep_ElMu_Event',    	      '( nmuon == 1 && nel == 1 && lep_flavour[0] == 11 )') )
vardb.registerCut( Cut('2Lep_ElEtaCut',               '( nel == 0 || ( nel > 0 && Max$( TMath::Abs(el_caloCluster_eta) ) < 1.37 ) )') )
vardb.registerCut( Cut('2Lep_ElTagEtaCut',            '( ( lep_tag_flavour[0] == 13 ) || ( lep_tag_flavour[0] == 11 && TMath::Abs(lep_tag_eta[0]) < 1.37 ) )') )

vardb.registerCut( Cut('2Lep_Zsidescut',  	      '( ( nmuon == 1 && nel == 1 ) || ( TMath::Abs( mll01 - 91.187e3 ) > 7.5e3 ) )' ) )   # Use this to require the 2 SF leptons to be outside Z peak
vardb.registerCut( Cut('2Lep_Zpeakcut',   	      '( ( nmuon == 1 && nel == 1 ) || ( TMath::Abs( mll01 - 91.187e3 ) < 30e3  ) )' ) )   # Use this to require the 2 SF leptons to be around Z peak
vardb.registerCut( Cut('2Lep_Zmincut',    	      '( ( nmuon == 1 && nel == 1 ) || ( mll01  > 40e3 ) )' ) )   # We don't have Sherpa Z+jets w/ low invariant mass (yet...)

vardb.registerCut( Cut('2Lep_LepTagTightTrigMatched',   '( lep_tag_isTightSelected[0] == 1 && lep_tag_isTrigMatched[0] == 1 && ( ( lep_tag_flavour[0] == 11 && lep_tag_pt[0] > 25e3 ) || ( lep_tag_flavour[0] == 13 && lep_tag_pt[0] > 21e3 ) ) )') )
vardb.registerCut( Cut('2Lep_ProbeTight',	        '( lep_probe_isTightSelected[0] == 1 )') )
vardb.registerCut( Cut('2Lep_ProbeAntiTight',	        '( lep_probe_isTightSelected[0] == 0 )') )
vardb.registerCut( Cut('2Lep_ProbeEl',  	        '( lep_probe_flavour[0] == 11 )') )
vardb.registerCut( Cut('2Lep_ProbeMu',  	        '( lep_probe_flavour[0] == 13 )') )
vardb.registerCut( Cut('2Lep_ElRealFakeRateCR',   	'( isProbeElEvent == 1 )') )
vardb.registerCut( Cut('2Lep_ElProbeTight',	  	'( el_probe_isTightSelected[0] == 1 )') )
vardb.registerCut( Cut('2Lep_ElProbeAntiTight',   	'( el_probe_isTightSelected[0] == 0 )') )
vardb.registerCut( Cut('2Lep_MuRealFakeRateCR',   	'( isProbeMuEvent == 1 )') )
vardb.registerCut( Cut('2Lep_MuProbeTight',	  	'( muon_probe_isTightSelected[0] == 1 )') )
vardb.registerCut( Cut('2Lep_MuProbeAntiTight',   	'( muon_probe_isTightSelected[0] == 0 )') )

"""
vardb.registerCut( Cut('TT',      '( Sum$( lep_isTightSelected ) == 2 )') )
vardb.registerCut( Cut('TL',      '( lep_isTightSelected[0] == 1 && lep_isTightSelected[1] == 0 )') )
vardb.registerCut( Cut('LT',      '( lep_isTightSelected[0] == 0 && lep_isTightSelected[1] == 1 )') )
vardb.registerCut( Cut('TL_LT',   '( ( lep_isTightSelected[0] == 1 && lep_isTightSelected[1] == 0 ) || ( lep_isTightSelected[0] == 0 && lep_isTightSelected[1] == 1 ) )') )
vardb.registerCut( Cut('LL',      '( Sum$( lep_isTightSelected ) == 0 )') )
vardb.registerCut( Cut('TelLmu',  '( ( lep_flavour[0] == 11 && lep_isTightSelected[0] == 1 ) && ( lep_flavour[1] == 13 && lep_isTightSelected[1] == 0 ) )') )
vardb.registerCut( Cut('LelTmu',  '( ( lep_flavour[0] == 11 && lep_isTightSelected[0] == 0 ) && ( lep_flavour[1] == 13 && lep_isTightSelected[1] == 1 ) )') )
vardb.registerCut( Cut('TmuLel',  '( ( lep_flavour[0] == 13 && lep_isTightSelected[0] == 1 ) && ( lep_flavour[1] == 11 && lep_isTightSelected[1] == 0 ) )') )
vardb.registerCut( Cut('LmuTel',  '( ( lep_flavour[0] == 13 && lep_isTightSelected[0] == 0 ) && ( lep_flavour[1] == 11 && lep_isTightSelected[1] == 1 ) )') )
"""
#"""
vardb.registerCut( Cut('TT',      '( is_T_T == 1 )') )
vardb.registerCut( Cut('TL',      '( is_T_AntiT == 1 )') )
vardb.registerCut( Cut('LT',      '( is_AntiT_T == 1 )') )
vardb.registerCut( Cut('TL_LT',   '( is_T_AntiT == 1 || is_AntiT_T == 1 )') )
vardb.registerCut( Cut('LL',      '( is_AntiT_AntiT == 1 )') )
vardb.registerCut( Cut('TelLmu',  '( is_Tel_AntiTmu == 1 )') )
vardb.registerCut( Cut('LelTmu',  '( is_Lel_Tmu == 1 )') )
vardb.registerCut( Cut('TmuLel',  '( is_Tmu_AntiTel == 1 )') )
vardb.registerCut( Cut('LmuTel',  '( is_Lmu_Tel == 1 )') )
#"""

# -------------------------------------------------------------------------------
# the following cuts must be used only on MC :
#
#   -) plot only prompt-matched MC to avoid double counting of non prompt
#      (as they are estimated via MM or FF)
#   -) in case there are electrons in the regions, also brem electrons must be taken into account
#   -) what is left must be charge flip
#
# Use these cuts to plot specific MC contaminations in:
#
#   -) SR and low njet CR
#     (here you want the pure prompt MC contribution
#      to be plotted in SR/low-jet CR if
#      fakes and charge flips are already taken into account
#      by rescaling the data)
#
#   -) CR for FAKE rate measurement
#      (here you want to plot whatever is NON-non-prompt,
#       which then you will subtract to data for the fake rate estimate)
#
# Values of MC truth matching flags ( i.e., 'truthType', 'truthOrigin' ) are defined in MCTruthClassifier:
#
# https://svnweb.cern.ch/trac/atlasoff/browser/PhysicsAnalysis/MCTruthClassifier/trunk/MCTruthClassifier/MCTruthClassifierDefs.h
#
# -------------------------------------------------------------------------------

# 1.
# event passes this cut if ALL leptons are prompt (MCTruthClassifier --> Iso), and none is charge flip
#
vardb.registerCut( Cut('2Lep_PurePromptEvent', '( isMC==0 || ( isMC==1 && ( lep_truthType[0] == 6 || lep_truthType[0] == 2 ) && ( lep_truthType[1] == 6 || lep_truthType[1] == 2 ) && ( lep_isChFlip[0] == 0 && lep_isChFlip[1] == 0 ) ) )') )
# 2.
# event passes this cut if AT LEAST ONE lepton is !prompt (MCTruthClassifier --> !Iso), and none is charge flip
# (i.e., the !prompt lepton will be ( HF lepton || photon conv || lepton from Dalitz decay || mis-reco jet...)
# --> USED FOR CLOSURE TEST
#
vardb.registerCut( Cut('2Lep_NonPromptEvent', '( isMC==0 || ( isMC==1 && ( ( lep_truthType[0] != 6 && lep_truthType[0] != 2 ) || ( lep_truthType[1] != 6 && lep_truthType[1] != 2 ) ) && ( lep_isChFlip[0] == 0 && lep_isChFlip[1] == 0 ) ) )') )
#
# 3.
# event passes this cut if AT LEAST ONE lepton is charge flip (does not distinguish trident VS charge-misId)
#
vardb.registerCut( Cut('2Lep_ChFlipEvent',   '( isMC==0 || ( isMC==1 && ( lep_isChFlip[0] == 1 || lep_isChFlip[1] == 1 ) ) )') )
#  CHECK: should rather use this ???? vardb.registerCut( Cut('2Lep_ChFlipEvent',   '( isMC==0 || ( isMC==1 && ( ( lep_isChFlip[0] == 1 && lep_isChFlip[1] == 0 ) || ( lep_isChFlip[1] == 1 && lep_isChFlip[0] == 0 ) ) ) )') )

# WRONG definition in group ntuples:
#vardb.registerCut( Cut('2Lep_ChFlipEvent',   '( isMC==0 || ( isMC==1 && ( (lep_truthType[0] == 4 && lep_truthOrigin[0] == 5) || ((lep_truthType[1] == 4 && lep_truthOrigin[1] == 5)) ) ) )') )

# 3a.
# event passes this cut if AT LEAST ONE lepton is (prompt and charge flip) (it will be a charge-misId charge flip)
#
vardb.registerCut( Cut('2Lep_ChFlipPromptEvent',  '( isMC==0 || ( isMC==1 && ( ( lep_isChFlip[0] == 1 && ( lep_truthType[0] == 6 || lep_truthType[0] == 2 ) ) || ( lep_isChFlip[1] == 1 && ( lep_truthType[1] == 6 || lep_truthType[1] == 2 ) ) ) ) )') )
# 3b.
# event passes this cut if AT LEAST ONE object is (!prompt and charge flip and from bremsstrahlung) (this will be a trident charge flip)
#
vardb.registerCut( Cut('2Lep_ChFlipBremEvent', '( isMC==0 || ( isMC==1 && ( ( lep_isChFlip[0] == 1 && lep_isBrem[0] == 1 && ( lep_truthType[0] != 6 && lep_truthType[0] != 2 ) ) || ( lep_isChFlip[1] == 1 && lep_isBrem[1] == 1 && ( lep_truthType[1] != 6 && lep_truthType[1] != 2 ) ) ) ) )') )
# 3c.
# event passes this cut if AT LEAST ONE lepton is (!prompt and charge flip)
#
vardb.registerCut( Cut('2Lep_ChFlipNonPromptEvent', '( isMC==0 || ( isMC==1 && ( ( lep_isChFlip[0] == 1 && ( lep_truthType[0] != 6 && lep_truthType[0] != 2 ) ) || ( lep_isChFlip[1] == 1 && ( lep_truthType[1] != 6 && lep_truthType[1] != 2 ) ) ) ) )') )
# 4.
# event passes this cut if NONE of the leptons is charge flip
vardb.registerCut( Cut('2Lep_ChFlipVeto',   '( isMC==0 || ( isMC==1 && ( lep_isChFlip[0] == 0 && lep_isChFlip[1] == 0 ) ) )') )

# ------------------
# ttW Control Region
# ------------------

# >= 2 btag
# <= 4 jets
# HT(jets) > 220 GeV in ee and emu
# Z peak [75,105] veto in ee
# MET > 50 GeV in ee


# ---------------------------
# A list of variables to plot
# ---------------------------

# Reconstructed pT of the Z
#
pT_Z = '( TMath::Sqrt( (lep_pt[0]*lep_pt[0]) + (lep_pt[1]*lep_pt[1]) + 2*lep_pt[0]*lep_pt[1]*(TMath::Cos( lep_phi[0] - lep_phi[1] )) ) )/1e3'

# Calculate DeltaR(lep0,lep1) in 2LepSS + 0 tau category
#
gROOT.LoadMacro("$ROOTCOREBIN/user_scripts/HTopMultilepAnalysis/ROOT_TTreeFormulas/deltaR.cxx+")
from ROOT import deltaR
delta_R_lep0lep1 = 'deltaR( lep_eta[0], lep_phi[0], lep_eta[1], lep_phi[1] )'


if doSR or doLowNJetCR:
    print ''
    #vardb.registerVar( Variable(shortname = 'NJets', latexname = 'Jet multiplicity', ntuplename = 'njets', bins = 10, minval = -0.5, maxval = 9.5) )
    if doSR:
        vardb.registerVar( Variable(shortname = 'NJets4j5j', latexname = 'Jet multiplicity', ntuplename = 'njets', bins = 4, minval = 1.5, maxval = 5.5) )
    elif doLowNJetCR:
        vardb.registerVar( Variable(shortname = 'NJets2j3j', latexname = 'Jet multiplicity', ntuplename = 'njets', bins = 3, minval = 1.5, maxval = 4.5) )
    #vardb.registerVar( Variable(shortname = 'NBJets', latexname = 'BJet multiplicity', ntuplename = 'njets_mv2c20_Fix77', bins = 4, minval = -0.5, maxval = 3.5) )
    #vardb.registerVar( Variable(shortname = 'Mll01_inc', latexname = 'm(l_{0}l_{1}) [GeV]', ntuplename = 'mll01/1e3', bins = 13, minval = 0.0, maxval = 260.0,) )
    #vardb.registerVar( Variable(shortname = 'Lep0Pt', latexname = 'p_{T}^{lead lep} [GeV]', ntuplename = 'lep_pt[0]/1e3', bins = 11, minval = 20.0, maxval = 240.0,) )
    #vardb.registerVar( Variable(shortname = 'Lep0Eta', latexname = '|#eta^{lead lep}|', ntuplename = 'TMath::Abs(lep_eta[0])', bins = 8, minval = 0.0, maxval = 2.6) )
    #vardb.registerVar( Variable(shortname = 'deltaRLep0Lep1', latexname = '#DeltaR(lep_{0},lep_{1})', ntuplename = delta_R_lep0lep1, bins = 10, minval = 0.0, maxval = 5.0) )

if doMMRates or doMMClosureRates:
    print ''
    #vardb.registerVar( Variable(shortname = 'NJets', latexname = 'Jet multiplicity', ntuplename = 'njets', bins = 10, minval = 0, maxval = 10) )
    #vardb.registerVar( Variable(shortname = 'NBJets', latexname = 'BJet multiplicity', ntuplename = 'njets_mv2c20_Fix77', bins = 4, minval = 0, maxval = 4) )

if doMMRatesLHFit:
    vardb.registerVar( Variable(shortname = 'LepPt', latexname = 'p_{T}^{lep} [GeV]', ntuplename = 'lep_pt/1e3', bins = 15, minval = 10.0, maxval = 150.0,) )

if doMMClosureTest:
    print ''
    vardb.registerVar( Variable(shortname = 'NJets', latexname = 'Jet multiplicity', ntuplename = 'njets', bins = 10, minval = -0.5, maxval = 9.5) )
    #vardb.registerVar( Variable(shortname = 'NBJets', latexname = 'BJet multiplicity', ntuplename = 'njets_mv2c20_Fix77', bins = 4, minval = -0.5, maxval = 3.5) )
    vardb.registerVar( Variable(shortname = 'Mll01_inc', latexname = 'm(l_{0}l_{1}) [GeV]', ntuplename = 'mll01/1e3', bins = 13, minval = 0.0, maxval = 260.0,) )
    vardb.registerVar( Variable(shortname = 'Lep0Pt', latexname = 'p_{T}^{lead lep} [GeV]', ntuplename = 'lep_pt[0]/1e3', bins = 15, minval = 25.0, maxval = 175.0,) )
    vardb.registerVar( Variable(shortname = 'deltaRLep0Lep1', latexname = '#DeltaR(lep_{0},lep_{1})', ntuplename = delta_R_lep0lep1, bins = 20, minval = 0.0, maxval = 5.0) )

if doZSSpeakCR:
    print ''
    vardb.registerVar( Variable(shortname = 'Mll01_NarrowPeak', latexname = 'm(l_{0}l_{1}) [GeV]', ntuplename = 'mll01/1e3', bins = 26, minval = 60.0, maxval = 125.0) )

if doStandardPlots:
    print ''
    #vardb.registerVar( Variable(shortname = 'Jet0Pt', latexname = 'p_{T}^{lead jet} [GeV]', ntuplename = 'jet_pt[0]/1e3', bins = 36, minval = 20.0, maxval = 200.0,) )
    #vardb.registerVar( Variable(shortname = 'Jet0Eta', latexname = '#eta^{lead jet}', ntuplename = 'jet_eta[0]', bins = 50, minval = -5.0, maxval = 5.0) )
    vardb.registerVar( Variable(shortname = 'NJets', latexname = 'Jet multiplicity', ntuplename = 'njets', bins = 10, minval = -0.5, maxval = 9.5) )
    #vardb.registerVar( Variable(shortname = 'NBJets', latexname = 'BJet multiplicity', ntuplename = 'njets_mv2c20_Fix77', bins = 4, minval = -0.5, maxval = 3.5) )
    #vardb.registerVar( Variable(shortname = 'NJetsPlus10NBJets', latexname = 'N_{Jets}+10*N_{BJets}', ntuplename = 'njets+10.0*njets_mv2c20_Fix77', bins = 40, minval = 0, maxval = 40, basecut = vardb.getCut('VetoLargeNBJet')) )
    #
    # Inclusive m(ll) plot
    #
    #vardb.registerVar( Variable(shortname = 'Mll01_inc', latexname = 'm(l_{0}l_{1}) [GeV]', ntuplename = 'mll01/1e3', bins = 40, minval = 40.0, maxval = 240.0,) )
    #
    # Z peak plot
    #
    #vardb.registerVar( Variable(shortname = 'Mll01_peak', latexname = 'm(l_{0}l_{1}) [GeV]', ntuplename = 'mll01/1e3', bins = 40, minval = 40.0, maxval = 120.0,) )
    #
    #vardb.registerVar( Variable(shortname = 'pT_Z', latexname = 'p_{T} Z (reco) [GeV]', ntuplename = pT_Z, bins = 100, minval = 0.0, maxval = 1000.0, logaxisX = True) )
    #
    #vardb.registerVar( Variable(shortname = 'Lep0Pt', latexname = 'p_{T}^{lead lep} [GeV]', ntuplename = 'lep_pt[0]/1e3', bins = 36, minval = 10.0, maxval = 190.0,) )
    #vardb.registerVar( Variable(shortname = 'Lep1Pt', latexname = 'p_{T}^{2nd lead lep} [GeV]', ntuplename = 'lep_pt[1]/1e3', bins = 20, minval = 10.0, maxval = 110.0,) )
    #vardb.registerVar( Variable(shortname = 'Lep0Eta', latexname = '#eta^{lead lep}', ntuplename = 'lep_eta[0]', bins = 16, minval = -2.6, maxval = 2.6) )
    #vardb.registerVar( Variable(shortname = 'Lep1Eta', latexname = '#eta^{2nd lead lep}', ntuplename = 'lep_eta[1]', bins = 16, minval = -2.6, maxval = 2.6) )
    #vardb.registerVar( Variable(shortname = 'deltaRLep0Lep1', latexname = '#DeltaR(lep_{0},lep_{1})', ntuplename = delta_R_lep0lep1, bins = 20, minval = 0.0, maxval = 5.0) )
    #vardb.registerVar( Variable(shortname = 'Mll12', latexname = 'm(l_{1}l_{2}) [GeV]', ntuplename = 'mll12/1e3', bins = 15, minval = 0.0, maxval = 300.0,) )
    #vardb.registerVar( Variable(shortname = 'avgint', latexname = 'Average Interactions Per Bunch Crossing', ntuplename = 'averageInteractionsPerCrossing', bins = 50, minval = 0, maxval = 50, typeval = TH1I) )
    #vardb.registerVar( Variable(shortname = 'MET_FinalClus', latexname = 'E_{T}^{miss} (FinalClus) [GeV]', ntuplename = 'metFinalClus/1e3', bins = 45, minval = 0.0, maxval = 180.0,))
    #vardb.registerVar( Variable(shortname = 'MET_FinalTrk', latexname = 'E_{T}^{miss} (FinalTrk) [GeV]', ntuplename = 'metFinalTrk/1e3', bins = 45, minval = 0.0, maxval = 180.0,) )
    #vardb.registerVar( Variable(shortname = 'MET_SoftClus', latexname = 'E_{T}^{miss} (SoftClus) [GeV]', ntuplename = 'metSoftClus/1e3', bins = 45, minval = 0.0, maxval = 180.0,) )
    #vardb.registerVar( Variable(shortname = 'MET_SoftTrk', latexname = 'E_{T}^{miss} (SoftTrk) [GeV]', ntuplename = 'metSoftTrk/1e3', bins = 45, minval = 0.0, maxval = 180.0,) )
    #vardb.registerVar( Variable(shortname = 'MET_Electrons', latexname = 'E_{T}^{miss} (Electrons) [GeV]', ntuplename = 'metEle/1e3', bins = 45, minval = 0.0, maxval = 180.0,) )
    #vardb.registerVar( Variable(shortname = 'MET_Muons', latexname = 'E_{T}^{miss} (Muons) [GeV]', ntuplename = 'metMuons/1e3', bins = 45, minval = 0.0, maxval = 180.0,) )
    #vardb.registerVar( Variable(shortname = 'MET_Jets', latexname = 'E_{T}^{miss} (Jets) [GeV]', ntuplename = 'metJet/1e3', bins = 45, minval = 0.0, maxval = 180.0,) )

    #vardb.registerVar( Variable(shortname = 'MT_Lep0MET', latexname = 'm_{T}(l_{0},MET) [GeV]', ntuplename = 'mT_lep0MET/1e3', bins = 40, minval = 0.0, maxval = 160.0,) )
    #vardb.registerVar( Variable(shortname = 'MT_Lep1MET', latexname = 'm_{T}(l_{1},MET) [GeV]', ntuplename = 'mT_lep1MET/1e3', bins = 40, minval = 0.0, maxval = 160.0,) )
    #vardb.registerVar( Variable(shortname = 'Tau0Pt', latexname = 'p_{T}^{lead tau} [GeV]', ntuplename = 'tau_pt[0]/1e3', bins = 30, minval = 25.0, maxval = 100.0,) )

    #vardb.registerVar( Variable(shortname = 'El0Pt', latexname = 'p_{T}^{lead e} [GeV]', ntuplename = 'el_pt[0]/1e3', bins = 36, minval = 10.0, maxval = 190.0,) )
    #vardb.registerVar( Variable(shortname = 'El1Pt', latexname = 'p_{T}^{2nd lead e} [GeV]', ntuplename = 'el_pt[1]/1e3', bins = 20, minval = 10.0, maxval = 110.0,) )
    #vardb.registerVar( Variable(shortname = 'El0Eta', latexname = '#eta^{lead e}', ntuplename = 'el_caloCluster_eta[0]', bins = 16, minval = -2.6, maxval = 2.6, manualbins = [ -2.6, -2.25, -2.0, -1.52, -1.37, -1.1, -0.8, -0.5, 0.0, 0.5, 0.8, 1.1, 1.37, 1.52, 2.0, 2.25, 2.6]) )
    #vardb.registerVar( Variable(shortname = 'El1Eta', latexname = '#eta^{2nd lead e}', ntuplename = 'el_caloCluster_eta[1]', bins = 16, minval = -2.6, maxval = 2.6, manualbins = [ -2.6, -2.25, -2.0, -1.52, -1.37, -1.1, -0.8, -0.5, 0.0, 0.5, 0.8, 1.1, 1.37, 1.52, 2.0, 2.25, 2.6]) )
    #vardb.registerVar( Variable(shortname = 'El0TopoEtCone20', latexname = 'topoetcone20^{lead e} [GeV]', ntuplename = 'el_topoetcone20[0]/1e3', bins = 40, minval = 0.0, maxval = 10.0, manualbins = [ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0] ) )
    #vardb.registerVar( Variable(shortname = 'El1TopoEtCone20', latexname = 'topoetcone20^{2nd lead e} [GeV]', ntuplename = 'el_topoetcone20[1]/1e3', bins = 40, minval = 0.0, maxval = 10.0, manualbins = [ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0] ) )
    #vardb.registerVar( Variable(shortname = 'El0PtVarCone20', latexname = 'ptvarcone20^{lead e} [GeV]', ntuplename = 'el_ptvarcone20[0]/1e3', bins = 40, minval = 1.0, maxval = 5.0) )
    #vardb.registerVar( Variable(shortname = 'El1PtVarCone20', latexname = 'ptvarcone20^{2nd lead e} [GeV]', ntuplename = 'el_ptvarcone20[1]/1e3', bins = 40, minval = 1.0, maxval = 5.0) )
    #vardb.registerVar( Variable(shortname = 'El0TopoEtCone20OverPt', latexname = 'topoetcone20/p_{T} lead e [GeV]', ntuplename = 'el_topoetcone20[0]/el_pt[0]', bins = 50, minval = -0.2, maxval = 0.8) )
    #vardb.registerVar( Variable(shortname = 'El1TopoEtCone20OverPt', latexname = 'topoetcone20/p_{T} 2nd lead e [GeV]', ntuplename = 'el_topoetcone20[1]/el_pt[1]', bins = 50, minval = -0.2, maxval = 0.8) )
    #vardb.registerVar( Variable(shortname = 'El0PtVarCone20OverPt', latexname = 'ptvarcone20/p_{T} lead e [GeV]', ntuplename = 'el_ptvarcone20[0]/el_pt[0]', bins = 50, minval = 0.0, maxval = 1.0) )
    #vardb.registerVar( Variable(shortname = 'El1PtVarCone20OverPt', latexname = 'ptvarcone20/p_{T} 2nd lead e [GeV]', ntuplename = 'el_ptvarcone20[1]/el_pt[1]', bins = 50, minval = 0.0, maxval = 1.0) )
    #vardb.registerVar( Variable(shortname = 'El0d0sig', latexname = '|d_{0}^{sig}| lead e', ntuplename = 'el_trkd0sig[0]', bins = 40, minval = 0.0, maxval = 10.0,) )
    #vardb.registerVar( Variable(shortname = 'El1d0sig', latexname = '|d_{0}^{sig}| 2nd lead e', ntuplename = 'el_trkd0sig[1]', bins = 40, minval = 0.0, maxval = 10.0,) )
    #vardb.registerVar( Variable(shortname = 'El0z0sintheta', latexname = 'z_{0}*sin(#theta) lead e [mm]', ntuplename = 'el_trkz0sintheta[0]', bins = 20, minval = -1.0, maxval = 1.0,) )
    #vardb.registerVar( Variable(shortname = 'El1z0sintheta', latexname = 'z_{0}*sin(#theta) 2nd lead e [mm]', ntuplename = 'el_trkz0sintheta[1]', bins = 20, minval = -1.0, maxval = 1.0,) )
    #vardb.registerVar( Variable(shortname = 'El0isProbe', latexname = 'lead e isProbe', ntuplename = '!el_isTag[0]', bins = 2, minval = -0.5, maxval = 1.5) )
    #vardb.registerVar( Variable(shortname = 'El1isProbe', latexname = '2nd lead e isProbe', ntuplename = '!el_isTag[1]', bins = 2, minval = -0.5, maxval = 1.5) )


    #vardb.registerVar( Variable(shortname = 'Mu0Pt', latexname = 'p_{T}^{lead #mu} [GeV]', ntuplename = 'muon_pt[0]/1e3', bins = 36, minval = 10.0, maxval = 190.0,) )
    #vardb.registerVar( Variable(shortname = 'Mu1Pt', latexname = 'p_{T}^{2nd lead #mu} [GeV]', ntuplename = 'muon_pt[1]/1e3', bins = 20, minval = 10.0, maxval = 110.0,) )
    #vardb.registerVar( Variable(shortname = 'Mu0Eta', latexname = '#eta^{lead #mu}', ntuplename = 'muon_eta[0]', bins = 16, minval = -2.5, maxval = 2.5, manualbins = [-2.5, -2.2, -1.9, -1.6, -1.3, -1.0, -0.7, -0.4, -0.1, 0.0 , 0.1 , 0.4 , 0.7, 1.0,  1.3 , 1.6 , 1.9, 2.2, 2.5 ]) )
    #vardb.registerVar( Variable(shortname = 'Mu1Eta', latexname = '#eta^{2nd lead #mu}', ntuplename = 'muon_eta[1]', bins = 16, minval = -2.5, maxval = 2.5, manualbins = [-2.5, -2.2, -1.9, -1.6, -1.3, -1.0, -0.7, -0.4, -0.1, 0.0 , 0.1 , 0.4 , 0.7, 1.0,  1.3 , 1.6 , 1.9, 2.2, 2.5 ]) )
    #vardb.registerVar( Variable(shortname = 'Mu0TopoEtCone20', latexname = 'topoetcone20^{lead #mu} [GeV]', ntuplename = 'muon_topoetcone20[0]/1e3', bins = 40, minval = 0.0, maxval = 10.0, manualbins = [ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0]) )
    #vardb.registerVar( Variable(shortname = 'Mu1TopoEtCone20', latexname = 'topoetcone20^{2nd lead #mu} [GeV]', ntuplename = 'muon_topoetcone20[1]/1e3', bins = 40, minval = 0.0, maxval = 10.0, manualbins = [ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0]) )
    #vardb.registerVar( Variable(shortname = 'Mu0PtVarCone30', latexname = 'ptvarcone20^{lead #mu} [GeV]', ntuplename = 'muon_ptvarcone30[0]/1e3', bins = 40, minval = 1.0, maxval = 5.0) )
    #vardb.registerVar( Variable(shortname = 'Mu1PtVarCone30', latexname = 'ptvarcone20^{2nd lead #mu} [GeV]', ntuplename = 'muon_ptvarcone30[1]/1e3', bins = 40, minval = 1.0, maxval = 5.0) )
    #vardb.registerVar( Variable(shortname = 'Mu0TopoEtCone20OverPt', latexname = 'topoetcone20/p_{T} lead #mu [GeV]', ntuplename = 'muon_topoetcone20[0]/muon_pt[0]', bins = 50, minval = -0.2, maxval = 0.8) )
    #vardb.registerVar( Variable(shortname = 'Mu1TopoEtCone20OverPt', latexname = 'topoetcone20/p_{T} 2nd lead #mu [GeV]', ntuplename = 'muon_topoetcone20[1]/muon_pt[1]', bins = 50, minval = -0.2, maxval = 0.8) )
    #vardb.registerVar( Variable(shortname = 'Mu0PtVarCone30OverPt', latexname = 'ptvarcone30/p_{T} lead #mu [GeV]', ntuplename = 'muon_ptvarcone30[0]/muon_pt[0]', bins = 50, minval = 0.0, maxval = 1.0) )
    #vardb.registerVar( Variable(shortname = 'Mu1PtVarCone30OverPt', latexname = 'ptvarcone30/p_{T} 2nd lead #mu [GeV]', ntuplename = 'muon_ptvarcone30[1]/muon_pt[1]', bins = 50, minval = 0.0, maxval = 1.0) )
    #vardb.registerVar( Variable(shortname = 'Mu0d0sig', latexname = '|d_{0}^{sig}| lead #mu', ntuplename = 'muon_trkd0sig[0]', bins = 40, minval = 0.0, maxval = 10.0,) )
    #vardb.registerVar( Variable(shortname = 'Mu1d0sig', latexname = '|d_{0}^{sig}| 2nd lead #mu', ntuplename = 'muon_trkd0sig[1]', bins = 40, minval = 0.0, maxval = 10.0,) )
    #vardb.registerVar( Variable(shortname = 'Mu0z0sintheta', latexname = 'z_{0}*sin(#theta) lead #mu [mm]', ntuplename = 'muon_trkz0sintheta[0]', bins = 20, minval = -1.0, maxval = 1.0,) )
    #vardb.registerVar( Variable(shortname = 'Mu1z0sintheta', latexname = 'z_{0}*sin(#theta) 2nd lead #mu [mm]', ntuplename = 'muon_trkz0sintheta[1]', bins = 20, minval = -1.0, maxval = 1.0,) )
    #vardb.registerVar( Variable(shortname = 'Mu0isProbe', latexname = 'lead #mu isProbe', ntuplename = '!muon_isTag[0]', bins = 2, minval = -0.5, maxval = 1.5) )
    #vardb.registerVar( Variable(shortname = 'Mu1isProbe', latexname = '2nd lead #mu isProbe', ntuplename = '!muon_isTag[1]', bins = 2, minval = -0.5, maxval = 1.5) )

# -------------------------------------------------
# Alterantive ranges and binning for the histograms
# -------------------------------------------------
midstatsbin = {
    'MMC': (25, 0., 250.),
    'mvis': (25, 0., 250.),
    'mT': (30, 0., 120.),
    'MET': (25, 0., 100.),
    'leppt': (30, 17., 77.),
    'taupt': (25, 20., 70.),
    'jetpt': (25, 25., 125.),
}
lowstatsbin = {
    'MMC': (12, 0., 240.),
    'mvis': (12, 0., 240.),
    'mT': (12, 0., 120.),
    'MET': (12, 0., 120.),
    'leppt': (12, 17., 77.),
    'taupt': (12, 20., 80.),
    'jetpt': (12, 25., 121.),
}

# ---------------------
# A list of systematics
# ---------------------
if args.doSyst:
    # if doTwoLepSR or doTwoLepLowNJetCR or dottWCR:
    #	 vardb.registerSystematics( Systematics(name='CFsys',	   eventweight='sys_weight_CF_') ) ## uncertainties on the kfactors used to normalize the various MC distributions
    if doMM:
        #vardb.registerSystematics( Systematics(name='MMrsys',	   eventweight='sys_weight_MMr_') )
    	#vardb.registerSystematics( Systematics(name='MMfsys',	   eventweight='sys_weight_MMf_') )
        vardb.registerSystematics( Systematics(name='MMrsys',	   eventweight='MMWeight') )
    	vardb.registerSystematics( Systematics(name='MMfsys',	   eventweight='MMWeight') )
    if doFF:
    	#vardb.registerSystematics( Systematics(name='FFsys',	   eventweight='sys_weight_FF_') )
    	vardb.registerSystematics( Systematics(name='FFsys',	   eventweight='FFWeight') )

    '''
    vardb.registerSystematics( Systematics(name='PU',		  eventweight='evtsel_sys_PU_rescaling_') )
    vardb.registerSystematics( Systematics(name='el_reco',	  eventweight='evtsel_sys_sf_el_reco_') )
    vardb.registerSystematics( Systematics(name='el_id',	  eventweight='evtsel_sys_sf_el_id_') )
    vardb.registerSystematics( Systematics(name='el_iso',	  eventweight='evtsel_sys_sf_el_iso_') )
    vardb.registerSystematics( Systematics(name='mu_id',	  eventweight='evtsel_sys_sf_mu_id_') )
    vardb.registerSystematics( Systematics(name='mu_iso',	  eventweight='evtsel_sys_sf_mu_iso_') )
    vardb.registerSystematics( Systematics(name='lep_trig',	  eventweight='evtsel_sys_sf_lep_trig_') )
    vardb.registerSystematics( Systematics(name='bjet_b',	  eventweight='evtsel_sys_sf_bjet_b_') )
    vardb.registerSystematics( Systematics(name='bjet_c',	  eventweight='evtsel_sys_sf_bjet_c_') )
    vardb.registerSystematics( Systematics(name='bjet_m',	  eventweight='evtsel_sys_sf_bjet_m_') )

    vardb.registerSystematics( Systematics(name='METSys',	  treename='METSys') )
    vardb.registerSystematics( Systematics(name='ElEnResSys',	  treename='ElEnResSys') )
    vardb.registerSystematics( Systematics(name='ElES_LowPt',	  treename='ElES_LowPt') )
    vardb.registerSystematics( Systematics(name='ElES_Zee',	  treename='ElES_Zee') )
    vardb.registerSystematics( Systematics(name='ElES_R12',	  treename='ElES_R12') )
    vardb.registerSystematics( Systematics(name='ElES_PS',	  treename='ElES_PS') )
    vardb.registerSystematics( Systematics(name='EESSys',	  treename='EESSys') )
    vardb.registerSystematics( Systematics(name='MuSys',	  treename='MuSys') )
    vardb.registerSystematics( Systematics(name='JES_Total',	  treename='JES_Total') )
    vardb.registerSystematics( Systematics(name='JER',  	  treename='JER') )
    '''
# -------------------------------------------------------------------
# Definition of the categories for which one wants produce histograms
# -------------------------------------------------------------------

# ------------
# SRs
# ------------
if doSR:

    if doTwoLepSR :

    	# when using MM or FF or ABCD for non-prompt bkg estimate, make sure you plot
    	# only pure prompt MC (to avoid double counting of fake background events!)
    	#
    	if ( args.fakeMethod == 'MM' or args.fakeMethod == 'FF' or args.fakeMethod == 'ABCD' ):
            # MuMu region
            #
            vardb.registerCategory( MyCategory('MuMuSS_SR_DataDriven',  cut = vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet_SR','2Lep_NLep','2Lep_SS','2Lep_MuMu_Event','TauVeto','2Lep_PurePromptEvent','2Lep_NJet_SR']) ) )
            # ElEl region
            #
            vardb.registerCategory( MyCategory('ElElSS_SR_DataDriven',  cut = vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet_SR','2Lep_NLep','2Lep_SS','2Lep_ElEl_Event','TauVeto','2Lep_PurePromptEvent','2Lep_ElEtaCut','2Lep_NJet_SR']) ) )
            # OF region
            #
            vardb.registerCategory( MyCategory('OFSS_SR_DataDriven',	cut = vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet_SR','2Lep_NLep','2Lep_SS','2Lep_OF_Event','TauVeto','2Lep_PurePromptEvent','2Lep_ElEtaCut','2Lep_NJet_SR']) ) )
    	else:
            vardb.registerCategory( MyCategory('MuMuSS_SR',		cut = vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet_SR','2Lep_NLep','2Lep_SS','2Lep_MuMu_Event','TauVeto','2Lep_NJet_SR']) ) )
            vardb.registerCategory( MyCategory('OFSS_SR',		cut = vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet_SR','2Lep_NLep','2Lep_SS','2Lep_OF_Event','TauVeto','2Lep_ElEtaCut','2Lep_NJet_SR']) ) )
            vardb.registerCategory( MyCategory('ElElSS_SR',		cut = vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet_SR','2Lep_NLep','2Lep_SS','2Lep_ElEl_Event','TauVeto','2Lep_ElEtaCut','2Lep_NJet_SR']) ) )
	    #
	    # Plot the inclusive flavour SR, by looking only at fake backgrounds (use truth)
	    #
            #vardb.registerCategory( MyCategory('InclusiveFlavourSS_SR',	 cut = vardb.getCuts(['TrigDec', '2Lep_TrigMatch', '2Lep_NBJet_SR', '2Lep_NLep', '2Lep_NonPromptEvent', '2Lep_SS', 'TauVeto', '2Lep_ElEtaCut', '2Lep_NJet_SR']) ) )
            #
            # 2lep+tau region
            #
            vardb.registerCategory( MyCategory('TwoLepSSTau_SR',	cut = vardb.getCuts(['TrigDec','2Lep1Tau_NLep','2Lep1Tau_pT','2Lep1Tau_TrigMatch','2Lep1Tau_SS','2Lep1Tau_1Tau','2Lep1Tau_Zsidescut','2Lep1Tau_NJet_SR','2Lep1Tau_NBJet']) ) )

    if doThreeLepSR:
    	vardb.registerCategory( MyCategory('ThreeLep_SR',    cut = vardb.getCuts(['TrigDec','3Lep_TrigMatch','3Lep_NLep','3Lep_Charge','3Lep_TightLeptons','3Lep_ZVeto','3Lep_MinZCut','3Lep_NJets']) ) )

    if doFourLepSR:
        vardb.registerCategory( MyCategory('FourLep_SR',     cut = vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet_SR','4Lep','4Lep_NJets']) ) )

# -------------
# low N-jet CRs
# -------------
if doTwoLepLowNJetCR :

    if ( args.fakeMethod == 'MM' or args.fakeMethod == 'FF' ):
        # MuMu region
        #
        vardb.registerCategory( MyCategory('MuMuSS_LowNJetCR_DataDriven',  cut = vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet','2Lep_NLep','2Lep_SS','2Lep_MuMu_Event','TauVeto','2Lep_PurePromptEvent','2Lep_NJet_CR']) ) )
        # OF region
        #
        vardb.registerCategory( MyCategory('OFSS_LowNJetCR_DataDriven',    cut = vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet','2Lep_NLep','2Lep_SS','2Lep_OF_Event','TauVeto','2Lep_PurePromptEvent','2Lep_ElEtaCut','2Lep_NJet_CR']) ) )
        # ElEl region
        #
        vardb.registerCategory( MyCategory('ElElSS_LowNJetCR_DataDriven',  cut = vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet','2Lep_NLep','2Lep_SS','2Lep_ElEl_Event','TauVeto','2Lep_PurePromptEvent','2Lep_ElEtaCut','2Lep_NJet_CR']) ) )

    elif ( args.fakeMethod == 'ABCD' ):
        # OF region
        #
        vardb.registerCategory( MyCategory('OFSS_LowNJetCR_DataDriven',    cut = vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet','2Lep_NLep','2Lep_SS','2Lep_OF_Event','TauVeto','2Lep_PurePromptEvent','2Lep_ElEtaCut','2Lep_NJet_CR']) ) )
    else:
    	vardb.registerCategory( MyCategory('MuMuSS_LowNJetCR',	     cut = vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet','2Lep_NLep','2Lep_SS','2Lep_MuMu_Event','TauVeto','2Lep_NJet_CR']) ) )
    	vardb.registerCategory( MyCategory('OFSS_LowNJetCR',	     cut = vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet','2Lep_NLep','2Lep_SS','2Lep_OF_Event','TauVeto','2Lep_ElEtaCut','2Lep_NJet_CR']) ) )
    	vardb.registerCategory( MyCategory('ElElSS_LowNJetCR',	     cut = vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet','2Lep_NLep','2Lep_SS','2Lep_ElEl_Event','TauVeto','2Lep_ElEtaCut','2Lep_NJet_CR']) ) )
        #
        # 2lep+tau region
        #
        vardb.registerCategory( MyCategory('TwoLepSSTau_LowNJetCR',  cut = vardb.getCuts(['TrigDec','2Lep1Tau_NLep','2Lep1Tau_pT','2Lep1Tau_TrigMatch','2Lep1Tau_SS','2Lep1Tau_1Tau','2Lep1Tau_Zsidescut','2Lep1Tau_NJet_CR','2Lep1Tau_NBJet']) ) )

if doThreeLepLowNJetCR:
    # take OS pairs
    #
    vardb.registerCategory( MyCategory('ThreeLep_LowNJetCR',   cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet','3Lep_NLep','2Lep_NJet_CR','ZOSsidescut','2Lep_OS']) ) ) )

# -------------
# other CRs
# -------------

if doWZonCR:
    vardb.registerCategory( MyCategory('WZonCR',      cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','BJetVeto','3Lep_NLep','2Lep_Zpeakcut','2Lep_OS']) ) ) )

if doWZoffCR:
    vardb.registerCategory( MyCategory('WZoffCR',     cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','BJetVeto','3Lep_NLep','2Lep_Zsidescut','2Lep_OS']) ) ) )

if doWZHFonCR:
    vardb.registerCategory( MyCategory('WZHFonCR',    cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet','3Lep_NLep','2Lep_Zpeakcut','2Lep_OS']) ) ) )

if doWZHFoffCR:
    vardb.registerCategory( MyCategory('WZHFoffCR',   cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet','3Lep_NLep','2Lep_Zsidescut','2Lep_OS']) ) ) )

if dottZCR:
    vardb.registerCategory( MyCategory('ttZCR',       cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet','3Lep_NLep','NJet3L','2Lep_Zpeakcut','2Lep_OS']) ) ) )

if dottWCR:
    #vardb.registerCategory( MyCategory('ttWCR',       cut =   vardb.getCuts(['TrigDec', '2Lep_TrigMatch', 'LargeNBjet', '2Lep_NLep',  '2Lep_NJet_CR',  'TauVeto', '2Lep_SS']) ) )
    vardb.registerCategory( MyCategory('ttWCR',       cut = vardb.getCuts(['TrigDec', '2Lep_TrigMatch', '2Lep_NJet_CR', 'LargeNBJet', '2Lep_NLep_Relaxed', 'TauVeto']) ) )
    vardb.registerCategory( MyCategory('ttWCR_TT',    cut = vardb.getCuts(['TrigDec', '2Lep_TrigMatch', '2Lep_NJet_CR', 'LargeNBJet', '2Lep_NLep_Relaxed', 'TauVeto', 'TT']) ) )

if doZSSpeakCR:
    vardb.registerCategory( MyCategory('ZSSpeakCR_ElEl',   cut = vardb.getCuts(['2Lep_NLep','TrigDec','2Lep_TrigMatch','TauVeto','2Lep_ElEl_Event','2Lep_SS','2Lep_Zpeakcut','2Lep_PurePromptEvent']) ) )
    vardb.registerCategory( MyCategory('ZSSpeakCR_MuMu',   cut = vardb.getCuts(['2Lep_NLep','TrigDec','2Lep_TrigMatch','TauVeto','2Lep_MuMu_Event','2Lep_SS','2Lep_Zpeakcut','2Lep_PurePromptEvent']) ) )

# ------------------------------------
# Special CR for Data/MC control plots
# ------------------------------------

if doDataMCCR:

    # ----------------------------------------------------
    # Inclusive OS dilepton (ee,mumu, emu)
    #
    vardb.registerCategory( MyCategory('DataMC_InclusiveOS_MuMu', cut = ( vardb.getCuts(['2Lep_NLep_Relaxed','TrigDec','2Lep_TrigMatchDataMC','2Lep_MuMu_Event','TT','2Lep_ElEtaCut','2Lep_Zmincut','2Lep_OS']) ) ) )
    vardb.registerCategory( MyCategory('DataMC_InclusiveOS_ElEl', cut = ( vardb.getCuts(['2Lep_NLep_Relaxed','TrigDec','2Lep_TrigMatchDataMC','2Lep_ElEl_Event','TT','2Lep_ElEtaCut','2Lep_Zmincut','2Lep_OS']) ) ) )
    vardb.registerCategory( MyCategory('DataMC_InclusiveOS_OF',   cut = ( vardb.getCuts(['2Lep_NLep_Relaxed','TrigDec','2Lep_TrigMatchDataMC','2Lep_OF_Event','TT','2Lep_ElEtaCut','2Lep_Zmincut','2Lep_OS']) ) ) )

    # ----------------------------------------------------
    # OS ttbar ( top dilepton) (ee,mumu,emu)
    #
    vardb.registerCategory( MyCategory('DataMC_OS_ttbar_MuMu',    cut = ( vardb.getCuts(['2Lep_NLep_Relaxed','TrigDec','2Lep_TrigMatchDataMC','2Lep_MuMu_Event','TT','2Lep_ElEtaCut','4Lep_NJets','2Lep_NBJet','2Lep_Zsidescut','2Lep_Zmincut','2Lep_OS']) ) ) )
    vardb.registerCategory( MyCategory('DataMC_OS_ttbar_ElEl',    cut = ( vardb.getCuts(['2Lep_NLep_Relaxed','TrigDec','2Lep_TrigMatchDataMC','2Lep_ElEl_Event','TT','2Lep_ElEtaCut','4Lep_NJets','2Lep_NBJet','2Lep_Zsidescut','2Lep_Zmincut','2Lep_OS']) ) ) )
    vardb.registerCategory( MyCategory('DataMC_OS_ttbar_OF',	  cut = ( vardb.getCuts(['2Lep_NLep_Relaxed','TrigDec','2Lep_TrigMatchDataMC','2Lep_OF_Event','TT','2Lep_ElEtaCut','4Lep_NJets','2Lep_NBJet','2Lep_Zsidescut','2Lep_Zmincut','2Lep_OS']) ) ) )

    # ----------------------------------------------------
    # SS ttbar (ee,mumu,emu)
    #
    vardb.registerCategory( MyCategory('DataMC_SS_ttbar', 	  cut = ( vardb.getCuts(['2Lep_NLep_Relaxed','TrigDec','2Lep_TrigMatchDataMC','2Lep_NJet_CR_SStt','TT','2Lep_ElEtaCut','OneBJet','2Lep_SS','2Lep_Zmincut']) ) ) )

    # this is the Real CR for MM
    #
    vardb.registerCategory( MyCategory('DataMC_MuMu_RealCR',	  cut = ( vardb.getCuts(['2Lep_NLep_MMRates','TrigDec','2Lep_TrigMatch','2Lep_MuMu_Event','2Lep_NBJet','2Lep_NJet_CR','2Lep_LepTagTightTrigMatched','2Lep_OS','2Lep_Zmincut']) ) ) )
    vardb.registerCategory( MyCategory('DataMC_ElEl_RealCR',	  cut = ( vardb.getCuts(['2Lep_NLep_MMRates','TrigDec','2Lep_TrigMatch','2Lep_ElEl_Event','2Lep_NBJet','2Lep_NJet_CR','2Lep_LepTagTightTrigMatched','2Lep_OS','2Lep_Zmincut']) ) ) )
    vardb.registerCategory( MyCategory('DataMC_OF_RealCR',        cut = ( vardb.getCuts(['2Lep_NLep_MMRates','TrigDec','2Lep_TrigMatch','2Lep_OF_Event','2Lep_NBJet','2Lep_NJet_CR','2Lep_LepTagTightTrigMatched','2Lep_OS']) ) ) )
    #
    # this is the Fake CR for MM
    #
    vardb.registerCategory( MyCategory('DataMC_MuMu_FakeCR',	  cut = vardb.getCuts(['2Lep_NLep_MMRates','TrigDec','2Lep_TrigMatch','2Lep_MuMu_Event','2Lep_NBJet','2Lep_NJet_CR','2Lep_LepTagTightTrigMatched','2Lep_SS','2Lep_Zsidescut','2Lep_ElTagEtaCut','2Lep_Zmincut'])  ) )
    vardb.registerCategory( MyCategory('DataMC_ElEl_FakeCR',	  cut = vardb.getCuts(['2Lep_NLep_MMRates','TrigDec','2Lep_TrigMatch','2Lep_ElEl_Event','2Lep_NBJet','2Lep_NJet_CR','2Lep_LepTagTightTrigMatched','2Lep_SS','2Lep_Zsidescut','2Lep_ElTagEtaCut','2Lep_Zmincut'])  ) )
    vardb.registerCategory( MyCategory('DataMC_OF_FakeCR',        cut = vardb.getCuts(['2Lep_NLep_MMRates','TrigDec','2Lep_TrigMatch','2Lep_OF_Event','2Lep_NBJet','2Lep_NJet_CR','2Lep_LepTagTightTrigMatched','2Lep_SS','2Lep_Zsidescut','2Lep_ElTagEtaCut'])  ) )

# --------------------------------------------
# Full breakdown of cuts for cutflow challenge
# --------------------------------------------

if doCFChallenge:

    vardb.registerCut( Cut('2Lep_JustNLep',  '( nlep == 2 )') )
    vardb.registerCut( Cut('2Lep_pT',        '( Min$( lep_pt ) > 25e3 )') )

    vardb.registerCut( Cut('3Lep_JustNLep',  '( nlep == 3 )') )
    vardb.registerCut( Cut('3Lep_pT',        '( lep_3lepClosestSS_pt[0] > 20e3 && lep_3lepOtherSS_pt[0] > 20e3 && lep_3lepOS_pt[0] > 10e3 )') )

    # CF Challenge for MM rates measurement
    #
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep_MMRates_NLep',		  cut = vardb.getCuts(['2Lep_JustNLep'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep_MMRates_TrigMatch',		  cut = vardb.getCuts(['2Lep_JustNLep','2Lep_TrigMatch'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep_MMRates_TauVeto',		  cut = vardb.getCuts(['2Lep_JustNLep','2Lep_TrigMatch','TauVeto'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep_MMRates_NJets',		  cut = vardb.getCuts(['2Lep_JustNLep','2Lep_TrigMatch','TauVeto','2Lep_NJet_CR'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep_MMRates_NBJets',		  cut = vardb.getCuts(['2Lep_JustNLep','2Lep_TrigMatch','TauVeto','2Lep_NJet_CR','2Lep_NBJet'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep_MMRates_LepTagTightTrigMatched', cut = vardb.getCuts(['2Lep_JustNLep','2Lep_TrigMatch','TauVeto','2Lep_NJet_CR','2Lep_NBJet','2Lep_LepTagTightTrigMatched'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep_MMRates_OS',  		  cut = vardb.getCuts(['2Lep_JustNLep','2Lep_TrigMatch','TauVeto','2Lep_NJet_CR','2Lep_NBJet','2Lep_LepTagTightTrigMatched','2Lep_OS'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep_MMRates_OS_ProbeT',		  cut = vardb.getCuts(['2Lep_JustNLep','2Lep_TrigMatch','TauVeto','2Lep_NJet_CR','2Lep_NBJet','2Lep_LepTagTightTrigMatched','2Lep_OS','2Lep_ProbeTight'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep_MMRates_OS_ProbeAntiT',	  cut = vardb.getCuts(['2Lep_JustNLep','2Lep_TrigMatch','TauVeto','2Lep_NJet_CR','2Lep_NBJet','2Lep_LepTagTightTrigMatched','2Lep_OS','2Lep_ProbeAntiTight'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep_MMRates_SS',  		  cut = vardb.getCuts(['2Lep_JustNLep','2Lep_TrigMatch','TauVeto','2Lep_NJet_CR','2Lep_NBJet','2Lep_LepTagTightTrigMatched','2Lep_SS'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep_MMRates_SS_Zveto',		  cut = vardb.getCuts(['2Lep_JustNLep','2Lep_TrigMatch','TauVeto','2Lep_NJet_CR','2Lep_NBJet','2Lep_LepTagTightTrigMatched','2Lep_SS','2Lep_Zsidescut'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep_MMRates_SS_ElEl',		  cut = vardb.getCuts(['2Lep_JustNLep','2Lep_TrigMatch','TauVeto','2Lep_NJet_CR','2Lep_NBJet','2Lep_LepTagTightTrigMatched','2Lep_SS','2Lep_Zsidescut', '2Lep_ElEl_Event'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep_MMRates_SS_MuMu',		  cut = vardb.getCuts(['2Lep_JustNLep','2Lep_TrigMatch','TauVeto','2Lep_NJet_CR','2Lep_NBJet','2Lep_LepTagTightTrigMatched','2Lep_SS','2Lep_Zsidescut', '2Lep_MuMu_Event'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep_MMRates_SS_ElTagMuProbe',	  cut = vardb.getCuts(['2Lep_JustNLep','2Lep_TrigMatch','TauVeto','2Lep_NJet_CR','2Lep_NBJet','2Lep_LepTagTightTrigMatched','2Lep_SS','2Lep_Zsidescut', '2Lep_OF_Event', '2Lep_ProbeMu'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep_MMRates_SS_MuTagElProbe',	  cut = vardb.getCuts(['2Lep_JustNLep','2Lep_TrigMatch','TauVeto','2Lep_NJet_CR','2Lep_NBJet','2Lep_LepTagTightTrigMatched','2Lep_SS','2Lep_Zsidescut', '2Lep_OF_Event', '2Lep_ProbeEl'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep_MMRates_SS_Eta',		  cut = vardb.getCuts(['2Lep_JustNLep','2Lep_TrigMatch','TauVeto','2Lep_NJet_CR','2Lep_NBJet','2Lep_LepTagTightTrigMatched','2Lep_SS','2Lep_Zsidescut','2Lep_ElTagEtaCut'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep_MMRates_SS_ProbeT',		  cut = vardb.getCuts(['2Lep_JustNLep','2Lep_TrigMatch','TauVeto','2Lep_NJet_CR','2Lep_NBJet','2Lep_LepTagTightTrigMatched','2Lep_SS','2Lep_Zsidescut','2Lep_ElTagEtaCut','2Lep_ProbeTight'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep_MMRates_SS_ProbeAntiT',	  cut = vardb.getCuts(['2Lep_JustNLep','2Lep_TrigMatch','TauVeto','2Lep_NJet_CR','2Lep_NBJet','2Lep_LepTagTightTrigMatched','2Lep_SS','2Lep_Zsidescut','2Lep_ElTagEtaCut','2Lep_ProbeAntiTight'])  ) )

    # 2lepSS + 0tau
    #
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep_NLep', 	 cut = vardb.getCuts(['2Lep_JustNLep'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep_TT',           cut = vardb.getCuts(['2Lep_JustNLep','TT'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep_pT', 	         cut = vardb.getCuts(['2Lep_JustNLep','TT','2Lep_pT'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep_eta', 	 cut = vardb.getCuts(['2Lep_JustNLep','TT','2Lep_pT','2Lep_ElEtaCut'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep_TrigMatch', 	 cut = vardb.getCuts(['2Lep_JustNLep','TT','2Lep_pT','2Lep_ElEtaCut','2Lep_TrigMatch'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep_SS',  	 cut = vardb.getCuts(['2Lep_JustNLep','TT','2Lep_pT','2Lep_ElEtaCut','2Lep_TrigMatch','2Lep_SS'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep_TauVeto',	 cut = vardb.getCuts(['2Lep_JustNLep','TT','2Lep_pT','2Lep_ElEtaCut','2Lep_TrigMatch','2Lep_SS','TauVeto'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep_NJets',	 cut = vardb.getCuts(['2Lep_JustNLep','TT','2Lep_pT','2Lep_ElEtaCut','2Lep_TrigMatch','2Lep_SS','TauVeto','2Lep_NJet_SR'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep_NBJets',	 cut = vardb.getCuts(['2Lep_JustNLep','TT','2Lep_pT','2Lep_ElEtaCut','2Lep_TrigMatch','2Lep_SS','TauVeto','2Lep_NJet_SR','2Lep_NBJet'])  ) )

    # 2lepSS + 1tau
    #
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep1Tau_NLep',	     cut = vardb.getCuts(['2Lep1Tau_NLep'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep1Tau_TightLeptons', cut = vardb.getCuts(['2Lep1Tau_NLep','2Lep1Tau_TightLeptons'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep1Tau_pT',	     cut = vardb.getCuts(['2Lep1Tau_NLep','2Lep1Tau_TightLeptons','2Lep1Tau_pT'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep1Tau_TrigMatch',    cut = vardb.getCuts(['2Lep1Tau_NLep','2Lep1Tau_TightLeptons','2Lep1Tau_pT','2Lep1Tau_TrigMatch'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep1Tau_SS',	     cut = vardb.getCuts(['2Lep1Tau_NLep','2Lep1Tau_TightLeptons','2Lep1Tau_pT','2Lep1Tau_TrigMatch','2Lep1Tau_SS'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep1Tau_1Tau',	     cut = vardb.getCuts(['2Lep1Tau_NLep','2Lep1Tau_TightLeptons','2Lep1Tau_pT','2Lep1Tau_TrigMatch','2Lep1Tau_SS','2Lep1Tau_1Tau'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep1Tau_Zsidescut',    cut = vardb.getCuts(['2Lep1Tau_NLep','2Lep1Tau_TightLeptons','2Lep1Tau_pT','2Lep1Tau_TrigMatch','2Lep1Tau_SS','2Lep1Tau_1Tau','2Lep1Tau_Zsidescut'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep1Tau_NJets',	     cut = vardb.getCuts(['2Lep1Tau_NLep','2Lep1Tau_TightLeptons','2Lep1Tau_pT','2Lep1Tau_TrigMatch','2Lep1Tau_SS','2Lep1Tau_1Tau','2Lep1Tau_Zsidescut','2Lep1Tau_NJet_SR'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_2Lep1Tau_NBJets',       cut = vardb.getCuts(['2Lep1Tau_NLep','2Lep1Tau_TightLeptons','2Lep1Tau_pT','2Lep1Tau_TrigMatch','2Lep1Tau_SS','2Lep1Tau_1Tau','2Lep1Tau_Zsidescut','2Lep1Tau_NJet_SR','2Lep1Tau_NBJet'])  ) )

    # 3 lep
    #
    #vardb.registerCategory( MyCategory('CFChallenge_3Lep_NLep',	 cut = vardb.getCuts(['3Lep_JustNLep'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_3Lep_Charge',	 cut = vardb.getCuts(['3Lep_JustNLep','3Lep_Charge'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_3Lep_TightLeptons', cut = vardb.getCuts(['3Lep_JustNLep','3Lep_Charge','3Lep_TightLeptons'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_3Lep_pT',  	 cut = vardb.getCuts(['3Lep_JustNLep','3Lep_Charge','3Lep_TightLeptons','3Lep_pT'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_3Lep_TrigMatch',	 cut = vardb.getCuts(['3Lep_JustNLep','3Lep_Charge','3Lep_TightLeptons','3Lep_pT','3Lep_TrigMatch'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_3Lep_ZVeto',	 cut = vardb.getCuts(['3Lep_JustNLep','3Lep_Charge','3Lep_TightLeptons','3Lep_pT','3Lep_TrigMatch','3Lep_ZVeto'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_3Lep_MinZCut',	 cut = vardb.getCuts(['3Lep_JustNLep','3Lep_Charge','3Lep_TightLeptons','3Lep_pT','3Lep_TrigMatch','3Lep_ZVeto','3Lep_MinZCut'])  ) )
    #vardb.registerCategory( MyCategory('CFChallenge_3Lep_NJets',	 cut = vardb.getCuts(['3Lep_JustNLep','3Lep_Charge','3Lep_TightLeptons','3Lep_pT','3Lep_TrigMatch','3Lep_ZVeto','3Lep_MinZCut','3Lep_NJets'])  ) )

# ----------------------------------------------
# CRs where r/f rates for MM method are measured
# ----------------------------------------------

# ---------------------------------------------------------------
#  electron(muon) REAL rate measurement region:
#
#  -) OS
#  -) elel(mumu) || OF (where the electron(muon) is the probe)
#
#  electron(muon) FAKE rate measurement region:
#  -) SS
#  -) elel(mumu) || OF (where the electron(muon) is the probe)
#
#   NB: for MC samples we want to plot only the prompt (i.e., the prompt ch-flip? needs reviewing!!) contribution
#      (i.e, what will be subtracted to measure FAKE rate)
#
# Side note:
#	in these events, there is by construction only one probe
#	lepton and one tag lepton.
#	The vectorial-component notation 'el_probe_*[0]' is kept only
#	for consistency with the other sections of the code.
#
# ---------------------------------------------------------------

if doMMRates or doMMClosureRates:

    # ---------------------------------------
    # Special plots for MM real/fake rate CRs
    # ---------------------------------------

    #vardb.registerVar( Variable(shortname = 'ElTagPt', latexname = 'p_{T}^{tag e} [GeV]', ntuplename = 'el_tag_pt[0]/1e3', bins = 40, minval = 10.0, maxval = 210.0,) )
    #vardb.registerVar( Variable(shortname = 'ElTagEta', latexname = '#eta^{tag e}', ntuplename = 'TMath::Abs( el_tag_eta[0] )',bins = 8, minval = 0.0,  maxval = 2.6, manualbins = [ 0.0 , 0.5 , 0.8 , 1.1 , 1.37 , 1.52 , 2.0 , 2.25 , 2.6]) )
    vardb.registerVar( Variable(shortname = 'ElProbePt', latexname = 'p_{T}^{probe e} [GeV]', ntuplename = 'el_probe_pt[0]/1e3', bins = 40, minval = 10.0, maxval = 210.0,) )
    vardb.registerVar( Variable(shortname = 'ElProbeEta', latexname = '#eta^{probe e}', ntuplename = 'TMath::Abs( el_probe_eta[0] )', bins = 8, minval = 0.0,  maxval = 2.6, manualbins = [ 0.0 , 0.5 , 0.8 , 1.1 , 1.37 , 1.52 , 2.0 , 2.25 , 2.6]) )
    #vardb.registerVar( Variable(shortname = 'ElProbeNJets', latexname = 'Jet multiplicity', ntuplename = 'njets', bins = 8, minval = 2, maxval = 10) )

    #vardb.registerVar( Variable(shortname = 'MuTagPt', latexname = 'p_{T}^{tag #mu} [GeV]', ntuplename = 'muon_tag_pt[0]/1e3', bins = 40, minval = 10.0, maxval = 210.0,) )
    #vardb.registerVar( Variable(shortname = 'MuTagEta', latexname = '#eta^{tag #mu}', ntuplename = 'TMath::Abs( muon_tag_eta[0] )', bins = 8,  minval = 0.0, maxval = 2.5, manualbins = [ 0.0 , 0.1 , 0.4 , 0.7, 1.0,  1.3 , 1.6 , 1.9, 2.2, 2.5 ]) )
    vardb.registerVar( Variable(shortname = 'MuProbePt', latexname = 'p_{T}^{probe #mu} [GeV]', ntuplename = 'muon_probe_pt[0]/1e3', bins = 40, minval = 10.0, maxval = 210.0) )
    vardb.registerVar( Variable(shortname = 'MuProbeEta', latexname = '#eta^{probe #mu}', ntuplename = 'TMath::Abs( muon_probe_eta[0] )', bins = 8, minval = 0.0, maxval = 2.5, manualbins = [ 0.0 , 0.1 , 0.4 , 0.7, 1.0,  1.3 , 1.6 , 1.9, 2.2, 2.5 ]) )
    #vardb.registerVar( Variable(shortname = 'MuProbeNJets', latexname = 'Jet multiplicity', ntuplename = 'njets', bins = 8, minval = 2, maxval = 10) )

    # -----------------------------------------------------------------------------------------------------------------
    # MC Prompt/ch-flip subtraction in Fake SS CR: make sure you only plot:
    #
    # ---> to account for ttV,VV, rare top: use '2Lep_PurePromptEvent' (NB: this does not take into account charge flips!)
    # ---> to account for charge flip background : just use the data-driven estimate
    #
    truth_sub_SS = vardb.getCut('2Lep_PurePromptEvent')
    truth_sub_OS = ( vardb.getCut('DummyCut') )

    # Only if looking at MC-based QMisID...
    #
    if args.useMCChFlip :
       truth_sub_SS = ( vardb.getCut('2Lep_PurePromptEvent') | vardb.getCut('2Lep_ChFlipEvent') )

    # Require at least one non-prompt or charge flip lepton if you want to see the MC bkg composition in the fake region
    #
    if args.MCCompRF:
       truth_sub_SS = ( vardb.getCut('2Lep_NonPromptEvent') | vardb.getCut('2Lep_ChFlipEvent') )

    # Use this when extracting FAKE rates from simulation, or if doing MMClosure (--> use MC events w/ at least 1 non-prompt, but veto charge flips)
    # Use this when extracting REAL rates from simulation, or if doing MMClosure (--> use MC events w/ only prompt leptons, and veto on charge flips)
    #
    if ( doMMClosureRates or args.ratesFromMC ) :

       truth_sub_OS = ( vardb.getCut('2Lep_PurePromptEvent') )
       truth_sub_SS = ( vardb.getCut('2Lep_NonPromptEvent') ) # --> if tag/probe assignment has been done w/ truth, then the probe will be automatically a !prompt and not charge flip (DEFAULT)

       if doMMClosureRates:
          print '*********************************\n Doing MMClosure : looking at ttbar only! \n*********************************'
       if args.ratesFromMC:
          print '*********************************\n Measuring rates from MC simulation! \n*********************************'

       if ( args.doChFlipRate ):
          print '*********************************\nMEASURING CHARGE FLIP RATE IN MC\n*********************************'
          truth_sub_SS = vardb.getCut('2Lep_ChFlipEvent')     # --> if tag/probe assignment has been done w/ truth, then the probe will be automatically a charge flip

    # electron/muon R/F region(s)
    #

    # combine OF + SF
    #
    if ( args.lepFlavComp == "INCLUSIVE" ):
        #"""
    	vardb.registerCategory( MyCategory('FakeCRMuL',    cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_SS','2Lep_NJet_CR','2Lep_MuRealFakeRateCR','2Lep_ElEtaCut','2Lep_MuProbeAntiTight','2Lep_Zmincut']) & truth_sub_SS ) ) )
    	vardb.registerCategory( MyCategory('FakeCRMuT',    cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_SS','2Lep_NJet_CR','2Lep_MuRealFakeRateCR','2Lep_ElEtaCut','2Lep_MuProbeTight','2Lep_Zmincut']) & truth_sub_SS ) ) )
    	vardb.registerCategory( MyCategory('RealCRMuL',    cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_OS','2Lep_NJet_CR','2Lep_MuRealFakeRateCR','2Lep_MuProbeAntiTight','2Lep_Zmincut']) & truth_sub_OS ) ) )
    	vardb.registerCategory( MyCategory('RealCRMuT',    cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_OS','2Lep_NJet_CR','2Lep_MuRealFakeRateCR','2Lep_MuProbeTight','2Lep_Zmincut']) & truth_sub_OS ) ) )

    	vardb.registerCategory( MyCategory('FakeCRElL',    cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_SS','2Lep_NJet_CR','2Lep_ElRealFakeRateCR','2Lep_ElEtaCut','2Lep_ElProbeAntiTight','2Lep_Zsidescut','2Lep_Zmincut']) & truth_sub_SS ) ) )
    	vardb.registerCategory( MyCategory('FakeCRElT',    cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_SS','2Lep_NJet_CR','2Lep_ElRealFakeRateCR','2Lep_ElEtaCut','2Lep_ElProbeTight','2Lep_Zsidescut','2Lep_Zmincut']) & truth_sub_SS ) ) )
    	vardb.registerCategory( MyCategory('RealCRElL',    cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_OS','2Lep_NJet_CR','2Lep_ElRealFakeRateCR','2Lep_ElProbeAntiTight','2Lep_Zmincut']) & truth_sub_OS ) ) )
    	vardb.registerCategory( MyCategory('RealCRElT',    cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_OS','2Lep_NJet_CR','2Lep_ElRealFakeRateCR','2Lep_ElProbeTight','2Lep_Zmincut']) & truth_sub_OS ) ) )
    	#"""
    	# temp: measure rates in inclusive jet multiplicity region
    	"""
    	vardb.registerCategory( MyCategory('FakeCRMuL',    cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_SS','2Lep_MuRealFakeRateCR','2Lep_ElEtaCut','2Lep_MuProbeAntiTight','2Lep_Zmincut']) & truth_sub_SS ) ) )
    	vardb.registerCategory( MyCategory('FakeCRMuT',    cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_SS','2Lep_MuRealFakeRateCR','2Lep_ElEtaCut','2Lep_MuProbeTight','2Lep_Zmincut']) & truth_sub_SS ) ) )
    	vardb.registerCategory( MyCategory('RealCRMuL',    cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_OS','2Lep_MuRealFakeRateCR','2Lep_MuProbeAntiTight','2Lep_Zmincut']) & truth_sub_OS ) ) )
    	vardb.registerCategory( MyCategory('RealCRMuT',    cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_OS','2Lep_MuRealFakeRateCR','2Lep_MuProbeTight','2Lep_Zmincut']) & truth_sub_OS ) ) )
    	#
    	vardb.registerCategory( MyCategory('FakeCRElL',    cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_SS','2Lep_ElRealFakeRateCR','2Lep_ElEtaCut','2Lep_ElProbeAntiTight','2Lep_Zsidescut','2Lep_Zmincut']) & truth_sub_SS ) ) )
    	vardb.registerCategory( MyCategory('FakeCRElT',    cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_SS','2Lep_ElRealFakeRateCR','2Lep_ElEtaCut','2Lep_ElProbeTight','2Lep_Zsidescut','2Lep_Zmincut']) & truth_sub_SS ) ) )
    	vardb.registerCategory( MyCategory('RealCRElL',    cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_OS','2Lep_ElRealFakeRateCR','2Lep_ElProbeAntiTight','2Lep_Zmincut']) & truth_sub_OS ) ) )
    	vardb.registerCategory( MyCategory('RealCRElT',    cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_OS','2Lep_ElRealFakeRateCR','2Lep_ElProbeTight','2Lep_Zmincut']) & truth_sub_OS ) ) )
    	"""
    # SF only
    #
    if ( args.lepFlavComp == "SF" ):

    	vardb.registerCategory( MyCategory('MuMuFakeCRMuL',  cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_SS','2Lep_NJet_CR','2Lep_MuRealFakeRateCR','2Lep_MuProbeAntiTight','2Lep_MuMu_Event','2Lep_Zmincut']) & truth_sub_SS ) ) )
    	vardb.registerCategory( MyCategory('MuMuFakeCRMuT',  cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_SS','2Lep_NJet_CR','2Lep_MuRealFakeRateCR','2Lep_MuProbeTight','2Lep_MuMu_Event','2Lep_Zmincut']) & truth_sub_SS ) ) )
    	vardb.registerCategory( MyCategory('MuMuRealCRMuL',  cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_OS','2Lep_NJet_CR','2Lep_MuRealFakeRateCR','2Lep_MuProbeAntiTight','2Lep_MuMu_Event','2Lep_Zmincut']) & truth_sub_OS ) ) )
    	vardb.registerCategory( MyCategory('MuMuRealCRMuT',  cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_OS','2Lep_NJet_CR','2Lep_MuRealFakeRateCR','2Lep_MuProbeTight','2Lep_MuMu_Event','2Lep_Zmincut']) & truth_sub_OS ) ) )

    	vardb.registerCategory( MyCategory('ElElFakeCRElL',  cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_SS','2Lep_NJet_CR','2Lep_ElRealFakeRateCR','2Lep_ElEtaCut','2Lep_ElProbeAntiTight', '2Lep_ElEl_Event','2Lep_Zsidescut','2Lep_Zmincut']) & truth_sub_SS ) ) )
    	vardb.registerCategory( MyCategory('ElElFakeCRElT',  cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_SS','2Lep_NJet_CR','2Lep_ElRealFakeRateCR','2Lep_ElEtaCut','2Lep_ElProbeTight','2Lep_ElEl_Event','2Lep_Zsidescut','2Lep_Zmincut']) & truth_sub_SS ) ) )
    	vardb.registerCategory( MyCategory('ElElRealCRElL',  cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_OS','2Lep_NJet_CR','2Lep_ElRealFakeRateCR','2Lep_ElProbeAntiTight', '2Lep_ElEl_Event','2Lep_Zmincut']) & truth_sub_OS ) ) )
    	vardb.registerCategory( MyCategory('ElElRealCRElT',  cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_OS','2Lep_NJet_CR','2Lep_ElRealFakeRateCR','2Lep_ElProbeTight','2Lep_ElEl_Event','2Lep_Zmincut']) & truth_sub_OS ) ) )

    # OF only
    #
    if ( args.lepFlavComp == "OF" ):

        vardb.registerCategory( MyCategory('OFFakeCRMuL',   cut =   vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_SS','2Lep_NJet_CR','2Lep_MuRealFakeRateCR','2Lep_ElEtaCut','2Lep_MuProbeAntiTight','2Lep_OF_Event']) & truth_sub_SS ) )
        vardb.registerCategory( MyCategory('OFFakeCRMuT',   cut =   vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_SS','2Lep_NJet_CR','2Lep_MuRealFakeRateCR','2Lep_ElEtaCut','2Lep_MuProbeTight','2Lep_OF_Event']) & truth_sub_SS ) )
        vardb.registerCategory( MyCategory('OFRealCRMuL',   cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_OS','2Lep_NJet_CR','2Lep_MuRealFakeRateCR','2Lep_MuProbeAntiTight','2Lep_OF_Event']) & truth_sub_OS ) ) )
        vardb.registerCategory( MyCategory('OFRealCRMuT',   cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_OS','2Lep_NJet_CR','2Lep_MuRealFakeRateCR','2Lep_MuProbeTight','2Lep_OF_Event']) & truth_sub_OS ) ) )

        vardb.registerCategory( MyCategory('OFFakeCRElL',   cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_SS','2Lep_NJet_CR','2Lep_ElRealFakeRateCR','2Lep_ElEtaCut','2Lep_ElProbeAntiTight','2Lep_OF_Event']) & truth_sub_SS ) ) )
        vardb.registerCategory( MyCategory('OFFakeCRElT',   cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_SS','2Lep_NJet_CR','2Lep_ElRealFakeRateCR','2Lep_ElEtaCut','2Lep_ElProbeTight','2Lep_OF_Event']) & truth_sub_SS ) ) )
        vardb.registerCategory( MyCategory('OFRealCRElL',   cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_OS','2Lep_NJet_CR','2Lep_ElRealFakeRateCR','2Lep_ElProbeAntiTight','2Lep_OF_Event']) & truth_sub_OS ) ) )
        vardb.registerCategory( MyCategory('OFRealCRElT',   cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_LepTagTightTrigMatched','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_OS','2Lep_NJet_CR','2Lep_ElRealFakeRateCR','2Lep_ElProbeTight','2Lep_OF_Event']) & truth_sub_OS ) ) )


# NB: do NOT apply truth cuts here! They will be implemented in the appropriate background classes...
#
if doMMClosureTest:
    print ''

    if ( args.fakeMethod == 'MM' or args.fakeMethod == 'FF' ):
        #
        # MuMu region
        #
        vardb.registerCategory( MyCategory('MuMuSS_SR_HighJet_DataDriven_Closure',   cut = vardb.getCuts(['2Lep_TrigMatch', 'TrigDec', '2Lep_NBJet', '2Lep_NLep', 'TauVeto', '2Lep_SS', '2Lep_MuMu_Event',  '2Lep_NJet_SR']) ) )
        #vardb.registerCategory( MyCategory('MuMuSS_SR_LowJet_DataDriven_Closure',    cut = vardb.getCuts(['2Lep_TrigMatch', 'TrigDec', '2Lep_NBJet', '2Lep_NLep', 'TauVeto', '2Lep_SS', '2Lep_MuMu_Event',  '2Lep_NJet_CR']) ) )
        #vardb.registerCategory( MyCategory('MuMuSS_SR_AllJet_DataDriven_Closure',    cut = vardb.getCuts(['2Lep_TrigMatch', 'TrigDec', '2Lep_NBJet', '2Lep_NLep', 'TauVeto', '2Lep_SS', '2Lep_MuMu_Event']) ) )
        #
	# OF region
	#
	vardb.registerCategory( MyCategory('OFSS_SR_HighJet_DataDriven_Closure',     cut = vardb.getCuts(['2Lep_TrigMatch', 'TrigDec', '2Lep_NBJet', '2Lep_NLep', 'TauVeto', '2Lep_SS', '2Lep_OF_Event', '2Lep_ElEtaCut',  '2Lep_NJet_SR']) ) )
	#vardb.registerCategory( MyCategory('OFSS_SR_LowJet_DataDriven_Closure',      cut = vardb.getCuts(['2Lep_TrigMatch', 'TrigDec', '2Lep_NBJet', '2Lep_NLep', 'TauVeto', '2Lep_SS', '2Lep_OF_Event', '2Lep_ElEtaCut',  '2Lep_NJet_CR']) ) )
	#vardb.registerCategory( MyCategory('OFSS_SR_AllJet_DataDriven_Closure',      cut = vardb.getCuts(['2Lep_TrigMatch', 'TrigDec', '2Lep_NBJet', '2Lep_NLep', 'TauVeto', '2Lep_SS', '2Lep_OF_Event', '2Lep_ElEtaCut']) ) )
	#
	# ElEl region
	#
	vardb.registerCategory( MyCategory('ElElSS_SR_HighJet_DataDriven_Closure',   cut = vardb.getCuts(['2Lep_TrigMatch', 'TrigDec', '2Lep_NBJet', '2Lep_NLep', 'TauVeto', '2Lep_SS', '2Lep_ElEl_Event', '2Lep_ElEtaCut', '2Lep_NJet_SR']) ) )
	#vardb.registerCategory( MyCategory('ElElSS_SR_LowJet_DataDriven_Closure',    cut = vardb.getCuts(['2Lep_TrigMatch', 'TrigDec', '2Lep_NBJet', '2Lep_NLep', 'TauVeto', '2Lep_SS', '2Lep_ElEl_Event', '2Lep_ElEtaCut', '2Lep_NJet_CR']) ) )
        #vardb.registerCategory( MyCategory('ElElSS_SR_AllJet_DataDriven_Closure',    cut = vardb.getCuts(['2Lep_TrigMatch', 'TrigDec', '2Lep_NBJet', '2Lep_NLep', 'TauVeto', '2Lep_SS', '2Lep_ElEl_Event', '2Lep_ElEtaCut']) ) )

    elif ( args.fakeMethod == 'ABCD' ):
        # NB: the closure test for ABCD is meaningful only in SR, b/c the [TT,low nr. jet] CR is already used to derive the theta factor
        #
        # MuMu region
        #
        vardb.registerCategory( MyCategory('MuMuSS_SR_HighJet_DataDriven_Closure',   cut = vardb.getCuts(['2Lep_TrigMatch', 'TrigDec', '2Lep_NBJet', '2Lep_NLep', 'TauVeto', '2Lep_SS', '2Lep_MuMu_Event',  '2Lep_NJet_SR']) ) )
        #
	# OF region
	#
	vardb.registerCategory( MyCategory('OFSS_SR_HighJet_DataDriven_Closure',     cut = vardb.getCuts(['2Lep_TrigMatch', 'TrigDec', '2Lep_NBJet', '2Lep_NLep', 'TauVeto', '2Lep_SS', '2Lep_OF_Event', '2Lep_ElEtaCut',  '2Lep_NJet_SR']) ) )
	#
	# ElEl region
	#
	vardb.registerCategory( MyCategory('ElElSS_SR_HighJet_DataDriven_Closure',   cut = vardb.getCuts(['2Lep_TrigMatch', 'TrigDec', '2Lep_NBJet', '2Lep_NLep', 'TauVeto', '2Lep_SS', '2Lep_ElEl_Event', '2Lep_ElEtaCut', '2Lep_NJet_SR']) ) )

if doMMRatesLHFit:

   truth_sub_SS = vardb.getCut('2Lep_PurePromptEvent')
   truth_sub_OS = ( vardb.getCut('DummyCut') )

   vardb.registerCategory( MyCategory('OS_ElEl_TT',    cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_OS','2Lep_NJet_CR','2Lep_ElEl_Event','TT']) & truth_sub_OS ) ) )
   vardb.registerCategory( MyCategory('OS_ElEl_TL',    cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_OS','2Lep_NJet_CR','2Lep_ElEl_Event','TL_LT']) & truth_sub_OS ) ) )
   vardb.registerCategory( MyCategory('OS_ElEl_LL',    cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_OS','2Lep_NJet_CR','2Lep_ElEl_Event','LL']) & truth_sub_OS ) ) )

   vardb.registerCategory( MyCategory('OS_MuMu_TT',    cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_OS','2Lep_NJet_CR','2Lep_MuMu_Event','TT']) & truth_sub_OS ) ) )
   vardb.registerCategory( MyCategory('OS_MuMu_TL',    cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_OS','2Lep_NJet_CR','2Lep_MuMu_Event','TL_LT']) & truth_sub_OS ) ) )
   vardb.registerCategory( MyCategory('OS_MuMu_LL',    cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_OS','2Lep_NJet_CR','2Lep_MuMu_Event','LL']) & truth_sub_OS ) ) )

   vardb.registerCategory( MyCategory('OS_OF_TT',      cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_OS','2Lep_NJet_CR','2Lep_OF_Event','TT']) & truth_sub_OS ) ) )
   vardb.registerCategory( MyCategory('OS_OF_TL',      cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_OS','2Lep_NJet_CR','2Lep_OF_Event','TL_LT']) & truth_sub_OS ) ) )
   vardb.registerCategory( MyCategory('OS_OF_LL',      cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_OS','2Lep_NJet_CR','2Lep_OF_Event','LL']) & truth_sub_OS ) ) )

   vardb.registerCategory( MyCategory('SS_ElEl_TT',    cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_SS','2Lep_NJet_CR','2Lep_ElEl_Event','TT','2Lep_Zsidescut','2Lep_ElEtaCut']) & truth_sub_SS ) ) )
   vardb.registerCategory( MyCategory('SS_ElEl_TL',    cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_SS','2Lep_NJet_CR','2Lep_ElEl_Event','TL_LT','2Lep_Zsidescut','2Lep_ElEtaCut']) & truth_sub_SS ) ) )
   vardb.registerCategory( MyCategory('SS_ElEl_LL',    cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_SS','2Lep_NJet_CR','2Lep_ElEl_Event','LL','2Lep_Zsidescut','2Lep_ElEtaCut']) & truth_sub_SS ) ) )

   vardb.registerCategory( MyCategory('SS_MuMu_TT',    cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_SS','2Lep_NJet_CR','2Lep_MuMu_Event','TT']) & truth_sub_SS ) ) )
   vardb.registerCategory( MyCategory('SS_MuMu_TL',    cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_SS','2Lep_NJet_CR','2Lep_MuMu_Event','TL_LT']) & truth_sub_SS ) ) )
   vardb.registerCategory( MyCategory('SS_MuMu_LL',    cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_SS','2Lep_NJet_CR','2Lep_MuMu_Event','LL']) & truth_sub_SS ) ) )

   vardb.registerCategory( MyCategory('SS_OF_TT',      cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_SS','2Lep_NJet_CR','2Lep_OF_Event','TT','2Lep_ElEtaCut']) & truth_sub_SS ) ) )
   vardb.registerCategory( MyCategory('SS_OF_TL',      cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_SS','2Lep_NJet_CR','2Lep_OF_Event','TL_LT','2Lep_ElEtaCut']) & truth_sub_SS ) ) )
   vardb.registerCategory( MyCategory('SS_OF_LL',      cut = ( vardb.getCuts(['TrigDec','2Lep_TrigMatch','2Lep_NBJet','2Lep_NLep_MMRates','TauVeto','2Lep_SS','2Lep_NJet_CR','2Lep_OF_Event','LL','2Lep_ElEtaCut']) & truth_sub_SS ) ) )

# ------------------------------------------------------------
# TTHBackgrounds2015 is the class used to manage each process:
#
#   Pass the input informations and the definitions and it
#   will perform the background estimationTTBarClosu
# ------------------------------------------------------------

ttH2015 = TTHBackgrounds2015(inputs, vardb)

# ------------------------------------
# Set the integrated luminosity (fb-1)
# ------------------------------------

# period D1-J6

#ttH2015.luminosity = 3.302 # GRL v71 - EOYE GRL
ttH2015.luminosity = 3.209 # GRL v73 - Moriond GRL
ttH2015.lumi_units = 'fb-1'

# for MM closure
if doMMClosureTest or doMMClosureRates:
	ttH2015.luminosity = 3.209
	ttH2015.lumi_units = 'fb-1'

# --------------------
# set the event weight
# --------------------

# MC generator event weight
#
weight_generator = 'mcEventWeight'

# PRW weight
#
weight_pileup = '1.0'
#if not ( args.noWeights or doMMClosureTest or doMMClosureRates ):
if not ( args.noWeights ):
    weight_pileup = 'weight_pileup'

weight_glob = str(weight_generator) + ' * ' + str(weight_pileup)

print ("Global event weight (apply to ALL categories to MC only) --> {0}\n".format( weight_glob ) )

ttH2015.eventweight = weight_glob

# ------------------------------------

ttH2015.useZCorrections = False

if doTwoLepSR or doTwoLepLowNJetCR or dottWCR or doMMClosureTest:
    ttH2015.channel = 'TwoLepSS'
elif doThreeLepSR or doThreeLepLowNJetCR or dottZCR or doWZonCR or doWZoffCR or doWZHFonCR or doWZHFoffCR:
    ttH2015.channel = 'ThreeLep'
elif doFourLepSR:
    ttH2015.channel = 'FourLep'
elif doDataMCCR or doZSSpeakCR or doMMRates or doMMRatesLHFit or doMMClosureRates or doCFChallenge:
    ttH2015.channel = 'TwoLepCR'

cut = None
systematics = None
systematicsdirection = None # 'UP', 'DOWN'

events = {}
hists  = {}
# --------------------------------------
# Dictionary with systematics histograms
# --------------------------------------
systs = {}

# ----------------------------------
# List of the backgrounds considered
# ----------------------------------

samplenames = { 'Observed':'observed',
                'TTBarH':'signal',
		'TTBarW':'ttbarwbkg',
		'TTBarZ':'ttbarzbkg',
		'Top':'topbkg',
		'TTBar':'ttbarbkg',
		'TTBarClosure':'ttbarbkg',
		'TopCF':'topcfbkg',
		'Diboson':'dibosonbkg',
		'DibosonCF':'dibosoncfbkg',
		'HtoZZ':'htozzbkg',
		'Zjets':'zjetsbkg',
		'Zeejets':'zeejetsbkg',
		'Zmumujets':'zmumujetsbkg',
		'Ztautaujets':'ztautaujetsbkg',
		'ZjetsHF':'zjetsbkg',
		'ZjetsLF':'zjetsbkg',
		'Wjets':'wjetsbkg',
		'Wenujets':'wenujets',
		'Wmunujets':'wmunujets',
		'Wtaunujets':'wtaunujets',
		'Prompt':'promptbkg',
		'ChargeFlip':'chargeflipbkg',
		'ChargeFlipInclusiveFlavRates':'chargeflipbkg',
		'ChargeFlipMC':'chargeflipbkg',
		'FakesMC':'fakesbkg',
		'FakesFF':'fakesbkg',
		'FakesMM':'fakesbkg',
		'FakesABCD':'fakesbkg',
                'FakesClosureMM':'fakesbgk',
                'FakesClosureABCD':'fakesbgk',
                'FakesClosureDataABCD':'fakesbgk',
	      }
#
# Override colours!
#
colours      = {'Observed':kBlack,
        	'TTBarH':kRed, # kBlack,
        	'TTBarW':kYellow-9, # kRed-4,
        	'TTBarZ':kAzure+1,  # kRed-7,
        	'Top':kGray, #kBlue+1,
        	'TTBar': kRed - 4,# kAzure+8,
		'TTBarClosure':kAzure+8,
        	'TopCF':kAzure-4,
        	'Diboson':kGreen-9, # kYellow-9,
        	'DibosonCF':kOrange-3,
        	'HtoZZ':kTeal+9,
        	'Zjets': kCyan -9, # kGreen,
		'Zeejets':kGreen-7,
		'Zmumujets':kTeal+2,
		'Ztautaujets':kTeal,
        	'ZjetsHF':kGreen+2,
        	'ZjetsLF':kGreen,
        	'Wjets':kWhite,
		'Wenujets':kGray,
		'Wmunujets':kGray+1,
		'Wtaunujets':kGray+2,
        	'Prompt':kOrange-3,
        	'ChargeFlip':kMagenta+1,# kAzure-4,
        	'ChargeFlipInclusiveFlavRates':kMagenta+1,# kAzure-4,
        	'ChargeFlipMC':kMagenta+1,#kAzure-4,
		'FakesMC':kMagenta-9,
        	'FakesMM':kMagenta-9,
        	'FakesFF':kMagenta-9, #kAzure-9,
        	'FakesMM':kMagenta-9, #kTeal-9,
		'FakesABCD':kMagenta-9, # kCyan-9,
                'FakesClosureMM':kTeal+1,
		'FakesClosureABCD':kCyan-9,
		'FakesClosureDataABCD': kMagenta-9, # kCyan-9,
              }

if ( doSR or doLowNJetCR ):

    ttH2015.signals     = ['TTBarH']
    ttH2015.observed    = ['Observed']

    if not doFourLepSR:

	if doMM:
            # ---> all the MC backgrounds use a truth req. of only prompt leptons in the event (and ch-flip veto) to avoid double counting with
            #      data-driven charge flip and fakes estimate

	    plotbackgrounds	= ['TTBarW','TTBarZ','Diboson','Top','FakesMM']
    	    ttH2015.backgrounds = ['TTBarW','TTBarZ','Diboson','Top','FakesMM']

	    if args.useMCChFlip:
	      plotbackgrounds.append('ChargeFlipMC')
	      ttH2015.backgrounds.append('ChargeFlipMC')
	    else:
	      plotbackgrounds.append('ChargeFlip')
	      ttH2015.backgrounds.append('ChargeFlip')

	elif doFF:

    	    plotbackgrounds	    = ['TTBarW','TTBarZ','Diboson','Top','FakesFF']
    	    ttH2015.backgrounds     = ['TTBarW','TTBarZ','Diboson','Top','FakesFF']
	    ttH2015.sub_backgrounds = ['TTBarW','TTBarZ','Diboson','Top']

	    if args.useMCChFlip:
	      plotbackgrounds.append('ChargeFlipMC')
	      ttH2015.backgrounds.append('ChargeFlipMC')
	      ttH2015.sub_backgrounds.append('ChargeFlipMC')
	    else:
	      plotbackgrounds.append('ChargeFlip')
	      ttH2015.backgrounds.append('ChargeFlip')
	      ttH2015.sub_backgrounds.append('ChargeFlip')

	elif doABCD:

      	    plotbackgrounds         = ['TTBarW','TTBarZ','Diboson','Top','FakesABCD']
    	    ttH2015.backgrounds     = ['TTBarW','TTBarZ','Diboson','Top','FakesABCD']
	    ttH2015.sub_backgrounds = ['TTBarW','TTBarZ','Diboson','Top']

            if doTwoLepLowNJetCR :
	        # Closure in OF, low njet w/ ttbar MC rewweighted by theta factors
    	        plotbackgrounds	    = ['TTBarW','TTBarZ','Diboson','Top','FakesClosureDataABCD']
    	        ttH2015.backgrounds = ['TTBarW','TTBarZ','Diboson','Top','FakesClosureDataABCD']

	    if args.useMCChFlip:
	      plotbackgrounds.append('ChargeFlipMC')
	      ttH2015.backgrounds.append('ChargeFlipMC')
	      ttH2015.sub_backgrounds.append('ChargeFlipMC')
	    else:
	      plotbackgrounds.append('ChargeFlip')
	      ttH2015.backgrounds.append('ChargeFlip')
	      ttH2015.sub_backgrounds.append('ChargeFlip')

	else:
    	    # MC based estimate of fakes (and charge flips) - make sure any truth cut is removed!!
    	    plotbackgrounds	= ['TTBarW','TTBarZ','Diboson','TTBar','Top','Zjets']
    	    ttH2015.backgrounds = ['TTBarW','TTBarZ','Diboson','TTBar','Top','Zjets']
    else:
        # no fakes in 4lep
        plotbackgrounds	    = ['TTBarW','TTBarZ','Diboson','TTBar','Top','Zjets']
        ttH2015.backgrounds = ['TTBarW','TTBarZ','Diboson','TTBar','Top','Zjets']

if doMMRates:

    ttH2015.signals     = ['TTBarH']
    ttH2015.observed    = ['Observed']
    if args.ratesFromMC:
        ttH2015.observed    = []

    if args.useMCChFlip:
       plotbackgrounds     = ['TTBar','Zjets','TTBarW','TTBarZ','Diboson','Wjets','Top']
       ttH2015.backgrounds = ['TTBar','Zjets','TTBarW','TTBarZ','Diboson','Wjets','Top']
    else:
       plotbackgrounds	    = ['TTBar','Zjets','TTBarW','TTBarZ','Diboson','Wjets','ChargeFlipInclusiveFlavRates','Top']
       ttH2015.backgrounds  = ['TTBar','Zjets','TTBarW','TTBarZ','Diboson','Wjets','ChargeFlipInclusiveFlavRates','Top']

if doMMRatesLHFit:

    ttH2015.signals     = ['TTBarH']
    ttH2015.observed    = ['Observed']
    plotbackgrounds	= ['TTBar','Zjets','TTBarW','TTBarZ','Diboson','Wjets','Top']
    ttH2015.backgrounds = ['TTBar','Zjets','TTBarW','TTBarZ','Diboson','Wjets','Top']

    if args.useMCChFlip:
      plotbackgrounds.append('ChargeFlipMC')
      ttH2015.backgrounds.append('ChargeFlipMC')
    else:
      plotbackgrounds.append('ChargeFlip')
      ttH2015.backgrounds.append('ChargeFlip')

if doDataMCCR:

    ttH2015.signals     = ['TTBarH']
    ttH2015.observed    = ['Observed']
    plotbackgrounds     = ['Zeejets','Zmumujets','Ztautaujets','TTBar','Diboson','TTBarW','TTBarZ','Top','Wjets']
    ttH2015.backgrounds = ['Zeejets','Zmumujets','Ztautaujets','TTBar','Diboson','TTBarW','TTBarZ','Top','Wjets']

if doZSSpeakCR:

    ttH2015.signals     = ['TTBarH']
    ttH2015.observed    = ['Observed']
    plotbackgrounds     = ['ChargeFlipMC','FakesMC','Prompt']
    ttH2015.backgrounds = ['ChargeFlipMC','FakesMC','Prompt']

if doCFChallenge:

    ttH2015.signals     = []
    ttH2015.observed    = ['Observed']
    plotbackgrounds     = ['TTBar'] #['TTBarH']
    ttH2015.backgrounds = ['TTBar'] #['TTBarH']

if doMMClosureRates:

      ttH2015.signals	  = []
      ttH2015.observed    = []
      plotbackgrounds	  = ['TTBar']
      ttH2015.backgrounds = ['TTBar']

if doMMClosureTest:

    if doMM:
        ttH2015.signals     = ['FakesClosureABCD']
        ttH2015.observed    = ['TTBarClosure'] # truth cuts done internally in TTBarClosure class
        plotbackgrounds	    = ['FakesClosureMM']
        ttH2015.backgrounds = ['FakesClosureMM'] # truth cuts done internally in FakesClosureMM class
    elif doFF:
        ttH2015.signals     = ['FakesClosureABCD']
        ttH2015.observed    = ['TTBar']
        plotbackgrounds	    = ['FakesFF']
        ttH2015.backgrounds = ['FakesFF']
    elif doABCD:
        ttH2015.signals     = []
        ttH2015.observed    = ['TTBarClosure'] # truth cuts done internally in TTBarClosure class
        #ttH2015.observed    = ['FakesClosureMM'] # truth cuts done internally in TTBarClosure class
        plotbackgrounds	    = ['FakesClosureABCD']
        ttH2015.backgrounds = ['FakesClosureABCD'] # truth cuts done internally in FakesClosureABCD class
    else:
        ttH2015.signals	    = []
        ttH2015.observed    = []
        #ttH2015.observed    = ['TTBar']
        plotbackgrounds	    = ['TTBar']
        ttH2015.backgrounds = ['TTBar']

if args.noSignal:
    ttH2015.signals = []

doShowRatio = True
# Make blinded plots in SR unless configured from input
#
if doSR and not args.doUnblinding:
    ttH2015.observed = []
    doShowRatio = False

# -------------------------------------------------------
# Filling histname with the name of the variables we want
#
# Override colours as well
# -------------------------------------------------------
histname   = {'Expected':'expected'}
histcolour = {'Expected':kBlack}
for samp in ttH2015.backgrounds:
    histname[samp]  = samplenames[samp]
    histcolour[samp] = colours[samp]
    #
    # Will override default colour based on the dictionary provided above
    #
    ttH2015.str_to_class(samp).colour = colours[samp]
for samp in ttH2015.observed:
    histname[samp]  = samplenames[samp]
    histcolour[samp] = colours[samp]
    #
    # Will override default colour based on the dictionary provided above
    #
    ttH2015.str_to_class(samp).colour = colours[samp]
for samp in ttH2015.signals:
    histname[samp]  = samplenames[samp]
    histcolour[samp] = colours[samp]
    #
    # Will override default colour based on the dictionary provided above
    #
    ttH2015.str_to_class(samp).colour = colours[samp]

print histname
print histcolour

# ---------------------------------
# Processing categories in sequence
# ---------------------------------
for category in vardb.categorylist:

    print ("\n*********************************************\n\nMaking plots in category: {0}\n".format( category.name ))
    if ( category.cut != None ):
        print ("\tdefined by cuts --> {0}\n".format( category.cut.cutname ))

    signalfactor = 1.0
    background = ttH2015

    # MMClosureTest: do not look at ABCD fakes, unless in SR
    #
    if doMMClosureTest:
       if doMM or doFF:
          if not ( "2Lep_NJet_SR" in category.cut.cutname ):
	     ttH2015.signals = []
	  elif ( ( "2Lep_NJet_SR" in category.cut.cutname ) and not ttH2015.signals ):
             ttH2015.signals = ['FakesClosureABCD']

    # NB: *must* initialise this to 1.0 !!
    #
    lepSF_weight  = '1.0'

    # ---> apply the lepton SFs to the event here!
    #
    #if not ( args.noWeights or doMMClosureTest or doMMClosureRates ):
    if not ( args.noWeights ):

        lepSF_weight = 'weight_lepton_trig_HTop[0] * weight_lepton_reco_HTop[0] * weight_lepton_iso_HTop[0] * weight_lepton_ID_HTop[0] * weight_lepton_TTVA_HTop[0]'

    print ("\tApplying lepton SFs (to MC only) --> {0}\n".format( lepSF_weight ))

    # ------------------------------
    # Processing different variables
    # ------------------------------
    for idx,var in enumerate(vardb.varlist, start=0):

        # NB: *must* initialise this to 1.0 !!
        #
        bjetSF_weight = '1.0'
        combined_SF_weight = '1.0'

        print ("\t\tNow plotting variable:\t{0}\n".format(var.shortname))

        # When looking at jet multiplicity distributions w/ bjets, BTag SF must be applied also to categories w/o any bjet cut
        #
        #if not ( args.noWeights or doMMClosureTest or doMMClosureRates ):
        if not ( args.noWeights ):
            if  ( ( ("BJet") in category.cut.cutname and not doRelaxedBJetCut ) or  ( ("BJetSR") in category.cut.cutname ) ) or ("BJet") in var.shortname:

                bjetSF_weight = 'weight_jet__MV2c20_SFFix77[0]'

                print ("\t\tCategory contains a cut on BJet multiplicity, or plotting variable \'Bjet\' : apply BTagging SF (to MC only) --> {0}\n".format( bjetSF_weight ))

        combined_SF_weight = str(lepSF_weight) + ' * ' + str(bjetSF_weight)

        print ("\t\t----------------------\n\t\tCombined SF weight --> {0}\n\t\t----------------------\n".format( combined_SF_weight ) )

	# Get event yields for *this* category. Do it only for the first variable in the list
        #
        if ( args.printEventYields and idx is 0 ):

            events[category.name] = background.events(cut=cut, eventweight=combined_SF_weight, category=category, hmass=['125'], systematics=systematics, systematicsdirection=systematicsdirection)

	#"""

        # --------------------------
        # Avoid making useless plots
        # --------------------------

        if ( ("MuMu") in category.name and ("El") in var.shortname ) or ( ("ElEl") in category.name and ("Mu") in var.shortname ):
            print ("\tSkipping variable: {0}\n".format( var.shortname ))
            continue
        if ( ( ("MuEl") in category.name or ("ElMu") in category.name or ("OF") in category.name ) and ( ("El1") in var.shortname or ("Mu1") in var.shortname ) ) :
            print ("\tSkipping variable: {0}\n".format( var.shortname ))
            continue

        if doMMRates:
            # if probe is a muon, do not plot ElProbe* stuff!
            if ( ( ("MuRealFakeRateCR") in category.cut.cutname ) and ( ("ElProbe") in var.shortname ) ):
                print ("\tSkipping variable: {0}\n".format( var.shortname ))
                continue
            # if probe is an electron, do not plot MuProbe* stuff!
            if ( ( ("ElRealFakeRateCR") in category.cut.cutname ) and ( ("MuProbe") in var.shortname ) ):
                print ("\tSkipping variable: {0}\n".format( var.shortname ))
                continue
            # be smart when looking at OF regions!
            if ( ("OF") in category.name and( ("MuRealFakeRateCR") in category.cut.cutname ) and ( ("MuTag") in var.shortname ) ):
                print ("\tSkipping variable: {0}\n".format( var.shortname ))
                continue
            if ( ("OF") in category.name and( ("ElRealFakeRateCR") in category.cut.cutname ) and ( ("ElTag") in var.shortname ) ):
                print ("\tSkipping variable: {0}\n".format( var.shortname ))
                continue

        # ---------------------------------------------------------
        # Creating a directory for the category if it doesn't exist
        # ---------------------------------------------------------
        fakeestimate=''
        if doMM:
            fakeestimate='_MM'
        if doFF:
            fakeestimate='_FF'
        if doABCD:
            fakeestimate='_ABCD'

        dirname =  'OutputPlots' + fakeestimate + '_' + args.outdirname + '/'

        # If specified a cut before entering the loop, it will be applied
        # Otherwise, all the cuts registered above for *this* category will be applied
        #
        if cut:
            dirname = dirname  + category.name + ' ' + cut.cutname
        else:
            dirname = dirname  + category.name
        if args.doLogScaleX:
            var.logaxisX = True # activate X-axis log scale in plot
            dirname += '_LOGX'
        if args.doLogScaleY:
            var.logaxis  = True # activate Y-axis log scale in plot
            dirname += '_LOGY'
        dirname = dirname.replace(' ', '_')

        try:
            os.makedirs(dirname)
        except:
            pass

        # -----------------------------------------------
        # Making a plot with ( category + variable ) name
        # -----------------------------------------------
        plotname = dirname + '/' + category.name + ' ' + var.shortname
        plotname = plotname.replace(' ', '_')

        if ( args.debug ):
            print ("\t\tPlotname: {0}\n".format( plotname ))

        wantooverflow = True

        list_formats = [ plotname + '.png' ] #, plotname + '_canvas.root' ]
        if args.doEPS:
            list_formats.append( plotname + '.eps' )

        # Here is where the plotting is actually performed!
        #
        hists[category.name + ' ' + var.shortname] = background.plot( var,
                                                                      cut=cut,
                                                                      eventweight=combined_SF_weight,
                                                                      category=category,
                                                                      signal='',#'125',
                                                                      signalfactor=signalfactor,
                                                                      overridebackground=plotbackgrounds,
                                                                      systematics=systematics,
                                                                      systematicsdirection=systematicsdirection,
                                                                      overflowbins=wantooverflow,
                                                                      showratio=doShowRatio,
                                                                      wait=False,
                                                                      save=list_formats,
                                                                      log=None,
                                                                      logx=None
                                                                      )

        # Creating a file with the observed and expected distributions and systematics.
        # We fit them for TES uncertainty studies
        #
        foutput = TFile(plotname + '.root','RECREATE')
        if ( 'Mll01' in var.shortname ) or ( 'NJets' in var.shortname ):
            outfile = open(plotname + '_yields.txt', 'w')

        if args.doSyst:

            # systematics go into a different folder
            #
            dirname = dirname.replace(' ', '_') + '_Syst'

            # loop on the defined systematics
            #
            total_syst      = 0.0
            total_syst_up   = 0.0
            total_syst_down = 0.0
            histograms_syst = {}

            for syst in vardb.systlist:
                try:
                    os.makedirs(dirname)
                except:
                    pass
                plotname = dirname + '/' + category.name + ' ' + var.shortname + ' ' + syst.name
                plotname = plotname.replace(' ', '_')
                #
                # plotSystematics is the function which takes care of the systematics
                #
                systs[category.name + ' ' + var.shortname] = background.plotSystematics( syst,
                                                                                         var=var,
                                                                                         cut=cut,
                                                                                         eventweight=combined_SF_weight,
                                                                                         category=category,
                                                                                         overflowbins=wantooverflow,
                                                                                         showratio=True,
                                                                                         wait=False,
                                                                                         save=[plotname+'.png']
                                                                                         )

                # Obtains the total MC histograms with a particular systematics shifted and saving it in the ROOT file
                #
		print ("\t\t\tPlotname: {0}\n".format( plotname ))
		print ("\t\t\tSystematic: {0}\n".format( syst.name ))

                systobs, systnom, systup, systdown, systlistup, systlistdown = systs[category.name + ' ' + var.shortname]

                histograms_syst['Expected_'+syst.name+'_up']=systup
                histograms_syst['Expected_'+syst.name+'_up'].SetNameTitle(histname['Expected']+'_'+syst.name+'_up','')
                histograms_syst['Expected_'+syst.name+'_up'].SetLineColor(histcolour['Expected'])
                histograms_syst['Expected_'+syst.name+'_up'].Write()
                histograms_syst['Expected_'+syst.name+'_down']=systdown
                histograms_syst['Expected_'+syst.name+'_down'].SetNameTitle(histname['Expected']+'_'+syst.name+'_down','')
                histograms_syst['Expected_'+syst.name+'_down'].SetLineColor(histcolour['Expected'])
                histograms_syst['Expected_'+syst.name+'_down'].Write()
                for samp in ttH2015.backgrounds:
                    histograms_syst[samp+'_'+syst.name+'_up'] = systlistup[samp]
                    histograms_syst[samp+'_'+syst.name+'_up'].SetNameTitle(histname[samp]+'_'+syst.name+'_up','')
                    histograms_syst[samp+'_'+syst.name+'_up'].SetLineColor(histcolour[samp])
                    histograms_syst[samp+'_'+syst.name+'_up'].Write()
                    histograms_syst[samp+'_'+syst.name+'_down'] = systlistdown[samp]
                    histograms_syst[samp+'_'+syst.name+'_down'].SetNameTitle(histname[samp]+'_'+syst.name+'_down','')
                    histograms_syst[samp+'_'+syst.name+'_down'].SetLineColor(histcolour[samp])
                    histograms_syst[samp+'_'+syst.name+'_down'].Write()
                #ACTUALLY THE CODE DOES NOT CONSIDER SYSTEMATICS FOR THE SIGNAL. PUT IT AMONG THE BACKGROUNDS IF YOU WANT SYST ON IT
                if ( 'Mll01' in var.shortname ) or ( 'NJets' in var.shortname ):
                    outfile.write('Integral syst: \n')
                    outfile.write('syst %s up:   delta_yields = %f \n' %(syst.name,(systup.Integral()-systnom.Integral())))
                    outfile.write('syst %s down: delta_yields = %f \n' %(syst.name,(systdown.Integral()-systnom.Integral())))
                    if ( args.debug ):
                        outfile.write('GetEntries syst: \n')
                        outfile.write('syst %s up:   delta_entries %f \n' %(syst.name,(systup.GetEntries()-systnom.GetEntries())))
                        outfile.write('syst %s down: delta_entries %f \n' %(syst.name,(systdown.GetEntries()-systnom.GetEntries())))
                total_syst = total_syst + (systup.Integral()-systdown.Integral())/2.0*(systup.Integral()-systdown.Integral())/2.0
                total_syst_up   += (systup.Integral()-systnom.Integral())*(systup.Integral()-systnom.Integral())
                total_syst_down += (systdown.Integral()-systnom.Integral())*(systdown.Integral()-systnom.Integral())

            total_syst      = math.sqrt(total_syst)
            total_syst_up   = math.sqrt(total_syst_up)
            total_syst_down = math.sqrt(total_syst_down)

            if ( 'Mll01' in var.shortname ) or ( 'NJets' in var.shortname ):
                outfile.write('yields total syst UP: %f \n' %(total_syst_up))
                outfile.write('yields total syst DN: %f \n' %(total_syst_down))
                outfile.write('yields total syst: %f \n' %(total_syst))

        # Obtains the histograms correctly normalized
        #
        mclist, expected, observed, signal, _ = hists[category.name + ' ' + var.shortname]
        histograms = {}

        for samp in ttH2015.observed:
            histograms[samp] = observed
        if ttH2015.backgrounds:
            histograms['Expected']=expected
            for samp in ttH2015.backgrounds:
                histograms[samp] = mclist[samp]
                #in case you have to add other histograms you maybe prefer to use the method clone:
                #histograms[samp] = mclist[samp].Clone(histname[samp])
        if ttH2015.signals:
            for samp in ttH2015.signals:
                histograms[samp] = signal

        #print histograms
        #
        for samp in histograms.keys():
            histograms[samp].SetNameTitle(histname[samp],'')
            histograms[samp].SetLineColor(histcolour[samp])

        if ( 'Mll01' in var.shortname ) or ( var.shortname == 'NJets' ):
	    print (" ")
            print ("\t\tCategory: {0} - Variable: {1}\n".format( category.name, var.shortname ))
	    print ("\t\tIntegral:\n")
            outfile.write('Category: %s \n' %(category.name))
            outfile.write('Variable: %s \n' %(var.shortname))
            outfile.write('Integral: \n')
            err=Double(0)  # integral error
            value=0        # integral value
            for samp in histograms.keys():
                value = histograms[samp].IntegralAndError(1,histograms[samp].GetNbinsX(),err)
                print ("\t\t{0}: {1} +- {2}".format( histname[samp], value, err ))
                outfile.write('yields %s: %f +- %f \n' %(histname[samp], value, err))
                if ( var.shortname == 'NJets' ):
                    for bin in range(1,histograms[samp].GetNbinsX()):
                        err_bin=Double(0)
                        value_bin=0
                        value_bin=histograms[samp].GetBinContent(bin)
                        err_bin=histograms[samp].GetBinError(bin)
                        print ("\t\t {0} - {1}-jets bin: {2} +- {3}".format( histname[samp], bin-1, value_bin, err_bin ))
                    # Get integral and error from njets=5 bin onwards in one go!
		    err_HJ   = Double(0)
                    value_HJ = histograms[samp].IntegralAndError (6,histograms[samp].GetNbinsX(),err_HJ)
                    print ("\t\t {0} - >=5-jets bin: {1} +- {2}".format( histname[samp], value_HJ, err_HJ ))
	    print ("\n\t\tGetEntries:\n")
            outfile.write('GetEntries: \n')
            for samp in histograms.keys():
                print ("\t\t{0}: {1}".format( histname[samp], histograms[samp].GetEntries() ))
                outfile.write('entries %s: %f \n' %(histname[samp], histograms[samp].GetEntries()))

        for samp in histograms.keys():
                    histograms[samp].Write()
        foutput.Close()

        if ( 'Mll01' in var.shortname ) or ( var.shortname == 'NJets' ):
            outfile.close()

	#"""
