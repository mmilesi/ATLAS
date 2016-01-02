 #!/usr/bin/python
 #
 # *********************************************************************************
 # A python plotting script for the main HTopMultilep RunII analysis
 # 
 #
 # Authors:
 #  Marco Milesi ( marco.milesi@cern.ch ), Francesco Nuti ( francesco.nuti@cern.ch )
 #
 # *********************************************************************************

import os 
import sys
import math

sys.path.append(os.path.abspath(os.path.curdir))

# -------------------------------
# Parser for command line options 
# -------------------------------
import argparse

parser = argparse.ArgumentParser(description='Plotting python macro for ttH to multileptons analysis.')

#***********************************
# positional arguments (compulsory!)
#***********************************
parser.add_argument('inputDir', metavar='inputDir',type=str,
                   help='path to the directory containing input files')
parser.add_argument('samplesCSV', metavar='samplesCSV',type=str,
                   help='path to the csv file containing the processes of interest with their cross sections and other metadata')
#*******************
# optional arguments
#*******************
parser.add_argument('--selection', dest='selection', action='store', default='', type=str,
                    help='the selection chosen (e.g. noTypeOriginMatch, noTruthMatch... )')
parser.add_argument('--channel', dest='channel', action='store', default='TwoLepSR', type=str,  
		    help='the channel chosen (e.g. TwoLepSR, ThreeLepSR... )')
parser.add_argument('--outdirname', dest='outdirname', action='store', default='', type=str,  
		    help='specify a name to append to the output directory')
parser.add_argument('--fakeMethod', dest='fakeMethod', action='store', default='MC', type=str,
		    help='the fake estimation method chosen ( MC,MM,FF,ABCD )')
parser.add_argument('--lepFlavComp', dest='lepFlavComp', action='store', default='', type=str,
                    help='flavour composition of the SS pair ( ee, mm, em )')
parser.add_argument('--doLogScaleX', dest='doLogScaleX', action='store_true',
                    help='use log scale on the X axis')
parser.add_argument('--doLogScaleY', dest='doLogScaleY', action='store_true',
                    help='use log scale on the Y axis')	
parser.add_argument('--doSyst', dest='doSyst', action='store_true',
                    help='run systematics')
parser.add_argument('--debug', dest='debug', action='store_true',
                    help='run in debug mode')
parser.add_argument('--noSignal', action='store_true', dest='noSignal',
                    help='exclude signal')
parser.add_argument('--noStandardPlots', action='store_true', dest='noStandardPlots',
                    help='exclude all standard plots')		    
parser.add_argument('--doEPS', action='store_true', dest='doEPS',
                    help='make a .eps output file (NB: heavy size!)')
parser.add_argument('--doChFlipRate', dest='doChFlipRate', action='store_true',
                    help='measure charge flip rate in MC (to be used with MMClosureRates channel')		    

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

# -------------
# Check channel
# -------------
list_available_channel = ['TwoLepSR','ThreeLepSR','FourLepSR','Mu_RFRate_CR', 'El_RFRate_CR', 
			  'TwoLepLowNJetCR', 'ThreeLepLowNJetCR', 
			  'WZonCR', 'WZoffCR', 'WZHFonCR', 'WZHFoffCR', 
			  'ttWCR', 'ttZCR', 
			  'ZSSpeakCR', 'DataMC', 'MMClosureTest', 'MMClosureRates']

try:
    args.channel in list_available_channel
except ValueError:
    sys.exit('the channel specified (', args.channel ,') is incorrect. Must be one of: ', list_available_channel)

doTwoLepSR      	= bool( args.channel == 'TwoLepSR' )
doThreeLepSR    	= bool( args.channel == 'ThreeLepSR' )
doFourLepSR     	= bool( args.channel == 'FourLepSR' )
doMuRFRateCR          	= bool( args.channel == 'Mu_RFRate_CR' )
doElRFRateCR          	= bool( args.channel == 'El_RFRate_CR' )
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
doOtherCR = (doWZonCR or doWZoffCR or doWZHFonCR or doWZHFoffCR or dottWCR or dottZCR or doZSSpeakCR or doDataMCCR or doMMClosureTest or doMMClosureRates) 

# ------------------------------------------------
# make standard plots unless differently specified
# ------------------------------------------------
doStandardPlots = False if (args.noStandardPlots) else (doSR or doLowNJetCR or doOtherCR)

# ----------------------------
# Check fake estimation method
# ----------------------------
list_available_fakemethod = ['MC','MM','FF']

try:
    args.fakeMethod in list_available_fakemethod
except ValueError:	
    sys.exit('the Fake Method specified (', args.fakeMethod ,') is incorrect. Must be one of ', list_available_fakemethod)

doMM   = bool( args.fakeMethod == 'MM' )
doFF   = bool( args.fakeMethod == 'FF' )
doABCD = bool( args.fakeMethod == 'ABCD' )

# -----------------------------------------------
# Check lepton flavour composition of the SS pair
# -----------------------------------------------
list_available_flavcomp = ['ee','mm','em', ''] 

try:
    args.lepFlavComp in list_available_flavcomp
except ValueError:	 
    sys.exit('ERROR: the lepton flavour composition of the SS pair is incorrect. Must be one of ', list_available_flavcomp)
else: 
    if not ( args.lepFlavComp ):	 
        print 'WARNING: lepton flavour composition of the SS pair unspecified ( i.e., empty string '' )... will not use any restriction on the lepton flavour'
	

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
                        inputdir    = args.inputDir + args.selection,
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

# -----------------------------------------------------
# The list of event-level cuts ( a-la TTree->Draw("") )
#
# WARNING:
# To avoid unexpected behaviour, 
# ALWAYS enclose the cut string in '()'!!!
#
# -----------------------------------------------------
vardb.registerCut( Cut('IsMC',        '( isMC == 1 )') )
vardb.registerCut( Cut('TrigDec',     '( passHLT == 1 && ( passedTriggers == \"HLT_e24_lhmedium_iloose_L1EM18VH\" || passedTriggers == \"HLT_e60_lhmedium\" || passedTriggers == \"HLT_e120_lhloose\" || passedTriggers == \"HLT_mu20_iloose_L1MU15\" || passedTriggers == \"HLT_mu50\" ) )') )
vardb.registerCut( Cut('TrigMatch',   '( ( lep_flavour[0] == 11 && lep_pt[0] > 28e3 ) || ( lep_flavour[0] == 13 && lep_pt[0] > 25e3 )  )') ) # for the time being, just require the leading lepton pT over a high-enough threshold
vardb.registerCut( Cut('NBJet',       '( njets_mv2c20_Fix70 > 0 )') )
vardb.registerCut( Cut('LargeNBJet',  '( njets_mv2c20_Fix70 > 1 )') )
vardb.registerCut( Cut('BJetVeto',    '( njets_mv2c20_Fix70 == 0 )') )
vardb.registerCut( Cut('TauVeto',     '( ntau == 0 )') )
vardb.registerCut( Cut('OneTau',      '( ntau == 1 )') )
vardb.registerCut( Cut('NJet2L',      '( njet >= 4 )') )
vardb.registerCut( Cut('NJet3L',      '( njet >= 4 || ( njets_mv2c20_Fix70 > 1 && njet == 3 ) )') )
vardb.registerCut( Cut('NJet4L',      '( njet >= 2 )') )
vardb.registerCut( Cut('LowJetCR',    '( njet > 0 && njet < 4 )') )
vardb.registerCut( Cut('2Lep',        '( nlep == 2 && ( lep_pt[0] > 20e3 && lep_pt[1] > 20e3 ) )') )
vardb.registerCut( Cut('2LepTau',     '( nlep == 2 && ntau > 0 && ( lep_charge[0] * tau_charge[0] ) < 0 && ( lep_pt[0] > 15e3 && lep_pt[1] > 15e3 ) )') )
vardb.registerCut( Cut('SS',          '( isSS01 == 1 )') )
vardb.registerCut( Cut('3Lep',        '( nlep == 3 && isSS12 == 1 && lep_isTightSelected[0] == 1 && TMath::Abs( lep_charge[0] + lep_charge[1] + lep_charge[2] ) == 1 && lep_pt[1] > 20e3 && lep_pt[2] > 20e3 && mll01 > 12e3 && mll02 > 12e3 )') )
vardb.registerCut( Cut('4Lep',        '( nlep == 4 && lep_pt[0] > 25e3 && lep_pt[1] > 15e3 && lep_isTightSelected[0] == 1 && lep_isTightSelected[1] == 1 && lep_isTightSelected[2] == 1 && lep_isTightSelected[3] == 1 && TMath::Abs( lep_charge[0] + lep_charge[1] + lep_charge[2] + lep_charge[3] ) == 0 && ( ( mJPsiCand_ee > 10e3 || mJPsiCand_ee < 0.0 ) && ( mJPsiCand_mm > 10e3 || mJPsiCand_mm < 0.0 ) ) )') )
vardb.registerCut( Cut('SF_Event',    '( nmuon == 2 || nel == 2 )') )
vardb.registerCut( Cut('MuMu_Event',  '( nmuon == 2 )') )
vardb.registerCut( Cut('ElEl_Event',  '( nel == 2 )') )
vardb.registerCut( Cut('OF_Event',    '( nmuon == 1 && nel == 1 )') )
vardb.registerCut( Cut('MuEl_Event',  '( nmuon == 1 && nel == 1 && lep_flavour[0] == 13 )') )
vardb.registerCut( Cut('ElMu_Event',  '( nmuon == 1 && nel == 1 && lep_flavour[0] == 11 )') )
# temporary, to be used for 2 lep SF, OS CR (Z peak): keep events where leading jet |eta| in the barrel
#
vardb.registerCut( Cut('Jet0FwdCut',  '( TMath::Abs(jet_eta[0]) < 2.5 )') )
# apply to reduce charge flip contamination
#
vardb.registerCut( Cut('OF_ElEtaCut', '( TMath::Abs(el_eta[0]) < 1.5 )') )
vardb.registerCut( Cut('SF_ElEtaCut', '( ( TMath::Abs(el_eta[0]) < 1.5 && TMath::Abs(el_eta[1]) < 1.5 ) )') )

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
# old definition: specify the !prompt lepton to be from a HF decay
#
#vardb.registerCut( Cut('2Lep_NonPromptEvent', '( isMC==0 || ( isMC==1 && ( lep_truthType[0] == 7 || lep_truthType[0] == 3 || lep_truthType[1] == 7 || lep_truthType[1] == 3 ) && ( lep_isChFlip[0] == 0 && lep_isChFlip[1] == 0 ) ) )') )
#
# 3.
# event passes this cut if ALL leptons are !non-prompt (MCTruthClassifier --> !(NonIso))
# ---> USED FOR PROMPT (AND CH-FLIP?) SUB IN FAKE CR ----> TO BE CHECKED
vardb.registerCut( Cut('2Lep_NonNonPromptEvent',  '( isMC==0 || ( isMC==1 && ( lep_truthType[0] != 7 && lep_truthType[0] != 3 && lep_truthType[1] != 7 && lep_truthType[1] != 3 ) ) )') )
# 4.
# event passes this cut if AT LEAST ONE lepton is charge flip (does not distinguish trident VS charge-misId)
#
vardb.registerCut( Cut('2Lep_ChFlipEvent',   '( isMC==0 || ( isMC==1 && ( lep_isChFlip[0] == 1 || lep_isChFlip[1] == 1 ) ) )') )
# 4a.
# event passes this cut if AT LEAST ONE lepton is (prompt and charge flip) (it will be a charge-misId charge flip)
#
vardb.registerCut( Cut('2Lep_ChFlipPromptEvent',  '( isMC==0 || ( isMC==1 && ( ( lep_isChFlip[0] == 1 && ( lep_truthType[0] == 6 || lep_truthType[0] == 2 ) ) || ( lep_isChFlip[1] == 1 && ( lep_truthType[1] == 6 || lep_truthType[1] == 2 ) ) ) ) )') )
# 4b.
# event passes this cut if AT LEAST ONE object is (!prompt and charge flip and from bremsstrahlung) (this will be a trident charge flip)
#
vardb.registerCut( Cut('2Lep_ChFlipBremEvent', '( isMC==0 || ( isMC==1 && ( ( lep_isChFlip[0] == 1 && lep_isBrem[0] == 1 && ( lep_truthType[0] != 6 && lep_truthType[0] != 2 ) ) || ( lep_isChFlip[1] == 1 && lep_isBrem[1] == 1 && ( lep_truthType[1] != 6 && lep_truthType[1] != 2 ) ) ) ) )') ) 
# 4c.
# event passes this cut if AT LEAST ONE lepton is (!prompt and charge flip) 
#
vardb.registerCut( Cut('2Lep_ChFlipNonPromptEvent', '( isMC==0 || ( isMC==1 && ( ( lep_isChFlip[0] == 1 && ( lep_truthType[0] != 6 && lep_truthType[0] != 2 ) ) || ( lep_isChFlip[1] == 1 && ( lep_truthType[1] != 6 && lep_truthType[1] != 2 ) ) ) ) )') )
# 5.
# event passes this cut if NONE of the leptons is charge flip
vardb.registerCut( Cut('2Lep_ChFlipVeto',   '( isMC==0 || ( isMC==1 && ( lep_isChFlip[0] == 0 && lep_isChFlip[1] == 0 ) ) )') )

# -------------------------------------------------------------------------------

# compute invariant masses on-the-fly

lep0_px = '( lep_pt[0] * TMath::Cos( lep_phi[0] ) )' 
lep0_py = '( lep_pt[0] * TMath::Sin( lep_phi[0] ) )' 
lep0_pz = '( lep_pt[0] * TMath::SinH( lep_eta[0] ) )' 
lep0_E  = '( lep_pt[0] * TMath::CosH( lep_eta[0] ) )' 

lep1_px = '( lep_pt[1] * TMath::Cos( lep_phi[1] ) )' 
lep1_py = '( lep_pt[1] * TMath::Sin( lep_phi[1] ) )' 
lep1_pz = '( lep_pt[1] * TMath::SinH( lep_eta[1] ) )' 
lep1_E  = '( lep_pt[1] * TMath::CosH( lep_eta[1] ) )' 

# neglecting lepton masses
m_ll = '( TMath::Sqrt( 2*( ' + lep0_E + '*' + lep1_E + '-' + lep0_px + '*' + lep1_px + '-' + lep0_py + '*' + lep1_py + '-' + lep0_pz + '*' + lep1_pz + ' ) ) )'

if args.debug:
    print 'string for invariant mass of a generic lepton pair: ', m_ll

mZCand_sidescut = '( TMath::Abs( ' + m_ll + ' - 91.187e3 ) > 30e3 )'
mZCand_peakcut  = '( TMath::Abs( ' + m_ll + ' - 91.187e3 ) < 30e3 )'

# Powheg Z+jets has a cut at 60 GeV
mZCand_mincut  = '( ' + m_ll + '  > 60e3 )'

if args.debug:
    print 'string for selecting events w/ 2 leptons outside Z peak: ',   mZCand_sidescut
    print 'string for selecting events w/ 2 leptons around Z peak: ',    mZCand_peakcut

vardb.registerCut( Cut('Zsidescut', mZCand_sidescut ) )   # use this to require the 2 leptons to be outside Z peak
vardb.registerCut( Cut('Zpeakcut',  mZCand_peakcut )  )   # use this to require the 2 leptons to be around Z peak
vardb.registerCut( Cut('Zmincut',   mZCand_mincut )   )   # Powheg Z+jets has a cut at 60 GeV

# reconstructed pT of the Z
#
pT_Z = '( TMath::Sqrt( (lep_pt[0]*lep_pt[0]) + (lep_pt[1]*lep_pt[1]) + 2*lep_pt[0]*lep_pt[1]*(TMath::Cos( lep_phi[0] - lep_phi[1] )) ) )/1e3'

# -------------------------------------------------------------------
#  Used for the fake estimate in the 2lepSS and 3lep channel
#
# 'isXY' refers to the SS pair in a 2Lep or trilep event. 
#    X is the lepton w/ highest pT of the pair, Y is the second.
#    X,Y = 'L' means that the lepton is ( loose && !tight )
#
#
#  NB: in SR, the TT cut is applied within the ttH2015_Background*.py classes
#      ( except for FakesMM, where the MM weight defines automatically the category )
#  NB: this cut will be used ONLY when looking at 2lep SS or 3lep events!!!
#
# -------------------------------------------------------------------
vardb.registerCut( Cut('TT',  '( isTT == 1 )') )
vardb.registerCut( Cut('TL',  '( isTL == 1 )') ) 
vardb.registerCut( Cut('LT',  '( isLT == 1 )') ) 
vardb.registerCut( Cut('LL',  '( isLL == 1 )') ) 

# ---------------------------
# A list of variables to plot
# ---------------------------

if doStandardPlots :
    print ''    
    if doMMClosureTest:
        vardb.registerVar( Variable(shortname = 'NJets',    	   latexname = 'Jet multiplicity',		    	     ntuplename = 'njets',			      bins = 10,  minval = 0,	 maxval = 10) )
    	vardb.registerVar( Variable(shortname = 'Mll01_inc',       latexname = 'm(l_{0}l_{1}) [GeV]',               	     ntuplename = 'mll01/1e3',  		      bins = 13,  minval = 0.0,  maxval = 260.0,) )

    #vardb.registerVar( Variable(shortname = 'Jet0Pt',   	   latexname = 'p_{T}^{lead jet} [GeV]',	    	     ntuplename = 'jet_pt[0]/1e3',		      bins = 36, minval = 20.0,  maxval = 200.0,) )
    #vardb.registerVar( Variable(shortname = 'Jet0Eta',		   latexname = '#eta^{lead jet}',  			     ntuplename = 'jet_eta[0]',	                      bins = 50,  minval = -5.0, maxval = 5.0) )
    #vardb.registerVar( Variable(shortname = 'NJets',    	   latexname = 'Jet multiplicity',		    	     ntuplename = 'njets',			      bins = 10,  minval = 0,	 maxval = 10) )
    #vardb.registerVar( Variable(shortname = 'NBJets',   	   latexname = 'BJet multiplicity',		    	     ntuplename = 'njets_mv2c20_Fix70',		      bins = 4,  minval = 0,	 maxval = 4) )
    #vardb.registerVar( Variable(shortname = 'NJetsPlus10NBJets',   latexname = 'N_{Jets}+10*N_{BJets}',  		     ntuplename = 'njets+10.0*njets_mv2c20_Fix70',    bins = 40,  minval = 0,	 maxval = 40) )
    #
    # Inclusive m(ll) plot
    #    
    #vardb.registerVar( Variable(shortname = 'Mll01_inc',    	   latexname = 'm(l_{0}l_{1}) [GeV]',               	     ntuplename = 'mll01/1e3',  		      bins = 40, minval = 60.0,  maxval = 260.0,) )
    #
    # Z peak plot
    #
    #vardb.registerVar( Variable(shortname = 'Mll01_peak',    	   latexname = 'm(l_{0}l_{1}) [GeV]',               	     ntuplename = 'mll01/1e3',  		      bins = 30, minval = 60.0,   maxval = 120.0,) )
    #
    #vardb.registerVar( Variable(shortname = 'pT_Z',    	           latexname = 'p_{T} Z (reco) [GeV]',               	     ntuplename = pT_Z,  		              bins = 100, minval = 0.0,   maxval = 1000.0, logaxisX = True) )
    #
    #vardb.registerVar( Variable(shortname = 'Lep0Pt',		   latexname = 'p_{T}^{lead lep} [GeV]',		     ntuplename = 'lep_pt[0]/1e3',		      bins = 11, minval = 20.0,  maxval = 240.0,) )
    #vardb.registerVar( Variable(shortname = 'Lep1Pt',		   latexname = 'p_{T}^{2nd lead lep} [GeV]',	             ntuplename = 'lep_pt[1]/1e3',		      bins = 7,  minval = 20.0,  maxval = 160.0,) )
    #vardb.registerVar( Variable(shortname = 'Lep0Eta',		   latexname = '#eta^{lead lep}',  			     ntuplename = 'TMath::Abs(lep_eta[0])',	      bins = 8,  minval = 0.0,	 maxval = 2.6) )
    #vardb.registerVar( Variable(shortname = 'Lep1Eta',		   latexname = '#eta^{2nd lead lep}',			     ntuplename = 'TMath::Abs(lep_eta[1])',	      bins = 8,  minval = 0.0,	 maxval = 2.6) )
    
    #vardb.registerVar( Variable(shortname = 'Mll12',              latexname = 'm(l_{1}l_{2}) [GeV]',	          	     ntuplename = 'mll12/1e3',  	              bins = 15, minval = 0.0,   maxval = 300.0,) )
    #vardb.registerVar( Variable(shortname = 'avgint',  	           latexname = 'Average Interactions Per Bunch Crossing',    ntuplename = 'averageInteractionsPerCrossing',   bins = 50, minval = 0,	 maxval = 50,  typeval = TH1I) )
    #vardb.registerVar( Variable(shortname = 'MET_FinalClus',       latexname = 'E_{T}^{miss} (FinalClus) [GeV]',	     ntuplename = 'metFinalClus/1e3',	              bins = 45, minval = 0.0,   maxval = 180.0,) )
    #vardb.registerVar( Variable(shortname = 'MET_FinalTrk',        latexname = 'E_{T}^{miss} (FinalTrk) [GeV]',		     ntuplename = 'metFinalTrk/1e3',	              bins = 45, minval = 0.0,   maxval = 180.0,) )
    #vardb.registerVar( Variable(shortname = 'MET_SoftClus',        latexname = 'E_{T}^{miss} (SoftClus) [GeV]',	             ntuplename = 'metSoftClus/1e3',	              bins = 45, minval = 0.0,   maxval = 180.0,) )
    #vardb.registerVar( Variable(shortname = 'MET_SoftTrk',         latexname = 'E_{T}^{miss} (SoftTrk) [GeV]',		     ntuplename = 'metSoftTrk/1e3',	              bins = 45, minval = 0.0,   maxval = 180.0,) )
    #vardb.registerVar( Variable(shortname = 'MET_Electrons',	   latexname = 'E_{T}^{miss} (Electrons) [GeV]',	     ntuplename = 'metEle/1e3', 		      bins = 45, minval = 0.0,   maxval = 180.0,) )
    #vardb.registerVar( Variable(shortname = 'MET_Muons',	   latexname = 'E_{T}^{miss} (Muons) [GeV]',		     ntuplename = 'metMuons/1e3',		      bins = 45, minval = 0.0,   maxval = 180.0,) )
    #vardb.registerVar( Variable(shortname = 'MET_Jets',	           latexname = 'E_{T}^{miss} (Jets) [GeV]',		     ntuplename = 'metJet/1e3', 		      bins = 45, minval = 0.0,   maxval = 180.0,) )

    #vardb.registerVar( Variable(shortname = 'MT_Lep0MET',     	   latexname = 'm_{T}(l_{0},MET) [GeV]',	  	     ntuplename = 'mT_lep0MET/1e3',		      bins = 40, minval = 0.0,   maxval = 160.0,) )
    #vardb.registerVar( Variable(shortname = 'MT_Lep1MET',     	   latexname = 'm_{T}(l_{1},MET) [GeV]',	  	     ntuplename = 'mT_lep1MET/1e3',		      bins = 40, minval = 0.0,   maxval = 160.0,) )
    #vardb.registerVar( Variable(shortname = 'Tau0Pt',  	   latexname = 'p_{T}^{lead tau} [GeV]',	  	     ntuplename = 'tau_pt[0]/1e3',	              bins = 30, minval = 25.0,  maxval = 100.0,) )
    
    #vardb.registerVar( Variable(shortname = 'El0Pt',		  latexname = 'p_{T}^{lead e} [GeV]',			    ntuplename = 'el_pt[0]/1e3',		     bins = 36, minval = 10.0,  maxval = 190.0,) )
    #vardb.registerVar( Variable(shortname = 'El1Pt',		  latexname = 'p_{T}^{2nd lead e} [GeV]',		    ntuplename = 'el_pt[1]/1e3',		     bins = 36, minval = 10.0,  maxval = 190.0,) )
    #vardb.registerVar( Variable(shortname = 'El0Eta',		  latexname = '#eta^{lead e}',  			    ntuplename = 'TMath::Abs(el_eta[0])',	     bins = 8,  minval = 0.0,	maxval = 2.6, manualbins = [ 0.0, 0.5, 0.8, 1.1, 1.37, 1.52, 2.0, 2.25, 2.6]) )
    #vardb.registerVar( Variable(shortname = 'El1Eta',		  latexname = '#eta^{2nd lead e}',			    ntuplename = 'TMath::Abs(el_eta[1])',	     bins = 8,  minval = 0.0,	maxval = 2.6, manualbins = [ 0.0, 0.5, 0.8, 1.1, 1.37, 1.52, 2.0, 2.25, 2.6]) )
    #vardb.registerVar( Variable(shortname = 'El0TopoEtCone20',     latexname = 'topoetcone20^{lead e} [GeV]',                ntuplename = 'el_topoetcone20[0]/1e3',           bins = 40, minval = 0.0,   maxval = 10.0, manualbins = [ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0] ) )
    #vardb.registerVar( Variable(shortname = 'El1TopoEtCone20',     latexname = 'topoetcone20^{2nd lead e} [GeV]',            ntuplename = 'el_topoetcone20[1]/1e3',           bins = 40, minval = 0.0,   maxval = 10.0, manualbins = [ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0] ) )
    #vardb.registerVar( Variable(shortname = 'El0PtVarCone20',      latexname = 'ptvarcone20^{lead e} [GeV]',                 ntuplename = 'el_ptvarcone20[0]/1e3',            bins = 40, minval = 1.0,   maxval = 5.0) )
    #vardb.registerVar( Variable(shortname = 'El1PtVarCone20',      latexname = 'ptvarcone20^{2nd lead e} [GeV]',             ntuplename = 'el_ptvarcone20[1]/1e3',            bins = 40, minval = 1.0,   maxval = 5.0) )
#    vardb.registerVar( Variable(shortname = 'El0TopoEtCone20OverPt', latexname = 'topoetcone20/p_{T} lead e [GeV]',	    ntuplename = 'el_topoetcone20[0]/el_pt[0]',  bins = 50, minval = -0.2,  maxval = 0.8) )
#    vardb.registerVar( Variable(shortname = 'El1TopoEtCone20OverPt', latexname = 'topoetcone20/p_{T} 2nd lead e [GeV]',    ntuplename = 'el_topoetcone20[1]/el_pt[1]',  bins = 50, minval = -0.2,  maxval = 0.8) )
#    vardb.registerVar( Variable(shortname = 'El0PtVarCone20OverPt',  latexname = 'ptvarcone20/p_{T} lead e [GeV]',	    ntuplename = 'el_ptvarcone20[0]/el_pt[0]',   bins = 50, minval = 0.0,  maxval = 1.0) )
#    vardb.registerVar( Variable(shortname = 'El1PtVarCone20OverPt',  latexname = 'ptvarcone20/p_{T} 2nd lead e [GeV]',     ntuplename = 'el_ptvarcone20[1]/el_pt[1]',   bins = 50, minval = 0.0,  maxval = 1.0) )
#    vardb.registerVar( Variable(shortname = 'El0d0sig',	    latexname = '|d_{0}^{sig}| lead e', 		      ntuplename = 'el_trkd0sig[0]',		       bins = 40, minval = 0.0,  maxval = 10.0,) )
#    vardb.registerVar( Variable(shortname = 'El1d0sig',	    latexname = '|d_{0}^{sig}| 2nd lead e',		      ntuplename = 'el_trkd0sig[1]',		       bins = 40, minval = 0.0,  maxval = 10.0,) )
#    vardb.registerVar( Variable(shortname = 'El0z0sintheta',	    latexname = 'z_{0}*sin(#theta) lead e [mm]',	      ntuplename = 'el_trkz0sintheta[0]',	       bins = 20, minval = -1.0,  maxval = 1.0,) )
#    vardb.registerVar( Variable(shortname = 'El1z0sintheta',	    latexname = 'z_{0}*sin(#theta) 2nd lead e [mm]',	      ntuplename = 'el_trkz0sintheta[1]',	       bins = 20, minval = -1.0,  maxval = 1.0,) )
#    vardb.registerVar( Variable(shortname = 'El0LHTight',	   latexname = 'lead e IsLHTight',			     ntuplename = 'el_LHTight[0]',		      bins = 2, minval = -0.5,  maxval = 1.5,) )
#    vardb.registerVar( Variable(shortname = 'El1LHTight',	   latexname = '2nd lead e IsLHTight',			     ntuplename = 'el_LHTight[1]',		      bins = 2, minval = -0.5,  maxval = 1.5,) )

#    vardb.registerVar( Variable(shortname = 'Mu0Pt',		   latexname = 'p_{T}^{lead #mu} [GeV]',		     ntuplename = 'muon_pt[0]/1e3',		      bins = 36, minval = 10.0,  maxval = 190.0,) )
#    vardb.registerVar( Variable(shortname = 'Mu1Pt',		   latexname = 'p_{T}^{2nd lead #mu} [GeV]',		     ntuplename = 'muon_pt[1]/1e3',		      bins = 36, minval = 10.0,  maxval = 190.0,) )
#    vardb.registerVar( Variable(shortname = 'Mu0Eta',  	   latexname = '#eta^{lead #mu}',			     ntuplename = 'muon_eta[0]',		      bins = 16,  minval = -2.6, maxval = 2.6,) )
#    vardb.registerVar( Variable(shortname = 'Mu1Eta',  	   latexname = '#eta^{2nd lead #mu}',			     ntuplename = 'muon_eta[1]',		      bins = 16,  minval = -2.6, maxval = 2.6,) )
    #vardb.registerVar( Variable(shortname = 'Mu0TopoEtCone20',    latexname = 'topoetcone20^{lead #mu} [GeV]', 	     ntuplename = 'muon_topoetcone20[0]/1e3',	      bins = 40, minval = 0.0,   maxval = 10.0, manualbins = [ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0]) )
    #vardb.registerVar( Variable(shortname = 'Mu1TopoEtCone20',    latexname = 'topoetcone20^{2nd lead #mu} [GeV]',	     ntuplename = 'muon_topoetcone20[1]/1e3',	      bins = 40, minval = 0.0,   maxval = 10.0, manualbins = [ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0]) )
    #vardb.registerVar( Variable(shortname = 'Mu0PtVarCone30',     latexname = 'ptvarcone20^{lead #mu} [GeV]',  	     ntuplename = 'muon_ptvarcone30[0]/1e3',	      bins = 40, minval = 1.0,   maxval = 5.0) )
    #vardb.registerVar( Variable(shortname = 'Mu1PtVarCone30',     latexname = 'ptvarcone20^{2nd lead #mu} [GeV]',	     ntuplename = 'muon_ptvarcone30[1]/1e3',	      bins = 40, minval = 1.0,   maxval = 5.0) )
#    vardb.registerVar( Variable(shortname = 'Mu0TopoEtCone20OverPt', latexname = 'topoetcone20/p_{T} lead #mu [GeV]',        ntuplename = 'muon_topoetcone20[0]/muon_pt[0]',	 bins = 50, minval = -0.2,  maxval = 0.8) )
#    vardb.registerVar( Variable(shortname = 'Mu1TopoEtCone20OverPt', latexname = 'topoetcone20/p_{T} 2nd lead #mu [GeV]',    ntuplename = 'muon_topoetcone20[1]/muon_pt[1]',	 bins = 50, minval = -0.2,  maxval = 0.8) )
#    vardb.registerVar( Variable(shortname = 'Mu0PtVarCone30OverPt',  latexname = 'ptvarcone30/p_{T} lead #mu [GeV]',	      ntuplename = 'muon_ptvarcone30[0]/muon_pt[0]',	 bins = 50, minval = 0.0,  maxval = 1.0) )
#    vardb.registerVar( Variable(shortname = 'Mu1PtVarCone30OverPt',  latexname = 'ptvarcone30/p_{T} 2nd lead #mu [GeV]',     ntuplename = 'muon_ptvarcone30[1]/muon_pt[1]',	 bins = 50, minval = 0.0,  maxval = 1.0) )
#    vardb.registerVar( Variable(shortname = 'Mu0d0sig',	   latexname = '|d_{0}^{sig}| lead #mu',		     ntuplename = 'muon_trkd0sig[0]',		      bins = 40, minval = 0.0,  maxval = 10.0,) )
#    vardb.registerVar( Variable(shortname = 'Mu1d0sig',	   latexname = '|d_{0}^{sig}| 2nd lead #mu',		     ntuplename = 'muon_trkd0sig[1]',		      bins = 40, minval = 0.0,  maxval = 10.0,) )
#    vardb.registerVar( Variable(shortname = 'Mu0z0sintheta',	   latexname = 'z_{0}*sin(#theta) lead #mu [mm]',	     ntuplename = 'muon_trkz0sintheta[0]',	      bins = 20, minval = -1.0,  maxval = 1.0,) )
#    vardb.registerVar( Variable(shortname = 'Mu1z0sintheta',	   latexname = 'z_{0}*sin(#theta) 2nd lead #mu [mm]',	     ntuplename = 'muon_trkz0sintheta[1]',	      bins = 20, minval = -1.0,  maxval = 1.0,) )

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
if doTwoLepSR :

    #
    # when using MM or FF for non-prompt bkg estimate, make sure you plot 
    # only pure prompt MC (to avoid double counting of fake background events!)
    #
    if ( args.fakeMethod == 'MM' or args.fakeMethod == 'FF' ):
    	 
	 #
	 # Preliminary selections
	 #
        
	 #vardb.registerCategory( MyCategory('Base',	            			     cut = vardb.getCuts(['NBJet']) ) )
	 #vardb.registerCategory( MyCategory('2LepBase',          			     cut = vardb.getCuts(['NBJet', '2Lep']) ) )
	 #vardb.registerCategory( MyCategory('2LepSSBase',        			     cut = vardb.getCuts(['NBJet', '2Lep','SS']) ) )

	 #
	 # MuMu region
	 #

	 #vardb.registerCategory( MyCategory('MuMuSSBase',                                    cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'MuMu_Event']) ) )
	 #vardb.registerCategory( MyCategory('MuMuTruthMSSBase',                              cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'MuMu_Event', '2Lep_PurePromptEvent']) ) )
	 #vardb.registerCategory( MyCategory('MuMuNoTauTruthMSSBase',                         cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'MuMu_Event', '2Lep_PurePromptEvent', 'TauVeto']) ) )
	 #vardb.registerCategory( MyCategory('MuMuTrigMatchNoTauTruthMSSBase',                cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'MuMu_Event', '2Lep_PurePromptEvent', 'TauVeto', 'TrigMatch']) ) )
	 vardb.registerCategory( MyCategory('MuMuHighJetTrigMatchNoTauTruthMSSBase',         cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'MuMu_Event', '2Lep_PurePromptEvent', 'TauVeto', 'TrigMatch', 'NJet2L']) ) )
	 vardb.registerCategory( MyCategory('MuMuLowJetTrigMatchNoTauTruthMSSBase',          cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'MuMu_Event', '2Lep_PurePromptEvent', 'TauVeto', 'TrigMatch', 'LowJetCR']) ) )

	 #
	 # OF region
	 #

	 #vardb.registerCategory( MyCategory('OFSSBase',                                    cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'OF_Event']) ) )
	 #vardb.registerCategory( MyCategory('OFTruthMSSBase',                              cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'OF_Event', '2Lep_PurePromptEvent']) ) )
	 #vardb.registerCategory( MyCategory('OFNoTauTruthMSSBase',                         cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'OF_Event', '2Lep_PurePromptEvent', 'TauVeto']) ) )
	 #vardb.registerCategory( MyCategory('OFTrigMatchNoTauTruthMSSBase',                cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'OF_Event', '2Lep_PurePromptEvent', 'TauVeto', 'TrigMatch']) ) )
	 #vardb.registerCategory( MyCategory('OFEtaCutTrigMatchNoTauTruthMSSBase',          cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'OF_Event', '2Lep_PurePromptEvent', 'TauVeto', 'TrigMatch', 'OF_ElEtaCut',]) ) )
	 vardb.registerCategory( MyCategory('OFHighJetEtaCutTrigMatchNoTauTruthMSSBase',   cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'OF_Event', '2Lep_PurePromptEvent', 'TauVeto', 'TrigMatch', 'OF_ElEtaCut', 'NJet2L']) ) )
	 vardb.registerCategory( MyCategory('OFLowJetEtaCutTrigMatchNoTauTruthMSSBase',    cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'OF_Event', '2Lep_PurePromptEvent', 'TauVeto', 'TrigMatch', 'OF_ElEtaCut', 'LowJetCR']) ) )

	 #
	 # ElEl region
	 #

	 #vardb.registerCategory( MyCategory('ElElSSBase',                                    cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'ElEl_Event']) ) )
	 #vardb.registerCategory( MyCategory('ElElTruthMSSBase',                              cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'ElEl_Event', '2Lep_PurePromptEvent']) ) )
	 #vardb.registerCategory( MyCategory('ElElNoTauTruthMSSBase',                         cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'ElEl_Event', '2Lep_PurePromptEvent', 'TauVeto']) ) )
	 #vardb.registerCategory( MyCategory('ElElTrigMatchNoTauTruthMSSBase',                cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'ElEl_Event', '2Lep_PurePromptEvent', 'TauVeto', 'TrigMatch']) ) )
	 #vardb.registerCategory( MyCategory('ElElEtaCutTrigMatchNoTauTruthMSSBase',          cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'ElEl_Event', '2Lep_PurePromptEvent', 'TauVeto', 'TrigMatch', 'SF_ElEtaCut',]) ) )
	 vardb.registerCategory( MyCategory('ElElHighJetEtaCutTrigMatchNoTauTruthMSSBase',   cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'ElEl_Event', '2Lep_PurePromptEvent', 'TauVeto', 'TrigMatch', 'SF_ElEtaCut', 'NJet2L']) ) )
	 vardb.registerCategory( MyCategory('ElElLowJetEtaCutTrigMatchNoTauTruthMSSBase',    cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'ElEl_Event', '2Lep_PurePromptEvent', 'TauVeto', 'TrigMatch', 'SF_ElEtaCut', 'LowJetCR']) ) )

    else:

         #vardb.registerCategory( MyCategory('MuMuSS',           			     cut = vardb.getCuts(['TrigMatch', 'NBJet', '2Lep', 'SS', 'MuMu_Event', 'TauVeto', 'TrigMatch', 'NJet2L']) ) )
         #vardb.registerCategory( MyCategory('OFSS',           			     cut = vardb.getCuts(['TrigMatch', 'NBJet', '2Lep', 'SS', 'OF_Event', 'TauVeto', 'TrigMatch', 'OF_ElEtaCut', 'NJet2L']) ) )
         #vardb.registerCategory( MyCategory('ElElSS',           			     cut = vardb.getCuts(['TrigMatch', 'NBJet', '2Lep', 'SS', 'ElEl_Event', 'TauVeto', 'TrigMatch', 'SF_ElEtaCut', 'NJet2L']) ) )
	 
         #
	 # 2lep+tau region
	 #
	 vardb.registerCategory( MyCategory('TwoLepSSTau',      			     cut = vardb.getCuts(['NBJet', '2LepTau', 'SS', 'NJet2L', 'OneTau', 'TrigMatch', 'Zsidescut']) ) )


if doThreeLepSR:
    vardb.registerCategory( MyCategory('ThreeLep',    cut = vardb.getCuts(['NBJet', '3Lep', 'TrigMatch', 'Zsidescut', 'NJet3L']) ) )

if doFourLepSR:
    vardb.registerCategory( MyCategory('FourLep',     cut = vardb.getCuts(['NBJet', '4Lep', 'TrigMatch', 'NJet4L']) ) )

# -------------
# low N-jet CRs
# -------------
if doTwoLepLowNJetCR :
    
    if ( args.fakeMethod == 'MM' or args.fakeMethod == 'FF' ):
	 #
	 # MuMu region
	 #
	 vardb.registerCategory( MyCategory('MuMuSSLowNJetCR_AllCuts',    cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'MuMu_Event', '2Lep_PurePromptEvent', 'TauVeto', 'TrigMatch', 'LowJetCR']) ) )
	 #
	 # OF region
	 #
	 vardb.registerCategory( MyCategory('OFSSLowNJetCR_AllCuts',    cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'OF_Event', '2Lep_PurePromptEvent', 'TauVeto', 'TrigMatch', 'OF_ElEtaCut', 'LowJetCR']) ) )
	 #
	 # ElEl region
	 #
	 vardb.registerCategory( MyCategory('ElElSSLowNJetCR_AllCuts',    cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'ElEl_Event', '2Lep_PurePromptEvent', 'TauVeto', 'TrigMatch', 'SF_ElEtaCut', 'LowJetCR']) ) )
    else:
    	 vardb.registerCategory( MyCategory('MuMuSSLowNJetCR',	cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'MuMu_Event', 'TauVeto', 'TrigMatch', 'LowJetCR']) ) )
    	 vardb.registerCategory( MyCategory('OFSSLowNJetCR',	cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'OF_Event', 'TauVeto', 'TrigMatch', 'OF_ElEtaCut', 'LowJetCR']) ) )
    	 vardb.registerCategory( MyCategory('ElElSSLowNJetCR',	cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'ElEl_Event', 'TauVeto', 'TrigMatch', 'SF_ElEtaCut', 'LowJetCR']) ) )
	 vardb.registerCategory( MyCategory('TwoLepSSTauLowNJetCR', cut = vardb.getCuts(['NBJet', '2LepTau', 'SS', 'LowJetCR', 'OneTau', 'TrigMatch', 'Zsidescut']) ) )

if doThreeLepLowNJetCR:
    #
    # take OS pairs
    #
    vardb.registerCategory( MyCategory('ThreeLepLowNJetCR',   cut = ( vardb.getCuts(['NBJet', '3Lep', 'TrigMatch', 'LowJetCR','ZOSsidescut']) & -vardb.getCut('SS') ) ) )  

# -------------
# other CRs
# -------------

if doWZonCR:
    vardb.registerCategory( MyCategory('WZonCR',      cut = ( vardb.getCuts(['BJetVeto',   '3Lep',  		 'TrigMatch',  'Zpeakcut'])  & -vardb.getCut('SS') ) ) ) 

if doWZoffCR:
    vardb.registerCategory( MyCategory('WZoffCR',     cut = ( vardb.getCuts(['BJetVeto',   '3Lep',  		 'TrigMatch',  'Zsidescut']) & -vardb.getCut('SS') ) ) ) 

if doWZHFonCR:
    vardb.registerCategory( MyCategory('WZHFonCR',    cut = ( vardb.getCuts(['NBJet',      '3Lep',  		 'TrigMatch',  'Zpeakcut'])  & -vardb.getCut('SS') ) ) )

if doWZHFoffCR:
    vardb.registerCategory( MyCategory('WZHFoffCR',   cut = ( vardb.getCuts(['NBJet',      '3Lep',  		 'TrigMatch',  'Zsidescut']) & -vardb.getCut('SS') ) ) )

if dottZCR:
    vardb.registerCategory( MyCategory('ttZCR',       cut = ( vardb.getCuts(['NBJet',      '3Lep',   'NJet3L',   'TrigMatch',  'Zpeakcut'])  & -vardb.getCut('SS') ) ) ) 

if dottWCR:
    vardb.registerCategory( MyCategory('ttWCR',       cut =   vardb.getCuts(['LargeNBjet', '2Lep',  'LowJetCR',  'TrigMatch', 'TauVeto', 'SS']) ) )

if doZSSpeakCR:
    vardb.registerCategory( MyCategory('ZSSpeakCR_ElEl',   cut = vardb.getCuts(['2Lep', 'SS', 'ElEl_Event', 'Zpeakcut',  'TrigMatch', 'SF_ElEtaCut']) ) )
    vardb.registerCategory( MyCategory('ZSSpeakCR_MuMu',   cut = vardb.getCuts(['2Lep', 'SS', 'MuMu_Event', 'Zpeakcut',  'TrigMatch']) ) )

# ------------------------------------
# Special CR for Data/MC control plots 
# ------------------------------------

vardb.registerCut( Cut('Mu0Tight',    '( muon_isIsolated_Tight[0] == 1 && muon_trkd0sig[0] < 10.0 && TMath::Abs(muon_trkz0sintheta[0]) < 2.0 )') ) 
vardb.registerCut( Cut('Mu1Tight',    '( muon_isIsolated_Tight[1] == 1 && muon_trkd0sig[1] < 10.0 && TMath::Abs(muon_trkz0sintheta[1]) < 2.0 )') ) 
vardb.registerCut( Cut('El0Tight',    '( el_isIsolated_Tight[0] == 1 && el_LHTight[0] == 1 && el_trkd0sig[0] < 10 && TMath::Abs(el_trkz0sintheta[0]) < 2.0 )') )  
vardb.registerCut( Cut('El1Tight',    '( el_isIsolated_Tight[1] == 1 && el_LHTight[1] == 1 && el_trkd0sig[1] < 10 && TMath::Abs(el_trkz0sintheta[1]) < 2.0 )') )  

# Run1 isolation
#
vardb.registerCut( Cut('Mu0Tight_Run1Iso',    '( muon_isIsolated_UserDefinedCut[0] == 1 && muon_trkd0sig[0] < 10.0 && TMath::Abs(muon_trkz0sintheta[0]) < 2.0 )') ) 
vardb.registerCut( Cut('Mu1Tight_Run1Iso',    '( muon_isIsolated_UserDefinedCut[1] == 1 && muon_trkd0sig[1] < 10.0 && TMath::Abs(muon_trkz0sintheta[1]) < 2.0 )') ) 
vardb.registerCut( Cut('El0Tight_Run1Iso',    '( el_isIsolated_UserDefinedCut[0] == 1 && el_LHTight[0] == 1 && el_trkd0sig[0] < 10 && TMath::Abs(el_trkz0sintheta[0]) < 2.0 )') ) 
vardb.registerCut( Cut('El1Tight_Run1Iso',    '( el_isIsolated_UserDefinedCut[1] == 1 && el_LHTight[1] == 1 && el_trkd0sig[1] < 10 && TMath::Abs(el_trkz0sintheta[1]) < 2.0 )') ) 

# no isolation
#
vardb.registerCut( Cut('Mu0Tight_NoIso',    '( muon_trkd0sig[0] < 10.0 && TMath::Abs(muon_trkz0sintheta[0]) < 2.0 )') ) 
vardb.registerCut( Cut('Mu1Tight_NoIso',    '( muon_trkd0sig[1] < 10.0 && TMath::Abs(muon_trkz0sintheta[1]) < 2.0 )') ) 
vardb.registerCut( Cut('El0Tight_NoIso',    '( el_LHTight[0] == 1 && el_trkd0sig[0] < 10 && TMath::Abs(el_trkz0sintheta[0]) < 2.0 )') )  
vardb.registerCut( Cut('El1Tight_NoIso',    '( el_LHTight[1] == 1 && el_trkd0sig[1] < 10 && TMath::Abs(el_trkz0sintheta[1]) < 2.0 )') )  


if doDataMCCR:
    #
    # Z+jets CR
    #
    # mumu
    #
    #vardb.registerCategory( MyCategory('DataMC_2Lep_TrigMatch_MuMu_Zmincut_OS_TT_NoIso',   cut = ( vardb.getCuts(['2Lep', 'TrigDec', 'TrigMatch', 'MuMu_Event',        'Zmincut',  'Mu0Tight_NoIso', 'Mu1Tight_NoIso']) & -vardb.getCut('SS') ) ) )   
    #vardb.registerCategory( MyCategory('DataMC_2Lep_TrigMatch_MuMu_Zpeakcut_OS_TT_NoIso',  cut = ( vardb.getCuts(['2Lep', 'TrigDec', 'TrigMatch', 'MuMu_Event',        'Zpeakcut', 'Mu0Tight_NoIso', 'Mu1Tight_NoIso']) & -vardb.getCut('SS') ) ) )   
    
    #vardb.registerCategory( MyCategory('DataMC_2Lep_TrigMatch_MuMu_Zmincut_TT',	   cut =   vardb.getCuts(['2Lep', 'TrigDec', 'TrigMatch', 'MuMu_Event',        'Zmincut',  'Mu0Tight', 'Mu1Tight']) ) )    
    vardb.registerCategory( MyCategory('DataMC_2Lep_TrigMatch_MuMu_Zmincut_OS_TT',   cut = ( vardb.getCuts(['2Lep', 'TrigDec', 'TrigMatch', 'MuMu_Event',        'Zmincut', 'Mu0Tight', 'Mu1Tight']) & -vardb.getCut('SS') ) ) ) 
    #vardb.registerCategory( MyCategory('DataMC_2Lep_TrigMatch_MuMu_Zpeak_OS_TT',   cut = ( vardb.getCuts(['2Lep', 'TrigDec', 'TrigMatch', 'MuMu_Event',        'Zpeakcut', 'Mu0Tight', 'Mu1Tight']) & -vardb.getCut('SS') ) ) ) 
    #vardb.registerCategory( MyCategory('DataMC_2Lep_TrigMatch_MuMu_Zpeak_JetBarrel_OS_TT',   cut = ( vardb.getCuts(['2Lep', 'TrigDec', 'TrigMatch', 'MuMu_Event',        'Jet0FwdCut', 'Zpeakcut', 'Mu0Tight', 'Mu1Tight']) & -vardb.getCut('SS') ) ) ) 
    #vardb.registerCategory( MyCategory('DataMC_2Lep_TrigMatch_MuMu_Zmincut_OS_TT_Run1Iso',   cut = ( vardb.getCuts(['2Lep', 'TrigDec', 'TrigMatch', 'MuMu_Event',        'Zmincut', 'Mu0Tight_Run1Iso', 'Mu1Tight_Run1Iso']) & -vardb.getCut('SS') ) ) )   
    #vardb.registerCategory( MyCategory('DataMC_2Lep_TrigMatch_MuMu_Zmincut_SS_TT',   cut =   vardb.getCuts(['2Lep', 'TrigDec', 'TrigMatch', 'MuMu_Event', 'SS', 'Zmincut', 'Mu0Tight', 'Mu1Tight']) ) )   
    #
    # elel
    #
    #vardb.registerCategory( MyCategory('DataMC_2Lep_TrigMatch_ElEl_Zmincut_OS_TT_NonIso',   cut = ( vardb.getCuts(['2Lep', 'TrigDec', 'TrigMatch', 'ElEl_Event',        'Zmincut', 'El0Tight_NoIso', 'El1Tight_NoIso']) & -vardb.getCut('SS') ) ) )   
    #vardb.registerCategory( MyCategory('DataMC_2Lep_TrigMatch_ElEl_Zpeakcut_OS_TT_NonIso',  cut = ( vardb.getCuts(['2Lep', 'TrigDec', 'TrigMatch', 'ElEl_Event',        'Zpeakcut','El0Tight_NoIso', 'El1Tight_NoIso']) & -vardb.getCut('SS') ) ) )   
    
    #vardb.registerCategory( MyCategory('DataMC_2Lep_TrigMatch_ElEl_Zmincut_TT', 	   cut =   vardb.getCuts(['2Lep', 'TrigDec', 'TrigMatch', 'ElEl_Event',         'Zmincut','El0Tight', 'El1Tight', 'SF_ElEtaCut']) ) )  
    vardb.registerCategory( MyCategory('DataMC_2Lep_TrigMatch_ElEl_Zmincut_OS_TT',   cut = ( vardb.getCuts(['2Lep', 'TrigDec', 'TrigMatch', 'ElEl_Event',        'Zmincut','El0Tight', 'El1Tight']) & -vardb.getCut('SS') ) ) )   
    #vardb.registerCategory( MyCategory('DataMC_2Lep_TrigMatch_ElEl_Zpeak_OS_TT',   cut = ( vardb.getCuts(['2Lep', 'TrigDec', 'TrigMatch', 'ElEl_Event',        'Zpeakcut','El0Tight', 'El1Tight']) & -vardb.getCut('SS') ) ) )   
    #vardb.registerCategory( MyCategory('DataMC_2Lep_TrigMatch_ElEl_Zmincut_OS_TT_Run1Iso',   cut = ( vardb.getCuts(['2Lep', 'TrigDec', 'TrigMatch', 'ElEl_Event',        'Zmincut','El0Tight_Run1Iso', 'El1Tight_Run1Iso']) & -vardb.getCut('SS') ) ) )   
    #vardb.registerCategory( MyCategory('DataMC_2Lep_TrigMatch_ElEl_Zmincut_SS_TT',   cut =   vardb.getCuts(['2Lep', 'TrigDec', 'TrigMatch', 'ElEl_Event', 'SS',  'Zmincut','El0Tight', 'El1Tight']) ) )   
    #
    # ttbar, Z tau tau (2Lep, OF) CR
    #    
    #vardb.registerCategory( MyCategory('DataMC_2Lep_TrigMatch_OF_Zmincut_TT_OS', 	   cut = vardb.getCuts(['2Lep', 'TrigDec', 'TrigMatch', 'OF_Event',       'Zmincut','Mu0Tight', 'El0Tight']) & -vardb.getCut('SS') ) ) 
    #vardb.registerCategory( MyCategory('DataMC_2Lep_TrigMatch_MuEl_Zmincut_TT_OS',	  cut = vardb.getCuts(['2Lep', 'TrigDec', 'TrigMatch', 'MuEl_Event',      'Zmincut','Mu0Tight', 'El0Tight']) & -vardb.getCut('SS') ) ) 
    #vardb.registerCategory( MyCategory('DataMC_2Lep_TrigMatch_ElMu_Zmincut_TT_OS',	  cut = vardb.getCuts(['2Lep', 'TrigDec', 'TrigMatch', 'ElMu_Event',      'Zmincut','Mu0Tight', 'El0Tight']) & -vardb.getCut('SS') ) ) 
    #vardb.registerCategory( MyCategory('DataMC_2Lep_TrigMatch_OF_Zmincut_TT_NBJets_OS',  cut = ( vardb.getCuts(['2Lep', 'TrigDec', 'TrigMatch', 'OF_Event',     'Zmincut','Mu0Tight', 'El0Tight', 'NBJet']) & -vardb.getCut('SS') ) ) )  
    #vardb.registerCategory( MyCategory('DataMC_2Lep_TrigMatch_MuEl_Zmincut_TT_NBJets_OS',  cut = ( vardb.getCuts(['2Lep', 'TrigDec', 'TrigMatch', 'MuEl_Event', 'Zmincut','Mu0Tight', 'El0Tight', 'NBJet']) & -vardb.getCut('SS') ) ) )  
    #vardb.registerCategory( MyCategory('DataMC_2Lep_TrigMatch_ElMu_Zmincut_TT_NBJets_OS',  cut = ( vardb.getCuts(['2Lep', 'TrigDec', 'TrigMatch', 'ElMu_Event', 'Zmincut','Mu0Tight', 'El0Tight', 'NBJet']) & -vardb.getCut('SS') ) ) )  
    

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
if doElRFRateCR :
    vardb.registerCut( Cut('ElRealFakeRateCR',	'( isProbeElEvent == 1 )') ) 
    vardb.registerCut( Cut('ElProbeTight',      '( el_probe_isTightSelected[0] == 1 )') )
    vardb.registerCut( Cut('ElProbeAntiTight',  '( el_probe_isTightSelected[0] == 0 )') )
if doMuRFRateCR :
    vardb.registerCut( Cut('MuRealFakeRateCR',	'( isProbeMuEvent == 1 )') ) 
    vardb.registerCut( Cut('MuProbeTight',      '( muon_probe_isTightSelected[0] == 1 )') )
    vardb.registerCut( Cut('MuProbeAntiTight',  '( muon_probe_isTightSelected[0] == 0 )') )
    
#
# define regions with at least one muon/electron
#
mu_region = ( vardb.getCut('MuMu_Event') | vardb.getCut('OF_Event') )
el_region = ( vardb.getCut('ElEl_Event') | vardb.getCut('OF_Event') )   

# ---------------------------------------
# Special plots for MM real/fake rate CRs 
# ---------------------------------------
if doElRFRateCR or doMMClosureRates:
    print ''
    #vardb.registerVar( Variable(shortname = 'ElTagPt',    latexname = 'p_{T}^{tag e} [GeV]',     ntuplename = 'el_tag_pt[0]/1e3',                bins = 36, minval = 10.0, maxval = 100.0,) )
    #vardb.registerVar( Variable(shortname = 'ElTagEta',   latexname = '#eta^{tag e}',            ntuplename = 'TMath::Abs( el_tag_eta[0] )',     bins = 8,  minval = 0.0,  maxval = 2.6, manualbins = [ 0.0 , 0.5 , 0.8 , 1.1 , 1.37 , 1.52 , 2.0 , 2.25 , 2.6],) )
    vardb.registerVar( Variable(shortname = 'ElProbePt',  latexname = 'p_{T}^{probe e} [GeV]',   ntuplename = 'el_probe_pt[0]/1e3',		 bins = 100, minval = 10.0, maxval = 310.0,) )
    vardb.registerVar( Variable(shortname = 'ElProbeEta', latexname = '#eta^{probe e}', 	 ntuplename = 'TMath::Abs( el_probe_eta[0] )',   bins = 8,  minval = 0.0,  maxval = 2.6, manualbins = [ 0.0 , 0.5 , 0.8 , 1.1 , 1.37 , 1.52 , 2.0 , 2.25 , 2.6],) )
    #
    # large bin for high eta
    #
    #vardb.registerVar( Variable(shortname = 'ElProbeEta', latexname = '#eta^{probe e}',          ntuplename = 'TMath::Abs( el_probe_eta[0] )',   bins = 8,  minval = 0.0,  maxval = 2.6, manualbins = [ 0.0 , 0.5 , 0.8 , 1.1 , 1.37 , 1.52 , 2.6],) )


if doMuRFRateCR or doMMClosureRates:
    print ''
    #vardb.registerVar( Variable(shortname = 'MuTagPt',    latexname = 'p_{T}^{tag #mu} [GeV]',   ntuplename = 'muon_tag_pt[0]/1e3',              bins = 36, minval = 10.0, maxval = 100.0,) )
    #vardb.registerVar( Variable(shortname = 'MuTagEta',   latexname = '#eta^{tag #mu}',          ntuplename = 'TMath::Abs( muon_tag_eta[0] )',   bins = 8,  minval = 0.0,  maxval = 2.6, manualbins = [ 0.0 , 0.5 , 0.8 , 1.1 , 1.37 , 1.52 , 2.0 , 2.25 , 2.6],) )
    vardb.registerVar( Variable(shortname = 'MuProbePt',  latexname = 'p_{T}^{probe #mu} [GeV]', ntuplename = 'muon_probe_pt[0]/1e3',		bins = 100, minval = 10.0, maxval = 310.0,) )
    vardb.registerVar( Variable(shortname = 'MuProbeEta', latexname = '#eta^{probe #mu}',	ntuplename = 'TMath::Abs( muon_probe_eta[0] )', bins = 8,  minval = 0.0,  maxval = 2.6, manualbins = [ 0.0 , 0.5 , 0.8 , 1.1 , 1.37 , 1.52 , 2.0 , 2.25 , 2.6],) )
    #
    # large bin for high eta
    #
    #vardb.registerVar( Variable(shortname = 'MuProbeEta', latexname = '#eta^{probe #mu}',	 ntuplename = 'TMath::Abs( muon_probe_eta[0] )', bins = 8,  minval = 0.0,  maxval = 2.6, manualbins = [ 0.0 , 0.5 , 0.8 , 1.1 , 1.37 , 1.52 , 2.6],) )

#
# do the non-non-prompt subtraction only in Fake CR:
# ---> i.e., plot ONLY the !(non-prompt) contamination
#
    
if doMuRFRateCR:
    #
    # muon R/F region(s)
    #    
    vardb.registerCategory( MyCategory('FakeCRMuL',	 cut = ( vardb.getCuts(['TrigMatch', 'NBJet', '2Lep', 'SS', 'LowJetCR', 'MuRealFakeRateCR', '2Lep_NonNonPromptEvent', 'MuProbeAntiTight']) & mu_region ) ) )
    vardb.registerCategory( MyCategory('FakeCRMuT',	 cut = ( vardb.getCuts(['TrigMatch', 'NBJet', '2Lep', 'SS', 'LowJetCR', 'MuRealFakeRateCR', '2Lep_NonNonPromptEvent', 'MuProbeTight'    ]) & mu_region ) ) )
    vardb.registerCategory( MyCategory('RealCRMuL',	 cut = ( vardb.getCuts(['TrigMatch', 'NBJet', '2Lep',       'LowJetCR', 'MuRealFakeRateCR',			      'MuProbeAntiTight']) & mu_region & -vardb.getCut('SS') ) ) )
    vardb.registerCategory( MyCategory('RealCRMuT',	 cut = ( vardb.getCuts(['TrigMatch', 'NBJet', '2Lep',       'LowJetCR', 'MuRealFakeRateCR',			      'MuProbeTight'    ]) & mu_region & -vardb.getCut('SS') ) ) )
    #
    # split into MuMu, OF contributions
    #
    """
    vardb.registerCategory( MyCategory('MuMuFakeCRMuL',  cut =   vardb.getCuts(['TrigMatch', 'NBJet', '2Lep', 'SS', 'LowJetCR', 'MuRealFakeRateCR', '2Lep_NonNonPromptEvent',  'MuProbeAntiTight', 'MuMu_Event']) ) )
    vardb.registerCategory( MyCategory('OFFakeCRMuL',  cut =   vardb.getCuts(['TrigMatch', 'NBJet', '2Lep', 'SS', 'LowJetCR', 'MuRealFakeRateCR', '2Lep_NonNonPromptEvent',  'MuProbeAntiTight', 'OF_Event']) ) )
    vardb.registerCategory( MyCategory('MuMuFakeCRMuT',  cut =   vardb.getCuts(['TrigMatch', 'NBJet', '2Lep', 'SS', 'LowJetCR', 'MuRealFakeRateCR', '2Lep_NonNonPromptEvent',  'MuProbeTight',     'MuMu_Event']) ) )
    vardb.registerCategory( MyCategory('OFFakeCRMuT',  cut =   vardb.getCuts(['TrigMatch', 'NBJet', '2Lep', 'SS', 'LowJetCR', 'MuRealFakeRateCR', '2Lep_NonNonPromptEvent',  'MuProbeTight',     'OF_Event']) ) )
    vardb.registerCategory( MyCategory('MuMuRealCRMuL',  cut = ( vardb.getCuts(['TrigMatch', 'NBJet', '2Lep',	    'LowJetCR', 'MuRealFakeRateCR',			       'MuProbeAntiTight', 'MuMu_Event']) & -vardb.getCut('SS') ) ) )
    vardb.registerCategory( MyCategory('OFRealCRMuL',  cut = ( vardb.getCuts(['TrigMatch', 'NBJet', '2Lep',	    'LowJetCR', 'MuRealFakeRateCR',			       'MuProbeAntiTight', 'OF_Event']) & -vardb.getCut('SS') ) ) )
    vardb.registerCategory( MyCategory('MuMuRealCRMuT',  cut = ( vardb.getCuts(['TrigMatch', 'NBJet', '2Lep',	    'LowJetCR', 'MuRealFakeRateCR',			       'MuProbeTight',	   'MuMu_Event']) & -vardb.getCut('SS') ) ) )
    vardb.registerCategory( MyCategory('OFRealCRMuT',  cut = ( vardb.getCuts(['TrigMatch', 'NBJet', '2Lep',	    'LowJetCR', 'MuRealFakeRateCR',			       'MuProbeTight',	   'OF_Event']) & -vardb.getCut('SS') ) ) )
    """
if doElRFRateCR:
    #
    # electron R/F region(s)
    #
    # NB: for FAKE region, additionally require the SS leptons to be outside Z peak 
    # to reduce charge flip contamination (in principle we want to have only non-prompt here!) 
    #
    vardb.registerCategory( MyCategory('FakeCRElL',	 cut = ( vardb.getCuts(['TrigMatch', 'NBJet', '2Lep', 'SS', 'LowJetCR', 'ElRealFakeRateCR', '2Lep_NonNonPromptEvent',  'ElProbeAntiTight', 'Zsidescut']) & el_region ) ) )    
    vardb.registerCategory( MyCategory('FakeCRElT',	 cut = ( vardb.getCuts(['TrigMatch', 'NBJet', '2Lep', 'SS', 'LowJetCR', 'ElRealFakeRateCR', '2Lep_NonNonPromptEvent',  'ElProbeTight',     'Zsidescut']) & el_region ) ) )
    vardb.registerCategory( MyCategory('RealCRElL',	 cut = ( vardb.getCuts(['TrigMatch', 'NBJet', '2Lep',	    'LowJetCR', 'ElRealFakeRateCR',			       'ElProbeAntiTight'             ])  & el_region & -vardb.getCut('SS') ) ) )
    vardb.registerCategory( MyCategory('RealCRElT',	 cut = ( vardb.getCuts(['TrigMatch', 'NBJet', '2Lep',	    'LowJetCR', 'ElRealFakeRateCR',			       'ElProbeTight'                 ])  & el_region & -vardb.getCut('SS') ) ) )
    #
    # split into ElEl, OF contributions
    #
    """
    vardb.registerCategory( MyCategory('ElElFakeCRElL',  cut =   vardb.getCuts(['TrigMatch', 'NBJet', '2Lep', 'SS', 'LowJetCR', 'ElRealFakeRateCR', '2Lep_NonNonPromptEvent',  'ElProbeAntiTight', 'ElEl_Event', 'Zsidescut']) ) )
    vardb.registerCategory( MyCategory('OFFakeCRElL',  cut =   vardb.getCuts(['TrigMatch', 'NBJet', '2Lep', 'SS', 'LowJetCR', 'ElRealFakeRateCR', '2Lep_NonNonPromptEvent',  'ElProbeAntiTight', 'OF_Event', 'Zsidescut']) ) )
    vardb.registerCategory( MyCategory('ElElFakeCRElT',  cut =   vardb.getCuts(['TrigMatch', 'NBJet', '2Lep', 'SS', 'LowJetCR', 'ElRealFakeRateCR', '2Lep_NonNonPromptEvent',  'ElProbeTight',     'ElEl_Event', 'Zsidescut']) ) )
    vardb.registerCategory( MyCategory('OFFakeCRElT',  cut =   vardb.getCuts(['TrigMatch', 'NBJet', '2Lep', 'SS', 'LowJetCR', 'ElRealFakeRateCR', '2Lep_NonNonPromptEvent',  'ElProbeTight',     'OF_Event', 'Zsidescut']) ) )
    vardb.registerCategory( MyCategory('ElElRealCRElL',  cut = ( vardb.getCuts(['TrigMatch', 'NBJet', '2Lep',	    'LowJetCR', 'ElRealFakeRateCR',			       'ElProbeAntiTight', 'ElEl_Event'      	    ]) & -vardb.getCut('SS') ) ) )
    vardb.registerCategory( MyCategory('OFRealCRElL',  cut = ( vardb.getCuts(['TrigMatch', 'NBJet', '2Lep',	    'LowJetCR', 'ElRealFakeRateCR',			       'ElProbeAntiTight', 'OF_Event'      	    ]) & -vardb.getCut('SS') ) ) )
    vardb.registerCategory( MyCategory('ElElRealCRElT',  cut = ( vardb.getCuts(['TrigMatch', 'NBJet', '2Lep',	    'LowJetCR', 'ElRealFakeRateCR',			       'ElProbeTight',     'ElEl_Event'      	    ]) & -vardb.getCut('SS') ) ) )
    vardb.registerCategory( MyCategory('OFRealCRElT',  cut = ( vardb.getCuts(['TrigMatch', 'NBJet', '2Lep',	    'LowJetCR', 'ElRealFakeRateCR',			       'ElProbeTight',     'OF_Event'      	    ]) & -vardb.getCut('SS') ) ) )
    """
    
# MM CLOSURE RATES: 
#
# use MC ttbar as data 
#
# --> FAKE region: SS leptons, at least one lepton must be !(prompt), and none charge flip (see 2. truth cut above)
# --> REAL region: OS leptons, all leptons must be prompt (and none charge flip)
    
if doMMClosureRates :
    #"""
    vardb.registerCut( Cut('ElRealFakeRateCR',	'( isProbeElEvent == 1 )') ) 
    #
    # Run1 isolation, LHTight
    #    
    vardb.registerCut( Cut('ElProbeTight',	'( el_probe_LHTight[0] == 1 && el_probe_isIsolated_UserDefinedCut[0] == 1 )') )
    vardb.registerCut( Cut('ElProbeAntiTight',  '( el_probe_LHTight[0] == 0 || el_probe_isIsolated_UserDefinedCut[0] == 0 )') )	 
    
    vardb.registerCut( Cut('MuRealFakeRateCR',	'( isProbeMuEvent == 1 )') ) 
    #
    # Run1 isolation, |d0sig| 
    #    
    vardb.registerCut( Cut('MuProbeTight',      '( muon_probe_trkd0sig[0] < 3.0  && muon_probe_isIsolated_UserDefinedCut[0] == 1 )') )
    vardb.registerCut( Cut('MuProbeAntiTight',  '( muon_probe_trkd0sig[0] >= 3.0 || muon_probe_isIsolated_UserDefinedCut[0] == 0 )') )
    #"""
    """    
    vardb.registerCut( Cut('ElRealFakeRateCR',	'( isProbeElEvent == 1 )') ) 
    #
    # Run1 isolated, non-isolated
    #    
    vardb.registerCut( Cut('ElProbeTight',	'( el_probe_LHTight[0] == 1 && el_probe_isIsolated_UserDefinedCut[0] == 1 )') )
    vardb.registerCut( Cut('ElProbeAntiTight',  '( el_probe_LHTight[0] == 1 && el_probe_isIsolated_UserDefinedCut[0] == 0 )') )	 
    
    vardb.registerCut( Cut('MuRealFakeRateCR',	'( isProbeMuEvent == 1 )') ) 
    #
    # Run1 isolated,non-isolated
    #    
    vardb.registerCut( Cut('MuProbeTight',      '( muon_probe_trkd0sig[0] < 3.0 && muon_probe_isIsolated_UserDefinedCut[0] == 1 )') )
    vardb.registerCut( Cut('MuProbeAntiTight',  '( muon_probe_trkd0sig[0] < 3.0 && muon_probe_isIsolated_UserDefinedCut[0] == 0 )') )
    """    
    #
    # TRUTH composition of event in R/F region
    # ( NB: in fake region, this choice automatically tells whether the probe will be a HF non-prompt lepton or a charge flip )
    #
    truth_event_cut_R = vardb.getCut('2Lep_PurePromptEvent')      
    if ( args.doChFlipRate ):
        print '*********************************\nMEASURING CHARGE FLIP RATE IN MC\n*********************************'
        truth_event_cut_F = vardb.getCut('2Lep_ChFlipEvent')     # --> probe will be automatically a charge flip
    else :
        truth_event_cut_F = vardb.getCut('2Lep_NonPromptEvent')  # --> probe will be automatically a !prompt and not charge flip (DEFAULT)
    
if doMMClosureRates:
    print ''
    #
    # muon R/F region(s)
    #
    
#    vardb.registerCategory( MyCategory('FakeCRMu',    cut =  ( vardb.getCuts(['TrigMatch', 'NBJet', 'SS', 'LowJetCR', 'MuRealFakeRateCR']) & mu_region & truth_event_cut_F ) ) )  
#    vardb.registerCategory( MyCategory('RealCRMu',    cut =  ( vardb.getCuts(['TrigMatch', 'NBJet',	  'LowJetCR', 'MuRealFakeRateCR']) & mu_region & truth_event_cut_R & -vardb.getCut('SS') ) ) )  
    
    vardb.registerCategory( MyCategory('FakeCRMuL',    cut =  ( vardb.getCuts(['TrigMatch', 'NBJet', 'SS', 'LowJetCR', 'MuRealFakeRateCR',  'MuProbeAntiTight']) & mu_region & truth_event_cut_F ) ) )  
    vardb.registerCategory( MyCategory('FakeCRMuT',    cut =  ( vardb.getCuts(['TrigMatch', 'NBJet', 'SS', 'LowJetCR', 'MuRealFakeRateCR',  'MuProbeTight'    ]) & mu_region & truth_event_cut_F ) ) )  
    vardb.registerCategory( MyCategory('RealCRMuL',    cut =  ( vardb.getCuts(['TrigMatch', 'NBJet',	  'LowJetCR', 'MuRealFakeRateCR',  'MuProbeAntiTight']) & mu_region & truth_event_cut_R & -vardb.getCut('SS') ) ) )  
    vardb.registerCategory( MyCategory('RealCRMuT',    cut =  ( vardb.getCuts(['TrigMatch', 'NBJet',	  'LowJetCR', 'MuRealFakeRateCR',  'MuProbeTight'    ]) & mu_region & truth_event_cut_R & -vardb.getCut('SS') ) ) )  
    #
    # electron R/F region(s)
    #

#    vardb.registerCategory( MyCategory('FakeCREl',    cut =  ( vardb.getCuts(['TrigMatch', 'NBJet', 'SS', 'LowJetCR', 'ElRealFakeRateCR']) & el_region & truth_event_cut_F ) ) )  
#    vardb.registerCategory( MyCategory('RealCREl',    cut =  ( vardb.getCuts(['TrigMatch', 'NBJet',	  'LowJetCR', 'ElRealFakeRateCR']) & el_region & truth_event_cut_R & -vardb.getCut('SS') ) ) )  

    vardb.registerCategory( MyCategory('FakeCRElL',    cut =  ( vardb.getCuts(['TrigMatch', 'NBJet', 'SS', 'LowJetCR', 'ElRealFakeRateCR',  'ElProbeAntiTight']) & el_region & truth_event_cut_F ) ) )  
    vardb.registerCategory( MyCategory('FakeCRElT',    cut =  ( vardb.getCuts(['TrigMatch', 'NBJet', 'SS', 'LowJetCR', 'ElRealFakeRateCR',  'ElProbeTight'    ]) & el_region & truth_event_cut_F ) ) )  
    vardb.registerCategory( MyCategory('RealCRElL',    cut =  ( vardb.getCuts(['TrigMatch', 'NBJet',	   'LowJetCR', 'ElRealFakeRateCR',  'ElProbeAntiTight']) & el_region & truth_event_cut_R & -vardb.getCut('SS') ) ) )  
    vardb.registerCategory( MyCategory('RealCRElT',    cut =  ( vardb.getCuts(['TrigMatch', 'NBJet',	   'LowJetCR', 'ElRealFakeRateCR',  'ElProbeTight'    ]) & el_region & truth_event_cut_R & -vardb.getCut('SS') ) ) )  
    
    
if doMMClosureTest:
    print ''   
    
    if ( args.fakeMethod == 'MM' or args.fakeMethod == 'FF' ):
        #
        # MuMu region
        #
        vardb.registerCategory( MyCategory('MuMuHighJetTrigMatchNoTauTruthMSSBase',         cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'MuMu_Event', 'TrigMatch', 'NJet2L']) ) )
        #vardb.registerCategory( MyCategory('MuMuLowJetTrigMatchNoTauTruthMSSBase',          cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'MuMu_Event', 'TrigMatch', 'LowJetCR']) ) )
        #vardb.registerCategory( MyCategory('MuMuAllJetTrigMatchNoTauTruthMSSBase',          cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'MuMu_Event', 'TrigMatch']) ) )
        #
	# OF region
	#
	vardb.registerCategory( MyCategory('OFHighJetEtaCutTrigMatchNoTauTruthMSSBase',   cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'OF_Event', 'TrigMatch', 'NJet2L']) ) )   # 'OF_ElEtaCut', 
	#vardb.registerCategory( MyCategory('OFLowJetEtaCutTrigMatchNoTauTruthMSSBase',    cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'OF_Event', 'TrigMatch', 'LowJetCR']) ) ) # 'OF_ElEtaCut', 
	#vardb.registerCategory( MyCategory('OFAllJetEtaCutTrigMatchNoTauTruthMSSBase',    cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'OF_Event', 'TrigMatch']) ) ) # 'OF_ElEtaCut', 
	#
	# ElEl region
	#
	vardb.registerCategory( MyCategory('ElElHighJetEtaCutTrigMatchNoTauTruthMSSBase',   cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'ElEl_Event', 'TrigMatch', 'NJet2L']) ) )   # 'SF_ElEtaCut', 
	#vardb.registerCategory( MyCategory('ElElLowJetEtaCutTrigMatchNoTauTruthMSSBase',    cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'ElEl_Event', 'TrigMatch', 'LowJetCR']) ) ) # 'SF_ElEtaCut', 
        #vardb.registerCategory( MyCategory('ElElAllJetEtaCutTrigMatchNoTauTruthMSSBase',    cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'ElEl_Event', 'TrigMatch']) ) ) # 'SF_ElEtaCut', 
    
    if ( args.fakeMethod == 'ABCD' ):  
        #
        # MuMu region
        #
        vardb.registerCategory( MyCategory('MuMuHighJetTrigMatchNoTauTruthMSSBase',         cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'MuMu_Event', 'TrigMatch', 'NJet2L']) ) )
        #
	# OF region
	#
	vardb.registerCategory( MyCategory('OFHighJetEtaCutTrigMatchNoTauTruthMSSBase',     cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'OF_Event', 'TrigMatch', 'NJet2L']) ) )    # 'OF_ElEtaCut', 
	#
	# ElEl region
	#
	vardb.registerCategory( MyCategory('ElElHighJetEtaCutTrigMatchNoTauTruthMSSBase',   cut = vardb.getCuts(['NBJet', '2Lep', 'SS', 'ElEl_Event', 'TrigMatch', 'NJet2L']) ) )  # 'SF_ElEtaCut', 
    
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

# period A3,A4,C2,C3,C4,C5 GRL (EPS)
#ttH2015.luminosity = 0.0849676 

# period D1-D6
#ttH2015.luminosity = 0.0803592

# period D-E, GRL v65
ttH2015.luminosity = 0.3224
#ttH2015.luminosity = 0.278979
ttH2015.lumi_units = 'pb-1'

# for MM closure
if doMMClosureTest or doMMClosureRates:
	ttH2015.luminosity = 5.4
	ttH2015.lumi_units = 'fb-1'
	
# --------------------
# set the event weight
# --------------------

if not ( doMMClosureTest or doMMClosureRates ):
	ttH2015.eventweight = 'weight_pileup'

print "\t Global eventweight (apply to ALL categories) - MC only --> ", ttH2015.eventweight

# ------------------------------------

ttH2015.useZCorrections = False

isblinded = False
#
# Make blinded plots in SR
#
if doSR and not doMMClosureTest:
        ttH2015.observed = []
        isblinded=True

if doTwoLepSR or doTwoLepLowNJetCR or dottWCR or doZSSpeakCR or doMMClosureTest :
    ttH2015.channel = 'TwoLepSS'
elif doThreeLepSR or doThreeLepLowNJetCR or dottZCR or doWZonCR or doWZoffCR or doWZHFonCR or doWZHFoffCR:
    ttH2015.channel = 'ThreeLep'
elif doFourLepSR:
    ttH2015.channel = 'FourLep'
elif doDataMCCR or doMuRFRateCR or doElRFRateCR or doMMClosureRates :
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
		'TTBarClosureMM':'ttbarbkg',		
		'TopCF':'topcfbkg', 
		'Diboson':'dibosonbkg', 
		'PowhegDiboson':'powhegdibosonbkg', 
		'PowhegDibosonWW':'powhegdibosonwwbkg', 
		'PowhegDibosonWZ':'powhegdibosonwzbkg', 
		'PowhegDibosonZZ':'powhegdibosonzzbkg', 
		'DibosonCF':'dibosoncfbkg', 
		'HtoZZ':'htozzbkg', 
		'Zjets':'zjetsbkg', 
		'Zeejets':'zeejetsbkg', 
		'Zmumujets':'zmumujetsbkg',  
		'Ztautaujets':'ztautaujetsbkg',
		'ZjetsLF':'zjetsbkg', 
		'SherpaZjets':'sherpazjetsbkg', 
		'SherpaZeejets':'sherpazeejetsbkg', 
		'SherpaZmumujets':'sherpazmumujetsbkg',  
		'SherpaZtautaujets':'sherpaztautaujetsbkg',	
		'SherpaZjetsBFilter':'sherpazjetsbfilter',
		'SherpaZjetsCFilterBVeto':'sherpazjetsbfiltercveto',	  
		'SherpaZjetsCVetoBVeto':'sherpazjetscvetobveto',		
		'SherpaZeejetsBFilter':'sherpazeejetsbfilter',
		'SherpaZeejetsCFilterBVeto':'sherpazeejetsbfiltercveto',	  
		'SherpaZeejetsCVetoBVeto':'sherpazjetscvetobveto',		
		'SherpaZmumujetsBFilter':'sherpazmumujetsbfilter',
		'SherpaZmumujetsCFilterBVeto':',sherpazmumujetsbfiltercveto',	  
		'SherpaZmumujetsCVetoBVeto':'sherpazmumujetscvetobveto',		
		'SherpaZtautaujetsBFilter':'sherpaztautaujetsbfilter',
		'SherpaZtautaujetsCFilterBVeto':',sherpaztautaujetsbfiltercveto',	  
		'SherpaZtautaujetsCVetoBVeto':',sherpaztautaujetscvetobveto',		
		'Wjets':'wjetsbkg', 
		'PowhegPythiaWjets':'powhegpythiawjets',
		'SherpaWjets':'sherpawjets',		
		'SherpaWenujets':'sherpawenujets',		
		'SherpaWmunujets':'sherpawmunujets',
		'SherpaWtaunujets':'sherpawtaunujets',		
		'SherpaWjetsBFilter':'sherpawjetsbfilter',
		'SherpaWjetsCFilterBVeto':'sherpawjetsbfiltercveto',
		'SherpaWjetsCVetoBVeto':'sherpawjetscvetobveto',
		'SherpaWenujetsBFilter':'sherpawenujetsbfilter',
		'SherpaWenujetsCFilterBVeto':'sherpawenujetsbfiltercveto',
		'SherpaWenujetsCVetoBVeto':'sherpawenujetscvetobveto',		
		'SherpaWmunujetsBFilter':'sherpawmunujetsbfilter',
		'SherpaWmunujetsCFilterBVeto':'sherpawmunujetsbfiltercveto',
		'SherpaWmunujetsCVetoBVeto':'sherpawmunujetscvetobveto',		
		'SherpaWtaunujetsBFilter':'sherpawtaunujetsbfilter',
		'SherpaWtaunujetsCFilterBVeto':'sherpawtaunujetsbfiltercveto',
		'SherpaWtaunujetsCVetoBVeto':'sherpawtaunujetscvetobveto',
		'Prompt':'promptbkg', 
		'ChargeFlip':'chargeflipbkg', 
		'FakesFF':'fakesbkg', 
		'FakesMM':'fakesbkg',
                'FakesClosureMM':'fakesbgk', 
                'FakesClosureABCD':'fakesbgk', 		
	      }
#	      
# Override colours!
#      
colours      = {'Observed':kBlack, 
        	'TTBarH':kBlack, 
        	'TTBarW':kRed-4, 
        	'TTBarZ':kRed-7, 
        	'Top':kAzure+1, 
        	'TTBar':kAzure+8,
		'TTBarClosureMM':kAzure+8,
        	'TopCF':kAzure-4, 
        	'Diboson':kYellow-9, 
		'PowhegDiboson':kYellow-9, 
		'PowhegDibosonWW':kYellow-7, 
		'PowhegDibosonWZ':kYellow-4, 
		'PowhegDibosonZZ':kYellow-3, 
        	'DibosonCF':kOrange-3, 
        	'HtoZZ':kTeal+9, 
        	'Zjets':kGreen, 
		'Zeejets':kGreen-7,
		'Zmumujets':kTeal+2, 
		'Ztautaujets':kTeal, 
        	'ZjetsLF':kGreen, 
		'SherpaZjets':kGreen, 
		'SherpaZeejets':kGreen-7, 
		'SherpaZmumujets':kTeal+2,  
		'SherpaZtautaujets':kTeal,	
		'SherpaZjetsBFilter':kGreen,
		'SherpaZjetsCFilterBVeto':kGreen,	  
		'SherpaZjetsCVetoBVeto':kGreen,		
		'SherpaZeejetsBFilter':kGreen-7,
		'SherpaZeejetsCFilterBVeto':kGreen-7,	  
		'SherpaZeejetsCVetoBVeto':kGreen-7,		
		'SherpaZmumujetsBFilter':kTeal+2,
		'SherpaZmumujetsCFilterBVeto':kTeal+2,	  
		'SherpaZmumujetsCVetoBVeto':kTeal+2,		
		'SherpaZtautaujetsBFilter':kTeal,
		'SherpaZtautaujetsCFilterBVeto':kTeal,	  
		'SherpaZtautaujetsCVetoBVeto':kTeal,		
        	'Wjets':kWhite, 
		'PowhegPythiaWjets':kWhite,
		'SherpaWjets':kWhite,	
		'SherpaWenujets':kGray,	
		'SherpaWmunujets':kGray+1,	
		'SherpaWtaunujets':kGray+2,	
		'SherpaWjetsBFilter':kWhite,
		'SherpaWjetsCFilterBVeto':kWhite,
		'SherpaWjetsCVetoBVeto':kWhite,
		'SherpaWenujetsBFilter':kGray,
		'SherpaWenujetsCFilterBVeto':kGray,
		'SherpaWenujetsCVetoBVeto':kGray,	  	  
		'SherpaWmunujetsBFilter':kGray+1,
		'SherpaWmunujetsCFilterBVeto':kGray+1,
		'SherpaWmunujetsCVetoBVeto':kGray+1,
		'SherpaWtaunujetsBFilter':kGray+2,
		'SherpaWtaunujetsCFilterBVeto':kGray+2,
		'SherpaWtaunujetsCVetoBVeto':kGray+2,
        	'Prompt':kOrange, 
        	'ChargeFlip':kAzure-4, 
        	'FakesFF':kAzure-9, 
        	'FakesMM':kTeal-9, 
                'FakesClosureMM':kTeal+1,
		'FakesClosureABCD':kCyan - 9, 
              }

if ( doSR or doLowNJetCR ):
    
    if not doFourLepSR: 
    
    	if doMM:
    	    plotbackgrounds	= [ 'TTBarW', 'TTBarZ', 'Top', 'TopCF', 'Diboson', 'DibosonCF', 'ZjetsLF', 'FakesMM']
    	    ttH2015.backgrounds = [ 'TTBarW', 'TTBarZ', 'Top', 'TopCF', 'Diboson', 'DibosonCF', 'ZjetsLF', 'FakesMM']
    	elif doFF:
    	    plotbackgrounds	= [ 'TTBarW', 'TTBarZ', 'Top', 'TopCF', 'Diboson', 'DibosonCF', 'ZjetsLF', 'FakesFF']
    	    ttH2015.backgrounds = [ 'TTBarW', 'TTBarZ', 'Top', 'TopCF', 'Diboson', 'DibosonCF', 'ZjetsLF', 'FakesFF']
    	else:
    	    # MC based estimate of fakes
    	    plotbackgrounds	= [ 'TTBarW', 'TTBarZ', 'Top', 'TopCF', 'Diboson', 'DibosonCF', 'ZjetsLF']
    	    ttH2015.backgrounds = [ 'TTBarW', 'TTBarZ', 'Top', 'TopCF', 'Diboson', 'DibosonCF', 'ZjetsLF'] 
    	
    else:
        # no fakes in 4lep
        plotbackgrounds	    = [ 'TTBarW', 'TTBarZ', 'Top', 'TopCF', 'Diboson', 'DibosonCF', 'ZjetsLF']
        ttH2015.backgrounds = [ 'TTBarW', 'TTBarZ', 'Top', 'TopCF', 'Diboson', 'DibosonCF', 'ZjetsLF'] 

if (doElRFRateCR or doMuRFRateCR):
    # here the fakes are the difference between data and MC
    #plotbackgrounds     = [ 'TTBarW', 'TTBarZ', 'Top', 'TopCF', 'Diboson', 'DibosonCF', 'ZjetsLF'] #, 'Wjets'] # the W+jets should be 0 in the truthmatching case and should be a rough estimate of fakes in the other case. Also ttbar should be bigger in the no truth matching case
    #ttH2015.backgrounds = [ 'TTBarW', 'TTBarZ', 'Top', 'TopCF', 'Diboson', 'DibosonCF', 'ZjetsLF'] #, 'Wjets'] # this is the list of backgrounds that will be calculated
    plotbackgrounds	= ['Top','TTBar','SherpaZeejets','SherpaZmumujets','SherpaZtautaujets','Diboson','SherpaWenujets','SherpaWmunujets','SherpaWtaunujets','TTBarW','TTBarZ']
    ttH2015.backgrounds = ['Top','TTBar','SherpaZeejets','SherpaZmumujets','SherpaZtautaujets','Diboson','SherpaWenujets','SherpaWmunujets','SherpaWtaunujets','TTBarW','TTBarZ'] 
    ttH2015.signals     = []

if doDataMCCR:
    #plotbackgrounds	 = ['Top','TTBar','Zeejets','Zmumujets','Ztautaujets','PowhegDibosonWW','PowhegDibosonWZ','PowhegDibosonZZ']
    #ttH2015.backgrounds = ['Top','TTBar','Zeejets','Zmumujets','Ztautaujets','PowhegDibosonWW','PowhegDibosonWZ','PowhegDibosonZZ'] 
     
    # PowhegPythia Z+jets
    # 
    #plotbackgrounds	= ['Top','TTBar','Zeejets','Zmumujets','Ztautaujets','Diboson','SherpaWjets','TTBarW','TTBarZ']
    #ttH2015.backgrounds = ['Top','TTBar','Zeejets','Zmumujets','Ztautaujets','Diboson','SherpaWjets','TTBarW','TTBarZ'] 
    
    # Sherpa Z+jets
    #     
    plotbackgrounds	= ['Top','TTBar','SherpaZeejets','SherpaZmumujets','SherpaZtautaujets','Diboson','SherpaWjets','TTBarW','TTBarZ']
    ttH2015.backgrounds = ['Top','TTBar','SherpaZeejets','SherpaZmumujets','SherpaZtautaujets','Diboson','SherpaWjets','TTBarW','TTBarZ']     
    
    #plotbackgrounds     = ['TTBar']
    #ttH2015.backgrounds = ['TTBar'] 
    
    ttH2015.signals     = []
    ttH2015.observed    = ['Observed']           
    

if doMMClosureRates:
      plotbackgrounds	  = ['TTBar'] 
      ttH2015.backgrounds = ['TTBar'] 
      ttH2015.signals	  = []
      ttH2015.observed    = []      
      #ttH2015.observed    = ['TTBar']  

if doMMClosureTest:
    if doMM:
      plotbackgrounds	  = ['FakesClosureMM']
      ttH2015.backgrounds = ['FakesClosureMM'] # truth cuts done internally in FakesClosureMM class
      ttH2015.signals     = ['FakesClosureABCD']
      ttH2015.observed    = ['TTBarClosureMM'] # truth cuts done internally in TTBarClosureMM class
      #ttH2015.observed    = ['TTBar']         
    elif doFF:
      plotbackgrounds	  = ['FakesFF']
      ttH2015.backgrounds = ['FakesFF'] 
      ttH2015.signals     = []
      ttH2015.observed    = ['TTBar']
    elif doABCD:      
      plotbackgrounds	  = ['FakesClosureABCD']
      ttH2015.backgrounds = ['FakesClosureABCD'] # truth cuts done internally in FakesClosureABCD class
      ttH2015.signals     = []
      ttH2015.observed    = ['TTBarClosureMM'] # truth cuts done internally in TTBarClosureMM class
      #ttH2015.observed    = ['FakesClosureMM'] # truth cuts done internally in TTBarClosureMM class
      
    else:
      plotbackgrounds	  = ['TTBar'] 
      ttH2015.backgrounds = ['TTBar'] 
      ttH2015.signals	  = []
      ttH2015.observed    = []      
      #ttH2015.observed    = ['TTBar']           


if args.noSignal:
    ttH2015.signals = []

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
    
    print "Making plots in category: {0}".format( category.name )
    if ( category.cut != None ): 
        print " defined by cuts --> {0}".format( category.cut.cutname )
	    
    signalfactor = 1.0
    background = ttH2015

    # NB: *must* initialise this to 1.0 !!
    #
    lepSF_weight  = '1.0'

    if not ( doMMClosureTest or doMMClosureRates ):
    	if ("ElEl_Event") in category.cut.cutname :
    	   lepSF_weight = 'weight_electron_RecoEff_SF[0] * weight_electron_PIDEff_SF_LHTight[0] * weight_electron_trig[0]'
    	   #
    	   # Do NOT apply LH SF weight ---> is it buggy?
    	   #
    	   #lepSF_weight = 'weight_electron_RecoEff_SF[0] * weight_electron_trig[0]'
    	   print "\t Category contains \'ElEl_Event\': apply this extra weight to events - MC only --> ", lepSF_weight
    
    	elif ("MuMu_Event") in category.cut.cutname :
    	   lepSF_weight = 'weight_muon_RecoEff_SF[0] * weight_muon_IsoEff_SF_Tight[0] * weight_muon_trig[0]'
    	   print "\t Category contains \'MuMu_Event\': apply this extra weight to events - MC only --> ", lepSF_weight
    
    	elif ("MuEl_Event") in category.cut.cutname :
    	   lepSF_weight = 'weight_electron_RecoEff_SF[0] * weight_electron_PIDEff_SF_LHTight[0] * weight_muon_RecoEff_SF[0] * weight_muon_IsoEff_SF_Tight[0] * weight_muon_trig[0]'
    	   #
    	   # Do NOT apply LH SF weight ---> is it buggy?
    	   #
    	   #lepSF_weight = 'weight_electron_RecoEff_SF[0] * weight_muon_RecoEff_SF[0] * weight_muon_IsoEff_SF_Tight[0] * weight_muon_trig[0]'
    	   
    	   print "\t Category contains \'MuEl_Event\': apply this extra weight to events - MC only --> ", lepSF_weight
    
    	elif ("ElMu_Event") in category.cut.cutname :
    	   lepSF_weight = 'weight_electron_RecoEff_SF[0] * weight_electron_PIDEff_SF_LHTight[0] * weight_muon_RecoEff_SF[0] * weight_muon_IsoEff_SF_Tight[0] * weight_electron_trig[0]'
    	   #
    	   # Do NOT apply LH SF weight ---> is it buggy?
    	   #
    	   #lepSF_weight = 'weight_electron_RecoEff_SF[0] * weight_muon_RecoEff_SF[0] * weight_muon_IsoEff_SF_Tight[0] * weight_electron_trig[0]'
    	   
    	   print "\t Category contains \'ElMu_Event\': apply this extra weight to events - MC only --> ", lepSF_weight
    

    # ------------------------------
    # Processing different variables
    # ------------------------------
    for idx,var in enumerate(vardb.varlist, start=0):
        
	# NB: *must* initialise this to 1.0 !!
        #
        bjetSF_weight = '1.0'
        combined_SF_weight = '1.0'

        print "\t now plotting variable: ", var.shortname, "\n"

        # When looking at jet multiplicity distributions w/ bjets, BTag SF must be applied also to categories w/o any bjet cut
        #
	if not ( doMMClosureTest or doMMClosureRates ):
           if ("NBJet") in category.cut.cutname or ("NBJet") in var.shortname:
               bjetSF_weight = 'weight_jet_MV2c20_SFFix70[0]'
               print "\t Category contains \'NBJet\', or plotting variable \'NBjet\' : apply this extra weight to events - MC only --> ", bjetSF_weight
     
        combined_SF_weight = str(lepSF_weight) + ' * ' + str(bjetSF_weight)
            
        print "**************************\n Combined SF weight for events in category:\n ", category.name ,"\n for variable:\n ", var.shortname ,"\n --> ", combined_SF_weight, "\n**************************\n"

        #if idx is 0:
        #    Get event yields for *this* category. Do it only for the 
        #    first variable in the list
        #    
        #    events[category.name] = background.events(cut=cut, eventweight=combined_SF_weight, category=category, hmass=['125'], systematics=systematics, systematicsdirection=systematicsdirection)

	# avoid making useless plots
	#
	if ( ("MuMu") in category.name and ("El") in var.shortname ) or ( ("ElEl") in category.name and ("Mu") in var.shortname ):
            print "\t skipping variable: ", var.shortname
	    continue
	if ( ( ("MuEl") in category.name or ("ElMu") in category.name or ("OF") in category.name ) and ( ("El1") in var.shortname or ("Mu1") in var.shortname ) ) :
            print "\t skipping variable: ", var.shortname
	    continue

	# ---------------------------------------------------------
        # Creating a directory for the category if it doesn't exist
	# ---------------------------------------------------------
        fakeestimate=''
        if doMM:
                fakeestimate='_MM'
        if doFF:
                fakeestimate='_FF'

        dirname =  'OutputPlots' + args.selection + fakeestimate + '_' + args.outdirname + '/'
        
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
	   print "plotname: ", plotname
	
	wantooverflow = True
	
	list_formats = [ plotname + '.png' ] #, plotname + '_canvas.root' ]
	if args.doEPS:
	    list_formats.append( plotname + '.eps' )
	
	doShowRatio = not isblinded
	  
	#
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
            
	    #
	    # systematics go into a different folder
            #
	    dirname = dirname.replace(' ', '_') + '_Syst'
            
	    #
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
		#								       
                # Obtains the total MC histograms with a particular systematics shifted and saving it in the ROOT file
                #
		
		print 'plotname: ', plotname
		systobs, systnom, systup, systdown, systlistup, systlistdown = systs[category.name + ' ' + var.shortname]
		
		print 'systematic: ', syst.name
		
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

        for samp in histograms.keys():
                histograms[samp].SetNameTitle(histname[samp],'')
                histograms[samp].SetLineColor(histcolour[samp])

	if ( 'Mll01' in var.shortname ) or ( 'NJets' in var.shortname ):
		print 'Category: ', category.name
		print 'Variable: ', var.shortname
		print 'Integral: '
		outfile.write('Category: %s \n' %(category.name))
		outfile.write('Variable: %s \n' %(var.shortname))
		outfile.write('Integral: \n')
		err=Double(0)  # integral error
		value=0        # integral value
                for samp in histograms.keys():
			print '%s: %f' %(histname[samp], histograms[samp].Integral())
			value=histograms[samp].IntegralAndError(1,histograms[samp].GetNbinsX(),err)
			outfile.write('yields %s: %f +- %f \n' %(histname[samp], value, err))
		print 'GetEntries: '
		outfile.write('GetEntries: \n')
                for samp in histograms.keys():
			print '%s: %f' %(histname[samp], histograms[samp].GetEntries())
			outfile.write('entries %s: %f \n' %(histname[samp], histograms[samp].GetEntries()))

        for samp in histograms.keys():
                histograms[samp].Write()
        foutput.Close()
	
	if ( 'Mll01' in var.shortname ) or ( 'NJets' in var.shortname ):
		outfile.close()

