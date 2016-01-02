#Run with the command MakePlotsttH_CR.py name_selection name_channel(ex.'TwoLepSR') FakeMethod dolog wantsys lep_flav(ee,em,mm) noSignal
import os, sys, math

sys.path.append(os.path.abspath(os.path.curdir))

from ROOT import gROOT
gROOT.SetBatch(True)#important to run without popups
from ROOT import TH1I, TMath, TFile, TAttFill, TColor, kBlack, kWhite, kBlue, kRed, kYellow, kAzure, kTeal, kSpring, kOrange, kGreen, Double

#importing all the tools and the definitions used to produce the plots 
from Plotter.BackgroundTools_ttH import loadSamples, Category, Background, Process, VariableDB, Variable, Cut, Systematics, Category
#importing the classes for the different processes. They contains many info on the normalization and treatment of the processes
from Plotter.ttH2013_Background_Run2Test import MyCategory, TTHBackgrounds2013

#OPTIONS
debug=True

selection = 'default' # other selection for example are: _noTypeOriginMatch, noTruthMatch
selection = sys.argv[1] #this is the type of reprocessed dataset to use as imput. The name specifies the selection applied
if selection == 'default':
	selection = ''


doSR=False
doStandardPlots=False
doTwoLepSR=False
doThreeLepSR=False
doFourLepSR=False
doMuCR=False
doElCR=False
doTwoLepFakeClosureCR=False
doThreeLepFakeClosureCR=False
doFakeClosureCR=False
doWZonCR=False
doWZoffCR=False
doWZHFonCR=False
doWZHFoffCR=False
dottWCR=False
dottZCR=False
channel = 'TwoLepSR' #this is the channel selected
if len(sys.argv)>2:
	channel = sys.argv[2]
	if(sys.argv[2] not in ['TwoLepSR','ThreeLepSR','FourLepSR','MuCR', 'ElCR', 'TwoLepFakeClosureCR', 'ThreeLepFakeClosureCR', 'WZonCR', 'WZoffCR', 'WZHFonCR', 'WZHFoffCR', 'ttWCR', 'ttZCR']):
		sys.exit('ERROR: the channel specified (', sys.argv[2] ,') is incorrect')
else:
	print 'WARNING: the channel is not specified, using the default 2l channel'
if channel == 'TwoLepSR':
        doTwoLepSR=True
elif channel == 'ThreeLepSR':
        doThreeLepSR=True
elif channel == 'FourLepSR':
        doFourLepSR=True
elif channel == 'MuCR':
        doMuCR=True
elif channel == 'ElCR':
        doElCR=True
elif channel == 'TwoLepFakeClosureCR':
        doTwoLepFakeClosureCR=True
elif channel == 'ThreeLepFakeClosureCR':
        doThreeLepFakeClosureCR=True
elif channel == 'WZonCR':
        doWZonCR=True
elif channel == 'WZoffCR':
        doWZoffCR=True
elif channel == 'WZHFonCR':
        doWZHFonCR=True
elif channel == 'WZHFoffCR':
        doWZHFoffCR=True
elif channel == 'ttWCR':
        dottWCR=True
elif channel == 'ttZCR':
        dottZCR=True

doSR = (doTwoLepSR or doThreeLepSR or doFourLepSR) 
doFakeClosureCR = (doTwoLepFakeClosureCR or doThreeLepFakeClosureCR)
doCR= (doWZonCR or doWZoffCR or doWZHFonCR or doWZHFoffCR or dottWCR or dottZCR)
doStandardPlots = (doSR or doFakeClosureCR or doCR)

doMM=False
doFF=False
method = 'MC' #this is the channel selected
if len(sys.argv)>3:
	method = sys.argv[3]
	if(sys.argv[3] not in ['MC','MM','FF']):
		sys.exit('ERROR: the Fake Method specified (', sys.argv[3] ,') is incorrect')
else:
	print 'WARNING: the channel is not specified, using the default MC method'
if method == 'MM':
        doMM=True
elif method == 'FF':
        doFF=True


dolog = False #
if len(sys.argv)>4:
	if sys.argv[4] == 'True':
		dolog = True

wantsysts = False
if len(sys.argv)>5:
	if sys.argv[5] == 'True':
		wantsysts = True

lep_flav = '' #this indicates the flavour of the SS pair of lepton selected
if len(sys.argv)>6:
	if(sys.argv[6] not in ['ee','mm','em']):
                print 'WARNING: the lep_flav is not specified, not using any restriction on the lepton flavour'
        else:
                lep_flav = sys.argv[6]

noSignal = False
if len(sys.argv)>7:
	if sys.argv[7] == 'True':
		noSignal = True

if debug:
        print 'selection ', selection, ' channel ', channel, ' method ', method, ' dolog ', dolog, ' wantsysts ', wantsysts

inputs = loadSamples(   
                        #path of the data to be processed
                        inputdir = '/home/data/ttH/NewReprocessed_Melb13_ttH_Run2Test'+selection,
                        samplescsv  = 'Files/samples2013_ttH_Run2Test.csv',#Files containing the processes of interest with their cross sections and other info
                        nomtree  = 'physics',# name of the TTree inside the files
                        #name of the trees that contains values for shifted systematics
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
# Here you include all names of variables to be plotted, with min, max, number of bins and ntuple name.
vardb = VariableDB()

# A list of cuts
#vardb.registerCut( Cut('HFOR',      'top_hfor_type!=4') )  #DUMMY VALUE BECAUSE WE DON'T HAVE THIS INFO
vardb.registerCut( Cut('HFOR',      '1!=4') )
vardb.registerCut( Cut('TrigMatch', 'evtsel_is_sltmatch') ) #FILLED WITH 1 IF AT LEAST ONE LEPTONS HAS PT>25
vardb.registerCut( Cut('NBJet',     'evtsel_bjets_num>=1') )
vardb.registerCut( Cut('TwoBJets',  'evtsel_bjets_num>=2') )
vardb.registerCut( Cut('BJetVeto',  'evtsel_bjets_num==0') )
vardb.registerCut( Cut('DiLep',     'evtsel_dilep_type') )
vardb.registerCut( Cut('DiLepTau',  'evtsel_dileptau_type') )
vardb.registerCut( Cut('TriLep',    'evtsel_trilep_type') )
vardb.registerCut( Cut('QuadLep',   'evtsel_fourlep_type') )
vardb.registerCut( Cut('NJet2L',    'evtsel_jets_num>=4') )
vardb.registerCut( Cut('NJet3L',    '(evtsel_jets_num>=4 || (evtsel_bjets_num>=2 && evtsel_jets_num==3))') )
vardb.registerCut( Cut('NJet4L',    'evtsel_jets_num>=2') )
vardb.registerCut( Cut('LowJetCR',  'evtsel_jets_num>=1 && evtsel_jets_num<=3') )
vardb.registerCut( Cut('JPsiRemovalCR',  '((evtsel_JPsiee_candidate_mass>12 || evtsel_JPsiee_candidate_mass<0 ) && (evtsel_JPsimm_candidate_mass>12 || evtsel_JPsimm_candidate_mass<0))') )
vardb.registerCut( Cut('NoTau',     'evtsel_tau_num==0') )
vardb.registerCut( Cut('OneTau',    'evtsel_tau_num==1') )
vardb.registerCut( Cut('ElEta',    '((fabs(Lep0PDG)!=11 || fabs(Lep0Eta)<1.5) && (fabs(Lep1PDG)!=11 || fabs(Lep1Eta)<1.5))') )
vardb.registerCut( Cut('ZSSmasscut',     'fabs(evtsel_ZSSee_candidate_mass-91)>20') )
vardb.registerCut( Cut('ZOSmasscut',     '(fabs(evtsel_Zmm_candidate_mass-91)>10 && fabs(evtsel_Zee_candidate_mass-91)>10)') )
vardb.registerCut( Cut('ZOSpeakcut',     '(fabs(evtsel_Zmm_candidate_mass-91)<10 || fabs(evtsel_Zee_candidate_mass-91)<10)') )
vardb.registerCut( Cut('2LepSS',    '(isSS01==1 && Lep0Pt>20 && Lep1Pt>20)') )
vardb.registerCut( Cut('2LepSSTau', '(Lep0Charge*Tau0Charge==-1 && isSS01==1 && Lep0Pt>15 && Lep1Pt>15)') )
vardb.registerCut( Cut('3Lep',      '(isSS12==1 && Lep0Tight==1 && fabs(Lep0Charge+Lep1Charge+Lep2Charge)==1 && Lep1Pt>20 && Lep2Pt>20 && Mll01>12 && Mll02>12)') )
vardb.registerCut( Cut('4Lep',      '(Lep0Pt>25 && Lep1Pt>15 && Lep0Tight==1 && Lep1Tight==1 && Lep2Tight==1 && Lep3Tight==1 && fabs(Lep0Charge+Lep1Charge+Lep2Charge+Lep3Charge)==0 && ((evtsel_JPsiee_candidate_mass>10 || evtsel_JPsiee_candidate_mass<0 ) && (evtsel_JPsimm_candidate_mass>10 || evtsel_JPsimm_candidate_mass<0)))') )
#used for the real and fake estimate in the 2lep SS channel
vardb.registerCut( Cut('2LepTT',    'isTT01==1') )
vardb.registerCut( Cut('2LepTL',    'isTL01==1') )#L means that the lepton is loose&&!tight
vardb.registerCut( Cut('2LepLL',    'isLL01==1') )#L means that the lepton is loose&&!tight
#used for the real and fake estimate in the 3lep channel
vardb.registerCut( Cut('3LepTT',    'isTT12==1') )
vardb.registerCut( Cut('3LepTL',    'isTL12==1') )#L means that the lepton is loose&&!tight
vardb.registerCut( Cut('3LepLL',    'isLL12==1') )#L means that the lepton is loose&&!tight

vardb.registerCut( Cut('ElRealOSCR',    '(isSS01==0 && isElRateEvt)') )
vardb.registerCut( Cut('ElFakeSSCR',    '(isSS01==1 && isElRateEvt)') )
vardb.registerCut( Cut('MuRealOSCR',    '(isSS01==0 && isMuRateEvt)') )
vardb.registerCut( Cut('MuFakeSSCR',    '(isSS01==1 && isMuRateEvt)') )
vardb.registerCut( Cut('ElTagEta',      '(ElTagEta==-5. || fabs(ElTagEta)<1.5)') )#to remove some charge flip
vardb.registerCut( Cut('ElProbeTight',    'ElProbeTight==1') )
vardb.registerCut( Cut('MuProbeTight',    'MuProbeTight==1') )
vardb.registerCut( Cut('ElProbeLoose',    'ElProbeTight==0') )
vardb.registerCut( Cut('MuProbeLoose',    'MuProbeTight==0') )
vardb.registerCut( Cut('MuMu_event',  'evtsel_is_dimu') )
vardb.registerCut( Cut('MuEl_event',  'evtsel_is_muel') )
vardb.registerCut( Cut('ElEl_event',  'evtsel_is_diel') )
## vardb.registerCut( Cut('PeriodAB',    '((evtsel_run_number<200000 && (evtsel_randomRunNumber>=200804 && evtsel_randomRunNumber<=205113)) || (evtsel_run_number>=200804 && evtsel_run_number<=205113))') )
## vardb.registerCut( Cut('PeriodCDE',    '((evtsel_run_number<200000 && (evtsel_randomRunNumber>=206248 && evtsel_randomRunNumber<=210308)) || (evtsel_run_number>=206248 && evtsel_run_number<=210308))') )

# A list of variables to plot
#vardb.registerVar( Variable(shortname = 'MET', latexname = 'E_{T}^{miss} [GeV]', ntuplename = 'evtsel_met_et', bins = 20, minval = 0., maxval = 200.,) )
#vardb.registerVar( Variable(shortname = 'Jet0Pt', latexname = 'p_{T}^{lead jet} [GeV]', ntuplename = 'Jet0Pt', bins = 18, minval = 20., maxval = 200.,) )
vardb.registerVar( Variable(shortname = 'NJets', latexname = 'Jet multiplicity', ntuplename = 'evtsel_jets_num', bins = 10, minval = 0, maxval = 10) )
vardb.registerVar( Variable(shortname = 'NBJets', latexname = 'BJet multiplicity', ntuplename = 'evtsel_bjets_num', bins = 10, minval = 0, maxval = 10) )
#vardb.registerVar( Variable(shortname = 'Tau0Pt', latexname = 'p_{T}^{lead tau} [GeV]', ntuplename = 'Tau0Pt', bins = 30, minval = 25., maxval = 100.,) )
vardb.registerVar( Variable(shortname = 'Mll01', latexname = 'm(l0^{#pm}l1^{#pm}) [GeV]', ntuplename = 'Mll01',  bins = 15, minval = 0., maxval = 300.,) )
#vardb.registerVar( Variable(shortname = 'Mll12', latexname = 'm(l1^{#pm}l2^{#pm}) [GeV]', ntuplename = 'Mll12',  bins = 15, minval = 0., maxval = 300.,) )
if doStandardPlots:
        vardb.registerVar( Variable(shortname = 'Lep0Pt', latexname = 'Lep0 p_{T} [GeV]', ntuplename = 'Lep0Pt', bins = 18, minval = 10., maxval = 100.,) )
        vardb.registerVar( Variable(shortname = 'Lep0Eta', latexname = 'Lep0 #eta', ntuplename = 'fabs(Lep0Eta)', bins = 8, minval = 0., maxval = 2.6, manualbins = [ 0.0 , 0.5 , 0.8 , 1.1 , 1.37 , 1.52 , 2.0 , 2.25 , 2.6],) )
        vardb.registerVar( Variable(shortname = 'Lep1Pt', latexname = 'Lep1 p_{T} [GeV]', ntuplename = 'Lep1Pt', bins = 18, minval = 10., maxval = 100.,) )
        vardb.registerVar( Variable(shortname = 'Lep1Eta', latexname = 'Lep1 #eta', ntuplename = 'fabs(Lep1Eta)', bins = 8, minval = 0., maxval = 2.6, manualbins = [ 0.0 , 0.5 , 0.8 , 1.1 , 1.37 , 1.52 , 2.0 , 2.25 , 2.6],) )
if doElCR:
        #vardb.registerVar( Variable(shortname = 'mTElTag', latexname = 'm_{T}^{lep,MET} [GeV]', ntuplename = 'mTElTag', bins = 32, minval = 0., maxval = 160.,) )
        #vardb.registerVar( Variable(shortname = 'mTMuTag', latexname = 'm_{T}^{lep,MET} [GeV]', ntuplename = 'mTMuTag', bins = 32, minval = 0., maxval = 160.,) )
        #vardb.registerVar( Variable(shortname = 'ElTagPt', latexname = 'El Tag p_{T} [GeV]', ntuplename = 'ElTagPt', bins = 36, minval = 10., maxval = 100.,) )
        #vardb.registerVar( Variable(shortname = 'ElTagEta', latexname = 'El Tag #eta', ntuplename = 'fabs(ElTagEta)', bins = 8, minval = 0., maxval = 2.6, manualbins = [ 0.0 , 0.5 , 0.8 , 1.1 , 1.37 , 1.52 , 2.0 , 2.25 , 2.6],) )
        #vardb.registerVar( Variable(shortname = 'MuTagPt', latexname = 'Mu Tag p_{T} [GeV]', ntuplename = 'MuTagPt', bins = 36, minval = 10., maxval = 100.,) )
        #vardb.registerVar( Variable(shortname = 'MuTagEta', latexname = 'Mu Tag #eta', ntuplename = 'fabs(MuTagEta)', bins = 8, minval = 0., maxval = 2.6, manualbins = [ 0.0 , 0.5 , 0.8 , 1.1 , 1.37 , 1.52 , 2.0 , 2.25 , 2.6],) )
        vardb.registerVar( Variable(shortname = 'mTProbe', latexname = 'm_{T}^{lep,MET} [GeV]', ntuplename = 'mTElProbe', bins = 32, minval = 0., maxval = 160.,) )
        #vardb.registerVar( Variable(shortname = 'ProbePt', latexname = 'Probe p_{T} [GeV]', ntuplename = 'ElProbePt', bins = 16, minval = 20., maxval = 2000., manualbins = [10.0, 15.0, 20.0 , 22.0 , 25.0 , 30.0 , 35.0 , 40.0 , 50.0 , 60.0 , 80.0 , 100.0 , 120.0 , 150.0 , 200.0 , 300.0 , 2000.0],) )
        vardb.registerVar( Variable(shortname = 'ProbePt', latexname = 'Probe p_{T} [GeV]', ntuplename = 'ElProbePt', bins = 36, minval = 10., maxval = 100.,) )
        vardb.registerVar( Variable(shortname = 'ProbeEta', latexname = 'Probe #eta', ntuplename = 'fabs(ElProbeEta)', bins = 8, minval = 0., maxval = 2.6, manualbins = [ 0.0 , 0.5 , 0.8 , 1.1 , 1.37 , 1.52 , 2.0 , 2.25 , 2.6],) )
if doMuCR:
        #vardb.registerVar( Variable(shortname = 'mTElTag', latexname = 'm_{T}^{lep,MET} [GeV]', ntuplename = 'mTElTag', bins = 32, minval = 0., maxval = 160.,) )
        #vardb.registerVar( Variable(shortname = 'mTMuTag', latexname = 'm_{T}^{lep,MET} [GeV]', ntuplename = 'mTMuTag', bins = 32, minval = 0., maxval = 160.,) )
        #vardb.registerVar( Variable(shortname = 'ElTagPt', latexname = 'El Tag p_{T} [GeV]', ntuplename = 'ElTagPt', bins = 36, minval = 10., maxval = 100.,) )
        #vardb.registerVar( Variable(shortname = 'ElTagEta', latexname = 'El Tag #eta', ntuplename = 'fabs(ElTagEta)', bins = 8, minval = 0., maxval = 2.6, manualbins = [ 0.0 , 0.5 , 0.8 , 1.1 , 1.37 , 1.52 , 2.0 , 2.25 , 2.6],) )
        #vardb.registerVar( Variable(shortname = 'MuTagPt', latexname = 'Mu Tag p_{T} [GeV]', ntuplename = 'MuTagPt', bins = 36, minval = 10., maxval = 100.,) )
        #vardb.registerVar( Variable(shortname = 'MuTagEta', latexname = 'Mu Tag #eta', ntuplename = 'fabs(MuTagEta)', bins = 8, minval = 0., maxval = 2.6, manualbins = [ 0.0 , 0.5 , 0.8 , 1.1 , 1.37 , 1.52 , 2.0 , 2.25 , 2.6],) )
        vardb.registerVar( Variable(shortname = 'mTProbe', latexname = 'm_{T}^{lep,MET} [GeV]', ntuplename = 'mTMuProbe', bins = 40, minval = 0., maxval = 200.,) )
        #vardb.registerVar( Variable(shortname = 'ProbePt', latexname = 'Probe p_{T} [GeV]', ntuplename = 'MuProbePt', bins = 16, minval = 20., maxval = 2000., manualbins = [10.0, 15.0, 20.0 , 22.0 , 25.0 , 30.0 , 35.0 , 40.0 , 50.0 , 60.0 , 80.0 , 100.0 , 120.0 , 150.0 , 200.0 , 300.0 , 2000.0],) )
        vardb.registerVar( Variable(shortname = 'ProbePt', latexname = 'Probe p_{T} [GeV]', ntuplename = 'MuProbePt', bins = 36, minval = 10., maxval = 100.,) )
        vardb.registerVar( Variable(shortname = 'ProbeEta', latexname = 'Probe #eta', ntuplename = 'fabs(MuProbeEta)', bins = 8, minval = 0., maxval = 2.6, manualbins = [ 0.0 , 0.5 , 0.8 , 1.1 , 1.37 , 1.52 , 2.0 , 2.25 , 2.6],) )


#vardb.registerVar( Variable(shortname = 'avgint', latexname = 'Average Interactions', ntuplename = 'averageIntPerXing', bins = 40, minval = 0, maxval = 40, typeval = TH1I) )

		
#alterantive ranges and binning for the histograms 
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


# A list of systematics
if wantsysts:
        ## #uncertainties on the kfactors used to normalize the various MC distributions 
        #if doTwoLepSR or doTwoLepFakeClosureCR or dottWCR:
        #        vardb.registerSystematics( Systematics(name='CFsys',             eventweight='sys_weight_CF_') )
        if doMM:
                vardb.registerSystematics( Systematics(name='MMrsys',           eventweight='sys_weight_MMr_') )
                vardb.registerSystematics( Systematics(name='MMfsys',           eventweight='sys_weight_MMf_') )
        if doFF:
                vardb.registerSystematics( Systematics(name='FFsys',           eventweight='sys_weight_FF_') )
	'''
	vardb.registerSystematics( Systematics(name='PU',             eventweight='evtsel_sys_PU_rescaling_') )
	vardb.registerSystematics( Systematics(name='el_reco',        eventweight='evtsel_sys_sf_el_reco_') )
	vardb.registerSystematics( Systematics(name='el_id',          eventweight='evtsel_sys_sf_el_id_') )
	#vardb.registerSystematics( Systematics(name='el_iso',         eventweight='evtsel_sys_sf_el_iso_') )
	vardb.registerSystematics( Systematics(name='mu_id',          eventweight='evtsel_sys_sf_mu_id_') )
	#vardb.registerSystematics( Systematics(name='mu_iso',         eventweight='evtsel_sys_sf_mu_iso_') )
        vardb.registerSystematics( Systematics(name='lep_trig',       eventweight='evtsel_sys_sf_lep_trig_') )
        vardb.registerSystematics( Systematics(name='bjet_b',         eventweight='evtsel_sys_sf_bjet_b_') )
        vardb.registerSystematics( Systematics(name='bjet_c',         eventweight='evtsel_sys_sf_bjet_c_') )
        vardb.registerSystematics( Systematics(name='bjet_m',         eventweight='evtsel_sys_sf_bjet_m_') )
	'''
        ## #vardb.registerSystematics( Systematics(name='METSys',        treename='METSys') )
        #vardb.registerSystematics( Systematics(name='ElEnResSys',     treename='ElEnResSys') )
        #vardb.registerSystematics( Systematics(name='ElES_LowPt',     treename='ElES_LowPt') )
        #vardb.registerSystematics( Systematics(name='ElES_Zee',       treename='ElES_Zee') )
        #vardb.registerSystematics( Systematics(name='ElES_R12',       treename='ElES_R12') )
        #vardb.registerSystematics( Systematics(name='ElES_PS',        treename='ElES_PS') )
        ## #vardb.registerSystematics( Systematics(name='EESSys',        treename='EESSys') )
        #vardb.registerSystematics( Systematics(name='MuSys',          treename='MuSys') )
        #vardb.registerSystematics( Systematics(name='JES_Total',      treename='JES_Total') )
        #vardb.registerSystematics( Systematics(name='JER',            treename='JER') )


#definition of the categories for which one want produce histograms
if doTwoLepSR:
        vardb.registerCategory( MyCategory('MuMuSS',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLep', 'NJet2L', 'NoTau', '2LepSS', 'MuMu_event',])) )
        vardb.registerCategory( MyCategory('MuElSS',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLep', 'NJet2L', 'NoTau', '2LepSS', 'MuEl_event', 'ElEta'])) )
        vardb.registerCategory( MyCategory('ElElSS',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLep', 'NJet2L', 'NoTau', '2LepSS', 'ElEl_event', 'ElEta'])) )
        vardb.registerCategory( MyCategory('TwoLepSSTau',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLepTau', 'NJet2L', 'OneTau', '2LepSSTau', 'ZSSmasscut',])) )
if doThreeLepSR:
        vardb.registerCategory( MyCategory('ThreeLep',     cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'TriLep', 'NJet3L', '3Lep', 'ZOSmasscut',])) )
if doFourLepSR:
        vardb.registerCategory( MyCategory('FourLep',      cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'QuadLep', 'NJet4L', '4Lep',])) )
if doTwoLepFakeClosureCR:
        vardb.registerCategory( MyCategory('MuMuSSFakeClosureCR',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLep', 'LowJetCR', 'NoTau', '2LepSS', 'MuMu_event',])) )
        vardb.registerCategory( MyCategory('MuElSSFakeClosureCR',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLep', 'LowJetCR', 'NoTau', '2LepSS', 'MuEl_event', 'ElEta'])) )
        vardb.registerCategory( MyCategory('ElElSSFakeClosureCR',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLep', 'LowJetCR', 'NoTau', '2LepSS', 'ElEl_event', 'ElEta'])) )
        vardb.registerCategory( MyCategory('TwoLepSSTauFakeClosureCR',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLepTau', 'LowJetCR', 'OneTau', '2LepSSTau', 'ZSSmasscut',])) )
if doThreeLepFakeClosureCR:
        vardb.registerCategory( MyCategory('ThreeLepFakeClosureCR',     cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'TriLep', 'LowJetCR', '3Lep', 'ZOSmasscut',])) )

if doWZonCR:
        vardb.registerCategory( MyCategory('WZonCR',     cut=vardb.getCuts(['HFOR', 'TrigMatch', 'BJetVeto', 'TriLep', '3Lep', 'ZOSpeakcut',])) )
if doWZoffCR:
        vardb.registerCategory( MyCategory('WZoffCR',     cut=vardb.getCuts(['HFOR', 'TrigMatch', 'BJetVeto', 'TriLep', '3Lep', 'ZOSmasscut',])) )
if doWZHFonCR:
        vardb.registerCategory( MyCategory('WZHFonCR',     cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'TriLep', '3Lep', 'ZOSpeakcut',])) )
if doWZHFoffCR:
        vardb.registerCategory( MyCategory('WZHFoffCR',     cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'TriLep', '3Lep', 'ZOSmasscut',])) )
if dottZCR:
        vardb.registerCategory( MyCategory('ttZCR',     cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'TriLep', 'NJet3L', '3Lep', 'ZOSpeakcut',])) )
if dottWCR:
        vardb.registerCategory( MyCategory('ttWCR',     cut=vardb.getCuts(['HFOR', 'TrigMatch', 'TwoBJets', 'DiLep', 'LowJetCR', 'NoTau', '2LepSS', 'ElEta',])) )

#loose means that is not tight and loose
if doMuCR:
        vardb.registerCategory( MyCategory('FakeCRMuL',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLep', 'LowJetCR', 'MuFakeSSCR', 'MuProbeLoose', 'ElTagEta',])) )
        vardb.registerCategory( MyCategory('MuMuFakeCRMuL',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLep', 'LowJetCR', 'MuFakeSSCR', 'MuProbeLoose', 'MuMu_event',])) )
        vardb.registerCategory( MyCategory('MuElFakeCRMuL',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLep', 'LowJetCR', 'MuFakeSSCR', 'MuProbeLoose', 'MuEl_event', 'ElTagEta',])) )
        vardb.registerCategory( MyCategory('FakeCRMuT',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLep', 'LowJetCR', 'MuFakeSSCR', 'MuProbeTight', 'ElTagEta',])) )
        vardb.registerCategory( MyCategory('MuMuFakeCRMuT',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLep', 'LowJetCR', 'MuFakeSSCR', 'MuProbeTight', 'MuMu_event',])) )
        vardb.registerCategory( MyCategory('MuElFakeCRMuT',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLep', 'LowJetCR', 'MuFakeSSCR', 'MuProbeTight', 'MuEl_event', 'ElTagEta',])) )
        vardb.registerCategory( MyCategory('RealCRMuL',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLep', 'LowJetCR', 'MuRealOSCR', 'JPsiRemovalCR', 'MuProbeLoose',])) )
        vardb.registerCategory( MyCategory('MuMuRealCRMuL',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLep', 'LowJetCR', 'MuRealOSCR', 'JPsiRemovalCR', 'MuProbeLoose', 'MuMu_event',])) )
        vardb.registerCategory( MyCategory('MuElRealCRMuL',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLep', 'LowJetCR', 'MuRealOSCR', 'JPsiRemovalCR', 'MuProbeLoose', 'MuEl_event',])) )
        vardb.registerCategory( MyCategory('RealCRMuT',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLep', 'LowJetCR', 'MuRealOSCR', 'JPsiRemovalCR', 'MuProbeTight',])) )
        vardb.registerCategory( MyCategory('MuMuRealCRMuT',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLep', 'LowJetCR', 'MuRealOSCR', 'JPsiRemovalCR', 'MuProbeTight', 'MuMu_event',])) )
        vardb.registerCategory( MyCategory('MuElRealCRMuT',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLep', 'LowJetCR', 'MuRealOSCR', 'JPsiRemovalCR', 'MuProbeTight', 'MuEl_event',])) )

if doElCR:
        vardb.registerCategory( MyCategory('FakeCRElL',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLep', 'LowJetCR', 'ElFakeSSCR', 'ZSSmasscut', 'ElProbeLoose',])) )
        vardb.registerCategory( MyCategory('ElElFakeCRElL',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLep', 'LowJetCR', 'ElFakeSSCR', 'ZSSmasscut', 'ElProbeLoose', 'ElEl_event',])) )
        vardb.registerCategory( MyCategory('MuElFakeCRElL',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLep', 'LowJetCR', 'ElFakeSSCR', 'ElProbeLoose', 'MuEl_event',])) )
        vardb.registerCategory( MyCategory('FakeCRElT',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLep', 'LowJetCR', 'ElFakeSSCR', 'ZSSmasscut', 'ElProbeTight',])) )
        vardb.registerCategory( MyCategory('ElElFakeCRElT',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLep', 'LowJetCR', 'ElFakeSSCR', 'ZSSmasscut', 'ElProbeTight', 'ElEl_event',])) )
        vardb.registerCategory( MyCategory('MuElFakeCRElT',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLep', 'LowJetCR', 'ElFakeSSCR', 'ElProbeTight', 'MuEl_event',])) )
        vardb.registerCategory( MyCategory('RealCRElL',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLep', 'LowJetCR', 'ElRealOSCR', 'JPsiRemovalCR', 'ElProbeLoose',])) )
        vardb.registerCategory( MyCategory('ElElRealCRElL',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLep', 'LowJetCR', 'ElRealOSCR', 'JPsiRemovalCR', 'ElProbeLoose', 'ElEl_event',])) )
        vardb.registerCategory( MyCategory('MuElRealCRElL',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLep', 'LowJetCR', 'ElRealOSCR', 'JPsiRemovalCR', 'ElProbeLoose', 'MuEl_event',])) )
        vardb.registerCategory( MyCategory('RealCRElT',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLep', 'LowJetCR', 'ElRealOSCR', 'JPsiRemovalCR', 'ElProbeTight',])) )
        vardb.registerCategory( MyCategory('ElElRealCRElT',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLep', 'LowJetCR', 'ElRealOSCR', 'JPsiRemovalCR', 'ElProbeTight', 'ElEl_event',])) )
        vardb.registerCategory( MyCategory('MuElRealCRElT',  cut=vardb.getCuts(['HFOR', 'TrigMatch', 'NBJet', 'DiLep', 'LowJetCR', 'ElRealOSCR', 'JPsiRemovalCR', 'ElProbeTight', 'MuEl_event',])) )

#TTHBackgrounds2013 is what is used to manage each process. One passes to it the input informations and the definitions and it calculates 
tth2013 = TTHBackgrounds2013(inputs, vardb)


## # Settings for 2012
## #tth2013.turnOffRQCD = True
## tth2013.useEmbedding = False

tth2013.luminosity = 5.4
#tth2013.luminosity = 13.0605
#tth2013.luminosity = 5.83521

tth2013.useZCorrections = False #because not using the alpgen pythia samples

isblinded=0
# Make blinded plots
#if doSR:
#        tth2013.observed = []
#        isblinded=1

if doTwoLepSR or doTwoLepFakeClosureCR or dottWCR:
	tth2013.channel = 'TwoLepSS'
elif doThreeLepSR or doThreeLepFakeClosureCR or dottZCR or doWZonCR or doWZoffCR or doWZHFonCR or doWZHFoffCR:
	tth2013.channel = 'ThreeLep'
elif doFourLepSR:
	tth2013.channel = 'FourLep'
elif doMuCR or doElCR:
	tth2013.channel = 'TwoLepCR'

#cut = vardb.getCut('InvDDR')
cut = None
systematics = None
systematicsdirection = None # 'UP', 'DOWN'
#hmass = [mass_value]
events = {}
hists = {}
#dict with systematics histograms
systs = {}
#list of the backgrounds considered
samplenames = { 'Observed':'observed', 'TTBarH':'signal', 'TTBarW':'ttbarwbkg', 'TTBarZ':'ttbarzbkg', 'Top':'topbkg', 'TopCF':'topcfbkg', 'Diboson':'dibosonbkg', 'DibosonCF':'dibosoncfbkg', 'HtoZZ':'htozzbkg', 'Zjets':'zjetsbkg', 'ZjetsLF':'zjetsbkg', 'Wjets':'wjetsbkg', 'Prompt':'promptbkg', 'ChargeFlip':'chargeflipbkg', 'FakesFF':'fakesbkg', 'FakesMM':'fakesbkg' }
colors = { 'Observed':kBlack, 'TTBarH':kBlack, 'TTBarW':kRed-4, 'TTBarZ':kRed-7, 'Top':kBlue-3, 'TopCF':kAzure-4, 'Diboson':kOrange, 'DibosonCF':kOrange-3, 'HtoZZ':kTeal+9, 'Zjets':kGreen, 'ZjetsLF':kGreen, 'Wjets':kYellow, 'Prompt':kOrange, 'ChargeFlip':kAzure-4, 'FakesFF':kAzure-9, 'FakesMM':kTeal-9 }

if doStandardPlots:
        if doMM:
                plotbackgrounds = [ 'TTBarW', 'TTBarZ', 'Top', 'TopCF', 'Diboson', 'DibosonCF', 'ZjetsLF', 'FakesMM']
                tth2013.backgrounds = [ 'TTBarW', 'TTBarZ', 'Top', 'TopCF', 'Diboson', 'DibosonCF', 'ZjetsLF', 'FakesMM']#this is the list of backgrounds that will be calculated
        elif doFF:
                plotbackgrounds = [ 'TTBarW', 'TTBarZ', 'Top', 'TopCF', 'Diboson', 'DibosonCF', 'ZjetsLF', 'FakesFF']
                tth2013.backgrounds = [ 'TTBarW', 'TTBarZ', 'Top', 'TopCF', 'Diboson', 'DibosonCF', 'ZjetsLF', 'FakesFF']#this is the list of backgrounds that will be calculated
        else:
                #MC based estimate of fakes
                plotbackgrounds = [ 'TTBarW', 'TTBarZ', 'Top', 'TopCF', 'Diboson', 'DibosonCF', 'ZjetsLF']
                tth2013.backgrounds = [ 'TTBarW', 'TTBarZ', 'Top', 'TopCF', 'Diboson', 'DibosonCF', 'ZjetsLF']#this is the list of backgrounds that will be calculated
	if doFourLepSR:
		#no fakes
                plotbackgrounds = [ 'TTBarW', 'TTBarZ', 'Top', 'TopCF', 'Diboson', 'DibosonCF', 'ZjetsLF']
                tth2013.backgrounds = [ 'TTBarW', 'TTBarZ', 'Top', 'TopCF', 'Diboson', 'DibosonCF', 'ZjetsLF']#this is the list of backgrounds that will be calculated

if doElCR or doMuCR:
        #here the fakes are the difference between data and mc
        plotbackgrounds = [ 'TTBarW', 'TTBarZ', 'Top', 'TopCF', 'Diboson', 'DibosonCF', 'ZjetsLF']#, 'Wjets']#The W+jets should be 0 in the truthmatching case and should be a rough estimate of fakes in the other case. Also ttbar should be bigger in the no truth matching case
        tth2013.backgrounds = [ 'TTBarW', 'TTBarZ', 'Top', 'TopCF', 'Diboson', 'DibosonCF', 'ZjetsLF']#, 'Wjets']#this is the list of backgrounds that will be calculated

if noSignal:
	tth2013.signals = []

#Filling histname with the name of the variables we want
histname = {'Expected':'expected'}
histcolor = {'Expected':kBlack}
for samp in tth2013.backgrounds:
        histname[samp] = samplenames[samp]
        histcolor[samp] = colors[samp]
for samp in tth2013.observed:
        histname[samp] = samplenames[samp]
        histcolor[samp] = colors[samp]
for samp in tth2013.signals:
        histname[samp] = samplenames[samp]
        histcolor[samp] = colors[samp]
print histname
print histcolor

#processing categories in sequence
for category in vardb.categorylist:
    print "Making plot in category:", category.name
    signalfactor = 1.
    background = tth2013

    events[category.name] = background.events(cut=cut, category=category, hmass=['125'],systematics=systematics, systematicsdirection=systematicsdirection)
    #processing different variables
    for var in vardb.varlist:
        #creating a directory for the category if it doesn't exist
        fakeestimate=''
        if doMM:
                fakeestimate='_MM'
        if doFF:
                fakeestimate='_FF'
                
        if cut:
            dirname = 'output/ttH/Run2Test_results'+selection+fakeestimate+'/'  + category.name + ' ' + cut.cutname
        else:
            dirname = 'output/ttH/Run2Test_results'+selection+fakeestimate+'/'  + category.name
	if dolog:
		dirname = dirname+'_LOG'
        dirname = dirname.replace(' ', '_')
        try:
            os.makedirs(dirname)
        except:
            pass
        #producing a plot with the name of the category and the variable
        plotname = dirname + '/' + category.name + ' ' + var.shortname
        plotname = plotname.replace(' ', '_')
	wantooverflow=True
        hists[category.name + ' ' + var.shortname] = background.plot(var, cut=cut, category=category, signal='125', signalfactor=signalfactor, overridebackground=plotbackgrounds, systematics=systematics, systematicsdirection=systematicsdirection, overflowbins=wantooverflow, showratio= not isblinded, wait=False, save=[plotname+'.png', plotname+'.eps', plotname+'_canvas.root'], log=dolog)

        # creating a file with the observed and expected distributions and systematics. We fit them for TES uncertainty studies
        foutput = TFile(plotname+'.root','RECREATE')
	if 'Mll01' in var.shortname:
		outfile = open(plotname+'_yields.txt', 'w')
        if wantsysts:
            #systematics go into a different folder
            dirname = dirname.replace(' ', '_') + '_Syst'
            #loop on the defined systematics
	    total_syst = 0
            histograms_syst = {}
            for syst in vardb.systlist:
                try:
                    os.makedirs(dirname)
                except:
                    pass
                plotname = dirname + '/' + category.name + ' ' + var.shortname + ' ' + syst.name
                plotname = plotname.replace(' ', '_')
                #plotSystematics is the function which takes care of the systematics
                systs[category.name + ' ' + var.shortname] = background.plotSystematics(syst, var=var, cut=cut, category=category, overflowbins=wantooverflow, showratio=True, wait=False, save=[plotname+'.png', plotname+'.eps'])
                #obtaining the total mc histograms with a particular systematics shifted and saving it in the root file
                systobs, systnom, systup, systdown, systlistup, systlistdown = systs[category.name + ' ' + var.shortname]
                histograms_syst['Expected_'+syst.name+'_up']=systup
                histograms_syst['Expected_'+syst.name+'_up'].SetNameTitle(histname['Expected']+'_'+syst.name+'_up','')
                histograms_syst['Expected_'+syst.name+'_up'].SetLineColor(histcolor['Expected'])
                histograms_syst['Expected_'+syst.name+'_up'].Write()
                histograms_syst['Expected_'+syst.name+'_down']=systdown
                histograms_syst['Expected_'+syst.name+'_down'].SetNameTitle(histname['Expected']+'_'+syst.name+'_down','')
                histograms_syst['Expected_'+syst.name+'_down'].SetLineColor(histcolor['Expected'])
                histograms_syst['Expected_'+syst.name+'_down'].Write()
                for samp in tth2013.backgrounds:
                        histograms_syst[samp+'_'+syst.name+'_up'] = systlistup[samp]
                        histograms_syst[samp+'_'+syst.name+'_up'].SetNameTitle(histname[samp]+'_'+syst.name+'_up','')
                        histograms_syst[samp+'_'+syst.name+'_up'].SetLineColor(histcolor[samp])
                        histograms_syst[samp+'_'+syst.name+'_up'].Write()
                        histograms_syst[samp+'_'+syst.name+'_down'] = systlistdown[samp]
                        histograms_syst[samp+'_'+syst.name+'_down'].SetNameTitle(histname[samp]+'_'+syst.name+'_down','')
                        histograms_syst[samp+'_'+syst.name+'_down'].SetLineColor(histcolor[samp])
                        histograms_syst[samp+'_'+syst.name+'_down'].Write()
                #ACTUALLY THE CODE DOES NOT CONSIDER SYSTEMATICS FOR THE SIGNAL. PUT IT AMONG THE BACKGROUNDS IF YOU WANT SYST ON IT        
		if 'Mll01' in var.shortname:
			outfile.write('Integral syst: \n')
			outfile.write('yields syst %s up: %f \n' %(syst.name,(systup.Integral()-systnom.Integral())))
			outfile.write('yields syst %s down: %f \n' %(syst.name,(systdown.Integral()-systnom.Integral())))
			outfile.write('GetEntries syst: \n')
			outfile.write('entries syst %s up: %f \n' %(syst.name,(systup.GetEntries()-systnom.GetEntries())))
			outfile.write('entries syst %s down: %f \n' %(syst.name,(systdown.GetEntries()-systnom.GetEntries())))
		total_syst = total_syst + (systup.Integral()-systdown.Integral())/2*(systup.Integral()-systdown.Integral())/2
	    total_syst = math.sqrt(total_syst)
	    if 'Mll01' in var.shortname:
		    outfile.write('yields total syst: %f \n' %(total_syst))


        #obtaining the histograms correctly normalized
        mclist, expected, observed, signal, _ = hists[category.name + ' ' + var.shortname]
        histograms = {}
        histograms['Expected']=expected
        for samp in tth2013.backgrounds:
                histograms[samp] = mclist[samp]
                #in case you have to add other histograms you maybe prefer to use the method clone:
                #histograms[samp] = mclist[samp].Clone(histname[samp])
        for samp in tth2013.observed:
                histograms[samp] = observed
        for samp in tth2013.signals:
                histograms[samp] = signal
        #print histograms

        for samp in histograms.keys():
                histograms[samp].SetNameTitle(histname[samp],'')
                histograms[samp].SetLineColor(histcolor[samp])

 	if 'Mll01' in var.shortname:
		print 'Category: ', category.name
		print 'Variable: ', var.shortname
		print 'Integral: '
		outfile.write('Category: %s \n' %(category.name))
		outfile.write('Variable: %s \n' %(var.shortname))
		outfile.write('Integral: \n')
		err=Double(0)#integral error
		value=0#integral value
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
 	if 'Mll01' in var.shortname:
		outfile.close()

