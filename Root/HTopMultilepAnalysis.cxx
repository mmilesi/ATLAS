/*********************************************
 *
 * The actual HTopMultilepAnalysis algorithm.
 * Here the user categorises events, and does 
 * the background estimation.
 *
 * M. Milesi (marco.milesi@cern.ch)
 * 
 *
 *********************************************/

// EL include(s):
#include <EventLoop/Job.h>
#include <EventLoop/Worker.h>
#include <EventLoop/OutputStream.h>

// EDM include(s): 
#include "xAODEventInfo/EventInfo.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODTau/TauJet.h"
#include "xAODTau/TauJetContainer.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/TruthParticle.h"
#include "xAODTruth/TruthVertex.h"

// package include(s):
#include "xAODAnaHelpers/HelperFunctions.h"
#include "xAODAnaHelpers/HelperClasses.h"
#include "xAODAnaHelpers/JetHists.h"
#include "xAODAnaHelpers/tools/ReturnCheck.h"
#include "xAODAnaHelpers/tools/ReturnCheckConfig.h"
#include "HTopMultilepAnalysis/HTopMultilepAnalysis.h"

// external tools include(s):
#include "TauAnalysisTools/TauSelectionTool.h"
#include "TauAnalysisTools/Enums.h"

// ROOT include(s):
#include "TFile.h"
#include "TTree.h"
#include "TEnv.h"
#include "TSystem.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

// this is needed to distribute the algorithm to the workers
ClassImp(HTopMultilepAnalysis)


HTopMultilepAnalysis :: HTopMultilepAnalysis () :
  m_cutflowHist(nullptr),
  m_cutflowHistW(nullptr),
  m_histEventCount(nullptr),
  /* For Francesco */
  m_totalEvents(nullptr),
  m_totalEventsW(nullptr),
  /* ************* */
  m_jetPlots(nullptr),
  m_TauSelTool(nullptr)
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
  
  Info("HTopMultilepAnalysis()", "Calling constructor");

  m_useCutFlow                = true; 
  
  m_inContainerName_Muons     = "";     
  m_inContainerName_Electrons = "";     
  m_inContainerName_Leptons   = "";
  m_inContainerName_Jets      = "";
  m_inContainerName_Taus      = "";      

  m_doLHPIDCut                = true;  
  m_doCutBasedPIDCut          = false;    
}

EL::StatusCode  HTopMultilepAnalysis :: configure ()
{
 
  if ( !getConfig().empty() ) {

    // read in user configuration from text file
    TEnv *config = new TEnv(getConfig(true).c_str());
    if ( !config ) {
      Error("BasicEventSelection()", "Failed to initialize reading of config file. Exiting." );
      return EL::StatusCode::FAILURE;
    }
    Info("configure()", "Configuing HTopMultilepAnalysis Interface. User configuration read from : %s", getConfig().c_str());
    
    // read debug flag from .config file
    m_debug	                 = config->GetValue("Debug" ,      m_debug);
    m_useCutFlow                 = config->GetValue("UseCutFlow",  m_useCutFlow);
    
    // input container to be read from TEvent or TStore
    m_inContainerName_Electrons	 = config->GetValue("InputContainerElectrons", m_inContainerName_Electrons.c_str());
    m_inContainerName_Muons	 = config->GetValue("InputContainerMuons",     m_inContainerName_Muons.c_str());
    m_inContainerName_Leptons    = config->GetValue("InputContainerLeptons",   m_inContainerName_Leptons.c_str());
    m_inContainerName_Jets	 = config->GetValue("InputContainerJets",      m_inContainerName_Jets.c_str());
    m_inContainerName_Taus	 = config->GetValue("InputContainerTaus",      m_inContainerName_Taus.c_str());
 
    // electron ID stuff - choose which one to define "Tight" electrons
    m_doLHPIDCut                 = config->GetValue("DoLHPIDCut"          ,  m_doLHPIDCut );      
    m_doCutBasedPIDCut           = config->GetValue("DoCutBasedPIDCut"    ,  m_doCutBasedPIDCut );
 
    config->Print();
  
    Info("configure()", "HTopMultilepAnalysis Interface succesfully configured!");
  
    delete config; config = nullptr;
  }

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode HTopMultilepAnalysis :: setupJob (EL::Job& job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.
  
  job.useXAOD();
  xAOD::Init( "HTopMultilepAnalysis" ).ignore(); // call before opening first file
  
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepAnalysis :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.
  Info("histInitialize()", "Calling histInitialize");

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepAnalysis :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepAnalysis :: changeInput (bool /*firstFile*/)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepAnalysis :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.
  
  m_event = wk()->xaodEvent();
  m_store = wk()->xaodStore();
  
  // count number of events
  m_eventCounter = 0;

  if ( this->configure() == EL::StatusCode::FAILURE ) {
    Error("initialize()", "Failed to properly configure. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  if ( m_useCutFlow ) {
    TFile *fileCF  = wk()->getOutputFile ("cutflow");  
    m_cutflowHist  = (TH1D*)fileCF->Get("cutflow");
    m_cutflowHistW = (TH1D*)fileCF->Get("cutflow_weighted");
    // label the bins for the cutflow
    m_cutflow_bin     = m_cutflowHist->GetXaxis()->FindBin(m_name.c_str());
    // do it again for the weighted cutflow hist
    m_cutflowHistW->GetXaxis()->FindBin(m_name.c_str());
  }

  TFile *fileMD = wk()->getOutputFile ("metadata");  
  m_histEventCount  = (TH1D*)fileMD->Get("MetaData_EventCount");
  if ( !m_histEventCount ) { 
    Error("histInitialize()", "Failed to retrieve MetaData histogram. Aborting");
    return EL::StatusCode::FAILURE;
  }

  /* for Francesco */
  m_totalEvents  = new TH1D("TotalEvents",  "TotalEvents",  2, 1, 3);    
  m_totalEventsW = new TH1D("TotalEventsW", "TotalEventsW", 2, 1, 3);  
  wk() -> addOutput(m_totalEvents);
  wk() -> addOutput(m_totalEventsW);
  
  m_jetPlots = new JetHists( "highPtJets", "clean" ); // second argument: "kinematic", "clean", "energy", "resolution"
  m_jetPlots -> initialize();
  m_jetPlots -> record( wk() );
  
  // initialise TauSelectionTool                                                                                                                        
  m_TauSelTool = new TauAnalysisTools::TauSelectionTool( "TauSelectionTool" );
  m_TauSelTool->msg().setLevel( MSG::INFO );
  m_TauSelTool->setProperty( "SelectionCuts", int(TauAnalysisTools::CutJetIDWP) );
  // requiring tau to pass a jet BDT working point                  
  TauAnalysisTools::e_JETID tauID = TauAnalysisTools::JETIDBDTTIGHT;
  m_TauSelTool->setProperty( "JetIDWP", static_cast<int>( tauID ) ); // need to cast enum to int!                                                        
  RETURN_CHECK( "HTopMultilepAnalysis::initialize()", m_TauSelTool->initialize(), "Failed to properly initialize TauSelectionTool." );

  Info("initialize()", "HTopMultilepAnalysis Interface succesfully initialized!" );
  
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode HTopMultilepAnalysis :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.
  if ( m_debug ) { Info("execute()", "Applying analysis selection"); }
  
  /* for Francesco */
  double n_init_evts(0.);   
  double n_init_evts_W(0.);
  // if MetaData is not empty, use it!
  if ( m_histEventCount->GetBinContent(1) > 0 && m_histEventCount->GetBinContent(3) > 0 ) {
    n_init_evts   =  m_histEventCount->GetBinContent(1);   //nEvents initial
    n_init_evts_W =  m_histEventCount->GetBinContent(3);  // sumOfWeights initial
  } 
  // ...else, retrieve it from cutflow
  else 
  {
    int init_evts_bin    =  m_cutflowHist->GetXaxis()->FindBin("all");
    n_init_evts          =  m_cutflowHist->GetBinContent( init_evts_bin );
    int init_evts_bin_W  =  m_cutflowHistW->GetXaxis()->FindBin("all");
    n_init_evts_W        =  m_cutflowHistW->GetBinContent( init_evts_bin_W );  
  }
  // set the value into both bins of histogram
  m_totalEvents->SetBinContent( 1, n_init_evts );
  m_totalEvents->SetBinContent( 2, n_init_evts );  
  m_totalEventsW->SetBinContent( 1, n_init_evts_W );
  m_totalEventsW->SetBinContent( 2, n_init_evts_W );   
  
  //---------------------------
  //***** Event information
  //--------------------------- 
  
  // retrieve event info
  const xAOD::EventInfo* eventInfo(nullptr);
  RETURN_CHECK("HTopMultilepAnalysis::execute()", HelperFunctions::retrieve(eventInfo, "EventInfo", m_event, m_store, m_debug) , "");

  // bool isMC = ( eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) );
  
  // retrieve vertices  
  const xAOD::VertexContainer* vertices(nullptr);
  RETURN_CHECK("HTopMultilepAnalysis::execute()", HelperFunctions::retrieve(vertices, "PrimaryVertices", m_event, m_store,  m_debug) , "");

  // retrieve selected objects
  const xAOD::ElectronContainer* signalElectrons(nullptr); 
  RETURN_CHECK("HTopMultilepAnalysis::execute()", HelperFunctions::retrieve(signalElectrons, m_inContainerName_Electrons, m_event, m_store,  m_debug) , "");
  const xAOD::MuonContainer*     signalMuons(nullptr);    
  RETURN_CHECK("HTopMultilepAnalysis::execute()", HelperFunctions::retrieve(signalMuons, m_inContainerName_Muons, m_event, m_store, m_debug ) , "");
  const xAOD::JetContainer*      signalJets(nullptr);   
  RETURN_CHECK("HTopMultilepAnalysis::execute()", HelperFunctions::retrieve(signalJets, m_inContainerName_Jets, m_event, m_store,  m_debug) , "");
  const xAOD::TauJetContainer*   signalTauJets(nullptr);
  RETURN_CHECK("HTopMultilepAnalysis::execute()", HelperFunctions::retrieve(signalTauJets,  m_inContainerName_Taus, m_event, m_store, m_debug) , "");
  ConstDataVector<xAOD::IParticleContainer>* leptonsCDV(nullptr);
  RETURN_CHECK("HTopMultilepAnalysis::execute()", HelperFunctions::retrieve(leptonsCDV, m_inContainerName_Leptons, m_event, m_store, m_debug) ,"");
  // Make a sorted version of the lepton container 
  // (this can be on the stack! Will not be pushed to the store...)
  const xAOD::IParticleContainer leptonsSorted = HelperFunctions::sort_container_pt( leptonsCDV->asDataVector() );

  unsigned int nSignalLeptons   = leptonsCDV->size();
  unsigned int nSignalJets      = signalJets->size();
  unsigned int nSignalTaus      = signalTauJets->size();
  static SG::AuxElement::Accessor< unsigned int > nBjetsMediumAcc("nBjetsMedium");
  unsigned int nBjetsMedium(-1);
  if ( nBjetsMediumAcc.isAvailable( *eventInfo ) ) {
     nBjetsMedium = nBjetsMediumAcc( *eventInfo );
  } else {
    Error("execute()"," nBjetsMedium is not available. Aborting" ); 
    return EL::StatusCode::FAILURE;
  }
    
  if ( m_debug ) {

    // retrieve preselected objects ( only for debugging purposes )
    const xAOD::ElectronContainer* preselElectrons(nullptr);
    RETURN_CHECK("HTopMultilepAnalysis::execute()", HelperFunctions::retrieve(preselElectrons, "Electrons_PreSelected", m_event, m_store, m_debug) , "");
    const xAOD::MuonContainer*     preselMuons(nullptr);     
    RETURN_CHECK("HTopMultilepAnalysis::execute()", HelperFunctions::retrieve(preselMuons, "Muons_PreSelected", m_event, m_store, m_debug) , "");
    const xAOD::JetContainer*      preselJets(nullptr); 
    RETURN_CHECK("HTopMultilepAnalysis::execute()", HelperFunctions::retrieve(preselJets, "AntiKt4EMTopoJets_Selected", m_event, m_store, m_debug ) , "");

    unsigned int nPreselElectrons = preselElectrons->size();
    unsigned int nPreselMuons     = preselMuons->size();
    unsigned int nPreselJets      = preselJets->size();
    unsigned int nSignalElectrons = signalElectrons->size();
    unsigned int nSignalMuons     = signalMuons->size();
    Info("execute()","Event %i ", static_cast<int>(m_eventCounter));
    Info("execute()"," Preselected vs Selected Signal Jets: \t %u \t %u "      , nPreselJets, nSignalJets );
    Info("execute()"," Preselected vs Selected Signal Muons: \t %u \t %u  "    , nPreselMuons, nSignalMuons ); 
    Info("execute()"," Preselected vs Selected Signal Electrons: \t %u \t %u " , nPreselElectrons, nSignalElectrons );     
    Info("execute()"," Selected Signal Leptons: \t %u " , nSignalLeptons );   
    Info("execute()"," Selected Signal Taus: \t %u " , nSignalTaus );   

  }

  ++m_eventCounter;
  
  //-------------------------------
  //***** Retrieve event weight
  //-------------------------------  
  
  // retrieve event weight from eventInfo (NB: will be always 1 for Data - see xAODAnaHelpers/root/BasicEventSelection.cxx)
  // for MC, it includes also PU weight. Need to multiply it by object SF weights later on! 
  float mcEvtWeight(1.0); 
  static SG::AuxElement::Accessor< float > mcEvtWeightAcc("mcEventWeight");
  if ( !mcEvtWeightAcc.isAvailable(*eventInfo) ) {
    Error("execute()", "event weight is not available. Aborting ");
    return EL::StatusCode::FAILURE;
  } 
  mcEvtWeight = mcEvtWeightAcc( *eventInfo );
   
  // multiply SFs to event weight 
  /*
  if ( isMC ) {
    // electron efficiency SF
    static SG::AuxElement::Accessor< float > elEffSFAcc("SF");
    for ( auto el_itr : *(signalElectrons) ) {  
      if ( !elEffSFAcc.isAvailable(*el_itr) ) { continue; }
      mcEvtWeight *= elEffSFAcc(*el_itr);
    }   

    // muon efficiency SF
    // bTagging efficiency SF
    // trigger efficiency SF
  }
  */

  if ( m_debug ) {
    Info("execute()", "event weight = %2f. ", mcEvtWeight );
  }

  //-------------------------------
  //***** Fill cutflow histogram
  //-------------------------------  

  
  if ( m_useCutFlow ) {
    m_cutflowHist ->Fill( m_cutflow_bin, 1 );
    m_cutflowHistW->Fill( m_cutflow_bin, mcEvtWeight);
  }
    
  //-------------------------------
  //***** decorate selected events
  //-------------------------------
  
  // declare static event decorators
  static SG::AuxElement::Decorator< float >  ystarDecor("ystar");

  // now decorate event!
  ystarDecor(*eventInfo) = ( signalJets->size() > 1 ) ? ( signalJets->at(0)->rapidity() - signalJets->at(1)->rapidity() ) / 2.0 : 999.0;

  // decorate taus
  static SG::AuxElement::Decorator< char > isTauBDTTightDecor("isTauBDTTight");  
  for ( auto tau_itr : *signalTauJets ) { 
    isTauBDTTightDecor( *tau_itr ) = 0;
    if ( m_TauSelTool->accept(*tau_itr) ) { isTauBDTTightDecor( *tau_itr ) = 1; }
  }

  // accessor to lepton isolation flag 
  static SG::AuxElement::Accessor< char > isIsoAcc("isIsolated");
  // accessor to likelihood PID for electrons
  static SG::AuxElement::Accessor< char > LHTightAcc("LHTight");
   // accessor to cut-based PID for electrons
  static SG::AuxElement::Accessor< char > EMTightAcc("Tight");

  // decorator for "Tight" leptons  
  static SG::AuxElement::Decorator< char > isTightDecor("isTight"); // electrons: isolated + tight PID (tightPP or TightLH, depending on user's choice)
  								    // muons: isolated + d0sig < 3.0

  for ( auto el_itr : *(signalElectrons) ) {

    // default
    isTightDecor( *el_itr ) = 0; 

    if ( !isIsoAcc.isAvailable( *el_itr ) ) {
      Error("execute()", "isIsolated attribute is not available for this electron. Aborting ");
      return EL::StatusCode::FAILURE;
    } 
    if ( !LHTightAcc.isAvailable( *el_itr ) ) {
      Error("execute()", "LHTight attribute is not available for this electron. Aborting ");
      return EL::StatusCode::FAILURE;
    } 
    if ( !EMTightAcc.isAvailable( *el_itr ) ) {
      Error("execute()", "Tight attribute is not available for this electron. Aborting ");
      return EL::StatusCode::FAILURE;
    } 

    // flag the "Tight" electrons
    if (  isIsoAcc( *el_itr ) ) {
        if      ( m_doLHPIDCut       && LHTightAcc( *el_itr )  ) { isTightDecor( *el_itr ) = 1; }
	else if ( m_doCutBasedPIDCut && EMTightAcc( *el_itr ) )  { isTightDecor( *el_itr ) = 1; }
    }
 
  }
  
  for ( auto mu_itr : *(signalMuons) ) {

      // default
      isTightDecor( *mu_itr ) =  0;

      if ( !isIsoAcc.isAvailable( *mu_itr ) ) {
	Error("execute()", "isIsolated attribute is not available for this muon. Aborting ");
	return EL::StatusCode::FAILURE;
      } 

      // flag the "Tight" muons
      if ( isIsoAcc( *mu_itr) ) {
          const xAOD::TrackParticle* tp  = mu_itr->primaryTrackParticle();
          float d0_significance = fabs( tp->d0() ) / sqrt( tp->definingParametersCovMatrix()(0,0) );	  
          if  ( d0_significance < 3.0 ) { isTightDecor( *mu_itr ) = 1; }
      }	
  }
  
  //-------------------------------- 
  //***** histogram filling
  //--------------------------------
  
  // fill plots for all preselected jets
  //m_jetPlots->FillHistograms( jets, mcEvtWeight );

  //  for( auto jet_itr : *(preselJets) ) {
  //    // fill plots for a select set of jets - those > 50 GeV
  //   if ( jet_itr->pt() > 50e3 ) {
  //      m_jetPlots->execute( jet_itr, mcEvtWeight );
  //    }
  //  }
 
  //------------------------------------------- 
  //***** set T&P variables in r/f rate meas CR
  //------------------------------------------- 
  
  if ( nSignalLeptons == 2 && ( nSignalJets >= 1 && nSignalJets <= 3 ) &&  nBjetsMedium >= 1 ) {
    this->applyTagAndProbeRFRateMeasurement( eventInfo, leptonsSorted );
  }
  
  //----------------------------------- 
  //***** Matrix Method event weighting
  //-----------------------------------

  // first, decorate event with specific variables ...
  this->addChannelDecorations( eventInfo, leptonsSorted );
  
  // ...then, decorate event with MM and FF weight!
  this->fakeWeightCalculator( eventInfo, leptonsSorted ); 
  
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode HTopMultilepAnalysis :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepAnalysis :: finalize ()
{
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.  This is different from histFinalize() in that it only
  // gets called on worker nodes that processed input events.

  Info("finalize()", "Deleting tool instances...");

  if ( m_TauSelTool ) { delete m_TauSelTool; m_TauSelTool = nullptr; }    
   
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepAnalysis :: histFinalize ()
{
  // This method is the mirror image of histInitialize(), meaning it
  // gets called after the last event has been processed on the worker
  // node and allows you to finish up any objects you created in
  // histInitialize() before they are written to disk.  This is
  // actually fairly rare, since this happens separately for each
  // worker node.  Most of the time you want to do your
  // post-processing on the submission node after all your histogram
  // outputs have been merged.  This is different from finalize() in
  // that it gets called on all worker nodes regardless of whether
  // they processed input events.
  return EL::StatusCode::SUCCESS;
}
  
//***************************************************
//
// Set Tag&Probe variables in r/f rate measurement CR 
//
//
EL::StatusCode HTopMultilepAnalysis :: applyTagAndProbeRFRateMeasurement( const xAOD::EventInfo* eventInfo, const xAOD::IParticleContainer& leptons )
{
  
  // ----------------------------------------
  // Categorise event based on lepton flavour
  // ----------------------------------------
  
  int prod_lep_charge(1);
  unsigned int count_el(0), count_mu(0);
  for ( auto lep_it : leptons ) {
      // get the lepton flavour
      xAOD::Type::ObjectType leptonFlavour = lep_it->type();
      if ( leptonFlavour == xAOD::Type::Electron ) { 
        ++count_el; 
        prod_lep_charge *= dynamic_cast<const xAOD::Electron*>( lep_it )->charge();
      } else if ( leptonFlavour == xAOD::Type::Muon ) { 
        ++count_mu; 
        prod_lep_charge *= dynamic_cast<const xAOD::Muon*>( lep_it )->charge();
      }
  }
  bool isElEl = ( count_el == 2 && count_mu == 0 );
  bool isElMu = ( count_el == 1 && count_mu == 1 );
  bool isMuMu = ( count_el == 0 && count_mu == 2 );
  
  bool isSS = (   prod_lep_charge > 0  );  
  bool isOS = ( !(prod_lep_charge > 0) ); 

  // --------------------------------
  // Declare accessors and decorators
  // --------------------------------
  
  // decorate lepton w/ is tag/probe
  static SG::AuxElement::Decorator< char > isTagDecor("isTag");
  // declare an event decoration in case there are no "good" leptons for the rate calculation
  static SG::AuxElement::Decorator< char > isNonTightEventDecor("isNonTightEvent");
  // declare event decorations for checking the probe type in event
  static SG::AuxElement::Decorator< char > isProbeElEventDecor( "isProbeElEvent" ); 
  static SG::AuxElement::Decorator< char > isProbeMuEventDecor( "isProbeMuEvent" ); 

  // accessor to tight leptons 
  static SG::AuxElement::Accessor< char > isTightAcc("isTight");
  // accessor to trigger-matched leptons 
  static SG::AuxElement::Accessor< char > isTrigMatchedAcc("isTrigMatched");
  // accessor to tag leptons 
  static SG::AuxElement::Accessor< char > isTagAcc("isTag");
    
  // decorate with default values
  for ( auto lep_it : leptons ) { isTagDecor( *lep_it ) = 0; }
  isNonTightEventDecor( *eventInfo ) = 0;
  isProbeElEventDecor( *eventInfo )  = 0;
  isProbeMuEventDecor( *eventInfo )  = 0;
    
  // -------------------------------------------
  // Now, take the leading and subleading lepton 
  // -------------------------------------------
    
  const xAOD::IParticle* leadingLepton           = *(leptons.begin());
  xAOD::Type::ObjectType leadingLeptonFlavour    = leadingLepton->type();
  bool                   isLeadingTight(false);          
  if ( isTightAcc.isAvailable( *leadingLepton ) ) {
    isLeadingTight = isTightAcc( *leadingLepton );
  } else {
    Warning("applyTagAndProbeRFRateMeasurement()","SG::AuxElement::Accessor('isTight') is not available for the leading lepton. Should not happen. Assigning 'isTight' = false" );
  }
  bool 			 isLeadingTrigMatched(false);  
  
  const xAOD::IParticle* subLeadingLepton        = *(std::next(leptons.begin(),1));
  xAOD::Type::ObjectType subLeadingLeptonFlavour = subLeadingLepton->type();
  bool                   isSubLeadingTight (false);
  if ( isTightAcc.isAvailable( *subLeadingLepton ) ) {
    isLeadingTight = isTightAcc( *subLeadingLepton );
  } else {
    Warning("applyTagAndProbeRFRateMeasurement()","SG::AuxElement::Accessor('isTight') is not available for the subleading lepton. Should not happen. Assigning 'isTight' = false" );
  }
  bool 			 isSubLeadingTrigMatched(false);  
    
  if ( m_debug ) { Info("applyTagAndProbeRFRateMeasurement()","Leading lepton: isTight: %i \t isTrigMatched: %i \t isTag: %i    ", isLeadingTight, isLeadingTrigMatched,  isTagAcc( *leadingLepton ) ); }
  if ( m_debug ) { Info("applyTagAndProbeRFRateMeasurement()","Subleading lepton: isTight: %i \t isTrigMatched: %i \t isTag: %i ", isSubLeadingTight, isSubLeadingTrigMatched, isTagAcc( *subLeadingLepton ) ); }
  
  // --------------------------------------------
  // Now decide who's the tag and who's the probe
  // --------------------------------------------
  
  // generate a uniformly distributed random number in [0,1]			  
  TRandom3 rndm(0);
  float unif = rndm.Uniform(); 
  
  if ( isElEl || isMuMu ) {  
    
    //In the real case (i.e. OS SF leptons) we ask randomly if either the leading lep or the subleading lep satisfy the tag requirement otherwise we have trigger pt treshold problem. 
    //In the fake case (i.e. SS SF leptons) if both leptons are matched, the one with smaller pt is the probe because it has been seen that the low pt lepton is the fake. 
			   
    // 1. leading & subleading lepton are both tight & trigger-matched 
    //    OR
    // 2. leading lepton is tight & trigger-matched, but subleading is not  
    if ( ( isLeadingTight && isSubLeadingTight /* && isLeadingTrigMatched && isSubLeadingTrigMatched */ ) || 
         ( isLeadingTight /* && isLeadingTrigMatched */ ) 
       ) {
      
      if ( isSS ) {
	
	// the probe will be the subleading lepton, no matter what
	isTagDecor( *leadingLepton ) = 1;
	
	// categorise event
	isProbeElEventDecor( *eventInfo ) = ( subLeadingLeptonFlavour == xAOD::Type::Electron );
	isProbeMuEventDecor( *eventInfo ) = ( subLeadingLeptonFlavour == xAOD::Type::Muon );
	
      
      } else if ( isOS ) {
        
	// assign randomly who's tag and who's probe 
	if ( unif >= 0.5 ) {
	    
	    isTagDecor( *leadingLepton )    = 1;
	    
	    // categorise event
	    isProbeElEventDecor( *eventInfo ) = ( subLeadingLeptonFlavour == xAOD::Type::Electron );
	    isProbeMuEventDecor( *eventInfo ) = ( subLeadingLeptonFlavour == xAOD::Type::Muon );
	    
	} else {
	    
	    isTagDecor( *subLeadingLepton ) = 1;
	    
	    // categorise event
	    isProbeElEventDecor( *eventInfo ) = ( leadingLeptonFlavour == xAOD::Type::Electron );
	    isProbeMuEventDecor( *eventInfo ) = ( leadingLeptonFlavour == xAOD::Type::Muon );
	
	}
	
      }
      
    }
    // 3. subleading lepton is tight & trigger-matched, but leading is not 
    else if ( isSubLeadingTight /* && isSubLeadingTrigMatched  */ ) {
      
      if ( isSS ) {
	
	// the probe will be the leading lepton, which is also !tight
	isTagDecor( *subLeadingLepton ) = 1;
	
	// categorise event
	isProbeElEventDecor( *eventInfo ) = ( leadingLeptonFlavour == xAOD::Type::Electron );
	isProbeMuEventDecor( *eventInfo ) = ( leadingLeptonFlavour == xAOD::Type::Muon );
     
      } else if ( isOS ) {
        
	// assign randomly who's tag and who's probe 
        if ( unif >= 0.5 ) {
	
	    isTagDecor( *leadingLepton )    = 1;
	    
	    // categorise event
	    isProbeElEventDecor( *eventInfo ) = ( subLeadingLeptonFlavour == xAOD::Type::Electron );
	    isProbeMuEventDecor( *eventInfo ) = ( subLeadingLeptonFlavour == xAOD::Type::Muon );
	    
	} else {
	    
	    isTagDecor( *subLeadingLepton ) = 1;
	    
	    // categorise event
	    isProbeElEventDecor( *eventInfo ) = ( leadingLeptonFlavour == xAOD::Type::Electron );
	    isProbeMuEventDecor( *eventInfo ) = ( leadingLeptonFlavour == xAOD::Type::Muon );
	    
	}
      
      }
    
    } 
     // 4. neither the leading, nor the subleading lepton are tight & trigger-matched  
    else {
       // take note that this event has no Tight leptons
       isNonTightEventDecor( *eventInfo ) = 1;
       
       // the probe will be the subleading lepton
       isTagDecor( *leadingLepton )    = 1;
	    
       // categorise event
       isProbeElEventDecor( *eventInfo ) = ( subLeadingLeptonFlavour == xAOD::Type::Electron );
       isProbeMuEventDecor( *eventInfo ) = ( subLeadingLeptonFlavour == xAOD::Type::Muon );
       
    }
     
  } // end check ( isElEl || isMuMu )
  else if ( isElMu ) 
  {
    
    // 1. leading & subleading lepton are both tight & trigger-matched 
    //   need to distinguish assignment btwn/ SS and OS  
    if ( isLeadingTight && isSubLeadingTight /* && isLeadingTrigMatched && isSubLeadingTrigMatched  */ ) {
        
	// In the fake case (i.e. SS leptons), if both leptons are matched, I must chose which lepton is the tag and which one is the probe 
        // because I want the probe to be the fake and the real the tag. 
        // I chose the low pt one as probe. 
        // In the real case (OS leptons), the leading is the tag and subleading is the probe, no matter what.
        
	// ---> no need to separate OS/SS here!  

        isTagDecor( *leadingLepton )    = 1;
	
        // categorise event
        isProbeElEventDecor( *eventInfo ) = ( subLeadingLeptonFlavour == xAOD::Type::Electron );
        isProbeMuEventDecor( *eventInfo ) = ( subLeadingLeptonFlavour == xAOD::Type::Muon );

    }
    // 2. leading lepton is tight & trigger-matched, but subleading is not 
    else if ( isLeadingTight /* && isLeadingTrigMatched */ ) {
        
	// the probe will be the subleading, which is also !tight
	isTagDecor( *leadingLepton )    = 1;
	
        // categorise event
        isProbeElEventDecor( *eventInfo ) = ( subLeadingLeptonFlavour == xAOD::Type::Electron );
        isProbeMuEventDecor( *eventInfo ) = ( subLeadingLeptonFlavour == xAOD::Type::Muon );
	
    }
    // 3. subleading lepton is tight & trigger-matched, but leading is not 
    else if ( isSubLeadingTight /* && isSubLeadingTrigMatched  */ ) {
        
	// the probe will be the leading, which is also !tight
        isTagDecor( *subLeadingLepton ) = 1;
	
        // categorise event
        isProbeElEventDecor( *eventInfo ) = ( leadingLeptonFlavour == xAOD::Type::Electron );
        isProbeMuEventDecor( *eventInfo ) = ( leadingLeptonFlavour == xAOD::Type::Muon );

    }
     // 4. neither the leading, nor the subleading lepton are tight & trigger-matched  
    else {
       // take note that this event has no Tight leptons
       isNonTightEventDecor( *eventInfo ) = 1;
       
       // the probe will be the subleading lepton
       isTagDecor( *leadingLepton )    = 1;
       
       // categorise event
       isProbeElEventDecor( *eventInfo ) = ( subLeadingLeptonFlavour == xAOD::Type::Electron );
       isProbeMuEventDecor( *eventInfo ) = ( subLeadingLeptonFlavour == xAOD::Type::Muon );
       
    }
    
  } // end check isElMu

  
  if ( m_debug ) {
    Info("execute()","Checking *isTag* lepton decoration"); 
    for ( auto lep_it : leptons ) {
      Info("execute()","\t lepton isTag: %i ", isTagAcc( *lep_it ) ); 
    }
  }
   
  return EL::StatusCode::SUCCESS;

}

//*****************************************************************************
//
// Channel variables
//
//
EL::StatusCode HTopMultilepAnalysis ::  addChannelDecorations(const xAOD::EventInfo* eventInfo, const xAOD::IParticleContainer& leptons)
{
  
  // declare new event decorations
  static SG::AuxElement::Decorator< char > isSS01Decor("isSS01");  // in the dilepton category, this will identify the events where the 2 leptons are SS
  static SG::AuxElement::Decorator< char > isSS12Decor("isSS12");  // in the trilepton category, this will identify the events where 2 of the leptons are SS (other will be OS by construction)
  // declare new lepton object decorations
  static SG::AuxElement::Decorator< char > isOSlepDecor("isOSlep");               // in the trilepton category, this will identify the lepton which is OS w.r.t. the other two
  static SG::AuxElement::Decorator< char > isClosestSSlepDecor("isClosestSSlep"); // in the trilepton category, this will identify the SS lepton which is closest in deltaR to the OS lepton
  
  // the invariant masses of the combined lepton systems 
  static SG::AuxElement::Decorator< float > mll01Decor("mll01"); // lead + 2nd lead  
  static SG::AuxElement::Decorator< float > mll02Decor("mll02"); // lead + 3rd lead 
  static SG::AuxElement::Decorator< float > mll12Decor("mll12"); // 2nd lead + 3rd lead 
  static SG::AuxElement::Decorator< float > mlll012Decor("mlll012");   

  // start decorating with default values 
  isSS01Decor(*eventInfo)  = 0; // false
  isSS12Decor(*eventInfo)  = 0; // false
  mll01Decor(*eventInfo)   = -1.0;   
  mll02Decor(*eventInfo)   = -1.0;   
  mll12Decor(*eventInfo)   = -1.0;   
  mlll012Decor(*eventInfo) = -1.0; 

  unsigned int nLeptons = leptons.size();

  if ( nLeptons > 0 ) {
    for( auto lep_it : leptons ) {
      isOSlepDecor(*(lep_it)) = 0; // false
      isSS12Decor(*(lep_it))  = 0; // false
    }
  }
  
  // initialize TLorentzVectors (for invariant mass computation)
  xAOD::IParticle::FourMom_t lepA_4mom, lepB_4mom, lepC_4mom;
 
  if ( nLeptons == 2 )
  {
     // retrieve lepA 
     const xAOD::IParticle* lepA = *(leptons.begin());
     // retrieve lepB
     const xAOD::IParticle* lepB = *(std::next( leptons.begin(), 1 ));     
     
     // compute invariant mass of the pair
     lepA_4mom.SetPtEtaPhiM( lepA->pt(), lepA->eta(), lepA->phi(), lepA->m() );
     lepB_4mom.SetPtEtaPhiM( lepB->pt(), lepB->eta(), lepB->phi(), lepB->m() );
     xAOD::IParticle::FourMom_t pair = lepA_4mom + lepB_4mom;
     mll01Decor(*eventInfo) = pair.M();
     
     // check if the two leptons are SS, and decorate event
     int prod_lep_charge(1);
     // need to cast to derived xAOD type as xAOD::IParticle does not have charge info...
     xAOD::Type::ObjectType lepAFlavour = lepA->type();
     xAOD::Type::ObjectType lepBFlavour = lepB->type();
     prod_lep_charge *=  ( lepAFlavour == xAOD::Type::Electron  ) ?  ( (dynamic_cast<const xAOD::Electron*>(lepA))->charge() ) : ( (dynamic_cast<const xAOD::Muon*>(lepA))->charge() ) ;
     prod_lep_charge *=  ( lepBFlavour == xAOD::Type::Electron  ) ?  ( (dynamic_cast<const xAOD::Electron*>(lepB))->charge() ) : ( (dynamic_cast<const xAOD::Muon*>(lepB))->charge() ) ;
     
     // now decorate event!
     isSS01Decor(*eventInfo) = ( prod_lep_charge > 0 );
     
     // FIXME : decorate w/ MET, mT	
     /*
     TLorentzVector met_4mom, lepton_4mom;
     met_4mom.SetPtEtaPhiM(ao->eventselections.met_et, 0, ao->eventselections.met_phi,0);
     double mt;// Transverse mass
     lepton_4mom.SetPtEtaPhiM(ao->channelvariables.lepAPt, ao->channelvariables.lepAEta, ao->channelvariables.lepAPhi, ao->channelvariables.lepAM);
     mt = sqrt( 2*lepton_4mom.Pt()*(met_4mom.Pt())*(1.-cos(lepton_4mom.DeltaPhi(met_4mom))) );
     ao->channelvariables.mTLep0 = mt;
     lepton_4mom.SetPtEtaPhiM(ao->channelvariables.Lep1Pt, ao->channelvariables.Lep1Eta, ao->channelvariables.Lep1Phi, ao->channelvariables.Lep1M);
     mt = sqrt( 2*lepton_4mom.Pt()*(met_4mom.Pt())*(1.-cos(lepton_4mom.DeltaPhi(met_4mom))) );
     ao->channelvariables.mTLep1 = mt;
     */
  }
  else if ( nLeptons == 3 )
  {
     // retrieve lepA 
     const xAOD::IParticle* lepA = *(leptons.begin());
     // retrieve lepB
     const xAOD::IParticle* lepB = *(std::next( leptons.begin(), 1 ));     
     // retrieve lepC
     const xAOD::IParticle* lepC = *(std::next( leptons.begin(), 2 ));   
     
     // compute invariant mass of all the possible pairs, and of the triplet as well
     lepA_4mom.SetPtEtaPhiM( lepA->pt(), lepA->eta(), lepA->phi(), lepA->m() );
     lepB_4mom.SetPtEtaPhiM( lepB->pt(), lepB->eta(), lepB->phi(), lepB->m() );
     lepC_4mom.SetPtEtaPhiM( lepC->pt(), lepC->eta(), lepC->phi(), lepC->m() );     
     xAOD::IParticle::FourMom_t pairAB    = lepA_4mom + lepB_4mom;
     mll01Decor(*eventInfo)   = pairAB.M();
     xAOD::IParticle::FourMom_t pairAC    = lepA_4mom + lepC_4mom;
     mll02Decor(*eventInfo)   = pairAC.M();  
     xAOD::IParticle::FourMom_t pairBC    = lepB_4mom + lepC_4mom;
     mll12Decor(*eventInfo)   = pairBC.M();       
     xAOD::IParticle::FourMom_t triplet   = lepA_4mom + lepB_4mom + lepC_4mom;
     mlll012Decor(*eventInfo) = triplet.M(); 
            
     // retrieve charge
     int lepAcharge(0), lepBcharge(0), lepCcharge(0);
     // need to dynamic_cast to derived xAOD type as xAOD::IParticle does not have charge info...
     xAOD::Type::ObjectType lepAFlavour = lepA->type();
     xAOD::Type::ObjectType lepBFlavour = lepB->type();
     xAOD::Type::ObjectType lepCFlavour = lepC->type();
     lepAcharge =  ( lepAFlavour == xAOD::Type::Electron  ) ?  ( (dynamic_cast<const xAOD::Electron*>(lepA))->charge() ) : ( (dynamic_cast<const xAOD::Muon*>(lepA))->charge() ) ;
     lepBcharge =  ( lepBFlavour == xAOD::Type::Electron  ) ?  ( (dynamic_cast<const xAOD::Electron*>(lepB))->charge() ) : ( (dynamic_cast<const xAOD::Muon*>(lepB))->charge() ) ;
     lepCcharge =  ( lepCFlavour == xAOD::Type::Electron  ) ?  ( (dynamic_cast<const xAOD::Electron*>(lepC))->charge() ) : ( (dynamic_cast<const xAOD::Muon*>(lepC))->charge() ) ;
     
     // look only at events where there is an OS pair
     if ( fabs( lepAcharge/fabs(lepAcharge) + lepBcharge/fabs(lepBcharge) + lepCcharge/fabs(lepCcharge)  ) != 3 )
     {
        // now decorate event!
	isSS12Decor(*eventInfo) = 1;
     
        // now find the OS lepton and decorate! 
        isOSlepDecor(*lepA) = ( (lepBcharge * lepCcharge) > 0 ); 
        isOSlepDecor(*lepB) = ( (lepAcharge * lepCcharge) > 0 ); 
        isOSlepDecor(*lepC) = ( (lepAcharge * lepBcharge) > 0 ); 
  
        // once OS lepton has been found, find the SS lep with min{ deltaR(lep0) } and decorate!
	float dR_i(-1), dR_j(-1);
	if ( isOSlepDecor(*lepA) ) {
	   dR_i = lepA_4mom.DeltaR(lepB_4mom ); 
	   dR_j = lepA_4mom.DeltaR(lepC_4mom );
	   isClosestSSlepDecor(*lepB) = ( dR_i <=  dR_j );
	   isClosestSSlepDecor(*lepC) = ( dR_j  <  dR_i );
	}
	else if ( isOSlepDecor(*lepB) ) {
	   dR_i = lepB_4mom.DeltaR(lepA_4mom ); 
	   dR_j = lepB_4mom.DeltaR(lepC_4mom );
	   isClosestSSlepDecor(*lepA) = ( dR_i <=  dR_j );
	   isClosestSSlepDecor(*lepC) = ( dR_j  <  dR_i );
	}
	else if ( isOSlepDecor(*lepC) ) {
	   dR_i = lepC_4mom.DeltaR(lepA_4mom ); 
	   dR_j = lepC_4mom.DeltaR(lepB_4mom );
	   isClosestSSlepDecor(*lepA) = ( dR_i <=  dR_j );
	   isClosestSSlepDecor(*lepB) = ( dR_j  <  dR_i );
	}
     }
  } 
  
  // FIXME: add Z/JPsi mass window 
  /*   
  const float Zmass(91187.6);     // in MeV
  const float JPsimass(3096.916); // in MeV
  
  mJPsiCand_ee, mJPsiCand_mm
  mZCand_ee, mZCand_mm, mZCand_ee_SS
  
  */
   
  return EL::StatusCode::SUCCESS;
}

//*****************************************************************************
//
// Matrix Method  & Fake Factor Method stuff
//
//

EL::StatusCode HTopMultilepAnalysis :: fakeWeightCalculator (const xAOD::EventInfo* eventInfo, const xAOD::IParticleContainer& leptons )
{
  // This method calculates the final weights to be applied to (for MM) data events
  // with leptons TT, TL, LT, LL in order to obtain the fake estimate
  // in the TT signal/control region
  // We manually plug in the lepton fake/real rates calculated via tag-and-probe analysis (FIXME) 
  //
  //NB: MM and FF Method are applied only to the 2lep and the 3lep case when there are 2 SS leptons
  
  unsigned int nLeptons = leptons.size();

  std::string region;
  // features of the two leptons that are considered for the MM
  float lepA_pt(-1.), lepA_eta(-999.), lepB_pt(-1.), lepB_eta(-999.);
  int lepA_flavour(0), lepB_flavour(0);
  
  // retrieve some previously applied event object decorations 
  static SG::AuxElement::Accessor< char > isSS01("isSS01");
  if ( !isSS01.isAvailable(*eventInfo) ) {
    Error("fakeWeightCalculator()", "isSS01 is not available. Aborting ");
    return EL::StatusCode::FAILURE;
  }   
  static SG::AuxElement::Accessor< char > isSS12("isSS12");
  if ( !isSS12.isAvailable(*eventInfo) ) {
    Error("fakeWeightCalculator()", "isSS12 is not available. Aborting ");
    return EL::StatusCode::FAILURE;
  }      
   
  // accessor to tight leptons 
  static SG::AuxElement::Accessor< char > isTightAcc("isTight");
   
  // event decorators to identify the TT,TL,LT,LL regions (looking at the two SS leptons: first is the leading, second is subleading) 
  static SG::AuxElement::Decorator< char > isTTDecor("isTT");
  static SG::AuxElement::Decorator< char > isTLDecor("isTL");
  static SG::AuxElement::Decorator< char > isLTDecor("isLT");
  static SG::AuxElement::Decorator< char > isLLDecor("isLL");
  isTTDecor( *eventInfo ) = 0;
  isTLDecor( *eventInfo ) = 0;
  isLTDecor( *eventInfo ) = 0;
  isLLDecor( *eventInfo ) = 0;
   
  // Assigning lepton kin and identifying the signal region (not taking jets into account)
  if ( nLeptons == 2 && isSS01(*eventInfo) )
  {
        // start from lepton container
	//
	// ordering criterion is simply based on pT
	// by construction, the first element of the DV is the leading lepton, the second (and last!) element is the subleading
	
	// retrieve lep0 : the leading lepton
	const xAOD::IParticle* lep0 = *(leptons.begin());
	// retrieve lep1: the subleading lepton
	const xAOD::IParticle* lep1 = *(std::next( leptons.begin(), 1));
	
	// set the region
	if      (  isTightAcc( *lep0 )    &&  isTightAcc( *lep1 )    ) { region = "TT"; isTTDecor( *eventInfo ) = 1; }
	else if (  isTightAcc( *lep0 )    &&  !(isTightAcc( *lep1 )) ) { region = "TL"; isTLDecor( *eventInfo ) = 1; }
	else if (  !(isTightAcc( *lep0 )) &&  isTightAcc( *lep1 )    ) { region = "LT"; isLTDecor( *eventInfo ) = 1; }
	else if (  !(isTightAcc( *lep0 )) &&  !(isTightAcc( *lep1 )) ) { region = "LL"; isLLDecor( *eventInfo ) = 1; }

  	if ( m_debug ) { Info("fakeWeightCalculator()", "Dilepton SS category. Region is %s ", region.c_str() ); }

        // set the properties of the two relevant leptons for future convenience
	lepA_pt  = lep0->pt();
	lepA_eta = lep0->eta();
	if ( lep0->type() == xAOD::Type::Electron ) {
	  lepA_flavour = 11;
	} else if ( lep0->type() == xAOD::Type::Muon ) {
	  lepA_flavour = 13;
	}
	
	lepB_pt  = lep1->pt();
	lepB_eta = lep1->eta();
	if ( lep1->type() == xAOD::Type::Electron ) {
	  lepB_flavour = 11;
	} else if ( lep1->type() == xAOD::Type::Muon ) {
	  lepB_flavour = 13;
	}	

  }
  else if ( nLeptons == 3 && isSS12(*eventInfo) )
  {        
        // start from lepton container
	//
	// for trilepton, ordering criterion is:
	// lep0: the OS lepton - lep1: the SS lepton with min{ deltaR(lep0) } - lep2: the other SS lepton 
        
	// retrieve some previously applied lepton object decorations 
  	static SG::AuxElement::Accessor< char > isOSlep("isOSlep");
  	static SG::AuxElement::Accessor< char > isClosestSSlep("isClosestSSlep");
	
	// need to declare these non-const pointers, and then do a const_cast.
	// this is bad practice in general, but then I do not need these pointers except for defining regions, so that's okay atm
	xAOD::IParticle* lep0(nullptr);
	xAOD::IParticle* lep1(nullptr);
	xAOD::IParticle* lep2(nullptr);
	
	for (auto lep_it :leptons ) {
	  
	  xAOD::IParticle* this_lep = const_cast< xAOD::IParticle* >(lep_it);
  	  
	  if ( !isOSlep.isAvailable(*this_lep) ) {
	    Error("fakeWeightCalculator()", "isOSlep lepton decoration is not available. Aborting ");
    	    return EL::StatusCode::FAILURE;
 	  } 	
	  
	  if ( isOSlep(*this_lep) ) {
	    
	    // retrieve lep0 : the OS lepton
	    lep0 = this_lep;
	    continue;
	    
	  } else {
	  
	    if ( !isClosestSSlep.isAvailable(*this_lep) ) {
	      Error("fakeWeightCalculator()", "isClosestSSlep lepton decoration is not available. Aborting");
    	      return EL::StatusCode::FAILURE;
 	    } 	
	    
	    // retrieve lep1 : the SS lepton with min{ deltaR(lep0) }
	    if ( isClosestSSlep(*this_lep) ) { 
	      lep1 = this_lep; 
	      continue; 
	    }
	    // retrieve lep2 : the other SS lepton 
	    lep2 = this_lep;
	  
	  }	
		
	} // close loop over lepton container
	
	// just a safety check
	if ( !( lep0 && lep1 && lep2 ) ) {
	  Error("fakeWeightCalculator()", "Trilepton region, but no lep1 and lep2 pointers! Aborting");
	  return EL::StatusCode::FAILURE;
	}
	
	// set the region
	if      (  isTightAcc( *lep1 )    &&   isTightAcc( *lep2 )   ) { region = "TT"; isTTDecor( *eventInfo ) = 1; }
	else if (  isTightAcc( *lep1 )    &&  !(isTightAcc( *lep2 )) ) { region = "TL"; isTLDecor( *eventInfo ) = 1; }
	else if (  !(isTightAcc( *lep1 )) &&  isTightAcc( *lep2 )    ) { region = "LT"; isLTDecor( *eventInfo ) = 1; }
	else if (  !(isTightAcc( *lep1 )) &&  !(isTightAcc( *lep2 )) ) { region = "LL"; isLLDecor( *eventInfo ) = 1; }
  	
	if ( m_debug ) { Info("fakeWeightCalculator()", "Trilepton 2SS+1OS category. Region (defined by the 2SS leptons) is %s ", region.c_str() ); }
	
	// set the properties of the two SS leptons for future convenience
	lepA_pt  = lep1->pt();
	lepA_eta = lep1->eta();
	if ( lep1->type() == xAOD::Type::Electron ) {
	  lepA_flavour = 11;
	} else if ( lep1->type() == xAOD::Type::Muon ) {
	  lepA_flavour = 13;
	}
	
	lepB_pt  = lep2->pt();
	lepB_eta = lep2->eta();
	if ( lep2->type() == xAOD::Type::Electron ) {
	  lepB_flavour = 11;
	} else if ( lep2->type() == xAOD::Type::Muon ) {
	  lepB_flavour = 13;
	}	
  }
  else
  {
    return EL::StatusCode::SUCCESS; //no weights in the other categories
  }


  if ( m_debug ) {
    Info("fakeWeightCalculator()", "Start calculating MM and FF weights... ");
  }

  // *******************************************
  // Now calculating MM real rate and fake rate. 
  //
  // NB: NO NEED TO CALCULATE FF RATE BECAUSE FF ARE NOW OBTAINED FROM THE MM FOR r1=r2=1
  
  // real and fake rates w/ syst variations
  std::vector<double> r1, r2, f1, f2;
  double r1up,r1dn, r2up, r2dn, f1up, f1dn, f2up, f2dn;
  
  if ( lepA_flavour == 11 )
  {
  	  r1 = calc_el_weights(lepA_pt, lepA_eta, true, false);  // first bool --> isDataDerived; second bool --> isFakeLep 
  	  f1 = calc_el_weights(lepA_pt, lepA_eta, true, true);
  }
  else if ( lepA_flavour == 13 )
  {
  	  r1 = calc_mu_weights(lepA_pt, lepA_eta, true, false);
  	  f1 = calc_mu_weights(lepA_pt, lepA_eta, true, true);
  }
  if ( lepB_flavour == 11 )
  {
  	  r2 = calc_el_weights(lepB_pt, lepB_eta, true, false);
  	  f2 = calc_el_weights(lepB_pt, lepB_eta, true, true);
  }
  else if ( lepB_flavour == 13 )
  {
  	  r2 = calc_mu_weights(lepB_pt, lepB_eta, true, false);
  	  f2 = calc_mu_weights(lepB_pt, lepB_eta, true, true);
  }
    
  if ( m_debug ) {
    Info("fakeWeightCalculator()", "Nominal real and fake rates: \n r1 = %f , r2 = %f , f1 = %f , f2 = %f ", r1.at(0), r2.at(0), f1.at(0), f2.at(0) );
  }

  // **********************************************************************************
  // declare MM and FF weight decorations
  //
  // decorate each event with a vector<double>: the first component will be the nominal,
  // the next components are w/ syst
  //
  // construct the vector with fixed number of components, and default weight 1.0 :
  //
  // nominal: MMWeightDecor( *eventInfo ).at(0);
  // rup:     MMWeightDecor( *eventInfo ).at(1);
  // rdn:     MMWeightDecor( *eventInfo ).at(2);
  // fup:     MMWeightDecor( *eventInfo ).at(3);	      
  // fdn:     MMWeightDecor( *eventInfo ).at(4);
  //
  // nominal: FFWeightDecor( *eventInfo ).at(0);
  // up:      FFWeightDecor( *eventInfo ).at(1);
  // dn:      FFWeightDecor( *eventInfo ).at(2);
  

  SG::AuxElement::Decorator< std::vector<double> > MMWeightDecor( "MMWeight" );
  if ( !MMWeightDecor.isAvailable( *eventInfo ) ) {
    MMWeightDecor( *eventInfo ) = std::vector<double>( 5, 1.0 );
  }
  SG::AuxElement::Decorator< std::vector<double> > FFWeightDecor( "FFWeight" );
  if ( !FFWeightDecor.isAvailable( *eventInfo ) ) {
    FFWeightDecor( *eventInfo ) = std::vector<double>( 3, 1.0 );
  }
  
  // *************************************
  // The Matrix Method: weight the events! 
  //
  //  Cannot be done also for the MC (even if MMweight is not used) because
  //  when it will be required to calculate the systematic sys_MMrweight_,
  //  it will then shift also MC and this is wrong. 
  //  Not a problem for FF method
  
  bool isMC = ( eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) );
  
  if ( !isMC )
  {
  	  // r cannot be 0 and has always to be more than f
  	  // WARNING! 
	  // WE SHOULD BE ALSO SURE THAT REAL EFFICIENCY IS ALWAYS < 1 FOR ANY PART OF THE PHASE SPACE. 
	  // YOU COULD HAVE SOME R(PT)*R(ETA)>1 AND THIS CANNOT HAPPEN
	  
	  double mm_weight(0.0);
	  
  	  if ( (r1.at(0) == 0) || (r2.at(0) == 0) || (r1.at(0) <= f1.at(0)) || (r2.at(0) <= f2.at(0)) ) {
	      // event will be decorated w/ null weight - will basically cancel out this event at plotting
	      if ( m_debug ) {
		Warning("fakeWeightCalculator()", "Warning! The Matrix Method cannot be applied in event %llu for run %i because : \n r1 = %f , r2 = %f , f1 = %f , f2 = %f \n ,given that pt1 = %f , eta1 = %f , pt2 = %f , eta2 = %f ", eventInfo->eventNumber(), eventInfo->runNumber(), r1.at(0), r2.at(0),  f1.at(0), f2.at(0), lepA_pt/1e3, lepA_eta, lepB_pt/1e3, lepB_eta);
		Warning("fakeWeightCalculator()", "applying MM weight = 0 ...");
	      }
	  } else {   
  	      //decorate event w/ calculated MM weight  
	      mm_weight      = calc_fake_weight( region, f1.at(0), f2.at(0), r1.at(0), r2.at(0) );
	  }
	  
	  // decorate with nominal MM weight
	  MMWeightDecor( *eventInfo ).at(0) = mm_weight;
	  
          if ( m_debug ) {
            Info("fakeWeightCalculator()", "MM final weight = %f ", mm_weight );
          }
	  
  	  if ( mm_weight != 0.0 )
  	  {
	        // decorate event w/ MM weight with systematics
	  
		 r1up = ( r1.at(1) > 1.0 ) ? 1.0 :  r1.at(1) ;
		 r2up = ( r2.at(1) > 1.0 ) ? 1.0 :  r2.at(1) ;
  	     	 r1dn = r1.at(2);
  	     	 r2dn = r2.at(2);

  	     	 f1up = f1.at(1);
  	     	 f2up = f2.at(1);
  	     	 f1dn = ( f1.at(2) < 0.0 ) ? 0.0 :  f1.at(2) ;
  	     	 f2dn = ( f2.at(2) < 0.0 ) ? 0.0 :  f2.at(2) ;

		 // rup syst
  	     	 MMWeightDecor( *eventInfo ).at(1) = ( calc_fake_weight( region, f1.at(0), f2.at(0), r1up, r2up ) / mm_weight );
		 // fdn syst
		 MMWeightDecor( *eventInfo ).at(4) = ( calc_fake_weight( region, f1dn, f2dn, r1.at(0), r2.at(0) ) / mm_weight );

  	     	 if ( (r1dn > f1.at(0)) && (r2dn > f2.at(0)) ) {
		     // rdn syst
		     MMWeightDecor( *eventInfo ).at(2) = ( calc_fake_weight(region, f1.at(0), f2.at(0), r1dn, r2dn) / mm_weight );
		 } else {
		     if ( m_debug ) {
		        Warning("fakeWeightCalculator()", "Warning! Systematic MMWeight_rdn cannot be calculated in event %llu for run %i because : \n r1dn = %f , r2dn = %f , f1 = %f , f2 = %f \n ,given that pt1 = %f , eta1 = %f , pt2 = %f , eta2 = %f ", eventInfo->eventNumber(), eventInfo->runNumber(), r1dn, r2dn,  f1.at(0), f2.at(0), lepA_pt/1e3, lepA_eta, lepB_pt/1e3, lepB_eta);
		     }
		 }
		 
  	     	 if ( (r1.at(0) > f1up) && (r2.at(0) > f2up) ) {
		     // fup syst
		     MMWeightDecor( *eventInfo ).at(3) = ( calc_fake_weight(region, f1up, f2up, r1.at(0),  r2.at(0)) / mm_weight );
		 } else {
		     if ( m_debug ) {
			Warning("fakeWeightCalculator()", "Warning! Systematic MMWeight_fup cannot be calculated in event %llu for run %i because : \n r1 = %f , r2 = %f , f1up = %f , f2up = %f \n ,given that pt1 = %f , eta1 = %f , pt2 = %f , eta2 = %f ", eventInfo->eventNumber(), eventInfo->runNumber(), r1.at(0), r2.at(0),  f1up, f2up, lepA_pt/1e3, lepA_eta, lepB_pt/1e3, lepB_eta);
  	     	     }
		 }
	  }
  	 
  } // close check on isMC
    
  // *****************************************
  // The Fake Factor Method: weight the events! 
  //
  // NB: it must hold: r = 1, f != 1 
  //  
  //  
  
  double ff_weight(0.0);
  
  if ( (f1.at(0) == 1.0) || (f2.at(0) == 1.0) )
  {
  	  if (m_debug) {
	    Warning("fakeWeightCalculator()", "Warning! The Fake Factor Method cannot be applied in event %llu for run %i because : \n f1 = %f , f2 = %f \n ,given that pt1 = %f , eta1 = %f , pt2 = %f , eta2 = %f ", eventInfo->eventNumber(), eventInfo->runNumber(), f1.at(0), f2.at(0), lepA_pt/1e3, lepA_eta, lepB_pt/1e3, lepB_eta);
	    Warning("fakeWeightCalculator()", "applying FF weight = 0 ...");
	  }
	  
	  // decorate event w/ (nominal) null weight
	  FFWeightDecor( *eventInfo ).at(0) = ff_weight;	      
  }
  else
  {	  
  	  //decorate event w/ calculated nominal FF weight  
	  ff_weight            =  calc_fake_weight( region, f1.at(0), f2.at(0) );
  	  FFWeightDecor( *eventInfo ).at(0) = ff_weight;
  	  
  	  if ( ff_weight != 0.0 )
  	  {
	        // decorate event w/ FF weight with systematics
  		
		f1up = f1.at(1);
  		f2up = f2.at(1);
		f1dn = ( f1.at(2) < 0.0 ) ? 0.0 : f1.at(2) ;
  	  	f2dn = ( f2.at(2) < 0.0 ) ? 0.0 : f2.at(2);
		
		// dn syst
  		FFWeightDecor( *eventInfo ).at(2) = ( calc_fake_weight( region, f1dn, f2dn ) / ff_weight );
  		  
  		if ( (f1up < 1) && (f2up < 1)) {
  		    // up syst
		    FFWeightDecor( *eventInfo ).at(1) = ( calc_fake_weight( region, f1up, f2up ) / ff_weight );
  		} else {
		    Warning("fakeWeightCalculator()", "Warning! Systematic FFWeight_up cannot be calculated in event %llu for run %i because : \n f1up = %f , f2up = %f \n ,given that pt1 = %f , eta1 = %f , pt2 = %f , eta2 = %f ", eventInfo->eventNumber(), eventInfo->runNumber(), f1up, f2up, lepA_pt/1e3, lepA_eta, lepB_pt/1e3, lepB_eta);
		}
	  }
  }

  if ( m_debug ) {
    Info("fakeWeightCalculator()", "FF final weight = %f ", ff_weight );
  }

  return EL::StatusCode::SUCCESS;

}

double  HTopMultilepAnalysis :: scaleFactorToRate( double val )
{
  if ( val < 0 ) val = 0.0;
  return (val /(val+1) );
}

std::vector<double>  HTopMultilepAnalysis :: calc_el_weights( float pt, float eta, bool isDataDerived, bool isFakeLep )
{
  
  // as a first thing, convert pT in GeV!
  pt = pt/1e3;
  
  int n_bins_pt_r(14);
  int n_bins_pt_f(6);
  int n_bins_eta(8);
  float bins_pt_r[] = { 10,15,20,25,30,35,40,45,50,55,60,70,80,90,2000 };
  float bins_pt_f[] = { 10,15,20,25,35,50,2000 };
  float bins_eta[]  = {  0.0 , 0.5 , 0.8 , 1.1 , 1.37 , 1.52 , 2.0 , 2.25 , 5 };

  // FIXME: for now, manually plug in Francesco's Run1 rates

  // THESE REAL/FAKE RATES ARE TAKEN FROM THE EMU, EE CHANNELS COMBINED TOGETHER
  const float el_ff_tot       = 0.077;
  const float el_ff_pt[]      = { 0.128, 0.096, 0.058, 0.0005, 0.0005, 0.0005 }; //the 0.0005 value are set by hand because the 3 digit rounding of the factor gives 0 which means that the vaule is smaller than 0.0005. So the upper value 0.0005 is used
  const float el_ff_pt_err[]  = { 0.012, 0.018, 0.026, 0.029, 0.057, 0.116 };
  const float el_ff_eta[]     = { 0.072, 0.066, 0.081, 0.065, 0.0, 0.118, 0.174, 0.103 };
  const float el_ff_eta_err[] = { 0.015, 0.02, 0.019, 0.022, 0.044, 0.031, 0.054, 0.058 };

  const float el_rf_tot       = 4.146;
  const float el_rf_pt[]      = { 1.095, 1.808, 2.634, 3.152, 4.086, 4.511, 4.991, 5.744, 6.822, 8.374, 8.781, 9.455, 13.333, 16.324 };
  const float el_rf_pt_err[]  = { 0.027, 0.043, 0.064, 0.074, 0.103, 0.121, 0.146, 0.186, 0.257, 0.371, 0.328, 0.448, 0.916, 0.797 };
  const float el_rf_eta[]     = { 4.232, 4.656, 4.02, 3.326, 2.218, 3.993, 5.329, 4.343 };
  const float el_rf_eta_err[] = { 0.058, 0.089, 0.079, 0.072, 0.232, 0.086, 0.201, 0.201 };
  
  // real rates derived from MC only
  const float el_rfmc_tot       = 5.154;
  const float el_rfmc_pt[]      = { 1.897, 2.462, 3.137, 3.926, 4.765, 5.453, 6.047, 6.67, 7.965, 8.498, 10.205, 11.451, 14.468, 18.656 };
  const float el_rfmc_pt_err[]  = { 0.041, 0.043, 0.054, 0.066, 0.084, 0.122, 0.137, 0.164, 0.201, 0.245, 0.247, 0.34, 0.592, 0.427 };
  const float el_rfmc_eta[]     = { 5.617, 6.056, 4.953, 4.027, 2.999, 4.334, 6.239, 5.406 };
  const float el_rfmc_eta_err[] = { 0.058, 0.09, 0.074, 0.06, 0.209, 0.065, 0.2, 0.218 };
  
  /*
  //THESE FAKE FACTORS ARE TAKEN FROM THE EMU ONLY!!!
    const double el_ff_tot = 0.083;
  //THESE FAKE FACTORS ARE TAKEN FROM THE EE ONLY!!!
    const double el_ff_tot = 0.073;
  */

  std::vector<double> weights(3,0.0); //initialized with zeroes
  
  weights.at(0) = 1.0;
  double error;
  
  for(int e = 0; e < n_bins_eta; e++)
  {
     if ( ( fabs(eta) >= bins_eta[e] ) && ( fabs(eta) < bins_eta[e+1] ) )
     {
     	if (isFakeLep)
     	{
     	   for(int p = 0; p < n_bins_pt_f; p++)
	   {
     	     if ( ( pt >= bins_pt_f[p] ) && ( pt<bins_pt_f[p+1] ) )
     	     {
     	       // these now calculating are still ScaleFactors
     	       weights.at(0) = ( el_ff_pt[p] * el_ff_eta[e] / el_ff_tot );
     	       error         = sqrt( (el_ff_pt[p]*el_ff_eta_err[e])*(el_ff_pt[p]*el_ff_eta_err[e]) + (el_ff_pt_err[p]*el_ff_eta[e])*(el_ff_pt_err[p]*el_ff_eta[e]) );
     	       // up syst
     	       weights.at(1) = ( el_ff_pt[p] * el_ff_eta[e] + error ) / el_ff_tot;
     	       // down syst
     	       if (el_ff_pt[p]*el_ff_eta[e] - error > 0) {
     	   	  weights.at(2) = ( el_ff_pt[p] * el_ff_eta[e] - error ) / el_ff_tot;
     	       } else {
     	   	  weights.at(2) = 0.0;
	       }	     
     	     }
	   } // close loop on pT bins: fake lepton   
     	}
     	else
     	{
     	   for(int p = 0; p < n_bins_pt_r; p++)
	   {
     	     if ( ( pt >= bins_pt_r[p] ) && ( pt < bins_pt_r[p+1] ) )
     	     {
     	   	 if (isDataDerived)
     	   	 {
     	   	    //these now calculating are still ScaleFactors
     	   	    weights.at(0) = ( el_rf_pt[p] * el_rf_eta[e] / el_rf_tot );
     	   	    error	  = sqrt( (el_rf_pt[p]*el_rf_eta_err[e])*(el_rf_pt[p]*el_rf_eta_err[e]) + (el_rf_pt_err[p]*el_rf_eta[e])*(el_rf_pt_err[p]*el_rf_eta[e]) );
     	   	    // up syst
     	   	    weights.at(1) = ( el_rf_pt[p] * el_rf_eta[e] + error ) / el_rf_tot;
     	   	    //down
     	   	    if (el_rf_pt[p]*el_rf_eta[e] - error > 0) {
     	   	       weights.at(2) = ( el_rf_pt[p] * el_rf_eta[e] - error ) / el_rf_tot;
     	   	    } else {
     	   	       weights.at(2) = 0.0;						    
     	   	    }
		 }
     	   	 else
     	   	 {
     	   	    //these now calculating are still ScaleFactors
     	   	    weights.at(0) = ( el_rfmc_pt[p] * el_rfmc_eta[e] / el_rfmc_tot );
     	   	    error         = sqrt( (el_rfmc_pt[p]*el_rfmc_eta_err[e])*(el_rfmc_pt[p]*el_rfmc_eta_err[e]) + (el_rfmc_pt_err[p]*el_rfmc_eta[e])*(el_rfmc_pt_err[p]*el_rfmc_eta[e]) );
     	   	    // up syst
     	   	    weights.at(1) = ( el_rfmc_pt[p] * el_rfmc_eta[e] + error ) / el_rfmc_tot;
     	   	    // down syst
     	   	    if (el_rfmc_pt[p]*el_rfmc_eta[e] - error > 0) {
     	   	       weights.at(2) = ( el_rfmc_pt[p] * el_rfmc_eta[e] - error ) / el_rfmc_tot;
     	   	    } else {
     	   	       weights.at(2) = 0.0;
     	   	    }
		 }
     	     } 
	   } // close loop on pT bins: real lepton				   
     	}
     }
  } // close loop on eta bins
  	  
  //Now converting to the rates for the MM
  if ( m_debug ) {
    Info("calc_el_weights()", " Electron SF = %f ( up = %f , dn = %f )", weights.at(0), weights.at(1), weights.at(2) );
  }
  weights.at(0)=scaleFactorToRate(weights.at(0));
  weights.at(1)=scaleFactorToRate(weights.at(1));
  weights.at(2)=scaleFactorToRate(weights.at(2));
  
  if ( m_debug ) {
    Info("calc_el_weights()", " Electron MM weight = %f ( up = %f , dn = %f )", weights.at(0), weights.at(1), weights.at(2) );
  }
  
  return weights;
}  

std::vector<double>  HTopMultilepAnalysis :: calc_mu_weights( float pt, float eta, bool isDataDerived, bool isFakeLep )
{  
  
  // as a first thing, convert pT in GeV!
  pt = pt/1e3;

  int n_bins_pt_r(14);
  int n_bins_pt_f(6);
  int n_bins_eta(8);
  float bins_pt_r[] = { 10,15,20,25,30,35,40,45,50,55,60,70,80,90,2000 };
  float bins_pt_f[] = { 10,15,20,25,35,50,2000 };
  float bins_eta[]  = {  0.0 , 0.5 , 0.8 , 1.1 , 1.37 , 1.52 , 2.0 , 2.25 , 5 };

  // FIXME: for now, manually plug in Francesco's Run1 rates

  // THESE REAL/FAKE RATES ARE TAKEN FROM THE EMU, MUMU CHANNELS COMBINED TOGETHER
  const float mu_ff_tot       = 0.095;
  const float mu_ff_pt[]      = { 0.108, 0.067, 0.028, 0.06, 0.006, 0.0005 }; //the 0.0005 value are set by hand because the 3 digit rounding of the factor gives 0 which means that the vaule is smaller than 0.0005. So the upper value 0.0005 is used
  const float mu_ff_pt_err[]  = { 0.008, 0.01, 0.01, 0.018, 0.03, 0.109 }; 
  const float mu_ff_eta[]     = { 0.069, 0.063, 0.056, 0.073, 0.067, 0.067, 0.155, 0.167 };
  const float mu_ff_eta_err[] = { 0.01, 0.012, 0.013, 0.017, 0.02, 0.011, 0.028, 0.031 };
  
  const float mu_rf_tot       = 7.837;
  const float mu_rf_pt[]      = { 1.448, 2.697, 4.356, 6.815, 8.906, 12.391, 15.068, 18.331, 23.531, 27.884, 30.334, 32.93, 34.282, 43.026 };
  const float mu_rf_pt_err[]  = { 0.031, 0.061, 0.103, 0.181, 0.264, 0.429, 0.591, 0.848, 1.331, 1.859, 1.732, 2.458, 3.316, 3.117 };
  const float mu_rf_eta[]     = { 7.19, 7.424, 8.202, 8.334, 8.121, 8.298, 8.438, 9.143 };
  const float mu_rf_eta_err[] = { 0.116, 0.155, 0.188, 0.225, 0.297, 0.192, 0.322, 0.41 };
  
  //real from MC only
  const float mu_rfmc_tot       = 10.34;
  const float mu_rfmc_pt[]      = { 2.936, 4.259, 5.719, 8.204, 10.501, 13.668, 17.139, 21.019, 23.871, 27.838, 31.526, 33.667, 37.315, 46.951 };
  const float mu_rfmc_pt_err[]  = { 0.059, 0.074, 0.102, 0.164, 0.26, 0.375, 0.519, 0.626, 0.735, 0.887, 0.907, 1.971, 1.547, 1.32 };
  const float mu_rfmc_eta[]     = { 9.201, 9.217, 11.054, 10.85, 11.636, 11.473, 11.565, 14.412 };
  const float mu_rfmc_eta_err[] = { 0.112, 0.139, 0.184, 0.231, 0.352, 0.23, 0.484, 0.558 };
  			   
  /*
  //THESE FAKE FACTORS ARE TAKEN FROM THE EMU ONLY!!!
  const double mu_ff_tot = 0.08;
  //THESE FAKE FACTORS ARE TAKEN FROM THE MUMU ONLY!!!
  const double mu_ff_tot = 0.107;
  */
 
  std::vector<double> weights(3,0.0); //initialized with zeroes
  
  weights.at(0) = 1.0;
  double error;
  
  for(int e = 0; e < n_bins_eta; e++)
  {
     if ( ( fabs(eta) >= bins_eta[e] ) && ( fabs(eta) < bins_eta[e+1] ) )
     {
     	if (isFakeLep)
     	{
     	   for(int p = 0; p < n_bins_pt_f; p++)
	   {
     	     if ( ( pt >= bins_pt_f[p] ) && ( pt<bins_pt_f[p+1] ) )
     	     {
     	       // these now calculating are still ScaleFactors
     	       weights.at(0) = ( mu_ff_pt[p] * mu_ff_eta[e] / mu_ff_tot );
     	       error         = sqrt( (mu_ff_pt[p]*mu_ff_eta_err[e])*(mu_ff_pt[p]*mu_ff_eta_err[e]) + (mu_ff_pt_err[p]*mu_ff_eta[e])*(mu_ff_pt_err[p]*mu_ff_eta[e]) );
     	       // up syst
     	       weights.at(1) = ( mu_ff_pt[p] * mu_ff_eta[e] + error ) / mu_ff_tot;
     	       // down syst
     	       if (mu_ff_pt[p]*mu_ff_eta[e] - error > 0) {
     	   	  weights.at(2) = ( mu_ff_pt[p] * mu_ff_eta[e] - error ) / mu_ff_tot;
     	       } else {
     	   	  weights.at(2) = 0.0;
	       }	     
     	     }
	   } // close loop on pT bins: fake lepton   
     	}
     	else
     	{
     	   for(int p = 0; p < n_bins_pt_r; p++)
	   {
     	     if ( ( pt >= bins_pt_r[p] ) && ( pt < bins_pt_r[p+1] ) )
     	     {
     	   	 if (isDataDerived)
     	   	 {
     	   	    //these now calculating are still ScaleFactors
     	   	    weights.at(0) = ( mu_rf_pt[p] * mu_rf_eta[e] / mu_rf_tot );
     	   	    error	  = sqrt( (mu_rf_pt[p]*mu_rf_eta_err[e])*(mu_rf_pt[p]*mu_rf_eta_err[e]) + (mu_rf_pt_err[p]*mu_rf_eta[e])*(mu_rf_pt_err[p]*mu_rf_eta[e]) );
     	   	    // up syst
     	   	    weights.at(1) = ( mu_rf_pt[p] * mu_rf_eta[e] + error ) / mu_rf_tot;
     	   	    //down
     	   	    if (mu_rf_pt[p]*mu_rf_eta[e] - error > 0) {
     	   	       weights.at(2) = ( mu_rf_pt[p] * mu_rf_eta[e] - error ) / mu_rf_tot;
     	   	    } else {
     	   	       weights.at(2) = 0.0;						    
     	   	    }
		 }
     	   	 else
     	   	 {
     	   	    //these now calculating are still ScaleFactors
     	   	    weights.at(0) = ( mu_rfmc_pt[p] * mu_rfmc_eta[e] / mu_rfmc_tot );
     	   	    error         = sqrt( (mu_rfmc_pt[p]*mu_rfmc_eta_err[e])*(mu_rfmc_pt[p]*mu_rfmc_eta_err[e]) + (mu_rfmc_pt_err[p]*mu_rfmc_eta[e])*(mu_rfmc_pt_err[p]*mu_rfmc_eta[e]) );
     	   	    // up syst
     	   	    weights.at(1) = ( mu_rfmc_pt[p] * mu_rfmc_eta[e] + error ) / mu_rfmc_tot;
     	   	    // down syst
     	   	    if (mu_rfmc_pt[p]*mu_rfmc_eta[e] - error > 0) {
     	   	       weights.at(2) = ( mu_rfmc_pt[p] * mu_rfmc_eta[e] - error ) / mu_rfmc_tot;
     	   	    } else {
     	   	       weights.at(2) = 0.0;
     	   	    }
		 }
     	     } 
	   } // close loop on pT bins: real lepton				   
     	}
     }
  } // close loop on eta bins
  	  
  //Now converting to the rates for the MM
  if ( m_debug ) {
    Info("calc_mu_weights()", " Muon SF = %f ( up = %f , dn = %f )", weights.at(0), weights.at(1), weights.at(2) );
  }
  weights.at(0)=scaleFactorToRate(weights.at(0));
  weights.at(1)=scaleFactorToRate(weights.at(1));
  weights.at(2)=scaleFactorToRate(weights.at(2));
  
  if ( m_debug ) {
    Info("calc_mu_weights()", " Muon MM weight = %f ( up = %f , dn = %f )", weights.at(0), weights.at(1), weights.at(2) );
  }
  
  return weights;

}  


double HTopMultilepAnalysis :: calc_fake_weight( std::string region, double f1, double f2, double r1, double r2 )
{
   //The Fake Factor Method weight is obtained under the hypothesis
   // r1=1 and r2=1 
   // So for the FF Method, you will need to pass just region, f1 and f2 (NB: make sure they are different from 0 and from 1!)
   // For the MM, you will need to pass also r1 and r2
   // region is defined as one between "TT", "TL", "LT", "LL"

   double weight = 1.0; 
   double alpha  = 1.0 / ( (r1-f1) * (r2-f2) );
   
   if      (region=="TT") { weight = 1 - ( r1 * r2 * (1-f1) * (1-f2) * alpha ); }
   else if (region=="TL") { weight = r1 * r2 * f2 * (1-f1) * alpha; }
   else if (region=="LT") { weight = r1 * r2 * f1 * (1-f2) * alpha; }
   else if (region=="LL") { weight = -1 * r1 * r2 * f1 * f2 * alpha; }
  
   if ( m_debug ) {
     Info("calc_fake_weight()", "In region %s : \n weight = %.15f , alpha = %.15f ", region.c_str(), weight, alpha);
   }

   return weight;
}

