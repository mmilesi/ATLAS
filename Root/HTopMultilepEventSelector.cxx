/**********************************************
 *
 * Module that performs an event selection for 
 * HTopMultilepAnalysis.
 *
 * M. Milesi (marco.milesi@cern.ch)
 * 
 *
 *********************************************/

// c++ include(s):
#include <iostream>
#include <typeinfo>
#include <sstream>

// EL include(s):
#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>

// EDM include(s):
#include "xAODEventInfo/EventInfo.h"
#include "xAODEgamma/Electron.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODMuon/Muon.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODTau/TauJet.h"
#include "xAODTau/TauJetContainer.h"
#include "xAODBase/IParticleContainer.h"
#include "xAODBase/IParticle.h"
#include "AthContainers/ConstDataVector.h"
#include "AthContainers/DataVector.h"

// package include(s):
#include "HTopMultilepAnalysis/HTopMultilepEventSelector.h"
#include "xAODAnaHelpers/HelperClasses.h"
#include "xAODAnaHelpers/HelperFunctions.h"
#include <xAODAnaHelpers/tools/ReturnCheck.h>
#include <xAODAnaHelpers/tools/ReturnCheckConfig.h>

// external tools include(s):
#include "TauAnalysisTools/TauSelectionTool.h"
#include "TauAnalysisTools/Enums.h"

// ROOT include(s):
#include "TEnv.h"
#include "TFile.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "TObjString.h"

// this is needed to distribute the algorithm to the workers
ClassImp(HTopMultilepEventSelector)


HTopMultilepEventSelector :: HTopMultilepEventSelector () :
  m_cutflowHist(nullptr),
  m_cutflowHistW(nullptr),
  m_TauSelTool(nullptr)
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().

  Info("HTopMultilepEventSelector()", "Calling constructor");

  m_useCutFlow           = true; 
  m_DC14                 = false;
  
  m_inContainerName_el   = "";     
  m_inContainerName_mu   = "";     
  m_inContainerName_jets = "";   
  m_inContainerName_tau  = "";    
  
  m_doMinObjCut = false;
  m_doMaxObjCut = false;
  m_n_leptons_min = 0;
  m_n_leptons_max = 100000;
  m_n_leptons_with_tau_min = 0;
  m_n_jets_min = 0;
  m_n_jets_max = 100000; 
  m_n_bjets_min = 0;
  m_n_taus_min = 0;

  m_BTag_WP = "FixedCutBEff_77";

  m_leptons_eta_max = 2.6;	    
  m_leading_lep_pT_min = 0.0;    
  m_subleading_lep_pT_min = 0.0; 
  
  m_passAuxDecorKeys = ""; 
  m_failAuxDecorKeys = ""; 
  
}

HTopMultilepEventSelector::~HTopMultilepEventSelector() {}

EL::StatusCode  HTopMultilepEventSelector :: configure ()
{

  if ( !getConfig().empty() ) {

    // read in user configuration from text file
    //
    TEnv *config = new TEnv(getConfig(true).c_str());
    if ( !config ) {
      Error("HTopMultilepEventSelection()", "Failed to initialize reading of config file. Exiting." );
      return EL::StatusCode::FAILURE;
    }
    Info("configure()", "Configuing HTopMultilepEventSelector Interface. User configuration read from : %s", getConfig().c_str());
 
    // read debug flag from .config file
    //
    m_debug	                 = config->GetValue("Debug" ,      m_debug);
    m_useCutFlow                 = config->GetValue("UseCutFlow",  m_useCutFlow);
    m_DC14                       = config->GetValue("DC14", m_DC14 );
    
    // input container to be read from TEvent or TStore
    //
    m_inContainerName_el	 = config->GetValue("InputContainerElectrons", m_inContainerName_el.c_str());
    m_inContainerName_mu	 = config->GetValue("InputContainerMuons",     m_inContainerName_mu.c_str());
    m_inContainerName_jets	 = config->GetValue("InputContainerJets",      m_inContainerName_jets.c_str());
    m_inContainerName_tau	 = config->GetValue("InputContainerTaus",      m_inContainerName_tau.c_str());
    
    // configurable cuts
    //
    m_doMinObjCut		 = config->GetValue("DoMinObjCut", m_doMinObjCut);
    m_doMaxObjCut		 = config->GetValue("DoMaxObjCut", m_doMaxObjCut);  
    m_n_leptons_min		 = config->GetValue("MinNLeptons", m_n_leptons_min);
    m_n_leptons_max		 = config->GetValue("MaxNLeptons", m_n_leptons_max);
    m_n_leptons_with_tau_min	 = config->GetValue("MinNLeptonsWithTau", m_n_leptons_with_tau_min);
    m_n_taus_min		 = config->GetValue("MinNTaus", m_n_taus_min);
    m_n_jets_min		 = config->GetValue("MinNJets", m_n_jets_min); 
    m_n_jets_max		 = config->GetValue("MaxNJets", m_n_jets_max); 
    m_n_bjets_min		 = config->GetValue("MinNBjets", m_n_bjets_min); 
    m_leading_lep_pT_min	 = config->GetValue("pTMinLeadingLepton",  m_leading_lep_pT_min);
    m_subleading_lep_pT_min	 = config->GetValue("pTMinSubLeadingLepton",  m_subleading_lep_pT_min);

    // BTag WP to count nbjets
    //
    m_BTag_WP                     = config->GetValue("BTagWP", m_BTag_WP.c_str());
     
    if ( m_inContainerName_el.empty() ) {
      Error("configure()", "InputContainerElectrons is empty!");
      return EL::StatusCode::FAILURE;
    }
    if ( m_inContainerName_mu.empty() ) {
      Error("configure()", "InputContainerMuons is empty!");
    return EL::StatusCode::FAILURE;
    }
    
    config->Print();
    
    Info("configure()", "HTopMultilepEventSelector Interface succesfully configured!");
    
    delete config; config = nullptr;
  }
  
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode HTopMultilepEventSelector :: setupJob (EL::Job& job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.

  Info("setupJob()", "Calling setupJob");

  job.useXAOD ();
  xAOD::Init( "HTopMultilepEventSelector" ).ignore(); // call before opening first file

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepEventSelector :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.

  Info("histInitialize()", "Calling histInitialize");

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepEventSelector :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed

  Info("fileExecute()", "Calling fileExecute");


  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepEventSelector :: changeInput (bool /*firstFile*/)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.

  Info("changeInput()", "Calling changeInput");

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepEventSelector :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.

  Info("initialize()", "Initializing HTopMultilepEventSelector Interface...");

  m_event = wk()->xaodEvent();
  m_store = wk()->xaodStore();

  if ( m_useCutFlow ) {
    TFile *file = wk()->getOutputFile("cutflow");
    m_cutflowHist  = (TH1D*)file->Get("cutflow");
    m_cutflowHistW = (TH1D*)file->Get("cutflow_weighted");
    m_cutflow_bin  = m_cutflowHist->GetXaxis()->FindBin(m_name.c_str());
    m_cutflowHistW->GetXaxis()->FindBin(m_name.c_str());
  }

  Info("initialize()", "Number of events: %lld ", m_event->getEntries() );

  if ( configure() == EL::StatusCode::FAILURE ) {
    Error("initialize()", "Failed to properly configure. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  m_numEvent      = 0;
  m_numObject     = 0;
  m_numEventPass  = 0;
  m_weightNumEventPass  = 0;
  m_numObjectPass = 0;

  // initialise TauSelectionTool 
  //
  m_TauSelTool = new TauAnalysisTools::TauSelectionTool( "TauSelectionTool" );
  m_TauSelTool->setProperty("ConfigPath", "$ROOTCOREBIN/data/HTopMultilepAnalysis/Taus/recommended_selection_mc15.conf");
  m_TauSelTool->msg().setLevel( MSG::INFO );
 
  RETURN_CHECK( "HTopMultilepEventSelector::initialize()", m_TauSelTool->initialize(), "Failed to properly initialize TauSelectionTool." );

  Info("initialize()", "HTopMultilepEventSelector Interface succesfully initialized!" );

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepEventSelector :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  if ( m_debug ) { Info("execute()", "Applying HTopMultilepAnalysis Event Selection... \n"); }

  // retrieve event
  //
  const xAOD::EventInfo* eventInfo(nullptr);
  RETURN_CHECK("HTopMultilepEventSelector::execute()", HelperFunctions::retrieve(eventInfo, "EventInfo", m_event, m_store, m_debug) , "");
 
  // MC event weight 
  //
  float mcEvtWeight(1.0); 
  static SG::AuxElement::Accessor< float > mcEvtWeightAcc("mcEventWeight");
  if ( !mcEvtWeightAcc.isAvailable(*eventInfo) ) {
    Info("execute()", "event weight is not available. Aborting ");
    return EL::StatusCode::FAILURE;
  
  } 
  mcEvtWeight = mcEvtWeightAcc(*eventInfo);

  m_numEvent++;

  //bool isMC = ( eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) );
  
  // this will be the collection processed - no matter what!!
  //
  const xAOD::ElectronContainer* inElectrons(nullptr);
  RETURN_CHECK("HTopMultilepEventSelector::execute()", HelperFunctions::retrieve(inElectrons, m_inContainerName_el, m_event, m_store, m_debug) , "");

  const xAOD::MuonContainer* inMuons(nullptr);
  RETURN_CHECK("HTopMultilepEventSelector::execute()", HelperFunctions::retrieve(inMuons, m_inContainerName_mu, m_event, m_store, m_debug) , "");
  
  const xAOD::JetContainer* inJets(nullptr);   
  RETURN_CHECK("HTopMultilepAnalysis::execute()", HelperFunctions::retrieve(inJets, m_inContainerName_jets, m_event, m_store,  m_debug) , "");

  const xAOD::TauJetContainer* inTauJets(nullptr);
  RETURN_CHECK("HTopMultilepEventSelector::execute()", HelperFunctions::retrieve(inTauJets, m_inContainerName_tau, m_event, m_store, m_debug) , "");
    
  if ( m_debug ) { Info("execute()"," inElectrons N = %lu , inMuons N = %lu , inJets N = %lu, inTaus N = %lu ", inElectrons->size(), inMuons->size(), inJets->size(), inTauJets->size() ); }

  
  // create a CDV for leptons ( electrons + muons )
  //
  ConstDataVector<xAOD::IParticleContainer>* leptonsCDV(nullptr);
  leptonsCDV =  new ConstDataVector<xAOD::IParticleContainer>(SG::VIEW_ELEMENTS);
  
  // fill this CDV with electrons and muons
  //
  for ( auto el_itr : *(inElectrons) ) {
    
    const xAOD::IParticle *lepton = el_itr;
    leptonsCDV->push_back( lepton );
    if ( m_debug ) { Info("execute()","pushing electron to leptonsCDV - pT = %2f ", el_itr->pt() / 1e3 ); }
  
  }
  for ( auto mu_itr : *(inMuons) ) {
      
    const xAOD::IParticle *lepton = mu_itr;
    leptonsCDV->push_back( lepton );
    if ( m_debug ) { Info("execute()","pushing muon to leptonsCDV - pT = %2f ", mu_itr->pt() / 1e3 ); }
      
  }
  // Make a sorted version of the container 
  // (this can be on the stack! Will not be the one that is pushed to the store...)
  //
  const xAOD::IParticleContainer leptonsSorted = HelperFunctions::sort_container_pt( leptonsCDV->asDataVector() );

  unsigned int nLeptons = leptonsCDV->size();
    
  if ( m_debug ) { 
    Info("execute()","inLeptons N = %u ", nLeptons ); 
    for ( auto lep_it : leptonsSorted ) {
      Info("execute()","\t after sorting - lepton pT = %2f ", lep_it->pt()/1e3 );
    }
  }

  float leading_lep_pt(-1.0), subleading_lep_pt(-1.0); 
  if ( nLeptons > 0 ) {
    const xAOD::IParticle* leadingLepton = *(leptonsSorted.begin());
    leading_lep_pt	 = leadingLepton->pt();
    if ( nLeptons > 1 ) {
      const xAOD::IParticle* subLeadingLepton = *(std::next(leptonsSorted.begin(),1));
      subleading_lep_pt  = subLeadingLepton->pt(); 
    }
  }
  

  // create selected tau container 
  //
  ConstDataVector<xAOD::TauJetContainer>* selectedTaus(nullptr); 
  selectedTaus =  new ConstDataVector<xAOD::TauJetContainer>(SG::VIEW_ELEMENTS);
    
  for ( auto tau_itr : *inTauJets ) { 
      
    if ( m_debug ) { Info("execute()","input tau pT = %2f ", tau_itr->pt() / 1e3 ); }

    if ( !m_TauSelTool->accept(*tau_itr) ) { continue; }
    
    selectedTaus->push_back( tau_itr );
  }
    
  unsigned int nSelectedTaus = selectedTaus->size();
    
  if ( m_debug ) { 
    Info("execute()"," inTaus = %lu -  selectedTaus = %u ", inTauJets->size(), nSelectedTaus ); 
    for ( auto tau_itr : *selectedTaus ) {
      Info("execute()","\t selected tau pT = %2f ", tau_itr->pt() / 1e3 ); 
    }
  }

  // ***************************** //
  // now make the event selection
  // ***************************** //

  bool passTwoLep(false), passLepTau(false), passnLepMax(false);

  if ( nLeptons >= static_cast<unsigned int>(m_n_leptons_min)		      &&
       leading_lep_pt > static_cast<unsigned int>(m_leading_lep_pT_min)       &&
       subleading_lep_pt > static_cast<unsigned int>(m_subleading_lep_pT_min)  
       ) { 
    if ( m_debug ) {
      Info("execute()","\t leading lepton pT = %2f    ", leading_lep_pt / 1e3 );
      Info("execute()","\t subleading lepton pT = %2f ", subleading_lep_pt / 1e3 );
    }    
    passTwoLep = true;
  }

  if ( nLeptons >= static_cast<unsigned int>(m_n_leptons_with_tau_min)	   &&
       nSelectedTaus >= static_cast<unsigned int>(m_n_taus_min)		   &&
       leading_lep_pt > static_cast<unsigned int>(m_leading_lep_pT_min)     
       ) { 
    if ( m_debug ) {
      for ( auto tau_it : *selectedTaus ) {
	Info("execute()", "selected tau pT: %2f ", tau_it->pt() );
      }
      Info("execute()","\t leading lepton pT = %2f ", leading_lep_pt / 1e3 );
    }   
    passLepTau = true;
  }
   
  if ( nLeptons <= static_cast<unsigned int>(m_n_leptons_max) ) {
    passnLepMax = true;
  } 
   
  // count number of LF jets 
  //
  unsigned int nJets = inJets->size();
  bool passnJetsMin(false), passnJetsMax(false); 
  passnJetsMin  = ( nJets >= static_cast<unsigned int>(m_n_jets_min) );
  passnJetsMax  = ( nJets <= static_cast<unsigned int>(m_n_jets_max) );

  // count number of bjets
  //
  unsigned int nBjets(0);
  //
  // 2015 data,MC use MV2c20 as default
  // Look into xAODAnaHelpers/Root/BJetEfficiencyCorrector.cxx for more info
  //
  static SG::AuxElement::ConstAccessor< int > isBTag("BTag_"+m_BTag_WP);
  for( auto jet_itr : *(inJets) ) {
    if ( isBTag.isAvailable(*jet_itr) ) {
      if ( isBTag(*jet_itr) == 1 ) ++nBjets;
    }
  } 

  bool passnBJetsMin(false); 
  passnBJetsMin  = ( nBjets >= static_cast<unsigned int>(m_n_bjets_min) );
    
  if ( m_debug ) {
    Info("execute()","***********************************");     
    Info("execute()","event passes TwoLep? %i    ", static_cast<int>( passTwoLep ) );
    Info("execute()","event passes LepTau? %i    ", static_cast<int>( passLepTau ) );  
    Info("execute()","event passes nLepMax? %i   ", static_cast<int>( passnLepMax ) );  
    Info("execute()","event passes nJetsMin? %i  ", static_cast<int>( passnJetsMin ) );  
    Info("execute()","event passes nJetsMax? %i  ", static_cast<int>( passnJetsMax ) );  
    Info("execute()","event passes nBJetsMin? %i ", static_cast<int>( passnBJetsMin ) );  
    Info("execute()","*********************************** \n");  
  }
  
  bool passMinObj(false), passMaxObj(false); 
  passMinObj = ( ( passTwoLep || passLepTau ) && passnJetsMin && passnBJetsMin );
  passMaxObj = ( passnJetsMax && passnLepMax );
  
  // decide whether to skip event or not
  //
  if ( m_doMinObjCut && !passMinObj ) {
    if ( m_debug ) { Info("execute()","event did not pass minObjCut. Reject it"); }
    wk()->skipEvent();
    return EL::StatusCode::SUCCESS;
  } 
  if ( m_doMaxObjCut && !passMaxObj ) {
    if ( m_debug ) { Info("execute()","event did not pass maxObjCut. Reject it"); }
    wk()->skipEvent();
    return EL::StatusCode::SUCCESS;
  }   
  
  // add ConstDataVector(s) to TStore
  // NB: don't store a sorted container to TStore and expect it will be still sorted at retrieval!!
  //
  RETURN_CHECK( "HTopMultilepEventSelector::execute()", m_store->record( leptonsCDV, "Leptons_Selected" ), "Failed to store const data container");
  RETURN_CHECK( "HTopMultilepEventSelector::execute()", m_store->record( selectedTaus, "Taus_Selected" ), "Failed to store const data container");
  
  m_numEventPass++;
  m_weightNumEventPass += mcEvtWeight;

  // add some decorations to the event
  //
  static SG::AuxElement::Decorator< unsigned int > nLeptonsDecor("nLeptons");
  static SG::AuxElement::Decorator< unsigned int > nBjets_Decor("nBjets_"+m_BTag_WP);
  static SG::AuxElement::Decorator< unsigned int > categoryFlagDecor("categoryFlag");
  
  // compact way to categorise event based on object counting (exploiting prime numbers)
  //
  unsigned int categoryFlag(0); 
  categoryFlag = pow(2.0,static_cast<float>(nLeptons))*pow(3.0,static_cast<float>(nJets))*pow(5.0,static_cast<float>(nBjets));
  
  nLeptonsDecor( *eventInfo )             = nLeptons;
  nBjets_Decor( *eventInfo )              = nBjets;
  categoryFlagDecor( *eventInfo )         = categoryFlag;

   
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode HTopMultilepEventSelector :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.

  if ( m_debug ) { Info("postExecute()", "Calling postExecute"); }

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepEventSelector :: finalize ()
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

  if ( m_useCutFlow ) {
    Info("histFinalize()", "Filling cutflow");
    m_cutflowHist ->SetBinContent( m_cutflow_bin, m_numEventPass        );
    m_cutflowHistW->SetBinContent( m_cutflow_bin, m_weightNumEventPass  );
  }

  Info("finalize()", "Deleting tool instances...");

  if ( m_TauSelTool ) { m_TauSelTool = nullptr; delete m_TauSelTool; }

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode HTopMultilepEventSelector :: histFinalize ()
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

  Info("histFinalize()", "Calling histFinalize");

  return EL::StatusCode::SUCCESS;
}



