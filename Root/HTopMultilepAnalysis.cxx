/************************************************
 *
 * The actual HTopMultilepAnalysis algorithm.
 * Here the user categorises events, and performs
 * the background estimation.
 *
 * M. Milesi (marco.milesi@cern.ch)
 *
 ** *********************************************/

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
#include "xAODMissingET/MissingETContainer.h"

// package include(s):
#include "xAODAnaHelpers/HelperFunctions.h"
#include "xAODAnaHelpers/HelperClasses.h"
#include "xAODAnaHelpers/JetHists.h"
#include "xAODAnaHelpers/tools/ReturnCheck.h"
#include "xAODAnaHelpers/tools/ReturnCheckConfig.h"
#include "HTopMultilepAnalysis/HTopMultilepAnalysis.h"
#include "HTopMultilepAnalysis/tools/HTopReturnCheck.h"

// external tools include(s):
#include "TauAnalysisTools/TauSelectionTool.h"
#include "TauAnalysisTools/Enums.h"

// ROOT include(s):
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TEnv.h"
#include "TSystem.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

// c++ include(s)
#include <stdexcept>

// a template function to get ROOT objects from a TFile
//
template<typename T>
T* get_object( TFile& file, const std::string& name ) {
  T* obj = dynamic_cast<T*>( file.Get(name.c_str()) );
  if ( !obj ) { throw std::runtime_error("object " + name + " not found"); }
  return obj;
}

// this is needed to distribute the algorithm to the workers
ClassImp(HTopMultilepAnalysis)


HTopMultilepAnalysis :: HTopMultilepAnalysis () :
  m_cutflowHist(nullptr),
  m_cutflowHistW(nullptr),
  m_histEventCount(nullptr),
  m_totalEvents(nullptr),
  m_totalEventsW(nullptr),
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

  m_inContainerName_PreSelectedElectrons = "Electrons";
  m_inContainerName_PreSelectedMuons     = "Muons";
  m_inContainerName_PreSelectedJets      = "AntiKt4EMTopoJets";

  // to define "Tight" leptons
  m_TightElectronPID_WP            = "LHTight";
  m_TightElectronIso_WP            = "isIsolated_Gradient";
  m_TightElectronD0sig_cut         = 5.0;
  m_TightElectronTrkz0sinTheta_cut = 0.5;

  m_TightMuonD0sig_cut         = -1.0;
  m_TightMuonTrkz0sinTheta_cut = 0.5;
  m_TightMuonIso_WP            = "isIsolated_Gradient";

  // to define "Tight" taus
  m_ConfigPathTightTaus       = "$ROOTCOREBIN/data/HTopMultilepAnalysis/Taus/recommended_selection_mc15_final_sel.conf";

  m_useLooseAsLoosest         = true;
  m_useMediumAsLoosest        = false;
  m_vetoMediumNonTight        = false;
  m_useMediumAsTightest       = false;
  m_useTightAsTightest        = true;

  m_useMCForTagAndProbe       = false;

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

  TFile *fileMD = wk()->getOutputFile ("metadata");
  m_histEventCount  = dynamic_cast<TH1D*>( fileMD->Get("MetaData_EventCount") );
  if ( !m_histEventCount ) {
    Error("initialize()", "Failed to retrieve MetaData histogram. Aborting");
    return EL::StatusCode::FAILURE;
  }

  m_totalEvents  = new TH1D("TotalEvents",  "TotalEvents",  2, 1, 3);
  m_totalEventsW = new TH1D("TotalEventsW", "TotalEventsW", 2, 1, 3);
  wk() -> addOutput(m_totalEvents);
  wk() -> addOutput(m_totalEventsW);

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode HTopMultilepAnalysis :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed

  // get the MetaData tree once a new file is opened, with
  TTree *MetaData = dynamic_cast<TTree*>( wk()->inputFile()->Get("MetaData") );
  if ( !MetaData ) {
    Error("fileExecute()", "MetaData not found! Exiting.");
    return EL::StatusCode::FAILURE;
  }
  MetaData->LoadTree(0);

  //check if file is from a DxAOD
  m_isDerivation = !MetaData->GetBranch("StreamAOD");

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

  // check if sample is MC
  //
  const xAOD::EventInfo* eventInfo(nullptr);
  RETURN_CHECK("HTopMultilepAnalysis::initialize()", HelperFunctions::retrieve(eventInfo, "EventInfo", m_event, m_store, m_debug) , "");
  m_isMC = ( eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) );

  if ( m_useCutFlow ) {
    TFile *fileCF  = wk()->getOutputFile ("cutflow");
    m_cutflowHist  = dynamic_cast<TH1D*>( fileCF->Get("cutflow") );
    m_cutflowHistW = dynamic_cast<TH1D*>( fileCF->Get("cutflow_weighted") );
    // label the bins for the cutflow
    m_cutflow_bin     = m_cutflowHist->GetXaxis()->FindBin(m_name.c_str());
    // do it again for the weighted cutflow hist
    m_cutflowHistW->GetXaxis()->FindBin(m_name.c_str());
  }

  /*
  m_jetPlots = new JetHists( "highPtJets", "clean" ); // second argument: "kinematic", "clean", "energy", "resolution"
  m_jetPlots -> initialize();
  m_jetPlots -> record( wk() );
  */

  // initialise TauSelectionTool
  //
  std::string tau_sel_tool_name = std::string("TauSelectionTool_") + m_name;

  if ( asg::ToolStore::contains<TauAnalysisTools::TauSelectionTool>(tau_sel_tool_name) ) {
    m_TauSelTool = asg::ToolStore::get<TauAnalysisTools::TauSelectionTool>(tau_sel_tool_name);
  } else {

    m_TauSelTool = new TauAnalysisTools::TauSelectionTool( tau_sel_tool_name );
    m_TauSelTool->msg().setLevel( MSG::INFO ); // VERBOSE, INFO, DEBUG

    RETURN_CHECK( "HTopMultilepAnalysis::initialize()", m_TauSelTool->setProperty("ConfigPath",m_ConfigPathTightTaus.c_str()), "Failed to set ConfigPath property");
    RETURN_CHECK( "HTopMultilepAnalysis::initialize()", m_TauSelTool->initialize(), "Failed to properly initialize TauSelectionTool_HTop" );
  }

  const std::string path("$ROOTCOREBIN/data/HTopMultilepAnalysis/External/");

  // Read QMisID rates from input ROOT histograms
  //
  EL_RETURN_CHECK("HTopMultilepAnalysis::initialize()", this->readQMisIDRates( path ) );

  // Read fake(real) rates from input ROOT histograms
  //
  EL_RETURN_CHECK("HTopMultilepAnalysis::initialize()", this->readFakeRates( path ) );

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

  ++m_eventCounter;

  //---------------------------
  //***** Event information
  //---------------------------

  // retrieve event info
  const xAOD::EventInfo* eventInfo(nullptr);
  RETURN_CHECK("HTopMultilepAnalysis::execute()", HelperFunctions::retrieve(eventInfo, "EventInfo", m_event, m_store, m_verbose), "");

  if ( m_debug ) {
    Info( "execute()", " ******************************** ");
    Info( "execute()", "eventNumber =  %llu \n", eventInfo->eventNumber() );
  }

  // retrieve vertices
  const xAOD::VertexContainer* vertices(nullptr);
  RETURN_CHECK("HTopMultilepAnalysis::execute()", HelperFunctions::retrieve(vertices, "PrimaryVertices", m_event, m_store,  m_verbose), "");

  // retrieve selected objects
  const xAOD::ElectronContainer* signalElectrons(nullptr);
  RETURN_CHECK("HTopMultilepAnalysis::execute()", HelperFunctions::retrieve(signalElectrons, m_inContainerName_Electrons, m_event, m_store,  m_verbose), "");
  const xAOD::MuonContainer*     signalMuons(nullptr);
  RETURN_CHECK("HTopMultilepAnalysis::execute()", HelperFunctions::retrieve(signalMuons, m_inContainerName_Muons, m_event, m_store, m_verbose), "");
  const xAOD::JetContainer*      signalJets(nullptr);
  RETURN_CHECK("HTopMultilepAnalysis::execute()", HelperFunctions::retrieve(signalJets, m_inContainerName_Jets, m_event, m_store,  m_verbose), "");
  const xAOD::TauJetContainer*   signalTauJets(nullptr);
  RETURN_CHECK("HTopMultilepAnalysis::execute()", HelperFunctions::retrieve(signalTauJets,  m_inContainerName_Taus, m_event, m_store, m_verbose), "");
  ConstDataVector<xAOD::IParticleContainer>* leptonsCDV(nullptr);
  RETURN_CHECK("HTopMultilepAnalysis::execute()", HelperFunctions::retrieve(leptonsCDV, m_inContainerName_Leptons, m_event, m_store, m_verbose),"");
  // Make a sorted version of the lepton container
  // (this can be on the stack! Will not be pushed to the store...)
  //
  const xAOD::IParticleContainer leptonsSorted = HelperFunctions::sort_container_pt( leptonsCDV->asDataVector() );

  unsigned int nSignalLeptons   = leptonsCDV->size();
  unsigned int nSignalJets      = signalJets->size();
  unsigned int nSignalTaus      = signalTauJets->size();

  if ( m_debug ) {

    // retrieve initial objects ( only for debugging purposes )
    // NB: the name of the initial containers is hard-coded as it is not expected to change at all!
    //
    const xAOD::ElectronContainer* inElectrons(nullptr);
    RETURN_CHECK("HTopMultilepAnalysis::execute()", HelperFunctions::retrieve(inElectrons, "Electrons", m_event, m_store, m_verbose), "");
    const xAOD::MuonContainer*     inMuons(nullptr);
    RETURN_CHECK("HTopMultilepAnalysis::execute()", HelperFunctions::retrieve(inMuons, "Muons", m_event, m_store, m_verbose), "");
    const xAOD::JetContainer*      inJets(nullptr);
    RETURN_CHECK("HTopMultilepAnalysis::execute()", HelperFunctions::retrieve(inJets, "AntiKt4EMTopoJets", m_event, m_store, m_verbose), "");

    // retrieve preselected objects ( only for debugging purposes )
    //
    const xAOD::ElectronContainer* preselElectrons(nullptr);
    RETURN_CHECK("HTopMultilepAnalysis::execute()", HelperFunctions::retrieve(preselElectrons, m_inContainerName_PreSelectedElectrons, m_event, m_store, m_verbose), "");
    const xAOD::MuonContainer*     preselMuons(nullptr);
    RETURN_CHECK("HTopMultilepAnalysis::execute()", HelperFunctions::retrieve(preselMuons, m_inContainerName_PreSelectedMuons, m_event, m_store, m_verbose), "");
    const xAOD::JetContainer*      preselJets(nullptr);
    RETURN_CHECK("HTopMultilepAnalysis::execute()", HelperFunctions::retrieve(preselJets, m_inContainerName_PreSelectedJets, m_event, m_store, m_verbose), "");

    unsigned int nInElectrons     = inElectrons->size();
    unsigned int nInMuons         = inMuons->size();
    unsigned int nInJets          = inJets->size();

    unsigned int nPreselElectrons = preselElectrons->size();
    unsigned int nPreselMuons     = preselMuons->size();
    unsigned int nPreselJets      = preselJets->size();

    unsigned int nSignalElectrons = signalElectrons->size();
    unsigned int nSignalMuons     = signalMuons->size();

    Info("execute()"," Initial vs Preselected vs Selected Signal Muons: \t %u \t %u \t %u  "    , nInMuons, nPreselMuons, nSignalMuons );
    Info("execute()"," Initial vs Preselected vs Selected Signal Electrons: %u \t %u \t %u "    , nInElectrons, nPreselElectrons, nSignalElectrons );
    Info("execute()"," Initial vs Preselected vs Selected Signal Jets: \t %u \t %u \t %u "      , nInJets, nPreselJets, nSignalJets );
    Info("execute()"," Selected Signal Leptons: \t %u " , nSignalLeptons );
    Info("execute()"," Selected Signal Taus: \t %u " , nSignalTaus );

  }
  /*
  if ( nSignalLeptons > 1 ) {
    Info( "execute()", " ******************************** ");
    Info( "execute()", "eventNumber =  %llu \n", eventInfo->eventNumber() );
    Info( "execute()"," Selected Signal Leptons: \t %u " , nSignalLeptons );
  }
  */
  //-------------------------------
  //***** Retrieve event weight
  //-------------------------------

  // retrieve event weight from eventInfo (NB: will be always 1 for Data - see xAODAnaHelpers/root/BasicEventSelection.cxx)
  // for MC, it includes also PU weight. Need to multiply it by object SF weights later on!
  //
  float mcEvtWeight(1.0);
  static SG::AuxElement::Accessor< float > mcEvtWeightAcc("mcEventWeight");
  if ( !mcEvtWeightAcc.isAvailable(*eventInfo) ) {
    Error("execute()", "event weight is not available. Aborting ");
    return EL::StatusCode::FAILURE;
  }
  mcEvtWeight = mcEvtWeightAcc( *eventInfo );

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

  //-------------------------------------------
  // definition of "Tight" and "Medium" leptons
  //-------------------------------------------

  static SG::AuxElement::Decorator< char > isTightDecor("isTight");
  static SG::AuxElement::Decorator< char > isMediumDecor("isMedium");

  static SG::AuxElement::Accessor< char >  TightElectronIsoAcc(m_TightElectronIso_WP);
  static SG::AuxElement::Accessor< char >  TightElectronIDAcc(m_TightElectronPID_WP);
  static SG::AuxElement::Accessor< char >  TightMuonIsoAcc(m_TightMuonIso_WP);
  static SG::AuxElement::Accessor< float > d0SigAcc ("d0sig");
  static SG::AuxElement::Accessor< float > z0sinthetaAcc ("z0sintheta");

  // -----------------------
  //        electrons
  // -----------------------

  // first: isolation
  // second: Electron ID
  //
  std::pair<std::string,std::string> tightness_def_el = std::make_pair(m_TightElectronIso_WP,m_TightElectronPID_WP);

  for ( auto el_itr : *(signalElectrons) ) {

    // set default decoration
    //
    isTightDecor( *el_itr ) = 0;
    isMediumDecor( *el_itr ) = 0;

    if ( !d0SigAcc.isAvailable( *el_itr ) ) {
      Error("execute()", "'d0sig' attribute is not available for this electron. Aborting " );
      return EL::StatusCode::FAILURE;
    }

    if ( !z0sinthetaAcc.isAvailable( *el_itr ) ) {
      Error("execute()", "'z0sintheta' attribute is not available for this electron. Aborting " );
      return EL::StatusCode::FAILURE;
    }

    // preliminary: tighten impact parameter cuts
    //
    if ( fabs( d0SigAcc( *el_itr ) ) < m_TightElectronD0sig_cut && fabs( z0sinthetaAcc( *el_itr ) ) < m_TightElectronTrkz0sinTheta_cut ) {

      // if using isolation...
      //
      if ( !tightness_def_el.first.empty() ) {

	if ( !TightElectronIsoAcc.isAvailable( *el_itr ) ) {
	  Error("execute()", "'%s' attribute is not available for this electron. Aborting ", m_TightElectronIso_WP.c_str() );
	  return EL::StatusCode::FAILURE;
	}

	// if using *also* Electron ID...
	//
	if ( !tightness_def_el.second.empty() ) {

	  if ( !TightElectronIDAcc.isAvailable( *el_itr ) ) {
	    Error("execute()", "'%s' attribute is not available for this electron. Aborting ", m_TightElectronPID_WP.c_str() );
	    return EL::StatusCode::FAILURE;
	  }

	  // "Tight"  ---> IP + ID + iso
	  // "Medium" ---> IP + ID (any iso)
	  //
	  if ( ( TightElectronIsoAcc( *el_itr ) == 1 ) && ( TightElectronIDAcc( *el_itr ) == 1 ) ) { isTightDecor( *el_itr ) = 1; }
	  if ( ( TightElectronIDAcc( *el_itr ) == 1 )  )                                           { isMediumDecor( *el_itr ) = 1; }

	} else {

	  // "Tight"  ---> IP + iso (any ID)
	  // "Medium" ---> IP (any iso, ID)
	  //
	  if ( TightElectronIsoAcc( *el_itr ) == 1 ) { isTightDecor( *el_itr ) = 1; }
	  isMediumDecor( *el_itr ) = 1;

	}

      }
      // if not using isolation, but using Electron ID...
      //
      else if ( !tightness_def_el.second.empty() ) {

	if ( !TightElectronIDAcc.isAvailable( *el_itr ) ) {
	  Error("execute()", "'%s' attribute is not available for this electron. Aborting ", m_TightElectronPID_WP.c_str() );
	  return EL::StatusCode::FAILURE;
	}

	// "Tight"  ---> IP + ID (any iso)
	// "Medium" ---> IP (any ID, iso)
	//
	if ( TightElectronIDAcc( *el_itr ) == 1 ) { isTightDecor( *el_itr ) = 1; }
	isMediumDecor( *el_itr ) = 1;

      }
      // if using neither isolation, nor Electron ID..
      //
      else {
	Error("execute()", "Need at least isolation or ElectronID requirement to define 'Tight' electrons. Aborting" );
	return EL::StatusCode::FAILURE;
      }
    }
  }

  // -----------------------
  //          muons
  // -----------------------

  // first: isolation
  // second: d0sig
  //
  std::pair<std::string,float> tightness_def_mu = std::make_pair(m_TightMuonIso_WP,m_TightMuonD0sig_cut);

  for ( auto mu_itr : *(signalMuons) ) {

    // set default decoration
    //
    isTightDecor( *mu_itr ) =  0;
    isMediumDecor( *mu_itr ) =  0;

    if ( !z0sinthetaAcc.isAvailable( *mu_itr ) ) {
      Error("execute()", "'z0sintheta' attribute is not available for this muon. Aborting " );
      return EL::StatusCode::FAILURE;
    }

    // preliminary: tighten impact parameter cuts
    //
    if ( fabs( z0sinthetaAcc( *mu_itr ) ) < m_TightMuonTrkz0sinTheta_cut ) {

      // if using isolation...
      //
      if ( !tightness_def_mu.first.empty() ) {

	if ( !TightMuonIsoAcc.isAvailable( *mu_itr ) ) {
	  Error("execute()", "'%s' attribute is not available for this muon. Aborting ", m_TightMuonIso_WP.c_str() );
	  return EL::StatusCode::FAILURE;
	}

	// if using *also* d0sig...
	//
	if ( tightness_def_mu.second > 0.0 ) {

	  if ( !d0SigAcc.isAvailable( *mu_itr ) ) {
	    Error("execute()", "'d0sig' attribute is not available for this muon. Aborting " );
	    return EL::StatusCode::FAILURE;
	  }

	  // "Tight"  ---> z0 + iso + d0sig
	  // "Medium" ---> z0 + d0sig (any iso)
	  //
	  if ( ( TightMuonIsoAcc( *mu_itr ) == 1 ) && ( fabs( d0SigAcc( *mu_itr ) ) < tightness_def_mu.second ) ) { isTightDecor( *mu_itr ) = 1; }
	  if ( fabs( d0SigAcc( *mu_itr ) ) < tightness_def_mu.second )                                            { isMediumDecor( *mu_itr ) = 1; }

	} else {

	  // "Tight"  ---> z0 + iso (any d0sig)
	  // "Medium" ---> z0 (any iso, d0sig)
	  //
	  if ( TightMuonIsoAcc( *mu_itr ) == 1 ) { isTightDecor( *mu_itr ) = 1; }
	  isMediumDecor( *mu_itr ) = 1;

	}

      }
      // if not using isolation, but using d0sig...
      //
      else if ( tightness_def_mu.second > 0.0 ) {

	if ( !d0SigAcc.isAvailable( *mu_itr ) ) {
	  Error("execute()", "'d0sig' attribute is not available for this muon. Aborting " );
	  return EL::StatusCode::FAILURE;
	}

	// "Tight"  ---> z0 + d0sig (any iso)
	// "Medium" ---> z0 (any d0sig, iso)
	//
	if ( fabs( d0SigAcc( *mu_itr ) ) < tightness_def_mu.second ) { isTightDecor( *mu_itr ) = 1; }
	isMediumDecor( *mu_itr ) = 1;

      }
      // if using neither isolation, nor d0sig..
      //
      else {
	Error("execute()", "Need at least isolation or d0sig requirement to define 'Tight' muons. Aborting" );
	return EL::StatusCode::FAILURE;
      }
    }
  }

  //-------------------------------------------
  // calculate lepton SFs for the event
  //-------------------------------------------

  if ( m_isMC ) {
    EL_RETURN_CHECK("HTopMultilepAnalysis::execute()", this->computeEventLepTrigSF( eventInfo, leptonsSorted ) );
    EL_RETURN_CHECK("HTopMultilepAnalysis::execute()", this->computeEventLepSF( eventInfo, leptonsSorted, SFType::RECO ) );
    EL_RETURN_CHECK("HTopMultilepAnalysis::execute()", this->computeEventLepSF( eventInfo, leptonsSorted, SFType::ISOLATION ) );
    EL_RETURN_CHECK("HTopMultilepAnalysis::execute()", this->computeEventLepSF( eventInfo, leptonsSorted, SFType::ID ) );
    EL_RETURN_CHECK("HTopMultilepAnalysis::execute()", this->computeEventLepSF( eventInfo, leptonsSorted, SFType::TTVA ) );
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

  if ( nSignalLeptons == 2 ) {
    if ( m_useMCForTagAndProbe && m_isMC ) { EL_RETURN_CHECK("HTopMultilepAnalysis::execute()", this->defineTagAndProbeRFRateVars_MC( eventInfo, leptonsSorted ) ); }
    else                                   { EL_RETURN_CHECK("HTopMultilepAnalysis::execute()", this->defineTagAndProbeRFRateVars( eventInfo, leptonsSorted ) ); }
  }

  //-----------------------------------
  //***** Matrix Method event weighting
  //-----------------------------------

  // first, decorate event with specific variables ...
  //
  EL_RETURN_CHECK("HTopMultilepAnalysis::execute()", this->addChannelDecorations( eventInfo, leptonsSorted ) );

  // ...then, decorate event with MM and FF weight!
  //
  EL_RETURN_CHECK("HTopMultilepAnalysis::execute()", this->fakeWeightCalculator( eventInfo, leptonsSorted ) );

  // ...and the QMisID weight!
  //
  EL_RETURN_CHECK("HTopMultilepAnalysis::execute()", this->QMisIDWeightCalculator( eventInfo, signalElectrons ) );

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

  if ( m_TauSelTool ) { m_TauSelTool = nullptr; delete m_TauSelTool; }

  /*
  // delete TH1D histograms
  //
  for ( auto itr : (m_el_hist_map) ) {
    itr.second = nullptr; delete itr.second;
  }
  for ( auto itr : (m_mu_hist_map) ) {
    itr.second = nullptr; delete itr.second;
  }
  */

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

  // fill histograms needed by Run1 merging script to normalise MC
  //
  double n_init_evts(0.0);
  double n_init_evts_W(0.0);

  // if MetaData is not empty, use it
  //
  if ( m_histEventCount->GetBinContent(1) > 0 && m_histEventCount->GetBinContent(3) > 0 ) {
    n_init_evts   =  m_histEventCount->GetBinContent(1);  // nEvents initial
    n_init_evts_W =  m_histEventCount->GetBinContent(3);  // sumOfWeights initial
  }
  // ...else, retrieve event count from cutflow
  else
  {
    if ( m_cutflowHist ) {
      int init_evts_bin    =  m_cutflowHist->GetXaxis()->FindBin("all");
      n_init_evts          =  m_cutflowHist->GetBinContent( init_evts_bin );
    }
    if ( m_cutflowHistW ) {
      int init_evts_bin_W  =  m_cutflowHistW->GetXaxis()->FindBin("all");
      n_init_evts_W        =  m_cutflowHistW->GetBinContent( init_evts_bin_W );
    }
  }
  // set the value into both bins of histogram
  m_totalEvents->SetBinContent( 1, n_init_evts );
  m_totalEvents->SetBinContent( 2, n_init_evts );
  m_totalEventsW->SetBinContent( 1, n_init_evts_W );
  m_totalEventsW->SetBinContent( 2, n_init_evts_W );

  return EL::StatusCode::SUCCESS;

}

EL::StatusCode HTopMultilepAnalysis :: readQMisIDRates ( const std::string& input_path )
{

  std::string path_AntiT = input_path + "QMisIDRates_Data_Loose.root";
  std::string path_T     = input_path + "QMisIDRates_Data_nominal_v4.root";

  TFile *file_AntiT = TFile::Open(path_AntiT.c_str());
  TFile *file_T     = TFile::Open(path_T.c_str());

  HTOP_RETURN_CHECK( "HTopMultilepAnalysis::readQMisIDRates()", file_AntiT->IsOpen(), "Failed to open ROOT file" );
  HTOP_RETURN_CHECK( "HTopMultilepAnalysis::readQMisIDRates()", file_T->IsOpen(), "Failed to open ROOT file" );

  Info("readQMisIDRates()", "Successfully opened ROOT files with QMisID rates from path:\n AntiT --> %s \n T --> %s", path_AntiT.c_str(), path_T.c_str() );

  TH2D *hist_QMisID_AntiT = get_object<TH2D>( *file_AntiT, "Rates" );
  TH2D *hist_QMisID_T     = get_object<TH2D>( *file_T, "Rates" );

  // fill a map for later usage
  //
  m_QMisID_hist_map["AntiT"] = hist_QMisID_AntiT;
  m_QMisID_hist_map["T"]     = hist_QMisID_T;

  return EL::StatusCode::SUCCESS;

}

EL::StatusCode HTopMultilepAnalysis :: readFakeRates ( const std::string& input_path )
{

  std::string path = ( !m_useMCForTagAndProbe ) ? ( input_path + "ObservedRates.root" ) : ( input_path + "ExpectedRates.root" );

  TFile *file = TFile::Open(path.c_str());

  HTOP_RETURN_CHECK( "HTopMultilepAnalysis::readFakeRates()", file->IsOpen() , "Failed to open ROOT file");

  Info("readFakeRates()", " Successfully opened ROOT file with r/f rates from path: %s ", path.c_str() );

  // 1) ELECTRONS
  //
  std::string histname_el_eta_rr   = "El_ProbeEta_Real_Rate_";
  std::string histname_el_eta_fr   = "El_ProbeEta_Fake_Rate_";
  std::string histname_el_pt_rr    = "El_ProbePt_Real_Rate_";
  std::string histname_el_pt_fr    = "El_ProbePt_Fake_Rate_";
  std::string histname_el_eta_r_T  = "El_ProbeEta_Real_T_";
  std::string histname_el_eta_r_L  = "El_ProbeEta_Real_L_";
  std::string histname_el_pt_r_T   = "El_ProbePt_Real_T_";
  std::string histname_el_pt_r_L   = "El_ProbePt_Real_L_";
  std::string histname_el_eta_f_T  = "El_ProbeEta_Fake_T_";
  std::string histname_el_eta_f_L  = "El_ProbeEta_Fake_L_";
  std::string histname_el_pt_f_T   = "El_ProbePt_Fake_T_";
  std::string histname_el_pt_f_L   = "El_ProbePt_Fake_L_";

  if ( !m_useMCForTagAndProbe ) {
    histname_el_eta_rr   += "observed";
    histname_el_eta_fr   += "observed";
    histname_el_pt_rr    += "observed";
    histname_el_pt_fr    += "observed";
    histname_el_eta_r_T  += "observed";
    histname_el_eta_r_L  += "observed";
    histname_el_pt_r_T   += "observed";
    histname_el_pt_r_L   += "observed";
    histname_el_eta_f_T  += "observed";
    histname_el_eta_f_L  += "observed";
    histname_el_pt_f_T   += "observed";
    histname_el_pt_f_L   += "observed";
  } else {
    histname_el_eta_rr   += "expected";
    histname_el_eta_fr   += "expected";
    histname_el_pt_rr    += "expected";
    histname_el_pt_fr    += "expected";
    histname_el_eta_r_T  += "expected";
    histname_el_eta_r_L  += "expected";
    histname_el_pt_r_T   += "expected";
    histname_el_pt_r_L   += "expected";
    histname_el_eta_f_T  += "expected";
    histname_el_eta_f_L  += "expected";
    histname_el_pt_f_T   += "expected";
    histname_el_pt_f_L   += "expected";
  }

  // get eta real/fake rate hist
  //
  TH1D *hist_el_eta_rr   = get_object<TH1D>( *file, histname_el_eta_rr );
  TH1D *hist_el_eta_fr   = get_object<TH1D>( *file, histname_el_eta_fr );
  TH1D *hist_el_pt_rr    = get_object<TH1D>( *file, histname_el_pt_rr );
  TH1D *hist_el_pt_fr    = get_object<TH1D>( *file, histname_el_pt_fr );
  TH1D *hist_el_eta_r_T  = get_object<TH1D>( *file, histname_el_eta_r_T );
  TH1D *hist_el_eta_r_L  = get_object<TH1D>( *file, histname_el_eta_r_L );
  TH1D *hist_el_pt_r_T   = get_object<TH1D>( *file, histname_el_pt_r_T );
  TH1D *hist_el_pt_r_L   = get_object<TH1D>( *file, histname_el_pt_r_L );
  TH1D *hist_el_eta_f_T  = get_object<TH1D>( *file, histname_el_eta_f_T );
  TH1D *hist_el_eta_f_L  = get_object<TH1D>( *file, histname_el_eta_f_L );
  TH1D *hist_el_pt_f_T   = get_object<TH1D>( *file, histname_el_pt_f_T );
  TH1D *hist_el_pt_f_L   = get_object<TH1D>( *file, histname_el_pt_f_L );

  // fill a map for later usage
  //
  m_el_hist_map["eta_rr"]   = hist_el_eta_rr;
  m_el_hist_map["eta_fr"]   = hist_el_eta_fr;
  m_el_hist_map["pt_rr"]    = hist_el_pt_rr;
  m_el_hist_map["pt_fr"]    = hist_el_pt_fr;
  m_el_hist_map["eta_r_T"]  = hist_el_eta_r_T;
  m_el_hist_map["eta_r_L"]  = hist_el_eta_r_L;
  m_el_hist_map["pt_r_T"]   = hist_el_pt_r_T;
  m_el_hist_map["pt_r_L"]   = hist_el_pt_r_L;
  m_el_hist_map["eta_f_T"]  = hist_el_eta_f_T;
  m_el_hist_map["eta_f_L"]  = hist_el_eta_f_L;
  m_el_hist_map["pt_f_T"]   = hist_el_pt_f_T;
  m_el_hist_map["pt_f_L"]   = hist_el_pt_f_L;

  // eta hist has same binning for r/f
  //
  m_n_el_bins_eta   =  hist_el_eta_rr->GetNbinsX();

  // pt hist has two different binning for r/f
  //
  m_n_el_bins_pt_rr =  hist_el_pt_rr->GetNbinsX();
  m_n_el_bins_pt_fr =  hist_el_pt_fr->GetNbinsX();

  // normalistaion factor is the same for eta and pt r/f histograms: use eta
  //
  m_el_rr_tot = ( hist_el_eta_r_T->Integral() ) / ( hist_el_eta_r_L->Integral() );
  m_el_fr_tot = ( hist_el_eta_f_T->Integral() ) / ( hist_el_eta_f_L->Integral() );

  // 2) MUONS
  //
  std::string histname_mu_eta_rr   = "Mu_ProbeEta_Real_Rate_";
  std::string histname_mu_eta_fr   = "Mu_ProbeEta_Fake_Rate_";
  std::string histname_mu_pt_rr    = "Mu_ProbePt_Real_Rate_";
  std::string histname_mu_pt_fr    = "Mu_ProbePt_Fake_Rate_";
  std::string histname_mu_eta_r_T  = "Mu_ProbeEta_Real_T_";
  std::string histname_mu_eta_r_L  = "Mu_ProbeEta_Real_L_";
  std::string histname_mu_pt_r_T   = "Mu_ProbePt_Real_T_";
  std::string histname_mu_pt_r_L   = "Mu_ProbePt_Real_L_";
  std::string histname_mu_eta_f_T  = "Mu_ProbeEta_Fake_T_";
  std::string histname_mu_eta_f_L  = "Mu_ProbeEta_Fake_L_";
  std::string histname_mu_pt_f_T   = "Mu_ProbePt_Fake_T_";
  std::string histname_mu_pt_f_L   = "Mu_ProbePt_Fake_L_";

  if ( !m_useMCForTagAndProbe ) {
    histname_mu_eta_rr   += "observed";
    histname_mu_eta_fr   += "observed";
    histname_mu_pt_rr    += "observed";
    histname_mu_pt_fr    += "observed";
    histname_mu_eta_r_T  += "observed";
    histname_mu_eta_r_L  += "observed";
    histname_mu_pt_r_T   += "observed";
    histname_mu_pt_r_L   += "observed";
    histname_mu_eta_f_T  += "observed";
    histname_mu_eta_f_L  += "observed";
    histname_mu_pt_f_T   += "observed";
    histname_mu_pt_f_L   += "observed";
  } else {
    histname_mu_eta_rr   += "expected";
    histname_mu_eta_fr   += "expected";
    histname_mu_pt_rr    += "expected";
    histname_mu_pt_fr    += "expected";
    histname_mu_eta_r_T  += "expected";
    histname_mu_eta_r_L  += "expected";
    histname_mu_pt_r_T   += "expected";
    histname_mu_pt_r_L   += "expected";
    histname_mu_eta_f_T  += "expected";
    histname_mu_eta_f_L  += "expected";
    histname_mu_pt_f_T   += "expected";
    histname_mu_pt_f_L   += "expected";
  }

  // get eta real/fake rate hist
  //
  TH1D *hist_mu_eta_rr	 = get_object<TH1D>( *file, histname_mu_eta_rr );
  TH1D *hist_mu_eta_fr	 = get_object<TH1D>( *file, histname_mu_eta_fr );
  TH1D *hist_mu_pt_rr	 = get_object<TH1D>( *file, histname_mu_pt_rr );
  TH1D *hist_mu_pt_fr	 = get_object<TH1D>( *file, histname_mu_pt_fr );
  TH1D *hist_mu_eta_r_T  = get_object<TH1D>( *file, histname_mu_eta_r_T );
  TH1D *hist_mu_eta_r_L  = get_object<TH1D>( *file, histname_mu_eta_r_L );
  TH1D *hist_mu_pt_r_T	 = get_object<TH1D>( *file, histname_mu_pt_r_T );
  TH1D *hist_mu_pt_r_L	 = get_object<TH1D>( *file, histname_mu_pt_r_L );
  TH1D *hist_mu_eta_f_T  = get_object<TH1D>( *file, histname_mu_eta_f_T );
  TH1D *hist_mu_eta_f_L  = get_object<TH1D>( *file, histname_mu_eta_f_L );
  TH1D *hist_mu_pt_f_T	 = get_object<TH1D>( *file, histname_mu_pt_f_T );
  TH1D *hist_mu_pt_f_L	 = get_object<TH1D>( *file, histname_mu_pt_f_L );

  // fill a map for later usage
  //
  m_mu_hist_map["eta_rr"]   = hist_mu_eta_rr;
  m_mu_hist_map["eta_fr"]   = hist_mu_eta_fr;
  m_mu_hist_map["pt_rr"]    = hist_mu_pt_rr;
  m_mu_hist_map["pt_fr"]    = hist_mu_pt_fr;
  m_mu_hist_map["eta_r_T"]  = hist_mu_eta_r_T;
  m_mu_hist_map["eta_r_L"]  = hist_mu_eta_r_L;
  m_mu_hist_map["pt_r_T"]   = hist_mu_pt_r_T;
  m_mu_hist_map["pt_r_L"]   = hist_mu_pt_r_L;
  m_mu_hist_map["eta_f_T"]  = hist_mu_eta_f_T;
  m_mu_hist_map["eta_f_L"]  = hist_mu_eta_f_L;
  m_mu_hist_map["pt_f_T"]   = hist_mu_pt_f_T;
  m_mu_hist_map["pt_f_L"]   = hist_mu_pt_f_L;

  // eta hist has same binning for r/f
  //
  m_n_mu_bins_eta   =  hist_mu_eta_rr->GetNbinsX();

  // pt hist has two different binning for r/f
  //
  m_n_mu_bins_pt_rr =  hist_mu_pt_rr->GetNbinsX();
  m_n_mu_bins_pt_fr =  hist_mu_pt_fr->GetNbinsX();

  // normalistaion factor is the same for eta and pt r/f histograms: use eta
  //
  m_mu_rr_tot = ( hist_mu_eta_r_T->Integral() ) / ( hist_mu_eta_r_L->Integral() );
  m_mu_fr_tot = ( hist_mu_eta_f_T->Integral() ) / ( hist_mu_eta_f_L->Integral() );

  // ***********************************************************

  return EL::StatusCode::SUCCESS;

}

//***************************************************
//
// Set Tag&Probe variables in r/f rate measurement CR
//
//
EL::StatusCode HTopMultilepAnalysis :: defineTagAndProbeRFRateVars( const xAOD::EventInfo* eventInfo, const xAOD::IParticleContainer& leptons )
{

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
  static SG::AuxElement::Accessor< char > isTrigMatchedLepAcc("isTrigMatchedLep");

  // decorate with default values
  for ( auto lep_it : leptons ) { isTagDecor( *lep_it ) = 0; }
  isNonTightEventDecor( *eventInfo ) = 0;
  isProbeElEventDecor( *eventInfo )  = 0;
  isProbeMuEventDecor( *eventInfo )  = 0;

  // -------------------------------------------
  // Now, take the leading and subleading lepton
  // -------------------------------------------

  const xAOD::IParticle* leadingLepton = leptons.at(0);

  // --------------------------------------------
  // Now decide who's the tag and who's the probe
  // --------------------------------------------

  static SG::AuxElement::Accessor< char > isChFlipAcc("isChFlip");
  static SG::AuxElement::ConstAccessor< int > truthTypeAcc("truthType");
  static SG::AuxElement::ConstAccessor< int > truthOriginAcc("truthOrigin");

  // The first lepton found that is tight && trigger-matched will be the tag, the other the probe
  // Note that the lepton container is sorted, so in case both are T & TM, the leading will be the tag and the subleading the probe
  // If none satisfies the above, choose the leading as tight and the subleading as probe, but flag the event as "bad"
  //
  // Use the same logic for OS and SS events

  bool found_tag(false);

  for ( auto lep_itr : leptons ) {
    if ( m_debug ) { Info("defineTagAndProbeRFRateVars()","lepton pT = %f", lep_itr->pt()/1e3 ); }
    if ( isTightAcc.isAvailable( *lep_itr ) && isTrigMatchedLepAcc.isAvailable( *lep_itr ) ) {
      if ( !found_tag && ( isTightAcc( *lep_itr ) == 1 && isTrigMatchedLepAcc( *lep_itr ) == 1 ) ) {
	isTagDecor( *lep_itr ) = 1;
	found_tag = true;
        if ( m_debug ) { Info("defineTagAndProbeRFRateVars()","\t ===> found tag!"); }
      }
    }
  }

  if ( !found_tag ) {
    isTagDecor( *leadingLepton ) = 1;
    // take note that this event should not be used
    isNonTightEventDecor( *eventInfo ) = 1;
    if ( m_debug ) { Info("defineTagAndProbeRFRateVars()","None lepton is T&TM - choose leading as tag (pT = %f)",leadingLepton->pt()/1e3 ); }
  }

  // Accessor to tag leptons
  //
  static SG::AuxElement::Accessor< char > isTagAcc("isTag");

  // Now loop over the leptons, and check whether the probe is an electron or a muon
  //
  for ( auto lep_itr : leptons ) {
    if ( !isTagAcc( *lep_itr ) ) {
      isProbeElEventDecor( *eventInfo ) = ( lep_itr->type() == xAOD::Type::Electron );
      isProbeMuEventDecor( *eventInfo ) = ( lep_itr->type() == xAOD::Type::Muon );
    }
  }

  if ( m_debug ) {
    Info("defineTagAndProbeRFRateVars()"," ********** Checking 'isTag' lepton decoration ********** ");
    for ( auto lep_itr : leptons ) {
      Info("defineTagAndProbeRFRateVars()","\t lepton \n \t isTag?: %i \n \t truthType(): %i \n \t truthOrigin(): %i \n ", isTagAcc( *lep_itr ), truthTypeAcc( *lep_itr ), truthOriginAcc( *lep_itr ) );
    }
    Info("defineTagAndProbeRFRateVars()"," ********** ");
  }

  return EL::StatusCode::SUCCESS;

}

// ******************************************
// Use this function only when doing
// MM estimate on pure MC (i.e, closure test)
//
EL::StatusCode HTopMultilepAnalysis :: defineTagAndProbeRFRateVars_MC( const xAOD::EventInfo* eventInfo, const xAOD::IParticleContainer& leptons )
{

  // -------------------------------
  // Categorise event based on OS/SS
  // -------------------------------

  int prod_lep_charge(1);
  for ( auto lep_it : leptons ) {
    // get the lepton flavour
    xAOD::Type::ObjectType leptonFlavour = lep_it->type();
    if ( leptonFlavour == xAOD::Type::Electron ) {
      prod_lep_charge *= dynamic_cast<const xAOD::Electron*>( lep_it )->charge();
    } else if ( leptonFlavour == xAOD::Type::Muon ) {
      prod_lep_charge *= dynamic_cast<const xAOD::Muon*>( lep_it )->charge();
    }
  }
  bool isSS = (  prod_lep_charge > 0  );

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

  // decorate with default values
  for ( auto lep_itr : leptons ) { isTagDecor( *lep_itr ) = 0; }
  isNonTightEventDecor( *eventInfo ) = 0;
  isProbeElEventDecor( *eventInfo )  = 0;
  isProbeMuEventDecor( *eventInfo )  = 0;

  // -------------------------------------------
  // Now, take the leading and subleading lepton
  // -------------------------------------------

  const xAOD::IParticle* leadingLepton           = leptons.at(0);

  // --------------------------------------------
  // Now decide who's the tag and who's the probe
  // --------------------------------------------

  static SG::AuxElement::Accessor< char > isChFlipAcc("isChFlip");
  static SG::AuxElement::ConstAccessor< int > truthTypeAcc("truthType");
  static SG::AuxElement::ConstAccessor< int > truthOriginAcc("truthOrigin");
  static SG::AuxElement::Accessor< char > isTightAcc("isTight");
  static SG::AuxElement::Accessor< char > isTrigMatchedLepAcc("isTrigMatchedLep");

  bool found_tag(false);
  if ( isSS ) {

    for ( auto lep_itr : leptons ) {

      // In the SS case, the tag will be the first prompt lepton in the event which is not charge flip, provided it's found.
      // The other will be the probe
      // See below the treatment for the case where all leptons in SS event are non prompt...
      //
      if ( truthTypeAcc.isAvailable( *lep_itr ) ) {

        if ( !found_tag && ( truthTypeAcc( *lep_itr ) == 2 || truthTypeAcc( *lep_itr ) == 6 ) && isChFlipAcc( *lep_itr ) != 1 ) {
          isTagDecor( *lep_itr ) = 1;
	  found_tag = true;
        }

      } else {
        Warning("defineTagAndProbeRFRateVars_MC()","SG::AuxElement::Accessor('truthType') is not available for this lepton. Should not happen. Tag will be the leading, probe the subleading" );
        break;
      }

    }

    if ( !found_tag ) {
      // the probe will be the subleading lepton
      //
      if ( m_debug ) { Info("defineTagAndProbeRFRateVars_MC()","There are no prompt leptons in this SS event. Tag will be the leading, probe the subleading" ); }
      isTagDecor( *leadingLepton ) = 1;
      // take note that this event should not be used
      isNonTightEventDecor( *eventInfo ) = 1;
    }

  } else {

    // OS case
    // The first lepton found that is tight && trigger-matched will be the tag, the other the probe
    // Note that the lepton container is sorted, so in case both are T & TM, the leading will be the tag and the subleading the probe
    // If none satisfies the above, choose the leading as tight and the subleading as probe, but flag the event as "bad"
    //

    bool found_tag(false);

    for ( auto lep_itr : leptons ) {
      if ( isTightAcc.isAvailable( *lep_itr ) && isTrigMatchedLepAcc.isAvailable( *lep_itr ) ) {
	if ( !found_tag && ( isTightAcc( *lep_itr ) == 1 && isTrigMatchedLepAcc( *lep_itr ) == 1 ) ) {
	  isTagDecor( *lep_itr ) = 1;
	  found_tag = true;
	}
      }
    }

    if ( !found_tag ) {
      isTagDecor( *leadingLepton ) = 1;
      // take note that this event should not be used
      isNonTightEventDecor( *eventInfo ) = 1;
    }

  }

  // Accessor to tag leptons
  //
  static SG::AuxElement::Accessor< char > isTagAcc("isTag");

  // Now loop over the leptons, and check whether the probe is an electron or a muon
  //
  for ( auto lep_itr : leptons ) {

    if ( !isTagAcc( *lep_itr ) ) {
      isProbeElEventDecor( *eventInfo ) = ( lep_itr->type() == xAOD::Type::Electron );
      isProbeMuEventDecor( *eventInfo ) = ( lep_itr->type() == xAOD::Type::Muon );
    }

  }

  if ( m_debug && isSS ) {
    Info("defineTagAndProbeRFRateVars_MC()"," ********** SS Event - Checking 'isTag' lepton decoration ********** ");
    for ( auto lep_itr : leptons ) {
      Info("defineTagAndProbeRFRateVars_MC()","\t lepton \n \t isTag?: %i \n \t truthType(): %i \n \t truthOrigin(): %i \n ", isTagAcc( *lep_itr ), truthTypeAcc( *lep_itr ), truthOriginAcc( *lep_itr ) );
    }
    Info("defineTagAndProbeRFRateVars_MC()"," ********** ");
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

  static SG::AuxElement::Decorator< float > mOSPair01Decor("mOSPair01"); // for events with 1 OS (lep0) and 2 SS (lep1,lep2), save m(lep0,lep1)
  static SG::AuxElement::Decorator< float > mOSPair02Decor("mOSPair02"); // for events with 1 OS (lep0) and 2 SS (lep1,lep2), save m(lep0,lep2)
  static SG::AuxElement::Decorator< char > isOSPairSF01Decor("isOSPairSF01"); // for events with 1 OS and 2 SS, check whether lep0,lep1 are same-flavour
  static SG::AuxElement::Decorator< char > isOSPairSF02Decor("isOSPairSF02"); // for events with 1 OS and 2 SS, check whether lep0,lep1 are same-flavour

  // start decorating with default values
  isSS01Decor(*eventInfo)  = 0;
  isSS12Decor(*eventInfo)  = 0;
  mll01Decor(*eventInfo)   = -1.0;
  mll02Decor(*eventInfo)   = -1.0;
  mll12Decor(*eventInfo)   = -1.0;
  mlll012Decor(*eventInfo) = -1.0;

  mOSPair01Decor(*eventInfo) = -1.0;
  mOSPair02Decor(*eventInfo) = -1.0;
  isOSPairSF01Decor(*eventInfo) = 0;
  isOSPairSF02Decor(*eventInfo) = 0;

  unsigned int nLeptons = leptons.size();

  // Set defaults
  if ( nLeptons > 0 ) {
    for( auto lep_it : leptons ) {
      isOSlepDecor(*lep_it) = 0;
      isClosestSSlepDecor(*lep_it) = 0;
      isSS12Decor(*lep_it)  = 0;
    }
  }

  // initialize TLorentzVectors (for invariant mass computation)
  xAOD::IParticle::FourMom_t lepA_4mom, lepB_4mom, lepC_4mom;

  if ( nLeptons == 2 ) {

      // retrieve lepA
     const xAOD::IParticle* lepA = leptons.at(0);
     // retrieve lepB
     const xAOD::IParticle* lepB = leptons.at(1);

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

     // mT( lep, MET )
     //
     const xAOD::MissingETContainer* inMETCont(nullptr);
     RETURN_CHECK("HTopMultilepAnalysis::execute()", HelperFunctions::retrieve(inMETCont, "RefFinal_HTopMultilep", m_event, m_store, m_verbose), "");

     static SG::AuxElement::Decorator<float> mTLep0METDecor("mT_lep0MET");
     static SG::AuxElement::Decorator<float> mTLep1METDecor("mT_lep1MET");

     const xAOD::MissingET* final = *inMETCont->find("FinalTrk"); // ("FinalClus" uses the calocluster-based soft terms, "FinalTrk" uses the track-based ones)
     TLorentzVector MET;
     MET.SetPtEtaPhiM( final->met(), 0, final->phi(), 0 );

     float mT(-1.0);
     mT = sqrt( 2.0 * lepA_4mom.Pt() * ( MET.Pt() ) * ( 1.0 - cos( lepA_4mom.DeltaPhi(MET) ) ) );
     mTLep0METDecor( *eventInfo ) = mT;
     mT = sqrt( 2.0 * lepB_4mom.Pt() * ( MET.Pt() ) * ( 1.0 - cos( lepB_4mom.DeltaPhi(MET) ) ) );
     mTLep1METDecor( *eventInfo ) = mT;

  } else if ( nLeptons == 3 ) {

      // retrieve lepA
     const xAOD::IParticle* lepA = leptons.at(0);
     // retrieve lepB
     const xAOD::IParticle* lepB = leptons.at(1);
     // retrieve lepC
     const xAOD::IParticle* lepC = leptons.at(2);

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

     // look only at events where there is (at least) one OS pair
     //
     if ( fabs( lepAcharge/fabs(lepAcharge) + lepBcharge/fabs(lepBcharge) + lepCcharge/fabs(lepCcharge)  ) != 3 ) {

	 // now decorate event!
	isSS12Decor(*eventInfo) = 1;

        // now find the OS lepton and decorate!
        isOSlepDecor(*lepA) = ( (lepBcharge * lepCcharge) > 0 );
        isOSlepDecor(*lepB) = ( (lepAcharge * lepCcharge) > 0 );
        isOSlepDecor(*lepC) = ( (lepAcharge * lepBcharge) > 0 );

        // once OS lepton has been found, find the SS lep with min{ deltaR(lep0) } and decorate!
        //
	// Save the invariant mass of all the lepton pairs containing the OS one, and whether the pair is
        // made up of same flavour leptons
        //
	float dR_i(-1), dR_j(-1);
	if ( isOSlepDecor(*lepA) ) {
	   dR_i = lepA_4mom.DeltaR(lepB_4mom );
	   dR_j = lepA_4mom.DeltaR(lepC_4mom );
	   if ( dR_i <=  dR_j ) {
	       isClosestSSlepDecor(*lepB) = 1;
	       mOSPair01Decor( *eventInfo ) = pairAB.M();
	       mOSPair02Decor( *eventInfo ) = pairAC.M();
	       isOSPairSF01Decor( *eventInfo ) = ( lepAFlavour == lepBFlavour );
	       isOSPairSF02Decor( *eventInfo ) = ( lepAFlavour == lepCFlavour );
	   } else {
	       isClosestSSlepDecor(*lepC) = 1;
	       mOSPair01Decor( *eventInfo ) = pairAC.M();
	       mOSPair02Decor( *eventInfo ) = pairAB.M();
	       isOSPairSF01Decor( *eventInfo ) = ( lepAFlavour == lepCFlavour );
	       isOSPairSF02Decor( *eventInfo ) = ( lepAFlavour == lepBFlavour );
	   }
	}
	else if ( isOSlepDecor(*lepB) ) {
	   dR_i = lepB_4mom.DeltaR(lepA_4mom );
	   dR_j = lepB_4mom.DeltaR(lepC_4mom );
	   if ( dR_i <=  dR_j ) {
	       isClosestSSlepDecor(*lepA) = 1;
	       mOSPair01Decor( *eventInfo ) = pairAB.M();
	       mOSPair02Decor( *eventInfo ) = pairBC.M();
	       isOSPairSF01Decor( *eventInfo ) = ( lepAFlavour == lepBFlavour );
	       isOSPairSF02Decor( *eventInfo ) = ( lepBFlavour == lepCFlavour );
	   } else {
	       isClosestSSlepDecor(*lepC) = 1;
	       mOSPair01Decor( *eventInfo ) = pairBC.M();
	       mOSPair02Decor( *eventInfo ) = pairAB.M();
	       isOSPairSF01Decor( *eventInfo ) = ( lepBFlavour == lepCFlavour );
	       isOSPairSF02Decor( *eventInfo ) = ( lepAFlavour == lepBFlavour );
	   }
	}
	else if ( isOSlepDecor(*lepC) ) {
	   dR_i = lepC_4mom.DeltaR(lepA_4mom );
	   dR_j = lepC_4mom.DeltaR(lepB_4mom );
	   if ( dR_i <=  dR_j ) {
	       isClosestSSlepDecor(*lepA) = 1;
	       mOSPair01Decor( *eventInfo ) = pairAC.M();
	       mOSPair02Decor( *eventInfo ) = pairBC.M();
	       isOSPairSF01Decor( *eventInfo ) = ( lepAFlavour == lepCFlavour );
	       isOSPairSF02Decor( *eventInfo ) = ( lepBFlavour == lepCFlavour );
	   } else {
	       isClosestSSlepDecor(*lepB) = 1;
	       mOSPair01Decor( *eventInfo ) = pairBC.M();
	       mOSPair02Decor( *eventInfo ) = pairAC.M();
	       isOSPairSF01Decor( *eventInfo ) = ( lepBFlavour == lepCFlavour );
	       isOSPairSF02Decor( *eventInfo ) = ( lepAFlavour == lepCFlavour );
	   }
	}
     }
  }

  return EL::StatusCode::SUCCESS;
}

//*****************************************************************************
//
// Matrix Method & Fake Factor Method stuff
//
//

EL::StatusCode HTopMultilepAnalysis :: fakeWeightCalculator (const xAOD::EventInfo* eventInfo, const xAOD::IParticleContainer& leptons )
{
  // This method calculates the final weights to be applied to (for MM) data events
  // with leptons TT, TL, LT, LL in order to obtain the fake estimate
  // in the TT signal/control region
  // We manually plug in the lepton fake/real rates calculated via tag-and-probe analysis (FIXME)
  //
  //NB: MM and FF Method are applied only to the 2lepSS and the 3lep case when there are 2 SS leptons

  unsigned int nLeptons = leptons.size();

  // Make sure nothing is done if there are less than 2 leptons...
  //
  if ( nLeptons < 2 ) { return EL::StatusCode::SUCCESS; }

  // retrieve some previously applied event object decorations
  //
  static SG::AuxElement::Accessor< char > isSS01("isSS01");
  static SG::AuxElement::Accessor< char > isSS12("isSS12");

  HTOP_RETURN_CHECK( "HTopMultilepAnalysis::fakeWeightCalculator()", isSS01.isAvailable(*eventInfo), "isSS01 event decoration is not available. Aborting");
  HTOP_RETURN_CHECK( "HTopMultilepAnalysis::fakeWeightCalculator()", isSS12.isAvailable(*eventInfo), "isSS12 event decoration is not available. Aborting");


  // **********************************************************************
  //
  // These will be the two leptons used for the fake estimate and to define
  // the "tightness" of the region, both in 2 lep SS and in 3 lep category
  //
  xAOD::IParticle* lepA(nullptr);
  xAOD::IParticle* lepB(nullptr);

  // The "real" lepton for 3 lep category
  //
  //xAOD::IParticle* lep0(nullptr);

  // Features of the two leptons that are considered for the fake estimate
  //
  float lepA_pt(-1.0), lepA_eta(-999.0), lepB_pt(-1.0), lepB_eta(-999.0);
  int lepA_flavour(0), lepB_flavour(0);

  if ( nLeptons == 2 ) {

    // start from lepton container
    //
    // ordering criterion is simply based on pT
    // by construction, the first element of the DV is the leading lepton, the second (and last!) element is the subleading

    // retrieve lepA : the leading lepton
    //
    lepA = const_cast<xAOD::IParticle*>(leptons.at(0));
    // retrieve lepB: the subleading lepton
    //
    lepB = const_cast<xAOD::IParticle*>(leptons.at(1));

  } else if ( nLeptons == 3 && isSS12( *eventInfo ) ) {

    // start from lepton container
    //
    // for trilepton, ordering criterion is:
    // lep0: the OS lepton (assume this is real!)
    // lepA: the SS lepton with min{ deltaR(lep0) }
    // lepB: the other SS lepton
    // --> lepA and lepB define the "tightness" of the region

    // retrieve some previously applied lepton object decorations
    //
    static SG::AuxElement::Accessor< char > isOSlep("isOSlep");
    static SG::AuxElement::Accessor< char > isClosestSSlep("isClosestSSlep");

    for (auto lep_it :leptons ) {

      xAOD::IParticle* this_lep = const_cast<xAOD::IParticle*>(lep_it);

      HTOP_RETURN_CHECK( "HTopMultilepAnalysis::fakeWeightCalculator()", isOSlep.isAvailable(*lep_it), "isOSlep lepton decoration is not available. Aborting");

      if ( isOSlep(*this_lep) ) {

    	// retrieve lep0 : the OS lepton
    	//
    	//lep0 = this_lep;
    	continue;

      } else {

	  HTOP_RETURN_CHECK( "HTopMultilepAnalysis::fakeWeightCalculator()", isClosestSSlep.isAvailable(*lep_it), "isClosestSSlep lepton decoration is not available. Aborting");

    	// retrieve lepA : the SS lepton with min{ deltaR(lep0) }
    	//
    	if ( isClosestSSlep(*this_lep) ) {
    	  lepA = this_lep;
    	  continue;
    	}
    	// retrieve lepB : the other SS lepton
    	//
    	lepB = this_lep;

      }

    } // close loop over lepton container

  } else {
    return EL::StatusCode::SUCCESS;
  }

  // just a safety check
  //
  HTOP_RETURN_CHECK( "HTopMultilepAnalysis::fakeWeightCalculator()", ( lepA && lepB ), "2Lep/3Lep(12SS) region, but no (lepA and lepB) pointers! Aborting" );

  // set the properties of the two relevant leptons for future convenience
  //
  lepA_pt  = lepA->pt();
  lepA_eta = lepA->eta();
  if ( lepA->type() == xAOD::Type::Electron )  { lepA_flavour = 11; }
  else if ( lepA->type() == xAOD::Type::Muon ) { lepA_flavour = 13; }

  lepB_pt  = lepB->pt();
  lepB_eta = lepB->eta();
  if ( lepB->type() == xAOD::Type::Electron )  { lepB_flavour = 11; }
  else if ( lepB->type() == xAOD::Type::Muon ) { lepB_flavour = 13; }

  bool OF = ( ( lepA_flavour == 11 && lepB_flavour == 13 ) || ( lepA_flavour == 13 && lepB_flavour == 11 ) );

  // ******************************************************************************************
  //
  // Set the region string for MM/FF based on tightness. This depends on the
  // choice of loosest and tightest definition made at configuration.
  //
  // Regardless of this choice, save an event flag with the T-M-L composition of the
  // dilepton (or dilepton+1lep) events
  //
  // For the theta method, specify the flavour of the two leptons in TL,LT regions

  // event decorators to identify the regions based on lepton tightness (lepton ordering criterion depends on the category)
  //
  static SG::AuxElement::Decorator< char > is_T_T_Decor("is_T_T");
  static SG::AuxElement::Decorator< char > is_T_AntiT_Decor("is_T_AntiT");
  static SG::AuxElement::Decorator< char > is_AntiT_T_Decor("is_AntiT_T");
  static SG::AuxElement::Decorator< char > is_AntiT_AntiT_Decor("is_AntiT_AntiT");

  is_T_T_Decor( *eventInfo ) = 0;
  is_T_AntiT_Decor( *eventInfo ) = 0;
  is_AntiT_T_Decor( *eventInfo ) = 0;
  is_AntiT_AntiT_Decor( *eventInfo ) = 0;

  static SG::AuxElement::Decorator< char > is_T_MAntiT_Decor("is_T_MAntiT");
  static SG::AuxElement::Decorator< char > is_MAntiT_T_Decor("is_MAntiT_T");
  static SG::AuxElement::Decorator< char > is_MAntiT_MAntiT_Decor("is_MAntiT_MAntiT");

  is_T_MAntiT_Decor( *eventInfo ) = 0;
  is_MAntiT_T_Decor( *eventInfo ) = 0;
  is_MAntiT_MAntiT_Decor( *eventInfo ) = 0;

  static SG::AuxElement::Decorator< char > is_M_M_Decor("is_M_M");
  static SG::AuxElement::Decorator< char > is_T_AntiM_Decor("is_T_AntiM");
  static SG::AuxElement::Decorator< char > is_AntiM_T_Decor("is_AntiM_T");
  static SG::AuxElement::Decorator< char > is_M_AntiM_Decor("is_M_AntiM");
  static SG::AuxElement::Decorator< char > is_AntiM_M_Decor("is_AntiM_M");
  static SG::AuxElement::Decorator< char > is_AntiM_AntiM_Decor("is_AntiM_AntiM");

  is_M_M_Decor( *eventInfo ) = 0;
  is_T_AntiM_Decor( *eventInfo ) = 0;
  is_AntiM_T_Decor( *eventInfo ) = 0;
  is_M_AntiM_Decor( *eventInfo ) = 0;
  is_AntiM_M_Decor( *eventInfo ) = 0;
  is_AntiM_AntiM_Decor( *eventInfo ) = 0;

  // Save a special set of flags for TL,LT OF events (used for "theta factior" ABCD method)
  //
  static SG::AuxElement::Decorator< char > is_Tel_AntiTmu_Decor("is_Tel_AntiTmu");
  static SG::AuxElement::Decorator< char > is_Tmu_AntiTel_Decor("is_Tmu_AntiTel");
  static SG::AuxElement::Decorator< char > is_AntiTel_Tmu_Decor("is_AntiTel_Tmu");
  static SG::AuxElement::Decorator< char > is_AntiTmu_Tel_Decor("is_AntiTmu_Tel");
  is_Tel_AntiTmu_Decor( *eventInfo ) = 0;
  is_Tmu_AntiTel_Decor( *eventInfo ) = 0;
  is_AntiTel_Tmu_Decor( *eventInfo ) = 0;
  is_AntiTmu_Tel_Decor( *eventInfo ) = 0;
  static SG::AuxElement::Decorator< char > is_Tel_MAntiTmu_Decor("is_Tel_MAntiTmu");
  static SG::AuxElement::Decorator< char > is_Tmu_MAntiTel_Decor("is_Tmu_MAntiTel");
  static SG::AuxElement::Decorator< char > is_MAntiTel_Tmu_Decor("is_MAntiTel_Tmu");
  static SG::AuxElement::Decorator< char > is_MAntiTmu_Tel_Decor("is_MAntiTmu_Tel");
  is_Tel_MAntiTmu_Decor( *eventInfo ) = 0;
  is_Tmu_MAntiTel_Decor( *eventInfo ) = 0;
  is_MAntiTel_Tmu_Decor( *eventInfo ) = 0;
  is_MAntiTmu_Tel_Decor( *eventInfo ) = 0;

  static SG::AuxElement::Decorator< char > is_Mel_AntiMmu_Decor("is_Mel_AntiMmu");
  static SG::AuxElement::Decorator< char > is_Mmu_AntiMel_Decor("is_Mmu_AntiMel");
  static SG::AuxElement::Decorator< char > is_AntiMel_Mmu_Decor("is_AntiMel_Mmu");
  static SG::AuxElement::Decorator< char > is_AntiMmu_Mel_Decor("is_AntiMmu_Mel");
  is_Mel_AntiMmu_Decor( *eventInfo ) = 0;
  is_Mmu_AntiMel_Decor( *eventInfo ) = 0;
  is_AntiMel_Mmu_Decor( *eventInfo ) = 0;
  is_AntiMmu_Mel_Decor( *eventInfo ) = 0;

  static SG::AuxElement::Decorator< char > is_Tel_AntiMmu_Decor("is_Tel_AntiMmu");
  static SG::AuxElement::Decorator< char > is_Tmu_AntiMel_Decor("is_Tmu_AntiMel");
  static SG::AuxElement::Decorator< char > is_AntiMel_Tmu_Decor("is_AntiMel_Tmu");
  static SG::AuxElement::Decorator< char > is_AntiMmu_Tel_Decor("is_AntiMmu_Tel");
  is_Tel_AntiMmu_Decor( *eventInfo ) = 0;
  is_Tmu_AntiMel_Decor( *eventInfo ) = 0;
  is_AntiMel_Tmu_Decor( *eventInfo ) = 0;
  is_AntiMmu_Tel_Decor( *eventInfo ) = 0;

  std::string region("");

  // Accessor to tight/medium selected leptons
  //
  static SG::AuxElement::Accessor< char > isTightAcc("isTight");
  static SG::AuxElement::Accessor< char > isMediumAcc("isMedium");

  // the TT region is defined in this way EXCEPT when "Medium" defines our tightest definition
  //
  if (  isTightAcc( *lepA ) &&  isTightAcc( *lepB ) ) {
      if ( !( m_useLooseAsLoosest && m_useMediumAsTightest ) ) { region = "TT"; }
      is_T_T_Decor( *eventInfo ) = 1;
  }

  // -) "Loose" defines the loosest definition, "Tight" defines the tightest
  //
  if (  isTightAcc( *lepA ) && !isTightAcc( *lepB ) ) {
      if ( m_useLooseAsLoosest && m_useTightAsTightest ) { region = "TL"; }
      is_T_AntiT_Decor( *eventInfo ) = 1;
      if ( OF ) {
	  // the tight is an electron, the loose is a muon
	  if (  lepA_flavour == 11 )     { is_Tel_AntiTmu_Decor( *eventInfo ) = 1; }
	  // the tight is a muon, the loose is an electron
	  else if ( lepA_flavour == 13 ) { is_Tmu_AntiTel_Decor( *eventInfo ) = 1; }
      }
  }
  if ( !isTightAcc( *lepA ) &&  isTightAcc( *lepB ) ) {
      if ( m_useLooseAsLoosest && m_useTightAsTightest ) { region = "LT"; }
      is_AntiT_T_Decor( *eventInfo ) = 1;
      if ( OF ) {
	  if (  lepA_flavour == 11 )     { is_AntiTel_Tmu_Decor( *eventInfo ) = 1; }
	  else if ( lepA_flavour == 13 ) { is_AntiTmu_Tel_Decor( *eventInfo ) = 1; }
      }
  }
  if ( !isTightAcc( *lepA ) && !isTightAcc( *lepB ) ) {
      if ( m_useLooseAsLoosest && m_useTightAsTightest ) { region = "LL"; }
      is_AntiT_AntiT_Decor( *eventInfo ) = 1;
  }

  // -) "Medium" defines the loosest definition, "Tight" defines the tightest
  //
  if (  isTightAcc( *lepA ) && ( isMediumAcc( *lepB ) && !isTightAcc( *lepB ) ) ) {
      if ( m_useMediumAsLoosest && m_useTightAsTightest ) { region = "TL"; }
      is_T_MAntiT_Decor( *eventInfo ) = 1;
      if ( OF ) {
	  if (  lepA_flavour == 11 )     { is_Tel_MAntiTmu_Decor( *eventInfo ) = 1; }
	  else if ( lepA_flavour == 13 ) { is_Tmu_MAntiTel_Decor( *eventInfo ) = 1; }
      }
  }
  if ( ( isMediumAcc( *lepA ) && !isTightAcc( *lepA ) ) && isTightAcc( *lepB ) ) {
      if ( m_useMediumAsLoosest && m_useTightAsTightest ) { region = "LT"; }
      is_MAntiT_T_Decor( *eventInfo ) = 1;
      if ( OF ) {
	  if (  lepA_flavour == 11 )     { is_MAntiTel_Tmu_Decor( *eventInfo ) = 1; }
	  else if ( lepA_flavour == 13 ) { is_MAntiTmu_Tel_Decor( *eventInfo ) = 1; }
      }
  }
  if ( ( isMediumAcc( *lepA ) && !isTightAcc( *lepA ) ) && ( isMediumAcc( *lepB ) && !isTightAcc( *lepB ) ) ) {
      if ( m_useMediumAsLoosest && m_useTightAsTightest ) { region = "LL"; }
      is_MAntiT_MAntiT_Decor( *eventInfo ) = 1;
  }

  // -) "Loose" defines the loosest definition, "Tight" defines the tightest, and exclude all the "Medium-!Tight" leptons from the loose set
  //
  if (  isTightAcc( *lepA ) && !isMediumAcc( *lepB )  ) {
      if ( ( m_useLooseAsLoosest && m_vetoMediumNonTight ) && m_useTightAsTightest ) { region = "TL"; }
      is_T_AntiM_Decor( *eventInfo ) = 1;
      if ( OF ) {
	  if (  lepA_flavour == 11 )     { is_Tel_AntiMmu_Decor( *eventInfo ) = 1; }
	  else if ( lepA_flavour == 13 ) { is_Tmu_AntiMel_Decor( *eventInfo ) = 1; }
      }
  }
  if (  !isMediumAcc( *lepA ) && isTightAcc( *lepB ) ) {
      if ( ( m_useLooseAsLoosest && m_vetoMediumNonTight ) && m_useTightAsTightest ) { region = "LT"; }
      is_AntiM_T_Decor( *eventInfo ) = 1;
      if ( OF ) {
	  if (  lepA_flavour == 11 )     { is_AntiMel_Tmu_Decor( *eventInfo ) = 1; }
	  else if ( lepA_flavour == 13 ) { is_AntiMmu_Tel_Decor( *eventInfo ) = 1; }
      }
  }

  // -) "Loose" defines the loosest definition, "Medium" defines the tightest
  //
  if (  isMediumAcc( *lepA )   &&  isMediumAcc( *lepB )   ) {
      if ( m_useLooseAsLoosest && m_useMediumAsTightest ) { region = "TT"; }
      is_M_M_Decor( *eventInfo ) = 1;
  }
  if (  isMediumAcc( *lepA )   &&  !isMediumAcc( *lepB )  ) {
      if ( m_useLooseAsLoosest && m_useMediumAsTightest ) { region = "TL"; }
      is_M_AntiM_Decor( *eventInfo ) = 1;
      if ( OF ) {
	  if (  lepA_flavour == 11 )     { is_Mel_AntiMmu_Decor( *eventInfo ) = 1; }
	  else if ( lepA_flavour == 13 ) { is_Mmu_AntiMel_Decor( *eventInfo ) = 1; }
      }
  }
  if (  !isMediumAcc( *lepA )  &&  isMediumAcc( *lepB )   ) {
      if ( m_useLooseAsLoosest && m_useMediumAsTightest ) { region = "LT"; }
      is_AntiM_M_Decor( *eventInfo ) = 1;
      if ( OF ) {
	  if (  lepA_flavour == 11 )     { is_AntiMel_Mmu_Decor( *eventInfo ) = 1; }
	  else if ( lepA_flavour == 13 ) { is_AntiMmu_Mel_Decor( *eventInfo ) = 1; }
      }
  }
  if (  !isMediumAcc( *lepA )  &&  !isMediumAcc( *lepB )  ) {
      if ( m_useLooseAsLoosest && m_useMediumAsTightest ) { region = "LL"; }
      is_AntiM_AntiM_Decor( *eventInfo ) = 1;
  }

  // ***********************************************

  bool is2lepSS = ( nLeptons == 2 && isSS01(*eventInfo) );
  bool is3lep   = ( nLeptons == 3 && isSS12(*eventInfo) );

  if ( !( is2lepSS || is3lep ) ) {
    return EL::StatusCode::SUCCESS; // no weights in the other categories: just return
  }

  if ( is2lepSS && m_debug ) { Info("fakeWeightCalculator()", "Dilepton SS category. Region is %s ", region.c_str() ); }
  if ( is3lep && m_debug )   { Info("fakeWeightCalculator()", "Trilepton 2SS+1OS category. Region (defined by the 2SS leptons) is %s ", region.c_str() ); }

  if ( m_debug ) { Info("fakeWeightCalculator()", "Start calculating MM and FF weights... "); }

  // *******************************************
  // Now calculating MM real rate and fake rate.
  //
  // NB: NO NEED TO CALCULATE FF RATE BECAUSE FF ARE NOW OBTAINED FROM THE MM FOR r1=r2=1

  // real and fake rates w/ syst variations
  //
  std::vector<double> r1, r2, f1, f2;
  double r1up,r1dn, r2up, r2dn, f1up, f1dn, f2up, f2dn;

  bool isFakeLep(true);

  if ( lepA_flavour == 11 ) {
     r1 = calc_weights( m_el_hist_map, lepA_pt, lepA_eta, !isFakeLep, m_n_el_bins_eta, m_n_el_bins_pt_fr, m_n_el_bins_pt_rr, m_el_fr_tot, m_el_rr_tot );
     f1 = calc_weights( m_el_hist_map, lepA_pt, lepA_eta, isFakeLep,  m_n_el_bins_eta, m_n_el_bins_pt_fr, m_n_el_bins_pt_rr, m_el_fr_tot, m_el_rr_tot );
  } else if ( lepA_flavour == 13 ) {
     r1 = calc_weights( m_mu_hist_map, lepA_pt, lepA_eta, !isFakeLep, m_n_mu_bins_eta, m_n_mu_bins_pt_fr, m_n_mu_bins_pt_rr, m_mu_fr_tot, m_mu_rr_tot );
     f1 = calc_weights( m_mu_hist_map, lepA_pt, lepA_eta, isFakeLep,  m_n_mu_bins_eta, m_n_mu_bins_pt_fr, m_n_mu_bins_pt_rr, m_mu_fr_tot, m_mu_rr_tot );
  }

  if ( lepB_flavour == 11 ) {
     r2 = calc_weights( m_el_hist_map, lepB_pt, lepB_eta, !isFakeLep, m_n_el_bins_eta, m_n_el_bins_pt_fr, m_n_el_bins_pt_rr, m_el_fr_tot, m_el_rr_tot );
     f2 = calc_weights( m_el_hist_map, lepB_pt, lepB_eta, isFakeLep,  m_n_el_bins_eta, m_n_el_bins_pt_fr, m_n_el_bins_pt_rr, m_el_fr_tot, m_el_rr_tot );
  } else if ( lepB_flavour == 13 ) {
     r2 = calc_weights( m_mu_hist_map, lepB_pt, lepB_eta, !isFakeLep, m_n_mu_bins_eta, m_n_mu_bins_pt_fr, m_n_mu_bins_pt_rr, m_mu_fr_tot, m_mu_rr_tot );
     f2 = calc_weights( m_mu_hist_map, lepB_pt, lepB_eta, isFakeLep,  m_n_mu_bins_eta, m_n_mu_bins_pt_fr, m_n_mu_bins_pt_rr, m_mu_fr_tot, m_mu_rr_tot );
  }

  if ( m_debug ) {
    Info("fakeWeightCalculator()", "\n Lepton 1 \n flavour: %i \n pT = %.2f \n eta = %.2f \n ****** \n Nominal real and fake rates: \n r1 = %f , f1 = %f ", lepA_flavour, lepA_pt, lepA_eta,  r1.at(0), f1.at(0) );
    Info("fakeWeightCalculator()", "\n Lepton 2 \n flavour: %i \n pT = %.2f \n eta = %.2f \n ****** \n Nominal real and fake rates: \n r2 = %f , f2 = %f ", lepB_flavour, lepB_pt, lepB_eta,  r2.at(0), f2.at(0) );
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

  if ( !m_isMC || m_useMCForTagAndProbe ) {

    // r cannot be 0 and has always to be more than f
    // WARNING!
    // WE SHOULD BE ALSO SURE THAT REAL EFFICIENCY IS ALWAYS < 1 FOR ANY PART OF THE PHASE SPACE.
    // YOU COULD HAVE SOME R(PT)*R(ETA)>1 AND THIS CANNOT HAPPEN

    double mm_weight(0.0);

    if ( (r1.at(0) == 0) || (r2.at(0) == 0) || (r1.at(0) <= f1.at(0)) || (r2.at(0) <= f2.at(0)) ) {
    	// event will be decorated w/ null weight - will basically cancel out this event at plotting
    	//
    	if ( m_debug ) {
	  Warning("fakeWeightCalculator()", "Warning! The Matrix Method cannot be applied in event %llu for run %i because : \n r1 = %f , r2 = %f , f1 = %f , f2 = %f \n ,given that pt1 = %f , eta1 = %f , pt2 = %f , eta2 = %f ", eventInfo->eventNumber(), eventInfo->runNumber(), r1.at(0), r2.at(0),  f1.at(0), f2.at(0), lepA_pt/1e3, lepA_eta, lepB_pt/1e3, lepB_eta);
	  Warning("fakeWeightCalculator()", "applying MM weight = 0 ...");
    	}
    } else {
    	//decorate event w/ calculated MM weight
    	//
    	mm_weight      = calc_final_event_weight( region, f1.at(0), f2.at(0), r1.at(0), r2.at(0) );
    }

    // decorate with nominal MM weight
    MMWeightDecor( *eventInfo ).at(0) = mm_weight;

    if ( m_debug ) { Info("fakeWeightCalculator()", "MM final weight = %f ", mm_weight ); }

    if ( mm_weight != 0.0 ) {

      // decorate event w/ MM weight with systematics
      //
      r1up = ( r1.at(1) > 1.0 ) ? 1.0 :  r1.at(1) ;
      r2up = ( r2.at(1) > 1.0 ) ? 1.0 :  r2.at(1) ;
      r1dn = r1.at(2);
      r2dn = r2.at(2);

      f1up = f1.at(1);
      f2up = f2.at(1);
      f1dn = ( f1.at(2) < 0.0 ) ? 0.0 :  f1.at(2) ;
      f2dn = ( f2.at(2) < 0.0 ) ? 0.0 :  f2.at(2) ;

      // rup syst
      //
      MMWeightDecor( *eventInfo ).at(1) = ( calc_final_event_weight( region, f1.at(0), f2.at(0), r1up, r2up ) / mm_weight );
      // fdn syst
      //
      MMWeightDecor( *eventInfo ).at(4) = ( calc_final_event_weight( region, f1dn, f2dn, r1.at(0), r2.at(0) ) / mm_weight );

      if ( (r1dn > f1.at(0)) && (r2dn > f2.at(0)) ) {
        // rdn syst
        //
        MMWeightDecor( *eventInfo ).at(2) = ( calc_final_event_weight(region, f1.at(0), f2.at(0), r1dn, r2dn) / mm_weight );
      } else {
        if ( m_debug ) {
	   Warning("fakeWeightCalculator()", "Warning! Systematic MMWeight_rdn cannot be calculated in event %llu for run %i because : \n r1dn = %f , r2dn = %f , f1 = %f , f2 = %f \n ,given that pt1 = %f , eta1 = %f , pt2 = %f , eta2 = %f ", eventInfo->eventNumber(), eventInfo->runNumber(), r1dn, r2dn,  f1.at(0), f2.at(0), lepA_pt/1e3, lepA_eta, lepB_pt/1e3, lepB_eta);
        }
      }

      if ( (r1.at(0) > f1up) && (r2.at(0) > f2up) ) {
        // fup syst
        //
        MMWeightDecor( *eventInfo ).at(3) = ( calc_final_event_weight(region, f1up, f2up, r1.at(0),  r2.at(0)) / mm_weight );
      } else {
        if ( m_debug ) {
	   Warning("fakeWeightCalculator()", "Warning! Systematic MMWeight_fup cannot be calculated in event %llu for run %i because : \n r1 = %f , r2 = %f , f1up = %f , f2up = %f \n ,given that pt1 = %f , eta1 = %f , pt2 = %f , eta2 = %f ", eventInfo->eventNumber(), eventInfo->runNumber(), r1.at(0), r2.at(0),  f1up, f2up, lepA_pt/1e3, lepA_eta, lepB_pt/1e3, lepB_eta);
        }
      }
    }

  } // close check on m_isMC

  // *****************************************
  // The Fake Factor Method: weight the events!
  //
  // NB: it must hold: r = 1, f != 1
  //
  //

  double ff_weight(0.0);

  if ( (f1.at(0) == 1.0) || (f2.at(0) == 1.0) ) {

    if (true) {
      Warning("fakeWeightCalculator()", "Warning! The Fake Factor Method cannot be applied in event %llu for run %i because : \n f1 = %f , f2 = %f \n ,given that pt1 = %f , eta1 = %f , pt2 = %f , eta2 = %f ", eventInfo->eventNumber(), eventInfo->runNumber(), f1.at(0), f2.at(0), lepA_pt/1e3, lepA_eta, lepB_pt/1e3, lepB_eta);
      Warning("fakeWeightCalculator()", "applying FF weight = 0 ...");
    }

    // decorate event w/ (nominal) null weight
    //
    FFWeightDecor( *eventInfo ).at(0) = ff_weight;

  } else {

    //decorate event w/ calculated nominal FF weight
    //
    ff_weight		 =  calc_final_event_weight( region, f1.at(0), f2.at(0) );
    FFWeightDecor( *eventInfo ).at(0) = ff_weight;

    if ( ff_weight != 0.0 ) {

      // decorate event w/ FF weight with systematics
      //
      f1up = f1.at(1);
      f2up = f2.at(1);
      f1dn = ( f1.at(2) < 0.0 ) ? 0.0 : f1.at(2) ;
      f2dn = ( f2.at(2) < 0.0 ) ? 0.0 : f2.at(2);

      // dn syst
      //
      FFWeightDecor( *eventInfo ).at(2) = ( calc_final_event_weight( region, f1dn, f2dn ) / ff_weight );

      if ( (f1up < 1) && (f2up < 1)) {
        // up syst
        //
        FFWeightDecor( *eventInfo ).at(1) = ( calc_final_event_weight( region, f1up, f2up ) / ff_weight );
      } else {
        Warning("fakeWeightCalculator()", "Warning! Systematic FFWeight_up cannot be calculated in event %llu for run %i because : \n f1up = %f , f2up = %f \n ,given that pt1 = %f , eta1 = %f , pt2 = %f , eta2 = %f ", eventInfo->eventNumber(), eventInfo->runNumber(), f1up, f2up, lepA_pt/1e3, lepA_eta, lepB_pt/1e3, lepB_eta);
      }

    }
  }

  if ( m_debug ) { Info("fakeWeightCalculator()", "FF final weight = %f ", ff_weight ); }

  return EL::StatusCode::SUCCESS;

}

double  HTopMultilepAnalysis :: scaleRateToEfficiency( double rate )
{
  if ( rate < 0 ) { rate = 0.0; }

  double eff = ( rate / (rate+1.0) );

  return eff;
}

std::vector<double>  HTopMultilepAnalysis :: calc_weights( std::map< std::string, TH1D* > &histograms,
							   float pt,
							   float eta,
							   bool isFakeLep,
							   int n_bins_eta,
							   int n_bins_pt_fr,
							   int n_bins_pt_rr,
							   float fr_tot,
							   float rr_tot
							 )
{

  // Read the real/fake rates from input histograms
  //
  // Will eventually convert these to real/fake FACTORS

  // As a first thing, convert pT in GeV!
  //
  pt = pt/1e3;

  std::vector<double> weights(3,0.0); //initialized with zeroes

  weights.at(0) = 1.0;
  double error(0.0);

  // loop over number of eta bins
  // do not consider underflow, i.e. 0th bin
  //
  for ( int e = 1; e <= n_bins_eta; e++ ) {

    // check whether the eta under question is in *this* eta range
    //
    if ( ( fabs(eta) >= (histograms.find("eta_rr")->second)->GetXaxis()->GetBinLowEdge(e) ) && ( fabs(eta) < (histograms.find("eta_rr")->second)->GetXaxis()->GetBinLowEdge(e+1) ) ) {

      // case 1) : lepton is fake: choose correct pt histogram
      //
      if ( isFakeLep ) {

	// loop over number of pt bins
        // do not consider underflow, i.e. 0th bin
        //
        for ( int p = 1; p <= n_bins_pt_fr; p++ ) {

     	  if ( ( pt >= (histograms.find("pt_fr")->second)->GetXaxis()->GetBinLowEdge(p) ) && ( pt < (histograms.find("pt_fr")->second)->GetXaxis()->GetBinLowEdge(p+1) ) ) {

	    // combine eta and pt rates
	    //
	    double fr_pt  = (histograms.find("pt_fr")->second)->GetBinContent(p);
	    double fr_eta = (histograms.find("eta_fr")->second)->GetBinContent(e);

	    double fr_pt_err  = (histograms.find("pt_fr")->second)->GetBinError(p);
	    double fr_eta_err = (histograms.find("eta_fr")->second)->GetBinError(e);

	    // nominal
	    // (since we take the 1Dx1D rates, we use a normalisation factor at the denominator)
	    //
     	    weights.at(0) = ( fr_pt * fr_eta ) / fr_tot;

	    // (assuming  fr_pt,fr_eta are independent) this is the error on the product
	    // ( the constant factor at denominator will be put back later in the def of weight...
	    //
	    error	  = sqrt( (fr_eta*fr_pt_err)*(fr_eta*fr_pt_err) + (fr_pt*fr_eta_err)*(fr_pt*fr_eta_err) );

	    // up syst
	    //
     	    weights.at(1) = ( (fr_pt * fr_eta) + error ) / fr_tot;

     	    // down syst
     	    //
	    if ( (fr_pt * fr_eta) - error > 0 ) { weights.at(2) = ( (fr_pt * fr_eta) - error ) / fr_tot;}
	    else                                { weights.at(2) = 0.0; }

     	  }

	} // close loop on pT bins: fake lepton

      // lepton is real: choose correct pt histogram
      //
      }	else {

	// loop over number of pt bins
        // do not consider underflow, i.e. 0th bin
        //
        for ( int p = 1; p <= n_bins_pt_rr; p++ ) {

     	  if ( ( pt >= (histograms.find("pt_rr")->second)->GetXaxis()->GetBinLowEdge(p) ) && ( pt < (histograms.find("pt_rr")->second)->GetXaxis()->GetBinLowEdge(p+1) ) ) {


	    // combine eta and pt rates
	    //
	    double rr_pt  = (histograms.find("pt_rr")->second)->GetBinContent(p);
	    double rr_eta = (histograms.find("eta_rr")->second)->GetBinContent(e);

	    double rr_pt_err  = (histograms.find("pt_rr")->second)->GetBinError(p);
	    double rr_eta_err = (histograms.find("eta_rr")->second)->GetBinError(e);

	    // nominal
	    // (since we take the 1Dx1D rates, we use a normalisation factor at the denominator)
	    //
	    weights.at(0) = ( rr_pt * rr_eta ) / rr_tot;

	    // (assuming  rr_pt,rr_eta are independent) this is the error on the product
	    // ( the constant factor at denominator will be put back in the def of weight...
	    //
	    error	  = sqrt( (rr_eta*rr_pt_err)*(rr_eta*rr_pt_err) + (rr_pt*rr_eta_err)*(rr_pt*rr_eta_err) );

     	    // up syst
	    //
     	    weights.at(1) = ( (rr_pt * rr_eta) + error ) / rr_tot;

	    // down syst
	    //
     	    if ( (rr_pt * rr_eta) - error > 0 ) { weights.at(2) = ( (rr_pt * rr_eta) - error ) / rr_tot; }
	    else                                { weights.at(2) = 0.0; }

     	  }
        } // close loop on pT bins: real lepton

      } // close check isFakeLep

    } // close check on eta bin

  } // close loop on eta bins

  // Now converting rates to the factors for the MM/FF
  //
  if ( m_debug ) { Info("calc_weights()", "Rates = %f ( up = %f , dn = %f )", weights.at(0), weights.at(1), weights.at(2) ); }

  weights.at(0) = scaleRateToEfficiency(weights.at(0));
  weights.at(1) = scaleRateToEfficiency(weights.at(1));
  weights.at(2) = scaleRateToEfficiency(weights.at(2));

  if ( m_debug ) { Info("calc_weights()", "MM/FF factor = %f ( up = %f , dn = %f )", weights.at(0), weights.at(1), weights.at(2) ); }

  return weights;
}


// This function saves the final MM/FF event weight, depending on the event type (TT,TL...)
//
// The Fake Factor Method weight is obtained under the hypothesis
// r1=1 and r2=1
// So for the FF Method, you will need to pass just region, f1 and f2 (NB: make sure they are different from 0 and from 1!)
// For the MM, you will need to pass also r1 and r2
//
double HTopMultilepAnalysis :: calc_final_event_weight( std::string region, double f1, double f2, double r1, double r2 )
{

   double weight = 1.0;
   double alpha  = 1.0 / ( (r1-f1) * (r2-f2) );

   if      ( region=="TT" ) { weight = 1 - ( r1 * r2 * (1-f1) * (1-f2) * alpha ); }
   else if ( region=="TL" ) { weight = r1 * r2 * f2 * (1-f1) * alpha;  }
   else if ( region=="LT" ) { weight = r1 * r2 * f1 * (1-f2) * alpha;  }
   else if ( region=="LL" ) { weight = -1 * r1 * r2 * f1 * f2 * alpha; }

   if ( m_debug ) { Info("calc_final_event_weight()", "In region %s : \n weight = %.15f , alpha = %.15f ", region.c_str(), weight, alpha); }

   return weight;
}

// Calculate QMisID weight
//
EL::StatusCode HTopMultilepAnalysis :: QMisIDWeightCalculator (const xAOD::EventInfo* eventInfo, const xAOD::ElectronContainer* electrons )
{

  SG::AuxElement::Decorator< std::vector<float> > QMisIDWeightDecor( "QMisIDWeight" );

  // Set default weight (w/ up and dn variations)
  //
  std::vector<float> weights(3,1.0);

  if ( m_debug ) { Info("QMisIDWeightCalculator()", "QMisID initial weight = %f ( up = %f, dn = %f )", weights.at(0), weights.at(1), weights.at(2) ); }

  // QMisID is a data weight only
  //
  if ( !m_isMC ) {

    const xAOD::Electron* elA(nullptr);
    const xAOD::Electron* elB(nullptr);

    unsigned int nel = electrons->size();

    if ( nel > 0 ) { elA = electrons->at(0); }
    if ( nel > 1 ) { elB = electrons->at(1); }

    EL_RETURN_CHECK("HTopMultilepAnalysis::QMisIDWeightCalculator()", this->calc_QMisID_weights( weights, elA, elB ) );

    QMisIDWeightDecor( *eventInfo ) = weights;
  }

  if ( m_debug ) { Info("QMisIDWeightCalculator()", "QMisID final weight = %f ( up = %f, dn = %f )", QMisIDWeightDecor( *eventInfo ).at(0), QMisIDWeightDecor( *eventInfo ).at(1), QMisIDWeightDecor( *eventInfo ).at(2) ); }

  return EL::StatusCode::SUCCESS;

}

EL::StatusCode HTopMultilepAnalysis :: calc_QMisID_weights( std::vector<float>& weights, const xAOD::Electron* elA, const xAOD::Electron* elB )
{

  // If there are no electrons, return
  //
  if ( !elA && !elB ) { return EL::StatusCode::SUCCESS; }

  // accessor to tight selected electrons
  //
  static SG::AuxElement::Accessor< char > isTightAcc("isTight");

  float elA_eta = elA->caloCluster()->etaBE(2);
  float elA_pt  = elA->pt()/1e3;
  bool  elA_isT = isTightAcc( *elA );

  float elB_eta = ( elB ) ? elB->caloCluster()->etaBE(2) : -999.0;
  float elB_pt  = ( elB ) ? elB->pt()/1e3 : -1.0;
  bool  elB_isT = ( elB ) ? isTightAcc( *elB ) : false;

  float rA(0.0), rA_up(0.0), rA_dn(0.0), rB(0.0), rB_up(0.0), rB_dn(0.0);

  // Get the 2D histogram from the map
  //
  TH2D* twoD_rates_T     = ( m_QMisID_hist_map.find("T")->second );
  TH2D* twoD_rates_AntiT = ( m_QMisID_hist_map.find("AntiT")->second );

  // Make X and Y projections of the twoD histogram with the rates
  //
  TH1D* proj_eta_T     = twoD_rates_T->ProjectionX("proj_eta_T");
  TH1D* proj_pt_T      = twoD_rates_T->ProjectionY("proj_pt_T");
  TH1D* proj_eta_AntiT = twoD_rates_AntiT->ProjectionX("proj_eta_AntiT");
  TH1D* proj_pt_AntiT  = twoD_rates_AntiT->ProjectionY("proj_pt_AntiT");

  // Look at elA first...
  //
  if ( elA_isT ) {
    EL_RETURN_CHECK("HTopMultilepAnalysis::calc_QMisID_weights()", this->readRatesAndError(twoD_rates_T, proj_eta_T, proj_pt_T, elA_eta, elA_pt, rA, rA_up, rA_dn) );
  } else {
    EL_RETURN_CHECK("HTopMultilepAnalysis::calc_QMisID_weights()", this->readRatesAndError(twoD_rates_AntiT, proj_eta_AntiT, proj_pt_AntiT, elA_eta, elA_pt, rA, rA_up, rA_dn) );
  }

  // .. and now at elB (if any...otherwise rB weights will be zero by default)
  //
  if ( elB ) {
    if ( elB_isT ) {
      EL_RETURN_CHECK("HTopMultilepAnalysis::calc_QMisID_weights()", this->readRatesAndError(twoD_rates_T, proj_eta_T, proj_pt_T, elB_eta, elB_pt, rB, rB_up, rB_dn) );
    } else {
      EL_RETURN_CHECK("HTopMultilepAnalysis::calc_QMisID_weights()", this->readRatesAndError(twoD_rates_AntiT, proj_eta_AntiT, proj_pt_AntiT, elB_eta, elB_pt, rB, rB_up, rB_dn) );
    }
  }

  // Finally, store the event weight + variations
  //
  weights.at(0) = ( rA + rB - 2.0 * rA * rB ) / ( 1.0 - rA - rB + 2.0 * rA * rB ) ;
  weights.at(1) = ( rA_up + rB_up - 2.0 * rA_up * rB_up ) / ( 1.0 - rA_up - rB_up + 2.0 * rA_up * rB_up );
  weights.at(2) = ( rA_dn + rB_dn - 2.0 * rA_dn * rB_dn ) / ( 1.0 - rA_dn - rB_dn + 2.0 * rA_dn * rB_dn );

  return EL::StatusCode::SUCCESS;

}

EL::StatusCode HTopMultilepAnalysis :: readRatesAndError(TH2D* rate_map, TH1D* proj_X, TH1D* proj_Y,
                                                         const float& x, const float& y,
							 float& r, float& r_up, float& r_dn )
{

  float this_low_edge(-999.0),this_up_edge(-999.0);

  int xbin_nr(-1), ybin_nr(-1);

  // Loop over the projections, and keep track of the bin number where (x,y) is found
  //
  for ( int xbin = 0; xbin < proj_X->GetNbinsX()+1; ++xbin  ) {

    this_low_edge = proj_X->GetXaxis()->GetBinLowEdge(xbin);
    this_up_edge  = proj_X->GetXaxis()->GetBinLowEdge(xbin+1);

    if ( fabs(x) >= this_low_edge && fabs(x) < this_up_edge ) {
      xbin_nr = proj_X->GetBin(xbin);
      break;
    }

  }
  for ( int ybin = 0; ybin < proj_Y->GetNbinsX()+1; ++ ybin ) {

    this_low_edge = proj_Y->GetXaxis()->GetBinLowEdge(ybin);
    this_up_edge  = proj_Y->GetXaxis()->GetBinLowEdge(ybin+1);

    if ( y >= this_low_edge && y < this_up_edge ) {
      ybin_nr = proj_Y->GetBin(ybin);
      break;
    }

  }

  // Now get the NOMINAL rate via global bin number (x,y)

  r = rate_map->GetBinContent( rate_map->GetBin( xbin_nr, ybin_nr ) );

  // Get the UP and DOWN variations
  //
  // QUESTION: Why the hell ROOT has GetBinErrorUp and GetBinErrorLow for TH2 ??
  // They seem to give always the same result...
  //
  r_up = r + rate_map->GetBinErrorUp( rate_map->GetBin( xbin_nr, ybin_nr ) );
  r_dn = r - rate_map->GetBinErrorUp( rate_map->GetBin( xbin_nr, ybin_nr ) );
  r_dn = ( r_dn > 0.0 ) ? r_dn : 0.0;

  return EL::StatusCode::SUCCESS;

}

// Calculate lepton trigger SF for the event
//
// The CP tools do not calculate the final event SF. Rather, they give back the SF and MC efficiency *per lepton*.
// With such ingredients in hand, the HTop way to compute the event SF is the following:
//
// eventSF = ( 1 - prod( 1 - SF(i)*eff(i) ) ) / ( 1 - prod ( 1 - eff(i) ) );
//
// where the productory is over the selected leptons in the event.
// The trick at the numerator is just to get the efficiency in data.
//
// The SF systematics are obtained by coherently varying the SF for each object (i.e. assume full correlation).
// The MC efficiency is assumed to have negligible uncertainty.

EL::StatusCode HTopMultilepAnalysis :: computeEventLepTrigSF( const xAOD::EventInfo* eventInfo,
                                                              const xAOD::IParticleContainer& leptons
		                                             )
{
  // Initialise product of SFs
  // (use a large size just to make sure...)
  //
  std::vector<float> numerator(10,1.0), denominator(10,1.0);

  for ( auto lep_itr : leptons ) {

    if ( m_debug ) { Info("computeEventLepTrigSF()", "lepton pT = %f \t lepton flavour = %i ", lep_itr->pt()/1e3, lep_itr->type()); }

    // check whether this lepton is T,L
    //
    bool lep_T = ( lep_itr->auxdecor< char >( "isTight" ) );

    std::string decor_name_SF(""), decor_name_eff("");

    std::vector < float > decor_SF, decor_eff;

    if ( lep_itr->type() == xAOD::Type::Electron ) {

      // hard-code this for now...FIXME
      //
      if ( lep_T ) {
        decor_name_SF  =  "ElectronEfficiencyCorrector_TrigSyst_LHLooseAndBLayer";
        decor_name_eff =  "ElectronEfficiencyCorrector_TrigMCEffSyst_LHLooseAndBLayer";
      } else {
        decor_name_SF  =  "ElectronEfficiencyCorrector_TrigSyst_LHTight";
        decor_name_eff =  "ElectronEfficiencyCorrector_TrigMCEffSyst_LHTight";
      }

    } else if ( lep_itr->type()  == xAOD::Type::Muon ) {

      if ( lep_T ) {
        decor_name_SF  = "MuonEfficiencyCorrector_TrigSyst_RecoLoose_IsoLoose";
        decor_name_eff = "MuonEfficiencyCorrector_TrigMCEff_RecoLoose_IsoLoose";
      } else {
        decor_name_SF  = "MuonEfficiencyCorrector_TrigSyst_RecoLoose_IsoFixedCutTightTrackOnly";
        decor_name_eff = "MuonEfficiencyCorrector_TrigMCEff_RecoLoose_IsoFixedCutTightTrackOnly";
      }

    }

    if ( !lep_itr->isAvailable< std::vector< float > >( decor_name_SF ) )  { Error("computeEventLepTrigSF()", "trigger SF not available for this lepton. Aborting."); return EL::StatusCode::FAILURE; }
    if ( !lep_itr->isAvailable< std::vector< float > >( decor_name_eff ) ) { Error("computeEventLepTrigSF()", "trigger MC eff. not available for this lepton. Aborting."); return EL::StatusCode::FAILURE; }

    decor_SF  = lep_itr->auxdecor< std::vector< float > >( decor_name_SF );
    decor_eff = lep_itr->auxdecor< std::vector< float > >( decor_name_eff );

    // check SF and efficiency vectors have same size
    //
    if ( decor_SF.size() != decor_eff.size() ) { Error("computeEventLepTrigSF()","trigger SF vector and trigger MC eff. vector have different size. Aborting."); return EL::StatusCode::FAILURE; }

    // consider syst variations only for SF (keep nominal eff.)
    //
    auto effSF      = decor_SF.begin();
    auto effNominal = decor_eff.at(0);

    // ---------- NUMERATOR ----------
    //

    std::vector<float> num_factor;

    int idx(0);
    for ( ; effSF != decor_SF.end(); ++effSF, ++idx ) {

      num_factor.push_back( 1.0 - *effSF  * effNominal );

      if ( m_debug ) {
        std::cout << "\t Numerator: " << std::endl;
        std::cout << "\t efficiency SF[" << idx << "] = " << *effSF	 << std::endl;
        std::cout << "\t efficiency[" << idx << "] = "    <<  effNominal << std::endl;
        std::cout << "\t num_factor[" << idx << "] = "    <<   1.0 - *effSF  * effNominal << std::endl;
      }

    }
    // this makes sure numerator has got as much elements as num_factor
    //
    numerator.resize(num_factor.size());

    // this updates the value of numerator to contain the element-wise product
    // of the previous numerator and the num_factor for *this* lepton
    //
    std::transform( numerator.begin(), numerator.end(), num_factor.begin(), numerator.begin(), std::multiplies<float>() );

    // ---------- DENOMINATOR ----------
    //

    // just use the nominal MC efficiency
    //
    std::vector<float> denom_factor( decor_eff.size(), 1.0 - effNominal );

    if ( m_debug ) {
      for ( auto eff = denom_factor.begin(); eff != denom_factor.end(); ++eff, ++idx ) {
        std::cout << "\t Denominator: " << std::endl;
        std::cout << "\t efficiency[" << idx << "] = "    <<  effNominal << std::endl;
        std::cout << "\t denom_factor[" << idx << "] = "  <<   1.0 - effNominal << std::endl;
      }
    }

    // this makes sure denominator has got as much elements as denom_factor
    //
    denominator.resize(denom_factor.size());

    // check numerator and denominator have same size
    //
    if ( numerator.size() != denominator.size() ) { Error("computeEventLepTrigSF()","Numerator and denominator in SF formula must have same size. Aborting."); return EL::StatusCode::FAILURE; }

    // this updates the value of denominator to contain the element-wise product
    // of the previous denominator and the denom_factor for *this* electron
    //
    std::transform( denominator.begin(), denominator.end(), denom_factor.begin(), denominator.begin(), std::multiplies<float>() );

  }

  // update numerator and denominator
  // Make sure the SF in the 0/0 case (i.e, when efficiency=0) will be set equal to 1
  //
  for ( auto &it : numerator )   { it = ( it != 1.0 ) ? (1.0 - it) : 1.0; }
  for ( auto &it : denominator ) { it = ( it != 1.0 ) ? (1.0 - it) : 1.0; }

  std::vector<float> lepTrigSF(numerator.size());

  // ...and finally get the global SF w/ all systematics!
  //
  std::transform( numerator.begin(), numerator.end(), denominator.begin(), lepTrigSF.begin(), std::divides<float>() );

  // decorator for event trigger SF
  //
  SG::AuxElement::Decorator< std::vector<float> > lepTrigSFDecor_GLOBAL("lepTrigEffSF_GLOBAL_HTop");
  lepTrigSFDecor_GLOBAL( *eventInfo ) = lepTrigSF;

  if ( m_debug ) {
    int idx(0);
    for ( const auto &effSF : lepTrigSF ) {
      Info( "computeEventLepTrigSF()", "===>>>");
      Info( "computeEventLepTrigSF()", " ");
      Info( "computeEventLepTrigSF()", "Trigger efficiency SF[%i] = %f", idx, effSF);
      ++idx;
    }
    Info( "computeEventLepTrigSF()", "--------------------------------------");

  }

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode HTopMultilepAnalysis :: computeEventLepSF( const xAOD::EventInfo* eventInfo,
                                                          const xAOD::IParticleContainer& leptons,
							  SFType TYPE
		                                        )
{
  // Initialise global per-event SF ( will be the product of each object SF)
  // (use a large size just to make sure...)
  //
  std::vector<float> SF(10,1.0);

  std::string type;

  for ( auto lep_itr : leptons ) {

    // in some cases we want the SF vector to be the default one (i.e, all 1.0's...)
    bool skip(false);

    std::string prefix(""), append_L(""), append_T("");
    bool lep_passWP_T(false); // for Isolation/ID

    bool isEl = ( lep_itr->type() == xAOD::Type::Electron );
    bool isMu = ( lep_itr->type() == xAOD::Type::Muon );

    if ( isEl )      { prefix = "ElectronEfficiencyCorrector"; }
    else if ( isMu ) { prefix = "MuonEfficiencyCorrector"; }

    // check whether this lepton is T,L
    //
    bool lep_T = ( lep_itr->auxdecor< char >( "isTight" ) );

    switch ( TYPE )
      {
      case SFType::RECO:
	type = "Reco";
	if ( isEl ) {
	  append_L = append_T = "RecoSyst";
	} else if ( isMu ) {
	  append_L = append_T = "RecoSyst_Loose";
	}
	break;
      case SFType::ISOLATION:
	type = "Iso";
	if ( isEl ) {
	  append_L = "IsoSyst_IsoLoose";
	  append_T = "IsoSyst_IsoFixedCutTight";
	  lep_passWP_T = ( lep_itr->auxdecor< char >( "isIsolated_FixedCutTight" ) );
	} else if ( isMu ) {
	  append_L = "IsoSyst_IsoLoose";
	  append_T = "IsoSyst_IsoFixedCutTightTrackOnly";
	  lep_passWP_T = ( lep_itr->auxdecor< char >( "isIsolated_FixedCutTightTrackOnly" ) );
	}
	break;
      case SFType::ID:
	type = "ID";
	if ( isEl ) {
	  append_L = "PIDSyst_LHLooseAndBLayer";
	  append_T = "PIDSyst_LHTight";
	  lep_passWP_T  = ( lep_itr->auxdecor< char >( "LHTight" ) );
	} else if ( isMu ) {
	  skip = true;
	}
	break;
      case SFType::TTVA:
	type = "TTVA";
	if ( isEl ) {
	  skip = true;
	} else if ( isMu ) {
	  // TTVA SF is to be applied only to tight muons!
	  if ( !lep_T ) { skip = true; }
	  append_T = "TTVASyst_TTVA";
	}
	break;
      default:
	Error ("computeEventLepSF()", "Unsupported SF type. Aborting." );
	return EL::StatusCode::FAILURE;
	break;
      }

    if ( m_debug ) {
      Info("computeEventLepSF()", "--------------------------------------");
      Info("computeEventLepSF()", "SF type = %s ", type.c_str() );
      Info("computeEventLepSF()", "lepton pT = %f \t lepton flavour = %i ", lep_itr->pt()/1e3, lep_itr->type());
    }

    if ( skip ) { continue; }

    std::string decor_name_SF("");

    std::vector < float > thisLepSF;

    // NB: a loose lepton could still have tight iso/ tight ID since we use ( iso & (ID) & IP) to discriminate T/L...
    //
    if ( lep_T || lep_passWP_T ) { decor_name_SF  = prefix + "_" + append_T; }
    else                         { decor_name_SF  = prefix + "_" + append_L; }

    if ( m_debug ) { Info("computeEventLepSF()", "Reading SF decoration name = %s ", decor_name_SF.c_str() ); }

    if ( !lep_itr->isAvailable< std::vector< float > >( decor_name_SF ) )  { Error("computeEventLepSF", "SF %s not available for this lepton. Aborting.", decor_name_SF.c_str() ); return EL::StatusCode::FAILURE; }

    thisLepSF  = lep_itr->auxdecor< std::vector< float > >( decor_name_SF );

    if ( m_debug ) {
      unsigned int idx(0);
      for ( const auto &effSF : thisLepSF ) {
	std::cout << "\t efficiency SF[" << idx << "] = " << effSF  << std::endl;
	 ++idx;
      }
    }

    // this makes sure the global SF has got as much elements as this per-object SF
    //
    SF.resize(thisLepSF.size());

    // this updates the value of the global SF to contain the element-wise product
    // of the previous global SF and the SF for *this* lepton
    //
    std::transform( SF.begin(), SF.end(), thisLepSF.begin(), SF.begin(), std::multiplies<float>() );

  }

  // decorator for event SF
  //
  const std::string decor_name = "lep" + type + "EffSF_GLOBAL_HTop";
  SG::AuxElement::Decorator< std::vector<float> > lepSFDecor_GLOBAL(decor_name);
  lepSFDecor_GLOBAL( *eventInfo ) = SF;

  if ( m_debug ) {
    unsigned int idx(0);
    for ( const auto &effSF : SF ) {
      Info( "computeEventLepSF()", "===>>>");
      Info( "computeEventLepSF()", " ");
      Info( "computeEventLepSF()", "%s SF[%i] = %f", type.c_str(), idx, effSF);
      ++idx;
    }
    Info( "computeEventLepSF()", "--------------------------------------");

  }

  return EL::StatusCode::SUCCESS;
}
