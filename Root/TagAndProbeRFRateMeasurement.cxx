/****************************************************************************************
 *
 * A fairly simple tag-and-probe selection to allow measuring lepton fake/real rates used 
 * later in the  Matrix Method / Fake Rate Method
 * 
 * The selection is aimed to target ttbar events enriched w/ real or fake leptons
 * 
 * 1 ) Select events w/ ( 2 leptons ) && ( 0 < jets < 4 ) && ( bjets > 0 ) 
       ( this is done by HTopMultilepEventSelector algo upstream. It also storesa CDv for the leptons )
 *
 * 2) 
 *
 * OS region : enriched with real leptons (targets dileptonic ttbar, some Z+jets)
 * SS region : enriched with fake leptons (targets semileptonic ttbar with one 
 * 					   non-prompt lepton, and W+jets with a 
 *					   fake lepton from one of the jets)
 * 
 * 
 * In OS region and SS region, we have different definitions of "tag" and "probe" lepton
 * See within the code
 *
 * Once tag and probe leptons are identified, the rates will be defined as:
 *
 * -) r = ( N probe T ) / ( N probe !T )  requiring OS
 * -) f = ( N probe T ) / ( N probe !T )  requiring SS
 *
 * "Tight" defintion:
 *  -) electrons: isolation (track and calo based), PID Tight/VeryTight (LH or cut-based, it's configurable)
 *  -) muons:     isolation (track and calo based), d0sig
 *
 *
 * These ratios will be calculated at plotting stage, using pT and eta distributions 
 *
 * M. Milesi (marco.milesi@cern.ch), F. Nuti (francesco.nuti@cern.ch)
 *
 ****************************************************************************************/

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
#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"
#include "xAODBase/IParticleContainer.h"
#include "xAODBase/IParticle.h"
#include "xAODBase/ObjectType.h"
#include "AthContainers/ConstDataVector.h"
#include "AthContainers/DataVector.h"

// package include(s):
#include "HTopMultilepAnalysis/TagAndProbeRFRateMeasurement.h"
#include "xAODAnaHelpers/HelperClasses.h"
#include "xAODAnaHelpers/HelperFunctions.h"
#include <xAODAnaHelpers/tools/ReturnCheck.h>
#include <xAODAnaHelpers/tools/ReturnCheckConfig.h>

// external tools include(s):

// ROOT include(s):
#include "TEnv.h"
#include "TFile.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TRandom3.h"

// external tools include(s):
#include "ElectronIsolationSelection/ElectronIsolationSelectionTool.h"
#include "ElectronPhotonSelectorTools/AsgElectronLikelihoodTool.h"
#include "ElectronPhotonSelectorTools/AsgElectronIsEMSelector.h"
#include "ElectronPhotonSelectorTools/egammaPIDdefs.h"

// this is needed to distribute the algorithm to the workers
ClassImp(TagAndProbeRFRateMeasurement)


TagAndProbeRFRateMeasurement :: TagAndProbeRFRateMeasurement () {
}

TagAndProbeRFRateMeasurement :: TagAndProbeRFRateMeasurement (std::string name, std::string configName) :
  Algorithm(),
  m_name(name),
  m_configName(configName),
  m_cutflowHist(nullptr),
  m_cutflowHistW(nullptr)
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().

  Info("TagAndProbeRFRateMeasurement()", "Calling constructor \n");
}

TagAndProbeRFRateMeasurement::~TagAndProbeRFRateMeasurement() {}

EL::StatusCode  TagAndProbeRFRateMeasurement :: configure ()
{
  Info("configure()", "Configuing TagAndProbeRFRateMeasurement Interface. User configuration read from : %s \n", m_configName.c_str());

  m_configName = gSystem->ExpandPathName( m_configName.c_str() );
  RETURN_CHECK_CONFIG( "TagAndProbeRFRateMeasurement::configure()", m_configName);

  TEnv* config = new TEnv(m_configName.c_str());

  // read debug flag from .config file
  m_debug         = config->GetValue("Debug" ,      false );
  m_useCutFlow    = config->GetValue("UseCutFlow",  true);

  // input container to be read from TEvent or TStore
  m_inContainerName_el         = config->GetValue("InputContainerElectrons",  "");
  m_inContainerName_mu         = config->GetValue("InputContainerMuons",  "");
  m_inContainerName_jet        = config->GetValue("InputContainerJets",  "");
  
  // configurable cuts
  m_doLHPIDCut                 = config->GetValue("DoLHPIDCut"         , true     );
  m_doCutBasedPIDCut           = config->GetValue("DoCutBasedPIDCut"   , false    );
  
  // isolation stuff
  m_useRelativeIso_Mu          = config->GetValue("UseRelativeIso_Mu"   , true      );
  m_CaloBasedIsoType_Mu        = config->GetValue("CaloBasedIsoType_Mu" , "topoetcone20");
  m_CaloBasedIsoCut_Mu         = config->GetValue("CaloBasedIsoCut_Mu"  , 0.05      );
  m_TrackBasedIsoType_Mu       = config->GetValue("TrackBasedIsoType_Mu", "ptcone20");
  m_TrackBasedIsoCut_Mu        = config->GetValue("TrackBasedIsoCut_Mu" , 0.05      );

  m_useRelativeIso_El          = config->GetValue("UseRelativeIso_El"   , true      );
  m_CaloBasedIsoType_El        = config->GetValue("CaloBasedIsoType_El" , "topoetcone20");
  m_CaloBasedIsoCut_El         = config->GetValue("CaloBasedIsoCut_El"  , 0.05      );
  m_TrackBasedIsoType_El       = config->GetValue("TrackBasedIsoType_El", "ptcone20");
  m_TrackBasedIsoCut_El        = config->GetValue("TrackBasedIsoCut_El" , 0.05      );

  // parse and split by comma
  std::string token;

  m_passAuxDecorKeys        = config->GetValue("PassDecorKeys", "");
  std::istringstream ss(m_passAuxDecorKeys);
  while ( std::getline(ss, token, ',') ) {
    m_passKeys.push_back(token);
  }

  m_failAuxDecorKeys        = config->GetValue("FailDecorKeys", "");
  ss.str(m_failAuxDecorKeys);
  while ( std::getline(ss, token, ',') ) {
    m_failKeys.push_back(token);
  }

  if ( m_inContainerName_el.empty() ) {
    Error("configure()", "InputContainerElectrons is empty!");
    return EL::StatusCode::FAILURE;
  }
  if ( m_inContainerName_mu.empty() ) {
    Error("configure()", "InputContainerMuons is empty!");
    return EL::StatusCode::FAILURE;
  }

  config->Print();
  Info("configure()", "TagAndProbeRFRateMeasurement Interface succesfully configured! \n");

  delete config; config = nullptr;
  
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode TagAndProbeRFRateMeasurement :: setupJob (EL::Job& job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.

  Info("setupJob()", "Calling setupJob \n");

  job.useXAOD ();
  xAOD::Init( "TagAndProbeRFRateMeasurement" ).ignore(); // call before opening first file

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TagAndProbeRFRateMeasurement :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.

  Info("histInitialize()", "Calling histInitialize \n");
  
  if ( m_useCutFlow ) {
    TFile *file = wk()->getOutputFile("cutflow");
    m_cutflowHist  = (TH1D*)file->Get("cutflow");
    m_cutflowHistW = (TH1D*)file->Get("cutflow_weighted");
    m_cutflow_bin  = m_cutflowHist->GetXaxis()->FindBin(m_name.c_str());
    m_cutflowHistW->GetXaxis()->FindBin(m_name.c_str());
  }

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TagAndProbeRFRateMeasurement :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed

  Info("fileExecute()", "Calling fileExecute \n");


  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TagAndProbeRFRateMeasurement :: changeInput (bool /*firstFile*/)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.

  Info("changeInput()", "Calling changeInput \n");

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode TagAndProbeRFRateMeasurement :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.

  Info("initialize()", "Initializing TagAndProbeRFRateMeasurement Interface... \n");

  m_event = wk()->xaodEvent();
  m_store = wk()->xaodStore();

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
  
  // initialise ElectronIsolationSelectionTool
  m_electronIsolationSelectionTool = new CP::ElectronIsolationSelectionTool( "ElectronIsolationSelectionTool" );
  m_electronIsolationSelectionTool->msg().setLevel( MSG::VERBOSE); // ERROR, VERBOSE, DEBUG, INFO
  // https://twiki.cern.ch/twiki/bin/view/AtlasProtected/ElectronIsolationSelectionTool
  HelperClasses::EnumParser<xAOD::Iso::IsolationType> isoParser;
  m_electronIsolationSelectionTool->configureCutBasedIsolation( isoParser.parseEnum(m_CaloBasedIsoType_El),   static_cast<double>(m_CaloBasedIsoCut_El),  m_useRelativeIso_El );
  m_electronIsolationSelectionTool->configureCutBasedIsolation( isoParser.parseEnum(m_TrackBasedIsoType_El),  static_cast<double>(m_TrackBasedIsoCut_El), m_useRelativeIso_El );
  RETURN_CHECK( "TagAndProbeRFRateMeasurement::initialize()", m_electronIsolationSelectionTool->initialize(), "Failed to properly initialize ElectronIsolationSelectionTool." );

  // needed by tools below
  std::string confDir = "ElectronPhotonSelectorTools/offline/mc15_20150224/";

  // initialise AsgElectronLikelihoodTool (likelihood-based PID)
  m_LHToolTight2012     = new AsgElectronLikelihoodTool ("LHToolTight2012");
  m_LHToolVeryTight2012 = new AsgElectronLikelihoodTool ("LHToolVeryTight2012");
  // initialize the primary vertex container for the tool to have access to the number of vertices used to adapt cuts based on the pileup
  m_LHToolTight2012     ->setProperty("primaryVertexContainer","PrimaryVertices");
  m_LHToolVeryTight2012 ->setProperty("primaryVertexContainer","PrimaryVertices");
  m_LHToolTight2012     ->setProperty("ConfigFile",confDir+"ElectronLikelihoodTightOfflineConfig2012.conf");
  m_LHToolVeryTight2012 ->setProperty("ConfigFile",confDir+"ElectronLikelihoodVeryTightOfflineConfig2012.conf");
  RETURN_CHECK( "TagAndProbeRFRateMeasurement::initialize()", m_LHToolTight2012->initialize(), "Failed to properly initialize AsgElectronLikelihoodTool w/ Tight WP." );
  RETURN_CHECK( "TagAndProbeRFRateMeasurement::initialize()", m_LHToolVeryTight2012->initialize(), "Failed to properly initialize AsgElectronLikelihoodTool w/ VeryTight WP." );

  // initialise AsgElectronIsEMSelector (cut-based PID)
  m_TightPP2012  = new AsgElectronIsEMSelector("TightPP2012");
  m_TightPP2012->setProperty("ConfigFile", confDir + "ElectronIsEMTightSelectorCutDefs2012.conf" );
  // bitmask to be set only for DC14 w/ 2012 configuration
  m_TightPP2012->setProperty("isEMMask", static_cast<unsigned int>(egammaPID::ElectronTightPP) );    
  RETURN_CHECK( "ElectronSelector::initialize()", m_TightPP2012->initialize(), "Failed to properly initialize AsgElectronIsEMSelector w/ Tight WP." );

  Info("initialize()", "TagAndProbeRFRateMeasurement Interface succesfully initialized!" );

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TagAndProbeRFRateMeasurement :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  if ( m_debug ) { Info("execute()", "Applying TagAndProbeRFRateMeasurement MM Rate Measurement... \n"); }

  // retrieve event 
  const xAOD::EventInfo* eventInfo(nullptr);
  RETURN_CHECK("TagAndProbeRFRateMeasurement::execute()", HelperFunctions::retrieve(eventInfo, "EventInfo", m_event, m_store, m_debug) , "");

  bool isMC = ( eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) );

  // MC event weight 
  float mcEvtWeight(1.0);
  static SG::AuxElement::Accessor< float > mcEvtWeightAcc("mcEventWeight");
  if ( ! mcEvtWeightAcc.isAvailable( *eventInfo ) ) {
    Error("execute()  ", "mcEventWeight is not available as decoration! Aborting" );
    return EL::StatusCode::FAILURE;
  }
  mcEvtWeight = mcEvtWeightAcc( *eventInfo );
  
  m_numEvent++;

  // this will be the collection processed - no matter what!!
  const xAOD::ElectronContainer* inElectrons(nullptr); 
  RETURN_CHECK("TagAndProbeRFRateMeasurement::execute()", HelperFunctions::retrieve(inElectrons, m_inContainerName_el, m_event, m_store, m_debug) , "");
  const xAOD::MuonContainer* inMuons(nullptr);  
  RETURN_CHECK("TagAndProbeRFRateMeasurement::execute()", HelperFunctions::retrieve(inMuons, m_inContainerName_mu, m_event, m_store, m_debug) , "");
  const xAOD::JetContainer* inJets(nullptr);
  RETURN_CHECK("TagAndProbeRFRateMeasurement::execute()", HelperFunctions::retrieve(inJets, m_inContainerName_jet, m_event, m_store, m_debug) , "");

  ConstDataVector<xAOD::IParticleContainer>* leptonsCDV(nullptr);
  RETURN_CHECK("HTopMultilepAnalysis::execute()", HelperFunctions::retrieve(leptonsCDV, "Leptons_Sel", m_event, m_store, m_debug) ,"");
  // Make a sorted version of the container 
  // (this can be on the stack! Will not be pushed to the store...)
  const xAOD::IParticleContainer leptonsSorted = HelperFunctions::sort_container_pt( leptonsCDV->asDataVector() );
  
  unsigned int nSignalLeptons = leptonsCDV->size();
  
  if ( m_debug ){ Info("execute()"," inElectrons: %lu , inMuons: %lu , inLeptons: %u ", inElectrons->size(), inMuons->size(), nSignalLeptons ); }

  // 2.
  // -------------------------------------
  // Categorise based on lepton flavour
  // -------------------------------------
  
  int prod_lep_charge(1);
  unsigned int count_el(0), count_mu(0);
  for ( auto lep_it : leptonsSorted ) {
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
  
  // declare some event decorations
  static SG::AuxElement::Decorator< char > isDileptonOSDecor("isDileptonOS");
  static SG::AuxElement::Decorator< char > isDileptonSSDecor("isDileptonSS");
  // decorate event
  isDileptonOSDecor( *eventInfo ) = isOS;
  isDileptonSSDecor( *eventInfo ) = isSS;

  // 3.
  // --------------------------------------------------
  // Check if leading & subleading pass *tight* selection
  //
  //  Electrons : ( (relative isolation (track & calo based, DeltaR < 0.2) < 0.05)  && ( LH/cut-based PID = tight ) )
  //  Muons:      ( (relative isolation (track & calo based, DeltaR < 0.2) < 0.05)  && ( d0sig < 3.0) ) 
  //
  // --------------------------------------------------
  
  const xAOD::IParticle* leadingLepton           = *(leptonsSorted.begin());
  xAOD::Type::ObjectType leadingLeptonFlavour    = leadingLepton->type();
  float                  leadingLeptonPt         = leadingLepton->pt();
  bool                   isLeadingTight(false);   
  bool 			 isLeadingTrigMatched(false);  
  const xAOD::IParticle* subLeadingLepton        = *(std::next(leptonsSorted.begin(),1));
  xAOD::Type::ObjectType subLeadingLeptonFlavour = subLeadingLepton->type();
  float                  subLeadingLeptonPt      = subLeadingLepton->pt();
  bool                   isSubLeadingTight(false); 
  bool 			 isSubLeadingTrigMatched(false);  
    
  // declare some lepton decorations
  static SG::AuxElement::Decorator< char > isTightDecor("isTight");
  static SG::AuxElement::Decorator< char > isTrigMatchedDecor("isTrigMatched");
  static SG::AuxElement::Decorator< char > isTagDecor("isTag");
  // declare an event decoration in case there are no "good" leptons for the rate calculation
  static SG::AuxElement::Decorator< char > isNonTightEventDecor("isNonTightEvent");
  // declare event decorations for checking the probe type in event
  static SG::AuxElement::Decorator< char > isProbeElEventDecor( "isProbeElEvent" ); 
  static SG::AuxElement::Decorator< char > isProbeMuEventDecor( "isProbeMuEvent" ); 
  
  // decorate with default values
  if ( m_debug ) { Info("execute()","Default lepton decorations:"); }
  for ( auto lep_it : leptonsSorted ) {
     isTightDecor( *lep_it )       = 0;
     isTrigMatchedDecor( *lep_it ) = 0;
     isTagDecor( *lep_it )         = 0;
     if ( m_debug ) { Info("execute()","isTight: %i \t isTrigMatched: %i \t isTag: %i ", lep_it->auxdata<char>("isTight"), lep_it->auxdata<char>("isTrigMatched"), lep_it->auxdata<char>("isTag") ); }
  }
  isNonTightEventDecor( *eventInfo ) = 0;
  isProbeElEventDecor( *eventInfo )  = 0;
  isProbeMuEventDecor( *eventInfo )  = 0;
    
  if ( m_debug ) { Info("execute()","Leading lepton: isTight: %i \t isTrigMatched: %i \t isTag: %i ", leadingLepton->auxdata<char>("isTight"),  leadingLepton->auxdata<char>("isTrigMatched"), leadingLepton->auxdata<char>("isTag") ); }
  if ( m_debug ) { Info("execute()","Subleading lepton: isTight: %i \t isTrigMatched: %i \t isTag: %i ", subLeadingLepton->auxdata<char>("isTight"),  subLeadingLepton->auxdata<char>("isTrigMatched"), subLeadingLepton->auxdata<char>("isTag") ); }
    
  HelperClasses::EnumParser<xAOD::Iso::IsolationType> isoParser;
  
  if ( isElEl ) {
    
    if ( m_debug ) { Info("execute()","ElEl event"); }

    if  ( m_doLHPIDCut ) { 
      isLeadingTight                     = ( m_LHToolTight2012->accept( *leadingLepton ) && m_electronIsolationSelectionTool->accept( *( dynamic_cast<const xAOD::Electron*>( leadingLepton ) ) ) );  
      isTightDecor( *leadingLepton  )    = isLeadingTight;
      isSubLeadingTight                  = ( m_LHToolTight2012->accept( *subLeadingLepton ) && m_electronIsolationSelectionTool->accept( *( dynamic_cast<const xAOD::Electron*>( subLeadingLepton ) ) ) ); 
      isTightDecor( *subLeadingLepton  ) = isSubLeadingTight;
    } else if ( m_doCutBasedPIDCut ) { 
      isLeadingTight                     = ( m_TightPP2012->accept( *leadingLepton ) && m_electronIsolationSelectionTool->accept( *( dynamic_cast<const xAOD::Electron*>( leadingLepton ) ) )  ) ; 
      isTightDecor( *leadingLepton  )    = isLeadingTight;
      isSubLeadingTight                  = ( m_TightPP2012->accept( *subLeadingLepton ) && m_electronIsolationSelectionTool->accept( *( dynamic_cast<const xAOD::Electron*>( subLeadingLepton ) ) ) ); 
      isTightDecor( *subLeadingLepton  ) = isSubLeadingTight;
    }

  } else if ( isMuMu ) { 

      if ( m_debug ) { Info("execute()","MuMu event"); }

      const xAOD::TrackParticle* tpLeading     = dynamic_cast<const xAOD::Muon*>(leadingLepton)->primaryTrackParticle();
      float d0SignificanceLeading              = fabs( tpLeading->d0() ) / sqrt(tpLeading->definingParametersCovMatrix()(0,0));
      bool  isLeadingd0SigSmall                = ( d0SignificanceLeading < 3.0 );

      const xAOD::TrackParticle* tpSubLeading  = dynamic_cast<const xAOD::Muon*>(subLeadingLepton)->primaryTrackParticle();
      float d0SignificanceSubLeading           = fabs( tpSubLeading->d0() ) / sqrt(tpSubLeading->definingParametersCovMatrix()(0,0));
      bool  isSubLeadingd0SigSmall             = ( d0SignificanceSubLeading < 3.0 );

      float ptcone_dr(-999.), etcone_dr(-999.);
      if ( dynamic_cast<const xAOD::Muon*>(leadingLepton)->isolation(ptcone_dr, isoParser.parseEnum(m_TrackBasedIsoType_Mu)) &&  
           dynamic_cast<const xAOD::Muon*>(leadingLepton)->isolation(etcone_dr,isoParser.parseEnum(m_CaloBasedIsoType_Mu))	  )
      {   
         bool isTrackIsoLeading = ( ptcone_dr / (leadingLeptonPt) <  m_TrackBasedIsoCut_Mu);
         bool isCaloIsoLeading  = ( etcone_dr / (leadingLeptonPt) <  m_CaloBasedIsoCut_Mu);
         isLeadingTight                   = ( isTrackIsoLeading && isCaloIsoLeading && isLeadingd0SigSmall );      
	 isTightDecor( *leadingLepton  )  = isLeadingTight;
      }
      ptcone_dr = -999., etcone_dr= -999.;
      if ( dynamic_cast<const xAOD::Muon*>(subLeadingLepton)->isolation(ptcone_dr, isoParser.parseEnum(m_TrackBasedIsoType_Mu)) &&  
           dynamic_cast<const xAOD::Muon*>(subLeadingLepton)->isolation(etcone_dr,isoParser.parseEnum(m_CaloBasedIsoType_Mu))      )
      {   
         bool isTrackIsoSubLeading = ( ptcone_dr / (subLeadingLeptonPt) <  m_TrackBasedIsoCut_Mu);
         bool isCaloIsoSubLeading  = ( etcone_dr / (subLeadingLeptonPt) <  m_CaloBasedIsoCut_Mu);
         isSubLeadingTight                  = ( isTrackIsoSubLeading && isCaloIsoSubLeading && isSubLeadingd0SigSmall );
         isTightDecor( *subLeadingLepton  ) = isSubLeadingTight;
      }

  } else if ( isElMu ) {
     
      if ( m_debug ) { Info("execute()","ElMu event"); }  
      
      if ( leadingLeptonFlavour == xAOD::Type::Electron ) {
     
        // leading is electron
  
        if ( m_doLHPIDCut ) { 
	  isLeadingTight                   = ( m_LHToolTight2012->accept( *leadingLepton ) && m_electronIsolationSelectionTool->accept( *( dynamic_cast<const xAOD::Electron*>( leadingLepton ) ) ) ); 
	  isTightDecor( *leadingLepton  )  = isLeadingTight;
	} else if ( m_doCutBasedPIDCut ) { 
	  isLeadingTight                   = ( m_TightPP2012->accept( *leadingLepton ) && m_electronIsolationSelectionTool->accept( *( dynamic_cast<const xAOD::Electron*>( leadingLepton ) ) ) ) ; 
	  isTightDecor( *leadingLepton  )  = isLeadingTight;
	}
	
        // subleading must be muon

        const xAOD::TrackParticle* tpSubLeading  = dynamic_cast<const xAOD::Muon*>(subLeadingLepton)->primaryTrackParticle();
        float d0SignificanceSubLeading           = fabs( tpSubLeading->d0() ) / sqrt(tpSubLeading->definingParametersCovMatrix()(0,0));
        bool  isSubLeadingd0SigSmall             = ( d0SignificanceSubLeading < 3.0 );
        
	float ptcone_dr(-999.), etcone_dr(-999.);
        if ( dynamic_cast<const xAOD::Muon*>(subLeadingLepton)->isolation(ptcone_dr, isoParser.parseEnum(m_TrackBasedIsoType_Mu)) &&  
             dynamic_cast<const xAOD::Muon*>(subLeadingLepton)->isolation(etcone_dr,isoParser.parseEnum(m_CaloBasedIsoType_Mu))      )
        {   
           bool isTrackIsoSubLeading = ( ptcone_dr / (subLeadingLeptonPt) <  m_TrackBasedIsoCut_Mu);
           bool isCaloIsoSubLeading  = ( etcone_dr / (subLeadingLeptonPt) <  m_CaloBasedIsoCut_Mu);
           isSubLeadingTight                  = ( isTrackIsoSubLeading && isCaloIsoSubLeading && isSubLeadingd0SigSmall );
           isTightDecor( *subLeadingLepton  ) = isSubLeadingTight;
        }

     } else {
        
	// leading is muon
        
	const xAOD::TrackParticle* tpLeading     = dynamic_cast<const xAOD::Muon*>(leadingLepton)->primaryTrackParticle();
        float d0SignificanceLeading              = fabs( tpLeading->d0() ) / sqrt(tpLeading->definingParametersCovMatrix()(0,0));
        bool  isLeadingd0SigSmall                = ( d0SignificanceLeading < 3.0 );

        float ptcone_dr(-999.), etcone_dr(-999.);
        if ( dynamic_cast<const xAOD::Muon*>(leadingLepton)->isolation(ptcone_dr, isoParser.parseEnum(m_TrackBasedIsoType_Mu)) &&  
             dynamic_cast<const xAOD::Muon*>(leadingLepton)->isolation(etcone_dr,isoParser.parseEnum(m_CaloBasedIsoType_Mu))      )
        {   
           bool isTrackIsoLeading = ( ptcone_dr / (leadingLeptonPt) <  m_TrackBasedIsoCut_Mu);
           bool isCaloIsoLeading  = ( etcone_dr / (leadingLeptonPt) <  m_CaloBasedIsoCut_Mu);
           isLeadingTight                   = ( isTrackIsoLeading && isCaloIsoLeading && isLeadingd0SigSmall );	  
	   isTightDecor( *leadingLepton  )  = isLeadingTight;
        }

        // subleading must be electron
        
	if ( m_doLHPIDCut ) { 
	  isSubLeadingTight                  = ( m_LHToolTight2012->accept( *subLeadingLepton ) && m_electronIsolationSelectionTool->accept( *( dynamic_cast<const xAOD::Electron*>( subLeadingLepton ) ) ) ); 
          isTightDecor( *subLeadingLepton  ) = isSubLeadingTight;
	} else if ( m_doCutBasedPIDCut ) { 
	  isSubLeadingTight = ( m_TightPP2012->accept( *subLeadingLepton ) && m_electronIsolationSelectionTool->accept( *( dynamic_cast<const xAOD::Electron*>( subLeadingLepton ) ) ) ); 
          isTightDecor( *subLeadingLepton  ) = isSubLeadingTight;
	}     
     
     }
  
  } 
  
  if ( m_debug ) {
    Info("execute()","Checking *isTight* lepton decoration"); 
    for ( auto lep_it : leptonsSorted ) {
      Info("execute()","isTight: %i ",  lep_it->auxdata<char>("isTight") ); 
    }
  }
  
  // 4.
  // ----------------------------------
  // Now do the lepton trigger matching
  // ----------------------------------

  // FIXME! 
  // No xAOD ASG trigger tools (outside Athena) in DC14 (maybe there are for DAOD mc14_13TeV ) 
  // 
  // TODO: do the trigger matching in separate algo, prior to this (as it is done w/ truth matching)


  // 5.
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
    for ( auto lep_it : leptonsSorted ) {
      Info("execute()","isTag: %i ",  lep_it->auxdata<char>("isTag") ); 
    }
  }
  
  m_numEventPass++;
  m_weightNumEventPass += mcEvtWeight;
   
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode TagAndProbeRFRateMeasurement :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.

  if ( m_debug ) { Info("postExecute()", "Calling postExecute \n"); }

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TagAndProbeRFRateMeasurement :: finalize ()
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

  Info("finalize()", "Deleting tool instances... \n");

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode TagAndProbeRFRateMeasurement :: histFinalize ()
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

  Info("histFinalize()", "Calling histFinalize \n");
  
  if ( m_useCutFlow ) {
    Info("histFinalize()", "Filling cutflow");
    m_cutflowHist ->SetBinContent( m_cutflow_bin, m_numEventPass        );
    m_cutflowHistW->SetBinContent( m_cutflow_bin, m_weightNumEventPass  );
  }
  
  return EL::StatusCode::SUCCESS;
}

