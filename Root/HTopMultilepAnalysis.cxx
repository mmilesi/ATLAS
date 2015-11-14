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
#include "xAODMissingET/MissingETContainer.h"

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
  m_TightElectronPID_WP       = "LHTight";
  m_TightElectronIso_WP       = "isIsolated_Gradient";
  m_TightMuonD0sig_cut        = -1.0;	      
  m_TightMuonIso_WP           = "isIsolated_Gradient";      
  
  // BTag WP to define nbjets
  m_BTag_WP                   = "MV2c20_Fix77";
  
  m_useMCForTagAndProbe       = false;
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
 
    m_inContainerName_PreSelectedElectrons = config->GetValue("InputContainerPreSelectedElectrons", m_inContainerName_PreSelectedElectrons.c_str());
    m_inContainerName_PreSelectedMuons     = config->GetValue("InputContainerPreSelectedMuons",     m_inContainerName_PreSelectedMuons.c_str());
    m_inContainerName_PreSelectedJets	   = config->GetValue("InputContainerPreSelectedJets",      m_inContainerName_PreSelectedJets.c_str());
    
    // to define "Tight" leptons
    m_TightElectronPID_WP         = config->GetValue("TightElectronPID_WP"  ,  m_TightElectronPID_WP.c_str());
    m_TightElectronIso_WP         = config->GetValue("TightElectronIso_WP"  ,  m_TightElectronIso_WP.c_str());
    m_TightMuonD0sig_cut          = config->GetValue("TightMuonD0sig_cut"   ,  m_TightMuonD0sig_cut );
    m_TightMuonIso_WP             = config->GetValue("TightMuonIso_WP"      ,  m_TightMuonIso_WP.c_str()); 

    // BTag WP to define nbjets
    m_BTag_WP                     = config->GetValue("BTagWP"               ,  m_BTag_WP.c_str());

    m_useMCForTagAndProbe         = config->GetValue("UseMCForTagAndProbe"  ,  m_useMCForTagAndProbe );

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

  // get the MetaData tree once a new file is opened, with
  TTree *MetaData = dynamic_cast<TTree*>(wk()->inputFile()->Get("MetaData"));
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

  if ( this->configure() == EL::StatusCode::FAILURE ) {
    Error("initialize()", "Failed to properly configure. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  // check if sample is MC
  //
  const xAOD::EventInfo* eventInfo(nullptr);
  RETURN_CHECK("HTopMultilepAnalysis::initialize()", HelperFunctions::retrieve(eventInfo, "EventInfo", m_event, m_store, m_debug) , "");
  m_isMC = ( eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) );

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
    Error("initialize()", "Failed to retrieve MetaData histogram. Aborting");
    return EL::StatusCode::FAILURE;
  }

  m_totalEvents  = new TH1D("TotalEvents",  "TotalEvents",  2, 1, 3);    
  m_totalEventsW = new TH1D("TotalEventsW", "TotalEventsW", 2, 1, 3);  
  wk() -> addOutput(m_totalEvents);
  wk() -> addOutput(m_totalEventsW);
  

  m_jetPlots = new JetHists( "highPtJets", "clean" ); // second argument: "kinematic", "clean", "energy", "resolution"
  m_jetPlots -> initialize();
  m_jetPlots -> record( wk() );
  
  // initialise TauSelectionTool  
  //
  if ( asg::ToolStore::contains<TauAnalysisTools::TauSelectionTool>("TauSelectionTool") ) {
    m_TauSelTool = asg::ToolStore::get<TauAnalysisTools::TauSelectionTool>("TauSelectionTool");
  } else {
    m_TauSelTool = new TauAnalysisTools::TauSelectionTool( "TauSelectionTool" );
  }  
  m_TauSelTool->setProperty("ConfigPath", "$ROOTCOREBIN/data/HTopMultilepAnalysis/Taus/recommended_selection_mc15_final_sel.conf");
                                                     
  RETURN_CHECK( "HTopMultilepAnalysis::initialize()", m_TauSelTool->initialize(), "Failed to properly initialize TauSelectionTool." );

  // ***********************************************************
  // For MM/FF: read r/f rates from input ROOT histograms
  //

  std::string path = "$ROOTCOREBIN/data/HTopMultilepAnalysis/External/";
  path = ( !m_useMCForTagAndProbe ) ? ( path +"ObservedRates.root" ) : ( path+"ExpectedRates.root" );

  TFile *file = TFile::Open(path.c_str());
  if ( !file->IsOpen() ) {
    Error("histInitialize()", "Failed to open ROOT file with r/f rates from path: %s . Aborting", path.c_str() );
    return EL::StatusCode::FAILURE;
  }    
  
  if ( m_debug ) { Info("initialize()", " Successfully opened ROOT file with r/f rates from path: %s ", path.c_str() ); }
  
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
  static SG::AuxElement::Accessor< unsigned int > nBjets_Acc("nBjets_"+m_BTag_WP);
  unsigned int nBjets(0);
  if ( nBjets_Acc.isAvailable( *eventInfo ) ) {
    nBjets = nBjets_Acc( *eventInfo );
  } else {
    Error("execute()"," 'nBjets_%s' is not available as decoration. Aborting", m_BTag_WP.c_str() ); 
    return EL::StatusCode::FAILURE;
  }
    
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
    
    Info("execute()","Event %i ", static_cast<int>(m_eventCounter));
    Info("execute()"," Initial vs Preselected vs Selected Signal Muons: \t %u \t %u \t %u  "    , nInMuons, nPreselMuons, nSignalMuons ); 
    Info("execute()"," Initial vs Preselected vs Selected Signal Electrons: %u \t %u \t %u "    , nInElectrons, nPreselElectrons, nSignalElectrons );     
    Info("execute()"," Initial vs Preselected vs Selected Signal Jets: \t %u \t %u \t %u "      , nInJets, nPreselJets, nSignalJets );
    Info("execute()"," Selected Signal Leptons: \t %u " , nSignalLeptons );   
    Info("execute()"," Selected Signal Taus: \t %u " , nSignalTaus );   

  }
  
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
  
  //-------------------------------------------
  // definition of "Tight" and "Medium" leptons
  //-------------------------------------------
 
  static SG::AuxElement::Decorator< char > isTightDecor("isTight"); 
  static SG::AuxElement::Decorator< char > isMediumDecor("isMedium"); 

  static SG::AuxElement::Accessor< char >  TightElectronIsoAcc(m_TightElectronIso_WP); 
  static SG::AuxElement::Accessor< char >  TightElectronIDAcc(m_TightElectronPID_WP);
  static SG::AuxElement::Accessor< char >  TightMuonIsoAcc(m_TightMuonIso_WP); 
  static SG::AuxElement::Accessor< float > d0SigAcc ("d0sig");  

  const xAOD::Vertex *primaryVertex = HelperFunctions::getPrimaryVertex(vertices);

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

    const xAOD::TrackParticle* trk = el_itr->trackParticle();
    if ( !trk ) {
      Error("execute()", "no track available for this electron. Aborting " );
      return EL::StatusCode::FAILURE;
    }

    float z0 =  trk->z0()  - ( primaryVertex->z() - trk->vz() ) ; // distance between z0 and zPV ( after referring the PV z coordinate to the beamspot position, given by vz() )
    // see https://twiki.cern.ch/twiki/bin/view/AtlasProtected/InDetTrackingDC14 for further reference
    float theta = trk->theta();

    // preliminary: tighten impact parameter cuts
    //
    if ( fabs( d0SigAcc( *el_itr ) ) < 5.0 && fabs( z0*sin(theta) ) < 0.5 ) {

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
         
	  // "Tight"  ---> iso + ID
	  // "Medium" ---> !iso + ID
	  // 
	  if ( ( TightElectronIsoAcc( *el_itr ) == 1 ) && ( TightElectronIDAcc( *el_itr ) == 1 ) ) { isTightDecor( *el_itr ) = 1; }
	  else if ( ( TightElectronIDAcc( *el_itr ) == 1 ) )                                       { isMediumDecor( *el_itr ) = 1; } 
	  
	} else {

	  // "Tight"  ---> iso
	  // "Medium" ---> !iso
	  // 
	  if ( TightElectronIsoAcc( *el_itr ) == 1 ) { isTightDecor( *el_itr ) = 1; }
	  else                                       { isMediumDecor( *el_itr ) = 1; } 
	  
	}
	
      } 
      // if not using isolation, but using Electron ID...
      //
      else if ( !tightness_def_el.second.empty() ) { 
    
	if ( !TightElectronIDAcc.isAvailable( *el_itr ) ) {
	  Error("execute()", "'%s' attribute is not available for this electron. Aborting ", m_TightElectronPID_WP.c_str() );
	  return EL::StatusCode::FAILURE;
	} 

	// "Tight"  ---> ID
	// "Medium" ---> !ID
	// 
	if ( TightElectronIDAcc( *el_itr ) == 1 ) { isTightDecor( *el_itr ) = 1; }
	else                                      { isMediumDecor( *el_itr ) = 1; } 
	
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

    const xAOD::TrackParticle* trk = mu_itr->primaryTrackParticle();
    if ( !trk ) {
      Error("execute()", "no track available for this muon. Aborting " );
      return EL::StatusCode::FAILURE;
    }

    float z0 =  trk->z0()  - ( primaryVertex->z() - trk->vz() ) ; // distance between z0 and zPV ( after referring the PV z coordinate to the beamspot position, given by vz() )
    // see https://twiki.cern.ch/twiki/bin/view/AtlasProtected/InDetTrackingDC14 for further reference
    float theta = trk->theta();

    // preliminary: tighten impact parameter cuts
    //
    if ( fabs( z0*sin(theta) ) < 0.5 ) {

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

	  // "Tight"  ---> iso + d0sig
	  // "Medium" ---> !iso + d0sig
	  // 
	  if ( ( TightMuonIsoAcc( *mu_itr ) == 1 ) && ( fabs( d0SigAcc( *mu_itr ) ) < tightness_def_mu.second ) ) { isTightDecor( *mu_itr ) = 1; }
	  else if ( fabs( d0SigAcc( *mu_itr ) ) < tightness_def_mu.second )                                       { isMediumDecor( *mu_itr ) = 1; }
	  
	} else {

	  // "Tight"  ---> iso 
	  // "Medium" ---> !iso
	  // 
	  if ( TightMuonIsoAcc( *mu_itr ) == 1 ) { isTightDecor( *mu_itr ) = 1; }
	  else                                   { isMediumDecor( *mu_itr ) = 1; }

	}
	
      } 
      // if not using isolation, but using d0sig...
      //
      else if ( tightness_def_mu.second > 0.0 ) { 
	
	if ( !d0SigAcc.isAvailable( *mu_itr ) ) {
	  Error("execute()", "'d0sig' attribute is not available for this muon. Aborting " );
	  return EL::StatusCode::FAILURE;
	}   

	// "Tight"  ---> d0sig
	// "Medium" ---> !d0sig
	// 
	if ( fabs( d0SigAcc( *mu_itr ) ) < tightness_def_mu.second ) { isTightDecor( *mu_itr ) = 1; }
	else                                                         { isMediumDecor( *mu_itr ) = 1; }

      } 
      // if using neither isolation, nor d0sig..
      //
      else {
	Error("execute()", "Need at least isolation or d0sig requirement to define 'Tight' muons. Aborting" );
	return EL::StatusCode::FAILURE;
      }
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
  
  if ( nSignalLeptons == 2 && ( nSignalJets >= 1 && nSignalJets <= 3 ) &&  nBjets >= 1 ) {

    if ( m_useMCForTagAndProbe && m_isMC ) { this->defineTagAndProbeRFRateVars_MC( eventInfo, leptonsSorted ); }
    else { this->defineTagAndProbeRFRateVars( eventInfo, leptonsSorted ); }
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

  // if MetaData is not empty, and file is a DAOD, use it (to take DAOD skimming into account!)
  //
  if ( m_isDerivation && m_histEventCount->GetBinContent(1) > 0 && m_histEventCount->GetBinContent(3) > 0 ) {
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
  
//***************************************************
//
// Set Tag&Probe variables in r/f rate measurement CR 
//
//
EL::StatusCode HTopMultilepAnalysis :: defineTagAndProbeRFRateVars( const xAOD::EventInfo* eventInfo, const xAOD::IParticleContainer& leptons )
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
  static SG::AuxElement::Accessor< char > isTrigMatchedLepAcc("isTrigMatchedLep");
    
  // decorate with default values
  for ( auto lep_it : leptons ) { isTagDecor( *lep_it ) = 0; }
  isNonTightEventDecor( *eventInfo ) = 0;
  isProbeElEventDecor( *eventInfo )  = 0;
  isProbeMuEventDecor( *eventInfo )  = 0;
    
  // -------------------------------------------
  // Now, take the leading and subleading lepton 
  // -------------------------------------------
    
  const xAOD::IParticle* leadingLepton           = leptons.at(0);
  xAOD::Type::ObjectType leadingLeptonFlavour    = leadingLepton->type();
  bool                   isLeadingTight(false);          
  if ( isTightAcc.isAvailable( *leadingLepton ) ) {
    isLeadingTight = isTightAcc( *leadingLepton );
  } else {
    Warning("defineTagAndProbeRFRateVars()","SG::AuxElement::Accessor('isTight') is not available for the leading lepton. Should not happen. Assigning 'isTight' = false" );
  }
  bool 			 isLeadingTrigMatched(false);  
  if ( isTrigMatchedLepAcc.isAvailable( *leadingLepton ) ) {
    isLeadingTrigMatched = isTrigMatchedLepAcc( *leadingLepton );
  } else {
    Warning("defineTagAndProbeRFRateVars()","SG::AuxElement::Accessor('isTrigMatchedLep') is not available for the leading lepton. Should not happen. Assigning 'isLeadingTrigMatched' = false" );
  }
  
  const xAOD::IParticle* subLeadingLepton        = leptons.at(1);
  xAOD::Type::ObjectType subLeadingLeptonFlavour = subLeadingLepton->type();
  bool                   isSubLeadingTight (false);
  if ( isTightAcc.isAvailable( *subLeadingLepton ) ) {
    isLeadingTight = isTightAcc( *subLeadingLepton );
  } else {
    Warning("defineTagAndProbeRFRateVars()","SG::AuxElement::Accessor('isTight') is not available for the subleading lepton. Should not happen. Assigning 'isTight' = false" );
  }
  bool 			 isSubLeadingTrigMatched(false);  
  if ( isTrigMatchedLepAcc.isAvailable( *subLeadingLepton ) ) {
    isSubLeadingTrigMatched = isTrigMatchedLepAcc( *subLeadingLepton );
  } else {
    Warning("defineTagAndProbeRFRateVars()","SG::AuxElement::Accessor('isTrigMatchedLep') is not available for the subleading lepton. Should not happen. Assigning 'isSubLeadingTrigMatched' = false" );
  }
    
  if ( m_debug ) { Info("defineTagAndProbeRFRateVars()","Leading lepton: isTight: %i \t isTrigMatched: %i ", isLeadingTight, isLeadingTrigMatched ); }
  if ( m_debug ) { Info("defineTagAndProbeRFRateVars()","Subleading lepton: isTight: %i \t isTrigMatched: %i ", isSubLeadingTight, isSubLeadingTrigMatched ); }
  
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
    if ( ( isLeadingTight && isSubLeadingTight && isLeadingTrigMatched && isSubLeadingTrigMatched ) || 
         ( isLeadingTight && isLeadingTrigMatched ) 
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
    else if ( isSubLeadingTight  && isSubLeadingTrigMatched ) {
      
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
    if ( isLeadingTight && isSubLeadingTight && isLeadingTrigMatched && isSubLeadingTrigMatched ) {
        
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
    else if ( isLeadingTight && isLeadingTrigMatched ) {
        
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

  // accessor to tag leptons 
  static SG::AuxElement::Accessor< char > isTagAcc("isTag");
  
  if ( m_debug ) {
    Info("execute()","Checking *isTag* lepton decoration"); 
    for ( auto lep_it : leptons ) {
      Info("execute()","\t lepton isTag: %i ", isTagAcc( *lep_it ) ); 
    }
  }
   
  return EL::StatusCode::SUCCESS;

}

// ******************************************
// Use this function only when doing
// MM estimate on pure MC (i.e, closure test)
//
EL::StatusCode HTopMultilepAnalysis :: defineTagAndProbeRFRateVars_MC( const xAOD::EventInfo* eventInfo, const xAOD::IParticleContainer& leptons )
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
  bool isSS = (  prod_lep_charge > 0  );  

  // --------------------------------
  // Declare accessors and decorators
  // --------------------------------
  
  // decorate lepton w/ is tag/probe
  static SG::AuxElement::Decorator< char > isTagDecor("isTag");
  // declare event decorations for checking the probe type in event
  static SG::AuxElement::Decorator< char > isProbeElEventDecor( "isProbeElEvent" ); 
  static SG::AuxElement::Decorator< char > isProbeMuEventDecor( "isProbeMuEvent" ); 
    
  // decorate with default values
  for ( auto lep_itr : leptons ) { isTagDecor( *lep_itr ) = 0; }
  isProbeElEventDecor( *eventInfo )  = 0;
  isProbeMuEventDecor( *eventInfo )  = 0;
    
  // -------------------------------------------
  // Now, take the leading and subleading lepton 
  // -------------------------------------------
    
  const xAOD::IParticle* leadingLepton           = leptons.at(0);
  const xAOD::IParticle* subLeadingLepton        = leptons.at(1);

  // --------------------------------------------
  // Now decide who's the tag and who's the probe
  // --------------------------------------------
  
  static SG::AuxElement::Accessor< char > isChFlipAcc("isChFlip");
  static SG::AuxElement::ConstAccessor< int > truthTypeAcc("truthType"); 
  static SG::AuxElement::ConstAccessor< int > truthOriginAcc("truthOrigin");
   
  bool found_a_prompt(false);
  if ( isSS ) {
  
    for ( auto lep_itr : leptons ) {
       
      // The tag will be the first prompt lepton in the event which is not charge flip, provided it's found.
      // The other will be the probe
      // See below the treatment for the case where all leptons in SS event are non prompt...
      //
      if ( truthTypeAcc.isAvailable( *lep_itr ) ) {
      
        if ( !found_a_prompt && ( truthTypeAcc( *lep_itr ) == 2 || truthTypeAcc( *lep_itr ) == 6 ) && isChFlipAcc( *lep_itr ) != 1 ) {
          isTagDecor( *lep_itr ) = 1;
	  found_a_prompt = true;
        }
      
      } else {
        Warning("defineTagAndProbeRFRateVars_MC()","SG::AuxElement::Accessor('truthType') is not available for this lepton. Should not happen. Tag will be the leading, probe the subleading" ); 
        break;
      }
    
    }
   
    if ( !found_a_prompt ) {
      // the probe will be the subleading lepton, no matter what
      //
      if ( m_debug ) { Info("defineTagAndProbeRFRateVars_MC()","There are no prompt leptons in this SS event. Tag will b\
e the leading, probe the subleading" ); }
      isTagDecor( *leadingLepton ) = 1;
    }
     
  } else {
    
    // generate a uniformly distributed random number in [0,1]			  
    //
    TRandom3 rndm(0);
    float unif = rndm.Uniform();
    
    // assign randomly who's tag and who's probe 
    //	
    if ( unif >= 0.5 ) {
      isTagDecor( *leadingLepton )    = 1;
    } else {
      isTagDecor( *subLeadingLepton ) = 1;
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

     const xAOD::MissingET* final = *inMETCont->find("FinalClus"); // ("FinalClus" uses the calocluster-based soft terms, "FinalTrk" uses the track-based ones)
     TLorentzVector MET;
     MET.SetPtEtaPhiM( final->met(), 0, final->phi(), 0 );

     float mT(-1.0);
     mT = sqrt( 2.0 * lepA_4mom.Pt() * ( MET.Pt() ) * ( 1.0 - cos( lepA_4mom.DeltaPhi(MET) ) ) );
     mTLep0METDecor( *eventInfo ) = mT;  
     mT = sqrt( 2.0 * lepB_4mom.Pt() * ( MET.Pt() ) * ( 1.0 - cos( lepB_4mom.DeltaPhi(MET) ) ) );
     mTLep1METDecor( *eventInfo ) = mT;

  }
  else if ( nLeptons == 3 )
  {
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

  // retrieve some previously applied event object decorations
  // 
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
  //
  static SG::AuxElement::Accessor< char > isTightAcc("isTight");
  static SG::AuxElement::Accessor< char > isMediumAcc("isMedium");
   
  // event decorators to identify the TT,TL,LT,LL regions (ordering depends on the category)
  //
  static SG::AuxElement::Decorator< char > isTTDecor("isTT");
  static SG::AuxElement::Decorator< char > isTLDecor("isTL");
  static SG::AuxElement::Decorator< char > isLTDecor("isLT");
  static SG::AuxElement::Decorator< char > isLLDecor("isLL");
  isTTDecor( *eventInfo ) = 0;
  isTLDecor( *eventInfo ) = 0;
  isLTDecor( *eventInfo ) = 0;
  isLLDecor( *eventInfo ) = 0;

  static SG::AuxElement::Decorator< char > isTMDecor("isTM");
  static SG::AuxElement::Decorator< char > isMTDecor("isMT");
  static SG::AuxElement::Decorator< char > isMMDecor("isMM");
  isTMDecor( *eventInfo ) = 0;
  isMTDecor( *eventInfo ) = 0;
  isMMDecor( *eventInfo ) = 0;
  
  // These will be the two leptons used for the fake estimate and to define the "tightness" of the region, both in 2 lep SS and in 3 lep category
  //
  xAOD::IParticle* lepA(nullptr);
  xAOD::IParticle* lepB(nullptr);
  
  // The "real" lepton for 3 lep category
  //
  xAOD::IParticle* lep0(nullptr);
  
  // Features of the two leptons that are considered for the fake estimate
  //
  float lepA_pt(-1.0), lepA_eta(-999.0), lepB_pt(-1.0), lepB_eta(-999.0);
  int lepA_flavour(0), lepB_flavour(0);

  std::string region;
   
  if ( nLeptons == 2 && isSS01(*eventInfo) ) {
        
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
	
	// set the region string
	//
	if      (  isTightAcc( *lepA )    &&  isTightAcc( *lepB )    ) { region = "TT"; isTTDecor( *eventInfo ) = 1; }
	else if (  isTightAcc( *lepA )    &&  ( !(isTightAcc( *lepB )) && !(isMediumAcc( *lepB )) ) ) { region = "TL"; isTLDecor( *eventInfo ) = 1; }
	else if (  ( !(isTightAcc( *lepA )) && !(isMediumAcc( *lepA )) ) &&  isTightAcc( *lepB )    ) { region = "LT"; isLTDecor( *eventInfo ) = 1; }
	else if (  ( !(isTightAcc( *lepA )) && !(isMediumAcc( *lepA )) ) &&  ( !(isTightAcc( *lepB )) && !(isMediumAcc( *lepB )) ) ) { region = "LL"; isLLDecor( *eventInfo ) = 1; }
	else if (  isTightAcc( *lepA )    &&  isMediumAcc( *lepB ) ) { region = "TM"; isTMDecor( *eventInfo ) = 1; }
	else if (  isMediumAcc( *lepA )   &&  isTightAcc( *lepB )  ) { region = "MT"; isMTDecor( *eventInfo ) = 1; }
	else if (  isMediumAcc( *lepA )   &&  isMediumAcc( *lepB ) ) { region = "MM"; isMMDecor( *eventInfo ) = 1; }

  	if ( m_debug ) { Info("fakeWeightCalculator()", "Dilepton SS category. Region is %s ", region.c_str() ); }

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
	
  } else if ( nLeptons == 3 && isSS12(*eventInfo) ) {        
        
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
	  
	  if ( !isOSlep.isAvailable(*lep_it) ) {
	    Error("fakeWeightCalculator()", "isOSlep lepton decoration is not available. Aborting ");
    	    return EL::StatusCode::FAILURE;
 	  } 	
	  
	  if ( isOSlep(*this_lep) ) {
	    
	    // retrieve lep0 : the OS lepton
	    //
	    lep0 = this_lep;
	    continue;
	    
	  } else {
	  
	    if ( !isClosestSSlep.isAvailable(*this_lep) ) {
	      Error("fakeWeightCalculator()", "isClosestSSlep lepton decoration is not available. Aborting");
    	      return EL::StatusCode::FAILURE;
 	    } 	
	    
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
	
	// just a safety check
	//
	if ( !( lep0 && lepA && lepB ) ) {
	  Error("fakeWeightCalculator()", "Trilepton region, but no (lep0 and lepA and lepB) pointers! Aborting");
	  return EL::StatusCode::FAILURE;
	}
	
	// set the region string
	//
	if      (  isTightAcc( *lepA )    &&  isTightAcc( *lepB )    ) { region = "TT"; isTTDecor( *eventInfo ) = 1; }
	else if (  isTightAcc( *lepA )    &&  ( !(isTightAcc( *lepB )) && !(isMediumAcc( *lepB )) ) ) { region = "TL"; isTLDecor( *eventInfo ) = 1; }
	else if (  ( !(isTightAcc( *lepA )) && !(isMediumAcc( *lepA )) ) &&  isTightAcc( *lepB )    ) { region = "LT"; isLTDecor( *eventInfo ) = 1; }
	else if (  ( !(isTightAcc( *lepA )) && !(isMediumAcc( *lepA )) ) &&  ( !(isTightAcc( *lepB )) && !(isMediumAcc( *lepB )) ) ) { region = "LL"; isLLDecor( *eventInfo ) = 1; }
	else if (  isTightAcc( *lepA )    &&  isMediumAcc( *lepB ) ) { region = "TM"; isTMDecor( *eventInfo ) = 1; }
	else if (  isMediumAcc( *lepA )   &&  isTightAcc( *lepB )  ) { region = "MT"; isMTDecor( *eventInfo ) = 1; }
	else if (  isMediumAcc( *lepA )   &&  isMediumAcc( *lepB ) ) { region = "MM"; isMMDecor( *eventInfo ) = 1; }
  	
	if ( m_debug ) { Info("fakeWeightCalculator()", "Trilepton 2SS+1OS category. Region (defined by the 2SS leptons) is %s ", region.c_str() ); }
	
	// set the properties of the two SS leptons for future convenience
	//
	lepA_pt  = lepA->pt();
	lepA_eta = lepA->eta();
	if ( lepA->type() == xAOD::Type::Electron )  { lepA_flavour = 11; } 
	else if ( lepA->type() == xAOD::Type::Muon ) { lepA_flavour = 13; }
	
	lepB_pt  = lepB->pt();
	lepB_eta = lepB->eta();
	if ( lepB->type() == xAOD::Type::Electron )  { lepB_flavour = 11; } 
	else if ( lepB->type() == xAOD::Type::Muon ) { lepB_flavour = 13; }	
  
  } else  {
    return EL::StatusCode::SUCCESS; //no weights in the other categories
  }

  // Save a special flag for TL,LT OF events (used for "theta factior" ABCD method)
  //
  static SG::AuxElement::Decorator< char > isTelLmuDecor("isTelLmu");
  static SG::AuxElement::Decorator< char > isTmuLelDecor("isTmuLel");
  static SG::AuxElement::Decorator< char > isLelTmuDecor("isLelTmu");
  static SG::AuxElement::Decorator< char > isLmuTelDecor("isLmuTel");
  isTelLmuDecor( *eventInfo ) = 0;
  isTmuLelDecor( *eventInfo ) = 0;
  isLelTmuDecor( *eventInfo ) = 0;
  isLmuTelDecor( *eventInfo ) = 0;  
  static SG::AuxElement::Decorator< char > isTelMmuDecor("isTelMmu");
  static SG::AuxElement::Decorator< char > isTmuMelDecor("isTmuMel");
  static SG::AuxElement::Decorator< char > isMelTmuDecor("isMelTmu");
  static SG::AuxElement::Decorator< char > isMmuTelDecor("isMmuTel");
  isTelMmuDecor( *eventInfo ) = 0;
  isTmuMelDecor( *eventInfo ) = 0;
  isMelTmuDecor( *eventInfo ) = 0;
  isMmuTelDecor( *eventInfo ) = 0;  

  bool OF = ( ( lepA_flavour == 11 && lepB_flavour == 13 ) || ( lepA_flavour == 13 && lepB_flavour == 11 ) );
  if ( OF ) {
    if ( region == "TL" ) {
      // the tight is an electron, the loose is a muon
      if (  lepA_flavour == 11 )     { isTelLmuDecor( *eventInfo ) = 1; }
      // the tight is a muon, the loose is an electron
      else if ( lepA_flavour == 13 ) { isTmuLelDecor( *eventInfo ) = 1; }
    } else if ( region == "LT" ) {
      // the loose is an electron, the tight is a muon
      if (  lepA_flavour == 11 )     { isLelTmuDecor( *eventInfo ) = 1; }
      // the loose is a muon, the tight is an electron
      else if ( lepA_flavour == 13 ) { isLmuTelDecor( *eventInfo ) = 1; }
    } else if ( region == "TM" ) {
      // the tight is an electron, the medium is a muon
      if (  lepA_flavour == 11 )     { isTelMmuDecor( *eventInfo ) = 1; }
      // the tight is a muon, the medium is an electron
      else if ( lepA_flavour == 13 ) { isTmuMelDecor( *eventInfo ) = 1; }
    } else if ( region == "MT" ) {
      // the medium is an electron, the tight is a muon
      if (  lepA_flavour == 11 )     { isMelTmuDecor( *eventInfo ) = 1; }
      // the medium is a muon, the tight is an electron
      else if ( lepA_flavour == 13 ) { isMmuTelDecor( *eventInfo ) = 1; }
    }
  }

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

double  HTopMultilepAnalysis :: scaleRateToFactor( double rate )
{
  if ( rate < 0 ) { rate = 0.0; }
  
  double factor = ( rate /(rate+1.0) );
  
  return factor;
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
  
  weights.at(0) = scaleRateToFactor(weights.at(0));
  weights.at(1) = scaleRateToFactor(weights.at(1));
  weights.at(2) = scaleRateToFactor(weights.at(2));
  
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

