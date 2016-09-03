// package include(s):
#include <HTopMultilepAnalysis/HTopMultilepNTupReprocesser.h>
#include "HTopMultilepAnalysis/tools/HTopReturnCheck.h"

// ASG status code check
#include <AsgTools/MessageCheck.h>

// ROOT include(s)
#include "TObjArray.h"

// C++ include(s)
#include <iomanip>

using namespace NTupReprocesser;

// this is needed to distribute the algorithm to the workers
ClassImp(HTopMultilepNTupReprocesser)

HTopMultilepNTupReprocesser :: HTopMultilepNTupReprocesser(std::string className) :
    Algorithm(className),
    m_inputNTuple(nullptr),
    m_outputNTuple(nullptr),
    m_isQMisIDBranchIn(false),
    m_isMMBranchIn(false)
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().

  Info("HTopMultilepNTupReprocesser()", "Calling constructor");

  m_inputBranches       = "";

  m_outputNTupName       = "physics";
  m_outputNTupStreamName = "output";

  m_weightToCalc         = "";
  m_doQMisIDWeighting    = false;
  m_doMMWeighting        = false;

  m_QMisIDRates_dir            = "";
  m_QMisIDRates_Filename_T     = "";
  m_QMisIDRates_Filename_AntiT = "";
  m_useTAntiTRates             = false;

  m_RR_dir                  = "";
  m_FR_dir                  = "";
  m_RRFR_YES_TM_dir         = "";
  m_RRFR_NO_TM_dir          = "";
  m_Efficiency_Filename     = "";
  m_doMMClosure             = false;
  m_useEtaParametrisation   = false;
  m_useTrigMatchingInfo     = false;
  m_useScaledFakeEfficiency = false;
  
}



EL::StatusCode HTopMultilepNTupReprocesser :: setupJob (EL::Job& job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.

  //ANA_CHECK_SET_TYPE (EL::StatusCode); // set type of return code you are expecting (add to top of each function once)

  Info("setupJob()", "Calling setupJob");

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepNTupReprocesser :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.

  //ANA_CHECK_SET_TYPE (EL::StatusCode);

  Info("histInitialize()", "Calling histInitialize");

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepNTupReprocesser :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed

  //ANA_CHECK_SET_TYPE (EL::StatusCode);

  Info("fileExecute()", "Calling fileExecute");

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepNTupReprocesser :: changeInput (bool firstFile)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.

  ANA_CHECK_SET_TYPE (EL::StatusCode);

  firstFile = firstFile;

  Info("changeInput()", "Calling changeInput. Now reading file : %s", wk()->inputFile()->GetName() );

  // Get the pointer to the main input TTree
  //
  m_inputNTuple = wk()->tree();

  // Check content of input tree
  //
  TObjArray* branches = m_inputNTuple->GetListOfBranches();
  int nbranches = branches->GetEntriesFast();
  for ( int idx(0); idx < nbranches; ++idx ) {
      if ( strcmp( branches->At(idx)->GetName(), "QMisIDWeight" ) == 0 ) {
	  m_isQMisIDBranchIn = true;
	  break;
      }
  }
  for ( int idx(0); idx < nbranches; ++idx ) {
      if ( strcmp( branches->At(idx)->GetName(), "MMWeight" ) == 0 ) {
	  m_isMMBranchIn = true;
	  break;
      }
  }

  ANA_CHECK( this->enableSelectedBranches() );

  // Connect the branches of the input tree to the algorithm members
  //
  m_inputNTuple->SetBranchAddress ("EventNumber",   			      &m_EventNumber);
  m_inputNTuple->SetBranchAddress ("RunNumber",   			      &m_RunNumber);
  m_inputNTuple->SetBranchAddress ("mc_channel_number",                       &m_mc_channel_number);
  m_inputNTuple->SetBranchAddress ("isSS01",                                  &m_isSS01);
  m_inputNTuple->SetBranchAddress ("dilep_type",  			      &m_dilep_type);
  m_inputNTuple->SetBranchAddress ("trilep_type",  			      &m_trilep_type);
  m_inputNTuple->SetBranchAddress ("is_T_T",				      &m_is_T_T);
  m_inputNTuple->SetBranchAddress ("is_T_AntiT",			      &m_is_T_AntiT);
  m_inputNTuple->SetBranchAddress ("is_AntiT_T",			      &m_is_AntiT_T);
  m_inputNTuple->SetBranchAddress ("is_AntiT_AntiT",			      &m_is_AntiT_AntiT);

  m_inputNTuple->SetBranchAddress ("lep_ID_0",   			      &m_lep_ID_0);
  m_inputNTuple->SetBranchAddress ("lep_Pt_0",  			      &m_lep_Pt_0);
  m_inputNTuple->SetBranchAddress ("lep_E_0",   			      &m_lep_E_0);
  m_inputNTuple->SetBranchAddress ("lep_Eta_0",  			      &m_lep_Eta_0);
  m_inputNTuple->SetBranchAddress ("lep_Phi_0",   			      &m_lep_Phi_0);
  m_inputNTuple->SetBranchAddress ("lep_EtaBE2_0",   			      &m_lep_EtaBE2_0);
  m_inputNTuple->SetBranchAddress ("lep_isTightSelected_0",   		      &m_lep_isTightSelected_0);
  m_inputNTuple->SetBranchAddress ("lep_isTrigMatch_0",   		      &m_lep_isTrigMatch_0);

  m_inputNTuple->SetBranchAddress ("lep_ID_1",   			      &m_lep_ID_1);
  m_inputNTuple->SetBranchAddress ("lep_Pt_1",  			      &m_lep_Pt_1);
  m_inputNTuple->SetBranchAddress ("lep_E_1",   			      &m_lep_E_1);
  m_inputNTuple->SetBranchAddress ("lep_Eta_1",  			      &m_lep_Eta_1);
  m_inputNTuple->SetBranchAddress ("lep_Phi_1",   			      &m_lep_Phi_1);
  m_inputNTuple->SetBranchAddress ("lep_EtaBE2_1",   			      &m_lep_EtaBE2_1);
  m_inputNTuple->SetBranchAddress ("lep_isTightSelected_1",   		      &m_lep_isTightSelected_1);
  m_inputNTuple->SetBranchAddress ("lep_isTrigMatch_1",   		      &m_lep_isTrigMatch_1);

  if ( m_isQMisIDBranchIn ) {  m_inputNTuple->SetBranchAddress ("QMisIDWeight",  &m_QMisIDWeight_in); }
  if ( m_isMMBranchIn )     {  m_inputNTuple->SetBranchAddress ("MMWeight",      &m_MMWeight_in);     }

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepNTupReprocesser :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.

  ANA_CHECK_SET_TYPE (EL::StatusCode);

  Info("initialize()", "Initialising HTopMultilepNTupReprocesser...");

  m_outputNTuple = EL::getNTupleSvc (wk(), m_outputNTupStreamName);

  if       ( m_weightToCalc.compare("QMisID") == 0 ) { m_doQMisIDWeighting = true; }
  else if  ( m_weightToCalc.compare("MM") == 0 )     { m_doMMWeighting = true; }
  else {
      Error("initialize()","Weight %s is not known. Aborting.", m_weightToCalc.c_str() );
      return EL::StatusCode::FAILURE;
  }

  // Set new branches for output TTree
  //
  if ( m_doQMisIDWeighting ) { m_outputNTuple->tree()->Branch("QMisIDWeight",  &m_QMisIDWeight_out); }
  if ( m_doMMWeighting )     { m_outputNTuple->tree()->Branch("MMWeight",      &m_MMWeight_out); }

  // ---------------------------------------------------------------------------------------------------------------

  // Initialise counter for input TTree entries
  //
  m_numEntry = 0;

  // Initialise counter for events where inf/nan is read
  //
  m_count_inf = 0;

  // ---------------------------------------------------------------------------------------------------------------

  m_outputNTuple->tree()->SetName( m_outputNTupName.c_str() );

  // ---------------------------------------------------------------------------------------------------------------

  // Copy input TTree weight to output TTree
  
  m_outputNTuple->tree()->SetWeight( m_inputNTuple->GetWeight() );

  // ---------------------------------------------------------------------------------------------------------------

  if ( m_doQMisIDWeighting ) {  
      Info("initialize()","Reading QMisID rates from ROOT file(s)..");
      ANA_CHECK( this->readQMisIDRates() );
  }
  if ( m_doMMWeighting ) {  
      Info("initialize()","Reading MM efficiencies from ROOT file(s)..");
      ANA_CHECK( this->readRFEfficiencies() );
  }  
  
  // ---------------------------------------------------------------------------------------------------------------

  Info("initialize()", "All good!");

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepNTupReprocesser :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  ANA_CHECK_SET_TYPE (EL::StatusCode);

  if ( m_numEntry == 0 ) { Info("execute()", "Processing input TTree : %s\n", m_inputNTuple->GetName() ); }

  m_inputNTuple->GetEntry( wk()->treeEntry() );

  if ( m_debug ) { Info("execute()", "===> Entry %u - EventNumber = %u ", static_cast<uint32_t>(m_numEntry), static_cast<uint32_t>(m_EventNumber) ); }

  ++m_numEntry;

  if ( m_numEntry > 0 && ( static_cast<int>(m_numEntry) % 20000 == 0 ) ) { Info("execute()","Processed %u entries", static_cast<uint32_t>(m_numEntry)); }

  // ------------------------------------------------------------------------

  // This call is crucial, otherwise you'll get no entries in the output tree!
  //
  m_outputNTuple->setFilterPassed();

  // ------------------------------------------------------------------------

  m_event = std::make_shared<eventObj>();

  // ------------------------------------------------------------------------

  auto lep0 = std::make_shared<leptonObj>();

  lep0.get()->pt            = m_lep_Pt_0;
  lep0.get()->eta           = m_lep_Eta_0;
  lep0.get()->etaBE2        = m_lep_EtaBE2_0;
  lep0.get()->ID            = m_lep_ID_0;
  lep0.get()->flavour       = abs(m_lep_ID_0);
  lep0.get()->charge        = m_lep_ID_0 / fabs(m_lep_ID_0);
  lep0.get()->tightselected = m_lep_isTightSelected_0;
  lep0.get()->trigmatched   = m_lep_isTrigMatch_0;

  m_leptons.push_back(lep0);

  auto lep1 = std::make_shared<leptonObj>();

  lep1.get()->pt            = m_lep_Pt_1;
  lep1.get()->eta           = m_lep_Eta_1;
  lep1.get()->etaBE2        = m_lep_EtaBE2_1;
  lep1.get()->ID            = m_lep_ID_1;
  lep1.get()->flavour       = abs(m_lep_ID_1);
  lep1.get()->charge        = m_lep_ID_1 / fabs(m_lep_ID_1);
  lep1.get()->tightselected = m_lep_isTightSelected_1;
  lep1.get()->trigmatched   = m_lep_isTrigMatch_1;

  m_leptons.push_back(lep1);

  m_event.get()->isMC   = ( m_mc_channel_number > 0 );
  m_event.get()->dilep  = ( m_dilep_type > 0 );
  m_event.get()->isSS01 = ( m_isSS01 );
  
  m_event.get()->TT         = m_is_T_T;	  
  m_event.get()->TAntiT     = m_is_T_AntiT;	  
  m_event.get()->AntiTT     = m_is_AntiT_T;	  
  m_event.get()->AntiTAntiT = m_is_AntiT_AntiT; 

  if ( m_debug ) {
      Info("execute()","lep0:\n pT = %.2f\n etaBE2 = %.2f\n eta = %.2f\n flavour = %i\n tight? %i\n trigmatched? %i", lep0.get()->pt/1e3, lep0.get()->etaBE2, lep0.get()->eta, lep0.get()->flavour, lep0.get()->tightselected, lep0.get()->trigmatched );
      Info("execute()","lep1:\n pT = %.2f\n etaBE2 = %.2f\n eta = %.2f\n flavour = %i\n tight? %i\n trigmatched? %i", lep1.get()->pt/1e3, lep1.get()->etaBE2, lep1.get()->eta, lep1.get()->flavour, lep1.get()->tightselected, lep1.get()->trigmatched );
      Info("execute()","event:\n TT ? %i, TAntiT ? %i, AntiTT ? %i, AntiTAntiT ? %i", m_event.get()->TT, m_event.get()->TAntiT, m_event.get()->AntiTT, m_event.get()->AntiTAntiT );
  }

  if ( m_debug ) {
      if ( m_doQMisIDWeighting ) {
	  unsigned int idx(0);
	  if ( !m_isQMisIDBranchIn ) {
	      for ( const auto& itr : m_event.get()->weight_QMisID ) {
		  if ( idx == 0 ) Info("execute()","\t\tDefault QMisIDWeight[%i] = %.3f", idx, itr );
		  else            Info("execute()","\t\tDefault QMisIDWeight[%i] ( QMisIDWeight[%i] * QMisIDWeight[0] ) = %.3f ( %.3f )", idx, idx, itr, ( itr * *(m_event.get()->weight_QMisID.begin()) )  );
		  ++idx;
	      }
	  } else {
	      for ( const auto& itr : *m_QMisIDWeight_in ) {
		  if ( idx == 0 ) Info("execute()","\t\tIN QMisIDWeight[%i] = %.3f", idx, itr );
		  else            Info("execute()","\t\tIN QMisIDWeight[%i] ( QMisIDWeight[%i] * QMisIDWeight[0] ) = %.3f ( %.3f )", idx, idx, itr, ( itr * *(m_QMisIDWeight_in->begin()) )  );
		  ++idx;
	      }
	  }
      } else if ( m_doMMWeighting ) {
	  unsigned int idx(0);
	  if ( !m_isMMBranchIn ) {
	      for ( const auto& itr : m_event.get()->weight_MM ) {
		  if ( idx == 0 ) Info("execute()","\t\tDefault MMWeight[%i] = %.3f", idx, itr );
		  else            Info("execute()","\t\tDefault MMWeight[%i] ( MMWeight[%i] * MMWeight[0] ) = %.3f ( %.3f )", idx, idx, itr, ( itr * *(m_event.get()->weight_MM.begin()) )  );
		  ++idx;
	      }
	  } else {
	      for ( const auto& itr : *m_MMWeight_in ) {
		  if ( idx == 0 ) Info("execute()","\t\tIN MMWeight[%i] = %.3f", idx, itr );
		  else            Info("execute()","\t\tIN MMWeight[%i] ( MMWeight[%i] * MMWeight[0] ) = %.3f ( %.3f )", idx, idx, itr, ( itr * *(m_MMWeight_in->begin()) )  );		  
		  ++idx;
	      }
	  }
      }
  }

  // ------------------------------------------------------------------------

  if ( m_doQMisIDWeighting ) {
      ANA_CHECK( this->calculateQMisIDWeights () );
  } 
  if ( m_doMMWeighting ) {
      ANA_CHECK( this->calculateMMWeights () );
  }

  // ------------------------------------------------------------------------

  ANA_CHECK( this->setOutputBranches() );

  // ------------------------------------------------------------------------

  if ( m_debug ) {
      if ( m_doQMisIDWeighting ) {
	  unsigned int idx(0);
	  for ( const auto& itr : m_QMisIDWeight_out ) {
	      if ( idx == 0 ) Info("execute()","\t\tOUT QMisIDWeight[%i] = %.3f", idx, itr );
              else            Info("execute()","\t\tOUT QMisIDWeight[%i] ( QMisIDWeight[%i] * QMisIDWeight[0] ) = %.3f ( %.3f )", idx, idx, itr, ( itr * *(m_QMisIDWeight_out.begin()) )  );
	      ++idx;
	  }
      } 
      if ( m_doMMWeighting ) {
	  unsigned int idx(0);
	  for ( const auto& itr : m_MMWeight_out ) {
	      if ( idx == 0 ) Info("execute()","\t\tOUT MMWeight[%i] = %.3f", idx, itr );
              else            Info("execute()","\t\tOUT MMWeight[%i] ( MMWeight[%i] * MMWeight[0] ) = %.3f ( %.3f )", idx, idx, itr, ( itr * *(m_MMWeight_out.begin()) )  );
	      ++idx;
	  }
      }
  }

  // ------------------------------------------------------------------------

  m_leptons.clear();

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepNTupReprocesser :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.

  //ANA_CHECK_SET_TYPE (EL::StatusCode);

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepNTupReprocesser :: finalize ()
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

  //ANA_CHECK_SET_TYPE (EL::StatusCode);

  Info("finalize()", "Finalising HTopMultilepNTupReprocesser...");

  Info("finalize()", "Events where inf/nan input was read: %u", m_count_inf );

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepNTupReprocesser :: histFinalize ()
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

  //ANA_CHECK_SET_TYPE (EL::StatusCode);

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode HTopMultilepNTupReprocesser :: enableSelectedBranches ()
{

  if ( m_inputBranches.empty() ) {
    Info("enableSelectedBranches()", "Keeping all branches enabled...");
    return EL::StatusCode::SUCCESS;
  }

  // Firstly, disable all branches
  //
  m_inputNTuple->SetBranchStatus ("*", 0);

  std::vector<std::string> branch_vec;

  // Parse input list, split by comma, and put into a vector
  //
  std::string token;
  std::istringstream ss( m_inputBranches );
  while ( std::getline(ss, token, ',') ) { branch_vec.push_back(token); }

  // Re-enable only the branches we are going to use
  //
  for ( const auto& branch : branch_vec ) {
    if ( !m_isQMisIDBranchIn && branch.compare("QMisIDWeight") == 0 ) { continue; }
    if ( !m_isMMBranchIn && branch.compare("MMWeight") == 0 )         { continue; }
    m_inputNTuple->SetBranchStatus (branch.c_str(), 1);
  }

  return EL::StatusCode::SUCCESS;

}


EL::StatusCode HTopMultilepNTupReprocesser :: setOutputBranches ()
{

  // Clear vector branches from previous event
  //
  ANA_CHECK( this->clearBranches() );

  if ( m_doQMisIDWeighting ) {
      m_QMisIDWeight_out = m_event.get()->weight_QMisID;
  } 
  if ( m_doMMWeighting ) {
      m_MMWeight_out = m_event.get()->weight_MM;
  }

  return EL::StatusCode::SUCCESS;

}

EL::StatusCode HTopMultilepNTupReprocesser :: clearBranches ()
{

  if ( m_doQMisIDWeighting ) { m_QMisIDWeight_out.clear(); }
  if ( m_doMMWeighting )     { m_MMWeight_out.clear(); }

  return EL::StatusCode::SUCCESS;

}

EL::StatusCode HTopMultilepNTupReprocesser ::  readQMisIDRates()
{
    if ( m_QMisIDRates_dir.back() != '/' ) { m_QMisIDRates_dir += "/"; }

    std::string path_AntiT = m_QMisIDRates_dir + m_QMisIDRates_Filename_AntiT;
    std::string path_T     = m_QMisIDRates_dir + m_QMisIDRates_Filename_T;

    TFile *file_AntiT = TFile::Open(path_AntiT.c_str());
    TFile *file_T     = TFile::Open(path_T.c_str());

    if ( !file_AntiT->IsOpen() ) {
	Error("readQMisIDRates()", "Failed to open ROOT file from path: %s . Aborting", path_AntiT.c_str() );
	return EL::StatusCode::FAILURE;
    }
    if ( !file_T->IsOpen() ) {
	Error("readQMisIDRates()", "Failed to open ROOT file from path: %s . Aborting", path_T.c_str() );
	return EL::StatusCode::FAILURE;
    }

    Info("readQMisIDRates()", "Successfully opened ROOT files with QMisID rates from path:\n AntiT --> %s \n T --> %s", path_AntiT.c_str(), path_T.c_str() );

    TH2D *hist_QMisID_AntiT = get_object<TH2D>( *file_AntiT, "Rates" );
    TH2D *hist_QMisID_T     = get_object<TH2D>( *file_T, "Rates" );

    hist_QMisID_AntiT->SetDirectory(0);
    hist_QMisID_T->SetDirectory(0);

    // fill a map for later usage
    //
    m_QMisID_hist_map["AntiT"] = hist_QMisID_AntiT;
    m_QMisID_hist_map["T"]     = hist_QMisID_T;

    return EL::StatusCode::SUCCESS;
}

EL::StatusCode HTopMultilepNTupReprocesser :: calculateQMisIDWeights ()
{
    ANA_CHECK_SET_TYPE (EL::StatusCode);

    // If is not a dileptonic event, return
    //
    if ( m_dilep_type <= 0 ) { return EL::StatusCode::SUCCESS; }

    // If there are no electrons, return
    //
    if ( m_dilep_type == 1 ) { return EL::StatusCode::SUCCESS; }

    std::shared_ptr<leptonObj> el0;
    std::shared_ptr<leptonObj> el1;

    if ( m_dilep_type == 2 ) { // OF events
	el0 = ( m_leptons.at(0).get()->flavour == 11 ) ? m_leptons.at(0) : m_leptons.at(1);
    } else if ( m_dilep_type == 3 ) { // ee events
	el0 = m_leptons.at(0);
	el1 = m_leptons.at(1);
    }

    // Just a precaution...
    //
    if ( el0 && !( fabs(el0.get()->eta) < 2.5 && el0.get()->pt >= 0.0 ) ) { return EL::StatusCode::SUCCESS; }
    if ( el1 && !( fabs(el1.get()->eta) < 2.5 && el1.get()->pt >= 0.0 ) ) { return EL::StatusCode::SUCCESS; }

    float r0(0.0), r0_up(0.0), r0_dn(0.0), r1(0.0), r1_up(0.0), r1_dn(0.0);

    if ( m_useTAntiTRates ) {

	if ( el0 && el1 ) { // ee events
	    if ( el0.get()->tightselected && el1.get()->tightselected ) {
		ANA_CHECK( this->getQMisIDRatesAndError( el0, r0, r0_up, r0_dn, "TIGHT" ) );
		ANA_CHECK( this->getQMisIDRatesAndError( el1, r1, r1_up, r1_dn, "TIGHT" ) );
	    } else {
		ANA_CHECK( this->getQMisIDRatesAndError( el0, r0, r0_up, r0_dn, "ANTI_TIGHT" ) );
		ANA_CHECK( this->getQMisIDRatesAndError( el1, r1, r1_up, r1_dn, "ANTI_TIGHT" ) );
	    }
	} else { // OF events
	    if ( el0.get()->tightselected ) {
		ANA_CHECK( this->getQMisIDRatesAndError( el0, r0, r0_up, r0_dn, "TIGHT" ) );
	    } else {
		ANA_CHECK( this->getQMisIDRatesAndError( el0, r0, r0_up, r0_dn, "ANTI_TIGHT" ) );
	    }
	}

    } else {

	// Look at el0 first...
	//
	if ( el0.get()->tightselected ) {
	    ANA_CHECK( this->getQMisIDRatesAndError( el0, r0, r0_up, r0_dn, "TIGHT" ) );
	} else {
	    ANA_CHECK( this->getQMisIDRatesAndError( el0, r0, r0_up, r0_dn, "ANTI_TIGHT" ) );
	}
	// .. and now at el1 (if any...otherwise r1 weights will be default)
	//
	if ( el1 ) {
	    if (  el1.get()->tightselected ) {
		ANA_CHECK( this->getQMisIDRatesAndError( el1, r1, r1_up, r1_dn, "TIGHT" ) );
	    } else {
		ANA_CHECK( this->getQMisIDRatesAndError( el1, r1, r1_up, r1_dn, "ANTI_TIGHT" ) );
	    }
	}
    }

    if ( m_debug ) {
	Info("calculateQMisIDWeights()","\t r0 = %f ( up = %f, dn = %f )", r0, r0_up, r0_dn );
	Info("calculateQMisIDWeights()","\t r1 = %f ( up = %f, dn = %f )", r1, r1_up, r1_dn );
    }

    // Finally, store the event weight + variations
    //
    if ( !( std::isnan(r0) ) && !( std::isnan(r1) ) && !( std::isinf(r0) ) && !( std::isinf(r1) ) ) {
	m_event.get()->weight_QMisID.at(0) = ( r0 + r1 - 2.0 * r0 * r1 ) / ( 1.0 - r0 - r1 + 2.0 * r0 * r1 ) ;
	m_event.get()->weight_QMisID.at(1) = ( r0_up + r1_up - 2.0 * r0_up * r1_up ) / ( 1.0 - r0_up - r1_up + 2.0 * r0_up * r1_up ) / m_event.get()->weight_QMisID.at(0);
	m_event.get()->weight_QMisID.at(2) = ( r0_dn + r1_dn - 2.0 * r0_dn * r1_dn ) / ( 1.0 - r0_dn - r1_dn + 2.0 * r0_dn * r1_dn ) / m_event.get()->weight_QMisID.at(0);
    } else {
      ++m_count_inf;
    }

    return EL::StatusCode::SUCCESS;
}


EL::StatusCode HTopMultilepNTupReprocesser :: getQMisIDRatesAndError( std::shared_ptr<leptonObj> lep,
		                                     	              float& r, float& r_up, float& r_dn,
								      const std::string& selection )
{

    // Get the 2D histogram from the map
    //
    TH2D* rates_2D(nullptr);
    std::string name_eta(""), name_pt("");

    if ( selection.compare("TIGHT") == 0 ) {
      rates_2D = ( m_QMisID_hist_map.find("T")->second );
      name_eta = "proj_eta_T";
      name_pt  = "proj_pt_T";
    } else if ( selection.compare("ANTI_TIGHT") == 0 ) {
      rates_2D = ( m_QMisID_hist_map.find("AntiT")->second );
      name_eta = "proj_eta_AntiT";
      name_pt  = "proj_pt_AntiT";
    }

    // Make (eta,pT) projections of the 2D histogram with the rates
    //
    TH1D* proj_eta = rates_2D->ProjectionX(name_eta.c_str());
    TH1D* proj_pt  = rates_2D->ProjectionY(name_pt.c_str());

    float this_low_edge(-999.0),this_up_edge(-999.0);

    int eta_bin_nr(-1), pt_bin_nr(-1);

    float eta = lep.get()->etaBE2;
    float pt  = lep.get()->pt;

    // Loop over the projections, and keep track of the bin number where (x,y) is found
    //
    for ( int eta_bin = 0; eta_bin < proj_eta->GetNbinsX()+1; ++eta_bin  ) {

	this_low_edge = proj_eta->GetXaxis()->GetBinLowEdge(eta_bin);
	this_up_edge  = proj_eta->GetXaxis()->GetBinLowEdge(eta_bin+1);

	if ( fabs(eta) >= this_low_edge && fabs(eta) < this_up_edge ) {

	    if ( m_debug ) { Info("getQMisIDRatesAndError()","\t\t eta = %.2f found in %i-th bin", eta, eta_bin ); }

	    eta_bin_nr = proj_eta->GetBin(eta_bin);

	    break;
	}

    }
    for ( int pt_bin = 0; pt_bin < proj_pt->GetNbinsX()+1; ++ pt_bin ) {

	this_low_edge = proj_pt->GetXaxis()->GetBinLowEdge(pt_bin);
	this_up_edge  = proj_pt->GetXaxis()->GetBinLowEdge(pt_bin+1);

	if ( pt/1e3 >= this_low_edge && pt/1e3 < this_up_edge ) {

	    if ( m_debug ) { Info("getQMisIDRatesAndError()","\t\t pT = %.2f found in %i-th bin", pt/1e3, pt_bin ); }

	    pt_bin_nr = proj_pt->GetBin(pt_bin);

	    break;
	}

    }

    if ( m_debug ) { Info("getQMisIDRatesAndError()","\t\t coordinates of efficiency bin = (%i,%i)", eta_bin_nr, pt_bin_nr ); }

    // Now get the NOMINAL rate via global bin number (x,y)

    r = rates_2D->GetBinContent( rates_2D->GetBin( eta_bin_nr, pt_bin_nr ) );

    if ( std::isnan(r) ) {
	Warning("getQMisIDRatesAndError()", "Rate value being read in is nan. Will assign default QMisIDWeight...");
	return EL::StatusCode::SUCCESS;
    }
    if ( std::isinf(r) ) {
	Warning("getQMisIDRatesAndError()", "Rate value being read in is inf. Will assign default QMisIDWeight...");
	return EL::StatusCode::SUCCESS;
    }

    // Get the UP and DOWN variations
    //
    // QUESTION: Why the hell ROOT has GetBinErrorUp and GetBinErrorLow for TH2 ??
    // They seem to give always the same result...
    //
    r_up = r + rates_2D->GetBinErrorUp( rates_2D->GetBin( eta_bin_nr, pt_bin_nr ) );
    r_dn = r - rates_2D->GetBinErrorUp( rates_2D->GetBin( eta_bin_nr, pt_bin_nr ) );
    r_dn = ( r_dn > 0.0 ) ? r_dn : 0.0;

    return EL::StatusCode::SUCCESS;

}

EL::StatusCode HTopMultilepNTupReprocesser :: readRFEfficiencies()
{

  std::string rate_type = ( !m_doMMClosure ) ? "observed" : "expected";

  if ( m_FR_dir.empty() ) { m_FR_dir = m_RR_dir; }

  if ( m_RR_dir.back() != '/' ) { m_RR_dir += "/"; }
  if ( m_FR_dir.back() != '/' ) { m_FR_dir += "/"; }

  // Histogram names - electrons
  //
  std::string histname_el_pt_rr    = "El_ProbePt_Real_Efficiency_"  + rate_type;
  std::string histname_el_eta_rr   = "El_ProbeEta_Real_Efficiency_" + rate_type;
  std::string histname_el_pt_r_T   = "El_ProbePt_Real_T_" + rate_type;
  std::string histname_el_pt_r_L   = "El_ProbePt_Real_L_" + rate_type;
  
  std::string histname_el_pt_fr(""), histname_el_eta_fr(""), histname_el_pt_f_T(""), histname_el_pt_f_L("");
  
  if ( m_useScaledFakeEfficiency && !m_doMMClosure ) {
    histname_el_pt_fr	 = "El_ProbePt_ScaledFake_Efficiency_"  + rate_type; 
    histname_el_eta_fr   = "El_ProbeEta_ScaledFake_Efficiency_" + rate_type;
    histname_el_pt_f_T   = "El_ProbePt_ScaledFake_T_" + rate_type;
    histname_el_pt_f_L   = "El_ProbePt_ScaledFake_L_" + rate_type;
  } else {
    histname_el_pt_fr    = "El_ProbePt_Fake_Efficiency_"  + rate_type;
    histname_el_eta_fr   = "El_ProbeEta_Fake_Efficiency_" + rate_type;
    histname_el_pt_f_T   = "El_ProbePt_Fake_T_" + rate_type;
    histname_el_pt_f_L   = "El_ProbePt_Fake_L_" + rate_type;
  }

  // Histogram names - muons
  //
  std::string histname_mu_pt_rr    = "Mu_ProbePt_Real_Efficiency_"  + rate_type;
  std::string histname_mu_eta_rr   = "Mu_ProbeEta_Real_Efficiency_" + rate_type;
  std::string histname_mu_pt_r_T   = "Mu_ProbePt_Real_T_" + rate_type;
  std::string histname_mu_pt_r_L   = "Mu_ProbePt_Real_L_" + rate_type;

  std::string histname_mu_pt_fr    = "Mu_ProbePt_Fake_Efficiency_"  + rate_type;
  std::string histname_mu_eta_fr   = "Mu_ProbeEta_Fake_Efficiency_" + rate_type;
  std::string histname_mu_pt_f_T   = "Mu_ProbePt_Fake_T_" + rate_type;
  std::string histname_mu_pt_f_L   = "Mu_ProbePt_Fake_L_" + rate_type;

  // Histograms - electrons

  TH1D *hist_el_pt_rr(nullptr);  
  TH1D *hist_el_eta_rr(nullptr);   
  TH1D *hist_el_pt_r_T(nullptr);   
  TH1D *hist_el_pt_r_L(nullptr);   

  TH1D *hist_el_pt_rr_YES_TM(nullptr);
  TH1D *hist_el_pt_rr_NO_TM(nullptr);

  TH1D *hist_el_pt_fr(nullptr);  
  TH1D *hist_el_eta_fr(nullptr); 
  TH1D *hist_el_pt_f_T(nullptr); 
  TH1D *hist_el_pt_f_L(nullptr); 
  
  TH1D *hist_el_pt_fr_YES_TM(nullptr);
  TH1D *hist_el_pt_fr_NO_TM(nullptr);

  // Histograms - muons

  TH1D *hist_mu_pt_rr(nullptr);    
  TH1D *hist_mu_eta_rr(nullptr);   
  TH1D *hist_mu_pt_r_T(nullptr);  
  TH1D *hist_mu_pt_r_L(nullptr);   

  TH1D *hist_mu_pt_rr_YES_TM(nullptr);
  TH1D *hist_mu_pt_rr_NO_TM(nullptr);

  TH1D *hist_mu_pt_fr(nullptr);  
  TH1D *hist_mu_eta_fr(nullptr);  
  TH1D *hist_mu_pt_f_T(nullptr);  
  TH1D *hist_mu_pt_f_L(nullptr);  
  
  TH1D *hist_mu_pt_fr_YES_TM(nullptr);
  TH1D *hist_mu_pt_fr_NO_TM(nullptr);
  
  // 1. 'REAL' efficiency

  Info("readRFEfficiencies()", "REAL efficiency from directory: %s ", m_RR_dir.c_str() );

  std::string path_R_el = m_RR_dir + m_Efficiency_Filename;
  TFile *file_R_el = TFile::Open(path_R_el.c_str());
  
  HTOP_RETURN_CHECK( "HTopMultilepNTupReprocesser::readRFEfficiencies()", file_R_el->IsOpen(), "Failed to open ROOT file" );
  Info("readRFEfficiencies()", "ELECTRON REAL efficiency: %s ", path_R_el.c_str() );
  
  std::string path_R_mu = m_RR_dir + m_Efficiency_Filename;
  TFile *file_R_mu = TFile::Open(path_R_mu.c_str());
  
  HTOP_RETURN_CHECK( "HTopMultilepNTupReprocesser::readRFEfficiencies()", file_R_mu->IsOpen(), "Failed to open ROOT file" );
  Info("readRFEfficiencies()", "MUON REAL efficiency: %s ", path_R_mu.c_str() );

  // Get real efficiency histograms
  //
  hist_el_pt_rr  = get_object<TH1D>( *file_R_el, histname_el_pt_rr );
  if( m_useEtaParametrisation ) hist_el_eta_rr = get_object<TH1D>( *file_R_el, histname_el_eta_rr );
  hist_el_pt_r_T = get_object<TH1D>( *file_R_el, histname_el_pt_r_T );
  hist_el_pt_r_L = get_object<TH1D>( *file_R_el, histname_el_pt_r_L );

  hist_mu_pt_rr  = get_object<TH1D>( *file_R_mu, histname_mu_pt_rr );
  if ( m_useEtaParametrisation ) hist_mu_eta_rr = get_object<TH1D>( *file_R_mu, histname_mu_eta_rr );
  hist_mu_pt_r_T = get_object<TH1D>( *file_R_mu, histname_mu_pt_r_T );
  hist_mu_pt_r_L = get_object<TH1D>( *file_R_mu, histname_mu_pt_r_L );

  // 2. FAKE efficiency

  if ( m_FR_dir.compare(m_RR_dir) != 0 ) {
     Warning("readRFEfficiencies()", "FAKE efficiency is going to be read from %s. Check whether it's really what you want...", m_FR_dir.c_str());
  } else {
     Info("readRFEfficiencies()", "FAKE efficiency from same directory as REAL" );
  }

  std::string path_F_el = m_FR_dir + m_Efficiency_Filename; 
  TFile *file_F_el = TFile::Open(path_F_el.c_str());
  
  HTOP_RETURN_CHECK( "HTopMultilepNTupReprocesser::readRFEfficiencies()", file_F_el->IsOpen(), "Failed to open ROOT file" );
  Info("readRFEfficiencies()", "ELECTRON FAKE efficiency: %s ", path_F_el.c_str() );

  std::string path_F_mu = m_FR_dir + m_Efficiency_Filename;
  TFile *file_F_mu = TFile::Open(path_F_mu.c_str());

  HTOP_RETURN_CHECK( "HTopMultilepNTupReprocesser::readRFEfficiencies()", file_F_mu->IsOpen(), "Failed to open ROOT file" );
  Info("readRFEfficiencies()", "MUON FAKE efficiency: %s ", path_F_mu.c_str() );

  // Get fake efficiency histograms
  //
  hist_el_pt_fr    = get_object<TH1D>( *file_F_el, histname_el_pt_fr );
  if ( m_useEtaParametrisation ) hist_el_eta_fr = get_object<TH1D>( *file_F_el, histname_el_eta_fr );
  hist_el_pt_f_T   = get_object<TH1D>( *file_F_el, histname_el_pt_f_T );
  hist_el_pt_f_L   = get_object<TH1D>( *file_F_el, histname_el_pt_f_L );

  hist_mu_pt_fr    = get_object<TH1D>( *file_F_mu, histname_mu_pt_fr );
  if ( m_useEtaParametrisation ) hist_mu_eta_fr = get_object<TH1D>( *file_F_mu, histname_mu_eta_fr );
  hist_mu_pt_f_T   = get_object<TH1D>( *file_F_mu, histname_mu_pt_f_T );
  hist_mu_pt_f_L   = get_object<TH1D>( *file_F_mu, histname_mu_pt_f_L );

  // ***********************************************************************
  
  if ( m_useTrigMatchingInfo && m_useEtaParametrisation ) {
      Error("readRFEfficiencies()", "As of today, it's not possible to use eta parametrisation when reading trigger-matching-dependent efficiencies. Aborting" );
      return EL::StatusCode::FAILURE;
  }

  if ( m_useTrigMatchingInfo ) {

    if ( m_RRFR_YES_TM_dir.back() != '/' ) { m_RRFR_YES_TM_dir += "/"; }
    if ( m_RRFR_NO_TM_dir.back() != '/' )  { m_RRFR_NO_TM_dir += "/"; }

    Info("readRFEfficiencies()", "REAL/FAKE efficiency (probe TRIGGER-MATCHED) from directory: %s ", m_RRFR_YES_TM_dir.c_str() );
    
    std::string path_YES_TM = m_RRFR_YES_TM_dir + m_Efficiency_Filename;
    TFile *file_YES_TM = TFile::Open(path_YES_TM.c_str());
    HTOP_RETURN_CHECK( "HTopMultilepNTupReprocesser::readRFEfficiencies()", file_YES_TM->IsOpen(), "Failed to open ROOT file" );
    Info("readRFEfficiencies()", "REAL/FAKE efficiency: %s ", path_YES_TM.c_str() );
    
    Info("readRFEfficiencies()", "REAL/FAKE efficiency (probe NOT TRIGGER-MATCHED) from directory: %s ", m_RRFR_NO_TM_dir.c_str() );
    
    std::string path_NO_TM = m_RRFR_NO_TM_dir + m_Efficiency_Filename;
    TFile *file_NO_TM = TFile::Open(path_NO_TM.c_str());
    HTOP_RETURN_CHECK( "HTopMultilepNTupReprocesser::readRFEfficiencies()", file_NO_TM->IsOpen(), "Failed to open ROOT file" );
    Info("readRFEfficiencies()", "REAL/FAKE efficiency: %s ", path_NO_TM.c_str() );

    // Get real efficiency histograms
    //
    hist_el_pt_rr_YES_TM = get_object<TH1D>( *file_YES_TM, histname_el_pt_rr );
    hist_el_pt_rr_NO_TM  = get_object<TH1D>( *file_NO_TM, histname_el_pt_rr );
    
    hist_mu_pt_rr_YES_TM = get_object<TH1D>( *file_YES_TM, histname_mu_pt_rr );
    hist_mu_pt_rr_NO_TM  = get_object<TH1D>( *file_NO_TM, histname_mu_pt_rr );

    // Get fake efficiency histograms
    //
    hist_el_pt_fr_YES_TM = get_object<TH1D>( *file_YES_TM, histname_el_pt_fr );
    hist_el_pt_fr_NO_TM  = get_object<TH1D>( *file_NO_TM, histname_el_pt_fr );

    hist_mu_pt_fr_YES_TM = get_object<TH1D>( *file_YES_TM, histname_mu_pt_fr );
    hist_mu_pt_fr_NO_TM  = get_object<TH1D>( *file_NO_TM, histname_mu_pt_fr );
   
  } 
  // ***********************************************************************

  // Fill a map for later usage
  //
  m_el_hist_map["pt_rr"]   = hist_el_pt_rr;
  m_mu_hist_map["pt_rr"]   = hist_mu_pt_rr;
  m_el_hist_map["pt_fr"]   = hist_el_pt_fr;
  m_mu_hist_map["pt_fr"]   = hist_mu_pt_fr;
  m_el_hist_map["pt_rr_YES_TM"] = hist_el_pt_rr_YES_TM;
  m_el_hist_map["pt_rr_NO_TM"]  = hist_el_pt_rr_NO_TM;
  m_mu_hist_map["pt_rr_YES_TM"] = hist_mu_pt_rr_YES_TM;
  m_mu_hist_map["pt_rr_NO_TM"]  = hist_mu_pt_rr_NO_TM;

  m_el_hist_map["eta_rr"]  = hist_el_eta_rr;
  m_mu_hist_map["eta_rr"]  = hist_mu_eta_rr;
  m_el_hist_map["eta_fr"]  = hist_el_eta_fr;
  m_mu_hist_map["eta_fr"]  = hist_mu_eta_fr;
  m_el_hist_map["pt_fr_YES_TM"] = hist_el_pt_fr_YES_TM;
  m_el_hist_map["pt_fr_NO_TM"]  = hist_el_pt_fr_NO_TM;
  m_mu_hist_map["pt_fr_YES_TM"] = hist_mu_pt_fr_YES_TM;
  m_mu_hist_map["pt_fr_NO_TM"]  = hist_mu_pt_fr_NO_TM;
  
  // pt hist has two different binning for r/f
  // Trigger-match-dependent efficiencies can have different binning wrt non trigger-matching dependent
  //
  m_n_el_bins_pt_rr =  hist_el_pt_rr->GetNbinsX()+1;
  m_n_el_bins_pt_fr =  hist_el_pt_fr->GetNbinsX()+1;
  m_n_mu_bins_pt_rr =  hist_mu_pt_rr->GetNbinsX()+1;
  m_n_mu_bins_pt_fr =  hist_mu_pt_fr->GetNbinsX()+1;

  if ( m_useTrigMatchingInfo ) {
      m_n_el_bins_pt_rr =  hist_el_pt_rr_YES_TM->GetNbinsX()+1;
      m_n_el_bins_pt_fr =  hist_el_pt_fr_YES_TM->GetNbinsX()+1;
      m_n_mu_bins_pt_rr =  hist_mu_pt_rr_YES_TM->GetNbinsX()+1;
      m_n_mu_bins_pt_fr =  hist_mu_pt_fr_YES_TM->GetNbinsX()+1;
  }
  
  // eta hist has same binning for r/f
  //
  if ( m_useEtaParametrisation ) {
      
      m_n_el_bins_eta   =  hist_el_eta_rr->GetNbinsX()+1;
      m_n_mu_bins_eta   =  hist_mu_eta_rr->GetNbinsX()+1;

      // Calculate normalisation factor for (pT * eta) 1D efficiencies case.
      //
      // This factor is the same for eta and pT r/f histograms (it's just Integral(N) / Integral(D) for the efficiency definition ): use pT
      //
      m_el_rr_tot = ( hist_el_pt_r_T->Integral(1,hist_el_pt_r_T->GetNbinsX()+1) ) / ( hist_el_pt_r_L->Integral(1,hist_el_pt_r_L->GetNbinsX()+1) );
      m_el_fr_tot = ( hist_el_pt_f_T->Integral(1,hist_el_pt_f_T->GetNbinsX()+1) ) / ( hist_el_pt_f_L->Integral(1,hist_el_pt_f_L->GetNbinsX()+1) );
      m_mu_rr_tot = ( hist_mu_pt_r_T->Integral(1,hist_mu_pt_r_T->GetNbinsX()+1) ) / ( hist_mu_pt_r_L->Integral(1,hist_mu_pt_r_L->GetNbinsX()+1) );
      m_mu_fr_tot = ( hist_mu_pt_f_T->Integral(1,hist_mu_pt_f_T->GetNbinsX()+1) ) / ( hist_mu_pt_f_L->Integral(1,hist_mu_pt_f_L->GetNbinsX()+1) );
  
  }
  
  std::cout << "\n" << std::endl;
  Info("readRFEfficiencies()", "MUON REAL efficiency - pT histogram name: %s ", histname_mu_pt_rr.c_str() );
  Info("readRFEfficiencies()", "MUON FAKE efficiency - pT histogram name: %s ", histname_mu_pt_fr.c_str() );
  if ( m_useEtaParametrisation ) {
    Info("readRFEfficiencies()", "MUON REAL efficiency - eta histogram name: %s ", histname_mu_eta_rr.c_str() );
    Info("readRFEfficiencies()", "MUON FAKE efficiency - eta histogram name: %s ", histname_mu_eta_fr.c_str() );
  }
  std::cout << "            --------------------------------------------" << std::endl;
  Info("readRFEfficiencies()", "ELECTRON REAL efficiency - pT histogram name: %s ", histname_el_pt_rr.c_str() );
  Info("readRFEfficiencies()", "ELECTRON FAKE efficiency - pT histogram name: %s ", histname_el_pt_fr.c_str() );
  if ( m_useEtaParametrisation ) {
    Info("readRFEfficiencies()", "ELECTRON REAL efficiency - eta histogram name: %s ", histname_el_eta_rr.c_str() );
    Info("readRFEfficiencies()", "ELECTRON FAKE efficiency - eta histogram name: %s ", histname_el_eta_fr.c_str() );
  }
  std::cout << "\n" << std::endl;

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode HTopMultilepNTupReprocesser :: getMMEfficiencyAndError( std::shared_ptr<leptonObj> lep, std::vector<float>& efficiency, const std::string& type )
{

    float error(0.0);

    float pt  = lep.get()->pt/1e3; // Must be in GeV!
    float eta = ( lep.get()->flavour == 13 ) ? lep.get()->eta : lep.get()->etaBE2; 

    float this_low_edge_pt(-1.0), this_up_edge_pt(-1.0);
    float this_low_edge_eta(-999.0), this_up_edge_eta(-999.0);

    int n_bins_pt_fr = ( lep.get()->flavour == 13 ) ? m_n_mu_bins_pt_fr : m_n_el_bins_pt_fr;
    int n_bins_pt_rr = ( lep.get()->flavour == 13 ) ? m_n_mu_bins_pt_rr : m_n_el_bins_pt_rr;
    int n_bins_eta   = ( lep.get()->flavour == 13 ) ? m_n_mu_bins_eta   : m_n_el_bins_eta;

    std::map< std::string, TH1D* > *histograms = ( lep.get()->flavour == 13 ) ? &m_mu_hist_map : &m_el_hist_map;

    TH1D *hist_pt(nullptr);
    TH1D *hist_eta(nullptr);

    //  1) Fake case: choose appropriate histogram
    //
    if ( type.compare("FAKE") == 0 ) {
    
    	if ( m_verbose ) { Info("getMMEfficiencyAndError()", "\tReading fake efficiency..."); }

        hist_pt = ( histograms->find("pt_fr")->second );
	if ( m_useTrigMatchingInfo ) {
	    hist_pt = ( lep.get()->trigmatched ) ? ( histograms->find("pt_fr_YES_TM")->second ) : ( histograms->find("pt_fr_NO_TM")->second );
    	}																					      
	
	if ( m_useEtaParametrisation ) {
	    hist_eta = ( histograms->find("eta_fr")->second );
	}
	
	// Loop over number of pt bins
    	// Do not consider underflow, i.e. 0th bin
    	//
    	for ( int p(1); p <= n_bins_pt_fr; ++p ) {

    	    this_low_edge_pt = hist_pt->GetXaxis()->GetBinLowEdge(p);
    	    this_up_edge_pt  = hist_pt->GetXaxis()->GetBinLowEdge(p+1);

    	    if ( m_verbose ) { Info("getMMEfficiencyAndError()","\t\tpT bin %i : [%.0f,%.0f] GeV", p, this_low_edge_pt, this_up_edge_pt ); }

    	    if ( pt >= this_low_edge_pt && pt < this_up_edge_pt ) {
            
    	        float fr_pt	= hist_pt->GetBinContent(p);
    	        float fr_pt_err = hist_pt->GetBinError(p);
		
      	        if ( m_verbose ) { Info("getMMEfficiencyAndError()", "\t\tLepton pT = %.3f GeV ==> Reading fake efficiency in pT bin [%.0f,%.0f] GeV: fr_pt = %.3f", pt, this_low_edge_pt, this_up_edge_pt, fr_pt ); }

    	        float fr_eta(1.0), fr_eta_err(0.0);

    	        if ( m_useEtaParametrisation ) {
    	            
	    	    // Loop over number of eta bins
    	            // Do not consider underflow, i.e. 0th bin
    	            //
    	            for ( int e(1); e <= n_bins_eta; ++e ) {
    	        	
	    	  	this_low_edge_eta = hist_eta->GetXaxis()->GetBinLowEdge(e);
    	        	this_up_edge_eta  = hist_eta->GetXaxis()->GetBinLowEdge(e+1);

    	        	if ( m_verbose ) { Info("getMMEfficiencyAndError()","\t\t|eta| bin %i : [%.3f,%.3f]", e, this_low_edge_eta, this_up_edge_eta ); }

    	        	if ( fabs(eta) >= this_low_edge_eta && fabs(eta) < this_up_edge_eta ) {

    	        	    fr_eta     = hist_eta->GetBinContent(e);
    	        	    fr_eta_err = hist_eta->GetBinError(e);

    	        	    if ( m_verbose ) {
    	        		Info("getMMEfficiencyAndError()", "\t\tLepton |eta| = %.3f ==> Reading fake efficiency in |eta| bin [%.3f,%.3f]: fr_eta = %.3f", fabs(eta), this_low_edge_eta, this_up_edge_eta, fr_eta );
    	        	    }

    	        	    break;
    	        	}
    	            }
    	        }

    	        // Nominal
    	        //
    	        efficiency.at(0) = fr_pt;
    	        error		 = fr_pt_err;

    	        // UP syst
    	        //
    	        efficiency.at(1) = ( fr_pt + error );

    	        // DN syst
    	        //
    	        if ( fr_pt - error > 0 ) { efficiency.at(2) = ( fr_pt - error ); }
    	        else			 { efficiency.at(2) = 0.0; }

    	        if ( m_useEtaParametrisation ) {
    	            
	    	    float fr_tot = ( lep.get()->flavour == 13 ) ? m_mu_fr_tot : m_el_fr_tot;
    	            if ( m_verbose ) {Info("getMMEfficiencyAndError()", "\t\t norm factor = %.3f", fr_tot ); }
	
	    	    efficiency.at(0) = ( fr_pt * fr_eta ) / fr_tot;
    	            
	    	    // Assuming  fr_pt,fr_eta are independent, this is the error on the product
    	            // ( the constant factor at denominator will be put back later in the def of Efficiency...)
    	            //
    	            error  = sqrt( (fr_eta*fr_pt_err)*(fr_eta*fr_pt_err) + (fr_pt*fr_eta_err)*(fr_pt*fr_eta_err) );
    	            efficiency.at(1) = ( (fr_pt * fr_eta) + error ) / fr_tot;
    	            if ( (fr_pt * fr_eta) - error > 0 ) { efficiency.at(2) = ( (fr_pt * fr_eta) - error ) / fr_tot;}
    	            else				{ efficiency.at(2) = 0.0; }
    	        }

    	        break;
    	    }

    	} // close loop on pT bins: fake case     
  
    }
    //  2) Real case: choose appropriate histogram
    //    
    else if ( type.compare("REAL") == 0 ) {
    
    	if ( m_verbose ) { Info("getMMEfficiencyAndError()", "\tReading real efficiency..."); }

        hist_pt = ( histograms->find("pt_rr")->second );
	if ( m_useTrigMatchingInfo ) {
	    hist_pt = ( lep.get()->trigmatched ) ? ( histograms->find("pt_rr_YES_TM")->second ) : ( histograms->find("pt_rr_NO_TM")->second );
    	}																					      
	
	if ( m_useEtaParametrisation ) {
	    hist_eta = ( histograms->find("eta_rr")->second );
	}

    	// Loop over number of pt bins
    	// Do not consider underflow, i.e. 0th bin
    	//
    	for ( int p(1); p <= n_bins_pt_rr; ++p ) {

    	    this_low_edge_pt = hist_pt->GetXaxis()->GetBinLowEdge(p);
    	    this_up_edge_pt  = hist_pt->GetXaxis()->GetBinLowEdge(p+1);

    	    if ( m_verbose ) { Info("getMMEfficiencyAndError()","\t\tpT bin %i : [%.0f,%.0f] GeV", p, this_low_edge_pt, this_up_edge_pt ); }

    	    if ( pt >= this_low_edge_pt && pt < this_up_edge_pt ) {
            
    	        float rr_pt	= hist_pt->GetBinContent(p);
    	        float rr_pt_err = hist_pt->GetBinError(p);

      	        if ( m_verbose ) { Info("getMMEfficiencyAndError()", "\t\tLepton pT = %.3f GeV ==> Reading real efficiency in pT bin [%.0f,%.0f] GeV: rr_pt = %.3f", pt, this_low_edge_pt, this_up_edge_pt, rr_pt ); }

    	        float rr_eta(1.0), rr_eta_err(0.0);

    	        if ( m_useEtaParametrisation ) {
    	            
	    	    // Loop over number of eta bins
    	            // Do not consider underflow, i.e. 0th bin
    	            //
    	            for ( int e(1); e <= n_bins_eta; ++e ) {
    	        	
	    	  	this_low_edge_eta = ( histograms->find("eta_rr")->second )->GetXaxis()->GetBinLowEdge(e);
    	        	this_up_edge_eta  = ( histograms->find("eta_rr")->second )->GetXaxis()->GetBinLowEdge(e+1);

    	        	if ( m_verbose ) { Info("getMMEfficiencyAndError()","\t\t|eta| bin %i : [%.3f,%.3f]", e, this_low_edge_eta, this_up_edge_eta ); }

    	        	if ( fabs(eta) >= this_low_edge_eta && fabs(eta) < this_up_edge_eta ) {

    	        	    rr_eta     = ( histograms->find("eta_rr")->second )->GetBinContent(e);
    	        	    rr_eta_err = ( histograms->find("eta_rr")->second )->GetBinError(e);

    	        	    if ( m_verbose ) {
    	        		Info("getMMEfficiencyAndError()", "\t\tLepton |eta| = %.3f ==> Reading real efficiency in |eta| bin [%.3f,%.3f]: rr_eta = %.3f", fabs(eta), this_low_edge_eta, this_up_edge_eta, rr_eta );
    	        	    }

    	        	    break;
    	        	}
    	            }
    	        }

    	        // Nominal
    	        //
    	        efficiency.at(0) = rr_pt;
    	        error		 = rr_pt_err;

    	        // UP syst
    	        //
    	        efficiency.at(1) = ( rr_pt + error );

    	        // DN syst
    	        //
    	        if ( rr_pt - error > 0 ) { efficiency.at(2) = ( rr_pt - error ); }
    	        else			 { efficiency.at(2) = 0.0; }

    	        if ( m_useEtaParametrisation ) {
    	            
	    	    float rr_tot = ( lep.get()->flavour == 13 ) ? m_mu_rr_tot : m_el_rr_tot;
    	            if ( m_verbose ) { Info("getMMEfficiencyAndError()", "\t\t norm factor = %.3f", rr_tot ); }
	
	    	    efficiency.at(0) = ( rr_pt * rr_eta ) / rr_tot;
    	            
	    	    // Assuming  rr_pt,rr_eta are independent, this is the error on the product
    	            // ( the constant factor at denominator will be put back later in the def of Efficiency...)
    	            //
    	            error  = sqrt( (rr_eta*rr_pt_err)*(rr_eta*rr_pt_err) + (rr_pt*rr_eta_err)*(rr_pt*rr_eta_err) );
    	            efficiency.at(1) = ( (rr_pt * rr_eta) + error ) / rr_tot;
    	            if ( (rr_pt * rr_eta) - error > 0 ) { efficiency.at(2) = ( (rr_pt * rr_eta) - error ) / rr_tot;}
    	            else				{ efficiency.at(2) = 0.0; }
    	        }

    	        break;
    	    }

    	} // close loop on pT bins: real case     
    
    }    
    
    if ( m_verbose ) {
        if ( type.compare("REAL") == 0 ) { Info("getMMEfficiencyAndError()", "\t\tEffective REAL efficiency ==> r = %.3f ( r_up = %.3f , r_dn = %.3f )", efficiency.at(0), efficiency.at(1), efficiency.at(2) ); }
        if ( type.compare("FAKE") == 0 ) { Info("getMMEfficiencyAndError()", "\t\tEffective FAKE efficiency ==> f = %.3f ( f_up = %.3f , f_dn = %.3f )", efficiency.at(0), efficiency.at(1), efficiency.at(2) ); }
    }    
    
    return EL::StatusCode::SUCCESS;
}

EL::StatusCode HTopMultilepNTupReprocesser :: getMMWeightAndError( std::vector<float>& mm_weight,
								   const std::vector<float>& r0, const std::vector<float>& r1, 
								   const std::vector<float>& f0, const std::vector<float>& f1 )
{

    if ( (r0.at(0) == 0) || (r1.at(0) == 0) || (r0.at(0) <= f0.at(0)) || (r1.at(0) <= f1.at(0)) ) {
	
    	if ( m_debug ) {
    	    Warning("getMMWeightAndError()", "Warning! The Matrix Method cannot be applied because : \nr0 = %.3f , r1 = %.3f, \nf0 = %.3f , f1 = %.3f", r0.at(0), r1.at(0),  f0.at(0), f1.at(0) );
    	    Warning("getMMWeightAndError()", "Setting MM weight = 0 ...");
    	}
        return EL::StatusCode::SUCCESS;
	
    } else {
    	
	// Calculate nominal MM weight
    	//
    	mm_weight.at(0) = matrix_equation( f0.at(0), f1.at(0), r0.at(0), r1.at(0) );
    
    	// Calculate MM weight with systematics
    	//
    	float r0up = ( r0.at(1) > 1.0 ) ? 1.0 :  r0.at(1) ;
    	float r1up = ( r1.at(1) > 1.0 ) ? 1.0 :  r1.at(1) ;
    	float r0dn = r0.at(2);
    	float r1dn = r1.at(2);

    	float f0up = f0.at(1);
    	float f1up = f1.at(1);
    	float f0dn = ( f0.at(2) < 0.0 ) ? 0.0 :  f0.at(2) ;
    	float f1dn = ( f1.at(2) < 0.0 ) ? 0.0 :  f1.at(2) ;

    	// rup syst (save relative weight wrt. nominal)
    	//
    	mm_weight.at(1) = ( matrix_equation( f0.at(0), f1.at(0), r0up, r1up ) / mm_weight.at(0) );
    	
	// fdn syst (save relative weight wrt. nominal)
    	//
    	mm_weight.at(4) = ( matrix_equation( f0dn, f1dn, r0.at(0), r1.at(0) ) / mm_weight.at(0) );

    	if ( (r0dn > f0.at(0)) && (r1dn > f1.at(0)) ) {
    	    
	    // rdn syst (save relative weight wrt. nominal)
    	    //
    	    mm_weight.at(2) = ( matrix_equation( f0.at(0), f1.at(0), r0dn, r1dn ) / mm_weight.at(0) );
    	
	} else {
    	    if ( m_debug ) {
    	        Warning("getMMWeightAndError()", "Warning! Systematic rdn cannot be calculated because : \nr0dn = %.3f , r1dn = %.3f, \nf0 = %.3f , f1 = %.3f", r0dn, r1dn,  f0.at(0), f1.at(0) );
    	    }
    	}

    	if ( (r0.at(0) > f0up) && (r1.at(0) > f1up) ) {
    	    
	    // fup syst (save relative weight wrt. nominal)
    	    //
    	    mm_weight.at(3) = ( matrix_equation( f0up, f1up, r0.at(0), r1.at(0) ) / mm_weight.at(0) );
    	
	} else {
    	    if ( m_debug ) {
    	        Warning("getMMWeightAndError()", "Warning! Systematic fup cannot be calculated because : \nr0dn = %.3f , r1dn = %.3f, \nf0 = %.3f , f1 = %.3f", r0dn, r1dn,  f0.at(0), f1.at(0) );
    	    }
    	}    
     
    }
    
    return EL::StatusCode::SUCCESS;

}

float HTopMultilepNTupReprocesser :: matrix_equation ( const float& f0, const float& f1, const float& r0, const float& r1 )
{

    float w      = 1.0;
    float alpha  = 1.0 / ( (r0-f0) * (r1-f1) );

    if ( m_event.get()->TT ) { 
        if ( m_verbose ) { Info("matrix_equation()", "In region TT:"); }
        w = 1.0 - ( r0 * r1 * ( 1.0 -f0 ) * ( 1.0 - f1 ) * alpha ); 
    } else if ( m_event.get()->TAntiT ) { 
        if ( m_verbose ) { Info("matrix_equation()", "In region TAntiT:"); }
        w = r0 * r1 * f1 * ( 1.0 - f0 ) * alpha;  
    } else if ( m_event.get()->AntiTT ) { 
        if ( m_verbose ) { Info("matrix_equation()", "In region AntiTT:"); }
        w = r0 * r1 * f0 * ( 1.0 - f1 ) * alpha;  
    } else if ( m_event.get()->AntiTAntiT ) { 
        if ( m_verbose ) { Info("matrix_equation()", "In region AntiTAntiT:"); }
        w = -1.0 * r0 * r1 * f0 * f1 * alpha; 
    }

    if ( m_verbose ) { Info("matrix_equation()", "\nr0 = %.3f, r1 = %.3f, f0 = %.3f, f1 = %.3f\nw = %.3f , alpha = %.3f ", r0, r1, f0, f1, w, alpha); }

    // The above formulas are equivalent to the following:
    //
    //float w2 = 1.0;
    //if      ( m_event.get()->TT	   ) { w2 = alpha * ( r0 * f1 * ( (f0 - 1) * (1 - r1) ) + r1 * f0 * ( (r0 - 1) * (1 - f1) ) + f0 * f1 * ( (1 - r0) * (1 - r1) ) ); }
    //else if ( m_event.get()->TAntiT	   ) { w2 = alpha * ( r0 * f1 * ( (1 - f0) * r1 ) + r1 * f0 * ( (1 - r0) * f1 ) + f0 * f1 * ( (r0 - 1) * r1 ) ); }
    //else if ( m_event.get()->AntiTT	   ) { w2 = alpha * ( r0 * f1 * ( (1 - r1) * f0 ) + r1 * f0 * ( (1 - f1) * r0 ) + f0 * f1 * ( (r1 - 1) * r0 ) ); }
    //else if ( m_event.get()->AntiTAntiT  ) { w2 = alpha * ( r0 * f1 * ( -1.0 * f0 * r1 ) + r1 * f0 * ( -1.0 * r0 * f1 ) + f0 * f1 * ( r0 * r1 ) ); }

    return w;
}

EL::StatusCode HTopMultilepNTupReprocesser :: calculateMMWeights()
{
    ANA_CHECK_SET_TYPE (EL::StatusCode);

    // If is not a dileptonic/trileptonic event, return
    //
    if ( m_dilep_type <= 0 && m_trilep_type <= 0 ) { return EL::StatusCode::SUCCESS; }

    std::shared_ptr<leptonObj> lep0 = m_leptons.at(0);
    std::shared_ptr<leptonObj> lep1 = m_leptons.at(1);
    
    // These are the "effective" r/f efficiencies for each lepton, obtained by reading the input r/f histogram(s) 
    
    std::vector<float> r0 = { 1.0, 0.0, 0.0 };
    std::vector<float> r1 = { 1.0, 0.0, 0.0 };
    std::vector<float> f0 = { 1.0, 0.0, 0.0 };
    std::vector<float> f1 = { 1.0, 0.0, 0.0 };
        
    ANA_CHECK( this->getMMEfficiencyAndError( lep0, r0, "REAL" ) );
    ANA_CHECK( this->getMMEfficiencyAndError( lep1, r1, "REAL" ) );
    ANA_CHECK( this->getMMEfficiencyAndError( lep0, f0, "FAKE" ) );
    ANA_CHECK( this->getMMEfficiencyAndError( lep1, f1, "FAKE" ) );
    
    if ( m_debug ) {
        std::cout << "" << std::endl;
	Info("calculateMMWeights()", "Lepton 0 - effective real eff. (nominal, up, dn): " );
	for ( unsigned int idx(0); idx < r0.size(); ++idx ) { std::cout << "r0[" << idx << "] = " << std::setprecision(3) << r0.at(idx) << std::endl; }
	Info("calculateMMWeights()", "Lepton 1 - effective real eff. (nominal, up, dn): " );
	for ( unsigned int idx(0); idx < r1.size(); ++idx ) { std::cout << "r1[" << idx << "] = " << std::setprecision(3) << r1.at(idx) << std::endl; }
	Info("calculateMMWeights()", "Lepton 0 - effective fake eff. (nominal, up, dn): " );
	for ( unsigned int idx(0); idx < f0.size(); ++idx ) { std::cout << "f0[" << idx << "] = " << std::setprecision(3) << f0.at(idx) << std::endl; }
	Info("calculateMMWeights()", "Lepton 1 - effective fake eff. (nominal, up, dn): " );
	for ( unsigned int idx(0); idx < f1.size(); ++idx ) { std::cout << "f1[" << idx << "] = " << std::setprecision(3) << f1.at(idx) << std::endl; }
	std::cout << "" << std::endl;
    }
    
    std::vector<float> mm_weight = { 0.0, 0.0, 0.0, 0.0, 0.0 };
    
    ANA_CHECK( this->getMMWeightAndError( mm_weight, r0, r1, f0, f1 ) );
        
    m_event.get()->weight_MM = mm_weight;
    	
    return EL::StatusCode::SUCCESS;
}
