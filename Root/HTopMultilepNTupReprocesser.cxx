#include <HTopMultilepAnalysis/HTopMultilepNTupReprocesser.h>

// ASG status code check
#include <AsgTools/MessageCheck.h>

// ROOT include(s)
#include "TObjArray.h"

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

  m_RR_dir                = "";
  m_FR_dir                = "";
  m_doMMClosure           = false;
  m_useEtaParemetrisation = false;
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

  m_inputNTuple->SetBranchAddress ("lep_ID_0",   			      &m_lep_ID_0);
  m_inputNTuple->SetBranchAddress ("lep_Pt_0",  			      &m_lep_Pt_0);
  m_inputNTuple->SetBranchAddress ("lep_E_0",   			      &m_lep_E_0);
  m_inputNTuple->SetBranchAddress ("lep_Eta_0",  			      &m_lep_Eta_0);
  m_inputNTuple->SetBranchAddress ("lep_Phi_0",   			      &m_lep_Phi_0);
  m_inputNTuple->SetBranchAddress ("lep_EtaBE2_0",   			      &m_lep_EtaBE2_0);
  m_inputNTuple->SetBranchAddress ("lep_isTightSelected_0",   		      &m_lep_isTightSelected_0);

  m_inputNTuple->SetBranchAddress ("lep_ID_1",   			      &m_lep_ID_1);
  m_inputNTuple->SetBranchAddress ("lep_Pt_1",  			      &m_lep_Pt_1);
  m_inputNTuple->SetBranchAddress ("lep_E_1",   			      &m_lep_E_1);
  m_inputNTuple->SetBranchAddress ("lep_Eta_1",  			      &m_lep_Eta_1);
  m_inputNTuple->SetBranchAddress ("lep_Phi_1",   			      &m_lep_Phi_1);
  m_inputNTuple->SetBranchAddress ("lep_EtaBE2_1",   			      &m_lep_EtaBE2_1);
  m_inputNTuple->SetBranchAddress ("lep_isTightSelected_1",   		      &m_lep_isTightSelected_1);


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

  if ( m_doQMisIDWeighting && m_isQMisIDBranchIn ) {
      Error("initialize()","Trying to add QMisIDWeight branch, which is already in input TTree. Aborting." );
      return EL::StatusCode::FAILURE;  
  }
  if ( m_doMMWeighting && m_isMMBranchIn ) {
      Error("initialize()","Trying to add MMWeight branch, which is already in input TTree. Aborting." );
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

  if ( m_doQMisIDWeighting ) {  
      Info("initialize()","Reading QMisID rates from ROOT file(s)..");
      ANA_CHECK( this->readQMisIDRates() );
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

  lep0.get()->pt          = m_lep_Pt_0;
  lep0.get()->eta         = m_lep_Eta_0;
  lep0.get()->etaBE2      = m_lep_EtaBE2_0;
  lep0.get()->ID          = m_lep_ID_0;
  lep0.get()->flavour     = abs(m_lep_ID_0);
  lep0.get()->charge      = m_lep_ID_0 / fabs(m_lep_ID_0);
  lep0.get()->tightselected = m_lep_isTightSelected_0;

  m_leptons.push_back(lep0);

  auto lep1 = std::make_shared<leptonObj>();

  lep1.get()->pt          = m_lep_Pt_1;
  lep1.get()->eta         = m_lep_Eta_1;
  lep1.get()->etaBE2      = m_lep_EtaBE2_1;
  lep1.get()->ID          = m_lep_ID_1;
  lep1.get()->flavour     = abs(m_lep_ID_1);
  lep1.get()->charge      = m_lep_ID_1 / fabs(m_lep_ID_1);
  lep1.get()->tightselected = m_lep_isTightSelected_1;

  m_leptons.push_back(lep1);

  if ( m_debug ) {
      Info("execute()","lep0 pT = %.2f, lep1 pT = %.2f", lep0.get()->pt/1e3, lep1.get()->pt/1e3 );
      Info("execute()","lep0 flavour = %i, lep1 flavour = %i", lep0.get()->flavour, lep1.get()->flavour );
      Info("execute()","lep0 etaBE2 = %.2f, lep1 etaBE2 = %.2f", lep0.get()->etaBE2, lep1.get()->etaBE2 );
  }

  m_event.get()->isMC   = ( m_mc_channel_number > 0 );
  m_event.get()->dilep  = ( m_dilep_type > 0 );
  m_event.get()->isSS01 = ( m_isSS01 );

  if ( m_debug ) {
      if ( m_doQMisIDWeighting ) {
	  unsigned int idx(0);
	  for ( const auto& itr : m_event.get()->weight_QMisID ) {
	      Info("execute()","\t\tDefault QMisIDWeight[%i] = %f", idx, itr );
	      ++idx;
	  }
      } 
      if ( m_doMMWeighting ) {
	  unsigned int idx(0);
	  for ( const auto& itr : m_event.get()->weight_MM ) {
	      Info("execute()","\t\tDefault MMWeight[%i] = %f", idx, itr );
	      ++idx;
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
	      Info("execute()","\t\tOUT QMisIDWeight[%i] = %f", idx, itr );
	      ++idx;
	  }
      } 
      if ( m_doMMWeighting ) {
	  unsigned int idx(0);
	  for ( const auto& itr : m_MMWeight_out ) {
	      Info("execute()","\t\tOUT MMWeight[%i] = %f", idx, itr );
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
	m_event.get()->weight_QMisID.at(1) = ( r0_up + r1_up - 2.0 * r0_up * r1_up ) / ( 1.0 - r0_up - r1_up + 2.0 * r0_up * r1_up );
	m_event.get()->weight_QMisID.at(2) = ( r0_dn + r1_dn - 2.0 * r0_dn * r1_dn ) / ( 1.0 - r0_dn - r1_dn + 2.0 * r0_dn * r1_dn );
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
