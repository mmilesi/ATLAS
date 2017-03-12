// package include(s):
#include <HTopMultilepAnalysis/HTopMultilepNTupReprocesser.h>
#include "HTopMultilepAnalysis/tools/HTopReturnCheck.h"

// ASG status code check
#include <AsgTools/MessageCheck.h>

// ROOT include(s)
#include "TObjArray.h"

// C++ include(s)
#include <iomanip>
#include <memory>

using namespace NTupReprocesser;

// this is needed to distribute the algorithm to the workers
ClassImp(HTopMultilepNTupReprocesser)

HTopMultilepNTupReprocesser :: HTopMultilepNTupReprocesser(std::string className) :
    Algorithm(className),
    m_inputNTuple(nullptr),
    m_outputNTuple(nullptr),
    m_isQMisIDBranchIn(false),
    m_isMMBranchIn(false),
    m_doQMisIDWeighting(false),
    m_doMMWeighting(false)
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().

  Info("HTopMultilepNTupReprocesser()", "Calling constructor");

  m_inputBranches        = "";

  m_outputNTupName       = "physics";
  m_outputNTupStreamName = "output";

  m_weightToCalc         = "";

  m_QMisIDRates_dir            = "";
  m_QMisIDRates_Filename_T     = "";
  m_QMisIDRates_Filename_AntiT = "";
  m_QMisIDRates_Histname_T = m_QMisIDRates_Histname_AntiT = "Rates";
  m_useTAntiTRates             = false;

  m_REFF_dir                = "";
  m_FEFF_dir                = "";
  m_EFF_YES_TM_dir          = "";
  m_EFF_NO_TM_dir           = "";
  m_Efficiency_Filename     = "";
  m_doMMClosure             = false;
  m_useEtaParametrisation   = false;
  m_useTrigMatchingInfo     = false;
  m_useScaledFakeEfficiency = false;
  m_useCutBasedLep          = false;
  m_useTEfficiency          = false;

  m_systematics_list        = "Stat";

  m_correlatedMMWeights     = false;
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

  job = job;

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepNTupReprocesser :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.

  ANA_CHECK_SET_TYPE (EL::StatusCode);

  Info("histInitialize()", "Calling histInitialize");

  // Parse input weight list, split by comma, and put into a vector

  std::vector<std::string> weights;
  ANA_CHECK( this->tokenize( ',', weights, m_weightToCalc ) );

  if ( std::find( weights.begin(), weights.end(), "QMisID" ) != weights.end() ) { m_doQMisIDWeighting = true; }
  if ( std::find( weights.begin(), weights.end(), "MM" ) != weights.end() )	{ m_doMMWeighting = true;     }

  if ( m_doQMisIDWeighting ) {
      Info("initialize()","Reading QMisID rates from ROOT file(s)..");
      ANA_CHECK( this->readQMisIDRates() );
  }

  if ( m_doMMWeighting ) {
      Info("initialize()","Reading MM efficiencies from ROOT file(s)..");
      ANA_CHECK( this->readRFEfficiencies() );
  }

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

  m_inputNTuple = wk()->tree();

  // Check content of input tree

  TObjArray* branches = m_inputNTuple->GetListOfBranches();
  int nbranches = branches->GetEntriesFast();
  for ( int idx(0); idx < nbranches; ++idx ) {
      std::string this_branch(branches->At(idx)->GetName());
      if ( this_branch.find("QMisIDWeight") != std::string::npos ) {
	  m_isQMisIDBranchIn = true;
	  break;
      }
  }
  for ( int idx(0); idx < nbranches; ++idx ) {
      std::string this_branch(branches->At(idx)->GetName());
      if ( this_branch.find("MMWeight") != std::string::npos ) {
	  m_isMMBranchIn = true;
	  break;
      }
  }

  ANA_CHECK( this->enableSelectedBranches() );

  // Connect the branches of the input tree to the algorithm members

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

  m_inputNTuple->SetBranchAddress ("is_TMVA_TMVA",			      &m_is_TMVA_TMVA);
  m_inputNTuple->SetBranchAddress ("is_TMVA_AntiTMVA",			      &m_is_TMVA_AntiTMVA);
  m_inputNTuple->SetBranchAddress ("is_AntiTMVA_TMVA",			      &m_is_AntiTMVA_TMVA);
  m_inputNTuple->SetBranchAddress ("is_AntiTMVA_AntiTMVA",		      &m_is_AntiTMVA_AntiTMVA);

  m_inputNTuple->SetBranchAddress ("lep_ID_0",   			      &m_lep_ID_0);
  m_inputNTuple->SetBranchAddress ("lep_Pt_0",  			      &m_lep_Pt_0);
  m_inputNTuple->SetBranchAddress ("lep_E_0",   			      &m_lep_E_0);
  m_inputNTuple->SetBranchAddress ("lep_Eta_0",  			      &m_lep_Eta_0);
  m_inputNTuple->SetBranchAddress ("lep_Phi_0",   			      &m_lep_Phi_0);
  m_inputNTuple->SetBranchAddress ("lep_EtaBE2_0",   			      &m_lep_EtaBE2_0);
  m_inputNTuple->SetBranchAddress ("lep_isTightSelected_0",   		      &m_lep_isTightSelected_0);
  m_inputNTuple->SetBranchAddress ("lep_isTightSelectedMVA_0", 		      &m_lep_isTightSelectedMVA_0);
  m_inputNTuple->SetBranchAddress ("lep_isTrigMatch_0",   		      &m_lep_isTrigMatch_0);

  m_inputNTuple->SetBranchAddress ("lep_ID_1",   			      &m_lep_ID_1);
  m_inputNTuple->SetBranchAddress ("lep_Pt_1",  			      &m_lep_Pt_1);
  m_inputNTuple->SetBranchAddress ("lep_E_1",   			      &m_lep_E_1);
  m_inputNTuple->SetBranchAddress ("lep_Eta_1",  			      &m_lep_Eta_1);
  m_inputNTuple->SetBranchAddress ("lep_Phi_1",   			      &m_lep_Phi_1);
  m_inputNTuple->SetBranchAddress ("lep_EtaBE2_1",   			      &m_lep_EtaBE2_1);
  m_inputNTuple->SetBranchAddress ("lep_isTightSelected_1",   		      &m_lep_isTightSelected_1);
  m_inputNTuple->SetBranchAddress ("lep_isTightSelectedMVA_1", 		      &m_lep_isTightSelectedMVA_1);
  m_inputNTuple->SetBranchAddress ("lep_isTrigMatch_1",   		      &m_lep_isTrigMatch_1);

  if ( m_isQMisIDBranchIn ) {
      m_inputNTuple->SetBranchAddress ("QMisIDWeight",     &m_QMisIDWeight_NOMINAL_in);
      m_inputNTuple->SetBranchAddress ("QMisIDWeight_up",  &m_QMisIDWeight_UP_in);
      m_inputNTuple->SetBranchAddress ("QMisIDWeight_dn",  &m_QMisIDWeight_DN_in);
  }

  if ( m_isMMBranchIn ) {

      m_inputNTuple->SetBranchAddress ("MMWeight", &m_MMWeight_NOMINAL_in);

      for ( const auto& sys : m_systematics ) {

	  if ( sys.find("Nominal") != std::string::npos ) { continue; }

	  m_MMWeight_in[sys] = std::vector<float>(2); // need up and dn for each systematic

	  std::string branchname_0 = "MMWeight_" + sys + "_up";
	  std::string branchname_1 = "MMWeight_" + sys + "_dn";

	  m_inputNTuple->SetBranchAddress( (branchname_0).c_str(), &(m_MMWeight_in[sys].at(0)) );
	  m_inputNTuple->SetBranchAddress( (branchname_1).c_str(), &(m_MMWeight_in[sys].at(1)) );

      }

  }

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode  HTopMultilepNTupReprocesser :: tokenize ( char separator, std::vector<std::string>& vec_tokens, const std::string& list ) {

  std::string token;
  std::istringstream ss( list );
  while ( std::getline(ss, token, separator) ) { vec_tokens.push_back(token); }

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

  //ANA_CHECK_SET_TYPE (EL::StatusCode);

  Info("initialize()", "Initialising HTopMultilepNTupReprocesser...");

  m_outputNTuple = EL::getNTupleSvc (wk(), m_outputNTupStreamName);

  // Set new branches for output TTree

  if ( m_doQMisIDWeighting ) {
      m_outputNTuple->tree()->Branch("QMisIDWeight",     &m_QMisIDWeight_NOMINAL_out,    "QMisIDWeight/F");
      m_outputNTuple->tree()->Branch("QMisIDWeight_up",  &m_QMisIDWeight_UP_out,         "QMisIDWeight_up/F");
      m_outputNTuple->tree()->Branch("QMisIDWeight_dn",  &m_QMisIDWeight_DN_out,         "QMisIDWeight_dn/F");
  }

  if ( m_doMMWeighting ) {

      // Set output branch for the nominal weight

      m_outputNTuple->tree()->Branch("MMWeight", &m_MMWeight_NOMINAL_out, "MMWeight/F");

      // Initialise the map containing the variations of MM weights for each input systematics.

      for ( const auto& sys : m_systematics ) {
	  m_MMWeight_out[sys] = std::vector<float>(2,1.0);
      }

      // Set output branches for the variations of MM weight for each systematic.

      for ( const auto& sys : m_systematics ) {

	  if ( sys.find("Nominal") != std::string::npos ) { continue; }

	  std::string branchname_0 = "MMWeight_" + sys + "_up";
	  std::string branchname_1 = "MMWeight_" + sys + "_dn";

	  m_outputNTuple->tree()->Branch( (branchname_0).c_str(), &(m_MMWeight_out[sys].at(0)) );
	  m_outputNTuple->tree()->Branch( (branchname_1).c_str(), &(m_MMWeight_out[sys].at(1)) );

      }

  }

  // ---------------------------------------------------------------------------------------------------------------

  m_outputNTuple->tree()->SetName( m_outputNTupName.c_str() );

  // ---------------------------------------------------------------------------------------------------------------

  // Copy input TTree weight to output TTree

  m_outputNTuple->tree()->SetWeight( m_inputNTuple->GetWeight() );

  // ---------------------------------------------------------------------------------------------------------------

  // Initialise counter for input TTree entries

  m_numEntry = -1;

  // Initialise counter for events where inf/nan is read

  m_count_inf = 0;

  m_effectiveTotEntries = m_inputNTuple->GetEntries();

  unsigned int maxEvents = static_cast<int>( wk()->metaData()->castDouble("nc_EventLoop_MaxEvents") );

  if ( maxEvents > 0 ) {
    m_effectiveTotEntries = maxEvents;
  }

  Info("initialize()", "Name of input TTree : %s", m_inputNTuple->GetName() );
  Info("initialize()", "Total events to run on: %u", m_effectiveTotEntries );

  // ---------------------------------------------------------------------------------------------------------------

  Info("initialize()", "All good!");

  return EL::StatusCode::SUCCESS;
}


void HTopMultilepNTupReprocesser :: printWeights ( const std::string& in_out ) {

    if ( in_out.compare("IN" ) == 0 ) {

	if ( m_doQMisIDWeighting ) {
	    if ( !m_isQMisIDBranchIn ) {
		Info("printWeights()","\t\tDefault QMisIDWeight = %.3f",                                        m_QMisIDWeight_NOMINAL_out );
		Info("printWeights()","\t\tDefault QMisIDWeight (up) * nominal = %.3f ( not rescaled = %.3f )", m_QMisIDWeight_NOMINAL_out * m_QMisIDWeight_UP_out, m_QMisIDWeight_UP_out );
		Info("printWeights()","\t\tDefault QMisIDWeight (dn) * nominal = %.3f ( not rescaled = %.3f )", m_QMisIDWeight_NOMINAL_out * m_QMisIDWeight_DN_out, m_QMisIDWeight_DN_out );
	    } else {
		Info("printWeights()","\t\tIN QMisIDWeight = %.3f",                                        m_QMisIDWeight_NOMINAL_in );
		Info("printWeights()","\t\tIN QMisIDWeight (up) * nominal = %.3f ( not rescaled = %.3f )", m_QMisIDWeight_NOMINAL_in * m_QMisIDWeight_UP_in, m_QMisIDWeight_UP_in  );
		Info("printWeights()","\t\tIN QMisIDWeight (dn) * nominal = %.3f ( not rescaled = %.3f )", m_QMisIDWeight_NOMINAL_in * m_QMisIDWeight_DN_in, m_QMisIDWeight_DN_in  );
	    }
	}

	if ( m_doMMWeighting ) {
	    if ( !m_isMMBranchIn ) {
		Info("printWeights()","\t\tNominal Default MMWeight = %.3f", m_MMWeight_NOMINAL_out );
		for ( const auto& sys : m_systematics ) {
		    if ( sys.find("Nominal") != std::string::npos ) { continue; }
		    Info("printWeights()","\t\tSys: %s ==> Default MMWeight (up) * nominal = %.3f ( not rescaled = %.3f )", sys.c_str(), m_MMWeight_NOMINAL_out * m_MMWeight_out[sys].at(0), m_MMWeight_out[sys].at(0) );
		    Info("printWeights()","\t\tSys: %s ==> Default MMWeight (dn) * nominal = %.3f ( not rescaled = %.3f )", sys.c_str(), m_MMWeight_NOMINAL_out * m_MMWeight_out[sys].at(1), m_MMWeight_out[sys].at(1) );
		}
	    } else {
		Info("printWeights()","\t\tNominal IN MMWeight = %.3f", m_MMWeight_NOMINAL_in );
		for ( const auto& sys : m_systematics ) {
		    if ( sys.find("Nominal") != std::string::npos ) { continue; }
		    Info("printWeights()","\t\tSys: %s ==> IN MMWeight (up) * nominal = %.3f ( not rescaled = %.3f )", sys.c_str(), m_MMWeight_NOMINAL_in * m_MMWeight_in[sys].at(0), m_MMWeight_in[sys].at(0) );
		    Info("printWeights()","\t\tSys: %s ==> IN MMWeight (dn) * nominal = %.3f ( not rescaled = %.3f )", sys.c_str(), m_MMWeight_NOMINAL_in * m_MMWeight_in[sys].at(1), m_MMWeight_in[sys].at(1) );
		}
	    }
	}

    } else if ( in_out.compare("OUT" ) == 0 ) {

	if ( m_doQMisIDWeighting ) {
	    Info("printWeights()","\t\tOUT QMisIDWeight = %.3f",      m_QMisIDWeight_NOMINAL_out );
	    Info("printWeights()","\t\tOUT QMisIDWeight (up) * nominal = %.3f ( not rescaled = %.3f )", m_QMisIDWeight_NOMINAL_out * m_QMisIDWeight_UP_out, m_QMisIDWeight_UP_out );
	    Info("printWeights()","\t\tOUT QMisIDWeight (dn) * nominal = %.3f ( not rescaled = %.3f )", m_QMisIDWeight_NOMINAL_out * m_QMisIDWeight_DN_out, m_QMisIDWeight_DN_out );
	}

	if ( m_doMMWeighting ) {
	    Info("printWeights()","\t\tNominal OUT MMWeight = %.3f", m_MMWeight_NOMINAL_out );
	    for ( const auto& sys : m_systematics ) {
		if ( sys.find("Nominal") != std::string::npos ) { continue; }
		Info("printWeights()","\t\tSys: %s ==> OUT MMWeight (up) * nominal = %.3f ( not rescaled = %.3f )", sys.c_str(), m_MMWeight_NOMINAL_out * m_MMWeight_out[sys].at(0), m_MMWeight_out[sys].at(0) );
		Info("printWeights()","\t\tSys: %s ==> OUT MMWeight (dn) * nominal = %.3f ( not rescaled = %.3f )", sys.c_str(), m_MMWeight_NOMINAL_out * m_MMWeight_out[sys].at(1), m_MMWeight_out[sys].at(1) );
	    }
	}

    }

}


EL::StatusCode HTopMultilepNTupReprocesser :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  ++m_numEntry;

  ANA_CHECK_SET_TYPE (EL::StatusCode);

  if ( m_numEntry == 0 ) { Info("execute()", "Processing input TTree : %s\n", m_inputNTuple->GetName() ); }

  m_inputNTuple->GetEntry( wk()->treeEntry() );

  if ( m_debug ) {
    std::cout << "" << std::endl;
    Info("execute()", "===> Entry %u - EventNumber = %u ", static_cast<uint32_t>(m_numEntry), static_cast<uint32_t>(m_EventNumber) );
  }

  if ( m_numEntry >= 10000 && m_numEntry % 10000 == 0 ) {
    std::cout << "Processed " << std::setprecision(3) << ( (float) m_numEntry / m_effectiveTotEntries ) * 1e2 << " % of total entries" << std::endl;
  }

  // ------------------------------------------------------------------------

  // Need to ensure all weight branches are reset to their default values
  // before getting new values for the event.

  ANA_CHECK( this->resetDefaultWeights() );

  // ------------------------------------------------------------------------

  // This call is crucial, otherwise you'll get no entries in the output tree!

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
  lep0.get()->tightselected = ( m_useCutBasedLep ) ? m_lep_isTightSelected_0 : m_lep_isTightSelectedMVA_0;
  lep0.get()->trigmatched   = m_lep_isTrigMatch_0; // SLT matching

  m_leptons.push_back(lep0);

  auto lep1 = std::make_shared<leptonObj>();

  lep1.get()->pt            = m_lep_Pt_1;
  lep1.get()->eta           = m_lep_Eta_1;
  lep1.get()->etaBE2        = m_lep_EtaBE2_1;
  lep1.get()->ID            = m_lep_ID_1;
  lep1.get()->flavour       = abs(m_lep_ID_1);
  lep1.get()->charge        = m_lep_ID_1 / fabs(m_lep_ID_1);
  lep1.get()->tightselected = ( m_useCutBasedLep ) ? m_lep_isTightSelected_1 : m_lep_isTightSelectedMVA_1;
  lep1.get()->trigmatched   = m_lep_isTrigMatch_1; // SLT matching

  m_leptons.push_back(lep1);

  m_event.get()->isMC        = ( m_mc_channel_number > 0 );
  m_event.get()->dilep_type  = m_dilep_type;
  m_event.get()->trilep_type = m_trilep_type;
  m_event.get()->isSS01      = m_isSS01;

  m_event.get()->TT         = ( m_useCutBasedLep ) ? m_is_T_T : m_is_TMVA_TMVA;
  m_event.get()->TAntiT     = ( m_useCutBasedLep ) ? m_is_T_AntiT : m_is_TMVA_AntiTMVA;
  m_event.get()->AntiTT     = ( m_useCutBasedLep ) ? m_is_AntiT_T : m_is_AntiTMVA_TMVA;
  m_event.get()->AntiTAntiT = ( m_useCutBasedLep ) ? m_is_AntiT_AntiT : m_is_AntiTMVA_AntiTMVA;

  // ------------------------------------------------------------------------

  if ( m_debug ) { this->printWeights( "IN" ); }

  // ------------------------------------------------------------------------

  // std::cout << "\nForcing job dump...\n" << std::endl;
  //  return EL::StatusCode::FAILURE;

  if ( m_doQMisIDWeighting ) { ANA_CHECK( this->calculateQMisIDWeights () ); }

  if ( m_doMMWeighting ) {

      // Apply MM only when pT(l0), pT(l1) > 10 GeV (--> it's the lowest pT for which r/f efficiencies are measured)

      if ( lep0.get()->pt >= 10e3 && lep1.get()->pt >= 10e3 ) {

	  if ( m_debug ) {
	      Info("execute()","Main properties of event:");
	      Info("execute()","lep[0]: pT = %.2f\t etaBE2 = %.2f\t eta = %.2f\t flavour = %i\t tight? %i\t trigmatched (SLT)? %i", lep0.get()->pt/1e3, lep0.get()->etaBE2, lep0.get()->eta, lep0.get()->flavour, lep0.get()->tightselected, lep0.get()->trigmatched );
	      Info("execute()","lep[1]: pT = %.2f\t etaBE2 = %.2f\t eta = %.2f\t flavour = %i\t tight? %i\t trigmatched (SLT)? %i", lep1.get()->pt/1e3, lep1.get()->etaBE2, lep1.get()->eta, lep1.get()->flavour, lep1.get()->tightselected, lep1.get()->trigmatched );
	      Info("execute()","dilep_type = %i, TT ? %i, TAntiT ? %i, AntiTT ? %i, AntiTAntiT ? %i\n", m_event.get()->dilep_type, m_event.get()->TT, m_event.get()->TAntiT, m_event.get()->AntiTT, m_event.get()->AntiTAntiT );
	  }

	  // Get a full set of MM weights with systematic variations
	  // (Make sure no other variables than pT are considered for systematic variations)

	  for ( const auto& sys : m_systematics ) {
	      if ( sys.find("Pt") == std::string::npos ) { continue; }
	      m_this_syst = sys;
	      ANA_CHECK( this->calculateMMWeights () );
	  }

      }

  }

  // ------------------------------------------------------------------------

  if ( m_debug ) { this->printWeights( "OUT" ); }

  // ------------------------------------------------------------------------

  ANA_CHECK( this->clearBranches() );

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

  m_inputNTuple->SetBranchStatus ("*", 0);

  // Parse input list, split by comma, and put into a vector

  std::vector<std::string> branch_vec;
  ANA_CHECK( this->tokenize( ',', branch_vec, m_inputBranches ) );

  // Re-enable only the branches we are going to use

  Info("enableSelectedBranches()", "Activating branches:\n");
  for ( const auto& branch : branch_vec ) {

    if ( !m_isQMisIDBranchIn && branch.find("QMisIDWeight") != std::string::npos ) { continue; }
    if ( !m_isMMBranchIn && branch.find("MMWeight") != std::string::npos )         { continue; }

    if ( m_debug ) { std::cout << "SetBranchStatus(" << branch << ", 1)" << std::endl; }

    m_inputNTuple->SetBranchStatus (branch.c_str(), 1);

  }

  return EL::StatusCode::SUCCESS;

}

EL::StatusCode HTopMultilepNTupReprocesser :: resetDefaultWeights ()
{

  // Need to reset to default values the output weights!
  // (in a previous implementation these branches were also members of the "eventObj" class, which
  // would get instantiated from scratch for every event, and re-initialise its data members through constructor call)

  m_QMisIDWeight_NOMINAL_out = 1.0;
  m_QMisIDWeight_UP_out      = 1.0;
  m_QMisIDWeight_DN_out      = 1.0;

  m_MMWeight_NOMINAL_out = 1.0;
  for ( auto& weight_sys : m_MMWeight_out ) {
    std::fill( weight_sys.second.begin(), weight_sys.second.end(), 1.0 );
  }

  return EL::StatusCode::SUCCESS;

}

EL::StatusCode HTopMultilepNTupReprocesser :: clearBranches ()
{

  // If you store container branches in the output tree (eg., std::vector),
  // remember to clear them after every event!

  return EL::StatusCode::SUCCESS;

}


EL::StatusCode HTopMultilepNTupReprocesser ::  readQMisIDRates()
{
    if ( m_QMisIDRates_dir.back() != '/' ) { m_QMisIDRates_dir += "/"; }

    std::string path_AntiT = m_QMisIDRates_dir + m_QMisIDRates_Filename_AntiT;
    std::string path_T     = m_QMisIDRates_dir + m_QMisIDRates_Filename_T;

    TFile *file_AntiT = TFile::Open(path_AntiT.c_str());
    TFile *file_T     = TFile::Open(path_T.c_str());

    HTOP_RETURN_CHECK( "HTopMultilepNTupReprocesser::readQMisIDRates()", file_AntiT->IsOpen(), "Failed to open ROOT file" );
    HTOP_RETURN_CHECK( "HTopMultilepNTupReprocesser::readQMisIDRates()", file_T->IsOpen(), "Failed to open ROOT file" );

    Info("readQMisIDRates()", "Successfully opened ROOT files with QMisID rates from path:\n AntiT --> %s \n T --> %s", path_AntiT.c_str(), path_T.c_str() );

    TH2D *hist_QMisID_AntiT = get_object<TH2D>( *file_AntiT, m_QMisIDRates_Histname_AntiT.c_str() );
    TH2D *hist_QMisID_T     = get_object<TH2D>( *file_T, m_QMisIDRates_Histname_T.c_str() );

    hist_QMisID_AntiT->SetDirectory(0);
    hist_QMisID_T->SetDirectory(0);

    // Fill a map for later usage

    m_QMisID_hist_map["AntiT"] = hist_QMisID_AntiT;
    m_QMisID_hist_map["T"]     = hist_QMisID_T;

    return EL::StatusCode::SUCCESS;
}

EL::StatusCode HTopMultilepNTupReprocesser :: calculateQMisIDWeights ()
{
    ANA_CHECK_SET_TYPE (EL::StatusCode);

    // If is not a dileptonic event, return

    if ( m_event.get()->dilep_type <= 0 ) { return EL::StatusCode::SUCCESS; }

    // If there are no electrons, return

    if ( m_event.get()->dilep_type == 1 ) { return EL::StatusCode::SUCCESS; }

    std::shared_ptr<leptonObj> el0;
    std::shared_ptr<leptonObj> el1;

    if ( m_event.get()->dilep_type == 2 ) { // OF events
	el0 = ( m_leptons.at(0).get()->flavour == 11 ) ? m_leptons.at(0) : m_leptons.at(1);
    } else if ( m_event.get()->dilep_type == 3 ) { // ee events
	el0 = m_leptons.at(0);
	el1 = m_leptons.at(1);
    }

    // Just a precaution...

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

	if ( el0.get()->tightselected ) {
	    ANA_CHECK( this->getQMisIDRatesAndError( el0, r0, r0_up, r0_dn, "TIGHT" ) );
	} else {
	    ANA_CHECK( this->getQMisIDRatesAndError( el0, r0, r0_up, r0_dn, "ANTI_TIGHT" ) );
	}

	// .. and now at el1 (if any...otherwise r1 weights will be default)

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

    // Finally, store the event weight + (relative) variations

    if ( !( std::isnan(r0) ) && !( std::isnan(r1) ) && !( std::isinf(r0) ) && !( std::isinf(r1) ) ) {

        float nominal = ( r0 + r1 - 2.0 * r0 * r1 ) / ( 1.0 - r0 - r1 + 2.0 * r0 * r1 );
        float up      = ( r0_up + r1_up - 2.0 * r0_up * r1_up ) / ( 1.0 - r0_up - r1_up + 2.0 * r0_up * r1_up );
        float dn      = ( r0_dn + r1_dn - 2.0 * r0_dn * r1_dn ) / ( 1.0 - r0_dn - r1_dn + 2.0 * r0_dn * r1_dn );

	m_QMisIDWeight_NOMINAL_out    = nominal;
	m_QMisIDWeight_UP_out = ( !std::isnan(up/nominal) && !std::isinf(up/nominal) ) ? up/nominal : 0.0;
	m_QMisIDWeight_DN_out = ( !std::isnan(dn/nominal) && !std::isinf(dn/nominal) ) ? dn/nominal : 0.0;

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

    TH1D* proj_eta = rates_2D->ProjectionX(name_eta.c_str());
    TH1D* proj_pt  = rates_2D->ProjectionY(name_pt.c_str());

    float this_low_edge(-999.0),this_up_edge(-999.0);

    int eta_bin_nr(-1), pt_bin_nr(-1);

    float eta = lep.get()->etaBE2;
    float pt  = lep.get()->pt;

    // Loop over the projections, and keep track of the bin number where (x,y) is found

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

    // Now get the NOMINAL rate from the TH2 map via global bin number (x,y)

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

    r_up = r + rates_2D->GetBinErrorUp( rates_2D->GetBin( eta_bin_nr, pt_bin_nr ) );
    r_dn = r - rates_2D->GetBinErrorUp( rates_2D->GetBin( eta_bin_nr, pt_bin_nr ) );
    r_dn = ( r_dn > 0.0 ) ? r_dn : 0.0;

    return EL::StatusCode::SUCCESS;

}


std::string HTopMultilepNTupReprocesser :: str_replace( const std::string& input_str, const std::string& old_substr, const std::string& new_substr )
{
  size_t start_pos = input_str.find(old_substr);

  if ( start_pos == std::string::npos ) {
    Warning("str_replace()", "Substring: %s  not found in input string: %s. Returning input string.", old_substr.c_str(), input_str.c_str() );
    return input_str;
  }

  std::string output_str = input_str;

  return output_str.replace( start_pos, old_substr.length(),new_substr );
}


EL::StatusCode HTopMultilepNTupReprocesser :: readRFEfficiencies()
{

  std::string process = ( !m_doMMClosure ) ? "observed_sub" : "expectedbkg";

  if ( m_FEFF_dir.empty() ) { m_FEFF_dir = m_REFF_dir; }

  if ( m_REFF_dir.back() != '/' ) { m_REFF_dir += "/"; }
  if ( m_FEFF_dir.back() != '/' ) { m_FEFF_dir += "/"; }

  if ( m_useTrigMatchingInfo && m_useEtaParametrisation ) {
      Error("readRFEfficiencies()", "As of today, the only supported parametrisation when reading trigger-dependent efficiencies is pT. Check your job configuration and retry. Aborting" );
      return EL::StatusCode::FAILURE;
  }

  TFile *file_YES_TM(nullptr), *file_NO_TM(nullptr);

  if ( m_useTrigMatchingInfo ) {

      if ( m_EFF_YES_TM_dir.back() != '/' ) { m_EFF_YES_TM_dir += "/"; }
      if ( m_EFF_NO_TM_dir.back() != '/' )  { m_EFF_NO_TM_dir += "/"; }

      Info("readRFEfficiencies()", "REAL/FAKE efficiency (probe SLT-TRIGGER-MATCHED) from directory: %s ", m_EFF_YES_TM_dir.c_str() );

      std::string path_YES_TM = m_EFF_YES_TM_dir + m_Efficiency_Filename;
      file_YES_TM = TFile::Open(path_YES_TM.c_str());
      HTOP_RETURN_CHECK( "HTopMultilepNTupReprocesser::readRFEfficiencies()", file_YES_TM->IsOpen(), "Failed to open ROOT file" );
      Info("readRFEfficiencies()", "REAL/FAKE efficiency: %s ", path_YES_TM.c_str() );

      Info("readRFEfficiencies()", "REAL/FAKE efficiency (probe NOT SLT-TRIGGER-MATCHED) from directory: %s ", m_EFF_NO_TM_dir.c_str() );

      std::string path_NO_TM = m_EFF_NO_TM_dir + m_Efficiency_Filename;
      file_NO_TM = TFile::Open(path_NO_TM.c_str());
      HTOP_RETURN_CHECK( "HTopMultilepNTupReprocesser::readRFEfficiencies()", file_NO_TM->IsOpen(), "Failed to open ROOT file" );
      Info("readRFEfficiencies()", "REAL/FAKE efficiency: %s ", path_NO_TM.c_str() );

  }

  std::vector<std::string> efficiencies = { "Real","Fake" };
  if ( m_useScaledFakeEfficiency ) efficiencies.push_back("ScaledFake");
  std::vector<std::string> leptons      = { "El","Mu" };
  std::vector<std::string> variables    = { "Pt" };
  if ( m_useEtaParametrisation ) variables.push_back("Eta");
  std::vector<std::string> sysdirections = { "up","dn" };

  // Parse input systematic groups, split by comma, and put into a vector

  std::vector<std::string> systematic_groups;
  ANA_CHECK( this->tokenize( ',', systematic_groups, m_systematics_list ) );

  int n_sysbins;

  for ( const auto& eff : efficiencies ) {

      std::string path;

      if ( eff.compare("Real") == 0 ) {

	  Info("readRFEfficiencies()", "REAL efficiency from directory: %s ", m_REFF_dir.c_str() );

	  path = m_REFF_dir + m_Efficiency_Filename;

      } else if ( eff.compare("Fake") == 0 ) {
	  if ( m_FEFF_dir.compare(m_REFF_dir) != 0 ) {
	      Warning("readRFEfficiencies()", "FAKE efficiency is going to be read from %s. Check whether it's really what you want...", m_FEFF_dir.c_str());
	  } else {
	      Info("readRFEfficiencies()", "FAKE efficiency from same directory as REAL" );
	  }
	  path = m_FEFF_dir + m_Efficiency_Filename;
      }

      TFile *file = TFile::Open(path.c_str());

      TH1D *hist(nullptr), *hist_avg(nullptr), *hist_YES_TM(nullptr), *hist_NO_TM(nullptr);
      TEfficiency *teff(nullptr);

      HTOP_RETURN_CHECK( "HTopMultilepNTupReprocesser::readRFEfficiencies()", file->IsOpen(), "Failed to open ROOT file" );

      for ( const auto& lep : leptons ) {

	  for ( const auto& var : variables ) {

	      // Do not look at variables other than pT for muons

	      if ( lep.compare("Mu") == 0 && var.compare("Pt") != 0 ) { continue; }

	      std::string sys_append;

	      bool isNominal(false), isStat(false), nominal_read(false);

	      std::string histname = eff + "_" + lep + "_" + var + "_Efficiency_"  + process;

	      n_sysbins = get_object<TH1D>( *file, histname )->GetNbinsX()+1;

	      for ( const auto& sysgroup : systematic_groups ) {

		  isNominal = ( sysgroup.compare("Nominal") == 0 );
		  isStat    = ( sysgroup.compare("Stat") == 0 );

		  // For variables other than pT, just consider the nominal case

		  if ( !isNominal && var.compare("Pt") != 0 ) { continue; }

		  // Make sure only relevant systematic sources for this efficiency,flavour are read in

		  if ( eff.compare("Real") == 0 ) {
		      if ( sysgroup.find("PromptSS") != std::string::npos ) { continue; }
		      if ( sysgroup.find("QMisID") != std::string::npos )   { continue; }
		  }
		  if ( eff.compare("Fake") == 0 ) {
		      if ( sysgroup.find("FakesOS") != std::string::npos ) { continue; }
		      if ( lep.compare("Mu") == 0 && sysgroup.find("QMisID") != std::string::npos )   { continue; }
		  }

		  std::cout << "" << std::endl;
		  Info("readRFEfficiencies()", "Reading inputs for systematic group: ===> %s", sysgroup.c_str() );
		  std::cout << "" << std::endl;

		  if ( !m_correlatedMMWeights ) {

		      for ( int bin(1);  bin <= n_sysbins; ++bin ) {

			  if ( isNominal && nominal_read ) { break; } // Do this only once for the nominal case

			  std::string sys = eff + "_" + lep + "_" + var + "_" + sysgroup;
			  if ( !isNominal ) { sys += "_" + std::to_string(bin); }

			  m_systematics.push_back(sys);

			  bool stat_read(false);

			  for ( const auto& dir : sysdirections ) {

			      if ( isNominal && nominal_read ) { break; } // Do this only once for the nominal case
			      if ( isStat    && stat_read    ) { break; } // Do this only once for the stat case

			      sys_append = ( isNominal || isStat ) ? "" : ( "_" + sysgroup + "_" + dir +  "_" + std::to_string(bin) );

			      histname  = eff + "_" + lep + "_" + var + "_Efficiency_"  + process + sys_append; // Name of efficiency TObject in input ROOT file

			      std::cout << "\t\t\t\t  " << lep << "," << eff << "," << var << " efficiency - input TH1 name: " << histname << std::endl;

			      hist  = get_object<TH1D>( *file,  histname );
			      hist->SetDirectory(0);

			      if ( m_useTEfficiency ) {
				  teff  = get_object<TEfficiency>( *file, this->str_replace( histname, "Efficiency", "TEfficiency" ).c_str() );
				  teff->SetDirectory(0);
			      }

			      if ( m_useTrigMatchingInfo ) {

				  hist_YES_TM = get_object<TH1D>( *file_YES_TM, histname );
				  hist_NO_TM  = get_object<TH1D>( *file_NO_TM, histname );

				  hist_YES_TM->SetDirectory(0);
				  hist_NO_TM->SetDirectory(0);
			      }

			      // Fill maps for later usage

			      std::string mapkey, mapkeyhist, mapkeyhist_yes_tm, mapkeyhist_no_tm;
			      if ( var.compare("Pt") == 0 ) {
				  if ( eff.compare("Real") == 0 ) {
				      mapkey     = "pt_reff";
				      mapkeyhist = "pt_reff_hist";
				  } else if ( eff.compare("Fake") == 0 ) {
				      mapkey     = "pt_feff";
				      mapkeyhist = "pt_feff_hist";
				  }
			      } else if ( var.compare("Eta") == 0 ) {
				  if ( eff.compare("Real") == 0 ) {
				      mapkey     = "eta_reff";
				      mapkeyhist = "eta_reff_hist";
				  } else if ( eff.compare("Fake") == 0 ) {
				      mapkey     = "eta_feff";
				      mapkeyhist = "eta_feff_hist";
				  }
			      }
			      mapkeyhist_yes_tm = mapkeyhist + "_YES_TM";
			      mapkeyhist_no_tm  = mapkeyhist + "_NO_TM";

			      std::string syskey("");
			      if      ( isNominal ){ syskey = sysgroup; }
			      else if ( isStat )   { syskey = sysgroup + "_" + std::to_string(bin); }
			      else                 { syskey = sysgroup + "_" + dir + "_" + std::to_string(bin); }

			      std::cout << "\t\t\t\t  Storing efficiency histogram in map w/ the following key: " << syskey << std::endl;

			      // Save in the histogram map a clone of the denominator histogram associated to the TEfficiency object in order to access the axis binning
			      // If we are not using TEfficiency, take the TH1 efficiency histogram itself (use the denominator "total" histogram by convention)
			      //
			      // NB: Calling GetCopyTotalHisto() transfer the ownership of the histogram pointer to the user. This introduces a memory leak in the code,
			      // as we don't explicitly call delete anywhere. However, this is harmless, since this is executed only once per job.

			      if ( lep.compare("El") == 0 ) {
				  m_el_teff_map[syskey][mapkey]  = teff;
				  m_el_hist_map[syskey][mapkeyhist] = ( m_useTEfficiency ) ? dynamic_cast<TH1D*>( teff->GetCopyTotalHisto() ) : hist;
				  if ( m_useTrigMatchingInfo ) {
				      m_el_hist_map[syskey][mapkeyhist_yes_tm] = hist_YES_TM;
				      m_el_hist_map[syskey][mapkeyhist_no_tm]  = hist_NO_TM;
				  }
			      } else if ( lep.compare("Mu") == 0 ) {
				  m_mu_hist_map[syskey][mapkeyhist] = ( m_useTEfficiency ) ? dynamic_cast<TH1D*>( teff->GetCopyTotalHisto() ) : hist;
				  m_mu_teff_map[syskey][mapkey] = teff;
				  if ( m_useTrigMatchingInfo ) {
				      m_mu_hist_map[syskey][mapkeyhist_yes_tm] = hist_YES_TM;
				      m_mu_hist_map[syskey][mapkeyhist_no_tm]  = hist_NO_TM;
				  }
			      }

			      // Calculate average efficiency to normalise (pT * eta) 1D efficiency (ignore this if using pT-only parametrisation).
			      //
			      // This factor is the same for eta and pT r/f histograms (it's just Integral(N) / Integral(D) for the efficiency definition ): use pT
			      // (If using TEfficiency, can get the TH1 objects that were used for measuring efficiency directly from the TEfficiency object.
			      // Otherwise, the average efficiency histogram must be already in the input file)

			      // Store this only for electrons, as for muons there is no other dependency than on pT.

			      if ( m_useEtaParametrisation && lep.compare("El") == 0 ) {

				  if ( var.compare("Pt") == 0 ) {

				      if ( !teff ) {
					  if ( eff.compare("Real") == 0 ) {
					      hist_avg = get_object<TH1D>( *file,  histname + "_AVG" );
					      hist_avg->SetDirectory(0);
					      m_el_reff_avg[syskey] = hist_avg->GetBinContent(1);
					  }
					  if ( eff.compare("Fake") == 0 ) {
					      hist_avg = get_object<TH1D>( *file,  histname + "_AVG" );
					      hist_avg->SetDirectory(0);
					      m_el_feff_avg[syskey] = hist_avg->GetBinContent(1);
					  }
				      } else {
					  if ( eff.compare("Real") == 0 ) {
					      m_el_reff_avg[syskey] = ( teff->GetPassedHistogram()->Integral(1,teff->GetPassedHistogram()->GetNbinsX()+1) ) / ( teff->GetTotalHistogram()->Integral(1,teff->GetTotalHistogram()->GetNbinsX()+1) );
					  }
					  if ( eff.compare("Fake") == 0 ) {
					      m_el_feff_avg[syskey] = ( teff->GetPassedHistogram()->Integral(1,teff->GetPassedHistogram()->GetNbinsX()+1) ) / ( teff->GetTotalHistogram()->Integral(1,teff->GetTotalHistogram()->GetNbinsX()+1) );
					  }
				      }
				  }
			      }

			      if ( isNominal ) { nominal_read = true; }
			      if ( isStat )    { stat_read = true; }

			  } // loop over sys directions

		      } // loop over sys bins

		  } else {

		      std::string sys = eff + "_" + lep + "_" + var + "_" + sysgroup;

		      m_systematics.push_back(sys);

		      bool stat_read(false);

		      for ( const auto& dir : sysdirections ) {

			  if ( isNominal && nominal_read ) { break; } // Do this only once for the nominal case
			  if ( isStat    && stat_read    ) { break; } // Do this only once for the stat case

			  sys_append = ( isNominal || isStat ) ? "" : ( "_" + sysgroup + "_" + dir ); // Check formatting of input efficiency file (need to store a shifted efficiency for all bins simultaneously )

			  histname  = eff + "_" + lep + "_" + var + "_Efficiency_"  + process + sys_append; // Name of efficiency TObject in input ROOT file

			  std::cout << "\t\t\t\t  " << lep << "," << eff << "," << var << " efficiency - input TH1 name: " << histname << std::endl;

			  hist  = get_object<TH1D>( *file,  histname );
			  hist->SetDirectory(0);

			  if ( m_useTEfficiency ) {
			      teff  = get_object<TEfficiency>( *file, this->str_replace( histname, "Efficiency", "TEfficiency" ).c_str() );
			      teff->SetDirectory(0);
			  }

			  if ( m_useTrigMatchingInfo ) {

			      hist_YES_TM = get_object<TH1D>( *file_YES_TM, histname );
			      hist_NO_TM  = get_object<TH1D>( *file_NO_TM, histname );

			      hist_YES_TM->SetDirectory(0);
			      hist_NO_TM->SetDirectory(0);
			  }

			  // Fill maps for later usage

			  std::string mapkey, mapkeyhist, mapkeyhist_yes_tm, mapkeyhist_no_tm;
			  if ( var.compare("Pt") == 0 ) {
			      if ( eff.compare("Real") == 0 ) {
				  mapkey     = "pt_reff";
				  mapkeyhist = "pt_reff_hist";
			      } else if ( eff.compare("Fake") == 0 ) {
				  mapkey     = "pt_feff";
				  mapkeyhist = "pt_feff_hist";
			      }
			  } else if ( var.compare("Eta") == 0 ) {
			      if ( eff.compare("Real") == 0 ) {
				  mapkey     = "eta_reff";
				  mapkeyhist = "eta_reff_hist";
			      } else if ( eff.compare("Fake") == 0 ) {
				  mapkey     = "eta_feff";
				  mapkeyhist = "eta_feff_hist";
			      }
			  }
			  mapkeyhist_yes_tm = mapkeyhist + "_YES_TM";
			  mapkeyhist_no_tm  = mapkeyhist + "_NO_TM";

			  std::string syskey = ( isNominal || isStat ) ? sysgroup : ( sysgroup + "_" + dir ); // No need to append bin idx here, as we consider efficiencies w/ all bins shifted up/dn simultaneously.

			  std::cout << "\t\t\t\t  Storing efficiency histogram in map w/ the following key: " << syskey << std::endl;

			  // Save in the histogram map a clone of the denominator histogram associated to the TEfficiency object in order to access the axis binning
			  // If we are not using TEfficiency, take the TH1 efficiency histogram itself (use the denominator "total" histogram by convention)
			  //
			  // NB: Calling GetCopyTotalHisto() transfer the ownership of the histogram pointer to the user. This introduces a memory leak in the code,
			  // as we don't explicitly call delete anywhere. However, this is harmless, since this is executed only once per job.

			  if ( lep.compare("El") == 0 ) {
			      m_el_teff_map[syskey][mapkey]  = teff;
			      m_el_hist_map[syskey][mapkeyhist] = ( m_useTEfficiency ) ? dynamic_cast<TH1D*>( teff->GetCopyTotalHisto() ) : hist;
			      if ( m_useTrigMatchingInfo ) {
				  m_el_hist_map[syskey][mapkeyhist_yes_tm] = hist_YES_TM;
				  m_el_hist_map[syskey][mapkeyhist_no_tm]  = hist_NO_TM;
			      }
			  } else if ( lep.compare("Mu") == 0 ) {
			      m_mu_hist_map[syskey][mapkeyhist] = ( m_useTEfficiency ) ? dynamic_cast<TH1D*>( teff->GetCopyTotalHisto() ) : hist;
			      m_mu_teff_map[syskey][mapkey] = teff;
			      if ( m_useTrigMatchingInfo ) {
				  m_mu_hist_map[syskey][mapkeyhist_yes_tm] = hist_YES_TM;
				  m_mu_hist_map[syskey][mapkeyhist_no_tm]  = hist_NO_TM;
			      }
			  }

			  // Calculate average efficiency to normalise (pT * eta) 1D efficiency (ignore this if using pT-only parametrisation).
			  //
			  // This factor is the same for eta and pT r/f histograms (it's just Integral(N) / Integral(D) for the efficiency definition ): use pT
			  // (If using TEfficiency, can get the TH1 objects that were used for measuring efficiency directly from the TEfficiency object.
			  // Otherwise, the average efficiency histogram must be already in the input file)

			  // Store this only for electrons, as for muons there is no other dependency than on pT.

			  if ( m_useEtaParametrisation && lep.compare("El") == 0 ) {

			      if ( var.compare("Pt") == 0 ) {

				  if ( !teff ) {
				      if ( eff.compare("Real") == 0 ) {
					  hist_avg = get_object<TH1D>( *file,  histname + "_AVG" );
					  hist_avg->SetDirectory(0);
					  m_el_reff_avg[syskey] = hist_avg->GetBinContent(1);
				      }
				      if ( eff.compare("Fake") == 0 ) {
					  hist_avg = get_object<TH1D>( *file,  histname + "_AVG" );
					  hist_avg->SetDirectory(0);
					  m_el_feff_avg[syskey] = hist_avg->GetBinContent(1);
				      }
				  } else {
				      if ( eff.compare("Real") == 0 ) {
					  m_el_reff_avg[syskey] = ( teff->GetPassedHistogram()->Integral(1,teff->GetPassedHistogram()->GetNbinsX()+1) ) / ( teff->GetTotalHistogram()->Integral(1,teff->GetTotalHistogram()->GetNbinsX()+1) );
				      }
				      if ( eff.compare("Fake") == 0 ) {
					  m_el_feff_avg[syskey] = ( teff->GetPassedHistogram()->Integral(1,teff->GetPassedHistogram()->GetNbinsX()+1) ) / ( teff->GetTotalHistogram()->Integral(1,teff->GetTotalHistogram()->GetNbinsX()+1) );
				      }
				  }
			      }
			  }

			  if ( isNominal ) { nominal_read = true; }
			  if ( isStat )    { stat_read = true; }

		      } // loop over sys directions

		  } // closes if ( m_correlatedMMWeights )

	      } // loop over systematic sources

	  } // loop over variables

      } // loop over leptons

  } // loop over efficieny types

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode HTopMultilepNTupReprocesser :: getMMEfficiencyAndError( std::shared_ptr<leptonObj> lep, std::vector<float>& efficiency, const std::string& type )
{

    float error_up(0.0), error_dn(0.0);

    float pt  = lep.get()->pt/1e3; // Must be in GeV!
    float eta = ( lep.get()->flavour == 13 ) ? lep.get()->eta : lep.get()->etaBE2;

    float this_low_edge_pt(-1.0), this_up_edge_pt(-1.0);
    float this_low_edge_eta(-999.0), this_up_edge_eta(-999.0);

    std::map< std::string, std::map< std::string, TH1D* > >        *histograms    = ( lep.get()->flavour == 13 ) ? &m_mu_hist_map : &m_el_hist_map;
    std::map< std::string, std::map< std::string, TEfficiency* > > *tefficiencies = ( lep.get()->flavour == 13 ) ? &m_mu_teff_map : &m_el_teff_map;

    std::string key_pt, key_eta, key_pt_teff, key_eta_teff;
    std::string key_pt_NO_TM, key_pt_YES_TM;

    if ( type.compare("REAL") == 0 ) {
	key_pt        = "pt_reff_hist";
	key_pt_teff   = "pt_reff";
	key_eta       = "eta_reff_hist";
	key_eta_teff  = "eta_reff";
	key_pt_NO_TM  = "pt_reff_hist_NO_TM";
	key_pt_YES_TM = "pt_reff_hist_YES_TM";
    } else if ( type.compare("FAKE") == 0 ) {
	key_pt        = "pt_feff_hist";
	key_pt_teff   = "pt_feff";
	key_eta       = "eta_feff_hist";
	key_eta_teff  = "eta_feff";
	key_pt_NO_TM  = "pt_feff_hist_NO_TM";
	key_pt_YES_TM = "pt_feff_hist_YES_TM";
    }

    if ( m_verbose ) { Info("getMMEfficiencyAndError()", "\tReading %s efficiency...", type.c_str() ); }

    std::string syskey_up_pt, syskey_dn_pt, syskey_up_eta, syskey_dn_eta;

    std::vector<std::string> tokens;

    TH1D* hist_nominal_pt  = (histograms->find("Nominal")->second).find(key_pt)->second;
    TH1D* hist_nominal_eta = (histograms->find("Nominal")->second).find(key_eta)->second;

    TH1D* hist_nominal_pt_YES_TM(nullptr);
    TH1D* hist_nominal_pt_NO_TM(nullptr);

    if ( m_useTrigMatchingInfo ) {

        hist_nominal_pt_YES_TM = (histograms->find("Nominal")->second).find(key_pt_YES_TM)->second;
        hist_nominal_pt_NO_TM  = (histograms->find("Nominal")->second).find(key_pt_NO_TM)->second;

    }

    // Get the number of bins. Use nominal

    int nbins_pt = hist_nominal_pt->GetNbinsX()+1;

    float overflow_pt_up_edge = hist_nominal_pt->GetXaxis()->GetBinUpEdge(nbins_pt);

    bool isBinOverflowPt(false);

    bool isNominal = ( m_this_syst.find("Nominal") != std::string::npos );
    bool isStat    = ( m_this_syst.find("Stat") != std::string::npos );

    // Loop over number of pt bins
    // Do not consider underflow, i.e. 0th bin

    for ( int p(1); p <= nbins_pt; ++p ) {

	isBinOverflowPt = hist_nominal_pt->IsBinOverflow(p);

	this_low_edge_pt = hist_nominal_pt->GetBinLowEdge(p);
	this_up_edge_pt  = hist_nominal_pt->GetBinLowEdge(p+1);

	if ( m_verbose ) { Info("getMMEfficiencyAndError()","\t\tpT bin %i : [%.0f,%.0f] GeV", p, this_low_edge_pt, this_up_edge_pt ); }

	if ( ( pt >= this_low_edge_pt && pt < this_up_edge_pt ) || ( isBinOverflowPt && pt >= overflow_pt_up_edge ) ) {

	    float eff_pt(1.0), eff_pt_err_up(0.0), eff_pt_err_dn(0.0);

	    // The central value for the efficiency will always be read from the nominal efficinecy histogram if:
	    //
	    // -) m_this_syst contains "Nominal"
	    // OR
	    // -) the bin in question does NOT correspond to the bin for *this* systematic varied histogram (only when looking at fully uncorrelated syst variations)
	    // OR
	    // -) lepton is el and m_this_syst starts w/ "Mu",
	    // OR
	    // -) lepton is mu and m_this_syst starts w/ "El"
	    // OR
	    // -) checking r eff and m_this_syst contains "Fake"
	    // OR
	    // -) checking f eff and m_this_syst contains "Real"
	    // OR
	    // -) m_this_syst does NOT contain Pt
	    //
	    // , need to read "Nominal" for both up and dn  --> syskey_up_pt = syskey_dn_pt = "Nominal"

	    tokens.clear();
	    ANA_CHECK( this->tokenize( '_', tokens, m_this_syst ) );

            bool readNominalPt(false);

	    if ( ( isNominal ) ||
	         ( ( lep.get()->flavour == 13 )  && ( m_this_syst.find("El_") != std::string::npos ) )   ||
		 ( ( lep.get()->flavour == 11 )  && ( m_this_syst.find("Mu_") != std::string::npos ) )   ||
		 ( ( type.compare("REAL") == 0 ) && ( m_this_syst.find("Fake_") != std::string::npos ) ) ||
		 ( ( type.compare("FAKE") == 0 ) && ( m_this_syst.find("Real_") != std::string::npos ) ) ||
		 ( m_this_syst.find("Pt") == std::string::npos )
		)
	    {
		syskey_up_pt = syskey_dn_pt = "Nominal";
		readNominalPt = true;
	    } else {

		// std::cout << "\t\ttokens for this syst: " << tokens.size() << std::endl;
		// unsigned int idx(0);
		// for ( const auto& t : tokens ) {
		//     std::cout << "\t\t\tt[" << idx << "] = " << t << std::endl;
		//     ++idx;
		// }
		// std::cout << "" << std::endl;

		std::string addon("");
		unsigned int endtoken = ( !m_correlatedMMWeights ) ? tokens.size() - 1 : tokens.size();
		for ( unsigned int idx(3); idx < endtoken; ++idx ) { addon = addon + tokens.at(idx) + "_"; }

		if ( !m_correlatedMMWeights ) {
		    if ( p != std::stoi(tokens.back()) ) {
			syskey_up_pt = syskey_dn_pt = "Nominal";
			readNominalPt = true;
		    } else {
			syskey_up_pt = ( isStat ) ? "Nominal" : addon + "up_" + tokens.back();
			syskey_dn_pt = ( isStat ) ? "Nominal" : addon + "dn_" + tokens.back();
		    }
		} else {
		    syskey_up_pt = ( isStat ) ? "Nominal" : addon + "up";
		    syskey_dn_pt = ( isStat ) ? "Nominal" : addon + "dn";
		}

	    }

	    if ( m_verbose ) { Info("getMMEfficiencyAndError()", "\t\t===> Retrieving nominal pT histogram and variations from map w/ key: %s (up), %s (dn)", syskey_up_pt.c_str(), syskey_dn_pt.c_str() ); }

	    // NB: "Stat" systematics need special treatment
	    // up/dn variations for *this* systematic bin are obtained by reading the stat uncertainty of the bin itself for the *nominal* hist, rather than a different histogram.
	    // If this bin is not corresponding to *this* systematic bin, then get an error of 0

	    if ( !m_useTEfficiency ) {

		eff_pt = hist_nominal_pt->GetBinContent(p);

		if ( m_useTrigMatchingInfo ) {
		  eff_pt = ( lep.get()->trigmatched ) ? hist_nominal_pt_YES_TM->GetBinContent(p) : hist_nominal_pt_NO_TM->GetBinContent(p);
		}

		if ( isStat ) {

		  eff_pt_err_up = ( readNominalPt ) ? 0 : hist_nominal_pt->GetBinError(p);
		  eff_pt_err_dn = ( readNominalPt ) ? 0 : hist_nominal_pt->GetBinError(p);

		  if ( m_useTrigMatchingInfo ) {

		    if ( readNominalPt ) {
		      eff_pt_err_up = eff_pt_err_dn = 0;
		    } else {
		      eff_pt_err_up = ( lep.get()->trigmatched ) ? hist_nominal_pt_YES_TM->GetBinError(p) : hist_nominal_pt_NO_TM->GetBinError(p);
		      eff_pt_err_dn = ( lep.get()->trigmatched ) ? hist_nominal_pt_YES_TM->GetBinError(p) : hist_nominal_pt_NO_TM->GetBinError(p);
		    }

		  }

		} else {

		  eff_pt_err_up = (histograms->find(syskey_up_pt)->second).find(key_pt)->second->GetBinContent(p);
		  eff_pt_err_dn = (histograms->find(syskey_dn_pt)->second).find(key_pt)->second->GetBinContent(p);

		  if ( m_useTrigMatchingInfo ) {
		    eff_pt_err_up = ( lep.get()->trigmatched ) ? (histograms->find(syskey_up_pt)->second).find(key_pt_YES_TM)->second->GetBinContent(p) : (histograms->find(syskey_up_pt)->second).find(key_pt_NO_TM)->second->GetBinContent(p);
		    eff_pt_err_dn = ( lep.get()->trigmatched ) ? (histograms->find(syskey_dn_pt)->second).find(key_pt_YES_TM)->second->GetBinContent(p) : (histograms->find(syskey_dn_pt)->second).find(key_pt_NO_TM)->second->GetBinContent(p);
		  }

		}

	    } else {
		eff_pt	      = (tefficiencies->find("Nominal")->second).find(key_pt_teff)->second->GetEfficiency(p);
		if ( isStat ) {
		  eff_pt_err_up = ( readNominalPt ) ? 0 : (tefficiencies->find("Nominal")->second).find(key_pt_teff)->second->GetEfficiencyErrorUp(p);
		  eff_pt_err_dn = ( readNominalPt ) ? 0 : (tefficiencies->find("Nominal")->second).find(key_pt_teff)->second->GetEfficiencyErrorLow(p);
		} else {
		  eff_pt_err_up = (tefficiencies->find(syskey_up_pt)->second).find(key_pt_teff)->second->GetEfficiency(p);
		  eff_pt_err_dn = (tefficiencies->find(syskey_dn_pt)->second).find(key_pt_teff)->second->GetEfficiency(p);
		}
	    }

	    if ( m_verbose ) {
		std::string trigmatch_str("");
		if ( m_useTrigMatchingInfo ) {
		    trigmatch_str = ( lep.get()->trigmatched ) ? ", TRIGMATCHED (SLT)" : ", NOT TRIGMATCHED (SLT)";
		}
		Info("getMMEfficiencyAndError()", "\t\t\tLepton pT = %.3f GeV, flavour = %i%s ==> Reading %s efficiency in pT bin [%.0f,%.0f] GeV: eff_pt = %.3f", pt, lep.get()->flavour, trigmatch_str.c_str(), type.c_str(), this_low_edge_pt, this_up_edge_pt, eff_pt );
	    }

	    // If looking at eta parametrisation, always take the nominal value, except when looking at "Stat".
	    // In that case, add the statistical error for the relevant eta bin to this pt bin stat error

	    float eff_eta(1.0), eff_eta_err_up(0.0), eff_eta_err_dn(0.0);

	    if ( m_useEtaParametrisation && lep.get()->flavour == 11 ) {

		// Get the number of bins. Use nominal

		int nbins_eta = hist_nominal_eta->GetNbinsX()+1;

		float overflow_eta_up_edge = hist_nominal_eta->GetXaxis()->GetBinUpEdge(nbins_eta);

		bool isBinOverflowEta(false);

		// Loop over number of eta bins
		// Do not consider underflow, i.e. 0th bin

		for ( int e(1); e <= nbins_eta; ++e ) {

		    isBinOverflowEta = hist_nominal_eta->IsBinOverflow(e);

		    this_low_edge_eta = hist_nominal_eta->GetBinLowEdge(e);
		    this_up_edge_eta  = hist_nominal_eta->GetBinLowEdge(e+1);

		    if ( m_verbose ) { Info("getMMEfficiencyAndError()","\t\t|eta| bin %i : [%.3f,%.3f]", e, this_low_edge_eta, this_up_edge_eta ); }

		    if ( ( fabs(eta) >= this_low_edge_eta && fabs(eta) < this_up_edge_eta ) || ( isBinOverflowEta && fabs(eta) >= overflow_eta_up_edge ) ) {

			syskey_up_eta = syskey_dn_eta = "Nominal";

                	bool readNominalEta(false);

	        	if ( ( !isStat ) || // Make sure to get the stat variation on eta eff only when this syst is "Stat". In all other cases, only the nominal eta eff value will be taken (error will be set to 0)
	        	     ( ( lep.get()->flavour == 11 )  && ( m_this_syst.find("Mu_") != std::string::npos ) )   ||
	        	     ( ( type.compare("REAL") == 0 ) && ( m_this_syst.find("Fake_") != std::string::npos ) ) ||
	        	     ( ( type.compare("FAKE") == 0 ) && ( m_this_syst.find("Real_") != std::string::npos ) ) )
	        	{
	        	    readNominalEta = true;
	        	} else if ( isStat ) {
			    if ( !m_correlatedMMWeights ) {
				if ( p != std::stoi(tokens.back()) ) {
				    readNominalEta = true;
				}
			    } else {
				readNominalEta = true;
			    }
			}

			if ( m_verbose ) { Info("getMMEfficiencyAndError()", "\t\t===> Retrieving nominal eta histogram and variations from map w/ key: %s (up), %s (dn)", syskey_up_eta.c_str(), syskey_dn_eta.c_str() ); }

	    		if ( !m_useTEfficiency ) {
	    		    eff_eta	   = hist_nominal_eta->GetBinContent(e);
			    eff_eta_err_up = ( readNominalEta ) ? 0 : hist_nominal_eta->GetBinError(e);
			    eff_eta_err_dn = ( readNominalEta ) ? 0 : hist_nominal_eta->GetBinError(e);
	    		} else {
	    		    eff_eta	   = (tefficiencies->find("Nominal")->second).find(key_eta_teff)->second->GetEfficiency(e);
			    eff_eta_err_up = ( readNominalEta ) ? 0 : (tefficiencies->find("Nominal")->second).find(key_eta_teff)->second->GetEfficiencyErrorUp(e);
			    eff_eta_err_dn = ( readNominalEta ) ? 0 : (tefficiencies->find("Nominal")->second).find(key_eta_teff)->second->GetEfficiencyErrorLow(e);
	                }

			if ( m_verbose ) {
			    Info("getMMEfficiencyAndError()", "\t\t\tLepton |eta| = %.3f, flavour = %i ==> Reading %s efficiency in |eta| bin [%.3f,%.3f]: eff_eta = %.3f", fabs(eta), lep.get()->flavour, type.c_str(), this_low_edge_eta, this_up_edge_eta, eff_eta );
			}

			break;
		    }
		} // close loop on eta bins
	    }

	    // Nominal

	    efficiency.at(0) = eff_pt;

	    eff_pt_err_up  = ( isStat ) ? eff_pt_err_up : fabs( eff_pt - eff_pt_err_up );
	    eff_pt_err_dn  = ( isStat ) ? eff_pt_err_dn : fabs( eff_pt - eff_pt_err_dn );

	    error_up	 = eff_pt_err_up;
	    error_dn	 = eff_pt_err_dn;

	    // UP syst

	    efficiency.at(1) = ( eff_pt + error_up );

	    // DN syst

	    if ( eff_pt - error_dn > 0 ) { efficiency.at(2) = ( eff_pt - error_dn ); }
	    else			 { efficiency.at(2) = 0.0; }

	    if ( m_useEtaParametrisation && lep.get()->flavour == 11 ) {

		float eff_avg(1.0);
		if      ( type.compare("REAL") == 0 ) { eff_avg = m_el_reff_avg["Nominal"]; }
		else if ( type.compare("FAKE") == 0 ) { eff_avg = m_el_feff_avg["Nominal"]; }

		if ( m_verbose ) { Info("getMMEfficiencyAndError()", "\t\tnormalisation factor (<eff>) = %.3f", eff_avg ); }

		efficiency.at(0) = ( eff_pt * eff_eta ) / eff_avg;

		// Assuming  eff_pt,eff_eta are independent, this is the error on the product
		// ( the constant factor at denominator will be put back later in the def of Efficiency...)

		error_up = ( isStat ) ? sqrt( (eff_eta*eff_pt_err_up)*(eff_eta*eff_pt_err_up) + (eff_pt*eff_eta_err_up)*(eff_pt*eff_eta_err_up) ) : error_up;
		error_dn = ( isStat ) ? sqrt( (eff_eta*eff_pt_err_dn)*(eff_eta*eff_pt_err_dn) + (eff_pt*eff_eta_err_dn)*(eff_pt*eff_eta_err_dn) ) : error_dn;

		efficiency.at(1) = ( (eff_pt * eff_eta) + error_up ) / eff_avg;
		if ( (eff_pt * eff_eta) - error_dn > 0 ) { efficiency.at(2) = ( (eff_pt * eff_eta) - error_dn ) / eff_avg; }
		else				         { efficiency.at(2) = 0.0; }
	    }

	    break;
	}

    } // close loop on pT bins

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

    if ( (r0.at(0) == 0) || (r1.at(0) == 0) ||
         (r0.at(0) <= f0.at(0))             ||
	 (r1.at(0) <= f1.at(0))
	)
    {

	if ( m_debug ) { // TEMP! This should NOT be enclosed in debug opt!
	    Warning("getMMWeightAndError()", "Warning! The Matrix Method cannot be applied to Entry: %u, EventNumber: %u because:\n", static_cast<uint32_t>(m_numEntry), static_cast<uint32_t>(m_EventNumber) );

	    if ( (r0.at(0) == 0) || (r1.at(0) == 0) ) {
		std::cout << "r0 = " << r0.at(0) << ", r1 = " << r1.at(0) << std::endl;
	    }
	    if ( r0.at(0) <= f0.at(0) ) {
		std::cout << "r0 = " << r0.at(0) << ", f0 = " << f0.at(0) <<  " ==> r0 <= f0 !! " << std::endl;
	    }
	    if ( r1.at(0) <= f1.at(0) ) {
		std::cout << "r1 = " << r1.at(0) << ", f1 = " << f1.at(0) <<  " ==> r1 <= f1 !! " << std::endl;
	    }
	    Warning("getMMWeightAndError()", "Assigning MMWeight (nominal) = 0, aka will remove the event ...");
	}

        return EL::StatusCode::SUCCESS;
    }

    // Calculate nominal MM weight

    mm_weight.at(0) = matrix_equation( f0.at(0), f1.at(0), r0.at(0), r1.at(0) );

    // Calculate MM weight with variations up/dn for r/f for *this* systematic

    float r0up = ( r0.at(1) > 1.0 ) ? 1.0 : r0.at(1);
    float r1up = ( r1.at(1) > 1.0 ) ? 1.0 : r1.at(1);
    float r0dn = ( r0.at(2) < 0.0 ) ? 0.0 : r0.at(2);
    float r1dn = ( r1.at(2) < 0.0 ) ? 0.0 : r1.at(2);

    float f0up = ( f0.at(1) > 1.0 ) ? 1.0 : f0.at(1);
    float f1up = ( f1.at(1) > 1.0 ) ? 1.0 : f1.at(1);
    float f0dn = ( f0.at(2) < 0.0 ) ? 0.0 : f0.at(2);
    float f1dn = ( f1.at(2) < 0.0 ) ? 0.0 : f1.at(2);

    mm_weight.at(1) = matrix_equation( f0up, f1up, r0up, r1up );
    mm_weight.at(1) = ( !std::isnan(mm_weight.at(1)/mm_weight.at(0)) && !std::isinf(mm_weight.at(1)/mm_weight.at(0)) ) ? mm_weight.at(1)/mm_weight.at(0) : 0.0;

    mm_weight.at(2) = matrix_equation( f0dn, f1dn, r0dn, r1dn );
    mm_weight.at(2) = ( !std::isnan(mm_weight.at(2)/mm_weight.at(0)) && !std::isinf(mm_weight.at(2)/mm_weight.at(0)) ) ? mm_weight.at(2)/mm_weight.at(0) : 0.0;

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

    if ( m_debug ) {
      std::cout << "" << std::endl;
      Info("calculateMMWeights()", "Now checking systematic variation: ===> %s", m_this_syst.c_str() );
      std::cout << "" << std::endl;
    }

    // If is not a dileptonic/trileptonic event, return

    if ( m_event.get()->dilep_type <= 0 && m_event.get()->trilep_type <= 0 ) { return EL::StatusCode::SUCCESS; }

    std::shared_ptr<leptonObj> lep0 = m_leptons.at(0);
    std::shared_ptr<leptonObj> lep1 = m_leptons.at(1);

    // These are the "effective" r/f efficiencies for each lepton, obtained by reading the input r/f histogram(s)

    std::vector<float> r0 = { 1.0, 1.0, 1.0 };
    std::vector<float> r1 = { 1.0, 1.0, 1.0 };
    std::vector<float> f0 = { 1.0, 1.0, 1.0 };
    std::vector<float> f1 = { 1.0, 1.0, 1.0 };

    ANA_CHECK( this->getMMEfficiencyAndError( lep0, r0, "REAL" ) );
    ANA_CHECK( this->getMMEfficiencyAndError( lep1, r1, "REAL" ) );
    ANA_CHECK( this->getMMEfficiencyAndError( lep0, f0, "FAKE" ) );
    ANA_CHECK( this->getMMEfficiencyAndError( lep1, f1, "FAKE" ) );

    if ( m_debug ) {
        std::cout << "" << std::endl;
	Info("calculateMMWeights()", "===> Systematic variation: %s", m_this_syst.c_str() );
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

    // The final MM event weight:
    //
    // 1) component is NOMINAL (default is 0 ==> remove event at plotting)
    // 2) component is THIS-SYS up
    // 3) component is THIS-SYS dn

    std::vector<float> mm_weight = { 0.0, 1.0, 1.0 };

    // For variations, save relative weight wrt. nominal

    ANA_CHECK( this->getMMWeightAndError( mm_weight, r0, r1, f0, f1 ) );

    m_MMWeight_NOMINAL_out = mm_weight.at(0);

    m_MMWeight_out[m_this_syst].at(0) = mm_weight.at(1);
    m_MMWeight_out[m_this_syst].at(1) = mm_weight.at(2);

    return EL::StatusCode::SUCCESS;
}

