#include <HTopMultilepAnalysis/HTopMultilepMiniNTupMaker.h>

// ASG status code check
#include <AsgTools/MessageCheck.h>

// ROOT include(s)
#include "TObjArray.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TRegexp.h"
#include "TString.h"

// C++ include(s)
#include<algorithm>

using namespace MiniNTupMaker;

// this is needed to distribute the algorithm to the workers
ClassImp(HTopMultilepMiniNTupMaker)


HTopMultilepMiniNTupMaker :: HTopMultilepMiniNTupMaker(std::string className) :
    Algorithm(className),
    m_inputNTuple(nullptr),
    m_sumWeightsTree(nullptr),
    m_sumGenEventsHist(nullptr),
    m_outputNTuple(nullptr)
{
    // Here you put any code for the base initialization of variables,
    // e.g. initialize all pointers to 0.  Note that you should only put
    // the most basic initialization here, since this method will be
    // called on both the submission and the worker node.  Most of your
    // initialization code will go into histInitialize() and
    // initialize().

    Info("HTopMultilepMiniNTupMaker()", "Calling constructor");

    m_debug                = false;
    m_outputNTupName       = "nominal";
    m_outputNTupStreamName = "output";
    m_inputBranches        = "";
    m_useAlgSelect         = false;
    m_addStreamEventsHist  = false;
    m_useTruthTP           = false;
    m_useNominalTP         = true;
    m_ambiSolvingCrit      = "OF";
    m_lepSelForTP          = "MVA";
    m_jetTruthMatching     = false;
}



EL::StatusCode HTopMultilepMiniNTupMaker :: setupJob (EL::Job& job)
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

    // Add output stream for total sum of weighted events (for MC normalisation)

    if ( m_addStreamEventsHist ) {
	EL::OutputStream outputStreamSumW("sum_weights");
	job.outputAdd (outputStreamSumW);
    }

    return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepMiniNTupMaker :: histInitialize ()
{
    // Here you do everything that needs to be done at the very
    // beginning on each worker node, e.g. create histograms and output
    // trees.  This method gets called before any input files are
    // connected.

    //ANA_CHECK_SET_TYPE (EL::StatusCode);

    Info("histInitialize()", "Calling histInitialize");

    TFile *fileSW = ( m_addStreamEventsHist ) ?  wk()->getOutputFile ("sum_weights") :  wk()->getOutputFile ("output");
    fileSW->cd();

    m_sumGenEventsHist = new TH1F("TotalEventsW", "Total number of generated events", 2, 0.0, 2.0);
    m_sumGenEventsHist->GetXaxis()->SetBinLabel(1, "sum gen. events (RAW)");
    m_sumGenEventsHist->GetXaxis()->SetBinLabel(2, "sum gen. events (WEIGHTED)");

    m_sumGenEvents         = 0.0;
    m_sumGenEventsWeighted = 0.0 ;

    return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepMiniNTupMaker :: fileExecute ()
{
    // Here you do everything that needs to be done exactly once for every
    // single file, e.g. collect a list of all lumi-blocks processed

    //ANA_CHECK_SET_TYPE (EL::StatusCode);

    Info("fileExecute()", "Calling fileExecute");

    return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepMiniNTupMaker :: changeInput (bool firstFile)
{
    // Here you do everything you need to do when we change input files,
    // e.g. resetting branch addresses on trees.  If you are using
    // D3PDReader or a similar service this method is not needed.

    ANA_CHECK_SET_TYPE (EL::StatusCode);

    firstFile = firstFile;

    Info("changeInput()", "Calling changeInput. Now reading file : %s", wk()->inputFile()->GetName() );

    // Get the pointer to the main input TTree

    m_inputNTuple = wk()->tree();

    ANA_CHECK( this->enableSelectedBranches() );

    // Connect the branches of the input tree to the algorithm members

    m_inputNTuple->SetBranchAddress ("EventNumber",   			      &m_EventNumber);
    m_inputNTuple->SetBranchAddress ("RunNumber",   			      &m_RunNumber);
    m_inputNTuple->SetBranchAddress ("RunYear",   			      &m_RunYear);
    m_inputNTuple->SetBranchAddress ("passEventCleaning", 		      &m_passEventCleaning);
    m_inputNTuple->SetBranchAddress ("mc_channel_number",                       &m_mc_channel_number);

    m_inputNTuple->SetBranchAddress ("mcWeightOrg",			      &m_mcWeightOrg);
    m_inputNTuple->SetBranchAddress ("pileupEventWeight_090",		      &m_pileupEventWeight_090);
    m_inputNTuple->SetBranchAddress ("MV2c10_70_EventWeight",		      &m_MV2c10_70_EventWeight);
    m_inputNTuple->SetBranchAddress ("JVT_EventWeight",			      &m_JVT_EventWeight);

    m_inputNTuple->SetBranchAddress ("dilep_type",  			      &m_dilep_type);
    m_inputNTuple->SetBranchAddress ("trilep_type",  			      &m_trilep_type);
    m_inputNTuple->SetBranchAddress ("quadlep_type",  			      &m_quadlep_type);

    m_inputNTuple->SetBranchAddress ("nJets_OR", 			      &m_nJets_OR);
    m_inputNTuple->SetBranchAddress ("nJets_OR_MV2c10_70", 		      &m_nJets_OR_MV2c10_70);
    m_inputNTuple->SetBranchAddress ("nJets_OR_T", 			      &m_nJets_OR_T);
    m_inputNTuple->SetBranchAddress ("nJets_OR_T_MV2c10_70", 		      &m_nJets_OR_T_MV2c10_70);

    std::string branchname, branchtype, key;

    for ( const auto& var : m_LEP_INPUT_VARS ) {

	key = var.substr( 0, var.length() - 2 );

	branchtype = var.substr( var.length() - 1 );

	for (int idx(0); idx < 4; ++idx ) {

	    branchname = "lep_" + key + "_" + std::to_string(idx);

	    switch( idx ) {
	    case 0:
		if ( branchtype.compare("F") == 0 )   { m_inputNTuple->SetBranchAddress ( branchname.c_str(),   &m_lep0_INPUT_branches[key].f ); }
		if ( branchtype.compare("B") == 0 )   { m_inputNTuple->SetBranchAddress ( branchname.c_str(),   &m_lep0_INPUT_branches[key].c ); }
		if ( branchtype.compare("I") == 0 )   { m_inputNTuple->SetBranchAddress ( branchname.c_str(),   &m_lep0_INPUT_branches[key].i ); }
		break;
	    case 1:
		if ( branchtype.compare("F") == 0 )   { m_inputNTuple->SetBranchAddress ( branchname.c_str(),   &m_lep1_INPUT_branches[key].f ); }
		if ( branchtype.compare("B") == 0 )   { m_inputNTuple->SetBranchAddress ( branchname.c_str(),   &m_lep1_INPUT_branches[key].c ); }
		if ( branchtype.compare("I") == 0 )   { m_inputNTuple->SetBranchAddress ( branchname.c_str(),   &m_lep1_INPUT_branches[key].i ); }
		break;
	    case 2:
		if ( branchtype.compare("F") == 0 )   { m_inputNTuple->SetBranchAddress ( branchname.c_str(),   &m_lep2_INPUT_branches[key].f ); }
		if ( branchtype.compare("B") == 0 )   { m_inputNTuple->SetBranchAddress ( branchname.c_str(),   &m_lep2_INPUT_branches[key].c ); }
		if ( branchtype.compare("I") == 0 )   { m_inputNTuple->SetBranchAddress ( branchname.c_str(),   &m_lep2_INPUT_branches[key].i ); }
		break;
	    case 3:
		if ( branchtype.compare("F") == 0 )   { m_inputNTuple->SetBranchAddress ( branchname.c_str(),   &m_lep3_INPUT_branches[key].f ); }
		if ( branchtype.compare("B") == 0 )   { m_inputNTuple->SetBranchAddress ( branchname.c_str(),   &m_lep3_INPUT_branches[key].c ); }
		if ( branchtype.compare("I") == 0 )   { m_inputNTuple->SetBranchAddress ( branchname.c_str(),   &m_lep3_INPUT_branches[key].i ); }
		break;
	    default:
		break;
	    }
	}
    }

    m_inputNTuple->SetBranchAddress ("electron_passOR",                               &m_electron_passOR);
    m_inputNTuple->SetBranchAddress ("electron_pt", 	                              &m_electron_pt);
    m_inputNTuple->SetBranchAddress ("electron_eta",	      		              &m_electron_eta);
    m_inputNTuple->SetBranchAddress ("electron_EtaBE2",	                              &m_electron_EtaBE2);
    m_inputNTuple->SetBranchAddress ("electron_phi",  	                              &m_electron_phi);
    m_inputNTuple->SetBranchAddress ("electron_E",	      		              &m_electron_E);
    m_inputNTuple->SetBranchAddress ("electron_ID", 	                              &m_electron_ID);
    m_inputNTuple->SetBranchAddress ("electron_sigd0PV",	                      &m_electron_sigd0PV);
    m_inputNTuple->SetBranchAddress ("electron_z0SinTheta",     		      &m_electron_z0SinTheta);
    m_inputNTuple->SetBranchAddress ("electron_topoetcone20",   		      &m_electron_topoetcone20);
    m_inputNTuple->SetBranchAddress ("electron_ptvarcone20",    		      &m_electron_ptvarcone20);
    m_inputNTuple->SetBranchAddress ("electron_truthType",      		      &m_electron_truthType);
    m_inputNTuple->SetBranchAddress ("electron_truthOrigin",    		      &m_electron_truthOrigin);
    m_inputNTuple->SetBranchAddress ("electron_PromptLeptonIso_TagWeight",    	      &m_electron_PromptLeptonIso_TagWeight);
    m_inputNTuple->SetBranchAddress ("electron_ChargeIDBDTLoose",    		      &m_electron_ChargeIDBDTLoose);
    m_inputNTuple->SetBranchAddress ("electron_ChargeIDBDTMedium",    		      &m_electron_ChargeIDBDTMedium);
    m_inputNTuple->SetBranchAddress ("electron_ChargeIDBDTTight",    		      &m_electron_ChargeIDBDTTight);
    m_inputNTuple->SetBranchAddress ("electron_match_HLT_e24_lhmedium_L1EM20VH",      &m_electron_match_HLT_e24_lhmedium_L1EM20VH);
    m_inputNTuple->SetBranchAddress ("electron_match_HLT_e26_lhtight_nod0_ivarloose", &m_electron_match_HLT_e26_lhtight_nod0_ivarloose);
    m_inputNTuple->SetBranchAddress ("electron_match_HLT_e60_lhmedium",		      &m_electron_match_HLT_e60_lhmedium);
    m_inputNTuple->SetBranchAddress ("electron_match_HLT_e60_lhmedium_nod0",	      &m_electron_match_HLT_e60_lhmedium_nod0);
    m_inputNTuple->SetBranchAddress ("electron_match_HLT_e120_lhloose",		      &m_electron_match_HLT_e120_lhloose);
    m_inputNTuple->SetBranchAddress ("electron_match_HLT_e140_lhloose_nod0",	      &m_electron_match_HLT_e140_lhloose_nod0);
    m_inputNTuple->SetBranchAddress ("electron_match_HLT_2e12_lhloose_L12EM10VH",     &m_electron_match_HLT_2e12_lhloose_L12EM10VH);
    m_inputNTuple->SetBranchAddress ("electron_match_HLT_2e17_lhvloose_nod0",	      &m_electron_match_HLT_2e17_lhvloose_nod0);
    m_inputNTuple->SetBranchAddress ("electron_match_HLT_e17_lhloose_mu14",           &m_electron_match_HLT_e17_lhloose_mu14);
    m_inputNTuple->SetBranchAddress ("electron_match_HLT_e17_lhloose_nod0_mu14",      &m_electron_match_HLT_e17_lhloose_nod0_mu14);
    m_inputNTuple->SetBranchAddress ("electron_match_HLT_e7_medium_mu24", 	      &m_electron_match_HLT_e7_medium_mu24);
    m_inputNTuple->SetBranchAddress ("electron_match_HLT_e7_lhmedium_mu24",	      &m_electron_match_HLT_e7_lhmedium_mu24);
    m_inputNTuple->SetBranchAddress ("electron_match_HLT_e24_medium_L1EM20VHI_mu8noL1", &m_electron_match_HLT_e24_medium_L1EM20VHI_mu8noL1);
    m_inputNTuple->SetBranchAddress ("muon_passOR",                                   &m_muon_passOR);
    m_inputNTuple->SetBranchAddress ("muon_pt", 	                              &m_muon_pt);
    m_inputNTuple->SetBranchAddress ("muon_eta",		                      &m_muon_eta);
    m_inputNTuple->SetBranchAddress ("muon_phi",  	                              &m_muon_phi);
    m_inputNTuple->SetBranchAddress ("muon_ID", 	                              &m_muon_ID);
    m_inputNTuple->SetBranchAddress ("muon_sigd0PV",	                              &m_muon_sigd0PV);
    m_inputNTuple->SetBranchAddress ("muon_z0SinTheta",  	                      &m_muon_z0SinTheta);
    m_inputNTuple->SetBranchAddress ("muon_ptvarcone30",  	                      &m_muon_ptvarcone30);
    m_inputNTuple->SetBranchAddress ("muon_truthType",  	                      &m_muon_truthType);
    m_inputNTuple->SetBranchAddress ("muon_truthOrigin",  	                      &m_muon_truthOrigin);
    m_inputNTuple->SetBranchAddress ("muon_PromptLeptonIso_TagWeight",    	      &m_muon_PromptLeptonIso_TagWeight);
    m_inputNTuple->SetBranchAddress ("muon_match_HLT_mu20_iloose_L1MU15", 	      &m_muon_match_HLT_mu20_iloose_L1MU15);
    m_inputNTuple->SetBranchAddress ("muon_match_HLT_mu26_ivarmedium",		      &m_muon_match_HLT_mu26_ivarmedium);
    m_inputNTuple->SetBranchAddress ("muon_match_HLT_mu50",			      &m_muon_match_HLT_mu50);
    m_inputNTuple->SetBranchAddress ("muon_match_HLT_mu22_mu8noL1",		      &m_muon_match_HLT_mu22_mu8noL1);
    m_inputNTuple->SetBranchAddress ("muon_match_HLT_e17_lhloose_mu14",		      &m_muon_match_HLT_e17_lhloose_mu14);
    m_inputNTuple->SetBranchAddress ("muon_match_HLT_e17_lhloose_nod0_mu14",	      &m_muon_match_HLT_e17_lhloose_nod0_mu14);
    m_inputNTuple->SetBranchAddress ("muon_match_HLT_mu18_mu8noL1",		      &m_muon_match_HLT_mu18_mu8noL1);
    m_inputNTuple->SetBranchAddress ("muon_match_HLT_e24_medium_L1EM20VHI_mu8noL1",   &m_muon_match_HLT_e24_medium_L1EM20VHI_mu8noL1);
    m_inputNTuple->SetBranchAddress ("muon_match_HLT_e7_medium_mu24",                 &m_muon_match_HLT_e7_medium_mu24);

    m_inputNTuple->SetBranchAddress ("m_jet_pt",  &m_jet_pt);
    m_inputNTuple->SetBranchAddress ("m_jet_eta", &m_jet_eta);
    m_inputNTuple->SetBranchAddress ("m_jet_phi", &m_jet_phi);
    m_inputNTuple->SetBranchAddress ("m_jet_E",   &m_jet_E);
    m_inputNTuple->SetBranchAddress ("m_jet_flavor_weight_MV2c10",     &m_jet_flavor_weight_MV2c10);
    m_inputNTuple->SetBranchAddress ("m_jet_flavor_truth_label",       &m_jet_flavor_truth_label);
    m_inputNTuple->SetBranchAddress ("m_jet_flavor_truth_label_ghost", &m_jet_flavor_truth_label_ghost);

    m_inputNTuple->SetBranchAddress ("selected_jets",   &m_selected_jets);
    m_inputNTuple->SetBranchAddress ("selected_jets_T", &m_selected_jets_T);

    m_inputNTuple->SetBranchAddress ("m_truth_jet_pt",  &m_truth_jet_pt);
    m_inputNTuple->SetBranchAddress ("m_truth_jet_eta", &m_truth_jet_eta);
    m_inputNTuple->SetBranchAddress ("m_truth_jet_phi", &m_truth_jet_phi);
    m_inputNTuple->SetBranchAddress ("m_truth_jet_e",   &m_truth_jet_e);

    // Get the pointer to the sumWeights TTree

    m_sumWeightsTree = dynamic_cast<TTree*>(wk()->inputFile()->Get("sumWeights"));

    m_sumWeightsTree->SetBranchAddress ("totalEvents", &m_totalEvents);
    m_sumWeightsTree->SetBranchAddress ("totalEventsWeighted", &m_totalEventsWeighted);


    return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepMiniNTupMaker :: initialize ()
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

    Info("initialize()", "Initialising HTopMultilepMiniNTupMaker...");

    m_outputNTuple = EL::getNTupleSvc (wk(), m_outputNTupStreamName);

    // Set new branches for output TTree

    m_outputNTuple->tree()->Branch("weight_event", 	    	&m_weight_event, "weight_event/F");
    m_outputNTuple->tree()->Branch("weight_event_trig_SLT", 	&m_weight_event_trig_SLT, "weight_event_trig_SLT/F");
    m_outputNTuple->tree()->Branch("weight_event_lep", 	    	&m_weight_event_lep, "weight_event_lep/F");

    m_outputNTuple->tree()->Branch("weight_lep_tag",		&m_weight_lep_tag,    "weight_lep_tag/F");
    m_outputNTuple->tree()->Branch("weight_trig_tag",		&m_weight_trig_tag,   "weight_trig_tag/F");
    m_outputNTuple->tree()->Branch("weight_lep_probe", 	    	&m_weight_lep_probe,  "weight_lep_probe/F");
    m_outputNTuple->tree()->Branch("weight_trig_probe", 	&m_weight_trig_probe, "weight_trig_probe/F");

    m_outputNTuple->tree()->Branch("isSS01",               	&m_isSS01, "isSS01/B");
    m_outputNTuple->tree()->Branch("isSS12",               	&m_isSS12, "isSS12/B");

    m_outputNTuple->tree()->Branch("is_T_T",               	&m_is_T_T, "is_T_T/B");
    m_outputNTuple->tree()->Branch("is_T_AntiT",               	&m_is_T_AntiT, "is_T_AntiT/B");
    m_outputNTuple->tree()->Branch("is_AntiT_T",               	&m_is_AntiT_T, "is_AntiT_T/B");
    m_outputNTuple->tree()->Branch("is_AntiT_AntiT",            &m_is_AntiT_AntiT, "is_AntiT_AntiT/B");
    m_outputNTuple->tree()->Branch("is_Tel_AntiTmu",      	&m_is_Tel_AntiTmu, "is_Tel_AntiTmu/B");
    m_outputNTuple->tree()->Branch("is_AntiTel_Tmu",      	&m_is_AntiTel_Tmu, "is_AntiTel_Tmu/B");
    m_outputNTuple->tree()->Branch("is_Tmu_AntiTel",      	&m_is_Tmu_AntiTel, "is_Tmu_AntiTel/B");
    m_outputNTuple->tree()->Branch("is_AntiTmu_Tel",      	&m_is_AntiTmu_Tel, "is_AntiTmu_Tel/B");

    m_outputNTuple->tree()->Branch("is_TMVA_TMVA",              &m_is_TMVA_TMVA, "is_TMVA_TMVA/B");
    m_outputNTuple->tree()->Branch("is_TMVA_AntiTMVA",          &m_is_TMVA_AntiTMVA, "is_TMVA_AntiTMVA/B");
    m_outputNTuple->tree()->Branch("is_AntiTMVA_TMVA",          &m_is_AntiTMVA_TMVA, "is_AntiTMVA_TMVA/B");
    m_outputNTuple->tree()->Branch("is_AntiTMVA_AntiTMVA",      &m_is_AntiTMVA_AntiTMVA, "is_AntiTMVA_AntiTMVA/B");
    m_outputNTuple->tree()->Branch("is_TMVAel_AntiTMVAmu",      &m_is_TMVAel_AntiTMVAmu, "is_TMVAel_AntiTMVAmu/B");
    m_outputNTuple->tree()->Branch("is_AntiTMVAel_TMVAmu",      &m_is_AntiTMVAel_TMVAmu, "is_AntiTMVAel_TMVAmu/B");
    m_outputNTuple->tree()->Branch("is_TMVAmu_AntiTMVAel",      &m_is_TMVAmu_AntiTMVAel, "is_TMVAmu_AntiTMVAel/B");
    m_outputNTuple->tree()->Branch("is_AntiTMVAmu_TMVAel",      &m_is_AntiTMVAmu_TMVAel, "is_AntiTMVAmu_TMVAel/B");

    m_outputNTuple->tree()->Branch("nmuons",               	&m_nmuons, "nmuons/I");
    m_outputNTuple->tree()->Branch("nelectrons",               	&m_nelectrons, "nelectrons/I");
    m_outputNTuple->tree()->Branch("nleptons",               	&m_nleptons, "nleptons/I");

    // Store some additional lepton flat branches

    std::string branchtype, branchname, branchnametype, key;

    for ( const auto& var : m_LEP_OUTPUT_VARS ) {

	key = var.substr( 0, var.length() - 2 );

	branchtype = var.substr( var.length() - 1 );

	for (int idx(0); idx < 3; ++idx ) {

	    branchname = "lep_" + key + "_" + std::to_string(idx);

	    switch( idx ) {
	    case 0:
		if ( branchtype.compare("F") == 0 ) { m_outputNTuple->tree()->Branch ( branchname.c_str(),   &m_lep0_OUTPUT_branches[key].f ); }
		if ( branchtype.compare("B") == 0 ) { m_outputNTuple->tree()->Branch ( branchname.c_str(),   &m_lep0_OUTPUT_branches[key].c ); }
		if ( branchtype.compare("I") == 0 ) { m_outputNTuple->tree()->Branch ( branchname.c_str(),   &m_lep0_OUTPUT_branches[key].i ); }
		break;
	    case 1:
		if ( branchtype.compare("F") == 0 ) { m_outputNTuple->tree()->Branch ( branchname.c_str(),   &m_lep1_OUTPUT_branches[key].f ); }
		if ( branchtype.compare("B") == 0 ) { m_outputNTuple->tree()->Branch ( branchname.c_str(),   &m_lep1_OUTPUT_branches[key].c ); }
		if ( branchtype.compare("I") == 0 ) { m_outputNTuple->tree()->Branch ( branchname.c_str(),   &m_lep1_OUTPUT_branches[key].i ); }
		break;
	    case 2:
		if ( branchtype.compare("F") == 0 ) { m_outputNTuple->tree()->Branch ( branchname.c_str(),   &m_lep2_OUTPUT_branches[key].f ); }
		if ( branchtype.compare("B") == 0 ) { m_outputNTuple->tree()->Branch ( branchname.c_str(),   &m_lep2_OUTPUT_branches[key].c ); }
		if ( branchtype.compare("I") == 0 ) { m_outputNTuple->tree()->Branch ( branchname.c_str(),   &m_lep2_OUTPUT_branches[key].i ); }
		break;
	    default:
		break;
	    }
	}
    }

    // Store some vector branches for electrons and muons passing OLR (pT-ordered)

    for ( const auto& var : m_EL_VEC_VARS ) {

	branchtype = var.substr( var.length() - 1 );
	branchname = "electron_" + var.substr( 0, var.length() - 2 );

	if ( branchtype.compare("F") == 0 ) { m_outputNTuple->tree()->Branch( branchname.c_str(), &( m_electron_OR_branches[branchname].vec_f ) ); }
	if ( branchtype.compare("B") == 0 ) { m_outputNTuple->tree()->Branch( branchname.c_str(), &( m_electron_OR_branches[branchname].vec_c ) ); }
	if ( branchtype.compare("I") == 0 ) { m_outputNTuple->tree()->Branch( branchname.c_str(), &( m_electron_OR_branches[branchname].vec_i ) ); }

    }

    for ( const auto& var : m_MU_VEC_VARS ) {

	branchtype = var.substr( var.length() - 1 );
	branchname = "muon_" + var.substr( 0, var.length() - 2 );

	if ( branchtype.compare("F") == 0 ) { m_outputNTuple->tree()->Branch( branchname.c_str(), &( m_muon_OR_branches[branchname].vec_f ) ); }
	if ( branchtype.compare("B") == 0 ) { m_outputNTuple->tree()->Branch( branchname.c_str(), &( m_muon_OR_branches[branchname].vec_c ) ); }
	if ( branchtype.compare("I") == 0 ) { m_outputNTuple->tree()->Branch( branchname.c_str(), &( m_muon_OR_branches[branchname].vec_i ) ); }

    }

    // From v23 onwards, no need to do this anymore --> already being done in GFW
    //
    //m_outputNTuple->tree()->Branch("lep_isTrigMatch_SLT_0",       &m_lep_isTrigMatch_SLT_0, "lep_isTrigMatch_SLT_0/B");
    //m_outputNTuple->tree()->Branch("lep_isTrigMatch_SLT_1",       &m_lep_isTrigMatch_SLT_1, "lep_isTrigMatch_SLT_1/B");
    //m_outputNTuple->tree()->Branch("lep_isTrigMatch_DLT_0",       &m_lep_isTrigMatch_DLT_0, "lep_isTrigMatch_DLT_0/B");
    //m_outputNTuple->tree()->Branch("lep_isTrigMatch_DLT_1",       &m_lep_isTrigMatch_DLT_1, "lep_isTrigMatch_DLT_1/B");

    m_outputNTuple->tree()->Branch("event_isTrigMatch_DLT",       &m_event_isTrigMatch_DLT, "event_isTrigMatch_DLT/B");

    m_outputNTuple->tree()->Branch("event_isBadTP_SLT", &m_isBadTPEvent_SLT, "event_isBadTP_SLT/B");

    std::string branchname_el, branchname_mu;

    for ( const auto& tp : m_TPS ) {

	for ( const auto& trig : m_TRIGS ) {

	    for ( const auto& var : m_TP_VARS ) {

		bool isvec = ( var.find("VEC") != std::string::npos );

		if ( isvec ) {
		    branchtype     = var.substr( var.length() - 4 );
		    branchname     = "lep_" + tp + "Vec_" + trig + "_" + var.substr( 0, var.length() - 5 );
		    branchname_el  = "electron_" + tp + "Vec_" + trig + "_" + var.substr( 0, var.length() - 5 );
		    branchname_mu  = "muon_" + tp + "Vec_" + trig + "_" + var.substr( 0, var.length() - 5 );
		} else {
		    branchtype     = var.substr( var.length() - 1 );
		    branchname     = "lep_" + tp + "_" + trig + "_" + var.substr( 0, var.length() - 2 );
		    branchnametype = var;
		}

		if      ( branchtype.compare("F") == 0 )    { m_outputNTuple->tree()->Branch( branchname.c_str(), &( m_TagProbe_branches[branchname].f ), branchnametype.c_str() ); }
		else if ( branchtype.compare("B") == 0 )    { m_outputNTuple->tree()->Branch( branchname.c_str(), &( m_TagProbe_branches[branchname].c ), branchnametype.c_str() ); }
		else if ( branchtype.compare("I") == 0 )    { m_outputNTuple->tree()->Branch( branchname.c_str(), &( m_TagProbe_branches[branchname].i ), branchnametype.c_str() ); }
		else if ( branchtype.compare("VECF") == 0 ) {
		    m_outputNTuple->tree()->Branch( branchname.c_str(),    &( m_TagProbe_branches[branchname].vec_f ) );
		    m_outputNTuple->tree()->Branch( branchname_el.c_str(), &( m_TagProbe_branches[branchname_el].vec_f ) );
		    m_outputNTuple->tree()->Branch( branchname_mu.c_str(), &( m_TagProbe_branches[branchname_mu].vec_f ) );
		} else if ( branchtype.compare("VECB") == 0 ) {
		    m_outputNTuple->tree()->Branch( branchname.c_str(),    &( m_TagProbe_branches[branchname].vec_c ) );
		    m_outputNTuple->tree()->Branch( branchname_el.c_str(), &( m_TagProbe_branches[branchname_el].vec_c ) );
		    m_outputNTuple->tree()->Branch( branchname_mu.c_str(), &( m_TagProbe_branches[branchname_mu].vec_c ) );
		} else if ( branchtype.compare("VECI") == 0 ) {
		    m_outputNTuple->tree()->Branch( branchname.c_str(),    &( m_TagProbe_branches[branchname].vec_i ) );
		    m_outputNTuple->tree()->Branch( branchname_el.c_str(), &( m_TagProbe_branches[branchname_el].vec_i ) );
		    m_outputNTuple->tree()->Branch( branchname_mu.c_str(), &( m_TagProbe_branches[branchname_mu].vec_i ) );
		}
	    }

	}

    }

    m_outputNTuple->tree()->Branch("jet_OR_Pt",		   	        &m_jet_OR_Pt);
    m_outputNTuple->tree()->Branch("jet_OR_Eta",			&m_jet_OR_Eta);
    m_outputNTuple->tree()->Branch("jet_OR_Phi",  	   	        &m_jet_OR_Phi);
    m_outputNTuple->tree()->Branch("jet_OR_E",		   	        &m_jet_OR_E);

    if ( m_jetTruthMatching ) {
	m_outputNTuple->tree()->Branch("jet_OR_truthMatch_Pt",		&m_jet_OR_truthMatch_Pt);
	m_outputNTuple->tree()->Branch("jet_OR_truthMatch_Eta",  	&m_jet_OR_truthMatch_Eta);
	m_outputNTuple->tree()->Branch("jet_OR_truthMatch_Phi",		&m_jet_OR_truthMatch_Phi);
	m_outputNTuple->tree()->Branch("jet_OR_truthMatch_E",		&m_jet_OR_truthMatch_E);
	m_outputNTuple->tree()->Branch("jet_OR_truthMatch_isBJet",  	&m_jet_OR_truthMatch_isBJet);
	m_outputNTuple->tree()->Branch("jet_OR_truthMatch_isCJet",	&m_jet_OR_truthMatch_isCJet);
	m_outputNTuple->tree()->Branch("jet_OR_truthMatch_isLFJet",	&m_jet_OR_truthMatch_isLFJet);
	m_outputNTuple->tree()->Branch("jet_OR_truthMatch_isGluonJet",  &m_jet_OR_truthMatch_isGluonJet);
    }

    // ---------------------------------------------------------------------------------------------------------------

    // Initialise counter for input TTree entries processed

    m_numEntry = -1;

    m_effectiveTotEntries = m_inputNTuple->GetEntries();

    unsigned int maxEvents = static_cast<int>( wk()->metaData()->castDouble("nc_EventLoop_MaxEvents") );

    if ( maxEvents > 0 ) {
	m_effectiveTotEntries = maxEvents;
    }

    Info("initialize()", "Name of input TTree : %s", m_inputNTuple->GetName() );
    Info("initialize()", "Total events to run on: %u", m_effectiveTotEntries );

    // ---------------------------------------------------------------------------------------------------------------

    // Print whether we are using an EL::AlgSelect algorithm in our job

    if ( m_useAlgSelect ) { Info("initialize()", "Will be using an EL::AlgSelect algorithm for event skimming" ); }

    // ---------------------------------------------------------------------------------------------------------------

    m_outputNTuple->tree()->SetName( m_outputNTupName.c_str() );

    // ---------------------------------------------------------------------------------------------------------------

    if ( m_useTruthTP ) {
	std::cout << "\n" << std::endl;
	Warning("initialize()", "Will define T&P leptons based on truth. Check that you're running on TTBar (nonallhad) only!!" );
	std::cout << "\n" << std::endl;
    }

    // ---------------------------------------------------------------------------------------------------------------

    m_rand = new TRandom3(0); // "0" guarantees unique seed in space & time

    Info("initialize()", "All good! Start processing now...");

    // ---------------------------------------------------------------------------------------------------------------

    return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepMiniNTupMaker :: execute ()
{
    // Here you do everything that needs to be done on every single
    // events, e.g. read input variables, apply cuts, and fill
    // histograms and trees.  This is where most of your actual analysis
    // code will go.

    ANA_CHECK_SET_TYPE (EL::StatusCode);

    ++m_numEntry;

    if ( m_numEntry < m_sumWeightsTree->GetEntries() ) {

	m_sumWeightsTree->GetEntry( m_numEntry );

	m_sumGenEvents         += m_totalEvents;
	m_sumGenEventsWeighted += m_totalEventsWeighted;

    }

    m_inputNTuple->GetEntry( wk()->treeEntry() );

    if ( m_debug ) {
	std::cout << "" << std::endl;
	Info("execute()", "===> Entry %u - EventNumber = %u - RunYear = %i", static_cast<uint32_t>(m_numEntry), static_cast<uint32_t>(m_EventNumber), m_RunYear );
	std::cout << "" << std::endl;
    }

    if ( m_numEntry >= 10000 && m_numEntry % 10000 == 0 ) {
	std::cout << "Processed " << std::setprecision(3) << ( (float) m_numEntry / m_effectiveTotEntries ) * 1e2 << " % of total entries" << std::endl;
    }

    // ------------------------------------------------------------------------

    if ( !m_useAlgSelect ) {

	m_outputNTuple->setFilterPassed();

	if ( !m_passEventCleaning )			      { wk()->skipEvent(); return EL::StatusCode::SUCCESS; }
	//if ( m_nJets_OR_T < 1 ) 			      { wk()->skipEvent(); return EL::StatusCode::SUCCESS; }
	if ( !( m_dilep_type > 0 || m_trilep_type > 0 ) ) { wk()->skipEvent(); return EL::StatusCode::SUCCESS; }

    }

    // ------------------------------------------------------------------------

    m_event = std::make_shared<eventObj>();

    m_event.get()->dilep_type   = m_dilep_type;
    m_event.get()->trilep_type  = m_trilep_type;
    m_event.get()->quadlep_type = m_quadlep_type;

    // This is crucial to avoid methods called later to crash: DO NOT REMOVE!

    if ( !m_event.get()->dilep_type && !m_event.get()->trilep_type && !m_event.get()->quadlep_type ) { return EL::StatusCode::SUCCESS; }

    // ------------------------------------------------------------------------

    auto lep0 = std::make_shared<leptonObj>();

    std::string key;

    for ( auto& property : lep0.get()->props ) {

	key  = property.first;

	if ( key.compare("Flavour") == 0 )        { property.second.i = std::abs(m_lep0_INPUT_branches["ID"].f); continue; }
	if ( key.compare("Charge") == 0 )         { property.second.f = m_lep0_INPUT_branches["ID"].f / std::abs(m_lep0_INPUT_branches["ID"].f); continue; }
	if ( key.compare("PID") == 0 )            { property.second.c = m_lep0_INPUT_branches["isTightLH"].c; continue; }
	if ( key.compare("LooseIsolated") == 0 )  { property.second.c = m_lep0_INPUT_branches["isolationLoose"].i; continue; }
	if ( key.compare("Isolated") == 0 )       { property.second.c = ( std::abs(m_lep0_INPUT_branches["ID"].f) == 13 ) ? m_lep0_INPUT_branches["isolationFixedCutTightTrackOnly"].i : m_lep0_INPUT_branches["isolationFixedCutTight"].i; continue; }
	if ( key.compare("TrackIsoOverPt") == 0 ) { property.second.f = ( std::abs(m_lep0_INPUT_branches["ID"].f) == 13 ) ? m_lep0_INPUT_branches["ptKeycone30"].f / m_lep0_INPUT_branches["Pt"].f : m_lep0_INPUT_branches["ptKeycone20"].f / m_lep0_INPUT_branches["Pt"].f; continue; }
	if ( key.compare("CaloIsoOverPt") == 0 )  { property.second.f = ( std::abs(m_lep0_INPUT_branches["ID"].f) == 13 ) ? -1.0 : m_lep0_INPUT_branches["topoEtcone20"].f / m_lep0_INPUT_branches["Pt"].f; continue; }
	if ( key.compare("isQMisID") == 0 )       { property.second.c = m_lep0_INPUT_branches[key].c; continue; }

	if ( m_lep0_INPUT_branches[key].f != -999.0 ) { property.second.f = m_lep0_INPUT_branches[key].f; }
	if ( m_lep0_INPUT_branches[key].c != -1 )     { property.second.c = m_lep0_INPUT_branches[key].c; }
	if ( m_lep0_INPUT_branches[key].i != -999 )   { property.second.i = m_lep0_INPUT_branches[key].i; }

    }

    // Check if this lepton passes tight selection
    ANA_CHECK( this->checkIsTightLep( lep0 ) );
    ANA_CHECK( this->checkIsTightLep( lep0, "MVA" ) );

    if ( m_debug ) {
	std::cout << "Lep0 properties: " << std::endl;
	size_t maxlength(40), nblanks;
	for ( const auto& property : lep0.get()->props ) {
	    nblanks = ( maxlength > property.first.size() ) ?  maxlength - property.first.size() : 0;
	    std::string blankstr(nblanks,' ');
	    std::cout << "\t" << property.first << " = " << blankstr <<  std::setw(20) << property.second.f << " (F), " <<  std::setw(15) << (int)property.second.c << " (B), " <<  std::setw(15) << property.second.i << " (I) " << std::endl;
	}
    }

    m_leptons.push_back(lep0);

    auto lep1 = std::make_shared<leptonObj>();

    for ( auto& property : lep1.get()->props ) {

	key  = property.first;

	if ( key.compare("Flavour") == 0 )        { property.second.i = std::abs(m_lep1_INPUT_branches["ID"].f); continue; }
	if ( key.compare("Charge") == 0 )         { property.second.f = m_lep1_INPUT_branches["ID"].f / std::abs(m_lep1_INPUT_branches["ID"].f); continue; }
	if ( key.compare("PID") == 0 )            { property.second.c = m_lep1_INPUT_branches["isTightLH"].c; continue; }
	if ( key.compare("LooseIsolated") == 0 )  { property.second.c = m_lep1_INPUT_branches["isolationLoose"].i; continue; }
	if ( key.compare("Isolated") == 0 )       { property.second.c = ( std::abs(m_lep1_INPUT_branches["ID"].f) == 13 ) ? m_lep1_INPUT_branches["isolationFixedCutTightTrackOnly"].i : m_lep1_INPUT_branches["isolationFixedCutTight"].i; continue; }
	if ( key.compare("TrackIsoOverPt") == 0 ) { property.second.f = ( std::abs(m_lep1_INPUT_branches["ID"].f) == 13 ) ? m_lep1_INPUT_branches["ptKeycone30"].f / m_lep1_INPUT_branches["Pt"].f : m_lep1_INPUT_branches["ptKeycone20"].f / m_lep1_INPUT_branches["Pt"].f; continue; }
	if ( key.compare("CaloIsoOverPt") == 0 )  { property.second.f = ( std::abs(m_lep1_INPUT_branches["ID"].f) == 13 ) ? -1.0 : m_lep1_INPUT_branches["topoEtcone20"].f / m_lep1_INPUT_branches["Pt"].f; continue; }
	if ( key.compare("isQMisID") == 0 )       { property.second.c = m_lep1_INPUT_branches[key].c; continue; }

	if ( m_lep1_INPUT_branches[key].f != -999.0 ) { property.second.f = m_lep1_INPUT_branches[key].f; }
	if ( m_lep1_INPUT_branches[key].c != -1 )     { property.second.c = m_lep1_INPUT_branches[key].c; }
	if ( m_lep1_INPUT_branches[key].i != -999 )     { property.second.i = m_lep1_INPUT_branches[key].i; }

    }

    // Check if this lepton passes tight selection
    ANA_CHECK( this->checkIsTightLep( lep1 ) );
    ANA_CHECK( this->checkIsTightLep( lep1, "MVA" ) );

    m_leptons.push_back(lep1);

    if ( m_debug ) {
	Info("execute()","lep0 pT = %.2f - flavour = %i", lep0.get()->props["Pt"].f/1e3, lep0.get()->props["Flavour"].i );
	Info("execute()","lep1 pT = %.2f - flavour = %i", lep1.get()->props["Pt"].f/1e3, lep1.get()->props["Flavour"].i );
    }

    if ( m_event.get()->trilep_type || m_event.get()->quadlep_type) {

	auto lep2 = std::make_shared<leptonObj>();

	for ( auto& property : lep2.get()->props ) {

	    key  = property.first;

	    if ( key.compare("Flavour") == 0 )        { property.second.i = std::abs(m_lep2_INPUT_branches["ID"].f); continue; }
	    if ( key.compare("Charge") == 0 )         { property.second.f = m_lep2_INPUT_branches["ID"].f / std::abs(m_lep2_INPUT_branches["ID"].f); continue; }
	    if ( key.compare("PID") == 0 )            { property.second.c = m_lep2_INPUT_branches["isTightLH"].c; continue; }
	    if ( key.compare("LooseIsolated") == 0 )  { property.second.c = m_lep2_INPUT_branches["isolationLoose"].i; continue; }
	    if ( key.compare("Isolated") == 0 )       { property.second.c = ( std::abs(m_lep2_INPUT_branches["ID"].f) == 13 ) ? m_lep2_INPUT_branches["isolationFixedCutTightTrackOnly"].i : m_lep2_INPUT_branches["isolationFixedCutTight"].i; continue; }
	    if ( key.compare("TrackIsoOverPt") == 0 ) { property.second.f = ( std::abs(m_lep2_INPUT_branches["ID"].f) == 13 ) ? m_lep2_INPUT_branches["ptKeycone30"].f / m_lep2_INPUT_branches["Pt"].f : m_lep2_INPUT_branches["ptKeycone20"].f / m_lep2_INPUT_branches["Pt"].f; continue; }
	    if ( key.compare("CaloIsoOverPt") == 0 )  { property.second.f = ( std::abs(m_lep2_INPUT_branches["ID"].f) == 13 ) ? -1.0 : m_lep2_INPUT_branches["topoEtcone20"].f / m_lep2_INPUT_branches["Pt"].f; continue; }
	    if ( key.compare("isQMisID") == 0 )       { property.second.c = m_lep2_INPUT_branches[key].c; continue; }

	    if ( m_lep2_INPUT_branches[key].f != -999.0 ) { property.second.f = m_lep2_INPUT_branches[key].f; }
	    if ( m_lep2_INPUT_branches[key].c != -1 )     { property.second.c = m_lep2_INPUT_branches[key].c; }
	    if ( m_lep2_INPUT_branches[key].i != -999 )   { property.second.i = m_lep2_INPUT_branches[key].i; }

	}

	// Check if this lepton passes tight selection

	ANA_CHECK( this->checkIsTightLep( lep2 ) );
	ANA_CHECK( this->checkIsTightLep( lep2, "MVA" ) );

	if ( m_debug ) { Info("execute()","lep2 pT = %.2f - flavour = %i", lep2.get()->props["Pt"].f/1e3, lep2.get()->props["Flavour"].i ); }

	m_leptons.push_back(lep2);

	if ( m_event.get()->quadlep_type ) {

	    auto lep3 = std::make_shared<leptonObj>();

	    for ( auto& property : lep3.get()->props ) {

		key  = property.first;

		if ( key.compare("Flavour") == 0 )        { property.second.i = std::abs(m_lep3_INPUT_branches["ID"].f); continue; }
		if ( key.compare("Charge") == 0 )         { property.second.f = m_lep3_INPUT_branches["ID"].f / std::abs(m_lep3_INPUT_branches["ID"].f); continue; }
		if ( key.compare("PID") == 0 )            { property.second.c = m_lep3_INPUT_branches["isTightLH"].c; continue; }
		if ( key.compare("LooseIsolated") == 0 )  { property.second.c = m_lep3_INPUT_branches["isolationLoose"].i; continue; }
		if ( key.compare("Isolated") == 0 )       { property.second.c = ( std::abs(m_lep3_INPUT_branches["ID"].f) == 13 ) ? m_lep3_INPUT_branches["isolationFixedCutTightTrackOnly"].i : m_lep3_INPUT_branches["isolationFixedCutTight"].i; continue; }
		if ( key.compare("TrackIsoOverPt") == 0 ) { property.second.f = ( std::abs(m_lep3_INPUT_branches["ID"].f) == 13 ) ? m_lep3_INPUT_branches["ptKeycone30"].f / m_lep3_INPUT_branches["Pt"].f : m_lep3_INPUT_branches["ptKeycone20"].f / m_lep3_INPUT_branches["Pt"].f; continue; }
		if ( key.compare("CaloIsoOverPt") == 0 )  { property.second.f = ( std::abs(m_lep3_INPUT_branches["ID"].f) == 13 ) ? -1.0 : m_lep3_INPUT_branches["topoEtcone20"].f / m_lep3_INPUT_branches["Pt"].f; continue; }
		if ( key.compare("isQMisID") == 0 )       { property.second.c = m_lep3_INPUT_branches[key].c; continue; }

		if ( m_lep3_INPUT_branches[key].f != -999.0 ) { property.second.f = m_lep3_INPUT_branches[key].f; }
		if ( m_lep3_INPUT_branches[key].c != -1 )     { property.second.c = m_lep3_INPUT_branches[key].c; }
		if ( m_lep3_INPUT_branches[key].i != -999 )   { property.second.i = m_lep3_INPUT_branches[key].i; }

	    }

	    // Check if this lepton passes tight selection

	    ANA_CHECK( this->checkIsTightLep( lep3 ) );
	    ANA_CHECK( this->checkIsTightLep( lep3, "MVA" ) );

	    if ( m_debug ) { Info("execute()","lep3 pT = %.3f - flavour = %i", lep3.get()->props["Pt"].f/1e3, lep3.get()->props["Flavour"].i ); }

	    m_leptons.push_back(lep3);
	}

    }

    // ------------------------------------------------------------------------

    ANA_CHECK( this->decorateEvent() );

    // ------------------------------------------------------------------------

    ANA_CHECK( this->defineTagAndProbe() );

    // ------------------------------------------------------------------------

    ANA_CHECK( this->decorateWeights() );

    // ------------------------------------------------------------------------

    ANA_CHECK( this->jetKinematics() );

    // ------------------------------------------------------------------------

    if ( m_jetTruthMatching ) { ANA_CHECK( this->jetTruthMatching() ); }

    // ------------------------------------------------------------------------

    ANA_CHECK( this->setOutputBranches() );

    // ------------------------------------------------------------------------

    // Don't forget to clear object containers before moving to next event!

    m_leptons.clear();
    m_jets.clear();

    return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepMiniNTupMaker :: postExecute ()
{
    // Here you do everything that needs to be done after the main event
    // processing.  This is typically very rare, particularly in user
    // code.  It is mainly used in implementing the NTupleSvc.

    //ANA_CHECK_SET_TYPE (EL::StatusCode);

    return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepMiniNTupMaker :: finalize ()
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

    Info("finalize()", "Finalising HTopMultilepMiniNTupMaker...");

    delete m_rand; m_rand = nullptr;

    return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepMiniNTupMaker :: histFinalize ()
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

    m_sumGenEventsHist->SetBinContent( 1, m_sumGenEvents );
    m_sumGenEventsHist->SetBinContent( 2, m_sumGenEventsWeighted );

    Info("histFinalize()", "sumGenEvents = %f", m_sumGenEvents );
    Info("histFinalize()", "sumGenEventsWeighted = %f", m_sumGenEventsWeighted );

    return EL::StatusCode::SUCCESS;
}

EL::StatusCode HTopMultilepMiniNTupMaker :: enableSelectedBranches ()
{

    if ( m_inputBranches.empty() ) {
	Info("enableSelectedBranches()", "Keeping all input branches enabled...");
	return EL::StatusCode::SUCCESS;
    }

    // Firstly, disable all branches

    m_inputNTuple->SetBranchStatus ("*", 0);

    std::vector<std::string> branch_vec;

    // Parse input list, split by comma, and put into a vector

    std::string token;
    std::istringstream ss( m_inputBranches );
    while ( std::getline(ss, token, ',') ) { branch_vec.push_back(token); }

    // Re-enable only the branches we are going to use

    TObjArray* arraybranches = m_inputNTuple->GetListOfBranches();

    for ( const auto& branch : branch_vec ) {

	TString rgxp_branch(branch);
	TRegexp pattern(rgxp_branch);

	// If this branch name's regexp matches a name in the list of branches, then activate the branch

	for ( int idx(0); idx < arraybranches->GetEntries(); ++idx ) {
	    TString str( arraybranches->At(idx)->GetName() );
	    Ssiz_t len(0);
	    if ( str.Index(pattern,&len) != kNPOS && len == str.Length() ) {
		m_inputNTuple->SetBranchStatus (str.Data(), 1);
	    }
	}
    }

    return EL::StatusCode::SUCCESS;

}

EL::StatusCode HTopMultilepMiniNTupMaker ::  checkIsTightLep( std::shared_ptr<leptonObj> lep, const std::string& strategy )
{

    bool useMVA = ( strategy.compare("MVA") == 0 );

    bool isTight(false);

    switch ( lep.get()->props["Flavour"].i )
    {
    case 11:
	if ( useMVA ) {
	    isTight = ( lep.get()->props["LooseIsolated"].c &&
			lep.get()->props["PID"].c &&
			std::abs(lep.get()->props["sigd0PV"].f) < 5.0 &&
			std::abs(lep.get()->props["Z0SinTheta"].f) < 0.5 &&
			lep.get()->props["chargeIDBDTTight"].f > 0.0670415 &&
			lep.get()->props["promptLeptonIso_TagWeight"].f < -0.50 );
	} else {
	    isTight = ( lep.get()->props["Isolated"].c &&
			lep.get()->props["PID"].c &&
			std::abs(lep.get()->props["sigd0PV"].f) < 5.0 &&
			std::abs(lep.get()->props["Z0SinTheta"].f) < 0.5 );
	}
	break;
    case 13:
	if ( useMVA ) {
	    isTight = ( lep.get()->props["LooseIsolated"].c &&
			std::abs(lep.get()->props["sigd0PV"].f) < 3.0 &&
			std::abs(lep.get()->props["Z0SinTheta"].f) < 0.5 &&
			lep.get()->props["promptLeptonIso_TagWeight"].f < -0.50 );
	} else {
	    isTight = ( lep.get()->props["Isolated"].c &&
			std::abs(lep.get()->props["sigd0PV"].f) < 3.0 &&
			std::abs(lep.get()->props["Z0SinTheta"].f) < 0.5 );
	}
	break;
    default:
	break;
    }

    if ( useMVA ) { lep.get()->props["isTightSelectedMVA"].c = isTight; }
    else          { lep.get()->props["isTightSelected"].c = isTight; }

    return EL::StatusCode::SUCCESS;
}


EL::StatusCode HTopMultilepMiniNTupMaker :: decorateEvent ( )
{

    if ( m_event.get()->dilep_type ) {

	int prod_lep_charge(1);
	for ( const auto& lep : m_leptons ) { prod_lep_charge *= lep.get()->props["Charge"].f; }

	m_event.get()->isSS01 = ( prod_lep_charge > 0 );

	// Fill flat branches with trigger matching decision split in SLT/DLT, and per-event (pair) matching for DLT

	// From v23 onwards, no need to do this anymore --> already being done in GFW
	//ANA_CHECK( this->triggerMatching() );

	m_event_isTrigMatch_DLT = ( m_leptons.at(0).get()->props["isTrigMatchDLT"].c && m_leptons.at(1).get()->props["isTrigMatchDLT"].c );

    } else if ( m_event.get()->trilep_type ) {

	m_event.get()->isSS12 = ( fabs( m_leptons.at(0).get()->props["Charge"].f + m_leptons.at(1).get()->props["Charge"].f + m_leptons.at(2).get()->props["Charge"].f ) != 3 );

    }

    // No need for the following since v27:
    //
    // Flatten lepton vars which are not flat yet in the GN
    //
    // ANA_CHECK( this->flatLepVars() );
    //
    // Set the lepton MVA tight WP
    //
    // for ( auto lep : m_leptons ) { ANA_CHECK( this->checkIsTightLep( lep, "MVA" ) ); }

    // Store jet objects w/ some properties (e.g., isbtag). For dilepton events, use jets after OLR w/ taus.

    unsigned int njets = ( m_event.get()->dilep_type ) ? m_selected_jets_T->size() : m_selected_jets->size();

    for ( unsigned int j_idx(0); j_idx < njets; ++j_idx ) {

	short j = ( m_event.get()->dilep_type ) ? m_selected_jets_T->at(j_idx) : m_selected_jets->at(j_idx);

	auto jet = std::make_shared<jetObj>();

	jet.get()->pt     = m_jet_pt->at(j);
	jet.get()->eta    = m_jet_eta->at(j);
	jet.get()->phi    = m_jet_phi->at(j);
	jet.get()->isbtag = !( m_jet_flavor_weight_MV2c10->at(j) < 0.8244273 ); // MV2c10, FixedCutBEff_70 : https://twiki.cern.ch/twiki/bin/view/AtlasProtected/BTaggingBenchmarks#MV2c20_tagger

	m_jets.push_back(jet);

    }

    m_event.get()->njets  = njets;
    m_event.get()->nbjets = ( m_event.get()->dilep_type ) ? m_nJets_OR_T_MV2c10_70 : m_nJets_OR_MV2c10_70;

    if ( m_debug ) { Info("decorateEvent()","Number of jets = %i - Number of bjets (MV2c10, 70 WP) = %i", m_event.get()->njets, m_event.get()->nbjets); }

    // Find closest jet/b-tagged jet to each lepton

    ANA_CHECK( this->findClosestJetLep() );
    ANA_CHECK( this->findClosestJetLep("bjets") );

    return EL::StatusCode::SUCCESS;

}

EL::StatusCode HTopMultilepMiniNTupMaker :: findClosestJetLep( const std::string& jetCollection ) {

    unsigned int idx_lep(0);

    int lep_flavour(0);
    float lep_pt(-1.0), lep_eta(-999.0), lep_phi(-999.0), lep_m(-1.0);

    TLorentzVector lepTLV, jetTLV, closestjetTLV, lepjetTLV;

    bool useBJets =  ( jetCollection.compare("bjets") == 0 );
    std::string jet_flag = ( useBJets ) ? "b" : "";

    for ( auto lep : m_leptons ) {

	lep_flavour = lep.get()->props["Flavour"].i;
	lep_pt      = lep.get()->props["Pt"].f;
	lep_eta     = ( lep_flavour == 11 ) ? lep.get()->props["EtaBE2"].f : lep.get()->props["Eta"].f;
	lep_phi     = lep.get()->props["Phi"].f;
	lep_m       = ( lep_flavour == 11 ) ? 0.511 : 105.65;

	lepTLV.SetPtEtaPhiM( lep_pt, lep_eta, lep_phi, lep_m );

	if ( m_debug ) { Info("findClosestJetLep()","Checking lepton[%i] w/ flavour = %i, pT = %.2f [GeV], eta = %.2f, phi = %.2f, m = %.3f [GeV]", idx_lep, lep_flavour, lep_pt/1e3, lep_eta, lep_phi, lep_m/1e3 ); }

	// Now find closest jet to this lepton

	unsigned int idx_j(0);

	int idx_closest_j(-1);

	float jet_pt(-.0),  jet_eta(-999.0), jet_phi(-999.0);
	float jet_m = ( useBJets ) ? 418.0 : 0.0; // use mass of b quark [MeV]

        float dist_to_closest_j(9e9), this_dist(-1.0);
	float mass_lep_closest_j(-1.0);

	for ( auto jet : m_jets ) {

	    if ( useBJets && !jet.get()->isbtag ) { continue; } // Make sure to look only at bjets if requested.

	    jet_pt  = jet.get()->pt;
	    jet_eta = jet.get()->eta;
	    jet_phi = jet.get()->phi;

	    jetTLV.SetPtEtaPhiM( jet_pt, jet_eta, jet_phi, jet_m );

	    this_dist = lepTLV.DeltaR(jetTLV);

	    if ( m_debug ) { Info("findClosestJetLep()","\t %sjet[%i] w/ pT = %.2f [GeV], eta = %.2f, phi = %.2f, DeltaR(lep[%i],%sjet[%i]) = %.3f", jet_flag.c_str(), idx_j, jet_pt/1e3, jet_eta, jet_phi, idx_lep, jet_flag.c_str(), idx_j, this_dist ); }

	    if ( this_dist < dist_to_closest_j ) {
		dist_to_closest_j = this_dist;
		idx_closest_j     = idx_j;
	    }

	    ++idx_j;

	}

	// Get invariant mass of lepton and closest jet

	if ( m_event.get()->njets > 0 && idx_closest_j > -1 ) {

	    closestjetTLV.SetPtEtaPhiM( m_jets.at(idx_closest_j).get()->pt, m_jets.at(idx_closest_j).get()->eta, m_jets.at(idx_closest_j).get()->phi, 418.0  );

	    lepjetTLV = lepTLV + closestjetTLV;

	    mass_lep_closest_j = lepjetTLV.M();

	}

	// Decorate lepton (Set a dummy negative decorator if there are no jets/bjets in event)

	if ( useBJets ) {
	    lep.get()->props["massClosestBJet"].f   = mass_lep_closest_j;
	    lep.get()->props["deltaRClosestBJet"].f = ( m_event.get()->nbjets > 0 ) ? dist_to_closest_j : -1.0;
	} else {
	    lep.get()->props["massClosestJet"].f   = mass_lep_closest_j;
	    lep.get()->props["deltaRClosestJet"].f = ( m_event.get()->njets > 0 ) ? dist_to_closest_j : -1.0;
	}

	if ( m_debug ) { Info("findClosestJetLep()","==> DeltaR(lep[%i],closest %sjet[%i]) = %.3f, M(lep[%i],closest %sjet[%i]) = %.3f [GeV]", idx_lep, jet_flag.c_str(), idx_closest_j, dist_to_closest_j, idx_lep, jet_flag.c_str(), idx_closest_j, mass_lep_closest_j/1e3 ); }

	++idx_lep;
    }

    return EL::StatusCode::SUCCESS;

}


EL::StatusCode  HTopMultilepMiniNTupMaker :: getPostOLRIndex( int& idx, const unsigned int& pos, const std::string& lep_type ) {

    unsigned int idxOLR(0);
    char passOLR(0);

    if ( m_debug ) { std::cout << "\nlepton type: " << lep_type << "\nlooking at lepton ranked: " << pos << "\n" << std::endl; }

    if ( lep_type.compare("electron") == 0 ) {

	for ( unsigned int e(0); e < m_electron_passOR->size(); ++e ) {

	    passOLR = m_electron_passOR->at(e);

	    if ( m_debug ) { std::cout << "\telectron - idx: " << e << ", pasOLR? " << static_cast<int>(passOLR) << std::endl; }

	    if ( passOLR ) {

		if ( idxOLR == pos ) {
		    idx = e;
		    if ( m_debug ) { std::cout << "\t==> found! " << pos << "-th electron which pass OLR is at position " << idx << " in vector" << std::endl; }
		    return EL::StatusCode::SUCCESS;
		}
		++idxOLR;

	    }

	}

    } else if ( lep_type.compare("muon") == 0 ) {

	for ( unsigned int m(0); m < m_muon_passOR->size(); ++m ) {

	    passOLR = m_muon_passOR->at(m);

	    if ( m_debug ) { std::cout << "\tmuon - idx: " << m << ", pasOLR? " << static_cast<int>(passOLR) << std::endl; }

	    if ( passOLR ) {

		if ( idxOLR == pos ) {
		    idx = m;
		    if ( m_debug ) { std::cout << "\t==> found! " << pos << "-th muon which pass OLR is at position " << idx << " in vector" << std::endl; }
		    return EL::StatusCode::SUCCESS;
		}
		++idxOLR;

	    }

	}

    }

    if ( idx < 0 ) {
	Error("getPostOLRIndex()","Index of %i-th %s passing OLR was not found. Aborting.", pos, lep_type.c_str() );
	return EL::StatusCode::FAILURE;
    }

    return EL::StatusCode::SUCCESS;

}

EL::StatusCode HTopMultilepMiniNTupMaker :: flatLepVars()
{

    if ( m_debug ) { Info("flatLepVars()","Getting lepton flat branches from electron/muon vector branches..." ); }

    // 0.
    // Create a vector of lepton objects

    std::vector< std::shared_ptr<MiniNTupMaker::leptonObj> > allleptons;

    char passOLR(0);

    // 1.
    // Fill the vector of lepton objects w/ info taken from the electron/muon vector branches.
    // Make sure only leptons that pass the OLR are used.
    // Store the pt (needed to sort afterwards!), and all the branches of interest that we wish to "flatten"

    if ( m_debug ) { std::cout << "Electrons:" << std::endl; }
    for ( unsigned int e(0); e < m_electron_passOR->size(); ++e ) { // electrons

	passOLR = m_electron_passOR->at(e);
	if ( !passOLR ) {
	    if ( m_debug ) { std::cout << "\tpassOLR[" << e << "] = " << static_cast<int>(passOLR) << " - skip..." << std::endl; }
	    continue;
	}

	auto lep = std::make_shared<leptonObj>();

	lep.get()->props["Pt"].f                        = m_electron_pt->at(e);
	lep.get()->props["EtaBE2"].f                    = m_electron_EtaBE2->at(e);
	lep.get()->props["Flavour"].i                   = fabs(m_electron_ID->at(e));
	lep.get()->props["Charge"].f                    = m_electron_ID->at(e) / fabs(m_electron_ID->at(e));
	lep.get()->props["promptLeptonIso_TagWeight"].f = m_electron_PromptLeptonIso_TagWeight->at(e);
	lep.get()->props["chargeIDBDTLoose"].f          = m_electron_ChargeIDBDTLoose->at(e);
	lep.get()->props["chargeIDBDTMedium"].f         = m_electron_ChargeIDBDTMedium->at(e);
	lep.get()->props["chargeIDBDTTight"].f          = m_electron_ChargeIDBDTTight->at(e);

	if ( m_debug ) { std::cout << "\tpT[" << e << "] = " <<  lep.get()->props["Pt"].f/1e3 << std::endl; }

	allleptons.push_back(lep);
    }
    if ( m_debug ) { std::cout << "Muons:" << std::endl; }
    for ( unsigned int m(0); m < m_muon_passOR->size(); ++m ) { // muons

	passOLR = m_muon_passOR->at(m);
	if ( !passOLR ) {
	    if ( m_debug ) { std::cout << "\tpassOLR[" << m << "] = " << static_cast<int>(passOLR) << " - skip..." << std::endl; }
	    continue;
	}

	auto lep = std::make_shared<leptonObj>();

	lep.get()->props["Pt"].f                        = m_muon_pt->at(m);
	lep.get()->props["Eta"].f                       = m_muon_eta->at(m);
	lep.get()->props["Flavour"].i                   = fabs(m_muon_ID->at(m));
	lep.get()->props["Charge"].f                    = m_muon_ID->at(m) / fabs(m_muon_ID->at(m));
	lep.get()->props["promptLeptonIso_TagWeight"].f = m_muon_PromptLeptonIso_TagWeight->at(m);

	if ( m_debug ) { std::cout << "\tpT[" << m << "] = " <<  lep.get()->props["Pt"].f/1e3 << std::endl; }

	allleptons.push_back(lep);
    }


    if ( m_event.get()->dilep_type || m_event.get()->quadlep_type ) {

	// 2.
	// Sort the vector of leptons by pT

	std::sort( allleptons.begin(), allleptons.end(), MiniNTupMaker::SorterPt() );

	if ( m_debug ) {
	    std::cout << "pT sorted lepton container:" << std::endl;
	    unsigned int idx(0);
	    for ( auto lep : allleptons ) {
		std::cout << "\tpT[" << idx << "] = " <<  lep.get()->props["Pt"].f/1e3 << " - flavour = " << lep.get()->props["Flavour"].i << std::endl;
		++idx;
	    }
	}

    } else if ( m_event.get()->trilep_type ) {

	// 2.
	// Set the 3L ranking, and then sort the vector of leptons by the 3L ranking itself.

	ANA_CHECK( this->classify3L(allleptons) );

	std::sort( allleptons.begin(), allleptons.end(), MiniNTupMaker::Sorter3L() );

	if ( m_debug ) {
	    std::cout << "3L ranking idx sorted lepton container:" << std::endl;
	    unsigned int idx(0);
	    for ( auto lep : allleptons ) {
		std::cout << "\trank = " << lep.get()->props["Rank3L"].i << " - pT[" << idx << "] = " <<  lep.get()->props["Pt"].f/1e3 << " - flavour = " << lep.get()->props["Flavour"].i << std::endl;
		++idx;
	    }
	}

    }

    // 3.
    // Finally, save the relevant information into the global m_lepton conatiner

    if ( m_debug ) { Info("flatLepVars()","Flattening additional variables:" ); }
    for ( unsigned int idx(0); idx < m_leptons.size(); ++idx ) {

	if ( m_debug ) {
	    std::cout << "\tpT["                        << idx << "] = " << allleptons.at(idx).get()->props["Pt"].f/1e3 << " - (NTup pT = " << m_leptons.at(idx).get()->props["Pt"].f/1e3 << " )" << std::endl;
	    std::cout << "\tpromptLeptonIso_TagWeight[" << idx << "] = " << allleptons.at(idx).get()->props["promptLeptonIso_TagWeight"].f << std::endl;
	    std::cout << "\tchargeIDBDTLoose["          << idx << "] = " << allleptons.at(idx).get()->props["chargeIDBDTLoose"].f          << std::endl;
	    std::cout << "\tchargeIDBDTMedium["         << idx << "] = " << allleptons.at(idx).get()->props["chargeIDBDTMedium"].f         << std::endl;
	    std::cout << "\tchargeIDBDTTight["          << idx << "] = " << allleptons.at(idx).get()->props["chargeIDBDTTight"].f          << std::endl;
	}

	m_leptons.at(idx).get()->props["promptLeptonIso_TagWeight"].f = allleptons.at(idx).get()->props["promptLeptonIso_TagWeight"].f;
	m_leptons.at(idx).get()->props["chargeIDBDTLoose"].f          = allleptons.at(idx).get()->props["chargeIDBDTLoose"].f;
	m_leptons.at(idx).get()->props["chargeIDBDTMedium"].f         = allleptons.at(idx).get()->props["chargeIDBDTMedium"].f;
	m_leptons.at(idx).get()->props["chargeIDBDTTight"].f          = allleptons.at(idx).get()->props["chargeIDBDTTight"].f;
    }

    allleptons.clear();

    return EL::StatusCode::SUCCESS;

}

EL::StatusCode HTopMultilepMiniNTupMaker :: classify3L( std::vector< std::shared_ptr<MiniNTupMaker::leptonObj> >& leptons )
{

    if ( leptons.size() != 3 ) {
	Warning("Sorter3L()","Trying to apply 3L lepton sorting to an event which has %i leptons. Will do nothing...", static_cast<int>(leptons.size()) );
	return EL::StatusCode::SUCCESS;
    }

    int sumcharge(0);
    for ( auto lep : leptons ) { sumcharge += lep.get()->props["Charge"].f; }

    if ( abs(sumcharge) != 1 ) {
	if ( m_debug ) { Info("Sorter3L()","3L event has total charge = %i. Returning...", sumcharge ); }
	return EL::StatusCode::SUCCESS;
    }

    // 1.
    // Find the OS lepton
    // NB. The total charge for 3L SR will be +-1, always. Hence lep0 is the one which is OS to the total charge sign.

    TLorentzVector lep0TLV, lepTLV;
    float eta0(-999.0), eta(-999.0);

    for ( auto& lep : leptons ) {
	if  ( lep.get()->props["Charge"].f * sumcharge < 0 ) {
	    lep.get()->props["Rank3L"].i = 0;
	    eta0 = ( lep.get()->props["Flavour"].i == 13 ) ?  lep.get()->props["Eta"].f :  lep.get()->props["EtaBE2"].f;
	    lep0TLV.SetPtEtaPhiM(lep.get()->props["Pt"].f,eta0,lep.get()->props["Phi"].f,0.0);
	    break;
	}
    }

    // 2.
    // Find the closest lepton to lep0

    float min_dRLep0(999.0), this_dRLep0(999.0);
    int idx(-1), idxLep1(-1);
    for ( auto& lep : leptons ) {
	++idx;
	if ( lep.get()->props["Rank3L"].i == 0 ) { continue; }

	eta = ( lep.get()->props["Flavour"].i == 13 ) ?  lep.get()->props["Eta"].f : lep.get()->props["EtaBE2"].f;
	lepTLV.SetPtEtaPhiM(lep.get()->props["Pt"].f,eta,lep.get()->props["Phi"].f,0.0);

	this_dRLep0 = lepTLV.DeltaR(lep0TLV);

	if ( this_dRLep0 < min_dRLep0 ) {
	    min_dRLep0 = this_dRLep0;
	    idxLep1 = idx;
	}
    }

    leptons.at(idxLep1).get()->props["Rank3L"].i = 1;

    // 3.
    // Flag the other lepton

    for ( auto& lep : leptons ) {
	if ( lep.get()->props["Rank3L"].i == 0 ) { continue; }
	if ( lep.get()->props["Rank3L"].i == 1 ) { continue; }
	lep.get()->props["Rank3L"].i = 2;
    }

    if ( m_debug ) {
	Info("Sorter3L()","3L event lepton ranking:" );
	for ( auto lep : leptons ) {
	    if ( lep.get()->props["Rank3L"].i == 0 ) { std::cout << "Lepton 0 - pT = " << lep.get()->props["Pt"].f/1e3 << " (NTuple pT = " << m_lep0_INPUT_branches["Pt"].f/1e3 << " )" << std::endl; }
	    if ( lep.get()->props["Rank3L"].i == 1 ) { std::cout << "Lepton 1 - pT = " << lep.get()->props["Pt"].f/1e3 << " (NTuple pT = " << m_lep1_INPUT_branches["Pt"].f/1e3 << " )" << std::endl; }
	    if ( lep.get()->props["Rank3L"].i == 2 ) { std::cout << "Lepton 2 - pT = " << lep.get()->props["Pt"].f/1e3 << " (NTuple pT = " << m_lep2_INPUT_branches["Pt"].f/1e3 << " )" << std::endl; }
	}
    }

    return EL::StatusCode::SUCCESS;

}

EL::StatusCode HTopMultilepMiniNTupMaker :: triggerMatching()
{

    if ( m_debug ) { Info("triggerMatching()","Performing trigger matching..." ); }

    auto lep0 = m_leptons.at(0);
    auto lep1 = m_leptons.at(1);

    // Initialise w/ dummy values

    m_lep0_OUTPUT_branches["isTrigMatch"].c = m_lep1_OUTPUT_branches["isTrigMatch"].c = m_lep0_OUTPUT_branches["isTrigMatchDLT"].c = m_lep1_OUTPUT_branches["isTrigMatchDLT"].c = -1;

    // Get the indexes of leading/subleading leptons which passed the OLR

    // NB: For the SLT case, we always consider the OR of the SLT chains. Hence the pT threshold cut should be just the one for the lowest pT chain

    int el0_idx(-1), el1_idx(-1), mu0_idx(-1), mu1_idx(-1);

    if ( m_dilep_type == 1 ) { // 1) mumu

	ANA_CHECK( this->getPostOLRIndex( mu0_idx, 0, "muon" ) );
	ANA_CHECK( this->getPostOLRIndex( mu1_idx, 1, "muon" ) );

	if ( m_RunYear == 2015 ) {

	    m_lep0_OUTPUT_branches["isTrigMatch"].c = ( lep0.get()->props["Pt"].f > 1.05*20e3 && ( m_muon_match_HLT_mu20_iloose_L1MU15->at(mu0_idx) || m_muon_match_HLT_mu50->at(mu0_idx) ) );
	    m_lep1_OUTPUT_branches["isTrigMatch"].c = ( lep1.get()->props["Pt"].f > 1.05*20e3 && ( m_muon_match_HLT_mu20_iloose_L1MU15->at(mu1_idx) || m_muon_match_HLT_mu50->at(mu1_idx) ) );
	    m_lep0_OUTPUT_branches["isTrigMatchDLT"].c = ( m_muon_match_HLT_mu18_mu8noL1->at(mu0_idx) && lep0.get()->props["Pt"].f > 1.05*18e3 );
	    m_lep1_OUTPUT_branches["isTrigMatchDLT"].c = ( m_muon_match_HLT_mu18_mu8noL1->at(mu1_idx) && lep1.get()->props["Pt"].f > 1.05*8e3 );

	} else if ( m_RunYear == 2016 ) {


	    m_lep0_OUTPUT_branches["isTrigMatch"].c = ( lep0.get()->props["Pt"].f > 1.05*26e3 && ( m_muon_match_HLT_mu26_ivarmedium->at(mu0_idx) || m_muon_match_HLT_mu50->at(mu0_idx) ) );
	    m_lep1_OUTPUT_branches["isTrigMatch"].c = ( lep1.get()->props["Pt"].f > 1.05*26e3 && ( m_muon_match_HLT_mu26_ivarmedium->at(mu1_idx) || m_muon_match_HLT_mu50->at(mu1_idx) ) );
	    m_lep0_OUTPUT_branches["isTrigMatchDLT"].c = ( m_muon_match_HLT_mu22_mu8noL1->at(mu0_idx) && lep0.get()->props["Pt"].f > 1.05*22e3 );
	    m_lep1_OUTPUT_branches["isTrigMatchDLT"].c = ( m_muon_match_HLT_mu22_mu8noL1->at(mu1_idx) && lep1.get()->props["Pt"].f > 1.05*8e3 );

	}

    } else if ( m_dilep_type == 2 ) { // 2) OF

	ANA_CHECK( this->getPostOLRIndex( el0_idx, 0, "electron" ) );
	ANA_CHECK( this->getPostOLRIndex( mu0_idx, 0, "muon" ) );

	/*
	  std::cout << "\n2015: \n" << std::endl;
	  std::cout << "m_electron_match_HLT_e24_lhmedium_L1EM20VH->size() = " << m_electron_match_HLT_e24_lhmedium_L1EM20VH->size() << std::endl;
	  std::cout << "m_electron_match_HLT_e60_lhmedium->size() = " << m_electron_match_HLT_e60_lhmedium->size() << std::endl;
	  std::cout << "m_electron_match_HLT_e120_lhloose->size() = " << m_electron_match_HLT_e120_lhloose->size() << std::endl;
	  std::cout << "m_muon_match_HLT_mu20_iloose_L1MU15->size() = " << m_muon_match_HLT_mu20_iloose_L1MU15->size() << std::endl;
	  std::cout << "m_muon_match_HLT_mu50->size() = " << m_muon_match_HLT_mu50->size() << std::endl;
	  std::cout << "m_electron_match_HLT_e24_medium_L1EM20VHI_mu8noL1->size() = " << m_electron_match_HLT_e24_medium_L1EM20VHI_mu8noL1->size() << std::endl;
	  std::cout << "m_electron_match_HLT_e7_medium_mu24->size() = " << m_electron_match_HLT_e7_medium_mu24->size() << std::endl;
	  std::cout << "m_muon_match_HLT_e24_medium_L1EM20VHI_mu8noL1->size() = " << m_muon_match_HLT_e24_medium_L1EM20VHI_mu8noL1->size() << std::endl;
	  std::cout << "m_muon_match_HLT_e7_medium_mu24->size() = " << m_muon_match_HLT_e7_medium_mu24->size() << std::endl;
	  std::cout << "\n2016: \n" << std::endl;
	  std::cout << "m_electron_match_HLT_e26_lhtight_nod0_ivarloose->size() = " << m_electron_match_HLT_e26_lhtight_nod0_ivarloose->size() << std::endl;
	  std::cout << "m_electron_match_HLT_e60_lhmedium_nod0->size() = " << m_electron_match_HLT_e60_lhmedium_nod0->size() << std::endl;
	  std::cout << "m_electron_match_HLT_e140_lhloose_nod0->size() = " << m_electron_match_HLT_e140_lhloose_nod0->size() << std::endl;
	  std::cout << "m_muon_match_HLT_mu26_ivarmedium->size() = " << m_muon_match_HLT_mu26_ivarmedium->size() << std::endl;
	  std::cout << "m_muon_match_HLT_mu50->size() = " << m_muon_match_HLT_mu50->size() << std::endl;
	  std::cout << "m_electron_match_HLT_e17_lhloose_mu14->size() = " << m_electron_match_HLT_e17_lhloose_mu14->size() << std::endl;
	  std::cout << "m_electron_match_HLT_e17_lhloose_nod0_mu14->size() = " << m_electron_match_HLT_e17_lhloose_nod0_mu14->size() << std::endl;
	  std::cout << "m_muon_match_HLT_e17_lhloose_mu14->size() = " << m_muon_match_HLT_e17_lhloose_mu14->size() << std::endl;
	  std::cout << "m_muon_match_HLT_e17_lhloose_nod0_mu14->size() = " << m_muon_match_HLT_e17_lhloose_nod0_mu14->size() << std::endl;
	*/

	// NB: make sure to always read the leading component (passing OLR) of the trigmatch bits vector in this case!

	if ( m_RunYear == 2015 ) {

	    m_lep0_OUTPUT_branches["isTrigMatch"].c = ( ( lep0.get()->props["Flavour"].i == 11 && lep0.get()->props["Pt"].f > 25e3 && ( m_electron_match_HLT_e24_lhmedium_L1EM20VH->at(el0_idx) || m_electron_match_HLT_e60_lhmedium->at(el0_idx) || m_electron_match_HLT_e120_lhloose->at(el0_idx) ) ) || ( lep0.get()->props["Flavour"].i == 13 && lep0.get()->props["Pt"].f > 1.05*20e3 && ( m_muon_match_HLT_mu20_iloose_L1MU15->at(mu0_idx) || m_muon_match_HLT_mu50->at(mu0_idx) ) ) );
	    m_lep1_OUTPUT_branches["isTrigMatch"].c = ( ( lep1.get()->props["Flavour"].i == 11 && lep1.get()->props["Pt"].f > 25e3 && ( m_electron_match_HLT_e24_lhmedium_L1EM20VH->at(el0_idx) || m_electron_match_HLT_e60_lhmedium->at(el0_idx) || m_electron_match_HLT_e120_lhloose->at(el0_idx) ) ) || ( lep1.get()->props["Flavour"].i == 13 && lep1.get()->props["Pt"].f > 1.05*20e3 && ( m_muon_match_HLT_mu20_iloose_L1MU15->at(mu0_idx) || m_muon_match_HLT_mu50->at(mu0_idx) ) ) );
	    m_lep0_OUTPUT_branches["isTrigMatchDLT"].c = ( ( lep0.get()->props["Flavour"].i == 11 && ( ( m_electron_match_HLT_e7_medium_mu24->at(el0_idx) && lep0.get()->props["Pt"].f > 8e3 ) || ( m_electron_match_HLT_e24_medium_L1EM20VHI_mu8noL1->at(el0_idx) && lep0.get()->props["Pt"].f > 25e3 ) ) ) || ( lep0.get()->props["Flavour"].i == 13 && ( ( m_muon_match_HLT_e7_medium_mu24->at(mu0_idx) && lep0.get()->props["Pt"].f > 1.05*24e3 ) || ( m_muon_match_HLT_e24_medium_L1EM20VHI_mu8noL1->at(mu0_idx) && lep0.get()->props["Pt"].f > 1.05*8e3 ) ) ) );
	    m_lep1_OUTPUT_branches["isTrigMatchDLT"].c = ( ( lep1.get()->props["Flavour"].i == 11 && ( ( m_electron_match_HLT_e7_medium_mu24->at(el0_idx) && lep1.get()->props["Pt"].f > 8e3 ) || ( m_electron_match_HLT_e24_medium_L1EM20VHI_mu8noL1->at(el0_idx) && lep1.get()->props["Pt"].f > 25e3 ) ) ) || ( lep1.get()->props["Flavour"].i == 13 && ( ( m_muon_match_HLT_e7_medium_mu24->at(mu0_idx) && lep1.get()->props["Pt"].f > 1.05*24e3 ) || ( m_muon_match_HLT_e24_medium_L1EM20VHI_mu8noL1->at(mu0_idx) && lep1.get()->props["Pt"].f > 1.05*8e3 ) ) ) );

	} else if ( m_RunYear == 2016 ) {

	    m_lep0_OUTPUT_branches["isTrigMatch"].c = ( ( lep0.get()->props["Flavour"].i == 11 && lep0.get()->props["Pt"].f > 27e3 && ( m_electron_match_HLT_e26_lhtight_nod0_ivarloose->at(el0_idx) || m_electron_match_HLT_e60_lhmedium_nod0->at(el0_idx) || m_electron_match_HLT_e140_lhloose_nod0->at(el0_idx) ) ) || ( lep0.get()->props["Flavour"].i == 13 && lep0.get()->props["Pt"].f > 1.05*26e3 && ( m_muon_match_HLT_mu26_ivarmedium->at(mu0_idx) || m_muon_match_HLT_mu50->at(mu0_idx) ) ) );
	    m_lep1_OUTPUT_branches["isTrigMatch"].c = ( ( lep1.get()->props["Flavour"].i == 11 && lep1.get()->props["Pt"].f > 27e3 && ( m_electron_match_HLT_e26_lhtight_nod0_ivarloose->at(el0_idx) || m_electron_match_HLT_e60_lhmedium_nod0->at(el0_idx) || m_electron_match_HLT_e140_lhloose_nod0->at(el0_idx) ) ) || ( lep1.get()->props["Flavour"].i == 13 && lep1.get()->props["Pt"].f > 1.05*26e3 && ( m_muon_match_HLT_mu26_ivarmedium->at(mu0_idx) || m_muon_match_HLT_mu50->at(mu0_idx) ) ) );
	    m_lep0_OUTPUT_branches["isTrigMatchDLT"].c = ( ( lep0.get()->props["Flavour"].i == 11 && ( ( m_electron_match_HLT_e17_lhloose_mu14->at(el0_idx) && lep0.get()->props["Pt"].f > 18e3 ) || ( m_electron_match_HLT_e17_lhloose_nod0_mu14->at(el0_idx) && lep0.get()->props["Pt"].f > 18e3 ) ) ) || ( lep0.get()->props["Flavour"].i == 13 && ( ( m_muon_match_HLT_e17_lhloose_mu14->at(mu0_idx) && lep0.get()->props["Pt"].f > 1.05*14e3 ) || ( m_muon_match_HLT_e17_lhloose_nod0_mu14->at(mu0_idx) && lep0.get()->props["Pt"].f > 1.05*14e3 ) ) ) );
	    m_lep1_OUTPUT_branches["isTrigMatchDLT"].c = ( ( lep1.get()->props["Flavour"].i == 11 && ( ( m_electron_match_HLT_e17_lhloose_mu14->at(el0_idx) && lep1.get()->props["Pt"].f > 18e3 ) || ( m_electron_match_HLT_e17_lhloose_nod0_mu14->at(el0_idx) && lep1.get()->props["Pt"].f > 18e3 ) ) ) || ( lep1.get()->props["Flavour"].i == 13 && ( ( m_muon_match_HLT_e17_lhloose_mu14->at(mu0_idx) && lep1.get()->props["Pt"].f > 1.05*14e3 ) || ( m_muon_match_HLT_e17_lhloose_nod0_mu14->at(mu0_idx) && lep1.get()->props["Pt"].f > 1.05*14e3 ) ) ) );

	}

    } else if ( m_dilep_type == 3 ) { // 3) ee

	ANA_CHECK( this->getPostOLRIndex( el0_idx, 0, "electron" ) );
	ANA_CHECK( this->getPostOLRIndex( el1_idx, 1, "electron" ) );

	if ( m_RunYear == 2015 ) {

	    m_lep0_OUTPUT_branches["isTrigMatch"].c = ( lep0.get()->props["Pt"].f > 25e3 && ( m_electron_match_HLT_e24_lhmedium_L1EM20VH->at(el0_idx) || m_electron_match_HLT_e60_lhmedium->at(el0_idx) || m_electron_match_HLT_e120_lhloose->at(el0_idx) ) );
	    m_lep1_OUTPUT_branches["isTrigMatch"].c = ( lep1.get()->props["Pt"].f > 25e3 && ( m_electron_match_HLT_e24_lhmedium_L1EM20VH->at(el1_idx) || m_electron_match_HLT_e60_lhmedium->at(el1_idx) || m_electron_match_HLT_e120_lhloose->at(el1_idx) ) );
	    m_lep0_OUTPUT_branches["isTrigMatchDLT"].c = ( m_electron_match_HLT_2e12_lhloose_L12EM10VH->at(el0_idx) && lep0.get()->props["Pt"].f > 13e3 );
	    m_lep1_OUTPUT_branches["isTrigMatchDLT"].c = ( m_electron_match_HLT_2e12_lhloose_L12EM10VH->at(el1_idx) && lep1.get()->props["Pt"].f > 13e3 );

	} else if ( m_RunYear == 2016 ) {

	    m_lep0_OUTPUT_branches["isTrigMatch"].c = ( lep0.get()->props["Pt"].f > 27e3 && ( m_electron_match_HLT_e26_lhtight_nod0_ivarloose->at(el0_idx) || m_electron_match_HLT_e60_lhmedium_nod0->at(el0_idx) || m_electron_match_HLT_e140_lhloose_nod0->at(el0_idx) ) );
	    m_lep1_OUTPUT_branches["isTrigMatch"].c = ( lep1.get()->props["Pt"].f > 27e3 && ( m_electron_match_HLT_e26_lhtight_nod0_ivarloose->at(el1_idx) || m_electron_match_HLT_e60_lhmedium_nod0->at(el1_idx) || m_electron_match_HLT_e140_lhloose_nod0->at(el1_idx) ) );
	    m_lep0_OUTPUT_branches["isTrigMatchDLT"].c = ( m_electron_match_HLT_2e17_lhvloose_nod0->at(el0_idx) && lep0.get()->props["Pt"].f > 18e3 );
	    m_lep1_OUTPUT_branches["isTrigMatchDLT"].c = ( m_electron_match_HLT_2e17_lhvloose_nod0->at(el1_idx) && lep1.get()->props["Pt"].f > 18e3 );

	}

    }

    m_event_isTrigMatch_DLT = ( m_lep0_OUTPUT_branches["isTrigMatchDLT"].c && m_lep1_OUTPUT_branches["isTrigMatchDLT"].c );

    lep0.get()->props["isTrigMatch"].c = m_lep0_OUTPUT_branches["isTrigMatch"].c;
    lep1.get()->props["isTrigMatch"].c = m_lep1_OUTPUT_branches["isTrigMatch"].c;

    lep0.get()->props["isTrigMatchDLT"].c = m_lep0_OUTPUT_branches["isTrigMatchDLT"].c;
    lep1.get()->props["isTrigMatchDLT"].c = m_lep1_OUTPUT_branches["isTrigMatchDLT"].c;

    return EL::StatusCode::SUCCESS;

}


EL::StatusCode HTopMultilepMiniNTupMaker :: decorateWeights ()
{

    // Compute the event weights

    m_event.get()->weight_event = m_mcWeightOrg * m_pileupEventWeight_090 * m_MV2c10_70_EventWeight * m_JVT_EventWeight;

    float weight_lep(1.0);
    for ( auto lep : m_leptons ) {
	if ( lep.get()->props["isTightSelected"].c ) weight_lep *= lep.get()->props["SFObjTight"].f;
	else                                         weight_lep *= lep.get()->props["SFObjLoose"].f;
    }
    m_event.get()->weight_event_lep = weight_lep;

    // Set the weights (trigger and lepton reco, iso...) for tag and probe leptons.
    // Do it for SLT-only selection.

    if ( m_event.get()->dilep_type && !m_event.get()->isBadTPEvent_SLT ) {

	auto lep0 = m_leptons.at(0);
	auto lep1 = m_leptons.at(1);

	// We assume events at this stage will have exactly 1 lepton flagged TAG and 1 flagged probe
	// (NB: in the code this is the case also for those ambiguous events where both leptons are effectively counted as tag candidates in the vector branches. The "tag" flag is chosen at random...).
	//
	// Since we always will require the tag to be tight and trigger-matched, the tag weight will account for the SFTrigTight and SFObjTight contributions.
	//
	// As for the probe, in the baseline case the user will need to apply *only* the "weight_lep_probe" to correct for the reco,ID... efficiency.
	// In fact, we don't mind about the trigger for the probe lepton.
	// This weight will be SFObjTight or SFObjLoose depending on which selection the probe passes.
	// We save also the "weight_trig_probe" weight for the special case when one wants to measure (SLT) trigger-dependent efficiencies. In that case, also the probe
	// will be explicitly requested to fire the trigger. If it does, the weight will be "SFTrigTight" or "SFTrigLoose", depending on the offline selection on the probe.

	if ( lep0.get()->props["isTagSLT"].c && !lep1.get()->props["isTagSLT"].c ) {

	    // Trigger efficiency

	    m_event.get()->weight_trig_tag   = lep0.get()->props["SFTrigTight"].f;
	    if ( lep1.get()->props["isTrigMatch"].c ) {
		m_event.get()->weight_trig_probe = ( lep1.get()->props["isTightSelected"].c ) ? lep1.get()->props["SFTrigTight"].f : lep1.get()->props["SFTrigLoose"].f;
	    }

	    // Lepton reco, iso, ID... efficiency

	    m_event.get()->weight_lep_tag    = lep0.get()->props["SFObjTight"].f;
	    m_event.get()->weight_lep_probe  = ( lep1.get()->props["isTightSelected"].c ) ? lep1.get()->props["SFObjTight"].f : lep1.get()->props["SFObjLoose"].f;

	} else if ( !lep0.get()->props["isTagSLT"].c && lep1.get()->props["isTagSLT"].c ) {

	    // Trigger efficiency

	    m_event.get()->weight_trig_tag   = lep1.get()->props["SFTrigTight"].f;
	    if ( lep0.get()->props["isTrigMatch"].c ) {
		m_event.get()->weight_trig_probe = ( lep0.get()->props["isTightSelected"].c ) ? lep0.get()->props["SFTrigTight"].f : lep0.get()->props["SFTrigLoose"].f;
	    }

	    // Lepton reco, iso, ID... efficiency

	    m_event.get()->weight_lep_tag    = lep1.get()->props["SFObjTight"].f;
	    m_event.get()->weight_lep_probe  = ( lep0.get()->props["isTightSelected"].c ) ? lep0.get()->props["SFObjTight"].f : lep0.get()->props["SFObjLoose"].f;

	} else if ( lep0.get()->props["isTagSLT"].c && lep1.get()->props["isTagSLT"].c ) {

	    Error("decorateWeights()", "Entry %u - EventNumber = %u - RunYear = %i - has TWO leptons flagged TAG. This shouldn't happen. Aborting...", static_cast<uint32_t>(m_numEntry), static_cast<uint32_t>(m_EventNumber), m_RunYear );
	    return EL::StatusCode::FAILURE;

	} else if ( !lep0.get()->props["isTagSLT"].c && !lep1.get()->props["isTagSLT"].c ) {

	    Error("decorateWeights()", "Entry %u - EventNumber = %u - RunYear = %i - has ZERO leptons flagged TAG. This shouldn't happen. Aborting...", static_cast<uint32_t>(m_numEntry), static_cast<uint32_t>(m_EventNumber), m_RunYear );
	    return EL::StatusCode::FAILURE;

	}

    }

    // Calculate *single* lepton trigger SF for the event.
    //
    // The following is based on the assumption that one OR the other lepton will fire the SLT:
    //
    // ( lep_isTrigMatch_0 || lep_isTrigMatch_1 )
    //
    // The CP tools do not calculate the final event SF. Rather, they give back the SF and MC efficiency *per lepton*.
    // With such ingredients in hand, the HTop way to compute the event SF (for SLT lep0 OR lep1 selection) is the following:
    //
    // eventSF = ( 1 - prod( 1 - SF(i)*eff(i) ) ) / ( 1 - prod ( 1 - eff(i) ) );
    //
    // where the productory is over the selected leptons in the event.
    // The trick at the numerator is just to get the efficiency in data.
    //
    // The SF systematics are obtained by coherently varying the SF for each object (i.e. assume full correlation).
    // The MC efficiency is assumed to have negligible uncertainty.

    // This should roughly correspond to the weights "lepSFTrigTight", "lepSFTrigLoose" in the Group NTup.

    float trig_weight_N(1.0), trig_weight_D(1.0);
    float this_SF(1.0), this_eff(0.0);
    for ( auto lep : m_leptons ) {

	this_eff = ( lep.get()->props["isTightSelected"].c ) ? lep.get()->props["EffTrigTight"].f : lep.get()->props["EffTrigLoose"].f;
	this_SF  = ( lep.get()->props["isTightSelected"].c ) ? lep.get()->props["SFTrigTight"].f  : lep.get()->props["SFTrigLoose"].f;

	trig_weight_N *= ( 1.0 - this_SF * this_eff );
	trig_weight_D *= ( 1.0 - this_eff );

	if ( m_debug ) {
	    Info("decorateWeights()", "\tthis_eff = %.2f - this_SF = %.2f", this_eff, this_SF );
	    Info("decorateWeights()", "\tN block = %.2f - D block = %.2f", trig_weight_N, trig_weight_D );
	}

    }

    // Update numerator and denominator
    // Make sure the SF in the 0/0 case (i.e, when efficiency=0) will be set equal to 1

    trig_weight_N = ( trig_weight_N != 1.0 ) ? trig_weight_N : 0.0;
    trig_weight_D = ( trig_weight_D != 1.0 ) ? trig_weight_D : 0.0;

    m_event.get()->weight_event_trig_SLT = ( 1.0 - trig_weight_N ) / ( 1.0 - trig_weight_D );

    if ( m_debug ) {
	Info("decorateWeights()", "N = %.2f - D = %.2f", trig_weight_N, trig_weight_D );
	Info("decorateWeights()", "per-event trigger weight = %.2f", m_event.get()->weight_event_trig_SLT );
    }

    return EL::StatusCode::SUCCESS;
}


EL::StatusCode HTopMultilepMiniNTupMaker :: defineTagAndProbe ()
{

    // Clear vector branches from previous event

    ANA_CHECK( this->clearBranches("tag_and_probe") );

    // Do this only for dilepton events

    if ( !m_event.get()->dilep_type ) { return EL::StatusCode::SUCCESS; }

    // Minimal requirement of pT = 10 GeV on all leptons

    for ( auto lep : m_leptons ) {  if ( lep.get()->props["Pt"].f < 10e3 ) { return EL::StatusCode::SUCCESS; } }

    if ( m_useTruthTP ) {

	// NB: the rationale assumes this will be done **only** on ttbar nonallhad!

	if ( m_event.get()->isSS01 ) {

	    // SS events:
	    //
	    // -) tag: the prompt lepton (still, must not be qmisid!)
	    // -) probe: the other lepton.
	    //
	    // If probe is (QMisID or prompt), or tag hasn't been found, flag the event as bad.
	    // This enforces the probe to be always ( non-prompt (HF lep) ||  misID-jet || photon conversion...) ==> a fake!
	    //
	    // In fact, we assume no SS prompt-prompt events ever exist in ttbar nonallhad, unless one of the two has flipped the charge due to mis-reco.

	    bool found_tag(false);
	    int tag_idx(0);
	    for ( auto lep : m_leptons ) {
		if ( ( lep.get()->props["isPrompt"].c == 1 || ( lep.get()->props["isBrems"].c == 1 && lep.get()->props["isQMisID"].c == 0 ) ) &&
		     ( lep.get()->props["isQMisID"].c == 0 ) ) {
		    lep.get()->props["isTagSLT"].c = 1;
		    found_tag = true;
		    break;
		}
		++tag_idx;
	    }
	    int probe_idx = ( tag_idx ) ? 0 : 1; // Our lepton vector has only 2 components ;-)

	    m_event.get()->isBadTPEvent_SLT = ( !found_tag ||
						m_leptons.at(probe_idx).get()->props["isQMisID"].c == 1 ||
						m_leptons.at(probe_idx).get()->props["isPrompt"].c == 1 ||
						( m_leptons.at(probe_idx).get()->props["isBrems"].c == 1 && m_leptons.at(probe_idx).get()->props["isQMisID"].c == 0 ) );

	} else {

	    // OS events:
	    //
	    // If event contains a !prompt or a qmisid, flag the event as bad and return.
	    // Else:
	    // -) tag: choose randomly
	    // -) probe: the other lepton.

	    m_event.get()->isBadTPEvent_SLT = 0;  // be optimistic!

	    for ( auto lep : m_leptons ) {
		if ( lep.get()->props["isPrompt"].c == 0 || lep.get()->props["isQMisID"].c == 1 ) {
		    m_event.get()->isBadTPEvent_SLT = 1;
		    return EL::StatusCode::SUCCESS;
		}
	    }

	    int tag_idx = ( m_rand->Rndm() > 0.5 ); // will pick index 0 or 1 in lepton vector randomly

	    m_leptons.at(tag_idx).get()->props["isTagSLT"].c = 1;

	}

	return EL::StatusCode::SUCCESS;

    }

    if ( m_useNominalTP ) {

	// Standard tag-and-probe method (a-la-SUSY SS analysis):
	//
	// In the Real OS CR, ambiguous events where both leptons are T & T.M. are double counted -->
	// both leptons are alternatively considered as the possible tag/probe
	// to avoid any bias in the choice of the tag and to increase the statistics
	//
	// We assume SLT only are being used in the CRs where the efficiencies will be measured.

	const std::string lepSelection = ( m_lepSelForTP.compare("CutBased") == 0 ) ? "isTightSelected" : "isTightSelectedMVA";

	m_event.get()->isBadTPEvent_SLT = 0; // be optimistic!

	int tag_candidate_counter_SLT(0);
	int idx(0);
	int tag_idx_SLT(-1), probe_idx_SLT(-1);

	// Count the number of tag lepton candidates

	if ( m_debug ) { std::cout << "Looking for tag candidates:\n" << std::endl; }

	for ( auto lep : m_leptons ) {

	    if ( m_debug ) { Info("defineTagAndProbe()","Checking lepton[%i] w/ pT = %.2f, flavour = %i", idx, lep.get()->props["Pt"].f/1e3, lep.get()->props["Flavour"].i ); }

	    if ( lep.get()->props[lepSelection].c && lep.get()->props["isTrigMatch"].c ) {
		++tag_candidate_counter_SLT;
		tag_idx_SLT = idx;
		if ( m_debug ) { Info("defineTagAndProbe()","\t ===> found a tag candidate (SLT matching)! pT[%i] = %.2f, flavour = %i", idx, lep.get()->props["Pt"].f/1e3, lep.get()->props["Flavour"].i ); }
	    }

	    ++idx;
	}

	if ( m_debug ) { std::cout << "" << std::endl; }

	// Different treatment for OS and SS events

	std::string key, tp_key, type;

	int this_idx(-1);

	if ( !m_event.get()->isSS01 ) {

	    // SLT

	    switch ( tag_candidate_counter_SLT ) {

	    case 0:

		m_event.get()->isBadTPEvent_SLT = 1;

		if ( m_debug ) { Info("defineTagAndProbe()","No lepton (T & TM) (SLT matching) was found - flag this event as bad" ); }

		break;

	    case 1: // this is the non-ambiguous case: 1 tag and 1 probe per event

		probe_idx_SLT = abs( tag_idx_SLT - 1 ); // NB: this works since we only look at dilepton events (tag_idx can be 0 or 1 )

		for ( const auto& tp : m_TPS ) {

		    if      ( tp.compare("Tag") == 0 )   { this_idx = tag_idx_SLT; }
		    else if ( tp.compare("Probe") == 0 ) { this_idx = probe_idx_SLT; }

		    for ( const auto& trig : m_TRIGS ) {

			for ( const auto& var : m_TP_VARS ) {

			    if ( var.find("VEC") == std::string::npos ) { continue; }

			    key  = var.substr( 0, var.length() - 5 );
			    type = var.substr( var.length() - 4 );

			    tp_key = "lep_" + tp + "Vec" + "_" + trig + "_" + key;

			    // Fill T&P vector branches for *leptons*

			    if ( type.compare("VECF") == 0 ) { m_TagProbe_branches[tp_key].vec_f.push_back( m_leptons.at(this_idx).get()->props[key].f ); }
			    if ( type.compare("VECB") == 0 ) { m_TagProbe_branches[tp_key].vec_c.push_back( m_leptons.at(this_idx).get()->props[key].c ); }
			    if ( type.compare("VECI") == 0 ) { m_TagProbe_branches[tp_key].vec_i.push_back( m_leptons.at(this_idx).get()->props[key].i ); }

			    if      ( m_leptons.at(this_idx).get()->props["Flavour"].i == 11 ) { tp_key = "electron_" + tp + "Vec" + "_" + trig + "_" + key; }
			    else if ( m_leptons.at(this_idx).get()->props["Flavour"].i == 13 ) { tp_key = "muon_" + tp + "Vec" + "_" + trig + "_" + key; }

			    // Fill T&P vector branches for *electrons* or *muons*

			    if ( type.compare("VECF") == 0 ) { m_TagProbe_branches[tp_key].vec_f.push_back( m_leptons.at(this_idx).get()->props[key].f ); }
			    if ( type.compare("VECB") == 0 ) { m_TagProbe_branches[tp_key].vec_c.push_back( m_leptons.at(this_idx).get()->props[key].c ); }
			    if ( type.compare("VECI") == 0 ) { m_TagProbe_branches[tp_key].vec_i.push_back( m_leptons.at(this_idx).get()->props[key].i ); }

			}

		    }

		}

		// This will set the flat branch as well for tag and probe later on (needed to define cuts in plotting tools)

		m_leptons.at(tag_idx_SLT).get()->props["isTagSLT"].c = 1;

		if ( m_debug ) { Info("defineTagAndProbe()","Unambiguous event (SLT matching) : tag lepton pT = %.2f, probe lepton pT = %.2f", m_leptons.at(tag_idx_SLT).get()->props["Pt"].f/1e3, m_leptons.at(probe_idx_SLT).get()->props["Pt"].f/1e3 ); }

		break;

	    case 2: // this is the ambiguous case: both T & TM

		// In OS, we consider both leptons as tag and probe (event will be double-counted). The reason is b/c we require
		// both leptons be real in this region (in data we will subtract backgrounds), so here there's really no distinction possible

		for ( auto lep : m_leptons ) {

		    for ( const auto& tp : m_TPS ) {

			for ( const auto& trig : m_TRIGS ) {

			    for ( const auto& var : m_TP_VARS ) {

				if ( var.find("VEC") == std::string::npos ) { continue; }

				key  = var.substr( 0, var.length() - 5 );
				type = var.substr( var.length() - 4 );

				tp_key = "lep_" + tp + "Vec" + "_" + trig + "_" + key;

				// Fill T&P vector branches for *leptons*. The same lepton property goes both in the tag and probe vectors.

				if ( type.compare("VECF") == 0 ) { m_TagProbe_branches[tp_key].vec_f.push_back( lep.get()->props[key].f ); }
				if ( type.compare("VECB") == 0 ) { m_TagProbe_branches[tp_key].vec_c.push_back( lep.get()->props[key].c ); }
				if ( type.compare("VECI") == 0 ) { m_TagProbe_branches[tp_key].vec_i.push_back( lep.get()->props[key].i ); }

				if      ( lep.get()->props["Flavour"].i == 11 ) { tp_key = "electron_" + tp + "Vec" + "_" + trig + "_" + key; }
				else if ( lep.get()->props["Flavour"].i == 13 ) { tp_key = "muon_" + tp + "Vec" + "_" + trig + "_" + key; }

				// Fill T&P vector branches for *electrons* or *muons*

				if ( type.compare("VECF") == 0 ) { m_TagProbe_branches[tp_key].vec_f.push_back( lep.get()->props[key].f ); }
				if ( type.compare("VECB") == 0 ) { m_TagProbe_branches[tp_key].vec_c.push_back( lep.get()->props[key].c ); }
				if ( type.compare("VECI") == 0 ) { m_TagProbe_branches[tp_key].vec_i.push_back( lep.get()->props[key].i ); }

			    }

			}

		    }

		}

		// This will allow to set the flat branch as well for tag and probe later on
		//
		// For convenience, choose the tag and probe randomly
		//
		// This has to be done b/c at plotting level we ask for the probe to be T/L to define the N and D for efficiency,
		// but in these events both are T and T.M., so it doesn't really matter which lepton we picked as tag/probe
		// In fact, in this case, when plotting the vector branch "lep_ProbeVec_*", both leptons will be considered effectively as the probe

		tag_idx_SLT = ( m_rand->Rndm() > 0.5 );

		m_leptons.at( tag_idx_SLT ).get()->props["isTagSLT"].c = 1;

		if ( m_debug ) { Info("defineTagAndProbe()","Ambiguous event (SLT matching): both leptons T & T.M. Event is OS --> will consider both as tag and probe"); }

		break;

	    default:

		Error("defineTagAndProbe()","Number of tag lepton candidates (SLT matching) is %i...This shouldn't happen. Aborting.", tag_candidate_counter_SLT );
		return EL::StatusCode::FAILURE;

	    } // close switch SLT

	} else { // closes OS case

	    // In SS, even in the ambiguous case we still need to tell which is tag and which is probe,
	    // since we really do want to take the fake as the probe! (unlike OS, where we require both leptons be real, in SS we have 1 real and 1 fake)

	    // *****************************************************************
	    // *
	    // * This is what SUSY actually does for th electron fake efficiency
	    // *
	    // *

	    if ( m_ambiSolvingCrit.compare("OF") == 0 && m_event.get()->dilep_type == 2 ) {

		// Look only at emu,mue events.
		// Select only events where the muon is T (& TM), and flag it as the tag
		// The electron will be the probe, reagrdless of any selection on it.
		//
		// This works b/c if the muon is T (& TM) is less likely to be a fake. The electron probe will then be completely unbiased.

		int muon_tag_candidate_counter_SLT(0);

		if ( m_debug ) { Info("defineTagAndProbe()","Looking at OF event..." ); }

		// SLT

		int i_SLT(0);
		for ( auto lep : m_leptons ) {

		    if ( m_debug ) { Info("defineTagAndProbe()","Checking lepton[%i] w/ pT = %.2f, flavour = %i", i_SLT, lep.get()->props["Pt"].f/1e3, lep.get()->props["Flavour"].i ); }

		    if ( lep.get()->props["Flavour"].i == 13 && lep.get()->props[lepSelection].c && lep.get()->props["isTrigMatch"].c ) {

			tag_idx_SLT = i_SLT;

			if ( m_debug ) { Info("defineTagAndProbe()","\t ===> found a muon tag candidate (SLT matching)! pT[%i] = %.2f", i_SLT, lep.get()->props["Pt"].f/1e3 ); }

			++muon_tag_candidate_counter_SLT;

			break;

		    }

		    ++i_SLT;
		}

		if ( !muon_tag_candidate_counter_SLT ) {

		    m_event.get()->isBadTPEvent_SLT = 1;

		    if ( m_debug ) { Info("defineTagAndProbe()","No muon (T & TM) (SLT matching) was found - flag this event as bad" ); }

		} else {

		    probe_idx_SLT = ( tag_idx_SLT ) ? 0 : 1; // Our lepton vector has only 2 components ;-)

		    for ( const auto& tp : m_TPS ) {

			if      ( tp.compare("Tag") == 0 )   { this_idx = tag_idx_SLT; }
			else if ( tp.compare("Probe") == 0 ) { this_idx = probe_idx_SLT; }

			for ( const auto& trig : m_TRIGS ) {

			    for ( const auto& var : m_TP_VARS ) {

				if ( var.find("VEC") == std::string::npos ) { continue; }

				key  = var.substr( 0, var.length() - 5 );
				type = var.substr( var.length() - 4 );

				tp_key = "lep_" + tp + "Vec" + "_" + trig + "_" + key;

				// Fill T&P vector branches for *leptons*

				if ( type.compare("VECF") == 0 ) { m_TagProbe_branches[tp_key].vec_f.push_back( m_leptons.at(this_idx).get()->props[key].f ); }
				if ( type.compare("VECB") == 0 ) { m_TagProbe_branches[tp_key].vec_c.push_back( m_leptons.at(this_idx).get()->props[key].c ); }
				if ( type.compare("VECI") == 0 ) { m_TagProbe_branches[tp_key].vec_i.push_back( m_leptons.at(this_idx).get()->props[key].i ); }

				if      ( m_leptons.at(this_idx).get()->props["Flavour"].i == 11 ) { tp_key = "electron_" + tp + "Vec" + "_" + trig + "_" + key; }
				else if ( m_leptons.at(this_idx).get()->props["Flavour"].i == 13 ) { tp_key = "muon_" + tp + "Vec" + "_" + trig + "_" + key; }

				// Fill T&P vector branches for *electrons* or *muons*

				if ( type.compare("VECF") == 0 ) { m_TagProbe_branches[tp_key].vec_f.push_back( m_leptons.at(this_idx).get()->props[key].f ); }
				if ( type.compare("VECB") == 0 ) { m_TagProbe_branches[tp_key].vec_c.push_back( m_leptons.at(this_idx).get()->props[key].c ); }
				if ( type.compare("VECI") == 0 ) { m_TagProbe_branches[tp_key].vec_i.push_back( m_leptons.at(this_idx).get()->props[key].i ); }

			    }

			}

		    }

		    // This will set the flat branch as well for tag and probe later on (needed to define cuts in plotting tools)

		    m_leptons.at(tag_idx_SLT).get()->props["isTagSLT"].c = 1;

		    if ( m_debug ) { Info("defineTagAndProbe()","Good OF event (SLT matching) : tag lepton pT = %.2f, flavour = %i, probe lepton pT = %.2f, flavour = %i", m_leptons.at(tag_idx_SLT).get()->props["Pt"].f/1e3, m_leptons.at(tag_idx_SLT).get()->props["Flavour"].i, m_leptons.at(probe_idx_SLT).get()->props["Pt"].f/1e3, m_leptons.at(probe_idx_SLT).get()->props["Flavour"].i ); }

		}

		// ----------------------------------------

		// All done w/ this event!

		return EL::StatusCode::SUCCESS;

	    } // close "OF" case where event is actually OF

	    // *********************************************************************************
	    // *
	    // * These are attempts at solving the ambiguity w/ least bias possible on the probe
	    // *
	    // *

	    // SLT

	    switch ( tag_candidate_counter_SLT ) {

	    case 0:

		m_event.get()->isBadTPEvent_SLT = 1;

		if ( m_debug ) { Info("defineTagAndProbe()","No lepton (T & TM) (SLT matching) was found - flag this event as bad" ); }

		break;

	    case 1: // this is the non-ambiguous case: 1 tag and 1 probe per event

		probe_idx_SLT = abs( tag_idx_SLT - 1 ); // NB: this works since we only look at dilepton events (tag_idx can be 0 or 1 )

		for ( const auto& tp : m_TPS ) {

		    if      ( tp.compare("Tag") == 0 )   { this_idx = tag_idx_SLT; }
		    else if ( tp.compare("Probe") == 0 ) { this_idx = probe_idx_SLT; }

		    for ( const auto& trig : m_TRIGS ) {

			for ( const auto& var : m_TP_VARS ) {

			    if ( var.find("VEC") == std::string::npos ) { continue; }

			    key  = var.substr( 0, var.length() - 5 );
			    type = var.substr( var.length() - 4 );

			    tp_key = "lep_" + tp + "Vec" + "_" + trig + "_" + key;

			    // Fill T&P vector branches for *leptons*

			    if ( type.compare("VECF") == 0 ) { m_TagProbe_branches[tp_key].vec_f.push_back( m_leptons.at(this_idx).get()->props[key].f ); }
			    if ( type.compare("VECB") == 0 ) { m_TagProbe_branches[tp_key].vec_c.push_back( m_leptons.at(this_idx).get()->props[key].c ); }
			    if ( type.compare("VECI") == 0 ) { m_TagProbe_branches[tp_key].vec_i.push_back( m_leptons.at(this_idx).get()->props[key].i ); }

			    if      ( m_leptons.at(this_idx).get()->props["Flavour"].i == 11 ) { tp_key = "electron_" + tp + "Vec" + "_" + trig + "_" + key; }
			    else if ( m_leptons.at(this_idx).get()->props["Flavour"].i == 13 ) { tp_key = "muon_" + tp + "Vec" + "_" + trig + "_" + key; }

			    // Fill T&P vector branches for *electrons* or *muons*

			    if ( type.compare("VECF") == 0 ) { m_TagProbe_branches[tp_key].vec_f.push_back( m_leptons.at(this_idx).get()->props[key].f ); }
			    if ( type.compare("VECB") == 0 ) { m_TagProbe_branches[tp_key].vec_c.push_back( m_leptons.at(this_idx).get()->props[key].c ); }
			    if ( type.compare("VECI") == 0 ) { m_TagProbe_branches[tp_key].vec_i.push_back( m_leptons.at(this_idx).get()->props[key].i ); }

			}

		    }

		}

		// This will set the flat branch as well for tag and probe later on (needed to define cuts in plotting tools)

		m_leptons.at(tag_idx_SLT).get()->props["isTagSLT"].c = 1;

		if ( m_debug ) { Info("defineTagAndProbe()","Unambiguous event (SLT matching) : tag lepton pT = %.2f, probe lepton pT = %.2f", m_leptons.at(tag_idx_SLT).get()->props["Pt"].f/1e3, m_leptons.at(probe_idx_SLT).get()->props["Pt"].f/1e3 ); }

		break;

	    case 2: // this is the ambiguous case: both T & TM

		// In SS, even in the ambiguous case we still need to tell which is tag and which is probe,
		// since we really do want to take the fake as the probe! (unlike OS, where we require both leptons be real, in SS we have 1 real and 1 fake)

		if ( m_ambiSolvingCrit.compare("Pt") == 0 ||  m_ambiSolvingCrit.compare("OF") == 0 ) { // this will be used also when the solving criterion is "OF", but the event is SF

		    // Make sure lepton container is sorted in descending order of lepton pt

		    std::sort( m_leptons.begin(), m_leptons.end(), MiniNTupMaker::SorterPt() );
		    if ( m_debug ) {
			std:: cout << "\nLepton container ( descending pT sorting ):\n" << std::endl;
			for ( unsigned int idx(0); idx < m_leptons.size(); ++idx ) {
			    std::cout << "lepton[" << idx << "] - DeltaR(lep, closest bjet) = " << m_leptons.at(idx).get()->props["deltaRClosestBJet"].f << ", pT = " << m_leptons.at(idx).get()->props["Pt"].f/1e3 << std::endl;
			}
			std::cout << "" << std::endl;
		    }

		    // Take the leading as the tag, the subleading as the probe

		    tag_idx_SLT   = 0;
		    probe_idx_SLT = 1;

		} else if ( m_ambiSolvingCrit.compare("deltaRClosestBJet") == 0 ) {

		    // Try an event topology-based approach to select the tag and probe:
		    //
		    // -) take the lepton that is further apart from the closest bjet as the tag, the other will be the probe

		    // Require *at least* one bjet, otherwise flag the event as bad

		    if ( m_event.get()->nbjets < 1 ) {

			m_event.get()->isBadTPEvent_SLT = 1;

			break;
		    }

		    // Sort lepton container in descending order of distance to closest bjet

		    std::sort( m_leptons.begin(), m_leptons.end(), MiniNTupMaker::SorterDistanceClosestBJet() );
		    if ( m_debug ) {
			std:: cout << "\nLepton container ( descending DeltaR(lep, closest bjet) sorting ):\n" << std::endl;
			for ( unsigned int idx(0); idx < m_leptons.size(); ++idx ) {
			    std::cout << "lepton[" << idx << "] - DeltaR(lep, closest bjet) = " << m_leptons.at(idx).get()->props["deltaRClosestBJet"].f << ", pT = " << m_leptons.at(idx).get()->props["Pt"].f/1e3 << std::endl;
			}
			std::cout << "" << std::endl;
		    }

		    tag_idx_SLT   = 0;
		    probe_idx_SLT = 1;

		} else if ( m_ambiSolvingCrit.compare("massClosestBJet") == 0  ) {

		    // Try an event topology-based approach to select the tag and probe:
		    //
		    // -) take the lepton forming the largest invariant mass w/ the closest bjet as the tag, the other will be the probe

		    // Require *at least* one bjet, otherwise flag the event as bad

		    if ( m_event.get()->nbjets < 1 ) {

			m_event.get()->isBadTPEvent_SLT = 1;

			break;
		    }

		    // Sort lepton container in descending order of mass w/ closest bjet

		    std::sort( m_leptons.begin(), m_leptons.end(), MiniNTupMaker::SorterMassClosestBJet() );
		    if ( m_debug ) {
			std:: cout << "\nLepton container ( descending M(lep, closest bjet) sorting ):\n" << std::endl;
			for ( unsigned int idx(0); idx < m_leptons.size(); ++idx ) {
			    std::cout << "lepton[" << idx << "] - M(lep, closest bjet) = " << m_leptons.at(idx).get()->props["massClosestBJet"].f/1e3 << ", pT = " << m_leptons.at(idx).get()->props["Pt"].f/1e3 << std::endl;
			}
			std::cout << "" << std::endl;
		    }

		    tag_idx_SLT   = 0;
		    probe_idx_SLT = 1;

		}

		// Flag the tag lepton (used afterwards to define flat branches)

		m_leptons.at( tag_idx_SLT ).get()->props["isTagSLT"].c = 1;

		for ( const auto& tp : m_TPS ) {

		    if      ( tp.compare("Tag") == 0 )   { this_idx = tag_idx_SLT; }
		    else if ( tp.compare("Probe") == 0 ) { this_idx = probe_idx_SLT; }

		    for ( const auto& trig : m_TRIGS ) {

			for ( const auto& var : m_TP_VARS ) {

			    if ( var.find("VEC") == std::string::npos ) { continue; }

			    key  = var.substr( 0, var.length() - 5 );
			    type = var.substr( var.length() - 4 );

			    tp_key = "lep_" + tp + "Vec" + "_" + trig + "_" + key;

			    // Fill T&P vector branches for *leptons*

			    if ( type.compare("VECF") == 0 ) { m_TagProbe_branches[tp_key].vec_f.push_back( m_leptons.at(this_idx).get()->props[key].f ); }
			    if ( type.compare("VECB") == 0 ) { m_TagProbe_branches[tp_key].vec_c.push_back( m_leptons.at(this_idx).get()->props[key].c ); }
			    if ( type.compare("VECI") == 0 ) { m_TagProbe_branches[tp_key].vec_i.push_back( m_leptons.at(this_idx).get()->props[key].i ); }

			    if      ( m_leptons.at(this_idx).get()->props["Flavour"].i == 11 ) { tp_key = "electron_" + tp + "Vec" + "_" + trig + "_" + key; }
			    else if ( m_leptons.at(this_idx).get()->props["Flavour"].i == 13 ) { tp_key = "muon_" + tp + "Vec" + "_" + trig + "_" + key; }

			    // Fill T&P vector branches for *electrons* or *muons*

			    if ( type.compare("VECF") == 0 ) { m_TagProbe_branches[tp_key].vec_f.push_back( m_leptons.at(this_idx).get()->props[key].f ); }
			    if ( type.compare("VECB") == 0 ) { m_TagProbe_branches[tp_key].vec_c.push_back( m_leptons.at(this_idx).get()->props[key].c ); }
			    if ( type.compare("VECI") == 0 ) { m_TagProbe_branches[tp_key].vec_i.push_back( m_leptons.at(this_idx).get()->props[key].i ); }

			}

		    }

		}

		if ( m_debug ) {

		    float tag_pt   = m_leptons.at(tag_idx_SLT).get()->props["Pt"].f;
		    float probe_pt = m_leptons.at(probe_idx_SLT).get()->props["Pt"].f;
		    float tag_deltaRClosestBJet   = m_leptons.at(tag_idx_SLT).get()->props["deltaRClosestBJet"].f;
		    float probe_deltaRClosestBJet = m_leptons.at(probe_idx_SLT).get()->props["deltaRClosestBJet"].f;
		    float tag_massClosestBJet   = m_leptons.at(tag_idx_SLT).get()->props["massClosestBJet"].f;
		    float probe_massClosestBJet = m_leptons.at(probe_idx_SLT).get()->props["massClosestBJet"].f;

		    if (m_ambiSolvingCrit.compare("Pt") == 0  )               { Info("defineTagAndProbe()","Ambiguous event (SLT matching): both leptons T & T.M. Event is SS --> tag lepton pT = %.2f, probe lepton pT = %.2f", tag_pt/1e3, probe_pt/1e3 ); }
		    if (m_ambiSolvingCrit.compare("deltaRClosestBJet") == 0 ) { Info("defineTagAndProbe()","Ambiguous event (SLT matching): both leptons T & T.M. Event is SS --> tag lepton pT = %.2f, DeltaR(lep, closest bjet) = %.3f, probe lepton pT = %.2f, DeltaR(lep, closest bjet) = %.3f", tag_pt/1e3, tag_deltaRClosestBJet, probe_pt/1e3, probe_deltaRClosestBJet  ); }
		    if (m_ambiSolvingCrit.compare("massClosestBJet") == 0 )   { Info("defineTagAndProbe()","Ambiguous event (SLT matching): both leptons T & T.M. Event is SS --> tag lepton pT = %.2f, M(lep, closest bjet) = %.3f, probe lepton pT = %.2f, M(lep, closest bjet) = %.3f", tag_pt/1e3, tag_massClosestBJet/1e3, probe_pt/1e3, probe_massClosestBJet/1e3 ); }

		}

		// Revert back to pT-sorting (which is the default)

		std::sort( m_leptons.begin(), m_leptons.end(), MiniNTupMaker::SorterPt() );

		break;

	    default:

		Error("defineTagAndProbe()","Number of tag lepton candidates (SLT matching) is %i...This shouldn't happen. Aborting.", tag_candidate_counter_SLT );
		return EL::StatusCode::FAILURE;

	    } // close switch SLT

	} // close OS-SS cases

    } // close standard TP

    return EL::StatusCode::SUCCESS;

}

EL::StatusCode HTopMultilepMiniNTupMaker :: fillTPFlatBranches ( std::shared_ptr<MiniNTupMaker::leptonObj> lep, const std::string& trig ) {

    bool isTag(false);
    bool isBadTPEvent(false); // be optimistic!

    std::string this_tp_trig("");
    std::vector<std::string> tp_trigs;

    if ( trig.compare("SLT") == 0 ) {
	isTag = lep.get()->props["isTagSLT"].c;
	this_tp_trig = ( isTag ) ?  "Tag_SLT" : "Probe_SLT";
	isBadTPEvent = m_event.get()->isBadTPEvent_SLT;
	if ( isBadTPEvent ) {
	    tp_trigs.push_back("Tag_SLT");
	    tp_trigs.push_back("Probe_SLT");
	} else {
	    tp_trigs.push_back(this_tp_trig);
	}
    } else {
	Error("fillTPFlatBranches()","Invalid parameter: %s for this function!", trig.c_str() );
	return EL::StatusCode::FAILURE;
    }

    // Keep the default branch values if the event is bad

    if ( isBadTPEvent ) { return EL::StatusCode::SUCCESS; }

    std::string key, tp_key, type;

    for ( const auto& tp_trig : tp_trigs ) {

	for ( const auto& var : m_TP_VARS ) {

	    if ( var.find("VEC") != std::string::npos ) { continue; } // only flat branches here!

	    key  = var.substr( 0, var.length() - 2 );
	    type = var.substr( var.length() - 1 );

	    tp_key = "lep_" + tp_trig + "_" + key;

	    if ( type.compare("F") == 0 ) { m_TagProbe_branches[tp_key].f = lep.get()->props[key].f; }
	    if ( type.compare("B") == 0 ) { m_TagProbe_branches[tp_key].c = lep.get()->props[key].c; }
	    if ( type.compare("I") == 0 ) { m_TagProbe_branches[tp_key].i = lep.get()->props[key].i; }

	}

    }

    return EL::StatusCode::SUCCESS;

}

EL::StatusCode HTopMultilepMiniNTupMaker :: storeLeptonBranches()
{

    // Clear vector branches from previous event

    ANA_CHECK( this->clearBranches("leptons") );

    m_nleptons = m_leptons.size();

    // Some flat branches

    std::string key;

    for ( const auto& var : m_LEP_OUTPUT_VARS ) {

	key = var.substr( 0, var.length() - 2 );

	for ( unsigned int idx(0); idx < m_leptons.size(); ++idx ) {

	    switch( idx ) {
	    case 0:
		m_lep0_OUTPUT_branches[key].f = m_leptons.at(idx).get()->props[key].f;
		m_lep0_OUTPUT_branches[key].c = m_leptons.at(idx).get()->props[key].c;
		m_lep0_OUTPUT_branches[key].i = m_leptons.at(idx).get()->props[key].i;
		break;
	    case 1:
		m_lep1_OUTPUT_branches[key].f = m_leptons.at(idx).get()->props[key].f;
		m_lep1_OUTPUT_branches[key].c = m_leptons.at(idx).get()->props[key].c;
		m_lep1_OUTPUT_branches[key].i = m_leptons.at(idx).get()->props[key].i;
		break;
	    case 2:
		m_lep2_OUTPUT_branches[key].f = m_leptons.at(idx).get()->props[key].f;
		m_lep2_OUTPUT_branches[key].c = m_leptons.at(idx).get()->props[key].c;
		m_lep2_OUTPUT_branches[key].i = m_leptons.at(idx).get()->props[key].i;
		break;
	    default:
		break;
	    }
	}
    }

    int lep_flavour;

    m_nmuons     = 0;
    m_nelectrons = 0;

    std::string electron_key, muon_key;

    for ( auto lep : m_leptons ) {

	lep_flavour = lep.get()->props["Flavour"].i;

	// Fill electron and muon post OLR vector branches

	if ( lep_flavour == 11 ) {

	    ++m_nelectrons;

	    for ( const auto& var : m_EL_VEC_VARS ) {

		key = var.substr( 0, var.length() - 2 );
		electron_key = "electron_" + key;

		m_electron_OR_branches[electron_key].vec_f.push_back(lep.get()->props[key].f);
		m_electron_OR_branches[electron_key].vec_c.push_back(lep.get()->props[key].c);
		m_electron_OR_branches[electron_key].vec_i.push_back(lep.get()->props[key].i);

	    }

	} else if ( lep_flavour == 13 ) {

	    ++m_nmuons;

	    for ( const auto& var : m_MU_VEC_VARS ) {

		key = var.substr( 0, var.length() - 2 );
		muon_key = "muon_" + key;

		m_muon_OR_branches[muon_key].vec_f.push_back(lep.get()->props[key].f);
		m_muon_OR_branches[muon_key].vec_c.push_back(lep.get()->props[key].c);
		m_muon_OR_branches[muon_key].vec_i.push_back(lep.get()->props[key].i);

	    }

	}

    }

    return EL::StatusCode::SUCCESS;

}

bool HTopMultilepMiniNTupMaker :: isQMisIDBDTLoose( std::shared_ptr<leptonObj> lepA, std::shared_ptr<leptonObj> lepB ) {

    if ( m_lepSelForTP.compare("MVA") != 0 ) { return true; } // Forget about QMisID BDT for ICHEP selection!
    if ( m_event.get()->trilep_type )        { return true; } // No QMisID at all in 3L
    if ( m_event.get()->dilep_type == 1 )    { return true; } // No QMisID at all in 2L mm
    if ( m_event.get()->dilep_type == 3 ) { // 2L ee
	if ( lepA.get()->props["chargeIDBDTLoose"].f > 0.10083 && lepB.get()->props["chargeIDBDTLoose"].f > 0.10083 ) { return true; }
    }
    if ( m_event.get()->dilep_type == 2 ) { // 2L OF
	if ( lepA.get()->props["Flavour"].i == 11 && lepA.get()->props["chargeIDBDTLoose"].f > 0.10083 ) { return true; }
	if ( lepB.get()->props["Flavour"].i == 11 && lepB.get()->props["chargeIDBDTLoose"].f > 0.10083 ) { return true; }
    }
    return false;
}


EL::StatusCode HTopMultilepMiniNTupMaker :: setOutputBranches ()
{

    m_isSS01 = m_event.get()->isSS01;
    m_isSS12 = m_event.get()->isSS12;

    // std::shared_ptr<MiniNTupMaker::leptonObj> lepA, lepB;
    auto lepA = m_leptons.at(0);
    auto lepB = m_leptons.at(1);

    if ( m_event.get()->trilep_type ) {
	lepA = m_leptons.at(1);
	lepB = m_leptons.at(2);
    }

    m_is_T_T	     = (  lepA.get()->props["isTightSelected"].c &&  lepB.get()->props["isTightSelected"].c );
    m_is_T_AntiT     = (  lepA.get()->props["isTightSelected"].c && !lepB.get()->props["isTightSelected"].c );
    m_is_AntiT_T     = ( !lepA.get()->props["isTightSelected"].c &&  lepB.get()->props["isTightSelected"].c );
    m_is_AntiT_AntiT = ( !lepA.get()->props["isTightSelected"].c && !lepB.get()->props["isTightSelected"].c );

    m_is_Tel_AntiTmu = ( ( lepA.get()->props["Flavour"].i == 11 && lepA.get()->props["isTightSelected"].c )  && ( lepB.get()->props["Flavour"].i == 13 && !lepB.get()->props["isTightSelected"].c ) );
    m_is_Tmu_AntiTel = ( ( lepA.get()->props["Flavour"].i == 13 && lepA.get()->props["isTightSelected"].c )  && ( lepB.get()->props["Flavour"].i == 11 && !lepB.get()->props["isTightSelected"].c ) );
    m_is_AntiTel_Tmu = ( ( lepA.get()->props["Flavour"].i == 11 && !lepA.get()->props["isTightSelected"].c ) && ( lepB.get()->props["Flavour"].i == 13 && lepB.get()->props["isTightSelected"].c ) );
    m_is_AntiTmu_Tel = ( ( lepA.get()->props["Flavour"].i == 13 && !lepA.get()->props["isTightSelected"].c ) && ( lepB.get()->props["Flavour"].i == 11 && lepB.get()->props["isTightSelected"].c ) );

    // For lepton MVA selection, apply a baseline requirement for *all* electrons to be passing QMisID BDT Loose cut

    bool isQMIsIDBDTLooseFlag = this->isQMisIDBDTLoose( lepA, lepB );

    m_is_TMVA_TMVA	   = (  lepA.get()->props["isTightSelectedMVA"].c &&  lepB.get()->props["isTightSelectedMVA"].c ) && isQMIsIDBDTLooseFlag;
    m_is_TMVA_AntiTMVA     = (  lepA.get()->props["isTightSelectedMVA"].c && !lepB.get()->props["isTightSelectedMVA"].c ) && isQMIsIDBDTLooseFlag;
    m_is_AntiTMVA_TMVA     = ( !lepA.get()->props["isTightSelectedMVA"].c &&  lepB.get()->props["isTightSelectedMVA"].c ) && isQMIsIDBDTLooseFlag;
    m_is_AntiTMVA_AntiTMVA = ( !lepA.get()->props["isTightSelectedMVA"].c && !lepB.get()->props["isTightSelectedMVA"].c ) && isQMIsIDBDTLooseFlag;

    m_is_TMVAel_AntiTMVAmu = ( ( lepA.get()->props["Flavour"].i == 11 && lepA.get()->props["isTightSelectedMVA"].c )  && ( lepB.get()->props["Flavour"].i == 13 && !lepB.get()->props["isTightSelectedMVA"].c ) ) && isQMIsIDBDTLooseFlag;
    m_is_TMVAmu_AntiTMVAel = ( ( lepA.get()->props["Flavour"].i == 13 && lepA.get()->props["isTightSelectedMVA"].c )  && ( lepB.get()->props["Flavour"].i == 11 && !lepB.get()->props["isTightSelectedMVA"].c ) ) && isQMIsIDBDTLooseFlag;
    m_is_AntiTMVAel_TMVAmu = ( ( lepA.get()->props["Flavour"].i == 11 && !lepA.get()->props["isTightSelectedMVA"].c ) && ( lepB.get()->props["Flavour"].i == 13 && lepB.get()->props["isTightSelectedMVA"].c ) ) && isQMIsIDBDTLooseFlag;
    m_is_AntiTMVAmu_TMVAel = ( ( lepA.get()->props["Flavour"].i == 13 && !lepA.get()->props["isTightSelectedMVA"].c ) && ( lepB.get()->props["Flavour"].i == 11 && lepB.get()->props["isTightSelectedMVA"].c ) ) && isQMIsIDBDTLooseFlag;

    ANA_CHECK( this->storeLeptonBranches() );

    m_isBadTPEvent_SLT = m_event.get()->isBadTPEvent_SLT;

    for ( const auto lep : m_leptons ) {

	// Fill flat variables for tag/probe leptons

	if ( m_event.get()->dilep_type ) {
	    ANA_CHECK( this->fillTPFlatBranches( lep, "SLT" ) );
	}

    }

    m_weight_event            = m_event.get()->weight_event;
    m_weight_event_lep        = m_event.get()->weight_event_lep;
    m_weight_event_trig_SLT   = m_event.get()->weight_event_trig_SLT;

    m_weight_lep_tag          = m_event.get()->weight_lep_tag;
    m_weight_trig_tag         = m_event.get()->weight_trig_tag;
    m_weight_lep_probe        = m_event.get()->weight_lep_probe;
    m_weight_trig_probe       = m_event.get()->weight_trig_probe;

    return EL::StatusCode::SUCCESS;

}

EL::StatusCode HTopMultilepMiniNTupMaker :: clearBranches ( const std::string& type )
{

    if ( type.compare("leptons") == 0 ) {

	for ( auto& elem : m_electron_OR_branches ) {
	    elem.second.vec_f.clear();
	    elem.second.vec_c.clear();
	    elem.second.vec_i.clear();
	}
	for ( auto& elem : m_muon_OR_branches ) {
	    elem.second.vec_f.clear();
	    elem.second.vec_c.clear();
	    elem.second.vec_i.clear();
	}

    } else if ( type.compare("tag_and_probe") == 0 ) {

	for ( auto& elem : m_TagProbe_branches ) {
	    elem.second.vec_f.clear();
	    elem.second.vec_c.clear();
	    elem.second.vec_i.clear();
	}

    } else if ( type.compare("jets_kin") == 0 ) {
	m_jet_OR_Pt.clear();
	m_jet_OR_Eta.clear();
	m_jet_OR_Phi.clear();
	m_jet_OR_E.clear();
    } else if ( type.compare("jets_truth") == 0 ) {
	m_jet_OR_truthMatch_Pt.clear();
	m_jet_OR_truthMatch_Eta.clear();
	m_jet_OR_truthMatch_Phi.clear();
	m_jet_OR_truthMatch_E.clear();
	m_jet_OR_truthMatch_isBJet.clear();
	m_jet_OR_truthMatch_isCJet.clear();
	m_jet_OR_truthMatch_isLFJet.clear();
	m_jet_OR_truthMatch_isGluonJet.clear();
    } else {
	Error("clearBranches()","Trying to clear branches of a particle type that is not known. Aborting.");
	return EL::StatusCode::FAILURE;
    }

    return EL::StatusCode::SUCCESS;

}

EL::StatusCode HTopMultilepMiniNTupMaker :: jetKinematics ()
{

    // Clear vector branches from previous event
    //
    ANA_CHECK( this->clearBranches("jets_kin") );

    for ( unsigned int j_idx(0); j_idx < m_selected_jets_T->size(); ++j_idx ) {

	short j = m_selected_jets_T->at(j_idx);

	// Read the jet branches from the vector *before* OR via the index

	if ( m_debug ) { Info("jetKinematics()","reco jet[%i], pT = %.2f, eta = %.2f, phi = %.2f, E = %.2f",j, m_jet_pt->at(j)/1e3, m_jet_eta->at(j), m_jet_phi->at(j), m_jet_E->at(j)/1e3 ); }

	m_jet_OR_Pt.push_back( m_jet_pt->at(j) );
	m_jet_OR_Eta.push_back( m_jet_eta->at(j) );
	m_jet_OR_Phi.push_back( m_jet_phi->at(j) );
	m_jet_OR_E.push_back( m_jet_E->at(j) );

    }

    return EL::StatusCode::SUCCESS;

}

EL::StatusCode HTopMultilepMiniNTupMaker :: jetTruthMatching ()
{

    // Clear vector branches from previous event
    //
    ANA_CHECK( this->clearBranches("jets_truth") );

    float DRCONE(0.4);

    // Find the truth-matching jet via min(deltaR) matching

    TVector3 truthjet, jet;

    for ( unsigned int j_idx(0); j_idx < m_selected_jets_T->size(); ++j_idx ) {

	short j = m_selected_jets_T->at(j_idx);
	// Read the jet branches from the vector *before* OR via the index
	//
	jet.SetPtEtaPhi( m_jet_pt->at(j), m_jet_eta->at(j), m_jet_phi->at(j) );

	if ( m_debug ) { Info("jetTruthMatching()","reco jet[%i], pT = %.2f, eta = %.2f, phi = %.2f",j, jet.Pt()/1e3, jet.Eta(), jet.Phi() ); }

	float minDR(999.0), thisDR(999.0);
	int best_tj(-1);

	for ( unsigned int tj(0); tj < m_truth_jet_pt->size() ; ++tj ) {

	    truthjet.SetPtEtaPhi( m_truth_jet_pt->at(tj), m_truth_jet_eta->at(tj), m_truth_jet_phi->at(tj) );

	    thisDR = jet.DeltaR(truthjet);

	    if ( m_debug ) { Info("jetTruthMatching()","\ttruth jet - pT = %.2f, eta = %.2f, phi = %.2f, DR = %.2f", truthjet.Pt()/1e3, truthjet.Eta(), truthjet.Phi(), thisDR ); }

	    if ( thisDR < DRCONE && thisDR < minDR ) {
		minDR = thisDR;
		best_tj = tj;
	    }

	}

	float match_tj_pt(-1.0), match_tj_eta(-999.0), match_tj_phi(-999.0), match_tj_e(-1.0);
	char match_tj_isbjet(-1), match_tj_iscjet(-1), match_tj_islfjet(-1), match_tj_isgluonjet(-1);

	if ( best_tj != -1 ) {

	    match_tj_pt  =  m_truth_jet_pt->at(best_tj);
	    match_tj_eta =  m_truth_jet_eta->at(best_tj);
	    match_tj_phi =  m_truth_jet_phi->at(best_tj);
	    match_tj_e   =  m_truth_jet_e->at(best_tj);

	    // Truth flavour classification already done w/ MCTruthClassifier in DF
	    //
	    match_tj_isbjet  = ( m_jet_flavor_truth_label_ghost->at(j) == 5 );
	    match_tj_iscjet  = ( m_jet_flavor_truth_label_ghost->at(j) == 4 );
	    match_tj_islfjet = ( m_jet_flavor_truth_label_ghost->at(j) != 5 && m_jet_flavor_truth_label_ghost->at(j) != 4 && m_jet_flavor_truth_label_ghost->at(j) != 21 );
	    match_tj_isgluonjet  = ( m_jet_flavor_truth_label_ghost->at(j) == 21 );

	    if ( m_debug ) { Info("jetTruthMatching()","matching truth jet found - idx = %i, pT = %.2f, eta = %.2f, phi = %.2f, DR = %.2f\nisBJet? = %i\nisCJet? = %i\nisLFJet? = %i\nisGluonJet? = %i", best_tj, match_tj_pt/1e3, match_tj_eta, match_tj_phi, minDR, match_tj_isbjet, match_tj_iscjet, match_tj_islfjet, match_tj_isgluonjet); }

	} else {
	    if ( m_debug ) { Info("jetTruthMatching()","matching reco jet NOT found! Setting dummy truth match branches"); }
	}

	m_jet_OR_truthMatch_Pt.push_back( match_tj_pt );
	m_jet_OR_truthMatch_Eta.push_back( match_tj_eta );
	m_jet_OR_truthMatch_Phi.push_back( match_tj_phi );
	m_jet_OR_truthMatch_E.push_back( match_tj_e );
	m_jet_OR_truthMatch_isBJet.push_back( match_tj_isbjet );
	m_jet_OR_truthMatch_isCJet.push_back( match_tj_iscjet );
	m_jet_OR_truthMatch_isLFJet.push_back( match_tj_islfjet );
	m_jet_OR_truthMatch_isGluonJet.push_back( match_tj_isgluonjet );

    }

    return EL::StatusCode::SUCCESS;

}
