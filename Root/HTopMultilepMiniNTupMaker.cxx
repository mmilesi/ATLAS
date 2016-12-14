#include <HTopMultilepAnalysis/HTopMultilepMiniNTupMaker.h>

// ASG status code check
#include <AsgTools/MessageCheck.h>

// ROOT include(s)
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

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
  m_useSUSYSSTP          = false;
  m_ambiSolvingCrit      = "Pt";
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
  //
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
  //
  m_inputNTuple = wk()->tree();

  ANA_CHECK( this->enableSelectedBranches() );

  // Connect the branches of the input tree to the algorithm members
  //
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

  m_inputNTuple->SetBranchAddress ("nJets_OR", 			              &m_nJets_OR);
  m_inputNTuple->SetBranchAddress ("nJets_OR_MV2c10_70", 		      &m_nJets_OR_MV2c10_70);
  m_inputNTuple->SetBranchAddress ("nJets_OR_T", 			      &m_nJets_OR_T);
  m_inputNTuple->SetBranchAddress ("nJets_OR_T_MV2c10_70", 		      &m_nJets_OR_T_MV2c10_70);

  m_inputNTuple->SetBranchAddress ("lep_ID_0",   			      &m_lep_ID_0);
  m_inputNTuple->SetBranchAddress ("lep_Pt_0",  			      &m_lep_Pt_0);
  m_inputNTuple->SetBranchAddress ("lep_E_0",   			      &m_lep_E_0);
  m_inputNTuple->SetBranchAddress ("lep_Eta_0",  			      &m_lep_Eta_0);
  m_inputNTuple->SetBranchAddress ("lep_Phi_0",   			      &m_lep_Phi_0);
  m_inputNTuple->SetBranchAddress ("lep_EtaBE2_0",   			      &m_lep_EtaBE2_0);
  m_inputNTuple->SetBranchAddress ("lep_sigd0PV_0",   			      &m_lep_sigd0PV_0);
  m_inputNTuple->SetBranchAddress ("lep_Z0SinTheta_0",   		      &m_lep_Z0SinTheta_0);
  m_inputNTuple->SetBranchAddress ("lep_isTightLH_0",   		      &m_lep_isTightLH_0);
  m_inputNTuple->SetBranchAddress ("lep_isMediumLH_0",   		      &m_lep_isMediumLH_0);
  m_inputNTuple->SetBranchAddress ("lep_isLooseLH_0",   		      &m_lep_isLooseLH_0);
  m_inputNTuple->SetBranchAddress ("lep_isTight_0",   			      &m_lep_isTight_0);
  m_inputNTuple->SetBranchAddress ("lep_isMedium_0",   			      &m_lep_isMedium_0);
  m_inputNTuple->SetBranchAddress ("lep_isLoose_0",   			      &m_lep_isLoose_0);
  m_inputNTuple->SetBranchAddress ("lep_isolationLooseTrackOnly_0",   	      &m_lep_isolationLooseTrackOnly_0);
  m_inputNTuple->SetBranchAddress ("lep_isolationLoose_0",   		      &m_lep_isolationLoose_0);
  m_inputNTuple->SetBranchAddress ("lep_isolationFixedCutTight_0",   	      &m_lep_isolationFixedCutTight_0);
  m_inputNTuple->SetBranchAddress ("lep_isolationFixedCutTightTrackOnly_0",   &m_lep_isolationFixedCutTightTrackOnly_0);
  m_inputNTuple->SetBranchAddress ("lep_isolationFixedCutLoose_0",   	      &m_lep_isolationFixedCutLoose_0);
  m_inputNTuple->SetBranchAddress ("lep_topoEtcone20_0",                      &m_lep_topoEtcone20_0);
  m_inputNTuple->SetBranchAddress ("lep_ptVarcone20_0",                       &m_lep_ptVarcone20_0);
  m_inputNTuple->SetBranchAddress ("lep_ptVarcone30_0",                       &m_lep_ptVarcone30_0);
  m_inputNTuple->SetBranchAddress ("lep_isTrigMatch_0",   		      &m_lep_isTrigMatch_0);
  m_inputNTuple->SetBranchAddress ("lep_isTrigMatchDLT_0",   		      &m_lep_isTrigMatchDLT_0);
  m_inputNTuple->SetBranchAddress ("lep_isPrompt_0",   			      &m_lep_isPrompt_0);
  m_inputNTuple->SetBranchAddress ("lep_isBrems_0",   		              &m_lep_isBrems_0);
  m_inputNTuple->SetBranchAddress ("lep_isFakeLep_0",   		      &m_lep_isFakeLep_0);
  m_inputNTuple->SetBranchAddress ("lep_isQMisID_0",   		              &m_lep_isQMisID_0);
  m_inputNTuple->SetBranchAddress ("lep_isConvPh_0",   		              &m_lep_isConvPh_0);
  m_inputNTuple->SetBranchAddress ("lep_truthType_0",   		      &m_lep_truthType_0);
  m_inputNTuple->SetBranchAddress ("lep_truthOrigin_0",   		      &m_lep_truthOrigin_0);
  m_inputNTuple->SetBranchAddress ("lep_SFIDLoose_0",   		      &m_lep_SFIDLoose_0);
  m_inputNTuple->SetBranchAddress ("lep_SFIDTight_0",   		      &m_lep_SFIDTight_0);
  m_inputNTuple->SetBranchAddress ("lep_SFTrigLoose_0",   		      &m_lep_SFTrigLoose_0);
  m_inputNTuple->SetBranchAddress ("lep_SFTrigTight_0",   		      &m_lep_SFTrigTight_0);
  m_inputNTuple->SetBranchAddress ("lep_EffTrigLoose_0",                      &m_lep_EffTrigLoose_0);
  m_inputNTuple->SetBranchAddress ("lep_EffTrigTight_0",                      &m_lep_EffTrigTight_0);
  m_inputNTuple->SetBranchAddress ("lep_SFIsoLoose_0",   		      &m_lep_SFIsoLoose_0);
  m_inputNTuple->SetBranchAddress ("lep_SFIsoTight_0",   		      &m_lep_SFIsoTight_0);
  m_inputNTuple->SetBranchAddress ("lep_SFReco_0",   			      &m_lep_SFReco_0);
  m_inputNTuple->SetBranchAddress ("lep_SFTTVA_0",   			      &m_lep_SFTTVA_0);
  m_inputNTuple->SetBranchAddress ("lep_SFObjLoose_0",   		      &m_lep_SFObjLoose_0);
  m_inputNTuple->SetBranchAddress ("lep_SFObjTight_0",   		      &m_lep_SFObjTight_0);

  m_inputNTuple->SetBranchAddress ("lep_ID_1",   			      &m_lep_ID_1);
  m_inputNTuple->SetBranchAddress ("lep_Pt_1",  			      &m_lep_Pt_1);
  m_inputNTuple->SetBranchAddress ("lep_E_1",   			      &m_lep_E_1);
  m_inputNTuple->SetBranchAddress ("lep_Eta_1",  			      &m_lep_Eta_1);
  m_inputNTuple->SetBranchAddress ("lep_Phi_1",   			      &m_lep_Phi_1);
  m_inputNTuple->SetBranchAddress ("lep_EtaBE2_1",   			      &m_lep_EtaBE2_1);
  m_inputNTuple->SetBranchAddress ("lep_sigd0PV_1",   			      &m_lep_sigd0PV_1);
  m_inputNTuple->SetBranchAddress ("lep_Z0SinTheta_1",   		      &m_lep_Z0SinTheta_1);
  m_inputNTuple->SetBranchAddress ("lep_isTightLH_1",   		      &m_lep_isTightLH_1);
  m_inputNTuple->SetBranchAddress ("lep_isMediumLH_1",   		      &m_lep_isMediumLH_1);
  m_inputNTuple->SetBranchAddress ("lep_isLooseLH_1",   		      &m_lep_isLooseLH_1);
  m_inputNTuple->SetBranchAddress ("lep_isTight_1",   			      &m_lep_isTight_1);
  m_inputNTuple->SetBranchAddress ("lep_isMedium_1",   			      &m_lep_isMedium_1);
  m_inputNTuple->SetBranchAddress ("lep_isLoose_1",   			      &m_lep_isLoose_1);
  m_inputNTuple->SetBranchAddress ("lep_isolationLooseTrackOnly_1",   	      &m_lep_isolationLooseTrackOnly_1);
  m_inputNTuple->SetBranchAddress ("lep_isolationLoose_1",   		      &m_lep_isolationLoose_1);
  m_inputNTuple->SetBranchAddress ("lep_isolationFixedCutTight_1",   	      &m_lep_isolationFixedCutTight_1);
  m_inputNTuple->SetBranchAddress ("lep_isolationFixedCutTightTrackOnly_1",   &m_lep_isolationFixedCutTightTrackOnly_1);
  m_inputNTuple->SetBranchAddress ("lep_isolationFixedCutLoose_1",   	      &m_lep_isolationFixedCutLoose_1);
  m_inputNTuple->SetBranchAddress ("lep_topoEtcone20_1",                      &m_lep_topoEtcone20_1);
  m_inputNTuple->SetBranchAddress ("lep_ptVarcone20_1",                       &m_lep_ptVarcone20_1);
  m_inputNTuple->SetBranchAddress ("lep_ptVarcone30_1",                       &m_lep_ptVarcone30_1);
  m_inputNTuple->SetBranchAddress ("lep_isTrigMatch_1",   		      &m_lep_isTrigMatch_1);
  m_inputNTuple->SetBranchAddress ("lep_isTrigMatchDLT_1",   		      &m_lep_isTrigMatchDLT_1);
  m_inputNTuple->SetBranchAddress ("lep_isPrompt_1",   			      &m_lep_isPrompt_1);
  m_inputNTuple->SetBranchAddress ("lep_isBrems_1",   		              &m_lep_isBrems_1);
  m_inputNTuple->SetBranchAddress ("lep_isQMisID_1",   		              &m_lep_isQMisID_1);
  m_inputNTuple->SetBranchAddress ("lep_isConvPh_1",   		              &m_lep_isConvPh_1);
  m_inputNTuple->SetBranchAddress ("lep_truthType_1",   		      &m_lep_truthType_1);
  m_inputNTuple->SetBranchAddress ("lep_truthOrigin_1",   		      &m_lep_truthOrigin_1);
  m_inputNTuple->SetBranchAddress ("lep_isFakeLep_1",   		      &m_lep_isFakeLep_1);
  m_inputNTuple->SetBranchAddress ("lep_SFIDLoose_1",   		      &m_lep_SFIDLoose_1);
  m_inputNTuple->SetBranchAddress ("lep_SFIDTight_1",   		      &m_lep_SFIDTight_1);
  m_inputNTuple->SetBranchAddress ("lep_SFTrigLoose_1",   		      &m_lep_SFTrigLoose_1);
  m_inputNTuple->SetBranchAddress ("lep_SFTrigTight_1",   		      &m_lep_SFTrigTight_1);
  m_inputNTuple->SetBranchAddress ("lep_EffTrigLoose_1",                      &m_lep_EffTrigLoose_1);
  m_inputNTuple->SetBranchAddress ("lep_EffTrigTight_1",                      &m_lep_EffTrigTight_1);
  m_inputNTuple->SetBranchAddress ("lep_SFIsoLoose_1",   		      &m_lep_SFIsoLoose_1);
  m_inputNTuple->SetBranchAddress ("lep_SFIsoTight_1",   		      &m_lep_SFIsoTight_1);
  m_inputNTuple->SetBranchAddress ("lep_SFReco_1",   			      &m_lep_SFReco_1);
  m_inputNTuple->SetBranchAddress ("lep_SFTTVA_1",   			      &m_lep_SFTTVA_1);
  m_inputNTuple->SetBranchAddress ("lep_SFObjLoose_1",   		      &m_lep_SFObjLoose_1);
  m_inputNTuple->SetBranchAddress ("lep_SFObjTight_1",   		      &m_lep_SFObjTight_1);

  m_inputNTuple->SetBranchAddress ("lep_ID_2",   			      &m_lep_ID_2);
  m_inputNTuple->SetBranchAddress ("lep_Pt_2",  			      &m_lep_Pt_2);
  m_inputNTuple->SetBranchAddress ("lep_E_2",   			      &m_lep_E_2);
  m_inputNTuple->SetBranchAddress ("lep_Eta_2",  			      &m_lep_Eta_2);
  m_inputNTuple->SetBranchAddress ("lep_Phi_2",   			      &m_lep_Phi_2);
  m_inputNTuple->SetBranchAddress ("lep_EtaBE2_2",   			      &m_lep_EtaBE2_2);
  m_inputNTuple->SetBranchAddress ("lep_sigd0PV_2",   			      &m_lep_sigd0PV_2);
  m_inputNTuple->SetBranchAddress ("lep_Z0SinTheta_2",   		      &m_lep_Z0SinTheta_2);
  m_inputNTuple->SetBranchAddress ("lep_isTightLH_2",   		      &m_lep_isTightLH_2);
  m_inputNTuple->SetBranchAddress ("lep_isMediumLH_2",   		      &m_lep_isMediumLH_2);
  m_inputNTuple->SetBranchAddress ("lep_isLooseLH_2",   		      &m_lep_isLooseLH_2);
  m_inputNTuple->SetBranchAddress ("lep_isTight_2",   			      &m_lep_isTight_2);
  m_inputNTuple->SetBranchAddress ("lep_isMedium_2",   			      &m_lep_isMedium_2);
  m_inputNTuple->SetBranchAddress ("lep_isLoose_2",   			      &m_lep_isLoose_2);
  m_inputNTuple->SetBranchAddress ("lep_isolationLooseTrackOnly_2",   	      &m_lep_isolationLooseTrackOnly_2);
  m_inputNTuple->SetBranchAddress ("lep_isolationLoose_2",   		      &m_lep_isolationLoose_2);
  m_inputNTuple->SetBranchAddress ("lep_isolationFixedCutTight_2",   	      &m_lep_isolationFixedCutTight_2);
  m_inputNTuple->SetBranchAddress ("lep_isolationFixedCutTightTrackOnly_2",   &m_lep_isolationFixedCutTightTrackOnly_2);
  m_inputNTuple->SetBranchAddress ("lep_isolationFixedCutLoose_2",   	      &m_lep_isolationFixedCutLoose_2);
  m_inputNTuple->SetBranchAddress ("lep_topoEtcone20_2",                      &m_lep_topoEtcone20_2);
  m_inputNTuple->SetBranchAddress ("lep_ptVarcone20_2",                       &m_lep_ptVarcone20_2);
  m_inputNTuple->SetBranchAddress ("lep_ptVarcone30_2",                       &m_lep_ptVarcone30_2);
  m_inputNTuple->SetBranchAddress ("lep_isTrigMatch_2",   		      &m_lep_isTrigMatch_2);
  m_inputNTuple->SetBranchAddress ("lep_isTrigMatchDLT_2",   		      &m_lep_isTrigMatchDLT_2);
  m_inputNTuple->SetBranchAddress ("lep_isPrompt_2",   			      &m_lep_isPrompt_2);
  m_inputNTuple->SetBranchAddress ("lep_isBrems_2",   		              &m_lep_isBrems_2);
  m_inputNTuple->SetBranchAddress ("lep_isFakeLep_2",   		      &m_lep_isFakeLep_2);
  m_inputNTuple->SetBranchAddress ("lep_isQMisID_2",   		              &m_lep_isQMisID_2);
  m_inputNTuple->SetBranchAddress ("lep_isConvPh_2",   		              &m_lep_isConvPh_2);
  m_inputNTuple->SetBranchAddress ("lep_truthType_2",   		      &m_lep_truthType_2);
  m_inputNTuple->SetBranchAddress ("lep_truthOrigin_2",   		      &m_lep_truthOrigin_2);
  m_inputNTuple->SetBranchAddress ("lep_SFIDLoose_2",   		      &m_lep_SFIDLoose_2);
  m_inputNTuple->SetBranchAddress ("lep_SFIDTight_2",   		      &m_lep_SFIDTight_2);
  m_inputNTuple->SetBranchAddress ("lep_SFTrigLoose_2",   		      &m_lep_SFTrigLoose_2);
  m_inputNTuple->SetBranchAddress ("lep_SFTrigTight_2",   		      &m_lep_SFTrigTight_2);
  m_inputNTuple->SetBranchAddress ("lep_EffTrigLoose_2",                      &m_lep_EffTrigLoose_2);
  m_inputNTuple->SetBranchAddress ("lep_EffTrigTight_2",                      &m_lep_EffTrigTight_2);
  m_inputNTuple->SetBranchAddress ("lep_SFIsoLoose_2",   		      &m_lep_SFIsoLoose_2);
  m_inputNTuple->SetBranchAddress ("lep_SFIsoTight_2",   		      &m_lep_SFIsoTight_2);
  m_inputNTuple->SetBranchAddress ("lep_SFReco_2",   			      &m_lep_SFReco_2);
  m_inputNTuple->SetBranchAddress ("lep_SFTTVA_2",   			      &m_lep_SFTTVA_2);
  m_inputNTuple->SetBranchAddress ("lep_SFObjLoose_2",   		      &m_lep_SFObjLoose_2);
  m_inputNTuple->SetBranchAddress ("lep_SFObjTight_2",   		      &m_lep_SFObjTight_2);

  m_inputNTuple->SetBranchAddress ("electron_match_HLT_e26_lhtight_nod0_ivarloose",   &m_electron_match_HLT_e26_lhtight_nod0_ivarloose);
  m_inputNTuple->SetBranchAddress ("electron_match_HLT_e24_lhmedium_L1EM20VH",        &m_electron_match_HLT_e24_lhmedium_L1EM20VH);
  m_inputNTuple->SetBranchAddress ("electron_match_HLT_e60_lhmedium",		      &m_electron_match_HLT_e60_lhmedium);
  m_inputNTuple->SetBranchAddress ("electron_match_HLT_e120_lhloose",		      &m_electron_match_HLT_e120_lhloose);
  m_inputNTuple->SetBranchAddress ("electron_match_HLT_2e12_lhloose_L12EM10VH",       &m_electron_match_HLT_2e12_lhloose_L12EM10VH);
  m_inputNTuple->SetBranchAddress ("electron_match_HLT_e24_medium_L1EM20VHI_mu8noL1", &m_electron_match_HLT_e24_medium_L1EM20VHI_mu8noL1);
  m_inputNTuple->SetBranchAddress ("electron_match_HLT_e7_medium_mu24", 	      &m_electron_match_HLT_e7_medium_mu24);
  m_inputNTuple->SetBranchAddress ("muon_match_HLT_mu20_iloose_L1MU15", 	      &m_muon_match_HLT_mu20_iloose_L1MU15);
  m_inputNTuple->SetBranchAddress ("muon_match_HLT_mu50",			      &m_muon_match_HLT_mu50);
  m_inputNTuple->SetBranchAddress ("muon_match_HLT_mu18_mu8noL1",		      &m_muon_match_HLT_mu18_mu8noL1);
  m_inputNTuple->SetBranchAddress ("muon_match_HLT_e24_medium_L1EM20VHI_mu8noL1",     &m_muon_match_HLT_e24_medium_L1EM20VHI_mu8noL1);
  m_inputNTuple->SetBranchAddress ("muon_match_HLT_e7_medium_mu24",                   &m_muon_match_HLT_e7_medium_mu24);

  m_inputNTuple->SetBranchAddress ("electron_match_HLT_e26_lhtight_nod0_ivarloose",   &m_electron_match_HLT_e26_lhtight_nod0_ivarloose);
  m_inputNTuple->SetBranchAddress ("electron_match_HLT_e60_lhmedium_nod0",	      &m_electron_match_HLT_e60_lhmedium_nod0);
  m_inputNTuple->SetBranchAddress ("electron_match_HLT_e140_lhloose_nod0",	      &m_electron_match_HLT_e140_lhloose_nod0);
  m_inputNTuple->SetBranchAddress ("electron_match_HLT_2e17_lhvloose_nod0",	      &m_electron_match_HLT_2e17_lhvloose_nod0);
  m_inputNTuple->SetBranchAddress ("electron_match_HLT_e17_lhloose_mu14",             &m_electron_match_HLT_e17_lhloose_mu14);
  m_inputNTuple->SetBranchAddress ("electron_match_HLT_e17_lhloose_nod0_mu14",        &m_electron_match_HLT_e17_lhloose_nod0_mu14);
  m_inputNTuple->SetBranchAddress ("electron_match_HLT_e7_lhmedium_mu24",	      &m_electron_match_HLT_e7_lhmedium_mu24);
  m_inputNTuple->SetBranchAddress ("muon_match_HLT_mu26_ivarmedium",		      &m_muon_match_HLT_mu26_ivarmedium);
  m_inputNTuple->SetBranchAddress ("muon_match_HLT_mu22_mu8noL1",		      &m_muon_match_HLT_mu22_mu8noL1);
  m_inputNTuple->SetBranchAddress ("muon_match_HLT_e17_lhloose_mu14",		      &m_muon_match_HLT_e17_lhloose_mu14);
  m_inputNTuple->SetBranchAddress ("muon_match_HLT_e17_lhloose_nod0_mu14",	      &m_muon_match_HLT_e17_lhloose_nod0_mu14);

  m_inputNTuple->SetBranchAddress ("electron_passOR",        &m_electron_passOR);
  m_inputNTuple->SetBranchAddress ("muon_passOR",            &m_muon_passOR);

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
  //
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
  //
  m_outputNTuple->tree()->Branch("isMC",               	        &m_isMC, "isMC/B");

  m_outputNTuple->tree()->Branch("weight_event", 	    	&m_weight_event, "weight_event/F");
  m_outputNTuple->tree()->Branch("weight_event_trig", 	    	&m_weight_event_trig, "weight_event_trig/F");
  m_outputNTuple->tree()->Branch("weight_event_lep", 	    	&m_weight_event_lep, "weight_event_lep/F");
  m_outputNTuple->tree()->Branch("weight_tag",		    	&m_weight_tag,   "weight_tag/F");
  m_outputNTuple->tree()->Branch("weight_probe", 	    	&m_weight_probe, "weight_probe/F");

  m_outputNTuple->tree()->Branch("isSS01",               	&m_isSS01, "isSS01/B");
  m_outputNTuple->tree()->Branch("isSS12",               	&m_isSS12, "isSS12/B");
  m_outputNTuple->tree()->Branch("is_T_T",               	&m_is_T_T, "is_T_T/B");
  m_outputNTuple->tree()->Branch("is_T_AntiT",               	&m_is_T_AntiT, "is_T_AntiT/B");
  m_outputNTuple->tree()->Branch("is_AntiT_T",               	&m_is_AntiT_T, "is_AntiT_T/B");
  m_outputNTuple->tree()->Branch("is_AntiT_AntiT",              &m_is_AntiT_AntiT, "is_AntiT_AntiT/B");
  m_outputNTuple->tree()->Branch("is_Tel_AntiTmu",      	&m_is_Tel_AntiTmu, "is_Tel_AntiTmu/B");
  m_outputNTuple->tree()->Branch("is_AntiTel_Tmu",      	&m_is_AntiTel_Tmu, "is_AntiTel_Tmu/B");
  m_outputNTuple->tree()->Branch("is_Tmu_AntiTel",      	&m_is_Tmu_AntiTel, "is_Tmu_AntiTel/B");
  m_outputNTuple->tree()->Branch("is_AntiTmu_Tel",      	&m_is_AntiTmu_Tel, "is_AntiTmu_Tel/B");

  m_outputNTuple->tree()->Branch("nmuons",               	&m_nmuons, "nmuons/I");
  m_outputNTuple->tree()->Branch("nelectrons",               	&m_nelectrons, "nelectrons/I");
  m_outputNTuple->tree()->Branch("nleptons",               	&m_nleptons, "nleptons/I");

  m_outputNTuple->tree()->Branch("electron_Pt_0", 	        &m_el_Pt_0, "electron_Pt_0/F");
  m_outputNTuple->tree()->Branch("electron_Pt_1",		&m_el_Pt_1, "electron_Pt_1/F");
  m_outputNTuple->tree()->Branch("muon_Pt_0",		        &m_mu_Pt_0, "muon_Pt_0/F");
  m_outputNTuple->tree()->Branch("muon_Pt_1",  	                &m_mu_Pt_1, "muon_Pt_1/F");

  m_outputNTuple->tree()->Branch("lep_isTightSelected_0",       &m_lep_isTightSelected_0, "lep_isTightSelected_0/B");
  m_outputNTuple->tree()->Branch("lep_isTightSelected_1",       &m_lep_isTightSelected_1, "lep_isTightSelected_1/B");
  m_outputNTuple->tree()->Branch("lep_isTightSelected_2",       &m_lep_isTightSelected_2, "lep_isTightSelected_2/B");

  // From v23 onwards, no need to do this anymore --> already being done in GFW
  //
  //m_outputNTuple->tree()->Branch("lep_isTrigMatch_SLT_0",       &m_lep_isTrigMatch_SLT_0, "lep_isTrigMatch_SLT_0/B");
  //m_outputNTuple->tree()->Branch("lep_isTrigMatch_SLT_1",       &m_lep_isTrigMatch_SLT_1, "lep_isTrigMatch_SLT_1/B");
  //m_outputNTuple->tree()->Branch("lep_isTrigMatch_DLT_0",       &m_lep_isTrigMatch_DLT_0, "lep_isTrigMatch_DLT_0/B");
  //m_outputNTuple->tree()->Branch("lep_isTrigMatch_DLT_1",       &m_lep_isTrigMatch_DLT_1, "lep_isTrigMatch_DLT_1/B");

  m_outputNTuple->tree()->Branch("event_isTrigMatch_DLT",       &m_event_isTrigMatch_DLT, "event_isTrigMatch_DLT/B");

  std::vector<std::string> TPS   = { "Tag", "Probe" };
  std::vector<std::string> TRIGS = { "SLT", "DLT" };
  std::vector<std::string> VARS  = { "Pt/F", "Eta/F", "EtaBE2/F", "ptVarcone20/F", "ptVarcone30/F", "topoEtcone20/F", "sigd0PV/F", "Z0SinTheta/F", "ID/F", "deltaRClosestBJet/F", "massClosestBJet/F", "isTrigMatch/B", "isTightSelected/B", "isPrompt/B", "isBrems/B", "isFakeLep/B", "isQMisID/B", "isConvPh/B", "truthType/I", "truthOrigin/I" };

  m_outputNTuple->tree()->Branch("event_isBadTP_SLT", &m_isBadTPEvent_SLT, "event_isBadTP_SLT/B");
  m_outputNTuple->tree()->Branch("event_isBadTP_DLT", &m_isBadTPEvent_DLT, "event_isBadTP_DLT/B");

  for ( const auto& tp : TPS ) {

    for ( const auto& trig : TRIGS ) {

      for ( const auto& var : VARS ) {

        std::string branchtype     = var.substr( var.length() - 1 );
        std::string branchname     = "lep_" + tp + "_" + trig + "_" + var.substr( 0, var.length() - 2 );
        std::string branchnametype = var;

        if ( branchtype.compare("F") == 0 ) {
          m_outputNTuple->tree()->Branch( branchname.c_str(), &( m_TagProbe_branches[branchname].f ), branchnametype.c_str() );
        } else if ( branchtype.compare("B") == 0 ) {
          m_outputNTuple->tree()->Branch( branchname.c_str(), &( m_TagProbe_branches[branchname].c ), branchnametype.c_str() );
        } else if ( branchtype.compare("I") == 0 ) {
          m_outputNTuple->tree()->Branch( branchname.c_str(), &( m_TagProbe_branches[branchname].i ), branchnametype.c_str() );
        }

      }

    }

  }

  m_outputNTuple->tree()->Branch("lep_Pt",		   	&m_lep_Pt);
  m_outputNTuple->tree()->Branch("lep_Eta",		   	&m_lep_Eta);
  m_outputNTuple->tree()->Branch("lep_EtaBE2",  	   	&m_lep_EtaBE2);

  m_outputNTuple->tree()->Branch("lep_TagVec_SLT_Pt",		&m_lep_TagVec_SLT_Pt);
  m_outputNTuple->tree()->Branch("lep_ProbeVec_SLT_Pt",		&m_lep_ProbeVec_SLT_Pt);
  m_outputNTuple->tree()->Branch("lep_TagVec_DLT_Pt",  	   	&m_lep_TagVec_DLT_Pt);
  m_outputNTuple->tree()->Branch("lep_ProbeVec_DLT_Pt",		&m_lep_ProbeVec_DLT_Pt);

  m_outputNTuple->tree()->Branch("electron_TagVec_SLT_Pt",	&m_el_TagVec_SLT_Pt);
  m_outputNTuple->tree()->Branch("electron_ProbeVec_SLT_Pt",	&m_el_ProbeVec_SLT_Pt);
  m_outputNTuple->tree()->Branch("electron_TagVec_DLT_Pt",  	&m_el_TagVec_DLT_Pt);
  m_outputNTuple->tree()->Branch("electron_ProbeVec_DLT_Pt",	&m_el_ProbeVec_DLT_Pt);

  m_outputNTuple->tree()->Branch("muon_TagVec_SLT_Pt",		&m_mu_TagVec_SLT_Pt);
  m_outputNTuple->tree()->Branch("muon_ProbeVec_SLT_Pt",	&m_mu_ProbeVec_SLT_Pt);
  m_outputNTuple->tree()->Branch("muon_TagVec_DLT_Pt",  	&m_mu_TagVec_DLT_Pt);
  m_outputNTuple->tree()->Branch("muon_ProbeVec_DLT_Pt",	&m_mu_ProbeVec_DLT_Pt);

  m_outputNTuple->tree()->Branch("jet_OR_Pt",		   	        &m_jet_OR_Pt);
  m_outputNTuple->tree()->Branch("jet_OR_Eta",		   	        &m_jet_OR_Eta);
  m_outputNTuple->tree()->Branch("jet_OR_Phi",  	   	        &m_jet_OR_Phi);
  m_outputNTuple->tree()->Branch("jet_OR_E",		   	        &m_jet_OR_E);
  m_outputNTuple->tree()->Branch("jet_OR_truthMatch_Pt",		&m_jet_OR_truthMatch_Pt);
  m_outputNTuple->tree()->Branch("jet_OR_truthMatch_Eta",  	   	&m_jet_OR_truthMatch_Eta);
  m_outputNTuple->tree()->Branch("jet_OR_truthMatch_Phi",		&m_jet_OR_truthMatch_Phi);
  m_outputNTuple->tree()->Branch("jet_OR_truthMatch_E",		        &m_jet_OR_truthMatch_E);
  m_outputNTuple->tree()->Branch("jet_OR_truthMatch_isBJet",  	        &m_jet_OR_truthMatch_isBJet);
  m_outputNTuple->tree()->Branch("jet_OR_truthMatch_isCJet",		&m_jet_OR_truthMatch_isCJet);
  m_outputNTuple->tree()->Branch("jet_OR_truthMatch_isLFJet",	        &m_jet_OR_truthMatch_isLFJet);
  m_outputNTuple->tree()->Branch("jet_OR_truthMatch_isGluonJet",        &m_jet_OR_truthMatch_isGluonJet);

  // ---------------------------------------------------------------------------------------------------------------

  // Initialise counter for input TTree entries processed
  //
  m_numEntry = 0;

  m_effectiveTotEntries = m_inputNTuple->GetEntries();

  unsigned int maxEvents = static_cast<int>( wk()->metaData()->castDouble("nc_EventLoop_MaxEvents") );

  if ( maxEvents > 0 ) {
    m_effectiveTotEntries = maxEvents;
  }

  Info("initialize()", "Name of input TTree : %s", m_inputNTuple->GetName() );
  Info("initialize()", "Total events to run on: %u", m_effectiveTotEntries );

  // ---------------------------------------------------------------------------------------------------------------

  // Print whether we are using an EL::AlgSelect algorithm in our job
  //
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

  if ( m_numEntry < m_sumWeightsTree->GetEntries() ) {

    m_sumWeightsTree->GetEntry( m_numEntry );

    m_sumGenEvents         += m_totalEvents;
    m_sumGenEventsWeighted += m_totalEventsWeighted;

  }

  m_inputNTuple->GetEntry( wk()->treeEntry() );

  if ( m_debug ) {
    std::cout << "" << std::endl;
    Info("execute()", "===> Entry %u - EventNumber = %u - RunYear = %i", static_cast<uint32_t>(m_numEntry), static_cast<uint32_t>(m_EventNumber), m_RunYear );
  }

  ++m_numEntry;

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

  // ------------------------------------------------------------------------

  auto lep0 = std::make_shared<leptonObj>();

  lep0.get()->pt          = m_lep_Pt_0;
  lep0.get()->eta         = m_lep_Eta_0;
  lep0.get()->etaBE2      = m_lep_EtaBE2_0;
  lep0.get()->phi         = m_lep_Phi_0;
  lep0.get()->ID          = m_lep_ID_0;
  lep0.get()->flavour     = fabs(m_lep_ID_0);
  lep0.get()->charge      = m_lep_ID_0 / fabs(m_lep_ID_0);
  lep0.get()->d0sig       = m_lep_sigd0PV_0;
  lep0.get()->z0sintheta  = m_lep_Z0SinTheta_0;
  lep0.get()->pid         = m_lep_isTightLH_0;
  lep0.get()->isolated       = ( fabs(m_lep_ID_0) == 13 ) ?  m_lep_isolationFixedCutTightTrackOnly_0 : m_lep_isolationFixedCutTight_0;
  lep0.get()->ptVarcone20    = m_lep_ptVarcone20_0;
  lep0.get()->ptVarcone30    = m_lep_ptVarcone30_0;
  lep0.get()->topoEtcone20   = m_lep_topoEtcone20_0;
  lep0.get()->trackisooverpt = ( fabs(m_lep_ID_0) == 13 ) ?  m_lep_ptVarcone30_0/m_lep_Pt_0 : m_lep_ptVarcone20_0/m_lep_Pt_0;
  lep0.get()->caloisooverpt  = ( fabs(m_lep_ID_0) == 13 ) ?  -1.0 : m_lep_topoEtcone20_0/m_lep_Pt_0;
  lep0.get()->trigmatched     = m_lep_isTrigMatch_0;
  lep0.get()->trigmatched_DLT = m_lep_isTrigMatchDLT_0;
  lep0.get()->prompt      = m_lep_isPrompt_0;
  lep0.get()->fake        = m_lep_isFakeLep_0;
  lep0.get()->brems       = m_lep_isBrems_0;
  lep0.get()->qmisid      = m_lep_isQMisID_0;
  lep0.get()->convph      = m_lep_isConvPh_0;
  lep0.get()->truthType   = m_lep_truthType_0;
  lep0.get()->truthOrigin = m_lep_truthOrigin_0;
  ANA_CHECK( this->checkIsTightLep( lep0 ) );

  lep0.get()->SFIDLoose    = m_lep_SFIDLoose_0;
  lep0.get()->SFIDTight    = m_lep_SFIDTight_0;
  lep0.get()->SFTrigLoose  = m_lep_SFTrigLoose_0;
  lep0.get()->SFTrigTight  = m_lep_SFTrigTight_0;
  lep0.get()->EffTrigLoose = m_lep_EffTrigLoose_0;
  lep0.get()->EffTrigTight = m_lep_EffTrigTight_0;
  lep0.get()->SFIsoLoose   = m_lep_SFIsoLoose_0;
  lep0.get()->SFIsoTight   = m_lep_SFIsoTight_0;
  lep0.get()->SFReco       = m_lep_SFReco_0;
  lep0.get()->SFTTVA       = m_lep_SFTTVA_0;
  lep0.get()->SFObjLoose   = m_lep_SFObjLoose_0;
  lep0.get()->SFObjTight   = m_lep_SFObjTight_0;

  m_leptons.push_back(lep0);

  auto lep1 = std::make_shared<leptonObj>();

  lep1.get()->pt          = m_lep_Pt_1;
  lep1.get()->eta         = m_lep_Eta_1;
  lep1.get()->etaBE2      = m_lep_EtaBE2_1;
  lep1.get()->phi         = m_lep_Phi_1;
  lep1.get()->ID          = m_lep_ID_1;
  lep1.get()->flavour     = fabs(m_lep_ID_1);
  lep1.get()->charge      = m_lep_ID_1 / fabs(m_lep_ID_1);
  lep1.get()->d0sig       = m_lep_sigd0PV_1;
  lep1.get()->z0sintheta  = m_lep_Z0SinTheta_1;
  lep1.get()->pid         = m_lep_isTightLH_1;
  lep1.get()->isolated       = ( fabs(m_lep_ID_1) == 13 ) ?  m_lep_isolationFixedCutTightTrackOnly_1 : m_lep_isolationFixedCutTight_1;
  lep1.get()->ptVarcone20    = m_lep_ptVarcone20_1;
  lep1.get()->ptVarcone30    = m_lep_ptVarcone30_1;
  lep1.get()->topoEtcone20   = m_lep_topoEtcone20_1;
  lep1.get()->trackisooverpt = ( fabs(m_lep_ID_1) == 13 ) ?  m_lep_ptVarcone30_1/m_lep_Pt_1 : m_lep_ptVarcone20_1/m_lep_Pt_1;
  lep1.get()->caloisooverpt  = ( fabs(m_lep_ID_1) == 13 ) ?  -1.0 : m_lep_topoEtcone20_1/m_lep_Pt_1;
  lep1.get()->trigmatched     = m_lep_isTrigMatch_1;
  lep1.get()->trigmatched_DLT = m_lep_isTrigMatchDLT_1;
  lep1.get()->prompt      = m_lep_isPrompt_1;
  lep1.get()->fake        = m_lep_isFakeLep_1;
  lep1.get()->brems       = m_lep_isBrems_1;
  lep1.get()->qmisid      = m_lep_isQMisID_1;
  lep1.get()->convph      = m_lep_isConvPh_1;
  lep1.get()->truthType   = m_lep_truthType_1;
  lep1.get()->truthOrigin = m_lep_truthOrigin_1;
  ANA_CHECK( this->checkIsTightLep( lep1 ) );

  lep1.get()->SFIDLoose    = m_lep_SFIDLoose_1;
  lep1.get()->SFIDTight    = m_lep_SFIDTight_1;
  lep1.get()->SFTrigLoose  = m_lep_SFTrigLoose_1;
  lep1.get()->SFTrigTight  = m_lep_SFTrigTight_1;
  lep1.get()->EffTrigLoose = m_lep_EffTrigLoose_1;
  lep1.get()->EffTrigTight = m_lep_EffTrigTight_1;
  lep1.get()->SFIsoLoose   = m_lep_SFIsoLoose_1;
  lep1.get()->SFIsoTight   = m_lep_SFIsoTight_1;
  lep1.get()->SFReco       = m_lep_SFReco_1;
  lep1.get()->SFTTVA       = m_lep_SFTTVA_1;
  lep1.get()->SFObjLoose   = m_lep_SFObjLoose_1;
  lep1.get()->SFObjTight   = m_lep_SFObjTight_1;

  m_leptons.push_back(lep1);

  if ( m_debug ) {
    Info("execute()","lep0 pT = %.2f", lep0.get()->pt/1e3 );
    Info("execute()","lep1 pT = %.2f", lep1.get()->pt/1e3 );
    Info("execute()","lep0 SFTrigTight = %.2f - SFTrigLoose = %.2f ", lep0.get()->SFTrigTight, lep0.get()->SFTrigLoose );
    Info("execute()","lep1 SFTrigTight = %.2f - SFTrigLoose = %.2f ", lep1.get()->SFTrigTight, lep1.get()->SFTrigLoose );
    Info("execute()","lep0 EffTrigTight = %.2f - EffTrigLoose = %.2f ", lep0.get()->EffTrigTight, lep0.get()->EffTrigLoose );
    Info("execute()","lep1 EffTrigTight = %.2f - EffTrigLoose = %.2f ", lep1.get()->EffTrigTight, lep1.get()->EffTrigLoose );
    Info("execute()","lep0 SFObjTight = %.2f - SFObjLoose = %.2f", lep0.get()->SFObjTight, lep0.get()->SFObjLoose );
    Info("execute()","lep1 SFObjTight = %.2f - SFObjLoose = %.2f", lep1.get()->SFObjTight, lep1.get()->SFObjLoose );
  }

  if ( m_trilep_type ) {

    auto lep2 = std::make_shared<leptonObj>();

    lep2.get()->pt	    = m_lep_Pt_2;
    lep2.get()->eta	    = m_lep_Eta_2;
    lep2.get()->etaBE2      = m_lep_EtaBE2_2;
    lep2.get()->phi         = m_lep_Phi_2;
    lep2.get()->ID          = m_lep_ID_2;
    lep2.get()->flavour     = fabs(m_lep_ID_2);
    lep2.get()->charge      = m_lep_ID_2 / fabs(m_lep_ID_2);
    lep2.get()->d0sig       = m_lep_sigd0PV_2;
    lep2.get()->z0sintheta  = m_lep_Z0SinTheta_2;
    lep2.get()->pid	    = m_lep_isTightLH_2;
    lep2.get()->isolated       = ( fabs(m_lep_ID_2) == 13 ) ?  m_lep_isolationFixedCutTightTrackOnly_2 : m_lep_isolationFixedCutTight_2;
    lep2.get()->ptVarcone20    = m_lep_ptVarcone20_2;
    lep2.get()->ptVarcone30    = m_lep_ptVarcone30_2;
    lep2.get()->topoEtcone20   = m_lep_topoEtcone20_2;
    lep2.get()->trackisooverpt = ( fabs(m_lep_ID_2) == 13 ) ?  m_lep_ptVarcone30_2/m_lep_Pt_2 : m_lep_ptVarcone20_2/m_lep_Pt_2;
    lep2.get()->caloisooverpt  = ( fabs(m_lep_ID_2) == 13 ) ?  -1.0 : m_lep_topoEtcone20_2/m_lep_Pt_2;
    lep2.get()->trigmatched     = m_lep_isTrigMatch_2;
    lep2.get()->trigmatched_DLT = m_lep_isTrigMatchDLT_2;
    lep2.get()->prompt      = m_lep_isPrompt_2;
    lep2.get()->fake	    = m_lep_isFakeLep_2;
    lep2.get()->brems       = m_lep_isBrems_2;
    lep2.get()->qmisid      = m_lep_isQMisID_2;
    lep2.get()->convph      = m_lep_isConvPh_2;
    lep2.get()->truthType   = m_lep_truthType_2;
    lep2.get()->truthOrigin = m_lep_truthOrigin_2;
    ANA_CHECK( this->checkIsTightLep( lep2 ) );

    lep2.get()->SFIDLoose    = m_lep_SFIDLoose_2;
    lep2.get()->SFIDTight    = m_lep_SFIDTight_2;
    lep2.get()->SFTrigLoose  = m_lep_SFTrigLoose_2;
    lep2.get()->SFTrigTight  = m_lep_SFTrigTight_2;
    lep2.get()->EffTrigLoose = m_lep_EffTrigLoose_2;
    lep2.get()->EffTrigTight = m_lep_EffTrigTight_2;
    lep2.get()->SFIsoLoose   = m_lep_SFIsoLoose_2;
    lep2.get()->SFIsoTight   = m_lep_SFIsoTight_2;
    lep2.get()->SFReco       = m_lep_SFReco_2;
    lep2.get()->SFTTVA       = m_lep_SFTTVA_2;
    lep2.get()->SFObjLoose   = m_lep_SFObjLoose_2;
    lep2.get()->SFObjTight   = m_lep_SFObjTight_2;

    if ( m_debug ) { Info("execute()","lep2 pT = %.2f", lep2.get()->pt/1e3 ); }

    m_leptons.push_back(lep2);
  }

  // ------------------------------------------------------------------------

  ANA_CHECK( this->decorateEvent() );
  // ------------------------------------------------------------------------

  ANA_CHECK( this->defineTagAndProbe() );
  //m_debug = false;

  // ------------------------------------------------------------------------

  ANA_CHECK( this->decorateWeights() );

  // ------------------------------------------------------------------------

  ANA_CHECK( this->jetTruthMatching() );

  // ------------------------------------------------------------------------

  ANA_CHECK( this->setOutputBranches() );

  // ------------------------------------------------------------------------

  // Don't forget to clear object containers before going to next event!

  m_leptons.clear();
  m_bjets.clear();

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
  //
  m_inputNTuple->SetBranchStatus ("*", 0);

  std::vector<std::string> branch_vec;

  // Parse input list, split by comma, and put into a vector
  //
  std::string token;
  std::istringstream ss( m_inputBranches );
  while ( std::getline(ss, token, ',') ) { branch_vec.push_back(token); }

  // Re-enable the branches we are going to use
  //
  for ( const auto& branch : branch_vec ) {

    size_t found_regexp = branch.find(".*");
    if (  found_regexp != std::string::npos ) {
      std::string wildcarded_branch(branch);
      wildcarded_branch =  wildcarded_branch.replace( found_regexp, 2, "*" ); // replace ".*" with "*"
      m_inputNTuple->SetBranchStatus (wildcarded_branch.c_str(), 1);
    } else {
      m_inputNTuple->SetBranchStatus (branch.c_str(), 1);
    }
  }

  return EL::StatusCode::SUCCESS;

}

EL::StatusCode HTopMultilepMiniNTupMaker ::  checkIsTightLep( std::shared_ptr<leptonObj> lep )
{

  bool isTight(false);

  if ( m_debug ) { Info("checkIsTightLep()", "Lepton flavour = %i", lep.get()->flavour ); }

  switch ( lep.get()->flavour )
  {
     case 11:
       if ( lep.get()->isolated && lep.get()->pid && fabs(lep.get()->d0sig) < 5.0 && fabs(lep.get()->z0sintheta) < 0.5 ) { isTight = true; }
       break;
     case 13:
       if ( lep.get()->isolated && fabs(lep.get()->d0sig) < 3.0 && fabs(lep.get()->z0sintheta) < 0.5 ) { isTight = true; }
       break;
     default:
       break;
  }

  lep.get()->tight = isTight;

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode HTopMultilepMiniNTupMaker :: decorateEvent ( )
{

  m_event.get()->isMC        = ( m_mc_channel_number > 0 );
  m_event.get()->dilep       = ( m_leptons.size() == 2 );
  m_event.get()->trilep      = ( m_leptons.size() == 3 );
  m_event.get()->dilep_type  = ( m_dilep_type );

  if ( m_event.get()->dilep ) {

     int prod_lep_charge(1);
     for ( const auto& lep : m_leptons ) { prod_lep_charge *= lep.get()->charge; }

     m_event.get()->isSS01 = ( prod_lep_charge > 0 );

     // Fill flat branches with trigger matching decision split in SLT/DLT, and per-event (pair) matching for DLT

     // From v23 onwards, no need to do this anymore --> already being done in GFW
     //ANA_CHECK( this->triggerMatching() );

     m_event_isTrigMatch_DLT = ( m_lep_isTrigMatchDLT_0 && m_lep_isTrigMatchDLT_1 );

  } else if ( m_event.get()->trilep ) {

     m_event.get()->isSS12 = ( fabs( m_leptons.at(0).get()->charge + m_leptons.at(1).get()->charge + m_leptons.at(2).get()->charge ) != 3 );

  }

  // Store b-tagged jets (for dilepton events, use jets after OLR w/ taus)

  unsigned int njets = ( m_event.get()->dilep ) ? m_selected_jets_T->size() : m_selected_jets->size();

  for ( unsigned int j_idx(0); j_idx < njets; ++j_idx ) {

      short j = ( m_event.get()->dilep ) ? m_selected_jets_T->at(j_idx) : m_selected_jets->at(j_idx);

      if ( m_jet_flavor_weight_MV2c10->at(j) < 0.8244273 ) { continue; } // MV2c10, FixedCutBEff_70 : https://twiki.cern.ch/twiki/bin/view/AtlasProtected/BTaggingBenchmarks#MV2c20_tagger

      auto bjet = std::make_shared<bjetObj>();

      bjet.get()->pt  = m_jet_pt->at(j);
      bjet.get()->eta = m_jet_eta->at(j);
      bjet.get()->phi = m_jet_phi->at(j);

      m_bjets.push_back(bjet);

  }

  m_event.get()->nbjets = ( m_event.get()->dilep ) ? m_nJets_OR_T_MV2c10_70 : m_nJets_OR_MV2c10_70;

  if ( m_debug ) { Info("decorateEvent()","Number of bjets (MV2c10, 70 WP) = %i", m_event.get()->nbjets); }

  // Find closest b-tagged jet to each lepton

  ANA_CHECK( this->findClosestBJetLep() );

  return EL::StatusCode::SUCCESS;

}

EL::StatusCode HTopMultilepMiniNTupMaker :: findClosestBJetLep() {

    if ( ! m_event.get()->dilep && ! m_event.get()->trilep ) { return EL::StatusCode::SUCCESS;  }

    unsigned int idx_lep(0);

    int lep_flavour(0);
    float lep_pt(-1.0), lep_eta(-999.0), lep_phi(-999.0), lep_m(-1.0);

    TLorentzVector lepTLV, bjetTLV, closestbjetTLV, lepbjetTLV;

    for ( auto lep : m_leptons ) {

	lep_flavour = lep.get()->flavour;
	lep_pt      = lep.get()->pt;
	lep_eta     = ( lep_flavour == 11 ) ? lep.get()->etaBE2 : lep.get()->eta;
	lep_phi     = lep.get()->phi;
	lep_m       = ( lep_flavour == 11 ) ? 0.511 : 105.65;

	lepTLV.SetPtEtaPhiM( lep_pt, lep_eta, lep_phi, lep_m );

	if ( m_debug ) { Info("findClosestBJetLep()","Checking lepton[%i] w/ flavour = %i, pT = %.2f [GeV], eta = %.2f, phi = %.2f, m = %.3f [GeV]", idx_lep, lep_flavour, lep_pt/1e3, lep_eta, lep_phi, lep_m/1e3 ); }

	// Now find closest bjet to this lepton

	unsigned int idx_bj(0);

	int idx_closest_bj(-1);

	float bjet_pt(-.0),  bjet_eta(-999.0), bjet_phi(-999.0);

        float dist_to_closest_bj(9e9), this_dist(-1.0);
	float mass_lep_closest_bj(-1.0);

	for ( auto bjet : m_bjets ) {

	    bjet_pt  = bjet.get()->pt;
	    bjet_eta = bjet.get()->eta;
	    bjet_phi = bjet.get()->phi;

	    bjetTLV.SetPtEtaPhiM( bjet_pt, bjet_eta, bjet_phi, 418.0 ); // use mass of b quark [MeV]

	    this_dist = lepTLV.DeltaR(bjetTLV);

	    if ( m_debug ) { Info("findClosestBJetLep()","\t bjet[%i] w/ pT = %.2f [GeV], eta = %.2f, phi = %.2f, DeltaR(lep[%i],bjet[%i]) = %.3f", idx_bj, bjet_pt/1e3, bjet_eta, bjet_phi, idx_lep, idx_bj, this_dist ); }

	    if ( this_dist < dist_to_closest_bj ) {
		dist_to_closest_bj = this_dist;
		idx_closest_bj     = idx_bj;
	    }

	    ++idx_bj;

	}

	// Get invariant mass of lepton and closest bjet

	if ( m_event.get()->nbjets > 0 && idx_closest_bj > -1 ) {

	    closestbjetTLV.SetPtEtaPhiM( m_bjets.at(idx_closest_bj).get()->pt, m_bjets.at(idx_closest_bj).get()->eta, m_bjets.at(idx_closest_bj).get()->phi, 418.0  );

	    lepbjetTLV = lepTLV + closestbjetTLV;

	    mass_lep_closest_bj = lepbjetTLV.M();

	}

	// Decorate lepton (set a dummy negative decorator if there are no bjets in event)

	lep.get()->massClosestBJet   = mass_lep_closest_bj;
	lep.get()->deltaRClosestBJet = ( m_event.get()->nbjets > 0 ) ? dist_to_closest_bj : -1.0;

	if ( m_debug ) { Info("findClosestBJetLep()","==> DeltaR(lep[%i],closest bjet[%i]) = %.3f, M(lep[%i],closest bjet[%i]) = %.3f [GeV]", idx_lep, idx_closest_bj, dist_to_closest_bj, idx_lep, idx_closest_bj, mass_lep_closest_bj/1e3 ); }

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


EL::StatusCode HTopMultilepMiniNTupMaker :: triggerMatching()
{

  if ( m_debug ) { Info("triggerMatching()","Performing trigger matching..." ); }

  auto lep0 = m_leptons.at(0);
  auto lep1 = m_leptons.at(1);

  // Initialise w/ dummy values

  m_lep_isTrigMatch_SLT_0 = m_lep_isTrigMatch_SLT_1 = m_lep_isTrigMatch_DLT_0 = m_lep_isTrigMatch_DLT_1 = -1;

  // Get the indexes of leading/subleading leptons which passed the OLR

  // NB: For the SLT case, we always consider the OR of the SLT chains. Hence the pT threshold cut should be just the one for the lowest pT chain

  int el0_idx(-1), el1_idx(-1), mu0_idx(-1), mu1_idx(-1);

  if ( m_dilep_type == 1 ) { // 1) mumu

    ANA_CHECK( this->getPostOLRIndex( mu0_idx, 0, "muon" ) );
    ANA_CHECK( this->getPostOLRIndex( mu1_idx, 1, "muon" ) );

    if ( m_RunYear == 2015 ) {

      m_lep_isTrigMatch_SLT_0 = ( lep0.get()->pt > 1.05*20e3 && ( m_muon_match_HLT_mu20_iloose_L1MU15->at(mu0_idx) || m_muon_match_HLT_mu50->at(mu0_idx) ) );
      m_lep_isTrigMatch_SLT_1 = ( lep1.get()->pt > 1.05*20e3 && ( m_muon_match_HLT_mu20_iloose_L1MU15->at(mu1_idx) || m_muon_match_HLT_mu50->at(mu1_idx) ) );
      m_lep_isTrigMatch_DLT_0 = ( m_muon_match_HLT_mu18_mu8noL1->at(mu0_idx) && lep0.get()->pt > 1.05*18e3 );
      m_lep_isTrigMatch_DLT_1 = ( m_muon_match_HLT_mu18_mu8noL1->at(mu1_idx) && lep1.get()->pt > 1.05*8e3 );

    } else if ( m_RunYear == 2016 ) {


      m_lep_isTrigMatch_SLT_0 = ( lep0.get()->pt > 1.05*26e3 && ( m_muon_match_HLT_mu26_ivarmedium->at(mu0_idx) || m_muon_match_HLT_mu50->at(mu0_idx) ) );
      m_lep_isTrigMatch_SLT_1 = ( lep1.get()->pt > 1.05*26e3 && ( m_muon_match_HLT_mu26_ivarmedium->at(mu1_idx) || m_muon_match_HLT_mu50->at(mu1_idx) ) );
      m_lep_isTrigMatch_DLT_0 = ( m_muon_match_HLT_mu22_mu8noL1->at(mu0_idx) && lep0.get()->pt > 1.05*22e3 );
      m_lep_isTrigMatch_DLT_1 = ( m_muon_match_HLT_mu22_mu8noL1->at(mu1_idx) && lep1.get()->pt > 1.05*8e3 );

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

      m_lep_isTrigMatch_SLT_0 = ( ( lep0.get()->flavour == 11 && lep0.get()->pt > 25e3 && ( m_electron_match_HLT_e24_lhmedium_L1EM20VH->at(el0_idx) || m_electron_match_HLT_e60_lhmedium->at(el0_idx) || m_electron_match_HLT_e120_lhloose->at(el0_idx) ) ) || ( lep0.get()->flavour == 13 && lep0.get()->pt > 1.05*20e3 && ( m_muon_match_HLT_mu20_iloose_L1MU15->at(mu0_idx) || m_muon_match_HLT_mu50->at(mu0_idx) ) ) );
      m_lep_isTrigMatch_SLT_1 = ( ( lep1.get()->flavour == 11 && lep1.get()->pt > 25e3 && ( m_electron_match_HLT_e24_lhmedium_L1EM20VH->at(el0_idx) || m_electron_match_HLT_e60_lhmedium->at(el0_idx) || m_electron_match_HLT_e120_lhloose->at(el0_idx) ) ) || ( lep1.get()->flavour == 13 && lep1.get()->pt > 1.05*20e3 && ( m_muon_match_HLT_mu20_iloose_L1MU15->at(mu0_idx) || m_muon_match_HLT_mu50->at(mu0_idx) ) ) );
      m_lep_isTrigMatch_DLT_0 = ( ( lep0.get()->flavour == 11 && ( ( m_electron_match_HLT_e7_medium_mu24->at(el0_idx) && lep0.get()->pt > 8e3 ) || ( m_electron_match_HLT_e24_medium_L1EM20VHI_mu8noL1->at(el0_idx) && lep0.get()->pt > 25e3 ) ) ) || ( lep0.get()->flavour == 13 && ( ( m_muon_match_HLT_e7_medium_mu24->at(mu0_idx) && lep0.get()->pt > 1.05*24e3 ) || ( m_muon_match_HLT_e24_medium_L1EM20VHI_mu8noL1->at(mu0_idx) && lep0.get()->pt > 1.05*8e3 ) ) ) );
      m_lep_isTrigMatch_DLT_1 = ( ( lep1.get()->flavour == 11 && ( ( m_electron_match_HLT_e7_medium_mu24->at(el0_idx) && lep1.get()->pt > 8e3 ) || ( m_electron_match_HLT_e24_medium_L1EM20VHI_mu8noL1->at(el0_idx) && lep1.get()->pt > 25e3 ) ) ) || ( lep1.get()->flavour == 13 && ( ( m_muon_match_HLT_e7_medium_mu24->at(mu0_idx) && lep1.get()->pt > 1.05*24e3 ) || ( m_muon_match_HLT_e24_medium_L1EM20VHI_mu8noL1->at(mu0_idx) && lep1.get()->pt > 1.05*8e3 ) ) ) );

    } else if ( m_RunYear == 2016 ) {

      m_lep_isTrigMatch_SLT_0 = ( ( lep0.get()->flavour == 11 && lep0.get()->pt > 27e3 && ( m_electron_match_HLT_e26_lhtight_nod0_ivarloose->at(el0_idx) || m_electron_match_HLT_e60_lhmedium_nod0->at(el0_idx) || m_electron_match_HLT_e140_lhloose_nod0->at(el0_idx) ) ) || ( lep0.get()->flavour == 13 && lep0.get()->pt > 1.05*26e3 && ( m_muon_match_HLT_mu26_ivarmedium->at(mu0_idx) || m_muon_match_HLT_mu50->at(mu0_idx) ) ) );
      m_lep_isTrigMatch_SLT_1 = ( ( lep1.get()->flavour == 11 && lep1.get()->pt > 27e3 && ( m_electron_match_HLT_e26_lhtight_nod0_ivarloose->at(el0_idx) || m_electron_match_HLT_e60_lhmedium_nod0->at(el0_idx) || m_electron_match_HLT_e140_lhloose_nod0->at(el0_idx) ) ) || ( lep1.get()->flavour == 13 && lep1.get()->pt > 1.05*26e3 && ( m_muon_match_HLT_mu26_ivarmedium->at(mu0_idx) || m_muon_match_HLT_mu50->at(mu0_idx) ) ) );
      m_lep_isTrigMatch_DLT_0 = ( ( lep0.get()->flavour == 11 && ( ( m_electron_match_HLT_e17_lhloose_mu14->at(el0_idx) && lep0.get()->pt > 18e3 ) || ( m_electron_match_HLT_e17_lhloose_nod0_mu14->at(el0_idx) && lep0.get()->pt > 18e3 ) ) ) || ( lep0.get()->flavour == 13 && ( ( m_muon_match_HLT_e17_lhloose_mu14->at(mu0_idx) && lep0.get()->pt > 1.05*14e3 ) || ( m_muon_match_HLT_e17_lhloose_nod0_mu14->at(mu0_idx) && lep0.get()->pt > 1.05*14e3 ) ) ) );
      m_lep_isTrigMatch_DLT_1 = ( ( lep1.get()->flavour == 11 && ( ( m_electron_match_HLT_e17_lhloose_mu14->at(el0_idx) && lep1.get()->pt > 18e3 ) || ( m_electron_match_HLT_e17_lhloose_nod0_mu14->at(el0_idx) && lep1.get()->pt > 18e3 ) ) ) || ( lep1.get()->flavour == 13 && ( ( m_muon_match_HLT_e17_lhloose_mu14->at(mu0_idx) && lep1.get()->pt > 1.05*14e3 ) || ( m_muon_match_HLT_e17_lhloose_nod0_mu14->at(mu0_idx) && lep1.get()->pt > 1.05*14e3 ) ) ) );

    }

  } else if ( m_dilep_type == 3 ) { // 3) ee

    ANA_CHECK( this->getPostOLRIndex( el0_idx, 0, "electron" ) );
    ANA_CHECK( this->getPostOLRIndex( el1_idx, 1, "electron" ) );

    if ( m_RunYear == 2015 ) {

      m_lep_isTrigMatch_SLT_0 = ( lep0.get()->pt > 25e3 && ( m_electron_match_HLT_e24_lhmedium_L1EM20VH->at(el0_idx) || m_electron_match_HLT_e60_lhmedium->at(el0_idx) || m_electron_match_HLT_e120_lhloose->at(el0_idx) ) );
      m_lep_isTrigMatch_SLT_1 = ( lep1.get()->pt > 25e3 && ( m_electron_match_HLT_e24_lhmedium_L1EM20VH->at(el1_idx) || m_electron_match_HLT_e60_lhmedium->at(el1_idx) || m_electron_match_HLT_e120_lhloose->at(el1_idx) ) );
      m_lep_isTrigMatch_DLT_0 = ( m_electron_match_HLT_2e12_lhloose_L12EM10VH->at(el0_idx) && lep0.get()->pt > 13e3 );
      m_lep_isTrigMatch_DLT_1 = ( m_electron_match_HLT_2e12_lhloose_L12EM10VH->at(el1_idx) && lep1.get()->pt > 13e3 );

    } else if ( m_RunYear == 2016 ) {

      m_lep_isTrigMatch_SLT_0 = ( lep0.get()->pt > 27e3 && ( m_electron_match_HLT_e26_lhtight_nod0_ivarloose->at(el0_idx) || m_electron_match_HLT_e60_lhmedium_nod0->at(el0_idx) || m_electron_match_HLT_e140_lhloose_nod0->at(el0_idx) ) );
      m_lep_isTrigMatch_SLT_1 = ( lep1.get()->pt > 27e3 && ( m_electron_match_HLT_e26_lhtight_nod0_ivarloose->at(el1_idx) || m_electron_match_HLT_e60_lhmedium_nod0->at(el1_idx) || m_electron_match_HLT_e140_lhloose_nod0->at(el1_idx) ) );
      m_lep_isTrigMatch_DLT_0 = ( m_electron_match_HLT_2e17_lhvloose_nod0->at(el0_idx) && lep0.get()->pt > 18e3 );
      m_lep_isTrigMatch_DLT_1 = ( m_electron_match_HLT_2e17_lhvloose_nod0->at(el1_idx) && lep1.get()->pt > 18e3 );

    }

  }

  m_event_isTrigMatch_DLT = ( m_lep_isTrigMatch_DLT_0 && m_lep_isTrigMatch_DLT_1 );

  lep0.get()->trigmatched_SLT = m_lep_isTrigMatch_SLT_0;
  lep1.get()->trigmatched_SLT = m_lep_isTrigMatch_SLT_1;

  lep0.get()->trigmatched_DLT = m_lep_isTrigMatch_DLT_0;
  lep1.get()->trigmatched_DLT = m_lep_isTrigMatch_DLT_1;

  return EL::StatusCode::SUCCESS;

}


EL::StatusCode HTopMultilepMiniNTupMaker :: decorateWeights ()
{

  // Compute the event weights
  //
  m_event.get()->weight_event = m_mcWeightOrg * m_pileupEventWeight_090 * m_MV2c10_70_EventWeight * m_JVT_EventWeight;

  float weight_lep(1.0);
  for ( auto lep : m_leptons ) {
    if ( lep.get()->tight ) weight_lep *= lep.get()->SFObjTight;
    else                    weight_lep *= lep.get()->SFObjLoose;
  }
  m_event.get()->weight_event_lep = weight_lep;

  if ( m_event.get()->dilep ) {

    auto lep0 = m_leptons.at(0);
    auto lep1 = m_leptons.at(1);

    if ( lep0.get()->tight && lep1.get()->tight ) {

      if ( lep0.get()->trigmatched && lep1.get()->trigmatched ) { // I believe this is wrong...

        m_event.get()->weight_tag   = lep0.get()->SFObjTight * lep0.get()->SFTrigTight;
        m_event.get()->weight_probe = lep1.get()->SFObjTight;

      } else if ( lep0.get()->trigmatched && !lep1.get()->trigmatched ) {

        m_event.get()->weight_tag   = lep0.get()->SFObjTight * lep0.get()->SFTrigTight;
        m_event.get()->weight_probe = lep1.get()->SFObjTight;

      } else if ( !lep0.get()->trigmatched && lep1.get()->trigmatched ) {

        m_event.get()->weight_tag   = lep1.get()->SFObjTight * lep1.get()->SFTrigTight;
        m_event.get()->weight_probe = lep0.get()->SFObjTight;

      }

    } else if ( lep0.get()->tight && !lep1.get()->tight ) { // I believe this is wrong...what if both are trigger matched?

      m_event.get()->weight_tag   = lep0.get()->SFObjTight * lep0.get()->SFTrigTight;
      m_event.get()->weight_probe = lep1.get()->SFObjLoose;


    } else if ( !lep0.get()->tight && lep1.get()->tight ) { // I believe this is wrong...what if both are trigger matched?

      m_event.get()->weight_tag   = lep1.get()->SFObjTight * lep1.get()->SFTrigTight;
      m_event.get()->weight_probe = lep0.get()->SFObjLoose;

    }

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

  float trig_weight_N(1.0), trig_weight_D(1.0);
  float this_SF(1.0), this_eff(0.0);
  for ( auto lep : m_leptons ) {

    this_eff = ( lep.get()->tight ) ? lep.get()->EffTrigTight : lep.get()->EffTrigLoose;
    this_SF  = ( lep.get()->tight ) ? lep.get()->SFTrigTight  : lep.get()->SFTrigLoose;

    trig_weight_N *= ( 1.0 - this_SF * this_eff );
    trig_weight_D *= ( 1.0 - this_eff );

    if ( m_debug ) {
      Info("decorateWeights()", "\tthis_eff = %.2f - this_SF = %.2f", this_eff, this_SF );
      Info("decorateWeights()", "\tN block = %.2f - D block = %.2f", trig_weight_N, trig_weight_D );
    }

  }

  // Update numerator and denominator
  // Make sure the SF in the 0/0 case (i.e, when efficiency=0) will be set equal to 1
  //
  trig_weight_N = ( trig_weight_N != 1.0 ) ? trig_weight_N : 0.0;
  trig_weight_D = ( trig_weight_D != 1.0 ) ? trig_weight_D : 0.0;

  m_event.get()->weight_event_trig = ( 1.0 - trig_weight_N ) / ( 1.0 - trig_weight_D );

  if ( m_debug ) {
    Info("decorateWeights()", "N = %.2f - D = %.2f", trig_weight_N, trig_weight_D );
    Info("decorateWeights()", "per-event trigger weight = %.2f", m_event.get()->weight_event_trig );
  }

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode HTopMultilepMiniNTupMaker :: defineTagAndProbe ()
{

  // Clear vector branches from previous event
  //
  ANA_CHECK( this->clearBranches("tag_and_probe") );

  // Do this only for dilepton events
  //
  if ( !m_event.get()->dilep ) { return EL::StatusCode::SUCCESS; }

  // Minimal requirement of pT = 10 GeV on all leptons
  //
  for ( auto lep : m_leptons ) {  if ( lep.get()->pt < 10e3 ) { return EL::StatusCode::SUCCESS; } }

  /*
  if ( m_event.get()->isSS01 ) {
    std::cout << "" << std::endl;
    Info("execute()", "===> Entry %u - EventNumber = %u - RunYear = %i", static_cast<uint32_t>(m_numEntry), static_cast<uint32_t>(m_EventNumber), m_RunYear );
    m_debug = true;
  }
  */

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
	  if ( ( lep.get()->prompt == 1 || ( lep.get()->brems == 1 && lep.get()->qmisid == 0 ) ) && ( lep.get()->qmisid == 0 ) ) {
	      lep.get()->tag_SLT = lep.get()->tag_DLT = 1;
	      found_tag = true;
	      break;
	  }
	  ++tag_idx;
      }
      int probe_idx = ( tag_idx ) ? 0 : 1; // Our lepton vector has only 2 components ;-)

      m_isBadTPEvent_SLT = m_isBadTPEvent_DLT = ( !found_tag || m_leptons.at(probe_idx).get()->qmisid == 1 || m_leptons.at(probe_idx).get()->prompt == 1 || ( m_leptons.at(probe_idx).get()->brems == 1 && m_leptons.at(probe_idx).get()->qmisid == 0 ) );

    } else {

      // OS events:
      //
      // If event contains a !prompt or a qmisid, flag the event as bad and return.
      // Else:
      // -) tag: choose randomly
      // -) probe: the other lepton.

      m_isBadTPEvent_SLT = m_isBadTPEvent_DLT = 0;  // be optimistic!

      for ( auto lep : m_leptons ) {
	  if ( lep.get()->prompt == 0 || lep.get()->qmisid == 1 ) {
	      m_isBadTPEvent_SLT = m_isBadTPEvent_DLT = 1;
	      return EL::StatusCode::SUCCESS;
	  }
      }

      int tag_idx = ( m_rand->Rndm() > 0.5 ); // will pick index 0 or 1 in lepton vector randomly

      m_leptons.at(tag_idx).get()->tag_SLT = m_leptons.at(tag_idx).get()->tag_DLT = 1;

    }

    return EL::StatusCode::SUCCESS;
  }

  if ( m_useSUSYSSTP ) {

    // Tag-and-probe a-la-SUSY SS analysis:
    //
    // In the Real OS CR, ambiguous events where both leptons are T & T.M. are double counted -->
    // both leptons are alternatively considered as the possible tag/probe
    // to avoid any bias in the choice of the tag and to increase the statistics

    m_isBadTPEvent_SLT = m_isBadTPEvent_DLT = 0; // be optimistic!

    int tag_candidate_counter_SLT(0), tag_candidate_counter_DLT(0);
    int idx(0);
    int tag_idx_SLT(-1), tag_idx_DLT(-1), probe_idx_SLT(-1), probe_idx_DLT(-1);

    // Count the number of tag lepton candidates, for both SLT and DLT

    if ( m_debug ) { std::cout << "Looking for tag candidates:\n" << std::endl; }

    for ( auto lep : m_leptons ) {

	if ( m_debug ) { Info("defineTagAndProbe()","Checking lepton[%i] w/ pT = %.2f, flavour = %i", idx, lep.get()->pt/1e3, lep.get()->flavour ); }

      if ( lep.get()->tight && lep.get()->trigmatched ) {
        ++tag_candidate_counter_SLT;
	tag_idx_SLT = idx;
        if ( m_debug ) { Info("defineTagAndProbe()","\t ===> found a tag candidate (SLT matching)! pT[%i] = %.2f, flavour = %i", idx, lep.get()->pt/1e3, lep.get()->flavour ); }
      }
      if ( lep.get()->tight && lep.get()->trigmatched_DLT ) {
        ++tag_candidate_counter_DLT;
	tag_idx_DLT = idx;
        if ( m_debug ) { Info("defineTagAndProbe()","\t ===> found a tag candidate (DLT matching)! pT[%i] = %.2f, flavour = %i", idx, lep.get()->pt/1e3, lep.get()->flavour ); }
      }

      ++idx;
    }
    if ( m_debug ) { std::cout << "" << std::endl; }

    float tag_pt(-1.0), probe_pt(-1.0);

    // Different treatment for OS and SS events

    if ( !m_event.get()->isSS01 ) {

	// SLT

	switch ( tag_candidate_counter_SLT ) {

	case 0:

	    m_isBadTPEvent_SLT = 1;

	    if ( m_debug ) { Info("defineTagAndProbe()","No lepton (T & TM) (SLT matching) was found - flag this event as bad" ); }

	    break;

	case 1: // this is the non-ambiguous case: 1 tag and 1 probe per event

	    probe_idx_SLT = abs( tag_idx_SLT - 1 ); // NB: this works since we only look at dilepton events (tag_idx can be 0 or 1 )

	    tag_pt   = m_leptons.at(tag_idx_SLT).get()->pt;
	    probe_pt = m_leptons.at(probe_idx_SLT).get()->pt;

	    m_lep_TagVec_SLT_Pt.push_back( tag_pt );
	    m_lep_ProbeVec_SLT_Pt.push_back( probe_pt );

	    if      ( m_leptons.at(tag_idx_SLT).get()->flavour == 11 ) { m_el_TagVec_SLT_Pt.push_back( tag_pt ); }
	    else if ( m_leptons.at(tag_idx_SLT).get()->flavour == 13 ) { m_mu_TagVec_SLT_Pt.push_back( tag_pt ); }

	    if      ( m_leptons.at(probe_idx_SLT).get()->flavour == 11 ) { m_el_ProbeVec_SLT_Pt.push_back( probe_pt ); }
	    else if ( m_leptons.at(probe_idx_SLT).get()->flavour == 13 ) { m_mu_ProbeVec_SLT_Pt.push_back( probe_pt ); }

	    // This will set the flat branch as well for tag and probe later on (needed to define cuts in plotting tools)

	    m_leptons.at(tag_idx_SLT).get()->tag_SLT = 1;

	    if ( m_debug ) { Info("defineTagAndProbe()","Unambiguous event (SLT matching) : tag lepton pT = %.2f, probe lepton pT = %.2f", m_leptons.at(tag_idx_SLT).get()->pt/1e3, m_leptons.at(probe_idx_SLT).get()->pt/1e3 ); }

	    break;

	case 2: // this is the ambiguous case: both T & TM

	    // In OS, we consider both leptons as tag and probe (event will be double-counted). The reason is b/c we require
	    // both leptons be real in this region (in data we will subtract backgrounds), so here there's really no distinction possible

	    for ( auto lep : m_leptons ) {

		tag_pt = probe_pt = lep.get()->pt;

		m_lep_TagVec_SLT_Pt.push_back( tag_pt );
		m_lep_ProbeVec_SLT_Pt.push_back( probe_pt );

		if ( lep.get()->flavour == 11 ) {

		    m_el_TagVec_SLT_Pt.push_back( tag_pt );
		    m_el_ProbeVec_SLT_Pt.push_back(probe_pt );

		} else if ( lep.get()->flavour == 13 ) {

		    m_mu_TagVec_SLT_Pt.push_back( tag_pt );
		    m_mu_ProbeVec_SLT_Pt.push_back(probe_pt );

		}

	    }

	    // This will allow to set the flat branch as well for tag and probe later on
	    //
	    // For convenience, choose the tag and probe randomly
	    //
	    // This has to be done b/c at plotting level we ask for the probe to be T/L/!T and/or TM/!TM to define the N and D for efficiency,
	    // but in these events both are T and T.M., so it doesn't really matter which lepton we picked as tag/probe
	    // In fact, in this case, when plotting the vector branch "lep_ProbeVec_*", both leptons will be considered effectively as the probe

	    tag_idx_SLT = ( m_rand->Rndm() > 0.5 );

	    m_leptons.at( tag_idx_SLT ).get()->tag_SLT = 1;

	    if ( m_debug ) { Info("defineTagAndProbe()","Ambiguous event (SLT matching): both leptons T & T.M. Event is OS --> will consider both as tag and probe"); }

	    break;

	default:

	    Error("defineTagAndProbe()","Number of tag lepton candidates (SLT matching) is %i...This shouldn't happen. Aborting.", tag_candidate_counter_SLT );
	    return EL::StatusCode::FAILURE;

	} // close switch SLT

	// DLT

	switch ( tag_candidate_counter_DLT ) {

	case 0:

	    m_isBadTPEvent_DLT = 1;

	    if ( m_debug ) { Info("defineTagAndProbe()","No lepton (T & TM) (DLT matching) was found - flag this event as bad" ); }

	    break;

	case 1: // this is the non-ambiguous case: 1 tag and 1 probe per event

	    probe_idx_DLT = abs( tag_idx_DLT - 1 ); // NB: this works since we only look at dilepton events (tag_idx can be 0 or 1 )

	    tag_pt   = m_leptons.at(tag_idx_DLT).get()->pt;
	    probe_pt = m_leptons.at(probe_idx_DLT).get()->pt;

	    m_lep_TagVec_DLT_Pt.push_back( tag_pt );
	    m_lep_ProbeVec_DLT_Pt.push_back( probe_pt );

	    if      ( m_leptons.at(tag_idx_DLT).get()->flavour == 11 ) { m_el_TagVec_DLT_Pt.push_back( tag_pt ); }
	    else if ( m_leptons.at(tag_idx_DLT).get()->flavour == 13 ) { m_mu_TagVec_DLT_Pt.push_back( tag_pt ); }

	    if      ( m_leptons.at(probe_idx_DLT).get()->flavour == 11 ) { m_el_ProbeVec_DLT_Pt.push_back( probe_pt ); }
	    else if ( m_leptons.at(probe_idx_DLT).get()->flavour == 13 ) { m_mu_ProbeVec_DLT_Pt.push_back( probe_pt ); }

	    // This will set the flat branch as well for tag and probe later on (needed to define cuts in plotting tools)

	    m_leptons.at(tag_idx_DLT).get()->tag_DLT = 1;

	    if ( m_debug ) { Info("defineTagAndProbe()","Unambiguous event (DLT matching) : tag lepton pT = %.2f, probe lepton pT = %.2f", m_leptons.at(tag_idx_DLT).get()->pt/1e3, m_leptons.at(probe_idx_DLT).get()->pt/1e3 ); }

	    break;

	case 2: // this is the ambiguous case: both T & TM

	    // In OS, we consider both leptons as tag and probe (event will be double-counted). The reason is b/c we require
	    // both leptons be real in this region (in data we will subtract backgrounds), so here there's really no distinction possible

	    for ( auto lep : m_leptons ) {

		tag_pt = probe_pt = lep.get()->pt;

		m_lep_TagVec_DLT_Pt.push_back( tag_pt );
		m_lep_ProbeVec_DLT_Pt.push_back( probe_pt );

		if ( lep.get()->flavour == 11 ) {

		    m_el_TagVec_DLT_Pt.push_back( tag_pt );
		    m_el_ProbeVec_DLT_Pt.push_back(probe_pt );

		} else if ( lep.get()->flavour == 13 ) {

		    m_mu_TagVec_DLT_Pt.push_back( tag_pt );
		    m_mu_ProbeVec_DLT_Pt.push_back(probe_pt );

		}

	    }

	    // This will allow to set the flat branch as well for tag and probe later on
	    //
	    // For convenience, choose the tag and probe randomly
	    //
	    // This has to be done b/c at plotting level we ask for the probe to be T/L/!T and/or TM/!TM to define the N and D for efficiency,
	    // but in these events both are T and T.M., so it doesn't really matter which lepton we picked as tag/probe
	    // In fact, in this case, when plotting the vector branch "lep_ProbeVec_*", both leptons will be considered effectively as the probe

	    tag_idx_DLT = ( m_rand->Rndm() > 0.5 );

	    m_leptons.at( tag_idx_DLT ).get()->tag_DLT = 1;

	    if ( m_debug ) { Info("defineTagAndProbe()","Ambiguous event (DLT matching): both leptons T & T.M. Event is OS --> will consider both as tag and probe"); }

	    break;

	default:

	    Error("defineTagAndProbe()","Number of tag lepton candidates (DLT matching) is %i...This shouldn't happen. Aborting.", tag_candidate_counter_DLT );
	    return EL::StatusCode::FAILURE;

	} // close switch DLT

    } else {

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

	    int muon_tag_candidate_counter_SLT(0), muon_tag_candidate_counter_DLT(0);

	    if ( m_debug ) { Info("defineTagAndProbe()","Looking at OF event..." ); }

	    // SLT

	    int i_SLT(0);
	    for ( auto lep : m_leptons ) {

		if ( m_debug ) { Info("defineTagAndProbe()","Checking lepton[%i] w/ pT = %.2f, flavour = %i", i_SLT, lep.get()->pt/1e3, lep.get()->flavour ); }

		if ( lep.get()->flavour == 13 && lep.get()->tight && lep.get()->trigmatched ) {

		    tag_idx_SLT = i_SLT;

		    if ( m_debug ) { Info("defineTagAndProbe()","\t ===> found a muon tag candidate (SLT matching)! pT[%i] = %.2f", i_SLT, lep.get()->pt/1e3 ); }

		    ++muon_tag_candidate_counter_SLT;

		    break;

		}

		++i_SLT;
	    }

	    if ( !muon_tag_candidate_counter_SLT ) {

		m_isBadTPEvent_SLT = 1;

		if ( m_debug ) { Info("defineTagAndProbe()","No muon (T & TM) (SLT matching) was found - flag this event as bad" ); }

	    } else {

		probe_idx_SLT = ( tag_idx_SLT ) ? 0 : 1; // Our lepton vector has only 2 components ;-)

		tag_pt   = m_leptons.at(tag_idx_SLT).get()->pt;
		probe_pt = m_leptons.at(probe_idx_SLT).get()->pt;

		m_lep_TagVec_SLT_Pt.push_back( tag_pt );
		m_lep_ProbeVec_SLT_Pt.push_back( probe_pt );

		if      ( m_leptons.at(tag_idx_SLT).get()->flavour == 11 ) { m_el_TagVec_SLT_Pt.push_back( tag_pt ); }
		else if ( m_leptons.at(tag_idx_SLT).get()->flavour == 13 ) { m_mu_TagVec_SLT_Pt.push_back( tag_pt ); }

		if      ( m_leptons.at(probe_idx_SLT).get()->flavour == 11 ) { m_el_ProbeVec_SLT_Pt.push_back( probe_pt ); }
		else if ( m_leptons.at(probe_idx_SLT).get()->flavour == 13 ) { m_mu_ProbeVec_SLT_Pt.push_back( probe_pt ); }

		// This will set the flat branch as well for tag and probe later on (needed to define cuts in plotting tools)

		m_leptons.at(tag_idx_SLT).get()->tag_SLT = 1;

		if ( m_debug ) { Info("defineTagAndProbe()","Good OF event (SLT matching) : tag lepton pT = %.2f, flavour = %i, probe lepton pT = %.2f, flavour = %i", m_leptons.at(tag_idx_SLT).get()->pt/1e3, m_leptons.at(tag_idx_SLT).get()->flavour, m_leptons.at(probe_idx_SLT).get()->pt/1e3, m_leptons.at(probe_idx_SLT).get()->flavour ); }

	    }

	    // DLT

	    int i_DLT(0);
	    for ( auto lep : m_leptons ) {

		if ( m_debug ) { Info("defineTagAndProbe()","Checking lepton[%i] w/ pT = %.2f, flavour = %i", i_DLT, lep.get()->pt/1e3, lep.get()->flavour ); }

		if ( lep.get()->flavour == 13 && lep.get()->tight && lep.get()->trigmatched_DLT ) {

		    tag_idx_DLT = i_DLT;

		    if ( m_debug ) { Info("defineTagAndProbe()","\t ===> found a muon tag candidate (DLT matching)! pT[%i] = %.2f", i_DLT, lep.get()->pt/1e3 ); }

		    ++muon_tag_candidate_counter_DLT;

		    break;

		}

		++i_DLT;
	    }

	    if ( !muon_tag_candidate_counter_DLT ) {

		m_isBadTPEvent_DLT = 1;

		if ( m_debug ) { Info("defineTagAndProbe()","No muon (T & TM) (DLT matching) was found - flag this event as bad" ); }

	    } else {

		probe_idx_DLT = ( tag_idx_DLT ) ? 0 : 1; // Our lepton vector has only 2 components ;-)

		tag_pt   = m_leptons.at(tag_idx_DLT).get()->pt;
		probe_pt = m_leptons.at(probe_idx_DLT).get()->pt;

		m_lep_TagVec_DLT_Pt.push_back( tag_pt );
		m_lep_ProbeVec_DLT_Pt.push_back( probe_pt );

		if      ( m_leptons.at(tag_idx_DLT).get()->flavour == 11 ) { m_el_TagVec_DLT_Pt.push_back( tag_pt ); }
		else if ( m_leptons.at(tag_idx_DLT).get()->flavour == 13 ) { m_mu_TagVec_DLT_Pt.push_back( tag_pt ); }

		if      ( m_leptons.at(probe_idx_DLT).get()->flavour == 11 ) { m_el_ProbeVec_DLT_Pt.push_back( probe_pt ); }
		else if ( m_leptons.at(probe_idx_DLT).get()->flavour == 13 ) { m_mu_ProbeVec_DLT_Pt.push_back( probe_pt ); }

		// This will set the flat branch as well for tag and probe later on (needed to define cuts in plotting tools)

		m_leptons.at(tag_idx_DLT).get()->tag_DLT = 1;

		if ( m_debug ) { Info("defineTagAndProbe()","Good OF event (DLT matching) : tag lepton pT = %.2f, flavour = %i, probe lepton pT = %.2f, flavour = %i", m_leptons.at(tag_idx_DLT).get()->pt/1e3, m_leptons.at(tag_idx_DLT).get()->flavour, m_leptons.at(probe_idx_DLT).get()->pt/1e3, m_leptons.at(probe_idx_DLT).get()->flavour ); }

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

	    m_isBadTPEvent_SLT = 1;

	    if ( m_debug ) { Info("defineTagAndProbe()","No lepton (T & TM) (SLT matching) was found - flag this event as bad" ); }

	    break;

	case 1: // this is the non-ambiguous case: 1 tag and 1 probe per event

	    probe_idx_SLT = abs( tag_idx_SLT - 1 ); // NB: this works since we only look at dilepton events (tag_idx can be 0 or 1 )

	    tag_pt   = m_leptons.at(tag_idx_SLT).get()->pt;
	    probe_pt = m_leptons.at(probe_idx_SLT).get()->pt;

	    m_lep_TagVec_SLT_Pt.push_back( tag_pt );
	    m_lep_ProbeVec_SLT_Pt.push_back( probe_pt );

	    if      ( m_leptons.at(tag_idx_SLT).get()->flavour == 11 ) { m_el_TagVec_SLT_Pt.push_back( tag_pt ); }
	    else if ( m_leptons.at(tag_idx_SLT).get()->flavour == 13 ) { m_mu_TagVec_SLT_Pt.push_back( tag_pt ); }

	    if      ( m_leptons.at(probe_idx_SLT).get()->flavour == 11 ) { m_el_ProbeVec_SLT_Pt.push_back( probe_pt ); }
	    else if ( m_leptons.at(probe_idx_SLT).get()->flavour == 13 ) { m_mu_ProbeVec_SLT_Pt.push_back( probe_pt ); }

	    // This will set the flat branch as well for tag and probe later on (needed to define cuts in plotting tools)

	    m_leptons.at(tag_idx_SLT).get()->tag_SLT = 1;

	    if ( m_debug ) { Info("defineTagAndProbe()","Unambiguous event (SLT matching) : tag lepton pT = %.2f, probe lepton pT = %.2f", m_leptons.at(tag_idx_SLT).get()->pt/1e3, m_leptons.at(probe_idx_SLT).get()->pt/1e3 ); }

	    break;

	case 2: // this is the ambiguous case: both T & TM

	    // In SS, even in the ambiguous case we still need to tell which is tag and which is probe,
	    // since we really do want to take the fake as the probe! (unlike OS, where we require both leptons be real, in SS we have 1 real and 1 fake)

	    if ( m_ambiSolvingCrit.compare("Pt") == 0 ||  m_ambiSolvingCrit.compare("OF") == 0 ) { // this will be used also when the solving criterion is OF, but the event is SF

		// Make sure lepton container is sorted in descending order of lepton pt

		std::sort( m_leptons.begin(), m_leptons.end(), MiniNTupMaker::SorterPt() );
		if ( m_debug ) {
		    std:: cout << "\nLepton container ( descending pT sorting ):\n" << std::endl;
		    for ( unsigned int idx(0); idx < m_leptons.size(); ++idx ) {
			std::cout << "lepton[" << idx << "] - DeltaR(lep, closest bjet) = " << m_leptons.at(idx).get()->deltaRClosestBJet << ", pT = " << m_leptons.at(idx).get()->pt/1e3 << std::endl;
		    }
		    std::cout << "" << std::endl;
		}

		// Take the leading as the tag, the subleading as the probe

		tag_idx_SLT   = 0;
		probe_idx_SLT = 1;

		/*
		// We take the following approach to select it:
		//
		// -) require *at least* one lepton of pT > 40 GeV, otherwise flag the event as bad
		// -) if the other lepton has pT < 40 GeV, that will be the probe, the other the tag
		// -) if both leptons have pT > 40 GeV, choose the tag/probe randomly (rather look at isolation?)

		if ( m_leptons.at(0).get()->pt < 40e3  && m_leptons.at(1).get()->pt < 40e3 ) {

		    m_isBadTPEvent_SLT = 1;

		    break;

		} else if ( m_leptons.at(0).get()->pt >= 40e3 && m_leptons.at(1).get()->pt < 40e3 ) {

		    tag_idx_SLT   = 0;
		    probe_idx_SLT = 1;

		} else if ( m_leptons.at(0).get()->pt >= 40e3 && m_leptons.at(1).get()->pt < 40e3 ) {

		    tag_idx_SLT   = 1;
		    probe_idx_SLT = 0;

		} else {

		    tag_idx_SLT   = ( m_rand->Rndm() > 0.5 );
		    probe_idx_SLT = ( tag_idx_SLT ) ? 0: 1; // Our lepton vector has only 2 components ;-)

		}
		*/

	    } else if ( m_ambiSolvingCrit.compare("deltaRClosestBJet") == 0 ) {

		// Try an event topology-based approach to select the tag and probe:
		//
		// -) take the lepton that is further apart from the closest bjet as the tag, the other will be the probe

		// Require *at least* one bjet, otherwise flag the event as bad

		if ( m_event.get()->nbjets < 1 ) {

		    m_isBadTPEvent_SLT = 1;

		    break;
		}

		// Sort lepton container in descending order of distance to closest bjet

		std::sort( m_leptons.begin(), m_leptons.end(), MiniNTupMaker::SorterDistanceClosestBJet() );
		if ( m_debug ) {
		    std:: cout << "\nLepton container ( descending DeltaR(lep, closest bjet) sorting ):\n" << std::endl;
		    for ( unsigned int idx(0); idx < m_leptons.size(); ++idx ) {
			std::cout << "lepton[" << idx << "] - DeltaR(lep, closest bjet) = " << m_leptons.at(idx).get()->deltaRClosestBJet << ", pT = " << m_leptons.at(idx).get()->pt/1e3 << std::endl;
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

		    m_isBadTPEvent_SLT = 1;

		    break;
		}

		// Sort lepton container in descending order of mass w/ closest bjet

		std::sort( m_leptons.begin(), m_leptons.end(), MiniNTupMaker::SorterMassClosestBJet() );
		if ( m_debug ) {
		    std:: cout << "\nLepton container ( descending M(lep, closest bjet) sorting ):\n" << std::endl;
		    for ( unsigned int idx(0); idx < m_leptons.size(); ++idx ) {
			std::cout << "lepton[" << idx << "] - M(lep, closest bjet) = " << m_leptons.at(idx).get()->massClosestBJet/1e3 << ", pT = " << m_leptons.at(idx).get()->pt/1e3 << std::endl;
		    }
		    std::cout << "" << std::endl;
		}

		tag_idx_SLT   = 0;
		probe_idx_SLT = 1;

	    }

	    // Flag the tag lepton (used afterwards to define flat branches)

	    m_leptons.at( tag_idx_SLT ).get()->tag_SLT = 1;

	    tag_pt   = m_leptons.at( tag_idx_SLT ).get()->pt;
	    probe_pt = m_leptons.at( probe_idx_SLT ).get()->pt;

	    m_lep_TagVec_SLT_Pt.push_back( tag_pt );
	    m_lep_ProbeVec_SLT_Pt.push_back( probe_pt );

	    if      ( m_leptons.at( tag_idx_SLT ).get()->flavour == 11 ) { m_el_TagVec_SLT_Pt.push_back( tag_pt ); }
	    else if ( m_leptons.at( tag_idx_SLT ).get()->flavour == 13 ) { m_mu_TagVec_SLT_Pt.push_back( tag_pt ); }

	    if	    ( m_leptons.at( probe_idx_SLT ).get()->flavour == 11 ) { m_el_ProbeVec_SLT_Pt.push_back( probe_pt ); }
	    else if ( m_leptons.at( probe_idx_SLT ).get()->flavour == 13 ) { m_mu_ProbeVec_SLT_Pt.push_back( probe_pt ); }

	    if ( m_debug && m_ambiSolvingCrit.compare("Pt") == 0  )               { Info("defineTagAndProbe()","Ambiguous event (SLT matching): both leptons T & T.M. Event is SS --> tag lepton pT = %.2f, probe lepton pT = %.2f", tag_pt/1e3, probe_pt/1e3 ); }
	    if ( m_debug && m_ambiSolvingCrit.compare("deltaRClosestBJet") == 0 ) { Info("defineTagAndProbe()","Ambiguous event (SLT matching): both leptons T & T.M. Event is SS --> tag lepton pT = %.2f, DeltaR(lep, closest bjet) = %.3f, probe lepton pT = %.2f, DeltaR(lep, closest bjet) = %.3f", tag_pt/1e3, m_leptons.at(tag_idx_SLT).get()->deltaRClosestBJet, probe_pt/1e3 , m_leptons.at(probe_idx_SLT).get()->deltaRClosestBJet ); }
	    if ( m_debug && m_ambiSolvingCrit.compare("massClosestBJet") == 0 )   { Info("defineTagAndProbe()","Ambiguous event (SLT matching): both leptons T & T.M. Event is SS --> tag lepton pT = %.2f, M(lep, closest bjet) = %.3f, probe lepton pT = %.2f, M(lep, closest bjet) = %.3f", tag_pt/1e3, m_leptons.at(tag_idx_SLT).get()->massClosestBJet/1e3, probe_pt/1e3 , m_leptons.at(probe_idx_SLT).get()->massClosestBJet/1e3 ); }

	    // Revert back to pT-sorting (which is the default)

	    std::sort( m_leptons.begin(), m_leptons.end(), MiniNTupMaker::SorterPt() );

	    break;

	default:

	    Error("defineTagAndProbe()","Number of tag lepton candidates (SLT matching) is %i...This shouldn't happen. Aborting.", tag_candidate_counter_SLT );
	    return EL::StatusCode::FAILURE;

	} // close switch SLT

	// DLT

	switch ( tag_candidate_counter_DLT ) {

	case 0:

	    m_isBadTPEvent_DLT = 1;

	    if ( m_debug ) { Info("defineTagAndProbe()","No lepton (T & TM) (DLT matching) was found - flag this event as bad" ); }

	    break;

	case 1: // this is the non-ambiguous case: 1 tag and 1 probe per event

	    probe_idx_DLT = abs( tag_idx_DLT - 1 ); // NB: this works since we only look at dilepton events (tag_idx can be 0 or 1 )

	    tag_pt   = m_leptons.at(tag_idx_DLT).get()->pt;
	    probe_pt = m_leptons.at(probe_idx_DLT).get()->pt;

	    m_lep_TagVec_DLT_Pt.push_back( tag_pt );
	    m_lep_ProbeVec_DLT_Pt.push_back( probe_pt );

	    if      ( m_leptons.at(tag_idx_DLT).get()->flavour == 11 ) { m_el_TagVec_DLT_Pt.push_back( tag_pt ); }
	    else if ( m_leptons.at(tag_idx_DLT).get()->flavour == 13 ) { m_mu_TagVec_DLT_Pt.push_back( tag_pt ); }

	    if      ( m_leptons.at(probe_idx_DLT).get()->flavour == 11 ) { m_el_ProbeVec_DLT_Pt.push_back( probe_pt ); }
	    else if ( m_leptons.at(probe_idx_DLT).get()->flavour == 13 ) { m_mu_ProbeVec_DLT_Pt.push_back( probe_pt ); }

	    // This will set the flat branch as well for tag and probe later on (needed to define cuts in plotting tools)

	    m_leptons.at(tag_idx_DLT).get()->tag_DLT = 1;

	    if ( m_debug ) { Info("defineTagAndProbe()","Unambiguous event (DLT matching) : tag lepton pT = %.2f, probe lepton pT = %.2f", m_leptons.at(tag_idx_DLT).get()->pt/1e3, m_leptons.at(probe_idx_DLT).get()->pt/1e3 ); }

	    break;

	case 2: // this is the ambiguous case: both T & TM

	    // In SS, even in the ambiguous case we still need to tell which is tag and which is probe,
	    // since we really do want to take the fake as the probe! (unlike OS, where we require both leptons be real, in SS we have 1 real and 1 fake)

	    if ( m_ambiSolvingCrit.compare("Pt") == 0 ||  m_ambiSolvingCrit.compare("OF") == 0 ) { // this will be used also when the solving criterion is OF, but the event is SF

		// Make sure lepton container is sorted in descending order of lepton pt

		std::sort( m_leptons.begin(), m_leptons.end(), MiniNTupMaker::SorterPt() );
		if ( m_debug ) {
		    std:: cout << "\nLepton container ( descending pT sorting ):\n" << std::endl;
		    for ( unsigned int idx(0); idx < m_leptons.size(); ++idx ) {
			std::cout << "lepton[" << idx << "] - DeltaR(lep, closest bjet) = " << m_leptons.at(idx).get()->deltaRClosestBJet << ", pT = " << m_leptons.at(idx).get()->pt/1e3 << std::endl;
		    }
		    std::cout << "" << std::endl;
		}

		// Take the leading as the tag, the subleading as the probe

		tag_idx_DLT   = 0;
		probe_idx_DLT = 1;

		/*
		// We take the following approach to select it:
		//
		// -) require *at least* one lepton of pT > 40 GeV, otherwise flag the event as bad
		// -) if the other lepton has pT < 40 GeV, that will be the probe, the other the tag
		// -) if both leptons have pT > 40 GeV, choose the tag/probe randomly (rather look at isolation?)

		if ( m_leptons.at(0).get()->pt < 40e3  && m_leptons.at(1).get()->pt < 40e3 ) {

		    m_isBadTPEvent_DLT = 1;

		    break;

		} else if ( m_leptons.at(0).get()->pt >= 40e3 && m_leptons.at(1).get()->pt < 40e3 ) {

		    tag_idx_DLT   = 0;
		    probe_idx_DLT = 1;

		} else if ( m_leptons.at(0).get()->pt >= 40e3 && m_leptons.at(1).get()->pt < 40e3 ) {

		    tag_idx_DLT   = 1;
		    probe_idx_DLT = 0;

		} else {

		    tag_idx_DLT   = ( m_rand->Rndm() > 0.5 );
		    probe_idx_DLT = ( tag_idx_DLT ) ? 0: 1; // Our lepton vector has only 2 components ;-)

		}
		*/

	    } else if ( m_ambiSolvingCrit.compare("deltaRClosestBJet") == 0 ) {

		// Try an event topology-based approach to select the tag and probe:
		//
		// -) take the lepton that is further apart from the closest bjet as the tag, the other will be the probe

		// Require *at least* one bjet, otherwise flag the event as bad

		if ( m_event.get()->nbjets < 1 ) {

		    m_isBadTPEvent_DLT = 1;

		    break;
		}

		// Sort lepton container in descending order of distance to closest bjet

		std::sort( m_leptons.begin(), m_leptons.end(), MiniNTupMaker::SorterDistanceClosestBJet() );
		if ( m_debug ) {
		    std:: cout << "\nLepton container ( descending DeltaR(lep, closest bjet) sorting ):\n" << std::endl;
		    for ( unsigned int idx(0); idx < m_leptons.size(); ++idx ) {
			std::cout << "lepton[" << idx << "] - DeltaR(lep, closest bjet) = " << m_leptons.at(idx).get()->deltaRClosestBJet << ", pT = " << m_leptons.at(idx).get()->pt/1e3 << std::endl;
		    }
		    std::cout << "" << std::endl;
		}

		tag_idx_DLT   = 0;
		probe_idx_DLT = 1;

	    } else if ( m_ambiSolvingCrit.compare("massClosestBJet") == 0  ) {

		// Try an event topology-based approach to select the tag and probe:
		//
		// -) take the lepton forming the largest invariant mass w/ the closest bjet as the tag, the other will be the probe

		// Require *at least* one bjet, otherwise flag the event as bad

		if ( m_event.get()->nbjets < 1 ) {

		    m_isBadTPEvent_DLT = 1;

		    break;
		}

		// Sort lepton container in descending order of mass w/ closest bjet

		std::sort( m_leptons.begin(), m_leptons.end(), MiniNTupMaker::SorterMassClosestBJet() );
		if ( m_debug ) {
		    std:: cout << "\nLepton container ( descending M(lep, closest bjet) sorting ):\n" << std::endl;
		    for ( unsigned int idx(0); idx < m_leptons.size(); ++idx ) {
			std::cout << "lepton[" << idx << "] - M(lep, closest bjet) = " << m_leptons.at(idx).get()->massClosestBJet/1e3 << ", pT = " << m_leptons.at(idx).get()->pt/1e3 << std::endl;
		    }
		    std::cout << "" << std::endl;
		}

		tag_idx_DLT   = 0;
		probe_idx_DLT = 1;

	    }

	    // Flag the tag lepton (used afterwards to define flat branches)

	    m_leptons.at( tag_idx_DLT ).get()->tag_DLT = 1;

	    tag_pt   = m_leptons.at( tag_idx_DLT ).get()->pt;
	    probe_pt = m_leptons.at( probe_idx_DLT ).get()->pt;

	    m_lep_TagVec_DLT_Pt.push_back( tag_pt );
	    m_lep_ProbeVec_DLT_Pt.push_back( probe_pt );

	    if      ( m_leptons.at( tag_idx_DLT ).get()->flavour == 11 ) { m_el_TagVec_DLT_Pt.push_back( tag_pt ); }
	    else if ( m_leptons.at( tag_idx_DLT ).get()->flavour == 13 ) { m_mu_TagVec_DLT_Pt.push_back( tag_pt ); }

	    if	    ( m_leptons.at( probe_idx_DLT ).get()->flavour == 11 ) { m_el_ProbeVec_DLT_Pt.push_back( probe_pt ); }
	    else if ( m_leptons.at( probe_idx_DLT ).get()->flavour == 13 ) { m_mu_ProbeVec_DLT_Pt.push_back( probe_pt ); }

	    if ( m_debug && m_ambiSolvingCrit.compare("Pt") == 0  )               { Info("defineTagAndProbe()","Ambiguous event (DLT matching): both leptons T & T.M. Event is SS --> tag lepton pT = %.2f, probe lepton pT = %.2f", tag_pt/1e3, probe_pt/1e3 ); }
	    if ( m_debug && m_ambiSolvingCrit.compare("deltaRClosestBJet") == 0 ) { Info("defineTagAndProbe()","Ambiguous event (DLT matching): both leptons T & T.M. Event is SS --> tag lepton pT = %.2f, DeltaR(lep, closest bjet) = %.3f, probe lepton pT = %.2f, DeltaR(lep, closest bjet) = %.3f", tag_pt/1e3, m_leptons.at(tag_idx_DLT).get()->deltaRClosestBJet, probe_pt/1e3 , m_leptons.at(probe_idx_DLT).get()->deltaRClosestBJet ); }
	    if ( m_debug && m_ambiSolvingCrit.compare("massClosestBJet") == 0 )   { Info("defineTagAndProbe()","Ambiguous event (DLT matching): both leptons T & T.M. Event is SS --> tag lepton pT = %.2f, M(lep, closest bjet) = %.3f, probe lepton pT = %.2f, M(lep, closest bjet) = %.3f", tag_pt/1e3, m_leptons.at(tag_idx_DLT).get()->massClosestBJet/1e3, probe_pt/1e3 , m_leptons.at(probe_idx_DLT).get()->massClosestBJet/1e3 ); }

	    // Revert back to pT-sorting (which is the default)

	    std::sort( m_leptons.begin(), m_leptons.end(), MiniNTupMaker::SorterPt() );

	    break;

	default:

	    Error("defineTagAndProbe()","Number of tag lepton candidates (DLT matching) is %i...This shouldn't happen. Aborting.", tag_candidate_counter_DLT );
	    return EL::StatusCode::FAILURE;

	} // close switch DLT

    } // close OS-SS cases

    return EL::StatusCode::SUCCESS;

  } // close SUSYTP

  bool found_tag_SLT(false), found_tag_DLT(false);

  // The first trigger-matched. tight lepton encountered in the event will be the tag (NB: the lepton container is pT-ordered).
  // Events where no such lepton is found are flagged and will be discarded

  // By construction, in the ambiguous case where both are (T & T.M.), the leading will be taken as the tag, the subleading will be the probe

  m_isBadTPEvent_SLT = m_isBadTPEvent_DLT = 0; // be optimistic!

  int idx(0);

  for ( auto lep : m_leptons ) {

    if ( m_debug ) { Info("defineTagAndProbe()","Checking lepton[%i] w/ pT = %.2f", idx, lep.get()->pt/1e3 ); }

    if ( !found_tag_SLT && ( lep.get()->tight && lep.get()->trigmatched ) ) {
      lep.get()->tag_SLT = 1;
      found_tag_SLT = true;
      if ( m_debug ) { Info("defineTagAndProbe()","\t ===> found a tag candidate (SLT matching)! pT[%i] = %.2f", idx, lep.get()->pt/1e3 ); }
    }
    if ( !found_tag_DLT && ( lep.get()->tight && lep.get()->trigmatched_DLT ) ) {
      lep.get()->tag_DLT = 1;
      found_tag_DLT = true;
      if ( m_debug ) { Info("defineTagAndProbe()","\t ===> found a tag candidate (DLT matching)! pT[%i] = %.2f", idx, lep.get()->pt/1e3 ); }
    }

    ++idx;
  }

  if ( !found_tag_SLT ) {
    m_leptons.at(0).get()->tag_SLT = 1;
    m_isBadTPEvent_SLT = 1;
    if ( m_debug ) { Info("defineTagAndProbe()","None lepton is T & TM (SLT matching) - choose leading as tag (pT = %.2f) and flag this event as bad", m_leptons.at(0).get()->pt/1e3 ); }
  }
  if ( !found_tag_DLT ) {
    m_leptons.at(0).get()->tag_DLT = 1;
    m_isBadTPEvent_DLT = 1;
    if ( m_debug ) { Info("defineTagAndProbe()","None lepton is T & TM (DLT matching) - choose leading as tag (pT = %.2f) and flag this event as bad", m_leptons.at(0).get()->pt/1e3 ); }
  }

  return EL::StatusCode::SUCCESS;

}

EL::StatusCode HTopMultilepMiniNTupMaker :: fillTPFlatBranches ( std::shared_ptr<MiniNTupMaker::leptonObj> lep, const std::string& trig ) {

    bool isTag(false);
    bool isBadTPEvent(false); // be optimistic!

    std::string this_tp_trig("");
    std::vector<std::string> tp_trigs;

    if ( trig.compare("SLT") == 0 ) {
	isTag = lep.get()->tag_SLT;
	this_tp_trig = ( isTag ) ?  "Tag_SLT" : "Probe_SLT";
	isBadTPEvent = m_isBadTPEvent_SLT;
	if ( isBadTPEvent ) {
	    tp_trigs.push_back("Tag_SLT");
	    tp_trigs.push_back("Probe_SLT");
	} else {
	    tp_trigs.push_back(this_tp_trig);
	}
    } else if ( trig.compare("DLT") == 0 ) {
	isTag = lep.get()->tag_DLT;
	this_tp_trig = ( isTag ) ?  "Tag_DLT" : "Probe_DLT";
	isBadTPEvent = m_isBadTPEvent_DLT;
	if ( isBadTPEvent ) {
	    tp_trigs.push_back("Tag_DLT");
	    tp_trigs.push_back("Probe_DLT");
	} else {
	    tp_trigs.push_back(this_tp_trig);
	}
    } else {
	Error("fillTPFlatBranches()","Invalid parameter: %s for this function!", trig.c_str() );
	return EL::StatusCode::FAILURE;
    }

    for ( const auto& tp_trig : tp_trigs ) {
	m_TagProbe_branches["lep_" + tp_trig + "_Pt"].f              = ( !isBadTPEvent ) ? lep.get()->pt : -1.0;
	m_TagProbe_branches["lep_" + tp_trig + "_Eta"].f             = ( !isBadTPEvent ) ? lep.get()->eta : -999.0;
	m_TagProbe_branches["lep_" + tp_trig + "_EtaBE2"].f          = ( !isBadTPEvent ) ? lep.get()->etaBE2 : -999.0;
	m_TagProbe_branches["lep_" + tp_trig + "_ptVarcone20"].f     = ( !isBadTPEvent ) ? lep.get()->ptVarcone20 : -1.0;
	m_TagProbe_branches["lep_" + tp_trig + "_ptVarcone30"].f     = ( !isBadTPEvent ) ? lep.get()->ptVarcone30 : -1.0;
	m_TagProbe_branches["lep_" + tp_trig + "_topoEtcone20"].f    = ( !isBadTPEvent ) ? lep.get()->topoEtcone20 : -1.0;
	m_TagProbe_branches["lep_" + tp_trig + "_sigd0PV"].f         = ( !isBadTPEvent ) ? lep.get()->d0sig : -1.0;
	m_TagProbe_branches["lep_" + tp_trig + "_Z0SinTheta"].f      = ( !isBadTPEvent ) ? lep.get()->z0sintheta : -1.0;
	m_TagProbe_branches["lep_" + tp_trig + "_deltaRClosestBJet"].f  = ( !isBadTPEvent ) ? lep.get()->deltaRClosestBJet : -1.0;
	m_TagProbe_branches["lep_" + tp_trig + "_massClosestBJet"].f    = ( !isBadTPEvent ) ? lep.get()->massClosestBJet : -1.0;
	m_TagProbe_branches["lep_" + tp_trig + "_ID"].f              = ( !isBadTPEvent ) ? lep.get()->ID : 0.0;
	m_TagProbe_branches["lep_" + tp_trig + "_isTrigMatch"].c     = ( !isBadTPEvent ) ? lep.get()->trigmatched : 0;
	m_TagProbe_branches["lep_" + tp_trig + "_isTightSelected"].c = ( !isBadTPEvent ) ? lep.get()->tight : 0;
	m_TagProbe_branches["lep_" + tp_trig + "_isPrompt"].c        = ( !isBadTPEvent ) ? lep.get()->prompt : 0;
	m_TagProbe_branches["lep_" + tp_trig + "_isBrems"].c         = ( !isBadTPEvent ) ? lep.get()->brems : 0;
	m_TagProbe_branches["lep_" + tp_trig + "_isFakeLep"].c       = ( !isBadTPEvent ) ? lep.get()->fake : 0;
	m_TagProbe_branches["lep_" + tp_trig + "_isQMisID"].c        = ( !isBadTPEvent ) ? lep.get()->qmisid : 0;
	m_TagProbe_branches["lep_" + tp_trig + "_isConvPh"].c        = ( !isBadTPEvent ) ? lep.get()->convph : 0;
	m_TagProbe_branches["lep_" + tp_trig + "_truthType"].i       = ( !isBadTPEvent ) ? lep.get()->truthType : 0;
	m_TagProbe_branches["lep_" + tp_trig + "_truthOrigin"].i     = ( !isBadTPEvent ) ? lep.get()->truthOrigin : 0;
    }

    return EL::StatusCode::SUCCESS;

}

EL::StatusCode HTopMultilepMiniNTupMaker :: setOutputBranches ()
{

  // Clear vector branches from previous event
  //
  ANA_CHECK( this->clearBranches("leptons") );

  m_isMC   = m_event.get()->isMC;
  m_isSS01 = m_event.get()->isSS01;
  m_isSS12 = m_event.get()->isSS12;

  m_is_T_T	   = (  m_leptons.at(0).get()->tight &&  m_leptons.at(1).get()->tight );
  m_is_T_AntiT     = (  m_leptons.at(0).get()->tight && !m_leptons.at(1).get()->tight );
  m_is_AntiT_T     = ( !m_leptons.at(0).get()->tight &&  m_leptons.at(1).get()->tight );
  m_is_AntiT_AntiT = ( !m_leptons.at(0).get()->tight && !m_leptons.at(1).get()->tight );

  m_is_Tel_AntiTmu = ( ( m_leptons.at(0).get()->flavour == 11 && m_leptons.at(0).get()->tight )  && ( m_leptons.at(1).get()->flavour == 13 && !m_leptons.at(1).get()->tight ) );
  m_is_Tmu_AntiTel = ( ( m_leptons.at(0).get()->flavour == 13 && m_leptons.at(0).get()->tight )  && ( m_leptons.at(1).get()->flavour == 11 && !m_leptons.at(1).get()->tight ) );
  m_is_AntiTel_Tmu = ( ( m_leptons.at(0).get()->flavour == 11 && !m_leptons.at(0).get()->tight ) && ( m_leptons.at(1).get()->flavour == 13 && m_leptons.at(1).get()->tight ) );
  m_is_AntiTmu_Tel = ( ( m_leptons.at(0).get()->flavour == 13 && !m_leptons.at(0).get()->tight ) && ( m_leptons.at(1).get()->flavour == 11 && m_leptons.at(1).get()->tight ) );

  m_nleptons = m_leptons.size();

  m_lep_isTightSelected_0 = m_leptons.at(0).get()->tight;
  m_lep_isTightSelected_1 = m_leptons.at(1).get()->tight;
  m_lep_isTightSelected_2 = ( m_event.get()->trilep ) ? m_leptons.at(2).get()->tight : -1;

  m_nmuons     = 0;
  m_nelectrons = 0;

  if ( m_event.get()->dilep ) {

    int idx_el(0), idx_mu(0);
    m_el_Pt_0 = m_el_Pt_1 = m_mu_Pt_0 = m_mu_Pt_1 = -1.0;

    // NB: m_leptons is pT sorted by construction
    //
    for ( const auto lep : m_leptons ) {
      if ( lep.get()->flavour == 13 ) {
        if ( idx_mu == 0 )      { m_mu_Pt_0 = lep.get()->pt; ++idx_mu; }
        else if ( idx_mu == 1 ) { m_mu_Pt_1 = lep.get()->pt; ++idx_mu; }
      } else if ( lep.get()->flavour == 11 ) {
        if ( idx_el == 0 )      { m_el_Pt_0 = lep.get()->pt; ++idx_el; }
        else if ( idx_el == 1 ) { m_el_Pt_1 = lep.get()->pt; ++idx_el; }
      }
    }
  }

  for ( const auto lep : m_leptons ) {

    if ( lep.get()->flavour == 13 ) ++m_nmuons;
    if ( lep.get()->flavour == 11 ) ++m_nelectrons;

    // Fill vector branches
    //
    m_lep_Pt.push_back(lep.get()->pt);
    m_lep_Eta.push_back(lep.get()->eta);
    m_lep_EtaBE2.push_back(lep.get()->etaBE2);

    // Fill flat variables for tag/probe leptons
    //
    if ( m_event.get()->dilep ) {
	ANA_CHECK( this->fillTPFlatBranches( lep, "SLT" ) );
	ANA_CHECK( this->fillTPFlatBranches( lep, "DLT" ) );
    }

  }

  m_weight_event      = m_event.get()->weight_event;
  m_weight_event_lep  = m_event.get()->weight_event_lep;
  m_weight_event_trig = m_event.get()->weight_event_trig;
  m_weight_tag        = m_event.get()->weight_tag;
  m_weight_probe      = m_event.get()->weight_probe;

  return EL::StatusCode::SUCCESS;

}

EL::StatusCode HTopMultilepMiniNTupMaker :: clearBranches ( const std::string& type )
{

  if ( type.compare("leptons") == 0 ) {
    m_lep_Pt.clear();
    m_lep_Eta.clear();
    m_lep_EtaBE2.clear();
  } else if ( type.compare("tag_and_probe") == 0 ) {
    m_lep_TagVec_SLT_Pt.clear();
    m_lep_ProbeVec_SLT_Pt.clear();
    m_lep_TagVec_DLT_Pt.clear();
    m_lep_ProbeVec_DLT_Pt.clear();
    m_el_TagVec_SLT_Pt.clear();
    m_el_ProbeVec_SLT_Pt.clear();
    m_el_TagVec_DLT_Pt.clear();
    m_el_ProbeVec_DLT_Pt.clear();
    m_mu_TagVec_SLT_Pt.clear();
    m_mu_ProbeVec_SLT_Pt.clear();
    m_mu_TagVec_DLT_Pt.clear();
    m_mu_ProbeVec_DLT_Pt.clear();
  } else if ( type.compare("jets") == 0 ) {
    m_jet_OR_Pt.clear();
    m_jet_OR_Eta.clear();
    m_jet_OR_Phi.clear();
    m_jet_OR_E.clear();
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

EL::StatusCode HTopMultilepMiniNTupMaker :: jetTruthMatching ()
{

  // Clear vector branches from previous event
  //
  ANA_CHECK( this->clearBranches("jets") );

  float DRCONE(0.4);

  // Find the truth-matching jet via min(deltaR) matching

  TVector3 truthjet, jet;

  for ( unsigned int j_idx(0); j_idx < m_selected_jets_T->size(); ++j_idx ) {

    short j = m_selected_jets_T->at(j_idx);
    // Read the jet branches from the vector *before* OR via the index
    //
    jet.SetPtEtaPhi( m_jet_pt->at(j), m_jet_eta->at(j), m_jet_phi->at(j) );

    if ( m_debug ) { Info("jetTruthMatching()","reco jet - idx = %i, pT = %.2f, eta = %.2f, phi = %.2f",j, jet.Pt()/1e3, jet.Eta(), jet.Phi() ); }

    m_jet_OR_Pt.push_back( m_jet_pt->at(j) );
    m_jet_OR_Eta.push_back( m_jet_eta->at(j) );
    m_jet_OR_Phi.push_back( m_jet_phi->at(j) );
    m_jet_OR_E.push_back( m_jet_E->at(j) );

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
