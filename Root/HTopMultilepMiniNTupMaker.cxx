#include <HTopMultilepAnalysis/HTopMultilepMiniNTupMaker.h>

// ASG status code check
#include <AsgTools/MessageCheck.h>

// this is needed to distribute the algorithm to the workers
ClassImp(HTopMultilepMiniNTupMaker)

HTopMultilepMiniNTupMaker :: HTopMultilepMiniNTupMaker(std::string className) :
    Algorithm(className),
    m_inputNTuple(nullptr),
    m_outputNTuple(nullptr)
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().

  Info("HTopMultilepMiniNTupMaker()", "Calling constructor");
  
  m_outputNTupStreamName = "output";
  m_inputBranches        = "";
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
  
  Info("changeInput()", "Calling changeInput");
  
  m_inputNTuple = wk()->tree();

  ANA_CHECK( this->enableSelectedBranches() );
  
  // Connect the branches of the input tree to the algorithm members
  //
  m_inputNTuple->SetBranchAddress ("EventNumber",   			      &m_EventNumber);
  m_inputNTuple->SetBranchAddress ("RunNumber",   			      &m_RunNumber);

  m_inputNTuple->SetBranchAddress ("dilep_type",  			      &m_dilep_type);
  m_inputNTuple->SetBranchAddress ("trilep_type",  			      &m_trilep_type);

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
  m_inputNTuple->SetBranchAddress ("lep_isTrigMatch_0",   		      &m_lep_isTrigMatch_0);
  m_inputNTuple->SetBranchAddress ("lep_isPrompt_0",   			      &m_lep_isPrompt_0);
  m_inputNTuple->SetBranchAddress ("lep_isBremsElec_0",   		      &m_lep_isBremsElec_0);
  m_inputNTuple->SetBranchAddress ("lep_isFakeLep_0",   		      &m_lep_isFakeLep_0);
  m_inputNTuple->SetBranchAddress ("lep_SFIDLoose_0",   		      &m_lep_SFIDLoose_0);
  m_inputNTuple->SetBranchAddress ("lep_SFIDTight_0",   		      &m_lep_SFIDTight_0);
  m_inputNTuple->SetBranchAddress ("lep_SFTrigLoose_0",   		      &m_lep_SFTrigLoose_0);
  m_inputNTuple->SetBranchAddress ("lep_SFTrigTight_0",   		      &m_lep_SFTrigTight_0);
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
  m_inputNTuple->SetBranchAddress ("lep_isTrigMatch_1",   		      &m_lep_isTrigMatch_1);
  m_inputNTuple->SetBranchAddress ("lep_isPrompt_1",   			      &m_lep_isPrompt_1);
  m_inputNTuple->SetBranchAddress ("lep_isBremsElec_1",   		      &m_lep_isBremsElec_1);
  m_inputNTuple->SetBranchAddress ("lep_isFakeLep_1",   		      &m_lep_isFakeLep_1);
  m_inputNTuple->SetBranchAddress ("lep_SFIDLoose_1",   		      &m_lep_SFIDLoose_1);
  m_inputNTuple->SetBranchAddress ("lep_SFIDTight_1",   		      &m_lep_SFIDTight_1);
  m_inputNTuple->SetBranchAddress ("lep_SFTrigLoose_1",   		      &m_lep_SFTrigLoose_1);
  m_inputNTuple->SetBranchAddress ("lep_SFTrigTight_1",   		      &m_lep_SFTrigTight_1);
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
  m_inputNTuple->SetBranchAddress ("lep_isTrigMatch_2",   		      &m_lep_isTrigMatch_2);
  m_inputNTuple->SetBranchAddress ("lep_isPrompt_2",   			      &m_lep_isPrompt_2);
  m_inputNTuple->SetBranchAddress ("lep_isBremsElec_2",   		      &m_lep_isBremsElec_2);
  m_inputNTuple->SetBranchAddress ("lep_isFakeLep_2",   		      &m_lep_isFakeLep_2);
  m_inputNTuple->SetBranchAddress ("lep_SFIDLoose_2",   		      &m_lep_SFIDLoose_2);
  m_inputNTuple->SetBranchAddress ("lep_SFIDTight_2",   		      &m_lep_SFIDTight_2);
  m_inputNTuple->SetBranchAddress ("lep_SFTrigLoose_2",   		      &m_lep_SFTrigLoose_2);
  m_inputNTuple->SetBranchAddress ("lep_SFTrigTight_2",   		      &m_lep_SFTrigTight_2);
  m_inputNTuple->SetBranchAddress ("lep_SFIsoLoose_2",   		      &m_lep_SFIsoLoose_2);
  m_inputNTuple->SetBranchAddress ("lep_SFIsoTight_2",   		      &m_lep_SFIsoTight_2);
  m_inputNTuple->SetBranchAddress ("lep_SFReco_2",   			      &m_lep_SFReco_2);
  m_inputNTuple->SetBranchAddress ("lep_SFTTVA_2",   			      &m_lep_SFTTVA_2);
  m_inputNTuple->SetBranchAddress ("lep_SFObjLoose_2",   		      &m_lep_SFObjLoose_2);
  m_inputNTuple->SetBranchAddress ("lep_SFObjTight_2",   		      &m_lep_SFObjTight_2);
 
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
  m_outputNTuple->tree()->Branch("m_isMC",               	&m_isMC, "m_isMC/B");
  m_outputNTuple->tree()->Branch("m_isSS01",               	&m_isSS01, "m_isSS01/B");
  m_outputNTuple->tree()->Branch("m_isSS12",               	&m_isSS12, "m_isSS12/B");
  m_outputNTuple->tree()->Branch("m_is_T_T",               	&m_is_T_T, "m_is_T_T/B");
  m_outputNTuple->tree()->Branch("m_is_T_AntiT",               	&m_is_T_AntiT, "m_is_T_AntiT/B");  
  m_outputNTuple->tree()->Branch("m_is_AntiT_T",               	&m_is_AntiT_T, "m_is_AntiT_T/B");
  m_outputNTuple->tree()->Branch("m_is_AntiT_AntiT",            &m_is_AntiT_AntiT, "m_is_AntiT_AntiT/B");  
  
  m_outputNTuple->tree()->Branch("m_nmuons",               	&m_nmuons, "m_nmuons/I");
  m_outputNTuple->tree()->Branch("m_nelectrons",               	&m_nelectrons, "m_nelectrons/I");
  m_outputNTuple->tree()->Branch("m_nleptons",               	&m_nleptons, "m_nleptons/I");
  
  m_outputNTuple->tree()->Branch("m_lep_isTightSelected_0",     &m_lep_isTightSelected_0, "m_lep_isTightSelected_0/B");
  m_outputNTuple->tree()->Branch("m_lep_isTightSelected_1",     &m_lep_isTightSelected_1, "m_lep_isTightSelected_1/B");
  m_outputNTuple->tree()->Branch("m_lep_isTightSelected_2",     &m_lep_isTightSelected_2, "m_lep_isTightSelected_2/B");
  
  m_outputNTuple->tree()->Branch("m_lep_Tag_Pt",               	&m_lep_Tag_Pt		   , "m_lep_Tag_Pt/F");
  m_outputNTuple->tree()->Branch("m_lep_Tag_Eta",               &m_lep_Tag_Eta  	   , "m_lep_Tag_Eta/F");
  m_outputNTuple->tree()->Branch("m_lep_Tag_EtaBE2",            &m_lep_Tag_EtaBE2	   , "m_lep_Tag_EtaBE2/F");
  m_outputNTuple->tree()->Branch("m_lep_Tag_sigd0PV",           &m_lep_Tag_sigd0PV	   , "m_lep_Tag_sigd0PV/F");
  m_outputNTuple->tree()->Branch("m_lep_Tag_Z0SinTheta",        &m_lep_Tag_Z0SinTheta	   , "m_lep_Tag_Z0SinTheta/F");
  m_outputNTuple->tree()->Branch("m_lep_Tag_ID",               	&m_lep_Tag_ID		   , "m_lep_Tag_ID/F");
  m_outputNTuple->tree()->Branch("m_lep_Tag_isTrigMatch",       &m_lep_Tag_isTrigMatch     , "m_lep_Tag_isTrigMatch/B");
  m_outputNTuple->tree()->Branch("m_lep_Tag_isTightSelected",   &m_lep_Tag_isTightSelected , "m_lep_Tag_isTightSelected/B");
  
  m_outputNTuple->tree()->Branch("m_lep_Probe_Pt",                &m_lep_Probe_Pt	   , "m_lep_Probe_Pt/F");
  m_outputNTuple->tree()->Branch("m_lep_Probe_Eta",               &m_lep_Probe_Eta  	   , "m_lep_Probe_Eta/F");
  m_outputNTuple->tree()->Branch("m_lep_Probe_EtaBE2",            &m_lep_Probe_EtaBE2	   , "m_lep_Probe_EtaBE2/F");
  m_outputNTuple->tree()->Branch("m_lep_Probe_sigd0PV",           &m_lep_Probe_sigd0PV	   , "m_lep_Probe_sigd0PV/F");
  m_outputNTuple->tree()->Branch("m_lep_Probe_Z0SinTheta",        &m_lep_Probe_Z0SinTheta  , "m_lep_Probe_Z0SinTheta/F");
  m_outputNTuple->tree()->Branch("m_lep_Probe_ID",                &m_lep_Probe_ID	   , "m_lep_Probe_ID/F");
  m_outputNTuple->tree()->Branch("m_lep_Probe_isTrigMatch",       &m_lep_Probe_isTrigMatch     , "m_lep_Probe_isTrigMatch/B");
  m_outputNTuple->tree()->Branch("m_lep_Probe_isTightSelected",   &m_lep_Probe_isTightSelected , "m_lep_Probe_isTightSelected/B");  
  
  // ---------------------------------------------------------------------------------------------------------------
  
  m_numEvent = 0;
  
  Info("initialize()", "All good!");
  
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepMiniNTupMaker :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.
 
  //ANA_CHECK_SET_TYPE (EL::StatusCode);
  
  ++m_numEvent;
  
  if ( m_numEvent == 1 ) { Info("execute()", "Processing input TTree : %s\n", m_inputNTuple->GetName() ); }
 
  if ( m_debug ) { Info("execute()", "===> Event %u - EventNumber = %u ", static_cast<uint32_t>(m_numEvent), static_cast<uint32_t>(m_EventNumber) ); }
 
  m_inputNTuple->GetEntry (wk()->treeEntry());
  
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
    m_inputNTuple->SetBranchStatus (branch.c_str(), 1);
  }

  return EL::StatusCode::SUCCESS;

}
