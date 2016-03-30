#include <HTopMultilepAnalysis/HTopMultilepMiniNTupMaker.h>

// this is needed to distribute the algorithm to the workers
ClassImp(HTopMultilepMiniNTupMaker)

HTopMultilepMiniNTupMaker :: HTopMultilepMiniNTupMaker(std::string className) :
    Algorithm(className)//,
    //m_inputNTuple(nullptr),
    //m_outputNTuple(nullptr)
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().

  Info("HTopMultilepMiniNTupMaker()", "Calling constructor");
  
  m_outputNTupleName = "output";
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
  
  Info("setupJob()", "Calling setupJob");

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepMiniNTupMaker :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.
  
  Info("histInitialize()", "Calling histInitialize");
  
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepMiniNTupMaker :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  
  Info("fileExecute()", "Calling fileExecute");

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepMiniNTupMaker :: changeInput (bool firstFile)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.

  Info("changeInput()", "Calling changeInput");
  
  m_inputNTuple = wk()->tree();
  
  // Firstly, disable all branches
  //
  m_inputNTuple->SetBranchStatus ("*", 0);

  // Re-enable the branches we are going to use
  //
  m_inputNTuple->SetBranchStatus ("EventNumber", 1);  
  m_inputNTuple->SetBranchStatus ("lep_Pt_0", 1);
  m_inputNTuple->SetBranchStatus ("lep_E_0", 1);
  m_inputNTuple->SetBranchStatus ("lep_Eta_0", 1);
  m_inputNTuple->SetBranchStatus ("lep_Phi_0", 1);
  m_inputNTuple->SetBranchStatus ("lep_EtaBE2_0", 1);
  
  // Connect the branches of the input tree to the algorithm members
  //
  m_inputNTuple->SetBranchAddress ("EventNumber",  &m_EventNumber);
  m_inputNTuple->SetBranchAddress ("lep_Pt_0",     &m_lep_Pt_0);
  m_inputNTuple->SetBranchAddress ("lep_E_0",      &m_lep_E_0);
  m_inputNTuple->SetBranchAddress ("lep_Eta_0",    &m_lep_Eta_0);
  m_inputNTuple->SetBranchAddress ("lep_Phi_0",    &m_lep_Phi_0);
  m_inputNTuple->SetBranchAddress ("lep_EtaBE2_0", &m_lep_EtaBE2_0);

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

  Info("initialize()", "Initialising HTopMultilepMiniNTupMaker...");
  
  m_outputNTuple = EL::getNTupleSvc (wk(), m_outputNTupleName);
 
  m_outputNTuple->tree()->Branch("lep_Pt_0_Squared",  &m_lep_Pt_0_Squared, "lep_Pt_0_Squared/F");
  
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
 
  ++m_numEvent;
  
  if ( m_numEvent == 1 ) { Info("execute()", "Processing input TTree : %s\n", m_inputNTuple->GetName() ); }
 
  if ( m_debug ) { Info("execute()", "===> Event %u - EventNumber = %u ", static_cast<uint32_t>(m_numEvent), static_cast<uint32_t>(m_EventNumber) ); }
 
  m_inputNTuple->GetEntry (wk()->treeEntry());
  
  m_lep_Pt_0_Squared = m_lep_Pt_0 * m_lep_Pt_0;
  
  if ( m_debug ) { Info("execute()", "\t lep_Pt_0_Squared = %.2f", m_lep_Pt_0_Squared ); }
  
  
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HTopMultilepMiniNTupMaker :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
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
  return EL::StatusCode::SUCCESS;
}
