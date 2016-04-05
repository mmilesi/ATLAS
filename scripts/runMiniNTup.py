import ROOT

import shutil, os

# load needed package inside the ROOT namespace
ROOT.gROOT.Macro('$ROOTCOREDIR/scripts/load_packages.C')

# create a new sample handler to access and find datasets
sh = ROOT.SH.SampleHandler()

indirname = "/coepp/cephfs/mel/mmilesi/ttH/multileptons_ntuple_run2/25ns_v6/Nominal"
ROOT.SH.scanDir(sh, indirname)

# specify the name of the input TTree
#sh.setMetaString(ROOT.SH.MetaFields.treeName, "nominal");
ROOT.SH.scanForTrees(sh);

# print the content of the sample handler
sh.printContent()

# create the job
job = ROOT.EL.Job()

# configure the job
job.sampleHandler(sh)
job.options().setDouble(ROOT.EL.Job.optMaxEvents, 1000)
job.options().setDouble (ROOT.EL.Job.optCacheSize, 10*1024*1024);  # set TTree cache size, needed for FAX

# instantiate and configure our algorithm
myalgo = ROOT.HTopMultilepMiniNTupMaker()
myalgo.m_name = "HTopMultilepMiniNTupMaker"
myalgo.m_debug = True
myalgo.m_outputNTupStreamName = "output"

# Set the output stream for the NTupleSvc,AlgSelect (must be done before adding the NTupleSvc,AlgSelect to the job)
outputstream = ROOT.EL.OutputStream(myalgo.m_outputNTupStreamName)
job.outputAdd(outputstream)

ntuplesvc = ROOT.EL.NTupleSvc(myalgo.m_outputNTupStreamName)
ntuplesvc.copyBranch("dilep_type")
ntuplesvc.copyBranch("lep_Pt_0")
ntuplesvc.copyBranch("lep_E_0")
ntuplesvc.copyBranch("lep_Eta_0")
ntuplesvc.copyBranch("lep_Phi_0")
ntuplesvc.copyBranch("lep_EtaBE2_0")
ntuplesvc.copyBranch("lep_Pt_1")

algskim = ROOT.EL.AlgSelect(myalgo.m_outputNTupStreamName)
algskim.addCut ("dilep_type>=1")
algskim.histName ("cut_flow")

# add the algorithm(s) to the job
job.algsAdd(ntuplesvc)
job.algsAdd(algskim)
job.algsAdd(myalgo)

# instantiate the driver, in this case we run locally
driver = ROOT.EL.DirectDriver()

# run the job through the driver
output_directory = "output_test"

if os.path.exists(output_directory):
    shutil.rmtree(output_directory)

driver.submit(job, output_directory)
