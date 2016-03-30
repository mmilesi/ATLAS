import ROOT

# load needed package inside the ROOT namespace
ROOT.gROOT.Macro('$ROOTCOREDIR/scripts/load_packages.C')

# create a new sample handler to access and find datasets
sh = ROOT.SH.SampleHandler()

indirname = "/coepp/cephfs/mel/mmilesi/ttH/multileptons_ntuple_run2/25ns_v6/Nominal/"
ROOT.SH.scanDir(sh, indirname)

# print the content of the sample handler
sh.printContent()

# create the job
job = ROOT.EL.Job()

# configure the job
job.sampleHandler(sh)
job.options().setDouble(ROOT.EL.Job.optMaxEvents, 1000) 
job.options().setDouble (ROOT.EL.Job.optCacheSize, 10*1024*1024);  # set TTree cache size, needed for FAX

# instantiate our algorithm
myAlg = ROOT.HTopMultilepMiniNTupMaker()

# configure algorithm
myAlg.m_name = "HTopMultilepMiniNTupMaker"
myAlg.m_debug = True
myAlg.m_outputNTupleName = "output"

output = ROOT.EL.OutputStream("output")    
job.outputAdd(output)

ntuplesvc = ROOT.EL.NTupleSvc("output");

ntuplesvc.copyBranch("lep_Pt_0")
ntuplesvc.copyBranch("lep_E_0")
ntuplesvc.copyBranch("lep_Eta_0")
ntuplesvc.copyBranch("lep_Phi_0")
ntuplesvc.copyBranch("lep_EtaBE2_0")

# add the algorithm(s) to the job
job.algsAdd(ntuplesvc)    
job.algsAdd(myAlg)

# instantiate the driver, in this case we run locally
driver = ROOT.EL.DirectDriver()

# run the job throught the driver
output_directory = "output_test"
driver.submit(job, output_directory)
