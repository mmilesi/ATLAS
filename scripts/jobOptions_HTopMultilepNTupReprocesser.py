import ROOT
from xAH_config import xAH_config
import sys, os

sys.path.insert(0, os.environ['ROOTCOREBIN']+"/user_scripts/HTopMultilepAnalysis/")

c = xAH_config()

event_branches = ["EventNumber","RunNumber","mc_channel_number","isSS01","dilep_type","QMisIDWeight","MMWeight"]
lep_branches   = ["lep_ID_0","lep_Pt_0","lep_Eta_0","lep_Phi_0","lep_EtaBE2_0","lep_isTightSelected_0",
                  "lep_ID_1","lep_Pt_1","lep_Eta_1","lep_Phi_1","lep_EtaBE2_1","lep_isTightSelected_1"]

branches_to_activate = event_branches + lep_branches

# Trick to pass the list as a comma-separated string to the C++ algorithm
#
branches_to_activate_str = ",".join(branches_to_activate)

# Instantiate the main algorithm
#
HTopMultilepNTupReprocesserDict = { "m_name"                       : "HTopMultilepNTupReprocesser",
                                    "m_debug"                      : True,
                                    "m_verbose"                    : True,
                                    "m_inputBranches"              : branches_to_activate_str,
                                    "m_weightToCalc"               : "QMisID",
                                    "m_QMisIDRates_dir"            : "$ROOTCOREBIN/data/HTopMultilepAnalysis/External/",
                                    "m_QMisIDRates_Filename_T"     : "QMisIDRates_Data_2016_T_25ns_v19.root",
                                    "m_QMisIDRates_Filename_AntiT" : "QMisIDRates_Data_2016_TanitiT_25ns_v19.root",
                                    "m_useTAntiTRates"             : True,
                                    "m_RR_dir"                     : "/afs/cern.ch/user/m/mmilesi/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMRates_25ns_v19_DDQMisID",
                                    "m_FR_dir"                     : "/afs/cern.ch/user/m/mmilesi/ttH/RUN2/HTopMultilepAnalysisCode/trunk/HTopMultilepAnalysis/PlotUtils/OutputPlots_MMRates_25ns_v19_DDQMisID",
                                    "m_doMMClosure"                : False,
                                    "m_useEtaParemetrisation"      : False,
                                  }

# Instantiate the NTupleSvc algorithm
#
ntuplesvc = ROOT.EL.NTupleSvc(HTopMultilepNTupReprocesserDict["m_outputNTupStreamName"])

# Copy ALL branches over from the input TTree
#
print("Copying branches from input TTree to output:")
ntuplesvc.copyBranch(".*")

# Add the algorithms to the job.
#
# Here order matters!
#
c._algorithms.append(ntuplesvc)
c.setalg("HTopMultilepNTupReprocesser", HTopMultilepNTupReprocesserDict)


