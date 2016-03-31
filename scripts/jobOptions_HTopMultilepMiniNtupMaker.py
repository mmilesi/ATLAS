import ROOT
from xAH_config import xAH_config
import sys, os

sys.path.insert(0, os.environ['ROOTCOREBIN']+"/user_scripts/HTopMultilepAnalysis/")

c = xAH_config()

HTopMultilepMiniNTupMakerDict = { "m_name" : "HTopMultilepMiniNTupMaker",
                                  "m_debug" : False,
                                  "m_outputNTupStreamName" : "output",
                                }

# Instantiate the NTupleSvc algorithm
ntuplesvc = ROOT.EL.NTupleSvc(HTopMultilepMiniNTupMakerDict["m_outputNTupStreamName"])

event_branches = ["mcWeightOrg","dilep_type"]
lep_branches   = ["lep_Pt_0","lep_E_0","lep_Eta_0","lep_Phi_0","lep_EtaBE2_0","lep_Pt_1"]

branches_to_copy = event_branches + lep_branches

# Set the branches to be copied over from the input TTree
for branch in branches_to_copy:
   ntuplesvc.copyBranch(branch)

# Instantiate the AlgSelect algorithm to skim the input ntuple
algskim = ROOT.EL.AlgSelect(HTopMultilepMiniNTupMakerDict["m_outputNTupStreamName"])
algskim.addCut ("dilep_type>=1")
algskim.histName ("cut_flow")

# Add the algorithms to the job.
#
# Here order matters! 
#
c._algorithms.append(ntuplesvc)
c._algorithms.append(algskim)
c.setalg("HTopMultilepMiniNTupMaker", HTopMultilepMiniNTupMakerDict)


