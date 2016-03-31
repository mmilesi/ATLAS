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

# List the branches to be copied over from the input TTree
eventweight_branches = ["mcWeightOrg","pileupEventWeight_090","MV2c20_77_EventWeight","JVT_EventWeight"]
event_branches       = ["EventNumber","RunNumber",
                        "onelep_type","dilep_type","trilep_type","quadlep_type","total_charge","total_leptons",
                        "Mll01","HT","HT_lep","HT_jets"
]
trigbits_branches    = ["HLT_e24_lhmedium_L1EM20VH","HLT_e24_lhmedium_L1EM18VH","HLT_e60_lhmedium","HLT_e120_lhloose",
                        "HLT_mu20_iloose_L1MU15","HLT_mu50"]

jet_branches   = ["lead_jetPt","lead_jetEta","lead_jetPhi",
                  "sublead_jetPt","sublead_jetEta","sublead_jetPhi"]



lep_branches   = ["lep_ID_0","lep_Index_0","lep_Pt_0","lep_E_0","lep_Eta_0","lep_Phi_0","lep_EtaBE2_0","lep_sigd0PV_0","lep_Z0SinTheta_0","lep_isTrigMatch_0","lep_isPrompt_0","lep_isBremsElec_0","lep_isFakeLep_0",
                  "lep_ID_1","lep_Index_1","lep_Pt_1","lep_E_1","lep_Eta_1","lep_Phi_1","lep_EtaBE2_1","lep_sigd0PV_1","lep_Z0SinTheta_1","lep_isTrigMatch_1","lep_isPrompt_1","lep_isBremsElec_1","lep_isFakeLep_1",
                  "lep_ID_2","lep_Index_2","lep_Pt_2","lep_E_2","lep_Eta_2","lep_Phi_2","lep_EtaBE2_2","lep_sigd0PV_2","lep_Z0SinTheta_2","lep_isTrigMatch_2","lep_isPrompt_2","lep_isBremsElec_2","lep_isFakeLep_2",
                  "lep_ID_3","lep_Index_3","lep_Pt_3","lep_E_3","lep_Eta_3","lep_Phi_3","lep_EtaBE2_3","lep_sigd0PV_3","lep_Z0SinTheta_3","lep_isTrigMatch_3","lep_isPrompt_3","lep_isBremsElec_3","lep_isFakeLep_3",
                  ]

branches_to_copy = eventweight_branches + event_branches + trigbits_branches + jet_branches + lep_branches

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


