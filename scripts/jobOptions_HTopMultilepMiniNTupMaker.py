import ROOT
from xAH_config import xAH_config
import sys, os

sys.path.insert(0, os.environ['ROOTCOREBIN']+"/user_scripts/HTopMultilepAnalysis/")

c = xAH_config()

# List the branches to be copied over from the input TTree
#
eventweight_branches = ["mcWeightOrg","SherpaNJetWeight","pileupEventWeight_090",
                        "MV2c10_70_EventWeight","JVT_EventWeight",
			"lepSFTrigLoose","lepSFTrigTight",
			"lepSFObjLoose","lepSFObjTight",
			"tauSFTight","tauSFLoose"]

event_branches       = ["EventNumber","RunNumber","RunYear","mc_channel_number","passEventCleaning",
                        "onelep_type","dilep_type","trilep_type","quadlep_type","total_charge","total_leptons",
                        "isQMisIDEvent","isFakeEvent","isLepFromPhEvent",
			"nJets_OR","nJets_OR_MV2c10_70",
                        "nJets_OR_T","nJets_OR_T_MV2c10_70",
			"nTaus_OR_Pt25",
			"Mll01","Ptll01","DRll01","Mlll012","Mll02","Ptll02","DRll02","Mll12","Ptll12","DRll12","HT","HT_lep","HT_jets",
			"MET_RefFinal_et","MET_RefFinal_phi",
			"isBlinded"]

trigbits_branches    = ["HLT_.*"]

jet_branches         = ["lead_jetPt","lead_jetEta","lead_jetPhi","sublead_jetPt","sublead_jetEta","sublead_jetPhi"]

lep_branches         = ["lep_ID_0",
                        "lep_Index_0",
			"lep_Pt_0",
			"lep_E_0",
			"lep_Eta_0",
			"lep_Phi_0",
			"lep_EtaBE2_0",
			"lep_sigd0PV_0",
			"lep_Z0SinTheta_0",
                        "lep_isTightLH_0",
			"lep_isMediumLH_0",
			"lep_isLooseLH_0",
			"lep_isTight_0",
			"lep_isMedium_0",
			"lep_isLoose_0",
			"lep_isolationLooseTrackOnly_0",
			"lep_isolationLoose_0",
			"lep_isolationFixedCutTight_0",
			"lep_isolationFixedCutTightTrackOnly_0",
			"lep_isolationFixedCutLoose_0",
                        "lep_isTrigMatch_0",
			"lep_isPrompt_0",
			"lep_isBrems_0",
			"lep_isFakeLep_0",
			"lep_isQMisID_0",
			"lep_isConvPh_0",
			"lep_isISR_FSR_Ph_0",
                        "lep_isTruthMatched_0",
			"lep_truthType_0",
			"lep_truthOrigin_0",
			"lep_truthPdgId_0",
			"lep_truthStatus_0",
                        "lep_truthParentType_0",
			"lep_truthParentOrigin_0",
			"lep_truthParentPdgId_0",
			"lep_truthParentStatus_0",
                        "lep_truthPt_0",
			"lep_truthEta_0",
			"lep_truthPhi_0",
			"lep_truthM_0",
			"lep_truthE_0",
			"lep_truthRapidity_0",
                        "lep_SFIDLoose_0",
			"lep_SFIDTight_0",
			"lep_SFTrigLoose_0",
			"lep_SFTrigTight_0",
			"lep_EffTrigLoose_0",
			"lep_EffTrigTight_0",
			"lep_SFIsoLoose_0",
			"lep_SFIsoTight_0",
			"lep_SFReco_0",
			"lep_SFTTVA_0",
			"lep_SFObjLoose_0",
			"lep_SFObjTight_0",
			#
                        "lep_ID_1",
                        "lep_Index_1",
			"lep_Pt_1",
			"lep_E_1",
			"lep_Eta_1",
			"lep_Phi_1",
			"lep_EtaBE2_1",
			"lep_sigd0PV_1",
			"lep_Z0SinTheta_1",
                        "lep_isTightLH_1",
			"lep_isMediumLH_1",
			"lep_isLooseLH_1",
			"lep_isTight_1",
			"lep_isMedium_1",
			"lep_isLoose_1",
			"lep_isolationLooseTrackOnly_1",
			"lep_isolationLoose_1",
			"lep_isolationFixedCutTight_1",
			"lep_isolationFixedCutTightTrackOnly_1",
			"lep_isolationFixedCutLoose_1",
                        "lep_isTrigMatch_1",
			"lep_isPrompt_1",
			"lep_isBrems_1",
			"lep_isFakeLep_1",
			"lep_isQMisID_1",
			"lep_isConvPh_1",
			"lep_isISR_FSR_Ph_1",
                        "lep_isTruthMatched_1",
			"lep_truthType_1",
			"lep_truthOrigin_1",
			"lep_truthPdgId_1",
			"lep_truthStatus_1",
                        "lep_truthParentType_1",
			"lep_truthParentOrigin_1",
			"lep_truthParentPdgId_1",
			"lep_truthParentStatus_1",
                        "lep_truthPt_1",
			"lep_truthEta_1",
			"lep_truthPhi_1",
			"lep_truthM_1",
			"lep_truthE_1",
			"lep_truthRapidity_1",
                        "lep_SFIDLoose_1",
			"lep_SFIDTight_1",
			"lep_SFTrigLoose_1",
			"lep_SFTrigTight_1",
			"lep_EffTrigLoose_1",
			"lep_EffTrigTight_1",
			"lep_SFIsoLoose_1",
			"lep_SFIsoTight_1",
			"lep_SFReco_1",
			"lep_SFTTVA_1",
			"lep_SFObjLoose_1",
			"lep_SFObjTight_1",			
			#
                        "lep_ID_2",
                        "lep_Index_2",
			"lep_Pt_2",
			"lep_E_2",
			"lep_Eta_2",
			"lep_Phi_2",
			"lep_EtaBE2_2",
			"lep_sigd0PV_2",
			"lep_Z0SinTheta_2",
                        "lep_isTightLH_2",
			"lep_isMediumLH_2",
			"lep_isLooseLH_2",
			"lep_isTight_2",
			"lep_isMedium_2",
			"lep_isLoose_2",
			"lep_isolationLooseTrackOnly_2",
			"lep_isolationLoose_2",
			"lep_isolationFixedCutTight_2",
			"lep_isolationFixedCutTightTrackOnly_2",
			"lep_isolationFixedCutLoose_2",
                        "lep_isTrigMatch_2",
			"lep_isPrompt_2",
			"lep_isBrems_2",
			"lep_isFakeLep_2",
			"lep_isQMisID_2",
			"lep_isConvPh_2",
			"lep_isISR_FSR_Ph_2",
                        "lep_isTruthMatched_2",
			"lep_truthType_2",
			"lep_truthOrigin_2",
			"lep_truthPdgId_2",
			"lep_truthStatus_2",
                        "lep_truthParentType_2",
			"lep_truthParentOrigin_2",
			"lep_truthParentPdgId_2",
			"lep_truthParentStatus_2",
                        "lep_truthPt_2",
			"lep_truthEta_2",
			"lep_truthPhi_2",
			"lep_truthM_2",
			"lep_truthE_2",
			"lep_truthRapidity_2",
                        "lep_SFIDLoose_2",
			"lep_SFIDTight_2",
			"lep_SFTrigLoose_2",
			"lep_SFTrigTight_2",
			"lep_EffTrigLoose_2",
			"lep_EffTrigTight_2",
			"lep_SFIsoLoose_2",
			"lep_SFIsoTight_2",
			"lep_SFReco_2",
			"lep_SFTTVA_2",
			"lep_SFObjLoose_2",
			"lep_SFObjTight_2",			
                      ]
tau_branches        = ["tau_pt_0","tau_eta_0","tau_phi_0","tau_charge_0","tau_BDTJetScore_0","tau_JetBDTSigLoose_0","tau_JetBDTSigMedium_0","tau_JetBDTSigTight_0","tau_numTrack_0","tau_SFTight_0","tau_SFLoose_0"]

MET_truth_branches  = ["MET_Truth_px","MET_Truth_py","MET_Truth_phi","MET_Truth_sumet"]

mc_truth_branches   = [] # ["m_truth_m","m_truth_pt","m_truth_eta","m_truth_phi","m_truth_e","m_truth_pdgId","m_truth_status","m_truth_barcode","m_truth_parents","m_truth_children"]

branches_to_copy = eventweight_branches + event_branches + trigbits_branches + jet_branches + lep_branches + tau_branches + mc_truth_branches + MET_truth_branches

# ---------------------------
#
# Add here branches that need to be used (hence activated), but not need to be copied over
#
jet_vec_branches       = ["selected_jets","selected_jets_T","m_jet_pt","m_jet_eta","m_jet_phi","m_jet_E","m_jet_flavor_truth_label","m_jet_flavor_truth_label_ghost"]
jet_truth_vec_branches = ["m_truth_jet_pt","m_truth_jet_eta","m_truth_jet_phi","m_truth_jet_e"]

branches_to_activate = branches_to_copy + jet_vec_branches + jet_truth_vec_branches

# Trick to pass the list as a comma-separated string to the C++ algorithm
#
branches_to_activate_str = ",".join(branches_to_activate)

print("Activating branches from input TTree:")
for branch in branches_to_activate:
   print("\t{0}".format(branch))

# Instantiate the main algorithm
#
HTopMultilepMiniNTupMakerDict = { "m_name"                 : "HTopMultilepMiniNTupMaker",
                                  "m_debug"                : False,
				  "m_outputNTupName"       : "physics",
                                  "m_outputNTupStreamName" : "output",
				  "m_inputBranches"        : branches_to_activate_str,
	                          "m_useAlgSelect"         : True,
				  "m_addStreamEventsHist"  : False,
                                }

# Instantiate the NTupleSvc algorithm
#
ntuplesvc = ROOT.EL.NTupleSvc(HTopMultilepMiniNTupMakerDict["m_outputNTupStreamName"])

# Set the branches to be copied over from the input TTree
#
print("Copying branches from input TTree to output:")
for branch in branches_to_copy:
   print("\t{0}".format(branch))
   ntuplesvc.copyBranch(branch)

# Instantiate the AlgSelect algorithm to skim the input ntuple
#
algskim = ROOT.EL.AlgSelect(HTopMultilepMiniNTupMakerDict["m_outputNTupStreamName"])
algskim.addCut("passEventCleaning==1")
algskim.addCut("dilep_type>0||trilep_type>0")
algskim.addCut("(dilep_type>0&&nJets_OR_T>=2&&nJets_OR_T_MV2c10_70>=1)||(trilep_type>0&&nJets_OR>=2&&nJets_OR_MV2c10_70>=1)")
algskim.histName("cutflow")

# Add the algorithms to the job.
#
# Here order matters!
#
c._algorithms.append(ntuplesvc)
if ( HTopMultilepMiniNTupMakerDict["m_useAlgSelect"] ):
  c._algorithms.append(algskim)
c.setalg("HTopMultilepMiniNTupMaker", HTopMultilepMiniNTupMakerDict)


