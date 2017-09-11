
import ROOT
from xAH_config import xAH_config
import sys, os

sys.path.insert(0, os.environ['ROOTCOREBIN']+"/user_scripts/HTopMultilepAnalysis/")

c = xAH_config()

# List the branches to be copied over from the input TTree

eventweight_branches = ["mcWeightOrg","SherpaNJetWeight","pileupEventWeight_090",
                        "bTagSF_weight_Continuous","MV2c10_Continuous_EventWeight",
                        "MV2c10_70_EventWeight","MV2c10_77_EventWeight","JVT_EventWeight",
                        "lepSFTrigLoose","lepSFTrigTight",
                        "lepSFTrigTightLoose","lepSFTrigLooseTight",
                        "lepEffTrigLoose","lepEffTrigTight",
                        "lepEffTrigTightLoose","lepEffTrigLooseTight",
                        "lepDataEffTrigLoose","lepDataEffTrigTight",
                        "lepDataEffTrigTightLoose","lepDataEffTrigLooseTight",
                        "lepSFObjLoose","lepSFObjTight",
			"tauSFTight","tauSFLoose",
                        "m_hasMEphoton.*",
                        "m_MEphoton_.*",
                        ]

event_branches       = ["EventNumber","RunNumber","RunYear","mc_channel_number","averageIntPerXing","passEventCleaning",
                        "onelep_type","dilep_type","trilep_type","quadlep_type","total_charge","total_leptons",
                        "isQMisIDEvent","isFakeEvent","isLepFromPhEvent",
			"nJets_OR","nJets_OR_MV2c10_70","nJets_OR_MV2c10_77",
                        "nJets_OR_T","nJets_OR_T_MV2c10_70","nJets_OR_T_MV2c10_77",
			"nTaus_OR_Pt25",
			"Mll01","Ptll01","DRll01","Mlll012","Mll02","Ptll02","DRll02","Mll12","Ptll12","DRll12","HT","HT_lep","HT_jets",
                        "matchDLTll.*",
			"MET_RefFinal_et","MET_RefFinal_phi",
			"isBlinded"]

trigbits_branches    = ["HLT.*[^(_PS)]"] # No need to store the prescales!

trigmatch_branches   = [
        		# 2015
        		"electron_match_HLT_e24_lhmedium_L1EM20VH",
        		"electron_match_HLT_e60_lhmedium",
        		"electron_match_HLT_e120_lhloose",
        		"electron_match_HLT_2e12_lhloose_L12EM10VH",
        		"electron_match_HLT_e24_medium_L1EM20VHI_mu8noL1",
        		"electron_match_HLT_e7_medium_mu24",
        		"muon_match_HLT_mu20_iloose_L1MU15",
        		"muon_match_HLT_mu18_mu8noL1",
			"muon_match_HLT_e24_medium_L1EM20VHI_mu8noL1",
			"muon_match_HLT_e7_medium_mu24",
        		# 2016
        		"electron_match_HLT_e26_lhtight_nod0_ivarloose",
        		"electron_match_HLT_e60_lhmedium_nod0",
        		"electron_match_HLT_e140_lhloose_nod0",
        		"electron_match_HLT_2e17_lhvloose_nod0",
			"electron_match_HLT_e17_lhloose_mu14",
        		"electron_match_HLT_e17_lhloose_nod0_mu14",
			"electron_match_HLT_e26_lhmedium_nod0_L1EM22VHI_mu8noL1",
        		"electron_match_HLT_e7_lhmedium_mu24",
        		"muon_match_HLT_mu26_ivarmedium",
        		"muon_match_HLT_mu22_mu8noL1",
			"muon_match_HLT_e17_lhloose_mu14",
			"muon_match_HLT_e17_lhloose_nod0_mu14",
			"muon_match_HLT_e26_lhmedium_nod0_L1EM22VHI_mu8noL1",
			"muon_match_HLT_e7_lhmedium_nod0_mu24",
        		 # 2015 & 2016
        		"muon_match_HLT_mu50",
                       ]

vec_lep_branches     = ["electron_passOR",
                        "electron_pt",
                        "electron_eta",
                        "electron_EtaBE2",
                        "electron_phi",
                        "electron_E",
                        "electron_ID",
                        "electron_sigd0PV",
                        "electron_z0SinTheta",
                        "electron_topoetcone20",
                        "electron_ptvarcone20",
                        "electron_truthType",
                        "electron_truthOrigin",
                        #
                        "muon_passOR",
                        "muon_pt",
                        "muon_eta",
                        "muon_phi",
                        "muon_ID",
                        "muon_sigd0PV",
                        "muon_z0SinTheta",
                        "muon_ptvarcone30",
                        "muon_truthType",
                        "muon_truthOrigin",
                       ]

base_lep_branches = [
   "lep_ID_",
   "lep_Index_",
   "lep_Pt_",
   "lep_E_",
   "lep_Eta_",
   "lep_Phi_",
   "lep_EtaBE2_",
   "lep_topoEtcone20_",
   "lep_topoEtcone30_",
   "lep_topoEtcone40_",
   "lep_ptVarcone20_",
   "lep_ptVarcone30_",
   "lep_ptVarcone40_",
   "lep_sigd0PV_",
   "lep_Z0SinTheta_",
   "lep_d0_",
   "lep_z0_",
   "lep_vz_",
   "lep_deltaz0_",
   "lep_isTightLH_",
   "lep_isMediumLH_",
   "lep_isLooseLH_",
   "lep_isTight_",
   "lep_isMedium_",
   "lep_isLoose_",
   "lep_isolationLooseTrackOnly_",
   "lep_isolationLoose_",
   "lep_isolationGradient_",
   "lep_isolationGradientLoose_",
   "lep_isolationFixedCutTight_",
   "lep_isolationFixedCutTightTrackOnly_",
   "lep_isolationFixedCutLoose_",
   "lep_isTrigMatch_",
   "lep_isTrigMatchDLT_",
   "lep_isPrompt_",
   "lep_isFakeLep_",
   "lep_isQMisID_",
   "lep_isConvPh_",
   "lep_isISR_FSR_Ph_",
   "lep_isBrems_",
   "lep_chargeIDBDTLoose_",
   "lep_chargeIDBDTMedium_",
   "lep_chargeIDBDTTight_",
   "lep_promptLeptonIso_TagWeight_",
   "lep_promptLeptonIso_sv1_jf_ntrkv_",
   "lep_promptLeptonIso_TrackJetNTrack_",
   "lep_promptLeptonIso_ip2_",
   "lep_promptLeptonIso_ip3_",
   "lep_promptLeptonIso_DRlj_",
   "lep_promptLeptonIso_LepJetPtFrac_",
   "lep_promptLepton_TagWeight_",
   "lep_promptLeptonNoIso_TagWeight_",
   "lep_jet_Pt_",
   "lep_jet_nTrk_",
   "lep_jet_sumPtTrk_",
   "lep_jet_mv2c10_",
   "lep_jet_deltaR_",
   "lep_jet_ptRel_",
   "lep_jet_ptOverMuPt_",
   "lep_jet_BDT_",
   "lep_isTruthMatched_",
   "lep_truthType_",
   "lep_truthOrigin_",
   "lep_truthPdgId_",
   "lep_truthStatus_",
   "lep_truthParentType_",
   "lep_truthParentOrigin_",
   "lep_truthParentPdgId_",
   "lep_truthParentStatus_",
   "lep_truthPt_",
   "lep_truthEta_",
   "lep_truthPhi_",
   "lep_truthM_",
   "lep_truthE_",
   "lep_truthRapidity_",
   "lep_SFIDLoose_",
   "lep_SFIDTight_",
   "lep_SFTrigLoose_",
   "lep_SFTrigTight_",
   "lep_EffTrigLoose_",
   "lep_EffTrigTight_",
   "lep_SFIsoLoose_",
   "lep_SFIsoTight_",
   "lep_SFReco_",
   "lep_SFTTVA_",
   "lep_SFObjLoose_",
   "lep_SFObjTight_",
   "lep_nInnerPix_",
   ]

lep_branches = []

for b in base_lep_branches:
   for idx in range(0,3):
      this_b = b + str(idx)
      lep_branches.append(this_b)

# lep_branches        = ["lep_.*[^4]"]

tau_branches        = ["tau_pt_0","tau_eta_0","tau_phi_0","tau_charge_0","tau_BDTJetScore_0","tau_JetBDTSigLoose_0","tau_JetBDTSigMedium_0","tau_JetBDTSigTight_0","tau_numTrack_0","tau_SFTight_0","tau_SFLoose_0"]

MET_truth_branches  = ["MET_Truth_px","MET_Truth_py","MET_Truth_phi","MET_Truth_sumet"]

mc_truth_branches   = [] # ["m_truth_m","m_truth_pt","m_truth_eta","m_truth_phi","m_truth_e","m_truth_pdgId","m_truth_status","m_truth_barcode","m_truth_parents","m_truth_children"]

branches_to_copy = eventweight_branches + event_branches + trigbits_branches + lep_branches + tau_branches + mc_truth_branches + MET_truth_branches

# ---------------------------

# Add here branches that need to be used (hence activated), but not need to be copied over

jet_vec_branches       = ["selected_jets","selected_jets_T","m_jet_pt","m_jet_eta","m_jet_phi","m_jet_E","m_jet_flavor_weight_MV2c10","m_jet_flavor_truth_label","m_jet_flavor_truth_label_ghost"]
jet_truth_vec_branches = ["m_truth_jet_pt","m_truth_jet_eta","m_truth_jet_phi","m_truth_jet_e"]

prompt_lep_branches = ["electron_PromptLeptonIso_TagWeight",
                       "electron_ChargeIDBDTLoose",
                       "electron_ChargeIDBDTMedium",
                       "electron_ChargeIDBDTTight",
                       "muon_PromptLeptonIso_TagWeight",
                       ]

branches_to_activate = branches_to_copy + jet_vec_branches + jet_truth_vec_branches + vec_lep_branches + prompt_lep_branches # + trigmatch_branches

# Trick to pass the list as a comma-separated string to the C++ algorithm

branches_to_activate_str = ",".join(branches_to_activate)

print("\nActivating branches from input TTree:\n")
print( ",".join("{0}".format(branch) for branch in branches_to_activate) + "\n" )

# Instantiate the main algorithm

HTopMultilepMiniNTupMakerDict = { "m_name"                 : "HTopMultilepMiniNTupMaker",
                                  "m_debug"                : False,
				  "m_outputNTupName"       : "physics",
                                  "m_outputNTupStreamName" : "output",
				  "m_inputBranches"        : branches_to_activate_str,
	                          "m_useAlgSelect"         : True,
				  "m_addStreamEventsHist"  : False,
				  "m_useTruthTP"           : True,
				  "m_useNominalTP"         : False,
                                  "m_do3LTP"               : False,
                                  "m_ambiSolvingCrit"      : "OF",
                                  #"m_ambiSolvingCrit"     : "Pt",
                                  #"m_ambiSolvingCrit"     : "deltaRClosestBJet",
                                  #"m_ambiSolvingCrit"     : "massClosestBJet",
                                  "m_lepSelForTP"          : "MVA",
                                  #"m_lepSelForTP"          : "CutBased",
                                  "m_jetTruthMatching"     : False,
                                }

# Instantiate the NTupleSvc algorithm

ntuplesvc = ROOT.EL.NTupleSvc(HTopMultilepMiniNTupMakerDict["m_outputNTupStreamName"])

# Set the branches to be copied over from the input TTree

print("Copying branches from input TTree to output:\n")
print( ",".join("{0}".format(branch) for branch in branches_to_copy) + "\n" )
for branch in branches_to_copy:
   ntuplesvc.copyBranch(branch)

# Instantiate the AlgSelect algorithm to skim the input ntuple

jet_skim_cut_2L = "( dilep_type > 0 && ( ( nJets_OR_T >= 2 && nJets_OR_T_MV2c10_70 >= 1 ) || ( nJets_OR_T >= 4 && nJets_OR_T_MV2c10_70 == 0 ) ) )"
jet_skim_cut_3L = "( trilep_type > 0 && nJets_OR_T >= 2 && nJets_OR_T_MV2c10_70 >= 1 )"
jet_skim_cut = jet_skim_cut_2L + " || " + jet_skim_cut_3L

tau_skim_cut = "nTaus_OR_Pt25 == 0"

algskim = ROOT.EL.AlgSelect(HTopMultilepMiniNTupMakerDict["m_outputNTupStreamName"])
algskim.addCut("RunYear==2015 || RunYear==2016")
algskim.addCut("passEventCleaning == 1")
algskim.addCut("( dilep_type > 0 && lep_Pt_1 > 10e3 ) || ( trilep_type > 0 && lep_Pt_0 > 10e3 && lep_Pt_1 > 10e3 && lep_Pt_2 > 10e3 )") # keep only dilepton,trilepton events, ensuring min(pT lep) > 10 GeV (NB: in 3L, flat branches are not pT-ranked, so we have to be explicit...)
algskim.addCut(jet_skim_cut)
algskim.addCut(tau_skim_cut)

if ( HTopMultilepMiniNTupMakerDict["m_do3LTP"] ):
   algskim.addCut("dilep_type > 0 && ( nJets_OR_T >= 2 && nJets_OR_T<= 3 && nJets_OR_T_MV2c10_70 >= 1 )") # Reduce the size as much as possible when doing 3L T&P
# algskim.addCut("trilep_type > 0 && lep_Pt_0 > 10e3 && lep_Pt_1 > 10e3 && lep_Pt_2 > 10e3") # Use this skimming cut only for Z-->mm+e T&P
algskim.histName("cutflow")

# Add the algorithms to the job.
# Here order matters!

c._algorithms.append(ntuplesvc)
if ( HTopMultilepMiniNTupMakerDict["m_useAlgSelect"] ):
  c._algorithms.append(algskim)
c.setalg("HTopMultilepMiniNTupMaker", HTopMultilepMiniNTupMakerDict)


