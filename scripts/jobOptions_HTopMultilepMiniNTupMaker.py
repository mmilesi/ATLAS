
import ROOT
from xAH_config import xAH_config
import sys, os

sys.path.insert(0, os.environ['ROOTCOREBIN']+"/user_scripts/HTopMultilepAnalysis/")

c = xAH_config()

# List the branches to be copied over from the input TTree

eventweight_branches = ["mcWeightOrg","SherpaNJetWeight","pileupEventWeight_090",
                        "MV2c10_70_EventWeight","MV2c10_77_EventWeight","JVT_EventWeight",
			"lepSFTrigLoose","lepSFTrigTight",
			"lepSFObjLoose","lepSFObjTight",
			"tauSFTight","tauSFLoose"]

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

lep_branches        = ["lep_.*[^4]"]

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
				  "m_useTruthTP"           : False,
				  "m_useNominalTP"         : True,
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
jet_skim_cut_3L = "( trilep_type > 0 && nJets_OR >= 2 )"
jet_skim_cut_4L = "( quadlep_type > 0 && nJets_OR >= 2 )"
jet_skim_cut = jet_skim_cut_2L + " || " + jet_skim_cut_3L + " || " + jet_skim_cut_4L

algskim = ROOT.EL.AlgSelect(HTopMultilepMiniNTupMakerDict["m_outputNTupStreamName"])
algskim.addCut("RunYear==2015 || RunYear==2016")
algskim.addCut("passEventCleaning == 1")
algskim.addCut("( dilep_type > 0 && lep_Pt_1 > 10e3 ) || ( trilep_type > 0 && lep_Pt_0 > 10e3 && lep_Pt_1 > 10e3 && lep_Pt_2 > 10e3 ) || ( quadlep_type > 0 && lep_Pt_3 > 10e3 )") # keep only dilepton,trilepton and quadlepton events, ensuring min(pT lep) > 10 GeV (NB: in 3L, flat branches are not pT-ranked, so we have to be explicit...)
algskim.addCut(jet_skim_cut)
algskim.histName("cutflow")

# Add the algorithms to the job.
# Here order matters!

c._algorithms.append(ntuplesvc)
if ( HTopMultilepMiniNTupMakerDict["m_useAlgSelect"] ):
  c._algorithms.append(algskim)
c.setalg("HTopMultilepMiniNTupMaker", HTopMultilepMiniNTupMakerDict)


