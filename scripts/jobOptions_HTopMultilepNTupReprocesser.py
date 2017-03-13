import ROOT
from xAH_config import xAH_config
import sys, os

sys.path.insert(0, os.environ['ROOTCOREBIN']+"/user_scripts/HTopMultilepAnalysis/")

c = xAH_config()

event_branches = ["EventNumber","RunNumber","mc_channel_number","isSS01","dilep_type","trilep_type",
                  "is_T_T","is_T_AntiT","is_AntiT_T","is_AntiT_AntiT",
                  "is_TMVA_TMVA","is_TMVA_AntiTMVA","is_AntiTMVA_TMVA","is_AntiTMVA_AntiTMVA",
                  ]

lep_branches   = ["lep_ID_0","lep_Pt_0","lep_Eta_0","lep_Phi_0","lep_EtaBE2_0","lep_isTightSelected_0","lep_isTightSelectedMVA_0","lep_isTrigMatch_0",
                  "lep_ID_1","lep_Pt_1","lep_Eta_1","lep_Phi_1","lep_EtaBE2_1","lep_isTightSelected_1","lep_isTightSelectedMVA_1","lep_isTrigMatch_1"]

branches_to_activate = event_branches + lep_branches

# Trick to pass the list as a comma-separated string to the C++ algorithm

branches_to_activate_str = ",".join(branches_to_activate)

# Instantiate the main algorithm

base_dir = "/imports/home/mmilesi/PhD/ttH_MultiLeptons/RUN2/HTopMultilepAnalysisCode/trunk/"
#base_dir = "/afs/cern.ch/user/m/mmilesi/ttH/RUN2/HTopMultilepAnalysisCode/trunk"

HTopMultilepNTupReprocesserDict = { "m_name"                       : "HTopMultilepNTupReprocesser",
                                    "m_debug"                      : False,
                                    "m_verbose"                    : False,
				    "m_outputNTupStreamName"       : "output",
                                    "m_inputBranches"              : branches_to_activate_str,
                                    "m_weightToCalc"               : "MM", #"QMisID",
                                    "m_QMisIDRates_dir"            : "$ROOTCOREBIN/data/HTopMultilepAnalysis/External/",
                                    "m_QMisIDRates_Filename_T"     : "Rates_3D_Tot_v26_Tight.root",
                                    "m_QMisIDRates_Filename_AntiT" : "Rates_3D_Tot_v26_Loose.root",
                                    "m_QMisIDRates_Histname_T"     : "LikelihoodEtaPtTight",
                                    "m_QMisIDRates_Histname_AntiT" : "LikelihoodEtaPtLoose",
                                    "m_useTAntiTRates"             : False, # --> set this option to True if QMisID rates have NOT been measured independently for T and AntiT electrons - Set True for v19
                                    #
                                    # ------------------------------------------------------------
                                    #
                                    # v24
                                    #
                                    #"m_REFF_dir" 		   : base_dir + "HTopMultilepAnalysis/PlotUtils/PLOTS_25ns_v24/MMClosure_v24_TruthTP/OutputPlots_MMClosureRates_TruthTagProbe_NoCorr_SLT_RealOFmuOFel_FakeSFmuOFel_25ns_v24", # Truth-based efficiencies
                                    #
                                    #"m_REFF_dir" 		   : base_dir + "HTopMultilepAnalysis/PlotUtils/PLOTS_25ns_v24/CombinedEfficiencies",
                                    #"m_EFF_YES_TM_dir" 	   : base_dir + "HTopMultilepAnalysis/PlotUtils/PLOTS_25ns_v24/CombinedEfficiencies_TRIGMATCH_EFF",
                                    #"m_EFF_NO_TM_dir" 		   : base_dir + "HTopMultilepAnalysis/PlotUtils/PLOTS_25ns_v24/CombinedEfficiencies_NOT_TRIGMATCH_EFF",
                                    #
                                    # Removed isolation on baseline electrons
                                    #
                                    #"m_REFF_dir" 		   : base_dir + "HTopMultilepAnalysis/PlotUtils/PLOTS_25ns_v24_ElNoIso/CombinedEfficiencies",
                                    #"m_EFF_YES_TM_dir" 	   : base_dir + "HTopMultilepAnalysis/PlotUtils/PLOTS_25ns_v24_ElNoIso/CombinedEfficiencies_TRIGMATCH_EFF",
                                    #"m_EFF_NO_TM_dir" 		   : base_dir + "HTopMultilepAnalysis/PlotUtils/PLOTS_25ns_v24_ElNoIso/CombinedEfficiencies_NOT_TRIGMATCH_EFF",
                                    #
                                    #"m_REFF_dir" 		   : base_dir + "HTopMultilepAnalysis/PlotUtils/PLOTS_25ns_v24_ElNoIso/CombinedEfficiencies_REBINNED_1",
                                    #"m_EFF_YES_TM_dir" 	   : base_dir + "HTopMultilepAnalysis/PlotUtils/PLOTS_25ns_v24_ElNoIso/CombinedEfficiencies_REBINNED_1_TRIGMATCH_EFF",
                                    #"m_EFF_NO_TM_dir" 		   : base_dir + "HTopMultilepAnalysis/PlotUtils/PLOTS_25ns_v24_ElNoIso/CombinedEfficiencies_REBINNED_1_NOT_TRIGMATCH_EFF",
                                    #
                                    #"m_REFF_dir" 		   : base_dir + "HTopMultilepAnalysis/PlotUtils/PLOTS_25ns_v24_ElNoIso/CombinedEfficiencies_REBINNED_2",
                                    #"m_EFF_YES_TM_dir" 	   : base_dir + "HTopMultilepAnalysis/PlotUtils/PLOTS_25ns_v24_ElNoIso/CombinedEfficiencies_REBINNED_2_TRIGMATCH_EFF",
                                    #"m_EFF_NO_TM_dir" 		   : base_dir + "HTopMultilepAnalysis/PlotUtils/PLOTS_25ns_v24_ElNoIso/CombinedEfficiencies_REBINNED_2_NOT_TRIGMATCH_EFF",
                                    #
                                    #"m_REFF_dir" 		   : base_dir + "HTopMultilepAnalysis/PlotUtils/PLOTS_25ns_v24_ElNoIso/CombinedEfficiencies_REBINNED_3",
                                    #"m_EFF_YES_TM_dir" 	   : base_dir + "HTopMultilepAnalysis/PlotUtils/PLOTS_25ns_v24_ElNoIso/CombinedEfficiencies_REBINNED_3_TRIGMATCH_EFF",
                                    #"m_EFF_NO_TM_dir" 		   : base_dir + "HTopMultilepAnalysis/PlotUtils/PLOTS_25ns_v24_ElNoIso/CombinedEfficiencies_REBINNED_3_NOT_TRIGMATCH_EFF",
                                    #
                                    # ------------------------------------------------------------
                                    #
                                    # v26
                                    #
                                    #"m_REFF_dir" 		   : base_dir + "HTopMultilepAnalysis/PlotUtils/PLOTS_25ns_v26/CombinedEfficiencies_LeptonMVA_410501",
                                    #"m_REFF_dir" 		   : base_dir + "HTopMultilepAnalysis/PlotUtils/PLOTS_25ns_v26/CombinedEfficiencies_LeptonMVA_410000",
                                    #
                                    #"m_REFF_dir" 		   : base_dir + "HTopMultilepAnalysis/PlotUtils/PLOTS_25ns_v26/CombinedEfficiencies_LeptonCutBased_410501",
                                    #"m_REFF_dir" 		   : base_dir + "HTopMultilepAnalysis/PlotUtils/PLOTS_25ns_v26/CombinedEfficiencies_LeptonCutBased_410000",
                                    #"m_useCutBasedLep"             : True,
                                    #
                                    #"m_REFF_dir" 		   : base_dir + "HTopMultilepAnalysis/PlotUtils/PLOTS_25ns_v26/MMRates_DATA/OutputPlots_MMRates_25ns_v26_LeptonMVA_DDQMisID",
                                    "m_REFF_dir" 		   : base_dir + "HTopMultilepAnalysis/PlotUtils/OutputPlots_MMRates_25ns_v26_LeptonMVA_MCQMisID",
                                    #
                                    # ------------------------------------------------------------
                                    #
				    "m_systematics_list"           : "Nominal:,Stat:UncorrBins,N_TTBarW:CorrBins,D_TTBarW:CorrBins,N_OtherPromptSS:CorrBins,D_OtherPromptSS:CorrBins,N_FakesOS:CorrBins,D_FakesOS:CorrBins,N_QMisID:UncorrBins,D_QMisID:UncorrBins",
				    #"m_systematics_list"           : "Nominal:Stat:UncorrBins",
                                    "m_correlatedMMWeights"        : False,
                                    #
				    "m_useTrigMatchingInfo"        : False,
				    #
				    "m_Efficiency_Filename"        : "LeptonEfficiencies.root",
				    #
                                    "m_doMMClosure"                : False,
                                    "m_useEtaParametrisation"      : False,
				    "m_useTEfficiency"             : False,
                                  }

# Instantiate the NTupleSvc algorithm

ntuplesvc = ROOT.EL.NTupleSvc(HTopMultilepNTupReprocesserDict["m_outputNTupStreamName"])

# Copy ALL branches over from the input TTree

print("Copying all branches from input TTree to output...")
ntuplesvc.copyBranch(".*")

# Add the algorithms to the job.
#
# Here order matters!

c._algorithms.append(ntuplesvc)
c.setalg("HTopMultilepNTupReprocesser", HTopMultilepNTupReprocesserDict)
