import ROOT
from xAH_config import xAH_config
import sys, os

sys.path.insert(0, os.environ['ROOTCOREBIN']+"/user_scripts/HTopMultilepAnalysis/")

c = xAH_config()

event_branches = ["EventNumber","RunNumber","mc_channel_number","isSS01","dilep_type","trilep_type",
                  "is_T_T","is_T_AntiT","is_AntiT_T","is_AntiT_AntiT",
                  "is_TMVA_TMVA","is_TMVA_AntiTMVA","is_AntiTMVA_TMVA","is_AntiTMVA_AntiTMVA",
                  "nJets_OR_T","nJets_OR_T_MV2c10_70",
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
                                    # "m_weightToCalc"               : "QMisID,MM",
                                    # "m_weightToCalc"               : "MM",
                                    "m_weightToCalc"               : "QMisID",
                                    "m_QMisIDRates_dir"            : "$ROOTCOREBIN/data/HTopMultilepAnalysis/External/",
                                    "m_QMisIDRates_Filename_T"     : "Rates_v27.root",
                                    "m_QMisIDRates_Filename_AntiT" : "Rates_v27_AntiT.root",
                                    # "m_QMisIDRates_Filename_T"     : "Rates_v27_For3l.root",
                                    # "m_QMisIDRates_Filename_AntiT" : "Rates_v27_For3l_AntiT.root",
                                    # "m_QMisIDRates_Filename_T"     : "Rates_Data_3l_2D_tight.root",
                                    # "m_QMisIDRates_Filename_AntiT" : "Rates_Data_3l_2D_Loose.root",
                                    "m_QMisIDRates_Histname_T"     : "LikelihoodEtaPtTight",
                                    "m_QMisIDRates_Histname_AntiT" : "LikelihoodEtaPtLoose",
                                    # "m_QMisIDRates_Histname_T"     : "LikelihoodEtaPt",
                                    # "m_QMisIDRates_Histname_AntiT" : "LikelihoodEtaPt",
                                    "m_useTAntiTRates"             : False, # --> set this option to True if QMisID rates have NOT been measured independently for T and AntiT electrons - Set True for v19
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
                                    #"m_REFF_dir" 		   : base_dir + "HTopMultilepAnalysis/PlotUtils/OutputPlots_MMClosureRates_25ns_v26_LeptonMVA", # Closure
                                    #
                                    # ------------------------------------------------------------
                                    #
                                    # v27
                                    #
                                    # "m_REFF_dir" 		   : base_dir + "HTopMultilepAnalysis/PlotUtils/PLOTS_25ns_v27/OutputPlots_MMRates_25ns_v27",
                                    # "m_REFF_dir" 		   : base_dir + "HTopMultilepAnalysis/PlotUtils/PLOTS_25ns_v27/OutputPlots_MMRates_TTWx2_25ns_v27",
                                    "m_REFF_dir" 		   : base_dir + "HTopMultilepAnalysis/PlotUtils/PLOTS_25ns_v27_v2/OutputPlots_MMClosureRates_25ns_v27", # Closure
                                    #
                                    # ------------------------------------------------------------
                                    #
                                    # v28
                                    #
                                    # "m_REFF_dir" 		   : base_dir + "HTopMultilepAnalysis/PlotUtils/PLOTS_25ns_v28/",
                                    # "m_REFF_dir" 		   : base_dir + "HTopMultilepAnalysis/PlotUtils/PLOTS_25ns_v28/", # Closure
                                    #
                                    # ------------------------------------------------------------
                                    #
                                    #"m_parametrisation_list"       : "Real_El:Pt,Real_Mu:Pt,Fake_El:Pt,Fake_Mu:Pt",
                                    #"m_parametrisation_list"       : "Real_El:Pt,Real_Mu:Pt,Fake_El:PtxEta,Fake_Mu:Pt",
                                    "m_parametrisation_list"       : "Real_El:Pt,Real_Mu:Pt,Fake_El:NBJets_VS_Pt,Fake_Mu:Pt",
				    # "m_systematics_list"         : "Nominal:,Stat:UncorrBins,ND_TTV:CorrBins,ND_VV:CorrBins,ND_OtherPromptSS:CorrBins,ND_FakesOS:CorrBins,N_QMisID:UncorrBins,D_QMisID:UncorrBins",
				    "m_systematics_list"           : "Nominal:,Stat:UncorrBins",
                                    "m_correlatedMMWeights"        : False,
                                    #
				    "m_useTrigMatchingInfo"        : False,
				    #
				    "m_Efficiency_Filename"        : "LeptonEfficiencies.root",
				    # "m_Efficiency_Filename"        : "LeptonEfficiencies_FakeElPtRescaledPhConv.root",
				    #
                                    "m_doMMClosure"                : False,
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
