import ROOT
from xAH_config import xAH_config
import sys, os

sys.path.insert(0, os.environ['ROOTCOREBIN']+"/user_scripts/HTopMultilepAnalysis/")

c = xAH_config()

event_branches = ["EventNumber","RunNumber","mc_channel_number","isSS01","dilep_type","trilep_type",
                  "is_T_T","is_T_AntiT","is_AntiT_T","is_AntiT_AntiT",
                  "QMisIDWeight","QMisIDWeight_up","QMisIDWeight_dn",
		  "MMWeight",
		  "MMWeight_lep0_r_Stat_up",
		  "MMWeight_lep0_r_Stat_dn",
		  "MMWeight_lep1_r_Stat_up",
		  "MMWeight_lep1_r_Stat_dn",
		  "MMWeight_lep0_f_Stat_up",
		  "MMWeight_lep0_f_Stat_dn",
		  "MMWeight_lep1_f_Stat_up",
		  "MMWeight_lep1_f_Stat_dn"]

lep_branches   = ["lep_ID_0","lep_Pt_0","lep_Eta_0","lep_Phi_0","lep_EtaBE2_0","lep_isTightSelected_0","lep_isTrigMatch_0",
                  "lep_ID_1","lep_Pt_1","lep_Eta_1","lep_Phi_1","lep_EtaBE2_1","lep_isTightSelected_1","lep_isTrigMatch_1"]

branches_to_activate = event_branches + lep_branches

# Trick to pass the list as a comma-separated string to the C++ algorithm
#
branches_to_activate_str = ",".join(branches_to_activate)

# Instantiate the main algorithm
#

base_dir = "/imports/home/mmilesi/PhD/ttH_MultiLeptons/RUN2/HTopMultilepAnalysisCode/trunk/"
#base_dir = "/afs/cern.ch/user/m/mmilesi/ttH/RUN2/HTopMultilepAnalysisCode/trunk"

HTopMultilepNTupReprocesserDict = { "m_name"                       : "HTopMultilepNTupReprocesser",
                                    "m_debug"                      : False,
                                    "m_verbose"                    : False,
				    "m_outputNTupStreamName"       : "output",
                                    "m_inputBranches"              : branches_to_activate_str,
                                    "m_weightToCalc"               : "QMisID", #,MM",
                                    "m_QMisIDRates_dir"            : "$ROOTCOREBIN/data/HTopMultilepAnalysis/External/",
                                    "m_QMisIDRates_Filename_T"     : "QMisIDRates_Data_2016_T_25ns_v21.root",
                                    "m_QMisIDRates_Filename_AntiT" : "QMisIDRates_Data_2016_antiT_25ns_v21.root",
                                    "m_QMisIDRates_Histname_T"     : "LikelihoodEtaPtTight",
                                    "m_QMisIDRates_Histname_AntiT" : "LikelihoodEtaPtLoose",
                                    #"m_useTAntiTRates"             : True, # for v19
                                    #
				    #"m_REFF_dir" 		   : base_dir + "HTopMultilepAnalysis/PlotUtils/PlotVault/PLOTS_25ns_v19/OutputPlots_MMRates_25ns_v19_DDQMisID",
				    #
				    "m_REFF_dir" 		   : base_dir + "HTopMultilepAnalysis/PlotUtils/PlotVault/PLOTS_25ns_v19/OutputPlots_MMRates_25ns_v19_DDQMisID_QMisIDSys",
				    #"m_REFF_dir" 		   : base_dir + "HTopMultilepAnalysis/PlotUtils/PlotVault/PLOTS_25ns_v19/OutputPlots_MMRates_25ns_v19_DDQMisID_QMisIDSys_WrongOverflow",
				    #
				    #"m_REFF_dir" 		   : base_dir + "HTopMultilepAnalysis/PlotUtils/PlotVault/PLOTS_25ns_v19/OutputPlots_MMRates_25ns_v19_DDQMisID_CPPM",
				    #"m_EFF_YES_TM_dir"  	   : base_dir + "HTopMultilepAnalysis/PlotUtils/PlotVault/PLOTS_25ns_v19/OutputPlots_MMRates_25ns_v19_DDQMisID_Probe_YES_TM",
				    #"m_EFF_NO_TM_dir"		   : base_dir + "HTopMultilepAnalysis/PlotUtils/PlotVault/PLOTS_25ns_v19/OutputPlots_MMRates_25ns_v19_DDQMisID_Probe_NO_TM",
                                    #
				    #"m_REFF_dir" 		   : base_dir + "HTopMultilepAnalysis/PlotUtils/PlotVault/PLOTS_25ns_v19/OutputPlots_MMClosureRates_NoCorrections_25ns_v19",
				    #"m_EFF_YES_TM_dir"	           : base_dir + "HTopMultilepAnalysis/PlotUtils/PlotVault/PLOTS_25ns_v19/OutputPlots_MMClosureRates_NoCorrections_25ns_v19_Probe_YES_TM",
				    #"m_EFF_NO_TM_dir" 	           : base_dir + "HTopMultilepAnalysis/PlotUtils/PlotVault/PLOTS_25ns_v19/OutputPlots_MMClosureRates_NoCorrections_25ns_v19_Probe_NO_TM",
				    #
				    #"m_REFF_dir"  		   : base_dir + "HTopMultilepAnalysis/PlotUtils/PlotVault/PLOTS_25ns_v20_02/OutputPlots_MMClosureRates_NoCorrections_TTBarNonAllHad_25ns_v20_02",
				    #"m_REFF_dir"  		   : base_dir + "HTopMultilepAnalysis/PlotUtils/PlotVault/PLOTS_25ns_v20_02/OutputPlots_MMClosureRates_NoCorrections_TTBarSemilep_25ns_v20_02",
                                    #
				    "m_systematics_list"           : "Nominal,Stat,numerator_QMisID,denominator_QMisID",
				    "m_useTrigMatchingInfo"        : False,
				    "m_Efficiency_Filename"        : "LeptonEfficiencies.root",
                                    "m_doMMClosure"                : False,
                                    "m_useEtaParametrisation"      : False,
				    "m_useTEfficiency"             : False,
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
