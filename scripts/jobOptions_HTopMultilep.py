import ROOT
from xAH_config import xAH_config
import sys, os

sys.path.insert(0, os.environ['ROOTCOREBIN']+"/user_scripts/HTopMultilepAnalysis/")

from HTopMultilep_config import *

c = xAH_config()

# Here order matters!
#
# NB: here users can update values in the dictionaries before setting them to the algorithm
#
c.setalg("BasicEventSelection", BasicEventSelectionDict)
c.setalg("JetCalibrator", JetCalibratorDict)
c.setalg("MuonCalibrator", MuonCalibratorDict)
c.setalg("ElectronCalibrator", ElectronCalibratorDict)
c.setalg("JetSelector", JetSelectorDict)
c.setalg("MuonSelector", MuonSelectorDict)
c.setalg("ElectronSelector", ElectronSelectorDict)
c.setalg("TauSelector", TauSelectorDict)
c.setalg("METConstructor", METConstructorDict)
c.setalg("OverlapRemover_HTopRun1", OverlapRemoverDict)
c.setalg("BJetEfficiencyCorrector", BJetEfficiencyCorrectorDict)
c.setalg("MuonEfficiencyCorrector", MuonEfficiencyCorrectorDict)
c.setalg("ElectronEfficiencyCorrector", ElectronEfficiencyCorrectorDict)
c.setalg("MuonEfficiencyCorrector", MuonEfficiencyCorrectorTightDict)
c.setalg("ElectronEfficiencyCorrector", ElectronEfficiencyCorrectorTightDict)
c.setalg("HTopMultilepEventSelector", HTopMultilepEventSelectorDict)
c.setalg("TruthMatchAlgo", TruthMatchAlgoDict)
c.setalg("HTopMultilepAnalysis", HTopMultilepAnalysisDict)
c.setalg("HTopMultilepTreeAlgo", HTopMultilepTreeAlgoDict)

