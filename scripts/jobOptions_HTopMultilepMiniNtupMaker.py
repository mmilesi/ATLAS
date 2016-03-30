import ROOT
from xAH_config import xAH_config
import sys, os

sys.path.insert(0, os.environ['ROOTCOREBIN']+"/user_scripts/HTopMultilepAnalysis/")

c = xAH_config()

HTopMultilepMiniNTupMakerDict = { "m_name" : "HTopMultilepMiniNTupMaker", 
                                  "m_debug" : True,
                                  "m_outputNTupleName" : "output",
                                }

# Here order matters!
#
# NB: here users can update values in the dictionaries before setting them to the algorithm
#
c.setalg("HTopMultilepMiniNTupMaker", HTopMultilepMiniNTupMakerDict)


