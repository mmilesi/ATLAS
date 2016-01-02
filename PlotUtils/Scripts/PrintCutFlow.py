"""
File    : Scripts/PrintCutFlow.py
Authors : KG <Kong.Guan.Tan@cern.ch>, Guilherme Hanninger <gnunes@cern.ch>, David Jennens <David.Thomas.Jennens@cern.ch>, Geng-Yuan Jeng <Geng-Yuan.Jeng@cern.ch>

Script to print the cutflow stored in an ntuple output
TODO: Show an example
"""


import os
import sys
import glob

defaults_suppress = {
        'NoControl' : ['AntiTau', 'OtherMT', 'HighMT', 'NonIsoLep', 'NoTrig', 'SS'],
        'NoSys': ['EESSysDown', 'EESSysUp', 'ElEnResSysDown', 'ElEnResSysUp', 'JERDown', 'JERUp',
                  'JES_BASEDown', 'JES_BASEUp', 'JES_FLAVDown', 'JES_FLAVUp', 'JES_FWDDown', 'JES_FWDUp',
                  'TESDown', 'TESUp', 'METResSysDown', 'METResSysUp', 'METScaleSysDown', 'METScaleSysUp', 'MuSysDown', 'MuSysUp'],
        'SSOnly': ['OS'],
        'MuOnly': ['El'],
        'ElOnly': ['Mu'],
        'VBFOnly' : ['0-Jet', '1-Jet', 'Boosted', 'NoCat'],
        '1JetOnly': ['0-Jet', 'VBF', 'Boosted', 'NoCat'],
        '0JetOnly': ['1-Jet', 'VBF', 'Boosted', 'NoCat'],
        'Boosted' : ['0-Jet', '1-Jet', 'VBF', 'NoCat'],
    }
defaults_dename = {
        'NoControl' : ['Tau', 'LowMT', 'IsoLep', 'OS'],
        'NoSys': [],
        'SSOnly': ['SS'],
        'MuOnly': ['Mu'],
        'ElOnly': ['El'],
        'VBFOnly' : ['VBF'],
        '1JetOnly': ['1-Jet'],
        '0JetOnly': ['0-Jet'],
        'Boosted' : ['Boosted'],
    }

defaults_suppress['SimpleEl'] = defaults_suppress['NoControl'] + defaults_suppress['NoSys'] + defaults_suppress['ElOnly']
defaults_suppress['SimpleMu'] = defaults_suppress['NoControl'] + defaults_suppress['NoSys'] + defaults_suppress['MuOnly']
defaults_dename['SimpleEl'] = defaults_dename['NoControl'] + defaults_dename['NoSys'] + defaults_dename['ElOnly']
defaults_dename['SimpleMu'] = defaults_dename['NoControl'] + defaults_dename['NoSys'] + defaults_dename['MuOnly']

def main():
    inputpath = 'output/ntuple*.root'

    sys.path.append(os.path.abspath(os.path.curdir))
    from Core import NTupleTools, parseInputArgs, listifyInputFiles

    options = parseInputArgs()
    try:
        inputs = listifyInputFiles(options.inFiles)
    except:
        inputs = listifyInputFiles(inputpath)
    cutflow = NTupleTools.getCutFlowFromHistogram(inputs)

    try:
        suppress = options.cutFlowSuppress.split(',')
    except:
        suppress = ['NoControl','NoSys']

    fullsuppress = []
    fulldename = []
    for s in suppress:
        if s in defaults_suppress:
            fullsuppress += defaults_suppress[s]
        if s in defaults_dename:
            fulldename += defaults_dename[s]
        else:
            print "ERROR: No such streamlet group:", s
            print "See Scripts/PrintCutFlow.py for all streamlet groups or to define new ones"
            print "Exiting script..."
            sys.exit()

    NTupleTools.printCutFlow(cutflow, suppressStreamlet=fullsuppress, denameStreamlet=fulldename)

if __name__ == "__main__":
    main()
