"""
Simple script for bookkeeping skims/ntuples. Modeled after Submit_Example_lhCN.py.

Authors:
 Alex Tuna

Instructions:
1. Place this file in Dev/ folder
2. Run it from the main folder, e.g.:
    > python Dev/Bookkeep.py --dir=./datasets/

Options:
    --input
        The directory of datasets.
    --metadata
        The name of the metadata histogram.
    --output
        The name of the output file where the table is dumped.

^.^

Created: 2012-11-27 <tuna>
Updated: 2013-01-16 <tuna>
"""

import optparse
import os
import sys
import ROOT
import time

sys.path.append(os.path.abspath(os.path.curdir))
from Core import DatasetManager

timestamp = time.strftime('%Y-%m-%d-%Hh%Mm%Ss')

# handle args
parser = optparse.OptionParser()
parser.add_option('-i', '--input',    dest='input',    type=str, default='./datasets/')
parser.add_option('-m', '--metadata', dest='metadata', type=str, default='TotalEvents')
parser.add_option('-o', '--output',   dest='output',   type=str, default='table.bookkeep.%s.txt' % timestamp)
ops,args = parser.parse_args()

# function for deriving dsid from d3pd
def dsid(d3pd):
    if d3pd.count('group.') + d3pd.count('user.') == 0:
        return d3pd.split('.')[1]
    if d3pd.count('group.') + d3pd.count('user.') == 1:
        return d3pd.split('.')[3]
    if d3pd.count('group.') + d3pd.count('user.') == 2:
        return d3pd.split('.')[5]
    sys.exit(' ERROR: Unable to parse dsid from %s. Exiting.' % d3pd)

# function for deriving name from d3pd
def name(d3pd):
    if not d3pd: return ''
    if 'data12' in d3pd: # temporary hack because mc12_8TeV was removed from d3pd output names.
        full_name = d3pd.split('.')[2] if not d3pd.startswith('user') and not d3pd.startswith('group') else d3pd.split('.')[4]
    else:
        full_name = d3pd.split('.')[1] if not d3pd.startswith('user') and not d3pd.startswith('group') else d3pd.split('.')[3]
    long_name = full_name
    # tunes
    long_name = long_name.replace('Auto'         , '')
    long_name = long_name.replace('A2CTEQ6L1'    , '')
    long_name = long_name.replace('AUET2BCTEQ6L1', '')
    long_name = long_name.replace('AUET2CTEQ6L1' , '')
    long_name = long_name.replace('AUET2CT10'    , '')
    long_name = long_name.replace('AU2CTEQ6L1'   , '')
    long_name = long_name.replace('AU2MSTW2008LO', '')
    long_name = long_name.replace('AU2'          , '')
    long_name = long_name.replace('CT10'         , '')
    # physics boundaries
    long_name = long_name.replace('_LeptonFilter', '')
    long_name = long_name.replace('_pt20'        , '')
    long_name = long_name.replace('_Mll150to250' , '_150M250')
    long_name = long_name.replace('_Mll250to400' , '_250M400')
    long_name = long_name.replace('_Mll400'      , '_M400')
    # misc
    long_name = long_name.replace('jetjet'    , '')
    long_name = long_name.replace('.merge.'   , '.')
    long_name = long_name.replace('physics_'  , '')
    long_name = long_name.replace('tautaulh'  , '')
    long_name = long_name.replace('unfiltered', '')
    long_name = long_name.replace('ATau', '')
    # punctuation
    long_name = long_name.replace('.', '')
    long_name = long_name.replace('_', '')
    return long_name

# function for twiki-formatting
def row(d3pd, obs, exp):
    agree = '%ICON{choice-yes}%' if obs == exp else '%ICON{choice-no}%'
    if exp != 0:
        fraction = obs / exp
    elif exp == obs == 0:
        fraction = 1.0
    else:
        fraction = 0.0
    return '| =%s= | %s | %i | %i | %s | %.3f |\n' % (name(d3pd), dsid(d3pd), obs, exp, agree, fraction)

# Load the dataset names
Data  = 'Files/Data_D3PD_datasets.txt'
MC    = 'Files/MC_D3PD_datasets.txt'
Embed = 'Files/Embedding_D3PD_datasets.txt'
dm = DatasetManager.DatasetManager([MC, Data, Embed])

# Dataset tags
version = "lhCNv00-00"
mc_tags = ['mc12a', 'p1130']
emb_tags = ['EXT0']
data_tags = ['data12', 'p1131'] # JetTauEtmiss, Muons, Egamma

# Note: Uncomment sample groups that you want to bookkeep
d3pds  = []
d3pds += dm.getList(tags=["Muons", "PeriodA"]+data_tags)

### Signal samples
#d3pds += dm.getList(tags=['Signal', 'VBF']+mc_tags)
#d3pds += dm.getList(tags=['Signal', 'ggF']+mc_tags)
#d3pds += dm.getList(tags=['Signal', 'WH']+mc_tags)
#d3pds += dm.getList(tags=['Signal', 'ZH']+mc_tags)

### Background samples
#d3pds += dm.getList(tags=['Background', 'Z+jets']+mc_tags)
#d3pds += dm.getList(tags=['Background', 'W+jets']+mc_tags)
#d3pds += dm.getList(tags=['Background', 'Top']+mc_tags)
#d3pds += dm.getList(tags=['Background', 'Diboson']+mc_tags)
#d3pds += dm.getList(tags=['Background', 'DYZ+jets']+mc_tags)

### Special samples
#d3pds += dm.getList(tags=['Background', 'Z+VBF']+mc_tags)
#d3pds += dm.getList(tags=['Background', 'DYZ+VBF']+mc_tags)
#d3pds += dm.getList(tags=['embedding12', 'test']+emb_tags)
#d3pds += dm.getList(tags=['embedding12', 'filtered', 'lephad']+emb_tags)

### Data samples
#d3pds += dm.getList(tags=["Egamma", "PeriodA"]+data_tags)
#d3pds += dm.getList(tags=["Egamma", "PeriodB"]+data_tags)
#d3pds += dm.getList(tags=["Egamma", "PeriodC"]+data_tags)
#d3pds += dm.getList(tags=["Egamma", "PeriodD"]+data_tags)
#d3pds += dm.getList(tags=["Egamma", "PeriodE"]+data_tags)
#d3pds += dm.getList(tags=["Muons", "PeriodA"]+data_tags)
#d3pds += dm.getList(tags=["Muons", "PeriodB"]+data_tags)
#d3pds += dm.getList(tags=["Muons", "PeriodC"]+data_tags)
#d3pds += dm.getList(tags=["Muons", "PeriodD"]+data_tags)
#d3pds += dm.getList(tags=["Muons", "PeriodE"]+data_tags)
#d3pds += dm.getList(tags=["JetTauEtmiss", "PeriodA"]+data_tags)
#d3pds += dm.getList(tags=["JetTauEtmiss", "PeriodB"]+data_tags)
#d3pds += dm.getList(tags=["JetTauEtmiss", "PeriodC"]+data_tags)
#d3pds += dm.getList(tags=["JetTauEtmiss", "PeriodD"]+data_tags)
#d3pds += dm.getList(tags=["JetTauEtmiss", "PeriodE"]+data_tags)

# remove duplicates
d3pds = list(set(d3pds))

# for string formatting
digits = len(str(len(d3pds)))
chars = max([len(name(d3pd)) for d3pd in d3pds])

# get AMI usr and pw
print
import pyAMI.pyAMI as AMI
import getpass as gp
usr = raw_input( ' AMI username: ')
pw  = gp.getpass(' AMI password: ')
print

# begin output
if ops.output:
    print ' Writing to %s' % ops.output
    print
    output = open(ops.output, 'w')
    output.write('| *sample* | *dsid* | *# observed events* | *# expected events* | *agree* | *obs/exp* |\n')

# begin loop
for d3pd in d3pds:

    # retrieve expected number of events from pyAMI
    client  = AMI.AMI(False)
    cmd  = []
    cmd += ['getDatasetInfo']
    cmd += ['logicalDatasetName=%s' % d3pd]
    cmd += ['AMIUser=%s'            %  usr]
    cmd += ['AMIPass=%s'            %   pw]
    exe = client.execute(cmd)
    outdict = exe.getDict()
    eventsExp = float(outdict['Element_Info']['row_1']['totalEvents'])

    # retrieve observed number of events from metadata
    datasetsRaw   = os.listdir(os.path.abspath(ops.input))
    datasetsMatch = filter(lambda d: dsid(d) == dsid(d3pd), datasetsRaw)
    datasetsAbs   = [os.path.join(ops.input, dataset) for dataset in datasetsMatch]
    eventsObs = 0.0
    for dataset in datasetsAbs:
        filesRel = filter(lambda f: '.root' in f, os.listdir(dataset))
        filesAbs = [os.path.join(dataset, fileRel) for fileRel in filesRel]
        for filename in filesAbs:
            file = ROOT.TFile.Open(filename)
            if file:
                hist = file.Get(ops.metadata)
                if hist:
                    eventsObs += float(hist.GetBinContent(1))
                file.Close()

    # squawk
    result = 'AGREE' if eventsExp == eventsObs else 'DISAGREE'
    print '%-*s of %-*s: %-*s %8s %10s %10s %8s' % (digits, d3pds.index(d3pd), digits, len(d3pds), chars, name(d3pd), dsid(d3pd), eventsExp, eventsObs, result)
    if output:
        output.write(row(d3pd, eventsObs, eventsExp))

# end output
if output:
    output.close()

print
print ' Done! ^.^'
print

