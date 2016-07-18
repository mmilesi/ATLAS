import os, subprocess, datetime
from subprocess import call
from argparse import ArgumentParser

date = datetime.datetime.now().strftime("%y%m%d_%H%M%S_")

parser = ArgumentParser()
parser.add_argument( '-w', '--workspace', help='Path to workspace',     default='/coepp/cephfs/mel/brianl/125.root' )
parser.add_argument( '-o', '--outputdir', help='Where to store output', default='test' )
result = parser.parse_args() 

nps = [
    "alpha_eg_reso_all",
    "alpha_eg_scale_all",
    "alpha_el_eff_id",
    "alpha_el_eff_iso",
    "alpha_el_eff_reco",
    "alpha_el_eff_slt",
    "alpha_el_eff_tlt",
    "alpha_jet_EtaIntercalib",
    "alpha_jet_JER_NP_0",
    "alpha_jet_JER_NP_1",
    "alpha_jet_JER_NP_2",
    "alpha_jet_JER_NP_3",
    "alpha_jet_JER_NP_4",
    "alpha_jet_JER_NP_5",
    "alpha_jet_JER_NP_6",
    "alpha_jet_JER_NP_7",
    "alpha_jet_JER_NP_8",
    "alpha_jet_JVT_eff",
    "alpha_jet_grouped_NP_1",
    "alpha_jet_grouped_NP_2",
    "alpha_jet_grouped_NP_3",
    "alpha_met_soft_scale",
    "alpha_mu_eff_isostat",
    "alpha_mu_eff_isosys",
    "alpha_mu_eff_slt_stat",
    "alpha_mu_eff_slt_sys",
    "alpha_mu_eff_stat",
    "alpha_mu_eff_sys",
    "alpha_mu_eff_tlt_stat",
    "alpha_mu_eff_tlt_sys",
    "alpha_muon_es_id",
    "alpha_muon_es_ms",
    "alpha_muon_es_scale",
    "alpha_tau_eff_eleolr",
    "alpha_tau_eff_id",
    "alpha_tau_eff_reco",
    "alpha_tau_el_tlt_stat_mc",
    "alpha_tau_el_tlt_sys",
    "alpha_tau_ele_tlt_stat_data",
    "alpha_tau_tes_detector",
    "alpha_tau_tes_insitu",
    "alpha_tau_tes_model",
]

n = 0
args = []

for np in nps:
    args.append({'Number' : n , 'NP' : np, 'Curdir' : os.path.abspath(os.curdir), 'WS' : result.workspace, 'Output' : result.outputdir})
    n += 1

template = """#!/bin/bash

#PBS -q long
#PBS -N nllscan_test_%(Number)s

setupATLAS
MYTMP=`mktemp -d $TMPDIR/brianl.XXXX`
echo `mktemp -d $TMPDIR/brianl.XXXX`
rsync -arvxSH /coepp/cephfs/mel/brianl/FitBox $MYTMP
cd $MYTMP/FitBox
rcSetup Base,2.3.48
FitBox -n %(NP)s -w %(WS)s -o %(Output)s
cp $MYTMP/FitBox/%(Output)s/*.root %(Curdir)s/%(Output)s/
"""

os.system('mkdir -p tmp/')

def submitQ( input ):
    f = open('tmp/'+date+str(input['Number'])+'.sh' , 'w')
    f.write(template % input)
    f.close()
    call(['chmod', '755', 'tmp/'+date+str(input['Number'])+'.sh'], shell=False)
    call(['qsub', 'tmp/'+date+str(input['Number'])+'.sh'], shell=False)

call(['mkdir', result.outputdir], shell=False)
map( submitQ, args ) 
