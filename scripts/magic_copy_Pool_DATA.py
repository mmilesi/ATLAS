#!/usr/bin/env python

""" magic_copy_Pool.py: parallelise xrdcp copy via multiprocessing.Pool """

__author__     = "Marco Milesi"
__email__      = "marco.milesi@cern.ch"
__maintainer__ = "Marco Milesi"

import glob, os, sys, subprocess, shutil

import multiprocessing

def listchunks(l, n):
    n = max(1, n)
    return [l[i:i + n] for i in range(0, len(l), n)]

def copy(sample):
    cmd = sample
    subprocess.call(cmd,shell=True)

if __name__ == '__main__':

    username = "mmilesi"

    version     = "25ns_v17"
    sample_type = "Data"

    basedir = "/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/25ns_v17/" + sample_type

    copylist = [
276262,
276329,
276336,
276416,
276511,
276689,
276778,
276790,
276952,
276954,
278880,
278912,
278968,
279169,
279259,
279279,
279284,
279345,
279515,
279598,
279685,
279813,
279867,
279928,
279932,
279984,
280231,
280273,
280319,
280368,
280423,
280464,
280500,
280520,
280614,
280673,
280753,
280853,
280862,
280950,
280977,
281070,
281074,
281075,
281317,
281385,
281411,
282625,
282631,
282712,
282784,
282992,
283074,
283155,
283270,
283429,
283608,
283780,
284006,
284154,
284213,
284285,
284420,
284427,
284484,
297730,
298595,
298609,
298633,
298687,
298690,
298771,
298773,
298862,
298967,
299055,
299144,
299147,
299184,
299243,
299584,
300279,
300345,
300415,
300418,
300487,
300540,
300571,
300600,
300655,
300687,
300784,
300800,
300863,
300908,
301912,
301915,
301918,
301932,
301973,
302053,
302137,
302265,
302269,
302300,
302347,
302380,
302391,
    ]

    cmdlist = []
    for sample in copylist:
        cmd = "cd " + basedir + " && mkdir -p 00" + str(sample) + " && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/" + version + "/" + sample_type + "/00" + str(sample) + ".root ."
        cmdlist.append(cmd)

    MAX_PARALLEL = 4

    #print listchunks(cmdlist,MAX_PARALLEL)

    for chunk in listchunks(cmdlist,MAX_PARALLEL):

        if not os.path.exists("/tmp/krb5cc_1016"):
	    print("Please get a Kerberos ticket first:")
	    krb_auth = "kinit " + username + "@CERN.CH"
	    subprocess.call(krb_auth,shell=True)
        subprocess.call("kinit -R",shell=True)

        print("Copying samples: ")
        print("\n".join("{0} - {1}".format(elem[0],elem[1].split()[5]) for elem in enumerate(chunk)))
        p = multiprocessing.Pool(MAX_PARALLEL)
        p.map(copy,chunk)
	p.close()
        p.join()

    os.chdir(basedir)
    print("Transfer finished!")

