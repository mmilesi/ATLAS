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

    version     = "25ns_v20"
    sample_type = "02/Nominal"

    basedir = "/coepp/cephfs/mel/mmilesi/ttH/GroupNTup/" + version + "/" + sample_type

    copylist = [
341177,
341270,
341271,
341988,
341989,
341990,
341992,
341995,
341997,
341998,
341999,
342001,
342004,
342170,
342171,
342172,
342284,
342285,
343266,
343267,
343268,
343365,
343366,
343367,
361063,
361064,
361065,
361066,
361067,
361068,
361069,
361070,
361071,
361072,
361073,
361074,
361075,
361076,
361077,
361078,
361079,
361080,
361081,
361082,
361083,
361084,
361085,
361086,
361087,
361088,
361089,
361090,
361091,
361092,
361093,
361094,
361095,
361096,
361097,
410000,
410001,
410002,
410007,
410009,
410011,
410012,
410013,
410014,
410015,
410016,
410025,
410026,
410049,
410050,
410066,
410067,
410068,
410073,
410074,
410075,
410080,
410081,
410111,
410112,
410113,
410114,
410115,
410116,
410120,
410121,
410142,
410143,
410144,
410155,
410156,
410157,
410159,
410187,
410188,
410189,
410215,
410218,
410219,
410220,
410500,
    ]

    cmdlist = []
    for sample in copylist:
        cmd = "cd " + basedir + " && mkdir -p " + str(sample) + " && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/" + version + "/" + sample_type + "/" + str(sample) + ".root ."
        cmdlist.append(cmd)

    MAX_PARALLEL = 6

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
