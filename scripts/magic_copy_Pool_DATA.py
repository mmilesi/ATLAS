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

    copylist = [
        "mkdir 00276262 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00276262.root . && cd ..",
        "mkdir 00276329 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00276329.root . && cd ..",
        "mkdir 00276336 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00276336.root . && cd ..",
        "mkdir 00276416 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00276416.root . && cd ..",
        "mkdir 00276511 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00276511.root . && cd ..",
        "mkdir 00276689 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00276689.root . && cd ..",
        "mkdir 00276778 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00276778.root . && cd ..",
        "mkdir 00276790 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00276790.root . && cd ..",
        "mkdir 00276952 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00276952.root . && cd ..",
        "mkdir 00276954 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00276954.root . && cd ..",
        "mkdir 00278880 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00278880.root . && cd ..",
        "mkdir 00278912 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00278912.root . && cd ..",
        "mkdir 00278968 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00278968.root . && cd ..",
        "mkdir 00279169 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00279169.root . && cd ..",
        "mkdir 00279259 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00279259.root . && cd ..",
        "mkdir 00279279 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00279279.root . && cd ..",
        "mkdir 00279284 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00279284.root . && cd ..",
        "mkdir 00279345 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00279345.root . && cd ..",
        "mkdir 00279515 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00279515.root . && cd ..",
        "mkdir 00279598 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00279598.root . && cd ..",
        "mkdir 00279685 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00279685.root . && cd ..",
        "mkdir 00279813 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00279813.root . && cd ..",
        "mkdir 00279867 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00279867.root . && cd ..",
        "mkdir 00279928 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00279928.root . && cd ..",
        "mkdir 00279932 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00279932.root . && cd ..",
        "mkdir 00279984 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00279984.root . && cd ..",
        "mkdir 00280231 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00280231.root . && cd ..",
        "mkdir 00280273 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00280273.root . && cd ..",
        "mkdir 00280319 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00280319.root . && cd ..",
        "mkdir 00280368 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00280368.root . && cd ..",
        "mkdir 00280423 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00280423.root . && cd ..",
        "mkdir 00280464 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00280464.root . && cd ..",
        "mkdir 00280500 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00280500.root . && cd ..",
        "mkdir 00280520 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00280520.root . && cd ..",
        "mkdir 00280614 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00280614.root . && cd ..",
        "mkdir 00280673 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00280673.root . && cd ..",
        "mkdir 00280753 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00280753.root . && cd ..",
        "mkdir 00280853 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00280853.root . && cd ..",
        "mkdir 00280950 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00280950.root . && cd ..",
        "mkdir 00280977 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00280977.root . && cd ..",
        "mkdir 00281070 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00281070.root . && cd ..",
        "mkdir 00281074 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00281074.root . && cd ..",
        "mkdir 00281075 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00281075.root . && cd ..",
        "mkdir 00281317 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00281317.root . && cd ..",
        "mkdir 00281385 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00281385.root . && cd ..",
        "mkdir 00281411 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00281411.root . && cd ..",
        "mkdir 00282625 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00282625.root . && cd ..",
        "mkdir 00282631 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00282631.root . && cd ..",
        "mkdir 00282712 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00282712.root . && cd ..",
        "mkdir 00282784 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00282784.root . && cd ..",
        "mkdir 00282992 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00282992.root . && cd ..",
        "mkdir 00283074 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00283074.root . && cd ..",
        "mkdir 00283155 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00283155.root . && cd ..",
        "mkdir 00283270 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00283270.root . && cd ..",
        "mkdir 00283429 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00283429.root . && cd ..",
        "mkdir 00283608 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00283608.root . && cd ..",
        "mkdir 00283780 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00283780.root . && cd ..",
        "mkdir 00284006 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00284006.root . && cd ..",
        "mkdir 00284154 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00284154.root . && cd ..",
        "mkdir 00284285 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00284285.root . && cd ..",
        "mkdir 00284420 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00284420.root . && cd ..",
        "mkdir 00284427 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00284427.root . && cd ..",
        "mkdir 00284484 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00284484.root . && cd ..",
        "mkdir 00297730 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00297730.root . && cd ..",
        "mkdir 00298595 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00298595.root . && cd ..",
        "mkdir 00298609 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00298609.root . && cd ..",
        "mkdir 00298633 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00298633.root . && cd ..",
        "mkdir 00298687 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00298687.root . && cd ..",
        "mkdir 00298690 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00298690.root . && cd ..",
        "mkdir 00298771 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00298771.root . && cd ..",
        "mkdir 00298773 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00298773.root . && cd ..",
        "mkdir 00298862 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00298862.root . && cd ..",
        "mkdir 00298967 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00298967.root . && cd ..",
        "mkdir 00299055 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00299055.root . && cd ..",
        "mkdir 00299144 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00299144.root . && cd ..",
        "mkdir 00299147 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00299147.root . && cd ..",
        "mkdir 00299184 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00299184.root . && cd ..",
        "mkdir 00299243 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00299243.root . && cd ..",
        "mkdir 00299584 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00299584.root . && cd ..",
        "mkdir 00300279 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00300279.root . && cd ..",
        "mkdir 00300345 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00300345.root . && cd ..",
        "mkdir 00300415 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00300415.root . && cd ..",
        "mkdir 00300418 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00300418.root . && cd ..",
        "mkdir 00300487 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00300487.root . && cd ..",
        "mkdir 00300540 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00300540.root . && cd ..",
        "mkdir 00300571 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00300571.root . && cd ..",
        "mkdir 00300600 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00300600.root . && cd ..",
        "mkdir 00300655 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00300655.root . && cd ..",
        "mkdir 00300687 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00300687.root . && cd ..",
        "mkdir 00300784 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00300784.root . && cd ..",
        "mkdir 00300800 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00300800.root . && cd ..",
        "mkdir 00300863 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00300863.root . && cd ..",
        "mkdir 00300908 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00300908.root . && cd ..",
        "mkdir 00301912 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00301912.root . && cd ..",
        "mkdir 00301915 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00301915.root . && cd ..",
        "mkdir 00301918 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00301918.root . && cd ..",
        "mkdir 00301932 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00301932.root . && cd ..",
        "mkdir 00301973 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00301973.root . && cd ..",
        "mkdir 00302053 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00302053.root . && cd ..",
        "mkdir 00302137 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00302137.root . && cd ..",
        "mkdir 00302265 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00302265.root . && cd ..",
        "mkdir 00302269 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00302269.root . && cd ..",
        "mkdir 00302300 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00302300.root . && cd ..",
        "mkdir 00302347 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00302347.root . && cd ..",
        "mkdir 00302380 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00302380.root . && cd ..",
        "mkdir 00302391 && cd $_ && xrdcp root://eospublic.cern.ch//eos/escience/UniTexas/HSG8/multileptons_ntuple_run2/25ns_v17/Data/00302391.root . && cd ..",
    ]

    MAX_PARALLEL=4

    print listchunks(copylist,MAX_PARALLEL)

    for chunk in listchunks(copylist,MAX_PARALLEL):

        print("Copying samples: ")
        print("\n".join("{0} - {1}".format(elem[0],elem[1].split()[1]) for elem in enumerate(chunk)))
        p = multiprocessing.Pool(MAX_PARALLEL)
        p.map(copy,chunk)
	p.close()
        p.join()
