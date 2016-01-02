"""
File    : Core/NTupleTools.py
Authors : KG <Kong.Guan.Tan@cern.ch>

Functions to parse/create/read ntuple inputs and outputs
TODO: Show examples of usage?
"""

import os, sys, math, glob
#from ROOT import TH1D, TFile, TObjString, TProofOutputFile, TObjArray
from ROOT import TH1D, TFile, TObjString, TObjArray

from Core import CodeGenerator, compileMinimal, listifyInputFiles, parseInputArgs

filedir = os.path.dirname(os.path.abspath(__file__))

outputFile = None
outputProofFile = None

def printBranchesFromFile(treename, inputpath):
    from ROOT import TChain
    tree = TChain(treename)
    tree.Add(inputpath)

    for l in tree.GetListOfLeaves():
        print "%60s , %20s" % ("'"+l.GetName()+"'", "'"+l.GetTypeName()+"'")

def setBranchAddresses(tree, allobjects):
    obs = CodeGenerator.obs
    if not obs:
        print "Need to setup with a Branch.py (or similar) file first!"
        return

    from ROOT import AddressOf
    tree.SetBranchStatus("*", 0)
    for name in sorted(obs.keys()):
        o = getattr(allobjects, name.lower())
        for bname in obs[name].elemsList():
            bname, baddress, btype, bdefault, bbranch = obs[name].findBranch(bname)
            if tree.GetBranch(baddress):
                tree.SetBranchStatus(baddress, 1)
                tree.SetBranchAddress(baddress, AddressOf(o, bname))
    return tree

def loadNTuple(treename, inputdir):
    obs = CodeGenerator.obs
    if not obs:
        print "Need to setup with a Branch.py (or similar) file first!"
        return

    location = listifyInputFiles(inputdir)

    from ROOT import TChain, AddressOf
    tree = TChain(treename)
    for l in location:
        tree.Add(l)
    return tree

def makeOutputFile(outputpath, recreate=True):
    global outputFile
    global outputProofFile
    if outputFile or outputProofFile:
        return
    options = parseInputArgs()
    if recreate:
        opt = "RECREATE"
    else:
        opt = "UPDATE"
    if options.noProof:
        outputFile = TFile.Open(outputpath, opt)
    else:
        dirname = os.path.dirname(outputpath)
        basename = '.temp.' + os.path.basename(outputpath)
        outputProofFile = TProofOutputFile(basename)
        outputFile = outputProofFile.OpenFile(opt)

def makeNTuple(treename, allobjects, dump):
    obs = CodeGenerator.obs
    if not obs:
        print "Need to setup with a Branch.py (or similar) file first!"
        return

    from ROOT import TTree, TFile, AddressOf, gROOT, gInterpreter
    if "/" in treename:
        if outputFile:
            dirname = treename.split('/')[0]
            if not outputFile.GetDirectory(dirname):
                outputFile.mkdir(dirname)
            outputFile.cd(dirname)
        treename = treename.split('/')[1]
    else:
        if outputFile:
            outputFile.cd()
    tree = TTree(treename, "NTuple output of AnalysisFramework")

    if not dump:
        objdumplist = []
    elif dump is True:
        objdumplist = obs.keys()
    elif type(dump) is str:
        objdumplist = [dump]
    else:
        objdumplist = dump

    for name in sorted(obs.keys()):
        if not name in objdumplist:
            continue
        o = getattr(allobjects, name.lower())
        for bname in obs[name].elemsList():
            bname, baddress, btype, bdefault, bbranch = obs[name].findBranch(bname)
            if not (obs[name].dumplist is True or bname in obs[name].dumplist): continue
            if bbranch:
                tree.Branch(baddress, AddressOf(o, bname), baddress+bbranch)
            else:
                tree.Branch(baddress, btype, getattr(o, bname))
    return tree

def getTotalEventsHistogram(inputdir):
    inputpath = listifyInputFiles(inputdir)

    if outputFile:
        outputFile.cd()
    totalEventsHistogram = TH1D("TotalEvents", "", 2, 1, 3)
    try:
        for d in inputpath:
            f = TFile.Open(d)
            htemp = f.Get("TotalEvents")
            if not htemp: htemp = f.Get("cutflow")
            else: totalEventsHistogram.SetTitle(htemp.GetTitle())
            for i in range(1, 3):
                totalEventsHistogram.SetBinContent(i, totalEventsHistogram.GetBinContent(i) + htemp.GetBinContent(i))
            f.Close()
        return totalEventsHistogram
    except:
        return None

def makeTotalEventsHistogram(flow):
    if outputFile:
        outputFile.cd()
    totalEventsHistogram = TH1D("TotalEvents", "", 2, 1, 3)
    totalEvents = flow.totalEvents[""]
    totalEventsHistogram.SetBinContent(1, totalEvents.passed)
    totalEventsHistogram.SetBinContent(2, totalEvents.passedW)
    return totalEventsHistogram

def getCutFlowFromHistogram(inputdir):
    try:
        from ROOT import AnalysisFramework
    except:
        compileMinimal()
        from ROOT import AnalysisFramework
    CutFlowHist = AnalysisFramework.CutFlows.CutFlowHist

    inputpath = listifyInputFiles(inputdir)

    htemp = CutFlowHist("CutFlow", "CutFlow output of AnalysisFramework", 400000, 0, 1)
    for d in inputpath:
        f = TFile.Open(d)
        heach = f.Get("CutFlow")
        #for i in range(heach.GetNbinsX()):
        #    if not heach.GetBinLabel(i+1):
        #        break
        #    print i+1, heach.GetBinLabel(i+1)
        col = TObjArray()
        col.Add(heach)
        htemp.Merge(col)
        f.Close()

    #xaxis = htemp.GetXaxis()
    temp = {}
    for i in range(htemp.GetNbinsX()):
        label = htemp.GetBinLabel(i+1)
        if not label:
            continue
        flownum = int(label.split('/')[0])
        isweighted = label.split('/')[1] == 'W'
        flowname = label.split('/')[2]
        streamlet = label.split('/')[3]
        if isweighted:
            raw = 0.
            weighted = htemp.GetBinContent(i+1)
        else:
            raw = htemp.GetBinContent(i+1)
            weighted = 0.
        if not flownum in temp:
            temp[flownum] = (flowname, {})
        flownametemp, numberstemp = temp[flownum]
        if not streamlet in numberstemp:
            numberstemp[streamlet] = (raw, weighted)
        else:
            rawtemp, weightedtemp = numberstemp[streamlet]
            numberstemp[streamlet] = (raw+rawtemp, weighted+weightedtemp)

    cutflow = []
    totalEvents = getTotalEventsHistogram(inputdir)
    if totalEvents:
        cutflow.append( ('OriginalTotalEvents', {'All': (totalEvents.GetBinContent(1), totalEvents.GetBinContent(2))}) )
    for i in sorted(temp.keys()):
        cutflow.append(temp[i])
    return cutflow

def makeCutFlowHistogram(cutflow):
    try:
        from ROOT import AnalysisFramework
    except:
        compileMinimal()
        from ROOT import AnalysisFramework
    CutFlowHist = AnalysisFramework.CutFlows.CutFlowHist

    if outputFile:
        outputFile.cd()
    totalsize = 0
    for flowname, numbers in cutflow:
        totalsize += len(numbers)
    cutFlowHist = CutFlowHist("CutFlow", "CutFlow output of AnalysisFramework", 400000, 0, 1)
    #xaxis = cutFlowHist.GetXaxis()
    flownum = 0
    index = 0
    for flowname, numbers in cutflow:
        for streamlet in sorted(numbers.keys()):
            raw, weighted = numbers[streamlet]
            cutFlowHist.SetBinContent(index+1, raw)
            cutFlowHist.SetBinContent(index+2, weighted)
            cutFlowHist.SetBinLabel(index+1, str(flownum) + '/R/' + flowname + '/' + streamlet)
            cutFlowHist.SetBinLabel(index+2, str(flownum) + '/W/' + flowname + '/' + streamlet)
            index += 2
        flownum += 1
    print 'Created CutFlow histogram with', index, 'entries.'
    return cutFlowHist

def copyLumi(inputdir):
    inputpath = listifyInputFiles(inputdir)

    lumidir = outputFile.mkdir("Lumi")
    for d in inputpath:
        f = TFile.Open(d)
        try:
            l  = f.GetDirectory("Lumi")
            keys = l.GetListOfKeys()
            for entry in range(keys.GetEntries()):
                objstr = l.Get(keys.At(entry).GetName() + ";" + str(keys.At(entry).GetCycle()))
                if objstr:
                    lumidir.cd()
                    objnew = TObjString(objstr.GetString().Data())
                    objnew.Write(keys.At(entry).GetName())
        except:
            pass
        f.Close()
    outputFile.cd()

def getSystematicNames(inputdir):
    inputpath = listifyInputFiles(inputdir)

    d = inputpath[0]
    f = TFile.Open(d)

    systlist = []
    l  = f.GetDirectory("SystematicsUP")
    if l:
        keys = l.GetListOfKeys()
        for entry in range(keys.GetEntries()):
            systlist.append('SystematicsUP/' + keys.At(entry).GetName())
    l  = f.GetDirectory("SystematicsDOWN")
    if l:
        keys = l.GetListOfKeys()
        for entry in range(keys.GetEntries()):
            systlist.append('SystematicsDOWN/' + keys.At(entry).GetName())
    return systlist

def printCutFlow(cutflow, suppressStreamlet=[], denameStreamlet=[], suppressEmpty=False):
    currentheadings = None
    currenttemplate = ""
    print
    for flowname, numbers in cutflow:
        if suppressStreamlet:
            newNumbers = {}
            for key in numbers:
                streamletslist = key.split('&')
                suppress = False
                for s in suppressStreamlet:
                    if s in streamletslist:
                        suppress = True
                        break
                if not suppress:
                    newNumbers[key] = numbers[key]
            numbers = newNumbers

        if denameStreamlet:
            newNumbers = {}
            for key in numbers:
                streamletslist = key.split('&')
                streamletslist = [s for s in streamletslist if not s in denameStreamlet]
                newNumbers['&'.join(streamletslist)] = numbers[key]
            numbers = newNumbers

        if suppressEmpty:
            newNumbers = {}
            for key in numbers:
                raw, weighted = numbers[key]
                if raw:
                    newNumbers[key] = numbers[key]
            numbers = newNumbers

        if not currentheadings == sorted(numbers.keys()):
            currentheadings = sorted(numbers.keys())
            currenttemplate = "%22s " + "%*s%s%*s" * len(numbers)
            args = ["FlowName"]
            
            for streamlet in sorted(numbers.keys()):
                args += [int(math.ceil((23-len(streamlet.replace('&', '')))/2.)), "", streamlet.replace('&', ''), int(math.floor((23-len(streamlet.replace('&', '')))/2.)), ""]
            print '-' * (21 + 23*len(numbers))
            print currenttemplate % tuple(args)
            print '-' * (21 + 23*len(numbers))
            currenttemplate = "%22s " + "   %7d %9.1f   " * len(numbers)
        args = [flowname]
        for streamlet in sorted(numbers.keys()):
            args += [round(numbers[streamlet][0], 1), round(numbers[streamlet][1], 1)]
        print currenttemplate % tuple(args)
    print
