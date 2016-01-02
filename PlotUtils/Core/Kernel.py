"""
File    : Core/Kernel.py
Authors : KG <Kong.Guan.Tan@cern.ch>

Sets up the cutflow chain, event loop, and runs it
TODO:   - List out not-so-obvious features for debugging
        - Some comments explaning what's going on
        - Fix experimental PROOF
"""

import sys, glob, math, time, os

from Core import makeVectorString, NTupleTools, checkInputs, parseInputArgs
from ROOT import TNamed, TH1D, TFile, TObjString, TProof, TPySelector, TDSet, TTree

class Kernel(TNamed):
    def __init__(self):
        TNamed.__init__(self, 'Kernel', 'Kernel')
        self.FlowList = []
        self.wb = None
        self.ao = None
        self.dataset = None
        self.treename = None
        self.samplename = None
        self.inputpath = None
        self.outputpath = None
        self.outputdump = None
        self.corrections = []
        self.systematics = []
        self.totalflow = None
        self.inputtree = None
        self.outputtree = None
        self._tempVec = []
        self._tempTree = {}

    def getDefaultOutputTree(self):
        if not self.outputtree:
            self.outputtree = NTupleTools.makeNTuple(self.treename, self.ao, self.outputdump)
        return self.outputtree

    def setupFlowChain(self, flowlist):
        from ROOT import AnalysisFramework

        systematicsPartners = {}
        systematicsList = []
        systematicsOff = []
        systematicsOn = []
        for syslist in self.systematics:
            if not type(syslist) in [list, tuple]:
                syslist = [syslist]
            else:
                friendflows = [flow for flow in flowlist if flow.name in syslist]
                for flow in friendflows:
                    systematicsPartners[flow.name] = [s for s in syslist if not s == flow.name]
                    sharedVec =  makeVectorString(systematicsPartners[flow.name])
                    flow.systematicsPartners = sharedVec
                    self._tempVec.append(sharedVec)
            for s in syslist:
                systematicsList.append(s)
        self.wb.systematics = makeVectorString(systematicsList)
        for flow in flowlist:
            if issubclass(flow.__class__, AnalysisFramework.CutFlows.SystematicsItem):
                if flow.name in systematicsList:
                    systematicsOn.append(flow.name)
                    if (self.wb.isData and flow.doData) or (self.wb.isMC and flow.doMC) or (self.wb.isEmbedding and flow.doEmbedding):
                        if flow.name in systematicsPartners:
                            treename = '_'.join(sorted([flow.name] + systematicsPartners[flow.name]))
                        else:
                            treename = flow.name
                        flow.treeDownDelayed = self._tempTree.get("SystematicsDOWN/"+treename, NTupleTools.makeNTuple("SystematicsDOWN/"+treename, self.ao, self.outputdump))
                        flow.treeUpDelayed = self._tempTree.get("SystematicsUP/"+treename, NTupleTools.makeNTuple("SystematicsUP/"+treename, self.ao, self.outputdump))
                        self._tempTree["SystematicsDOWN/"+treename] = flow.treeDownDelayed
                        self._tempTree["SystematicsUP/"+treename] = flow.treeUpDelayed
                else:
                    if not flow.doNominalWhenOff:
                        flow.doNominalInternal = False
                        flow.doNominalExternal = False
                    systematicsOff.append(flow.name)
        print '@@@@@ The following systematics are enabled:'
        for syst in systematicsOn:
            if syst in systematicsPartners: print '    ', syst, 'with', ', '.join(systematicsPartners[syst])
            else: print '    ', syst
        if not systematicsOn:
            print '     None!'
        print "@@@@@ The following systematics are disabled:"
        for syst in systematicsOff:
            if syst in systematicsPartners: print '    ', syst, 'with', ', '.join(systematicsPartners[syst])
            else: print '    ', syst
        if not systematicsOff:
            print '     None!'
        for syst in systematicsList:
            if not syst in systematicsOff and not syst in systematicsOn:
                if syst in systematicsPartners: print '     WARNING: Unrecognised systematic:', syst, 'with', ', '.join(systematicsPartners[syst])
                else: print '     WARNING: Unrecognised systematic:', syst

        toRemove = []
        correctionsOff = []
        correctionsOn = []
        for flow in flowlist:
            if issubclass(flow.__class__, AnalysisFramework.CutFlows.CorrectionItem):
                if flow.name in self.corrections:
                    correctionsOn.append(flow.name)
                else:
                    correctionsOff.append(flow.name)
                    toRemove.append(flow)
        for flow in toRemove:
            flow.wb = self.wb
            flow.ao = self.ao
            flowlist.remove(flow)
        print "@@@@@ The following corrections are enabled:"
        for corr in correctionsOn:
            print '    ', corr
        if not correctionsOn:
            print '     None!'
        print "@@@@@ The following corrections are disabled:"
        for corr in correctionsOff:
            print '    ', corr
        if not correctionsOff:
            print '     None!'
        for corr in self.corrections:
            if not corr in correctionsOff and not corr in correctionsOn:
                print "     WARNING: Unrecognised correction:", corr
        self.wb.corrections = makeVectorString(correctionsOn)

        for flow in flowlist:
            if flow.name == "TotalEvents":
                self.totalflow = flow
                break
        if not self.totalflow:
            print "@@@@@ ERROR: Must have a flow named 'TotalEvents' to run!"
            exit()

        for i in range(len(flowlist)):
            flow = flowlist[i]
            try:
                flows = list(flow.nextList)
                flow.wb = self.wb
                flow.ao = self.ao
            except:
                flows = [flow]
            for flow in flows:
                if i==0:
                    flow.hasPrevious = False
                else:
                    prevflow = flowlist[i-1]
                    flow.previous = prevflow
                    flow.hasPrevious = True
                if i+1<len(flowlist):
                    nextflow = flowlist[i+1]
                    flow.next = nextflow
                    flow.hasNext = True
                else:
                    flow.hasNext = False
                flow.wb = self.wb
                flow.ao = self.ao
        return flowlist[0].process

    def getCutFlowFromFlowChain(self, flowlist):
        from ROOT import AnalysisFramework

        cutflow = []
        for i in range(len(flowlist)):
            flow = flowlist[i]
            try:
                flows = list(flow.nextList)
            except:
                flows = [flow]
            atleastone = [flow for flow in flows if issubclass(flow.__class__, AnalysisFramework.CutFlows.CutItem)]
            if not atleastone: continue
            numbers = {}
            for flow in flows:
                for streamlet in flow.streamlets:
                    totalEvents = flow.totalEvents[streamlet]
                    if not streamlet:
                        streamlet = "All"
                    numbers[streamlet] = (totalEvents.passed, totalEvents.passedW)
            if not numbers:
                numbers["All"] = (0, 0.0)
            cutflow.append((flow.name, numbers))
        return cutflow

    def run(self, entryNumbers = ()):
        if not self.FlowList:
            raise Exception("Nothing defined in CutFlow!")
        firstflow_process = self.setupFlowChain(self.FlowList)

        self.inputtree = NTupleTools.loadNTuple(self.treename, self.inputpath)
        NTupleTools.setBranchAddresses(self.inputtree, self.ao)
        num = self.inputtree.GetEntries()
        if entryNumbers:
            start, end = entryNumbers
            if end > num:
                end = num
            allentries = xrange(start, end)
        else:
            allentries = xrange(num)
        printmod = int(round(num/1000))
        wb_resetToDefaultValues = self.wb.resetToDefaultValues
        ao_resetToDefaultValues = self.ao.resetToDefaultValues
        ao_forceResetToDefaultValues = self.ao.forceResetToDefaultValues
        inputtree_GetEntry = self.inputtree.GetEntry
        inputtree_GetTreeNumber = self.inputtree.GetTreeNumber
        ao_resetSelected = self.ao.resetSelected
        outputdump = self.outputdump
        ao_removeUnselected = self.ao.removeUnselected
        outputtree_Fill = self.outputtree.Fill
        sys_stdout_flush = sys.stdout.flush

        print
        print "Number of input entries:", num
        try:
            current_tree_num = 0
            for i in allentries:
                wb_resetToDefaultValues()
                ao_resetToDefaultValues()
                if not inputtree_GetTreeNumber() == current_tree_num:
                    current_tree_num = inputtree_GetTreeNumber()
                    print "Switching to tree number %i in TChain..." % (current_tree_num)
                    NTupleTools.setBranchAddresses(self.inputtree, self.ao)
                    ao_forceResetToDefaultValues()
                inputtree_GetEntry(i)
                ao_resetSelected()
                if printmod and i % printmod == 0:
                    print "\r%.1f%%" % (i*100./num),
                    sys_stdout_flush()
                passflow = firstflow_process()
        except Exception, e:
            crashedflow = None
            for fl in range(len(self.FlowList)):
                flow = self.FlowList[fl]
                try:
                    flows = [flow] + list(flow.nextList)
                except:
                    flows = [flow]
                for flow in flows:
                    if flow.running:
                        crashedflow = flow
                        break
            print "-" * 80
            if crashedflow:
                print "@@@@@ Crash detected at entry", i, "in", crashedflow.name, crashedflow
                print "        Passport :", ' > '.join([fl.name for fl in crashedflow.runningPassport.flowHistory])
                print "      Streamlets :", ' > '.join(list(crashedflow.runningPassport.streamletHistory))
            else:
                print "@@@@@ Crash detected at entry", i, "but unable to determine which flow"
            print "    ErrorMessage :", e
            print "@@@@@ Starting interactive session so you can inspect your variables...use Ctrl-D to exit"
            print
            import code
            vars = globals().copy()
            vars.update(locals())
            shell = code.InteractiveConsole(vars)
            shell.interact()
            sys.exit(1)

        self.outputtree.SetWeight(self.inputtree.GetWeight())
        print " Done!"
        print
        print "Setting tree weight to", self.outputtree.GetWeight()
        print "Number of entries to write in default tree:", self.outputtree.GetEntries()

        if self.treename.startswith('SystematicsUP/') or self.treename.startswith('SystematicsDOWN/'):
            # Do the write operation (might take a while if the tree is huge)
            self.outputtree.Write()
            return

        if self._tempTree:
            print "Number of entries to write in other trees:"
        for t in sorted(self._tempTree.keys()):
            self._tempTree[t].SetWeight(self.inputtree.GetWeight())
            print '     %30s : ' % (t), self._tempTree[t].GetEntries()

        # Make (or copy if exist) the total events histogram
        totalEventsHistogram = NTupleTools.getTotalEventsHistogram(self.inputpath)
        if totalEventsHistogram:
            print "TotalEvents histogram already exist. Using it instead."
        else:
            totalEventsHistogram = NTupleTools.makeTotalEventsHistogram(self.totalflow)
            print "Creating TotalEvents histogram. Bin 1 = raw, Bin 2 = weighted."
        print "Sample name (" + self.samplename + ") stored in TotalEvents.GetTitle()"
        totalEventsHistogram.SetTitle(self.samplename)

        # Save the cutflow
        print "Creating CutFlow histogram. Bins labeled as FlowNumber/FlowName/Streamlet."
        cutflow = self.getCutFlowFromFlowChain(self.FlowList)
        cutFlowHist = NTupleTools.makeCutFlowHistogram(cutflow)

        # Copy the Lumi XML strings across to the new ntuple
        NTupleTools.copyLumi(self.inputpath)

        # Do the write operation (might take a while if the tree is huge)
        NTupleTools.outputFile.Write()

        # Move the output file to the output location
        options = parseInputArgs()
        if not options.noProof:
            dirname = os.path.dirname(self.outputpath)
            basename = '.temp.' + os.path.basename(self.outputpath)
            if dirname and not os.path.isdir(dirname):
                os.makedirs(dirname)
            os.rename(basename, self.outputpath)


#def RunProof(treename, inputpath, outputpath, entries = (), cores = None):
#    inputpath = checkInputs(inputpath)
#    if not cores:
#        proof = TProof.Open('')
#    else:
#        proof = TProof.Open('workers=%d' %(cores) )
#    proof.Exec ('TPython::Exec("%s");' % \
#        ("import sys; sys.path.insert(0,'"+os.path.dirname(os.path.abspath("Core"))+"')"))
#    TProof.AddEnvVar ("PYTHONHOME", sys.prefix+":"+sys.exec_prefix)
#    tree = TDSet('TTree', treename)
#    for l in inputpath:
#        tree.Add(l)
#    time.sleep(1)
#    if entries:
#        tree.Process( 'TPySelector', 'Core.Kernel',  entries[1]-entries[0], entries[0])
#    else:
#        tree.Process( 'TPySelector', 'Core.Kernel')
#
#    # Move the output file to the output location
#    options = parseInputArgs()
#    if not options.noProof:
#        dirname = os.path.dirname(outputpath)
#        basename = '.temp.' + os.path.basename(outputpath)
#        if dirname and not os.path.isdir(dirname):
#            os.makedirs(dirname)
#        os.rename(basename, outputpath)
#
#    # Reopen output files to write the total events and copy the lumi files
#    NTupleTools.outputFile = TFile.Open(outputpath)
#    outputtree = NTupleTools.outputFile.Get(treename)
#    print
#    print "Number of entries to write in default tree:", outputtree.GetEntries()
#    NTupleTools.outputFile = TFile.Open(outputpath, 'UPDATE')
#
#    # Make (or copy if exist) the total events histogram
#    totalEventsHistogram = NTupleTools.getTotalEventsHistogram(inputpath)
#    if totalEventsHistogram:
#        print "TotalEvents histogram already exist. Using it instead."
#    else:
#        print "Creating TotalEvents histogram. Bin 1 = raw, Bin 2 = weighted."
#    print "Creating CutFlow histogram. Bins labeled as FlowNumber/FlowName/Streamlet."
##    print "Sample name (" + self.samplename + ") stored in TotalEvents.GetTitle()"
##    totalEventsHistogram.SetTitle(self.samplename)
#
#    # Copy the Lumi XML strings across to the new ntuple
#    NTupleTools.copyLumi(inputpath)
#
#    # Do the write operation (might take a while if the tree is huge)
#    NTupleTools.outputFile.Write()
#
#class MyPySelector(TPySelector):
#    def Begin( self ):
#        print '@@@@@ PROOF: Begin'
#
#    def SlaveBegin( self, tree ):
#        print '@@@@@ PROOF: Slave Begin'
#        import DefaultCutFlow_Example
#        self.kernel = DefaultCutFlow_Example.ConfigureJob()
#        if not self.kernel.FlowList:
#            raise Exception("Nothing defined in CutFlow!")
#        self.firstflow_process = self.kernel.setupFlowChain(self.kernel.FlowList)
#        self.ao_resetToDefaultValues = self.kernel.ao.resetToDefaultValues
#        self.wb_resetToDefaultValues = self.kernel.wb.resetToDefaultValues
#        self.ao_resetSelected = self.kernel.ao.resetSelected
#
#    def Init( self, tree):
#        NTupleTools.setBranchAddresses(self.fChain, self.kernel.ao)
#        self.GetEntry = self.fChain.GetEntry
#
#    def Process( self, entry ):
#        self.ao_resetToDefaultValues()
#        self.wb_resetToDefaultValues()
#        if self.GetEntry( entry ) <= 0: return 0
#        self.ao_resetSelected()
#        self.firstflow_process()
#        return 1
#
#    def SlaveTerminate( self ):
#        print '@@@@@ PROOF: Slave Terminating'
#        print " Done!"
#        print
#        print "Number of entries to write in default tree:", self.kernel.outputtree.GetEntries()
#        if self.kernel._tempTree:
#            print "Number of entries to write in other trees:"
#        for t in self.kernel._tempTree:
#            print '     %20s : ' % (t), self.kernel._tempTree[t].GetEntries()
#        print
#
#        # Make (or copy if exist) the total events histogram
#        totalEventsHistogram = NTupleTools.getTotalEventsHistogram(self.kernel.inputpath)
#        if totalEventsHistogram:
#            print "TotalEvents histogram already exist. Using it instead."
#            del totalEventsHistogram
#        else:
#            totalEventsHistogram = NTupleTools.makeTotalEventsHistogram(self.kernel.totalflow, self.samplename)
#            print "Creating TotalEvents histogram. Bin 1 = raw, Bin 2 = weighted."
#        print "Sample name (" + self.kernel.samplename + ") stored in TotalEvents.GetTitle()"
#        totalEventsHistogram.SetTitle(self.kernel.samplename)
#
#        # Save the cutflow
#        cutflow = self.kernel.getCutFlowFromFlowChain(self.kernel.FlowList)
#        cutFlowHist = NTupleTools.makeCutFlowHistogram(cutflow)
#        print "Creating CutFlow histogram. Bins labeled as FlowNumber/FlowName/Streamlet."
#
#        NTupleTools.printCutFlow(cutflow, suppressStreamlet=['SS', 'OtherMT', 'HighMT'])
#        NTupleTools.outputFile.Write()
#        NTupleTools.outputFile.Close()
#        NTupleTools.outputProofFile.Print()
#        self.GetOutputList().Add(NTupleTools.outputProofFile)
#
#    def Terminate( self ):
#        print '@@@@@ PROOF: Terminating'
