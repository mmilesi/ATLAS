"""
File    : Core/__init__.py
Authors : KG <Kong.Guan.Tan@cern.ch>

Misc core functions of the framework
TODO: Show example of joboptions file setup here?
"""

print "@@@@@ Initialising the analysis framework"
import os, sys, copy, glob, optparse

filedir = os.path.dirname(os.path.abspath(__file__))
compiledir = 'libs'

def compileC(filename, target=None):
    global compiledir
    if not target is None:
        originaldir = compiledir
        compiledir = target
    from ROOT import gSystem
    if not gSystem.CompileMacro(filename, "k-", "", filedir+"/../"+compiledir):
        print "@@@@@ ERROR: Failed to compile", filename
        sys.exit(1)
    if not target is None:
        compiledir = originaldir

def makeVectorString(strlist):
    from ROOT import std
    vec = std.vector('string')()
    for s in strlist:
        vec.push_back(s)
    return vec

def parseInputArgs():
    parser = optparse.OptionParser(description='AnalysisFramework job configuration.')
    parser.add_option('-e', '--entries', default=None,
                      help='Start to end entries to run i.e. --entries=Start,End')
    parser.add_option('-o', '--outFile', default=None,
                      help='Name of output file')
    parser.add_option('-i', '--inFiles', default=None,
                      help='List of comma-separated input files')
    parser.add_option('-s', '--sampleName', default=None,
                      help='Name of the sample')
    parser.add_option('-p', '--useProof', default=None,
                      help='Use PROOF with number of cores (default or 0 = all cores on machine)')
    parser.add_option('-n', '--noProof', default=False, action='store_true',
                      help='Override all configuration and disable PROOF')
    parser.add_option('-f', '--stdinFiles', default=False, action='store_true',
                      help='Read input files from stdin')
    parser.add_option('-c', '--cutFlowSuppress', default=None,
                      help='Suppress comma separated streamlet groups (see Scripts/PrintCutFlow.py)')
    parser.add_option('-d', '--debugMode', default=False, action='store_true',
                      help='Turns debug mode on')

    (options, args) = parser.parse_args()
    return options

def listifyInputFiles(inputdir):
    if not type(inputdir) in [list, tuple]:
        inputdir = [inputdir]
    inputpath = []
    for i in inputdir:
        if '*' in i:
            inputpath += glob.glob(i)
        else:
            inputpath += [i]
    return inputpath

def checkInputs(inputpath):
    # Override the input files from commandline arguments if exists
    options = parseInputArgs()
    if options.inFiles:
        inputdir = options.inFiles.split(',')
    else:
        inputdir = inputpath

    # If stdin, then override (again)
    if options.stdinFiles:
        inputdir = sys.stdin.readline()[:-1].split(',') ## When running on grid!

    # Check if the input files exist and compile a list
    inputpath = listifyInputFiles(inputdir)
    if not inputpath:
        print "@@@@@ ERROR: No input files found!"
        sys.exit(1)
    else:
        inputpath = sorted(inputpath)
        print "@@@@@ Reading from the following file(s):"
        for i in inputpath:
            print '    ', i

    return inputpath

def checkAllInputs(outputpath, inputpath, samplename, datasets, entries, systematics, debugmode):
    inputpath = checkInputs(inputpath)

    # Override the sample name from commandlnie arguments if exists
    options = parseInputArgs()
    if options.sampleName:
        sampleinput = options.sampleName
    else:
        sampleinput = samplename
    samplename = sampleinput

    # Check sample name against the dataset manager to see if it's recognised
    from DatasetManager import DatasetManager
    datasets_paths = [ filedir+'/../Files/D3PD_datasets/' + d + '.txt' for d in datasets ]
    dm = DatasetManager(datasets_paths)
    # HACK: For now, this is manually set. It shouldn't!
    #dataset = DatasetManager([filedir+'/../Files/MC_D3PD_datasets.txt', filedir+'/../Files/Data_D3PD_datasets.txt', filedir+'/../Files/Embedding_D3PD_datasets.txt', filedir+'/../Files/MC_Win13_Melb_D3PD_datasets.txt'])

    if not dm.contains(samplename):
        print "@@@@@ ERROR: Sample name not recognised:", samplename
        sys.exit(1)
    else:
        print "@@@@@ Using sample name:", samplename
        print "     Tags:", dm.getTags(samplename)

    if options.outFile:
        outputpath = options.outFile

    if options.entries and len(options.entries.split(',')) == 2:
        entries = int(options.entries.split(',')[0]), int(options.entries.split(',')[1])

    if options.debugMode:
        debugmode = options.debugMode

    return outputpath, inputpath, samplename, datasets, entries, systematics, debugmode

def compileAll(
        branches            = [],
        treename            = 'tau',
        inputpath           = '',
        outputdump          = False,
        removeunselected    = False,
        cutflows            = [],
        toolset             = 'LepHad2011',
    ):
    global compiledir
    compiledir = 'libs/' + toolset

    from ROOT import gSystem
    gSystem.SetAclicMode(gSystem.kOpt)
    gSystem.AddIncludePath('-I' + filedir + '/../External/TauSpinner/extern/include')
    gSystem.AddIncludePath('-I' + filedir + '/../External/TauSpinner/extern/include/Tauola')
    compileC(filedir+"/Loader.C")
    compileC(filedir+"/CutFlowHist.C")

    print "@@@@@ Setting up the branches, kernel, whiteboard and external tools..."
    if type(branches) is str:
        branches = [branches]
    b = []
    for branch in branches:
        if not hasattr(__import__('Branches.'+branch), branch):
            raise Exception("Can't find Branches: " + branch)
        b.append(getattr(__import__('Branches.'+branch), branch))

    import CodeGenerator, NTupleTools
    if not CodeGenerator.setupBranches(b, treename, inputpath, outputdump):
        raise Exception("Setup branches failed!")

    compileC(filedir+"/WhiteBoard.C")
    compileC(filedir+"/CutFlow_Base.C")
    CodeGenerator.setupForkItem(removeunselected, outputdump)

    import External
    External.initialiseAllTools(toolset)

    import CutFlows
    CutFlows.load(cutflows)
    compileC(filedir+"/PostLoader.C")

def compileMinimal(compiletarget = 'libs/minimal'):
    global compiledir
    compiledir = compiletarget

    from ROOT import gSystem
    gSystem.SetAclicMode(gSystem.kOpt)
    compileC(filedir+"/CutFlowHist.C")

def setupApp(
        treename            = 'tau',
        samplename          = '',
        datasets            = [],
        inputpath           = '',
        outputpath          = 'ntuple.root',
        outputdump          = False,
        removeunselected    = False,
        corrections         = [],
        systematics         = [],
        debugmode           = False,
        
    ):

    # Make output file
    if treename.startswith('SystematicsUP/') or treename.startswith('SystematicsDOWN/'):
        recreate = False
    else:
        recreate = True
    NTupleTools.makeOutputFile(outputpath, recreate)

    # Check sample name against the dataset manager to see if it's recognised
    from DatasetManager import DatasetManager
    # HACK: See above HACK
    datasets_paths = [ filedir+'/../Files/D3PD_datasets/' + d + '.txt' for d in datasets ]
    dm = DatasetManager(datasets_paths)
#    dataset = DatasetManager([filedir+'/../Files/MC_D3PD_datasets.txt', filedir+'/../Files/Data_D3PD_datasets.txt', filedir+'/../Files/Embedding_D3PD_datasets.txt', filedir+'/../Files/MC_Win13_Melb_D3PD_datasets.txt'])
    if not dm.contains(samplename):
        print "@@@@@ ERROR: Sample name not recognised:", samplename
        sys.exit(1)

    from Kernel import Kernel
    kernel = Kernel()
    kernel.treename = treename
    kernel.inputpath = inputpath
    kernel.outputpath = outputpath
    kernel.outputdump = outputdump
    kernel.samplename = samplename

    from ROOT import AnalysisFramework
    kernel.ao = CodeGenerator.AllObjects()
    kernel.dataset = dm
    kernel.corrections = corrections
    kernel.systematics = systematics

    kernel.wb = AnalysisFramework.CutFlows.WhiteBoard()
    kernel.wb.debugMode = debugmode
    kernel.wb.doDump = not outputdump == False
    kernel.wb.doRemoveUnselected = removeunselected

    kernel.wb.isMC11a = dm.check(samplename, 'mc11a')
    kernel.wb.isMC11b = dm.check(samplename, 'mc11b')
    kernel.wb.isMC11c = dm.check(samplename, 'mc11c')
    kernel.wb.isMC12a = dm.check(samplename, 'mc12a')
    kernel.wb.isMC15  = dm.check(samplename, 'mc15')
    kernel.wb.isAFII_SS = dm.check(samplename, 'a188')
    kernel.wb.isAFII_fid  = dm.check(samplename, 'Fiducial')
    kernel.wb.isAFII_HSG4  = dm.check(samplename, 'FastSim')
    kernel.wb.isAFII  = kernel.wb.isAFII_SS or kernel.wb.isAFII_fid or kernel.wb.isAFII_HSG4
    #print 'isAFII', kernel.wb.isAFII
    kernel.wb.isMC = kernel.wb.isMC11a or kernel.wb.isMC11b or kernel.wb.isMC11c or kernel.wb.isMC12a or kernel.wb.isMC15 # or kernel.wb.isAFII
    kernel.wb.isEmbedding11 = dm.check(samplename, 'embedding')
    kernel.wb.isEmbedding12 = dm.check(samplename, 'embedding12')
    kernel.wb.isEmbedding13 = dm.check(samplename, 'embedding13')
    kernel.wb.isEmbedding   = kernel.wb.isEmbedding11 or kernel.wb.isEmbedding12 or kernel.wb.isEmbedding13
    kernel.wb.isVBFFiltered = dm.check(samplename, 'Z+VBF') or dm.check(samplename, 'Z+TightVBF')
    kernel.wb.isZDY = dm.check(samplename, 'Z+jets') or dm.check(samplename, 'DYZ+jets') or dm.check(samplename, 'AlpgenPythiaAutoZ+jets') or dm.check(samplename, 'DileptonFilteredLowmassDY') or kernel.wb.isVBFFiltered
    kernel.wb.isWAlpgenPythia = dm.check(samplename, 'AlpgenPythiaAutoW+jets')
    kernel.wb.isggF = dm.check(samplename, 'ggF')
    if dm.check(samplename, 'Signal'):
        kernel.wb.higgsMass = dm.getHiggsMass(samplename)
    kernel.wb.isData11 = dm.check(samplename, 'data')
    kernel.wb.isData12 = dm.check(samplename, 'data12')
    kernel.wb.isData15 = dm.check(samplename, 'data15')
    kernel.wb.isPeriodC = dm.check(samplename, 'PeriodC')
    kernel.wb.isPeriodD = dm.check(samplename, 'PeriodD')
    kernel.wb.isData = kernel.wb.isData11 or kernel.wb.isData12 or kernel.wb.isData15
    kernel.wb.isMuData = dm.check(samplename, 'Muons') and kernel.wb.isData
    kernel.wb.isElData = dm.check(samplename, 'Egamma') and kernel.wb.isData
    kernel.wb.isTauData = dm.check(samplename, 'JetTauEtmiss') and kernel.wb.isData

    return kernel
