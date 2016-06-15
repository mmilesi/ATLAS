"""
File    : MergeOutput.py
Authors : KG <Kong.Guan.Tan@cern.ch>

Script to merge output ntuples produced by the framework. Example usage:

    > python Scripts/MergeOutput.py - -s Files/samples*.csv -r /path/to/dq2downloadsdir -d /path/to/mergeddir

You can also provide a '-c' option to tell the script to merge the CutFlow histograms, as by default it doesn't to save time.

Additionally, you can provide a '-w' option to tell the script to use TTree.SetWeight to apply the cross-section weights listed in the samples.csv file, while also normalising it to the total number of events.
"""

from ROOT import gROOT
gROOT.SetBatch(True) # important to run without popups
from ROOT import TFile, TH1, TH1F, TH1D, TList, TObjString, TTree, TChain, TObjArray, TDirectoryFile
import sys, glob, os, optparse

sys.path.append(os.path.abspath(os.path.curdir))
from Core import NTupleTools, DatasetManager, listifyInputFiles

datasets = DatasetManager.DatasetManager()
samples = None

def parseInputArgs():
    parser = optparse.OptionParser(description='MergeOutput script configuration.')
    parser.add_option('-i', '--inFiles', default=None,
                      help='List of comma-separated input files')
    parser.add_option('-o', '--outFile', default=None,
                      help='Name of output file (use with --inFiles)')
    parser.add_option('-r', '--inRunDir', default=None,
                      help='Directory that contains downloaded samples (overrides --inFiles)')
    parser.add_option('-d', '--outRunDir', default=None,
                      help='Directory that contains combined samples (use with --inRunDir)')
    parser.add_option('-s', '--samplecsv', default='Files/samples2012.csv',
                      help='Specify the samples.csv file to use')
    parser.add_option('-c', '--cutflow', default=False, action='store_true',
                      help='Do merging of the CutFlow histograms, which can be a HUGE performance overhead')
    parser.add_option('-w', '--weight', default=False, action='store_true',
                      help='Apply cross-section weights to all TTree and normalise to total events')

    (options, args) = parser.parse_args()
    return options

def main():
    global samples

    # Defaults
    outputpath  = 'output/combined.root'
    inputpath   = 'output/ntuple*.root'

    from Core import compileMinimal
    compileMinimal()

    options = parseInputArgs()
    samples = datasets.getListSamples(options.samplecsv)
    if options.inRunDir:
        if not os.path.isdir(options.inRunDir):
            print("ERROR: input directory does not exist or is not a directory")
            return
        else:
            inputdir = options.inRunDir
        if not options.outRunDir:
            print("ERROR: invalid output directory (set with --outRunDir)")
            return
        elif os.path.isdir(options.outRunDir):
            print("ERROR:\noutput directory: \n{0} \nalready exists...to avoid inconsistencies, please remove it first".format(options.outRunDir))
            return
        else:
            outputdir = options.outRunDir
            os.makedirs(outputdir)
            logfile = open(outputdir + '/merge.log', 'w')

        for s in samples:
            sampledir = outputdir + '/' + s['group']
            if not os.path.isdir(sampledir):
                os.makedirs(sampledir)
            if not s['category'] == 'Data' and not s['group'] == 'Embedding':
                inputpath = inputdir + '/*' + s['ID'] + '*/*.root*'
            else:
                inputpath = inputdir + '/*' + s['name'] + '*/*.root*'

            separator = '.'
            if not s['ID']:
                separator = ''
            outputpath = sampledir + '/' + s['ID'] + separator + s['name'] + '.root'

            weight = None
            try:
                if options.weight:
                    weight = float(s['xsection']) * float(s['efficiency']) * float(s['kfactor']) * 1e3 # to get the weight in fb (the Xsec is in pb)
            except:
                pass
            mergeOne(inputpath, outputpath, logfile, weight, options.cutflow)
    else:
        if options.inFiles:
            inputpath = options.inFiles

        if options.outFile:
            outputpath = options.outFile

        mergeOne(inputpath, outputpath, cutflow=options.cutflow)

def mergeOne(inputpath, outputpath, logfile=None, weight=None, cutflow=True):
    print("Merging {0} ...\n".format(inputpath)),
    original_inputpath = inputpath
    inputpath = listifyInputFiles(inputpath)
    if not inputpath:
        print("ERROR: No inputs here specified!")
        if logfile:
            logfile.write("ERROR: No inputs found for " + original_inputpath + "\n")
        return
    else:
        missingfiles = []
        for i in inputpath:
            if not os.path.isfile(i):
                missingfiles.append(i)
        if missingfiles:
            print("ERROR: File(s) not found:", ', '.join(missingfiles))
            if logfile:
                logfile.write("ERROR: Missing input files for " + original_inputpath + ":\n")
                for m in missingfiles:
                    logfile.write('    ' + m + '\n')
            return

    target = TFile.Open(outputpath, "RECREATE")
    path = target.GetPath()
    path = path[path.index(':')+2:]

    cache={'TOTALLUMI':0}
    errorfiles = []

    chain = TChain("physics")
    for i in inputpath:
        f = TFile.Open(i)
        if not f or f.IsZombie():
            errorfiles.append(i)
            continue
	# Take tree in file and add it to TChain
	chain.Add(i)
        print("\nMerging input file: \n{0} \nto target file...".format(i))
        recursiveMerge(target, f, path, cache, cutflow)
        f.Close()

    if errorfiles:
        print("ERROR in opening the following files:")
        for e in errorfiles:
            print "    ", e
        if logfile:
            logfile.write("ERROR: Cannot open input files for " + original_inputpath + ":\n")
            for e in errorfiles:
                logfile.write('    ' + e + '\n')

    if weight:
        totalevents = None
        for key in cache:
            obj = cache[key]
            if ( type(obj) == TH1D or type(obj) == TH1F ) and obj.GetName() == 'TotalEventsW':
                totalevents = obj.GetBinContent(2)
                break

        if totalevents:
            for key in cache:
                obj = cache[key]
                if type(obj) == TTree:
                    obj.SetWeight(weight/totalevents)

            print("Applying weight - w = ( xsec[fb]*kfactor*filter_eff) / (TOT EVTS W ) = {0}/{1} = {2} [fb]\n".format(weight,totalevents,(weight/totalevents)))

    target.Write()

    del cache
    print("Merged {0} files into {1}".format(len(inputpath), outputpath))

    merged_tree = target.Get("physics")

    if ( merged_tree and chain.GetEntries() != merged_tree.GetEntries() ):
        print("******** WARNING! ********\n entries in chain: {0} \n entries in merged TTree: {1} \n*************************".format(chain.GetEntries(), merged_tree.GetEntries()))

    target.Close()

def recursiveMerge(target, infile, path='', cache={'TOTALLUMI':0}, cutflow=True):
    l = infile.GetDirectory(path)
    keys = l.GetListOfKeys()
    cycles = {}

    #print("keys in input file: \n\n{0}\n\n".format(keys.ls()))

    for entry in range(keys.GetEntries()):
    	name = keys.At(entry).GetName() + ";" + str(keys.At(entry).GetCycle())
    	if path:
    	    cachename = path + "/" + name
    	else:
    	    cachename = name
    	obj = l.Get(name)

    	if type(obj) == TDirectoryFile:
    	    #print("TDirectory obj name: {0}".format(obj.GetName()))
    	    targetpath = keys.At(entry).GetName()
    	    if not target.Get(targetpath):
    		target.mkdir(targetpath)
    	    recursiveMerge(target, infile, path + "/" + obj.GetName(), cache)
    	elif type(obj) == TTree:
    	    #print("TTree obj name: {0} - cachename: {1} ".format(obj.GetName(), cachename))
    	    cyclename, cyclenumber = cachename.split(';')
    	    if cyclename in cycles: continue
    	    #print("cyclename: {0} - cyclenumber: {1}".format(cyclename, cyclenumber))
    	    cycles[cyclename] = cyclenumber
    	    if not cyclename in cache:
    		#print("adding cyclename {0} to cache (via TTree::CloneTree())".format(cyclename))
    		target.cd(path)
    		cache[cyclename] = obj.CloneTree()
    	    else:
    		objcached = cache[cyclename]
    		col = TObjArray()
    		col.Add(obj)
    		#print("merging TTree obj to cached object")
    		objcached.Merge(col)
    	elif issubclass(obj.__class__, TH1):
    	    #print("TH1 obj name: {0}".format(obj.GetName()))
    	    if not cutflow and keys.At(entry).GetName() == "CutFlow":
    		continue
    	    if not cachename in cache:
    		target.cd(path)
    		cache[cachename] = obj.Clone()
    	    else:
    		objcached = cache[cachename]
    		col = TObjArray()
    		col.Add(obj)
    		objcached.Merge(col)
    	elif type(obj) == TObjString:
    	    #print("TObjString obj name: {0}".format(obj.GetName()))
    	    if obj:
    		target.cd(path)
    		objnew = TObjString(obj.GetString().Data())
    		objnew.Write(keys.At(entry).GetName())
    		cache['TOTALLUMI'] += 1
    	elif issubclass(obj.__class__, TList):
    	    #print("TList obj name: {0}".format(obj.GetName()))
    	    if obj:
    		target.cd(path)
    		objnew = TList(obj)
    		objnew.Write(keys.At(entry).GetName()) # not working...
    	else:
    	    print "UNKNOWN OBJECT", name, "OF TYPE", type(obj)

main()
