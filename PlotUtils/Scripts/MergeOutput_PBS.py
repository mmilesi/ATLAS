"""
File    : MergeOutput.py
Authors : KG <Kong.Guan.Tan@cern.ch>

Script to merge output ntuples produced by the framework. Example usage:

    > python Scripts/MergeOutput.py - -s Files/samples2012.csv -r /path/to/dq2downloadsdir -d /path/to/mergeddir

You can also provide a '-c' option to tell the script to merge the CutFlow histograms, as by default it doesn't to save time.

Additionally, you can provide a '-w' option to tell the script to use TTree.SetWeight to apply the cross-section weights listed in the samples.csv file, while also normalising it to the total number of events.
Attention with the --inFiles option in place of --inRunDir the output tree will not be weighted!
"""

from ROOT import TFile, TH1, TH1D, TObjString, TTree, TChain, TObjArray, TDirectoryFile
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
    parser.add_option('-s', '--samplecsv', default='Files/samples.csv',
                      help='Specify the samples.csv file to use')
    parser.add_option('-t', '--singleOutput', default='physics_Muons',
                      help='Specify the singleOutput file to use')
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
    inputpath    = 'output/ntuple*.root'

    from Core import compileMinimal
    compileMinimal()

    options = parseInputArgs()
    samples = datasets.getListSamples(options.samplecsv)
    
    print "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
    print "Running over singleOutput: %s" % (options.singleOutput)
    print "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
    
    #Here checks if input and output folders are valid
    if options.inRunDir:
        if not os.path.isdir(options.inRunDir):
            print "ERROR: input directory does not exist or is not a directory"
            return
        else:
            inputdir = options.inRunDir
        if not options.outRunDir:
            print "ERROR: invalid output directory (set with --outRunDir)"
            return
        elif os.path.isdir(options.outRunDir):
            print "ERROR: output directory already exists...to avoid inconsistencies, please remove it first"
            return
        else:
            outputdir = options.outRunDir
            os.makedirs(outputdir)
            logfile = open(outputdir + '/merge.log', 'w')

		#Here creates the structure in subdirectories with the name of the group
        for s in samples:
            sampledir = outputdir + '/' + s['group']
            if not os.path.isdir(sampledir):
                os.makedirs(sampledir)
            if not s['category'] == 'Data' and not s['group'] == 'Embedding':
				#For MC samples the number and the name are used 
                inputpath = inputdir + '/*' + s['ID'] + '*' + s['name'] + '*/*.root*'
            else:
				#For data and embedding only the name
                inputpath = inputdir + '/*' + s['name'] + '*/*.root*'
            outputpath = sampledir + '/' + s['name'] + '.root'
            weight = None
            try:
                if options.weight:
                    weight = float(s['xsection']) * float(s['efficiency']) * float(s['kfactor']) * 1.0e6
            except:
                pass
                weight = None
			#Here start the difference with the standard tools. It merges the file only if it is specified options.singleOutput and coincides with the output ntuple name otherwise it does nothing
            if options.singleOutput in outputpath.split('/')[-1]:
               print "MERGING: %s" % (outputpath.split('/')[-1])
#               print "MERGING: %s" % (outputpath)
               mergeOne(inputpath, outputpath, logfile, weight, options.cutflow)
            else:
              print "SKIPPING: %s" % (outputpath.split('/')[-1])
#              print "SKIPPING: %s" % (outputpath)
    else:
		#This part of code should not be used when merging in pbs
        if options.inFiles:
            inputpath = options.inFiles

        if options.outFile:
            outputpath = options.outFile
        if options.singleOutput in outputpath.split('/')[-1]:
          print "MERGING: %s" % (outputpath.split('/')[-1])
#          print "MERGING: %s" % (outputpath)
		  #Why the weight is not specified in this case?
          mergeOne(inputpath, outputpath, cutflow=options.cutflow)
        else:
          print "SKIPPING: %s" % (outputpath.split('/')[-1])
#          print "SKIPPING: %s" % (outputpath)

def mergeOne(inputpath, outputpath, logfile=None, weight=None, cutflow=True):
    print "Merging", inputpath, "...",
    original_inputpath = inputpath
    inputpath = listifyInputFiles(inputpath)#it solves the wildcards of inputpath in a list of files
    if not inputpath:
        print "ERROR: No inputs here specified!"
        if logfile:
            logfile.write("ERROR: No inputs found for " + original_inputpath + "\n")
        return
    else:
        missingfiles = []
        for i in inputpath:
            if not os.path.isfile(i):
                missingfiles.append(i)
        if missingfiles:
            print "ERROR: File(s) not found:", ', '.join(missingfiles)
            if logfile:
                logfile.write("ERROR: Missing input files for " + original_inputpath + ":\n")
                for m in missingfiles:
                    logfile.write('    ' + m + '\n')
            return

    target = TFile.Open(outputpath, "RECREATE")#opening output file
    path = target.GetPath()
    path = path[path.index(':')+2:]

    cache={'TOTALLUMI':0}
    errorfiles = []
	#loop over input files
    for i in inputpath:
        f = TFile.Open(i)
        if not f or f.IsZombie():
            errorfiles.append(i)
            continue
        #print i
        recursiveMerge(target, f, path, cache, cutflow)#function used to merge the files
        f.Close()
    if errorfiles:
        print "ERROR in opening the following files:"
        for e in errorfiles:
            print "    ", e
        if logfile:
            logfile.write("ERROR: Cannot open input files for " + original_inputpath + ":\n")
            for e in errorfiles:
                logfile.write('    ' + e + '\n')

	#Setting the ttree weight
    if weight:
        totalevents = None
        for key in cache:
            obj = cache[key]
            if type(obj) == TH1D and obj.GetName() == 'TotalEvents':
                totalevents = obj.GetBinContent(2)
                break

        if totalevents:
            for key in cache:
                obj = cache[key]
                if type(obj) == TTree:
                    obj.SetWeight(weight / totalevents)

    target.Write()
    target.Close()

    del cache
    print "Merged", len(inputpath), "files into", outputpath

def recursiveMerge(target, infile, path='', cache={'TOTALLUMI':0}, cutflow=True):
        l = infile.GetDirectory(path)
        keys = l.GetListOfKeys()
        cycles = {}
        for entry in range(keys.GetEntries()):
            name = keys.At(entry).GetName() + ";" + str(keys.At(entry).GetCycle())
            if path:
                cachename = path + "/" + name
            else:
                cachename = name
            obj = l.Get(name)
            if type(obj) == TDirectoryFile:
                #print obj, "DIRECTORY"
                targetpath = keys.At(entry).GetName()
                if not target.Get(targetpath):
                    target.mkdir(targetpath)
                recursiveMerge(target, infile, path + "/" + obj.GetName(), cache)
            elif type(obj) == TTree:
#                print obj, cachename, "TTree"
                cyclename, cyclenumber = cachename.split(';')
                if cyclename in cycles: continue
#                print cachename, "Used!"
                cycles[cyclename] = cyclenumber
                if not cyclename in cache:
                    target.cd(path)
                    cache[cyclename] = obj.CloneTree()
                else:
                    objcached = cache[cyclename]
                    col = TObjArray()
                    col.Add(obj)
                    objcached.Merge(col)
            elif issubclass(obj.__class__, TH1):
                #print obj, "TH1"
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
                #print type(obj), name, "TObjString"
                if obj:
                    target.cd(path)
                    objnew = TObjString(obj.GetString().Data())
                    objnew.Write(keys.At(entry).GetName())
                    cache['TOTALLUMI'] += 1
            else:
                print "UNKNOWN OBJECT", name, "OF TYPE", type(obj)

main()
