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
    parser.add_option('-r', '--inRunDir', default=None,
                      help='Directory that contains downloaded samples (overrides --inFiles)')
    parser.add_option('-s', '--samplecsv', default='Files/samples2012.csv',
                      help='Specify the samples.csv file to use')
 
    (options, args) = parser.parse_args()
    return options

def main():

    # Defaults
    inputpath   = 'output/ntuple*.root'

    from Core import compileMinimal
    compileMinimal()

    options = parseInputArgs()

    samples = datasets.getListSamples(options.samplecsv)

    if options.inRunDir:
        if not os.path.isdir(options.inRunDir):
            print "ERROR: input directory does not exist or is not a directory"
            return
        else:
            inputdir = options.inRunDir
            
	inputpath = inputdir + '/*/*.root*'

	checkDuplicates(inputpath, samples)
	    
    else:
        if options.inFiles:
            inputpath = options.inFiles

	checkDuplicates(inputpath, samples)

def checkDuplicates(inputpath, samples):
    
    original_inputpath = inputpath
    inputpath = listifyInputFiles(inputpath)
    
    if not inputpath:
        print "ERROR: No inputs here specified!"
        return
    else:
        missingfiles = []
        for i in inputpath:
            if not os.path.isfile(i):
                missingfiles.append(i)
        if missingfiles:
            print "ERROR: File(s) not found:", ', '.join(missingfiles)
            return
    
    cache={'TOTALLUMI':0}
    errorfiles = []

    for i in inputpath:
        f = TFile.Open(i)
        if not f or f.IsZombie():
            errorfiles.append(i)
            continue
        
	print "Checking file ", i, "..."
	
	printSampleID(i, samples)
	
	totalevents_EL, totalevents_MetaData, totalweights_MetaData = check(f)
	
	if (totalevents_EL != totalevents_MetaData):
	   print "\t WARNING!\n\t total events from Metadata = ", totalevents_MetaData, "\n\t total sumOfWeights from Metadata = ", totalweights_MetaData,"\n\t total events from EL = ", totalevents_EL, "\n"
	else:
	   print "All good with this sample!\n"
	
	f.Close()

    if errorfiles:
        print "ERROR in opening the following files:"
        for e in errorfiles:
            print "    ", e
    
def printSampleID(infilename, samples):
	
	actualfilename = (infilename.split('/'))[-1]
	actualfilename = actualfilename.replace('.root','')
	
	for s in samples:
	   if actualfilename == s['name']:
		print "\t Sample ID = ", s['ID']	
            
    
def check(infile):
        l = infile.GetDirectory("")
        keys = l.GetListOfKeys()
	
        #print 'keys in input file: \n', keys.ls(), '\n'
        
	for entry in range(keys.GetEntries()):
            name = keys.At(entry).GetName() + ";" + str(keys.At(entry).GetCycle())
            cachename = name
            obj = l.Get(name)

            if issubclass(obj.__class__, TH1):
                
		#print '\t\tTH1 obj name: ', obj.GetName()
                
		if obj.GetName() == 'EventLoop_EventCount':
		    totalevents_EL = obj.GetBinContent(1)
		
		if obj.GetName() == 'TotalEvents':
		    totalevents_MetaData = obj.GetBinContent(1)
		
		if obj.GetName() == 'TotalEventsW':
		    totalweights_MetaData = obj.GetBinContent(1)		
		
        return totalevents_EL, totalevents_MetaData, totalweights_MetaData

main()
