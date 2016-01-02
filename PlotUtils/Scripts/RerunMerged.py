import os, sys, glob, subprocess, optparse

sys.path.append(os.path.abspath(os.path.curdir))
from Core import NTupleTools, parseInputArgs, listifyInputFiles, DatasetManager

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
    parser.add_option('-j', '--joboptions', default='DefaultCutFlow_Rerun.py',
                      help='Specify which job options to use for the re-run')

    (options, args) = parser.parse_args()
    return options


def main():
    global samples

    # Defaults
    outputpath  = 'output/combined.root'
    inputpath   = 'output/ntuple*.root'
    #joboptions = options.joboptions
    samplename = "data12_8TeV.periodA.physics_Muons.PhysCont.NTUP_EMBLHIM.grp14_v02_r4697_p1462"
    
    from Core import compileMinimal
    compileMinimal()

    options = parseInputArgs()
    samples = datasets.getListSamples(options.samplecsv)
    joboptions  = options.joboptions
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

        for s in samples:
            sampledir = outputdir + '/' + s['group']
            if not os.path.isdir(sampledir):
                os.makedirs(sampledir)
            inputpath = inputdir + '/' + s['group'] + '/' + s['name'] + '.root'
            outputpath = sampledir + '/' + s['name'] + '.root'
            rerunOne(inputpath, outputpath, samplename, joboptions)
    else:
        if options.inFiles:
            inputpath = options.inFiles

        if options.outFile:
            outputpath = options.outFile

        rerunOne(inputpath, outputpath, samplename, joboptions)

def rerunOne(inputpath, outputpath, samplename, joboptions):
    inputpath = [os.path.abspath(p) for p in listifyInputFiles(inputpath)]
    inputpath = '"' + ','.join(inputpath) + '"'
    command = ' '.join(['python', joboptions, '-', '-n', '-s', samplename, '-i', inputpath, '-o', outputpath])
    #subprocess.call(['python', joboptions, '-', '-i', inputpath, '-o', outputpath])
    os.system(command)

main()

