""" BackgroundsTools.py: base classes to define variables, cuts, processes operations and load input samples """

__author__     = "KG, Marco Milesi, Francesco Nuti"
__email__      = "Kong.Guan.Tan@cern.ch, marco.milesi@cern.ch, francesco.nuti@cern.ch"
__maintainer__ = "Marco Milesi"

import sys, glob, os, array, inspect, math

from ROOT import TFile, TH1, TH1D, TH1F, TH1I, TH2, TH2D, TH2F, TH2I, TObjString, TTree, TChain, TObjArray, TDirectoryFile, TNamed, TObject
from ROOT import gROOT, gPad, THStack, TColor, TCanvas, TPad, TLine, TLegend, kBlack, kWhite, kRed, kGray, kBlue, TMath, TGraphAsymmErrors, TLatex, gStyle

sys.path.append(os.path.abspath(os.path.curdir))
from Core import NTupleTools, DatasetManager, listifyInputFiles

gROOT.Reset()
gROOT.LoadMacro(os.path.abspath(os.path.curdir)+"/Plotter/AtlasStyle.C")
# gROOT.LoadMacro("$HOME/RootUtils/AtlasStyle.C")
from ROOT import SetAtlasStyle
SetAtlasStyle()

class Inputs:

    def __init__(self):

        self.alltrees = {}
        self.sampleids = {}
        self.nomtree = 'physics'
        self.friendtrees = []
        self.friendfile_extension = None
        self.systrees = []
        self.sysweights = []
        self.readGFW2 = False

    def setFriendTree(self, friendtrees=[], friendfile_extension=None):

        self.friendtrees = friendtrees
        self.friendfile_extension = friendfile_extension


    def registerTree(self, filegroup, nomtree = 'physics', systrees=[], ismc=True, isembedding=False, isdata=False, sample={}, resetTreeWeight=False):

	# This function add a tree to TChain contained in self.alltrees = {}.
	# In this dictionary each group and subgroup is separated and has its proper TChain.
	# This function create the TChain for the group specified if not exist and if already exist add the tree founded to this chain.
	# I'm not very sure of the role of the sampleid variable. Can do this for many samples (all those in filegroup) and for many trees (for the nominal and for all those in systrees)

	sampleid = sample["ID"]
	group    = sample["group"]
	subgroup = sample["subgroup"]

        self.nomtree = nomtree
        self.systrees = systrees # These are the names of the trees with systematics
        syslist = []
        for t in systrees:
            syslist.append('SystematicsUP/'+t)
            syslist.append('SystematicsDOWN/'+t)
        treelist = [nomtree] + syslist

        prefix = ''
        if ismc:
            prefix = '$ISMC$'
        elif isembedding:
            prefix = '$ISEMBED$'
        elif isdata:
            prefix = '$ISDATA$'

        filelist = glob.glob(filegroup)
        if not filelist:
            print "WARNING: file", filegroup, "cannot be found during tree registration"
            return False

        # If not already done, set the Xsec weight to each (MC) TTree

        if not self.readGFW2:
            for filepath in filelist:
                f = TFile.Open(filepath,"UPDATE")
                for treename in treelist:
                    t = f.Get(treename)
                    if not t:
                        print ("WARNING: tree {0} in file {1} cannot be found during tree registration".format(treename,filepath))
                        return False
                    if ismc and ( t.GetWeight() == 1.0 or resetTreeWeight ):
                        print("Weighting tree w/ Xsec weight...")
                        weight = float(sample['xsection']) * float(sample['efficiency']) * float(sample['kfactor']) * 1e3 # To get the weight in fb (assuming the Xsec is in pb)
                        h = f.Get("TotalEventsW")
                        if not h:
                            print ("WARNING: histogram named TotalEventsW in file {1} couldn't be found!".format(filepath))
                            return False
                        weight /= h.GetBinContent(2)
                        t.SetWeight(weight)
                        t.Write(t.GetName(),t.kOverwrite)
                f.Close()

        # Add the TTrees into a TChain

        for filepath in filelist:
	    for treename in treelist:
                if not treename in self.alltrees:
                    self.alltrees[treename] = {}
                processes = self.alltrees[treename]
                if not group in processes:
                    processes[group] = {}
                if not subgroup in processes[group]:
                    processes[group][subgroup] = TChain(treename)
                    processes[group][subgroup].SetTitle(prefix+treename+group+subgroup)
                processes[group][subgroup].Add(filepath)
                # Add friend trees to *this* tree (if any)
                if self.friendtrees:
                    for friendname in self.friendtrees:
                        friendfilepath = filepath + self.friendfile_extension
                        processes[group][subgroup].AddFriend(friendname,friendfilepath)
                if sampleid:
                    processes[group][subgroup].SetTitle(processes[group][subgroup].GetTitle()+'_'+sampleid)
                    if not treename in self.sampleids:
                        self.sampleids[treename] = {}
                    self.sampleids[treename][sampleid] = (group, subgroup) # Seems that using sampleid you can then recover the info about the group and the subgroup...

        return True

    # Load a tree from the list of all trees

    def getTree(self, treename='physics', group='', subgroup='', sampleid=None):

        if self.readGFW2:
            treename = 'nominal'

	if sampleid:
            group, subgroup = self.sampleids[self.nomtree][sampleid]

        if treename.startswith('SystematicsUP/') or treename.startswith('SystematicsDOWN/'):
            if self.getTree(self.nomtree, group, subgroup).GetTitle().startswith('$ISDATA$'): # In case of data the tree to be considered is the nominal
                treename = self.nomtree

        try:
            tree = self.alltrees[treename][group][subgroup]
        except:
            tree = None
            print("WARNING: Could not reach tree {0} in group {1}, subgroup {2}".format(treename, group, subgroup))

	#print("\nTree: {0} - Xsec weight = {1}".format(tree.GetName(),tree.GetWeight()))
        return tree

    def getTrees(self, treename='physics', grouplist=[]):

        # Group list is a list of tuple with two elements (strings),
        # e.g. [ ('groupname1', 'subgroupname1'),('groupname2', '*'),] accepts also wildcards.
        # The entire list of trees is returned, one tree for each tuple.

        if self.readGFW2:
            treename = 'nominal'

	newGroupList = []
        for group, subgroup in grouplist:
            if group == '*' and subgroup == '*': # If using wildcard, the functions getGroupList and/or getSubGroupList are called to parse the '*' symbol
                for g in self.getGroupList():
                    for sg in self.getSubGroupList(g):
                        newGroupList.append( (g, sg) )
            elif subgroup == '*':
                for sg in self.getSubGroupList(group):
                    newGroupList.append( (group, sg) )
            elif group == '*':
                for g in self.getGroupList():
                    for sg in self.getSubGroupList(g):
                        if sg == subgroup:
                            newGroupList.append( (g, sg) )
            else:
                newGroupList.append( (group, subgroup) )

        treelist = []
        for g, sg in newGroupList:
            t = self.getTree(treename, g, sg)
            if t and not t in treelist:
                treelist.append(t)
        return treelist

    def getTreenameList(self):
        return self.alltrees.keys()

    def getGroupList(self):
        return self.alltrees[self.nomtree].keys()

    def getSubGroupList(self, group=''):
        return self.alltrees[self.nomtree][group].keys()

    def getSysIndexes(self, sampleID, branchID):

        if sampleID:
            if sampleID == "fakes_mm":
                sampletree = self.getTree(group='Data',subgroup='fakes_mm')
            else:
                sampletree = self.getTree(sampleid=sampleID)
        else:
            sampletree = self.getTree(group='Data',subgroup='physics_Main')

        branches_array = sampletree.GetListOfBranches()

        idxlist = []
        for b in range(0,branches_array.GetEntries()):
            branchname = branches_array.At(b).GetName()
            if not branchID in branchname or "Grouped" in branchname: continue
            tokens = branchname.split("_")
            # The bin index is supposed to be the second to last token in the branch name
            idx = tokens[-2]
            try:
                int(idx)
            except ValueError:
                print("\nERROR! String token: {0} cannot be converted to an integer\n".format(idx))
                raise
            if not idx in idxlist:
                idxlist.append(idx)
        return idxlist

class Variable:

    def __init__(self, **kw):
        self.latexname          = kw.get('latexname',   "Variable_{Plot}^{Name} [Units]")
        self.latexnameX         = kw.get('latexnameX',  self.latexname)
        self.latexnameY         = kw.get('latexnameY',  self.latexname)
        self.plainname          = kw.get('plainname',   self.latexname)
        self.shortname          = kw.get('shortname',   self.plainname)
        self.ntuplename         = kw.get('ntuplename',  "evtsel_name_of_variable")
        self.bins               = kw.get('bins',        40)
        self.binsX              = kw.get('binsX',	self.bins)
        self.binsY              = kw.get('binsY',	self.bins)
        self.minval             = kw.get('minval',      0.)
        self.maxval             = kw.get('maxval',      400.)
        self.minvalX            = kw.get('minvalX',	self.minval)
        self.maxvalX            = kw.get('maxvalX',	self.maxval)
        self.minvalY            = kw.get('minvalY',	self.minval)
        self.maxvalY            = kw.get('maxvalY',	self.maxval)
        self.typeval            = kw.get('typeval',     TH1D)
        self.manualbins         = kw.get('manualbins',  None)
        self.manualbinsX        = kw.get('manualbinsX', None)
        self.manualbinsY        = kw.get('manualbinsY', None)
        self.binlabelsX         = kw.get('binlabelsX',  {})
        self.logaxis            = kw.get('logaxis',     False)
        self.logaxisX           = kw.get('logaxisX',    False)
        self.basecut            = kw.get('basecut',     None)
        self.weight             = kw.get('weight',      None)
        self.sysvar             = kw.get('sysvar',      False)
        self.drawOpt2D          = kw.get('drawOpt2D',   "COLZ1") # The "1" does not plot empty bins even if there are bins w/ negative weights (NB: ROOT by default forces empty bins to be painted if there are negative weights!). See https://root.cern.ch/doc/master/classTHistPainter.html

        if '[' in self.latexname and ']' in self.latexname:
            self.unit = self.latexname[self.latexname.index('[')+1:self.latexname.index(']')]
        else:
            self.unit = None
        self.binarray = None
        if self.manualbins:
            self.binarray = array.array('d', self.manualbins)
        self.binarrayX = None
        if self.manualbinsX:
            self.binarrayX = array.array('d', self.manualbinsX)
        self.binarrayY = None
        if self.manualbinsY:
            self.binarrayY = array.array('d', self.manualbinsY)

    def ytitle(self, manualbins=None):
        if (manualbins and type(manualbins) is list) or self.typeval is TH1I:
            return 'Events / bin'
	else:
            if manualbins:
                bins, minval, maxval = manualbins
            else:
                bins, minval, maxval = self.bins, self.minval, self.maxval
            binwidth = (maxval - minval) / bins
            if self.unit:
                return 'Events / %.2g %s' % (binwidth, self.unit)
            else:
                return 'Events / bin'

    def makeHist(self, name=None, title=None, category=None):
        if name is None:
            name = self.plainname
        if title is None:
            title = self.latexname

        if any( self.typeval is t for t in [TH2D,TH2F,TH2I] ):

            if not self.manualbinsX and not self.manualbinsY:
                h = self.typeval(name, title, self.binsX, self.minvalX, self.maxvalX, self.binsY, self.minvalY, self.maxvalY)
            elif self.manualbinsX and not self.manualbinsY:
                h = self.typeval(name, title, len(self.manualbinsX)-1, self.binarrayX, self.binsY, self.minvalY, self.maxvalY)
            elif not self.manualbinsX and self.manualbinsY:
                h = self.typeval(name, title, self.binsX, self.minvalX, self.maxvalX, len(self.manualbinsY)-1, self.binarrayY)
            else:
                h = self.typeval(name, title, len(self.manualbinsX)-1, self.binarrayX, len(self.manualbinsY)-1, self.binarrayY)

            h.GetXaxis().SetTitle(self.latexnameX)
            h.GetYaxis().SetTitle(self.latexnameY)

        else :
	    if category and category.overridebins and self.shortname in category.overridebins:
            	manualbins = category.overridebins[self.shortname]
            	if type(manualbins) is list:
            	    binarray = array.array('d', manualbins)
            	else:
            	    binarray = None
            else:
            	manualbins = self.manualbins
            	binarray = self.binarray

            if manualbins:
            	if type(manualbins) is list:
            	    h = self.typeval(name, title, len(self.manualbins)-1, self.binarray)
            	else:
            	    h = self.typeval(name, title, manualbins[0], manualbins[1], manualbins[2])
            else:
            	h = self.typeval(name, title, self.bins, self.minval, self.maxval)
            h.SetXTitle(self.latexname)
            h.SetYTitle(self.ytitle(manualbins=manualbins))
            if self.binlabelsX:
                for idx,label in self.binlabelsX.iteritems():
                    h.GetXaxis().SetBinLabel(idx,label)
                    h.GetYaxis().SetLabelSize(h.GetXaxis().GetLabelSize())
        h.Sumw2()
        return h

class Cut:

    # A Cut is defined by a name and a set of rules defined in the cut string,
    # but can be also the composition of several cuts specified in the cut list parameter.

    def __init__(self, cutname, cutstr, cutlist=None):
        self.cutname    = cutname
        self.cutstr     = cutstr
        if cutlist is None:
            cutlist = [self]
        self.cutlist     = cutlist
        self.cutnamelist = [c.cutname for c in cutlist]

    # Removes a cut, provided it's found in the list of cuts

    def removeCut(self, cut):
        newlistnames = []
        newliststr = []
        newlist = []
        for c in self.cutlist:
            if c.cutname == cut.cutname: continue
            newlistnames.append(c.cutname)
            newliststr.append(c.cutstr)
            newlist.append(c)
        newname = ' AND '.join(newlistnames)
        newstr  = ' && '.join(newliststr)
        return Cut(newname, newstr, newlist)

    # Substitute a cut w/ another cut

    def swapCut(self, cutremove, cutadd):
        cut = self.removeCut(cutremove)
        return cut & cutadd

    # Bitwise AND of cuts. Should be called in this way: cut1 & cut2
    # This works also for concatenating lists of cuts.

    def __and__(self, othercut):
        if not othercut:
            newname = self.cutname
            newstr = self.cutstr
            newlist = self.cutlist
        else:
            newlistnames = []
            newliststr = []
            newlist = []
            for c in sorted(self.cutlist + othercut.cutlist, key=lambda x: x.cutname):
                if c.cutname in newlistnames: continue # Don't use the same cut more than once!
                newlistnames.append(c.cutname)
                if not any( s == c.cutstr for s in ['( 1 )','1'] ): # Add only non-dummy cuts to the TTreeFormula
                    newliststr.append(c.cutstr)
                newlist.append(c)
            newname = ' AND '.join(newlistnames)
            newstr  = ' && '.join(newliststr)
        return Cut(newname, newstr, newlist)

    # Bitwise OR of two cuts. Should be called in this way: cut1 | cut2

    def __or__(self, othercut):
        if not othercut:
            newname = self.cutname
            newstr = self.cutstr
            newlist = self.cutlist
        else:
            cut1, cut2 = sorted([self, othercut], key=lambda x: x.cutname)
            newname = '( ( %s ) OR ( %s ) )' % (cut1.cutname, cut2.cutname)
            newstr  = '( ( %s ) || ( %s ) )' % (cut1.cutstr, cut2.cutstr)
            newlist = None
        return Cut(newname, newstr, newlist)

    # Negation of a cut. Should be called in this way: -cut

    def __neg__(self):
        newname = 'NOT ( %s )' % (self.cutname)
        newstr = '!( %s )' % (self.cutstr)
        return Cut(newname, newstr)

class Systematics:
    def __init__(self, name, treename=None, eventweight=None, process=[], categorytokens=None):
        self.name           = name
        self.treename       = treename
        self.eventweight    = eventweight
        self.process        = process
        self.categorytokens = categorytokens

class Category:
    def __init__(self, name, overridebins = {}):
        self.name   = name
        self.tokens = name.split(' ')
        self.overridebins = overridebins

class VariableDB:

    # This class contains the list of all the variables, cuts, systematics and categories to be used to analyse our data.
    # Info are registered using the register methods.
    # Info stored in this DB can be retrieved at any time using the getter methods...

    def __init__(self):
        self.vardb = {}
        self.cutdb = {}
        self.systdb = {}
        self.categorydb = {}
        self.systlist = []
        self.varlist = []
        self.cutlist = []
        self.categorylist = []

    def registerVar(self, var):
        self.vardb[var.shortname] = var
        self.varlist.append(var)

    def registerCut(self, cut):
        self.cutdb[cut.cutname] = cut
        self.cutlist.append(cut)

    def registerSystematics(self, syst):
        self.systdb[syst.name] = syst
        self.systlist.append(syst)

    def printSystematics(self):
        print("\nRegistered systematics:\n")
        for syst in self.systlist:
            print("{0}".format(syst.name))
        print("")

    def registerCategory(self, category):
        self.categorydb[category.name] = category
        self.categorylist.append(category)

    def getVar(self, name):
        return self.vardb[name]

    def getCut(self, cutname):
        return self.cutdb[cutname]

    def getCuts(self, cutlist):
        cut = None
        for cutname in cutlist:
            if not cutname: continue
            if not self.cutdb[cutname]: continue
            if not cut:
                cut = self.cutdb[cutname]
            else:
                cut = cut & self.cutdb[cutname]
        return cut

    def getSyst(self, name):
        return self.systdb[name]

    def getCategory(self, name):
        return self.categorydb[name]

class SubProcess:

    histcache = {}
    numcache = {}

    def __init__(self, tree, basecut=None, baseweight=1.0, eventweight=None, debug=False):
        treename = tree.GetName()
        self.name = 'SubProcess:'+tree.GetTitle()+'_EVTWGT_'+str(eventweight)
        self.tree = tree
        if type(basecut) is str:
            basecut = Cut(basecut, basecut)
        self.basecut = basecut
        self.baseweight = baseweight
        self.eventweight = eventweight
        self.debug = debug

    def __str__(self):
        if self.basecut:
            return self.name + '_BASEWGT_' + str(self.baseweight) + '_BASECUT_' + self.basecut.cutname
        else:
            return self.name + '_BASEWGT_' + str(self.baseweight) + '_NOBASECUT'

    def __repr__(self):
        return self.__str__()

    def __add__(self, other):
        return OperatorProcess(self, '+', other)

    def __sub__(self, other):
        return OperatorProcess(self, '-', other)

    def __mul__(self, other):
        return OperatorProcess(self, '*', other)

    def __div__(self, other):
        return OperatorProcess(self, '/', other)

    def rootNameFriendly(self, name):
        charlist = '[]()/\\ '
        for c in charlist:
            name = name.replace(c, '_')
        return name

    def subprocess(self, tree=None, cut=None, weight=1.0, eventweight=None, clearbasecut=False, debug=False):
        if clearbasecut:
            self.basecut = None
        if tree is None:
            tree = self.tree
        if self.eventweight and eventweight:
            eventweight = '%s * %s' % (self.eventweight, eventweight)
        elif self.eventweight:
            eventweight = self.eventweight
        if cut and self.basecut:
            cut = self.basecut & cut
        elif self.basecut and not cut:
            cut = self.basecut

        sp = SubProcess(tree=tree, basecut=cut, baseweight=self.baseweight*weight, eventweight=eventweight, debug=debug)
        return sp

    def numberstats(self, cut=None, weight=None, eventweight=None, category=None):

        if weight is None:
            weight = 1.0

        if type(cut) is str:
            cut = Cut(cut, cut)
        if cut and self.basecut:
            cut = self.basecut & cut
        elif self.basecut and not cut:
            cut = self.basecut

        if cut:
            cutstr = cut.cutstr
            cachename = self.name+'_CUT_'+cut.cutname
        else:
            cutstr = '1'
            cachename = self.name+'_NOCUT'
        cachename = self.rootNameFriendly(cachename)

        if cachename in self.numcache:
            num, stat = self.numcache[cachename]
            return num * self.baseweight * weight, stat * self.baseweight * weight

        h = TH1D('NUM'+cachename, 'NUM'+cachename, 1, 0., 2.)
        h.Sumw2()

        if self.eventweight and eventweight:
            self.tree.Project('NUM'+cachename, '1.0', '%s * %s * (%s)' % (self.eventweight, eventweight, cutstr))
            if self.debug:
                print("\nApplying TTreeFormula string:\n\n{0}->Project(\"1\",\"{1} * {2} * {3} * ( {4} )\")\n".format(self.tree.GetName(), self.baseweight, self.eventweight, eventweight, cutstr))
        elif self.eventweight :
            self.tree.Project('NUM'+cachename, '1.0', '%s * (%s)' % (self.eventweight, cutstr))
            if self.debug:
                print("\nApplying TTreeFormula string:\n\n{0}->Project(\"1\",\"{1} * {2} * ( {3} )\")\n".format(self.tree.GetName(), self.baseweight, self.eventweight, cutstr))
        else:
            self.tree.Project('NUM'+cachename, '1.0', '%s' % (cutstr))
            if self.debug:
                print("\nApplying TTreeFormula string:\n\n{0}->Project(\"1\",\"{1} * ( {2} )\")\n".format(self.tree.GetName(), self.baseweight, cutstr))
        self.numcache[cachename] = h.GetBinContent(1), h.GetBinError(1)
        del h

        num, stat = self.numcache[cachename]
        return num * self.baseweight * weight, stat * self.baseweight * weight

    def number(self, cut = None, weight = None, eventweight = None, category = None):
        return self.numberstats(cut, weight, eventweight, category)[0]

    def stats(self, cut = None, weight = None, eventweight = None, category = None):
        return self.numberstats(cut, weight, eventweight, category)[1]

    def hist(self, var, cut = None, weight = None, category = None):

        if weight is None:
            weight = 1.0

        if type(cut) is str:
            cut = Cut(cut, cut)
        if cut and self.basecut:
            cut = self.basecut & cut
        elif self.basecut and not cut:
            cut = self.basecut
        if cut and var.basecut:
            cut = cut & var.basecut
        elif var.basecut and not cut:
            cut = var.basecut
        if not "$ISDATA$" in self.name:
	    if type(var.weight) is str:
                if not var.weight in self.eventweight:
            	    self = self.subprocess(eventweight=var.weight)
            elif type(var.weight) is float:
            	weight *= var.weight

        if category and category.overridebins and var.shortname in category.overridebins:
            binname = str(category.overridebins[var.shortname])
        else:
            binname = 'Default'

        if cut:
            cutstr = cut.cutstr
            cachename = self.name+'_VAR_'+var.shortname+'_CUT_'+cut.cutname+'_BIN_'+binname
        else:
            cutstr = '1'
            cachename = self.name+'_VAR_'+var.shortname+'_NOCUT'+'_BIN_'+binname
        cachename = self.rootNameFriendly(cachename)

        if cachename in self.histcache:
            h = self.histcache[cachename].Clone()
            h.SetName(self.histcache[cachename].GetName()+str(weight))
            h.SetTitle(self.histcache[cachename].GetTitle()+str(weight))
            h.__imul__(self.baseweight * weight)
            return h

        self.histcache[cachename] = var.makeHist('HIST'+cachename, 'HIST'+cachename, category)

        if self.eventweight:
            self.tree.Project('HIST'+cachename, var.ntuplename, '%s * (%s)' % (self.eventweight, cutstr))
            if self.debug:
                print("\nApplying TTreeFormula string:\n\n{0}->Project(\"{1}\",\"{2} * {3} * ( {4} )\")\n".format(self.tree.GetName(), var.ntuplename, self.baseweight, self.eventweight, cutstr))
	else:
            self.tree.Project('HIST'+cachename, var.ntuplename, '%s' % (cutstr))
            if self.debug:
                print("\nApplying TTreeFormula string:\n\n{0}->Project(\"{1}\",\"{2} * ( {3} )\")\n".format(self.tree.GetName(), var.ntuplename, self.baseweight, cutstr))

        h = self.histcache[cachename].Clone()
        h.SetName(self.histcache[cachename].GetName()+str(weight))
        h.SetTitle(self.histcache[cachename].GetTitle()+str(weight))
        h.__imul__(self.baseweight * weight)

        return h

class OperatorProcess(SubProcess):

    def __init__(self, left, operator, right):

	leftname = left.name
        if type(right) is int:
            right = float(right)
        if type(right) is float:
            rightname = str(right)
        else:
            rightname = right.name
        self.basecut = left.basecut
        self.baseweight = left.baseweight
        self.eventweight = left.eventweight
        self.name = '('+leftname+operator+rightname+')'
        self.left = left
        self.operator = operator
        self.right = right
        self.debug = left.debug

    def __str__(self):
        return 'OperatorProcess:(' + self.left.__str__() + ' ' + self.operator + ' ' + self.right.__str__() + ')'

    def subprocess(self, tree=None, cut=None, weight=1.0, eventweight=None, clearbasecut=False, debug=False):
    	if hasattr(self.left, 'subprocess'):
    	    left = self.left.subprocess(tree, cut, weight, eventweight, clearbasecut, debug)
    	else:
    	    left = self.left
    	if hasattr(self.right, 'subprocess'):
    	    right = self.right.subprocess(tree, cut, weight, eventweight, clearbasecut, debug)
    	else:
    	    right = self.right

    	sp = OperatorProcess(left, self.operator, right)
    	return sp

    def numberstats(self, cut = None, weight = None, eventweight = None, category = None):
        number = self.number(cut, weight, eventweight, category)
        stats = self.stats(cut, weight, eventweight, category)
        return number, stats

    def number(self, cut = None, weight = None, eventweight = None, category = None):
        leftnum = self.left.number(cut, weight, eventweight, category)
        if type(self.right) is float:
            rightnum = self.right
        else:
            rightnum = self.right.number(cut, weight, eventweight, category)
        if self.operator == '+':
            return leftnum + rightnum
        if self.operator == '-':
            return leftnum - rightnum
        if self.operator == '*':
            return leftnum * rightnum
        if self.operator == '/':
            return leftnum / rightnum
        print "ERROR: Invalid operator", self.operator

    def stats(self, cut = None, weight = None, eventweight = None, category = None):
        leftnum, leftstats = self.left.numberstats(cut, weight, eventweight, category)
        if type(self.right) is float:
            rightnum, rightstats = self.right, 0.
        else:
            rightnum, rightstats = self.right.numberstats(cut, weight, eventweight, category)
        if self.operator in ['+', '-']:
            return math.sqrt(leftstats**2. + rightstats**2.)
        if self.operator in ['*', '/']:
            errfrac = 0.
            if leftnum:
                errfrac += (leftstats/leftnum)**2.
            if rightnum:
                errfrac += (rightstats/rightnum)**2.
            errfrac = math.sqrt(errfrac)
            if self.operator == '*':
                return errfrac * leftnum*rightnum
            else:
                return errfrac * leftnum/rightnum
        print "ERROR: Invalid operator", self.operator

    def hist(self, var, cut = None, weight = None, category = None):

        if weight is None:
            weight = 1.0

        if type(cut) is str:
            cut = Cut(cut, cut)

        lefthist = self.left.hist(var, cut, weight, category)
        if type(self.right) is float:
            righthist = self.right
        else:
            righthist = self.right.hist(var, cut, weight, category)

        name = self.name + var.shortname
        if cut: name += cut.cutname
        if weight: name += str(weight)
        if issubclass(type(lefthist), TH1):
            h = lefthist.Clone()
            h.SetName(name)
            h.SetTitle(name)
        elif issubclass(type(righthist), TH1):
            h = righthist.Clone()
            h.SetName(name)
            h.SetTitle(name)
        else:
            h = None
        if self.operator == '+':
            if type(lefthist) is float and type(righthist) is float:
                h = lefthist + righthist
            else:
                h.Add(righthist)
        elif self.operator == '-':
            if type(lefthist) is float and type(righthist) is float:
                h = lefthist - righthist
            else:
                h.Add(righthist, -1.0)
        elif self.operator == '*':
            if type(lefthist) is float and type(righthist) is float:
                h = lefthist * righthist
            elif type(righthist) is float:
                h.__imul__(righthist)
            elif type(lefthist) is float:
                h.__imul__(lefthist)
            else:
                h.Multiply(righthist)
        elif self.operator == '/':
            if type(lefthist) is float and type(righthist) is float:
                h = lefthist / righthist
            elif type(righthist) is float:
                h.__imul__(1.0/righthist)
            elif type(lefthist) is float:
                h.__imul__(1.0/lefthist)
            else:
                h.Divide(righthist)
        else:
            print "ERROR: Invalid operator", self.operator
            return
        return h

class ScaleFactors(SubProcess):
    def number(self, cut = None, weight = None, eventweight = None, category = None):
        return 1.0

    def hist(self, var, cut = None, weight = None, eventweight = None, category = None):
        return self.number(cut, weight, eventweight, category)

class Process:

    name = 'Default'

    def __init__(self, inputs, vardb, parent):
        self.inputs = inputs
        self.vardb  = vardb
        self.parent = parent

    def __call__(self, treename='physics', category=None, options={}):
        pass

    def subprocess(self, trees=[], basecut=None, baseweight=1.0, eventweight=None):

        if self.parent.eventweight and eventweight:
            eventweight = '%s * %s' % (self.parent.eventweight, eventweight)
        elif self.parent.eventweight:
            eventweight = self.parent.eventweight
        elif not self.parent.eventweight and not eventweight:
            eventweight = 1.0
        sp = None
        for tree in trees:
            weight = baseweight
	    if tree.GetTitle().startswith('$ISDATA$'): # Make sure no eventweight is applied if looking at data...
	        eventweight = None
            if tree.GetTitle().startswith('$ISMC$') or tree.GetTitle().startswith('$ISEMBED$'):
                weight *= self.parent.luminosity
		if self.parent.rescaleXsecAndLumi:
		    weight /= tree.GetWeight()
		    weight /= self.parent.luminosity

            s = SubProcess(tree=tree, basecut=basecut, baseweight=weight, eventweight=eventweight)
            if sp: sp = sp + s
            else: sp = s
        return sp

class Background:

    backgrounds = []
    signals     = []
    observed    = []
    luminosity  = 1.0
    eventweight = '1.0'
    rescaleXsecAndLumi = False
    style = {}
    readGFW2 = False

    def __init__(self, inputs, vardb):
        self.inputs      = inputs
        self.vardb       = vardb
        self.procmap     = {}
        self.debugprocs  = []
        self.colourcache = self.colours()
        self.var         = Variable(shortname = "Integral", latexname = "", ntuplename = "0.5", bins = 1, minval = 0.0, maxval = 1.0)

        for key in dir(self):
            attr = getattr(self, key)
            if inspect.isclass(attr) and issubclass(attr, Process):
                self.procmap[attr.__name__] = attr(inputs, vardb, self)

    def labels(self, legs, showratio):
        mid = int(len(legs)/2)
        high = math.ceil(len(legs)/2)
        lower = 0.92 - 0.04*high
        leg1 = TLegend(0.60,lower,0.75,0.92)
        leg2 = TLegend(0.75,lower,0.90,0.92)
        for leg in [leg1, leg2]:
            leg.SetFillColor(0)
            leg.SetFillStyle(0)
            leg.SetLineColor(10)
            leg.SetShadowColor(kWhite)
            leg.SetTextSize(0.03)
            leg.SetBorderSize(0)

        for l in legs[:mid]:
            leg1.AddEntry(l[0], l[1], l[2])
        for l in legs[mid:]:
            leg2.AddEntry(l[0], l[1], l[2])
        leg1.Draw()
        leg2.Draw()

        lumtext = drawText(text="  \\int L dt = %.2g fb^{-1}"%(self.luminosity), x=.2, y=.87, size=0.03)
        cmetext = drawText(text="         #sqrt{s} = 7TeV", x=.2, y=.82, size=0.03)
        atlastext = drawText(text="#bf{#it{ATLAS}} Preliminary", x=.2, y=.77, size=0.03)
        return lower, locals()

    def colours(self):
        cache = []
        cache.append(TColor(1040, 82/255., 195/255., 229/255.))     # blue
        cache.append(TColor(1041, 66/255., 156/255., 183/255.))     # dark blue
        cache.append(TColor(1030, 229/255., 229/255., 121/255.))    # yellow
        cache.append(TColor(1020, 220/255., 87/255., 60/255.))      # red
        cache.append(TColor(1010, 103/255., 73/255., 130/255.))     # purple
        cache.append(TColor(1000, 108/255., 178/255., 81/255.))     # green
        return cache

    def getProcess(self, name, category=None, systematics=None, systematicsdirection=None, options={}):
        treename = 'physics'
        eventweight = None
        if systematics:
            matchtoken = True
            if systematics.categorytokens:
                matchtoken = False
                for c in systematics.categorytokens:
                    if c in category.tokens:
                        matchtoken = True
                        break
            if matchtoken:
                if systematics.treename:
                    if systematicsdirection == 'UP':
                        treename = 'SystematicsUP/' + systematics.treename
                    elif systematicsdirection == 'DOWN':
                        treename = 'SystematicsDOWN/' + systematics.treename
                if systematics.eventweight and name in systematics.process:
		    if systematicsdirection == 'UP':
                        if type(systematics.eventweight) is str:
                            eventweight = ( systematics.eventweight + 'up' ) if not self.readGFW2 else ( systematics.eventweight + 'Up' )
                        else:
                            eventweight = 1.0 + systematics.eventweight
                    elif systematicsdirection == 'DOWN':
                        if type(systematics.eventweight) is str:
                            eventweight = ( systematics.eventweight + 'dn' ) if not self.readGFW2 else ( systematics.eventweight + 'Dn' )
                        else:
                            eventweight = 1.0 - systematics.eventweight

        if systematics and (name in systematics.process or not systematics.process):
            options['systematics'] = systematics
            options['systematicsdirection'] = systematicsdirection
        process = self.procmap[name](treename=treename, category=category, options=options)
        if eventweight:
            process = process.subprocess(eventweight=eventweight)

        return process

    def parseArguments(self, cut, category, systematics, overridebackground):
        if cut is None:
            cut =Cut('NoCut', '1')
        elif type(cut) is str:
            cut = self.vardb.getCut(cut)
        elif type(cut) is list:
            cut = self.vardb.getCuts(cut)
        if category is None:
            category = self.vardb.getCategory('VBF')
        elif type(category) is str:
            category = self.vardb.getCategory(category)
        if type(systematics) is str:
            systematics = self.vardb.getSyst(systematics)
        if overridebackground is None:
            overridebackground = self.backgrounds

        return cut, category, systematics, overridebackground

    def events(self, cut = None, eventweight=None, category = None, hmass = ['125'], systematics = None, systematicsdirection = None, overridebackground = None, show = True):
        cut, category, systematics, overridebackground = self.parseArguments(cut, category, systematics, overridebackground)

        if show:
	    print "\nNumber of events with", cut.cutname, "in", category.name, "..."

        hbkg = {}
        bkg = (0., 0.)
        for name in overridebackground:
            process = self.getProcess(name, category=category, systematics=systematics, systematicsdirection=systematicsdirection)
            process = process.subprocess(debug=any( proc in self.debugprocs for proc in [name,"ALL"]))
            events, stats = process.numberstats(cut=cut, eventweight=eventweight)
            hbkg[name] = (events, stats)
            bkg = bkg[0] + events, math.sqrt(bkg[1]**2. + stats**2.)
            if show: print "%40s : %.2f +- %.2f (stat)" % (name, round(events, 2), round(stats, 2))
        if show: print "%40s : %.2f +- %.2f (stat)\n" % ('TOTAL BACKGROUND', round(bkg[0], 2), round(bkg[1], 2))
        hbkg['TOTAL'] = bkg

        hobs = {}
        obs = (0., 0.)
        for name in self.observed:
            process = self.getProcess(name, category=category, systematics=systematics, systematicsdirection=systematicsdirection)
            process = process.subprocess(debug=any( proc in self.debugprocs for proc in [name,"ALL"]))
            events, stats = process.numberstats(cut=cut) # No extra event weight for data!
            hobs[name] = (events, stats)
            obs = obs[0] + events, math.sqrt(obs[1]**2. + stats**2.)
            if show: print "%40s : %.2f +- %.2f (stat)" % (name, round(events, 2), round(stats, 2))
        hobs['TOTAL'] = obs
        if show: print "%40s : %.1f +- %.1f (stat)\n" % ('TOTAL OBSERVED', round(obs[0], 1), round(obs[1], 1))

        hsig = {}
        for h in hmass:
            sig = (0., 0.)
            for name in self.signals:
                options = {'hmass': h}
                process = self.getProcess(name, category=category, systematics=systematics, systematicsdirection=systematicsdirection, options=options)
                process = process.subprocess(debug=any( proc in self.debugprocs for proc in [name,"ALL"]))
                events, stats = process.numberstats(cut=cut, eventweight=eventweight)
                sig = sig[0] + events, math.sqrt(sig[2]**2. + stats**2.)
                if show: print "%36s %3s : %.2f +- %.2f (stat)" % (name, h, round(events, 2), round(stats, 2))
            hsig[h] = sig
            if show: print "%40s : %.2f +- %.2f (stat)\n" % ('TOTAL SIGNAL', round(sig[0], 2), round(sig[1], 2))

        return hbkg, hobs, hsig

    def sumhist(self, var, processes = [], cut = None, eventweight = None, category = None, systematics = None, systematicsdirection = None, scale = 1.0, overflowbins = False, options={}):

	tSum = None
        histlist = []

        for name in processes:

            process = self.getProcess(name, category=category, systematics=systematics, systematicsdirection=systematicsdirection, options=options) * scale

            if eventweight:
                weight = eventweight
                if "$ISDATA$" in process.name: # Make sure no eventweight is applied if looking at data...
                    weight = None
                process = process.subprocess(eventweight=weight)

            # Check whether debug flag should be activated for *this* process

            process = process.subprocess(debug=any( proc in self.debugprocs for proc in [name,"ALL"]))

            h = process.hist(var, cut=cut, category=category)

            if overflowbins and not isinstance(h,TH2):
                lastbin = h.GetSize()-2
                over, overerror = h.GetBinContent(lastbin+1), h.GetBinError(lastbin+1)
                under, undererror = h.GetBinContent(0), h.GetBinError(0)
                if '_ext' in var.shortname: # Convention for "extended" histograms: plot overflow, but not underflow
                    under, undererror = 0., 0.
                h.SetBinContent(lastbin, h.GetBinContent(lastbin) + over)
                h.SetBinError(lastbin, math.sqrt( h.GetBinError(lastbin)**2. + overerror**2.))
                h.SetBinContent(1, h.GetBinContent(1) + under)
                h.SetBinError(1, math.sqrt( h.GetBinError(1)**2. + undererror**2. ))
            if not tSum:
                tSum = h.Clone()
                if cut: cutname = cut.cutname
                else: cutname = 'None'
                if category: catname = category.name
                else: catname = 'None'
                if systematics: systname = systematics.name + systematicsdirection
                else: systname = 'None'
                tSum.SetName('SumHist:'+var.shortname+'_PROC_'+'+'.join(processes)+'_CUT_'+cutname+'_CAT_'+catname+'_EVTWGT_'+str(eventweight)+'_SYST_'+systname)
                tSum.SetTitle('SumHist:'+var.shortname+'_PROC_'+'+'.join(processes)+'_CUT_'+cutname+'_CAT_'+catname+'_EVTWGT_'+str(eventweight)+'_SYST_'+systname)
            else:
                tSum.Add(h)
            histlist.append( (h, self.procmap[name]) )

        # Sort histlist, from smallest yield to largest

        histlist.sort( key=lambda x: x[0].Integral() )

        return tSum, histlist

    def plot(self, var, cut = None, eventweight=None, category = None, signal = '125', signalfactor = 1., systematics = None, systematicsdirection = None, overridebackground = None, overflowbins = False, showratio = True, wait = False, save = ['.eps'], options = {}, normalise = False, log=False, logx=False, showyields=False, nolegs=False):

	if not wait:
            gROOT.SetBatch(True)
        cut, category, systematics, overridebackground = self.parseArguments(cut, category, systematics, overridebackground)
        if type(var) is str:
            var = self.vardb.getVar(var)

        self.var = var

        isvar2D = ( var.typeval in [TH2,TH2I,TH2F,TH2D] )

        if not isvar2D:
            c = TCanvas("c1","Temp",50,50,600,600)
	else:
	    c = TCanvas("c1","Temp",50,50,800,600)

        legs = []

        obs, obslist = self.sumhist(var, processes=self.observed, cut=cut, eventweight=eventweight, category=category, systematics=systematics, systematicsdirection=systematicsdirection, overflowbins=overflowbins)

	if obs:
            process = obslist[0][1]
            datagr = None
            if not ( "$ISDATA$" in obslist[0][0].GetName() ):
                if not isvar2D:
                    datagr = TH1D(obs) # Equivalent to: datagr = makeMCErrors(obs)
                    datagr.SetLineStyle(2)
                else:
                    datagr = TH2D(obs)
            else :
                if not isvar2D:
                    datagr = makePoissonErrors(obs)
                    datagr.SetMarkerSize(0.8) # (1.2)
                    datagr.SetLineColor(self.style.get('ObservedLineColour', 1))
                    datagr.SetLineWidth(self.style.get('ObservedLineWidth', 1))
                    datagr.SetMarkerStyle(self.style.get('ObservedMarkerStyle', 20))
                else:
                    datagr = TH2D(obs)
            legobs = process.latexname + " ({0:.1f})".format(integrate(obs,overflowbins)) if showyields else process.latexname
            legs.append([datagr, legobs, "p"])

        tSum, bkglist = self.sumhist(var, processes=overridebackground, cut=cut, eventweight=eventweight, category=category, systematics=systematics, systematicsdirection=systematicsdirection, overflowbins=overflowbins, options=options)

        if obs and bkglist and normalise:
            if not isvar2D:
	    	num_data = obs.GetEntries()
            	num_mc = 0.0
            	for b, bname in bkglist:
            	    for i in range(b.GetSize()):
            		num_mc += b.GetBinContent(i)
            	if num_mc:
            	    normratio = num_data / num_mc
            	    for b, bname in bkglist:
            		b *= normratio
            	    tSum *= normratio

        bkg = {}
	if bkglist:
	    stack = THStack('Stack'+tSum.GetName(), self.__class__.__name__+';'+tSum.GetXaxis().GetTitle()+';'+tSum.GetYaxis().GetTitle())
	    for h, process in reversed(bkglist):
	    	h.Draw()
	    	h.SetLineWidth(self.style.get('BackgroundLineWidth', 1))
	    	h.SetLineStyle(self.style.get('BackgroundLineStyle', 1))
	    	pname = process.__class__.__name__
	    	h.SetLineColor(self.style.get(pname+'LineColour', self.style.get('BackgroundLineColour', 1)))
	    	h.SetFillColor(self.style.get(pname+'FillColour', process.colour))
	    	h.SetFillStyle(self.style.get(pname+'FillStyle', 1001))
	    	stack.Add(h)
                legbkg = process.latexname + " ({0:.1f})".format(integrate(h,overflowbins)) if showyields else process.latexname
	    	legs.append([h, legbkg, 'f'])
	    	bkg[pname] = h
	else:
	    stack = None
	    print "No Stack"

        if bkg:
            gStyle.SetHatchesLineWidth(1)
            gStyle.SetHatchesSpacing(0.8)#(0.4)
            tSum.SetFillColor(self.style.get('SumErrorFillColour', kBlue))#(self.style.get('SumErrorFillColour', kGray+1))
            tSum.SetLineColor(self.style.get('SumErrorLineColour', kBlue))#(self.style.get('SumErrorLineColour', kGray+1))
            tSum.SetFillStyle(self.style.get('SumErrorFillStyle', 3356))
            tSum.SetMarkerSize(0)
            legs.insert(0, (tSum,"#bf{Stat. unc.}","F"))
        else:
            print "No background processes are plotted!"

	options['hmass'] = signal
        sig, siglist = self.sumhist(var, processes=self.signals, cut=cut, eventweight=eventweight, category=category, systematics=systematics, systematicsdirection=systematicsdirection, overflowbins=overflowbins, scale=signalfactor, options=options)

        if sig:
            process = siglist[0][1]
            sig.SetFillColor(self.style.get('SignalFillColour', 10))
            sig.SetFillStyle(self.style.get('SignalFillStyle', 1001))
            sig.SetLineWidth(self.style.get('SignalLineWidth', 1))
            sig.SetLineColor(self.style.get('SignalLineColour', 2))
            sig.SetLineStyle(self.style.get('SignalLineStyle', 2))
            if not bkglist:
                stack = THStack('Stack'+sig.GetName(), self.__class__.__name__+';'+sig.GetXaxis().GetTitle()+';'+sig.GetYaxis().GetTitle())
            stack.Add(sig)
            h_name = process.latexname+signal
            if signalfactor != 1.:
               h_name += " [#times"+str(int(signalfactor))+']'
            legsig = h_name + " ({0:.1f})".format(integrate(sig,overflowbins)) if showyields else process.latexname
            legs.append([sig, legsig, 'f'])

        # Never show ratio if variable is 2D

        if isvar2D:
            showratio = False

        if showratio and obs and bkg:
            pad1 = TPad("pad1", "", 0, 0.25, 1, 1)
            pad2 = TPad("pad2", "", 0, 0, 1, 0.25)
            pad1.SetBottomMargin(0.02)
            pad2.SetBottomMargin(0.4)
            #pad2.SetTopMargin(0)
            pad2.SetGridy(1)
            if log or var.logaxis:
                pad1.SetLogy()
                stack.SetMinimum(0.1)
            if logx or var.logaxisX:
                pad1.SetLogx()
                pad2.SetLogx()
            pad1.Draw()
            pad2.Draw()
	if not showratio:
            if log or var.logaxis:
                gPad.SetLogy()
                stack.SetMinimum(0.1)
            if logx or var.logaxisX:
                gPad.SetLogx()

        if showratio and obs and bkg:
            if var.typeval in [TH1D,TH1F]:
                gStyle.SetHatchesLineWidth(1)
                gStyle.SetHatchesSpacing(0.4)
                ratiomc_props = {
                    "TitleY": "" if not ( "$ISDATA$" in obslist[0][0].GetName() ) else "Data/Exp.",
                    "TitleX": var.latexname,
                    "TitleSizeX":0.15,
                    "TitleSizeY":0.15,
                    "TitleOffsetX":0.90,
                    "TitleOffsetY":0.35,
                    "LabelSizeX":0.15,
                    "LabelSizeY": 0.12,
                    "NdivisionsY":505,
                    "FillColor":kBlue,#kGray+1,
                    "FillStyle":self.style.get('SumErrorFillStyle', 3356),
                    "LineColor":kBlue,#kGray+1,
                    "MarkerSize":0,
                }
                ratiomc = SelfDivide( tSum, "RatioMC", ratiomc_props )
                ratiodata = obs.Clone("RatioData")
            elif var.typeval is TH1I:
                ratiomc = TH1D("RatioMC", "RatioMC", tSum.GetNbinsX(), tSum.GetBinLowEdge(1), tSum.GetBinLowEdge(tSum.GetNbinsX()+1))
                ratiomc.GetXaxis().SetNdivisions(tSum.GetNbinsX())
                ratiomc.GetXaxis().CenterLabels(True)
                ratiodata = TH1D("RatioData", "RatioData", tSum.GetNbinsX(), tSum.GetBinLowEdge(1), tSum.GetBinLowEdge(tSum.GetNbinsX()+1))
                ratiodata.SetLineWidth(2)
                for i in range(1, tSum.GetNbinsX()+1):
                    ratiomc.SetBinContent(i, tSum.GetBinContent(i))
                    ratiodata.SetBinContent(i, obs.GetBinContent(i))

            ratiodata.SetMarkerSize(0.8) # (0.3)
	    if not ( "$ISDATA$" in obslist[0][0].GetName() ):
	         ratiodata.SetLineStyle(2)
            ratiodata.SetLineWidth(1)
            ratiodata.Divide(tSum)

            soverb = None
            if showratio and sig and bkg and not  isvar2D:
                soverb = sig.Clone("SoverB")
                soverb.SetLineStyle(1)
                if signalfactor != 1.0:
                    soverb.Scale(1.0/signalfactor)
                soverb.Add(tSum)
                soverb.Divide(tSum)

            valYmin =  99.
            valYmax = -99.
            for i in range(1, tSum.GetXaxis().GetNbins()+1):
                if tSum.GetBinContent(i):
                    ratiomc.SetBinError(i, tSum.GetBinError(i) / tSum.GetBinContent(i))
                    ratiodata.SetBinError(i, obs.GetBinError(i) / tSum.GetBinContent(i))
                else:
                    ratiomc.SetBinError(i, 0.)
                    ratiodata.SetBinError(i, 0.)
                ratiomc.SetBinContent(i, 1.)

                if (ratiodata.GetBinContent(i) - ratiodata.GetBinError(i)) != 0. and ratiodata.GetBinContent(i) < valYmin:
                    valYmin = ratiodata.GetBinContent(i) - ratiodata.GetBinError(i)
                    if valYmin < 0.:
                        valYmin = 0.
                if (ratiodata.GetBinContent(i) + ratiodata.GetBinError(i)) != 0. and ratiodata.GetBinContent(i) > valYmax:
                    valYmax = ratiodata.GetBinContent(i) + ratiodata.GetBinError(i)

                if (ratiomc.GetBinContent(i) - ratiomc.GetBinError(i)) != 0. and ratiomc.GetBinContent(i) < valYmin:
                    valYmin = ratiomc.GetBinContent(i) - ratiomc.GetBinError(i)
                    if valYmin < 0.:
                        valYmin = 0.
                if (ratiomc.GetBinContent(i) + ratiomc.GetBinError(i)) != 0. and ratiomc.GetBinContent(i) > valYmax:
                    valYmax = ratiomc.GetBinContent(i) + ratiomc.GetBinError(i)

            if abs(1-valYmin) > abs(valYmax-1):
                val = abs(1-valYmin)
            else:
                val = abs(valYmax-1)

            # ratiomc.SetMinimum((1-val)-0.1)
            # ratiomc.SetMaximum((1+val)+0.1)
            ratiomc.SetMinimum(0)
            ratiomc.SetMaximum(2)

            if type(showratio) is tuple:
                if showratio[0] == "MIN" and type(showratio[1]) is float:
                    ratiomc.SetMinimum(valYmin - 0.1)
                    ratiomc.SetMaximum(showratio[1])
                elif showratio[1] == "MAX" and type(showratio[0]) is float:
                    ratiomc.SetMinimum(showratio[0])
                    ratiomc.SetMaximum(valYmax + 0.1)
                else:
                    ratiomc.SetMinimum(showratio[0])
                    ratiomc.SetMaximum(showratio[1])

            pad2.cd()
            ratiomc.Draw("E2")
            if soverb:
                # gStyle.SetPaintTextFormat(".2f S/B")
                # soverb.SetMarkerSize(3.8)
                # soverb.SetMarkerColor(self.style.get('SignalMarkerColour', 2))
                # soverb.Draw("HIST SAME TEXT0")
                soverb.SetFillStyle(0)
                # soverb.Draw("HIST SAME")
                reflsoverb = TLine(soverb.GetBinLowEdge(1), 1.15, soverb.GetBinLowEdge(soverb.GetNbinsX()+1), 1.15)
                reflsoverb.SetLineStyle(2)
                reflsoverb.SetLineColor(kRed)
                # reflsoverb.Draw("SAME")
            refl = TLine(ratiomc.GetBinLowEdge(1), 1., ratiomc.GetBinLowEdge(ratiomc.GetNbinsX()+1), 1.)
            refl.SetLineStyle(2)
            refl.SetLineColor(kBlack)
            refl.Draw("SAME")
            ratiodata.Draw("PE1 SAME")
            pad1.cd()

        legs.reverse()
        if not nolegs:
            lower, labels = self.labels(legs, showratio and obs and bkg)

        # Trick to rescale:

	if stack:
	   if not isvar2D:
	      ymax_new = stack.GetMaximum()
	      if obs and obs.GetMaximum() > ymax_new:
	          ymax_new = obs.GetMaximum()
	      if stack and stack.GetMaximum() > ymax_new:
	          ymax_new = stack.GetMaximum()
	      if showratio and bkg and obs:
	          stack.SetMaximum(ymax_new*(2.-lower+0.075))
	      else:
	          stack.SetMaximum(ymax_new*(2.-lower+0.15))
	      if log or var.logaxis:
	          #stack.SetMaximum(stack.GetMaximum() * 10**(1.5))
	          stack.SetMaximum(stack.GetMaximum() * 3*10**(2))

	      stack.Draw('HIST')

	      if showratio and obs and bkg:
	          stack.GetHistogram().GetXaxis().SetLabelOffset(999)
	          stack.GetHistogram().GetXaxis().SetLabelSize(0)
	          stack.GetHistogram().GetYaxis().SetTitleSize(stack.GetHistogram().GetYaxis().GetTitleSize() * 1.2)
	          stack.GetHistogram().GetYaxis().SetTitleOffset(1.20)
	          #ratiomc.GetXaxis().SetNdivisions(8)
	          ratiomc.GetXaxis().SetTitleSize(ratiomc.GetXaxis().GetTitleSize() * 1.2)
	          ratiomc.GetYaxis().SetTitleSize(ratiomc.GetYaxis().GetTitleSize() * 1.2)
	      else:
	          stack.GetHistogram().GetXaxis().SetLabelSize(stack.GetHistogram().GetXaxis().GetLabelSize() * 0.75)
	          stack.GetHistogram().GetYaxis().SetLabelSize(stack.GetHistogram().GetYaxis().GetLabelSize() * 0.75)
	          if var.typeval is TH1I:
	              stack.GetHistogram().GetXaxis().SetNdivisions(tSum.GetNbinsX())
	              stack.GetHistogram().GetXaxis().CenterLabels(True)
	   else:
	      set_fancy_2D_style(57)
              gPad.SetRightMargin(0.2)
	      stack.Draw(var.drawOpt2D)
              if False and var.binsX == var.binsY and not var.binsX == var.bins:
                  xA = stack.GetXaxis().GetBinLowEdge(1)
                  yA = stack.GetYaxis().GetBinLowEdge(1)
                  xB = stack.GetXaxis().GetBinLowEdge(stack.GetXaxis().GetNbins()+1)
                  yB = stack.GetYaxis().GetBinLowEdge(stack.GetYaxis().GetNbins()+1)
                  diagonal = TLine( xA, yA, xB, yB )
                  diagonal.SetLineStyle(2)
                  diagonal.SetLineColor(kBlack)
                  diagonal.SetLineWidth(2)
                  diagonal.Draw("SAME")

              # if var.shortname == "BDTGScore_ttH_ttbarDD_VS_BDTGScore_ttH_ttV":
              #     xA_BDT_0 = -1
              #     yA_BDT_0 = 1
              #     xB_BDT_0 = 1
              #     yB_BDT_0 = -1
              #     diag_BDT_0 = TLine( xA_BDT_0, yA_BDT_0, xB_BDT_0, yB_BDT_0 )
              #     diag_BDT_0.SetLineStyle(2)
              #     diag_BDT_0.SetLineColor(kWhite)
              #     diag_BDT_0.SetLineWidth(2)
              #     diag_BDT_0.Draw("SAME")
              #     leg_BDT_0 = TLatex()
              #     leg_BDT_0.SetTextSize(0.037)
              #     leg_BDT_0.SetNDC()
              #     leg_BDT_0.SetTextColor(kWhite)
              #     leg_BDT_0.SetTextAngle(317.33)
              #     leg_BDT_0.DrawLatex(0.298, 0.808, "#bf{BDTG=0.0}")
              #     xA_BDT_05 = -0.5
              #     yA_BDT_05 = 1
              #     xB_BDT_05 = 1
              #     yB_BDT_05 = -0.5
              #     diag_BDT_05 = TLine( xA_BDT_05, yA_BDT_05, xB_BDT_05, yB_BDT_05 )
              #     diag_BDT_05.SetLineStyle(2)
              #     diag_BDT_05.SetLineColor(kWhite)
              #     diag_BDT_05.SetLineWidth(2)
              #     diag_BDT_05.Draw("SAME")
              #     leg_BDT_05 = TLatex()
              #     leg_BDT_05.SetTextSize(0.037)
              #     leg_BDT_05.SetNDC()
              #     leg_BDT_05.SetTextAngle(317.33)
              #     leg_BDT_05.SetTextColor(kWhite)
              #     leg_BDT_05.DrawLatex(0.367, 0.916, "#bf{BDTG=0.5}")
              #     xA_BDT_m05 = -1
              #     yA_BDT_m05 = 0.5
              #     xB_BDT_m05 = 0.5
              #     yB_BDT_m05 = -1
              #     diag_BDT_m05 = TLine( xA_BDT_m05, yA_BDT_m05, xB_BDT_m05, yB_BDT_m05 )
              #     diag_BDT_m05.SetLineStyle(2)
              #     diag_BDT_m05.SetLineColor(kWhite)
              #     diag_BDT_m05.SetLineWidth(2)
              #     diag_BDT_m05.Draw("SAME")
              #     leg_BDT_m05 = TLatex()
              #     leg_BDT_m05.SetTextSize(0.037)
              #     leg_BDT_m05.SetNDC()
              #     leg_BDT_m05.SetTextAngle(317.33)
              #     leg_BDT_m05.SetTextColor(kWhite)
              #     leg_BDT_m05.DrawLatex(0.228, 0.695, "#bf{BDTG=-0.5}")

              #     # from ROOT import TPaletteAxis
              #     # gPad.Update()
              #     # pal = stack.GetHistogram().GetListOfFunctions().FindObject("palette")
              #     # print type(pal)
              #     # pal.GetAxis().SetLabelSize(0.02)

        if bkg and not isvar2D:
            tSum.Draw("E2 SAME")

        if obs:
            if not isvar2D:
                if stack:
                    datagr.Draw("P SAME")
                else:
                    if logx or var.logaxisX:
                        gPad.SetLogx()
                    if log or var.logaxis:
                        gPad.SetLogy()
                    datagr.GetXaxis().SetTitle(var.latexname)
                    binwidth = (var.maxval - var.minval) / var.bins
                    ytitle = 'Events / %.2g %s' % (binwidth, var.unit) if var.unit else 'Events / bin'
                    datagr.GetYaxis().SetTitle(ytitle)
                    datagr.Draw('AP')
            else:
                set_fancy_2D_style()
                gPad.SetRightMargin(0.2)
                datagr.Draw(var.drawOpt2D)

        if not nolegs:
            lower, labels = self.labels(legs, showratio and obs and bkg)
        #gPad.RedrawAxis()
        c.Update()

        if wait: raw_input('Hit enter to continue...')
        if type(save) is str:
            save = [save]
        for filepath in save:
            c.SaveAs(filepath)
        if not wait:
            gROOT.SetBatch(False)

        return bkg, tSum, obs, sig, stack

    def plotSystematics(self, systematics, var = 'MMC', cut = None, eventweight=None, category = None, overridebackground = None, overflowbins = False, showratio = True, wait = False, save = ['.eps'], log=False, logx=False):

	if not wait:
            gROOT.SetBatch(True)

        cut, category, systematics, overridebackground = self.parseArguments(cut, category, systematics, overridebackground)

        if type(var) is str:
            var = self.vardb.getVar(var)

        isvar2D = ( var.typeval in [TH2,TH2I,TH2F,TH2D] )

        nom, nomlist    = self.sumhist(var, processes=overridebackground, cut=cut, eventweight=eventweight, category=category, systematics=None, systematicsdirection=None, overflowbins=overflowbins)
        up, uplist      = self.sumhist(var, processes=overridebackground, cut=cut, eventweight=eventweight, category=category, systematics=systematics, systematicsdirection='UP', overflowbins=overflowbins)
        down, downlist  = self.sumhist(var, processes=overridebackground, cut=cut, eventweight=eventweight, category=category, systematics=systematics, systematicsdirection='DOWN', overflowbins=overflowbins)
        obs, obslist    = self.sumhist(var, processes=self.observed,      cut=cut, eventweight=eventweight, category=category, systematics=None, systematicsdirection=None, overflowbins=overflowbins)
	sig, siglist    = self.sumhist(var, processes=self.signals,       cut=cut, eventweight=eventweight, category=category, systematics=None, systematicsdirection=None, overflowbins=overflowbins)

        nom_stat_err = nom.Clone("nom_stat_err")
        nom_stat_err.SetFillColor(self.style.get('SumErrorFillColour', kGray+3))
        nom_stat_err.SetLineColor(self.style.get('SumErrorLineColour', 10))
        nom_stat_err.SetFillStyle(self.style.get('SumErrorFillStyle', 3004))
        nom_stat_err.SetMarkerSize(0)

        legs = [
            [nom, 'Nominal' + " ({0:.1f})".format(integrate(nom)), "F"],
            [up, systematics.name + ' + 1#sigma' + " ({0:.1f})".format(integrate(up)), "F"],
            [down, systematics.name + ' - 1#sigma' + " ({0:.1f})".format(integrate(down)), "F"],
        ]

        bkguplist = {}
        bkgdownlist = {}
        for h, process in reversed(uplist):
            h.SetLineWidth(self.style.get('BackgroundLineWidth', 3))
            h.SetLineStyle(self.style.get('BackgroundLineStyle', 1))
            pname = process.__class__.__name__
            h.SetLineColor(self.style.get(pname+'LineColour', self.style.get('BackgroundLineColour', 1)))
            h.SetFillColor(self.style.get(pname+'FillColour', process.colour))
            h.SetFillStyle(self.style.get(pname+'FillStyle', 1001))
            bkguplist[pname] = h
        for h, process in reversed(downlist):
            h.SetLineWidth(self.style.get('BackgroundLineWidth', 3))
            h.SetLineStyle(self.style.get('BackgroundLineStyle', 1))
            pname = process.__class__.__name__
            h.SetLineColor(self.style.get(pname+'LineColour', self.style.get('BackgroundLineColour', 1)))
            h.SetFillColor(self.style.get(pname+'FillColour', process.colour))
            h.SetFillStyle(self.style.get(pname+'FillStyle', 1001))
            bkgdownlist[pname] = h

        if obs:
            process = obslist[0][1]
	    if not ( "$ISDATA$" in obslist[0][0].GetName() ):
                if not isvar2D:
                    datagr = TH1D(obs) # Equivalent to: datagr = makeMCErrors(obs)
                    datagr.SetLineStyle(2)
                else:
                    datagr = TH2D(obs)
	    else :
                if not isvar2D:
                    datagr = makePoissonErrors(obs)
                    datagr.SetMarkerSize(1.2)
                    datagr.SetMarkerStyle(20)
                    datagr.SetMarkerColor(1)
                    datagr.SetLineColor(1)
                else:
                    datagr = TH2D(obs)
            legs.insert(0, (datagr, process.latexname + " ({0:.1f})".format(integrate(obs)), "P"))

        c = TCanvas("c1","Temp",50,50,600,600)
        color_up = 46
        color_down = 36

        if isvar2D:
            showratio = False

        if showratio:
            pad1 = TPad("pad1", "", 0, 0.25, 1, 1)
            pad2 = TPad("pad2", "", 0, 0,   1, 0.25)
            pad1.SetBottomMargin(0.02)
            pad2.SetBottomMargin(0.4)
            #pad2.SetTopMargin(0)
            pad2.SetGridy(1)
            if log or var.logaxis:
                pad1.SetLogy()
            if logx or var.logaxisX:
                pad1.SetLogx()
                pad2.SetLogx()
            pad1.Draw()
            pad2.Draw()
	if not showratio:
            if log or var.logaxis:
                gPad.SetLogy()
            if logx or var.logaxisX:
                gPad.SetLogx()

        for h in nom, down, up:
            h.SetLineStyle(1)
            h.SetFillStyle(0)
            h.SetLineWidth(2)

        down.SetLineColor(color_down)
        up.SetLineColor(color_up)
        nom.SetLineColor(1)
        nom.SetLineWidth(2)
        nom.SetLineStyle(2)

        if showratio:
            if var.typeval is TH1D:
                ratioup = up.Clone("RatioUP")
                ratiodown = down.Clone("RatioDOWN")
                if obs:
                    ratioobs = obs.Clone("RatioOBS")
            elif var.typeval is TH1I:
                ratioup = TH1D("RatioUP", "RatioUP", up.GetNbinsX(), up.GetBinLowEdge(1), up.GetBinLowEdge(up.GetNbinsX()+1))
                ratioup.GetXaxis().SetNdivisions(up.GetNbinsX())
                ratioup.GetXaxis().CenterLabels(True)
                ratiodown = TH1D("RatioDOWN", "RatioDOWN", down.GetNbinsX(), down.GetBinLowEdge(1), down.GetBinLowEdge(down.GetNbinsX()+1))
                if obs:
                    ratioobs = TH1D("RatioOBS", "RatioOBS", obs.GetNbinsX(), obs.GetBinLowEdge(1), obs.GetBinLowEdge(obs.GetNbinsX()+1))
                for i in range(1, up.GetNbinsX()+1):
                    ratioup.SetBinContent(i, up.GetBinContent(i))
                    ratiodown.SetBinContent(i, down.GetBinContent(i))
                    if obs:
                        ratioobs.SetBinContent(i, obs.GetBinContent(i))

            ratioup.SetXTitle(var.latexname)
            ratioup.SetYTitle("Syst/Nom")
            ratioup.GetXaxis().SetTitleSize(0.15)
            ratioup.GetYaxis().SetTitleSize(0.15)
            ratioup.GetXaxis().SetTitleOffset(0.90)
            ratioup.GetYaxis().SetTitleOffset(0.35)
            ratioup.GetXaxis().SetLabelSize(0.15)
            ratioup.GetYaxis().SetLabelSize(0.12)
            ratioup.GetYaxis().SetNdivisions(5)
            ratioup.SetMarkerColor(color_up)
            ratioup.SetMarkerSize(1.0)
            ratioup.SetLineWidth(2)
            ratioup.Divide(nom)

            ratiodown.SetMarkerColor(color_down)
            ratiodown.SetMarkerSize(1.)
            ratiodown.SetLineWidth(2)
            ratiodown.Divide(nom)

            if obs:
                ratioobs.SetMarkerColor(1)
                ratioobs.SetMarkerSize(1.0)
                ratioobs.SetLineWidth(1)
                ratioobs.Divide(nom)

            # Histogram for statistical error bars

	    nom_stat_err_ratio = nom_stat_err.Clone("nom_stat_err_ratio")
	    nom_stat_err_ratio.Divide(nom_stat_err)
            legs.insert(len(legs), (nom_stat_err,"#bf{Stat. unc.}","F"))

            valYmin =  99.
            valYmax = -99.
            for i in range(1, nom.GetXaxis().GetNbins()+1):
                if (ratioup.GetBinContent(i) - ratioup.GetBinError(i)) != 0. and ratioup.GetBinContent(i) < valYmin:
                    valYmin = ratioup.GetBinContent(i) - ratioup.GetBinError(i)
                    if valYmin < 0.:
                        valYmin = 0.
                if (ratioup.GetBinContent(i) + ratioup.GetBinError(i)) != 0. and ratioup.GetBinContent(i) > valYmax:
                    valYmax = ratioup.GetBinContent(i) + ratioup.GetBinError(i)

                if (ratiodown.GetBinContent(i) - ratiodown.GetBinError(i)) != 0. and ratiodown.GetBinContent(i) < valYmin:
                    valYmin = ratiodown.GetBinContent(i) - ratiodown.GetBinError(i)
                    if valYmin < 0.:
                        valYmin = 0.
                if (ratiodown.GetBinContent(i) + ratiodown.GetBinError(i)) != 0. and ratiodown.GetBinContent(i) > valYmax:
                    valYmax = ratiodown.GetBinContent(i) + ratiodown.GetBinError(i)

                if (nom_stat_err_ratio.GetBinContent(i) - nom_stat_err_ratio.GetBinError(i)) != 0. and nom_stat_err_ratio.GetBinContent(i) < valYmin:
                    valYmin = nom_stat_err_ratio.GetBinContent(i) - nom_stat_err_ratio.GetBinError(i)
                    if valYmin < 0.:
                        valYmin = 0.
                if (nom_stat_err_ratio.GetBinContent(i) + nom_stat_err_ratio.GetBinError(i)) != 0. and nom_stat_err_ratio.GetBinContent(i) > valYmax:
                    valYmax = nom_stat_err_ratio.GetBinContent(i) + nom_stat_err_ratio.GetBinError(i)

                if obs:
                    if (ratioobs.GetBinContent(i) - ratioobs.GetBinError(i)) != 0. and ratioobs.GetBinContent(i) < valYmin:
                        valYmin = ratioobs.GetBinContent(i) - ratioobs.GetBinError(i)
                        if valYmin < 0.:
                            valYmin = 0.
                    if (ratioobs.GetBinContent(i) + ratioobs.GetBinError(i)) != 0. and ratioobs.GetBinContent(i) > valYmax:
                        valYmax = ratioobs.GetBinContent(i) + ratioobs.GetBinError(i)

            if abs(1-valYmin) > abs(valYmax-1):
                val = abs(1-valYmin)
            else:
                val = abs(valYmax-1)

            #ratioup.SetMinimum((1-val)-0.1)
            #ratioup.SetMaximum((1+val)+0.1)
            ratioup.SetMinimum(0.0)
            ratioup.SetMaximum(2.0)

            if type(showratio) is tuple:
                if showratio[0] == "MIN" and type(showratio[1]) is float:
                    ratioup.SetMinimum(valYmin - 0.1)
                    ratioup.SetMaximum(showratio[1])
                elif showratio[1] == "MAX" and type(showratio[0]) is float:
                    ratioup.SetMinimum(showratio[0])
                    ratioup.SetMaximum(valYmax + 0.1)
                else:
                    ratioup.SetMinimum(showratio[0])
                    ratioup.SetMaximum(showratio[1])

            pad2.cd()
            ratioup.Draw("HIST")        # Do not draw error bars
            ratiodown.Draw("HIST SAME")	# Do not draw error bars
	    nom_stat_err_ratio.Draw("E2 SAME")
            # ratioobs.Draw("SAME")
            refl = TLine(ratioup.GetBinLowEdge(1), 1., ratioup.GetBinLowEdge(ratioup.GetNbinsX()+1), 1.)
            refl.SetLineStyle(2)
            refl.SetLineColor(kRed)
            refl.Draw("SAME")
            pad1.cd()

        lower, labels = self.labels(legs, showratio)

        # Trick to rescale:

        ymax_new = up.GetMaximum()
        if down and down.GetMaximum() > ymax_new:
            ymax_new = down.GetMaximum()
        if nom and nom.GetMaximum() > ymax_new:
            ymax_new = nom.GetMaximum()
        if obs and obs.GetMaximum() > ymax_new:
            ymax_new = obs.GetMaximum()
        if showratio:
            up.SetMaximum(ymax_new*(2.-lower+0.075))
        else:
            up.SetMaximum(ymax_new*(2.-lower+0.15))
	if log or var.logaxis:
	    up.SetMaximum(up.GetMaximum() * 3*10**(2))

        if showratio:
            up.GetXaxis().SetLabelOffset(999)
            up.GetXaxis().SetLabelSize(0)
            up.GetYaxis().SetTitleSize(up.GetYaxis().GetTitleSize() * 1.2)
            up.GetYaxis().SetTitleOffset(1.20)
            ratioup.GetXaxis().SetTitleSize(ratioup.GetXaxis().GetTitleSize() * 1.2)
            ratioup.GetYaxis().SetTitleSize(ratioup.GetYaxis().GetTitleSize() * 1.2)
        else:
            up.GetXaxis().SetLabelSize(up.GetXaxis().GetLabelSize() * 0.75)
            up.GetYaxis().SetLabelSize(up.GetYaxis().GetLabelSize() * 0.75)
            if var.typeval is TH1I:
                up.GetXaxis().SetNdivisions(up.GetNbinsX())
                up.GetXaxis().CenterLabels(True)

        up.Draw("HIST")
        down.Draw("HIST SAME")
        nom.Draw("HIST SAME")
	nom_stat_err.Draw("E2 SAME")
        if obs:
           datagr.Draw("PE SAME")

        lower, labels = self.labels(legs, showratio)

        if wait: raw_input('Hit enter to continue...')
        if type(save) is str:
            save = [save]
        for filepath in save:
            c.SaveAs(filepath)
        if not wait:
            gROOT.SetBatch(False)

        return obs, nom, up, down, bkguplist, bkgdownlist


def drawText(text, x, y, size=0.05, colour=1):
    l = TLatex()
    l.SetTextSize(size)
    l.SetNDC()
    l.SetTextColor(colour)
    l.DrawLatex(x, y, text)
    return l

def makePoissonErrors(hist):
    removezerobin=False
    npoint=0
    for ibin in range(0,hist.GetNbinsX()+1):
        if(hist.GetBinContent(ibin)>0): npoint += 1

    graph = TGraphAsymmErrors(npoint)
    graph.SetLineWidth(hist.GetLineWidth())
    ipoint=0
    for ibin in range(0,hist.GetNbinsX()+1):
        bincontent=int(hist.GetBinContent(ibin))
        EYlow=bincontent-0.5*TMath.ChisquareQuantile(0.1586555,2.*bincontent)
        EYhigh=0.5*TMath.ChisquareQuantile(1-0.1586555,2.*(bincontent+1))-bincontent
        EX=hist.GetBinWidth(ibin)/2.
        if(bincontent!=0):
            graph.SetPoint(ipoint,hist.GetBinCenter(ibin),bincontent)

            graph.SetPointEXlow(ipoint,EX)
            graph.SetPointEXhigh(ipoint,EX)
            graph.SetPointEYlow(ipoint,EYlow)
            graph.SetPointEYhigh(ipoint,EYhigh)
            ipoint += 1

    return graph

def makeMCErrors(hist):
    removezerobin=False
    npoint=0
    for ibin in range(0,hist.GetNbinsX()+1):
        if(hist.GetBinContent(ibin)>0): npoint += 1

    graph = TGraphAsymmErrors(npoint)
    graph.SetLineWidth(hist.GetLineWidth())
    ipoint=0
    for ibin in range(0,hist.GetNbinsX()+1):
        bincontent= hist.GetBinContent(ibin)
        if(bincontent!=0):
            EY = hist.GetBinError(ibin)
            EX = 0.5 * (hist.GetBinWidth(ibin))
            graph.SetPoint(ipoint,hist.GetBinCenter(ibin),bincontent)
            graph.SetPointError(ipoint, EX, EX, EY, EY)
            ipoint += 1

    return graph

def integrate(hist, mergeOverflow = False):

    # Get effective integral of input histogram, taking u/oflow into account

    integral = -1
    if isinstance(hist,TH2):
        integral = hist.Integral(0,hist.GetXaxis().GetNbins()+1,0,hist.GetYaxis().GetNbins()+1)
    else:
        integral =  hist.Integral(0,hist.GetNbinsX()+1)

    if mergeOverflow:
        integral -= hist.GetBinContent( hist.GetSize() )

    return integral


def SelfDivide( h_self, name = "Ratio", h_props = None ):

    # Divide an histogram by itself, and set the errors correctly

    h = h_self.Clone(name)
    h.Divide(h_self)
    for ibin in range(0,h.GetSize()):
        if h_self.GetBinContent(ibin) > 0:
            h.SetBinError(ibin, h_self.GetBinError(ibin)/h_self.GetBinContent(ibin))

    if "TitleX" in h_props: h.SetXTitle(h_props["TitleX"])
    if "TitleY" in h_props: h.SetYTitle(h_props["TitleY"])
    if "TitleSizeX" in h_props: h.GetXaxis().SetTitleSize(h_props["TitleSizeX"])
    if "TitleSizeY" in h_props: h.GetYaxis().SetTitleSize(h_props["TitleSizeY"])
    if "TitleOffsetX" in h_props: h.GetXaxis().SetTitleOffset(h_props["TitleOffsetX"])
    if "TitleOffsetY" in h_props: h.GetYaxis().SetTitleOffset(h_props["TitleOffsetY"])
    if "LabelSizeX" in h_props: h.GetXaxis().SetLabelSize(h_props["LabelSizeX"])
    if "LabelSizeY" in h_props: h.GetYaxis().SetLabelSize(h_props["LabelSizeY"])
    if "NdivisionsY" in h_props: h.GetYaxis().SetNdivisions(h_props["NdivisionsY"])
    if "FillColor" in h_props: h.SetFillColor(h_props["FillColor"])
    if "FillStyle" in h_props: h.SetFillStyle(h_props["FillStyle"])
    if "LineStyle" in h_props: h.SetLineStyle(h_props["LineStyle"])
    if "LineColor" in h_props: h.SetLineColor(h_props["LineColor"])
    if "LineWidth" in h_props: h.SetLineWidth(h_props["LineWidth"])
    if "MarkerSize" in h_props: h.SetMarkerSize(h_props["MarkerSize"])
    if "MarkerColor" in h_props: h.SetMarkerColor(h_props["MarkerColor"])
    if "MarkerStyle" in h_props: h.SetMarkerStyle(h_props["MarkerStyle"])

    return h


# This function loads the samples metadata from the .csv file info

def loadSamples(inputdir, samplescsv='Files/samples.csv', nomtree='physics', friendtrees=[], friendfile_extension=None, systrees=[], readGFW2=False):

    # The datasat manager takes care of parsing the sample.csv files.

    datasets = DatasetManager.DatasetManager()

    # This returns a list of dictionaries which contain the features for each sample:
    # ID,category,xsection,kfactor,efficiency,name,group,subgroup

    samples = datasets.getListSamples(samplesfile=samplescsv)

    inputs = Inputs()

    if friendtrees:
        inputs.setFriendTree(friendtrees,friendfile_extension)

    inputs.readGFW2 = readGFW2

    for s in samples:

	sampleid = s['ID']
        name     = s['name']
        category = s['category']
        group    = s['group']
        subgroup = s['subgroup']

        separator = '.'
        if not sampleid:
          separator = ''

        filename = inputdir + '/' + group + '/' + sampleid + separator + name + '.root'

	ismc        = ( not category == 'Data' and not group == 'Embedding' )
        isembedding = ( group == 'Embedding' )
        isdata      = ( category == 'Data' )

	inputs.registerTree(filename, nomtree, systrees, ismc, isembedding, isdata, s)

    return inputs

# Function to set nice feats to 2D histograms
# By default it defines a dark rainbow color palette.
# User can pass an enum index to set one of the predefined ROOT palettes
# See:
# https://root.cern.ch/doc/master/classTColor.html#C05
# for a list of available enums.

def set_fancy_2D_style( palette_enum = -1 ):

    icol = 0
    gStyle.SetFrameBorderMode(icol);
    gStyle.SetFrameFillColor(icol);
    gStyle.SetCanvasBorderMode(icol);
    gStyle.SetCanvasColor(icol);
    gStyle.SetPadBorderMode(icol);
    gStyle.SetPadColor(icol);
    gStyle.SetStatColor(icol);
    gStyle.SetOptTitle(0);
    gStyle.SetOptStat(0);
    gStyle.SetOptFit(0);

    if palette_enum == -1:
        ncontours=999
        s = array.array('d', [0.00, 0.34, 0.61, 0.84, 1.00])
        r = array.array('d', [0.00, 0.00, 0.87, 1.00, 0.51])
        g = array.array('d', [0.00, 0.81, 1.00, 0.20, 0.00])
        b = array.array('d', [0.51, 1.00, 0.12, 0.00, 0.00])
        npoints = len(s)
        TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
        gStyle.SetNumberContours(ncontours)
    else:
        gStyle.SetPalette(palette_enum)
