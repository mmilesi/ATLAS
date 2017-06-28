""" Backgrounds_HTopMultilep.py: a set of classes to manage the various processes """

__author__     = "KG, Marco Milesi, Francesco Nuti"
__email__      = "Kong.Guan.Tan@cern.ch, marco.milesi@cern.ch, francesco.nuti@cern.ch"
__maintainer__ = "Marco Milesi"

import os, sys, math, types

sys.path.append(os.path.abspath(os.path.curdir))

from Plotter.BackgroundTools import loadSamples, drawText, Category, Background, Process, VariableDB, Variable, Cut, Systematics, Category

from ROOT import TColor, kBlack, kWhite, kGray, kBlue, kRed, kYellow, kGreen, kAzure, kTeal, kSpring, kOrange, kCyan, TLegend, TLatex, TCanvas, TH1I, TFile

class MyCategory(Category):

    def __init__(self, name, cut=None, controlcut=None, weight=None, overridebins=None, ratiolims=None):
        self.name = name
        self.tokens = name.split(' ')
        self.cut = cut
        self.controlcut = controlcut
	self.weight = weight
        self.overridebins = overridebins
        self.stream = '*'
        self.streamname = '*'
        self.jet = name
        self.met = 'HMET'
        self.ratiolims = ratiolims

        if len(self.tokens)>=3:
            self.stream = self.tokens[0]
            if   self.stream == 'El': self.streamname = 'Egamma'
            elif self.stream == 'Mu': self.streamname = 'Muons'
            self.jet = self.tokens[1]
            self.met = self.tokens[2]

class TTHBackgrounds(Background):

    backgrounds     = []
    sub_backgrounds = []
    signals         = []
    observed        = []
    luminosity  = 1.0
    lumi_units  = 'fb-1'
    norm_factor = 1.0           # Might be needed to correct units for Xsec weight
    rescaleXsecAndLumi = False  # Set to "True" if you don't want to take into account the Xsec*lumi weight
    channel     = '2LSS'        # Can be one of ['2LSS' '2LSS_CR', '3L']
    eventweight = 1.0
    noWeights = False
    useEmbedding    = False
    useZCorrections = False
    useSherpaNNPDF30NNLO = False

    RQCD = {
        'El': (1.00, 0.05),
        'Mu': (1.11, 0.08),
    }

    theta = {
        'El': (999.0, 0.0),
        'Mu': (999.0, 0.0),
    }

    theta_MC = {
        'El': (999.0, 0.0),
        'Mu': (999.0, 0.0),
    }

    def getFirstSimulatedProc(self, category, mybackgrounds = None):

        # Returns the name and position in bkg list of the first background process which
        # is *not* data-driven

        backgrounds = self.backgrounds
        if mybackgrounds:
            backgrounds = mybackgrounds
	for idx, samplename in enumerate(mybackgrounds):
            subprocess = self.getProcess(samplename, category)
            if not ("$ISDATA$") in subprocess.name:
	        return idx, samplename
	return (-1,None)

    def str_to_class(self, field):
        try:
            identifier = getattr(self, field)
        except AttributeError:
            raise NameError("%s doesn't exist." % field)
        if isinstance(identifier, (types.ClassType, types.TypeType)):
            return identifier
        raise TypeError("%s is not a class." % field)


    def applyKfactor(self, sp, category, kfactor, options):

        # Multiply a sample by a kfactor treating properly the systematics if needed

        systematics = options.get('systematics', None)
        systematicsdirection = options.get('systematicsdirection', None)

        kfactor = dict(kfactor)
        if systematics:
            if systematicsdirection == 'UP':
                systdir = 1.0
            elif systematicsdirection == 'DOWN':
                systdir = -1.0
            for key in kfactor:
                kfactor[key] = (kfactor[key][0] + kfactor[key][1]*systdir, kfactor[key][1])

        for key in kfactor:
            ssp = sp.subprocess()
            ssp *= kfactor[key][0]
        return ssp


    def applyRQCD(self, sp, category, options):
        return self.applyKfactor(sp, category, self.RQCD, options)

    def applyTightSF(self, sp):

	# Uses the new scale factors for BDT tight instead of the old BDT medium. ATTENTION the correction is done only with the nominal coefficients.
	# To have the right coefficients for systematics shifted correct the corrections code and produce again the ntuples

        sp1ploweta = sp.subprocess(cut=self.vardb.getCuts(['TauEta00to15', 'OneProng'])) * 0.941 / 0.992
        sp1phigheta = sp.subprocess(cut=self.vardb.getCuts(['TauEta15to25', 'OneProng'])) * 0.89 / 0.952
        sp3ploweta = sp.subprocess(cut=self.vardb.getCuts(['TauEta00to15', 'ThreeProng'])) * 1.006 / 1.073
        sp3phigheta = sp.subprocess(cut=self.vardb.getCuts(['TauEta15to25', 'ThreeProng'])) * 1.082 / 0.99
        sp_corr=sp1ploweta+sp1phigheta+sp3ploweta+sp3phigheta

	return sp_corr

    def labels(self, legs, showratio):
        scale = 1.
        if not showratio:
            scale = 0.75

        if len(legs) < 5:
            scale *= 1.4
            mid = len(legs)
            high = len(legs)
            lower = 0.92 - 0.04*high
            leg1 = TLegend(0.60,lower,0.90,0.92)
            leg2 = None
        else:
            #scale *= 1.2
            mid = int(len(legs)/2)
            high = math.ceil(len(legs)/2)
            lower = 0.92 - 0.04*high
            leg1 = TLegend(0.45,lower,0.65,0.92)
            #leg2 = TLegend(0.70,lower,0.90,0.92)
            leg2 = TLegend(0.65,lower,0.85,0.92)
        for leg in [leg1, leg2]:
            if not leg: continue
            leg.SetFillColor(0)
            leg.SetFillStyle(0)
            leg.SetLineColor(10)
            leg.SetShadowColor(kWhite)
            leg.SetTextSize(0.03 * scale)
            #leg.SetTextSize(0.03)
            leg.SetBorderSize(0)
            #leg.SetEntrySeparation(1.0)

        for l in legs[:mid]:
            leg1.AddEntry(l[0], l[1], l[2])
        leg1.Draw()
        if leg2:
            for l in legs[mid:]:
                leg2.AddEntry(l[0], l[1], l[2])
            leg2.Draw()

        # for O(fb-1) luminosity

        if self.lumi_units == 'fb-1':
            if self.readGFW2:
                lumtext = drawText(text="  #int L dt = %.1f fb^{-1}"%(self.luminosity/1e3), x=.2, y=.87, size=0.03 * scale)
            else:
                lumtext = drawText(text="  #int L dt = %.1f fb^{-1}"%(self.luminosity), x=.2, y=.87, size=0.03 * scale)

        # for O(pb-1) luminosity

        elif self.lumi_units == 'pb-1':
            lumtext = drawText(text="  #int L dt = %.2f pb^{-1}"%(self.luminosity*1000), x=.2, y=.87, size=0.04 * scale)

        cmetext   = drawText(text="         #sqrt{s} = 13 TeV", x=.2, y=.82, size=0.03 * scale)
        atlastext = drawText(text="#bf{#it{ATLAS}} Work In Progress", x=.2, y=.77, size=0.03 * scale)

        return lower, locals()

    def colours(self):
        cache = []
        cache.append(TColor(1040, 82/255., 195/255., 229/255.))     # light blue
        cache.append(TColor(1041, 66/255., 156/255., 183/255.))     # blue
        cache.append(TColor(1042, 44/255., 96/255., 125/255.))      # dark blue
        cache.append(TColor(1030, 229/255., 229/255., 121/255.))    # yellow
        cache.append(TColor(1020, 220/255., 87/255., 60/255.))      # red
        cache.append(TColor(1010, 103/255., 73/255., 130/255.))     # purple
        cache.append(TColor(1000, 108/255., 178/255., 81/255.))     # green
        return cache

    class Observed(Process):

        latexname = 'Data'

        # This method contains the instruction of which tree to load.
	# NB: this is not automatically executed when the instance of the class is created. In fact, it is executed via __call__

        def base(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

            inputgroup = [
                    ('Data', 'physics_Main'),
                ]

            if self.parent.readGFW2:
                del inputgroup[:]
                inputgroup = [ ('Data', 'fakes_mm') ]

            print("\n{0}:\n".format(self.__class__.__name__))
	    print("\n".join("{0} - {1} : {2}".format(idx,sample[0],sample[1]) for idx, sample in enumerate(inputgroup)))

            trees = self.inputs.getTrees(treename, inputgroup)
            return self.subprocess(trees=trees)

        # This method makes the class instance "callable" as a function:
	# E.g.:
	#
	# >>> obs = Observed(myproc)
	# >>> obs(mytreename,mycategory,myoptions)
	# ---> __call__ is executed
	#
        # Here, in addition to the tree loading, the selection, application of SFs and systematics are applied.

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

            sp = self.base(treename, category, options)

            weight = None
            TTcut  = ('','TT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]

            # Clean up from any truth cut

	    basecut = category.cut

	    for CUT in basecut.cutlist:
	    	if ("TRUTH") in CUT.cutname:
	    	    basecut = basecut.removeCut(CUT)

            # basecut = basecut & Cut('2LSS0TauFlag','is2LSS0Tau')

            sp = sp.subprocess(cut=basecut,eventweight=weight)

            if TTcut:
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            # logistic_NN_top  = "output_top"
            # logistic_RNN_top = "RNN_output_top"
            logistic_NN_top  = "1/(1+TMath::Exp(-600*(output_top-0.984176)))"
            logistic_RNN_top = "1/(1+TMath::Exp(-100*(RNN_output_top-0.968906)))"
            dist_NN  = "( TMath::Sqrt(2.0) - TMath::Sqrt( (1.0-{0})*(1.0-{0}) + (1.0-output_ttV)*(1.0-output_ttV) ) )/TMath::Sqrt(2.0)".format(logistic_NN_top)
            dist_RNN = "( TMath::Sqrt(2.0) - TMath::Sqrt( (1.0-{0})*(1.0-{0}) + (1.0-RNN_output_ttV)*(1.0-RNN_output_ttV) ) )/TMath::Sqrt(2.0)".format(logistic_RNN_top)

            if   self.parent.var.shortname == "NN_ttV": sp = sp.subprocess(cut=Cut("NN_Cut","( output_ttV < 0.5 )"))
            if   self.parent.var.shortname == "NN_top": sp = sp.subprocess(cut=Cut("NN_Cut","( {0} < 0.9 )".format(logistic_NN_top)))
            if   self.parent.var.shortname == "RNN_ttV": sp = sp.subprocess(cut=Cut("NN_Cut","( RNN_output_ttV < 0.5 )"))
            if   self.parent.var.shortname == "RNN_top": sp = sp.subprocess(cut=Cut("NN_Cut","( {0} < 0.9 )".format(logistic_RNN_top)))
            if   self.parent.var.shortname == "NN_Rebinned": sp = sp.subprocess(cut=Cut("NN_Cut","( output_bin <= 3 )"))
            if   self.parent.var.shortname == "RNN_Rebinned": sp = sp.subprocess(cut=Cut("NN_Cut","( RNN_output_bin <= 3 )"))
            if   self.parent.var.shortname == "NNComb": sp = sp.subprocess(cut=Cut("NN_Cut","( {0} < 0.7 )".format(dist_NN)))
            if   self.parent.var.shortname == "RNNComb": sp = sp.subprocess(cut=Cut("NN_Cut","( {0}  < 0.7 )".format(dist_RNN)))

            print("\n{0} - cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp.basecut.cutnamelist, weight))

            return sp


    class TTBarH(Process):

        latexname = 't#bar{t} H'
        colour = kBlack

        def base(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    # hmass = 300 is the default value if hmass is not in options.
	    # hmass can be specified in the option "signal" passed to the main plot function in the plotting script

            # hmass = options.get('hmass', '125')

            inputgroup = [
                ('ttH', 'ttH_dil_Pythia8'),
                ('ttH', 'ttH_semilep_Pythia8'),
                ('ttH', 'ttH_allhad_Pythia8'),
                        ]

            print("\n{0}:\n".format(self.__class__.__name__))
	    print("\n".join("{0} - {1} : {2}".format(idx,sample[0],sample[1]) for idx, sample in enumerate(inputgroup)))

            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) * self.parent.norm_factor
            return sp

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)

            weight = None
            TTcut  = ('','TT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

            if TTcut:
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\n{0} - cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp.basecut.cutnamelist, weight))

            return sp


    class TTBarHSemilep(Process):

        latexname = 't#bar{t} H (semilep)'
        colour = kBlack

        def base(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    # hmass = 300 is the default value if hmass is not in options.
	    # hmass can be specified in the option "signal" passed to the main plot function in the plotting script

            # hmass = options.get('hmass', '125')

            inputgroup = [
                ('ttH', 'ttH_semilep_Pythia8'),
                        ]

            print("\n{0}:\n".format(self.__class__.__name__))
	    print("\n".join("{0} - {1} : {2}".format(idx,sample[0],sample[1]) for idx, sample in enumerate(inputgroup)))

            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) * self.parent.norm_factor
            return sp

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)

            weight = None
            TTcut  = ('','TT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

            if TTcut:
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\n{0} - cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp.basecut.cutnamelist, weight))

            return sp


    class TTBarHDilep(Process):

        latexname = 't#bar{t} H (dilep)'
        colour = kBlack

        def base(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

            inputgroup = [
                ('ttH', 'ttH_dil_Pythia8'),
                         ]

            print("\n{0}:\n".format(self.__class__.__name__))
	    print("\n".join("{0} - {1} : {2}".format(idx,sample[0],sample[1]) for idx, sample in enumerate(inputgroup)))

            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) * self.parent.norm_factor
            return sp

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)

            weight = None
            TTcut  = ('','TT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

            if TTcut:
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\n{0} - cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp.basecut.cutnamelist, weight))

            return sp


    class THbj(Process):

        latexname = 'tHbj'
        colour = kBlack

        def base(self, treename='physics', category=None, options={}):

            inputgroup = [
                ('tHbj', '*'),
                        ]

            print("\n{0}:\n".format(self.__class__.__name__))
	    print("\n".join("{0} - {1} : {2}".format(idx,sample[0],sample[1]) for idx, sample in enumerate(inputgroup)))

            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) * self.parent.norm_factor
            return sp

        def __call__(self, treename='physics', category=None, options={}):

            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)

            basecut = category.cut

            if all( c in category.cut.cutname for c in ["2Lep_TRUTH_ProbePromptEvent","2Lep_TRUTH_QMisIDVeto"] ):
                basecut = basecut.swapCut(self.vardb.getCut('2Lep_TRUTH_ProbePromptEvent'), self.vardb.getCut('2Lep_TRUTH_PurePromptEvent'))
                basecut = basecut.removeCut(self.vardb.getCut("2Lep_TRUTH_QMisIDVeto"))

            weight = None
            TTcut  = ('','TT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]

            sp = sp.subprocess(cut=basecut,eventweight=weight)

            if TTcut:
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\n{0} - cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp.basecut.cutnamelist, weight))

            return sp

    class WtH(Process):

        latexname = 'WtH'
        colour = kBlack

        def base(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

            inputgroup = [
                ('tWH', '*'),
                        ]

            print("\n{0}:\n".format(self.__class__.__name__))
	    print("\n".join("{0} - {1} : {2}".format(idx,sample[0],sample[1]) for idx, sample in enumerate(inputgroup)))

            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) * self.parent.norm_factor
            return sp

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)

            basecut = category.cut

            if all( c in category.cut.cutname for c in ["2Lep_TRUTH_ProbePromptEvent","2Lep_TRUTH_QMisIDVeto"] ):
                basecut = basecut.swapCut(self.vardb.getCut('2Lep_TRUTH_ProbePromptEvent'), self.vardb.getCut('2Lep_TRUTH_PurePromptEvent'))
                basecut = basecut.removeCut(self.vardb.getCut("2Lep_TRUTH_QMisIDVeto"))

            weight = None
            TTcut  = ('','TT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]

            sp = sp.subprocess(cut=basecut,eventweight=weight)

            if TTcut:
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\n{0} - cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp.basecut.cutnamelist, weight))

            return sp

    class TTBarW(Process):

        latexname = 't#bar{t} W'
        colour = kRed - 4

        def base(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    inputgroup = [
		('tops', 'ttW_aMcAtNlo'),
                         ]

            print("\n{0}:\n".format(self.__class__.__name__))
	    print("\n".join("{0} - {1} : {2}".format(idx,sample[0],sample[1]) for idx, sample in enumerate(inputgroup)))

            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) * self.parent.norm_factor
            return sp

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)

            basecut = category.cut

            if all( c in category.cut.cutname for c in ["2Lep_TRUTH_ProbePromptEvent","2Lep_TRUTH_QMisIDVeto"] ):
                basecut = basecut.swapCut(self.vardb.getCut('2Lep_TRUTH_ProbePromptEvent'), self.vardb.getCut('2Lep_TRUTH_PurePromptEvent'))
                basecut = basecut.removeCut(self.vardb.getCut("2Lep_TRUTH_QMisIDVeto"))

            weight = None
            TTcut  = ('','TT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]

            sp = sp.subprocess(cut=basecut,eventweight=weight)

            if TTcut:
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\n{0} - cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp.basecut.cutnamelist, weight))

            return sp


    class TTBarZ(Process):

        latexname = 't#bar{t} Z'
        colour = kRed - 7

        def base(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    inputgroup = [
                ('tops', 'ttee_aMcAtNlo'),
                ('tops', 'ttmumu_aMcAtNlo'),
                ('tops', 'tttautau_aMcAtNlo'),
                ('tops', 'ttZnunu_aMcAtNlo'),
                ('tops', 'ttZqq_aMcAtNlo'),
                         ]

            print("\n{0}:\n".format(self.__class__.__name__))
	    print("\n".join("{0} - {1} : {2}".format(idx,sample[0],sample[1]) for idx, sample in enumerate(inputgroup)))

            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) * self.parent.norm_factor
            return sp

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)

            basecut = category.cut

            if all( c in category.cut.cutname for c in ["2Lep_TRUTH_ProbePromptEvent","2Lep_TRUTH_QMisIDVeto"] ):
                basecut = basecut.swapCut(self.vardb.getCut('2Lep_TRUTH_ProbePromptEvent'), self.vardb.getCut('2Lep_TRUTH_PurePromptEvent'))
                basecut = basecut.removeCut(self.vardb.getCut("2Lep_TRUTH_QMisIDVeto"))

            weight = None
            TTcut  = ('','TT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]

            sp = sp.subprocess(cut=basecut,eventweight=weight)

            if TTcut:
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\n{0} - cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp.basecut.cutnamelist, weight))

            return sp

    class TTBarGamma(Process):

        latexname = 't#bar{t}#gamma'
        colour = kRed - 6

        def base(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    inputgroup = [
		('tops', 'ttgamma'),
                         ]

            print("\n{0}:\n".format(self.__class__.__name__))
	    print("\n".join("{0} - {1} : {2}".format(idx,sample[0],sample[1]) for idx, sample in enumerate(inputgroup)))

            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) * self.parent.norm_factor
            return sp

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)

            weight = None
            TTcut  = ('','TT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

            if TTcut:
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            tt_ttgamma_OLR_cut = Cut("TTBar_TTGamma_OLR","( ( mc_channel_number == 410082 && m_MEphoton_OLtty_cat1 ) || ( mc_channel_number != 410082 && m_MEphoton_OLtty_keepEvent && !m_hasMEphoton_DRgt02_nonhad ) )")
            sp = sp.subprocess(cut=tt_ttgamma_OLR_cut)

            print("\n{0} - cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp.basecut.cutnamelist, weight))

            return sp


    class Zeejets(Process):

        latexname = 'Z/#gamma*#rightarrow#it{ee}+jets'
        colour = kGreen-7

        def base(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    inputgroup = [
                           ('Z+jetsLowMllBVeto', 'ee'),
                           ('Z+jetsLowMllBFilter', 'ee'),
                           ('Z+jets_HighZPt', 'ee'),
                         ]

            if self.parent.useSherpaNNPDF30NNLO:
                # Sherpa NNPDF30NNLO
                print("\nUsing Sherpa NNPDF30NNLO Z+jets! Jet reweighting needed...\n")
                inputgroup += [ ('Z+jetsCVetoBVeto_NNPDF30NNLO', 'ee'), ('Z+jetsCFilterBVeto_NNPDF30NNLO', 'ee'), ('Z+jetsBFilter_NNPDF30NNLO', 'ee') ]
            else:
                # Sherpa CT10
                inputgroup += [ ('Z+jetsCVetoBVeto', 'ee'), ('Z+jetsCFilterBVeto', 'ee'), ('Z+jetsBFilter', 'ee') ]

            print("\n{0}:\n".format(self.__class__.__name__))
	    print("\n".join("{0} - {1} : {2}".format(idx,sample[0],sample[1]) for idx, sample in enumerate(inputgroup)))

            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) * self.parent.norm_factor
            if self.parent.useZCorrections:
                sp = sp*0.94
            return sp

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)

            basecut = category.cut

            if all( c in category.cut.cutname for c in ["2Lep_TRUTH_ProbePromptEvent","2Lep_TRUTH_QMisIDVeto"] ):
                basecut = basecut.swapCut(self.vardb.getCut('2Lep_TRUTH_ProbePromptEvent'), self.vardb.getCut('2Lep_TRUTH_PurePromptEvent'))
                basecut = basecut.removeCut(self.vardb.getCut("2Lep_TRUTH_QMisIDVeto"))

            weight = None
            TTcut  = ('','TT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]

            if self.parent.useSherpaNNPDF30NNLO:
                weight = 'SherpaNJetWeight'

            sp = sp.subprocess(cut=basecut,eventweight=weight)

            if TTcut:
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\n{0} - cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp.basecut.cutnamelist, weight))

            return sp


    class Zmumujets(Process):

        latexname = 'Z/#gamma*#rightarrow#mu#mu+jets'
        colour = kTeal+2

        def base(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    inputgroup = [
                           ('Z+jetsLowMllBVeto', 'mumu'),
                           ('Z+jetsLowMllBFilter', 'mumu'),
                           ('Z+jets_HighZPt', 'mumu'),
                         ]

            if self.parent.useSherpaNNPDF30NNLO:
                # Sherpa NNPDF30NNLO
                print("\nUsing Sherpa NNPDF30NNLO Z+jets! Jet reweighting needed...\n")
                inputgroup += [ ('Z+jetsCVetoBVeto_NNPDF30NNLO', 'mumu'), ('Z+jetsCFilterBVeto_NNPDF30NNLO', 'mumu'), ('Z+jetsBFilter_NNPDF30NNLO', 'mumu') ]
            else:
                # Sherpa CT10
                inputgroup += [ ('Z+jetsCVetoBVeto', 'mumu'), ('Z+jetsCFilterBVeto', 'mumu'), ('Z+jetsBFilter', 'mumu') ]

            print("\n{0}:\n".format(self.__class__.__name__))
	    print("\n".join("{0} - {1} : {2}".format(idx,sample[0],sample[1]) for idx, sample in enumerate(inputgroup)))

            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) * self.parent.norm_factor
            if self.parent.useZCorrections:
                sp = sp*0.94
            return sp

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)

            basecut = category.cut

            if all( c in category.cut.cutname for c in ["2Lep_TRUTH_ProbePromptEvent","2Lep_TRUTH_QMisIDVeto"] ):
                basecut = basecut.swapCut(self.vardb.getCut('2Lep_TRUTH_ProbePromptEvent'), self.vardb.getCut('2Lep_TRUTH_PurePromptEvent'))
                basecut = basecut.removeCut(self.vardb.getCut("2Lep_TRUTH_QMisIDVeto"))

            weight = None
            TTcut  = ('','TT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]

            if self.parent.useSherpaNNPDF30NNLO:
                weight = 'SherpaNJetWeight'

            sp = sp.subprocess(cut=basecut,eventweight=weight)

            if TTcut:
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\n{0} - cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp.basecut.cutnamelist, weight))

            return sp

    class Ztautaujets(Process):

        latexname = 'Z/#gamma*#rightarrow#tau#tau+jets'
        colour = kTeal

        def base(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    inputgroup = [
                           ('Z+jetsLowMllBVeto', 'tautau'),
                           ('Z+jetsLowMllBFilter', 'tautau'),
                           ('Z+jets_HighZPt', 'tautau'),
                         ]

            if self.parent.useSherpaNNPDF30NNLO:
                # Sherpa NNPDF30NNLO
                print("\nUsing Sherpa NNPDF30NNLO Z+jets! Jet reweighting needed...\n")
                inputgroup += [ ('Z+jetsCVetoBVeto_NNPDF30NNLO', 'tautau'), ('Z+jetsCFilterBVeto_NNPDF30NNLO', 'tautau'), ('Z+jetsBFilter_NNPDF30NNLO', 'tautau') ]
            else:
                # Sherpa CT10
                inputgroup += [ ('Z+jetsCVetoBVeto', 'tautau'), ('Z+jetsCFilterBVeto', 'tautau'), ('Z+jetsBFilter', 'tautau') ]

            print("\n{0}:\n".format(self.__class__.__name__))
	    print("\n".join("{0} - {1} : {2}".format(idx,sample[0],sample[1]) for idx, sample in enumerate(inputgroup)))

            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) * self.parent.norm_factor
            if self.parent.useZCorrections:
                sp = sp*0.94
            return sp

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)

            basecut = category.cut

            if all( c in category.cut.cutname for c in ["2Lep_TRUTH_ProbePromptEvent","2Lep_TRUTH_QMisIDVeto"] ):
                basecut = basecut.swapCut(self.vardb.getCut('2Lep_TRUTH_ProbePromptEvent'), self.vardb.getCut('2Lep_TRUTH_PurePromptEvent'))
                basecut = basecut.removeCut(self.vardb.getCut("2Lep_TRUTH_QMisIDVeto"))

            weight = None
            TTcut  = ('','TT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]

            if self.parent.useSherpaNNPDF30NNLO:
                weight = 'SherpaNJetWeight'

            sp = sp.subprocess(cut=basecut,eventweight=weight)

            if TTcut:
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\n{0} - cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp.basecut.cutnamelist, weight))

            return sp


    class Zjets(Process):

        latexname = 'Z/#gamma* + jets'
        colour = kGreen

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

            print("\n{0}:\n".format(self.__class__.__name__))

            zee     = self.parent.procmap['Zeejets'](treename, category, options)
            zmumu   = self.parent.procmap['Zmumujets'](treename, category, options)
            ztautau = self.parent.procmap['Ztautaujets'](treename, category, options)

            return zee + zmumu + ztautau


    class ZjetsCF(Process):

        latexname = 'Z/#gamma*+jets (QMisID only)'
        colour = kAzure + 10

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

            print("\n{0}:\n".format(self.__class__.__name__))

            # Pick the subprocesses from the DB

            zeejets     = self.parent.procmap['Zeejets'](treename, category, options)
            zmumujets   = self.parent.procmap['Zmumujets'](treename, category, options)
            ztautaujets = self.parent.procmap['Ztautaujets'](treename, category, options)

            basecut = category.cut

            # Plot only events where at least one lepton is charge flip.
            # Do it only for SS events: for other events, just apply a weight = 0 to kill the process

            basecut = category.cut
            weight  = None
            #if ("2Lep_SS") in basecut.cutname:
            if True:
                for CUT in basecut.cutlist:
                    if ("TRUTH") in CUT.cutname:
        	        basecut = basecut.removeCut(CUT)
            else:
                weight = '0.0'

            updatedcut = basecut & self.vardb.getCut('2Lep_TRUTH_QMisIDEvent')

            # Remember to add TT cut (if needed) as it's not in basecut!

            TTcut  = ('','TT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]
            if TTcut:
                updatedcut = updatedcut & self.vardb.getCut(TTcut)

            # Reset the cut for the subprocesses, and apply the new one

            zeejets     = zeejets.subprocess(cut=updatedcut, clearbasecut=True, eventweight=weight)
            zmumujets   = zmumujets.subprocess(cut=updatedcut, clearbasecut=True, eventweight=weight)
            ztautaujets = ztautaujets.subprocess(cut=updatedcut, clearbasecut=True, eventweight=weight)

            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("Zeejets",zeejets.basecut.cutnamelist, zeejets.eventweight))
            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("Zmumujets",zmumujets.basecut.cutnamelist, zmumujets.eventweight))
            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("Ztautaujets",ztautaujets.basecut.cutnamelist, ztautaujets.eventweight))

            return zeejets + zmumujets + ztautaujets


    class Wenujets(Process):

        latexname = 'W#rightarrow#it{e}#nu+jets'
        colour = kYellow

        def base(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    inputgroup = [
                ('W+jetsBFilter', 'enu'),
                ('W+jetsCFilterBVeto', 'enu'),
                ('W+jetsCVetoBVeto', 'enu'),
                ('W+jets_HighWPt', 'enu')
                         ]

            print("\n{0}:\n".format(self.__class__.__name__))
	    print("\n".join("{0} - {1} : {2}".format(idx,sample[0],sample[1]) for idx, sample in enumerate(inputgroup)))

            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) * self.parent.norm_factor
            return sp

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)

            basecut = category.cut

            if all( c in category.cut.cutname for c in ["2Lep_TRUTH_ProbePromptEvent","2Lep_TRUTH_QMisIDVeto"] ):
                basecut = basecut.swapCut(self.vardb.getCut('2Lep_TRUTH_ProbePromptEvent'), self.vardb.getCut('2Lep_TRUTH_PurePromptEvent'))
                basecut = basecut.removeCut(self.vardb.getCut("2Lep_TRUTH_QMisIDVeto"))

            weight = None
            TTcut  = ('','TT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]

            sp = sp.subprocess(cut=basecut,eventweight=weight)

            if TTcut:
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\n{0} - cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp.basecut.cutnamelist, weight))

            return sp


    class Wmunujets(Process):

        latexname = 'W#rightarrow#mu#nu+jets'
        colour = kYellow

        def base(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    inputgroup = [
                ('W+jetsBFilter', 'munu'),
                ('W+jetsCFilterBVeto', 'munu'),
                ('W+jetsCVetoBVeto', 'munu'),
                ('W+jets_HighWPt', 'munu')
                         ]

            print("\n{0}:\n".format(self.__class__.__name__))
	    print("\n".join("{0} - {1} : {2}".format(idx,sample[0],sample[1]) for idx, sample in enumerate(inputgroup)))

            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) * self.parent.norm_factor
            return sp

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)

            basecut = category.cut

            if all( c in category.cut.cutname for c in ["2Lep_TRUTH_ProbePromptEvent","2Lep_TRUTH_QMisIDVeto"] ):
                basecut = basecut.swapCut(self.vardb.getCut('2Lep_TRUTH_ProbePromptEvent'), self.vardb.getCut('2Lep_TRUTH_PurePromptEvent'))
                basecut = basecut.removeCut(self.vardb.getCut("2Lep_TRUTH_QMisIDVeto"))

            weight = None
            TTcut  = ('','TT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]

            sp = sp.subprocess(cut=basecut,eventweight=weight)

            if TTcut:
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\n{0} - cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp.basecut.cutnamelist, weight))

            return sp


    class Wtaunujets(Process):

        latexname = 'W#rightarrow#tau#nu+jets'
        colour = kYellow

        def base(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    inputgroup = [
                ('W+jetsBFilter', 'taunu'),
                ('W+jetsCFilterBVeto', 'taunu'),
                ('W+jetsCVetoBVeto', 'taunu'),
                ('W+jets_HighWPt', 'taunu')
                         ]

            print("\n{0}:\n".format(self.__class__.__name__))
	    print("\n".join("{0} - {1} : {2}".format(idx,sample[0],sample[1]) for idx, sample in enumerate(inputgroup)))

            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) * self.parent.norm_factor
            return sp

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)

            basecut = category.cut

            if all( c in category.cut.cutname for c in ["2Lep_TRUTH_ProbePromptEvent","2Lep_TRUTH_QMisIDVeto"] ):
                basecut = basecut.swapCut(self.vardb.getCut('2Lep_TRUTH_ProbePromptEvent'), self.vardb.getCut('2Lep_TRUTH_PurePromptEvent'))
                basecut = basecut.removeCut(self.vardb.getCut("2Lep_TRUTH_QMisIDVeto"))

            weight = None
            TTcut  = ('','TT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]

            sp = sp.subprocess(cut=basecut,eventweight=weight)

            if TTcut:
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\n{0} - cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp.basecut.cutnamelist, weight))

            return sp

    class Wjets(Process):

        latexname = 'W+jets'
        colour = kYellow

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

            print("\n{0}:\n".format(self.__class__.__name__))

            wenu   = self.parent.procmap['Wenujets'](treename, category, options)
            wmunu  = self.parent.procmap['Wmunujets'](treename, category, options)
            wtaunu = self.parent.procmap['Wtaunujets'](treename, category, options)

            return (wenu + wmunu + wtaunu)


    class RareTop(Process):

        latexname = 'tZ, 4t, ttWW, tWZ'
        #latexname = 'rare top'

        colour = kAzure + 1

        def base(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    inputgroup = [
                ('tops', 'tZ'),
                ('tops', '3top'),
                ('tops', '4top'),
                ('tops', 'ttWW'),
                ('tops', 'tWZDR'),
                         ]

            print("\n{0}:\n".format(self.__class__.__name__))
	    print("\n".join("{0} - {1} : {2}".format(idx,sample[0],sample[1]) for idx, sample in enumerate(inputgroup)))

            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) * self.parent.norm_factor
            return sp

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)

            basecut = category.cut

            if all( c in category.cut.cutname for c in ["2Lep_TRUTH_ProbePromptEvent","2Lep_TRUTH_QMisIDVeto"] ):
                basecut = basecut.swapCut(self.vardb.getCut('2Lep_TRUTH_ProbePromptEvent'), self.vardb.getCut('2Lep_TRUTH_PurePromptEvent'))
                basecut = basecut.removeCut(self.vardb.getCut("2Lep_TRUTH_QMisIDVeto"))

            weight = None
            TTcut  = ('','TT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]

            sp = sp.subprocess(cut=basecut,eventweight=weight)

            if TTcut:
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\n{0} - cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp.basecut.cutnamelist, weight))

            return sp


    class SingleTop(Process):

        latexname = 'single t, tW'
        colour = kAzure + 1

        def base(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    inputgroup = [
                ('tops', 'singlet'),
                ('tops', 'SingleTopSchan_noAllHad'),
                ('tops', 'tW'),
                         ]

            print("\n{0}:\n".format(self.__class__.__name__))
	    print("\n".join("{0} - {1} : {2}".format(idx,sample[0],sample[1]) for idx, sample in enumerate(inputgroup)))

            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) * self.parent.norm_factor
            return sp

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)

            basecut = category.cut

            if all( c in category.cut.cutname for c in ["2Lep_TRUTH_ProbePromptEvent","2Lep_TRUTH_QMisIDVeto"] ):
                basecut = basecut.swapCut(self.vardb.getCut('2Lep_TRUTH_ProbePromptEvent'), self.vardb.getCut('2Lep_TRUTH_PurePromptEvent'))
                basecut = basecut.removeCut(self.vardb.getCut("2Lep_TRUTH_QMisIDVeto"))

            weight = None
            TTcut  = ('','TT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]

            sp = sp.subprocess(cut=basecut,eventweight=weight)

            if TTcut:
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\n{0} - cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp.basecut.cutnamelist, weight))

            return sp

    class Rare(Process):

        latexname = 'Others'
        colour = kGray

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

            print("\n{0}:\n".format(self.__class__.__name__))

            rare_top = self.parent.procmap['RareTop'](treename, category, options)
            triboson = self.parent.procmap['Triboson'](treename, category, options)
            tHbj     = self.parent.procmap['THbj'](treename, category, options)
            WtH      = self.parent.procmap['WtH'](treename, category, options)

            return (rare_top + triboson + tHbj + WtH)

    class TTBar(Process):

        latexname = 't#bar{t}, t#bar{t}#gamma'
        colour = kAzure + 8

        def base(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    inputgroup = [
                ('tops', 'ttbar_nonallhad_Pythia8'),
		('tops', 'ttgamma'),
                #('tops', 'ttbar_nonallhad'),
                #('tops', 'ttbar_dilep'),
                #('tops', 'ttbar_SingleLeptonP_MEPS_NLO'),
                #('tops', 'ttbar_SingleLeptonM_MEPS_NLO'),
                #('tops', 'ttbar_dilepton_MEPS_NLO'),
                         ]

            if self.parent.readGFW2:
                del inputgroup[:]
                inputgroup = [
                    ('tops_MMClosure', 'ttbar_nonallhad_Pythia8'),
                    ('tops_MMClosure', 'ttgamma'),
                    #('tops_MMClosure', 'ttgamma_vPP6'),
                    #('tops_MMClosure', 'ttgamma_vSherpa'),
                    #('tops_MMClosure', 'ttbar_nonallhad'),
                    #('tops_MMClosure', 'ttbar_SingleLeptonP_MEPS_NLO'),
                    #('tops_MMClosure', 'ttbar_SingleLeptonM_MEPS_NLO'),
                    #('tops_MMClosure', 'ttbar_dilepton_MEPS_NLO'),
                ]

            print("\n{0}:\n".format(self.__class__.__name__))
	    print("\n".join("{0} - {1} : {2}".format(idx,sample[0],sample[1]) for idx, sample in enumerate(inputgroup)))

            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) * self.parent.norm_factor

	    return sp

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)

            weight = None
            TTcut  = ('','TT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

            if TTcut:
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            tt_ttgamma_OLR_cut = Cut("TTBar_TTGamma_OLR","( ( mc_channel_number == 410082 && m_MEphoton_OLtty_cat1 ) || ( mc_channel_number != 410082 && m_MEphoton_OLtty_keepEvent && !m_hasMEphoton_DRgt02_nonhad ) )")
            # tt_ttgamma_OLR_cut = Cut("TTBar_TTGamma_OLR","( 1 )")
            sp = sp.subprocess(cut=tt_ttgamma_OLR_cut)

            print("\n{0} - cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp.basecut.cutnamelist, weight))

            return sp


    class TTBarNoTTBarGamma(Process):

        latexname = 't#bar{t}'
        colour = kAzure + 8

        def base(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    inputgroup = [
                ('tops', 'ttbar_nonallhad_Pythia8'),
                #('tops', 'ttbar_nonallhad'),
                #('tops', 'ttbar_dilep'),
                #('tops', 'ttbar_SingleLeptonP_MEPS_NLO'),
                #('tops', 'ttbar_SingleLeptonM_MEPS_NLO'),
                         ]

            print("\n{0}:\n".format(self.__class__.__name__))
	    print("\n".join("{0} - {1} : {2}".format(idx,sample[0],sample[1]) for idx, sample in enumerate(inputgroup)))

            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) * self.parent.norm_factor

	    return sp

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)

            weight = None
            TTcut  = ('','TT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

            if TTcut:
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\n{0} - cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp.basecut.cutnamelist, weight))

            return sp


    class TopCF(Process):

        latexname = 'tops (QMisID only)'
        colour = kAzure - 4

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

            print("\n{0}:\n".format(self.__class__.__name__))

            # Pick the subprocesses from the DB

            rare_top  = self.parent.procmap['RareTop'](treename, category, options)
            ttbar     = self.parent.procmap['TTBar'](treename, category, options)
            # ttbar     = self.parent.procmap['TTBarNoTTBarGamma'](treename, category, options)
            singletop = self.parent.procmap['SingleTop'](treename, category, options)

            basecut = category.cut

            # Plot only events where at least one lepton is charge flip.
            # Do it only for SS events: for other events, just apply a weight = 0 to kill the process

            basecut = category.cut
            weight  = None
            #if ("2Lep_SS") in basecut.cutname:
            if True:
                for CUT in basecut.cutlist:
                    if ("TRUTH") in CUT.cutname:
        	        basecut = basecut.removeCut(CUT)
            else:
                weight = '0.0'

            updatedcut = basecut & self.vardb.getCut('2Lep_TRUTH_QMisIDEvent')

            # Remember to add TT cut (if needed) as it's not in basecut!

            TTcut  = ('','TT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]
            if TTcut:
                updatedcut = updatedcut & self.vardb.getCut(TTcut)

            # Reset the cut for the subprocesses, and apply the new one

            rare_top  = rare_top.subprocess(cut=updatedcut, clearbasecut=True, eventweight=weight)
            ttbar     = ttbar.subprocess(cut=updatedcut, clearbasecut=True, eventweight=weight)
            singletop = singletop.subprocess(cut=updatedcut, clearbasecut=True, eventweight=weight)

            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("RareTop",rare_top.basecut.cutnamelist, rare_top.eventweight))
            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("TTBar",ttbar.basecut.cutnamelist, ttbar.eventweight))
            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("SingleTop",singletop.basecut.cutnamelist, singletop.eventweight))

            return rare_top + ttbar + singletop


    class Diboson(Process):

        latexname = 'WZ, ZZ, WW'
        colour = kYellow - 9

        def base(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    inputgroup = [
                ('Diboson', 'llll'),
                ('Diboson', 'lllvSFMinus'),
                ('Diboson', 'lllvOFMinus'),
                ('Diboson', 'lllvSFPlus'),
                ('Diboson', 'lllvOFPlus'),
                ('Diboson', 'llvv'),
                ('Diboson', 'llvvjj_ss_EW4'),
                ('Diboson', 'EWllssnunujj'),
                ('Diboson', 'lllvjj_EW6'),
                ('Diboson', 'lllljj_EW6'),
                ('Diboson', 'ggllll'),
                ('Diboson', 'ggllvv'),
                ('Diboson', 'WW_SHv21_improved'),
                ('Diboson', 'WZ_SHv21_improved'),
                ('Diboson', 'ZZ_SHv21_improved'),
                         ]

            print("\n{0}:\n".format(self.__class__.__name__))
	    print("\n".join("{0} - {1} : {2}".format(idx,sample[0],sample[1]) for idx, sample in enumerate(inputgroup)))

            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) * self.parent.norm_factor
            return sp

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)

            basecut = category.cut

            if all( c in category.cut.cutname for c in ["2Lep_TRUTH_ProbePromptEvent","2Lep_TRUTH_QMisIDVeto"] ):
                basecut = basecut.swapCut(self.vardb.getCut('2Lep_TRUTH_ProbePromptEvent'), self.vardb.getCut('2Lep_TRUTH_PurePromptEvent'))
                basecut = basecut.removeCut(self.vardb.getCut("2Lep_TRUTH_QMisIDVeto"))

            weight = None
            TTcut  = ('','TT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]

            sp = sp.subprocess(cut=basecut,eventweight=weight)

            if TTcut:
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\n{0} - cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp.basecut.cutnamelist, weight))

            return sp


    class DibosonCF(Process):

        latexname = 'WZ, ZZ, WW (QMisID only)'
        colour = kAzure - 9

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

            print("\n{0}:\n".format(self.__class__.__name__))

            # Pick the subprocesses from the DB

            diboson  = self.parent.procmap['Diboson'](treename, category, options)

            basecut = category.cut

            # Plot only events where at least one lepton is charge flip.
            # Do it only for SS events: for other events, just apply a weight = 0 to kill the process

            basecut = category.cut
            weight  = None
            #if ("2Lep_SS") in basecut.cutname:
            if True:
                for CUT in basecut.cutlist:
                    if ("TRUTH") in CUT.cutname:
        	        basecut = basecut.removeCut(CUT)
            else:
                weight = '0.0'

            updatedcut = basecut & self.vardb.getCut('2Lep_TRUTH_QMisIDEvent')

            # Remember to add TT cut (if needed) as it's not in basecut!

            TTcut  = ('','TT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]
            if TTcut:
                updatedcut = updatedcut & self.vardb.getCut(TTcut)

            # Reset the cut for the subprocesses, and apply the new one

            diboson = diboson.subprocess(cut=updatedcut, clearbasecut=True, eventweight=weight)

            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("Diboson",diboson.basecut.cutnamelist, diboson.eventweight))

            return diboson


    class Triboson(Process):

        latexname = 'VVV'
        colour = kYellow - 9

        def base(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    inputgroup = [
                ('Triboson', '*'),
                         ]

            print("\n{0}:\n".format(self.__class__.__name__))
	    print("\n".join("{0} - {1} : {2}".format(idx,sample[0],sample[1]) for idx, sample in enumerate(inputgroup)))

            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) * self.parent.norm_factor
            return sp

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)

            basecut = category.cut

            if all( c in category.cut.cutname for c in ["2Lep_TRUTH_ProbePromptEvent","2Lep_TRUTH_QMisIDVeto"] ):
                basecut = basecut.swapCut(self.vardb.getCut('2Lep_TRUTH_ProbePromptEvent'), self.vardb.getCut('2Lep_TRUTH_PurePromptEvent'))
                basecut = basecut.removeCut(self.vardb.getCut("2Lep_TRUTH_QMisIDVeto"))

            weight = None
            TTcut  = ('','TT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]

            sp = sp.subprocess(cut=basecut,eventweight=weight)

            if TTcut:
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\n{0} - cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp.basecut.cutnamelist, weight))

            return sp


    class Prompt(Process):

        latexname = 'Prompt'
        colour = kYellow - 9

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

            print("\n{0}:\n".format(self.__class__.__name__))

            # Pick the subprocesses from the DB

            diboson   = self.parent.procmap['Diboson'](treename, category, options)
            triboson  = self.parent.procmap['Triboson'](treename, category, options)
            rare_top  = self.parent.procmap['RareTop'](treename, category, options)
            ttbarw    = self.parent.procmap['TTBarW'](treename, category, options)
            ttbarz    = self.parent.procmap['TTBarZ'](treename, category, options)
            tHbj      = self.parent.procmap['THbj'](treename, category, options)
            WtH       = self.parent.procmap['WtH'](treename, category, options)
            #
            ttbar     = self.parent.procmap['TTBar'](treename, category, options)
            singletop = self.parent.procmap['SingleTop'](treename, category, options)
            zjets     = self.parent.procmap['Zjets'](treename, category, options)
            wjets     = self.parent.procmap['Wjets'](treename, category, options)

            # Clean up from any truth cut

	    basecut = category.cut

	    for CUT in basecut.cutlist:
	    	if ("TRUTH") in CUT.cutname:
	    	    basecut = basecut.removeCut(CUT)

            # Require all lepton be prompt, and veto on QMisID

            updatedcut = basecut & self.vardb.getCut('2Lep_TRUTH_PurePromptEvent')

            # Remember to add TT cut (if needed) as it's not in basecut!

            TTcut  = ('','TT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]
            if TTcut:
                updatedcut = updatedcut & self.vardb.getCut(TTcut)

            # Reset the cut for the subprocesses, and apply the new one

            diboson   = diboson.subprocess(cut=updatedcut, clearbasecut=True)
            triboson  = triboson.subprocess(cut=updatedcut, clearbasecut=True)
            rare_top  = rare_top.subprocess(cut=updatedcut, clearbasecut=True)
            ttbarw    = ttbarw.subprocess(cut=updatedcut, clearbasecut=True)
            ttbarz    = ttbarz.subprocess(cut=updatedcut, clearbasecut=True)
            tHbj      = tHbj.subprocess(cut=updatedcut, clearbasecut=True)
            WtH       = WtH.subprocess(cut=updatedcut, clearbasecut=True)
            ttbar     = ttbar.subprocess(cut=updatedcut, clearbasecut=True)
            singletop = singletop.subprocess(cut=updatedcut, clearbasecut=True)
            zjets     = zjets.subprocess(cut=updatedcut, clearbasecut=True)
            wjets     = wjets.subprocess(cut=updatedcut, clearbasecut=True)

            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("Diboson",diboson.basecut.cutnamelist, diboson.eventweight))
            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("Triboson",triboson.basecut.cutnamelist, triboson.eventweight))
            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("RareTop",rare_top.basecut.cutnamelist, rare_top.eventweight))
            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("THbj",tHbj.basecut.cutnamelist, tHbj.eventweight))
            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("WtH",WtH.basecut.cutnamelist, WtH.eventweight))
            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("TTBarW",ttbarw.basecut.cutnamelist, ttbarw.eventweight))
            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("TTBarZ",ttbarz.basecut.cutnamelist, ttbarz.eventweight))
            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("TTBar",ttbar.basecut.cutnamelist, ttbar.eventweight))
            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("SingleTop",singletop.basecut.cutnamelist, singletop.eventweight))
            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("Zjets",zjets.basecut.cutnamelist, zjets.eventweight))
            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("Wjets",wjets.basecut.cutnamelist, wjets.eventweight))

            return diboson + rare_top + ttbarw + ttbarz + ttbar + singletop + zjets + wjets + triboson + tHbj + WtH


    class OtherPrompt(Process):

        #latexname = 'Prompt (non t#bar{t} W)'
        latexname = 'Other Prompt'
        colour = kYellow - 10

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

            print("\n{0}:\n".format(self.__class__.__name__))

            # Pick the subprocesses from the DB

            diboson   = self.parent.procmap['Diboson'](treename, category, options)
            triboson  = self.parent.procmap['Triboson'](treename, category, options)
            rare_top  = self.parent.procmap['RareTop'](treename, category, options)
            ttbarz    = self.parent.procmap['TTBarZ'](treename, category, options)
            tHbj      = self.parent.procmap['THbj'](treename, category, options)
            WtH       = self.parent.procmap['WtH'](treename, category, options)
            #
            ttbar     = self.parent.procmap['TTBar'](treename, category, options)
            singletop = self.parent.procmap['SingleTop'](treename, category, options)
            zjets     = self.parent.procmap['Zjets'](treename, category, options)
            wjets     = self.parent.procmap['Wjets'](treename, category, options)

            # Clean up from any truth cut

	    basecut = category.cut

	    for CUT in basecut.cutlist:
	    	if ("TRUTH") in CUT.cutname:
	    	    basecut = basecut.removeCut(CUT)

            # Require all lepton be prompt, and veto on QMisID

            updatedcut = basecut & self.vardb.getCut('2Lep_TRUTH_PurePromptEvent')

            # Remember to add TT cut (if needed) as it's not in basecut!

            TTcut  = ('','TT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]
            if TTcut:
                updatedcut = updatedcut & self.vardb.getCut(TTcut)

            # Reset the cut for the subprocesses, and apply the new one

            diboson   = diboson.subprocess(cut=updatedcut, clearbasecut=True)
            triboson  = triboson.subprocess(cut=updatedcut, clearbasecut=True)
            rare_top  = rare_top.subprocess(cut=updatedcut, clearbasecut=True)
            ttbarz    = ttbarz.subprocess(cut=updatedcut, clearbasecut=True)
            tHbj      = tHbj.subprocess(cut=updatedcut, clearbasecut=True)
            WtH       = WtH.subprocess(cut=updatedcut, clearbasecut=True)
            ttbar     = ttbar.subprocess(cut=updatedcut, clearbasecut=True)
            singletop = singletop.subprocess(cut=updatedcut, clearbasecut=True)
            zjets     = zjets.subprocess(cut=updatedcut, clearbasecut=True)
            wjets     = wjets.subprocess(cut=updatedcut, clearbasecut=True)

            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("Diboson",diboson.basecut.cutnamelist, diboson.eventweight))
            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("Triboson",triboson.basecut.cutnamelist, triboson.eventweight))
            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("RareTop",rare_top.basecut.cutnamelist, rare_top.eventweight))
            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("THbj",tHbj.basecut.cutnamelist, tHbj.eventweight))
            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("WtH",WtH.basecut.cutnamelist, WtH.eventweight))
            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("TTBarZ",ttbarz.basecut.cutnamelist, ttbarz.eventweight))
            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("TTBar",ttbar.basecut.cutnamelist, ttbar.eventweight))
            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("SingleTop",singletop.basecut.cutnamelist, singletop.eventweight))
            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("Zjets",zjets.basecut.cutnamelist, zjets.eventweight))
            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("Wjets",wjets.basecut.cutnamelist, wjets.eventweight))

            return diboson + rare_top + ttbarz + ttbar + singletop + zjets + wjets + triboson + tHbj + WtH


    class QMisIDMC(Process):

        latexname = 'QMisID (MC)'
        colour = kAzure - 4

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

            print("\n{0}:\n".format(self.__class__.__name__))

            dibosoncf = self.parent.procmap['DibosonCF'](treename, category, options)
            topcf     = self.parent.procmap['TopCF'](treename, category, options)
            zjetscf   = self.parent.procmap['ZjetsCF'](treename, category, options)

            weight = 1.0

	    return ( dibosoncf.subprocess(eventweight=weight) + topcf.subprocess(eventweight=weight) + zjetscf.subprocess(eventweight=weight))


    class QMisID(Process):

        latexname = 'QMisID'
        colour = kAzure - 4

        def base(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    inputgroup = [
                    ('Data', 'physics_Main'),
                         ]

            if self.parent.readGFW2:
                del inputgroup[:]
                inputgroup = [ ('Data', 'fakes_mm') ]

            print("\n{0}:\n".format(self.__class__.__name__))
	    print("\n".join("{0} - {1} : {2}".format(idx,sample[0],sample[1]) for idx, sample in enumerate(inputgroup)))

            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees)
            return sp

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)

            TTcut  = ''
            weight = '0.0' # NB: do not change this!

            # Categories w/ any of these cuts must get the QMisID weight...

            cut_w_el  = ['2Lep_ElEl_Event','2Lep_MuEl_Event','2Lep_ElMu_Event','2Lep_OF_Event']

            # ... and these must not!

            cut_wo_el = ['2Lep_MuMu_Event']

            # Take the category cut (NB: cache it in a temporary cut, since we do NOT want to modify the original one!)

            basecut = category.cut

            print("\nQMisID initial cut: {0}".format(basecut.cutnamelist))

            if ( self.parent.channel=='2LSS' ) or ( self.parent.channel=='3L' ):
                TTcut = 'TT'

            # Clean up from any truth cut

	    for CUT in basecut.cutlist:
	    	if ("TRUTH") in CUT.cutname:
	    	    basecut = basecut.removeCut(CUT)

            if TTcut != '':
                basecut = basecut & self.vardb.getCut(TTcut)

            # If SS cut in category, need to switch to OS before applying QMisID weight.
            # Otherwise, just forget about QMisID, at all (NB: the weight is 0 by default, so the process will be suppressed).

            if ("2Lep_SS") in category.cut.cutname:
                basecut = basecut.swapCut(self.vardb.getCut('2Lep_SS'), -self.vardb.getCut('2Lep_SS'))
            else:
                sp = sp.subprocess(cut=basecut,eventweight=weight)
                print("\n{0} - cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp.basecut.cutnamelist, weight))
                return sp

            if bool([ cut for cut in cut_w_el if cut in category.cut.cutname ]):

                weight = 'QMisIDWeight' if not self.parent.readGFW2 else 'QMisID_MM_EventWeight'

                sp = sp.subprocess(cut=basecut,eventweight=weight)
                print("\n{0} - cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp.basecut.cutnamelist, weight))

            elif bool([ cut for cut in cut_wo_el if cut in category.cut.cutname ]):

                weight = '0.0'
                sp = sp.subprocess(cut=basecut,eventweight=weight)
                print("\n{0} - cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp.basecut.cutnamelist, weight))

            else:

                weight_SF = weight_OF = '0.0'
                sp_SF = sp_OF = None

                if ( 'FakeCREl' ) in category.name:

                    cut_SF = basecut & self.vardb.getCut('2Lep_ElEl_Event')
                    cut_OF = basecut & self.vardb.getCut('2Lep_OF_Event')

                    weight_SF = weight_OF = 'QMisIDWeight' if not self.parent.readGFW2 else 'QMisID_MM_EventWeight'

                    sp_SF = sp.subprocess(cut=cut_SF,eventweight=weight_SF)
                    sp_OF = sp.subprocess(cut=cut_OF,eventweight=weight_OF)

                    print("\nQMisID sp_SF - cuts: {0}, process weight: {1}".format(sp_SF.basecut.cutnamelist, weight_SF))
                    print("\nQMisID sp_OF - cuts: {0}, process weight: {1}".format(sp_OF.basecut.cutnamelist, weight_OF))

                elif ( 'FakeCRMu' ) in category.name:

                    cut_SF = basecut & self.vardb.getCut('2Lep_MuMu_Event')
                    cut_OF = basecut & self.vardb.getCut('2Lep_OF_Event')

                    weight_SF = '0.0' # mu-mu region MUST get zero QMisID weight!
                    weight_OF = 'QMisIDWeight' if not self.parent.readGFW2 else 'QMisID_MM_EventWeight'

                    sp_SF = sp.subprocess(cut=cut_SF,eventweight=weight_SF)
                    sp_OF = sp.subprocess(cut=cut_OF,eventweight=weight_OF)

                    print("\nQMisID sp_SF - cuts: {0}, process weight: {1}".format(sp_SF.basecut.cutnamelist, weight_SF))
                    print("\nQMisID sp_OF - cuts: {0}, process weight: {1}".format(sp_OF.basecut.cutnamelist, weight_OF))

                else:

                    cut_elel = basecut & self.vardb.getCut('2Lep_ElEl_Event')
                    cut_mumu = basecut & self.vardb.getCut('2Lep_MuMu_Event')
                    cut_OF   = basecut & self.vardb.getCut('2Lep_OF_Event')

                    weight_mumu = '0.0' # mu-mu region MUST get zero QMisID weight!
                    weight_elel = weight_OF = 'QMisIDWeight' if not self.parent.readGFW2 else 'QMisID_MM_EventWeight'

                    sp_elel = sp.subprocess(cut=cut_elel,eventweight=weight_elel)
                    sp_mumu = sp.subprocess(cut=cut_mumu,eventweight=weight_mumu)

                    sp_SF   = sp_elel + sp_mumu
                    sp_OF   = sp.subprocess(cut=cut_OF,eventweight=weight_OF)

                    print("\nQMisID sp_elel - cuts: {0}, process weight: {1}".format(sp_elel.basecut.cutnamelist, weight_elel))
                    print("\nQMisID sp_mumu - cuts: {0}, process weight: {1}".format(sp_mumu.basecut.cutnamelist, weight_mumu))
                    print("\nQMisID sp_OF - cuts: {0}, process weight: {1}".format(sp_OF.basecut.cutnamelist, weight_OF))

                sp = sp_SF + sp_OF

            return sp

    class ZpeakSidebandBkg(Process):

        latexname = 'Background (from Z sidebands)'
        colour = kBlue

        def base(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

            inputgroup = [
                    ('Data', 'physics_Main'),
                         ]

            print("\n{0}:\n".format(self.__class__.__name__))
	    print("\n".join("{0} - {1} : {2}".format(idx,sample[0],sample[1]) for idx, sample in enumerate(inputgroup)))

            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees)
            return sp

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)

            # WORK IN PROGRESS...

            return sp


    class FakesMC(Process):

        latexname = 'Fakes (MC)'
        colour = kCyan -9

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

            print("\n{0}:\n".format(self.__class__.__name__))

            # Pick the subprocesses from the DB

            diboson   = self.parent.procmap['Diboson'](treename, category, options)
            rare_top  = self.parent.procmap['RareTop'](treename, category, options)
            ttbarw    = self.parent.procmap['TTBarW'](treename, category, options)
            ttbarz    = self.parent.procmap['TTBarZ'](treename, category, options)
            # ttbargamma = self.parent.procmap['TTBarGamma'](treename, category, options)
            ttbar     = self.parent.procmap['TTBar'](treename, category, options)
            singletop = self.parent.procmap['SingleTop'](treename, category, options)
            zjets     = self.parent.procmap['Zjets'](treename, category, options)
            wjets     = self.parent.procmap['Wjets'](treename, category, options)

            # Clean up from any truth cut

	    basecut = category.cut

	    for CUT in basecut.cutlist:
	    	if ("TRUTH") in CUT.cutname:
	    	    basecut = basecut.removeCut(CUT)

            # Require at least one !prompt lepton, and veto on QMisID

            updatedcut = basecut & self.vardb.getCut('2Lep_TRUTH_NonPromptEvent')

            # Remember to add TT cut (if needed) as it's not in basecut!

            TTcut  = ('','TT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]
            if TTcut:
                updatedcut = updatedcut & self.vardb.getCut(TTcut)

            # Reset the cut to the subprocesses

            diboson   = diboson.subprocess(cut=updatedcut, clearbasecut=True)
            rare_top  = rare_top.subprocess(cut=updatedcut, clearbasecut=True)
            ttbarw    = ttbarw.subprocess(cut=updatedcut, clearbasecut=True)
            ttbarz    = ttbarz.subprocess(cut=updatedcut, clearbasecut=True)
            ttbar     = ttbar.subprocess(cut=updatedcut, clearbasecut=True)
            singletop = singletop.subprocess(cut=updatedcut, clearbasecut=True)
            zjets     = zjets.subprocess(cut=updatedcut, clearbasecut=True)
            wjets     = wjets.subprocess(cut=updatedcut, clearbasecut=True)

            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("Diboson",diboson.basecut.cutnamelist, diboson.eventweight))
            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("RareTop",rare_top.basecut.cutnamelist, rare_top.eventweight))
            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("TTBarW",ttbarw.basecut.cutnamelist, ttbarw.eventweight))
            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("TTBarZ",ttbarz.basecut.cutnamelist, ttbarz.eventweight))
            # print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("TTBarGamma",ttbargamma.basecut.cutnamelist, ttbarw.eventweight))
            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("TTBar",ttbar.basecut.cutnamelist, ttbar.eventweight))
            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("SingleTop",singletop.basecut.cutnamelist, singletop.eventweight))
            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("Zjets",zjets.basecut.cutnamelist, zjets.eventweight))
            print("\n{0} - UPDATED cuts: {1}, process weight: {2}".format("Wjets",wjets.basecut.cutnamelist, wjets.eventweight))

            return diboson + rare_top + ttbarw + ttbarz + ttbar + singletop + zjets + wjets # + ttbargamma


    class FakesFF(Process):

        latexname = 'Fakes FF'
        colour = kAzure - 9

        def base(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    inputgroup = [
                    ('Data', 'physics_Main'),
                         ]

            print("\n{0}:\n".format(self.__class__.__name__))
	    print("\n".join("{0} - {1} : {2}".format(idx,sample[0],sample[1]) for idx, sample in enumerate(inputgroup)))

            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees)
            return sp

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)

            TLcut=''
            LTcut=''
            LLcut=''
            weight = None

            if self.parent.channel=='2LSS' or self.parent.channel=='3L':
                TLcut='TL'
                LTcut='LT'
                LLcut='LL'
                weight='FFWeight'

            print 'FakesFF - weight name : ', weight, ', TLcut : ', TLcut, ', LTcut : ', LTcut, ', LLcut : ', LLcut

            # Now doing prompt background subtractions from fakes. sublist must contain both prompt and QMisIDs

            sublist = [ item for item in self.parent.sub_backgrounds ]
            for sample in sublist:
                sp = sp - self.parent.procmap[sample].base(treename, category, options) # Here it is important to have used base otherwise at the sub sample
                                                                                        # there would be already applied the selection specified in the category __call__ of that sample
                                                                                        # i.e. for example also the iso-iso cut which is orthogonal to the cuts done in this sample which are fake iso or fake fake.
            #print sp
            sp_TL = sp.subprocess(cut=self.vardb.getCut(TLcut))
            sp_LT = sp.subprocess(cut=self.vardb.getCut(LTcut))
            sp_LL = sp.subprocess(cut=self.vardb.getCut(LLcut))

            #IN THIS VERSION OF THE CODE THE SIGN MINUS FOR FF EVENTS IS ALREADY IN THE WEIGHT SO THE SUBRACTION OF THE DOUBLE FAKE COUNTING IS DONE ADDING SP_FF
            sp = sp_TL + sp_LT +sp_LL
            #OLD WAY subtracting the double counting
            #sp = sp_if - sp_ff
            sp = sp.subprocess(cut=category.cut,eventweight=weight)
            #print sp
            return sp


    class FakesMM(Process):

        latexname = 'Fakes MM'
        colour = kTeal - 9

        def base(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    inputgroup = [
                    ('Data', 'physics_Main'),
                         ]

            if self.parent.readGFW2:
                del inputgroup[:]
                inputgroup = [ ('Data', 'fakes_mm') ]

            print("\n{0}:\n".format(self.__class__.__name__))
	    print("\n".join("{0} - {1} : {2}".format(idx,sample[0],sample[1]) for idx, sample in enumerate(inputgroup)))

            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees)
            return sp

        def getFakesMM(self, treename='physics', category=None, options={}, sp=None, inputcut=None):

            if self.parent.readGFW2:
                treename = 'nominal'

            debugflag = any( proc in self.parent.debugprocs for proc in [self.__class__.__name__,"ALL"])

            TTcut  = ('','FakesSideband_TT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]
            TLcut  = ('','FakesSideband_TL')[any( c == self.parent.channel for c in ['2LSS','3L'] )]
            LTcut  = ('','FakesSideband_LT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]
            LLcut  = ('','FakesSideband_LL')[any( c == self.parent.channel for c in ['2LSS','3L'] )]
            weight = (None,'MMWeight')[any( c == self.parent.channel for c in ['2LSS','3L'] )]

            if weight and self.parent.readGFW2:
                weight = "MM_EventWeight"

            # Clean up from any truth cut

	    basecut = inputcut

	    for CUT in basecut.cutlist:
	    	if ("TRUTH") in CUT.cutname:
	    	    basecut = basecut.removeCut(CUT)

            sp_TT  = sp.subprocess(cut=basecut & self.vardb.getCut(TTcut), eventweight=weight)
            sp_TL  = sp.subprocess(cut=basecut & self.vardb.getCut(TLcut), eventweight=weight)
            sp_LT  = sp.subprocess(cut=basecut & self.vardb.getCut(LTcut), eventweight=weight)
            sp_LL  = sp.subprocess(cut=basecut & self.vardb.getCut(LLcut), eventweight=weight)

            if debugflag:
                print(" ")
                print("{0} - TT cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp_TT.basecut.cutnamelist, weight))
                print("{0} - TL cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp_TL.basecut.cutnamelist, weight))
                print("{0} - LT cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp_LT.basecut.cutnamelist, weight))
                print("{0} - LL cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp_LL.basecut.cutnamelist, weight))

            # Perform QMisID per-event subtraction taking advantage of the linearity of MM equations
            # NB: as it is now, this works only for DD QMisID! If we wanted to use MC QMiSID, we would have to
            #     rewweight w/ the DD MM weight also the TTBar ntuple, and then subtract the weighted TTBar contribution
            #     (w/ the truth QMisID cut applied to TTBar)

	    sublist = [ item for item in self.parent.sub_backgrounds if item == "QMisID" ]
            sp_OS_TT = sp_OS_TL = sp_OS_LT = sp_OS_LL = None
            sp_QMisID_TT = sp_QMisID_TL = sp_QMisID_LT = sp_QMisID_LL = None

            for process in sublist:

                if ("2Lep_MuMu_Event") in inputcut.cutname:
                    print("NO QMisID subtraction for MuMu events!!")
                    continue

	        QMisIDcut = basecut.swapCut(self.vardb.getCut('2Lep_SS'), -self.vardb.getCut('2Lep_SS'))

	        OSweight  = 'QMisIDWeight * MMWeight' if not self.parent.readGFW2 else 'QMisID_MM_EventWeight * MM_EventWeight'
                SSweight  = 'QMisIDWeight' if not self.parent.readGFW2 else 'QMisID_MM_EventWeight'

                # These are the actual QMisID events in the sidebands

                sp_QMisID_TT = (self.parent.procmap[process].base(treename,category,options)).subprocess(cut=QMisIDcut & self.vardb.getCut(TTcut), eventweight=SSweight)
                sp_QMisID_TL = (self.parent.procmap[process].base(treename,category,options)).subprocess(cut=QMisIDcut & self.vardb.getCut(TLcut), eventweight=SSweight)
                sp_QMisID_LT = (self.parent.procmap[process].base(treename,category,options)).subprocess(cut=QMisIDcut & self.vardb.getCut(LTcut), eventweight=SSweight)
                sp_QMisID_LL = (self.parent.procmap[process].base(treename,category,options)).subprocess(cut=QMisIDcut & self.vardb.getCut(LLcut), eventweight=SSweight)

                # These are the OS weighted events that need to be subtracted

                sp_OS_TT = (self.parent.procmap[process].base(treename,category,options)).subprocess(cut=QMisIDcut & self.vardb.getCut(TTcut), eventweight=OSweight)
                sp_OS_TL = (self.parent.procmap[process].base(treename,category,options)).subprocess(cut=QMisIDcut & self.vardb.getCut(TLcut), eventweight=OSweight)
                sp_OS_LT = (self.parent.procmap[process].base(treename,category,options)).subprocess(cut=QMisIDcut & self.vardb.getCut(LTcut), eventweight=OSweight)
                sp_OS_LL = (self.parent.procmap[process].base(treename,category,options)).subprocess(cut=QMisIDcut & self.vardb.getCut(LLcut), eventweight=OSweight)

                if debugflag:
                    print(" ")
                    print("Subtracting OS TT sp: {0}, weight: {1}".format(sp_OS_TT.basecut.cutnamelist, OSweight))
                    print("{0} - TT :\nFakes (before OS subtraction) = {1:.2f},\nQMisID = {2:.2f},\nOS = {3:.2f}".format(self.__class__.__name__,(sp_TT.numberstats())[0],(sp_QMisID_TT.numberstats())[0],(sp_OS_TT.numberstats())[0]))
                    print("Subtracting OS TL sp: {0}, weight: {1}".format(sp_OS_TL.basecut.cutnamelist, OSweight))
                    print("{0} - TL :\nFakes (before OS subtraction) = {1:.2f},\nQMisID = {2:.2f},\nOS = {3:.2f}".format(self.__class__.__name__,(sp_TL.numberstats())[0],(sp_QMisID_TL.numberstats())[0],(sp_OS_TL.numberstats())[0]))
                    print("Subtracting OS LT sp: {0}, weight: {1}".format(sp_OS_LT.basecut.cutnamelist, OSweight))
                    print("{0} - LT :\nFakes (before OS subtraction) = {1:.2f},\nQMisID = {2:.2f},\nOS = {3:.2f}".format(self.__class__.__name__,(sp_LT.numberstats())[0],(sp_QMisID_LT.numberstats())[0],(sp_OS_LT.numberstats())[0]))
                    print("Subtracting OS LL sp: {0}, weight: {1}".format(sp_OS_LL.basecut.cutnamelist, OSweight))
                    print("{0} - LL :\nFakes (before OS subtraction) = {1:.2f},\nQMisID = {2:.2f},\nOS = {3:.2f}".format(self.__class__.__name__,(sp_LL.numberstats())[0],(sp_QMisID_LL.numberstats())[0],(sp_OS_LL.numberstats())[0]))

                sp_TT  = sp_TT - sp_OS_TT
                sp_TL  = sp_TL - sp_OS_TL
                sp_LT  = sp_LT - sp_OS_LT
                sp_LL  = sp_LL - sp_OS_LL

                if debugflag:
                    print(" ")
                    print("{0} - TT : Fakes (AFTER OS subtraction) = {1:.2f}".format(self.__class__.__name__,(sp_TT.numberstats())[0]))
                    print("{0} - TL : Fakes (AFTER OS subtraction) = {1:.2f}".format(self.__class__.__name__,(sp_TL.numberstats())[0]))
                    print("{0} - LT : Fakes (AFTER OS subtraction) = {1:.2f}".format(self.__class__.__name__,(sp_LT.numberstats())[0]))
                    print("{0} - LL : Fakes (AFTER OS subtraction) = {1:.2f}".format(self.__class__.__name__,(sp_LL.numberstats())[0]))

            sp = sp_TT + sp_TL + sp_LT + sp_LL

            # Option A: do nothing
            # Option B: rescale fakes by the non-closure SF in each flavour category (do nothing for mm)
            # Option C: use a different SF for ee (since there's already a rescaling in the electron fake rate)

            optA = False
            optB = False
            optC = True

            SF = 1.0
            if not optA:
                if "2Lep_ElEl_Event" in inputcut.cutnamelist:
                    if optB:
                        SF = 1.428 if ("2Lep_NJet_SR" in inputcut.cutnamelist) else 1.417
                    elif optC:
                        SF = 1.130 if ("2Lep_NJet_SR" in inputcut.cutnamelist) else 1.122
                elif "2Lep_OF_Event" in inputcut.cutnamelist:
                    SF = 1.222 if ("2Lep_NJet_SR" in inputcut.cutnamelist) else 1.093
                print("\nApply SF: {0}\n".format(SF))
            sp = sp * SF

            if debugflag:
                print(" ")
                print("{0} ---> Total Fakes = {1:.2f} +- {2:.2f}".format(self.__class__.__name__,sp.numberstats()[0],sp.numberstats()[1]))

	    return sp

        def __call__(self, treename='physics', category=None, options={}):

            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction

            sp = self.base(treename, category, options)

            if any( cutname in category.cut.cutnamelist for cutname in ["2Lep_MuMu_Event","2Lep_ElEL_Event","2Lep_OF_Event"] ):
                sp = self.getFakesMM(treename, category, options, sp, category.cut)
            else:
                # If the input cut is inclusive in ee,mm,OF, we need to split the 3 subprocesses, otherwise the QMisID subtraction
                # will screw up results for mm events...

                cut_ee = category.cut & self.vardb.getCut('2Lep_ElEl_Event')
                cut_mm = category.cut & self.vardb.getCut('2Lep_MuMu_Event')
                cut_OF = category.cut & self.vardb.getCut('2Lep_OF_Event')

                sp_ee = self.getFakesMM(treename, category, options, sp, cut_ee)
                sp_mm = self.getFakesMM(treename, category, options, sp, cut_mm)
                sp_OF = self.getFakesMM(treename, category, options, sp, cut_OF)

                sp = sp_ee + sp_mm + sp_OF

            return sp


    class FakesTHETA(Process):

        # Consider the following regions:
        #
        # Region A: TT, njet >= 5 (SR)
        # Region B (in DATA): LT,TL njet >= 5
        # Region C (in DATA): TT, njet = [2,3,4]
        # Region D (in DATA): LT,TL  njet = [2,3,4]
        #
        # After subtracting prompt MC and charge flips from data in B, C and D, the estimate of fakes in A is:
        #
        # (C/D) * B = A
        #
        # where theta_e = (C/D)(ee), theta_mu = (C/D)(mumu)
        #
        # This corrects the normalisation (we assume the shape in B and A is the same!) for the fake probability as measured in data.
	#
	# This method is eventually just an "inclusive" fake factor method
        #
        # Hence we return B (ee, mumu, OF), scaled by the ratio of TT/LT events (ee, mumu) in low jet multiplicity region.
        # The OF case is obtained via:
        #
        # theta_e * B(L_e,T_mu) + theta_mu * B (L_mu,T_e) = A (OF)

        latexname = 'Fakes #theta method'
        colour = kCyan - 9

    	theta = {
    	    'El': (999.0, 0.0),
    	    'Mu': (999.0, 0.0),
    	}

        def base(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    inputgroup = [
                    ('Data', 'physics_Main'),
                ]

            print("\n{0}:\n".format(self.__class__.__name__))
	    print("\n".join("{0} - {1} : {2}".format(idx,sample[0],sample[1]) for idx, sample in enumerate(inputgroup)))

            trees = self.inputs.getTrees(treename, inputgroup)
            return self.subprocess(trees=trees)

        def clearTheta(self):
            self.theta['El'] = (999.0, 0.0)
            self.theta['Mu'] = (999.0, 0.0)

        def calcTheta(self, sp_num, sp_denom, stream=None, options={}):

            if not sp_denom.number():
                print ("ERROR: Cannot calculate theta transfer factor! Denominator = 0")

            theta = (sp_num/sp_denom).numberstats()

            print ("\n***********************************************************************\n")
            print ("Calculated theta transfer factor for stream {0}: theta = {1:.3f} +- {2:.3f}".format(stream,theta[0],theta[1]))
            print ("\n***********************************************************************\n")

            return theta

        def applyThetaFactor(self, sp, theta, options={}):
            systematics = options.get('systematics', None)
            systematicsdirection = options.get('systematicsdirection', None)

            if systematics:
                if systematicsdirection == 'UP':
                    systdir = 1.0
                elif systematicsdirection == 'DOWN':
                    systdir = -1.0
                theta = (theta[0] + theta[1]*systdir, theta[1])

            ssp = sp.subprocess()
            ssp *= theta[0]
            return ssp

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

            print("\n{0}\n".format(self.__class__.__name__))

            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction

            TTCut =''
            TLCut =''
            LTCut =''
            TelLmuCut =''
            LelTmuCut =''
            TmuLelCut =''
            LmuTelCut =''
            weight   = '1.0'
            weightMC = (category.weight,"1.0")[bool(TTHBackgrounds.noWeights)]

            if ( self.parent.channel=='2LSS' ):
                TTCut  = self.vardb.getCut('FakesSideband_TT')
                TLCut  = self.vardb.getCut('FakesSideband_TL')
                LTCut  = self.vardb.getCut('FakesSideband_LT')
                TelLmuCut = self.vardb.getCut('FakesSideband_TelLmu')
                LelTmuCut = self.vardb.getCut('FakesSideband_LelTmu')
                TmuLelCut = self.vardb.getCut('FakesSideband_TmuLel')
                LmuTelCut = self.vardb.getCut('FakesSideband_LmuTel')

            TL_LT_Cut = (TLCut | LTCut)

            # Take the base suprocess (DATA)

            sp = self.base(treename, category, options)

            # Cache the category cut in a base cut which can be modified at no risk

            basecut = category.cut

            print("base sp: {0}".format(basecut.cutnamelist))

            # Remove the cuts defining the flavour composition (this is made just for calculating thetas...)

            basecut = basecut.removeCut(self.vardb.getCut('2Lep_ElEl_Event'))
            basecut = basecut.removeCut(self.vardb.getCut('2Lep_MuMu_Event'))
            basecut = basecut.removeCut(self.vardb.getCut('2Lep_OF_Event'))

            # Remove the cuts defining the jet multiplicity

            basecut = basecut.removeCut(self.vardb.getCut('2Lep_NJet_SR'))
            basecut = basecut.removeCut(self.vardb.getCut('2Lep_NJet_CR'))

            self.clearTheta()

            if self.theta['El'][0] == 999.0 :

                print ("Calculating theta_el from data in regions C,D...")

                # Define selection for region C and D (ee)

                cut_sp_C_el = basecut & self.vardb.getCut('2Lep_NJet_CR') & self.vardb.getCut('2Lep_ElEl_Event') & TTCut & self.vardb.getCut('2Lep_Zsidescut') & self.vardb.getCut('2Lep_ElEtaCut')
                cut_sp_D_el = basecut & self.vardb.getCut('2Lep_NJet_CR') & self.vardb.getCut('2Lep_ElEl_Event') & TL_LT_Cut & self.vardb.getCut('2Lep_Zsidescut') & self.vardb.getCut('2Lep_ElEtaCut')

                sp_C_el = sp.subprocess( cut = cut_sp_C_el, eventweight=weight )
                sp_D_el = sp.subprocess( cut = cut_sp_D_el, eventweight=weight )

                print("Region C (el) sp: {0}".format(cut_sp_C_el.cutnamelist))
                print("Region D (el) sp: {0}".format(cut_sp_D_el.cutnamelist))

                # Get a list of stuff to subtract for region C and D (i.e, prompt MC and charge flips)

                sublist = [ item for item in self.parent.sub_backgrounds ]

                # ... and now subtract!

                for sample in sublist:

                    print ("Subtracting {0} from data in regions C,D...".format(sample))

                    this_cut_sp_C_el = cut_sp_C_el
                    this_cut_sp_D_el = cut_sp_D_el

                    this_weight = weightMC

                    # NB: here it is crucial to call .base() on the subprocess, otherwise the subprocess would have the cuts
                    # defined in its own __call__ method already applied, whcih in general is not what we want
                    # (e.g., it might have a TT selection applied, when we want to consider TL events instead...)

                    sub_sample_C_el = self.parent.procmap[sample].base(treename,category,options)
                    sub_sample_D_el = self.parent.procmap[sample].base(treename,category,options)

                    if sample == "QMisIDMC":
                        this_cut_sp_C_el = this_cut_sp_C_el.swapCut(self.vardb.getCut('2Lep_TRUTH_PurePromptEvent'),self.vardb.getCut('2Lep_TRUTH_QMisIDEvent'))
                        this_cut_sp_D_el = this_cut_sp_D_el.swapCut(self.vardb.getCut('2Lep_TRUTH_PurePromptEvent'),self.vardb.getCut('2Lep_TRUTH_QMisIDEvent'))
                    if sample == "QMisID":
                        this_cut_sp_C_el = this_cut_sp_C_el.removeCut(self.vardb.getCut('2Lep_TRUTH_PurePromptEvent'))
                        this_cut_sp_C_el = this_cut_sp_C_el.swapCut(self.vardb.getCut('2Lep_SS'), -self.vardb.getCut('2Lep_SS'))
                        this_cut_sp_D_el = this_cut_sp_D_el.removeCut(self.vardb.getCut('2Lep_TRUTH_PurePromptEvent'))
                        this_cut_sp_D_el = this_cut_sp_D_el.swapCut(self.vardb.getCut('2Lep_SS'), -self.vardb.getCut('2Lep_SS'))
                        this_weight = 'QMisIDWeight'

                    sub_sample_C_el = sub_sample_C_el.subprocess( cut = this_cut_sp_C_el, eventweight=this_weight )
                    sub_sample_D_el = sub_sample_D_el.subprocess( cut = this_cut_sp_D_el, eventweight=this_weight )

                    print ("sub_sample_C_el ={0} ".format(sub_sample_C_el.basecut.cutnamelist))

                    print ("C (el) - yields data: {0:.2f}".format(sp_C_el.numberstats()[0]))
                    print ("C (el) - yields bkg: {0:.2f}".format(sub_sample_C_el.numberstats()[0]))
                    sp_C_el = sp_C_el - sub_sample_C_el
                    print ("C (el) ==> yields after sub: {0:.2f}".format(sp_C_el.numberstats()[0]))

                    # *********************************************

                    print ("sub_sample_D_el ={0} ".format(sub_sample_D_el.basecut.cutnamelist))

                    print ("D (el) - yields data: {0:.2f}".format(sp_D_el.numberstats()[0]))
                    print ("D (el) - yields bkg: {0:.2f}".format(sub_sample_D_el.numberstats()[0]))
                    sp_D_el = sp_D_el - sub_sample_D_el
                    print ("D (el) ==> yields after sub: {0:.2f}".format(sp_D_el.numberstats()[0]))


                print ("---------------------------------------------------------------------\n")
                print ("C (el) - data yields after prompt/ch-flip subtraction: {0:.2f} +- {1:.2f}".format(sp_C_el.numberstats()[0], sp_C_el.numberstats()[1]))
                print ("D (el) - data yields after prompt/ch-flip subtraction: {0:.2f} +- {1:.2f}".format(sp_D_el.numberstats()[0], sp_D_el.numberstats()[1]))
                print ("\n---------------------------------------------------------------------")

                # Derive theta factors for el

                self.theta['El'] = self.calcTheta(sp_C_el,sp_D_el,stream='El')

            else :
                print ("Reading theta(el) value: {0} +- {1}".format(self.theta['El'][0], self.theta['El'][1]))


            if self.theta['Mu'][0] == 999.0 :

                print ("Calculating theta_mu from data in regions C,D...")

                # Define selection for region C and D (mumu)

                cut_sp_C_mu = basecut & self.vardb.getCut('2Lep_NJet_CR') & self.vardb.getCut('2Lep_MuMu_Event') & TTCut
                cut_sp_D_mu = basecut & self.vardb.getCut('2Lep_NJet_CR') & self.vardb.getCut('2Lep_MuMu_Event') & TL_LT_Cut

                sp_C_mu = sp.subprocess( cut = cut_sp_C_mu, eventweight=weight )
                sp_D_mu = sp.subprocess( cut = cut_sp_D_mu, eventweight=weight )

                print("Region C (mu) sp: {0}".format(cut_sp_C_mu.cutnamelist))
                print("Region D (mu) sp: {0}".format(cut_sp_D_mu.cutnamelist))

                # Get a list of stuff to subtract for region C and D (i.e, prompt MC and charge flips)

                sublist = [ item for item in self.parent.sub_backgrounds ]

                # ... and now subtract!

                for sample in sublist:

                    print ("Subtracting {0} from data in regions C,D...".format(sample))

                    this_cut_sp_C_mu = cut_sp_C_mu
                    this_cut_sp_D_mu = cut_sp_D_mu

                    this_weight = weightMC

                    # NB: here it is crucial to call .base() on the subprocess, otherwise the subprocess would have the cuts
                    # defined in its own __call__ method already applied, whcih in general is not what we want
                    # (e.g., it might have a TT selection applied, when we want to consider TL events instead...)

                    sub_sample_C_mu = self.parent.procmap[sample].base(treename,category,options)
                    sub_sample_D_mu = self.parent.procmap[sample].base(treename,category,options)

                    if ( sample == "QMisIDMC" ) or ( sample == "QMisID" ):
                        print("NO QMisID subtraction for muons!!")
                        continue

                    sub_sample_C_mu = sub_sample_C_mu.subprocess( cut = this_cut_sp_C_mu, eventweight=this_weight )
                    sub_sample_D_mu = sub_sample_D_mu.subprocess( cut = this_cut_sp_D_mu, eventweight=this_weight )

                    print ("sub_sample_C_mu ={0} ".format(sub_sample_C_mu.basecut.cutnamelist))

                    print ("C (mu) - yields data: {0:.2f}".format(sp_C_mu.numberstats()[0]))
                    print ("C (mu) - yields bkg: {0:.2f}".format(sub_sample_C_mu.numberstats()[0]))
                    sp_C_mu = sp_C_mu - sub_sample_C_mu
                    print ("C (mu) ==> yields after sub: {0:.2f}".format(sp_C_mu.numberstats()[0]))

                    # *********************************************

                    print ("sub_sample_D_mu ={0} ".format(sub_sample_D_mu.basecut.cutnamelist))

                    print ("D (mu) - yields data: {0:.2f}".format(sp_D_mu.numberstats()[0]))
                    print ("D (mu) - yields bkg: {0:.2f}".format(sub_sample_D_mu.numberstats()[0]))
                    sp_D_mu = sp_D_mu - sub_sample_D_mu
                    print ("D (mu) ==> yields after sub: {0:.2f}".format(sp_D_mu.numberstats()[0]))


                print ("---------------------------------------------------------------------\n")
                print ("C (mu) - data yields after prompt/ch-flip subtraction: {0:.2f} +- {1:.2f}".format(sp_C_mu.numberstats()[0],sp_C_mu.numberstats()[1]))
                print ("D (mu) - data yields after prompt/ch-flip subtraction: {0:.2f} +- {1:.2f}".format(sp_D_mu.numberstats()[0],sp_D_mu.numberstats()[1]))
                print ("\n---------------------------------------------------------------------")

                # Derive theta factors for mu

                self.theta['Mu'] = self.calcTheta(sp_C_mu,sp_D_mu,stream='Mu')

            else :
                print ("Reading theta(mu) value: {0:.3f} +- {1:.3f}".format(self.theta['Mu'][0], self.theta['Mu'][1]))


            # Define Region B,  depending on which flavour composition we are looking at:
	    # (NB now we must look at the original category cut!)

            cut_sp_B_SF     = category.cut & TL_LT_Cut
            cut_sp_B_OF_Lel = category.cut & (LelTmuCut | TmuLelCut)
            cut_sp_B_OF_Lmu = category.cut & (TelLmuCut | LmuTelCut)

            if not ("2Lep_OF_Event") in category.cut.cutname:

                sp_B = sp.subprocess( cut = cut_sp_B_SF, eventweight=weight )

                sublist = [ item for item in self.parent.sub_backgrounds ]
                for sample in sublist:

                    print ("Subtracting {0} from data in region B...".format(sample))

                    this_cut_sp_B_SF = cut_sp_B_SF
                    this_weight = weightMC

                    if ( ( sample == "QMisIDMC" ) or ( sample == "QMisID" ) ):

                        if ("2Lep_MuMu_Event") in category.cut.cutname:
                            print("NO QMisID subtraction for muons!!")
                            continue

                    if sample == "QMisIDMC":
                        this_cut_sp_B_SF = this_cut_sp_B_SF.swapCut(self.vardb.getCut('2Lep_TRUTH_PurePromptEvent'),self.vardb.getCut('2Lep_TRUTH_QMisIDEvent'))
                    if sample == "QMisID":
                        this_weight = 'QMisIDWeight'
                        this_cut_sp_B_SF = this_cut_sp_B_SF.removeCut(self.vardb.getCut('2Lep_TRUTH_PurePromptEvent'))
                        this_cut_sp_B_SF = this_cut_sp_B_SF.swapCut(self.vardb.getCut('2Lep_SS'), -self.vardb.getCut('2Lep_SS'))

                    sub_sample_B = self.parent.procmap[sample].base(treename,category,options)
                    sub_sample_B = sub_sample_B.subprocess( cut = this_cut_sp_B_SF, eventweight=this_weight )

                    print ("B (mumu,ee) - yields data: {0:.2f}".format(sp_B.numberstats()[0]))
                    print ("B (mumu,ee) - yields bkg: {0:.2f}".format(sub_sample_B.numberstats()[0]))

		    sp_B = sp_B - sub_sample_B

		    print ("B (mumu,ee) ==> yields after sub: {0:.2f}".format(sp_B.numberstats()[0]))

                print ("Region B (mumu,ee), yield: {0:.2f}\n".format(sp_B.numberstats()[0]))

                if ("2Lep_ElEl_Event") in category.cut.cutname:
                    sp_B = self.applyThetaFactor(sp_B,self.theta['El'])
                elif ("2Lep_MuMu_Event") in category.cut.cutname:
                    sp_B = self.applyThetaFactor(sp_B,self.theta['Mu'])
            else:

                sp_B_OF_Lel = sp.subprocess(cut=cut_sp_B_OF_Lel, eventweight=weight )
                sp_B_OF_Lmu = sp.subprocess(cut=cut_sp_B_OF_Lmu, eventweight=weight )

                sublist = [ item for item in self.parent.sub_backgrounds ]
                for sample in sublist:

                    print ("Subtracting {0} from data in region B...".format(sample))

                    this_cut_sp_B_OF_Lel = cut_sp_B_OF_Lel
                    this_cut_sp_B_OF_Lmu = cut_sp_B_OF_Lmu
                    this_weight = weightMC

                    if sample == "QMisIDMC":
                        this_cut_sp_B_OF_Lel = this_cut_sp_B_OF_Lel.swapCut(self.vardb.getCut('2Lep_TRUTH_PurePromptEvent'), self.vardb.getCut('2Lep_TRUTH_QMisIDEvent'))
                        this_cut_sp_B_OF_Lmu = this_cut_sp_B_OF_Lmu.swapCut(self.vardb.getCut('2Lep_TRUTH_PurePromptEvent'), self.vardb.getCut('2Lep_TRUTH_QMisIDEvent'))
		    if sample == "QMisID":
                        this_weight = 'QMisIDWeight'
                        this_cut_sp_B_OF_Lel = this_cut_sp_B_OF_Lel.removeCut(self.vardb.getCut('2Lep_TRUTH_PurePromptEvent'))
                        this_cut_sp_B_OF_Lel = this_cut_sp_B_OF_Lel.swapCut(self.vardb.getCut('2Lep_SS'), -self.vardb.getCut('2Lep_SS'))
                        this_cut_sp_B_OF_Lmu = this_cut_sp_B_OF_Lmu.removeCut(self.vardb.getCut('2Lep_TRUTH_PurePromptEvent'))
                        this_cut_sp_B_OF_Lmu = this_cut_sp_B_OF_Lmu.swapCut(self.vardb.getCut('2Lep_SS'), -self.vardb.getCut('2Lep_SS'))

                    sub_sample_B_OF_Lel = self.parent.procmap[sample].base(treename,category,options)
                    sub_sample_B_OF_Lel = sub_sample_B_OF_Lel.subprocess( cut = this_cut_sp_B_OF_Lel, eventweight=this_weight )
                    sub_sample_B_OF_Lmu = self.parent.procmap[sample].base(treename,category,options)
                    sub_sample_B_OF_Lmu = sub_sample_B_OF_Lmu.subprocess( cut = this_cut_sp_B_OF_Lmu, eventweight=this_weight )

		    print ("B (emu,mue, loose el) - yields data: {0:.2f}".format(sp_B_OF_Lel.numberstats()[0]))
                    print ("B (emu,mue, loose el) - yields bkg: {0:.2f}".format(sub_sample_B_OF_Lel.numberstats()[0]))
		    print ("B (emu,mue, loose mu) - yields data: {0:.2f}".format(sp_B_OF_Lmu.numberstats()[0]))
                    print ("B (emu,mue, loose mu) - yields bkg: {0:.2f}".format(sub_sample_B_OF_Lmu.numberstats()[0]))

                    sp_B_OF_Lel = sp_B_OF_Lel - sub_sample_B_OF_Lel
                    sp_B_OF_Lmu = sp_B_OF_Lmu - sub_sample_B_OF_Lmu

		    print ("B (emu,mue, loose el) ==> yields after sub: {0:.2f}".format(sp_B_OF_Lel.numberstats()[0]))
                    print ("B (emu,mue, loose mu) ==> yields after sub: {0:.2f}".format(sp_B_OF_Lmu.numberstats()[0]))

                print ("Region B (emu,mue) \n yield (loose el): {0:.2f}\n yield (loose mu): {1:.2f}\n".format(sp_B_OF_Lel.numberstats()[0],sp_B_OF_Lmu.numberstats()[0]))

                sp_B = self.applyThetaFactor(sp_B_OF_Lmu,self.theta['Mu']) + self.applyThetaFactor(sp_B_OF_Lel,self.theta['El'])

            print ("\n=========> Final fakes yield: {0:.2f} +- {1:.2f}\n".format(sp_B.numberstats()[0],sp_B.numberstats()[1]))

            return sp_B


    class FakesClosureMM(Process):

        latexname = 'Fakes MM - t#bar{t},t#bar{t}#gamma'
        colour = kTeal -9

        name = 'FakesClosureMM'

        def base(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    inputgroup = [
                ('tops', 'ttbar_nonallhad_Pythia8'),
                ('tops', 'ttgamma'),
                #('tops', 'ttbar_nonallhad'),
                #('tops', 'ttbar_dilep'),
                #('tops', 'ttbar_SingleLeptonP_MEPS_NLO'),
                #('tops', 'ttbar_SingleLeptonM_MEPS_NLO'),
                #('tops', 'ttbar_dilepton_MEPS_NLO'),
                ]

            if self.parent.readGFW2:
                del inputgroup[:]
                inputgroup = [
                    ('tops_MMClosure', 'ttbar_nonallhad_Pythia8'),
                    ('tops_MMClosure', 'ttgamma'),
                    #('tops_MMClosure', 'ttgamma_vPP6'),
                    #('tops_MMClosure', 'ttgamma_vSherpa'),
                    #('tops_MMClosure', 'ttbar_nonallhad'),
                    #('tops_MMClosure', 'ttbar_SingleLeptonP_MEPS_NLO'),
                    #('tops_MMClosure', 'ttbar_SingleLeptonM_MEPS_NLO'),
                    #('tops_MMClosure', 'ttbar_dilepton_MEPS_NLO'),
                ]

            print("\n{0}:\n".format(self.__class__.__name__))
	    print("\n".join("{0} - {1} : {2}".format(idx,sample[0],sample[1]) for idx, sample in enumerate(inputgroup)))

            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) * self.parent.norm_factor
            return sp

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

            debugflag = any( proc in self.parent.debugprocs for proc in [self.__class__.__name__,"ALL"])

            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)

            tt_ttgamma_OLR_cut = Cut("TTBar_TTGamma_OLR","( ( mc_channel_number == 410082 && m_MEphoton_OLtty_cat1 ) || ( mc_channel_number != 410082 && m_MEphoton_OLtty_keepEvent && !m_hasMEphoton_DRgt02_nonhad ) )")

            TTcut  = ('','FakesSideband_TT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]
            TLcut  = ('','FakesSideband_TL')[any( c == self.parent.channel for c in ['2LSS','3L'] )]
            LTcut  = ('','FakesSideband_LT')[any( c == self.parent.channel for c in ['2LSS','3L'] )]
            LLcut  = ('','FakesSideband_LL')[any( c == self.parent.channel for c in ['2LSS','3L'] )]
            weight = (None,'MMWeight')[any( c == self.parent.channel for c in ['2LSS','3L'] )]

            if weight and self.parent.readGFW2:
                weight = "MM_EventWeight"

            if debugflag:

                sp_TT_preweight = sp.subprocess(cut=category.cut & self.vardb.getCut(TTcut) & tt_ttgamma_OLR_cut)
                sp_TL_preweight = sp.subprocess(cut=category.cut & self.vardb.getCut(TLcut) & tt_ttgamma_OLR_cut)
                sp_LT_preweight = sp.subprocess(cut=category.cut & self.vardb.getCut(LTcut) & tt_ttgamma_OLR_cut)
                sp_LL_preweight = sp.subprocess(cut=category.cut & self.vardb.getCut(LLcut) & tt_ttgamma_OLR_cut)

                print(" ")
                print("{0} - TT : TTBar Fakes (pre-weighting) = {1:.2f} +- {2:.2f}".format(self.__class__.__name__,sp_TT_preweight.numberstats()[0],sp_TT_preweight.numberstats()[1]))
                print("{0} - TL : TTBar Fakes (pre-weighting) = {1:.2f} +- {2:.2f}".format(self.__class__.__name__,sp_TL_preweight.numberstats()[0],sp_TL_preweight.numberstats()[1]))
                print("{0} - LT : TTBar Fakes (pre-weighting) = {1:.2f} +- {2:.2f}".format(self.__class__.__name__,sp_LT_preweight.numberstats()[0],sp_LT_preweight.numberstats()[1]))
                print("{0} - LL : TTBar Fakes (pre-weighting) = {1:.2f} +- {2:.2f}".format(self.__class__.__name__,sp_LL_preweight.numberstats()[0],sp_LL_preweight.numberstats()[1]))

            sp_TT  = sp.subprocess(cut=category.cut & self.vardb.getCut(TTcut) & tt_ttgamma_OLR_cut, eventweight=weight)
            sp_TL  = sp.subprocess(cut=category.cut & self.vardb.getCut(TLcut) & tt_ttgamma_OLR_cut, eventweight=weight)
            sp_LT  = sp.subprocess(cut=category.cut & self.vardb.getCut(LTcut) & tt_ttgamma_OLR_cut, eventweight=weight)
            sp_LL  = sp.subprocess(cut=category.cut & self.vardb.getCut(LLcut) & tt_ttgamma_OLR_cut, eventweight=weight)

            if debugflag:

                print(" ")
                print("{0} - TT cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp_TT.basecut.cutnamelist, weight))
                print("{0} - TL cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp_TL.basecut.cutnamelist, weight))
                print("{0} - LT cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp_LT.basecut.cutnamelist, weight))
                print("{0} - LL cuts: {1}, process weight: {2}".format(self.__class__.__name__,sp_LL.basecut.cutnamelist, weight))
                print(" ")
                print("{0} - TT : TTBar Fakes = {1:.2f} +- {2:.2f}".format(self.__class__.__name__,sp_TT.numberstats()[0],sp_TT.numberstats()[1]))
                print("{0} - TL : TTBar Fakes = {1:.2f} +- {2:.2f}".format(self.__class__.__name__,sp_TL.numberstats()[0],sp_TL.numberstats()[1]))
                print("{0} - LT : TTBar Fakes = {1:.2f} +- {2:.2f}".format(self.__class__.__name__,sp_LT.numberstats()[0],sp_LT.numberstats()[1]))
                print("{0} - LL : TTBar Fakes = {1:.2f} +- {2:.2f}".format(self.__class__.__name__,sp_LL.numberstats()[0],sp_LL.numberstats()[1]))
                print(" ")

            sp = sp_TT + sp_TL + sp_LT + sp_LL

            return sp


    class FakesClosureTHETA(Process):

        latexname = 'Fakes #theta method - t#bar{t}'
        colour = kCyan - 9

        def base(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

	    inputgroup = [
                    ('tops', 'ttbar_nonallhad_Pythia8'),
                    #('tops', 'ttbar_nonallhad'),
                ]

            print("\n{0}:\n".format(self.__class__.__name__))
	    print("\n".join("{0} - {1} : {2}".format(idx,sample[0],sample[1]) for idx, sample in enumerate(inputgroup)))

            trees = self.inputs.getTrees(treename, inputgroup)
            return self.subprocess(trees=trees) * self.parent.norm_factor

        def calcTheta(self, sp_num, sp_denom, stream=None, options={}):

            if not sp_denom.number():
                print ("ERROR: Cannot calculate theta transfer factor! Denominator = 0")

            #print("N: ", sp_num.numberstats())
            #print("D: ", sp_denom.numberstats())

            theta = (sp_num/sp_denom).numberstats()

            print ("***********************************************************************\n")
            print ("Calculated theta transfer factor for stream {0}: theta = {1} +- {2}".format(stream,theta[0],theta[1]))
            print ("\n***********************************************************************")

            return theta

        def applyThetaFactor(self, sp, theta, options={}):
            systematics = options.get('systematics', None)
            systematicsdirection = options.get('systematicsdirection', None)

            if systematics:
                if systematicsdirection == 'UP':
                    systdir = 1.0
                elif systematicsdirection == 'DOWN':
                    systdir = -1.0
                theta = (theta[0] + theta[1]*systdir, theta[1])

            ssp = sp.subprocess()
            ssp *= theta[0]
            return ssp

        def __call__(self, treename='physics', category=None, options={}):

            if self.parent.readGFW2:
                treename = 'nominal'

            print("\nFakesClosureTHETA\n")

            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction

            TTCut =''
            TTCut =''
            TLCut =''
            LTCut =''
            TelLmuCut =''
            LelTmuCut =''
            TmuLelCut =''
            LmuTelCut =''
            weight = None
            if self.parent.channel=='2LSS':
                TTCut  = self.vardb.getCut('FakesSideband_TT')
                TLCut  = self.vardb.getCut('FakesSideband_TL')
                LTCut  = self.vardb.getCut('FakesSideband_LT')
                TelLmuCut = self.vardb.getCut('FakesSideband_TelLmu')
                LelTmuCut = self.vardb.getCut('FakesSideband_LelTmu')
                TmuLelCut = self.vardb.getCut('FakesSideband_TmuLel')
                LmuTelCut = self.vardb.getCut('FakesSideband_LmuTel')

            TL_LT_Cut = (TLCut | LTCut)

            # Take the base suprocess (i.e, TTBar)

            sp = self.base(treename, category, options)

            # Cache the category cut in a base cut which can be modified at no risk

            basecut = category.cut

            print("base sp: {0}".format(basecut.cutnamelist))

            # Remove the cuts defining the flavour composition (this is made just for calculating thetas...)

            basecut = basecut.removeCut(self.vardb.getCut('2Lep_ElEl_Event'))
            basecut = basecut.removeCut(self.vardb.getCut('2Lep_MuMu_Event'))
            basecut = basecut.removeCut(self.vardb.getCut('2Lep_OF_Event'))

            # Remove the cuts defiining the jet multiplicity

            basecut = basecut.removeCut(self.vardb.getCut('2Lep_NJet_SR'))
            basecut = basecut.removeCut(self.vardb.getCut('2Lep_NJet_CR'))

            if TTHBackgrounds.theta_MC['El'][0] == 999.0 :

                print ("Calculating theta_el from TTBar in regions C,D...")

                # Define selection for region C and D (ee)

                cut_sp_C_el = basecut & self.vardb.getCut('2Lep_NJet_CR') & self.vardb.getCut('2Lep_ElEl_Event') & TTCut & self.vardb.getCut('2Lep_Zsidescut') & self.vardb.getCut('2Lep_ElEtaCut')
                cut_sp_D_el = basecut & self.vardb.getCut('2Lep_NJet_CR') & self.vardb.getCut('2Lep_ElEl_Event') & TL_LT_Cut & self.vardb.getCut('2Lep_Zsidescut') & self.vardb.getCut('2Lep_ElEtaCut')

                sp_C_el = sp.subprocess( cut = cut_sp_C_el, eventweight=weight )
                sp_D_el = sp.subprocess( cut = cut_sp_D_el, eventweight=weight )

                print("Region C (el) sp: {0}".format(cut_sp_C_el.cutnamelist))
                print("Region D (el) sp: {0}".format(cut_sp_D_el.cutnamelist))

                # No subtraction in closure test!

                print ("---------------------------------------------------------------------\n")
                print ("C (el) - TTBar yields: ", sp_C_el.numberstats())
                print ("D (el) - TTBar yields: ", sp_D_el.numberstats())
                print ("\n---------------------------------------------------------------------")

                # Derive theta factor for el

                TTHBackgrounds.theta_MC['El'] = self.calcTheta(sp_C_el,sp_D_el,stream='El')

            else :
                print ("Reading theta(el) value: {0} +- {1}".format(TTHBackgrounds.theta_MC['El'][0], TTHBackgrounds.theta_MC['El'][1]))


            if TTHBackgrounds.theta_MC['Mu'][0] == 999.0 :

                print ("Calculating theta_mu from TTBar in regions C,D...")

                # Define selection for region C and D (mumu)

                cut_sp_C_mu = basecut & self.vardb.getCut('2Lep_NJet_CR') & self.vardb.getCut('2Lep_MuMu_Event') & TTCut
                cut_sp_D_mu = basecut & self.vardb.getCut('2Lep_NJet_CR') & self.vardb.getCut('2Lep_MuMu_Event') & TL_LT_Cut

                sp_C_mu = sp.subprocess( cut = cut_sp_C_mu, eventweight=weight )
                sp_D_mu = sp.subprocess( cut = cut_sp_D_mu, eventweight=weight )

                print("Region C (mu) sp: {0}".format(cut_sp_C_mu.cutnamelist))
                print("Region D (mu) sp: {0}".format(cut_sp_D_mu.cutnamelist))

                # No subtraction in closure test!

                print ("---------------------------------------------------------------------\n")
                print ("C (mu) - TTBar yields: ", sp_C_mu.numberstats())
                print ("D (mu) - TTBar yields: ", sp_D_mu.numberstats())
                print ("\n---------------------------------------------------------------------")

                # Derive theta factor for mu

                TTHBackgrounds.theta_MC['Mu'] = self.calcTheta(sp_C_mu,sp_D_mu,stream='Mu')

            else :
                print ("Reading theta(mu) value: {0} +- {1}".format(TTHBackgrounds.theta_MC['Mu'][0], TTHBackgrounds.theta_MC['Mu'][1]))


            # Define Region B,  depending on which flavour composition we are looking at:
            # Take TTbar MC events with fakes, vetoing all prompts and charge flips, and reweight it by the theta factors measured in ttbar MC

            cut_sp_B_SF     = category.cut.swapCut(self.vardb.getCut('2Lep_NJet_CR'),self.vardb.getCut('2Lep_NJet_SR')) & TL_LT_Cut
            cut_sp_B_OF_Lel = category.cut.swapCut(self.vardb.getCut('2Lep_NJet_CR'),self.vardb.getCut('2Lep_NJet_SR')) & (LelTmuCut | TmuLelCut)
            cut_sp_B_OF_Lmu = category.cut.swapCut(self.vardb.getCut('2Lep_NJet_CR'),self.vardb.getCut('2Lep_NJet_SR')) & (TelLmuCut | LmuTelCut)

            if not ("2Lep_OF_Event") in category.cut.cutname:
                sp_B = self.parent.procmap['TTBarClosure'].base(treename,category,options)
                sp_B = sp_B.subprocess(cut=cut_sp_B_SF,eventweight=weight)
                if ("2Lep_ElEl_Event") in category.cut.cutname:
                    sp_B = self.applyThetaFactor(sp_B,TTHBackgrounds.theta_MC['El'])
                elif ("2Lep_MuMu_Event") in category.cut.cutname:
                    sp_B = self.applyThetaFactor(sp_B,TTHBackgrounds.theta_MC['Mu'])
            else:
                sp_B_Lel = self.parent.procmap['TTBarClosure'].base(treename,category,options)
                sp_B_Lmu = self.parent.procmap['TTBarClosure'].base(treename,category,options)
                sp_B_Lel = sp_B_Lel.subprocess(cut=cut_sp_B_OF_Lel, eventweight=weight)
                sp_B_Lmu = sp_B_Lmu.subprocess(cut=cut_sp_B_OF_Lmu, eventweight=weight)
                sp_B = self.applyThetaFactor(sp_B_Lmu,TTHBackgrounds.theta_MC['Mu']) + self.applyThetaFactor(sp_B_Lel,TTHBackgrounds.theta_MC['El'])

            print ("=================>\n")
            print ("Region B sp: {0} \n".format(sp_B.basecut.cutnamelist))
            print ("fakes yield: ", sp_B.numberstats() ,"\n")

            return sp_B

    # Check Data vs. TTBar MC in low njet, TL,LT, e-mu region
    # Use the theta factors derived in data to reweight TTBar MC

    class FakesClosureDataTHETA(Process):

        latexname = 'Fakes #theta method'
        colour = kCyan - 9

        def base(self, treename='physics', category=None, options={}):

	    inputgroup = [
                    ('tops', 'ttbar_nonallhad_Pythia8'),
                    #('tops', 'ttbar_nonallhad'),
                ]

            print("\n{0}:\n".format(self.__class__.__name__))
	    print("\n".join("{0} - {1} : {2}".format(idx,sample[0],sample[1]) for idx, sample in enumerate(inputgroup)))

            trees = self.inputs.getTrees(treename, inputgroup)
            return self.subprocess(trees=trees) * self.parent.norm_factor

        def __call__(self, treename='physics', category=None, options={}):

            print 'FakesClosureDataTHETA \n '

            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction

            TTCut =''
            TTCut =''
            TLCut =''
            LTCut =''
            TelLmuCut =''
            LelTmuCut =''
            TmuLelCut =''
            LmuTelCut =''
            weight = None
            weightMC = 'weight_event_trig * weight_event_lep * tauSFTight * JVT_EventWeight * MV2c10_70_EventWeight'

            if self.parent.channel=='2LSS':
                TTCut  = self.vardb.getCut('FakesSideband_TT')
                TLCut  = self.vardb.getCut('FakesSideband_TL')
                LTCut  = self.vardb.getCut('FakesSideband_LT')
                TelLmuCut = self.vardb.getCut('FakesSideband_TelLmu')
                LelTmuCut = self.vardb.getCut('FakesSideband_LelTmu')
                TmuLelCut = self.vardb.getCut('FakesSideband_TmuLel')
                LmuTelCut = self.vardb.getCut('FakesSideband_LmuTel')

            TL_LT_Cut = (TLCut | LTCut)

            # Take the base suprocess (i.e, TTBar)

            sp = self.base(treename, category, options)

            print("base sp: {0}".format(category.cut.cutnamelist))

            # Define closure region (OF, Tl,LT, low njet)
            #
            # Veto all prompt and charge flips (for TTBar, basically it's a charge flip veto)

            cut_sp_OF_Lel = category.cut & (LelTmuCut | TmuLelCut)
            cut_sp_OF_Lmu = category.cut & (TelLmuCut | LmuTelCut)
            cut_sp_OF_Lel = cut_sp_OF_Lel.swapCut(self.vardb.getCut('2Lep_TRUTH_PurePromptEvent'),self.vardb.getCut('2Lep_TRUTH_NonPromptEvent'))
            cut_sp_OF_Lmu = cut_sp_OF_Lmu.swapCut(self.vardb.getCut('2Lep_TRUTH_PurePromptEvent'),self.vardb.getCut('2Lep_TRUTH_NonPromptEvent'))

            # Plug in the theta factors by hand...CHANGE ME!

            sp_Lel = sp.subprocess(cut=cut_sp_OF_Lel, eventweight=weightMC )
            sp_Lmu = sp.subprocess(cut=cut_sp_OF_Lmu, eventweight=weightMC )
            sp_final = ( sp_Lel * TTHBackgrounds.theta['El'][0] ) + ( sp_Lmu * TTHBackgrounds.theta['Mu'][0] )

            print ("=================>\n")
            print ("Region closure sp: {0}".format(sp_final.basecut.cutnamelist))
            print ("yield: ", sp_final.numberstats())

            return sp_final
