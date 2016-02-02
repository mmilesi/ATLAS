import os, sys, math, types

sys.path.append(os.path.abspath(os.path.curdir))
from Plotter.BackgroundTools_ttH import loadSamples, drawText, Category, Background, Process, VariableDB, Variable, Cut, Systematics, Category
from ROOT import TColor, kBlack, kWhite, kGray, kBlue, kRed, kYellow, kGreen, kAzure, kTeal, kSpring, kOrange, kCyan, TLegend, TLatex, TCanvas, TH1I, TFile

class MyCategory(Category):
    def __init__(self, name, cut=None, controlcut=None, overridebins=None):
        self.name = name
        self.tokens = name.split(' ')
        self.cut = cut
        self.controlcut = controlcut
        self.overridebins = overridebins
        self.stream = '*'
        self.streamname = '*'
        self.jet = name
        #self.met = 'HMET'

        if len(self.tokens)>=2:
            self.stream = self.tokens[0]
            if self.stream == 'El': self.streamname = 'Egamma'
            elif self.stream == 'Mu': self.streamname = 'Muons'
            self.jet = self.tokens[1]
            #self.met = self.tokens[2]

class TTHBackgrounds2015(Background):
    backgrounds = [ 'TTBarW', 'TTBarZ', 'Top', 'TopCF', 'Diboson', 'DibosonCF', 'HtoZZ', 'ZjetsLF', 'Zjets', 'Wjets', 'Prompt', 'ChargeFlipMC', 'FakesFF', 'FakesMM', 'FakesABCD']
    sub_backgrounds = [ 'TTBarW', 'TTBarZ', 'Top', 'TTBar', 'Diboson', 'Zjets', 'Wjets', 'ChargeFlipMC']
    signals     = ['TTBarH']
    observed    = ['Observed']
    luminosity  = 1.0
    lumi_units  = 'fb-1'
    channel     = 'TwoLepSS' # can be one of ['TwoLepSS' 'TwoLepCR', 'ThreeLep', 'FourLep']
    eventweight = 1.0
    #eventweight = 'weight_pileup*weight_muon_trig*weight_electron_trig'
    #eventweight = 'evtsel_weight*evtsel_weight_el*evtsel_weight_mu*evtsel_bjet_weight*evtsel_weight_lep_trigger*weight_CF'
    #eventweight = 'evtsel_weight*evtsel_weight_el*evtsel_weight_mu*evtsel_bjet_weight*evtsel_weight_lep_trigger*weight_MM*weight_CF'
    #eventweight = 'evtsel_weight*evtsel_weight_el*evtsel_weight_mu*evtsel_bjet_weight*evtsel_weight_lep_trigger*weight_FF*weight_CF'
    useEmbedding    = False
    useZCorrections = False
    RQCD = {
        #'El': (1.000, 0.051),
        #'Mu': (1.177, 0.066),
        'El': (1.00, 0.05),
        'Mu': (1.11, 0.08),
    }

    theta = {
    	'El': (999.0, 0.0),
    	'Mu': (999.0, 0.0),
    }

    def str_to_class(self, field):
        try:
            identifier = getattr(self, field)
        except AttributeError:
            raise NameError("%s doesn't exist." % field)
        if isinstance(identifier, (types.ClassType, types.TypeType)):
            return identifier
        raise TypeError("%s is not a class." % field)


    def applyKfactor(self, sp, category, kfactor, options):
        #function used to multiply a sample for a kfactor treating properly the systematics if needed
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
        #uses the new scale factors for BDT tight instead of the old BDT medium. ATTENTION the correction is done only with the nominal coefficients. To have the right coefficients for systematics shifted correct the corrections code and produce again the ntuples
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
            leg2 = TLegend(0.70,lower,0.90,0.92)
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
	#
        if self.lumi_units == 'fb-1':
	    lumtext = drawText(text="  #int L dt = %.1f fb^{-1}"%(self.luminosity), x=.2, y=.87, size=0.03 * scale)
	# for O(pb-1) luminosity
	#
        elif self.lumi_units == 'pb-1':
	    lumtext = drawText(text="  #int L dt = %.2f pb^{-1}"%(self.luminosity*1000), x=.2, y=.87, size=0.04 * scale)

        cmetext = drawText(text="         #sqrt{s} = 13 TeV", x=.2, y=.82, size=0.03 * scale)
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

        latexname = 'Data 2015'

        def base(self, treename='physics', category=None, options={}):
            #Contains the instuction of which tree load and eventually correct the units of the cross-setion (see the division /1000.). Note it is not automatically executed when the class is called. It is executed trought __call__
            inputgroup = [
                    ('Data', 'physics_Main'),
                ]
            trees = self.inputs.getTrees(treename, inputgroup)
            return self.subprocess(trees=trees)

        def __call__(self, treename='physics', category=None, options={}):

	    # This is the default procedure which is applied to the class when it is called in the code.
	    # Here in addition to the load of the tree are operated also the selection, application of factors and systematics
            #
	    sp = self.base(treename, category, options)
            TTcut=''
            weight=1.0
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

	    sp = sp.subprocess(cut=category.cut,eventweight=weight)

	    if TTcut is not '':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

	    print("Observed sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

            return sp


    class TTBarH(Process):

	latexname = 't#bar{t} H'
        colour = kBlack

        def base(self, treename='physics', category=None, options={}):

	    #hmass = options.get('hmass', '125') #300 is the default value if hmass is not in options, hmass is specified in one option passed to the plot function in the plotting script
            inputgroup = [
                #('ttH', hmass),
                #('ttH', '125'),
                ('ttH', 'ttH_inc_dil'),
                ('ttH', 'ttH_inc_semil'),
                ('ttH', 'ttH_inc_allhad'),
                         ]
            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) / 1000. # the division by 1000. is used to have the correct cross section
            return sp

        def __call__(self, treename='physics', category=None, options={}):

            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)
            TTcut=''
            weight=1.0
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

            if TTcut is not '':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

	    print("TTBarH sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

            return sp

    class TTBarW(Process):

	latexname = 't#bar{t} W'
        colour = kRed - 4

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
		('tops', 'ttW'),
		#('tops', 'Sherpa_ttW'),
                         ]
            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) / 1000.
            return sp

        def __call__(self, treename='physics', category=None, options={}):
            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)
            TTcut=''
            weight=1.0
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

            if TTcut is not '':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

	    print("TTBarW sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

            return sp


    class TTBarZ(Process):

	latexname = 't#bar{t} Z'
        colour = kRed - 7

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
		('tops', 'ttZnnqq'),
                ('tops', 'ttee'),
                ('tops', 'ttmumu'),
		('tops', 'tttautau'),
                #('tops', 'Sherpa_ttZnnqq'),
                #('tops', 'Sherpa_ttll'),
                         ]
            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) / 1000.
            return sp

        def __call__(self, treename='physics', category=None, options={}):
            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)
            TTcut=''
            weight=1.0
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

            if TTcut is not '':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

	    print("TTBarZ sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

            return sp


    class Zeejets(Process):

	latexname = 'Z/#gamma*#rightarrow#it{ee}+jets'
        colour = kGreen-7

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
		           #('Z+jets', 'ee'),
                           #('MadGraphZ+jets', 'ee'),
                           ('Z+jetsCVetoBVeto', 'ee'),
			   ('Z+jetsCFilterBVeto', 'ee'),
                           ('Z+jetsBFilter', 'ee'),
		           #('DYZ+jets', 'ee'),
		         ]
            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) / 1000.
            if self.parent.useZCorrections:
                sp = sp*0.94
            return sp

        def __call__(self, treename='physics', category=None, options={}):
            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)
            TTcut= ''
            weight= 1.0
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

	    if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

	    print("Zeejets sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

            return sp


    class Zmumujets(Process):

	latexname = 'Z/#gamma*#rightarrow#mu#mu+jets'
        colour = kTeal+2

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
		           #('Z+jets', 'mumu'),
                           #('MadGraphZ+jets', 'mumu'),
                           ('Z+jetsCVetoBVeto', 'mumu'),
			   ('Z+jetsCFilterBVeto', 'mumu'),
                           ('Z+jetsBFilter', 'mumu'),
		           #('DYZ+jets', 'mumu'),
		         ]
            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) / 1000.
            if self.parent.useZCorrections:
                sp = sp*0.94
            return sp

        def __call__(self, treename='physics', category=None, options={}):
            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)
            TTcut=''
            weight=1.0
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

	    if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

	    print("Zmumujets sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

	    return sp

    class Ztautaujets(Process):

	latexname = 'Z/#gamma*#rightarrow#tau#tau+jets'
        colour = kTeal

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
		           #('Z+jets', 'tautau'),
                           #('MadGraphZ+jets', 'tautau'),
                           ('Z+jetsCVetoBVeto', 'tautau'),
			   ('Z+jetsCFilterBVeto', 'tautau'),
                           ('Z+jetsBFilter', 'tautau'),
		           #('DYZ+jets', 'tautau'),
		         ]
            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) / 1000.
            if self.parent.useZCorrections:
                sp = sp*0.94
            return sp

        def __call__(self, treename='physics', category=None, options={}):
            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)
            TTcut=''
            weight=1.0
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

	    if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

	    print("Ztautaujets sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

            return sp


    class Zjets(Process):

	latexname = 'Z/#gamma* + jets'
        colour = kGreen

        def base(self, treename='physics', category=None, options={}):

	    zee = self.parent.procmap['Zeejets'].base(treename, category, options)
            zmumu = self.parent.procmap['Zmumujets'].base(treename, category, options)
            ztautau = self.parent.procmap['Ztautaujets'].base(treename, category, options)

            return (zee + zmumu + ztautau)

        def __call__(self, treename='physics', category=None, options={}):

	    zee = self.parent.procmap['Zeejets'](treename, category, options)
            zmumu = self.parent.procmap['Zmumujets'](treename, category, options)
            ztautau = self.parent.procmap['Ztautaujets'](treename, category, options)

            return (zee + zmumu + ztautau)


    class ZjetsLF(Process):

	latexname = 'Z/#gamma* + LF jets'
        colour = kGreen

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
		           ('Z+jetsCVetoBVeto', '*'),
		         ]
            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) / 1000.
            if self.parent.useZCorrections:
                sp = sp*0.94
            return sp

        def __call__(self, treename='physics', category=None, options={}):
            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)
            TTcut=''
            weight=1.0
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

	    if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

	    print("ZjetsLF sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

            return sp

    class ZjetsHF(Process):

        latexname = 'Z/#gamma* + HF jets'
        colour = kGreen

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
		           ('Z+jetsCFilterBVeto', '*'),
		           ('Z+jetsBFilter', '*'),
                         ]
            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) / 1000.
            if self.parent.useZCorrections:
                sp = sp*1.5
            return sp

        def __call__(self, treename='physics', category=None, options={}):
            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)
            TTcut=''
            weight=1.0
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

	    if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

	    print("ZjetsLF sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

            return sp


    class ZjetsCF(Process):

	latexname = 'Z/#gamma*+jets (ChFlip only)'
        colour = kAzure + 10

	def base(self, treename='physics', category=None, options={}):
            z = self.parent.procmap['Zjets'].base(treename, category, options)
            return z

        def __call__(self, treename='physics', category=None, options={}):

            z = self.parent.procmap['Zjets'](treename, category, options)

	    # plot only events where at least one lepton is charge flip. Remove req. where both lep must be prompt
	    #
	    truth_cut = category.cut.swapCut(self.vardb.getCut('2Lep_PurePromptEvent'),self.vardb.getCut('2Lep_ChFlipEvent'))

            z = z.subprocess(cut=truth_cut)

	    print("Zjets CF sp: {0}".format(z.basecut.cutnamelist))

            return z


    class Wenujets(Process):

	latexname = 'W#rightarrow#it{e}#nu+jets'
        colour = kYellow

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                #('PowhegPythiaW+jets', 'Wplusenu'),
                #('PowhegPythiaW+jets', 'Wminusenu'),
                #('MadGraphW+jets', 'enu'),
                ('W+jetsBFilter', 'enu'),
                ('W+jetsCFilterBVeto', 'enu'),
                ('W+jetsCVetoBVeto', 'enu'),
                         ]
            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) / 1000.
            return sp

        def __call__(self, treename='physics', category=None, options={}):
            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)
            TTcut=''
            weight=1.0
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

	    if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

	    print("Wenujets sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

            return sp

    class Wmunujets(Process):

	latexname = 'W#rightarrow#mu#nu+jets'
        colour = kYellow

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                #('PowhegPythiaW+jets', 'Wplusmunu'),
                #('PowhegPythiaW+jets', 'Wminusmunu'),
                #('MadGraphW+jets', 'munu'),
                ('W+jetsBFilter', 'munu'),
                ('W+jetsCFilterBVeto', 'munu'),
                ('W+jetsCVetoBVeto', 'munu'),
                         ]
            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) / 1000.
            return sp

        def __call__(self, treename='physics', category=None, options={}):
            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)
            TTcut=''
            weight=1.0
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

	    if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

	    print("Wmunujets sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

            return sp

    class Wtaunujets(Process):

	latexname = 'W#rightarrow#tau#nu+jets'
        colour = kYellow

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                #('PowhegPythiaW+jets', 'Wplustaunu'),
                #('PowhegPythiaW+jets', 'Wminustaunu'),
                #('MadGraphW+jets', 'taunu'),
                ('W+jetsBFilter', 'taunu'),
                ('W+jetsCFilterBVeto', 'taunu'),
                ('W+jetsCVetoBVeto', 'taunu'),
                         ]
            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) / 1000.
            return sp

        def __call__(self, treename='physics', category=None, options={}):
            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)
            TTcut=''
            weight=1.0
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

	    if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

	    print("Wtaunujets sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

            return sp

    class Wjets(Process):

        latexname = 'W+jets'
        colour = kYellow

        def base(self, treename='physics', category=None, options={}):

	    wenu = self.parent.procmap['Wenujets'].base(treename, category, options)
            wmunu = self.parent.procmap['Wmunujets'].base(treename, category, options)
            wtaunu = self.parent.procmap['Wtaunujets'].base(treename, category, options)

            return (wenu + wmunu + wtaunu)

        def __call__(self, treename='physics', category=None, options={}):

	    wenu = self.parent.procmap['Wenujets'](treename, category, options)
            wmunu = self.parent.procmap['Wmunujets'](treename, category, options)
            wtaunu = self.parent.procmap['Wtaunujets'](treename, category, options)

	    return (wenu + wmunu + wtaunu)

    class Top(Process):

        #latexname = 'tW, tZ, tt#it{l}{l}, ttWW, 4t, single t'
        latexname = 'rare top'

        colour = kAzure + 1

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                ('tops', 'tZ'),
                ('tops', 'tW'),
                ('tops', 'singlet'),
                ('tops', '4top'),
                ('tops', 'ttWW'),
                         ]
            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) / 1000.
            return sp

        def __call__(self, treename='physics', category=None, options={}):
            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)
            TTcut=''
            weight=1.0

            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

	    if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

	    print("Top sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

            return sp


    class TTBar(Process):

        latexname = 't#bar{t}'
        colour = kAzure + 8

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                ('tops', 'ttbar_nonallhad'),
                #('tops', 'ttbar_dilep'),
                         ]
            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) / 1000.
            return sp

        def __call__(self, treename='physics', category=None, options={}):
            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)
            TTcut  =''
            weight=1.0

            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

            if TTcut is not '':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

	    print("TTBar sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

            return sp

    # the only difference with the above class is
    # that we plot only events with at least one !prompt lepton, and veto charge flip
    #
    class TTBarClosureMM(Process):

        latexname = 't#bar{t}'
        colour = kAzure + 8

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                ('tops', 'ttbar_nonallhad'),
                #('tops', 'ttbar_dilep'),
                         ]
            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) / 1000.
            return sp

        def __call__(self, treename='physics', category=None, options={}):

	    systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)
            TTcut  =''
            weight=1.0

            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

            if TTcut is not '':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

	    # plot only events where at least one lepton is !prompt, and none is charge flip
	    #
            sp = sp.subprocess(cut=self.vardb.getCut('2Lep_NonPromptEvent'), eventweight=weight)

	    print("TTBarClosureMM sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

            return sp


    class TopCF(Process):

        latexname = 'tops (ChFlip only)'
        colour = kAzure - 4

        def base(self, treename='physics', category=None, options={}):

	    raretop = self.parent.procmap['Top'].base(treename, category, options)
            ttbar   = self.parent.procmap['TTBar'].base(treename, category, options)

            return (raretop + ttbar)

        def __call__(self, treename='physics', category=None, options={}):

            sp = self.base(treename, category, options)

	    raretop = self.parent.procmap['Top'](treename, category, options)
            ttbar   = self.parent.procmap['TTBar'](treename, category, options)

	    topcf = ttbar + raretop

	    # plot only events where at least one lepton is charge flip. Remove req. where both lep must be prompt
	    #
            truthcut = category.cut.swapCut(self.vardb.getCut('2Lep_PurePromptEvent'),self.vardb.getCut('2Lep_ChFlipEvent'))

	    topcf = topcf.subprocess(cut=truthcut)

	    print("TopCF sp: {0}".format(topcf.basecut.cutnamelist))

            return topcf


    class Diboson(Process):

        latexname = 'WZ, ZZ, WW'
        colour = kYellow - 9

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                ('Diboson', '*'),
                         ]
            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) / 1000.
            return sp

        def __call__(self, treename='physics', category=None, options={}):

	    systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)
            TTcut=''
            weight=1.0

            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

            if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

	    print("Diboson sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

            return sp


    class DibosonCF(Process):

        latexname = 'WW, W#gamma, ZZ#rightarrow ll#nu#nu (ChFlip only)'
        colour = kAzure - 9

        def base(self, treename='physics', category=None, options={}):

	    diboson = self.parent.procmap['Diboson'].base(treename, category, options)

            return diboson

        def __call__(self, treename='physics', category=None, options={}):

	    dibosoncf = self.parent.procmap['Diboson'](treename, category, options)

	    # plot only events where at least one lepton is charge flip. Remove req. where both lep must be prompt
	    #
            truthcut = category.cut.swapCut(self.vardb.getCut('2Lep_PurePromptEvent'),self.vardb.getCut('2Lep_ChFlipEvent'))

	    dibosoncf = dibosoncf.subprocess(cut=truthcut)

	    print("DibosonCF sp: {0}".format(dibosoncf.basecut.cutnamelist))

            return dibosoncf


    class HtoZZ(Process):

        latexname = 'H #rightarrow ZZ'
        colour = kTeal + 9

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                    ('HtoZZ', '125'),
                ]
            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) / 1000.
            return sp

        def __call__(self, treename='physics', category=None, options={}):
            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)
            TTcut=''
            weight=1.0

            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

            if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

	    print("HtoZZ sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

            return sp

    class Prompt(Process):

        latexname = 'Prompt'
        colour = kYellow - 9

        def base(self, treename='physics', category=None, options={}):

	    diboson = self.parent.procmap['Diboson'].base(treename, category, options)
            top = self.parent.procmap['Top'].base(treename, category, options)
            ttbarw = self.parent.procmap['TTBarW'].base(treename, category, options)
            ttbarz = self.parent.procmap['TTBarZ'].base(treename, category, options)
            return (diboson + top + ttbarw + ttbarz)

        def __call__(self, treename='physics', category=None, options={}):

            diboson = self.parent.procmap['Diboson'](treename, category, options)
            top = self.parent.procmap['Top'](treename, category, options)
            ttbarw = self.parent.procmap['TTBarW'](treename, category, options)
            ttbarz = self.parent.procmap['TTBarZ'](treename, category, options)

	    return (diboson + top + ttbarw + ttbarz)


    class ChargeFlipMC(Process):

	latexname = 'Charge flip (MC)'
        colour = kAzure -4

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                ('tops', 'tZ'),
                ('tops', 'tW'),
                ('tops', 'singlet'),
                ('tops', '4top'),
                ('tops', 'ttWW'),
                ('tops', 'ttbar_nonallhad'),
                ('Diboson','*'),
                ('Z+jetsCVetoBVeto', '*'),
                ('Z+jetsCFilterBVeto', '*'),
                ('Z+jetsBFilter', '*'),
                         ]
            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) / 1000.
            return sp

        def __call__(self, treename='physics', category=None, options={}):

	    systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)
            TTcut=''
            weight=1.0

            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'
            if TTcut != '':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

	    # plot only events where at least one lepton is charge flip. Remove req. where both lep must be prompt
	    #
            truthcut = category.cut.swapCut(self.vardb.getCut('2Lep_PurePromptEvent'),self.vardb.getCut('2Lep_ChFlipEvent'))

	    sp = sp.subprocess(cut=truthcut)

	    print("ChargeFlipMC sp: {0}".format(sp.basecut.cutnamelist))

            return sp

    """
    class ChargeFlipMC(Process):

	latexname = 'Charge flip (MC)'
        colour = kAzure - 4

        def base(self, treename='physics', category=None, options={}):

	    dibosoncf = self.parent.procmap['DibosonCF'].base(treename, category, options)
            topcf = self.parent.procmap['TopCF'].base(treename, category, options)
            zjetscf = self.parent.procmap['ZjetsCF'].base(treename, category, options)

	    return (dibosoncf + topcf + zjetscf)

        def __call__(self, treename='physics', category=None, options={}):

            dibosoncf = self.parent.procmap['DibosonCF'](treename, category, options)
            topcf = self.parent.procmap['TopCF'](treename, category, options)
            zjetscf = self.parent.procmap['ZjetsCF'](treename, category, options)

	    return (dibosoncf + topcf + zjetscf)

    """


    class FakesFF(Process):

	latexname = 'FakesFF'
        colour = kAzure - 9

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                    ('Data', 'physics_Main'),
                         ]
            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees)
            return sp

        def __call__(self, treename='physics', category=None, options={}):
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
            weight=1.0
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TLcut='TL'
                LTcut='LT'
                LLcut='LL'
                weight='FFWeight[0]'

	    print 'FakesFF - weight name : ', weight, ', TLcut : ', TLcut, ', LTcut : ', LTcut, ', LLcut : ', LLcut

            #Now doing prompt background subtractions from fakes. sublist must contain both prompt and chargeflips
	    #
            sublist = [ item for item in self.parent.sub_backgrounds ]
            for sample in sublist:
                sp = sp - self.parent.procmap[sample].base(treename, category, options) # here it is important to have used base otherwise at the sub sample there would be already applied the selection specified in the category __call__ of that sample i.e. for example also the iso-iso cut which is orthogonal to the cuts done in this sample which are fake iso or fake fake.
            #print sp
            sp_TL = sp.subprocess(cut=self.vardb.getCut(TLcut))
            sp_LT = sp.subprocess(cut=self.vardb.getCut(LTcut))
            sp_LL = sp.subprocess(cut=self.vardb.getCut(LLcut))

            #IN THIS VERSION OF THE CODE THE SIGN MINUS FOR FF EVENTS IS ALREADY IN THE WEIGHT SO THE SUBRACTION OF THE DOUBLE FAKE COUNTING IS DOE ADDING SP_FF
            sp = sp_TL + sp_LT +sp_LL
            #OLD WAY subtracting the double counting
            #sp = sp_if - sp_ff
            sp = sp.subprocess(cut=category.cut,eventweight=weight)
            #print sp
            return sp


    class FakesMM(Process):

	latexname = 'FakesMM'
        colour = kTeal - 9

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                    ('Data', 'physics_Main'),
                         ]
            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees)
            return sp

        def __call__(self, treename='physics', category=None, options={}):

	    systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)
            TTcut=''
            TLcut=''
            LTcut=''
            LLcut=''
            weight=1.0
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut  = 'TT'
                TLcut  = 'TL'
                LTcut  = 'LT'
                LLcut  = 'LL'
                weight = 'MMWeight[0]'

            # There is no prompt subtraction in the MM!

            sp_TT  = sp.subprocess(cut=category.cut & self.vardb.getCut(TTcut), eventweight=weight)
            sp_TL  = sp.subprocess(cut=category.cut & self.vardb.getCut(TLcut), eventweight=weight)
            sp_LT  = sp.subprocess(cut=category.cut & self.vardb.getCut(LTcut), eventweight=weight)
            sp_LL  = sp.subprocess(cut=category.cut & self.vardb.getCut(LLcut), eventweight=weight)

	    print("FakesMM - TT sp: {0}, weight: {1}".format(sp_TT.basecut.cutnamelist, weight))
	    print("FakesMM - TL sp: {0}, weight: {1}".format(sp_TL.basecut.cutnamelist, weight))
	    print("FakesMM - LT sp: {0}, weight: {1}".format(sp_LT.basecut.cutnamelist, weight))
	    print("FakesMM - LL sp: {0}, weight: {1}".format(sp_LL.basecut.cutnamelist, weight))

	    sp = sp_TT + sp_TL + sp_LT + sp_LL

            # THE ABOVE IS EQUIVALENT TO :
            #
            #sp = sp.subprocess(cut=category.cut,eventweight=weight)
            #
            # In fact, the MM weight for each event already contains the info on the event type: TT,TL,LT,LL :
            # no need to add any TT,LT ... cut!
            # (final estimate will be the sum of each contribution)
            #

            return sp


    class FakesABCD(Process):

        # Consider the following regions:
	#
        # Region A: TT, njet >= 4 (SR)
	# Region B (in MC): LT,TL njet >= 4
	# Region C (in DATA): TT, njet = [2,3]
	# Region D (in DATA): LT,TL  njet = [2,3]
	#
	# After subtracting prompt MC and charge flips from data in C and D, the estimate of fakes in A is:
	#
	# (C/D) * B = A
        #
	# where theta_e = (C/D)(ee), theta_mu = (C/D)(mumu)
	#
	# This corrects the MC normalisation (we assume the shape in B and A is the same!) for the fake probability as measured in data.
	#
	# Hence we return B (ee, mumu, OF), scaled by the ratio of TT/LT events (ee, mumu) in low jet multiplicity region.
	# The OF case is obtained via:
	#
	# theta_e * B(L_e,T_mu) + theta_mu * B (L_mu,T_e) = A (OF)

	latexname = 'Fakes #theta method'
        colour = kCyan - 9

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                    ('Data', 'physics_Main'),
                ]
            trees = self.inputs.getTrees(treename, inputgroup)
            return self.subprocess(trees=trees)


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

            print 'FakesABCD \n '

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
            weight=1.0
            if self.parent.channel=='TwoLepSS':
                TTCut  = self.vardb.getCut('TT')
                TLCut  = self.vardb.getCut('TL')
                LTCut  = self.vardb.getCut('LT')
                TelLmuCut = Cut('TelLmu',  '( isTelLmu == 1 )')
                LelTmuCut = Cut('LelTmu',  '( isLelTmu == 1 )')
                TmuLelCut = Cut('TmuLel',  '( isTmuLel == 1 )')
                LmuTelCut = Cut('LmuTel',  '( isLmuTel == 1 )')

            TL_LT_Cut = (TLCut | LTCut)

	    # take the base suprocess (DATA)
	    #
	    sp = self.base(treename, category, options)

	    print("base sp: {0}".format(category.cut.cutnamelist))

	    # Remove the cuts defining the flavour composition (this is made just for calculating thetas...)
	    #
            basecut = category.cut.removeCut(self.vardb.getCut('ElEl_Event'))
            basecut = basecut.removeCut(self.vardb.getCut('MuMu_Event'))
            basecut = basecut.removeCut(self.vardb.getCut('OF_Event'))


	    if TTHBackgrounds2015.theta['El'][0] == 999.0 :

                print ("Calculating theta_el from data in regions C,D...")

                # define selection for region C and D (ee)
                #
                cut_sp_C_el = basecut.swapCut(self.vardb.getCut('NJet2L'),self.vardb.getCut('LowJetCR')) & self.vardb.getCut('ElEl_Event') & TTCut & self.vardb.getCut('Zsidescut') & self.vardb.getCut('ElEtaCut')
                cut_sp_D_el = basecut.swapCut(self.vardb.getCut('NJet2L'),self.vardb.getCut('LowJetCR')) & self.vardb.getCut('ElEl_Event') & TL_LT_Cut & self.vardb.getCut('Zsidescut') & self.vardb.getCut('ElEtaCut')

                sp_C_el = sp.subprocess( cut = cut_sp_C_el, eventweight=weight )
                sp_D_el = sp.subprocess( cut = cut_sp_D_el, eventweight=weight )

                print("Region C (el) sp: {0}".format(cut_sp_C_el.cutnamelist))
                print("Region D (el) sp: {0}".format(cut_sp_D_el.cutnamelist))

                # get a list of stuff to subtract for region C and D (i.e, prompt MC and charge flips)
                #
                sublist = [ item for item in self.parent.sub_backgrounds ]

                # ... and now subtract!
                #
                # NB: here it is crucial to call .base() on the subprocess, otherwise the subprocess would have the cuts
                # defined in its own __call__ method already applied, whcih in general is not what we want
                # (e.g., it might have a TT selection applied, when we want to consider TL events instead...)
                #
                for sample in sublist:

                    print ("Subtracting {0} from data in regions C,D...".format(sample))

                    this_cut_sp_C_el = cut_sp_C_el
                    sub_sample_C_el = self.parent.procmap[sample].base(treename,category,options)

                    if sample == "ChargeFlipMC":
                        this_cut_sp_C_el = this_cut_sp_C_el.swapCut(self.vardb.getCut('2Lep_PurePromptEvent'),self.vardb.getCut('2Lep_ChFlipEvent'))

                    sub_sample_C_el = sub_sample_C_el.subprocess( cut = this_cut_sp_C_el, eventweight=weight )

                    #print ("sub_sample_C_el ={0} ".format(sub_sample_C_el.basecut.cutnamelist))

                    #print ("C (el) - yields data: ", sp_C_el.numberstats())
                    #print ("C (el) - yields MC: ", sub_sample_C_el.numberstats())
                    sp_C_el = sp_C_el - sub_sample_C_el
                    #print ("C (el) - yields after sub: ", sp_C_el.numberstats())

                    # *********************************************

                    this_cut_sp_D_el = cut_sp_D_el
                    sub_sample_D_el = self.parent.procmap[sample].base(treename,category,options)

                    if sample == "ChargeFlipMC":
                        this_cut_sp_D_el = this_cut_sp_D_el.swapCut(self.vardb.getCut('2Lep_PurePromptEvent'),self.vardb.getCut('2Lep_ChFlipEvent'))

                    sub_sample_D_el = sub_sample_D_el.subprocess( cut = this_cut_sp_D_el, eventweight='weight_lepton_trig_HTop[0] * weight_lepton_reco_HTop[0] * weight_lepton_iso_HTop[0] * weight_lepton_ID_HTop[0] * weight_lepton_TTVA_HTop[0] * weight_jet__MV2c20_SFFix77[0]' )

                    #print ("sub_sample_D_el ={0} ".format(sub_sample_D_el.basecut.cutnamelist))

                    #print ("D (el) - yields data: ", sp_D_el.numberstats())
                    #print ("D (el) - yields MC: ", sub_sample_D_el.numberstats())
                    sp_D_el = sp_D_el - sub_sample_D_el
                    #print ("D (el) - yields after sub: ", sp_D_el.numberstats())

                print ("---------------------------------------------------------------------\n")
		print ("C (el) - data yields after prompt/ch-flip subtraction: ", sp_C_el.numberstats())
                print ("D (el) - data yields after prompt/ch-flip subtraction: ", sp_D_el.numberstats())
                print ("\n---------------------------------------------------------------------")

                # derive theta factors for el
                #
                TTHBackgrounds2015.theta['El'] = self.calcTheta(sp_C_el,sp_D_el,stream='El')

	    else :
                print ("Reading theta(el) value: {0} +- {1}".format(TTHBackgrounds2015.theta['El'][0], TTHBackgrounds2015.theta['El'][1]))


	    if TTHBackgrounds2015.theta['Mu'][0] == 999.0 :

                print ("Calculating theta_mu from data in regions C,D...")

                # define selection for region C and D (mumu)
                #
                cut_sp_C_mu = basecut.swapCut(self.vardb.getCut('NJet2L'),self.vardb.getCut('LowJetCR')) & self.vardb.getCut('MuMu_Event') & TTCut
                cut_sp_D_mu = basecut.swapCut(self.vardb.getCut('NJet2L'),self.vardb.getCut('LowJetCR')) & self.vardb.getCut('MuMu_Event') & TL_LT_Cut

                sp_C_mu = sp.subprocess( cut = cut_sp_C_mu, eventweight=weight )
                sp_D_mu = sp.subprocess( cut = cut_sp_D_mu, eventweight=weight )

                print("Region C (mu) sp: {0}".format(cut_sp_C_mu.cutnamelist))
                print("Region D (mu) sp: {0}".format(cut_sp_D_mu.cutnamelist))

                # get a list of stuff to subtract for region C and D (i.e, prompt MC and charge flips)
                #
                sublist = [ item for item in self.parent.sub_backgrounds ]

                # ... and now subtract!
                #
                # NB: here it is crucial to call .base() on the subprocess, otherwise the subprocess would have the cuts
                # defined in its own __call__ method already applied, whcih in general is not what we want
                # (e.g., it might have a TT selection applied, when we want to consider TL events instead...)
                #
                for sample in sublist:

                    print ("Subtracting {0} from data in regions C,D...".format(sample))

                    this_cut_sp_C_mu = cut_sp_C_mu
                    sub_sample_C_mu = self.parent.procmap[sample].base(treename,category,options)

                    if sample == "ChargeFlipMC":
                        this_cut_sp_C_mu = this_cut_sp_C_mu.swapCut(self.vardb.getCut('2Lep_PurePromptEvent'),self.vardb.getCut('2Lep_ChFlipEvent'))

                    sub_sample_C_mu = sub_sample_C_mu.subprocess( cut = this_cut_sp_C_mu, eventweight=weight )

                    #print ("sub_sample_C_mu ={0} ".format(sub_sample_C_mu.basecut.cutnamelist))

                    #print ("C (mu) - yields data: ", sp_C_mu.numberstats())
                    #print ("C (mu) - yields MC: ", sub_sample_C_mu.numberstats())
                    sp_C_mu = sp_C_mu - sub_sample_C_mu
                    #print ("C (mu) - yields after sub: ", sp_C_mu.numberstats())

                    # *********************************************

                    this_cut_sp_D_mu = cut_sp_D_mu
                    sub_sample_D_mu = self.parent.procmap[sample].base(treename,category,options)

                    if sample == "ChargeFlipMC":
                        this_cut_sp_D_mu = this_cut_sp_D_mu.swapCut(self.vardb.getCut('2Lep_PurePromptEvent'),self.vardb.getCut('2Lep_ChFlipEvent'))

                    sub_sample_D_mu = sub_sample_D_mu.subprocess( cut = this_cut_sp_D_mu, eventweight='weight_lepton_trig_HTop[0] * weight_lepton_reco_HTop[0] * weight_lepton_iso_HTop[0] * weight_lepton_ID_HTop[0] * weight_lepton_TTVA_HTop[0] * weight_jet__MV2c20_SFFix77[0]' )

                    #print ("sub_sample_D_mu ={0} ".format(sub_sample_D_mu.basecut.cutnamelist))

                    #print ("D (mu) - yields data: ", sp_D_mu.numberstats())
                    #print ("D (mu) - yields MC: ", sub_sample_D_mu.numberstats())
                    sp_D_mu = sp_D_mu - sub_sample_D_mu
                    #print ("D (mu) - yields after sub: ", sp_D_mu.numberstats())


                print ("---------------------------------------------------------------------\n")
		print ("C (mu) - data yields after prompt/ch-flip subtraction: ", sp_C_mu.numberstats())
                print ("D (mu) - data yields after prompt/ch-flip subtraction: ", sp_D_mu.numberstats())
                print ("\n---------------------------------------------------------------------")

                # derive theta factors for mu
                #
                TTHBackgrounds2015.theta['Mu'] = self.calcTheta(sp_C_mu,sp_D_mu,stream='Mu')

	    else :
                print ("Reading theta(mu) value: {0} +- {1}".format(TTHBackgrounds2015.theta['Mu'][0], TTHBackgrounds2015.theta['Mu'][1]))


	    # define region B,  depending on which flavour composition we are looking at
	    #
            # Take TTbar MC, excluding all prompts and charge flips, and reweight it by the theta factors measured in data
            #
	    cut_sp_B_SF     = category.cut.swapCut(self.vardb.getCut('2Lep_PurePromptEvent'),self.vardb.getCut('2Lep_NonPromptEvent'))
            cut_sp_B_SF     = cut_sp_B_SF.swapCut(self.vardb.getCut('LowJetCR'),self.vardb.getCut('NJet2L')) & TL_LT_Cut
	    cut_sp_B_OF_Lel = category.cut.swapCut(self.vardb.getCut('2Lep_PurePromptEvent'),self.vardb.getCut('2Lep_NonPromptEvent'))
	    cut_sp_B_OF_Lel = cut_sp_B_OF_Lel.swapCut(self.vardb.getCut('LowJetCR'),self.vardb.getCut('NJet2L')) & (LelTmuCut | TmuLelCut)
            cut_sp_B_OF_Lmu = category.cut.swapCut(self.vardb.getCut('2Lep_PurePromptEvent'),self.vardb.getCut('2Lep_NonPromptEvent'))
	    cut_sp_B_OF_Lmu = cut_sp_B_OF_Lmu.swapCut(self.vardb.getCut('LowJetCR'),self.vardb.getCut('NJet2L')) & (TelLmuCut | LmuTelCut)

	    if not ("OF_Event") in category.cut.cutname:
                sp_B = self.parent.procmap['TTBar'].base(treename,category,options)
                sp_B = sp_B.subprocess(cut=cut_sp_B_SF,eventweight='weight_lepton_trig_HTop[0] * weight_lepton_reco_HTop[0] * weight_lepton_iso_HTop[0] * weight_lepton_ID_HTop[0] * weight_lepton_TTVA_HTop[0] * weight_jet__MV2c20_SFFix77[0]')
                if ("ElEl_Event") in category.cut.cutname:
                    sp_B = self.applyThetaFactor(sp_B,TTHBackgrounds2015.theta['El'])
                elif ("MuMu_Event") in category.cut.cutname:
                    sp_B = self.applyThetaFactor(sp_B,TTHBackgrounds2015.theta['Mu'])
            else:
                sp_B_Lel = self.parent.procmap['TTBar'].base(treename,category,options)
                sp_B_Lmu = self.parent.procmap['TTBar'].base(treename,category,options)
                sp_B_Lel = sp_B_Lel.subprocess(cut=cut_sp_B_OF_Lel, eventweight='weight_lepton_trig_HTop[0] * weight_lepton_reco_HTop[0] * weight_lepton_iso_HTop[0] * weight_lepton_ID_HTop[0] * weight_lepton_TTVA_HTop[0] * weight_jet__MV2c20_SFFix77[0]')
                sp_B_Lmu = sp_B_Lmu.subprocess(cut=cut_sp_B_OF_Lmu, eventweight='weight_lepton_trig_HTop[0] * weight_lepton_reco_HTop[0] * weight_lepton_iso_HTop[0] * weight_lepton_ID_HTop[0] * weight_lepton_TTVA_HTop[0] * weight_jet__MV2c20_SFFix77[0]')
                sp_B = self.applyThetaFactor(sp_B_Lmu,TTHBackgrounds2015.theta['Mu']) + self.applyThetaFactor(sp_B_Lel,TTHBackgrounds2015.theta['El'])

	    print ("=================>\n")
	    print ("Region B sp: {0} \n".format(sp_B.basecut.cutnamelist))
            print ("fakes yield: ", sp_B.numberstats() ,"\n")

            return sp_B


    # Check Data vs. TTBar MC in low njet, TL,LT, e-mu region
    # Use the theta factors derived in data to rewight TTBar MC
    #
    class FakesClosureDataABCD(Process):

	latexname = 'Fakes #theta method (t\bar{t} reweighted)'
        colour = kCyan - 9

        theta_el = 0.1467
        theta_mu = 0.6856

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
		    ('tops', 'ttbar_nonallhad'),
                ]
            trees = self.inputs.getTrees(treename, inputgroup)
            return self.subprocess(trees=trees) / 1000

        def __call__(self, treename='physics', category=None, options={}):

            print 'FakesClosureDataABCD \n '

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
            weight=1.0
            if self.parent.channel=='TwoLepSS':
                TTCut  = self.vardb.getCut('TT')
                TLCut  = self.vardb.getCut('TL')
                LTCut  = self.vardb.getCut('LT')
                TelLmuCut = Cut('TelLmu',  '( isTelLmu == 1 )')
                LelTmuCut = Cut('LelTmu',  '( isLelTmu == 1 )')
                TmuLelCut = Cut('TmuLel',  '( isTmuLel == 1 )')
                LmuTelCut = Cut('LmuTel',  '( isLmuTel == 1 )')

            TL_LT_Cut = (TLCut | LTCut)

	    # take the base suprocess (i.e, TTBar)
	    #
	    sp = self.base(treename, category, options)

	    print("base sp: {0}".format(category.cut.cutnamelist))

	    # define closure region (OF, Tl,LT, low njet)
	    #
            # Veto all prompt and charge flips (for TTBar, basically it's a charge flip veto)
	    cut_sp_OF_Lel = category.cut & (LelTmuCut | TmuLelCut)
	    cut_sp_OF_Lmu = category.cut & (TelLmuCut | LmuTelCut)
            cut_sp_OF_Lel = cut_sp_OF_Lel.swapCut(self.vardb.getCut('2Lep_PurePromptEvent'),self.vardb.getCut('2Lep_NonPromptEvent'))
	    cut_sp_OF_Lmu = cut_sp_OF_Lmu.swapCut(self.vardb.getCut('2Lep_PurePromptEvent'),self.vardb.getCut('2Lep_NonPromptEvent'))

            # plug in the theta factors by hand...CHANGE ME!
            #
            sp_Lel = sp.subprocess(cut=cut_sp_OF_Lel, eventweight='weight_lepton_trig_HTop[0] * weight_lepton_reco_HTop[0] * weight_lepton_iso_HTop[0] * weight_lepton_ID_HTop[0] * weight_lepton_TTVA_HTop[0] * weight_jet__MV2c20_SFFix77[0]')
            sp_Lmu = sp.subprocess(cut=cut_sp_OF_Lmu, eventweight='weight_lepton_trig_HTop[0] * weight_lepton_reco_HTop[0] * weight_lepton_iso_HTop[0] * weight_lepton_ID_HTop[0] * weight_lepton_TTVA_HTop[0] * weight_jet__MV2c20_SFFix77[0]')
            sp_final = ( sp_Lel * self.theta_el ) + ( sp_Lmu * self.theta_mu )

	    print ("=================>\n")
	    print ("Region closure sp: {0}".format(sp_final.basecut.cutnamelist))
            print ("yield: ", sp_final.numberstats())

            return sp_final

    class FakesClosureMM(Process):

	latexname = 'FakesMM - t#bar{t}'
	colour = kTeal -9

        name = 'FakesClosureMM'

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
		    ('tops', 'ttbar_nonallhad'),
		    #('tops', 'ttbar_dilep'),
                         ]
            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) / 1000.
            return sp

        def __call__(self, treename='physics', category=None, options={}):

            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)
            TTcut=''
            TLcut=''
            LTcut=''
            LLcut=''
            weight=1.0
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut  = 'TT'
                TLcut  = 'TL'
                LTcut  = 'LT'
                LLcut  = 'LL'
                weight = 'MMWeight[0]'

	    # plot only events where at least one lepton is non-prompt, and none of the leptons is charge flip
	    # ---> NB: truth req. *only* for closure test!
	    #
            sp_TT  = sp.subprocess(cut=category.cut & self.vardb.getCut('2Lep_NonPromptEvent') & self.vardb.getCut(TTcut), eventweight=weight)
            sp_TL  = sp.subprocess(cut=category.cut & self.vardb.getCut('2Lep_NonPromptEvent') & self.vardb.getCut(TLcut), eventweight=weight)
            sp_LT  = sp.subprocess(cut=category.cut & self.vardb.getCut('2Lep_NonPromptEvent') & self.vardb.getCut(LTcut), eventweight=weight)
            sp_LL  = sp.subprocess(cut=category.cut & self.vardb.getCut('2Lep_NonPromptEvent') & self.vardb.getCut(LLcut), eventweight=weight)

	    print("FakesClosureMM - TT sp: {0}, weight: {1}".format(sp_TT.basecut.cutnamelist, weight))
	    print("FakesClosureMM - TL sp: {0}, weight: {1}".format(sp_TL.basecut.cutnamelist, weight))
	    print("FakesClosureMM - LT sp: {0}, weight: {1}".format(sp_LT.basecut.cutnamelist, weight))
	    print("FakesClosureMM - LL sp: {0}, weight: {1}".format(sp_LL.basecut.cutnamelist, weight))

            sp = sp_TT + sp_TL + sp_LT + sp_LL

            return sp

    class FakesClosureABCD(Process):

	latexname = 'FakesABCD - t#bar{t}'
        colour = kCyan - 9

        theta_el = {
            'El': (1.0, 0.0),
        }
        theta_mu = {
            'Mu': (1.0, 0.0),
        }

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
		    ('tops', 'ttbar_nonallhad'),
		    #('tops', 'ttbar_dilep'),
                         ]
            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) / 1000.
            return sp

        def calcTheta(self, sp, cut_num, cut_denom, stream=None, options={}):

	    numerator   = sp.subprocess(cut=cut_num)
	    denominator = sp.subprocess(cut=cut_denom)

            if not denominator.number():
                print "ERROR: Cannot calculate theta transfer factor!"

            theta = (numerator/denominator).numberstats()
            print "Calculated theta transfer factor for stream ", stream, " - theta :", theta[0], '+-', theta[1]

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

            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)
            TTcut =''
            TTcut =''
            TLcut =''
            LTcut =''
            LLcut =''
            TelLmucut =''
            LelTmucut =''
            TmuLelcut =''
            LmuTelcut =''
            weight=1.0
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut  = self.vardb.getCut('TT')
                TLcut  = self.vardb.getCut('TL')
                LTcut  = self.vardb.getCut('LT')
                LLcut  = self.vardb.getCut('LL')
                TelLmucut = self.vardb.getCut('OF_Event') & Cut('TelLmu',  '( isTelLmu == 1 )')
                LelTmucut = self.vardb.getCut('OF_Event') & Cut('LelTmu',  '( isLelTmu == 1 )')
                TmuLelcut = self.vardb.getCut('OF_Event') & Cut('TmuLel',  '( isTmuLel == 1 )')
                LmuTelcut = self.vardb.getCut('OF_Event') & Cut('LmuTel',  '( isLmuTel == 1 )')

            print 'FakesClosureABCD \n '

	    # define a base subprocess to derive theta factors
	    #
	    # --> low njet region, do not specify for base sp neither the flavour composition, nor the tightness combination
	    #     (but require truth cut)
	    #
	    basecut            = self.vardb.getCuts(['NBJet', '2Lep', 'SS', 'TrigMatch', 'LowJetCR', '2Lep_NonPromptEvent'])
	    #
	    # now characterise!
	    #
	    numerator_cut_el   = basecut & self.vardb.getCut('ElEl_Event') & TTcut
	    numerator_cut_mu   = basecut & self.vardb.getCut('MuMu_Event') & TTcut
	    denominator_cut_el = basecut & self.vardb.getCut('ElEl_Event') & ( TLcut | LTcut )
	    denominator_cut_mu = basecut & self.vardb.getCut('MuMu_Event') & ( TLcut | LTcut )
	    #
	    theta_sp = self.base(treename) / 1000

	    # update default dictionaries for theta (do it only once to change the default value!)
	    #

	    #print '(self.theta_el[\'El\'])[0] = ', (self.theta_el['El'])[0], ' - (self.theta_el[\'El\'])[1] = ', (self.theta_el['El'])[1]
	    #print '(self.theta_mu[\'Mu\'])[0] = ', (self.theta_mu['Mu'])[0], ' - (self.theta_mu[\'Mu\'])[1] = ', (self.theta_mu['Mu'])[1]

	    #if ( (self.theta_el['El'])[0] == 1.0 and (self.theta_el['El'])[1] == 0.0 ):
	    self.theta_el = self.calcTheta(theta_sp,numerator_cut_el,denominator_cut_el,stream='El')
	    #if ( (self.theta_mu['Mu'])[0] == 1.0 and (self.theta_mu['Mu'])[1] == 0.0 ):
	    self.theta_mu = self.calcTheta(theta_sp,numerator_cut_mu,denominator_cut_mu,stream='Mu')

	    # look only at events where at least one lepton is !prompt, and none is charge flip
	    #
            sp = sp.subprocess(cut=category.cut & self.vardb.getCut('2Lep_NonPromptEvent'), eventweight=weight)
            print ' \ttruth cut : ', (self.vardb.getCut('2Lep_NonPromptEvent')).cutstr

	    if ('NJet2L') in category.cut.cutname:

	      if ("ElEl_Event") in category.cut.cutname:

		sp_elel = sp.subprocess(cut=category.cut & ( TLcut | LTcut ), eventweight=weight)
	        sp_elel = self.applyThetaFactor(sp_elel,self.theta_el)
		return sp_elel

	      elif ("MuMu_Event") in category.cut.cutname:

		sp_mumu = sp.subprocess(cut=category.cut & ( TLcut | LTcut ), eventweight=weight)
		sp_mumu = self.applyThetaFactor(sp_mumu,self.theta_mu)
	        return sp_mumu

	      elif ("OF_Event") in category.cut.cutname:

	        sp_OF_loose_is_el = sp.subprocess(cut=category.cut & ( LelTmucut | TmuLelcut ), eventweight=weight)
		sp_OF_loose_is_mu = sp.subprocess(cut=category.cut & ( TelLmucut | LmuTelcut ), eventweight=weight)
		sp_OF = self.applyThetaFactor(sp_OF_loose_is_mu,self.theta_mu) + self.applyThetaFactor(sp_OF_loose_is_el,self.theta_el)
		return sp_OF

	      else:
                 print ' \tcould not find ( \'ElEl_Event\' || \'MuMu_Event\' || \'OF_Event\' ) in category name - just doing standard TTBarClosureMM'
	         return sp

            print ' \tcould not find \'NJet2L\' in category name - just doing standard TTBarClosureMM'
	    return sp

'''
    class Prompt2(Process):
        latexname = 'Prompt2'
        colour = kYellow - 9

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                    ('Prompt', 'WZ'),
                    ('Prompt', 'ZZ'),
                    #('Prompt', 'WWW'),
                    #('Prompt', 'ZWW'),
                    #('Prompt', 'ZZZ'),
                    ('Prompt', 'ttW'),
                    ('Prompt', 'ttZ'),
                    ('Prompt', 'WWjj'),
                ]
            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) / 1000.
            return sp

        def __call__(self, treename='physics', category=None, options={}):
            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            kf_opts = {}#this dictionary tell if you have to apply the syst and in which direction
            if systematics and systematics.name == 'Prompt_kf':
                kf_opts['systematics'] = True
                kf_opts['systematicsdirection'] = direction
            kfactors = {'Norm_corr': (1.0, 0.1),}#this dictionary contain the value of the kfactor and its error

            sp = self.base(treename, category, options)
            TTcut=''
            weight=1.0
            if self.parent.channel=='TwoLepSS':
                TTcut='TT'
                kfactors = {'Norm_corr': (1.0, 0.1),}
            elif self.parent.channel=='ThreeLep':
                TTcut='TT'
                kfactors = {'Norm_corr': (1.0, 0.1),}
            sp = self.parent.applyKfactor(sp, category, kfactors, kf_opts)#this command is used to apply the kfactor
            if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))
            sp = sp.subprocess(cut=category.cut,eventweight=weight)
            return sp

    class DPI(Process):
        latexname = 'DPI'
        colour = kRed

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                    ('Prompt', 'DPIWW'),
                    ('Prompt', 'DPIWZ'),
                    ('Prompt', 'DPIZZ'),
                ]
            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) / 1000.
            return sp

        def __call__(self, treename='physics', category=None, options={}):
            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)
            TTcut=''
            weight=1.0
            if self.parent.channel=='TwoLepSS':
                TTcut='TT'
            elif self.parent.channel=='ThreeLep':
                TTcut='TT'
            if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))
            sp = sp.subprocess(cut=category.cut,eventweight=weight)
            return sp
'''
