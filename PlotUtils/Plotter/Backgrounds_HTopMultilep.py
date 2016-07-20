""" Backgrounds_HTopMultilep.py: a set of classes to manage the various processes """

__author__     = "KG, Marco Milesi, Francesco Nuti"
__email__      = "Kong.Guan.Tan@cern.ch, marco.milesi@cern.ch, francesco.nuti@cern.ch"
__maintainer__ = "Marco Milesi"

import os, sys, math, types

sys.path.append(os.path.abspath(os.path.curdir))

from Plotter.BackgroundTools_HTopMultilep import loadSamples, drawText, Category, Background, Process, VariableDB, Variable, Cut, Systematics, Category

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

class TTHBackgrounds(Background):

    backgrounds     = [ 'TTBarW', 'TTBarZ', 'Top', 'TopCF', 'Diboson', 'DibosonCF', 'HtoZZ', 'ZjetsLF', 'Zjets', 'Wjets', 'Prompt', 'ChargeFlipMC', 'ChargeFlip', 'FakesFF', 'FakesMM', 'FakesABCD']
    sub_backgrounds = [ 'TTBarW', 'TTBarZ', 'Top', 'TTBar', 'Diboson', 'Zjets', 'Wjets', 'ChargeFlipMC']
    signals     = ['TTBarH']
    observed    = ['Observed']
    luminosity  = 1.0
    lumi_units  = 'fb-1'
    norm_factor = 1.0           # Might be needed to correct units for Xsec weight
    rescaleXsecAndLumi = False  # Set to "True" if you don't want to take into account the Xsec*lumi weight
    channel     = 'TwoLepSS'    # Can be one of ['TwoLepSS' 'TwoLepCR', 'ThreeLep', 'FourLep']
    eventweight = 1.0
    #eventweight = 'weight_pileup*weight_muon_trig*weight_electron_trig'
    #eventweight = 'evtsel_weight*evtsel_weight_el*evtsel_weight_mu*evtsel_bjet_weight*evtsel_weight_lep_trigger*weight_CF'
    #eventweight = 'evtsel_weight*evtsel_weight_el*evtsel_weight_mu*evtsel_bjet_weight*evtsel_weight_lep_trigger*weight_MM*weight_CF'
    #eventweight = 'evtsel_weight*evtsel_weight_el*evtsel_weight_mu*evtsel_bjet_weight*evtsel_weight_lep_trigger*weight_FF*weight_CF'
    useEmbedding    = False
    useZCorrections = False
    useSherpaNNPDF30NNLO = False
    RQCD = {
        #'El': (1.000, 0.051),
        #'Mu': (1.177, 0.066),
        'El': (1.00, 0.05),
        'Mu': (1.11, 0.08),
    }

    theta = {
        'El': (999.0, 0.0),
        'Mu': (999.0, 0.0),
        #'El': (0.163, 0.036), # v027
        #'Mu': (0.159, 0.022), # v027
        #'El': (0.390, 0.156), # v15
        #'Mu': (0.602, 0.120), # v15
    }

    theta_MC = {
        'El': (999.0, 0.0),
        'Mu': (999.0, 0.0),
        #'El': (0.194, 0.037), # v027
        #'Mu': (0.153, 0.022), # v027
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

        latexname = 'Data'

        def base(self, treename='physics', category=None, options={}):
            #Contains the instuction of which tree load and eventually correct the units of the cross-setion (see the division /1000.). Note it is not automatically executed when the class is called. It is executed trought __call__
            inputgroup = [
                    ('Data', 'physics_Main'),
                ]

            print inputgroup

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

            print("\nObserved sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

            return sp


    class TTBarH(Process):

        latexname = 't#bar{t} H'
        colour = kBlack

        def base(self, treename='physics', category=None, options={}):

            #hmass = options.get('hmass', '125') #300 is the default value if hmass is not in options, hmass is specified in one option passed to the plot function in the plotting script
            inputgroup = [
                #('ttH', 'ttH_dil'),
                #('ttH', 'ttH_semilep'),
                #('ttH', 'ttH_allhad'),
                ('ttH', 'ttH_dil_Pythia8'),
                ('ttH', 'ttH_semilep_Pythia8'),
                ('ttH', 'ttH_allhad_Pythia8'),
                        ]

            print inputgroup

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
            TTcut=''
            weight=1.0
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

            if TTcut is not '':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\nTTBarH sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

            return sp

    class TTBarHDilep(Process):

        latexname = 't#bar{t} H (dilep)'
        colour = kBlack

        def base(self, treename='physics', category=None, options={}):

            #hmass = options.get('hmass', '125') #300 is the default value if hmass is not in options, hmass is specified in one option passed to the plot function in the plotting script
            inputgroup = [
                #('ttH', 'ttH_dil'),
                ('ttH', 'ttH_dil_Pythia8'),
                         ]

            print inputgroup

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
            TTcut=''
            weight=1.0
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

            if TTcut is not '':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\nTTBarHDilep sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

            return sp

    class TTBarW(Process):

        latexname = 't#bar{t} W'
        colour = kRed - 4

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                #('tops', 'ttW'),
		('tops', 'ttW_aMcAtNlo'),
                #('tops', 'Sherpa_ttW'),
                         ]

            print inputgroup

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
            TTcut=''
            weight=1.0
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

            if TTcut is not '':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\nTTBarW sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

            return sp


    class TTBarZ(Process):

        latexname = 't#bar{t} Z'
        colour = kRed - 7

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                #('tops', 'ttZnnqq'),
                #('tops', 'ttee'),
                #('tops', 'ttmumu'),
                #('tops', 'tttautau'),
                ('tops', 'ttee_aMcAtNlo'),
                ('tops', 'ttmumu_aMcAtNlo'),
                ('tops', 'tttautau_aMcAtNlo'),
                #('tops', 'Sherpa_ttZnnqq'),
                #('tops', 'Sherpa_ttll'),
                         ]

            print inputgroup

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
            TTcut=''
            weight=1.0
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

            if TTcut is not '':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\nTTBarZ sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

            return sp


    class Zeejets(Process):

        latexname = 'Z/#gamma*#rightarrow#it{ee}+jets'
        colour = kGreen-7

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                           #('Z+jets', 'ee'),
                           #('MadGraphZ+jets', 'ee'),
                           ('Z+jetsLowMllBVeto', 'ee'),
                           ('Z+jetsLowMllBFilter', 'ee'),
                           #('DYZ+jets', 'ee'),
                         ]

            if self.parent.useSherpaNNPDF30NNLO:
                # Sherpa NNPDF30NNLO
                #
                print("\nUsing Sherpa NNPDF30NNLO Z+jets! Jet reweighting needed...\n")
                inputgroup += [ ('Z+jetsCVetoBVeto_NNPDF30NNLO', 'ee'), ('Z+jetsCFilterBVeto_NNPDF30NNLO', 'ee'), ('Z+jetsBFilter_NNPDF30NNLO', 'ee') ]
            else:
                # Sherpa CT10
                #
                inputgroup += [ ('Z+jetsCVetoBVeto', 'ee'), ('Z+jetsCFilterBVeto', 'ee'), ('Z+jetsBFilter', 'ee') ]

            print inputgroup

            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) * self.parent.norm_factor
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

            if self.parent.useSherpaNNPDF30NNLO:
                weight = 'SherpaNJetWeight'

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

            if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\nZeejets sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

            return sp


    class Zmumujets(Process):

        latexname = 'Z/#gamma*#rightarrow#mu#mu+jets'
        colour = kTeal+2

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                           #('Z+jets', 'mumu'),
                           #('MadGraphZ+jets', 'mumu'),
                           ('Z+jetsLowMllBVeto', 'mumu'),
                           ('Z+jetsLowMllBFilter', 'mumu'),
                           #('DYZ+jets', 'mumu'),
                         ]

            if self.parent.useSherpaNNPDF30NNLO:
                # Sherpa NNPDF30NNLO
                #
                print("\nUsing Sherpa NNPDF30NNLO Z+jets! Jet reweighting needed...\n")
                inputgroup += [ ('Z+jetsCVetoBVeto_NNPDF30NNLO', 'mumu'), ('Z+jetsCFilterBVeto_NNPDF30NNLO', 'mumu'), ('Z+jetsBFilter_NNPDF30NNLO', 'mumu') ]
            else:
                # Sherpa CT10
                #
                inputgroup += [ ('Z+jetsCVetoBVeto', 'mumu'), ('Z+jetsCFilterBVeto', 'mumu'), ('Z+jetsBFilter', 'mumu') ]

            print inputgroup

            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) * self.parent.norm_factor
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

            if self.parent.useSherpaNNPDF30NNLO:
                weight = 'SherpaNJetWeight'

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

            if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\nZmumujets sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

            return sp

    class Ztautaujets(Process):

        latexname = 'Z/#gamma*#rightarrow#tau#tau+jets'
        colour = kTeal

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                           #('Z+jets', 'tautau'),
                           #('MadGraphZ+jets', 'tautau'),
                           ('Z+jetsLowMllBVeto', 'tautau'),
                           ('Z+jetsLowMllBFilter', 'tautau'),
                           #('DYZ+jets', 'tautau'),
                         ]

            if self.parent.useSherpaNNPDF30NNLO:
                # Sherpa NNPDF30NNLO
                #
                print("\nUsing Sherpa NNPDF30NNLO Z+jets! Jet reweighting needed...\n")
                inputgroup += [ ('Z+jetsCVetoBVeto_NNPDF30NNLO', 'tautau'), ('Z+jetsCFilterBVeto_NNPDF30NNLO', 'tautau'), ('Z+jetsBFilter_NNPDF30NNLO', 'tautau') ]
            else:
                # Sherpa CT10
                #
                inputgroup += [ ('Z+jetsCVetoBVeto', 'tautau'), ('Z+jetsCFilterBVeto', 'tautau'), ('Z+jetsBFilter', 'tautau') ]

            print inputgroup

            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) * self.parent.norm_factor
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

            if self.parent.useSherpaNNPDF30NNLO:
                weight = 'SherpaNJetWeight'

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

            if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\nZtautaujets sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

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
                           ('Z+jetsLowMllBVeto', '*'),
                         ]
            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) * self.parent.norm_factor
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

            print("\nZjetsLF sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

            return sp

    class ZjetsHF(Process):

        latexname = 'Z/#gamma* + HF jets'
        colour = kGreen

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                           ('Z+jetsCFilterBVeto', '*'),
                           ('Z+jetsBFilter', '*'),
                           ('Z+jetsLowMllBFilter', '*'),
                         ]
            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) * self.parent.norm_factor
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

            print("\nZjetsLF sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

            return sp


    class ZjetsCF(Process):

        latexname = 'Z/#gamma*+jets (QMisID only)'
        colour = kAzure + 10

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                           #('Z+jets', '*'),
                           #('MadGraphZ+jets', '*'),
                           ('Z+jetsLowMllBVeto', '*'),
                           ('Z+jetsLowMllBFilter', '*'),
                           #('DYZ+jets', '*'),
                         ]

            if self.parent.useSherpaNNPDF30NNLO:
                # Sherpa NNPDF30NNLO
                #
                print("\nUsing Sherpa NNPDF30NNLO Z+jets! Jet reweighting needed...\n")
                inputgroup += [ ('Z+jetsCVetoBVeto_NNPDF30NNLO', '*'), ('Z+jetsCFilterBVeto_NNPDF30NNLO', '*'), ('Z+jetsBFilter_NNPDF30NNLO', '*') ]
            else:
                # Sherpa CT10
                #
                inputgroup += [ ('Z+jetsCVetoBVeto', '*'), ('Z+jetsCFilterBVeto', '*'), ('Z+jetsBFilter', '*') ]

            print inputgroup

            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) * self.parent.norm_factor
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

            if self.parent.useSherpaNNPDF30NNLO:
                weight = 'SherpaNJetWeight'

            # plot only events where at least one lepton is charge flip. Remove req. where both lep must be prompt
            #
            #truthcut = category.cut.swapCut(self.vardb.getCut('2Lep_PurePromptEvent'),self.vardb.getCut('2Lep_QMisIDEvent'))
            truthcut = category.cut.removeCut(self.vardb.getCut('2Lep_ProbePromptEvent'))
            truthcut = truthcut.removeCut(self.vardb.getCut('2Lep_QMisIDVeto'))
            truthcut = truthcut & self.vardb.getCut('2Lep_ProbeQMisIDEvent')

            sp = sp.subprocess(cut=truthcut)

            print("\nZjets CF sp: {0}".format(sp.basecut.cutnamelist))

            return sp


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

            print inputgroup

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
            TTcut=''
            weight=1.0
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

            if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\nWenujets sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

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

            print inputgroup

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
            TTcut=''
            weight=1.0
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

            if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\nWmunujets sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

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

            print inputgroup

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
            TTcut=''
            weight=1.0
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

            if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\nWtaunujets sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

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

        latexname = 't, tW, tZ, tWZ, ttWW, 4t'
        #latexname = 'tZ, tWZ, ttWW, 4t'
        #latexname = 'rare top'

        colour = kAzure + 1

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                ('tops', 'tZ'),
                ('tops', 'tW'),
                ('tops', 'singlet'),
                ('tops', '4top'),
                ('tops', 'ttWW'),
                ('tops', 'tWZDR'),
                         ]

            print inputgroup

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
            TTcut=''
            weight=1.0

            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

            if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\nTop sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

            return sp


    class TTBar(Process):

        latexname = 't#bar{t}'
        colour = kAzure + 8

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                ('tops', 'ttbar_nonallhad'),
                #('tops', 'ttbar_dilep'),
                         ]

            print inputgroup

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
            TTcut  =''
            weight=1.0

            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

            if TTcut is not '':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\nTTBar sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

            return sp

    # The only difference with the above class is
    # that we plot only events with at least one !prompt lepton, and veto charge flip
    #
    class TTBarClosure(Process):

        latexname = 't#bar{t}'
        colour = kAzure + 8

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                ('tops', 'ttbar_nonallhad'),
                #('tops', 'ttbar_dilep'),
                         ]

            print inputgroup

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

            print("\nTTBarClosure sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

            return sp


    class TopCF(Process):

        latexname = 'tops (QMisID only)'
        colour = kAzure - 4

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                ('tops', 'ttbar_nonallhad'),
                #('tops', 'ttbar_dilep'),
                ('tops', 'tZ'),
                ('tops', 'tW'),
                ('tops', 'singlet'),
                ('tops', '4top'),
                ('tops', 'ttWW'),
                         ]

            print inputgroup

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
            TTcut=''
            weight=1.0

            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'
            if TTcut != '':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            # plot only events where at least one lepton is charge flip. Remove req. where both lep must be prompt
            #
            #truthcut = category.cut.swapCut(self.vardb.getCut('2Lep_PurePromptEvent'),self.vardb.getCut('2Lep_QMisIDEvent'))
            truthcut = category.cut.removeCut(self.vardb.getCut('2Lep_ProbePromptEvent'))
            truthcut = truthcut.removeCut(self.vardb.getCut('2Lep_QMisIDVeto'))
            truthcut = truthcut & self.vardb.getCut('2Lep_ProbeQMisIDEvent')

            sp = sp.subprocess(cut=truthcut)

            print("\nTopCF sp: {0}".format(sp.basecut.cutnamelist))

            return sp


    class Diboson(Process):

        latexname = 'WZ, ZZ, WW'
        colour = kYellow - 9

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                #('Diboson', '*'),
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
                #('Diboson', 'WW'),
                #('Diboson', 'WZ'),
                #('Diboson', 'ZZ'),
                ('Diboson', 'WW_SHv21_improved'),
                ('Diboson', 'WZ_SHv21_improved'),
                ('Diboson', 'ZZ_SHv21_improved'),
                         ]

            print inputgroup

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
            TTcut=''
            weight=1.0

            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

            if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\nDiboson sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

            return sp


    class DibosonCF(Process):

        latexname = 'WW, W#gamma, ZZ#rightarrow ll#nu#nu (QMisID only)'
        colour = kAzure - 9

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                #('Diboson', '*'),
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
                ('Diboson', 'WW'),
                ('Diboson', 'WZ'),
                ('Diboson', 'ZZ'),
                         ]

            print inputgroup

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
            TTcut=''
            weight=1.0

            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'
            if TTcut != '':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            # plot only events where at least one lepton is charge flip. Remove req. where both lep must be prompt
            #
            #truthcut = category.cut.swapCut(self.vardb.getCut('2Lep_PurePromptEvent'),self.vardb.getCut('2Lep_QMisIDEvent'))
            truthcut = category.cut.removeCut(self.vardb.getCut('2Lep_ProbePromptEvent'))
            truthcut = truthcut.removeCut(self.vardb.getCut('2Lep_QMisIDVeto'))
            truthcut = truthcut & self.vardb.getCut('2Lep_ProbeQMisIDEvent')

            sp = sp.subprocess(cut=truthcut)

            print("\nDibosonCF sp: {0}".format(sp.basecut.cutnamelist))

            return sp


    class HtoZZ(Process):

        latexname = 'H #rightarrow ZZ'
        colour = kTeal + 9

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                    ('HtoZZ', '125'),
                ]
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
            TTcut=''
            weight=1.0

            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

            sp = sp.subprocess(cut=category.cut,eventweight=weight)

            if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            print("\nHtoZZ sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

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

        latexname = 'QMisID (MC)'
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


    class ChargeFlip(Process):

        latexname = 'QMisID'
        colour = kAzure - 4

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                    ('Data', 'physics_Main'),
                         ]

            print inputgroup

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

            TTcut  = ''
            weight = '0.0'

            # Categories w/ any of these cuts must get the QMisID weight...
            #
            cut_w_el  = ['2Lep_ElEl_Event','2Lep_MuEl_Event','2Lep_ElMu_Event']
            # ... and these must not!
            #
            cut_wo_el = ['2Lep_MuMu_Event']

            # Take the category cut
            #
            basecut = category.cut

            print("\nChargeFlip initial cut: {0}".format(basecut.cutnamelist))

            if ( self.parent.channel=='TwoLepSS' ) or ( self.parent.channel=='ThreeLep' ):
                TTcut = 'TT'

            # Remove any truth cut
            #
            basecut = basecut.removeCut(self.vardb.getCut('2Lep_PurePromptEvent'))

            if TTcut != '':
                basecut = basecut & self.vardb.getCut(TTcut)

            # If SS cut in category, need to switch to OS before applying QMisID weight.
            # Otherwise, just forget about QMisID, at all.
            #
            if ("2Lep_SS") in category.cut.cutname:
                basecut = basecut.swapCut(self.vardb.getCut('2Lep_SS'), -self.vardb.getCut('2Lep_SS'))
            else:
                sp = sp.subprocess(cut=basecut,eventweight=weight)
                print("\nChargeFlip sp: {0}".format(sp.basecut.cutnamelist))
                return sp

            if bool([ cut for cut in cut_w_el if cut in category.cut.cutname ]):

                weight = 'QMisIDWeight[0]'
                sp = sp.subprocess(cut=basecut,eventweight=weight)
                print("\nChargeFlip sp: {0}".format(sp.basecut.cutnamelist))

            elif bool([ cut for cut in cut_wo_el if cut in category.cut.cutname ]):

                weight = '0.0'
                sp = sp.subprocess(cut=basecut,eventweight=weight)
                print("\nChargeFlip sp: {0}".format(sp.basecut.cutnamelist))

            else:

                weight_SF = '0.0'
                weight_OF = '0.0'

                sp_SF = None
                sp_OF = None

                if ( 'ElRealFakeRateCR' ) in category.cut.cutname:

                    cut_SF = basecut & self.vardb.getCut('2Lep_ElEl_Event')
                    cut_OF = basecut & self.vardb.getCut('2Lep_OF_Event')

                    weight_SF = weight_OF = 'QMisIDWeight[0]'

                    sp_SF = sp.subprocess(cut=cut_SF,eventweight=weight_SF)
                    sp_OF = sp.subprocess(cut=cut_OF,eventweight=weight_OF)

                    print("\nChargeFlip sp_SF: {0}, weight: {1}".format(sp_SF.basecut.cutnamelist, weight_SF))
                    print("\nChargeFlip sp_OF: {0}, weight: {1}".format(sp_OF.basecut.cutnamelist, weight_OF))

                elif ( 'MuRealFakeRateCR' ) in category.cut.cutname:

                    cut_SF = basecut & self.vardb.getCut('2Lep_MuMu_Event')
                    cut_OF = basecut & self.vardb.getCut('2Lep_OF_Event')

                    weight_SF = '0.0' # mu-mu region MUST get zero QMisID weight!
                    weight_OF = 'QMisIDWeight[0]'

                    sp_SF = sp.subprocess(cut=cut_SF,eventweight=weight_SF)
                    sp_OF = sp.subprocess(cut=cut_OF,eventweight=weight_OF)

                    print("\nChargeFlip sp_SF: {0}, weight: {1}".format(sp_SF.basecut.cutnamelist, weight_SF))
                    print("\nChargeFlip sp_OF: {0}, weight: {1}".format(sp_OF.basecut.cutnamelist, weight_OF))

                else:

                    cut_elel = basecut & self.vardb.getCut('2Lep_ElEl_Event')
                    cut_mumu = basecut & self.vardb.getCut('2Lep_MuMu_Event')
                    cut_OF   = basecut & self.vardb.getCut('2Lep_OF_Event')

                    weight_mumu = '0.0' # mu-mu region MUST get zero QMisID weight!
                    weight_elel = weight_OF = 'QMisIDWeight[0]'

                    sp_elel = sp.subprocess(cut=cut_elel,eventweight=weight_elel)
                    sp_mumu = sp.subprocess(cut=cut_mumu,eventweight=weight_mumu)
                    sp_OF = sp.subprocess(cut=cut_OF,eventweight=weight_OF)

                    print("\nChargeFlip sp_elel: {0}, weight: {1}".format(sp_elel.basecut.cutnamelist, weight_elel))
                    print("\nChargeFlip sp_mumu: {0}, weight: {1}".format(sp_mumu.basecut.cutnamelist, weight_mumu))
                    print("\nChargeFlip sp_OF: {0}, weight: {1}".format(sp_OF.basecut.cutnamelist, weight_OF))

                    sp_SF = sp_elel + sp_mumu

                sp = sp_SF + sp_OF

            return sp

    class ZpeakSidebandBkg(Process):

        latexname = 'Background (from Z sidebands)'
        colour = kBlue

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                    ('Data', 'physics_Main'),
                         ]

            print inputgroup

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

            weight=1.0

            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut  = 'TT'
                TLcut  = 'TL'
                LTcut  = 'LT'
                LLcut  = 'LL'
                weight = 'MMWeight[0]'

            # Remove any truth cut
            #
            basecut = category.cut.removeCut(self.vardb.getCut('2Lep_PurePromptEvent'))

            # QMisID subtraction? ---> CHECK IN PROGRESS

            sp_TT  = sp.subprocess(cut=basecut & self.vardb.getCut(TTcut), eventweight=weight)
            sp_TL  = sp.subprocess(cut=basecut & self.vardb.getCut(TLcut), eventweight=weight)
            sp_LT  = sp.subprocess(cut=basecut & self.vardb.getCut(LTcut), eventweight=weight)
            sp_LL  = sp.subprocess(cut=basecut & self.vardb.getCut(LLcut), eventweight=weight)

            print(" ")
            print("ZpeakSidebandBkg - TT sp: {0}, weight: {1}".format(sp.basecut.cutnamelist, weight))

            sp = sp_TT + sp_TL + sp_LT + sp_LL

            return sp

    class FakesFF(Process):

        latexname = 'FakesFF'
        colour = kAzure - 9

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                    ('Data', 'physics_Main'),
                         ]

            print inputgroup

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

            print inputgroup

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

            # Remove any truth cut (we use data!)
            #
            basecut = category.cut.removeCut(self.vardb.getCut('2Lep_PurePromptEvent'))

            # Subtract QMisID from final MM yield in TT, since MM estimates both

	    QMISID_SUB = True

            if QMISID_SUB:
	    	sublist = [ item for item in self.parent.sub_backgrounds ]
            	sp_QMisID_TT = None
            	for process in sublist:
            	    if not ( process == "ChargeFlip" or process == "ChargeFlipMC" ): continue

                    if ("2Lep_MuMu_Event") in category.cut.cutname:
                        print("NO QMisID subtraction for muons!!")
                        continue

	    	    QMisIDcut	 = basecut.swapCut(self.vardb.getCut('2Lep_SS'), -self.vardb.getCut('2Lep_SS'))
	    	    QMisIDweight = 'QMisIDWeight[0]'

            	    sp_QMisID_TT = (self.parent.procmap[process].base(treename,category,options)).subprocess(cut=QMisIDcut & self.vardb.getCut(TTcut), eventweight=QMisIDweight)
            	    print(" ")
            	    print("FakesMM - Subtracting QMisID TT sp: {0}, weight: {1}".format(sp_QMisID_TT.basecut.cutnamelist, QMisIDweight))
            	    print(" ")
            	    print("FakesMM - QMisID - TT : {0}".format(sp_QMisID_TT.numberstats()))

            sp_TT  = sp.subprocess(cut=basecut & self.vardb.getCut(TTcut), eventweight=weight)
            sp_TL  = sp.subprocess(cut=basecut & self.vardb.getCut(TLcut), eventweight=weight)
            sp_LT  = sp.subprocess(cut=basecut & self.vardb.getCut(LTcut), eventweight=weight)
            sp_LL  = sp.subprocess(cut=basecut & self.vardb.getCut(LLcut), eventweight=weight)

            print(" ")
            print("FakesMM - TT sp: {0}, weight: {1}".format(sp_TT.basecut.cutnamelist, weight))
            print("FakesMM - TL sp: {0}, weight: {1}".format(sp_TL.basecut.cutnamelist, weight))
            print("FakesMM - LT sp: {0}, weight: {1}".format(sp_LT.basecut.cutnamelist, weight))
            print("FakesMM - LL sp: {0}, weight: {1}".format(sp_LL.basecut.cutnamelist, weight))

            sp = sp_TT + sp_TL + sp_LT + sp_LL

            if QMISID_SUB and not ("2Lep_MuMu_Event") in category.cut.cutname:
            	print(" ")
            	print("FakesMM - Before QMisID subtraction : {0}".format(sp.numberstats()))
                sp = sp - sp_QMisID_TT
            	print("FakesMM - After QMisID subtraction : {0}".format(sp.numberstats()))

            return sp


    class FakesABCD(Process):

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

            print inputgroup

            trees = self.inputs.getTrees(treename, inputgroup)
            return self.subprocess(trees=trees)

        def calcTheta(self, sp_num, sp_denom, stream=None, options={}):

            if not sp_denom.number():
                print ("ERROR: Cannot calculate theta transfer factor! Denominator = 0")

            print("N: ", sp_num.numberstats())
            print("D: ", sp_denom.numberstats())

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

            print("\nFakesABCD\n")

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
            weight   = 1.0
            #weightMC ='weight_lepton_trig_HTop[0] * weight_lepton_reco_HTop[0] * weight_lepton_iso_HTop[0] * weight_lepton_ID_HTop[0] * weight_lepton_TTVA_HTop[0] * weight_jet__MV2c20_SFFix77[0]'
            weightMC = 'weight_event_trig * weight_event_lep * tauSFTight * JVT_EventWeight * MV2c10_70_EventWeight'

            if ( self.parent.channel=='TwoLepSS' ):
                TTCut  = self.vardb.getCut('TT')
                TLCut  = self.vardb.getCut('TL')
                LTCut  = self.vardb.getCut('LT')
                TelLmuCut = self.vardb.getCut('TelLmu')
                LelTmuCut = self.vardb.getCut('LelTmu')
                TmuLelCut = self.vardb.getCut('TmuLel')
                LmuTelCut = self.vardb.getCut('LmuTel')

            TL_LT_Cut = (TLCut | LTCut)

            # take the base suprocess (DATA)
            #
            sp = self.base(treename, category, options)

            print("base sp: {0}".format(category.cut.cutnamelist))

            # Remove the cuts defining the flavour composition (this is made just for calculating thetas...)
            #
            basecut = category.cut.removeCut(self.vardb.getCut('2Lep_ElEl_Event'))
            basecut = basecut.removeCut(self.vardb.getCut('2Lep_MuMu_Event'))
            basecut = basecut.removeCut(self.vardb.getCut('2Lep_OF_Event'))

            # Remove the cuts defining the jet multiplicity
            #
            basecut = basecut.removeCut(self.vardb.getCut('2Lep_NJet_SR'))
            basecut = basecut.removeCut(self.vardb.getCut('2Lep_NJet_CR'))

            if TTHBackgrounds.theta['El'][0] == 999.0 :

                print ("Calculating theta_el from data in regions C,D...")

                # define selection for region C and D (ee)
                #
                cut_sp_C_el = basecut & self.vardb.getCut('2Lep_NJet_CR') & self.vardb.getCut('2Lep_ElEl_Event') & TTCut & self.vardb.getCut('2Lep_Zsidescut') & self.vardb.getCut('2Lep_ElEtaCut')
                cut_sp_D_el = basecut & self.vardb.getCut('2Lep_NJet_CR') & self.vardb.getCut('2Lep_ElEl_Event') & TL_LT_Cut & self.vardb.getCut('2Lep_Zsidescut') & self.vardb.getCut('2Lep_ElEtaCut')

                # Lower pT threshold of subleading lepton to enrich in fakes
                #
                #cut_sp_C_el     = cut_sp_C_el.swapCut(self.vardb.getCut('2Lep_NLep'),self.vardb.getCut('2Lep_NLep_Relaxed'))
                #cut_sp_D_el     = cut_sp_D_el.swapCut(self.vardb.getCut('2Lep_NLep'),self.vardb.getCut('2Lep_NLep_Relaxed'))

                sp_C_el = sp.subprocess( cut = cut_sp_C_el, eventweight=weight )
                sp_D_el = sp.subprocess( cut = cut_sp_D_el, eventweight=weight )

                print("Region C (el) sp: {0}".format(cut_sp_C_el.cutnamelist))
                print("Region D (el) sp: {0}".format(cut_sp_D_el.cutnamelist))

                # get a list of stuff to subtract for region C and D (i.e, prompt MC and charge flips)
                #
                sublist = [ item for item in self.parent.sub_backgrounds ]

                # ... and now subtract!
                #
                for sample in sublist:

                    print ("Subtracting {0} from data in regions C,D...".format(sample))

                    this_cut_sp_C_el = cut_sp_C_el
                    this_cut_sp_D_el = cut_sp_D_el

                    this_weight = weightMC

                    # NB: here it is crucial to call .base() on the subprocess, otherwise the subprocess would have the cuts
                    # defined in its own __call__ method already applied, whcih in general is not what we want
                    # (e.g., it might have a TT selection applied, when we want to consider TL events instead...)
                    #
                    sub_sample_C_el = self.parent.procmap[sample].base(treename,category,options)
                    sub_sample_D_el = self.parent.procmap[sample].base(treename,category,options)

                    if sample == "ChargeFlipMC":
                        this_cut_sp_C_el = this_cut_sp_C_el.swapCut(self.vardb.getCut('2Lep_PurePromptEvent'),self.vardb.getCut('2Lep_QMisIDEvent'))
                        this_cut_sp_D_el = this_cut_sp_D_el.swapCut(self.vardb.getCut('2Lep_PurePromptEvent'),self.vardb.getCut('2Lep_QMisIDEvent'))
                    if sample == "ChargeFlip":
                        this_cut_sp_C_el = this_cut_sp_C_el.removeCut(self.vardb.getCut('2Lep_PurePromptEvent'))
                        this_cut_sp_C_el = this_cut_sp_C_el.swapCut(self.vardb.getCut('2Lep_SS'), -self.vardb.getCut('2Lep_SS'))
                        this_cut_sp_D_el = this_cut_sp_D_el.removeCut(self.vardb.getCut('2Lep_PurePromptEvent'))
                        this_cut_sp_D_el = this_cut_sp_D_el.swapCut(self.vardb.getCut('2Lep_SS'), -self.vardb.getCut('2Lep_SS'))
                        this_weight = 'QMisIDWeight[0]'

                    sub_sample_C_el = sub_sample_C_el.subprocess( cut = this_cut_sp_C_el, eventweight=this_weight )
                    sub_sample_D_el = sub_sample_D_el.subprocess( cut = this_cut_sp_D_el, eventweight=this_weight )

                    print ("sub_sample_C_el ={0} ".format(sub_sample_C_el.basecut.cutnamelist))

                    print ("C (el) - yields data: ", sp_C_el.numberstats())
                    print ("C (el) - yields bkg: ", sub_sample_C_el.numberstats())
                    sp_C_el = sp_C_el - sub_sample_C_el
                    print ("C (el) - yields after sub: ", sp_C_el.numberstats())

                    # *********************************************

                    print ("sub_sample_D_el ={0} ".format(sub_sample_D_el.basecut.cutnamelist))

                    print ("D (el) - yields data: ", sp_D_el.numberstats())
                    print ("D (el) - yields bkg: ", sub_sample_D_el.numberstats())
                    sp_D_el = sp_D_el - sub_sample_D_el
                    print ("D (el) - yields after sub: ", sp_D_el.numberstats())


                print ("---------------------------------------------------------------------\n")
                print ("C (el) - data yields after prompt/ch-flip subtraction: ", sp_C_el.numberstats())
                print ("D (el) - data yields after prompt/ch-flip subtraction: ", sp_D_el.numberstats())
                print ("\n---------------------------------------------------------------------")

                # derive theta factors for el
                #
                TTHBackgrounds.theta['El'] = self.calcTheta(sp_C_el,sp_D_el,stream='El')

            else :
                print ("Reading theta(el) value: {0} +- {1}".format(TTHBackgrounds.theta['El'][0], TTHBackgrounds.theta['El'][1]))


            if TTHBackgrounds.theta['Mu'][0] == 999.0 :

                print ("Calculating theta_mu from data in regions C,D...")

                # define selection for region C and D (mumu)
                #
                cut_sp_C_mu = basecut & self.vardb.getCut('2Lep_NJet_CR') & self.vardb.getCut('2Lep_MuMu_Event') & TTCut
                cut_sp_D_mu = basecut & self.vardb.getCut('2Lep_NJet_CR') & self.vardb.getCut('2Lep_MuMu_Event') & TL_LT_Cut

                # Lower pT threshold of subleading lepton to enrich in fakes
                #
                #cut_sp_C_mu     = cut_sp_C_mu.swapCut(self.vardb.getCut('2Lep_NLep'),self.vardb.getCut('2Lep_NLep_Relaxed'))
                #cut_sp_D_mu     = cut_sp_D_mu.swapCut(self.vardb.getCut('2Lep_NLep'),self.vardb.getCut('2Lep_NLep_Relaxed'))

                sp_C_mu = sp.subprocess( cut = cut_sp_C_mu, eventweight=weight )
                sp_D_mu = sp.subprocess( cut = cut_sp_D_mu, eventweight=weight )

                print("Region C (mu) sp: {0}".format(cut_sp_C_mu.cutnamelist))
                print("Region D (mu) sp: {0}".format(cut_sp_D_mu.cutnamelist))

                # get a list of stuff to subtract for region C and D (i.e, prompt MC and charge flips)
                #
                sublist = [ item for item in self.parent.sub_backgrounds ]

                # ... and now subtract!
                #
                for sample in sublist:

                    print ("Subtracting {0} from data in regions C,D...".format(sample))

                    this_cut_sp_C_mu = cut_sp_C_mu
                    this_cut_sp_D_mu = cut_sp_D_mu

                    this_weight = weightMC

                    # NB: here it is crucial to call .base() on the subprocess, otherwise the subprocess would have the cuts
                    # defined in its own __call__ method already applied, whcih in general is not what we want
                    # (e.g., it might have a TT selection applied, when we want to consider TL events instead...)
                    #
                    sub_sample_C_mu = self.parent.procmap[sample].base(treename,category,options)
                    sub_sample_D_mu = self.parent.procmap[sample].base(treename,category,options)

                    if ( sample == "ChargeFlipMC" ) or ( sample == "ChargeFlip" ):
                        print("NO QMisID subtraction for muons!!")
                        continue

                    sub_sample_C_mu = sub_sample_C_mu.subprocess( cut = this_cut_sp_C_mu, eventweight=this_weight )
                    sub_sample_D_mu = sub_sample_D_mu.subprocess( cut = this_cut_sp_D_mu, eventweight=this_weight )

                    print ("sub_sample_C_mu ={0} ".format(sub_sample_C_mu.basecut.cutnamelist))

                    print ("C (mu) - yields data: ", sp_C_mu.numberstats())
                    print ("C (mu) - yields bkg: ", sub_sample_C_mu.numberstats())
                    sp_C_mu = sp_C_mu - sub_sample_C_mu
                    print ("C (mu) - yields after sub: ", sp_C_mu.numberstats())

                    # *********************************************

                    print ("sub_sample_D_mu ={0} ".format(sub_sample_D_mu.basecut.cutnamelist))

                    print ("D (mu) - yields data: ", sp_D_mu.numberstats())
                    print ("D (mu) - yields bkg: ", sub_sample_D_mu.numberstats())
                    sp_D_mu = sp_D_mu - sub_sample_D_mu
                    print ("D (mu) - yields after sub: ", sp_D_mu.numberstats())


                print ("---------------------------------------------------------------------\n")
                print ("C (mu) - data yields after prompt/ch-flip subtraction: ", sp_C_mu.numberstats())
                print ("D (mu) - data yields after prompt/ch-flip subtraction: ", sp_D_mu.numberstats())
                print ("\n---------------------------------------------------------------------")

                # derive theta factors for mu
                #
                TTHBackgrounds.theta['Mu'] = self.calcTheta(sp_C_mu,sp_D_mu,stream='Mu')

            else :
                print ("Reading theta(mu) value: {0} +- {1}".format(TTHBackgrounds.theta['Mu'][0], TTHBackgrounds.theta['Mu'][1]))


            # Define Region B,  depending on which flavour composition we are looking at:
            #
            cut_sp_B_SF     = category.cut.swapCut(self.vardb.getCut('2Lep_NJet_CR'),self.vardb.getCut('2Lep_NJet_SR')) & TL_LT_Cut
            cut_sp_B_OF_Lel = category.cut.swapCut(self.vardb.getCut('2Lep_NJet_CR'),self.vardb.getCut('2Lep_NJet_SR')) & (LelTmuCut | TmuLelCut)
            cut_sp_B_OF_Lmu = category.cut.swapCut(self.vardb.getCut('2Lep_NJet_CR'),self.vardb.getCut('2Lep_NJet_SR')) & (TelLmuCut | LmuTelCut)

            if not ("2Lep_OF_Event") in category.cut.cutname:

                sp_B = sp.subprocess( cut = cut_sp_B_SF, eventweight=weight )

                sublist = [ item for item in self.parent.sub_backgrounds ]
                for sample in sublist:

                    print ("Subtracting {0} from data in region B...".format(sample))

                    this_cut_sp_B_SF = cut_sp_B_SF
                    this_weight = weightMC

                    if ( ( sample == "ChargeFlipMC" ) or ( sample == "ChargeFlip" ) ):

                        if ("2Lep_MuMu_Event") in category.cut.cutname:
                            print("NO QMisID subtraction for muons!!")
                            continue

                        this_weight = 'QMisIDWeight[0]'
                        this_cut_sp_B_SF = this_cut_sp_B_SF.removeCut(self.vardb.getCut('2Lep_PurePromptEvent'))
                        this_cut_sp_B_SF = this_cut_sp_B_SF.swapCut(self.vardb.getCut('2Lep_SS'), -self.vardb.getCut('2Lep_SS'))

                    sub_sample_B = self.parent.procmap[sample].base(treename,category,options)
                    sub_sample_B = sub_sample_B.subprocess( cut = this_cut_sp_B_SF, eventweight=this_weight )

                    sp_B = sp_B - sub_sample_B

                print ("Region B (mumu,ee), yield: ", sp_B.numberstats() , "\n")

                if ("2Lep_ElEl_Event") in category.cut.cutname:
                    sp_B = self.applyThetaFactor(sp_B,TTHBackgrounds.theta['El'])
                elif ("2Lep_MuMu_Event") in category.cut.cutname:
                    sp_B = self.applyThetaFactor(sp_B,TTHBackgrounds.theta['Mu'])
            else:

                sp_B_OF_Lel = sp.subprocess(cut=cut_sp_B_OF_Lel, eventweight=weight )
                sp_B_OF_Lmu = sp.subprocess(cut=cut_sp_B_OF_Lmu, eventweight=weight )

                sublist = [ item for item in self.parent.sub_backgrounds ]
                for sample in sublist:

                    print ("Subtracting {0} from data in region B...".format(sample))

                    this_cut_sp_B_OF_Lel = cut_sp_B_OF_Lel
                    this_cut_sp_B_OF_Lmu = cut_sp_B_OF_Lmu
                    this_weight = weightMC

                    if ( ( sample == "ChargeFlipMC" ) or ( sample == "ChargeFlip" ) ):

                        this_weight = 'QMisIDWeight[0]'
                        this_cut_sp_B_OF_Lel = this_cut_sp_B_OF_Lel.removeCut(self.vardb.getCut('2Lep_PurePromptEvent'))
                        this_cut_sp_B_OF_Lel = this_cut_sp_B_OF_Lel.swapCut(self.vardb.getCut('2Lep_SS'), -self.vardb.getCut('2Lep_SS'))
                        this_cut_sp_B_OF_Lmu = this_cut_sp_B_OF_Lmu.removeCut(self.vardb.getCut('2Lep_PurePromptEvent'))
                        this_cut_sp_B_OF_Lmu = this_cut_sp_B_OF_Lmu.swapCut(self.vardb.getCut('2Lep_SS'), -self.vardb.getCut('2Lep_SS'))

                    sub_sample_B_OF_Lel = self.parent.procmap[sample].base(treename,category,options)
                    sub_sample_B_OF_Lel = sub_sample_B_OF_Lel.subprocess( cut = this_cut_sp_B_OF_Lel, eventweight=this_weight )
                    sub_sample_B_OF_Lmu = self.parent.procmap[sample].base(treename,category,options)
                    sub_sample_B_OF_Lmu = sub_sample_B_OF_Lmu.subprocess( cut = this_cut_sp_B_OF_Lel, eventweight=this_weight )

                    sp_B_OF_Lel = sp_B_OF_Lel - sub_sample_B_OF_Lel
                    sp_B_OF_Lmu = sp_B_OF_Lmu - sub_sample_B_OF_Lmu

                print ("Region B (emu)","\n", "yield (loose el): ", sp_B_OF_Lel.numberstats() , "\n", "yield (loose mu): ", sp_B_OF_Lmu.numberstats(), "\n")

                sp_B = self.applyThetaFactor(sp_B_OF_Lmu,TTHBackgrounds.theta['Mu']) + self.applyThetaFactor(sp_B_OF_Lel,TTHBackgrounds.theta['El'])

            print ("=================>\n")
            print ("Final fakes yield: ", sp_B.numberstats() ,"\n")

            return sp_B


    class FakesClosureMM(Process):

        latexname = 'FakesMM - t#bar{t}'
        colour = kTeal -9

        name = 'FakesClosureMM'

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                    ('tops', 'ttbar_nonallhad'),
                    #('tops', 'ttbar_dilep'),
                         ]

            print inputgroup

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

            TTcut=''
            TLcut=''
            LTcut=''
            LLcut=''
            weight=1.0

            if ( self.parent.channel=='TwoLepSS' ) or ( self.parent.channel=='ThreeLep' ):
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

            print(" ")
            print("FakesClosureMM - TT sp: {0}, weight: {1}".format(sp_TT.basecut.cutnamelist, weight))
            print("FakesClosureMM - TL sp: {0}, weight: {1}".format(sp_TL.basecut.cutnamelist, weight))
            print("FakesClosureMM - LT sp: {0}, weight: {1}".format(sp_LT.basecut.cutnamelist, weight))
            print("FakesClosureMM - LL sp: {0}, weight: {1}".format(sp_LL.basecut.cutnamelist, weight))

            sp = sp_TT + sp_TL + sp_LT + sp_LL

            return sp

    class FakesClosureABCD(Process):

        latexname = 'Fakes #theta method - t#bar{t}'
        colour = kCyan - 9

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                    ('tops', 'ttbar_nonallhad'),
                ]

            print inputgroup

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

            print("\nFakesClosureABCD\n")

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
                TelLmuCut = self.vardb.getCut('TelLmu')
                LelTmuCut = self.vardb.getCut('LelTmu')
                TmuLelCut = self.vardb.getCut('TmuLel')
                LmuTelCut = self.vardb.getCut('LmuTel')

            TL_LT_Cut = (TLCut | LTCut)

            # take the base suprocess (i.e, TTBar)
            #
            sp = self.base(treename, category, options)

            print("base sp: {0}".format(category.cut.cutnamelist))

            # Remove the cuts defining the flavour composition (this is made just for calculating thetas...)
            #
            basecut = category.cut.removeCut(self.vardb.getCut('2Lep_ElEl_Event'))
            basecut = basecut.removeCut(self.vardb.getCut('2Lep_MuMu_Event'))
            basecut = basecut.removeCut(self.vardb.getCut('2Lep_OF_Event'))

            # Remove the cuts defiining the jet multiplicity
            basecut = basecut.removeCut(self.vardb.getCut('2Lep_NJet_SR'))
            basecut = basecut.removeCut(self.vardb.getCut('2Lep_NJet_CR'))

            # For closure test, consider only TTBar non-prompt, vetoing all charge flips!
            #
            basecut = basecut & self.vardb.getCut('2Lep_NonPromptEvent')

            if TTHBackgrounds.theta_MC['El'][0] == 999.0 :

                print ("Calculating theta_el from TTBar in regions C,D...")

                # define selection for region C and D (ee)
                #
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

                # derive theta factors for el
                #
                TTHBackgrounds.theta_MC['El'] = self.calcTheta(sp_C_el,sp_D_el,stream='El')

            else :
                print ("Reading theta(el) value: {0} +- {1}".format(TTHBackgrounds.theta_MC['El'][0], TTHBackgrounds.theta_MC['El'][1]))


            if TTHBackgrounds.theta_MC['Mu'][0] == 999.0 :

                print ("Calculating theta_mu from TTBar in regions C,D...")

                # define selection for region C and D (mumu)
                #
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

                # derive theta factors for mu
                #
                TTHBackgrounds.theta_MC['Mu'] = self.calcTheta(sp_C_mu,sp_D_mu,stream='Mu')

            else :
                print ("Reading theta(mu) value: {0} +- {1}".format(TTHBackgrounds.theta_MC['Mu'][0], TTHBackgrounds.theta_MC['Mu'][1]))


            # Define Region B,  depending on which flavour composition we are looking at:
            # take TTbar MC events with fakes, vetoing all prompts and charge flips, and reweight it by the theta factors measured in ttbar MC
            #
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
    #
    class FakesClosureDataABCD(Process):

        latexname = 'Fakes #theta method'
        colour = kCyan - 9

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                    ('tops', 'ttbar_nonallhad'),
                ]

            print inputgroup

            trees = self.inputs.getTrees(treename, inputgroup)
            return self.subprocess(trees=trees) * self.parent.norm_factor

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
            #weightMC='weight_lepton_trig_HTop[0] * weight_lepton_reco_HTop[0] * weight_lepton_iso_HTop[0] * weight_lepton_ID_HTop[0] * weight_lepton_TTVA_HTop[0] * weight_jet__MV2c20_SFFix77[0]'
            weightMC = 'weight_event_trig * weight_event_lep * tauSFTight * JVT_EventWeight * MV2c10_70_EventWeight'

            if self.parent.channel=='TwoLepSS':
                TTCut  = self.vardb.getCut('TT')
                TLCut  = self.vardb.getCut('TL')
                LTCut  = self.vardb.getCut('LT')
                TelLmuCut = self.vardb.getCut('TelLmu')
                LelTmuCut = self.vardb.getCut('LelTmu')
                TmuLelCut = self.vardb.getCut('TmuLel')
                LmuTelCut = self.vardb.getCut('LmuTel')

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
            sp_Lel = sp.subprocess(cut=cut_sp_OF_Lel, eventweight=weightMC )
            sp_Lmu = sp.subprocess(cut=cut_sp_OF_Lmu, eventweight=weightMC )
            sp_final = ( sp_Lel * TTHBackgrounds.theta['El'][0] ) + ( sp_Lmu * TTHBackgrounds.theta['Mu'][0] )

            print ("=================>\n")
            print ("Region closure sp: {0}".format(sp_final.basecut.cutnamelist))
            print ("yield: ", sp_final.numberstats())

            return sp_final

    class FakesMC(Process):

        latexname = 'Fakes (MC)'
        colour = kCyan -9

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                ('tops', 'ttbar_nonallhad'),
                ('Diboson','*'),
                ('W+jetsBFilter', '*'),
                ('W+jetsCFilterBVeto', '*'),
                ('W+jetsCVetoBVeto', '*'),
                         ]

            print inputgroup

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
            TTcut=''
            weight=1.0

            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'
            if TTcut != '':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            # plot only events where at least one lepton is non-prompt (+ch-flip veto). Remove any req. where both lep must be prompt
            #
            truthcut = category.cut.swapCut(self.vardb.getCut('2Lep_PurePromptEvent'),self.vardb.getCut('2Lep_NonPromptEvent'))

            sp = sp.subprocess(cut=truthcut)

            print("\nFakesMC sp: {0}".format(sp.basecut.cutnamelist))

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
            sp = self.subprocess(trees=trees) * self.parent.norm_factor
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
