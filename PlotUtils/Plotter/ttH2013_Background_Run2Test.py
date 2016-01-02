import os, sys, math

sys.path.append(os.path.abspath(os.path.curdir))
from Plotter.BackgroundTools_ttH import loadSamples, drawText, Category, Background, Process, VariableDB, Variable, Cut, Systematics, Category
from ROOT import TColor, kBlack, kWhite, kBlue, kRed, kYellow, kAzure, kTeal, kSpring, TLegend, TLatex, TCanvas, TH1I, TFile

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

class TTHBackgrounds2013(Background):
    backgrounds = [ 'TTBarW', 'TTBarZ', 'Top', 'TopCF', 'Diboson', 'DibosonCF', 'HtoZZ', 'ZjetsLF', 'Zjets', 'Wjets', 'Prompt', 'ChargeFlip', 'FakesFF', 'FakesMM']#this is the list of backgrounds that will be calculated
    signals = ['TTBarH']
    observed = ['Observed']
    #luminosity = 13.0258
    luminosity = 5.4
    channel= 'TwoLepSS' # 'TwoLepCR', 'ThreeLep', 'FourLep'
    eventweight = 'weight_CF'
    #eventweight = 1.
    #eventweight = 'evtsel_weight*evtsel_weight_el*evtsel_weight_mu*evtsel_bjet_weight*evtsel_weight_lep_trigger*weight_CF'
    #eventweight = 'evtsel_weight*evtsel_weight_el*evtsel_weight_mu*evtsel_bjet_weight*evtsel_weight_lep_trigger*weight_MM*weight_CF'
    #eventweight = 'evtsel_weight*evtsel_weight_el*evtsel_weight_mu*evtsel_bjet_weight*evtsel_weight_lep_trigger*weight_FF*weight_CF'
    #eventweight = 1.
    useEmbedding = False
    useZCorrections = False
    RQCD = {
        #'El': (1.000, 0.051),
        #'Mu': (1.177, 0.066),
        'El': (1.00, 0.05),
        'Mu': (1.11, 0.08),
    }

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
            leg.SetBorderSize(0)

        for l in legs[:mid]:
            leg1.AddEntry(l[0], l[1], l[2])
        leg1.Draw()
        if leg2:
            for l in legs[mid:]:
                leg2.AddEntry(l[0], l[1], l[2])
            leg2.Draw()

        lumtext = drawText(text="  \\int L dt = %.2g fb^{-1}"%(self.luminosity), x=.2, y=.87, size=0.04 * scale)
        cmetext = drawText(text="         #sqrt{s} = 8TeV", x=.2, y=.82, size=0.04 * scale)
        #atlastext = drawText(text="#bf{#it{ATLAS}} Internal", x=.2, y=.77, size=0.04 * scale)
        atlastext = drawText(text="#bf{#it{ATLAS}} Work In Progress", x=.2, y=.77, size=0.04 * scale)
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
        latexname = 'Data 2012'

        def base(self, treename='physics', category=None, options={}):
            #Contains the instuction of which tree load and eventually correct the units of the cross-setion (see the division /1000.). Note it is not automatically executed when the class is called. It is executed trought __call__ 
            inputgroup = [
                    #('Data', '*'),
                    ('Data', 'Egamma'),
                    ('Data', 'Muons'),
                ]
            trees = self.inputs.getTrees(treename, inputgroup)
            return self.subprocess(trees=trees)

        def __call__(self, treename='physics', category=None, options={}):
            #this is the default procedure which is applied to the class when it is called in the code. Here in addition to the load of the tree are operated also the selection, application of factors and systematics
            sp = self.base(treename, category, options)
            isocut=''
            weight=''
            if self.parent.channel=='TwoLepSS':
                isocut='2LepTT'
            elif self.parent.channel=='ThreeLep':
                isocut='3LepTT'
            #print 'weight ', weight, ' isocut ', isocut
            if isocut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(isocut))
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
            return sp

    class TTBarH(Process):
        latexname = 't\\bar{t} H'
        colour = kBlack

        def base(self, treename='physics', category=None, options={}):
            #hmass = options.get('hmass', '125')#300 is the default value if hmass is not in options, hmass is specified in one option passed to the plot function in the plotting script
            inputgroup = [
                #('ttH', hmass),
                ('ttH', '125'),
                ]
            trees = self.inputs.getTrees(treename, inputgroup)
            sp = self.subprocess(trees=trees) / 1000.
            #the division for /1000. is used to have the correct cross section
            return sp

        def __call__(self, treename='physics', category=None, options={}):
            #this part of code is used to fill systname_opts list which is then used to apply properly a kfactor systematics. See Prompt class for further details.
            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)
            isocut=''
            weight=''
            if self.parent.channel=='TwoLepSS':
                isocut='2LepTT'
            elif self.parent.channel=='ThreeLep':
                isocut='3LepTT'
            if isocut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(isocut))
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
            return sp

    class TTBarW(Process):
        latexname = 't\\bar{t} W'
        colour = kRed - 4

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                ('ttW', 'ttW'),
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
            isocut=''
            weight=''
            if self.parent.channel=='TwoLepSS':
                isocut='2LepTT'
            elif self.parent.channel=='ThreeLep':
                isocut='3LepTT'
            if isocut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(isocut))
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
            return sp

    class TTBarZ(Process):
        latexname = 't\\bar{t} Z'
        colour = kRed - 7

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                ('ttZ', 'ttZ'),
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
            isocut=''
            weight=''
            if self.parent.channel=='TwoLepSS':
                isocut='2LepTT'
            elif self.parent.channel=='ThreeLep':
                isocut='3LepTT'
            if isocut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(isocut))
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
            return sp

    class FakesFF(Process):
        latexname = 'FakesFF'
        colour = kAzure - 9

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                    ('Data', 'Egamma'),
                    ('Data', 'Muons'),
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
            notisocut=''
            notisocut2=''
            weight=''
            if self.parent.channel=='TwoLepSS':
                notisocut='2LepTL'
                notisocut2='2LepLL'
                weight='weight_FF'
            elif self.parent.channel=='ThreeLep':
                notisocut='3LepTL'
                notisocut2='3LepLL'
                weight='weight_FF'

            #Now doing prompt background subtractions from fakes. promptlist must contain either prompt and chargeflips
            promptlist = [ item for item in self.parent.backgrounds if 'Fakes' not in item ]
            #print 'prompt list after ',promptlist
            for sample in promptlist:
                sp = sp - self.parent.procmap[sample].base(treename, category, options) # here it is important to have used base otherwise at the prompt sample there would be already applied the selection specified in the category __call__ of that sample i.e. for example also the iso-iso cut which is orthogonal to the cuts done in this sample which are fake iso or fake fake.
            #print sp
            sp_if = sp.subprocess(cut=self.vardb.getCut(notisocut))
            sp_ff = sp.subprocess(cut=self.vardb.getCut(notisocut2))
            #print sp_if
            #print sp_ff
            #IN THIS VERSION OF THE CODE THE SIGN MINUS FOR FF EVENTS IS ALREADY IN THE WEIGHT SO THE SUBRACTION OF THE DOUBLE FAKE COUNTING IS DOE ADDING SP_FF
            sp = sp_if + sp_ff
            #OLD WAY subtracting the double counting
            #sp = sp_if - sp_ff
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
            #print sp
            return sp

    class FakesMM(Process):
        latexname = 'FakesMM'
        colour = kTeal - 9

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                    ('Data', 'Egamma'),
                    ('Data', 'Muons'),
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
            isocut=''
            notisocut=''
            notisocut2=''
            weight=''
            if self.parent.channel=='TwoLepSS':
                isocut='2LepTT'
                notisocut='2LepTL'
                notisocut2='2LepLL'
                weight='weight_MM'
            elif self.parent.channel=='ThreeLep':
                isocut='3LepTT'
                notisocut='3LepTL'
                notisocut2='3LepLL'
                weight='weight_MM'

            #THERE IS NO PROMPT SUBTRACTION IN THE MATRIX METHOD
            #promptlist = [ item for item in self.parent.backgrounds if 'Fakes' not in item ]
            #for sample in promptlist:
            #    sp = sp - self.parent.procmap[sample].base(treename, category, options) # here it is important to have used base otherwise at the prompt sample there would be already applied the selection specified in the category __call__ of that sample i.e. for example also the iso-iso cut which is orthogonal to the cuts done in this sample which are fake iso or fake fake.
            #SINCE ALL THE PIECES HAVE TO BE ADDED TOGETHER THERE IS NO NEED TO SPLIT II IF FF AND THEN ADD AGAIN THEM TOGHETER
            #sp_ii = sp.subprocess(cut=self.vardb.getCut(isocut))
            #sp_if = sp.subprocess(cut=self.vardb.getCut(notisocut))
            #sp_ff = sp.subprocess(cut=self.vardb.getCut(notisocut2))
            #print sp_if
            #print sp_ff
            #sp = sp_ii + sp_if + sp_ff
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
            #print sp
            return sp

    class ZjetsLF(Process):
        latexname = 'Z/\\gamma* + LF jets'
        colour = kWhite

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                    #('Z+jets', '*'),
                    ('SherpaZ+jets', '*'),
                    #('DYZ+jets', '*'),
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
            isocut=''
            weight=''
            if self.parent.channel=='TwoLepSS':
                isocut='2LepTT'
            elif self.parent.channel=='ThreeLep':
                isocut='3LepTT'
            #print 'weight ', weight, ' isocut ', isocut
            if isocut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(isocut))
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
            return sp
 
    class ZjetsHF(Process):
        latexname = 'Z/\\gamma* + HF jets'
        colour = kWhite

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                    ('Z+bb', '*'),
                    ('Z+cc', '*'),
                    ('DYZ+bb', '*'),
                    ('DYZ+cc', '*'),
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
            isocut=''
            weight=''
            if self.parent.channel=='TwoLepSS':
                isocut='2LepTT'
            elif self.parent.channel=='ThreeLep':
                isocut='3LepTT'
            #print 'weight ', weight, ' isocut ', isocut
            if isocut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(isocut))
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
            return sp
 
    class Zjets(Process):
        latexname = 'Z/\\gamma* + jets'
        colour = kWhite

        def base(self, treename='physics', category=None, options={}):
            zlf = self.parent.procmap['ZjetsLF'].base(treename, category, options)
            zhf = self.parent.procmap['ZjetsHF'].base(treename, category, options)
            return zlf+zhf

        def __call__(self, treename='physics', category=None, options={}):
            zlf = self.parent.procmap['ZjetsLF'](treename, category, options)
            zhf = self.parent.procmap['ZjetsHF'](treename, category, options)
            return zlf+zhf

    class Wjets(Process):
        latexname = 'W + jets'
        colour = kYellow

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                    ('W+jets', '*'),
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
            isocut=''
            weight=''
            if self.parent.channel=='TwoLepSS':
                isocut='2LepTT'
            elif self.parent.channel=='ThreeLep':
                isocut='3LepTT'
            #print 'weight ', weight, ' isocut ', isocut
            if isocut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(isocut))
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
            return sp
 
    class Top(Process):
        latexname = 'tZ, 4top, t\\bar{t}WW'
        colour = kBlue - 3

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                    #('tops', 'tZ'),
                    #('tops', '4top'),
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
            isocut=''
            weight=''
            if self.parent.channel=='TwoLepSS':
                isocut='2LepTT'
            elif self.parent.channel=='ThreeLep':
                isocut='3LepTT'
            if isocut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(isocut))
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
            return sp
        
    class TopCF(Process):
        latexname = 't\\bar{t}, single t'
        colour = kAzure - 4

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                    ('tops', 'ttbar'),
                    ('tops', 'singlet'),
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
            isocut=''
            weight=''
            if self.parent.channel=='TwoLepSS':
                isocut='2LepTT'
            elif self.parent.channel=='ThreeLep':
                isocut='3LepTT'
            if isocut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(isocut))
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
            return sp
        
    class Diboson(Process):
        latexname = 'WZ, ZZ, SSWW'
        colour = kYellow - 9

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                    ('Diboson', 'WZ'),
                    ('Diboson', 'ZZ'),
                    #('Diboson', 'Wgammastar'),
                    #('Diboson', 'SSWW'),
                    ('Zgamma',  'Zgammastar'),
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
            isocut=''
            weight=''
            if self.parent.channel=='TwoLepSS':
                isocut='2LepTT'
            elif self.parent.channel=='ThreeLep':
                isocut='3LepTT'
            if isocut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(isocut))
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
            return sp

    class DibosonCF(Process):
        latexname = 'WW, W\\gamma, ZZ\\rightarrow ll\\nu\\nu'
        colour = kSpring + 1

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                    #('Diboson', 'WW'),
                    ('Diboson', 'Wgamma'),
                    ('Diboson', 'ZZnunu'),
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
            isocut=''
            weight=''
            if self.parent.channel=='TwoLepSS':
                isocut='2LepTT'
            elif self.parent.channel=='ThreeLep':
                isocut='3LepTT'
            if isocut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(isocut))
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
            return sp

    class HtoZZ(Process):
        latexname = 'H \\rightarrow ZZ'
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
            isocut=''
            weight=''
            if self.parent.channel=='TwoLepSS':
                isocut='2LepTT'
            elif self.parent.channel=='ThreeLep':
                isocut='3LepTT'
            if isocut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(isocut))
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
            return sp

    class Prompt(Process):
        latexname = 'Prompt'
        colour = kYellow - 9

        def base(self, treename='physics', category=None, options={}):
            diboson = self.parent.procmap['Diboson'].base(treename, category, options)
            top = self.parent.procmap['Top'].base(treename, category, options)
            htozz = self.parent.procmap['HtoZZ'].base(treename, category, options)
            ttbarw = self.parent.procmap['TTBarW'].base(treename, category, options)
            ttbarz = self.parent.procmap['TTBarZ'].base(treename, category, options)
            return diboson+top+htozz+ttbarw+ttbarz

        def __call__(self, treename='physics', category=None, options={}):
            diboson = self.parent.procmap['Diboson'](treename, category, options)
            top = self.parent.procmap['Top'](treename, category, options)
            htozz = self.parent.procmap['HtoZZ'](treename, category, options)
            ttbarw = self.parent.procmap['TTBarW'](treename, category, options)
            ttbarz = self.parent.procmap['TTBarZ'](treename, category, options)
            return diboson+top+htozz+ttbarw+ttbarz

    class ChargeFlip(Process):
        latexname = 'Charge Flip'
        colour = kAzure - 4

        def base(self, treename='physics', category=None, options={}):
            diboson = self.parent.procmap['DibosonCF'].base(treename, category, options)
            top = self.parent.procmap['TopCF'].base(treename, category, options)
            zjets = self.parent.procmap['Zjets'].base(treename, category, options)
            return diboson+top+zjets

        def __call__(self, treename='physics', category=None, options={}):
            diboson = self.parent.procmap['DibosonCF'](treename, category, options)
            top = self.parent.procmap['TopCF'](treename, category, options)
            zjets = self.parent.procmap['Zjets'](treename, category, options)
            return diboson+top+zjets

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
            isocut=''
            weight=''
            if self.parent.channel=='TwoLepSS':
                isocut='2LepTT'
                kfactors = {'Norm_corr': (1.0, 0.1),}
            elif self.parent.channel=='ThreeLep':
                isocut='3LepTT'
                kfactors = {'Norm_corr': (1.0, 0.1),}
            sp = self.parent.applyKfactor(sp, category, kfactors, kf_opts)#this command is used to apply the kfactor
            if isocut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(isocut))
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
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
            isocut=''
            weight=''
            if self.parent.channel=='TwoLepSS':
                isocut='2LepTT'
            elif self.parent.channel=='ThreeLep':
                isocut='3LepTT'
            if isocut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(isocut))
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
            return sp
'''
