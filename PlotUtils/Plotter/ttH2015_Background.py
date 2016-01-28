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
    backgrounds = [ 'TTBarW', 'TTBarZ', 'Top', 'TopCF', 'Diboson', 'DibosonCF', 'HtoZZ', 'ZjetsLF', 'Zjets', 'Wjets', 'Prompt', 'ChargeFlip', 'FakesFF', 'FakesMM']
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
            #this is the default procedure which is applied to the class when it is called in the code. Here in addition to the load of the tree are operated also the selection, application of factors and systematics
            sp = self.base(treename, category, options)
            TTcut=''
            weight=''
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

	    print 'Observed - weight name : ', weight, ', TTcut :', TTcut

	    if TTcut is not '':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
            return sp


    class TTBarH(Process):

	latexname = 't#bar{t} H'
        colour = kBlack

        def base(self, treename='physics', category=None, options={}):
            #hmass = options.get('hmass', '125')#300 is the default value if hmass is not in options, hmass is specified in one option passed to the plot function in the plotting script
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
            #this part of code is used to fill systname_opts list which is then used to apply properly a kfactor systematics. See Prompt class for further details.
            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
            sp = self.base(treename, category, options)
            TTcut=''
            weight=''
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

	    print 'TTBarH - weight name : ', weight, ', TTcut : ', TTcut

            if TTcut is not '':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
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
            weight=''
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

	    print 'TTBarW - weight name : ', weight, ', TTcut : ', TTcut

            if TTcut is not '':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
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
            weight=''
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

	    print 'TTBarZ - weight name : ', weight, ', TTcut : ', TTcut

            if TTcut is not '':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
            return sp

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
            weight=''
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TLcut='TL'
                LTcut='LT'
                LLcut='LL'
                weight='FFWeight[0]'

	    print 'FakesFF - weight name : ', weight, ', TLcut : ', TLcut, ', LTcut : ', LTcut, ', LLcut : ', LLcut

            #Now doing prompt background subtractions from fakes. promptlist must contain both prompt and chargeflips
	    #
            promptlist = [ item for item in self.parent.backgrounds if 'Fakes' not in item ]
            #print 'prompt list after ',promptlist
            for sample in promptlist:
                sp = sp - self.parent.procmap[sample].base(treename, category, options) # here it is important to have used base otherwise at the prompt sample there would be already applied the selection specified in the category __call__ of that sample i.e. for example also the iso-iso cut which is orthogonal to the cuts done in this sample which are fake iso or fake fake.
            #print sp
            sp_TL = sp.subprocess(cut=self.vardb.getCut(TLcut))
            sp_LT = sp.subprocess(cut=self.vardb.getCut(LTcut))
            sp_LL = sp.subprocess(cut=self.vardb.getCut(LLcut))

            #IN THIS VERSION OF THE CODE THE SIGN MINUS FOR FF EVENTS IS ALREADY IN THE WEIGHT SO THE SUBRACTION OF THE DOUBLE FAKE COUNTING IS DOE ADDING SP_FF
            sp = sp_TL + sp_LT +sp_LL
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
            weight=''
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut  = 'TT'
                TLcut  = 'TL'
                LTcut  = 'LT'
                LLcut  = 'LL'
                weight = 'MMWeight[0]'

            print 'FakesMM \n \tweight name : ', weight, ',\n \tTTcut : ', TTcut,',\n \tTLcut : ', TLcut, ',\n  \tLTcut : ', LTcut, ',\n \tLLcut : ', LLcut

            # THERE IS NO PROMPT SUBTRACTION IN THE MATRIX METHOD

            sp_TT  = sp.subprocess(cut=category.cut & self.vardb.getCut(TTcut), eventweight=weight)
            sp_TL  = sp.subprocess(cut=category.cut & self.vardb.getCut(TLcut), eventweight=weight)
            sp_LT  = sp.subprocess(cut=category.cut & self.vardb.getCut(LTcut), eventweight=weight)
            sp_LL  = sp.subprocess(cut=category.cut & self.vardb.getCut(LLcut), eventweight=weight)
            sp = sp_TT + sp_TL + sp_LT + sp_LL

            # THE ABOVE IS EQUIVALENT TO :
            #
            #sp = sp.subprocess(cut=category.cut, eventweight=weight)
            #
            # In fact, the MM weight for each event already contains the info on the event type: TT,TL,LT,LL :
            # no need to add any TT,LT ... cut!
            # (final estimate will be the sum of each contribution)
            #

            return sp

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
            weight=''
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut  = 'TT'
                TLcut  = 'TL'
                LTcut  = 'LT'
                LLcut  = 'LL'
                weight = 'MMWeight[0]'

            print 'FakesClosureMM \n \tweight name : ', weight, ',\n \tTTcut : ', TTcut,',\n \tTLcut : ', TLcut, ',\n  \tLTcut : ', LTcut, ',\n \tLLcut : ', LLcut

	    # plot only events where at least one lepton is non-prompt, and none of the leptons is charge flip
	    # ---> NB: truth req. *only* for closure test!
	    #
            sp_TT  = sp.subprocess(cut=category.cut & self.vardb.getCut('2Lep_NonPromptEvent') & self.vardb.getCut(TTcut), eventweight=weight)
            sp_TL  = sp.subprocess(cut=category.cut & self.vardb.getCut('2Lep_NonPromptEvent') & self.vardb.getCut(TLcut), eventweight=weight)
            sp_LT  = sp.subprocess(cut=category.cut & self.vardb.getCut('2Lep_NonPromptEvent') & self.vardb.getCut(LTcut), eventweight=weight)
            sp_LL  = sp.subprocess(cut=category.cut & self.vardb.getCut('2Lep_NonPromptEvent') & self.vardb.getCut(LLcut), eventweight=weight)
            sp = sp_TT + sp_TL + sp_LT + sp_LL

            print ' \ttruth cut : ', (self.vardb.getCut('2Lep_NonPromptEvent')).cutstr

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
            weight =''
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
            weight= ''
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

	    print 'Zeejets - weight : ', weight, ' TTcut : ', TTcut

	    if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
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
            weight=''
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

	    print 'Zmumujets - weight : ', weight, ' TTcut : ', TTcut

	    if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
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
            weight=''
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

	    print 'Ztautaujets - weight : ', weight, ' TTcut : ', TTcut

	    if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
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
            weight=''
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

	    print 'ZjetsLF - weight : ', weight, ' TTcut : ', TTcut

	    if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
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
            weight=''
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

	    print 'ZjetsHF - weight : ', weight, ' TTcut : ', TTcut

	    if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
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
	    truth_cut = category.cut.removeCut(self.vardb.getCut('2Lep_PurePromptEvent'))
	    truth_cut = truth_cut & self.vardb.getCut('2Lep_ChFlipEvent')

            z = z.subprocess(cut=truth_cut)
            print 'Zjets CF - \ttruth cut : ', self.vardb.getCut('2Lep_ChFlipEvent').cutstr

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
            weight=''
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

	    print 'Wenujets - weight : ', weight, ' TTcut : ', TTcut

	    if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
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
            weight=''
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

	    print 'Wmunujets - weight : ', weight, ' TTcut : ', TTcut

	    if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
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
            weight=''
            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

	    print 'Wtaunujets - weight : ', weight, ' TTcut : ', TTcut

	    if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
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

            #sp = self.base(treename, category, options)
            #return sp
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
            weight=''

	    print 'Top - weight : ', weight, ' TTcut : ', TTcut

            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'
            if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
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
            weight =''

            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

            print 'TTBar - weight name : ', weight, ', TTcut : ', TTcut

            if TTcut is not '':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
            return sp

    # the only difference with the above class is
    # that we plot only events with at least one !prompt lepton, and none charge flip
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
            weight =''

            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'

            print 'TTBarClosureMM \n \tweight name : ', weight, ',\n \tTTcut : ', TTcut, ','

            if TTcut is not '':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

	    # plot only events where at least one lepton is !prompt, and none is charge flip
	    #
            sp = sp.subprocess(cut=category.cut & self.vardb.getCut('2Lep_NonPromptEvent'), eventweight=weight)
            print ' \ttruth cut : ', (self.vardb.getCut('2Lep_NonPromptEvent')).cutstr

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
	    truth_cut = category.cut.removeCut(self.vardb.getCut('2Lep_PurePromptEvent'))
	    truth_cut = truth_cut & self.vardb.getCut('2Lep_ChFlipEvent')

            topcf = topcf.subprocess(cut=truth_cut)
            print 'Top CF - \ttruth cut : ', self.vardb.getCut('2Lep_ChFlipEvent').cutstr

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
            weight=''

	    print 'Diboson - weight : ', weight, ' TTcut : ', TTcut

            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'
            if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
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
	    truth_cut = category.cut.removeCut(self.vardb.getCut('2Lep_PurePromptEvent'))
	    truth_cut = truth_cut & self.vardb.getCut('2Lep_ChFlipEvent')

            dibosoncf = dibosoncf.subprocess(cut=truth_cut)
            print 'DibosonCF - \ttruth cut : ', self.vardb.getCut('2Lep_ChFlipEvent').cutstr

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
            weight=''

	    print 'HtoZZ - weight : ', weight, ' TTcut : ', TTcut

            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'
            if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))
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
            return (diboson + top + htozz + ttbarw+ttbarz)

        def __call__(self, treename='physics', category=None, options={}):

            diboson = self.parent.procmap['Diboson'](treename, category, options)
            top = self.parent.procmap['Top'](treename, category, options)
            htozz = self.parent.procmap['HtoZZ'](treename, category, options)
            ttbarw = self.parent.procmap['TTBarW'](treename, category, options)
            ttbarz = self.parent.procmap['TTBarZ'](treename, category, options)
            return (diboson + top + htozz + ttbarw+ttbarz)


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
            weight=''

            print 'ChargeFlipMC - weight : ', weight, ' TTcut : ', TTcut

            if self.parent.channel=='TwoLepSS' or self.parent.channel=='ThreeLep':
                TTcut='TT'
            if TTcut != '':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))

            # plot only events where at least one lepton is charge flip. Remove req. where both lep must be prompt
            #
            truth_cut = category.cut.removeCut(self.vardb.getCut('2Lep_PurePromptEvent'))
            truth_cut = truth_cut & self.vardb.getCut('2Lep_ChFlipEvent')

            sp = sp.subprocess(cut=truth_cut, eventweight=weight)
            print ' \ttruth cut : ', self.vardb.getCut('2Lep_ChFlipEvent').cutstr

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
            weight=''
            if self.parent.channel=='TwoLepSS':
                TTcut='TT'
                kfactors = {'Norm_corr': (1.0, 0.1),}
            elif self.parent.channel=='ThreeLep':
                TTcut='TT'
                kfactors = {'Norm_corr': (1.0, 0.1),}
            sp = self.parent.applyKfactor(sp, category, kfactors, kf_opts)#this command is used to apply the kfactor
            if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))
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
            TTcut=''
            weight=''
            if self.parent.channel=='TwoLepSS':
                TTcut='TT'
            elif self.parent.channel=='ThreeLep':
                TTcut='TT'
            if TTcut!='':
                sp = sp.subprocess(cut=self.vardb.getCut(TTcut))
            sp = sp.subprocess(cut=category.cut, eventweight=weight)
            return sp
'''
