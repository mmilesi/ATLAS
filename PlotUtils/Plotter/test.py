
    class FakesABCD(Process):

	latexname = 'Fakes \theta method'
        colour = kCyan - 9

        theta_el = {
            'El': (1.0, 0.0),
        }
        theta_mu = {
            'Mu': (1.0, 0.0),
        }

        def base(self, treename='physics', category=None, options={}):
            inputgroup = [
                    ('Data', 'physics_Main'),
                ]
            trees = self.inputs.getTrees(treename, inputgroup)
            return self.subprocess(trees=trees)	    


        def calcTheta(self, sp_num, sp_denom, stream=None, options={}):

            if not sp_denom.number():
                print "ERROR: Cannot calculate theta transfer factor!"

            theta = (sp_num/sp_denom).numberstats()
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

            print 'FakesABCD \n '

            systematics = options.get('systematics', None)
            direction = options.get('systematicsdirection', 'UP')
            systname_opts = {}
            if systematics and systematics.name == 'SystName':
                systname_opts['systematics'] = True
                systname_opts['systematicsdirection'] = direction
	    
	    TTcut =''
            TTcut =''
            TLcut =''
            LTcut =''
            TelLmucut =''
            LelTmucut =''
            TmuLelcut =''
            LmuTelcut =''
            weight =''
            if self.parent.channel=='TwoLepSS':
                TTcut  = self.vardb.getCut('TT')
                TLcut  = self.vardb.getCut('TL')
                LTcut  = self.vardb.getCut('LT')
                TelLmucut = self.vardb.getCut('OF_Event') & Cut('TelLmu',  '( isTelLmu == 1 )')
                LelTmucut = self.vardb.getCut('OF_Event') & Cut('LelTmu',  '( isLelTmu == 1 )')
                TmuLelcut = self.vardb.getCut('OF_Event') & Cut('TmuLel',  '( isTmuLel == 1 )')
                LmuTelcut = self.vardb.getCut('OF_Event') & Cut('LmuTel',  '( isLmuTel == 1 )')

	    # take the base suprocess (i.e, data)
	    #
	    sp = self.base(treename, category, options)
	    
	    # define region B
	    #
	    sp_B = sp.subprocess(cut=( ( TLCut | LTCut) & self.vardb.getCut('NJet2L') ))
	    
	    # start defining region C and D
	    #
	    sp_C = sp.subprocess(cut=( TTCut & self.vardb.getCut('LowJetCR'))) 
	    sp_D = sp.subprocess(cut=( ( TLCut | LTCut) & self.vardb.getCut('LowJetCR') ))
    
            # get a list of stuff to subtract for region C and D (i.e, prompt MC and charge flips)
	    #
            sublist = [ item for item in self.parent.sub_backgrounds ]
            
	    # ... and now subtract!
	    #
            for sample in sublist:
                sp_C = sp_C - self.parent.procmap[sample].base(treename, category, options) # here it is important to have used base otherwise at the prompt sample there would be already applied the selection specified in the category __call__ of that sample i.e. for example also the iso-iso cut which is orthogonal to the cuts done in this sample which are fake iso or fake fake.
                sp_D = sp_D - self.parent.procmap[sample].base(treename, category, options) # here it is important to have used base otherwise at the prompt sample there would be already applied the selection specified in the category __call__ of that sample i.e. for example also the iso-iso cut which is orthogonal to the cuts done in this sample which are fake iso or fake fake.
   
            # derive theta factors for e and mu
	    #
	    if ( ("ElEl_Event") in category.cut.cutname or ("OF_Event") in category.cut.cutname ):
	      sp_C_el = sp_C.subprocess(cut=self.vardb.getCut('ElEl_Event'))
	      sp_D_el = sp_D.subprocess(cut=self.vardb.getCut('ElEl_Event'))
 	      self.theta_el = self.calcTheta(sp_C_el,sp_D_e,stream='El')
	    if ( ("MuMu_Event") in category.cut.cutname or ("OF_Event") in category.cut.cutname  ):
	      sp_C_mu = sp_C.subprocess(cut=self.vardb.getCut('MuMu_muvent'))
	      sp_D_mu = sp_D.subprocess(cut=self.vardb.getCut('MuMu_muvent'))
 	      self.theta_mu = self.calcTheta(sp_C_mu,sp_D_mu,stream='Mu')
	    
	    print ('(self.theta_el[\'El\'])[0] = {0} - (self.theta_el[\'El\'])[1] = {1}'.format( (self.theta_el['El'])[0], (self.theta_el['El'])[1] ) )
	    print ('(self.theta_mu[\'Mu\'])[0] = {0} - (self.theta_mu[\'Mu\'])[1] = {1}'.format( (self.theta_mu['Mu'])[0], (self.theta_mu['Mu'])[1] ) )
    
	    # characterise region B, depending on which flavour comp we are looking at
	    #
	    if ("ElEl_Event") in category.cut.cutname:
	      sp_B = self.applyThetaFactor(sp_B,self.theta_el)
	    elif ("MuMu_Event") in category.cut.cutname:
	      sp_B = self.applyThetaFactor(sp_B,self.theta_mu)
            elif ("OF_Event") in category.cut.cutname:
	      sp_B_Lel = sp.subprocess(cut=( LelTmucut | TmuLelcut ))
	      sp_B_Lmu = sp.subprocess(cut=( TelLmucut | LmuTelcut ))
	      sp_B = self.applyThetaFactor(sp_B_Lmu,self.theta_mu) + self.applyThetaFactor(sp_B_Lel,self.theta_el)

            return sp_B

