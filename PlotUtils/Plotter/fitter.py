# encoding: utf-8
'''
fracfitter.py

description:

'''
__author__    = "Will Davey"
__email__     = "will.davey@cern.ch"
__created__   = "2013-02-27"
__copyright__ = "Copyright 2013 Will Davey"
__license__   = "GPL http://www.gnu.org/licenses/gpl.html"



## modules
import ROOT
from ROOT import RooFit
#from ROOT.RooFit import RooArgSet, RooArgList, Minos, Extended, Save, Components, LineColor, LineStyle
from math import sqrt
import utils 
from array import array
import numpy

#grandom = ROOT.TRandom3(1000)
grandom = ROOT.TRandom3(0)
##ROOT.gROOT.SetBatch()

# - - - - - - - - - - - class defs  - - - - - - - - - - - - #
#------------------------------------------------------------
class Model():
    '''
    description of BaseModel 
    '''
    #____________________________________________________________
    def __init__(self,
            xvar = None,
            h_data = None,
            h_sig = None,
            h_sig_1p = None,
            h_sig_mp = None,
            h_bkg = None,
            h_anti = None,
            n_ev_min  = 0,
            n_ev_max = 1000000,
            tag = 'noid',
            color_model = ROOT.kBlack,
            color_sig = ROOT.kRed,
            color_sig_1p = ROOT.kRed,
            color_sig_mp = ROOT.kRed,
            color_bkg = ROOT.kBlue,
            color_anti = ROOT.kOrange+1,
            style_model = ROOT.kSolid,
            style_sig = ROOT.kDashed,
            style_sig_1p = ROOT.kSolid,
            style_sig_mp = ROOT.kDashed,
            style_bkg = ROOT.kDashed,
            style_anti = ROOT.kDashed,
            bkg_const_model  = None, ## to set fake-factors
            set_const_anti = True,
            stat_limit = 100.,
            fix_r1p3p = False,
            fout = None,
            quiet = False,
            fit_min_bin = None,
            fit_max_bin = None,
            draw_true_hists = False,
            plot_toy_fits = False,
            ):
        self.xvar      = xvar
        self.n_ev_min  = n_ev_min
        self.n_ev_max  = n_ev_max
        self.tag       = tag  # incorporate this as var in names

        ## histograms
        self.h_orig_data = h_data
        self.orig_templates = {}
        self.orig_templates['sig']    = h_sig
        self.orig_templates['sig_1p'] = h_sig_1p
        self.orig_templates['sig_mp'] = h_sig_mp
        self.orig_templates['bkg']    = h_bkg
        self.orig_templates['anti']   = h_anti
        self.h_orig_total = None

        ## additional hists derived from input hists
        ## used for fitting / pseudo experiments 
        ## (filled in initialise)
        self.h_curr_data = None 
        self.orig_raw_templates = {}
        self.orig_wei_templates = {}
        self.curr_raw_templates = {}

        self.color_model  = color_model
        self.color_sig    = color_sig
        self.color_sig_1p = color_sig_1p
        self.color_sig_mp = color_sig_mp
        self.color_bkg    = color_bkg
        self.color_anti   = color_anti
        self.style_model  = style_model
        self.style_sig    = style_sig
        self.style_sig_1p = style_sig_1p
        self.style_sig_mp = style_sig_mp
        self.style_bkg    = style_bkg
        self.style_anti   = style_anti

        self.bkg_const_model = bkg_const_model
        self.set_const_anti = set_const_anti

        self.stat_limit = stat_limit
        self.fix_r1p3p = fix_r1p3p
        self.fout = fout
        self.quiet = quiet

        self.fit_min = 1.
        self.fit_max = 20.
        self.fit_min_bin = fit_min_bin
        self.fit_max_bin = fit_max_bin
        self.draw_true_hists = draw_true_hists
        self.plot_toy_fits = plot_toy_fits

        ## fit results
        self.fitter       = None
        self.fit_pars     = {}
        self.fit_results  = {}
        self.fit_hists    = {}
        self.nfit_tot     = None
      

        ## corrections from toy MC
        self.toy_results = {}



        assert h_sig or (h_sig_1p and h_sig_mp), "ERROR - must config either incl or 1p/mp signal"



    #____________________________________________________________
    def initialise(self):

        if self.orig_templates['sig']:
            print 'Original total Ztautau events: %.3f'%(integral(self.orig_templates['sig']))
        if self.orig_templates['sig_1p']:
            print 'Original total Ztautau [1p] events: %.3f'%(integral(self.orig_templates['sig_1p']))
        if self.orig_templates['sig_mp']:
            print 'Original total Ztautau [Mp] events: %.3f'%(integral(self.orig_templates['sig_mp']))

        ## create merged signal template if fixing 1p3p ratio
        if self.orig_templates['sig_1p'] and self.fix_r1p3p:
            h_orig_sig = self.orig_templates['sig_1p'].Clone('h_sig')
            h_orig_sig.Add(self.orig_templates['sig_mp'])
            self.orig_templates['sig'] = h_orig_sig
        
        ## get list of active samples
        samples = []
        if self.orig_templates['sig']: samples += ['sig']
        else:                          samples += ['sig_1p','sig_mp']
        samples += ['bkg','anti']
        #samples += ['bkg'] 
        self.samples = samples

        ## make total hist
        for s in samples:
            h = self.orig_templates[s]
            if not self.h_orig_total: self.h_orig_total = h.Clone('h_orig_total')
            else: self.h_orig_total.Add(h)

        ## split templates into raw and weight components
        for s, h in self.orig_templates.items():
            h_raw = h_wei = None
            if h: h_raw, h_wei = split_hist_raw_weight_components(h)
            self.orig_raw_templates[s] = h_raw
            self.orig_wei_templates[s] = h_wei

        ## make copies of raw templates to pass to fit
        ## (to allow fluctuating in pseudo experiments)
        for s, h in self.orig_raw_templates.items():
            h_clone = None
            if h: h_clone = h.Clone('%s_curr'%h.GetName())
            self.curr_raw_templates[s] = h_clone

        ## make copy of data to pass to fit
        ## (to allow creation of pseudo-data)
        self.h_curr_data = self.h_orig_data.Clone('%s_curr'%self.h_orig_data.GetName())
        
        ## determine fit range
        self.determine_fit_range()

    #____________________________________________________________
    def determine_fit_range(self):
        print 'setting fit range...'

        if self.fit_min_bin != None:
            min_bin = self.fit_min_bin
            max_bin = self.fit_max_bin
        else:
            #print 'determining fit range'
            max_bin = self.h_curr_data.GetNbinsX()
            min_bin = 0
            for i in range(1,self.h_curr_data.GetNbinsX()-1):
                n = self.h_curr_data.GetBinContent(i)
                print 'bin%d: %.1f' % (i,n)
                if n > 0 and min_bin == 0:
                    min_bin = i
                    break
            for i in range(min_bin+1,self.h_curr_data.GetNbinsX()):
                n = self.h_curr_data.GetBinContent(i)
                print 'bin%d: %.1f' % (i,n)
                if n >= 0 and n < self.stat_limit:
                    max_bin = i-1
                    if n == 0 and max_bin-min_bin >= 2:
                        max_bin -= 1
                    break
            #max_bin = 12
            if max_bin == min_bin:
                print 'unable to find good fit range! EXIT'
                exit()
            print 'determined fit bin range: ', [min_bin,max_bin]

        self.fit_min = self.h_curr_data.GetXaxis().GetBinLowEdge(min_bin)
        fit_max_low = self.h_curr_data.GetXaxis().GetBinLowEdge(max_bin)
        self.fit_max = self.h_curr_data.GetXaxis().GetBinUpEdge(max_bin)
        self.fit_min_bin = min_bin 
        self.fit_max_bin = max_bin
        print 'determined fit range: [%.2f, %.2f - %.2f]' %(self.fit_min,fit_max_low,self.fit_max)

    
    #____________________________________________________________
    def init_fitter(self,toy=False):
        if self.fitter: del(self.fitter)

        ## create object array to pass templates to TFractionFitter
        self.mc_arr = mc_arr = ROOT.TObjArray()
        for s in self.samples:
            self.fit_pars[s] = mc_arr.GetEntries() 
            mc_arr.Add(self.curr_raw_templates[s])

        ## instantiate fitter
        opt = 'Q' if toy or self.quiet else ''
        self.fitter = fitter = ROOT.TFractionFitter(self.h_curr_data,self.mc_arr,opt)
        for s in self.samples:
            fitter.SetWeight(self.fit_pars[s],self.orig_wei_templates[s])

        self.fitter.SetRangeX(self.fit_min_bin,self.fit_max_bin)
        self.fitter.ExcludeBin(1)

        print 'initialised fitter'
        print 'templates:'
        for s in self.samples:
            print '%10s: %d'%(s,self.fit_pars[s])

    #____________________________________________________________
    def config_fit_params(self,toy=None):
        if toy==None: toy=False
        if not self.quiet: print 'configuring fit parameters...'
        fitter = self.fitter
        f = fitter.GetFitter()

        ## update fluctuated hists
        if toy:
            fitter.SetData(self.h_curr_data)
            for s in self.samples:
                fitter.SetMC(self.fit_pars[s],self.curr_raw_templates[s])
                fitter.SetWeight(self.fit_pars[s],self.orig_wei_templates[s])


        ## determine fractions
        n_data = self.ndata_curr()
        fractions = {}
        if self.orig_templates['sig_1p'] and not self.fix_r1p3p: 
            fractions['sig_1p'] = self.nsamp_curr('sig_1p') / n_data
            fractions['sig_mp'] = self.nsamp_curr('sig_mp') / n_data
        else:
            fractions['sig'] = self.nsamp_curr('sig') / n_data
        if self.bkg_const_model and not toy:
            """
            constrain bkg normalisation from other region
            using fake-factor.
            """
            fake_factor = integral(self.orig_templates['bkg']) / integral(self.bkg_const_model.orig_templates['bkg'])
            nbkg_fit = self.bkg_const_model.nsamp_fit('bkg')
            ebkg_fit = self.bkg_const_model.ensamp_fit('bkg')
            nbkg = nbkg_fit * fake_factor
            ebkg = ebkg_fit * fake_factor
            fractions['bkg'] = nbkg / n_data
            fractions['ebkg'] = ebkg / n_data

        else:
            fractions['bkg'] = self.nsamp_curr('bkg') / n_data
        if 'anti' in self.samples:
            # prevent setting to 0, b/c it breaks fit
            fractions['anti'] = max(0.001, self.nsamp_curr('anti') / n_data )


        ## release parameters
        for s in self.samples: f.ReleaseParameter(self.fit_pars[s])

        ## set parameters
        default_step_size = 0.001
        for s in self.samples:
            fname = 'f%s'%s
            frac = fractions[s]
            ipar = self.fit_pars[s]
            f.SetParameter(ipar,fname,frac,default_step_size,0,0)

        ## constrain parameters
        for s in self.samples:
            fitter.Constrain(self.fit_pars[s],0.,1.)
            if not self.quiet: print 'constrain %s frac: [0.,1.], init val: %f' %(s,fractions[s]) 

        ## fix parameters
        if self.bkg_const_model:
            #f_min = fractions['bkg']-fractions['ebkg']/2.
            #f_max = fractions['bkg']+fractions['ebkg']/2.
            #fitter.Constrain(self.fit_pars['bkg'],f_min,f_max)
            #if not self.quiet: print 'constrain bkg frac: %f +- %f'%(fractions['bkg'],fractions['ebkg']/2.)
            f.FixParameter(self.fit_pars['bkg'])
            if not self.quiet: print 'fixed bkg frac: ', fractions['bkg']
        if 'anti' in self.samples and self.set_const_anti:
            f.FixParameter(self.fit_pars['anti'])
            if not self.quiet: print 'fixed anti frac: ', fractions['anti']





    #____________________________________________________________
    def fit(self,save=True):
        """
        real fit to data
        """
        ## ensure current hists are original ones
        ## not fluctuated ones from toys
        self.reset_hists()
        self.init_fitter(toy=False)
        self.config_fit_params(toy=False) 
        self.perform_fit(save=True)


    #____________________________________________________________
    def toy_fit(self,save=True):
        """
        toy fit
        """
        self.config_fit_params(toy=True) 
        self.perform_fit(save=False)

    #____________________________________________________________
    def perform_fit(self,save=True):
      
        if not self.quiet:
            print 
            print 
            print '-------------->>>>>>> starting fit <<<<<<-----------------'
            print
            print


        if not self.quiet: print 'clearing previous fit results...'
        ## clear previous fit results
        self.chi2        = None
        self.ndf         = None
        self.nfit_tot    = None 
        self.fit_results = {}
        
        ## clear previous fit hists
        for s in self.samples: self.fit_hists[s] = None
        self.fit_hists['model']  = None

        if not self.quiet: print 'performing minimisation...'
        fitter = self.fitter 
        self.fit_status = fitter.Fit()
        if not self.fit_status == 0: 
            print 'ERROR - failure in fit'
            return

        ## fit stats
        self.chi2 = fitter.GetChisquare()
        self.ndf = fitter.GetNDF()

        ## save results manually since shitty
        ## TFractionFitter affected by next fit
        for i in xrange(self.mc_arr.GetEntries()):
            val = ROOT.Double()
            err = ROOT.Double()
            fitter.GetResult(i,val,err)
            self.fit_results[i] = [val,err]
            #print 'fit res par%d: %.4f +- %.4f' % (i,val,err)


        if not self.quiet:
            fsum = sum( d[0] for d in self.fit_results.values() )
            print 'fsum: ', fsum

        self.fit_hists['model'] = fitter.GetPlot()
        self.nfit_tot  = self.fit_hists['model'].Integral()

        for s in self.samples:
            htemp = fitter.GetMCPrediction(self.fit_pars[s])
            htemp.Multiply(self.orig_wei_templates[s])
            self.fit_hists[s] = htemp 

        if not self.quiet:
            print 
            print 
            print '-------------->>>>>>> fit SUCCESS! <<<<<<-----------------'
            print
            print
        
        

    #____________________________________________________________
    def plot(self):

        line_width = 3

        cname = 'c_%s'%self.tag
        c = ROOT.TCanvas(cname,cname,700,700)
        c.cd()
        fr = c.DrawFrame(0.,0.,20.,1.2*self.h_curr_data.GetMaximum(),';Number of Tracks;Events')
        
        ## construct legend
        leg = ROOT.TLegend(0.6,0.2,0.9,0.5)
        leg.SetBorderSize(0)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)

        ## get contributions
        h_bkg = self.fit_hists['bkg']
        h_bkg.Scale(self.nsamp_fit('bkg',corr=False) / h_bkg.Integral())
        if 'anti' in self.samples:
            h_anti = self.fit_hists['anti']
            h_anti.Scale(self.nsamp_fit('anti',corr=False) / h_anti.Integral())
        if self.orig_templates['sig_1p']:
            if self.fix_r1p3p: 
                h_sig_1p = self.orig_templates['sig_1p'].Clone()
                h_sig_mp = self.orig_templates['sig_mp'].Clone()
            else:
                h_sig_1p = self.fit_hists['sig_1p'] 
                h_sig_mp = self.fit_hists['sig_mp']

            h_sig_1p.Scale(self.nsamp_fit('sig_1p',corr=False) / h_sig_1p.Integral())
            h_sig_mp.Scale(self.nsamp_fit('sig_mp',corr=False) / h_sig_mp.Integral())
        else: 
            h_sig = self.fit_hists['sig']
            h_sig.Scale(self.nsamp_fit('sig',corr=False) / h_sig.Integral())

        ## draw error on total model
        ## TODO: fix this up to use orig hists
        h_total = None
        for s in self.samples:
            if not h_total: 
                h_total = self.orig_templates[s].Clone('h_total')
                h_total.Scale( self.nsamp_fit(s) / h_total.Integral() )
            else:
                htemp = self.orig_templates[s]
                h_total.Add( htemp, self.nsamp_fit(s) / htemp.Integral() )
        
        h_total.SetFillColor(ROOT.kRed)
        h_total.SetLineColor(ROOT.kRed)
        h_total.SetLineStyle(1)
        h_total.SetLineWidth(0)
        h_total.SetMarkerSize(0)
        h_total.Draw("SAME,E2")
            
        
        fsum = 0.
        for s in self.samples: fsum += self.fsamp_fit(s)
        print 'fsum: ', fsum

            
        ## total model central value from fit
        #self.fit_hists['model'].Scale(self.ndata() / self.nhist(self.fit_hists['model']))
        #self.fit_hists['model'].Scale(self.ndata_curr() / self.nhist(self.fit_hists['model']))
        self.fit_hists['model'].SetLineWidth(line_width)
        self.fit_hists['model'].Draw("SAME")
        leg.AddEntry(self.h_curr_data,'Data','PL')
        leg.AddEntry(self.fit_hists['model'],'Model','L')
        leg.AddEntry(h_total,'Model (stat.)','F')

        ## draw data
        self.h_curr_data.Draw("SAME")
        
        print 'nfit: ',       self.ntot_fit()
        print 'nfit(corr): ', self.ntot_fit()
        print 'h_mod: ', self.fit_hists['model'].Integral()
        print 'h_tot: ', h_total.Integral()


        # draw bkg
        h_bkg.SetLineColor(self.color_bkg)
        h_bkg.SetLineStyle(self.style_bkg)
        h_bkg.SetLineWidth(line_width)
        h_bkg.Draw("SAME,HIST")
        leg.AddEntry(h_bkg,'Jet','L')
        if self.draw_true_hists: self.orig_templates['bkg'].Draw("SAME,HIST")


        # draw anti
        if 'anti' in self.samples:
            h_anti.SetLineColor(self.color_anti)
            h_anti.SetLineStyle(self.style_anti)
            h_anti.SetLineWidth(line_width)
            h_anti.Draw("SAME,HIST")
            leg.AddEntry(h_anti,'Lep','L')
            if self.draw_true_hists: self.orig_templates['anti'].Draw("SAME,HIST")

        ## 1p3p split signal
        if self.orig_templates['sig_1p']:
            h_sig_1p.SetLineColor(self.color_sig_1p)
            h_sig_1p.SetLineStyle(self.style_sig_1p)
            h_sig_1p.SetLineWidth(line_width)
            
            h_sig_mp.SetLineColor(self.color_sig_mp)
            h_sig_mp.SetLineStyle(self.style_sig_mp)
            h_sig_mp.SetLineWidth(line_width)
            
            h_sig_1p.Draw("SAME,HIST")
            h_sig_mp.Draw("SAME,HIST")
            leg.AddEntry(h_sig_1p,'Tau (1p)','L')
            leg.AddEntry(h_sig_mp,'Tau (mp)','L')
            if self.draw_true_hists: self.orig_templates['sig_1p'].Draw("SAME,HIST")
            if self.draw_true_hists: self.orig_templates['sig_mp'].Draw("SAME,HIST")
        else: 
            h_sig = self.fit_hists['sig']
            h_sig.SetLineColor(self.color_sig)
            h_sig.SetLineStyle(self.style_sig)
            h_sig.SetLineWidth(line_width)
            h_sig.Draw("SAME,HIST")
            leg.AddEntry(h_sig,'Tau','L')
            if self.draw_true_hists: self.orig_templates['sig'].Draw("SAME,HIST")

        
        leg.Draw()

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextFont(42)
        latex.DrawLatex(0.4,0.85,self.tag)
        latex.DrawLatex(0.4,0.75,'#chi^{2}/NDF = %.1f / %d' % (self.chi2,self.ndf) ) 

        c.Update()
        c.SaveAs('fit_%s.eps'%self.tag)

        if self.fout: utils.save_object(c,self.fout)

    
    #____________________________________________________________
    def finalise(self):
        #self.fitter.Delete()
        del(self.fitter)
        #if self.fout: 
        #    utils.save_object(self.fit_result,self.fout)


    #____________________________________________________________
    def summary(self):
        print
        print 'fit stats'
        print 'chi2: ', self.chi2
        print 'NDF: ', self.ndf
        if not self.orig_templates['sig_1p']:
            print '%s fractions:' % self.tag
            print 'sig: %.3f+-%.3f' % (self.fsamp_fit('sig'),self.efsamp_fit('sig'))
            print 'bkg: %.3f+-%.3f' % (self.fsamp_fit('bkg'),self.efsamp_fit('bkg'))
            if 'anti' in self.samples:
                print 'anti: %.3f+-%.3f' % (self.fsamp_fit('anti'),self.efsamp_fit('anti'))


            print '%s corr Nmc: %.1f+-%.1f, Ndata: %.1f+-%.1f (%.3f%%)' % (self.tag,
                                                                  self.nsamp_curr('sig'),
                                                                  self.ensamp_curr('sig'),
                                                                  self.nsamp_fit('sig'),
                                                                  self.ensamp_fit('sig'),
                                                                  self.ensamp_fit('sig')/self.nsamp_fit('sig')*100.,
                                                                  )

            print '%s uncorr Nmc: %.1f+-%.1f, Ndata: %.1f+-%.1f (%.3f%%)' % (self.tag,
                                                                  self.nsamp_curr('sig'),
                                                                  self.ensamp_curr('sig'),
                                                                  self.nsamp_fit('sig',corr=False),
                                                                  self.ensamp_fit('sig',corr=False),
                                                                  self.ensamp_fit('sig',corr=False)/self.nsamp_fit('sig',corr=False)*100.,
                                                                  )

        else:
            print '%s corr 1p Nmc: %.1f+-%.1f, Ndata: %.1f+-%.1f (%.3f%%)' % (self.tag,
                                                                  self.nsamp_curr('sig_1p'),
                                                                  self.ensamp_curr('sig_1p'),
                                                                  self.nsamp_fit('sig_1p'),
                                                                  self.ensamp_fit('sig_1p'),
                                                                  self.ensamp_fit('sig_1p')/self.nsamp_fit('sig_1p')*100.,
                                                                  )
            print '%s corr mp Nmc: %.1f+-%.1f, Ndata: %.1f+-%.1f (%.3f%%)' % (self.tag,
                                                                  self.nsamp_curr('sig_mp'),
                                                                  self.ensamp_curr('sig_mp'),
                                                                  self.nsamp_fit('sig_mp'),
                                                                  self.ensamp_fit('sig_mp'),
                                                                  self.ensamp_fit('sig_mp')/self.nsamp_fit('sig_mp')*100.,
                                                                  )


            print '%s uncorr 1p Nmc: %.1f+-%.1f, Ndata: %.1f+-%.1f (%.3f%%)' % (self.tag,
                                                                  self.nsamp_curr('sig_1p'),
                                                                  self.ensamp_curr('sig_1p'),
                                                                  self.nsamp_fit('sig_1p',corr=False),
                                                                  self.ensamp_fit('sig_1p',corr=False),
                                                                  self.ensamp_fit('sig_1p',corr=False)/self.nsamp_fit('sig_1p',corr=False)*100.,
                                                                  )
            print '%s uncorr mp Nmc: %.1f+-%.1f, Ndata: %.1f+-%.1f (%.3f%%)' % (self.tag,
                                                                  self.nsamp_curr('sig_mp'),
                                                                  self.ensamp_curr('sig_mp'),
                                                                  self.nsamp_fit('sig_mp',corr=False),
                                                                  self.ensamp_fit('sig_mp',corr=False),
                                                                  self.ensamp_fit('sig_mp',corr=False)/self.nsamp_fit('sig_mp',corr=False)*100.,
                                                                  )


    
    #____________________________________________________________
    def randomise_hists(self):
        ## fluctuate templates and data
        for s in self.samples:
            fluctuate_hist_gaus(self.orig_raw_templates[s],self.curr_raw_templates[s])
            #if self.curr_raw_templates[s].Integral()<=0.:
            #    print 'NULL template: %s'%s
            #print '%s sum: %.1f'%(s,self.curr_raw_templates[s].Integral())
        fluctuate_hist_pois(self.h_orig_total,self.h_curr_data)

    #____________________________________________________________
    def prepare_toy_hists(self):
        ## setup hists for toys
        for s in self.samples:
            reset_hist(self.orig_raw_templates[s],self.curr_raw_templates[s])
        reset_hist(self.h_orig_total,self.h_curr_data)

    #____________________________________________________________
    def reset_hists(self):
        for s in self.samples: 
            reset_hist(self.orig_raw_templates[s],self.curr_raw_templates[s])
        reset_hist(self.h_orig_data,self.h_curr_data)       

    #____________________________________________________________
    def mc_study(self):
        """
        perform psedudo-experiments by generating random 
        data and models (using statistical bin uncertainty)

        pseudo data: 
            - sum all model hists
            - then randomise
        model: 
            - randomise each input hist individually
        """
        samples = self.samples

        ntrials = 10000
        if self.plot_toy_fits: ntrials = 10
        
        ## save options before mc study
        tag = self.tag 
        quiet = self.quiet
        self.quiet = True

        toy_arrays = {} 
        for s in samples: 
            if not s in toy_arrays: toy_arrays[s] = {}
            toy_arrays[s]['mean']  = []
            toy_arrays[s]['error'] = []
            toy_arrays[s]['pull']  = []
            toy_arrays[s]['diff']  = []
            toy_arrays[s]['mc']  = []
        data_array = []

        ## initialise toy fitter
        ## VERY important to prepare hists first
        ## before initialising, so TFractionFitter
        ## is not initialise with the real data
        self.prepare_toy_hists()
        self.init_fitter(toy=True)

        for i in xrange(ntrials):
            if i%100==0: print 'trial ',i
            
            #print
            #print 'trial%d summary: ' % i

            self.tag = '%s_trial%d'% (tag,i)
            self.randomise_hists()
            self.toy_fit()
            if not self.fit_status == 0: continue
            if self.plot_toy_fits: self.plot()

           
            temp_means = {}
            temp_errors = {}
            temp_pulls = {} 
            temp_diffs = {} 
            temp_mcs = {} 
            has_zero = False
            n_tot_mc = self.ntot_orig()
            #print 'n_tot_mc: %.1f, h_orig_tot.int: %.1f'% (n_tot_mc,self.h_orig_total.Integral())

            n_tot_fit = self.ntot_fit()
            #print 'n_tot_fit: %.1f, n_data: %.1f' % (n_tot_fit,self.ndata_curr())
            for s in samples: 
                n_mc   = self.nsamp_orig(s)
                n_fit  = self.nsamp_fit(s)
                en_fit = self.ensamp_fit(s)
                pull   = (n_fit-n_mc)/en_fit if en_fit else 0.0
                diff   = n_fit - n_mc
                temp_means[s] = n_fit
                temp_errors[s] = en_fit
                temp_pulls[s] = pull
                temp_diffs[s] = diff
                temp_mcs[s] = self.nsamp_curr(s)
                f_fit = self.fsamp_fit(s)
                #if f_fit == 0.: has_zero = True 
                if f_fit < 0.0000001: has_zero = True            
                f_mc = n_mc / n_tot_mc if n_tot_mc else 0.0
                f_fit = n_fit / n_tot_fit if n_tot_fit else 0.0
                #print '%s, mc: %.1f, fit: %.1f, fmc: %.4f, ffit: %.4f' % (s,n_mc,n_fit,f_mc,f_fit)

            #print 'data, mc: %.1f, fit: %.1f' % (n_tot_mc,n_tot_fit)
            ## remove cases where any component is fit to 0
            ## argument is that we would not take this 
            ## result if we got it in data
            ## probably should try to do something 
            ## better in future
            if not has_zero:      
                for s in samples: 
                    toy_arrays[s]['mean'].append(temp_means[s])
                    toy_arrays[s]['error'].append(temp_errors[s])
                    toy_arrays[s]['pull'].append(temp_pulls[s])
                    toy_arrays[s]['diff'].append(temp_diffs[s])
                    toy_arrays[s]['mc'].append(temp_mcs[s])
                    data_array.append(self.ndata_curr())
            else:
                print 'ERROR - component fit to zero'


        ## restore to original state before toys 
        self.reset_hists()
        self.tag = tag
        self.quiet = quiet

        ## set corrections from toy study
        filename = 'toy_%s.root'%(self.tag)
        for s in samples: 
            a_mean  = toy_arrays[s]['mean']
            a_error = toy_arrays[s]['error']
            a_pull  = toy_arrays[s]['pull']
            a_diff  = toy_arrays[s]['diff']
            a_mc    = toy_arrays[s]['mc']
 
            if not s in self.toy_results: self.toy_results[s] = {}
            self.toy_results[s]['meanm']  = numpy.mean(a_mean)
            self.toy_results[s]['meane']  = numpy.std(a_mean)
            self.toy_results[s]['errorm']  = numpy.mean(a_error)
            self.toy_results[s]['errore']  = numpy.std(a_error)
            self.toy_results[s]['pullm']  = numpy.mean(a_pull)
            self.toy_results[s]['pulle']  = numpy.std(a_pull)
            
            ## create plots
            h_mean  = create_mean_hist(s)
            h_error = create_error_hist(s)
            h_pull  = create_pull_hist(s)
            h_diff  = create_diff_hist(s)
            h_mc    = create_mc_hist(s)
            for v in a_mean: h_mean.Fill(v)
            for v in a_error: h_error.Fill(v)
            for v in a_pull: h_pull.Fill(v)
            for v in a_diff: h_diff.Fill(v)
            for v in a_mc  : h_mc.Fill(v)
            utils.save_object(h_mean,filename)
            utils.save_object(h_error,filename)
            utils.save_object(h_pull,filename)
            utils.save_object(h_diff,filename)
            utils.save_object(h_mc,filename)

        h_mc_data    = create_mc_hist('data')
        for v in data_array: h_mc_data.Fill(v)
        utils.save_object(h_mc_data,filename)


        for isamp in xrange(len(samples)):
            s1 = samples[isamp]
            for isamp2 in xrange(len(samples)):
                if not isamp2 < isamp: continue
                s2 = samples[isamp2]
                h = create_2d_mean_hist(s1,s2)
                for ns1,ns2 in zip(toy_arrays[s1]['mean'],
                                   toy_arrays[s2]['mean']):
                    h.Fill(ns1,ns2)
                utils.save_object(h,filename)
       
        for s in samples:
            h = create_2d_mean_hist('%s_mc'%s,'%s_fit'%s)
            for ns1,ns2 in zip(toy_arrays[s]['mc'],
                               toy_arrays[s]['mean']):
                h.Fill(ns1,ns2)
            utils.save_object(h,filename)



        f = utils.open_file(filename)
        f.Close()

    #____________________________________________________________
    def nhist(self,h):
        return integral(h)
    #____________________________________________________________
    def enhist(self,h):
        return integral_error(h)
    #____________________________________________________________
    def nsamp_orig(self,sample):
        """
        wrapper to hists for member samples
        """
        h = self.orig_templates[sample]
        return self.nhist(h)
    #____________________________________________________________
    def ensamp_orig(self,sample):
        """
        wrapper to hists for member samples
        """
        h = self.orig_templates[sample]
        return self.enhist(h)
    
    #____________________________________________________________
    def nsamp_curr(self,sample):
        """
        wrapper to hists for member samples
        """
        h_raw = self.curr_raw_templates[sample]
        h_wei = self.orig_wei_templates[sample]
        n = 0.
        for i in xrange(0,h_raw.GetNbinsX()+1):
            n+= h_raw.GetBinContent(i) * h_wei.GetBinContent(i)
        return n

    #____________________________________________________________
    def ensamp_curr(self,sample):
        """
        wrapper to hists for member samples
        """
        h_raw = self.curr_raw_templates[sample]
        h_wei = self.orig_wei_templates[sample]
        en2 = 0.
        for i in xrange(0,h_raw.GetNbinsX()+1):
            en2 += pow(h_raw.GetBinError(i) * h_wei.GetBinContent(i),2)
        return sqrt(en2)
   
    #____________________________________________________________
    def ndata_orig(self):
        return self.nhist(self.h_orig_data)

    #____________________________________________________________
    def endata_orig(self):
        return self.enhist(slef.h_orig_data)

    #____________________________________________________________
    def ndata_curr(self):
        return self.nhist(self.h_curr_data)

    #____________________________________________________________
    def endata_curr(self):
        return self.enhist(slef.h_curr_data)

    #____________________________________________________________
    def ntot_orig(self):
        return sum([self.nsamp_orig(s) for s in self.samples])

    #____________________________________________________________
    def entot_orig(self):
        en2 = sum([pow(self.ensamp_orig(s),2) for s in self.samples])
        return sqrt(en2)

    #____________________________________________________________
    def ntot_curr(self):
        return sum([self.nsamp_curr(s) for s in self.samples])

    #____________________________________________________________
    def entot_curr(self):
        en2 = sum([pow(self.ensamp_curr(s),2) for s in self.samples])
        return sqrt(en2)

    #____________________________________________________________
    def fsamp_true(self,sample):
        ntot = self.ntot_orig()
        nsamp = nsamp_orig(sample)
        return nsamp / ntot if ntot else 0.0
    



    #____________________________________________________________
    def fit_frac_and_error(self,ipar):
        return self.fit_results[ipar]
    
    #____________________________________________________________
    def fsamp_fit(self,sample,corr=None):
        if corr is None: corr = True
        if   sample == 'sig_1p': return self.fsig_1p_fit(corr)
        elif sample == 'sig_mp': return self.fsig_mp_fit(corr)
        ipar = self.fit_pars[sample]
        f, ef = self.fit_frac_and_error(ipar)
        if corr and self.toy_results:
            pullm = self.toy_results[sample]['pullm']
            f = f - ef * pullm
        return f
               
    #____________________________________________________________
    def efsamp_fit(self,sample,corr=None): 
        if corr is None: corr = True
        if   sample == 'sig_1p': return self.efsig_1p_fit(corr)
        elif sample == 'sig_mp': return self.efsig_mp_fit(corr)
        ipar = self.fit_pars[sample]
        ef = self.fit_frac_and_error(self.fit_pars[sample])[1]
        if corr and self.toy_results:
            ef *= self.toy_results[sample]['pulle']
        return ef
   
    #____________________________________________________________
    def fsig_1p_fit(self,corr=None):
        if corr is None: corr = True
        if not self.fix_r1p3p:
            f, ef = self.fit_frac_and_error(self.fit_pars['sig_1p'])
            if corr and self.toy_results:
                f = f - ef * self.toy_results['sig_1p']['pullm']
            return f
        else:
            fsig = self.fsamp_fit('sig',corr=corr)
            n1p = self.nsamp_orig('sig_1p') 
            ntot = self.nsamp_orig('sig')
            return fsig * n1p / ntot if ntot else 0.0
    #____________________________________________________________
    def efsig_1p_fit(self,corr=None):
        if corr is None: corr = True
        if not self.fix_r1p3p:
            ef = self.fit_frac_and_error(self.fit_pars['sig_1p'])[1]
            if corr and self.toy_results:
                ef *= self.toy_results['sig_1p']['pulle']
            return ef
        else:
            efsig = self.efsamp_fit('sig',corr=corr)
            n1p = self.nsamp_orig('sig_1p') 
            ntot = self.nsamp_orig('sig')
            return efsig * n1p / ntot if ntot else 0.0
    #____________________________________________________________
    def fsig_mp_fit(self,corr=None):
        if corr is None: corr = True
        if not self.fix_r1p3p:
            f, ef = self.fit_frac_and_error(self.fit_pars['sig_mp'])
            if corr and self.toy_results:
                f = f - ef  * self.toy_results['sig_mp']['pullm']
            return f
        else:
            fsig = self.fsamp_fit('sig')
            nmp = self.nsamp_orig('sig_mp') 
            ntot = self.nsamp_orig('sig')
            return fsig * nmp / ntot if ntot else 0.0
    #____________________________________________________________
    def efsig_mp_fit(self,corr=None):
        if corr is None: corr = True
        if not self.fix_r1p3p:
            ef = self.fit_frac_and_error(self.fit_pars['sig_mp'])[1]
            if corr and self.toy_results:
                ef *= self.toy_results['sig_mp']['pulle']
            return ef
        else:
            efsig = self.efsamp_fit('sig')
            nmp = self.nsamp_orig('sig_mp') 
            ntot = self.nsamp_orig('sig')
            return efsig * nmp / ntot if ntot else 0.0

   
    #____________________________________________________________
    def ntot_fit(self): 
        return self.nfit_tot
    #____________________________________________________________
    def nsamp_fit(self,sample,corr=None): 
        return self.fsamp_fit(sample,corr=corr) * self.ntot_fit() 
        #return self.fsamp_fit(sample,corr=corr) * self.ndata_curr()
    #____________________________________________________________
    def ensamp_fit(self,sample,corr=None): 
        return self.efsamp_fit(sample,corr=corr) * self.ntot_fit()
        #return self.efsamp_fit(sample,corr=corr) * self.ndata_curr()

    #____________________________________________________________
    def ensamp_frac_fit(self,sample,corr=None):
        f = self.fsamp_fit(sample,corr=corr)
        ef = self.efsamp_fit(sample,corr=corr)
        return ef / f if f else 0.0
    

#------------------------------------------------------------
class EffCalculator():
    '''
    description of EffCalculator
    '''
    #____________________________________________________________
    def __init__(self,
            name         = None,
            pass_models = [],
            fail_models = [],
            ):
        self.name = name
        self.pass_models = pass_models
        self.fail_models = fail_models
    
    #____________________________________________________________
    def nsamp_fit_pass(self,sample):
        return sum([m.nsamp_fit(sample) for m in self.pass_models])
    #____________________________________________________________
    def ensamp_fit_pass(self,sample):
        return sqrt(sum(pow(m.ensamp_fit(sample),2) for m in self.pass_models))
    #____________________________________________________________
    def nsamp_fit_fail(self,sample):
        return sum([m.nsamp_fit(sample) for m in self.fail_models])
    #____________________________________________________________
    def ensamp_fit_fail(self,sample):
        return sqrt(sum(pow(m.ensamp_fit(sample),2) for m in self.fail_models))
    #____________________________________________________________
    def eff_samp_fit(self,sample):
        Np = self.nsamp_fit_pass(sample)
        Nf = self.nsamp_fit_fail(sample)
        Ntot = Np + Nf
        eff = Np / Ntot if Ntot else 0.0
        return eff 
    #____________________________________________________________
    def eeff_samp_fit(self,sample):
        Np = self.nsamp_fit_pass(sample)
        Nf = self.nsamp_fit_fail(sample)
        eNp = self.ensamp_fit_pass(sample)
        eNf = self.ensamp_fit_fail(sample)
        return eff_err(Np,Nf,eNp,eNf)

    #____________________________________________________________
##     def eff_samp_mc(self,sample):
##         Np = self.nsamp_fit_pass(sample)
##         Nf = self.nsamp_fit_fail(sample)
##         Ntot = Np + Nf
##         eff = Np / Ntot if Ntot else 0.0
##         return eff 
##     #____________________________________________________________
##     def eeff_sample_mc(self,sample):
##         Np = self.nsamp_fit_pass(sample)
##         Nf = self.nsamp_fit_fail(sample)
##         eNp = self.ensamp_fit_pass(sample)
##         eNf = self.ensamp_fit_fail(sample)
##         return eff_err(Np,Nf,eNp,eNf)
    #____________________________________________________________
    def nsamp_mc_pass(self,sample):
        return sum([m.nsamp_orig(sample) for m in self.pass_models])
    #____________________________________________________________
    def ensamp_mc_pass(self,sample):
        return sqrt(sum(pow(m.ensamp_orig(sample),2) for m in self.pass_models))
    #____________________________________________________________
    def nsamp_mc_fail(self,sample):
        return sum([m.nsamp_orig(sample) for m in self.fail_models])
    #____________________________________________________________
    def ensamp_mc_fail(self,sample):
        return sqrt(sum(pow(m.ensamp_orig(sample),2) for m in self.pass_models))
    #____________________________________________________________
    def eff_samp_mc(self,sample):
        Np = self.nsamp_mc_pass(sample)
        Nf = self.nsamp_mc_fail(sample)
        Ntot = Np + Nf
        eff = Np / Ntot if Ntot else 0.0
        return eff 
    #____________________________________________________________
    def eeff_samp_mc(self,sample):
        Np = self.nsamp_mc_pass(sample)
        Nf = self.nsamp_mc_fail(sample)
        eNp = self.ensamp_mc_pass(sample)
        eNf = self.ensamp_mc_fail(sample)
        return eff_err(Np,Nf,eNp,eNf)
    #____________________________________________________________
    def eff_samp_mc(self,sample):
        Np = self.nsamp_mc_pass(sample)
        Nf = self.nsamp_mc_fail(sample)
        Ntot = Np + Nf
        eff = Np / Ntot if Ntot else 0.0
        return eff 
    #____________________________________________________________
    def eeff_samp_mc(self,sample):
        Np = self.nsamp_mc_pass(sample)
        Nf = self.nsamp_mc_fail(sample)
        eNp = self.ensamp_mc_pass(sample)
        eNf = self.ensamp_mc_fail(sample)
        return eff_err(Np,Nf,eNp,eNf)
    #____________________________________________________________
    def sf_samp(self,sample):
        eff_fit = self.eff_samp_fit(sample)
        eff_mc  = self.eff_samp_mc(sample)
        return eff_fit / eff_mc if eff_mc else 0.0

    #____________________________________________________________
    def esf_samp(self,sample):
        eff_fit = self.eff_samp_fit(sample)
        eff_mc  = self.eff_samp_mc(sample)
        eeff_fit = self.eeff_samp_fit(sample)
        eeff_mc  = self.eeff_samp_mc(sample)
        return sqrt(pow(eeff_fit/eff_fit,2)+pow(eeff_mc/eff_mc,2))
        
    

# - - - - - - - - - - function defs - - - - - - - - - - - - #
#____________________________________________________________
def integral_and_error(h,xmin=None,xmax=None):
    '''
    get integral and error in specified range
    '''
    min_bin = 1
    max_bin = h.GetNbinsX()
    if xmin: min_bin = h.GetXaxis().FindBin(xmin)
    if xmax: max_bin = h.GetXaxis().FindBin(xmax)
    error = ROOT.Double()
    total = h.IntegralAndError(min_bin,max_bin,error)
    return total,error

#____________________________________________________________
def integral(h,xmin=None,xmax=None):
    return integral_and_error(h,xmin,xmax)[0]

#____________________________________________________________
def integral_error(h,xmin=None,xmax=None):
    return integral_and_error(h,xmin,xmax)[1]



#____________________________________________________________
def generate_random_mc_hist(hin,hout):
    global grandom
    for i in xrange(hin.GetNbinsX()+1):
        n = hin.GetBinContent(i)
        en = hin.GetBinError(i)

        n = max(0.,grandom.Gaus(n,en))
        hout.SetBinContent(i,n)
        hout.SetBinError(i,en)


#____________________________________________________________
def fluctuate_hist_gaus(hin,hout):
    global grandom
    for i in xrange(hin.GetNbinsX()+1):
        n = hin.GetBinContent(i)
        en = sqrt(n) 
        n = max(0.001,grandom.Gaus(n,en))
        hout.SetBinContent(i,n)
        hout.SetBinError(i,en)


#____________________________________________________________
def fluctuate_hist_pois(hin,hout):
    global grandom
    #hout.Reset()
    for i in xrange(hin.GetNbinsX()+1):
        n = hin.GetBinContent(i)
        ##print 'n[%d] = %.2f'%(i,n)
        en = sqrt(n) if n>=0 else 0.0
        n = grandom.Poisson(n)
        hout.SetBinContent(i,n)
        hout.SetBinError(i,en)
    
    #ntot = grandom.Poisson(hin.Integral())
    #for i in xrange(ntot): hout.Fill(hin.GetRandom())

#____________________________________________________________
def fluctuate_hist_sumw2(hin,hout):
    global grandom
    for i in xrange(hin.GetNbinsX()+1):
        n = hin.GetBinContent(i)
        en = hin.GetBinError(i)
        n = max(0.,grandom.Gaus(n,en))
        hout.SetBinContent(i,n)
        hout.SetBinError(i,en)

#____________________________________________________________
def sum_hists(hins,hout):
    for i in xrange(hout.GetNbinsX()+1):
        n = 0.
        en2 = 0.
        for h in hins: 
            n+=h.GetBinContent(i)
            en2+=pow(h.GetBinError(i),2)
        hout.SetBinContent(i,n)
        #hout.SetBinError(i,sqrt(en2))
        hout.SetBinError(i,sqrt(n))




#____________________________________________________________
def generate_sum_hist(mc_hists,hout,random=True):
    global grandom
    # method 1 
    hout.SetEntries(0)
    for i in xrange(hout.GetNbinsX()+1):
        mctot = 0.
        for h in mc_hists: mctot+=h.GetBinContent(i)
        mctot = max(0.,mctot)
        if random: n = float(grandom.Poisson(mctot))
        else: n = float(int(mctot))
        en = sqrt(n)
        hout.SetBinContent(i,n)
        hout.SetBinError(i,en)

    # method 2
    '''
    clear_hist(hout)
    for h in mc_hists: 
        for i in xrange(grandom.Poisson(h.Integral())):
            hout.Fill(h.GetRandom())
    '''
    


#____________________________________________________________
def reset_hist(hin,hout):
    hout.Reset()
    for i in xrange(hin.GetNbinsX()+1):
        hout.SetBinContent(i,hin.GetBinContent(i))
        hout.SetBinError(i,hin.GetBinError(i))
    hout.SetEntries(hin.GetEntries())

#____________________________________________________________
def clear_hist(h):
    for i in xrange(h.GetNbinsX()+1):
        h.SetBinContent(i,0)
        h.SetBinError(i,0)

    h.SetEntries(0)




#____________________________________________________________
def eff_err(Np,Nf,eNp,eNf):

        # e = Np / (Np + Nf) = Np * (Np+Nf)^-1
        # 
        # de^2 = (de/dNp)^2(eNp)^2 + (de/dNf)^2(eNf)
        #
        # de/dNp = 1.(Np+Nf)^-1 + Np * -1(Np+Nf)^-2*1
        #   = (Np+Nf)^-1 - Np (Np+Nf)^-2
        #
        # de/dNf = -1 * Np * (Np+Nf)^-2
        #   = -Np(Np+Nf)^-2
        #
        #
        dedNp = pow(Np+Nf,-1) - Np*pow(Np+Nf,-2)
        dedNf = -Np*pow(Np+Nf,-2)
        return sqrt(pow(dedNp*eNp,2) + pow(dedNf*eNf,2)) 


#____________________________________________________________
def create_mean_hist(sample):
    h = ROOT.TH1F("h_%s_mean"%sample,";N %s Fit;Trials"%sample,100,999999,-999999)
    h.SetBit(ROOT.TH1.kCanRebin)
    return h

#____________________________________________________________
def create_error_hist(sample):
    h = ROOT.TH1F("h_%s_error"%sample,";Error %s Fit;Trials"%sample,100,999999,-999999)
    h.SetBit(ROOT.TH1.kCanRebin)
    return h

#____________________________________________________________
def create_pull_hist(sample):
    #h = ROOT.TH1F("h_%s_pull"%sample,";Pull %s Fit;Trials"%sample,100,999999,-999999)
    h = ROOT.TH1F("h_%s_pull"%sample,";Pull %s Fit;Trials"%sample,100,-5.,5.)
    #h.SetBit(ROOT.TH1.kCanRebin)
    return h

#____________________________________________________________
def create_diff_hist(sample):
    h = ROOT.TH1F("h_%s_diff"%sample,";N(Fit) - N(MC) %s;Trials"%sample,100,999999,-999999)
    h.SetBit(ROOT.TH1.kCanRebin)
    return h


#____________________________________________________________
def create_mc_hist(sample):
    h = ROOT.TH1F("h_%s_mc"%sample,";N %s MC;Trials"%sample,100,999999,-999999)
    h.SetBit(ROOT.TH1.kCanRebin)
    return h


#____________________________________________________________
def create_2d_mean_hist(sample1,sample2):
    h = ROOT.TH2F("h2_%s_%s_mean"%(sample1,sample2),";N %s Fit; N %s Fit;Trials"%(sample1,sample2),100,999999,-999999,100,999999,-999999)
    h.SetBit(ROOT.TH1.kCanRebin)
    return h


#____________________________________________________________
def split_hist_raw_weight_components(h):
    h_entries = h.Clone('%s_raw'%h.GetName())
    h_weights = h.Clone('%s_wei'%h.GetName())
    print 'splitting hist: ', h.GetName()
    for i in xrange(0,h.GetNbinsX()+1):
        n = h.GetBinContent(i)
        en = h.GetBinError(i)
        entries = round(pow(n/en,2)) if en else 0.0
        weight  = n / entries if entries else 1.0
        h_entries.SetBinContent(i,entries)
        h_entries.SetBinError(i,sqrt(entries))
        h_weights.SetBinContent(i,weight)
        h_weights.SetBinError(i,0.)
        print 'bin%d, orig: %.1f, raw ent: %.1f, wei: %.4f' % (i,n,entries,weight)

    return h_entries,h_weights





## EOF
