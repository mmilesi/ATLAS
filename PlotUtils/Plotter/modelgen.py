# encoding: utf-8
'''
modelgen.py

description:

'''
__author__    = "Will Davey"
__email__     = "will.davey@cern.ch"
__created__   = "2013-02-28"
__copyright__ = "Copyright 2013 Will Davey"
__license__   = "GPL http://www.gnu.org/licenses/gpl.html"



## modules
import os
import ROOT
from utils import get_object
from math import sqrt

# - - - - - - - - - - - class defs  - - - - - - - - - - - - #
#------------------------------------------------------------
class TemplateGenerator():
    '''
    description of TemplateGeneator 
    '''
    #____________________________________________________________
    def __init__(self,
            base_dir         = None, 
            input_abcd       = None,
            input_wcontrol   = None,
            input_sig        = None,
            input_sig_sys    = None,
            bkg_method       = 'ss',
            mcwjets_shape    = False,
            binned_rosss     = True,
            config           = 'will',
            ):

        self.base_dir         = base_dir 
        self.input_abcd       = input_abcd
        self.input_wcontrol   = input_wcontrol
        self.input_sig        = input_sig
        self.input_sig_sys    = input_sig_sys
        self.bkg_method       = bkg_method
        self.mcwjets_shape    = mcwjets_shape 
        self.binned_rosss     = binned_rosss
        
        # set the config
        self.config = config
        if config == 'fed':
            # defaults
            if not base_dir: self.base_dir = \
                    '/home/wedavey/analysis/FedericoTauIDEff/InputFiles'
            if not input_abcd: self.input_abcd = \
                    "PlotsAll_ABCD_%s_p1344_%s_tau_total_kt_tracks_D04_dR06.root"
            if not input_wcontrol: self.input_wcontrol = \
                    "WScale_PlotsAll_ABCD_%s_p1344_%s_tau_total_kt_tracks_D04_dR06.root"
            if not input_sig: self.input_sig = \
                    "PlotsAll_SingSamp_%s_p1344_%s_tau_total_kt_tracks_D04_dR06_%s.root"
            # methods
            self.get_nominal_hists = self.get_nominal_hists_fed
            self.filename_getter_abcd = self.filename_getter_abcd_fed
            self.filename_getter_wcontrol = self.filename_getter_wcontrol_fed
            self.filename_getter_sig = self.filename_getter_sig_fed
            self.hist_getter_abcd = self.hist_getter_abcd_fed
            self.hist_getter_wcontrol = self.hist_getter_wcontrol_fed
            self.hist_getter_sig = self.hist_getter_sig_fed
            self.get_kin_bin = self.get_kin_bin_fed
        elif config == 'will':
            # defaults
            if not base_dir: self.base_dir = \
                    '/home/wedavey/analysis/FedericoTauIDEff/CoEPPPyTools2012/run/FitPlots03_Ntrack_All'
            if not input_abcd: self.input_abcd = \
                    "ABCD_Plots_%s_%s_%s%s.root"
            if not input_wcontrol: self.input_wcontrol = \
                    "WScale_%s_%s_%s%s.root"
            if not input_sig: self.input_sig = \
                    "Signal_%s_%s_%s_%s%s.root"
            if not input_sig_sys: self.input_sig_sys = \
                    "SignalSysReweight.root"

            self.get_nominal_hists = self.get_nominal_hists_fed
            self.filename_getter_abcd = self.filename_getter_abcd_will
            self.filename_getter_wcontrol = self.filename_getter_wcontrol_will
            self.filename_getter_sig = self.filename_getter_sig_will
            self.filename_getter_sig_sys = self.filename_getter_sig_sys_will
            self.hist_getter_abcd = self.hist_getter_abcd_fed
            self.hist_getter_wcontrol = self.hist_getter_wcontrol_fed
            self.hist_getter_sig = self.hist_getter_sig_will
            self.get_kin_bin = self.get_kin_bin_will

        elif config == 'usyd':
            self.get_nominal_hists = self.get_nominal_hists_fed
            self.hist_getter_abcd = self.hist_getter_abcd_usyd
            self.hist_getter_wcontrol = self.hist_getter_wcontrol_usyd
            self.hist_getter_sig = self.hist_getter_sig_usyd
            self.get_kin_bin = self.get_kin_bin_will

        else:
            print 'ERROR - invalid config: ', config
            exit(1)






        self.kin_bin = [
            "PtEtaAll",
            "pt2027",
            "pt2735",
            "pt35plus",
            "Bar",
            "End",
            ] 

        self.prongs = [
            "inc",
            "1P",
            "mP",
            ]

        self.abcd_sys = [
            "iso03",
            "iso05",
            ] 

        self.w_sys = [
            "scdp_n100",
            "scdp_n25",
            ] 
        self.w_sys += self.abcd_sys

        self.anti_sys = [
            'antiUP',
            'antiDN',
            ]

        self.sig_sys = [
            "Alpgen",
            "FTFP_BERT",
            "QGSP",
            "ExtraMaterial",
            "A2Tune",
            "QGSP"
            ]

        

        '''
        self.idmap = {
                'noid' : 'Base',
                'llhl' : 'LLHl',
                'llhm' : 'LLHm',
                'llht' : 'LLHt',
                'bdtl' : 'BDTl',
                'bdtm' : 'BDTm',
                'bdtt' : 'BDTt',
                }
        '''

        self.idmap = {
                'llhnl'  : 'LLHnl',
                'llhlnm' : 'LLHlnm',
                'llhmnt' : 'LLHmnt',
                'llht'   : 'LLHt',
                'bdtnl'  : 'BDTnl',
                'bdtlnm' : 'BDTlnm',
                'bdtmnt' : 'BDTmnt',
                'bdtt'   : 'BDTt',
                }


        ## Note:
        ## this map was designed to calculate exclusive 
        ## regions by subtracting tighter ones from 
        ## looser ones. this was done using the func
        ## hist_subtract_correlated, which calculates
        ## the uncertainty in each bin as:
        ##  en = sqrt(en1*en1 - en2*en2) 
        ## but there were problems since in some cases
        ## en2 > en1 (which should not be possible)
        ## so this method was abandoned, and instead
        ## the exclusive templates are now generated 
        ## explicitly in the plotting code.
        self.idexclmap_old = {
                'bdt!l'  : ['noid','bdtl'],
                'bdtl!m' : ['bdtl','bdtm'],
                'bdtm!t' : ['bdtm','bdtt'],
                'bdtt'   : ['bdtt',None],
                'llh!l'  : ['noid','llhl'],
                'llhl!m' : ['llhl','llhm'],
                'llhm!t' : ['llhm','llht'],
                'llht'   : ['llht',None],
                }



    #____________________________________________________________
    def id_type_will(self,ID=None):
        if not ID: return 'Base'
        if ID.count('bdtnl') or ID.count('llhnl'): return 'FailID'
        if ID.lower().count('bdt'): return 'BDT'
        if ID.lower().count('llh'): return 'LLH'
        return 'Base'

    #____________________________________________________________
    def get_kin_bin_fed(self,pt_bin=None,eta_bin=None):
        if pt_bin is None and eta_bin is None: 
            return 'PtEtaAll'
        assert not (pt_bin and eta_bin), "ERROR - cant set both pt and eta bin in Fed Setup"
        if pt_bin: return pt_bin
        return eta_bin

    #____________________________________________________________
    def get_kin_bin_will(self,pt_bin=None,eta_bin=None):
        if pt_bin is None: pt_bin = 'ptAll'
        if eta_bin is None: eta_bin = 'etaAll'
        return '%s_%s'%(pt_bin,eta_bin) 

    #____________________________________________________________
    def coepp_id_tag(self,ID=None):
        if ID is None: return 'Base'
        return self.idmap[ID]

    #____________________________________________________________
    def filename_getter_abcd_fed(self,sys=None,pt_bin=None,eta_bin=None,ID=None):
        if sys is None: sys = 'Nom'
        kin_bin = self.get_kin_bin(pt_bin,eta_bin)
        fname = self.input_abcd % (sys,kin_bin)
        fpath = os.path.join(self.base_dir,fname)
        return fpath

    #____________________________________________________________
    def filename_getter_sig_fed(self,sys=None,pt_bin=None,eta_bin=None,prongs=None,ID=None):
        if sys is None: sys = 'Nom'
        if prongs is None: prongs = 'inc'
        kin_bin = self.get_kin_bin(pt_bin,eta_bin)
        fname = self.input_sig % (sys,kin_bin,prongs)
        fpath = os.path.join(self.base_dir,fname)
        return fpath
       
    #____________________________________________________________
    def filename_getter_wcontrol_fed(self,sys=None,pt_bin=None,eta_bin=None,ID=None):
        if sys is None: sys = 'Nom'
        kin_bin = self.get_kin_bin(pt_bin,eta_bin)
        fname = self.input_wcontrol % (sys,kin_bin)
        fpath = os.path.join(self.base_dir,fname)
        return fpath

    #____________________________________________________________
    def filename_getter_abcd_will(self,pt_bin=None,eta_bin=None,ID=None,sys=None):
        ## Note: not currently using 'sys'
        if sys is None: sysstr = ''
        else:           sysstr = '_%s'% sys
        if pt_bin is None: pt_bin = 'ptAll'
        if eta_bin is None: eta_bin = 'etaAll'
        IDtype = self.id_type_will(ID)

        fname = self.input_abcd % (pt_bin,eta_bin,IDtype,sysstr)
        fpath = os.path.join(self.base_dir,fname)
        return fpath
    
    #____________________________________________________________
    def filename_getter_sig_will(self,sys=None,pt_bin=None,eta_bin=None,prongs=None,ID=None):
        if pt_bin is None: pt_bin = 'ptAll'
        if eta_bin is None: eta_bin = 'etaAll'
        if sys in ['pileupUP','pileupDN']: sysstr = '_%s'%sys
        else: sysstr = ''
        ## prongs must be either 1P or mP (inc is just sum of two)

        IDtype = self.id_type_will(ID)
        fname = self.input_sig % (prongs,pt_bin,eta_bin,IDtype,sysstr)
        fpath = os.path.join(self.base_dir,fname)
        return fpath

    #____________________________________________________________
    def filename_getter_sig_sys_will(self):
        return self.input_sig_sys
 
    #____________________________________________________________
    def filename_getter_wcontrol_will(self,pt_bin=None,eta_bin=None,sys=None,ID=None):
        if sys is None: sysstr = ''
        else:           sysstr = '_%s'%sys
        ## Note: not currently using 'sys'
        if pt_bin is None: pt_bin = 'ptAll'
        if eta_bin is None: eta_bin = 'etaAll'
        IDtype = self.id_type_will(ID)
         
        fname = self.input_wcontrol % (pt_bin,eta_bin,IDtype,sysstr)
        fpath = os.path.join(self.base_dir,fname)
        return fpath

    #____________________________________________________________
    def hist_getter_abcd_fed(self,sample=None,region=None,ID=None,sys=None,pt_bin=None,eta_bin=None,var=None,newname=None, antiID=None):
        if sample is None: sample = 'data'
        if region is None: region = 'regA'
        if var is None: var = 'evtsel_tau_nTrackTotal'
        IDtag = self.coepp_id_tag(ID)
        if sys not in self.abcd_sys: sys = None

        fname = self.filename_getter_abcd(sys=sys,pt_bin=pt_bin,eta_bin=eta_bin,ID=ID)
        hname = 'h_%s_%s_%s' % (var,region,sample)
        hpath = '%s/%s/hists/%s' % (IDtag,region,hname)
        h = get_object(hpath,fname)
        if newname: 
            h = h.Clone(newname)

        if antiID: 
            hanti = self.hist_getter_abcd_fed(sample=sample,region=region,ID=antiID,sys=sys,pt_bin=pt_bin,eta_bin=eta_bin,var=var)
            h = hist_subtract_correlated(h,hanti)
        return h

    #____________________________________________________________
    def hist_getter_sig_fed(self,sample=None,ID=None,sys=None,pt_bin=None,eta_bin=None,prongs=None,var=None,newname=None,antiID=None):
        if sample is None: sample = 'ZtautauTruth'
        if var is None: var = 'evtsel_tau_nTrackTotal'
        if prongs is None: prongs = 'inc'
        kin_bin = self.get_kin_bin(pt_bin=pt_bin,eta_bin=eta_bin) 
        IDtag = self.coepp_id_tag(ID)
        fname = self.filename_getter_sig(sys=sys,pt_bin=pt_bin,eta_bin=eta_bin,prongs=prongs,ID=ID)
        hname = 'h_%s_%s' % (var,sample)
        hpath = '%s_%s_%s/hists/%s' % (IDtag,prongs,kin_bin,hname)
        h = get_object(hpath,fname)
        if newname: 
            h = h.Clone(newname)
        if antiID: 
            hanti = self.hist_getter_sig_fed(sample=sample,ID=antiID,sys=sys,pt_bin=pt_bin,eta_bin=eta_bin,prongs=prongs,var=var)
            h = hist_subtract_correlated(h,hanti)
        return h


    #____________________________________________________________
    def hist_getter_sig_will(self,sample=None,ID=None,sys=None,pt_bin=None,eta_bin=None,prongs=None,var=None,newname=None,antiID=None):
        if sample is None: sample = 'ZtautauTruth'
        if var is None: var = 'evtsel_tau_nTrackTotal'
        if prongs is None: prongs = 'inc'
        
        if prongs == 'inc': 
            h1P = self.hist_getter_sig_will(sample=sample,ID=ID,sys=sys,pt_bin=pt_bin,eta_bin=eta_bin,prongs='1P',var=var)
            hmP = self.hist_getter_sig_will(sample=sample,ID=ID,sys=sys,pt_bin=pt_bin,eta_bin=eta_bin,prongs='mP',var=var)
            if newname: h = h1P.Clone(newname)
            else: h = h1P.Clone()
            h.Add(hmP)
            return h

        IDtag = self.coepp_id_tag(ID)
        fname = self.filename_getter_sig(sys=sys,pt_bin=pt_bin,eta_bin=eta_bin,prongs=prongs,ID=ID)
        hname = 'h_%s_%s' % (var,sample)
        hpath = '%s/hists/%s' % (IDtag,hname)
        h = get_object(hpath,fname)
        if newname: 
            h = h.Clone(newname)
        if antiID: 
            hanti = self.hist_getter_sig_will(sample=sample,ID=antiID,sys=sys,pt_bin=pt_bin,eta_bin=eta_bin,prongs=prongs,var=var)
            h = hist_subtract_correlated(h,hanti)

        
        ## signal systematics
        print 'in hist getter, sys: ', sys
        if sys in self.sig_sys:
            hbits = ['h']
            if   prongs == '1P': hbits += ['1p']
            elif prongs == 'mP': hbits += ['mp'] 
            hbits += ['total_tracks',
                     IDtag,
                     'Ztautau',
                     sys,
                     pt_bin,
                     eta_bin,
                     'ratio'
                     ]
            hname = '_'.join(hbits)
            fname = os.path.join(self.base_dir,self.filename_getter_sig_sys())
            hsys = get_object(hname,fname)
            for i in xrange(h.GetNbinsX()+1):
                c = hsys.GetBinContent(i)
                ec = hsys.GetBinError(i)
                if not 0.5 < c < 1.5: 
                    c = 1.
                    ec = 0.2
                n = h.GetBinContent(i)
                en = h.GetBinError(i)
                if n:
                    nnew = n*c
                    ennew = sqrt(pow(ec/c,2)+pow(en/n,2)) * nnew
                    #print 'bin%d, n: %.1f+-%.1f, c: %.3f+-%.3f, nnew: %.1f+-%.1f'%(i,n,en,c,ec,nnew,ennew)
                    h.SetBinContent(i,nnew)
                    h.SetBinError(i,ennew)
            #h.Multiply(hsys)

        return h


    #____________________________________________________________
    def hist_getter_sig_usyd(self,sample=None,ID=None,sys=None,pt_bin=None,eta_bin=None,prongs=None,var=None,newname=None,antiID=None):
        if sample is None: sample = 'ZtautauTruth'
        if var is None: var = 'evtsel_tau_nTrackTotal'
        if prongs is None: prongs = 'inc'
        kin_bin = self.get_kin_bin(pt_bin=pt_bin,eta_bin=eta_bin) 
        IDtag = self.coepp_id_tag(ID)
        fname = self.filename_getter_sig(sys=sys,pt_bin=pt_bin,eta_bin=eta_bin,prongs=prongs,ID=ID)
        hname = 'h_%s_%s' % (var,sample)
        hpath = '%s_%s_%s/hists/%s' % (IDtag,prongs,kin_bin,hname)
        h = get_object(hpath,fname)
        if newname: 
            h = h.Clone(newname)
        if antiID: 
            hanti = self.hist_getter_sig_fed(sample=sample,ID=antiID,sys=sys,pt_bin=pt_bin,eta_bin=eta_bin,prongs=prongs,var=var)
            h = hist_subtract_correlated(h,hanti)
        return h

    #____________________________________________________________
    def hist_getter_wcontrol_fed(self,sample=None,sign=None,ID=None,sys=None,pt_bin=None,eta_bin=None,var=None,newname=None,antiID=None):
        if sample is None: sample = 'data'
        if sign is None: sign = 'OS'
        if var is None: var = 'evtsel_tau_nTrackTotal'
        IDtag = self.coepp_id_tag(ID)
        
        w_sys = sys if sys in self.w_sys else None

        fname = self.filename_getter_wcontrol(sys=w_sys,pt_bin=pt_bin,eta_bin=eta_bin,ID=ID)
        hname = 'h_%s_%s' % (var,sample)
        hpath = '%s_%s/hists/%s' % (IDtag,sign,hname)
        h = get_object(hpath,fname)
        if newname: 
            h = h.Clone(newname)
        if antiID: 
            hanti = self.hist_getter_wcontrol_fed(sample=sample,sign=sign,ID=antiID,sys=sys,pt_bin=pt_bin,eta_bin=eta_bin,var=var)
            h = hist_subtract_correlated(h,hanti)
        return h




    #____________________________________________________________
    def get_nominal_hists_fed(self, 
            var = None, 
            pt_bin = None,
            eta_bin = None,
            idmap = None,
            sys = None,
            ):
        if idmap == None: idmap = self.idmap
        ## set hist getters
        hist_abcd = self.hist_getter_abcd
        hist_sig  = self.hist_getter_sig

        ## map for hists
        ## hists[sample][id]
        ## id: noid, bdtl, bdtm, bdtt, llhl, llhm, llht
        hists = {}
        hists['data']   = {}
        hists['sig']    = {}
        hists['sig_1p'] = {}
        hists['sig_mp'] = {}
        hists['anti']   = {}
        if self.bkg_method == 'wjets_qcd_osss':
            hists['wjets']  = {}
            hists['qcd']    = {}
        else:
            hists['bkg']    = {}
       
        
        for ID in idmap:

            ## data
            hists['data'][ID] = hist_abcd(sample='data',
                                          region='regA',
                                          ID = ID,
                                          pt_bin=pt_bin,
                                          eta_bin=eta_bin,
                                          var = var,
                                          newname = 'h_data_%s' % ID
                                          )
        
            ## signal
            # inc
            hists['sig'][ID] = hist_sig(ID=ID,
                                        var=var,
                                        pt_bin=pt_bin,
                                        eta_bin=eta_bin,
                                        sys=sys,
                                        newname = 'h_sig_%s' % ID
                                        )
            # 1p
            hists['sig_1p'][ID] = hist_sig(ID=ID,
                                           var=var,
                                           pt_bin=pt_bin,
                                           eta_bin=eta_bin,
                                           prongs='1P',
                                           sys=sys,
                                           newname = 'h_sig_1p_%s' % ID
                                           )
            # mp
            hists['sig_mp'][ID] = hist_sig(ID=ID,
                                           var=var,
                                           pt_bin=pt_bin,
                                           eta_bin=eta_bin,
                                           prongs='mP',
                                           sys=sys,
                                           newname = 'h_sig_mp_%s' % ID
                                           )
           
            ## bkg
            if self.bkg_method == 'ss':
                hists['bkg'][ID] = self.build_bkg_ss(ID=ID,
                                                     pt_bin=pt_bin,
                                                     eta_bin=eta_bin,
                                                     var=var,
                                                     newname = 'h_bkg_%s' % ID,
                                                     sys=sys,
                                                     )
            elif self.bkg_method == 'osss':
                h_bkg, control_hists = self.build_bkg_osss(ID=ID,
                                                     pt_bin=pt_bin,
                                                     eta_bin=eta_bin,
                                                     var=var,
                                                     newname = 'h_bkg_%s' % ID,
                                                     sys=sys,
                                                     )
                hists['bkg'][ID] = h_bkg
                if not 'control' in hists: hists['control'] = {}
                if not ID in hists['control']: hists['control'][ID] = []
                hists['control'][ID] += control_hists 


            elif self.bkg_method == 'wjets_qcd_osss':
                h_wjets,h_qcd,control_hists = self.build_wjets_qcd_osss(ID=ID,
                                                     pt_bin=pt_bin,
                                                     eta_bin=eta_bin,
                                                     var=var,
                                                     sys=sys,
                                                     )
                hists['wjets'][ID] = h_wjets.Clone('h_wjets_%s'%ID)
                hists['qcd'][ID]   = h_qcd.Clone('h_qcd_%s'%ID)
                if not 'control' in hists: hists['control'] = {}
                if not ID in hists['control']: hists['control'][ID] = []
                hists['control'][ID] += control_hists 


            ## anti
            h_anti = hist_abcd(sample='ZtautauAntiTruth',
                               region='regA',
                               ID = ID,
                               pt_bin = pt_bin,
                               eta_bin = eta_bin,
                               var = var,
                               newname = 'h_anti_%s' % ID
                               )
            h_anti.Add(hist_abcd(sample='Zmumu',
                                region='regA',
                                ID = ID,
                                pt_bin = pt_bin,
                                eta_bin = eta_bin,
                                var = var,
                                ))
            if sys == 'antiUP': 
                h_anti.Scale(1.5)
            elif sys == 'antiDN': 
                h_anti.Scale(0.50)

            hists['anti'][ID] = h_anti

            
            for sname in hists: 
                if sname in ['data','control']: continue
                fix_neg_bins(hists[sname][ID])


        return hists    




    #____________________________________________________________
    def build_bkg_ss(self,ID = None,sys=None,pt_bin=None,eta_bin=None,var = None, newname = None,antiID=None):
        hist_abcd = self.hist_getter_abcd
        return hist_abcd(sample='data',
                         region='regB',
                         ID = ID,
                         antiID = antiID,
                         sys = sys,
                         pt_bin = pt_bin,
                         eta_bin = eta_bin,
                         var = var,
                         newname = newname,
                         )
            

    #____________________________________________________________
    def build_bkg_osss(self,ID = None, sys=None,pt_bin=None,eta_bin=None, var = None, newname = None,antiID=None):

        h_est_Wjets_os, h_est_qcd_os, control_plots = self.build_wjets_qcd_osss(ID=ID,sys=sys,pt_bin=pt_bin,eta_bin=eta_bin,var=var,antiID=antiID)

        if newname: 
            h = h_est_qcd_os.Clone(newname)
        else:
            h = h_est_qcd_os

        ## full template estimation
        h.Add(h_est_Wjets_os)
        return h, control_plots


    #____________________________________________________________
    def build_wjets_osss(self,ID = None, sys=None,pt_bin=None,eta_bin=None, var = None, newname = None,antiID=None):

        hist_abcd = self.hist_getter_abcd
        hist_wcr  = self.hist_getter_wcontrol
        ## wcontrol region
        samples = ['data','Wmunu','WtaunuAlpgen','Zmumu','ZtautauTruth','ZtautauAntiTruth']
        hists_wcr_os = {}
        hists_wcr_ss = {}

        ## use default var in WCR if taking SR shape from MC
        ## (in this case you only need the overall kW factor 
        ## from the WCR, not the shape)
        wcr_var = None
        if not self.mcwjets_shape: wcr_var = var

        for s in samples: 
            hists_wcr_os[s] = hist_wcr(sample=s,
                                     sign='OS',
                                     ID=ID,
                                     antiID=antiID,
                                     sys=sys,
                                     pt_bin=pt_bin,
                                     eta_bin=eta_bin,
                                     var=wcr_var,
                                     )
            hists_wcr_ss[s] = hist_wcr(sample=s,
                                     sign='SS',
                                     ID=ID,
                                     antiID=antiID,
                                     sys=sys,
                                     pt_bin=pt_bin,
                                     eta_bin=eta_bin,
                                     var=wcr_var,
                                     )
        h_wcr_data_os = hists_wcr_os['data'].Clone('h_wcr_data_os_%s'%ID)
        h_wcr_data_os.Add(hists_wcr_os['Zmumu'], -1)
        h_wcr_data_os.Add(hists_wcr_os['ZtautauTruth'], -1)
        h_wcr_data_os.Add(hists_wcr_os['ZtautauAntiTruth'], -1)

        h_wcr_data_ss = hists_wcr_ss['data'].Clone('h_wcr_data_ss_%s'%ID)
        h_wcr_data_ss.Add(hists_wcr_ss['Zmumu'], -1)
        h_wcr_data_ss.Add(hists_wcr_ss['ZtautauTruth'], -1)
        h_wcr_data_ss.Add(hists_wcr_ss['ZtautauAntiTruth'], -1)

        h_wcr_Wjets_os = hists_wcr_os['Wmunu'].Clone('h_wcr_Wjets_%s'%ID)
        h_wcr_Wjets_os.Add(hists_wcr_os['WtaunuAlpgen'])
        h_wcr_Wjets_ss = hists_wcr_ss['Wmunu'].Clone('h_wcr_Wjets_%s'%ID)
        h_wcr_Wjets_ss.Add(hists_wcr_ss['WtaunuAlpgen'])

        ## signal region
        h_sr_Wjets_os = hist_abcd(sample='Wmunu',
                      region='regA',
                      ID = ID,
                      antiID = antiID,
                      sys = sys,
                      pt_bin = pt_bin,
                      eta_bin = eta_bin,
                      var = var,
                      newname = 'h_sr_Wjets_os_%s'%ID,
                      )
        h_sr_Wjets_os.Add(hist_abcd(sample='WtaunuAlpgen',
                      region='regA',
                      ID = ID,
                      antiID = antiID,
                      sys = sys,
                      pt_bin = pt_bin,
                      eta_bin = eta_bin,
                      var = var,
                      ))

        h_sr_Wjets_ss = hist_abcd(sample='Wmunu',
                      region='regB',
                      ID = ID,
                      antiID = antiID,
                      sys = sys,
                      pt_bin = pt_bin,
                      eta_bin = eta_bin,
                      var = var,
                      newname = 'h_sr_Wjets_ss_%s'%ID,
                      )
        h_sr_Wjets_ss.Add(hist_abcd(sample='WtaunuAlpgen',
                      region='regB',
                      ID = ID,
                      antiID = antiID,
                      sys = sys,
                      pt_bin = pt_bin,
                      eta_bin = eta_bin,
                      var = var,
                      ))


        ## Estimate Wjets contribution from data
        ##  Note: this is conceptually the nicest 
        ##        way it could be done. But unfortunately
        ##        in these ntuples we have the kW 
        ##        factor applied in the signal region 
        ##        but not in the WCR.
        ##        So instead we use an equivalent method
        ##        that is less intuitive
        ## Ideal 1
        """
        # bin-by-bin MC WCR->SR transfer factors
        h_trans_os = h_sr_Wjets_os.Clone('h_trans_os')
        h_trans_os.Divide(h_wcr_Wjets_os)
        h_trans_ss = h_sr_Wjets_ss.Clone('h_trans_ss')
        h_trans_ss.Divide(h_wcr_Wjets_ss)
        # estimated W contribution
        h_est_Wjets_os = h_wcr_data_os.Clone('h_est_Wjets_os')
        h_est_Wjets_os.Multiply(h_trans_os)
        h_est_Wjets_ss = h_wcr_data_ss.Clone('h_est_Wjets_ss')
        h_est_Wjets_ss.Multiply(h_trans_ss)
        """
        ## Attempt 1
        """
        h_wcr_Wjets_os.Scale(1./h_wcr_Wjets_os.Integral())
        h_wcr_Wjets_ss.Scale(1./h_wcr_Wjets_ss.Integral())
        h_wcr_data_os.Scale(1./h_wcr_data_os.Integral())
        h_wcr_data_ss.Scale(1./h_wcr_data_ss.Integral())
        h_trans_os = h_wcr_data_os.Clone('h_trans_os_%s'%ID) 
        h_trans_os.Divide(h_wcr_Wjets_os)
        h_trans_ss = h_wcr_data_ss.Clone('h_trans_ss_%s'%ID) 
        h_trans_ss.Divide(h_wcr_Wjets_ss)
        h_est_Wjets_os = h_sr_Wjets_os.Clone('h_est_Wjets_os_%s'%ID)
        h_est_Wjets_os.Multiply(h_trans_os)
        h_est_Wjets_ss = h_sr_Wjets_ss.Clone('h_est_Wjets_ss_%s'%ID)
        h_est_Wjets_ss.Multiply(h_trans_ss)
        """
        ## Attempt 3
        '''
        take shape directly from data in WCR. 
        scale by single valued constant back to SR, 
        ie. dont do bin-by-bin reweight. 
        this should be ok, since we do a 5 GeV pt 
        binning now anyway.
        '''
        """
        n_wcr_Wjets_os = h_wcr_Wjets_os.Integral()
        n_wcr_data_os = h_wcr_data_os.Integral()
        kW_os = n_wcr_data_os / n_wcr_Wjets_os
        n_sig_Wjets_os = h_sr_Wjets_os.Integral() / kW_os
        transfer_wcr_os = n_sig_Wjets_os / n_wcr_Wjets_os
        h_est_Wjets_os = h_wcr_data_os.Clone('h_est_Wjets_os_%s'%ID)
        h_est_Wjets_os.Scale(transfer_wcr_os)

        n_wcr_Wjets_ss = h_wcr_Wjets_ss.Integral()
        n_wcr_data_ss = h_wcr_data_ss.Integral()
        kW_ss = n_wcr_data_ss / n_wcr_Wjets_ss
        n_sig_Wjets_ss = h_sr_Wjets_ss.Integral() / kW_ss
        transfer_wcr_ss = n_sig_Wjets_ss / n_wcr_Wjets_ss
        h_est_Wjets_ss = h_wcr_data_ss.Clone('h_est_Wjets_ss_%s'%ID)
        h_est_Wjets_ss.Scale(transfer_wcr_ss)
        """
        ## Attempt 4
        '''
        now we have removed kW from signal region
        so we can do a direct MC->MC transfer factor
        '''

        ## MC shape
        if self.mcwjets_shape:
            """
            W(i) = W_MC_SR(i) * kW_WCR
            ie. W in the signal region is estimated
            by taking the MC W contribution in the SR, 
            and scaling it by the kW factor measured in 
            the W control region.
            """
            n_wcr_Wjets_os = h_wcr_Wjets_os.Integral()
            n_wcr_data_os = h_wcr_data_os.Integral()
            kW_os = n_wcr_data_os / n_wcr_Wjets_os if n_wcr_Wjets_os else 0.0
            h_est_Wjets_os = h_sr_Wjets_os.Clone('h_est_Wjets_os_%s'%ID)
            h_est_Wjets_os.Scale(kW_os)
            #print 'ID: %s, kW_os: %s' % (ID,kW_os)

            n_wcr_Wjets_ss = h_wcr_Wjets_ss.Integral()
            n_wcr_data_ss = h_wcr_data_ss.Integral()
            kW_ss = n_wcr_data_ss / n_wcr_Wjets_ss if n_wcr_Wjets_ss else 0.0
            h_est_Wjets_ss = h_sr_Wjets_ss.Clone('h_est_Wjets_ss_%s'%ID)
            h_est_Wjets_ss.Scale(kW_ss)

        
        ## WCR data shape
        else:
            """
            W(i) = W_DATA_WCR(i) * W_MC_SR / W_MC_CR
            ie. W in the signal region is estimated
            by taking the data in the WCR and scaling
            it by the WCR->SR transfer factor measured
            from MC.

            Note: contamination from processes other than 
            Wmunu and Wtaunu in the WCR is subtracted from 
            the data using MC.
            """
            n_wcr_Wjets_os = h_wcr_Wjets_os.Integral()
            n_sig_Wjets_os = h_sr_Wjets_os.Integral()
            transfer_wcr_os = n_sig_Wjets_os / n_wcr_Wjets_os
            h_est_Wjets_os = h_wcr_data_os.Clone('h_est_Wjets_os_%s'%ID)
            h_est_Wjets_os.Scale(transfer_wcr_os)

            n_wcr_Wjets_ss = h_wcr_Wjets_ss.Integral()
            n_sig_Wjets_ss = h_sr_Wjets_ss.Integral()
            transfer_wcr_ss = n_sig_Wjets_ss / n_wcr_Wjets_ss
            h_est_Wjets_ss = h_wcr_data_ss.Clone('h_est_Wjets_ss_%s'%ID)
            h_est_Wjets_ss.Scale(transfer_wcr_ss)

        return h_est_Wjets_os, h_est_Wjets_ss

    #____________________________________________________________
    def build_wjets_qcd_osss(self,ID = None, sys=None,pt_bin=None,eta_bin=None, var = None, antiID=None):

        h_est_Wjets_os, h_est_Wjets_ss = self.build_wjets_osss(ID=ID,sys=sys,pt_bin=pt_bin,eta_bin=eta_bin,var=var,antiID=antiID)
        hist_abcd = self.hist_getter_abcd
        control_plots = [h_est_Wjets_ss]

        ## ss data
        h_est_qcd_ss = hist_abcd(sample='data',
                      region='regB',
                      ID = ID,
                      antiID=antiID,
                      sys = sys,
                      pt_bin = pt_bin,
                      eta_bin = eta_bin,
                      var = var,
                      newname = 'h_est_qcd_ss_%s'%(ID),
                      )

        ## anti-isolated region 
        anti_var = None
        if self.binned_rosss: anti_var = var
        h_data_regC = hist_abcd(sample='data',
                      region='regC',
                      ID = ID,
                      antiID = antiID,
                      sys = sys,
                      pt_bin = pt_bin,
                      eta_bin = eta_bin,
                      var = anti_var,
                      )
        h_data_regD = hist_abcd(sample='data',
                      region='regD',
                      ID = ID,
                      antiID = antiID,
                      sys = sys,
                      pt_bin = pt_bin,
                      eta_bin = eta_bin,
                      var = anti_var,
                      )

        ## full template estimation
        h_est_qcd_ss.Add( h_est_Wjets_ss, -1.)
        control_plots.append(h_est_qcd_ss)
        h_est_qcd_os = h_est_qcd_ss.Clone('h_est_qcd_os_%s'%ID)

        ## bin-by-bin Ros/ss factors
        if self.binned_rosss:
            h_rosss = h_data_regC.Clone('h_rosss_%s'%ID)
            h_rosss.Divide(h_data_regD)
            h_est_qcd_os.Multiply(h_rosss)
            control_plots.append(h_rosss)
        else: 
            nos = h_data_regC.Integral()
            nss = h_data_regD.Integral()
            rosss = nos / nss if nss else 0.0
            h_est_qcd_os.Scale(rosss)

        return h_est_Wjets_os, h_est_qcd_os, control_plots





# - - - - - - - - - - function defs - - - - - - - - - - - - #
#____________________________________________________________
def fix_neg_bins(h):
    for i in range(1,h.GetNbinsX()+1):
        n = h.GetBinContent(i)
        en = h.GetBinError(i)
        if n<0: 
            #h.SetBinContent(i,abs(en/2.))
            h.SetBinContent(i,0.)


#____________________________________________________________
def hist_subtract_correlated(h1,h2):
    print 'hist subtract corr. %s - %s' % (h1.GetName(),h2.GetName())
    for i in xrange(0,h1.GetNbinsX()+1):
        n1 = h1.GetBinContent(i)
        en1 = h1.GetBinError(i)
        n2 = h2.GetBinContent(i)
        en2 = h2.GetBinError(i)
        n = n1-n2
        print 'bin%d %10.1f+-%5.1f - %10.1f+-%5.1f = %10d' % (i,n1,en1,n2,en2,n)
        en = sqrt(en1*en1-en2*en2)
        h1.SetBinContent(i,n)
        h1.SetBinError(i,en)


    return h1






## EOF
