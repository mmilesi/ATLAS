#THIS SCRIPT IS MADE TO COMPARE THE NJET DISTRIBUTIONS FROM MC12 TTW EVENTS AND FROM MC14 NJET EVENTS. IT ALSO CALCULATE THE REAL AND THE FAKE RATES IN TTBAR EVENTS. IT USES THE NOMINAL TTBARMC SELECTED REQUIRING TRUTH MATCHING FOR REAL RATES (OS REGION) AND THE FAKE EVENTS (AT LEAST ONE NON-PROMPT LEPTON) FOR FAKE RATES IN THE SS REGION
from math import sqrt, pow
from ROOT import TCanvas, TFile, TGraph, TGraphErrors, TColor, TAttFill, TStyle, TLegend, TH1, TH2, TH1D, gROOT, TF1, TTree, TH1I, gDirectory, TChain, TH2D
import array
import copy
import csv

gROOT.SetBatch(True)

TH1.SetDefaultSumw2()
TH2.SetDefaultSumw2()

#path_data12 = '/mel_data/ttH/NewReprocessed_Melb13_ttH_v1.1_noBjetcut/ttW/'
#path_dc14 = '/mel_data/ttH/NewReprocessed_Melb13_ttH_Run2Test18_DxAOD_nobjetcut/ttW/'
path_data12 = '/mel_data/ttH/NewReprocessed_Melb13_ttH_v1.1/ttW/'
path_dc14 = '/mel_data/ttH/NewReprocessed_Melb13_ttH_Run2Test22_DxAOD/ttW/'

import os

dists=['evtsel_jets_num','Jet0Pt']
channels=['ee','em','mm']
#this selection has no bjet requirements
selection={'ee':'evtsel_is_sltmatch && evtsel_dilep_type && evtsel_tau_num==0 && (isSS01==1 && Lep0Pt>20 && Lep1Pt>20) && isTT01==1 && evtsel_is_diel && ((fabs(Lep0PDG)!=11 || fabs(Lep0Eta)<1.5) && (fabs(Lep1PDG)!=11 || fabs(Lep1Eta)<1.5))', 'em':'evtsel_is_sltmatch && evtsel_dilep_type && evtsel_tau_num==0 && (isSS01==1 && Lep0Pt>20 && Lep1Pt>20) && isTT01==1 && evtsel_is_muel && ((fabs(Lep0PDG)!=11 || fabs(Lep0Eta)<1.5) && (fabs(Lep1PDG)!=11 || fabs(Lep1Eta)<1.5))', 'mm':'evtsel_is_sltmatch && evtsel_dilep_type && evtsel_tau_num==0 && (isSS01==1 && Lep0Pt>20 && Lep1Pt>20) && isTT01==1 && evtsel_is_dimu'}
#selection={'ee':'evtsel_is_sltmatch && evtsel_bjets_num>=1 && evtsel_dilep_type && evtsel_jets_num>=1 && evtsel_tau_num==0 && (isSS01==1 && Lep0Pt>20 && Lep1Pt>20) && isTT01==1 && evtsel_is_diel && ((fabs(Lep0PDG)!=11 || fabs(Lep0Eta)<1.5) && (fabs(Lep1PDG)!=11 || fabs(Lep1Eta)<1.5))', 'em':'evtsel_is_sltmatch && evtsel_bjets_num>=1 && evtsel_dilep_type && evtsel_jets_num>=1 && evtsel_tau_num==0 && (isSS01==1 && Lep0Pt>20 && Lep1Pt>20) && isTT01==1 && evtsel_is_muel && ((fabs(Lep0PDG)!=11 || fabs(Lep0Eta)<1.5) && (fabs(Lep1PDG)!=11 || fabs(Lep1Eta)<1.5))', 'mm':'evtsel_is_sltmatch && evtsel_bjets_num>=1 && evtsel_dilep_type && evtsel_jets_num>=1 && evtsel_tau_num==0 && (isSS01==1 && Lep0Pt>20 && Lep1Pt>20) && isTT01==1 && evtsel_is_dimu'}

group_list = os.listdir(path_data12)
group_list = group_list[:]
for group in group_list:
    f_data12 = TFile.Open(path_data12+group)
    f_dc14 = TFile.Open(path_dc14+group)
    t_data12 = f_data12.Get('physics')
    t_dc14 = f_dc14.Get('physics')
    #t_data12.SetWeight(1)
    #t_dc14.SetWeight(1)
        
    c = TCanvas("c","",50,50,600,600)
    for ch in range(0,len(channels)):
        for di in range(0,len(dists)):        
            #creating histos
            if dists[di] == 'evtsel_jets_num':
                nBIN = 10
                xbins = [ -0.5 , 0.5 , 1.5 , 2.5 , 3.5 , 4.5 , 5.5 , 6.5 , 7.5 , 8.5 , 9.5 ]
            if dists[di] == 'Jet0Pt':
                nBIN = 20
                BINlimit=200
            vxbins=array.array('d', xbins)

            histtitle='data12_'+channels[ch]+'_'+dists[di]
            if dists[di] == 'evtsel_jets_num':
                h_data12   = TH1D('h_data12',histtitle,nBIN,vxbins)                
            if dists[di] == 'Jet0Pt':
                h_data12   = TH1D('h_data12',histtitle,nBIN,0,BINlimit)

            histtitle='dc14_'+channels[ch]+'_'+dists[di]
            if dists[di] == 'evtsel_jets_num':
                h_dc14   = TH1D('h_dc14',histtitle,nBIN,vxbins)
            if dists[di] == 'Jet0Pt':
                h_dc14   = TH1D('h_dc14',histtitle,nBIN,0,BINlimit)

            t_data12.Project('h_data12', dists[di], selection[channels[ch]])
            t_dc14.Project('h_dc14', dists[di], selection[channels[ch]])
            h_data12.SetLineColor(1)
            h_dc14.SetLineColor(2)
            h_data12.SetLineWidth(2)
            h_dc14.SetLineWidth(2)

            h_data12.Draw()
            h_dc14.Draw('same')
            
            c.SaveAs('Comparison_test22/Comparison_'+group[:-5]+'_'+channels[ch]+'_'+dists[di]+'.eps')

#2D rates
path1 = '/mel_data/ttH/NewReprocessed_Melb13_ttH_Run2Test22_DxAOD_ttbarfaketest/tops/' #this is the fake rate path
path2 = '/mel_data/ttH/NewReprocessed_Melb13_ttH_Run2Test22_DxAOD/tops/' #this is the real rate path

dists=['evtsel_jets_num','Jet0Pt']
channels=['ee','em','mm']
selection = {
    'FakeMuL':'evtsel_is_sltmatch && evtsel_bjets_num>=1 && evtsel_dilep_type && evtsel_jets_num>=1  && evtsel_jets_num<=3 && (ElTagEta==-5. || fabs(ElTagEta)<1.5) && (isSS01==1 && isMuRateEvt) && MuProbeTight==0', 
    'FakeMuT':'evtsel_is_sltmatch && evtsel_bjets_num>=1 && evtsel_dilep_type && evtsel_jets_num>=1  && evtsel_jets_num<=3 && (ElTagEta==-5. || fabs(ElTagEta)<1.5) && (isSS01==1 && isMuRateEvt) && MuProbeTight==1', 
    'FakeElL':'evtsel_is_sltmatch && evtsel_bjets_num>=1 && evtsel_dilep_type && evtsel_jets_num>=1  && evtsel_jets_num<=3 && fabs(evtsel_ZSSee_candidate_mass-91)>20 && (isSS01==1 && isElRateEvt) && ElProbeTight==0', 
    'FakeElT':'evtsel_is_sltmatch && evtsel_bjets_num>=1 && evtsel_dilep_type && evtsel_jets_num>=1  && evtsel_jets_num<=3 && fabs(evtsel_ZSSee_candidate_mass-91)>20 && (isSS01==1 && isElRateEvt) && ElProbeTight==1', 
    'RealMuL':'evtsel_is_sltmatch && evtsel_bjets_num>=1 && evtsel_dilep_type && evtsel_jets_num>=1  && evtsel_jets_num<=3 && ((evtsel_JPsiee_candidate_mass>12 || evtsel_JPsiee_candidate_mass<0 ) && (evtsel_JPsimm_candidate_mass>12 || evtsel_JPsimm_candidate_mass<0)) && (isSS01==0 && isMuRateEvt) && MuProbeTight==0', 
    'RealMuT':'evtsel_is_sltmatch && evtsel_bjets_num>=1 && evtsel_dilep_type && evtsel_jets_num>=1  && evtsel_jets_num<=3 && ((evtsel_JPsiee_candidate_mass>12 || evtsel_JPsiee_candidate_mass<0 ) && (evtsel_JPsimm_candidate_mass>12 || evtsel_JPsimm_candidate_mass<0)) && (isSS01==0 && isMuRateEvt) && MuProbeTight==1', 
    'RealElL':'evtsel_is_sltmatch && evtsel_bjets_num>=1 && evtsel_dilep_type && evtsel_jets_num>=1  && evtsel_jets_num<=3 && ((evtsel_JPsiee_candidate_mass>12 || evtsel_JPsiee_candidate_mass<0 ) && (evtsel_JPsimm_candidate_mass>12 || evtsel_JPsimm_candidate_mass<0)) && (isSS01==0 && isElRateEvt) && ElProbeTight==0', 
    'RealElT':'evtsel_is_sltmatch && evtsel_bjets_num>=1 && evtsel_dilep_type && evtsel_jets_num>=1  && evtsel_jets_num<=3 && ((evtsel_JPsiee_candidate_mass>12 || evtsel_JPsiee_candidate_mass<0 ) && (evtsel_JPsimm_candidate_mass>12 || evtsel_JPsimm_candidate_mass<0)) && (isSS01==0 && isElRateEvt) && ElProbeTight==1',
}

group_list = os.listdir(path1)
group_list = group_list[:]
for group in group_list:
    f1 = TFile.Open(path1+group)
    t1 = f1.Get('physics')
    #t1.SetWeight(1)
    f2 = TFile.Open(path2+group)
    t2 = f2.Get('physics')
    #t2.SetWeight(1)

    outfile = open('Comparison_test22/Rates2D'+group+'.txt', 'w')                    
    outfile.write('Rates for Fake Factors \n')
    foutname='Comparison_test22/Rates2D'+group+'.root'
    fout=TFile(foutname,'RECREATE')
    fout.cd()

    for ch in ['FakeMu','FakeEl','RealMu','RealEl']:
        #eta bins
        nBINy = 8
        ybins = [ 0.0 , 0.5 , 0.8 , 1.1 , 1.37, 1.52 , 2.0 , 2.25 , 2.6]
        #pt bins
        if 'Fake' in ch:
            nBINx = 6
            xbins = [10,15,20,25,35,50,100]
        else:
            nBINx = 14
            xbins = [10,15,20,25,30,35,40,45,50,55,60,70,80,90,100]
        #nBINx = 6
        #xbins = [10,15,20,25,35,50,100]
        vxbins=array.array('d', xbins)
        vybins=array.array('d', ybins)

        if 'Fake' in ch:
            histname=ch+'T'
            h_T = TH2D(histname,'', nBINx, vxbins, nBINy, vybins)
            if 'El' in ch:
                t1.Project(histname, 'fabs(ElProbeEta):ElProbePt', selection[ch+'T'])
            else:
                t1.Project(histname, 'fabs(MuProbeEta):MuProbePt', selection[ch+'T'])
            histname=ch+'L'
            h_L = TH2D(histname,'', nBINx, vxbins, nBINy, vybins)
            if 'El' in ch:
                t1.Project(histname, 'fabs(ElProbeEta):ElProbePt', selection[ch+'L'])
            else:
                t1.Project(histname, 'fabs(MuProbeEta):MuProbePt', selection[ch+'L'])
        else:
            #real case
            histname=ch+'T'
            h_T = TH2D(histname,'', nBINx, vxbins, nBINy, vybins)
            if 'El' in ch:
                t2.Project(histname, 'fabs(ElProbeEta):ElProbePt', selection[ch+'T'])
            else:
                t2.Project(histname, 'fabs(MuProbeEta):MuProbePt', selection[ch+'T'])
            histname=ch+'L'
            h_L = TH2D(histname,'', nBINx, vxbins, nBINy, vybins)
            if 'El' in ch:
                t2.Project(histname, 'fabs(ElProbeEta):ElProbePt', selection[ch+'L'])
            else:
                t2.Project(histname, 'fabs(MuProbeEta):MuProbePt', selection[ch+'L'])


        h_T.SetLineColor(1)
        h_L.SetLineColor(2)
        h_T.SetLineWidth(2)
        h_L.SetLineWidth(2)

        values=[]        
        errors=[]
        histname=ch+'Ratio'
        if h_L.Integral() > 0:
            h_R = h_T.Clone(histname)
            h_R.Divide(h_L)
            histname=ch+'Error'
            h_E = h_R.Clone(histname) # contain the error on the ratio
            h_E.Reset()
            for xbin in range(1, h_R.GetNbinsX()+1):
                for ybin in range(1, h_R.GetNbinsY()+1):
                    h_E.SetBinContent(xbin,ybin,h_R.GetBinError(xbin,ybin))
                    values.append(h_R.GetBinContent(xbin,ybin))
                    errors.append(h_R.GetBinError(xbin,ybin))
            h_R.Write()
            h_E.Write()
            outfile.write('%s \n' %(ch) )
            outfile.write('FFval = { %s }; \n' %(', '.join(str(e) for e in[ round(elem, 3) for elem in values ])) )
            outfile.write('FFerr = { %s }; \n' %(', '.join(str(e) for e in[ round(elem, 3) for elem in errors ])) )

    outfile.close()
    fout.Write()
    fout.Close()

    #1D rates
    outfile = open('Comparison_test22/Rates'+group+'.txt', 'w')
    outfile.write('Rates for Fake Factors \n')
    foutname='Comparison_test22/Rates'+group+'.root'
    fout=TFile(foutname,'RECREATE')
    fout.cd()

    for ch in ['FakeMu','FakeEl','RealMu','RealEl']:
        #eta bins
        nBINy = 8
        ybins = [ 0.0 , 0.5 , 0.8 , 1.1 , 1.37, 1.52 , 2.0 , 2.25 , 2.6]
        #pt bins
        if 'Fake' in ch:
            nBINx = 6
            xbins = [10,15,20,25,35,50,100]
        else:
            nBINx = 14
            xbins = [10,15,20,25,30,35,40,45,50,55,60,70,80,90,100]
        #nBINx = 6
        #xbins = [10,15,20,25,35,50,100]
        vxbins=array.array('d', xbins)
        vybins=array.array('d', ybins)
        
        for var in ['ElProbeEta','ElProbePt','MuProbeEta','MuProbePt']:
            if ('El' in ch and 'Mu' in var) or ('Mu' in ch and 'El' in var):
                continue
            if 'Fake' in ch:
                if 'Pt' in var:
                    histname=ch+'T'+var
                    h_T = TH1D(histname,'', nBINx, vxbins)
                    t1.Project(histname, var, selection[ch+'T'])
                    histname=ch+'L'+var
                    h_L = TH1D(histname,'', nBINx, vxbins)
                    t1.Project(histname, var, selection[ch+'L'])
                elif 'Eta' in var:
                    histname=ch+'T'+var
                    h_T = TH1D(histname,'', nBINy, vybins)
                    t1.Project(histname, 'fabs('+var+')', selection[ch+'T'])
                    histname=ch+'L'+var
                    h_L = TH1D(histname,'', nBINy, vybins)
                    t1.Project(histname, 'fabs('+var+')', selection[ch+'L'])
            else:
                #real case
                if 'Pt' in var:
                    histname=ch+'T'+var
                    h_T = TH1D(histname,'', nBINx, vxbins)
                    t2.Project(histname, var, selection[ch+'T'])
                    histname=ch+'L'+var
                    h_L = TH1D(histname,'', nBINx, vxbins)
                    t2.Project(histname, var, selection[ch+'L'])
                elif 'Eta' in var:
                    histname=ch+'T'+var
                    h_T = TH1D(histname,'', nBINy, vybins)
                    t2.Project(histname, 'fabs('+var+')', selection[ch+'T'])
                    histname=ch+'L'+var
                    h_L = TH1D(histname,'', nBINy, vybins)
                    t2.Project(histname, 'fabs('+var+')', selection[ch+'L'])

            h_T.SetLineColor(1)
            h_L.SetLineColor(2)
            h_T.SetLineWidth(2)
            h_L.SetLineWidth(2)

            values=[]        
            errors=[]
            histname=ch+'Ratio'+var
            if h_L.Integral() > 0:
                h_R = h_T.Clone(histname)
                h_R.Divide(h_L)
                histname=ch+'Error'+var
                h_E = h_R.Clone(histname) # contain the error on the ratio
                h_E.Reset()
                for xbin in range(1, h_R.GetNbinsX()+1):
                    h_E.SetBinContent(xbin,h_R.GetBinError(xbin))
                    values.append(h_R.GetBinContent(xbin))
                    errors.append(h_R.GetBinError(xbin))
                h_R.Write()
                h_E.Write()
                FFtot=h_T.Integral()/h_L.Integral()
                outfile.write('%s %s \n' %(ch,var) )
                outfile.write('FFtot = %s \n' %(round(FFtot, 3)))
                outfile.write('FFval = { %s }; \n' %(', '.join(str(e) for e in[ round(elem, 3) for elem in values ])) )
                outfile.write('FFerr = { %s }; \n' %(', '.join(str(e) for e in[ round(elem, 3) for elem in errors ])) )

    outfile.close()
    fout.Write()
    fout.Close()
    f1.Close()
    f2.Close()
