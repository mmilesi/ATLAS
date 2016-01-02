#THIS PROGRAM TAKE DISTRIBUTIONS FOR LOOSE AND TIGHT SELECTION AND DIVIDES THEM IN ORDER TO GET EFFICIENCIES OR CONVERSION FACTORS TO BE USED IN THE MATRIX OR THE FF METHOD
from ROOT import TH1, TFile, Double
import array, os, sys
TH1.SetDefaultSumw2()

input_path='/Users/fnuti/workarea/ATLASAnalyses/ttH/plots_from_common_ntuple_melbourne/ttH/v1_1_results/'
channels=['ElEl', 'MuEl', '', 'MuEl', 'MuMu', '']#the '' are the channels where ee and emu are together for el rates and emu,mumu for muon rates
leptons=['El', 'El', 'El', 'Mu', 'Mu', 'Mu']
lep_types=['Fake','Real']
dists=['ProbeEta', 'ProbePt']
samples=['expected', 'observed']
selections=['T','L']
out_samples=['factor', 'factorbkgsub' , 'rate', 'ratebkgsub']

#if you want to subtract bkgrounds to estimate the fake rates
doBkgSub=True
saveOnlyRates=True

def FFtoMM(v):
    v=float(v)
    if v<0:
        v=0.
    return v/(v+1)

hists={}
yields={}
fin=[]
for ch in range(0,len(channels)):
    for di in range(0,len(dists)):
        for ty in range(0,len(lep_types)):
            for se in range(0,len(selections)):
                fname=input_path+channels[ch]+lep_types[ty]+'CR'+leptons[ch]+selections[se]+'/'+channels[ch]+lep_types[ty]+'CR'+leptons[ch]+selections[se]+'_'+dists[di]+'.root'
                fin.append(TFile(fname))
                #cloning all the files
                for sam in range(0,len(samples)):
                    histname=leptons[ch]+'_'+dists[di]+'_'+lep_types[ty]+'_'+channels[ch]+'_'+selections[se]+'_'+samples[sam]
                    #print 'histname ', histname
                    htmp=fin[-1].Get(samples[sam])
                    if htmp:
                        if dists[di] == 'ProbeEta':                        
                            hists[histname]=htmp.Clone(histname)
                            yields[histname]=hists[histname].Integral()
                        if dists[di] == 'ProbePt':                        
                            #rebinning
                            if lep_types[ty] == 'Fake':
                                nBIN = 6
                                xbins = [10,15,20,25,35,50,100]
                            else:
                                nBIN = 14
                                xbins = [10,15,20,25,30,35,40,45,50,55,60,70,80,90,100]
                            vxbins=array.array('d', xbins)
                            #print 'vxbins ',vxbins
                            #the rebinning method automatically create a new histogram
                            hists[histname]=htmp.Rebin(nBIN, histname, vxbins)
                            yields[histname]=hists[histname].Integral()
                    else:
                        print 'error files ', fin[-1], ' or histogram htmp ', htmp, ' do not exist'
                    #print 'hists = ', hists
                if doBkgSub:
                    #Subtracting prompt bkg from data in fakes histos                    
                    if lep_types[ty] == 'Fake':
                        hists[leptons[ch]+'_'+dists[di]+'_'+lep_types[ty]+'_'+channels[ch]+'_'+selections[se]+'_observed'].Add(hists[leptons[ch]+'_'+dists[di]+'_'+lep_types[ty]+'_'+channels[ch]+'_'+selections[se]+'_expected'],-1)
                    for ibin in range(0, hists[leptons[ch]+'_'+dists[di]+'_'+lep_types[ty]+'_'+channels[ch]+'_'+selections[se]+'_observed'].GetNbinsX()+2):
                        if hists[leptons[ch]+'_'+dists[di]+'_'+lep_types[ty]+'_'+channels[ch]+'_'+selections[se]+'_observed'].GetBinContent(ibin) < 0:
                            hists[leptons[ch]+'_'+dists[di]+'_'+lep_types[ty]+'_'+channels[ch]+'_'+selections[se]+'_observed'].SetBinContent(ibin,0)
                            yields[histname]=hists[leptons[ch]+'_'+dists[di]+'_'+lep_types[ty]+'_'+channels[ch]+'_'+selections[se]+'_observed'].Integral()
            #NOW CALCULATING THE RATIO
            histname=leptons[ch]+'_'+dists[di]+'_'+lep_types[ty]+'_'+channels[ch]
            hists[histname+'_R_observed']=hists[histname+'_T_observed'].Clone(histname+'_R_observed')
            hists[histname+'_R_observed'].Divide(hists[histname+'_L_observed'])
            yields[histname+'_R_observed']=yields[histname+'_T_observed']/yields[histname+'_L_observed']
            #also calculating rates using mc only for real leptons
            if lep_types[ty] =='Real':
                hists[histname+'_R_expected']=hists[histname+'_T_expected'].Clone(histname+'_R_expected')
                hists[histname+'_R_expected'].Divide(hists[histname+'_L_expected'])
                yields[histname+'_R_expected']=yields[histname+'_T_expected']/yields[histname+'_L_expected']
                
outfile = open(input_path+'Rates.txt', 'w')                    
outfile.write('Rates for FF amd Matrix Method \n')
foutname=input_path+'Rates.root'
fout=TFile(foutname,'RECREATE')
fout.cd()
for h in sorted(hists.keys()):
    if not saveOnlyRates or '_R_' in h:
        hists[h].Write()
        FFval=[]
        FFerr=[]
        MMval=[]
        MMerrup=[]
        MMerrdn=[]
        if '_R_' in h:
            FFtot=yields[h]
            MMtot=FFtoMM(yields[h])
            for ibin in range(1, hists[h].GetNbinsX()+1):
                FFval.append(hists[h].GetBinContent(ibin))
                #FFerr.append(hists[h].GetBinError(ibin)/hists[h].GetBinContent(ibin))
                FFerr.append(hists[h].GetBinError(ibin))
                MMval.append(FFtoMM(hists[h].GetBinContent(ibin)))
                MMerrup.append(FFtoMM(hists[h].GetBinContent(ibin)+hists[h].GetBinError(ibin))-FFtoMM(hists[h].GetBinContent(ibin)))
                MMerrdn.append(FFtoMM(hists[h].GetBinContent(ibin))-FFtoMM(hists[h].GetBinContent(ibin)-hists[h].GetBinError(ibin)))
                #MMerrup.append(FFtoMM(hists[h].GetBinContent(ibin)+hists[h].GetBinError(ibin))/FFtoMM(hists[h].GetBinContent(ibin)))
                #MMerrdn.append(FFtoMM(hists[h].GetBinContent(ibin)-hists[h].GetBinError(ibin))/FFtoMM(hists[h].GetBinContent(ibin)))
            print 'Fake factors and rate for ', h
            print 'FFtot = ', round(FFtot, 3) 
            print 'FFval = ', [ round(elem, 3) for elem in FFval ] 
            print 'FFerr = ', [ round(elem, 3) for elem in FFerr ]
            print 'MMtot = ', round(MMtot, 3) 
            print 'MMval = ', [ round(elem, 3) for elem in MMval ]
            print 'MMerrdn = ', [ round(elem, 3) for elem in MMerrdn ]
            print 'MMerrup = ', [ round(elem, 3) for elem in MMerrup ]
            print ''
            outfile.write('%s \n' %(h) )
            outfile.write('FFtot = %s \n' %(round(FFtot, 3)))
            outfile.write('FFval = { %s }; \n' %(', '.join(str(e) for e in[ round(elem, 3) for elem in FFval ])) )
            outfile.write('FFerr = { %s }; \n' %(', '.join(str(e) for e in[ round(elem, 3) for elem in FFerr ])) )


outfile.close()
fout.Write()
fout.Close()
for f in range(0,len(fin)):
    fin[f].Close()
