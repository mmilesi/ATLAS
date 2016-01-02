 #!/usr/bin/python

import array
import os
import sys

sys.path.append(os.path.abspath(os.path.curdir))

# -------------------------------
# Parser for command line options 
# -------------------------------
import argparse

parser = argparse.ArgumentParser(description='Plotting python macro for deriving real/fake lepton rates for MM.')

#***********************************
# positional arguments (compulsory!)
#***********************************
parser.add_argument('inputDir', metavar='inputDir',type=str,
		  help='path to the directory containing subdirs w/ input files')
#*******************
# optional arguments
#*******************
parser.add_argument('--myChannel', metavar='CHANNEL', dest='myChannel', default='', type=str,
		  help='Flavour composition of the two leptons in CR (*empty_string*, ElEl, MuEl, MuMu) - default is *empty_string*')
parser.add_argument('--usePrediction', metavar='DATA_TYPE', dest='usePrediction', action='store', default='DATA', type=str,
		    help='use Monte-Carlo (MC) or data (DATA) to derive rates - default is DATA')
parser.add_argument('--debug', dest='debug', action='store_true',
		    help='run in debug mode')
parser.add_argument('--doBkgSub', dest='doBkgSub', action='store_true',
		    help='subtract prompt MC in fake region (only to data)')		    
parser.add_argument('--saveOnlyRates', dest='saveOnlyRates', action='store_true',
		    help='save only rates')

args = parser.parse_args()

if not args.inputDir.endswith('/'):
   args.inputDir += '/'

# ---------------------------------------------------------------------------------
# NB: the following string definitions are to be chosen according to the 
#     output name of the files (produced by the MakePlots_HTopMultilep.py scripts)
#     used for comuting the T/!T(=L) ratio to derive rates.
#
#     An example name could be:
#  	  
#        MuElFakeCRMuL_ElProbePt
#
#        channels[h] + lep_types[i] + CR + leptons[j] + selections[k] + leptons[l] + variables[m]
# ---------------------------------------------------------------------------------

# -------------------------------------------------------
# for each channel, store which leptons can be associated
# -------------------------------------------------------
dict_channels_lep = {
 	             'ElEl': ['El'],
                     'MuEl': ['El','Mu'],
                     ''    : ['El','Mu'],
                     'MuMu': ['Mu'] 
		    }

list_lep         = dict_channels_lep[args.myChannel]
list_types       = ['Fake','Real']
list_variables   = ['ProbeEta','ProbePt']
list_selections  = ['T','L']
list_prediction  = ['expected', 'observed']   # expected --> use MC distribution for probe lepton to derive the rate ( used only as a cross check for REAL rate and in closure test) 
                                              # observed --> use DATA distribution for probe lepton to derive the rate
list_out_samples = ['factor','factorbkgsub','rate','ratebkgsub']

# ----------------------------------------------------------------
# if you want to subtract backgrounds to estimate the fake rates
#
# ( prompt subtraction, charge flip subtraction )
# ----------------------------------------------------------------

from ROOT import gROOT, TH1, TFile, TGraphAsymmErrors, Double

gROOT.SetBatch(True) 

TH1.SetDefaultSumw2()

def RatesToFactors( v ):
    v = float(v)
    if v < 0:
        v = 0.
    return v/(v+1)

hists  = {}
graphs = {}
yields = {}
fin    = []
		
channel_str = ''
if ( args.myChannel is not '' ):
   channel_str = '_' + args.myChannel
   if ( args.debug ) :
      print ' channel string: ', channel_str 
		     
for iLep in list_lep:
    
   print 'looking at lepton of flavour: ' , iLep
   
   for iVar in list_variables:
      
      print '\t looking at variable: ' , iVar
      
      for iType in list_types:
     
          print '\t\t looking at rate of type: ', iType
          
          for iSel in list_selections:
      

              print '\t\t\t object  is: ', iSel
              
              fname = args.inputDir + args.myChannel + iType + 'CR' + iLep + iSel + '/' + args.myChannel + iType + 'CR' + iLep + iSel + '_' + iLep + iVar + '.root'

	      if ( args.debug ) :
	          print '\t\t\t input filename: ', fname
	       
              fin.append( TFile(fname) )
              
	      for iPred in list_prediction:

                  print '\t\t\t\t checking prediction from: ', iPred

                  # L[-1] can be used to access the last item in a list
          	  htmp = fin[-1].Get( iPred )
          	        	
		  #
		  # let the macro know the name of the histogram in the input ROOT file
		  #   
		  # do not even bother for 'observed' hist if usePrediction = MC
		  #
		  if ( args.usePrediction is 'MC' ) and ( iPred is 'observed' ) :
		     continue
          	  
	  	  if htmp:

                     histname = iLep + '_' + iVar + '_' + iType + channel_str + '_' + iSel + '_' + iPred
                     
                     if (args.debug ) :
                        print '\t\t\t\t\t output histname (Num,Denom): ', histname  
          	      
                     if iVar is 'ProbeEta':			      
          	  	  
                        hists[histname]  = htmp.Clone( histname )
                        yields[histname] = hists[histname].Integral()
          	      
                     # ----------------
		     # rebinning for pT
		     # ----------------
		      
                     if iVar is 'ProbePt':		     
	  	  	  
                        if iType is 'Fake':
			   # Francesco's binning
                           #nBIN = 6
                           #xbins = [10,15,20,25,35,50,100]
                           
			   # TEST_13F
			   # 
                           # use the following for non prompt rate
                           #nBIN  = 9
                           #xbins = [10,19,28,37,46,55,64,82,100,265] # the last number in the array is the lower edge of the overflow bin
                           #
			   # use the following for ChFlip rate
                           #nBIN  = 12
                           #xbins = [10,19,28,37,46,55,64,82,100,136,172,208,265] # the last number in the array is the lower edge of the overflow bin

			   # TEST_13F - 2
			   # 
                           # use the following for non prompt rate
                           #nBIN  = 10
                           #xbins = [10,15,20,25,35,50,60,70,100,150,310] # the last number in the array is the lower edge of the overflow bin
			   # use the following for ChFlip rate
                           nBIN  = 12
                           xbins = [10,15,20,25,35,50,60,70,100,130,160,200,310]
                        else:
                           # Francesco's binning
                           #nBIN = 14
                           #xbins = [10,15,20,25,30,35,40,45,50,55,60,70,80,90,100]
			   
			   # TEST_13F
			   # 			   
                           nBIN  = 28
                           xbins = [10,13,16,19,22,25,28,31,34,37,40,49,58,67,76,85,94,103,112,121,130,139,148,157,166,175,193,229,265]

			   # TEST_13F - 2
			   #
                           nBIN  = 18
                           xbins = [10,15,20,25,30,35,40,45,50,55,60,70,80,90,100,130,160,200,310]

                        vxbins = array.array('d', xbins)
	  	  	
                        print '\t\t\t\t\t vxbins: ',vxbins
	  	  	  
                        # the rebinning method automatically creates a new histogram
                        #
                        hists[histname]  = htmp.Rebin( nBIN, histname, vxbins )
                        yields[histname] = hists[histname].Integral()
          	  
                  else:
          	      
                     print '\t\t\t\t\t error: histogram htmp does not exist'
	  	  
	      # -------------------------------------------------
              # subtract prompt bkg from data in fakes histograms	
	      # -------------------------------------------------
	  	    		
	      if args.doBkgSub :
	      
             	 if ( iType is 'Fake' ) and ( args.usePrediction is 'DATA' ):
             	     
		     hists[ iLep + '_' +iVar + '_' + iType + channel_str + '_' + iSel + '_' + list_prediction[1] ].Add( hists[ iLep + '_' + iVar + '_' + iType +'_' + args.myChannel + '_' + iSel + '_' + list_prediction[0] ], -1 )
             	 
		 # ------------------------------------------------------------
		 # set bin content to zero if subtraction gives negative yields
		 # 
		 # compute the yield for this histogram
		 # ------------------------------------------------------------
		 
		 for iBin in range( 0, hists[ iLep + '_'+ iVar +'_'+ iType + channel_str +'_'+  iSel +  list_prediction[1] ].GetNbinsX()+2):
             	     
		     if hists[ iLep+'_'+ iVar +'_'+ iType + channel_str +'_'+ iSel +  list_prediction[1] ].GetBinContent( iBin ) < 0:
                        hists[ iLep+'_'+ iVar +'_'+ iType + channel_str +'_'+ iSel +  list_prediction[1] ].SetBinContent( iBin, 0 )
                
                 yields[ iLep+'_'+ iVar +'_'+ iType + channel_str +'_'+ iSel +  list_prediction[1] ] = hists[ iLep+'_'+ iVar +'_'+ iType + channel_str +'_'+ iSel +  list_prediction[1] ].Integral()

	  # -----------------
	  # compute the ratio:
	  # T / !T
	  #
	  # and the efficiency:
	  # T / ( !T + T )
	  # -----------------
          
	  histname =  iLep + '_' + iVar +'_'+ iType + channel_str
	 
          print '\t\t histogram to be used in the ratio: ', histname

          # -----------------------------------------------------
	  # use also MC based estimate, but only for real leptons
	  # or if usePrediction is MC
          # -----------------------------------------------------
	   
          #if ( args.usePrediction is 'MC' ): # or iType  is 'Real':

          hists[histname + '_Rate_expected']  = hists[histname + '_T_expected'].Clone(histname + '_Rate_expected')
          hists[histname + '_Rate_expected'].Divide(hists[histname+'_L_expected'])
          yields[histname + '_Rate_expected'] = yields[histname + '_T_expected'] / yields[histname + '_L_expected']
	  
          #
          # For efficiency, make sure the errors are computed correctly
          # (NB: in this case, numerator and denominator are not independent sets of events!)
          #
          hist_pass = hists[histname + '_T_expected']
          hist_tot  = hists[histname + '_T_expected'] + hists[histname+'_L_expected']
          hist_eff  = hist_pass.Clone(histname + '_Efficiency_expected')
          hist_eff.Divide(hist_pass, hist_tot,1.0,1.0,"B")
          hists[histname + '_Efficiency_expected'] = hist_eff
          #
          # even better, use TGraphAsymmErrors
          # (handles the errors in case eff is 100%)
          #
          g_efficiency = TGraphAsymmErrors(hist_eff)
          g_efficiency.Divide(hist_pass,hist_tot,"cl=0.683 b(1,1) mode")
          graphs[histname + '_Efficiency_expected_graph'] = g_efficiency

          print '\t\t --> ratio hist name: ', histname + '_Rate_expected'
              
	  if ( args.usePrediction is 'DATA' ):
              hists[histname + '_Rate_observed'] = hists[histname + '_T_observed'].Clone( histname + '_Rate_observed' )
              hists[histname + '_Rate_observed'].Divide(hists[histname + '_L_observed'])
              yields[histname+'_Rate_observed']  = yields[histname + '_T_observed'] / yields[histname + '_L_observed']
              #
              # For efficiency, make sure the errors are computed correctly
              # (NB: in this case, numerator and denominator are not independent sets of events!)
              #
              hist_pass = hists[histname + '_T_observed']
              hist_tot  = hists[histname + '_T_observed'] + hists[histname+'_L_observed']
	      hist_eff  = hist_pass.Clone(histname + '_Efficiency_observed')
	      hist_eff.Divide(hist_pass, hist_tot,1.0,1.0,"B")
	      hists[histname + '_Efficiency_observed'] = hist_eff
              #
              # even better, use TGraphAsymmErrors
              # (handles the errors in case eff is 100%)
              #
              g_efficiency = TGraphAsymmErrors(hist_eff)
              g_efficiency.Divide(hist_pass,hist_tot,"cl=0.683 b(1,1) mode")
              graphs[histname + '_Efficiency_observed_graph'] = g_efficiency
	  	 
              print '\t\t --> ratio hist name: ', histname + '_Rate_observed'


outfile = open( args.inputDir + 'Rates.txt', 'w')                    
outfile.write( 'Rates for FF amd Matrix Method \n')

foutname = args.inputDir + 'Rates.root'
fout = TFile( foutname,'RECREATE' )
fout.cd()

for g in sorted( graphs.keys() ):
   

   string_to_replace = graphs[g].GetName() 
   
   print 'old graph name: ', string_to_replace
   
   string_to_replace = string_to_replace.replace("expected", "observed")
   
   print 'new graph name: ', string_to_replace
   
   graphs[g].SetName(string_to_replace)

   print "saving graph: ", graphs[g].GetName() 
   graphs[g].Write(string_to_replace)

for h in sorted( hists.keys() ):
   
   string_to_replace = hists[h].GetName() 
   
   print 'old hist name: ', string_to_replace
   
   string_to_replace = string_to_replace.replace("expected", "observed")
   
   print 'new hist name: ', string_to_replace
   
   hists[h].SetName(string_to_replace)
   
   if (not args.saveOnlyRates) or ('_Rate_' in h)  :
        
       print "saving histogram: ", hists[h].GetName() 

       hists[h].Write(string_to_replace)
       Rval=[]
       Rerr=[]
       Fval=[]
       Ferrup=[]
       Ferrdn=[]
        
       if '_Rate_' in h:
            
	  if (args.debug ) :  
             print "saving factors from rates for histogram: ", hists[h].GetName() 	    
	    
          Rtot=yields[h]
          Ftot=RatesToFactors(yields[h])
            
          for ibin in range( 1, hists[h].GetNbinsX()+1 ):
                
             Rval.append(hists[h].GetBinContent(ibin))
             #Rerr.append(hists[h].GetBinError(ibin)/hists[h].GetBinContent(ibin))
             Rerr.append(hists[h].GetBinError(ibin))
             Fval.append(RatesToFactors(hists[h].GetBinContent(ibin)))
             Ferrup.append(RatesToFactors(hists[h].GetBinContent(ibin)+hists[h].GetBinError(ibin))-RatesToFactors(hists[h].GetBinContent(ibin)))
             Ferrdn.append(RatesToFactors(hists[h].GetBinContent(ibin))-RatesToFactors(hists[h].GetBinContent(ibin)-hists[h].GetBinError(ibin)))
             #Ferrup.append(RatesToFactors(hists[h].GetBinContent(ibin)+hists[h].GetBinError(ibin))/RatesToFactors(hists[h].GetBinContent(ibin)))
             #Ferrdn.append(RatesToFactors(hists[h].GetBinContent(ibin)-hists[h].GetBinError(ibin))/RatesToFactors(hists[h].GetBinContent(ibin)))
            
	     if (args.debug ) :
             	print 'Rtot = ', round(Rtot, 3) 
             	print 'Rval = ', [ round(elem, 3) for elem in Rval ] 
             	print 'Rerr = ', [ round(elem, 3) for elem in Rerr ]
             	print 'Ftot = ', round(Ftot, 3) 
             	print 'Fval = ', [ round(elem, 3) for elem in Fval ]
             	print 'Ferrdn = ', [ round(elem, 3) for elem in Ferrdn ]
             	print 'Ferrup = ', [ round(elem, 3) for elem in Ferrup ]
             	print ''
	     
             outfile.write('%s \n' %(h) )
             outfile.write('Rtot = %s \n' %(round(Rtot, 3)))
             outfile.write('Rval = { %s }; \n' %(', '.join(str(e) for e in[ round(elem, 3) for elem in Rval ])) )
             outfile.write('Rerr = { %s }; \n' %(', '.join(str(e) for e in[ round(elem, 3) for elem in Rerr ])) )


outfile.close()
fout.Write()
fout.Close()

for f in range( 0,len(fin) ):
   fin[f].Close()
