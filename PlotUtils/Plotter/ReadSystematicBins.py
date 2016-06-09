 #!/usr/bin/python

import os
import sys

from ROOT import gROOT, TH1D, TFile, Double

gROOT.SetBatch(True)

def get_yields(nominal, up=None, down=None):

  #print("\t\tGetNbinsX() = {0}".format(nominal.GetNbinsX()))
  #print("\t\trange(0,GetNbinsX()+1) = {0}".format(range(0,nominal.GetNbinsX()+1)))

  # pick also O-Flow bin
  for bin in range(0,nominal.GetNbinsX()+2):

     nextbin = bin
     #if nominal.IsBinOverflow(bin+1):
     #  nextbin = bin + 1

     stat_error    = Double(0)
     value_nominal = nominal.IntegralAndError(bin,nextbin,stat_error)

     if ( up and down ):
       dummy_err_up  = Double(0)
       value_up = up.IntegralAndError(bin,nextbin,dummy_err_up)
       dummy_err_down  = Double(0)
       value_down = down.IntegralAndError(bin,nextbin,dummy_err_down)

       delta_up   = value_up - value_nominal
       delta_down = value_down - value_nominal

       sys_up   = 0.0
       sys_down = 0.0

       if delta_up >= 0.0:
         sys_up = abs( delta_up )
       if delta_down <= 0.0:
         sys_down = abs( delta_down )

       # Take the sys error as max( abs(value(up,down) - value_nominal) )
       #
       #list_variations = []
       #list_variations.append(abs( value_up - value_nominal ))
       #list_variations.append(abs( value_down - value_nominal ))
       #
       #sys_error  = max(list_variations)

       if nominal.IsBinOverflow(bin):
         #print ("\t\t{0}-th bin (O-FLOW): integral = {1} +- {2} (stat) +- {3} (syst)".format( bin, value_nominal, stat_error, sys_error ))
         print ("\t\t{0}-th bin (O-FLOW): integral = {1} +- {2} (stat) (+ {3}, - {4}) (syst)".format( bin, value_nominal, stat_error, sys_up, sys_down ))
       else:
         #print ("\t\t{0}-th bin: integral = {1} +- {2} (stat) +- {3} (syst)".format( bin, value_nominal, stat_error, sys_error ))
         print ("\t\t{0}-th bin: integral = {1} +- {2} (stat) (+ {3}, - {4}) (syst)".format( bin, value_nominal, stat_error, sys_up, sys_down ))

     else:
       if nominal.IsBinOverflow(bin):
         print ("\t\t{0}-th bin (O-FLOW): integral = {1} +- {2} (stat)".format( bin, value_nominal, stat_error ))
       else:
         print ("\t\t{0}-th bin: integral = {1} +- {2} (stat)".format( bin, value_nominal, stat_error ))

  integral_stat_error = Double(0)
  integral_nominal    = nominal.IntegralAndError(0,nominal.GetNbinsX()+1,integral_stat_error)

  print ("\t\t--------------------")
  if ( up and down ):

    delta_up   = up.Integral(0,up.GetNbinsX()+1) - integral_nominal
    delta_down = down.Integral(0,down.GetNbinsX()+1) - integral_nominal

    integral_sys_up   = 0.0
    integral_sys_down = 0.0

    if delta_up >= 0.0:
      integral_sys_up = abs( delta_up )
    if delta_down <= 0.0:
      integral_sys_down = abs( delta_down )

    #list_variations = []
    #list_variations.append(abs( up.Integral(0,up.GetNbinsX()+1) - integral_nominal ))
    #list_variations.append(abs( down.Integral(0,down.GetNbinsX()+1) - integral_nominal ))
    #
    #integral_sys_error  = max(list_variations)

    #print ("\t\tIntegral = {0} +- {1} (stat) +- {2} (syst)".format(integral_nominal, integral_stat_error, integral_sys_error ))
    print ("\t\tIntegral = {0} +- {1} (stat) (+ {2}, - {3}) (syst)".format(integral_nominal, integral_stat_error, integral_sys_up, integral_sys_down ))
  else:
    print ("\t\tIntegral = {0} +- {1} (stat)".format(integral_nominal, integral_stat_error ))


def main():

  doClosure = False

  if not doClosure:

    # -----------
    # DATA-DRIVEN
    #------------

    #inputpath = "./OutputPlots_MM_TwoLepLowNJetCR_25ns_v7_FinalSelection_NominalBinning/"
    #region    = "SS_LowNJetCR_DataDriven"
    #var_name  = "NJets2j3j4j"

    #inputpath = "./OutputPlots_MM_TwoLepSR_25ns_v7_FinalSelection_NominalBinning/"
    #region    = "SS_SR_DataDriven"
    #var_name  = "NJets5j"

    #inputpath = "./OutputPlots_MM_TwoLepSR_25ns_v7_FinalSelection_NominalBinning_AvgMuFake/"
    #region    = "SS_SR_DataDriven"
    #var_name  = "NJets5j"

    # -----------------------------------------

    #inputpath = "./OutputPlots_MM_TwoLepLowNJetCR_25ns_v7_FinalSelection_DDQmisID/"
    #region    = "SS_SR_DataDriven"
    #var_name  = "NJets2j3j4j"

    inputpath = "./OutputPlots_MM_TwoLepSR_25ns_v7_FinalSelection_DDQMisID/"
    region    = "SS_SR_DataDriven"
    var_name  = "NJets5j"

  else:

    # -------------
    # TTBAR CLOSURE
    #--------------

    #inputpath = "./OutputPlots_MM_MMClosureTest_v029_Baseline_Mllgt40GeV_AllElEtaCut/"
    #region    = "SS_SR_HighJet_DataDriven_Closure"
    #var_name  = "NJets4j5j"
    ##region	= "SS_SR_LowJet_DataDriven_Closure"
    ##var_name  = "NJets2j3j"

    inputpath = "./OutputPlots_MM_MMClosureTest_v029_Baseline_Mllgt40GeV_AllElEtaCut_EffNoLepIso/"
    region    = "SS_SR_HighJet_DataDriven_Closure"
    var_name  = "NJets4j5j"
    ##region	= "SS_SR_LowJet_DataDriven_Closure"
    ##var_name  = "NJets2j3j"

    #inputpath = "./OutputPlots_MM_MMClosureTest_v029_NoLepIso_Mllgt40GeV_AllElEtaCut/"
    #region    = "SS_SR_HighJet_DataDriven_Closure"
    #var_name  = "NJets4j5j"
    ##region	= "SS_SR_LowJet_DataDriven_Closure"
    ##var_name  = "NJets2j3j"

    #inputpath = "./OutputPlots_MM_MMClosureTest_v029_NoLepIP_Mllgt40GeV_AllElEtaCut/"
    #region    = "SS_SR_HighJet_DataDriven_Closure"
    #var_name  = "NJets4j5j"
    ##region	= "SS_SR_LowJet_DataDriven_Closure"
    ##var_name  = "NJets2j3j"

  flavour_list = ["ElEl", "MuMu", "OF"]

  for flav in flavour_list:

    print ("\nFlavour region: {0}\n".format(flav))

    filename = inputpath + flav + region + "/" + flav + region + "_" + var_name + ".root"
    myfile = TFile(filename)

    print("Looking at file: {0}".format(filename))

    fakes_nominal = myfile.Get("fakesbkg")
    fakes_up	  = myfile.Get("fakesbkg_MMfsys_up")
    fakes_down    = myfile.Get("fakesbkg_MMfsys_down")

    print ("\n\tFakes: \n")
    get_yields(fakes_nominal,fakes_up,fakes_down)

    if doClosure:
    	ttbar_nominal = myfile.Get("ttbarbkg")

    	print ("\n\tTTbar: \n")
    	get_yields(ttbar_nominal)

    expected_nominal = myfile.Get("expected")
    expected_up	     = myfile.Get("expected_MMfsys_up")
    expected_down    = myfile.Get("expected_MMfsys_down")

    print ("\n\tExpected: \n")
    get_yields(expected_nominal,expected_up,expected_down)

    if not doClosure:
    	observed = myfile.Get("observed")
        if observed:
    	   print ("\n\tObserved: \n")
    	   get_yields(observed)

    	signal = myfile.Get("signal")

    	print ("\n\tSignal: \n")
    	get_yields(signal)

# -----------------------------------------------------------------------------------


main()


