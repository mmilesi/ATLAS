 #!/usr/bin/python

import os, sys, math, argparse

parser = argparse.ArgumentParser(description='Get yields and systematics for MM')

list_channel = ['HIGHNJ','LOWNJ','ALLNJ']

parser.add_argument('inputDir', metavar='inputDir',type=str,
                   help='Path to the directory containing input histograms')
parser.add_argument('--channel', dest='channel', action='store', default='HIGHNJ', type=str, nargs='+',
                    help='The channel chosen. Full list of available options:\n{0}'.format(list_channel))
parser.add_argument('--doClosure', dest='doClosure', action='store_true', default=False,
                    help="Check yields for MC closure test")

args = parser.parse_args()

from ROOT import gROOT, TH1D, TFile, Double

gROOT.SetBatch(True)

def get_yields(nominal, up=None, down=None):

    #print("\t\tGetNbinsX() = {0}".format(nominal.GetNbinsX()))
    #print("\t\trange(0,GetNbinsX()+1) = {0}".format(range(0,nominal.GetNbinsX()+1)))

    # Pick also O-Flow bin (NB: the last bin contains also the OFlow!)
    #
    for bin in range(0,nominal.GetNbinsX()+1):#+2):

        nextbin = bin

        stat_error    = Double(0)
        value_nominal = nominal.IntegralAndError(bin,nextbin,stat_error)

        bincenter = nominal.GetBinCenter(bin)

        if ( up and down ):

            dummy_err_up  = Double(0)
            value_up = up.IntegralAndError(bin,nextbin,dummy_err_up)
            dummy_err_down  = Double(0)
            value_down = down.IntegralAndError(bin,nextbin,dummy_err_down)

            delta_up   = value_up - value_nominal
            delta_down = value_down - value_nominal

            sys_up   = 0.0
            sys_dn = 0.0

            if delta_up >= 0.0:
                sys_up = abs( delta_up )
            if delta_down <= 0.0:
                sys_dn = abs( delta_down )

            if nominal.IsBinOverflow(bin):
                print ("\t\t{0}-jets bin (O-FLOW): integral = {1:.2f} +- {2:.2f} (stat) (+ {3:.2f}, - {4:.2f}) (syst)".format( bincenter, value_nominal, stat_error, sys_up, sys_dn ))
            else:
                print ("\t\t{0}-jets bin: integral = {1:.2f} +- {2:.2f} (stat) (+ {3:.2f}, - {4:.2f}) (syst)".format( bincenter, value_nominal, stat_error, sys_up, sys_dn ))

        else:
            if nominal.IsBinOverflow(bin):
                print ("\t\t{0}-jets bin (O-FLOW): integral = {1:.2f} +- {2:.2f} (stat)".format( bincenter, value_nominal, stat_error ))
            else:
                print ("\t\t{0}-jets bin: integral = {1:.2f} +- {2:.2f} (stat)".format( bincenter, value_nominal, stat_error ))

    integral_stat_error = Double(0)
    integral_nominal    = nominal.IntegralAndError(0,nominal.GetNbinsX(),integral_stat_error)

    integral_total_error = integral_stat_error

    print ("\t\t--------------------")
    if ( up and down ):

        delta_up   = up.Integral(0,up.GetNbinsX()+1) - integral_nominal
        delta_down = down.Integral(0,down.GetNbinsX()+1) - integral_nominal

        integral_sys_up   = 0.0
        integral_sys_dn = 0.0

        if delta_up >= 0.0:
            integral_sys_up = abs( delta_up )
        if delta_down <= 0.0:
            integral_sys_dn = abs( delta_down )

        # Total uncertainty
        #
        max_integral_sys   = max([integral_sys_up, integral_sys_dn])
        integral_tot_error = math.sqrt( ( integral_stat_error * integral_stat_error ) + ( max_integral_sys * max_integral_sys ) )

        print ("\t\tIntegral = {0:.2f} +- {1:.2f} (stat) (+ {2:.2f}, - {3:.2f}) (syst)".format(integral_nominal, integral_stat_error, integral_sys_up, integral_sys_dn ))
        print ("\t\t         = {0:.2f} +- {1:.2f} (TOTAL UNCERTAINTY)".format(integral_nominal, integral_tot_error ))

    else:
        print ("\t\tIntegral = {0:.2f} +- {1:.2f} (stat)".format(integral_nominal, integral_stat_error ))

    return integral_nominal, integral_total_error

if __name__ == '__main__':

    region = var_name = None

    inputpath = args.inputDir

    if not args.doClosure:

        # -----------
        # DATA-DRIVEN
        #------------

        if ( "HIGHNJ" in args.channel ):

            region    = "SS_SR_DataDriven"
            var_name  = "NJets5j"

        elif ( "LOWNJ" in args.channel ):

            region    = "SS_LowNJetCR_DataDriven"
            var_name  = "NJets2j3j4j"

    else:

        # -------------
        # TTBAR CLOSURE
        #--------------

        if ( "LOWNJ" in args.channel ):

            region    = "SS_SR_LowJet_DataDriven_Closure"
            var_name  = "NJets2j3j4j"

        elif ( "HIGHNJ" in args.channel ):

            region    = "SS_SR_HighJet_DataDriven_Closure"
            var_name  = "NJets5j"

        elif ( "ALLNJ" in args.channel ):

            region    = "SS_SR_AllJet_DataDriven_Closure"
            var_name  = "NJets"

    flavour_list = ["ElEl", "MuMu", "OF"]

    for flav in flavour_list:

        print ("\nFlavour region: {0}\n".format(flav))

        filename = inputpath + flav + region + "/" + flav + region + "_" + var_name + ".root"
        myfile = TFile(filename)

        print("Looking at file: {0}".format(filename))

        fakes_nominal = myfile.Get("fakesbkg")
        fakes_up      = myfile.Get("fakesbkg_MMfsys_up")
        fakes_down    = myfile.Get("fakesbkg_MMfsys_dn")

        print ("\n\tFakes: \n")
        fakes, fakes_err = get_yields(fakes_nominal,fakes_up,fakes_down)

        if args.doClosure:
            ttbar_nominal = myfile.Get("ttbarbkg")

            print ("\n\tTTbar: \n")
            ttbar, ttbar_err = get_yields(ttbar_nominal)

            closure     = ( (fakes - ttbar) / (ttbar) ) * 100
            closure_err = math.sqrt( ( ( math.pow(fakes,2.0) / math.pow(ttbar,4.0) ) * math.pow(ttbar_err,2.0) ) + ( math.pow(fakes_err,2.0) / math.pow(ttbar,2.0) ) ) * 100

            print("\nNon-closure ((fakes-ttbar)/ttbar) = {0:.2f} [%] +- {1:.2f} [%]".format(closure,closure_err))

        expected_nominal = myfile.Get("expected")
        expected_up      = myfile.Get("expected_MMfsys_up")
        expected_down    = myfile.Get("expected_MMfsys_dn")

        print ("\n\tExpected: \n")
        exp, exp_err = get_yields(expected_nominal,expected_up,expected_down)

        if not args.doClosure:
            observed = myfile.Get("observed")
            if observed:
                print ("\n\tObserved: \n")
                obs, obs_err = get_yields(observed)

            signal = myfile.Get("signal")

            print ("\n\tSignal: \n")
            sig, sig_err = get_yields(signal)
