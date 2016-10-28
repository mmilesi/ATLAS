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

# Store sys integral for each systematic name
g_sys_dict = {}

# Store list of sys integrals for each sys group 
g_sysgroup_dict = {}

def get_yields(nominal, up=None, down=None, sysname=None, sysgroup=None):

    # Pick also O-Flow bin (NB: the last bin contains also the OFlow! Thus, stop at bin before oflow)
    #
    for bin in range(0,nominal.GetNbinsX()+1):

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

            #if nominal.IsBinOverflow(bin):
            #    print ("\t\t{0}-jets bin (O-FLOW): integral = {1:.2f} +- {2:.2f} (stat) (+ {3:.2f}, - {4:.2f}) (syst: {5})".format( bincenter, value_nominal, stat_error, sys_up, sys_dn, sysname ))
            #else:
            #    print ("\t\t{0}-jets bin: integral = {1:.2f} +- {2:.2f} (stat) (+ {3:.2f}, - {4:.2f}) (syst: {5})".format( bincenter, value_nominal, stat_error, sys_up, sys_dn, sysname ))

        else:
            if nominal.IsBinOverflow(bin):
                print ("\t\t{0}-jets bin (O-FLOW): integral = {1:.2f} +- {2:.2f} (stat)".format( bincenter, value_nominal, stat_error ))
            else:
                print ("\t\t{0}-jets bin: integral = {1:.2f} +- {2:.2f} (stat)".format( bincenter, value_nominal, stat_error ))

    integral_stat_error = Double(0)
    integral_nominal    = nominal.IntegralAndError(0,nominal.GetNbinsX(),integral_stat_error)

    integral_total_error = integral_stat_error

    if ( up and down ):

        integral_sys_up = abs( up.Integral(0,up.GetNbinsX()) - integral_nominal )
        integral_sys_dn = abs( integral_nominal - down.Integral(0,down.GetNbinsX()) )

	# Symmetrised systematic uncertainty
	#
	simm_sys_unc = abs( integral_sys_up + integral_sys_dn ) / 2.0

        g_sys_dict[sysname] = simm_sys_unc
        
	if not g_sysgroup_dict.get(sysgroup):
	    g_sysgroup_dict[sysgroup] = [simm_sys_unc]
        else:
	    g_sysgroup_dict[sysgroup].append(simm_sys_unc)
	
	# Total uncertainty
        #
        max_integral_sys   = max([integral_sys_up, integral_sys_dn])
        integral_tot_error = math.sqrt( ( integral_stat_error * integral_stat_error ) + ( max_integral_sys * max_integral_sys ) )

        # This will print the yield w/ stat and syst error per each bin
	
        #print ("\t\tIntegral = {0:.2f} +- {1:.2f} (stat) ( +{2:.2f}, -{3:.2f} --> +- {4:.2f}) (syst: {5})".format(integral_nominal, integral_stat_error, integral_sys_up, integral_sys_dn, simm_sys_unc,  sysname ))

    else:
        print ("\t\tIntegral = {0:.2f} +- {1:.2f} (stat)".format(integral_nominal, integral_stat_error ))

    return integral_nominal, integral_stat_error


def sum_quad ( inlist ):
    sq = 0
    for elem in inlist:
    	sq +=  pow(elem,2.0)
    return math.sqrt(sq)


def printTotFakeUncertainty( nominal, stat, flav ):

    if ( "HIGHNJ" in args.channel ):
        if flav == "ElEl" : non_closure = 0.29 * nominal
        if flav == "OF"   : non_closure = 0.27 * nominal
        if flav == "MuMu" : non_closure = 0.20 * nominal
    elif ( "LOWNJ" in args.channel ):
        if flav == "ElEl" : non_closure = 0.34 * nominal
        if flav == "OF"   : non_closure = 0.22 * nominal
        if flav == "MuMu" : non_closure = 0.22 * nominal
    
    if args.doClosure:
        non_closure = 0.0
    
    # This prints out sorting systematics from smaller to larger
    #
    print ("\t\tIntegral = {0:.2f}\n\t\t+- {1:.2f} [{2:.2f} %] (stat)\n\t\t+-".format(nominal, stat, (stat/nominal)*100) + "\t\t+-".join( " {0:.4f} [{1:.4f} %] ({2}) \n".format( g_sys_dict[key], (g_sys_dict[key]/nominal)*100, key ) for key in sorted( g_sys_dict, key=g_sys_dict.get ) ) + "\t\t+- {0:.2f} [{1:.2f} %] (non-closure)".format(non_closure, (non_closure/nominal)*100) )
   
    print("")
    
    # Print sum in quadrature of syst for each syst group
    
    print ("\t\tIntegral = {0:.2f} +- {1:.2f} [{2:.2f} %] (stat)".format(nominal, stat, (stat/nominal)*100) )
    for sg, values in g_sysgroup_dict.iteritems():
    
        print("\t\t+- {0:.2f} [{1:.2f} %] ({2})".format( sum_quad(values), (sum_quad(values)/nominal)*100, sg ) )
    
    print("")

    # Print sum in quadrature of ALL syst + stat
    
    toterrlist = list(g_sys_dict.values())
    toterrlist.extend([stat,non_closure])
    sq = sum_quad( toterrlist )

    print ("\t\tIntegral = {0:.2f} +- {1:.2f} [{2:.2f} %] (TOTAL UNCERTAINTY)".format(nominal, sq, (sq/nominal)*100))


def clearDicts():
    g_sys_dict.clear()
    g_sysgroup_dict .clear()

if __name__ == '__main__':

    region = var_name = None

    inputpath = args.inputDir

    if not inputpath.endswith('/'):
         inputpath += '/'

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
        
	clearDicts()
	
	print ("\nFlavour region: {0}\n".format(flav))

        filename = inputpath + flav + region + "/" + flav + region + "_" + var_name + ".root"
        myfile = TFile(filename)

        print("Looking at file: {0}".format(filename))

        fakes_nominal = myfile.Get("fakesbkg")
        fakes_syst = {}
	for key in myfile.GetListOfKeys():
	    keyname = key.GetName()
	    if not ( "fakesbkg_" in keyname ): continue
            keyname = keyname.replace("fakesbkg_","")
	    keyname = keyname.replace("_dn","")
	    keyname = keyname.replace("_up","")
	    if "Stat" in keyname:
	        value = "Stat"
	    if "numerator_QMisID" in keyname:
	        value = "numerator_QMisID"
	    if "denominator_QMisID" in keyname:
	        value = "denominator_QMisID"
            if not fakes_syst.get(keyname):
	        fakes_syst[keyname] = value
	
        print ("\n\tFakes: \n")

	for sys, sysgroup in fakes_syst.iteritems():
	   fakes_up   = myfile.Get( "fakesbkg_" + sys + "_up")
	   fakes_down = myfile.Get( "fakesbkg_" + sys + "_dn")
	   #print(" ==> sys: {0}, sysgroup: {1}\n".format(sys, sysgroup))
           fakes, fakes_err = get_yields(fakes_nominal,fakes_up,fakes_down, sys, sysgroup)
        printTotFakeUncertainty( fakes, fakes_err, flav )

        if args.doClosure:
            ttbar_nominal = myfile.Get("ttbarbkg")

            print ("\n\tTTbar: \n")
            ttbar, ttbar_err = get_yields(ttbar_nominal)

            closure     = ( (fakes - ttbar) / (ttbar) ) * 100
            closure_err = math.sqrt( ( ( math.pow(fakes,2.0) / math.pow(ttbar,4.0) ) * math.pow(ttbar_err,2.0) ) + ( math.pow(fakes_err,2.0) / math.pow(ttbar,2.0) ) ) * 100

            print("\nNon-closure ((fakes-ttbar)/ttbar) = {0:.2f} [%] +- {1:.2f} [%]".format(closure,closure_err))

        expected_nominal = myfile.Get("expectedbkg")

        print ("\n\tExpected: \n")
        exp, exp_err = get_yields(expected_nominal)

        if not args.doClosure:

            chargemisid = myfile.Get("qmisidbkg")
            if chargemisid:
                print ("\n\tQMisID: \n")
                qmisid, qmisid_err = get_yields(chargemisid)

            observed = myfile.Get("observed")
            if observed:
                print ("\n\tObserved: \n")
                obs, obs_err = get_yields(observed)

            signal = myfile.Get("signal")

            print ("\n\tSignal: \n")
            sig, sig_err = get_yields(signal)
