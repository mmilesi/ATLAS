 #!/usr/bin/python

import os, sys, math, argparse

parser = argparse.ArgumentParser(description='Get yields and systematics for MM')

list_channel = ['HIGHNJ','LOWNJ','ALLNJ']

g_luminosities = { "Moriond 2016 GRL":3.209,            # March 2016
                 "ICHEP 2015+2016 DS":13.20768,       # August 2016
                 "POST-ICHEP 2015+2016 DS":22.07036,  # October 2016
                 "FULL 2015+2016 DS":36.4702          # December 2016
                 }

parser.add_argument('inputDir', metavar='inputDir',type=str,
                    help='Path to the directory containing input histograms')
parser.add_argument('--channel', dest='channel', action='store', default='HIGHNJ', type=str, nargs='+',
                    help='The channel chosen. Full list of available options:\n{0}'.format(list_channel))
parser.add_argument('--variables', dest='variables', action='store', type=str, nargs='*',
                    help='List of variables to be considered. Use a space-separated list. If unspecified, will consider Njets only.')
parser.add_argument('--closure', dest='closure', action='store_true', default=False,
                    help="Check yields for MC closure test")
parser.add_argument("--lumi", dest="lumi", action="store", type=float, default=g_luminosities["FULL 2015+2016 DS"],
                    help="The luminosity of the dataset. Pick one of these values: ==> " + ",".join( "{0} ({1})".format( lumi, tag ) for tag, lumi in g_luminosities.iteritems() ) )

args = parser.parse_args()

from ROOT import gROOT, TCanvas, TPad, TH1D, TFile, TLegend, TLatex, TLine, Double, kTeal, kGray

gROOT.Reset()
gROOT.LoadMacro("$HOME/RootUtils/AtlasStyle.C")
from ROOT import SetAtlasStyle
SetAtlasStyle()


# Store sys integral for each systematic name
g_sys_dict = {}

# Store for each histogram bin a list with the uncertainties for each source
#
# { bin_idx : [unc0, unc1, unc2,...] }
#
g_unc_bins = {}
# Store for each histogram bin the sum in quadrature of all uncertainties
#
# { bin_idx : sq }
#
g_sq_unc_bins = {}

# Store list of sys integrals for each sys group
g_sysgroup_dict = {}

def get_yields(nominal, up=None, down=None, sysname=None, sysgroup=None):

    for bin in range(0,nominal.GetNbinsX()+1):
        if not g_unc_bins.get(bin):
            g_unc_bins[bin] = []

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

            sys_up = sys_dn = 0.0

            if delta_up >= 0.0:
                sys_up = abs( delta_up )
            if delta_down <= 0.0:
                sys_dn = abs( delta_down )

	    # Symmetrised systematic uncertainty
	    #
            simm_sys_unc =  abs( sys_up + sys_dn ) / 2.0

            # Print yield, stat and syst uncertainty for each bin, for this systematic
            #
            if False:
                if nominal.IsBinOverflow(bin): # this condition should never be matched: here just for safety
                    print ("\t\t{0}-jets bin (O-FLOW): integral = {1:.3f} +- {2:.3f} (stat) (+ {3:.3f}, - {4:.3f} --> +- {5:.3f}) (syst: {6})".format( bincenter, value_nominal, stat_error, sys_up, sys_dn, simm_sys_unc, sysname ))
                else:
                    print ("\t\t{0}-jets bin: integral = {1:.3f} +- {2:.3f} (stat) (+ {3:.3f}, - {4:.3f} --> +- {5:.3f}) (syst: {6})".format( bincenter, value_nominal, stat_error, sys_up, sys_dn, simm_sys_unc, sysname ))

            # Store list of uncertainties for each bin
	    #
            if ( float("{0:.3f}".format(stat_error)) ) and not stat_error in g_unc_bins[bin]:
                g_unc_bins[bin].append(stat_error)
            if ( float("{0:.3f}".format(simm_sys_unc)) ) and not simm_sys_unc in g_unc_bins[bin]:
                g_unc_bins[bin].append(simm_sys_unc)

        else:
            if nominal.IsBinOverflow(bin): # this condition should never be matched: here just for safety
                print ("\t\t{0}-jets bin (O-FLOW): integral = {1:.3f} +- {2:.3f} (stat)".format( bincenter, value_nominal, stat_error ))
            else:
                print ("\t\t{0}-jets bin: integral = {1:.3f} +- {2:.3f} (stat)".format( bincenter, value_nominal, stat_error ))

    integral_stat_error = Double(0)
    integral_nominal    = nominal.IntegralAndError(0,nominal.GetNbinsX(),integral_stat_error)

    integral_total_error = integral_stat_error

    if ( up and down ):

        integral_sys_up = abs( up.Integral(0,up.GetNbinsX()) - integral_nominal )
        integral_sys_dn = abs( integral_nominal - down.Integral(0,down.GetNbinsX()) )

	# Symmetrised systematic uncertainty
	#
	integral_simm_sys_unc = abs( integral_sys_up + integral_sys_dn ) / 2.0

        g_sys_dict[sysname] = integral_simm_sys_unc

	# Store the total syst uncertainty for later use
	#
	if not g_sysgroup_dict.get(sysgroup):
	    g_sysgroup_dict[sysgroup] = [integral_simm_sys_unc]
        else:
	    g_sysgroup_dict[sysgroup].append(integral_simm_sys_unc)

	# Total uncertainty
        #
        max_integral_sys   = max([integral_sys_up, integral_sys_dn])
        integral_tot_error = math.sqrt( ( integral_stat_error * integral_stat_error ) + ( max_integral_sys * max_integral_sys ) )

        # This will print the total yield w/ stat and syst error per each sys
        #
        if False:
            print ("\t\tIntegral = {0:.3f} +- {1:.3f} (stat) ( +{2:.3f}, -{3:.3f} --> +- {4:.3f}) (syst: {5})".format(integral_nominal, integral_stat_error, integral_sys_up, integral_sys_dn, integral_simm_sys_unc,  sysname ))

        # This will sum in quadrature the stat & systematic uncertainties for each bin
        #
        for bin, list_unc in g_unc_bins.iteritems():
            if False:
                print("\t\tbin {0}".format(bin) + " list of uncertainties: [" + ",".join( "{0:.3f}".format(x) for x in list_unc ) + "]" )
            g_sq_unc_bins[bin] = sum_quad(list_unc)

        if False:
            for bin, sq in g_sq_unc_bins.iteritems():
                print("\t\tbin {0}".format(bin) + " sum quad: {0:.3f}".format(sq) )

    else:
        print ("\t\tIntegral = {0:.3f} +- {1:.3f} (stat)".format(integral_nominal, integral_stat_error ))

    return integral_nominal, integral_stat_error


def sum_quad ( inlist ):
    sq = 0
    for elem in inlist:
    	sq +=  pow(elem,2.0)
    return math.sqrt(sq)


def getTotFakeUncertainty( nominal, stat, flav ):

    if ( "HIGHNJ" in args.channel ):
        if flav == "ElEl" : non_closure = 0.29 * nominal
        if flav == "OF"   : non_closure = 0.27 * nominal
        if flav == "MuMu" : non_closure = 0.20 * nominal
    elif ( "LOWNJ" in args.channel ):
        if flav == "ElEl" : non_closure = 0.34 * nominal
        if flav == "OF"   : non_closure = 0.22 * nominal
        if flav == "MuMu" : non_closure = 0.22 * nominal

    if args.closure:
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

    return sq


def clearDicts():
    g_sys_dict.clear()
    g_sysgroup_dict .clear()
    g_unc_bins.clear()
    g_sq_unc_bins.clear()


def makeSysPlots( flav, var, MC_hist, MM_hist ):

    gROOT.SetBatch(True)

    c = TCanvas("c1","Temp",50,50,800,600)

    pad1 = TPad("pad1", "", 0, 0.25, 1, 1)
    pad2 = TPad("pad2", "", 0, 0,   1, 0.25)
    pad1.SetBottomMargin(0.02)
    pad2.SetBottomMargin(0.4)
    pad1.Draw()
    pad2.Draw()

    new_MM_hist = MM_hist.Clone(MM_hist.GetName())

    # For MM hist, set the bin error as the sum quad of stat+syst
    #
    for ibin in range(1,new_MM_hist.GetNbinsX()+1):
        new_MM_hist.SetBinError(ibin, g_sq_unc_bins[ibin] )

    MC_hist.SetLineStyle(2)
    MC_hist.SetMarkerSize(0.8)
    MC_hist.SetLineColor(1)
    MC_hist.SetMarkerStyle(20)
    MC_hist.SetLineWidth(1)

    new_MM_hist.SetLineWidth(3)
    new_MM_hist.SetLineStyle(1)
    new_MM_hist.SetLineColor(1)
    new_MM_hist.SetFillColor(kTeal+1)
    new_MM_hist.SetFillStyle(1001)

    err = new_MM_hist.Clone("tot_uncertainty")
    err.SetFillColor(kGray+3)
    err.SetLineColor(10)
    err.SetFillStyle(3004)
    err.SetMarkerSize(0)

    # Trick to rescale:
    #
    ymax      = new_MM_hist.GetMaximum()
    binmax    = new_MM_hist.GetMaximumBin()
    binmaxerr = new_MM_hist.GetBinError(binmax)

    #print("FLAV: {0} - MM max = {1:.2f} - ttbar max = {2:.2f}".format(flav,new_MM_hist.GetMaximum(),MC_hist.GetMaximum()))

    if MC_hist.GetMaximum() > ymax:
        ymax      = MC_hist.GetMaximum()
        binmax    = MC_hist.GetMaximumBin()
        binmaxerr = MC_hist.GetBinError(binmax)

    #print("FLAV: {0} - max = {1:.2f} - binmax = {2} - binmaxerr = {3:.2f}".format(flav,ymax,binmax,binmaxerr))

    new_MM_hist.SetMaximum( ymax + binmaxerr *1.3 )
    new_MM_hist.SetMinimum(0)

    # --------------------------

    # Stat +sys error on MM (ratio)

    ratio_err = err.Clone("RatioErr")
    ratio_err.SetXTitle(MC_hist.GetXaxis().GetTitle())
    ratio_err.SetYTitle("")
    ratio_err.GetXaxis().SetTitleSize(0.15)
    ratio_err.GetYaxis().SetTitleSize(0.15)
    ratio_err.GetXaxis().SetTitleOffset(0.90)
    ratio_err.GetYaxis().SetTitleOffset(0.35)
    ratio_err.GetXaxis().SetLabelSize(0.15)
    ratio_err.GetYaxis().SetLabelSize(0.12)
    ratio_err.GetYaxis().SetNdivisions(505) #(5)
    ratio_err.SetFillColor(kGray+3)
    ratio_err.SetLineColor(10)
    ratio_err.SetFillStyle(3004)
    ratio_err.SetMarkerSize(0)

    ratio_err.Divide(err)

    # MC ttbar / MM ttbar

    ratio_MC_MM = MC_hist.Clone("RatioMCMM")
    ratio_MC_MM.SetYTitle("")
    ratio_MC_MM.SetLineStyle(2)
    ratio_MC_MM.SetMarkerSize(0.8)
    ratio_MC_MM.SetLineColor(1)
    ratio_MC_MM.SetMarkerStyle(20)
    ratio_MC_MM.SetLineWidth(1)
    ratio_MC_MM.Divide(new_MM_hist)

    # Trick to rescale

    ymax_ratio = ratio_MC_MM.GetMaximum()
    ymin_ratio = ratio_MC_MM.GetMinimum()

    #print("FLAV: {0} - MC/MM max = {1:.2f} - MC/MM min = {2:.2f}".format(flav,ymax_ratio,ymin_ratio))

    # For ratio error, get the bin with maximum uncertainty and the uncertainty value

    max_unc = 0
    for ibin in range(1, ratio_err.GetNbinsX()+1):
        this_unc = ratio_err.GetBinError(ibin)
      	if this_unc > max_unc:
	    max_unc     = this_unc

    #print("FLAV: {0} - max unc./2 ratio = {1:.2f}".format(flav,max_unc/2.0))

    if ( 1 + max_unc/2.0 ) > ymax_ratio:
    	ymax_ratio = 1 + max_unc/2.0
    if ( 1 - max_unc/2.0 ) < ymin_ratio:
    	ymin_ratio = 1 - max_unc/2.0

    if ymin_ratio < 0: ymin_ratio = 0

    #print("FLAV: {0} - actual max = {1:.2f} - actual min = {2:.2f}".format(flav,ymax_ratio,ymin_ratio))

    ratio_err.SetMaximum( ymax_ratio * 1.2 )
    ratio_err.SetMinimum( ymin_ratio * 0.8 )

    # --------------------------

    pad1.cd()

    new_MM_hist.GetXaxis().SetLabelSize(0)
    new_MM_hist.GetXaxis().SetLabelOffset(999)

    new_MM_hist.Draw("HIST")
    MC_hist.Draw("PE SAME")
    err.Draw("E2 SAME")

    legend = TLegend(0.6,0.7,0.8,0.9) # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
    legend.SetBorderSize(0)	# no border
    legend.SetFillStyle(0)	# Legend transparent background
    legend.SetTextSize(0.03)	# Increase entry font size!
    legend.SetTextFont(42)	# Helvetica

    legend.AddEntry(new_MM_hist, "FakesMM - t#bar{t}", "F")
    legend.AddEntry(MC_hist, "t#bar{t}", "P")
    legend.AddEntry(err, "Stat.+Sys. Unc.", "F")

    leg_ATLAS = TLatex()
    leg_lumi  = TLatex()
    leg_ATLAS.SetTextSize(0.03)
    leg_ATLAS.SetNDC()
    leg_lumi.SetTextSize(0.03)
    leg_lumi.SetNDC()

    legend.Draw()
    leg_ATLAS.DrawLatex(0.6,0.55,"#bf{#it{ATLAS}} Work In Progress")
    leg_lumi.DrawLatex(0.6,0.45,"#sqrt{{s}} = 13 TeV, #int L dt = {0:.1f} fb^{{-1}}".format(args.lumi))

    pad2.cd()
    #ratio_err.GetYaxis().SetRangeUser(0.0, 2.0)
    ratio_err.GetYaxis().SetRangeUser(0.5, 1.5)
    ratio_err.Draw("E2")
    ratio_MC_MM.Draw("PE SAME")

    refl = TLine(ratio_err.GetBinLowEdge(1), 1.0, ratio_err.GetBinLowEdge(ratio_err.GetNbinsX()+1), 1.0)
    refl.SetLineStyle(2)
    refl.SetLineWidth(2)
    refl.Draw("SAME")

    # --------------------------

    outpath = args.inputDir
    if outpath[-1] == '/':
      outpath = outpath[:-1]

    c.SaveAs( outpath + "/" + flav + "_" + var + "_Sys.png" )

# -------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    region = var_name = None

    inputpath = args.inputDir

    if not inputpath.endswith('/'):
         inputpath += '/'

    var_list = []
    if args.variables:
        var_list = args.variables

    if not args.closure:

        # -----------
        # DATA-DRIVEN
        #------------

        if ( "HIGHNJ" in args.channel ):

            region    = "SS_SR_DataDriven"
	    if not "NJets5j" in var_list:
	        var_list.append("NJets5j")

        elif ( "LOWNJ" in args.channel ):

            region    = "SS_LowNJetCR_DataDriven"
	    if not "NJets2j3j4j" in var_list:
  	        var_list.append("NJets2j3j4j")

    else:

        # -------------
        # TTBAR CLOSURE
        #--------------

        if ( "LOWNJ" in args.channel ):

            region    = "SS_SR_LowJet_DataDriven_Closure"
	    if not "NJets2j3j4j" in var_list:
  	        var_list.append("NJets2j3j4j")

        elif ( "HIGHNJ" in args.channel ):

            region    = "SS_SR_HighJet_DataDriven_Closure"
	    if not "NJets5j" in var_list:
	        var_list.append("NJets5j")

        elif ( "ALLNJ" in args.channel ):

            region    = "SS_SR_AllJet_DataDriven_Closure"
	    if not "NJets" in var_list:
	        var_list.append("NJets")

    flavour_list = ["ElEl", "MuMu", "OF"]
    #flavour_list = ["ElEl"]

    print("Looking at variables : [" + ",".join( "{0}".format(v) for v in var_list ) + "]")

    for var_name in var_list:

	print ("\nVariable: {0}\n".format(var_name))

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
    	       fakes, fakes_stat_err = get_yields(fakes_nominal,fakes_up,fakes_down, sys, sysgroup)

    	    fakes_tot_err = getTotFakeUncertainty( fakes, fakes_stat_err, flav )

    	    if args.closure:
    		ttbar_nominal = myfile.Get("ttbarbkg")

    		print ("\n\tTTbar: \n")
    		ttbar, ttbar_err = get_yields(ttbar_nominal)

    		closure     = ( (fakes - ttbar) / (ttbar) ) * 100
    		closure_err = math.sqrt( ( ( math.pow(fakes,2.0) / math.pow(ttbar,4.0) ) * math.pow(ttbar_err,2.0) ) + ( math.pow(fakes_tot_err,2.0) / math.pow(ttbar,2.0) ) ) * 100

    		makeSysPlots(flav,var_name,ttbar_nominal,fakes_nominal)

    		print("\nNon-closure ((fakes-ttbar)/ttbar) = {0:.2f} [%] +- {1:.2f} [%]".format(closure,closure_err))

    	    expected_nominal = myfile.Get("expectedbkg")

    	    print ("\n\tExpected: \n")
    	    exp, exp_err = get_yields(expected_nominal)

    	    if not args.closure:

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
