 #!/usr/bin/python

import os, sys, math, argparse

sys.path.append(os.path.abspath(os.path.curdir))

parser = argparse.ArgumentParser(description='Get yields and systematics for MM')

list_channel = ['HIGHNJ','LOWNJ','ALLNJ']

categories   = ["ALL","ElEl","MuMu","OF","Inclusive"]

luminosities = { "Moriond 2016 GRL":3.209,            # March 2016
                 "ICHEP 2015+2016 DS":13.20768,       # August 2016
                 "POST-ICHEP 2015+2016 DS":22.07036,  # October 2016
                 "FULL 2015+2016 DS":36.0746          # December 2016
               }

parametrisations_fake_el = ["Pt","NBJets_VS_Pt"]
parametrisations_fake_mu = ["Pt","DistanceClosestJet_VS_Pt"]

parser.add_argument('inputDir', metavar='inputDir',type=str,
                    help='Path to the directory containing input histograms')
parser.add_argument('--channel', dest='channel', action='store', default='HIGHNJ', type=str, nargs='+',
                    help='The channel chosen. Full list of available options:\n{0}'.format(list_channel))
parser.add_argument('--category', dest='category', action='store', default=categories[0], type=str, nargs='+', choices=categories,
                    help='The category chosen. Can pass multiple space-separated arguments to this command-line option (picking amonge the above list). If this option is not specified, default will be \'{0}\' (NB: \'Inclusive\' is not included by default, no pun intended:))'.format(categories[0]))
parser.add_argument('--variables', dest='variables', action='store', type=str, nargs='*',
                    help='List of variables to be considered. Use a space-separated list.')
parser.add_argument('--paramFakeEl', dest='paramFakeEl', action='store', default=parametrisations_fake_el[1], const=parametrisations_fake_el[1], type=str, nargs='?', choices=parametrisations_fake_el,
                    help='The parametrisation for the electron fake rate. Option is needed to make sure the correct non-closure uncertainty is used. If this option is not specified, default will be \'{0}\''.format(parametrisations_fake_el[1]))
parser.add_argument('--paramFakeMu', dest='paramFakeMu', action='store', default=parametrisations_fake_mu[1], const=parametrisations_fake_mu[1], type=str, nargs='?', choices=parametrisations_fake_mu,
                    help='The parametrisation for the muon fake rate. Option is needed to make sure the correct non-closure uncertainty is used. If this option is not specified, default will be \'{0}\''.format(parametrisations_fake_mu[1]))
parser.add_argument('--closure', dest='closure', action='store_true', default=False,
                    help="Check yields for MC closure test")
parser.add_argument('--doKinematicsComparison', dest='doKinematicsComparison', action='store_true', default=False,
                    help="Compare DD Fakes Vs. MC Fakes distributions.")
parser.add_argument("--lumi", dest="lumi", action="store", type=float, default=luminosities["FULL 2015+2016 DS"],
                    help="The luminosity of the dataset. Pick one of these values: ==> " + ",".join( "{0} ({1})".format( lumi, tag ) for tag, lumi in luminosities.iteritems() ) + ". Default is {0}".format(luminosities["FULL 2015+2016 DS"] ) )
parser.add_argument("--mergeOverflow", dest="mergeOverflow", action="store_true", default=False,
                    help="If this option is used, the script will assume the overflow has been merged already w/ the last visible bin. Default is False.")
parser.add_argument('--doLogScaleY', dest='doLogScaleY', action='store_true', default=False,
                    help='Use log scale on the Y axis')
parser.add_argument('--debug', dest='debug', action='store_true', default=False,
                    help='Run in debug mode')

args = parser.parse_args()

from ROOT import gROOT, gStyle, gPad, TCanvas, TPad, TH1, TH1D, THStack, TFile, TLegend, TLatex, TLine, Double, kTeal, kGray, kBlack, kBlue, kOrange, kWhite, kViolet, kRed

gROOT.Reset()
gROOT.LoadMacro(os.path.abspath(os.path.curdir)+"/Plotter/AtlasStyle.C")
# gROOT.LoadMacro("$HOME/RootUtils/AtlasStyle.C")
from ROOT import SetAtlasStyle
SetAtlasStyle()

from Plotter.BackgroundTools import makePoissonErrors, integrate, SelfDivide

# Store sys integral for each systematic name

g_sys_dict = {}

# Store each histogram bin's content

g_content_bins = {}

# Store for each histogram bin a list with the uncertainties for each source
#
# { bin_idx : [ (unc0,sysname,sysgroup), (unc1,sysname,sysgroup), (unc2,sysname,sysgroup),...] }

g_unc_bins = {}
g_unc_bins_NO_STAT = {}

# Store for each histogram bin the sum in quadrature of all uncertainties
# { bin_idx : sq }

g_sq_unc_bins = {}
g_sq_unc_bins_NO_STAT = {}

# Store list of sys integrals for each sys group

g_sysgroup_dict = {}

def getYields(nominal, up=None, dn=None, sysname=None, sysgroup=None, debug=False):

    lower = 0
    upper = nominal.GetNbinsX()+2
    if args.mergeOverflow:
        upper = nominal.GetNbinsX()+1

    for bin in range(lower,upper):
        if not g_unc_bins.get(bin):
            g_unc_bins[bin] = []
            g_unc_bins_NO_STAT[bin] = []

    for bin in range(lower,upper):

        nextbin = bin

        stat_error    = Double(0)
        value_nominal = nominal.IntegralAndError(bin,nextbin,stat_error)

        g_content_bins[bin] = value_nominal

        bincenter = nominal.GetBinCenter(bin)

        if ( up and dn ):

            dummy_err_up = Double(0)
            value_up     = up.IntegralAndError(bin,nextbin,dummy_err_up)
            dummy_err_dn = Double(0)
            value_dn     = dn.IntegralAndError(bin,nextbin,dummy_err_dn)

            delta_up = value_up - value_nominal
            delta_dn = value_dn - value_nominal

            sys_up = sys_dn = 0.0

            if abs(delta_up) != abs(value_nominal) and delta_up >= 0.0:
                sys_up = abs( delta_up )
            if abs(delta_dn) != abs(value_nominal) and delta_dn <= 0.0:
                sys_dn = abs( delta_dn )

	    # Symmetrised systematic uncertainty for this bin

            simm_sys_unc =  abs( sys_up + sys_dn ) / 2.0

            # Print yield, stat and syst uncertainty for this bin, for this systematic

            if debug:
                if nominal.IsBinOverflow(bin):
                    if args.mergeOverflow:
                        print ("\nWARNING! Now checking the overflow bin content. This should not happen since you said in the script opt config that the last bin should already contain also the OFlow...check your input histograms.\n")
                    print ("\t\t{0}-th bin (O-FLOW), bincenter {1} : integral = {2:.3f} +- {3:.3f} (stat) +- {4:.3f} (syst: {5}) [ +{6:.3f}, -{7:.3f} ]".format( bin, bincenter, value_nominal, stat_error, simm_sys_unc, sysname, sys_up, sys_dn ))
                else:
                    print ("\t\t{0}-th bin, bincenter {1} : integral = {2:.3f} +- {3:.3f} (stat) +- {4:.3f} (syst: {5}) [ +{6:.3f}, -{7:.3f} ]".format( bin, bincenter, value_nominal, stat_error, simm_sys_unc, sysname, sys_up, sys_dn ))

            # Store list of uncertainties for each bin

            if ( float("{0:.3f}".format(stat_error)) ) and not (stat_error,"Statistical","Statistical") in g_unc_bins[bin]:
                g_unc_bins[bin].append( (stat_error,"Statistical","Statistical") )
            if ( float("{0:.3f}".format(simm_sys_unc)) ) and not (simm_sys_unc,sysname,sysgroup) in g_unc_bins[bin]:
                g_unc_bins[bin].append( (simm_sys_unc,sysname,sysgroup) )
                g_unc_bins_NO_STAT[bin].append( (simm_sys_unc,sysname,sysgroup) )
        else:
            if debug:
                if nominal.IsBinOverflow(bin):
                    if args.mergeOverflow:
                        print ("\nWARNING! Now checking the overflow bin content. This should not happen since you said in the script opt config that the last bin should already contain also the OFlow...check your input histograms.\n")
                    print ("\t\t{0}-th bin (O-FLOW), bincenter {1} : integral = {2:.3f} +- {3:.3f} (stat)".format( bin, bincenter, value_nominal, stat_error ))
                else:
                    print ("\t\t{0}-th bin, bincenter {1} : integral = {2:.3f} +- {3:.3f} (stat)".format( bin, bincenter, value_nominal, stat_error ))

    integral_stat_error = Double(0)
    integral_nominal    = nominal.IntegralAndError(lower,upper-1,integral_stat_error)

    integral_total_error = integral_stat_error

    if ( up and dn ):

        integral_sys_up = integral_sys_dn = 0.0

        if up.Integral(lower,upper-1):
            integral_sys_up = abs( up.Integral(lower,upper-1) - integral_nominal )
        if dn.Integral(lower,upper-1):
            integral_sys_dn = abs( integral_nominal - dn.Integral(lower,upper-1) )

	# Symmetrised total systematic uncertainty on yield for this sys source

	integral_simm_sys_unc = abs( integral_sys_up + integral_sys_dn ) / 2.0

        g_sys_dict[sysname] = integral_simm_sys_unc

	# Store the total syst uncertainty for later use
        # NB: multiple systematics may belong to the same group, e.g. "Real_Mu_Pt_Stat", "Fake_El_Pt_Stat" are part of group "Stat"...

	if not g_sysgroup_dict.get(sysgroup):
	    g_sysgroup_dict[sysgroup] = [integral_simm_sys_unc]
        else:
	    g_sysgroup_dict[sysgroup].append(integral_simm_sys_unc)

	# Total uncertainty

        max_integral_sys   = max([integral_sys_up, integral_sys_dn])
        integral_tot_error = math.sqrt( ( integral_stat_error * integral_stat_error ) + ( max_integral_sys * max_integral_sys ) )

        # This will print the total yield w/ stat and syst error per each sys

        if debug:
            print ("\n\t\tIntegral = {0:.3f} +- {1:.3f} (stat) +- {2:.3f} (syst: {3}) [ +{4:.3f}, -{5:.3f} ]\n".format(integral_nominal, integral_stat_error, integral_simm_sys_unc, sysname, integral_sys_up, integral_sys_dn) )

        # This will sum in quadrature the stat & systematic uncertainties for each bin

        for bin, list_unc in g_unc_bins.iteritems():
            g_sq_unc_bins[bin] = sumQuadrature([x[0] for x in list_unc])
            if debug:
                print("\t\t{0}-th bin,".format(bin) + " list of uncertainties: [" + ",".join( "{0:.3f} ({1},{2})".format(x[0],x[1],x[2]) for x in list_unc ) + "]" + " - tot. uncertainty : {0:.3f}".format(g_sq_unc_bins[bin]) )

        # This will sum in quadrature the systematic uncertainties (excluding the statistical) for each bin

        for bin, list_unc in g_unc_bins_NO_STAT.iteritems():
            g_sq_unc_bins_NO_STAT[bin] = sumQuadrature([x[0] for x in list_unc])
            if debug:
                print("\t\t{0}-th bin,".format(bin) + " list of uncertainties (excluded stat): [" + ",".join( "{0:.3f} ({1},{2})".format(x[0],x[1],x[2]) for x in list_unc ) + "]" + " - tot. uncertainty (excluded stat): {0:.3f}".format(g_sq_unc_bins_NO_STAT[bin]) )
        if debug:
            print("")

    else:
        print ("\t\tIntegral = {0:.3f} +- {1:.3f} (stat)".format(integral_nominal, integral_stat_error ))

    return integral_nominal, integral_stat_error


def sumQuadrature ( inlist ):
    sq = 0
    for elem in inlist:
    	sq +=  pow(elem,2.0)
    return math.sqrt(sq)


def getTotFakeUncertainty( nominal, stat, flav, debug=False ):

    non_closure_vals = {}

    if "HIGHNJ" in args.channel:

        # non_closure_vals["MuMu"] = 0.0461 # Updated on v27
        # non_closure_vals["OF"]   = 0.1316 # Updated on v27
        # non_closure_vals["ElEl"] = 0.1711 # Updated on v27

        non_closure_vals["MuMu"] = 0.003 # Updated on v28
        if args.paramFakeEl == "Pt":
            non_closure_vals["OF"]   = 0.2146 # Updated on v28
            non_closure_vals["ElEl"] = 0.3425 # Updated on v28
        elif args.paramFakeEl == "NBJets_VS_Pt":
            non_closure_vals["OF"]   = 0.1784 # Updated on v28
            non_closure_vals["ElEl"] = 0.2712 # Updated on v28

    elif "LOWNJ" in args.channel:

        # non_closure_vals["MuMu"] = 0.0281 # Updated on v27
        # non_closure_vals["OF"]   = 0.0481 # Updated on v27
        # non_closure_vals["ElEl"] = 0.0208 # Updated on v27

        non_closure_vals["MuMu"] = 0.0289 # Updated on v28
        if args.paramFakeEl == "Pt":
            non_closure_vals["OF"]   = 0.0568 # Updated on v28
            non_closure_vals["ElEl"] = 0.2730 # Updated on v28
        elif args.paramFakeEl == "NBJets_VS_Pt":
            non_closure_vals["OF"]   = 0.0499 # Updated on v28
            non_closure_vals["ElEl"] = 0.2657 # Updated on v28


    if flav == "Inclusive" :
        non_closure ={"MuMu":non_closure_vals["MuMu"] * nominal, "OF":non_closure_vals["OF"] * nominal, "ElEl":non_closure_vals["ElEl"] * nominal}
    else:
        non_closure = non_closure_vals[flav] * nominal

    # If you are doing closure test, do not consider closure syst!

    if args.closure:
        non_closure = 0.0

    g_sys_dict["Non_Closure"]      = non_closure
    g_sysgroup_dict["Non_Closure"] = [non_closure]

    # This prints out all systematic uncertainties, sorting systematics from smallest to largest

    if debug:
        print("\t\tTot. yields w/ systematics (ungrouped):\n")
        print g_sys_dict
        # print ("\t\tIntegral = {0:.2f}\n\t\t+- {1:.2f} [{2:.2f} %] (Sidebands Stat)\n\t\t+-".format(nominal, stat, (stat/nominal)*100) + "\t\t+-".join( " {0:.2f} [{1:.2f} %] ({2}) \n".format( g_sys_dict[key], (g_sys_dict[key]/nominal)*100, key ) for key in sorted( g_sys_dict, key=g_sys_dict.get ) and not type(g_sys_dict[key]) is dict ) )
        print("")

    # Print sum in quadrature of syst for each syst group from smallest to largest

    sq_list = []
    for sg, values in g_sysgroup_dict.iteritems():
        if debug:
            print "\t", sg
            print  "\t", values
        if any( type(t) is dict for t in values ): continue
        tup = ( sg, sumQuadrature(values), (sumQuadrature(values)/nominal)*100 )
        sq_list.append(tup)

    print("\t\tTot. yields w/ systematics (grouped):\n")
    print ("\t\tIntegral = {0:.2f} +- {1:.2f} [{2:.2f} %] (Sidebands Stat)".format(nominal, stat, (stat/nominal)*100) )
    for s in sorted( sq_list, key = lambda sq : sq[2] ):
        print("\t\t+- {0:.2f} [{1:.2f} %] ({2})".format( s[1], s[2], s[0] ) )
    print("")

    # Print sum in quadrature of ALL syst + stat

    toterrlist = list(g_sys_dict.values())
    toterrlist.extend([stat])
    toterrlist = [ e for e in toterrlist if not type(e) is dict ]
    sq = sumQuadrature( toterrlist )

    toterrlist_NO_STAT = list(g_sys_dict.values())
    toterrlist_NO_STAT = [ e for e in toterrlist_NO_STAT if not type(e) is dict ]
    sq_NO_STAT = sumQuadrature( toterrlist_NO_STAT )

    print ("\t\tIntegral = {0:.2f} +- {1:.2f} [{2:.2f} %] (TOTAL UNCERTAINTY)".format(nominal, sq, (sq/nominal)*100))
    print ("\t\t                   +- {0:.2f} [{1:.2f} %] (TOTAL SYST. UNCERTAINTY)".format(sq_NO_STAT, (sq_NO_STAT/nominal)*100))

    # Add in quadrature the non-closure uncertainty to each bin's total uncertainty

    if not flav == "Inclusive":
        print("")
        for bin, list_unc in g_unc_bins.iteritems():
            g_sq_unc_bins[bin] = sumQuadrature( [x[0] for x in list_unc] + [ non_closure_vals[flav] * g_content_bins[bin] ] )
            if debug:
                print("\t\t{0}-th bin,".format(bin) + " list of uncertainties (INCLUDING stat, including non-closure): [" + ",".join( "{0:.3f}".format(x[0]) for x in list_unc ) + ",{0:3f}]".format(non_closure_vals[flav] * g_content_bins[bin]) + " --> tot. uncertainty = {0:.3f}".format(g_sq_unc_bins[bin]) )
        print("")
        for bin, list_unc in g_unc_bins_NO_STAT.iteritems():
            g_sq_unc_bins_NO_STAT[bin] = sumQuadrature( [x[0] for x in list_unc] + [ non_closure_vals[flav] * g_content_bins[bin] ] )
            if debug:
                print("\t\t{0}-th bin,".format(bin) + " list of uncertainties (EXCLUDING stat, including non-closure): [" + ",".join( "{0:.3f}".format(x[0]) for x in list_unc ) + ",{0:3f}]".format(non_closure_vals[flav] * g_content_bins[bin]) + " --> tot. uncertainty = {0:.3f}".format(g_sq_unc_bins_NO_STAT[bin]) )
        print("")
    else:
        if debug:
            print("")
            print "bin[1] = ", g_content_bins[1]
            print "bin[2] = ", g_content_bins[2]
            print "bin[3] = ", g_content_bins[3]
        g_sq_unc_bins[1] = sumQuadrature( [x[0] for x in g_unc_bins[1]] + [ non_closure_vals["MuMu"] * g_content_bins[1] ] )
        g_sq_unc_bins[2] = sumQuadrature( [x[0] for x in g_unc_bins[2]] + [ non_closure_vals["OF"]   * g_content_bins[2] ] )
        g_sq_unc_bins[3] = sumQuadrature( [x[0] for x in g_unc_bins[3]] + [ non_closure_vals["ElEl"] * g_content_bins[3] ] )
        if debug:
            print("1-th bin (MuMu), list of uncertainties (INCLUDING stat, including non-closure): [" + ",".join( "{0:.3f}".format(x[0]) for x in g_unc_bins[1] ) + ",{0:3f}]".format(non_closure_vals["MuMu"] * g_content_bins[1]))
            print("1-th bin (OF), list of uncertainties (INCLUDING stat, including non-closure): ["   + ",".join( "{0:.3f}".format(x[0]) for x in g_unc_bins[2] ) + ",{0:3f}]".format(non_closure_vals["OF"]   * g_content_bins[2]))
            print("1-th bin (ElEl), list of uncertainties (INCLUDING stat, including non-closure): [" + ",".join( "{0:.3f}".format(x[0]) for x in g_unc_bins[3] ) + ",{0:3f}]".format(non_closure_vals["ElEl"] * g_content_bins[3]))
            print("")
        g_sq_unc_bins_NO_STAT[1] = sumQuadrature( [x[0] for x in g_unc_bins_NO_STAT[1]] + [ non_closure_vals["MuMu"] * g_content_bins[1] ] )
        g_sq_unc_bins_NO_STAT[2] = sumQuadrature( [x[0] for x in g_unc_bins_NO_STAT[2]] + [ non_closure_vals["OF"]   * g_content_bins[2] ] )
        g_sq_unc_bins_NO_STAT[3] = sumQuadrature( [x[0] for x in g_unc_bins_NO_STAT[3]] + [ non_closure_vals["ElEl"] * g_content_bins[3] ] )
        if debug:
            print("1-th bin (MuMu), list of uncertainties (EXCLUDING stat, including non-closure): [" + ",".join( "{0:.3f}".format(x[0]) for x in g_unc_bins_NO_STAT[1] ) + ",{0:3f}]".format(non_closure_vals["MuMu"] * g_content_bins[1]))
            print("1-th bin (OF), list of uncertainties (EXCLUDING stat, including non-closure): ["   + ",".join( "{0:.3f}".format(x[0]) for x in g_unc_bins_NO_STAT[2] ) + ",{0:3f}]".format(non_closure_vals["OF"]   * g_content_bins[2]))
            print("1-th bin (ElEl), list of uncertainties (EXCLUDING stat, including non-closure): [" + ",".join( "{0:.3f}".format(x[0]) for x in g_unc_bins_NO_STAT[3] ) + ",{0:3f}]".format(non_closure_vals["ElEl"] * g_content_bins[3]))

    # Print out results in LaTeX friendly format!

    print("")
    print("\\begin{{table}}\n\\begin{{center}}\n\\begin{{tabular}}{{ll}}\n\\toprule\n\\multicolumn{{2}}{{c}}{{MY_FLAVOUR - $2\\leq N_{{jets}} \\leq 3, N_{{b-tags}} \\geq 1$ VR, SR}} \\\\ \n\\midrule\nFakes (MM) = & {0:.2f} $\\pm$ \\\\\n & {1:.2f} $[{2:.2f}\%]$ (Sidebands stat.) $\pm$ \\\\".format(nominal, stat, (stat/nominal)*100))
    for s in sorted( sq_list, key = lambda sq : sq[2] ):
        proc = "?Unknown Process?"
        if s[0] == "Non_Closure" : proc = "Non closure"
        if s[0] == "ND_FakesOS"  : proc = "Fakes OS sub."
        if s[0] == "ND_TTV"      : proc = "$t\\bar{t}W,t\\bar{t}Z$ sub."
        if s[0] == "ND_VV"       : proc = "$WW,WZ,ZZ$ sub."
        if s[0] == "ND_OtherPromptSS" : proc = "Other Prompt SS sub."
        if s[0] == "D_QMisID"   : proc = "QMisID sub., $T\\bar{T}$"
        if s[0] == "N_QMisID"   : proc = "QMisID sub., $TT$"
        if s[0] == "Stat"       : proc = "$\\varepsilon$ stat. unc."
        print(" & {0:.2f} $[{1:.2f}\%]$ ({2}) $\pm$ \\\\".format( s[1], s[2], proc ) )
    print("\\midrule\nFakes (MM) = &{0:.2f} $\pm$ \\\\\n & {1:.2f} $[{2:.2f}\%]$ (Tot. uncertainty) \\\\\n\\bottomrule\n\\end{{tabular}}\n\\end{{center}}\n\\end{{table}}".format(nominal, sq, (sq/nominal)*100))
    print("")

    return sq, sq_NO_STAT


def clearDicts():

    g_sys_dict.clear()
    g_sysgroup_dict .clear()
    g_content_bins.clear()
    g_unc_bins.clear()
    g_unc_bins_NO_STAT.clear()
    g_sq_unc_bins.clear()
    g_sq_unc_bins_NO_STAT.clear()

def saveSystHistogram( flav, var, nominalhist, tot_syst ):

    gROOT.SetBatch(True)

    nbins = nominalhist.GetSize()-2
    lowest_edge  = nominalhist.GetXaxis().GetBinLowEdge(0)
    highest_edge = nominalhist.GetXaxis().GetBinUpEdge(nominalhist.GetNbinsX())

    if nominalhist.GetSize() != len(g_unc_bins_NO_STAT):
        sys.exit("ERROR: nominal hist nbins: {0}, g_unc_bins_NO_STAT size: {1}".format(nbins,len(g_unc_bins_NO_STAT)))

    sysstack = THStack("fakessys_stack_"+flav+"_"+var,"AllSys;"+nominalhist.GetXaxis().GetTitle()+";Sys. [%]")

    legend = TLegend(0.22,0.7,0.45,0.9) # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
    legend.SetBorderSize(0)     # no border
    legend.SetFillStyle(0)      # Legend transparent background
    legend.SetTextSize(0.025)   # Increase entry font size!
    #legend.SetTextFont(42)      # Helvetica

    # By construction, each bin has got the same number of systematic sources contributing.
    # We need to create an histogram for each one of them.
    # Pick the 1st bin by convention, and store the histograms that we're going to stack up.

    systhists = []
    for idx, elem in enumerate( sorted( g_unc_bins_NO_STAT[1], key = lambda elem : elem[0], reverse = True ) ):
        sysvalue = elem[0]
        sysname  = elem[1]
        sysgroup = elem[2]
        #print("\tsysvalue: {0:.3f}, sysname: {1}, sysgroup: {2}".format(sysvalue,sysname,sysgroup))
        hist = TH1D(sysname,sysname,nbins,lowest_edge,highest_edge)
        hist.SetLineColor(kBlack)
        hist.SetFillColor(idx+2)
        hist.SetFillStyle(1001)
        hist.SetDirectory(0)
        legend.AddEntry(hist,"{0}".format(sysname), "F")
        systhists.append(hist)

    # Loop over the histograms just created, and assign to every bin the corresponding systematic value

    for hist in systhists:
        thissystname = hist.GetName()
        #print("\tsysname: {0}:".format(thissystname))
        for bin in range(1,hist.GetSize()):
            thisbin_syslist = g_unc_bins_NO_STAT[bin]
            thistuple = [ tup for tup in thisbin_syslist if tup[1] == thissystname ]
            #print("\t\tthistuple: {0}".format(thistuple))
            if thistuple:
                #print("\t\t{0}-th bin, sysvalue: {1:.3f}".format(bin,thistuple[0][0]))
                if thistuple[0][0] < 0 or nominalhist.Integral(bin,bin) <= 0:
                    hist.SetBinContent(bin,0)
                else:
                    hist.SetBinContent(bin,(thistuple[0][0]/nominalhist.Integral(bin,bin))*100)
            else:
                #print("\t\t{0}-th bin, sysvalue: 0".format(bin))
                hist.SetBinContent(bin,0)

    # Ok, now add the histograms to the stack

    for hist in systhists:
        sysstack.Add(hist)

    # Now draw the stack, and on the same canvas draw the tot. systematic (i.e., the sum in quadrature of all systematics)

    c = TCanvas("cstack","Systematic Stack",50,50,600,600)

    sysstack.Draw("HIST")

    syssumquadhist = TH1D("sumquadsyshist","sumquadsyshist",nbins,lowest_edge,highest_edge)
    syssumquadhist.SetLineWidth(3)
    syssumquadhist.SetLineColor(1)
    syssumquadhist.SetLineStyle(2)
    for bin, sq in g_sq_unc_bins_NO_STAT.iteritems():
        if bin == 0 or nominalhist.Integral(bin,bin) <= 0:
            syssumquadhist.SetBinContent(bin,0)
        else:
            syssumquadhist.SetBinContent(bin,sq/nominalhist.Integral(bin,bin))
    legend.AddEntry(syssumquadhist,"TOT. SYST. (sum. quad.)", "L")
    syssumquadhist.Draw("HIST SAME")

    # Draw also the total systematic for the integral

    syssumquadhist_integral = TH1D("sumquadsyshist_integral","sumquadsyshist_integral",1,lowest_edge,highest_edge)
    syssumquadhist_integral.SetLineWidth(3)
    syssumquadhist_integral.SetLineColor(kOrange+7)
    syssumquadhist_integral.SetLineStyle(2)
    syssumquadhist_integral.SetBinContent(1,tot_syst/nominalhist.Integral(0,nominalhist.GetNbinsX()+1))
    legend.AddEntry(syssumquadhist_integral,"TOT. SYST. (sum. quad.) - Integral", "L")
    syssumquadhist_integral.Draw("HIST SAME")

    legend.Draw()

    outpath = args.inputDir
    if outpath[-1] == '/':
      outpath = outpath[:-1]

    outpath += "/" + flav + "_FakesSysSplit"
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    extensions = [(".png","PNG"),(".pdf","PDF"),(".root","ROOT")]
    for ext in extensions:
        finalpath = outpath+"/"+ext[1]
        if not os.path.exists(finalpath):
            os.makedirs(finalpath)
        outname = finalpath + "/" + flav + "_" + var + "_FakesSysSplit"
        c.SaveAs( outname + ext[0] )

    foutputpath = outpath + "/ROOT_HISTS"
    if not os.path.exists(foutputpath):
        os.makedirs(foutputpath)

    foutput = TFile(foutputpath+"/"+flav+"_"+var+"_SysUncertHists"+".root","RECREATE")
    for hist in systhists:
        hist.Write()
    foutput.Close()


def makeSysPlots( flav, var, observedhist, expectedhist, fakeshist, debug=False ):

    doLogY = args.doLogScaleY

    if debug:
        print("\nmakeSysPlots()\n")

    gROOT.SetBatch(True)

    c = TCanvas("c1","Temp",50,50,600,600)

    pad1 = TPad("pad1", "", 0, 0.25, 1, 1)
    pad2 = TPad("pad2", "", 0, 0,   1, 0.25)
    pad1.SetBottomMargin(0.02)
    pad2.SetBottomMargin(0.4)
    pad1.Draw()
    pad2.Draw()

    # For expected hist, set the bin error as the sum in quadrature of stat (on expected) + syst (on fakes)

    new_expectedhist = expectedhist.Clone(expectedhist.GetName())
    new_expectedhist.SetDirectory(0)

    if debug:
        print("NEW expectedhist:")
    for ibin in range(1,new_expectedhist.GetSize()):
        uncertlist = [ g_sq_unc_bins_NO_STAT[ibin], new_expectedhist.GetBinError(ibin) ]
        new_expectedhist.SetBinError(ibin, sumQuadrature(uncertlist) )
        if debug:
            print("\t\t{0}-th bin,".format(ibin) + " tot. uncertainty (stat+syst): {0:.3f}".format( sumQuadrature(uncertlist) ) )
        if new_expectedhist.Integral(ibin,ibin) > 0:
            if debug:
                print("\t\tbin[{0}] = {1:.3f} +- {2:.3f} ({3:.3f} [%])".format( ibin, new_expectedhist.GetBinContent(ibin), new_expectedhist.GetBinError(ibin) , ( new_expectedhist.GetBinError(ibin) / new_expectedhist.Integral(ibin,ibin) ) * 100 ) )

    new_expectedhist.SetLineWidth(2)
    new_expectedhist.SetLineStyle(1)
    new_expectedhist.SetLineColor(kBlue-4)
    new_expectedhist.SetFillColor(kWhite)
    new_expectedhist.SetFillStyle(1001)

    # For fakes hist, set the bin error as the sum in quadrature of stat + syst

    new_fakeshist = fakeshist.Clone(fakeshist.GetName())
    new_fakeshist.SetDirectory(0)

    if debug:
        print("NEW fakeshist:")
    for ibin in range(1,new_fakeshist.GetSize()):
        uncertlist = [ g_sq_unc_bins[ibin] ]
        new_fakeshist.SetBinError(ibin, sumQuadrature(uncertlist) )
        if debug:
            print("\t\t{0}-th bin,".format(ibin) + " tot. uncertainty (stat+syst): {0:.3f}".format( sumQuadrature(uncertlist) ) )
        if new_fakeshist.Integral(ibin,ibin) > 0:
            if debug:
                print("\t\tbin[{0}] = {1:.3f} +- {2:.3f} ({3:.3f} [%])".format( ibin, new_fakeshist.GetBinContent(ibin), new_fakeshist.GetBinError(ibin) , ( new_fakeshist.GetBinError(ibin) / new_fakeshist.Integral(ibin,ibin) ) * 100 ) )

    new_fakeshist.SetLineWidth(2)
    new_fakeshist.SetLineStyle(1)
    new_fakeshist.SetLineColor(kViolet-4)
    new_fakeshist.SetFillColor(kWhite)
    new_fakeshist.SetFillStyle(1001)

    if observedhist:
        observedgr = makePoissonErrors(observedhist)
        observedgr.SetMarkerSize(0.8) # (1.2)
        observedgr.SetLineColor(1)
        observedgr.SetLineWidth(2)
        observedgr.SetMarkerStyle(20)
        observedgr.SetLineStyle(1)

    err = new_expectedhist.Clone("tot_uncertainty")
    err.SetFillColor(kOrange)
    err.SetLineColor(10)
    err.SetFillStyle(3356)
    gStyle.SetHatchesLineWidth(2)
    gStyle.SetHatchesSpacing(0.8)
    err.SetMarkerSize(0)

    # Trick to rescale:

    ymax      = err.GetMaximum()
    binmax    = err.GetMaximumBin()
    binmaxerr = err.GetBinError(binmax)

    if new_expectedhist.GetMaximum() > ymax:
        ymax      = new_expectedhist.GetMaximum()
        binmax    = new_expectedhist.GetMaximumBin()
        binmaxerr = new_expectedhist.GetBinError(binmax)

    #print("FLAV: {0} - MM max = {1:.2f} - expected max = {2:.2f}".format(flav,new_expectedhist.GetMaximum(),observedhist.GetMaximum()))

    if observedhist:
        if observedhist.GetMaximum() > ymax:
            ymax      = observedhist.GetMaximum()
            binmax    = observedhist.GetMaximumBin()
            binmaxerr = observedhist.GetBinError(binmax)

    #print("FLAV: {0} - max = {1:.2f} - binmax = {2} - binmaxerr = {3:.2f}".format(flav,ymax,binmax,binmaxerr))

    maxfactor = 1.3
    minimum    = 0
    if doLogY:
        maxfactor = 1e2 * 1.3
        minimum = 0.01
    err.SetMaximum( ( ymax + binmaxerr ) * maxfactor )
    err.SetMinimum(minimum)

    # --------------------------

    # Stat + sys error on expected (ratio)

    # ratio_err = new_expectedhist.Clone("RatioErr")
    # ratio_err.SetXTitle(new_expectedhist.GetXaxis().GetTitle())
    # ratio_err.SetYTitle("Data/Exp.")
    # ratio_err.GetXaxis().SetTitleSize(0.15)
    # ratio_err.GetYaxis().SetTitleSize(0.15)
    # ratio_err.GetXaxis().SetTitleOffset(0.90)
    # ratio_err.GetYaxis().SetTitleOffset(0.35)
    # ratio_err.GetXaxis().SetLabelSize(0.15)
    # ratio_err.GetYaxis().SetLabelSize(0.12)
    # ratio_err.GetYaxis().SetNdivisions(505) #(5)
    # ratio_err.SetFillColor(kOrange)
    # ratio_err.SetLineColor(kOrange)#(10)
    # ratio_err.SetFillStyle(3356)
    # gStyle.SetHatchesLineWidth(2)
    # gStyle.SetHatchesSpacing(0.8)
    # ratio_err.SetMarkerSize(0)

    # # print("")
    # # for ibin in range(1,ratio_err.GetSize()):
    # #     print("\t\tratio_err[{0}]  = {1:.3f} +- {2:.4f}".format( ibin, ratio_err.GetBinContent(ibin), ratio_err.GetBinError(ibin) ))

    # ratio_err.Divide(new_expectedhist)
    # for ibin in range(1,ratio_err.GetSize()):
    #     if new_expectedhist.GetBinContent(ibin) > 0:
    #         ratio_err.SetBinError(ibin, new_expectedhist.GetBinError(ibin)/new_expectedhist.GetBinContent(ibin))

    # # for ibin in range(1,ratio_err.GetSize()):
    # #     print("\t\tratio_err[{0}]  = {1:.3f} +- {2:.4f}".format( ibin, ratio_err.GetBinContent(ibin), ratio_err.GetBinError(ibin) ))

    ratio_err_props = {
        "TitleX":new_expectedhist.GetXaxis().GetTitle(),
        "TitleY":"Data/Exp.",
        "TitleSizeX":0.15,
        "TitleSizeY":0.15,
        "TitleOffsetX":0.90,
        "TitleOffsetY":0.35,
        "LabelSizeX":0.15,
        "LabelSizeY": 0.12,
        "NdivisionsY":505,
        "FillColor":kOrange,
        "FillStyle":3356,
        "LineColor":kOrange,
        "MarkerSize":0,
    }

    ratio_err = SelfDivide( new_expectedhist, "RatioErr", ratio_err_props )

    # obs / exp

    ymax_ratio = -999
    ymin_ratio = 999

    if observedhist:
        ratio_obs_exp = observedhist.Clone("RatioObsExp")
        ratio_obs_exp.SetYTitle("Data/Exp.")
        ratio_obs_exp.SetLineStyle(1)
        ratio_obs_exp.SetMarkerSize(0.8)
        ratio_obs_exp.SetLineColor(1)
        ratio_obs_exp.SetMarkerStyle(20)
        ratio_obs_exp.SetLineWidth(1)
        ratio_obs_exp.Divide(new_expectedhist)

    # --------------------------

    pad1.cd()

    if doLogY:
        gPad.Update()
        gPad.SetLogy()

    err.GetXaxis().SetLabelSize(0)
    err.GetXaxis().SetLabelOffset(999)

    err.Draw("E2")
    new_expectedhist.Draw("HIST SAME")
    if observedhist:
        observedgr.Draw("P SAME")

    legend = TLegend(0.55,0.7,0.75,0.9) # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
    legend.SetBorderSize(0)	 # no border
    legend.SetFillStyle(0)	 # Legend transparent background
    legend.SetTextSize(0.035)    # Increase entry font size!
    #legend.SetTextFont(42)	 # Helvetica

    legend.AddEntry(new_expectedhist, "Tot. Expected ({0:.1f})".format(integrate(new_expectedhist)), "F")
    if observedhist:
        legend.AddEntry(observedhist, "Data ({0:.0f})".format(integrate(observedhist)), "P")
    legend.AddEntry(err, "Stat. + Sys. Unc. (Fakes)", "F")

    leg_ATLAS = TLatex()
    leg_lumi  = TLatex()
    leg_ATLAS.SetTextSize(0.03)
    leg_ATLAS.SetNDC()
    leg_lumi.SetTextSize(0.03)
    leg_lumi.SetNDC()

    legend.Draw()
    leg_ATLAS.DrawLatex(0.23,0.85,"#bf{#it{ATLAS}} Work In Progress")
    leg_lumi.DrawLatex(0.23,0.75,"#sqrt{{s}} = 13 TeV, #int L dt = {0:.1f} fb^{{-1}}".format(args.lumi))

    pad2.cd()

    # Hardcoded axis limits
    ratio_err.GetYaxis().SetRangeUser(0.0, 2.0)
    pad2.SetGridy(1)

    ratio_err.Draw("E2")
    if observedhist:
        ratio_obs_exp.Draw("PE SAME")

    refl = TLine(ratio_err.GetBinLowEdge(1), 1.0, ratio_err.GetBinLowEdge(ratio_err.GetNbinsX()+1), 1.0)
    refl.SetLineStyle(2)
    refl.SetLineColor(kRed)
    refl.SetLineWidth(2)
    refl.Draw("SAME")

    # --------------------------

    outpath = args.inputDir
    if outpath[-1] == '/':
      outpath = outpath[:-1]

    outpath += "/" + flav + "_DataVSExp_FakesSys"
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    extensions = [(".png","PNG"),(".pdf","PDF"),(".root","ROOT")]
    for ext in extensions:
        finalpath = outpath+"/"+ext[1]
        if not os.path.exists(finalpath):
            os.makedirs(finalpath)
        outname = finalpath + "/" + flav + "_" + var + "_DataVSExp_FakesSys"
        if doLogY:
            outname += "_LOGY"
        c.SaveAs( outname + ext[0] )

    return new_expectedhist, new_fakeshist


def makeSysPlotsClosure( flav, var, MC_hist, MM_hist ):

    gROOT.SetBatch(True)

    c = TCanvas("c1","Temp",50,50,600,600)

    pad1 = TPad("pad1", "", 0, 0.25, 1, 1)
    pad2 = TPad("pad2", "", 0, 0,   1, 0.25)
    pad1.SetBottomMargin(0.02)
    pad2.SetBottomMargin(0.4)
    pad1.Draw()
    pad2.Draw()

    new_MM_hist = MM_hist.Clone(MM_hist.GetName())

    # For MM closure hist, set the bin error as the sum in quadrature of stat+syst
    # If looking at the SR, take only the part of error which is not correlated to the MC stats error
    # That is, take only the syst error from the Fake CR size

    for ibin in range(1,new_MM_hist.GetNbinsX()+1):
        if "HIGHNJ" in args.channel:
            print("\nbin: {0}\tstat error: {1:.2f} - sys error: {2:.2f} - stat+sys error: {3:.2f}".format(ibin,new_MM_hist.GetBinError(ibin),g_sq_unc_bins_NO_STAT[ibin],g_sq_unc_bins[ibin]))
            new_MM_hist.SetBinError(ibin, g_sq_unc_bins_NO_STAT[ibin] )
        else:
            new_MM_hist.SetBinError(ibin, g_sq_unc_bins[ibin] )

    MC_hist.SetLineStyle(2)
    MC_hist.SetLineWidth(3)
    MC_hist.SetLineColor(kViolet-4)
    MC_hist.SetMarkerSize(0.8)
    MC_hist.SetMarkerColor(kViolet-4) # (1)
    MC_hist.SetMarkerStyle(20)

    new_MM_hist.SetLineWidth(2)
    new_MM_hist.SetLineStyle(1)
    new_MM_hist.SetLineColor(kViolet-4)
    new_MM_hist.SetFillColor(kWhite)
    new_MM_hist.SetFillStyle(1001)

    gStyle.SetHatchesLineWidth(2)
    gStyle.SetHatchesSpacing(0.8)

    err = new_MM_hist.Clone("tot_uncertainty")
    err.SetFillColor(kGray)
    err.SetLineColor(kGray)
    err.SetFillStyle(3356)
    err.SetMarkerSize(0)

    # Trick to rescale histograms in first pad

    ymax      = new_MM_hist.GetMaximum()
    binmax    = new_MM_hist.GetMaximumBin()
    binmaxerr = new_MM_hist.GetBinError(binmax)

    #print("FLAV: {0} - MM max = {1:.2f} - ttbar max = {2:.2f}".format(flav,new_MM_hist.GetMaximum(),MC_hist.GetMaximum()))

    if MC_hist.GetMaximum() > ymax:
        ymax      = MC_hist.GetMaximum()
        binmax    = MC_hist.GetMaximumBin()
        binmaxerr = MC_hist.GetBinError(binmax)

    #print("FLAV: {0} - max = {1:.2f} - binmax = {2} - binmaxerr = {3:.2f}".format(flav,ymax,binmax,binmaxerr))

    # Set this to the first histogram that will be plotted in the pad

    err.SetMaximum( (ymax + binmaxerr) *1.3 )
    err.SetMinimum(0)

    # --------------------------

    # Make the actual non closure plot

    ratio_MC_MM = new_MM_hist.Clone("RatioMCMM")
    ratio_MC_MM = ratio_MC_MM - MC_hist
    ratio_MC_MM.Divide(MC_hist) # Do not care about the errors: this histogram will be drawn w/o error bars
    ratio_MC_MM.SetYTitle("Non-closure")
    ratio_MC_MM.SetLineStyle(1)
    ratio_MC_MM.SetLineWidth(2)
    ratio_MC_MM.SetLineColor(kRed)
    ratio_MC_MM.SetMarkerSize(0.8)
    ratio_MC_MM.SetMarkerColor(1)
    ratio_MC_MM.SetMarkerStyle(20) # (22)

    ratio_MC_MM_err = ratio_MC_MM.Clone("RatioErr")

    # Calculate the error on the non-closure using the standard error propagation bin by bin
    # On the MM estimate, the error is mostly from the Fake CR size, which is orthogonal to all the MM sidebands.
    # Therefore we neglect correlations in the errror computation, which should be accounted for b/c the MM prediction
    # in the closure test actually contains the events of the pure MC prediction.

    # Dummy histogram to show the non-closure
    effective_NC = ratio_MC_MM.Clone("Effective")
    effective_NC.SetLineColor(kBlack)
    # effective_NC.SetLineStyle(2)

    for ibin in range(1,ratio_MC_MM_err.GetSize()):

        iMM    = new_MM_hist.GetBinContent(ibin)
        iMMerr = new_MM_hist.GetBinError(ibin)
        iMC    = MC_hist.GetBinContent(ibin)
        iMCerr = MC_hist.GetBinError(ibin)

        print("\nbin: {0}\tiMM = {1:.2f} - iMMerr = {2:.2f}\tiMC = {3:.2f} - iMCerr = {4:.2f}".format(ibin,iMM,iMMerr,iMC,iMCerr))
        iclosure_err = 0.0
        if iMC and iMCerr:
            iclosure_err = math.sqrt( (iMMerr*iMMerr) / (iMC*iMC) + (iMM*iMM) * (iMCerr*iMCerr) / (iMC*iMC*iMC*iMC) )

        ratio_MC_MM_err.SetBinError(ibin, iclosure_err)
        print("NON CLOSURE: {0:.2f} +- {1:.2f} [%]".format(ratio_MC_MM_err.GetBinContent(ibin)*1e2,iclosure_err*1e2))

        # Calculate the effective non-closure, and store it as errro for the dummy histogram

        effectiveNC = abs(ratio_MC_MM_err.GetBinContent(ibin)) - iclosure_err
        if effectiveNC <= 0: effectiveNC = 0
        print("EFFECTIVE NON CLOSURE: {0:.2f} [%]".format(effectiveNC))
        effective_NC.SetBinContent(ibin,effectiveNC)

    ratio_MC_MM_err.SetYTitle("#frac{MM-t#bar{t}}{t#bar{t}}")
    ratio_MC_MM_err.GetXaxis().SetTitleSize(0.15)
    ratio_MC_MM_err.GetYaxis().SetTitleSize(0.12)
    ratio_MC_MM_err.GetXaxis().SetTitleOffset(0.90)
    ratio_MC_MM_err.GetYaxis().SetTitleOffset(0.45)
    ratio_MC_MM_err.GetXaxis().SetLabelSize(0.15)
    ratio_MC_MM_err.GetYaxis().SetLabelSize(0.12)
    ratio_MC_MM_err.GetYaxis().SetNdivisions(505) #(5)
    ratio_MC_MM_err.SetFillColor(kOrange)
    ratio_MC_MM_err.SetLineColor(kOrange)
    ratio_MC_MM_err.SetFillStyle(3356)
    ratio_MC_MM_err.SetMarkerSize(0)

    # # Trick to rescale ratio pad

    # ymax_ratio = ratio_MC_MM.GetMaximum()
    # ymin_ratio = ratio_MC_MM.GetMinimum()

    # #print("FLAV: {0} - MC/MM max = {1:.2f} - MC/MM min = {2:.2f}".format(flav,ymax_ratio,ymin_ratio))

    # # For ratio error, get the bin with maximum uncertainty and the uncertainty value

    # max_unc = 0
    # for ibin in range(1, ratio_err.GetNbinsX()+1):
    #     this_unc = ratio_err.GetBinError(ibin)
    #   	if this_unc > max_unc:
    #         max_unc = this_unc

    # #print("FLAV: {0} - max unc./2 ratio = {1:.2f}".format(flav,max_unc/2.0))

    # if ( 1 + max_unc/2.0 ) > ymax_ratio:
    # 	ymax_ratio = 1 + max_unc/2.0
    # if ( 1 - max_unc/2.0 ) < ymin_ratio:
    # 	ymin_ratio = 1 - max_unc/2.0

    # if ymin_ratio < 0: ymin_ratio = 0

    # #print("FLAV: {0} - actual max = {1:.2f} - actual min = {2:.2f}".format(flav,ymax_ratio,ymin_ratio))

    # ratio_err.SetMaximum( ymax_ratio * 1.2 )
    # ratio_err.SetMinimum( ymin_ratio * 0.8 )

    # --------------------------

    pad1.cd()

    err.GetXaxis().SetLabelSize(0)
    err.GetXaxis().SetLabelOffset(999)

    err.Draw("E2")
    new_MM_hist.Draw("HIST SAME")
    MC_hist.Draw("PE SAME")

    legend = TLegend(0.55,0.7,0.75,0.9) # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
    legend.SetBorderSize(0)	# no border
    legend.SetFillStyle(0)	# Legend transparent background
    legend.SetTextSize(0.03)	# Increase entry font size!
    #legend.SetTextFont(42)	# Helvetica

    legend.AddEntry(new_MM_hist,"Fakes MM (t#bar{{t}}, t#bar{{t}}#gamma inputs) ({0:.1f})".format(integrate(new_MM_hist)),"F")
    legend.AddEntry(MC_hist, "t#bar{{t}}, t#bar{{t}}#gamma ({0:.1f})".format(integrate(MC_hist)), "P")
    err_type = "Stat.+Sys. Unc." if not  "HIGHNJ" in args.channel else "Sys. Unc."
    legend.AddEntry(err, "{0}".format(err_type), "F")

    leg_ATLAS = TLatex()
    leg_lumi  = TLatex()
    leg_ATLAS.SetTextSize(0.03)
    leg_ATLAS.SetNDC()
    leg_lumi.SetTextSize(0.03)
    leg_lumi.SetNDC()

    legend.Draw()
    leg_ATLAS.DrawLatex(0.19,0.85,"#bf{#it{ATLAS}} Work In Progress")
    leg_lumi.DrawLatex(0.19,0.75,"#sqrt{{s}} = 13 TeV, #int L dt = {0:.1f} fb^{{-1}}".format(args.lumi))

    pad2.cd()

    ratio_MC_MM_err.GetYaxis().SetRangeUser(-1.0, 1.0)
    ratio_MC_MM_err.Draw("E2")

    # Write the actual non-closure per bin on the histogram in selected cases

    if any( v == var for v in ["Integral","BDTGScore","NJets2j3j","NJets4j","NBJets"] ):
        gStyle.SetPaintTextFormat(".2f")
        ratio_MC_MM.SetMarkerSize(4.1)
        ratio_MC_MM.SetMarkerColor(kRed)
        ratio_MC_MM.Draw("HIST SAME TEXT0")
        if "HIGHNJ" in args.channel:
            effective_NC.SetMarkerSize(4.1)
            effective_NC.SetMarkerColor(kBlack)
            effective_NC.Draw("HIST SAME TEXT0")
    else:
        ratio_MC_MM.Draw("HIST SAME")

    refl = TLine(ratio_MC_MM.GetBinLowEdge(1), 0.0, ratio_MC_MM.GetBinLowEdge(ratio_MC_MM.GetNbinsX()+1), 0.0)
    refl.SetLineStyle(2)
    refl.SetLineColor(kBlack)
    refl.SetLineWidth(2)
    refl.Draw("SAME")

    # --------------------------

    outpath = args.inputDir
    if outpath[-1] == '/':
      outpath = outpath[:-1]

    outpath += "/" + flav + "_MM_vs_MC_Sys"
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    extensions = [(".png","PNG"),(".pdf","PDF"),(".root","ROOT")]
    for ext in extensions:
        finalpath = outpath+"/"+ext[1]
        if not os.path.exists(finalpath):
            os.makedirs(finalpath)
        outname = finalpath + "/" + flav + "_" + var + "_MM_vs_MC_Sys"
        c.SaveAs( outname + ext[0] )

# -------------------------------------------------------------------------------------------------------------------

def compareFakesKinematics(flav, var, MM_hist, infilename):

    doLogY = args.doLogScaleY

    doMCScaleNorm = True

    c = TCanvas("c1_x","Temp_x",50,50,600,600)

    pad1 = TPad("pad1_x", "", 0, 0.25, 1, 1)
    pad2 = TPad("pad2_x", "", 0, 0,   1, 0.25)
    pad1.SetBottomMargin(0.02)
    pad2.SetBottomMargin(0.4)
    pad1.Draw()
    pad2.Draw()

    filename_MC = infilename.replace("_MM_","_MC_")
    filename_MC = filename_MC.replace("_DataDriven","")

    file_MC = TFile(filename_MC)
    if not file_MC:
        sys.exit("ERROR: file:\n {0}\ndoes not exist!".format(filename_MC))

    MC_hist = file_MC.Get("fakesbkg")
    MC_hist.SetDirectory(0)

    MC_hist.SetLineStyle(2)
    MC_hist.SetLineWidth(3)
    MC_hist.SetLineColor(kViolet-4)
    MC_hist.SetMarkerSize(0.8)
    MC_hist.SetMarkerColor(1)
    MC_hist.SetMarkerStyle(20)
    # MC_hist.SetMarkerSize(0.8)
    # MC_hist.SetMarkerColor(kViolet-4)
    # MC_hist.SetMarkerStyle(24)

    # Fix MC normalisation to the DD norm

    norm_MC = integrate(MM_hist)/integrate(MC_hist)
    if doMCScaleNorm:
        MC_hist.Scale(norm_MC)

    err = MM_hist.Clone("tot_fakesMM_uncertainty")
    err.SetFillColor(kOrange)
    err.SetLineColor(10)
    err.SetFillStyle(3356)
    gStyle.SetHatchesLineWidth(2)
    gStyle.SetHatchesSpacing(0.8)
    err.SetMarkerSize(0)

    maxfactor = 1.3
    minimum = 0
    if doLogY:
        maxfactor = 1e2 * 1.3
        minimum   = 0.01

    err.SetMaximum( ( err.GetMaximum() + err.GetBinError(err.GetMaximumBin()) ) * maxfactor )
    err.SetMinimum(minimum)

    ratio_err = MM_hist.Clone("RatioErr")
    ratio_err.SetXTitle(MM_hist.GetXaxis().GetTitle())
    ratio_err.SetYTitle("MC/MM")
    ratio_err.GetXaxis().SetTitleSize(0.15)
    ratio_err.GetYaxis().SetTitleSize(0.15)
    ratio_err.GetXaxis().SetTitleOffset(0.90)
    ratio_err.GetYaxis().SetTitleOffset(0.35)
    ratio_err.GetXaxis().SetLabelSize(0.15)
    ratio_err.GetYaxis().SetLabelSize(0.12)
    ratio_err.GetYaxis().SetNdivisions(505) #(5)
    ratio_err.SetFillColor(kOrange)
    ratio_err.SetLineColor(kOrange)#(10)
    ratio_err.SetFillStyle(3356)
    gStyle.SetHatchesLineWidth(2)
    gStyle.SetHatchesSpacing(0.8)
    ratio_err.SetMarkerSize(0)

    ratio_err.Divide(MM_hist)
    for ibin in range(1,ratio_err.GetSize()):
        if MM_hist.GetBinContent(ibin) > 0:
            ratio_err.SetBinError(ibin, MM_hist.GetBinError(ibin)/MM_hist.GetBinContent(ibin))

    ratio_MC_MM = MC_hist.Clone("RatioMCMM")
    ratio_MC_MM.SetYTitle("MC/MM")
    ratio_MC_MM.SetLineStyle(2)
    ratio_MC_MM.SetLineWidth(3)
    ratio_MC_MM.SetLineColor(kViolet-4)
    ratio_MC_MM.SetMarkerSize(0.8)
    ratio_MC_MM.SetMarkerColor(1)
    ratio_MC_MM.SetMarkerStyle(20)
    # ratio_MC_MM.SetMarkerSize(0.8)
    # ratio_MC_MM.SetMarkerColor(kViolet-4)
    # ratio_MC_MM.SetMarkerStyle(24)

    ratio_MC_MM.Divide(MM_hist)

    # --------------------------

    pad1.cd()

    if doLogY:
        pad1.Update()
        pad1.SetLogy()

    err.GetXaxis().SetLabelSize(0)
    err.GetXaxis().SetLabelOffset(999)

    err.Draw("E2")
    MM_hist.Draw("HIST SAME")
    MC_hist.Draw("PE SAME")

    legend = TLegend(0.55,0.7,0.75,0.9) # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
    legend.SetBorderSize(0)	 # no border
    legend.SetFillStyle(0)	 # Legend transparent background
    legend.SetTextSize(0.035)    # Increase entry font size!
    #legend.SetTextFont(42)	 # Helvetica

    legend.AddEntry(MM_hist, "Fakes MM ({0:.2f})".format(integrate(MM_hist)), "F")
    if doMCScaleNorm:
        legend.AddEntry(MC_hist, "Fakes MC #times {0:.2f} ({1:.2f})".format(norm_MC,integrate(MC_hist)), "P")
    else:
        legend.AddEntry(MC_hist, "Fakes MC ({0:.2f})".format(integrate(MC_hist)), "P")
    legend.AddEntry(err, "Stat. + Sys. Unc. (Fakes MM)", "F")

    leg_ATLAS = TLatex()
    leg_lumi  = TLatex()
    leg_ATLAS.SetTextSize(0.03)
    leg_ATLAS.SetNDC()
    leg_lumi.SetTextSize(0.03)
    leg_lumi.SetNDC()

    legend.Draw()
    leg_ATLAS.DrawLatex(0.23,0.85,"#bf{#it{ATLAS}} Work In Progress")
    leg_lumi.DrawLatex(0.23,0.75,"#sqrt{{s}} = 13 TeV, #int L dt = {0:.1f} fb^{{-1}}".format(args.lumi))

    pad2.cd()

    # Hardcoded axis limits
    ratio_err.GetYaxis().SetRangeUser(0.0, 2.0)
    pad2.SetGridy(1)

    ratio_err.Draw("E2")
    ratio_MC_MM.Draw("PE SAME")

    refl = TLine(ratio_err.GetBinLowEdge(1), 1.0, ratio_err.GetBinLowEdge(ratio_err.GetNbinsX()+1), 1.0)
    refl.SetLineStyle(2)
    refl.SetLineColor(kRed)
    refl.SetLineWidth(2)
    refl.Draw("SAME")

    c.Update()

    # --------------------------

    outpath = args.inputDir
    if outpath[-1] == '/':
      outpath = outpath[:-1]

    outpath += "/" + flav + "_MC_vs_MM_CompareFakesKinematics"
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    extensions = [(".png","PNG"),(".pdf","PDF"),(".root","ROOT")]
    for ext in extensions:
        finalpath = outpath+"/"+ext[1]
        if not os.path.exists(finalpath):
            os.makedirs(finalpath)
        outname = finalpath + "/" + flav + "_" + var + "_MC_vs_MM_CompareFakesKinematics"
        if doLogY:
            outname += "_LOGY"
        c.SaveAs( outname + ext[0] )

# -------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    TH1.SetDefaultSumw2()

    region = var_name = None

    inputpath = args.inputDir

    if not inputpath.endswith('/'):
         inputpath += '/'

    # var_list = []
    var_list = [
        "Integral",
        "deltaRLep0Lep1",
        "deltaPhiLep0Lep1",
        "MET_FinalTrk",
        "Mll01_inc",
        # "TotLepCharge",
        "NBJets",
        "NJets2j3j",
        "NJets4j",
        #
        "BDTGScore",
        "Lep0Pt",
        "Lep1Pt",
        "Lep0Eta",
        "Lep1Eta",
        "Lep0EtaBE2",
        "Lep1EtaBE2",
        #
        # "Mu0Pt",
        # "Mu1Pt",
        # "Mu0Eta",
        # "Mu1Eta",
        # "Mu0DeltaRClosestJet",
        # "Mu1DeltaRClosestJet",
        # "El0Pt",
        # "El1Pt",
        # "El0Eta",
        # "El1Eta",
        # "El0DeltaRClosestJet",
        # "El1DeltaRClosestJet",
        #
        # "NN_Rebinned",
        # "RNN_Rebinned",
        # "NN_ttV",
        # "NN_top",
        # "RNN_ttV",
        # "RNN_top",
        # "NNComb",
        # "RNNComb",
        ]

    if args.variables:
        var_list = args.variables

    if not args.closure:

        # -----------
        # DATA-DRIVEN
        #------------

        if   "HIGHNJ" in args.channel: region = "SS_SR_DataDriven"
        elif "LOWNJ" in args.channel:  region = "SS_LowNJetCR_DataDriven"

    else:

        # -------------
        # TTBAR CLOSURE
        #--------------

        if "LOWNJ" in args.channel:
            region    = "SS_SR_LowJet_DataDriven_Closure"
        elif "HIGHNJ" in args.channel:
            region    = "SS_SR_HighJet_DataDriven_Closure"
        elif "ALLNJ" in args.channel:
            region    = "SS_SR_AllJet_DataDriven_Closure"

    flavour_list = []
    if "ALL" in args.category:
        flavour_list.extend(["ElEl", "MuMu", "OF"])
    else:
        flavour_list.extend(args.category)

    print("Looking at variables : [" + ",".join( "{0}".format(v) for v in var_list ) + "]")

    for var_name in var_list:

	print ("\n\tVariable: {0}\n".format(var_name))

        if var_name == "NJets2j3j" and "HIGHNJ" in args.channel: continue
        if var_name == "NJets4j" and "LOWNJ" in args.channel: continue

    	for flav in flavour_list:

    	    clearDicts()

            print("\n\t**********************************************************************\n")
    	    print ("\tFlavour region: {0}\n".format(flav))

            if flav == "ElEl" and "Mu" in var_name: continue
            if flav == "MuMu" and "El" in var_name: continue
            if flav == "OF" and any ( f in var_name for f in ["El1","Mu1"] ): continue
            if "HIGHNJ" in args.channel and "NJets2j3j" in var_name : continue
            if "LOWNJ" in args.channel and "NJets4j" in var_name : continue

            # if flav == "Inclusive" and var_name != "LepFlavours": continue
            if var_name == "LepFlavours" and flav != "Inclusive": continue

    	    filename = inputpath + flav + region + "/" + flav + region + "_" + var_name + ".root"

    	    myfile = TFile(filename)
            if not myfile:
                os.sys.exit("ERROR! file:\n{0}\ndoes not exist!".format(filename))

    	    print("\tLooking at file: {0}".format(filename))

    	    fakes_nominal = myfile.Get("fakesbkg")
    	    fakes_syst = {}
    	    for key in myfile.GetListOfKeys():

    		keyname = key.GetName()
    		if not ( "fakesbkg_" in keyname ): continue

                #if not any( k in keyname for k in ["TTV"] ): continue

    		keyname = keyname.replace("fakesbkg_","")
    		keyname = keyname.replace("_dn","")
    		keyname = keyname.replace("_up","")

    		if "Stat" in keyname:
    		    value = "Stat"
    		if "N_QMisID" in keyname:
    		    value = "N_QMisID"
    		if "D_QMisID" in keyname:
    		    value = "D_QMisID"
    		if "ND_TTV" in keyname:
    		    value = "ND_TTV"
    		if "ND_VV" in keyname:
    		    value = "ND_VV"
    		if "ND_OtherPromptSS" in keyname:
    		    value = "ND_OtherPromptSS"
    		if "ND_FakesOS" in keyname:
    		    value = "ND_FakesOS"

    		if not fakes_syst.get(keyname):
    		    fakes_syst[keyname] = value

    	    print ("\n\tFakes:\n")

    	    for sys, sysgroup in fakes_syst.iteritems():
                fakes_up = myfile.Get( "fakesbkg_" + sys + "_up")
                fakes_dn = myfile.Get( "fakesbkg_" + sys + "_dn")
                if args.debug:
                    print("\tsys: {0}, sysgroup: {1}\n".format(sys, sysgroup))
                fakes, fakes_stat_err = getYields(fakes_nominal,fakes_up,fakes_dn, sys, sysgroup, debug=args.debug)

    	    fakes_tot_err, fakes_tot_err_NO_STAT = getTotFakeUncertainty( fakes, fakes_stat_err, flav, debug=args.debug )

            saveSystHistogram(flav,var_name,fakes_nominal,fakes_tot_err_NO_STAT)

    	    expected_nominal = myfile.Get("expectedbkg")
    	    print ("\n\tExpected: \n")
    	    exp, exp_err = getYields(expected_nominal)

    	    if args.closure:

    		ttbar_nominal = myfile.Get("ttbarbkg")
    		print ("\n\tTTbar: \n")
    		ttbar, ttbar_err = getYields(ttbar_nominal)

    		closure     = ( (fakes - ttbar) / (ttbar) ) * 100
    		closure_err = math.sqrt( ( ( math.pow(fakes,2.0) / math.pow(ttbar,4.0) ) * math.pow(ttbar_err,2.0) ) + ( math.pow(fakes_tot_err,2.0) / math.pow(ttbar,2.0) ) ) * 100

    		makeSysPlotsClosure(flav,var_name,ttbar_nominal,fakes_nominal)

                print("\n\t======================================================================\n")
                print("\tCategory : {0} - Channel : {1}\n".format(flav,args.channel[0]))
    		print("\tNon-closure ((fakes-ttbar)/ttbar) = {0:.2f} [%] +- {1:.2f} [%]".format(closure,closure_err))
                if "HIGHNJ" in args.channel:
                    closure_err_uncorr = math.sqrt( ( ( math.pow(fakes,2.0) / math.pow(ttbar,4.0) ) * math.pow(ttbar_err,2.0) ) + ( math.pow(fakes_tot_err_NO_STAT,2.0) / math.pow(ttbar,2.0) ) ) * 100
                    effective_closure = abs(closure) - closure_err_uncorr
                    if effective_closure <= 0: effective_closure = 0
                    print("\tEFFECTIVE Non-closure ((fakes-ttbar)/ttbar) = {0:.2f} [%]".format(effective_closure))
                print("\n\t======================================================================\n")

            else:

    		chargemisid = myfile.Get("qmisidbkg")
    		if chargemisid:
    		    print ("\n\tQMisID: \n")
    		    qmisid, qmisid_err = getYields(chargemisid)

    		observed = myfile.Get("observed")
    		if observed:
    		    print ("\n\tObserved: \n")
    		    obs, obs_err = getYields(observed)

    		signal = myfile.Get("signal")
    		print ("\n\tSignal: \n")
    		sig, sig_err = getYields(signal)

                # Make Data/Expected plots w/ systematic bars for fakes included, and return expected, fakesMM histograms w/ updated errors

                expected_hist_syst, fakes_hist_syst = makeSysPlots(flav,var_name,observed,expected_nominal,fakes_nominal, debug=args.debug)

                if args.doKinematicsComparison:
                    # Compare distributions for DD fakes VS. MC fakes (scaled norm. to DD)
                    compareFakesKinematics(flav,var_name,fakes_hist_syst,filename)

