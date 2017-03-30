 #!/usr/bin/python

import os, sys, math, argparse

sys.path.append(os.path.abspath(os.path.curdir))

parser = argparse.ArgumentParser(description='Get yields and systematics for MM')

list_channel = ['HIGHNJ','LOWNJ','ALLNJ']

categories   = ["ALL","ElEl","MuMu","OF"]

luminosities = { "Moriond 2016 GRL":3.209,            # March 2016
                 "ICHEP 2015+2016 DS":13.20768,       # August 2016
                 "POST-ICHEP 2015+2016 DS":22.07036,  # October 2016
                 "FULL 2015+2016 DS":36.0746          # December 2016
               }

parser.add_argument('inputDir', metavar='inputDir',type=str,
                    help='Path to the directory containing input histograms')
parser.add_argument('--channel', dest='channel', action='store', default='HIGHNJ', type=str, nargs='+',
                    help='The channel chosen. Full list of available options:\n{0}'.format(list_channel))
parser.add_argument('--category', dest='category', action='store', default=categories[0], type=str, nargs='+', choices=categories,
                    help='The category chosen. Can pass multiple space-separated arguments to this command-line option (picking amonge the above list). If this option is not specified, default will be \'{0}\''.format(categories[0]))
parser.add_argument('--variables', dest='variables', action='store', type=str, nargs='*',
                    help='List of variables to be considered. Use a space-separated list.')
parser.add_argument('--closure', dest='closure', action='store_true', default=False,
                    help="Check yields for MC closure test")
parser.add_argument("--lumi", dest="lumi", action="store", type=float, default=luminosities["FULL 2015+2016 DS"],
                    help="The luminosity of the dataset. Pick one of these values: ==> " + ",".join( "{0} ({1})".format( lumi, tag ) for tag, lumi in luminosities.iteritems() ) + ". Default is {0}".format(luminosities["FULL 2015+2016 DS"] ) )
parser.add_argument("--mergeOverflow", dest="mergeOverflow", action="store_true", default=False,
                    help="If this option is used, the script will assume the overflow has been merged already w/ the last visible bin. Default is False.")

args = parser.parse_args()

from ROOT import gROOT, gStyle, TCanvas, TPad, TH1D, THStack, TFile, TLegend, TLatex, TLine, Double, kTeal, kGray, kBlack, kBlue, kOrange, kWhite, kViolet, kRed

gROOT.Reset()
gROOT.LoadMacro("$HOME/RootUtils/AtlasStyle.C")
from ROOT import SetAtlasStyle
SetAtlasStyle()

from Plotter.BackgroundTools import makePoissonErrors, integrate

# Store sys integral for each systematic name

g_sys_dict = {}

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


def getTotFakeUncertainty( nominal, stat, flav ):

    if ( "HIGHNJ" in args.channel ):
        if flav == "ElEl" : non_closure = 0.1711 * nominal # Updated on v27
        if flav == "OF"   : non_closure = 0.1316 * nominal # Updated on v27
        if flav == "MuMu" : non_closure = 0.0461 * nominal # Updated on v27
    elif ( "LOWNJ" in args.channel ):
        if flav == "ElEl" : non_closure = 0.0208 * nominal # Updated on v27
        if flav == "OF"   : non_closure = 0.0481 * nominal # Updated on v27
        if flav == "MuMu" : non_closure = 0.0281 * nominal # Updated on v27

    # If you are doing closure, do not consider closure syst!

    if args.closure:
        non_closure = 0.0

    g_sys_dict["Non_Closure"]      = non_closure
    g_sysgroup_dict["Non_Closure"] = [non_closure]

    # This prints out sorting systematics from smallest to largest

    print("\t\tTot. yields w/ systematics (ungrouped):\n")
    print ("\t\tIntegral = {0:.2f}\n\t\t+- {1:.2f} [{2:.2f} %] (Sidebands Stat)\n\t\t+-".format(nominal, stat, (stat/nominal)*100) + "\t\t+-".join( " {0:.2f} [{1:.2f} %] ({2}) \n".format( g_sys_dict[key], (g_sys_dict[key]/nominal)*100, key ) for key in sorted( g_sys_dict, key=g_sys_dict.get ) ) )

    print("")

    # Print sum in quadrature of syst for each syst group from smallest to largest

    sq_list = []
    for sg, values in g_sysgroup_dict.iteritems():
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
    sq = sumQuadrature( toterrlist )

    sq_NO_STAT = sumQuadrature( list(g_sys_dict.values()) )

    print ("\t\tIntegral = {0:.2f} +- {1:.2f} [{2:.2f} %] (TOTAL UNCERTAINTY)".format(nominal, sq, (sq/nominal)*100))
    print ("\t\t                   +- {0:.2f} [{1:.2f} %] (TOTAL SYST. UNCERTAINTY)".format(sq_NO_STAT, (sq_NO_STAT/nominal)*100))

    # Print out results in LaTeX friendly format!

    print("\\begin{{table}}\n\\begin{{center}}\n\\begin{{tabular}}{{ll}}\n\\toprule\n\\multicolumn{{2}}{{c}}{{MY_FLAVOUR - $2\\leq N_{{jets}} \\leq 3$ VR, SR}} \\\\ \n\\midrule\nFakes (MM) = & {0:.2f} $\\pm$ \\\\\n & {1:.2f} $[{2:.2f}\%]$ (Sidebands stat.) $\pm$ \\\\".format(nominal, stat, (stat/nominal)*100))
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

    return sq, sq_NO_STAT


def clearDicts():

    g_sys_dict.clear()
    g_sysgroup_dict .clear()
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

    extensions = [".png",".eps"]
    for ext in extensions:
        c.SaveAs( outpath + "/" + flav + "_" + var + "_SysUncertStackHist" + ext )

    foutput = TFile(outpath+"/"+flav+"_"+var+"_SysUncertHists"+".root","RECREATE")
    for hist in systhists:
        hist.Write()
    foutput.Close()


def makeSysPlots( flav, var, observedhist, expectedhist ):

    gROOT.SetBatch(True)

    c = TCanvas("c1","Temp",50,50,600,600)

    pad1 = TPad("pad1", "", 0, 0.25, 1, 1)
    pad2 = TPad("pad2", "", 0, 0,   1, 0.25)
    pad1.SetBottomMargin(0.02)
    pad2.SetBottomMargin(0.4)
    pad1.Draw()
    pad2.Draw()

    new_expectedhist = expectedhist.Clone(expectedhist.GetName())

    # For expected hist, set the bin error as the sum in quadrature of stat (on expected) + syst (on fakes)

    for ibin in range(1,new_expectedhist.GetNbinsX()+2):
        uncertlist = [ g_sq_unc_bins_NO_STAT[ibin], new_expectedhist.GetBinError(ibin) ]
        new_expectedhist.SetBinError(ibin, sumQuadrature(uncertlist) )

    if observedhist:
        observedgr = makePoissonErrors(observedhist)
        observedgr.SetMarkerSize(0.8) # (1.2)
        observedgr.SetLineColor(1)
        observedgr.SetLineWidth(2)
        observedgr.SetMarkerStyle(20)
        observedgr.SetLineStyle(1)

    new_expectedhist.SetLineWidth(2)
    new_expectedhist.SetLineStyle(1)
    new_expectedhist.SetLineColor(kBlue-4)
    new_expectedhist.SetFillColor(kWhite)
    new_expectedhist.SetFillStyle(1001)

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

    err.SetMaximum( ( ymax + binmaxerr ) * 1.5 )
    err.SetMinimum(0)

    # --------------------------

    # Stat + sys error on expected (ratio)

    ratio_err = err.Clone("RatioErr")
    ratio_err.SetXTitle(new_expectedhist.GetXaxis().GetTitle())
    ratio_err.SetYTitle("Data/Exp.")
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

    ratio_err.Divide(err)

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

        # Trick to rescale

        ymax_ratio = ratio_obs_exp.GetMaximum()
        ymin_ratio = ratio_obs_exp.GetMinimum()

    #print("FLAV: {0} - MC/MM max = {1:.2f} - MC/MM min = {2:.2f}".format(flav,ymax_ratio,ymin_ratio))

    # For ratio error, get the bin with maximum uncertainty and the uncertainty value

    max_unc = 0
    for ibin in range(1, ratio_err.GetNbinsX()+1):
        this_unc = ratio_err.GetBinError(ibin)
      	if this_unc > max_unc:
	    max_unc = this_unc

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
    err.GetXaxis().SetLabelSize(0)
    err.GetXaxis().SetLabelOffset(999)

    err.Draw("E2")
    new_expectedhist.Draw("HIST SAME")
    if observedhist:
        observedgr.Draw("P SAME")

    legend = TLegend(0.6,0.7,0.8,0.9) # (x1,y1 (--> bottom left corner), x2, y2 (--> top right corner) )
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
    leg_ATLAS.DrawLatex(0.19,0.85,"#bf{#it{ATLAS}} Work In Progress")
    leg_lumi.DrawLatex(0.19,0.75,"#sqrt{{s}} = 13 TeV, #int L dt = {0:.1f} fb^{{-1}}".format(args.lumi))

    pad2.cd()
    ratio_err.GetYaxis().SetRangeUser(0.0, 2.0)
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

    extensions = [".png",".eps",".root"]
    for ext in extensions:
        c.SaveAs( outpath + "/" + flav + "_" + var + "_Sys" + ext )


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

    # For MM hist, set the bin error as the sum in quadrature of stat+syst

    for ibin in range(1,new_MM_hist.GetNbinsX()+1):
        new_MM_hist.SetBinError(ibin, g_sq_unc_bins[ibin] )

    MC_hist.SetLineStyle(2)
    MC_hist.SetMarkerSize(0.8)
    MC_hist.SetLineColor(1)
    MC_hist.SetMarkerStyle(20)
    MC_hist.SetLineWidth(1)

    new_MM_hist.SetLineWidth(2)
    new_MM_hist.SetLineStyle(1)
    new_MM_hist.SetLineColor(kViolet-4)
    new_MM_hist.SetFillColor(kWhite)
    new_MM_hist.SetFillStyle(1001)

    err = new_MM_hist.Clone("tot_uncertainty")
    err.SetFillColor(kGray+3)
    err.SetLineColor(10)
    err.SetFillStyle(3004)
    err.SetMarkerSize(0)

    # Trick to rescale:

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
    leg_ATLAS.DrawLatex(0.19,0.85,"#bf{#it{ATLAS}} Work In Progress")
    leg_lumi.DrawLatex(0.19,0.75,"#sqrt{{s}} = 13 TeV, #int L dt = {0:.1f} fb^{{-1}}".format(args.lumi))

    pad2.cd()
    #ratio_err.GetYaxis().SetRangeUser(0.0, 2.0)
    ratio_err.GetYaxis().SetRangeUser(0.5, 1.5)
    ratio_err.Draw("E2")
    ratio_MC_MM.Draw("PE SAME")

    refl = TLine(ratio_err.GetBinLowEdge(1), 1.0, ratio_err.GetBinLowEdge(ratio_err.GetNbinsX()+1), 1.0)
    refl.SetLineStyle(2)
    refl.SetLineColor(kRed)
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

        if ( "HIGHNJ" in args.channel ):  region = "SS_SR_DataDriven"
        elif ( "LOWNJ" in args.channel ): region = "SS_LowNJetCR_DataDriven"

    else:

        # -------------
        # TTBAR CLOSURE
        #--------------

        if ( "LOWNJ" in args.channel ):
            region    = "SS_SR_LowJet_DataDriven_Closure"
        elif ( "HIGHNJ" in args.channel ):
            region    = "SS_SR_HighJet_DataDriven_Closure"
        elif ( "ALLNJ" in args.channel ):
            region    = "SS_SR_AllJet_DataDriven_Closure"

    flavour_list = []
    if "ALL" in args.category:
        flavour_list.extend(["ElEl", "MuMu", "OF"])
    else:
        flavour_list.extend(args.category)

    print("Looking at variables : [" + ",".join( "{0}".format(v) for v in var_list ) + "]")

    for var_name in var_list:

	print ("\n\tVariable: {0}\n".format(var_name))

    	for flav in flavour_list:

    	    clearDicts()

    	    print ("\tFlavour region: {0}\n".format(flav))

            if flav == "ElEl" and "Mu" in var_name: continue
            if flav == "MuMu" and "El" in var_name: continue
            if flav == "OF" and any ( f in var_name for f in ["El1","Mu1"] ): continue

    	    filename = inputpath + flav + region + "/" + flav + region + "_" + var_name + ".root"
    	    myfile = TFile(filename)

    	    print("\tLooking at file: {0}".format(filename))

    	    fakes_nominal = myfile.Get("fakesbkg")
    	    fakes_syst = {}
    	    for key in myfile.GetListOfKeys():

    		keyname = key.GetName()
    		if not ( "fakesbkg_" in keyname ): continue

                # if not any( k == keyname for k in ["fakesbkg_MMsys_Fake_El_Pt_Stat_up","fakesbkg_MMsys_Fake_El_Pt_Stat_dn","fakesbkg_MMsys_Real_El_Pt_Stat_up","fakesbkg_MMsys_Real_El_Pt_Stat_dn"] ): continue

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

    	    print ("\tFakes:\n")

    	    for sys, sysgroup in fakes_syst.iteritems():
                debugflag = False
                fakes_up = myfile.Get( "fakesbkg_" + sys + "_up")
                fakes_dn = myfile.Get( "fakesbkg_" + sys + "_dn")
                if debugflag: print("\tsys: {0}, sysgroup: {1}\n".format(sys, sysgroup))
                fakes, fakes_stat_err = getYields(fakes_nominal,fakes_up,fakes_dn, sys, sysgroup, debug=debugflag)

    	    fakes_tot_err, fakes_tot_err_NO_STAT = getTotFakeUncertainty( fakes, fakes_stat_err, flav )

            saveSystHistogram(flav,var_name,fakes_nominal,fakes_tot_err_NO_STAT)

    	    if args.closure:
    		ttbar_nominal = myfile.Get("ttbarbkg")

    		print ("\n\tTTbar: \n")
    		ttbar, ttbar_err = getYields(ttbar_nominal)

    		closure     = ( (fakes - ttbar) / (ttbar) ) * 100
    		closure_err = math.sqrt( ( ( math.pow(fakes,2.0) / math.pow(ttbar,4.0) ) * math.pow(ttbar_err,2.0) ) + ( math.pow(fakes_tot_err,2.0) / math.pow(ttbar,2.0) ) ) * 100

    		makeSysPlotsClosure(flav,var_name,ttbar_nominal,fakes_nominal)

    		print("\n\tNon-closure ((fakes-ttbar)/ttbar) = {0:.2f} [%] +- {1:.2f} [%]".format(closure,closure_err))

    	    expected_nominal = myfile.Get("expectedbkg")

    	    print ("\n\tExpected: \n")
    	    exp, exp_err = getYields(expected_nominal)

    	    if not args.closure:

    		chargemisid = myfile.Get("qmisidbkg")
    		if chargemisid:
    		    print ("\n\tQMisID: \n")
    		    qmisid, qmisid_err = getYields(chargemisid)

    		observed = myfile.Get("observed")
    		if observed:
    		    print ("\n\tObserved: \n")
    		    obs, obs_err = getYields(observed)

                makeSysPlots(flav,var_name,observed,expected_nominal)

    		signal = myfile.Get("signal")

    		print ("\n\tSignal: \n")
    		sig, sig_err = getYields(signal)
