import math

NOMINAL_EFF_CONV = 0.057
NOMINAL_EFF_NP   = 0.009

def __alpha_approx__(f, f_CR):

    # This works w/ the following approximation:
    # eff_conv >> eff_NP

    alpha = (f - f_CR)/f_CR

    return alpha


def __alpha__(f, f_CR, eff_np=NOMINAL_EFF_NP, eff_conv=NOMINAL_EFF_CONV):

    # This works w/o the above approximation

    alpha = ( ( ( 1 - f ) * eff_np + f * eff_conv ) / ( ( 1 - f_CR ) * eff_np + f_CR * eff_conv ) ) - 1

    return alpha

def __alpha_unc__(f, f_CR, key="", variation=[], symmetrise=False ):

    alpha_nominal = __alpha__(f, f_CR)

    alpha_var_up = alpha_var_dn = alpha_nominal
    if "frac_conv" in key:
        alpha_var_up =  __alpha__(f+variation[0], f_CR+variation[1])
        alpha_var_dn =  __alpha__(f-variation[0], f_CR-variation[1])
    elif key == "eff_np":
        alpha_var_up =   __alpha__(f, f_CR, eff_np=NOMINAL_EFF_NP+variation[0])
        alpha_var_dn =   __alpha__(f, f_CR, eff_np=NOMINAL_EFF_NP-variation[0])
    elif key == "eff_conv":
        alpha_var_up =   __alpha__(f, f_CR, eff_conv=NOMINAL_EFF_CONV+variation[0])
        alpha_var_dn =   __alpha__(f, f_CR, eff_conv=NOMINAL_EFF_CONV-variation[0])

    # print("key: {0}".format(key))
    # print("alpha nominal = {0:.2f} (sys var UP : {1:.2f}, sys var DN : {2:.2f})".format(alpha_nominal,alpha_var_up,alpha_var_dn))

    sigma_up = abs( alpha_var_up - alpha_nominal )
    sigma_dn = abs( alpha_var_dn - alpha_nominal )

    # print("alpha nominal = {0:.2f} (sigma UP : {1:.2f}, sigma DN : {2:.2f})".format(alpha_nominal,sigma_up,sigma_dn))

    if alpha_var_up < alpha_nominal and alpha_var_dn > alpha_nominal:
        sigma_up, sigma_dn = sigma_dn, sigma_up

    if symmetrise:
        sigma_up = sigma_dn = ( sigma_up + sigma_dn ) / 2.0

    return (sigma_up,sigma_dn)

def __alpha_unc_tot__(unclist_up=[], unclist_dn=[], symmetrise=False ):

    tot_up = tot_dn = 0

    for unc_up in unclist_up:
        tot_up += pow(unc_up,2)
    for unc_dn in unclist_dn:
        tot_dn += pow(unc_dn,2)

    tot_up = math.sqrt(tot_up)
    tot_dn = math.sqrt(tot_dn)

    if symmetrise:
        tot_up = tot_dn = ( tot_up + tot_dn ) / 2.0

    return (tot_up,tot_dn)


if __name__ == "__main__":

    # v28
    #
    # Values used for Higgs Approval on 20/07/17
    # Contained a bug in the ee fractions
    # f_CR = 0.28
    # f_ee_PreMVA = 0.43
    # f_em_PreMVA = 0.28
    # f_me_PreMVA = 0.39
    # f_OF_PreMVA = 0.33
    # f_ee_LJ = 0.36
    #
    # These contain a fix for ee conversion fractions
    # f_CR = 0.28
    # f_ee_PreMVA = 0.53
    # f_em_PreMVA = 0.28
    # f_me_PreMVA = 0.39
    # f_OF_PreMVA = 0.33
    # f_ee_LJ = 0.49

    # # v29 - PP8
    # # Prod. 25_07_17, pT(l) > 15 GeV
    # # f_CR = 0.22 # <--- WRONG! don't use this
    # f_CR = 0.25 # --> uses me events
    # # v29 - PP6
    # # f_CR = 0.20
    # f_ee_PreMVA = 0.42
    # f_em_PreMVA = 1
    # f_me_PreMVA = 1
    # f_OF_PreMVA = 0.24
    # f_ee_LJ = 0.33
    # f_Xee_PreMVA = 0.496
    # f_XOF_PreMVA = 0.278

    # v29 - PP8
    # Prod. 26_07_17, pT(l) > 20 GeV
    f_CR = 0.27 # --> uses me events
    f_ee_PreMVA = 0.46
    f_em_PreMVA = 1
    f_me_PreMVA = 1
    f_OF_PreMVA = 0.30
    f_ee_LJ = 0.43
    f_Xee_PreMVA = 0.496 # <--- ask Ximo for new vals
    f_XOF_PreMVA = 0.278 # <--- ask Ximo for new vals

    # 40% uncertainty on all f's (based on ID material study and ttgamma Xsec measurement at 8 TeV)
    # 50% uncertainty on eff_np
    # 50% uncertainty on eff_conv

    unc_dict= {
        "frac_conv_ee_PreMVA" :[0.4*f_ee_PreMVA,0.4*f_CR],
        "frac_conv_ee_LJ"     :[0.4*f_ee_LJ,0.4*f_CR],
        "frac_conv_OF_PreMVA" :[0.4*f_OF_PreMVA,0.4*f_CR],
        "eff_np"              :[0.5*NOMINAL_EFF_NP],
        "eff_conv"            :[0.5*NOMINAL_EFF_CONV]
    }

    # ---------------------------------

    alpha_ee_PreMVA              = __alpha__(f_ee_PreMVA, f_CR)
    alpha_ee_PreMVA_unc_f        = __alpha_unc__(f_ee_PreMVA, f_CR, key="frac_conv_ee_PreMVA", variation=unc_dict["frac_conv_ee_PreMVA"], symmetrise=False )
    alpha_ee_PreMVA_unc_eff_np   = __alpha_unc__(f_ee_PreMVA, f_CR, key="eff_np", variation=unc_dict["eff_np"], symmetrise=False )
    alpha_ee_PreMVA_unc_eff_conv = __alpha_unc__(f_ee_PreMVA, f_CR, key="eff_conv", variation=unc_dict["eff_conv"], symmetrise=False )
    alpha_ee_PreMVA_tot_unc      = __alpha_unc_tot__(unclist_up=[alpha_ee_PreMVA_unc_f[0],alpha_ee_PreMVA_unc_eff_np[0],alpha_ee_PreMVA_unc_eff_conv[0]],
                                                     unclist_dn=[alpha_ee_PreMVA_unc_f[1],alpha_ee_PreMVA_unc_eff_np[1],alpha_ee_PreMVA_unc_eff_conv[1]],
                                                     symmetrise=True )
    print("Pre-MVA extrapolation:\n")
    print("ee, Pre-MVA : alpha = {0:.2f} +-\n\t\t+{1:.2f},-{2:.2f} (f)\n\t\t+{3:.2f},-{4:.2f} (eff_np)\n\t\t+{5:.2f},-{6:.2f} (eff_conv)\n\t\t+{7:.2f},-{8:.2f} ({9:.2f}%,-{10:.2f}%) (TOT.)".format(alpha_ee_PreMVA,alpha_ee_PreMVA_unc_f[0],alpha_ee_PreMVA_unc_f[1],alpha_ee_PreMVA_unc_eff_np[0],alpha_ee_PreMVA_unc_eff_np[1],alpha_ee_PreMVA_unc_eff_conv[0],alpha_ee_PreMVA_unc_eff_conv[1],alpha_ee_PreMVA_tot_unc[0],alpha_ee_PreMVA_tot_unc[1], alpha_ee_PreMVA_tot_unc[0]/alpha_ee_PreMVA*1e2,alpha_ee_PreMVA_tot_unc[1]/alpha_ee_PreMVA*1e2) )

    # ---------------------------------

    alpha_OF_PreMVA              = __alpha__(f_OF_PreMVA, f_CR)
    alpha_OF_PreMVA_unc_f        = __alpha_unc__(f_OF_PreMVA, f_CR, key="frac_conv_OF_PreMVA", variation=unc_dict["frac_conv_OF_PreMVA"], symmetrise=False )
    alpha_OF_PreMVA_unc_eff_np   = __alpha_unc__(f_OF_PreMVA, f_CR, key="eff_np", variation=unc_dict["eff_np"], symmetrise=False )
    alpha_OF_PreMVA_unc_eff_conv = __alpha_unc__(f_OF_PreMVA, f_CR, key="eff_conv", variation=unc_dict["eff_conv"], symmetrise=False )
    alpha_OF_PreMVA_tot_unc      = __alpha_unc_tot__(unclist_up=[alpha_OF_PreMVA_unc_f[0],alpha_OF_PreMVA_unc_eff_np[0],alpha_OF_PreMVA_unc_eff_conv[0]],
                                                     unclist_dn=[alpha_OF_PreMVA_unc_f[1],alpha_OF_PreMVA_unc_eff_np[1],alpha_OF_PreMVA_unc_eff_conv[1]],
                                                     symmetrise=True )
    print("Pre-MVA extrapolation:\n")
    print("OF, Pre-MVA : alpha = {0:.2f} +-\n\t\t+{1:.2f},-{2:.2f} (f)\n\t\t+{3:.2f},-{4:.2f} (eff_np)\n\t\t+{5:.2f},-{6:.2f} (eff_conv)\n\t\t+{7:.2f},-{8:.2f} ({9:.2f}%,-{10:.2f}%) (TOT.)".format(alpha_OF_PreMVA,alpha_OF_PreMVA_unc_f[0],alpha_OF_PreMVA_unc_f[1],alpha_OF_PreMVA_unc_eff_np[0],alpha_OF_PreMVA_unc_eff_np[1],alpha_OF_PreMVA_unc_eff_conv[0],alpha_OF_PreMVA_unc_eff_conv[1],alpha_OF_PreMVA_tot_unc[0],alpha_OF_PreMVA_tot_unc[1], alpha_OF_PreMVA_tot_unc[0]/alpha_OF_PreMVA*1e2,alpha_OF_PreMVA_tot_unc[1]/alpha_OF_PreMVA*1e2) )

    # ---------------------------------

    alpha_ee_LJ              = __alpha__(f_ee_LJ, f_CR)
    alpha_ee_LJ_unc_f        = __alpha_unc__(f_ee_LJ, f_CR, key="frac_conv_ee_LJ", variation=unc_dict["frac_conv_ee_LJ"], symmetrise=False )
    alpha_ee_LJ_unc_eff_np   = __alpha_unc__(f_ee_LJ, f_CR, key="eff_np", variation=unc_dict["eff_np"], symmetrise=False )
    alpha_ee_LJ_unc_eff_conv = __alpha_unc__(f_ee_LJ, f_CR, key="eff_conv", variation=unc_dict["eff_conv"], symmetrise=False )
    alpha_ee_LJ_tot_unc      = __alpha_unc_tot__(unclist_up=[alpha_ee_LJ_unc_f[0],alpha_ee_LJ_unc_eff_np[0],alpha_ee_LJ_unc_eff_conv[0]],
                                                     unclist_dn=[alpha_ee_LJ_unc_f[1],alpha_ee_LJ_unc_eff_np[1],alpha_ee_LJ_unc_eff_conv[1]],
                                                     symmetrise=True )
    print("LJ extrapolation:\n")
    print("ee, LJ : alpha = {0:.2f} +-\n\t\t+{1:.2f},-{2:.2f} (f)\n\t\t+{3:.2f},-{4:.2f} (eff_np)\n\t\t+{5:.2f},-{6:.2f} (eff_conv)\n\t\t+{7:.2f},-{8:.2f} ({9:.2f}%,-{10:.2f}%) (TOT.)".format(alpha_ee_LJ,alpha_ee_LJ_unc_f[0],alpha_ee_LJ_unc_f[1],alpha_ee_LJ_unc_eff_np[0],alpha_ee_LJ_unc_eff_np[1],alpha_ee_LJ_unc_eff_conv[0],alpha_ee_LJ_unc_eff_conv[1],alpha_ee_LJ_tot_unc[0],alpha_ee_LJ_tot_unc[1], alpha_ee_LJ_tot_unc[0]/alpha_ee_LJ*1e2,alpha_ee_LJ_tot_unc[1]/alpha_ee_LJ*1e2) )

    # ---------------------------------

    # # Alpha(pT)

    # regions = []

    # # pT(e) = [15,20,30,210+] GeV
    # # f_CR_pt        = [(1, 0.219), (2, 0.284), (3, 0.237), (4, 0.0)] # me, LJ CR
    # # f_ee_PreMVA_pt = [(1, 0.337), (2, 0.335), (3, 0.625), (4, 0.0)]
    # # regions.append((f_ee_PreMVA_pt,"ee, SR"))
    # # f_OF_PreMVA_pt = [(1, 0.095), (2, 0.203), (3, 0.38), (4, 1.0)]
    # # regions.append((f_OF_PreMVA_pt,"OF, SR"))
    # # f_ee_LJ_pt     = [(1, 0.143), (2, 0.337), (3, 0.546), (4, 1.0)]
    # # regions.append((f_ee_LJ_pt,"ee, LJ"))

    # # pT(e) = [15,30,210+] GeV
    # f_CR_pt        = [(1, 0.251), (2, 0.237), (3, 0.0)] # me, LJ CR
    # f_ee_PreMVA_pt = [(1, 0.336), (2, 0.625), (3, 0.0)]
    # regions.append((f_ee_PreMVA_pt,"ee, SR"))
    # f_OF_PreMVA_pt = [(1, 0.148), (2, 0.38), (3, 1.0)]
    # regions.append((f_OF_PreMVA_pt,"OF, SR"))
    # f_ee_LJ_pt     = [(1, 0.248), (2, 0.546), (3, 1.0)]
    # regions.append((f_ee_LJ_pt,"ee, LJ"))

    # print("\nAlpha(pT):\n")

    # for r in regions:
    #     rname = r[1]
    #     l =  r[0]
    #     print("Region: {0}\n".format(rname))
    #     for i,bin in enumerate(l):
    #         bin_idx = bin[0]
    #         bin_value = bin[1]
    #         print("\talpha[{0}] = {1:.2f} +- {2:.2f}".format(bin_idx,__alpha__(bin_value, f_CR_pt[i][1])[0],__alpha__(bin_value, f_CR_pt[i][1])[1]))
