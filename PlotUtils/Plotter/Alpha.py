import math

def __alpha_approx__(f, f_CR):

    # This works w/ the following approximation:
    # eff_conv >> eff_NP

    alpha = (f - f_CR)/f_CR

    # Assume 30% uncertainty on all f's
    sigma_f_CR = 0.30 * f_CR
    sigma_f = 0.30 * f
    #sigma_alpha = math.sqrt( (1/(f_CR*f_CR)) * (sigma_f*sigma_f) + ((f*f)/(f_CR*f_CR*f_CR*f_CR)) * (sigma_f_CR*sigma_f_CR) )
    sigma_alpha = alpha # Decided to be super-conservative and take the full rescaling effect as a 100% uncertainty

    return (alpha,sigma_alpha)


def __alpha__(f, f_CR):

    # This works w/o the above approximation

    eff_conv = 0.057
    eff_np   = 0.009

    alpha = ( ( ( 1 - f ) * eff_np + f * eff_conv ) / ( ( 1 - f_CR ) * eff_np + f_CR * eff_conv ) ) - 1
    sigma_alpha = alpha # Decided to be super-conservative and take the full rescaling effect as a 100% uncertainty

    return (alpha,sigma_alpha)

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

    # v29 - PP8
    # Prod. 25_07_17, pT(l) > 15 GeV
    # f_CR = 0.22 # <--- WRONG! don't use this
    f_CR = 0.25 # --> uses me events
    # v29 - PP6
    # f_CR = 0.20
    f_ee_PreMVA = 0.42
    f_em_PreMVA = 1
    f_me_PreMVA = 1
    f_OF_PreMVA = 0.24
    f_ee_LJ = 0.33
    f_Xee_PreMVA = 0.496
    f_XOF_PreMVA = 0.278

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

    print("Pre-MVA extrapolation:\n")
    # print("ee : alpha_approx = {0:.2f} +- {1:.2f}".format(__alpha_approx__(f_ee_PreMVA,f_CR)[0],__alpha_approx__(f_ee_PreMVA,f_CR)[1]))
    # print("em : alpha_approx = {0:.2f} +- {1:.2f}".format(__alpha_approx__(f_em_PreMVA,f_CR)[0],__alpha_approx__(f_em_PreMVA,f_CR)[1]))
    # print("me : alpha_approx = {0:.2f} +- {1:.2f}".format(__alpha_approx__(f_me_PreMVA,f_CR)[0],__alpha_approx__(f_me_PreMVA,f_CR)[1]))
    # print("OF : alpha_approx = {0:.2f} +- {1:.2f}".format(__alpha_approx__(f_OF_PreMVA,f_CR)[0],__alpha_approx__(f_OF_PreMVA,f_CR)[1]))
    print("")
    print("ee, SR : alpha = {0:.2f} +- {1:.2f}".format(__alpha__(f_ee_PreMVA, f_CR)[0],__alpha__(f_ee_PreMVA,f_CR)[1]))
    #print("em : alpha = {0:.2f} +- {1:.2f}".format(__alpha__(f_em_PreMVA,f_CR)[0],__alpha__(f_em_PreMVA,f_CR)[1]))
    #print("me : alpha = {0:.2f} +- {1:.2f}".format(__alpha__(f_me_PreMVA,f_CR)[0],__alpha__(f_me_PreMVA,f_CR)[1]))
    print("OF, SR : alpha = {0:.2f} +- {1:.2f}".format(__alpha__(f_OF_PreMVA,f_CR)[0],__alpha__(f_OF_PreMVA,f_CR)[1]))
    print("Xee, SR : alpha = {0:.2f} +- {1:.2f}".format(__alpha__(f_Xee_PreMVA,f_CR)[0],__alpha__(f_Xee_PreMVA,f_CR)[1]))
    print("XOF, SR : alpha = {0:.2f} +- {1:.2f}".format(__alpha__(f_XOF_PreMVA,f_CR)[0],__alpha__(f_XOF_PreMVA,f_CR)[1]))
    print("\nNon-OF, LJ extrapolation:\n")
    #print("ee : alpha_approx = {0:.2f} +- {1:.2f}".format(__alpha_approx__(f_ee_LJ,f_CR)[0],__alpha_approx__(f_ee_LJ,f_CR)[1]))
    print("")
    print("ee, LJ : alpha = {0:.2f} +- {1:.2f}".format(__alpha__(f_ee_LJ,f_CR)[0],__alpha__(f_ee_LJ,f_CR)[1]))

    # Alpha(pT)

    regions = []

    # pT(e) = [15,20,30,210+] GeV
    # f_CR_pt        = [(1, 0.219), (2, 0.284), (3, 0.237), (4, 0.0)] # me, LJ CR
    # f_ee_PreMVA_pt = [(1, 0.337), (2, 0.335), (3, 0.625), (4, 0.0)]
    # regions.append((f_ee_PreMVA_pt,"ee, SR"))
    # f_OF_PreMVA_pt = [(1, 0.095), (2, 0.203), (3, 0.38), (4, 1.0)]
    # regions.append((f_OF_PreMVA_pt,"OF, SR"))
    # f_ee_LJ_pt     = [(1, 0.143), (2, 0.337), (3, 0.546), (4, 1.0)]
    # regions.append((f_ee_LJ_pt,"ee, LJ"))

    # pT(e) = [15,30,210+] GeV
    f_CR_pt        = [(1, 0.251), (2, 0.237), (3, 0.0)] # me, LJ CR
    f_ee_PreMVA_pt = [(1, 0.336), (2, 0.625), (3, 0.0)]
    regions.append((f_ee_PreMVA_pt,"ee, SR"))
    f_OF_PreMVA_pt = [(1, 0.148), (2, 0.38), (3, 1.0)]
    regions.append((f_OF_PreMVA_pt,"OF, SR"))
    f_ee_LJ_pt     = [(1, 0.248), (2, 0.546), (3, 1.0)]
    regions.append((f_ee_LJ_pt,"ee, LJ"))

    print("\nAlpha(pT):\n")

    for r in regions:
        rname = r[1]
        l =  r[0]
        print("Region: {0}\n".format(rname))
        for i,bin in enumerate(l):
            bin_idx = bin[0]
            bin_value = bin[1]
            print("\talpha[{0}] = {1:.2f} +- {2:.2f}".format(bin_idx,__alpha__(bin_value, f_CR_pt[i][1])[0],__alpha__(bin_value, f_CR_pt[i][1])[1]))
