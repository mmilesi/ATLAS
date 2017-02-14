#!/usr/bin/env python

if __name__ == "__main__":

    from ROOT import gROOT

    gROOT.Reset()
    gROOT.LoadMacro("$HOME/RootUtils/AtlasStyle.C")

    from ROOT import SetAtlasStyle
    SetAtlasStyle()

    gROOT.SetBatch(True)

    import EfficiencyPlotter_v24_ElNoIso

    # Use this one for v24 (has the right binning)
    if False:
        EfficiencyPlotter_v24_ElNoIso.plotRealElectron()
        EfficiencyPlotter_v24_ElNoIso.plotRealMuon()
        EfficiencyPlotter_v24_ElNoIso.plotFakeElectron()
        EfficiencyPlotter_v24_ElNoIso.plotFakeMuon()

    import EfficiencyPlotter_v26

    if False:

        EfficiencyPlotter_v26.plotRealElectron()
        EfficiencyPlotter_v26.plotRealMuon()
        EfficiencyPlotter_v26.plotFakeElectron()
        EfficiencyPlotter_v26.plotFakeMuon()

        EfficiencyPlotter_v26.plotRealElectron_CutBased()
        EfficiencyPlotter_v26.plotRealMuon_CutBased()
        EfficiencyPlotter_v26.plotFakeElectron_CutBased()
        EfficiencyPlotter_v26.plotFakeMuon_CutBased()

        # EfficiencyPlotter_v26.plotProbeElectronAssignEff_MVA()
        # EfficiencyPlotter_v26.plotProbeElectronAssignEff_CutBased()

    import EfficiencyPlotter_v26_VS_v24_ElNoIso

    if True:

        # EfficiencyPlotter_v26_VS_v24_ElNoIso.plotRealElectron()
        # EfficiencyPlotter_v26_VS_v24_ElNoIso.plotRealMuon()
        # EfficiencyPlotter_v26_VS_v24_ElNoIso.plotFakeElectron()
        # EfficiencyPlotter_v26_VS_v24_ElNoIso.plotFakeMuon()
        # EfficiencyPlotter_v26_VS_v24_ElNoIso.plotFakeElectron_NO_TRUTH()
        EfficiencyPlotter_v26_VS_v24_ElNoIso.plotFakeElectron_EtaBarrel()

    import TypeAndOriginPlots

    if False:

        samples       = ["ttbarbkg","wjetsbkg"]
        #lepSelections = ["FakeCRElT","FakeCRElL"]
        flavour       = "Mu"
        lepSelections = ["FakeCRMuT","FakeCRMuL"]
        jetSelections = ["ALLNJ","LOWNJ","HIGHNJ"]
        prodIDs       = ["25ns_v24","25ns_v24_ElNoIso"]
        normFactor    = 0 # 1.0

        for s in samples:
            print("\nLooking at sample: {0}".format(s))
            for ls in lepSelections:
                print("\n\tlepton selection: {0}".format(ls))
                for js in jetSelections:
                    print("\n\tjet selection: {0}".format(js))
                    for pid in prodIDs:
                        print("\n\t\tprodID: {0}\n".format(pid))
                        kwargs = {"flavour":flavour,"prodID":pid, "sample":s, "jetSelection":js, "lepSelection":ls}
                        TypeAndOriginPlots.plotTypeVSOrigin(normFactor=normFactor, **kwargs)
                        TypeAndOriginPlots.plotTypeVSNjets(normFactor=normFactor, **kwargs)
                        TypeAndOriginPlots.plotOriginVSNjets(normFactor=normFactor, **kwargs)

