
#!/usr/bin/env python

if __name__ == "__main__":

    from ROOT import gROOT

    gROOT.Reset()
    gROOT.LoadMacro("$HOME/RootUtils/AtlasStyle.C")

    from ROOT import SetAtlasStyle
    SetAtlasStyle()

    gROOT.SetBatch(True)

    import EfficiencyPlotter

    if True:

        EfficiencyPlotter.plotRealElectron()
        # EfficiencyPlotter.plotRealMuon()
        # EfficiencyPlotter.plotFakeElectron()
        # EfficiencyPlotter.plotFakeMuon()

        # EfficiencyPlotter.plotRealElectron_CutBased()
        # EfficiencyPlotter.plotRealMuon_CutBased()
        # EfficiencyPlotter.plotFakeElectron_CutBased()
        # EfficiencyPlotter.plotFakeMuon_CutBased()

        # EfficiencyPlotter.plotProbeElectronAssignEff_MVA()
        # EfficiencyPlotter.plotProbeElectronAssignEff_CutBased()

        EfficiencyPlotter.plotFakeElectron_NonPromptVSPhotonConv()

    import TypeAndOriginPlots

    if False:

        samples       = ["ttbarbkg","wjetsbkg"]
        #flavour       = "El"
        #lepSelections = ["FakeCRElAntiT"]# ["FakeCRElT","FakeCRElL"]
        flavour       = "Mu"
        lepSelections = ["FakeCRMuAntiT"] #["FakeCRMuT","FakeCRMuL"]
        jetSelections = ["ALLNJ","LOWNJ","HIGHNJ"]
        #prodIDs       = ["25ns_v24","25ns_v24_ElNoIso"]
        prodIDs       = ["25ns_v26"]
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


    if False:

        samples       = ["ttbarbkg"]
        flavour       = "El"
        prodIDs       = ["25ns_v27"]

        for s in samples:
            print("\nLooking at sample: {0}".format(s))
            for pid in prodIDs:
                print("\n\t\tprodID: {0}\n".format(pid))
                kwargs = {"flavour":flavour,"prodID":pid, "sample":s}
                TypeAndOriginPlots.plotFakeOriginFrac2L3L(**kwargs)


