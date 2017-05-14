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

        # EfficiencyPlotter.plotRealElectron()
        # EfficiencyPlotter.plotRealMuon()
        # EfficiencyPlotter.plotFakeElectron()
        # EfficiencyPlotter.plotFakeMuon()

        # EfficiencyPlotter.plotRealElectron_CutBased()
        # EfficiencyPlotter.plotRealMuon_CutBased()
        # EfficiencyPlotter.plotFakeElectron_CutBased()
        # EfficiencyPlotter.plotFakeMuon_CutBased()

        # EfficiencyPlotter.plotProbeElectronAssignEff_MVA()
        # EfficiencyPlotter.plotProbeElectronAssignEff_CutBased()

        EfficiencyPlotter.plotFakeElectron_NonPromptVSPhotonConv() # MC TTbar-based
        # EfficiencyPlotter.plotFakeElectron_NonPromptAndPhotonConvVSPhotonConv() # DD

    import EfficiencyPlotter_NBJetsGt2

    if False:

        EfficiencyPlotter_NBJetsGt2.plotRealElectron()
        EfficiencyPlotter_NBJetsGt2.plotRealMuon()
        EfficiencyPlotter_NBJetsGt2.plotFakeElectron()
        EfficiencyPlotter_NBJetsGt2.plotFakeMuon()

    import TypeAndOriginPlots

    if False:

        samples       = ["ttbarbkg"] # ["ttbarbkg","wjetsbkg"]
        flavours      = ["El"] # ["El","Mu"]
        # lepSelections = ["FakeCRElAntiT"]# ["FakeCRElT","FakeCRElL"]
        # lepSelections = ["FakeCRMuAntiT"] #["FakeCRMuT","FakeCRMuL"]
        #prodIDs       = ["25ns_v24","25ns_v24_ElNoIso"]
        prodIDs       = ["25ns_v27"]
        normFactor    = 0 # 1.0

        for s in samples:
            print("\nLooking at sample: {0}".format(s))
            for pid in prodIDs:
                print("\n\t\tprodID: {0}\n".format(pid))
                for flavour in flavours:
                    print("\n\t\t\tflavour: {0}\n".format(flavour))
                    kwargs = {"flavour":flavour,"prodID":pid, "sample":s}
                    # TypeAndOriginPlots.plotFakeOriginFrac2L3L(**kwargs)
                    # TypeAndOriginPlots.plotOriginVSNjets(normFactor=normFactor, **kwargs)
                    # TypeAndOriginPlots.plotOriginVSNBjets(normFactor=normFactor, **kwargs)
                    # TypeAndOriginPlots.plotOriginVSDistanceOtherLep(normFactor=normFactor, **kwargs)
                    TypeAndOriginPlots.plotOriginVSPt(normFactor=normFactor, **kwargs)
                    # TypeAndOriginPlots.plotOriginVSDistanceLepClosestJet(normFactor=normFactor, **kwargs)
                    # if flavour == "Mu":
                    #    TypeAndOriginPlots.plotOriginVSEta(normFactor=normFactor, **kwargs)
                    # elif flavour == "El":
                    #    TypeAndOriginPlots.plotOriginVSEtaBE2(normFactor=normFactor, **kwargs)
