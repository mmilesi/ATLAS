#!/usr/bin/env python

if __name__ == "__main__":

    from ROOT import gROOT

    gROOT.Reset()
    gROOT.LoadMacro("$HOME/RootUtils/AtlasStyle.C")

    from ROOT import SetAtlasStyle
    SetAtlasStyle()

    gROOT.SetBatch(True)

    import EfficiencyPlotter, EfficiencyPlotter_ElNoIso, EfficiencyPlotter_ElNoIso_Rebinned_1, EfficiencyPlotter_ElNoIso_Rebinned_2, EfficiencyPlotter_ElNoIso_Rebinned_3

    if False:
        EfficiencyPlotter.plotFakeElectron()
        EfficiencyPlotter.plotFakeElectron_anyProbe()
        EfficiencyPlotter.plotFakeMuon()
        EfficiencyPlotter.plotRealElectron()
        EfficiencyPlotter.plotRealMuon()
        EfficiencyPlotter.plotFakeElectron_diffTagSel()
        EfficiencyPlotter.plotFakeElectronAssignEff()

    if False:
        EfficiencyPlotter_ElNoIso.plotFakeElectron()
        EfficiencyPlotter_ElNoIso.plotFakeElectron_anyProbe()
        EfficiencyPlotter_ElNoIso.plotFakeElectron_BaselineElIso()
        EfficiencyPlotter_ElNoIso.plotRealElectron()
        EfficiencyPlotter_ElNoIso.plotRealElectron_BaselineElIso()

    if False:
        EfficiencyPlotter_ElNoIso_Rebinned_1.plotFakeElectron()
        EfficiencyPlotter_ElNoIso_Rebinned_1.plotFakeElectron_anyProbe()
        EfficiencyPlotter_ElNoIso_Rebinned_1.plotFakeElectron_BaselineElIso()

    # Use this one (has the right binning)
    if False:
        EfficiencyPlotter_ElNoIso_Rebinned_2.plotRealElectron()
        EfficiencyPlotter_ElNoIso_Rebinned_2.plotRealMuon()
        EfficiencyPlotter_ElNoIso_Rebinned_2.plotFakeElectron()
        EfficiencyPlotter_ElNoIso_Rebinned_2.plotFakeMuon()

    if False:
        EfficiencyPlotter_ElNoIso_Rebinned_3.plotFakeElectron()

    import TypeAndOriginPlots

    if True:

        samples       = ["ttbarbkg","wjetsbkg"]
        lepSelections = ["FakeCRElT"]#,"FakeCRElL"]
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
                        kwargs = {"prodID":pid, "sample":s, "jetSelection":js, "lepSelection":ls}
                        TypeAndOriginPlots.plotTypeVSOrigin(normFactor=normFactor, **kwargs)
                        TypeAndOriginPlots.plotTypeVSNjets(normFactor=normFactor, **kwargs)
                        TypeAndOriginPlots.plotOriginVSNjets(normFactor=normFactor, **kwargs)

