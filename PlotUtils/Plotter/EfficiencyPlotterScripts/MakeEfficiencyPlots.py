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

    if False:
        EfficiencyPlotter_ElNoIso_Rebinned_2.plotFakeElectron()

    if True:
        EfficiencyPlotter_ElNoIso_Rebinned_3.plotFakeElectron()

