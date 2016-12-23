#!/usr/bin/env python

if __name__ == "__main__":

    from ROOT import gROOT

    gROOT.Reset()
    gROOT.LoadMacro("$HOME/RootUtils/AtlasStyle.C")

    from ROOT import SetAtlasStyle
    SetAtlasStyle()

    gROOT.SetBatch(True)

    import EfficiencyPlotter, EfficiencyPlotter_ElNoIso, EfficiencyPlotter_ElNoIso_Rebinned

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

    if True:
        EfficiencyPlotter_ElNoIso_Rebinned.plotFakeElectron()
        EfficiencyPlotter_ElNoIso_Rebinned.plotFakeElectron_anyProbe()
        EfficiencyPlotter_ElNoIso_Rebinned.plotFakeElectron_BaselineElIso()
