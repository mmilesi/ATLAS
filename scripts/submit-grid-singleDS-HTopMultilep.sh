#######################################
# DAOD HIGG8D1 mc14_8TeV
# ----------------------
#test-HTopMultilep --outDir outputtest_grid_DxAOD-DC14-8TeV --grid 1 --gridUser mmilesi -inGridDSName mc14_8TeV.161305.Pythia8_AU2CTEQ6L1_ttH125_WWinclusive.merge.DAOD_HIGG8D1.e1530_s1933_s1911_r5591_r5625_p1854_tid04985328_00 -inDSType DxAOD-DC14-8TeV --outDSID HTopMultilep.test21_SINGLEDS_DxAOD
#######################################
# DAOD HIGG8D1 data12_8TeV egamma
# -------------------------------
#test-HTopMultilep --outDir outputtest_grid_DxAOD-DC14-8TeV --grid 1 --gridUser mmilesi -inGridDSName data12_8TeV.00205113.physics_Egamma.merge.DAOD_HIGG8D1.r5723_p1751_p2309_p1851_tid04984115_00 -inDSType DxAOD-DC14-8TeV --outDSID HTopMultilep.test21_SINGLEDS_DxAOD
#######################################
# DAOD HIGG8D1 mc14_13TeV
# ----------------------
test-HTopMultilep --outDir outputtest_grid_DxAOD-2015-13TeV --grid 1 --gridUser mmilesi -inGridDSName  mc15_13TeV.410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.merge.DAOD_HIGG8D1.e3698_s2608_s2183_r6630_r6264_p2352_tid05526323_00 -inDSType DxAOD-2015-13TeV --outDSID HTopMultilep.001_SINGLEDS_DxAOD

