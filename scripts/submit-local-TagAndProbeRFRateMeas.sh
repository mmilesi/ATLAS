#######################################
# full xAOD - mc14 13 TeV  
# -----------------------
#test-TagAndProbeRFRateMeas --inDSType xAOD-mc14_13TeV --inDir $HOME/work/xAOD_DC1413TeV/ --inDSName mc14_13TeV.110401.PowhegPythia_P2012_ttbar_nonallhad.merge.AOD.e2928_s1982_s2008_r5787_r5853_tid01597983_00 --outDir out_xAOD_TagAndProbeRFRateMeas_mc --maxEvents 10000
#######################################
# mc14 8TeV DAOD HIGG8D1 
# ----------------------
test-TagAndProbeRFRateMeas --inDSType DxAOD-DC14-8TeV --inDir $HOME/work/HIGG8D1_19.1.4.7/ --inDSName mc14_8TeV.117050.PowhegPythia_P2011C_ttbar.merge.DAOD_HIGG8D1.e1727_s1933_s1911_r5591_r5625_p1854_tid04985147_00 --outDir out_DxAOD_TagAndProbeRFRateMeas_mc --maxEvents 10000
#######################################
# data12 DAOD HIGG8D1
# -------------------
#test-TagAndProbeRFRateMeas --inDSType DxAOD-DC14-8TeV --inDir $HOME/work/HIGG8D1_19.1.4.9/ --inDSName data12_8TeV.00204158.physics_Muons.merge.DAOD_HIGG8D1.r5724_p1751_p2309_p1871_tid05250348_00 --outDir out_DxAOD_TagAndProbeRFRateMeas_data --maxEvents 10000
