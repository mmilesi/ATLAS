#######################################
# ----------------------
# full xAOD - mc15 13TeV  
# ----------------------
#
# ------------
# 25 ns
# ------------
#
# ttbar non all had - 410000
#
#test-HTopMultilep --inDSType DxAOD-2015-13TeV --inDir /data/mmilesi/HTopMultileptonsTestSamples/MC15/r6282/ --inDSName mc15_13TeV.410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.merge.AOD.e3698_s2608_s2183_r6765_r6282 --outDir out_xAOD_HTopMultilep_MC --maxEvents 1000
#
# ------------
# 50 ns
# ------------
#
# ttbar non all had - 410000
#
#test-HTopMultilep --inDSType DxAOD-2015-13TeV --inDir /data/mmilesi/HTopMultileptonsTestSamples/MC15/r6264/ --inDSName mc15_13TeV.410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.merge.AOD.e3698_s2608_s2183_r6630_r6264 --outDir out_xAOD_HTopMultilep_MC --maxEvents 1000
#
# Z-->mumu PowhegPythia inclusive - 361107
#
#test-HTopMultilep --inDSType DxAOD-2015-13TeV --inDir /data/mmilesi/HTopMultileptonsTestSamples/MC15/r6264/ --inDSName mc15_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.merge.AOD.e3601_s2576_s2132_r6630_r6264 --outDir out_xAOD_HTopMultilep_MC --maxEvents 1000
#
# Z-->ee PowhegPythia inclusive - 361106
#
#test-HTopMultilep --inDSType DxAOD-2015-13TeV --inDir /data/mmilesi/HTopMultileptonsTestSamples/MC15/r6264/ --inDSName mc15_13TeV.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee.merge.AOD.e3601_s2576_s2132_r6630_r6264 --outDir out_xAOD_HTopMultilep_MC --maxEvents -1
#
# Z-->ee Sherpa 0 < pT(Z) < 70 GeV - 361372
#
#test-HTopMultilep --inDSType DxAOD-2015-13TeV --inDir /data/mmilesi/HTopMultileptonsTestSamples/MC15/r6264/ --inDSName mc15_13TeV.361372.Sherpa_CT10_Zee_Pt0_70_CVetoBVeto.merge.AOD.e3651_s2586_s2174_r6793_r6264 --outDir out_xAOD_HTopMultilep_MC --maxEvents -1
#
#######################################
# ----------------------
# mc15 13TeV DAOD HIGG8D1 
# ----------------------
#
# --------------
# 20.1.5.1
# --------------
#
#test-HTopMultilep --inDSType DxAOD-2015-13TeV --inDir /data/mmilesi/HTopMultileptonsTestSamples/HIGG8D1_20.1.5.1/ --inDSName mc15_13TeV.410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.merge.DAOD_HIGG8D1.e3698_s2608_s2183_r6630_r6264_p2363 --outDir out_DxAOD_HTopMultilep_MC --maxEvents 2000
#test-HTopMultilep --inDSType DxAOD-2015-13TeV --inDir /data/mmilesi/HTopMultileptonsTestSamples/HIGG8D1_20.1.5.1/ --inDSName mc15_13TeV.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee.merge.DAOD_HIGG8D1.e3601_s2576_s2132_r6630_r6264_p2363 --outDir out_DxAOD_HTopMultilep_MC --maxEvents 10000
#
# --------------
# 20.1.5.4
# --------------
#
#test-HTopMultilep --inDSType DxAOD-2015-13TeV --inDir /data/mmilesi/HTopMultileptonsTestSamples/HIGG8D1_20.1.5.4/ --inDSName mc15_13TeV.410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.merge.DAOD_HIGG8D1.e3698_s2608_s2183_r6630_r6264_p2377 --outDir out_DxAOD_HTopMultilep_MC --maxEvents 1000
#
# ----------------
# 20.1.6.3 - 25 ns
#
# ----------------
#
#
test-HTopMultilep --inDSType DxAOD-2015-13TeV --inDir /data/mmilesi/HTopMultileptonsTestSamples/HIGG8D1_20.1.6.3/ --inDSName mc15_13TeV.410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.merge.DAOD_HIGG8D1.e3698_s2608_s2183_r6765_r6282_p2413 --outDir out_DxAOD_HTopMultilep_MC --maxEvents 1000
#
#test-HTopMultilep --inDSType DxAOD-2015-13TeV --inDir /data/mmilesi/HTopMultileptonsTestSamples/HIGG8D1_20.1.6.3/ --inDSName mc15_13TeV.361068.Sherpa_CT10_llvv.merge.DAOD_HIGG8D1.e3836_s2608_s2183_r6869_r6282_p2413 --outDir out_DxAOD_HTopMultilep_MC --maxEvents 1000
#
#test-HTopMultilep --inDSType DxAOD-2015-13TeV --inDir /data/mmilesi/HTopMultileptonsTestSamples/testing/ --inDSName ttbar_2 --outDir test_trigmatching --maxEvents -1
#######################################
# --------------------------
# data15 13TeV DAOD 
# --------------------------
#
# --------------
# 20.1.5.1
# --------------
#test-HTopMultilep --inDSType DxAOD-2015-13TeV --inDir /data/mmilesi/HTopMultileptonsTestSamples/HIGG8D1_20.1.5.1/ --inDSName data15_13TeV.00267639.physics_Main.merge.DAOD_HIGG8D1.f598_m1441_p2361 --outDir out_DxAOD_HTopMultilep_DATA --maxEvents -1
# --------------
# 20.1.5.4
# --------------
#test-HTopMultilep --inDSType DxAOD-2015-13TeV --inDir /data/mmilesi/HTopMultileptonsTestSamples/HIGG8D1_20.1.5.4/ --inDSName data15_13TeV.00271421.physics_Main.merge.DAOD_HIGG8D1.f611_m1463_p2375 --outDir out_DxAOD_HTopMultilep_DATA --maxEvents 1000
#
# ----------------
# 20.1.6.3 - 25 ns
# ----------------
#
#test-HTopMultilep --inDSType DxAOD-2015-13TeV --inDir /data/mmilesi/HTopMultileptonsTestSamples/HIGG8D1_20.1.6.3/ --inDSName data15_13TeV.00276329.physics_Main.merge.DAOD_HIGG8D1.f620_m1480_p2411 --outDir out_DxAOD_HTopMultilep_DATA --maxEvents 1000
#
#######################################
