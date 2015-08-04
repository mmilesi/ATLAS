#######################################
# full xAOD - mc15 13TeV  
# ----------------------
#
# ttbar non all had - 410000
#
#test-HTopMultilep --inDSType DxAOD-2015-13TeV --inDir /data/mmilesi/HTopMultileptonsTestSamples/MC15/ --inDSName mc15_13TeV.410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.merge.AOD.e3698_s2608_s2183_r6630_r6264 --outDir out_xAOD_HTopMultilep_mc --maxEvents -1
#
# Z-->mumu PowhegPythia inclusive - 361107
#
#test-HTopMultilep --inDSType DxAOD-2015-13TeV --inDir /data/mmilesi/HTopMultileptonsTestSamples/MC15/ --inDSName mc15_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.merge.AOD.e3601_s2576_s2132_r6630_r6264 --outDir out_xAOD_HTopMultilep_mc --maxEvents 1000
#
# Z-->ee PowhegPythia inclusive - 361106
#
#test-HTopMultilep --inDSType DxAOD-2015-13TeV --inDir /data/mmilesi/HTopMultileptonsTestSamples/MC15/ --inDSName mc15_13TeV.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee.merge.AOD.e3601_s2576_s2132_r6630_r6264 --outDir out_xAOD_HTopMultilep_mc --maxEvents -1
#
# Z-->ee Sherpa 0 < pT(Z) < 70 GeV - 361372
#
#test-HTopMultilep --inDSType DxAOD-2015-13TeV --inDir /data/mmilesi/HTopMultileptonsTestSamples/MC15/ --inDSName mc15_13TeV.361372.Sherpa_CT10_Zee_Pt0_70_CVetoBVeto.merge.AOD.e3651_s2586_s2174_r6793_r6264 --outDir out_xAOD_HTopMultilep_mc --maxEvents -1
#######################################
# data15 13TeV DAOD 
# --------------------------
#
# --------------
# 20.1.5.1
# --------------
#test-HTopMultilep --inDSType DxAOD-2015-13TeV --inDir /data/mmilesi/HTopMultileptonsTestSamples/HIGG8D1_20.1.5.1/ --inDSName data15_13TeV.00267639.physics_Main.merge.DAOD_HIGG8D1.f598_m1441_p2361 --outDir out_DxAOD_HTopMultilep_data --maxEvents -1
# --------------
# 20.1.5.4
# --------------
test-HTopMultilep --inDSType DxAOD-2015-13TeV --inDir /data/mmilesi/HTopMultileptonsTestSamples/HIGG8D1_20.1.5.4/ --inDSName data15_13TeV.00271421.physics_Main.merge.DAOD_HIGG8D1.f611_m1463_p2375 --outDir out_DxAOD_HTopMultilep_data --maxEvents -1
#######################################
# mc15 13TeV DAOD HIGG8D1 
# ----------------------
#
# --------------
# 20.1.5.1
# --------------
#test-HTopMultilep --inDSType DxAOD-2015-13TeV --inDir /data/mmilesi/HTopMultileptonsTestSamples/HIGG8D1_20.1.5.1/ --inDSName mc15_13TeV.410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.merge.DAOD_HIGG8D1.e3698_s2608_s2183_r6630_r6264_p2363 --outDir out_DxAOD_HTopMultilep_mc --maxEvents 2000
#test-HTopMultilep --inDSType DxAOD-2015-13TeV --inDir /data/mmilesi/HTopMultileptonsTestSamples/HIGG8D1_20.1.5.1/ --inDSName mc15_13TeV.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee.merge.DAOD_HIGG8D1.e3601_s2576_s2132_r6630_r6264_p2363 --outDir out_DxAOD_HTopMultilep_mc --maxEvents 10000
# --------------
# 20.1.5.4
# --------------
#test-HTopMultilep --inDSType DxAOD-2015-13TeV --inDir /data/mmilesi/HTopMultileptonsTestSamples/HIGG8D1_20.1.5.4/ --inDSName mc15_13TeV.410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.merge.DAOD_HIGG8D1.e3698_s2608_s2183_r6630_r6264_p2377 --outDir out_DxAOD_HTopMultilep_mc --maxEvents -1 
#######################################
