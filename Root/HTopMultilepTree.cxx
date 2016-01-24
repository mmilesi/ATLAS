#include "HTopMultilepAnalysis/HTopMultilepTree.h"

HTopMultilepTree :: HTopMultilepTree( TTree* tree, TFile* file, xAOD::TEvent* event, xAOD::TStore* store, const float units, bool debug, bool DC14 ) :
  HelpTreeBase(tree, file, event, store, units, debug, DC14 )
{
  Info("HTopMultilepTree", "Creating output TTree");
}

HTopMultilepTree :: ~HTopMultilepTree() {}

void HTopMultilepTree::AddEventUser(const std::string detailStrUser)
{

  // to get rid of warning when detailString is not used
  if ( m_debug ) { Info("AddEventUser()", "Adding branches w/ detail: %s", detailStrUser.c_str()); }

  // event variables
  m_tree->Branch("isMC",              &m_is_mc, "isMC/I");
  m_tree->Branch("ystar",             &m_ystar, "ystar/F");
  m_tree->Branch("categoryFlag",      &m_categoryFlag, "categoryFlag/i");
  m_tree->Branch("isSS01",            &m_isSS01, "isSS01/I");
  m_tree->Branch("isSS12",            &m_isSS12, "isSS12/I");
  m_tree->Branch("MMWeight",          &m_MMWeight);
  m_tree->Branch("FFWeight",          &m_FFWeight);
  m_tree->Branch("mll01",             &m_mll01, "mll01/F");
  m_tree->Branch("mll02",             &m_mll02, "mll02/F");
  m_tree->Branch("mll12",             &m_mll12, "mll12/F");
  m_tree->Branch("mlll012",           &m_mlll012, "mlll012/F");
  m_tree->Branch("mT_lep0MET",        &m_mT_lep0MET, "mT_lep0MET/F");
  m_tree->Branch("mT_lep1MET",        &m_mT_lep1MET, "mT_lep1MET/F");
  m_tree->Branch("isTT",              &m_isTT, "isTT/I");
  m_tree->Branch("isTL",              &m_isTL, "isTL/I");
  m_tree->Branch("isLT",              &m_isLT, "isLT/I");
  m_tree->Branch("isLL",              &m_isLL, "isLL/I");
  m_tree->Branch("isTM",              &m_isTM, "isTM/I");
  m_tree->Branch("isMT",              &m_isMT, "isMT/I");
  m_tree->Branch("isMM",              &m_isMM, "isMM/I");
  m_tree->Branch("isTelLmu",          &m_isTelLmu, "isTelLmu/I");
  m_tree->Branch("isLelTmu",          &m_isLelTmu, "isLelTmu/I");
  m_tree->Branch("isTmuLel",          &m_isTmuLel, "isTmuLel/I");
  m_tree->Branch("isLmuTel",          &m_isLmuTel, "isLmuTel/I");
  m_tree->Branch("isTelMmu",          &m_isTelMmu, "isTelMmu/I");
  m_tree->Branch("isMelTmu",          &m_isMelTmu, "isMelTmu/I");
  m_tree->Branch("isTmuMel",          &m_isTmuMel, "isTmuMel/I");
  m_tree->Branch("isMmuTel",          &m_isMmuTel, "isMmuTel/I");
  m_tree->Branch("isNonTightEvent",   &m_isNonTightEvent,    "isNonTightEvent/I");
  m_tree->Branch("isProbeElEvent",    &m_isProbeElEvent,     "isProbeElEvent/I");
  m_tree->Branch("isProbeMuEvent",    &m_isProbeMuEvent,     "isProbeMuEvent/I");

  if ( m_isMC ) {
    m_tree->Branch("weight_lepton_trig_HTop", &m_weight_lepton_trig_HTop);
    m_tree->Branch("weight_lepton_reco_HTop", &m_weight_lepton_reco_HTop);
    m_tree->Branch("weight_lepton_iso_HTop",  &m_weight_lepton_iso_HTop);
    m_tree->Branch("weight_lepton_ID_HTop",   &m_weight_lepton_ID_HTop);
    m_tree->Branch("weight_lepton_TTVA_HTop",   &m_weight_lepton_TTVA_HTop);
  }

}

/*
void HTopMultilepTree::AddTriggerUser(const std::string detailStrUser)
{

  // to get rid of warning when detailString is not used
  if ( m_debug ) { Info("AddTriggerUser()", "Adding branches w/ detail: %s", detailStrUser.c_str()); }

  // trigger variables
}
*/

void HTopMultilepTree::AddJetsUser(const std::string detailStrUser, const std::string jetName )
{

  // to get rid of warning when detailString is not used
  if ( m_debug ) { Info("AddJetsUser()", "Adding branches w/ detail: %s - Jet name: %s", detailStrUser.c_str(), jetName.c_str()); }

  // jet variables
  m_tree->Branch("jet_m",     &m_jet_m);
}

void HTopMultilepTree::AddMuonsUser(const std::string detailStrUser)
{

  // to get rid of warning when detailString is not used
  if ( m_debug ) { Info("AddMuonsUser()", "Adding branches w/ detail: %s", detailStrUser.c_str()); }

  // muon variables
  m_tree->Branch("muon_isTightSelected",               &m_muon_isTight);
  m_tree->Branch("muon_isMediumSelected",              &m_muon_isMedium);
  m_tree->Branch("muon_isTag",                         &m_muon_isTag);
  m_tree->Branch("muon_isOS",                          &m_muon_isOS);
  m_tree->Branch("muon_isClosestSS",                   &m_muon_isClosestSS);
  m_tree->Branch("muon_isTruthMatchedToMuon",          &m_muon_isTruthMatched);
  m_tree->Branch("muon_isChFlip",	               &m_muon_isChFlip);
  m_tree->Branch("muon_isBrem",	                       &m_muon_isBrem);
  m_tree->Branch("muon_truthType",                     &m_muon_truthType);
  m_tree->Branch("muon_truthPdgId",                    &m_muon_truthPdgId);
  m_tree->Branch("muon_truthOrigin",                   &m_muon_truthOrigin);
  m_tree->Branch("muon_truthStatus",                   &m_muon_truthStatus);
  m_tree->Branch("muon_ancestorTruthType",             &m_muon_ancestorTruthType);
  m_tree->Branch("muon_ancestorTruthPdgId",            &m_muon_ancestorTruthPdgId);
  m_tree->Branch("muon_ancestorTruthOrigin",           &m_muon_ancestorTruthOrigin);
  m_tree->Branch("muon_ancestorTruthStatus",           &m_muon_ancestorTruthStatus);
  // muon TAG variables
  m_tree->Branch("muon_tag_pt",	                       &m_muon_tag_pt);
  m_tree->Branch("muon_tag_eta",	               &m_muon_tag_eta);
  m_tree->Branch("muon_tag_trkd0sig",                  &m_muon_tag_trkd0sig);
  m_tree->Branch("muon_tag_ptvarcone30",               &m_muon_tag_ptvarcone30);
  m_tree->Branch("muon_tag_topoetcone20",              &m_muon_tag_topoetcone20);
  m_tree->Branch("muon_tag_isIsolated_Loose",	       &m_muon_tag_isIsolated_Loose);
  m_tree->Branch("muon_tag_isIsolated_FixedCutTightTrackOnly", &m_muon_tag_isIsolated_FixedCutTightTrackOnly);
  m_tree->Branch("muon_tag_isTightSelected",	       &m_muon_tag_isTight);
  m_tree->Branch("muon_tag_isMediumSelected",	       &m_muon_tag_isMedium);
  m_tree->Branch("muon_tag_isTruthMatchedToMuon",      &m_muon_tag_isTruthMatched);
  m_tree->Branch("muon_tag_isChFlip",	               &m_muon_tag_isChFlip);
  m_tree->Branch("muon_tag_isBrem",	               &m_muon_tag_isBrem);
  m_tree->Branch("muon_tag_truthType",                 &m_muon_tag_truthType);
  m_tree->Branch("muon_tag_truthPdgId",                &m_muon_tag_truthPdgId);
  m_tree->Branch("muon_tag_truthOrigin",               &m_muon_tag_truthOrigin);
  m_tree->Branch("muon_tag_truthStatus",               &m_muon_tag_truthStatus);
  m_tree->Branch("muon_tag_ancestorTruthType",         &m_muon_tag_ancestorTruthType);
  m_tree->Branch("muon_tag_ancestorTruthPdgId",        &m_muon_tag_ancestorTruthPdgId);
  m_tree->Branch("muon_tag_ancestorTruthOrigin",       &m_muon_tag_ancestorTruthOrigin);
  m_tree->Branch("muon_tag_ancestorTruthStatus",       &m_muon_tag_ancestorTruthStatus);
  // muon PROBE variables
  m_tree->Branch("muon_probe_pt",	               &m_muon_probe_pt);
  m_tree->Branch("muon_probe_eta",	               &m_muon_probe_eta);
  m_tree->Branch("muon_probe_trkd0sig",                &m_muon_probe_trkd0sig);
  m_tree->Branch("muon_probe_ptvarcone30",               &m_muon_probe_ptvarcone30);
  m_tree->Branch("muon_probe_topoetcone20",              &m_muon_probe_topoetcone20);
  m_tree->Branch("muon_probe_isIsolated_Loose",	         &m_muon_probe_isIsolated_Loose);
  m_tree->Branch("muon_probe_isIsolated_FixedCutTightTrackOnly", &m_muon_probe_isIsolated_FixedCutTightTrackOnly);
  m_tree->Branch("muon_probe_isTightSelected",         &m_muon_probe_isTight);
  m_tree->Branch("muon_probe_isMediumSelected",         &m_muon_probe_isMedium);
  m_tree->Branch("muon_probe_isTruthMatchedToMuon",    &m_muon_probe_isTruthMatched);
  m_tree->Branch("muon_probe_isChFlip",	               &m_muon_probe_isChFlip);
  m_tree->Branch("muon_probe_isBrem",	               &m_muon_probe_isBrem);
  m_tree->Branch("muon_probe_truthType",               &m_muon_probe_truthType);
  m_tree->Branch("muon_probe_truthPdgId",              &m_muon_probe_truthPdgId);
  m_tree->Branch("muon_probe_truthOrigin",             &m_muon_probe_truthOrigin);
  m_tree->Branch("muon_probe_truthStatus",             &m_muon_probe_truthStatus);
  m_tree->Branch("muon_probe_ancestorTruthType",       &m_muon_probe_ancestorTruthType);
  m_tree->Branch("muon_probe_ancestorTruthPdgId",      &m_muon_probe_ancestorTruthPdgId);
  m_tree->Branch("muon_probe_ancestorTruthOrigin",     &m_muon_probe_ancestorTruthOrigin);
  m_tree->Branch("muon_probe_ancestorTruthStatus",     &m_muon_probe_ancestorTruthStatus);
}

void HTopMultilepTree::AddElectronsUser(const std::string detailStrUser)
{

  // to get rid of warning when detailString is not used
  if ( m_debug ) { Info("AddElectronsUser()", "Adding branches w/ detail: %s", detailStrUser.c_str()); }

  // electron variables
  m_tree->Branch("el_crack",                           &m_electron_crack);
  m_tree->Branch("el_isTightSelected",                 &m_electron_isTight);
  m_tree->Branch("el_isMediumSelected",                &m_electron_isMedium);
  m_tree->Branch("el_isTag",                           &m_electron_isTag);
  m_tree->Branch("el_isOS",	                       &m_electron_isOS);
  m_tree->Branch("el_isClosestSS",                     &m_electron_isClosestSS);
  m_tree->Branch("el_isTruthMatchedToElectron",        &m_electron_isTruthMatched);
  m_tree->Branch("el_isChFlip",                        &m_electron_isChFlip);
  m_tree->Branch("el_isBrem",                          &m_electron_isBrem);
  m_tree->Branch("el_truthType",	               &m_electron_truthType);
  m_tree->Branch("el_truthPdgId",	               &m_electron_truthPdgId);
  m_tree->Branch("el_truthOrigin",                     &m_electron_truthOrigin);
  m_tree->Branch("el_truthStatus",                     &m_electron_truthStatus);
  m_tree->Branch("el_ancestorTruthType",	       &m_electron_ancestorTruthType);
  m_tree->Branch("el_ancestorTruthPdgId",	       &m_electron_ancestorTruthPdgId);
  m_tree->Branch("el_ancestorTruthOrigin",	       &m_electron_ancestorTruthOrigin);
  m_tree->Branch("el_ancestorTruthStatus",	       &m_electron_ancestorTruthStatus);
  // electron TAG variables
  m_tree->Branch("el_tag_pt",	                       &m_electron_tag_pt);
  m_tree->Branch("el_tag_eta",	                       &m_electron_tag_eta);
  m_tree->Branch("el_tag_caloCluster_eta",	       &m_electron_tag_caloCluster_eta);
  m_tree->Branch("el_tag_LHLoose",                     &m_electron_tag_LHLoose);
  m_tree->Branch("el_tag_LHMedium",                    &m_electron_tag_LHMedium);
  m_tree->Branch("el_tag_LHTight",                     &m_electron_tag_LHTight);
  m_tree->Branch("el_tag_IsEMLoose",                   &m_electron_tag_IsEMLoose);
  m_tree->Branch("el_tag_IsEMMedium",                  &m_electron_tag_IsEMMedium);
  m_tree->Branch("el_tag_IsEMTight",                   &m_electron_tag_IsEMTight);
  m_tree->Branch("el_tag_ptvarcone20",                 &m_electron_tag_ptvarcone20);
  m_tree->Branch("el_tag_topoetcone20",                &m_electron_tag_topoetcone20);
  m_tree->Branch("el_tag_isIsolated_Loose",	       &m_electron_tag_isIsolated_Loose);
  m_tree->Branch("el_tag_isIsolated_FixedCutTight",    &m_electron_tag_isIsolated_FixedCutTight);
  m_tree->Branch("el_tag_isTightSelected",             &m_electron_tag_isTight);
  m_tree->Branch("el_tag_isMediumSelected",            &m_electron_tag_isMedium);
  m_tree->Branch("el_tag_isTruthMatchedToElectron",    &m_electron_tag_isTruthMatched);
  m_tree->Branch("el_tag_isChFlip",	               &m_electron_tag_isChFlip);
  m_tree->Branch("el_tag_isBrem",	               &m_electron_tag_isBrem);
  m_tree->Branch("el_tag_truthPdgId",                  &m_electron_tag_truthPdgId);
  m_tree->Branch("el_tag_truthType",                   &m_electron_tag_truthType);
  m_tree->Branch("el_tag_truthOrigin",                 &m_electron_tag_truthOrigin);
  m_tree->Branch("el_tag_truthStatus",                 &m_electron_tag_truthStatus);
  m_tree->Branch("el_tag_ancestorTruthType",	       &m_electron_tag_ancestorTruthType);
  m_tree->Branch("el_tag_ancestorTruthPdgId",	       &m_electron_tag_ancestorTruthPdgId);
  m_tree->Branch("el_tag_ancestorTruthOrigin",         &m_electron_tag_ancestorTruthOrigin);
  m_tree->Branch("el_tag_ancestorTruthStatus",         &m_electron_tag_ancestorTruthStatus);
  // electron PROBE variables
  m_tree->Branch("el_probe_pt",	                       &m_electron_probe_pt);
  m_tree->Branch("el_probe_eta",	               &m_electron_probe_eta);
  m_tree->Branch("el_probe_caloCluster_eta",	       &m_electron_probe_caloCluster_eta);
  m_tree->Branch("el_probe_LHLoose",                   &m_electron_probe_LHLoose);
  m_tree->Branch("el_probe_LHMedium",                  &m_electron_probe_LHMedium);
  m_tree->Branch("el_probe_LHTight",                   &m_electron_probe_LHTight);
  m_tree->Branch("el_probe_IsEMLoose",                 &m_electron_probe_IsEMLoose);
  m_tree->Branch("el_probe_IsEMMedium",                &m_electron_probe_IsEMMedium);
  m_tree->Branch("el_probe_IsEMTight",                 &m_electron_probe_IsEMTight);
  m_tree->Branch("el_probe_ptvarcone20",               &m_electron_probe_ptvarcone20);
  m_tree->Branch("el_probe_topoetcone20",              &m_electron_probe_topoetcone20);
  m_tree->Branch("el_probe_isIsolated_Loose",	       &m_electron_probe_isIsolated_Loose);
  m_tree->Branch("el_probe_isIsolated_FixedCutTight",  &m_electron_probe_isIsolated_FixedCutTight);
  m_tree->Branch("el_probe_isTightSelected",	       &m_electron_probe_isTight);
  m_tree->Branch("el_probe_isMediumSelected",          &m_electron_probe_isMedium);
  m_tree->Branch("el_probe_isTruthMatchedToElectron",  &m_electron_probe_isTruthMatched);
  m_tree->Branch("el_probe_isChFlip",	               &m_electron_probe_isChFlip);
  m_tree->Branch("el_probe_isBrem",	               &m_electron_probe_isBrem);
  m_tree->Branch("el_probe_truthPdgId",                &m_electron_probe_truthPdgId);
  m_tree->Branch("el_probe_truthType",                 &m_electron_probe_truthType);
  m_tree->Branch("el_probe_truthOrigin",               &m_electron_probe_truthOrigin);
  m_tree->Branch("el_probe_truthStatus",               &m_electron_probe_truthStatus);
  m_tree->Branch("el_probe_ancestorTruthType",	       &m_electron_probe_ancestorTruthType);
  m_tree->Branch("el_probe_ancestorTruthPdgId",	       &m_electron_probe_ancestorTruthPdgId);
  m_tree->Branch("el_probe_ancestorTruthOrigin",       &m_electron_probe_ancestorTruthOrigin);
  m_tree->Branch("el_probe_ancestorTruthStatus",       &m_electron_probe_ancestorTruthStatus);
}


void HTopMultilepTree::AddLeptons()
{
  // always
  m_tree->Branch("nlep",   &m_nlep, "nlep/I");

  m_tree->Branch("lep_pt",          		       &m_lepton_pt);
  m_tree->Branch("lep_phi",         		       &m_lepton_phi);
  m_tree->Branch("lep_eta",         		       &m_lepton_eta);
  m_tree->Branch("lep_m",           		       &m_lepton_m);
  m_tree->Branch("lep_charge",      		       &m_lepton_charge);
  m_tree->Branch("lep_flavour",     		       &m_lepton_flavour);
  m_tree->Branch("lep_isTrigMatched",                  &m_lepton_isTrigMatched);
  m_tree->Branch("lep_isTightSelected",                &m_lepton_isTight);
  m_tree->Branch("lep_isMediumSelected",               &m_lepton_isMedium);
  m_tree->Branch("lep_isTag",                          &m_lepton_isTag);
  m_tree->Branch("lep_isOS",	                       &m_lepton_isOS);
  m_tree->Branch("lep_isClosestSS",                    &m_lepton_isClosestSS);
  m_tree->Branch("lep_isTruthMatchedToLepton",         &m_lepton_isTruthMatched);
  m_tree->Branch("lep_isChFlip",                       &m_lepton_isChFlip);
  m_tree->Branch("lep_isBrem",                         &m_lepton_isBrem);
  m_tree->Branch("lep_truthType",                      &m_lepton_truthType);
  m_tree->Branch("lep_truthPdgId",                     &m_lepton_truthPdgId);
  m_tree->Branch("lep_truthOrigin",                    &m_lepton_truthOrigin);
  m_tree->Branch("lep_truthStatus",                    &m_lepton_truthStatus);
  m_tree->Branch("lep_ancestorTruthType",              &m_lepton_ancestorTruthType);
  m_tree->Branch("lep_ancestorTruthPdgId",             &m_lepton_ancestorTruthPdgId);
  m_tree->Branch("lep_ancestorTruthOrigin",            &m_lepton_ancestorTruthOrigin);
  m_tree->Branch("lep_ancestorTruthStatus",            &m_lepton_ancestorTruthStatus);
  // lepton TAG variables
  m_tree->Branch("lep_tag_pt",	                       &m_lepton_tag_pt);
  m_tree->Branch("lep_tag_eta",	                       &m_lepton_tag_eta);
  m_tree->Branch("lep_tag_flavour",	               &m_lepton_tag_flavour);
  m_tree->Branch("lep_tag_charge",	               &m_lepton_tag_charge);
  m_tree->Branch("lep_tag_isTrigMatched",              &m_lepton_tag_isTrigMatched);
  m_tree->Branch("lep_tag_isTightSelected",            &m_lepton_tag_isTight);
  m_tree->Branch("lep_tag_isMediumSelected",           &m_lepton_tag_isMedium);
  m_tree->Branch("lep_tag_isTruthMatchedToLepton",     &m_lepton_tag_isTruthMatched);
  m_tree->Branch("lep_tag_isChFlip",	               &m_lepton_tag_isChFlip);
  m_tree->Branch("lep_tag_isBrem",	               &m_lepton_tag_isBrem);
  m_tree->Branch("lep_tag_truthType",                  &m_lepton_tag_truthType);
  m_tree->Branch("lep_tag_truthPdgId",                 &m_lepton_tag_truthPdgId);
  m_tree->Branch("lep_tag_truthOrigin",                &m_lepton_tag_truthOrigin);
  m_tree->Branch("lep_tag_truthStatus",                &m_lepton_tag_truthStatus);
  m_tree->Branch("lep_tag_ancestorTruthType",	       &m_lepton_tag_ancestorTruthType);
  m_tree->Branch("lep_tag_ancestorTruthPdgId",         &m_lepton_tag_ancestorTruthPdgId);
  m_tree->Branch("lep_tag_ancestorTruthOrigin",        &m_lepton_tag_ancestorTruthOrigin);
  m_tree->Branch("lep_tag_ancestorTruthStatus",        &m_lepton_tag_ancestorTruthStatus);
  // lepton PROBE variables
  m_tree->Branch("lep_probe_pt",	               &m_lepton_probe_pt);
  m_tree->Branch("lep_probe_eta",	               &m_lepton_probe_eta);
  m_tree->Branch("lep_probe_flavour",	               &m_lepton_probe_flavour);
  m_tree->Branch("lep_probe_charge",	               &m_lepton_probe_charge);
  m_tree->Branch("lep_probe_isTrigMatched",            &m_lepton_probe_isTrigMatched);
  m_tree->Branch("lep_probe_isTightSelected",          &m_lepton_probe_isTight);
  m_tree->Branch("lep_probe_isMediumSelected",         &m_lepton_probe_isMedium);
  m_tree->Branch("lep_probe_isTruthMatchedToLepton",   &m_lepton_probe_isTruthMatched);
  m_tree->Branch("lep_probe_isChFlip",	               &m_lepton_probe_isChFlip);
  m_tree->Branch("lep_probe_isBrem",	               &m_lepton_probe_isBrem);
  m_tree->Branch("lep_probe_truthType",                &m_lepton_probe_truthType);
  m_tree->Branch("lep_probe_truthPdgId",               &m_lepton_probe_truthPdgId);
  m_tree->Branch("lep_probe_truthOrigin",              &m_lepton_probe_truthOrigin);
  m_tree->Branch("lep_probe_truthStatus",              &m_lepton_probe_truthStatus);
  m_tree->Branch("lep_probe_ancestorTruthType",	       &m_lepton_probe_ancestorTruthType);
  m_tree->Branch("lep_probe_ancestorTruthPdgId",       &m_lepton_probe_ancestorTruthPdgId);
  m_tree->Branch("lep_probe_ancestorTruthOrigin",      &m_lepton_probe_ancestorTruthOrigin);
  m_tree->Branch("lep_probe_ancestorTruthStatus",      &m_lepton_probe_ancestorTruthStatus);
}

void HTopMultilepTree::AddTausUser(const std::string detailStrUser)
{

  // to get rid of warning when detailString is not used
  if ( m_debug ) { Info("AddTausUser()", "Adding branches w/ detail: %s", detailStrUser.c_str()); }

  m_tree->Branch("tau_isBDTTight",	                &m_tau_isBDTTight );
}

/*
void HTopMultilepTree::AddMETUser(const std::string detailStrUser)
{

  // to get rid of warning when detailString is not used
  if ( m_debug ) { Info("AddMETUser()", "Adding branches w/ detail: %s", detailStrUser.c_str()); }

}
*/

void HTopMultilepTree::ClearEventUser()
{
  m_MMWeight.clear();
  m_FFWeight.clear();
  if ( m_isMC ) {
    m_weight_lepton_trig_HTop.clear();
    m_weight_lepton_trig_HTop.clear();
    m_weight_lepton_reco_HTop.clear();
    m_weight_lepton_iso_HTop.clear();
    m_weight_lepton_ID_HTop.clear();
    m_weight_lepton_TTVA_HTop.clear();
  }
}

/*
void HTopMultilepTree::ClearTriggerUser() {}
*/

void HTopMultilepTree::ClearMuonsUser()
{
  // muon variables
  m_muon_isTight.clear();
  m_muon_isMedium.clear();
  m_muon_isOS.clear();
  m_muon_isClosestSS.clear();
  m_muon_isTag.clear();
  m_muon_isTruthMatched.clear();
  m_muon_isChFlip.clear();
  m_muon_isBrem.clear();
  m_muon_truthType.clear();
  m_muon_truthPdgId.clear();
  m_muon_truthOrigin.clear();
  m_muon_truthStatus.clear();
  m_muon_ancestorTruthType.clear();
  m_muon_ancestorTruthPdgId.clear();
  m_muon_ancestorTruthOrigin.clear();
  m_muon_ancestorTruthStatus.clear();
  m_muon_tag_pt.clear();
  m_muon_tag_eta.clear();
  m_muon_tag_trkd0sig.clear();
  m_muon_tag_ptvarcone30.clear();
  m_muon_tag_topoetcone20.clear();
  m_muon_tag_isIsolated_Loose.clear();
  m_muon_tag_isIsolated_FixedCutTightTrackOnly.clear();
  m_muon_tag_isTight.clear();
  m_muon_tag_isMedium.clear();
  m_muon_tag_isTruthMatched.clear();
  m_muon_tag_isChFlip.clear();
  m_muon_tag_isBrem.clear();
  m_muon_tag_truthType.clear();
  m_muon_tag_truthPdgId.clear();
  m_muon_tag_truthOrigin.clear();
  m_muon_tag_truthStatus.clear();
  m_muon_tag_ancestorTruthType.clear();
  m_muon_tag_ancestorTruthPdgId.clear();
  m_muon_tag_ancestorTruthOrigin.clear();
  m_muon_tag_ancestorTruthStatus.clear();
  m_muon_probe_pt.clear();
  m_muon_probe_eta.clear();
  m_muon_probe_trkd0sig.clear();
  m_muon_probe_ptvarcone30.clear();
  m_muon_probe_topoetcone20.clear();
  m_muon_probe_isIsolated_Loose.clear();
  m_muon_probe_isIsolated_FixedCutTightTrackOnly.clear();
  m_muon_probe_isTight.clear();
  m_muon_probe_isMedium.clear();
  m_muon_probe_isTruthMatched.clear();
  m_muon_probe_isChFlip.clear();
  m_muon_probe_isBrem.clear();
  m_muon_probe_truthType.clear();
  m_muon_probe_truthPdgId.clear();
  m_muon_probe_truthOrigin.clear();
  m_muon_probe_truthStatus.clear();
  m_muon_probe_ancestorTruthType.clear();
  m_muon_probe_ancestorTruthPdgId.clear();
  m_muon_probe_ancestorTruthOrigin.clear();
  m_muon_probe_ancestorTruthStatus.clear();
}

void HTopMultilepTree::ClearElectronsUser()
{
  // electron variables
  m_electron_crack.clear();
  m_electron_isTight.clear();
  m_electron_isMedium.clear();
  m_electron_isOS.clear();
  m_electron_isClosestSS.clear();
  m_electron_isTag.clear();
  m_electron_isTruthMatched.clear();
  m_electron_isChFlip.clear();
  m_electron_isBrem.clear();
  m_electron_truthType.clear();
  m_electron_truthPdgId.clear();
  m_electron_truthOrigin.clear();
  m_electron_truthStatus.clear();
  m_electron_ancestorTruthType.clear();
  m_electron_ancestorTruthPdgId.clear();
  m_electron_ancestorTruthOrigin.clear();
  m_electron_ancestorTruthStatus.clear();
  m_electron_tag_pt.clear();
  m_electron_tag_eta.clear();
  m_electron_tag_caloCluster_eta.clear();
  m_electron_tag_LHLoose.clear();
  m_electron_tag_LHMedium.clear();
  m_electron_tag_LHTight.clear();
  m_electron_tag_IsEMLoose.clear();
  m_electron_tag_IsEMMedium.clear();
  m_electron_tag_IsEMTight.clear();
  m_electron_tag_ptvarcone20.clear();
  m_electron_tag_topoetcone20.clear();
  m_electron_tag_isIsolated_Loose.clear();
  m_electron_tag_isIsolated_FixedCutTight.clear();
  m_electron_tag_isTight.clear();
  m_electron_tag_isMedium.clear();
  m_electron_tag_isTruthMatched.clear();
  m_electron_tag_isChFlip.clear();
  m_electron_tag_isBrem.clear();
  m_electron_tag_truthType.clear();
  m_electron_tag_truthPdgId.clear();
  m_electron_tag_truthOrigin.clear();
  m_electron_tag_truthStatus.clear();
  m_electron_tag_ancestorTruthType.clear();
  m_electron_tag_ancestorTruthPdgId.clear();
  m_electron_tag_ancestorTruthOrigin.clear();
  m_electron_tag_ancestorTruthStatus.clear();
  m_electron_probe_pt.clear();
  m_electron_probe_eta.clear();
  m_electron_probe_caloCluster_eta.clear();
  m_electron_probe_LHLoose.clear();
  m_electron_probe_LHMedium.clear();
  m_electron_probe_LHTight.clear();
  m_electron_probe_IsEMLoose.clear();
  m_electron_probe_IsEMMedium.clear();
  m_electron_probe_IsEMTight.clear();
  m_electron_probe_ptvarcone20.clear();
  m_electron_probe_topoetcone20.clear();
  m_electron_probe_isIsolated_Loose.clear();
  m_electron_probe_isIsolated_FixedCutTight.clear();
  m_electron_probe_isTight.clear();
  m_electron_probe_isMedium.clear();
  m_electron_probe_isTruthMatched.clear();
  m_electron_probe_isChFlip.clear();
  m_electron_probe_isBrem.clear();
  m_electron_probe_truthType.clear();
  m_electron_probe_truthPdgId.clear();
  m_electron_probe_truthOrigin.clear();
  m_electron_probe_truthStatus.clear();
  m_electron_probe_ancestorTruthType.clear();
  m_electron_probe_ancestorTruthPdgId.clear();
  m_electron_probe_ancestorTruthOrigin.clear();
  m_electron_probe_ancestorTruthStatus.clear();
}

void HTopMultilepTree::ClearJetsUser( const std::string jetName )
{
  if ( m_debug ) { Info("ClearJetsUser()", "Clearing jet branches - Jet name: %s", jetName.c_str()); }

  // jet variables
  m_jet_m.clear();
}


void HTopMultilepTree::ClearLeptons()
{
  m_nlep = 0;
  m_lepton_pt.clear();
  m_lepton_phi.clear();
  m_lepton_eta.clear();
  m_lepton_m.clear();
  m_lepton_charge.clear();
  m_lepton_flavour.clear();
  m_lepton_isTrigMatched.clear();
  m_lepton_isMedium.clear();
  m_lepton_isTight.clear();
  m_lepton_isOS.clear();
  m_lepton_isClosestSS.clear();
  m_lepton_isTag.clear();
  m_lepton_isTruthMatched.clear();
  m_lepton_isChFlip.clear();
  m_lepton_isBrem.clear();
  m_lepton_truthType.clear();
  m_lepton_truthPdgId.clear();
  m_lepton_truthOrigin.clear();
  m_lepton_truthStatus.clear();
  m_lepton_ancestorTruthType.clear();
  m_lepton_ancestorTruthPdgId.clear();
  m_lepton_ancestorTruthOrigin.clear();
  m_lepton_ancestorTruthStatus.clear();
  m_lepton_tag_pt.clear();
  m_lepton_tag_eta.clear();
  m_lepton_tag_flavour.clear();
  m_lepton_tag_charge.clear();
  m_lepton_tag_isTrigMatched.clear();
  m_lepton_tag_isMedium.clear();
  m_lepton_tag_isTight.clear();
  m_lepton_tag_isTruthMatched.clear();
  m_lepton_tag_isChFlip.clear();
  m_lepton_tag_isBrem.clear();
  m_lepton_tag_truthType.clear();
  m_lepton_tag_truthPdgId.clear();
  m_lepton_tag_truthOrigin.clear();
  m_lepton_tag_truthStatus.clear();
  m_lepton_tag_ancestorTruthType.clear();
  m_lepton_tag_ancestorTruthPdgId.clear();
  m_lepton_tag_ancestorTruthOrigin.clear();
  m_lepton_tag_ancestorTruthStatus.clear();
  m_lepton_probe_pt.clear();
  m_lepton_probe_eta.clear();
  m_lepton_probe_flavour.clear();
  m_lepton_probe_charge.clear();
  m_lepton_probe_isTrigMatched.clear();
  m_lepton_probe_isMedium.clear();
  m_lepton_probe_isTight.clear();
  m_lepton_probe_isTruthMatched.clear();
  m_lepton_probe_isChFlip.clear();
  m_lepton_probe_isBrem.clear();
  m_lepton_probe_truthType.clear();
  m_lepton_probe_truthPdgId.clear();
  m_lepton_probe_truthOrigin.clear();
  m_lepton_probe_truthStatus.clear();
  m_lepton_probe_ancestorTruthType.clear();
  m_lepton_probe_ancestorTruthPdgId.clear();
  m_lepton_probe_ancestorTruthOrigin.clear();
  m_lepton_probe_ancestorTruthStatus.clear();
}

void HTopMultilepTree::ClearTausUser()
{
  m_tau_isBDTTight.clear();
}

/*
void HTopMultilepTree::ClearMETUser() {}
*/

void HTopMultilepTree::FillEventUser( const xAOD::EventInfo* eventInfo )
{

  m_is_mc                =  ( eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) );
  m_ystar                =  ( eventInfo->isAvailable< float >( "ystar" ) )                      ?  eventInfo->auxdecor< float >( "ystar" )                       :   999.0;
  m_categoryFlag         =  ( eventInfo->isAvailable< unsigned int >( "categoryFlag" ) )        ?  eventInfo->auxdecor< unsigned int >( "categoryFlag" )         :   -1;
  m_isSS01   	         =  ( eventInfo->isAvailable< char >( "isSS01" ) )		        ?  eventInfo->auxdecor< char >( "isSS01" )		         :   -1;
  m_isSS12   	         =  ( eventInfo->isAvailable< char >( "isSS12" ) )		        ?  eventInfo->auxdecor< char >( "isSS12" )		     	 :   -1;
  m_MMWeight 	         =  ( eventInfo->isAvailable< std::vector<double> >( "MMWeight" ) )     ?  eventInfo->auxdecor< std::vector<double> >( "MMWeight" )  	 :  std::vector<double>( 5, 1.0 );
  m_FFWeight 	         =  ( eventInfo->isAvailable< std::vector<double> >( "FFWeight" ) )     ?  eventInfo->auxdecor< std::vector<double> >( "FFWeight" )  	 :  std::vector<double>( 5, 1.0 );
  m_mll01    	         =  ( eventInfo->isAvailable< float >( "mll01" ) )		        ?  ( eventInfo->auxdecor< float >( "mll01" ) / m_units )    	 : -1.0;
  m_mll02    	         =  ( eventInfo->isAvailable< float >( "mll02" ) )		        ?  ( eventInfo->auxdecor< float >( "mll02" ) / m_units )    	 : -1.0;
  m_mll12    	         =  ( eventInfo->isAvailable< float >( "mll12" ) )		        ?  ( eventInfo->auxdecor< float >( "mll12" ) / m_units )    	 : -1.0;
  m_mlll012  	         =  ( eventInfo->isAvailable< float >( "mlll012" ) )		        ?  ( eventInfo->auxdecor< float >( "mlll012" ) / m_units )       : -1.0;
  m_mT_lep0MET           =  ( eventInfo->isAvailable< float >( "mT_lep0MET" ) ) 	        ?  ( eventInfo->auxdecor< float >( "mT_lep0MET" ) / m_units )	 : -1.0;
  m_mT_lep1MET           =  ( eventInfo->isAvailable< float >( "mT_lep1MET" ) ) 	        ?  ( eventInfo->auxdecor< float >( "mT_lep1MET" ) / m_units )	 : -1.0;
  m_isTT   	         =  ( eventInfo->isAvailable< char >( "isTT" ) )		        ?  eventInfo->auxdecor< char >( "isTT" )		   	 :   -1;
  m_isTL   	         =  ( eventInfo->isAvailable< char >( "isTL" ) )		        ?  eventInfo->auxdecor< char >( "isTL" )		   	 :   -1;
  m_isLT   	         =  ( eventInfo->isAvailable< char >( "isLT" ) )		        ?  eventInfo->auxdecor< char >( "isLT" )		   	 :   -1;
  m_isLL   	         =  ( eventInfo->isAvailable< char >( "isLL" ) )		        ?  eventInfo->auxdecor< char >( "isLL" )		   	 :   -1;
  m_isTM   	         =  ( eventInfo->isAvailable< char >( "isTM" ) )		        ?  eventInfo->auxdecor< char >( "isTM" )		   	 :   -1;
  m_isMT   	         =  ( eventInfo->isAvailable< char >( "isMT" ) )		        ?  eventInfo->auxdecor< char >( "isMT" )		   	 :   -1;
  m_isMM   	         =  ( eventInfo->isAvailable< char >( "isMM" ) )		        ?  eventInfo->auxdecor< char >( "isMM" )		   	 :   -1;
  m_isTelLmu		 =  ( eventInfo->isAvailable< char >( "isTelLmu" ) )			?  eventInfo->auxdecor< char >( "isTelLmu" )			 :   -1;
  m_isLelTmu		 =  ( eventInfo->isAvailable< char >( "isLelTmu" ) )			?  eventInfo->auxdecor< char >( "isLelTmu" )			 :   -1;
  m_isTmuLel		 =  ( eventInfo->isAvailable< char >( "isTmuLel" ) )			?  eventInfo->auxdecor< char >( "isTmuLel" )			 :   -1;
  m_isLmuTel		 =  ( eventInfo->isAvailable< char >( "isLmuTel" ) )			?  eventInfo->auxdecor< char >( "isLmuTel" )			 :   -1;
  m_isTelMmu		 =  ( eventInfo->isAvailable< char >( "isTelMmu" ) )			?  eventInfo->auxdecor< char >( "isTelMmu" )			 :   -1;
  m_isMelTmu		 =  ( eventInfo->isAvailable< char >( "isMelTmu" ) )			?  eventInfo->auxdecor< char >( "isMelTmu" )			 :   -1;
  m_isTmuMel		 =  ( eventInfo->isAvailable< char >( "isTmuMel" ) )			?  eventInfo->auxdecor< char >( "isTmuMel" )			 :   -1;
  m_isMmuTel		 =  ( eventInfo->isAvailable< char >( "isMmuTel" ) )			?  eventInfo->auxdecor< char >( "isMmuTel" )			 :   -1;
  m_isNonTightEvent      =  ( eventInfo->isAvailable< char >( "isNonTightEvent" ) )	        ?  eventInfo->auxdecor< char >( "isNonTightEvent" )	   	 :  -1;
  m_isProbeElEvent       =  ( eventInfo->isAvailable< char >( "isProbeElEvent" ) )	        ?  eventInfo->auxdecor< char >( "isProbeElEvent" )	   	 :  -1;
  m_isProbeMuEvent       =  ( eventInfo->isAvailable< char >( "isProbeMuEvent" ) )	        ?  eventInfo->auxdecor< char >( "isProbeMuEvent" )	   	 :  -1;

  if ( m_isMC ) {
    std::vector<float> junk(1,1.0);
    static SG::AuxElement::ConstAccessor< std::vector<float> >  accLepTrigSF_GLOBAL("lepTrigEffSF_GLOBAL_HTop");
    static SG::AuxElement::ConstAccessor< std::vector<float> >  accLepRecoSF_GLOBAL("lepRecoEffSF_GLOBAL_HTop");
    static SG::AuxElement::ConstAccessor< std::vector<float> >  accLepIsoSF_GLOBAL("lepIsoEffSF_GLOBAL_HTop");
    static SG::AuxElement::ConstAccessor< std::vector<float> >  accLepIDSF_GLOBAL("lepIDEffSF_GLOBAL_HTop");
    static SG::AuxElement::ConstAccessor< std::vector<float> >  accLepTTVASF_GLOBAL("lepTTVAEffSF_GLOBAL_HTop");

    if ( accLepTrigSF_GLOBAL.isAvailable( *eventInfo ) ) { m_weight_lepton_trig_HTop = accLepTrigSF_GLOBAL( *eventInfo ); } else { m_weight_lepton_trig_HTop = junk; }
    if ( accLepRecoSF_GLOBAL.isAvailable( *eventInfo ) ) { m_weight_lepton_reco_HTop = accLepRecoSF_GLOBAL( *eventInfo ); } else { m_weight_lepton_reco_HTop = junk; }
    if ( accLepIsoSF_GLOBAL.isAvailable( *eventInfo ) )  { m_weight_lepton_iso_HTop = accLepIsoSF_GLOBAL( *eventInfo ); } else { m_weight_lepton_iso_HTop = junk; }
    if ( accLepIDSF_GLOBAL.isAvailable( *eventInfo ) )   { m_weight_lepton_ID_HTop = accLepIDSF_GLOBAL( *eventInfo ); } else { m_weight_lepton_ID_HTop = junk; }
    if ( accLepTTVASF_GLOBAL.isAvailable( *eventInfo ) )   { m_weight_lepton_TTVA_HTop = accLepTTVASF_GLOBAL( *eventInfo ); } else { m_weight_lepton_TTVA_HTop = junk; }

  }
}

/*
void HTopMultilepTree::FillTriggerUser( const xAOD::EventInfo* eventInfo ) { }
*/

void HTopMultilepTree::FillJetsUser( const xAOD::Jet* jet, const std::string jetName )
{
  if ( m_debug ) { Info("FillJetsUser()", "Filling jets - Jet name: %s", jetName.c_str()); }

  m_jet_m.push_back( jet->m() );
}

void HTopMultilepTree::FillMuonsUser( const xAOD::Muon* muon )
{

  // access this info only to fill tag/probe branches
  static SG::AuxElement::Accessor<float> d0SigAcc ("d0sig");
  float d0_significance =  ( d0SigAcc.isAvailable( *muon ) ) ? fabs( d0SigAcc( *muon ) ) : -9999.0;
  static SG::AuxElement::Accessor<char> isIsoLooseAcc ("isIsolated_Loose");
  static SG::AuxElement::Accessor<char> isIsoFixedCutTightTrackOnlyAcc ("isIsolated_FixedCutTightTrackOnly");

  static SG::AuxElement::Accessor< char > isTightAcc("isTight");
  static SG::AuxElement::Accessor< char > isMediumAcc("isMedium");
  static SG::AuxElement::Accessor< char > isOSlepAcc("isOSlep");
  static SG::AuxElement::Accessor< char > isClosestSSlepAcc("isClosestSSlep");
  static SG::AuxElement::Accessor< char > isTagAcc("isTag");

  static SG::AuxElement::Accessor< char > isTruthMatchedAcc("isTruthMatched");
  static SG::AuxElement::Accessor< char > isChFlipAcc("isChFlip");
  static SG::AuxElement::Accessor< char > isBremAcc("isBrem");
  static SG::AuxElement::ConstAccessor< int >  truthTypeAcc("truthType");
  static SG::AuxElement::Accessor< int >  truthPdgIdAcc("truthPdgId");
  static SG::AuxElement::ConstAccessor< int >  truthOriginAcc("truthOrigin");
  static SG::AuxElement::Accessor< int >  truthStatusAcc("truthStatus");
  static SG::AuxElement::Accessor< int >  ancestorTruthTypeAcc("ancestorTruthType");
  static SG::AuxElement::Accessor< int >  ancestorTruthPdgIdAcc("ancestorTruthPdgId");
  static SG::AuxElement::Accessor< int >  ancestorTruthOriginAcc("ancestorTruthOrigin");
  static SG::AuxElement::Accessor< int >  ancestorTruthStatusAcc("ancestorTruthStatus");

  if (  isTightAcc.isAvailable( *muon ) )         {  m_muon_isTight.push_back( isTightAcc( *muon ) ); }
  else {  m_muon_isTight.push_back(-1); }
  if (  isMediumAcc.isAvailable( *muon ) )         {  m_muon_isMedium.push_back( isMediumAcc( *muon ) ); }
  else {  m_muon_isMedium.push_back(-1); }
  if (  isOSlepAcc.isAvailable( *muon ) )         { m_muon_isOS.push_back( isOSlepAcc( *muon ) ); }
  else  { m_muon_isOS.push_back(-1); }
  if (  isClosestSSlepAcc.isAvailable( *muon ) )  { m_muon_isClosestSS.push_back( isClosestSSlepAcc( *muon ) ); }
  else  { m_muon_isClosestSS.push_back(-1); }
  if ( isTruthMatchedAcc.isAvailable( *muon ) )   { m_muon_isTruthMatched.push_back( isTruthMatchedAcc( *muon ) ); }
  else   { m_muon_isTruthMatched.push_back(-1); }
  if ( isChFlipAcc.isAvailable( *muon ) )         { m_muon_isChFlip.push_back( isChFlipAcc( *muon ) ); }
  else   { m_muon_isChFlip.push_back(-1); }
  if ( isBremAcc.isAvailable( *muon ) )           { m_muon_isBrem.push_back( isBremAcc( *muon ) ); }
  else   { m_muon_isBrem.push_back(-1); }
   if ( truthPdgIdAcc.isAvailable( *muon ) )      { m_muon_truthPdgId.push_back( truthPdgIdAcc( *muon ) ); }
  else   { m_muon_truthPdgId.push_back(0); }
  if ( truthTypeAcc.isAvailable( *muon ) )        { m_muon_truthType.push_back( truthTypeAcc( *muon ) ); }
  else   { m_muon_truthType.push_back(-1); }
  if ( truthOriginAcc.isAvailable( *muon ) )      { m_muon_truthOrigin.push_back( truthOriginAcc( *muon ) ); }
  else   { m_muon_truthOrigin.push_back(0); }
  if ( truthStatusAcc.isAvailable( *muon ) )      { m_muon_truthStatus.push_back( truthStatusAcc( *muon ) ); }
  else   { m_muon_truthStatus.push_back(0); }
  if ( ancestorTruthTypeAcc.isAvailable( *muon ) )      { m_muon_ancestorTruthType.push_back( ancestorTruthTypeAcc( *muon ) ); }
  else   { m_muon_ancestorTruthType.push_back(0); }
  if ( ancestorTruthPdgIdAcc.isAvailable( *muon ) )     { m_muon_ancestorTruthPdgId.push_back( ancestorTruthPdgIdAcc( *muon ) ); }
  else   { m_muon_ancestorTruthPdgId.push_back(0); }
  if ( ancestorTruthOriginAcc.isAvailable( *muon ) )    { m_muon_ancestorTruthOrigin.push_back( ancestorTruthOriginAcc( *muon ) ); }
  else   { m_muon_ancestorTruthOrigin.push_back(0); }
  if ( ancestorTruthStatusAcc.isAvailable( *muon ) )    { m_muon_ancestorTruthStatus.push_back( ancestorTruthStatusAcc( *muon ) ); }
  else   { m_muon_ancestorTruthStatus.push_back(0); }

  if ( isTagAcc.isAvailable( *muon ) ) {

    m_muon_isTag.push_back( isTagAcc( *muon ) );

    // fill muon TAG variables
    if ( isTagAcc( *muon ) ){

      m_muon_tag_pt.push_back(  muon->pt() );
      m_muon_tag_eta.push_back(  muon->eta() );
      m_muon_tag_trkd0sig.push_back( d0_significance );
      m_muon_tag_ptvarcone30.push_back( muon->isolation( xAOD::Iso::ptvarcone30 ) );
      m_muon_tag_topoetcone20.push_back( muon->isolation( xAOD::Iso::topoetcone20 ) );
      if ( isIsoLooseAcc.isAvailable( *muon ) )              { m_muon_tag_isIsolated_Loose.push_back( isIsoLooseAcc( *muon ) ); } else { m_muon_tag_isIsolated_Loose.push_back( -1 ); }
      if ( isIsoFixedCutTightTrackOnlyAcc.isAvailable( *muon ) ) { m_muon_tag_isIsolated_FixedCutTightTrackOnly.push_back( isIsoFixedCutTightTrackOnlyAcc( *muon ) ); } else { m_muon_tag_isIsolated_FixedCutTightTrackOnly.push_back( -1 ); }
      if ( isTightAcc.isAvailable( *muon ) )                 { m_muon_tag_isTight.push_back(  isTightAcc( *muon ) ); }
      else { m_muon_tag_isTight .push_back(  -1 ); }
      if ( isMediumAcc.isAvailable( *muon ) )                { m_muon_tag_isMedium.push_back(  isMediumAcc( *muon ) ); }
      else { m_muon_tag_isMedium .push_back(  -1 ); }
      if ( isTruthMatchedAcc.isAvailable( *muon ) )          { m_muon_tag_isTruthMatched.push_back( isTruthMatchedAcc( *muon ) ); }
      else {  m_muon_tag_isTruthMatched .push_back( -1 ); }
      if ( isChFlipAcc.isAvailable( *muon ) )                { m_muon_tag_isChFlip.push_back(  isChFlipAcc( *muon ) ); }
      else { m_muon_tag_isChFlip.push_back( -1 ); }
      if ( isBremAcc.isAvailable( *muon ) )                  { m_muon_tag_isBrem.push_back(  isBremAcc( *muon ) ); }
      else { m_muon_tag_isBrem.push_back( -1 ); }
      if ( truthPdgIdAcc.isAvailable( *muon ) )              { m_muon_tag_truthPdgId.push_back(  truthPdgIdAcc( *muon ) ); }
      else { m_muon_tag_truthPdgId.push_back( -1 ); }
      if ( truthTypeAcc.isAvailable( *muon ) )               { m_muon_tag_truthType.push_back(  truthTypeAcc( *muon ) ); }
      else { m_muon_tag_truthType.push_back( -1 ); }
      if ( truthOriginAcc.isAvailable( *muon ) )             { m_muon_tag_truthOrigin.push_back(  truthOriginAcc( *muon ) ); }
      else { m_muon_tag_truthOrigin.push_back( -1 ); }
      if ( truthStatusAcc.isAvailable( *muon ) )             { m_muon_tag_truthStatus.push_back(  truthStatusAcc( *muon ) ); }
      else { m_muon_tag_truthStatus.push_back( -1 ); }
      if ( ancestorTruthTypeAcc.isAvailable( *muon ) )       { m_muon_tag_ancestorTruthType.push_back( ancestorTruthTypeAcc( *muon ) ); }
      else   { m_muon_tag_ancestorTruthType.push_back(0); }
      if ( ancestorTruthPdgIdAcc.isAvailable( *muon ) )      { m_muon_tag_ancestorTruthPdgId.push_back( ancestorTruthPdgIdAcc( *muon ) ); }
      else   { m_muon_tag_ancestorTruthPdgId.push_back(0); }
      if ( ancestorTruthOriginAcc.isAvailable( *muon ) )     { m_muon_tag_ancestorTruthOrigin.push_back( ancestorTruthOriginAcc( *muon ) ); }
      else   { m_muon_tag_ancestorTruthOrigin.push_back(0); }
      if ( ancestorTruthStatusAcc.isAvailable( *muon ) )     { m_muon_tag_ancestorTruthStatus.push_back( ancestorTruthStatusAcc( *muon ) ); }
      else   { m_muon_tag_ancestorTruthStatus.push_back(0); }

    } else {
    // fill muon PROBE variables

      m_muon_probe_pt.push_back(  muon->pt() );
      m_muon_probe_eta.push_back(  muon->eta() );
      m_muon_probe_trkd0sig.push_back( d0_significance );
      m_muon_probe_ptvarcone30.push_back( muon->isolation( xAOD::Iso::ptvarcone30 ) );
      m_muon_probe_topoetcone20.push_back( muon->isolation( xAOD::Iso::topoetcone20 ) );
      if ( isIsoLooseAcc.isAvailable( *muon ) )          { m_muon_probe_isIsolated_Loose.push_back( isIsoLooseAcc( *muon ) ); } else { m_muon_probe_isIsolated_Loose.push_back( -1 ); }
      if ( isIsoFixedCutTightTrackOnlyAcc.isAvailable( *muon ) )  { m_muon_probe_isIsolated_FixedCutTightTrackOnly.push_back( isIsoFixedCutTightTrackOnlyAcc( *muon ) ); } else { m_muon_probe_isIsolated_FixedCutTightTrackOnly.push_back( -1 ); }
      if ( isTightAcc.isAvailable( *muon ) )                 { m_muon_probe_isTight.push_back(  isTightAcc( *muon ) ); }
      else { m_muon_probe_isTight.push_back(  -1 ); }
      if ( isMediumAcc.isAvailable( *muon ) )                 { m_muon_probe_isMedium.push_back(  isMediumAcc( *muon ) ); }
      else { m_muon_probe_isMedium.push_back(  -1 ); }
      if ( isTruthMatchedAcc.isAvailable( *muon ) )          { m_muon_probe_isTruthMatched.push_back( isTruthMatchedAcc( *muon ) ); }
      else {  m_muon_probe_isTruthMatched.push_back( -1 ); }
      if ( isChFlipAcc.isAvailable( *muon ) )                { m_muon_probe_isChFlip.push_back(  isChFlipAcc( *muon ) ); }
      else { m_muon_probe_isChFlip.push_back( -1 ); }
      if ( isBremAcc.isAvailable( *muon ) )                  { m_muon_probe_isBrem.push_back(  isBremAcc( *muon ) ); }
      else { m_muon_probe_isBrem.push_back( -1 ); }
      if ( truthPdgIdAcc.isAvailable( *muon ) )              { m_muon_probe_truthPdgId.push_back(  truthPdgIdAcc( *muon ) ); }
      else { m_muon_probe_truthPdgId.push_back( -1 ); }
      if ( truthTypeAcc.isAvailable( *muon ) )               { m_muon_probe_truthType.push_back(  truthTypeAcc( *muon ) ); }
      else { m_muon_probe_truthType.push_back( -1 ); }
      if ( truthOriginAcc.isAvailable( *muon ) )             { m_muon_probe_truthOrigin.push_back(  truthOriginAcc( *muon ) ); }
      else { m_muon_probe_truthOrigin.push_back( -1 ); }
      if ( truthStatusAcc.isAvailable( *muon ) )             { m_muon_probe_truthStatus.push_back(  truthStatusAcc( *muon ) ); }
      else { m_muon_probe_truthStatus.push_back( -1 ); }
      if ( ancestorTruthTypeAcc.isAvailable( *muon ) )       { m_muon_probe_ancestorTruthType.push_back( ancestorTruthTypeAcc( *muon ) ); }
      else   { m_muon_probe_ancestorTruthType.push_back(0); }
      if ( ancestorTruthPdgIdAcc.isAvailable( *muon ) )      { m_muon_probe_ancestorTruthPdgId.push_back( ancestorTruthPdgIdAcc( *muon ) ); }
      else   { m_muon_probe_ancestorTruthPdgId.push_back(0); }
      if ( ancestorTruthOriginAcc.isAvailable( *muon ) )     { m_muon_probe_ancestorTruthOrigin.push_back( ancestorTruthOriginAcc( *muon ) ); }
      else   { m_muon_probe_ancestorTruthOrigin.push_back(0); }
      if ( ancestorTruthStatusAcc.isAvailable( *muon ) )     { m_muon_probe_ancestorTruthStatus.push_back( ancestorTruthStatusAcc( *muon ) ); }
      else   { m_muon_probe_ancestorTruthStatus.push_back(0); }

    }

  }
  else   { m_muon_isTag.push_back(-1); }

}

void HTopMultilepTree::FillElectronsUser( const xAOD::Electron* electron )
{

  // access this info only to fill tag/probe branches

  static SG::AuxElement::Accessor<char> LHLooseAcc ("LHLoose");
  static SG::AuxElement::Accessor<char> LHMediumAcc ("LHMedium");
  static SG::AuxElement::Accessor<char> LHTightAcc ("LHTight");
  static SG::AuxElement::Accessor<char> EMLooseAcc ("Loose");
  static SG::AuxElement::Accessor<char> EMMediumAcc ("Medium");
  static SG::AuxElement::Accessor<char> EMTightAcc ("Tight");
  static SG::AuxElement::Accessor<char> isIsoLooseAcc ("isIsolated_Loose");
  static SG::AuxElement::Accessor<char> isIsoFixedCutTightAcc ("isIsolated_FixedCutTight");

  static SG::AuxElement::Accessor< char > isTightAcc("isTight");
  static SG::AuxElement::Accessor< char > isMediumAcc("isMedium");
  static SG::AuxElement::Accessor< char > isOSlepAcc("isOSlep");
  static SG::AuxElement::Accessor< char > isClosestSSlepAcc("isClosestSSlep");
  static SG::AuxElement::Accessor< char > isTagAcc("isTag");

  static SG::AuxElement::Accessor< char > isTruthMatchedAcc("isTruthMatched");
  static SG::AuxElement::Accessor< char > isChFlipAcc("isChFlip");
  static SG::AuxElement::Accessor< char > isBremAcc("isBrem");
  static SG::AuxElement::ConstAccessor< int >  truthTypeAcc("truthType");
  static SG::AuxElement::Accessor< int >  truthPdgIdAcc("truthPdgId");
  static SG::AuxElement::ConstAccessor< int >  truthOriginAcc("truthOrigin");
  static SG::AuxElement::Accessor< int >  truthStatusAcc("truthStatus");
  static SG::AuxElement::Accessor< int >  ancestorTruthTypeAcc("ancestorTruthType");
  static SG::AuxElement::Accessor< int >  ancestorTruthPdgIdAcc("ancestorTruthPdgId");
  static SG::AuxElement::Accessor< int >  ancestorTruthOriginAcc("ancestorTruthOrigin");
  static SG::AuxElement::Accessor< int >  ancestorTruthStatusAcc("ancestorTruthStatus");

  float calo_eta = ( electron->caloCluster() ) ? electron->caloCluster()->etaBE(2) : -999.0;

  if ( electron->caloCluster() ) {
    if ( fabs( calo_eta ) > 1.37 && fabs( calo_eta ) < 1.52 ) {
      m_electron_crack.push_back(1);
    } else {
      m_electron_crack.push_back(0);
    }
  } else {
    m_electron_crack.push_back(-1);
  }

  if (  isTightAcc.isAvailable( *electron ) )          {  m_electron_isTight.push_back( isTightAcc( *electron ) ); }
  else {  m_electron_isTight.push_back(-1); }
  if (  isMediumAcc.isAvailable( *electron ) )         {  m_electron_isMedium.push_back( isMediumAcc( *electron ) ); }
  else {  m_electron_isMedium.push_back(-1); }
  if (  isOSlepAcc.isAvailable( *electron ) )          { m_electron_isOS.push_back( isOSlepAcc( *electron ) ); }
  else  { m_electron_isOS.push_back(-1); }
  if (  isClosestSSlepAcc.isAvailable( *electron ) )   { m_electron_isClosestSS.push_back( isClosestSSlepAcc( *electron ) ); }
  else  { m_electron_isClosestSS.push_back(-1); }
  if ( isTruthMatchedAcc.isAvailable( *electron ) )    { m_electron_isTruthMatched.push_back( isTruthMatchedAcc( *electron ) ); }
  else   { m_electron_isTruthMatched.push_back(-1); }
  if ( isChFlipAcc.isAvailable( *electron ) )          { m_electron_isChFlip.push_back( isChFlipAcc( *electron ) ); }
  else   { m_electron_isChFlip.push_back(-1); }
  if ( isBremAcc.isAvailable( *electron ) )            { m_electron_isBrem.push_back( isBremAcc( *electron ) ); }
  else   { m_electron_isBrem.push_back(-1); }
  if ( truthPdgIdAcc.isAvailable( *electron ) )        { m_electron_truthPdgId.push_back( truthPdgIdAcc( *electron ) ); }
  else   { m_electron_truthPdgId.push_back(0); }
  if ( truthTypeAcc.isAvailable( *electron ) )         { m_electron_truthType.push_back( truthTypeAcc( *electron ) ); }
  else   { m_electron_truthType.push_back(-1); }
  if ( truthOriginAcc.isAvailable( *electron ) )       { m_electron_truthOrigin.push_back( truthOriginAcc( *electron ) ); }
  else   { m_electron_truthOrigin.push_back(0); }
  if ( truthStatusAcc.isAvailable( *electron ) )       { m_electron_truthStatus.push_back( truthStatusAcc( *electron ) ); }
  else   { m_electron_truthStatus.push_back(0); }
  if ( ancestorTruthTypeAcc.isAvailable( *electron ) )      { m_electron_ancestorTruthType.push_back( ancestorTruthTypeAcc( *electron ) ); }
  else   { m_electron_ancestorTruthType.push_back(0); }
  if ( ancestorTruthPdgIdAcc.isAvailable( *electron ) )     { m_electron_ancestorTruthPdgId.push_back( ancestorTruthPdgIdAcc( *electron ) ); }
  else   { m_electron_ancestorTruthPdgId.push_back(0); }
  if ( ancestorTruthOriginAcc.isAvailable( *electron ) )    { m_electron_ancestorTruthOrigin.push_back( ancestorTruthOriginAcc( *electron ) ); }
  else   { m_electron_ancestorTruthOrigin.push_back(0); }
  if ( ancestorTruthStatusAcc.isAvailable( *electron ) )    { m_electron_ancestorTruthStatus.push_back( ancestorTruthStatusAcc( *electron ) ); }
  else   { m_electron_ancestorTruthStatus.push_back(0); }

  if ( isTagAcc.isAvailable( *electron ) ) {

    m_electron_isTag.push_back( isTagAcc( *electron ) );

    // fill electron TAG variables
    if ( isTagAcc( *electron ) ){

      m_electron_tag_pt.push_back(  electron->pt() );
      m_electron_tag_eta.push_back(  electron->eta() );
      m_electron_tag_caloCluster_eta.push_back( calo_eta );

      if ( LHLooseAcc.isAvailable( *electron ) )     { m_electron_tag_LHLoose.push_back( LHLooseAcc( *electron ) );         } else { m_electron_tag_LHLoose.push_back( -1 ); }
      if ( LHMediumAcc.isAvailable( *electron ) )    { m_electron_tag_LHMedium.push_back( LHMediumAcc( *electron ) );       } else { m_electron_tag_LHMedium.push_back( -1 ); }
      if ( LHTightAcc.isAvailable( *electron ) )     { m_electron_tag_LHTight.push_back( LHTightAcc( *electron ) );         } else { m_electron_tag_LHTight.push_back( -1 ); }
      if ( EMLooseAcc.isAvailable( *electron ) )         { m_electron_tag_IsEMLoose.push_back( EMLooseAcc( *electron ) );   } else { m_electron_tag_IsEMLoose.push_back( -1 ); }
      if ( EMMediumAcc.isAvailable( *electron ) )        { m_electron_tag_IsEMMedium.push_back( EMMediumAcc( *electron ) ); } else { m_electron_tag_IsEMMedium.push_back( -1 ); }
      if ( EMTightAcc.isAvailable( *electron ) )         { m_electron_tag_IsEMTight.push_back( EMTightAcc( *electron ) );   } else { m_electron_tag_IsEMTight.push_back( -1 ); }
      m_electron_tag_ptvarcone20.push_back( electron->isolation( xAOD::Iso::ptvarcone20 ) );
      m_electron_tag_topoetcone20.push_back( electron->isolation( xAOD::Iso::topoetcone20 ) );
      if ( isIsoLooseAcc.isAvailable( *electron ) )              { m_electron_tag_isIsolated_Loose.push_back( isIsoLooseAcc( *electron ) ); } else { m_electron_tag_isIsolated_Loose.push_back( -1 ); }
      if ( isIsoFixedCutTightAcc.isAvailable( *electron ) )      { m_electron_tag_isIsolated_FixedCutTight.push_back( isIsoFixedCutTightAcc( *electron ) ); } else { m_electron_tag_isIsolated_FixedCutTight.push_back( -1 ); }
      if ( isTightAcc.isAvailable( *electron ) )                 { m_electron_tag_isTight.push_back(  isTightAcc( *electron ) ); }
      else { m_electron_tag_isMedium .push_back(  -1 ); }
      if ( isMediumAcc.isAvailable( *electron ) )                 { m_electron_tag_isMedium.push_back(  isMediumAcc( *electron ) ); }
      else { m_electron_tag_isMedium .push_back(  -1 ); }
      if ( isTruthMatchedAcc.isAvailable( *electron ) )          { m_electron_tag_isTruthMatched.push_back( isTruthMatchedAcc( *electron ) ); }
      else {  m_electron_tag_isTruthMatched .push_back( -1 ); }
      if ( isChFlipAcc.isAvailable( *electron ) )                { m_electron_tag_isChFlip.push_back(  isChFlipAcc( *electron ) ); }
      else { m_electron_tag_isChFlip.push_back( -1 ); }
      if ( isBremAcc.isAvailable( *electron ) )                  { m_electron_tag_isBrem.push_back(  isBremAcc( *electron ) ); }
      else { m_electron_tag_isBrem.push_back( -1 ); }
      if ( truthPdgIdAcc.isAvailable( *electron ) )              { m_electron_tag_truthPdgId.push_back(  truthPdgIdAcc( *electron ) ); }
      else { m_electron_tag_truthPdgId.push_back( -1 ); }
      if ( truthTypeAcc.isAvailable( *electron ) )               { m_electron_tag_truthType.push_back(  truthTypeAcc( *electron ) ); }
      else { m_electron_tag_truthType.push_back( -1 ); }
      if ( truthOriginAcc.isAvailable( *electron ) )             { m_electron_tag_truthOrigin.push_back(  truthOriginAcc( *electron ) ); }
      else { m_electron_tag_truthOrigin.push_back( -1 ); }
      if ( truthStatusAcc.isAvailable( *electron ) )             { m_electron_tag_truthStatus.push_back(  truthStatusAcc( *electron ) ); }
      else { m_electron_tag_truthStatus.push_back( -1 ); }
      if ( ancestorTruthTypeAcc.isAvailable( *electron ) )       { m_electron_tag_ancestorTruthType.push_back( ancestorTruthTypeAcc( *electron ) ); }
      else   { m_electron_tag_ancestorTruthType.push_back(0); }
      if ( ancestorTruthPdgIdAcc.isAvailable( *electron ) )      { m_electron_tag_ancestorTruthPdgId.push_back( ancestorTruthPdgIdAcc( *electron ) ); }
      else   { m_electron_tag_ancestorTruthPdgId.push_back(0); }
      if ( ancestorTruthOriginAcc.isAvailable( *electron ) )     { m_electron_tag_ancestorTruthOrigin.push_back( ancestorTruthOriginAcc( *electron ) ); }
      else   { m_electron_tag_ancestorTruthOrigin.push_back(0); }
      if ( ancestorTruthStatusAcc.isAvailable( *electron ) )     { m_electron_tag_ancestorTruthStatus.push_back( ancestorTruthStatusAcc( *electron ) ); }
      else   { m_electron_tag_ancestorTruthStatus.push_back(0); }

    } else {
    // fill electron PROBE variables

      m_electron_probe_pt.push_back(  electron->pt() );
      m_electron_probe_eta.push_back(  electron->eta() );
      m_electron_probe_caloCluster_eta.push_back( calo_eta );

      if ( LHLooseAcc.isAvailable( *electron ) )     { m_electron_probe_LHLoose.push_back( LHLooseAcc( *electron ) );         } else { m_electron_probe_LHLoose.push_back( -1 ); }
      if ( LHMediumAcc.isAvailable( *electron ) )    { m_electron_probe_LHMedium.push_back( LHMediumAcc( *electron ) );       } else { m_electron_probe_LHMedium.push_back( -1 ); }
      if ( LHTightAcc.isAvailable( *electron ) )     { m_electron_probe_LHTight.push_back( LHTightAcc( *electron ) );         } else { m_electron_probe_LHTight.push_back( -1 ); }
      if ( EMLooseAcc.isAvailable( *electron ) )         { m_electron_probe_IsEMLoose.push_back( EMLooseAcc( *electron ) );   } else { m_electron_probe_IsEMLoose.push_back( -1 ); }
      if ( EMMediumAcc.isAvailable( *electron ) )        { m_electron_probe_IsEMMedium.push_back( EMMediumAcc( *electron ) ); } else { m_electron_probe_IsEMMedium.push_back( -1 ); }
      if ( EMTightAcc.isAvailable( *electron ) )         { m_electron_probe_IsEMTight.push_back( EMTightAcc( *electron ) );   } else { m_electron_probe_IsEMTight.push_back( -1 ); }
      m_electron_probe_ptvarcone20.push_back( electron->isolation( xAOD::Iso::ptvarcone20 ) );
      m_electron_probe_topoetcone20.push_back( electron->isolation( xAOD::Iso::topoetcone20 ) );
      if ( isIsoLooseAcc.isAvailable( *electron ) )              { m_electron_probe_isIsolated_Loose.push_back( isIsoLooseAcc( *electron ) ); } else { m_electron_probe_isIsolated_Loose.push_back( -1 ); }
      if ( isIsoFixedCutTightAcc.isAvailable( *electron ) )      { m_electron_probe_isIsolated_FixedCutTight.push_back( isIsoFixedCutTightAcc( *electron ) ); } else { m_electron_probe_isIsolated_FixedCutTight.push_back( -1 ); }
      if ( isTightAcc.isAvailable( *electron ) )                 { m_electron_probe_isTight.push_back(  isTightAcc( *electron ) ); }
      else { m_electron_probe_isTight.push_back(  -1 ); }
      if ( isMediumAcc.isAvailable( *electron ) )                { m_electron_probe_isMedium.push_back(  isMediumAcc( *electron ) ); }
      else { m_electron_probe_isMedium.push_back(  -1 ); }
      if ( isTruthMatchedAcc.isAvailable( *electron ) )          { m_electron_probe_isTruthMatched.push_back( isTruthMatchedAcc( *electron ) ); }
      else {  m_electron_probe_isTruthMatched.push_back( -1 ); }
      if ( isChFlipAcc.isAvailable( *electron ) )                { m_electron_probe_isChFlip.push_back(  isChFlipAcc( *electron ) ); }
      else { m_electron_probe_isChFlip.push_back( -1 ); }
      if ( isBremAcc.isAvailable( *electron ) )                  { m_electron_probe_isBrem.push_back(  isBremAcc( *electron ) ); }
      else { m_electron_probe_isBrem.push_back( -1 ); }
      if ( truthPdgIdAcc.isAvailable( *electron ) )              { m_electron_probe_truthPdgId.push_back(  truthPdgIdAcc( *electron ) ); }
      else { m_electron_probe_truthPdgId.push_back( -1 ); }
      if ( truthTypeAcc.isAvailable( *electron ) )               { m_electron_probe_truthType.push_back(  truthTypeAcc( *electron ) ); }
      else { m_electron_probe_truthType.push_back( -1 ); }
      if ( truthOriginAcc.isAvailable( *electron ) )             { m_electron_probe_truthOrigin.push_back(  truthOriginAcc( *electron ) ); }
      else { m_electron_probe_truthOrigin.push_back( -1 ); }
      if ( truthStatusAcc.isAvailable( *electron ) )             { m_electron_probe_truthStatus.push_back(  truthStatusAcc( *electron ) ); }
      else { m_electron_probe_truthStatus.push_back( -1 ); }
      if ( ancestorTruthTypeAcc.isAvailable( *electron ) )       { m_electron_probe_ancestorTruthType.push_back( ancestorTruthTypeAcc( *electron ) ); }
      else   { m_electron_probe_ancestorTruthType.push_back(0); }
      if ( ancestorTruthPdgIdAcc.isAvailable( *electron ) )      { m_electron_probe_ancestorTruthPdgId.push_back( ancestorTruthPdgIdAcc( *electron ) ); }
      else   { m_electron_probe_ancestorTruthPdgId.push_back(0); }
      if ( ancestorTruthOriginAcc.isAvailable( *electron ) )     { m_electron_probe_ancestorTruthOrigin.push_back( ancestorTruthOriginAcc( *electron ) ); }
      else   { m_electron_probe_ancestorTruthOrigin.push_back(0); }
      if ( ancestorTruthStatusAcc.isAvailable( *electron ) )     { m_electron_probe_ancestorTruthStatus.push_back( ancestorTruthStatusAcc( *electron ) ); }
      else   { m_electron_probe_ancestorTruthStatus.push_back(0); }

    }

  }
  else   { m_electron_isTag.push_back(-1); }

}

void HTopMultilepTree::FillTausUser( const xAOD::TauJet* tau )
{

  static SG::AuxElement::Accessor< char > isTauBDTTightAcc("isTauBDTTight");

  if ( isTauBDTTightAcc.isAvailable( *tau ) ) { m_tau_isBDTTight.push_back( isTauBDTTightAcc( *tau ) ); }
  else { m_tau_isBDTTight.push_back(-1); }

}

/*
void HTopMultilepTree::FillMETUser( const xAOD::MissingETContainer* met ) {}
*/

void HTopMultilepTree::FillLeptons( const xAOD::IParticleContainer* leptons )
{

  this->ClearLeptons();

  m_nlep = 0;

  static SG::AuxElement::Accessor< char > isTrigMatchedLepAcc("isTrigMatchedLep");
  static SG::AuxElement::Accessor< char > isTightAcc("isTight");
  static SG::AuxElement::Accessor< char > isMediumAcc("isMedium");
  static SG::AuxElement::Accessor< char > isOSlepAcc("isOSlep");
  static SG::AuxElement::Accessor< char > isClosestSSlepAcc("isClosestSSlep");
  static SG::AuxElement::Accessor< char > isTagAcc("isTag");

  static SG::AuxElement::Accessor< char > isTruthMatchedAcc("isTruthMatched");
  static SG::AuxElement::Accessor< char > isChFlipAcc("isChFlip");
  static SG::AuxElement::Accessor< char > isBremAcc("isBrem");
  static SG::AuxElement::ConstAccessor< int >  truthTypeAcc("truthType");
  static SG::AuxElement::Accessor< int >  truthPdgIdAcc("truthPdgId");
  static SG::AuxElement::ConstAccessor< int >  truthOriginAcc("truthOrigin");
  static SG::AuxElement::Accessor< int >  truthStatusAcc("truthStatus");
  static SG::AuxElement::Accessor< int >  ancestorTruthTypeAcc("ancestorTruthType");
  static SG::AuxElement::Accessor< int >  ancestorTruthPdgIdAcc("ancestorTruthPdgId");
  static SG::AuxElement::Accessor< int >  ancestorTruthOriginAcc("ancestorTruthOrigin");
  static SG::AuxElement::Accessor< int >  ancestorTruthStatusAcc("ancestorTruthStatus");

  for ( auto lep_itr : *(leptons) ) {

      if ( m_debug ) { Info("HTopMultilepTree::FillLeptons()", "Filling lepton w/ pT = %2f", lep_itr->pt() / m_units ); }

      m_lepton_pt.push_back ( lep_itr->pt() / m_units  );
      m_lepton_eta.push_back( lep_itr->eta() );
      m_lepton_phi.push_back( lep_itr->phi() );
      m_lepton_m.push_back  ( lep_itr->m() / m_units  );

      int lep_flav(0);
      float lep_charge( -999.0 );
      if ( lep_itr->type() == xAOD::Type::Electron ) {
        lep_flav = 11;
	const xAOD::TrackParticle* trk = dynamic_cast<const xAOD::Electron*>(lep_itr)->trackParticle();
        if ( trk ) { lep_charge = trk->charge() / fabs(trk->charge()); }
      }
      else if ( lep_itr->type() == xAOD::Type::Muon ) {
        lep_flav = 13;
        const xAOD::TrackParticle* trk = dynamic_cast<const xAOD::Muon*>(lep_itr)->primaryTrackParticle();
        if ( trk ) { lep_charge = trk->charge() / fabs(trk->charge()); }
      }
      m_lepton_flavour.push_back( lep_flav );
      m_lepton_charge.push_back( lep_charge );

      if (  isTrigMatchedLepAcc.isAvailable( *lep_itr ) )  {  m_lepton_isTrigMatched.push_back( isTrigMatchedLepAcc( *lep_itr ) ); }
      else {  m_lepton_isTrigMatched.push_back(-1); }
      if (  isTightAcc.isAvailable( *lep_itr ) )        {  m_lepton_isTight.push_back( isTightAcc( *lep_itr ) ); }
      else {  m_lepton_isTight.push_back(-1); }
      if (  isMediumAcc.isAvailable( *lep_itr ) )        {  m_lepton_isMedium.push_back( isMediumAcc( *lep_itr ) ); }
      else {  m_lepton_isMedium.push_back(-1); }
      if (  isOSlepAcc.isAvailable( *lep_itr ) )        { m_lepton_isOS.push_back( isOSlepAcc( *lep_itr ) ); }
      else  { m_lepton_isOS.push_back(-1); }
      if (  isClosestSSlepAcc.isAvailable( *lep_itr ) ) { m_lepton_isClosestSS.push_back( isClosestSSlepAcc( *lep_itr ) ); }
      else  { m_lepton_isClosestSS.push_back(-1); }
      if ( isTruthMatchedAcc.isAvailable( *lep_itr ) )  { m_lepton_isTruthMatched.push_back( isTruthMatchedAcc( *lep_itr ) ); }
      else   { m_lepton_isTruthMatched.push_back(-1); }
      if ( isChFlipAcc.isAvailable( *lep_itr ) )        { m_lepton_isChFlip.push_back( isChFlipAcc( *lep_itr ) ); }
      else   { m_lepton_isChFlip.push_back(-1); }
      if ( isBremAcc.isAvailable( *lep_itr ) )          { m_lepton_isBrem.push_back( isBremAcc( *lep_itr ) ); }
      else   { m_lepton_isBrem.push_back(-1); }
      if ( truthPdgIdAcc.isAvailable( *lep_itr ) )      { m_lepton_truthPdgId.push_back( truthPdgIdAcc( *lep_itr ) ); }
      else   { m_lepton_truthPdgId.push_back(0); }
      if ( truthTypeAcc.isAvailable( *lep_itr ) )       { m_lepton_truthType.push_back( truthTypeAcc( *lep_itr ) ); }
      else   { m_lepton_truthType.push_back(-1); }
      if ( truthOriginAcc.isAvailable( *lep_itr ) )     { m_lepton_truthOrigin.push_back( truthOriginAcc( *lep_itr ) ); }
      else   { m_lepton_truthOrigin.push_back(0); }
      if ( truthStatusAcc.isAvailable( *lep_itr ) )     { m_lepton_truthStatus.push_back( truthStatusAcc( *lep_itr ) ); }
      else   { m_lepton_truthStatus.push_back(0); }
      if ( ancestorTruthTypeAcc.isAvailable( *lep_itr ) )	 { m_lepton_ancestorTruthType.push_back( ancestorTruthTypeAcc( *lep_itr ) ); }
      else   { m_lepton_ancestorTruthType.push_back(0); }
      if ( ancestorTruthPdgIdAcc.isAvailable( *lep_itr ) )	 { m_lepton_ancestorTruthPdgId.push_back( ancestorTruthPdgIdAcc( *lep_itr ) ); }
      else   { m_lepton_ancestorTruthPdgId.push_back(0); }
      if ( ancestorTruthOriginAcc.isAvailable( *lep_itr ) )	 { m_lepton_ancestorTruthOrigin.push_back( ancestorTruthOriginAcc( *lep_itr ) ); }
      else   { m_lepton_ancestorTruthOrigin.push_back(0); }
      if ( ancestorTruthStatusAcc.isAvailable( *lep_itr ) )	 { m_lepton_ancestorTruthStatus.push_back( ancestorTruthStatusAcc( *lep_itr ) ); }
      else   { m_lepton_ancestorTruthStatus.push_back(0); }

      if ( isTagAcc.isAvailable( *lep_itr ) ) {

    	m_lepton_isTag.push_back( isTagAcc( *lep_itr ) );

    	// fill lepton TAG variables
    	if ( isTagAcc( *lep_itr ) ){

	  m_lepton_tag_pt.push_back(  lep_itr->pt() );
    	  m_lepton_tag_eta.push_back(  lep_itr->eta() );
    	  m_lepton_tag_flavour.push_back( lep_flav );
	  m_lepton_tag_charge.push_back( lep_charge );

	  if ( isTrigMatchedLepAcc.isAvailable( *lep_itr ) )           { m_lepton_tag_isTrigMatched.push_back(  isTrigMatchedLepAcc( *lep_itr ) ); }
    	  else { m_lepton_tag_isTrigMatched .push_back(  -1 ); }
	  if ( isTightAcc.isAvailable( *lep_itr ) )		    { m_lepton_tag_isTight.push_back(  isTightAcc( *lep_itr ) ); }
    	  else { m_lepton_tag_isTight .push_back(  -1 ); }
	  if ( isMediumAcc.isAvailable( *lep_itr ) )		    { m_lepton_tag_isMedium.push_back(  isMediumAcc( *lep_itr ) ); }
    	  else { m_lepton_tag_isMedium .push_back(  -1 ); }
    	  if ( isTruthMatchedAcc.isAvailable( *lep_itr ) )	    { m_lepton_tag_isTruthMatched.push_back( isTruthMatchedAcc( *lep_itr ) ); }
    	  else {  m_lepton_tag_isTruthMatched .push_back( -1 ); }
    	  if ( isChFlipAcc.isAvailable( *lep_itr ) )		    { m_lepton_tag_isChFlip.push_back(  isChFlipAcc( *lep_itr ) ); }
    	  else { m_lepton_tag_isChFlip.push_back( -1 ); }
    	  if ( isBremAcc.isAvailable( *lep_itr ) )		    { m_lepton_tag_isBrem.push_back(  isBremAcc( *lep_itr ) ); }
    	  else { m_lepton_tag_isBrem.push_back( -1 ); }
   	  if ( truthPdgIdAcc.isAvailable( *lep_itr ) )	            { m_lepton_tag_truthPdgId.push_back(  truthPdgIdAcc( *lep_itr ) ); }
    	  else { m_lepton_tag_truthPdgId.push_back( -1 ); }
    	  if ( truthTypeAcc.isAvailable( *lep_itr ) )	            { m_lepton_tag_truthType.push_back(  truthTypeAcc( *lep_itr ) ); }
    	  else { m_lepton_tag_truthType.push_back( -1 ); }
    	  if ( truthOriginAcc.isAvailable( *lep_itr ) )	            { m_lepton_tag_truthOrigin.push_back(  truthOriginAcc( *lep_itr ) ); }
    	  else { m_lepton_tag_truthOrigin.push_back( -1 ); }
    	  if ( truthStatusAcc.isAvailable( *lep_itr ) )	            { m_lepton_tag_truthStatus.push_back(  truthStatusAcc( *lep_itr ) ); }
    	  else { m_lepton_tag_truthStatus.push_back( -1 ); }
          if ( ancestorTruthTypeAcc.isAvailable( *lep_itr ) )       { m_lepton_tag_ancestorTruthType.push_back( ancestorTruthTypeAcc( *lep_itr ) ); }
          else   { m_lepton_tag_ancestorTruthType.push_back(0); }
          if ( ancestorTruthPdgIdAcc.isAvailable( *lep_itr ) )      { m_lepton_tag_ancestorTruthPdgId.push_back( ancestorTruthPdgIdAcc( *lep_itr ) ); }
          else   { m_lepton_tag_ancestorTruthPdgId.push_back(0); }
          if ( ancestorTruthOriginAcc.isAvailable( *lep_itr ) )     { m_lepton_tag_ancestorTruthOrigin.push_back( ancestorTruthOriginAcc( *lep_itr ) ); }
          else   { m_lepton_tag_ancestorTruthOrigin.push_back(0); }
          if ( ancestorTruthStatusAcc.isAvailable( *lep_itr ) )     { m_lepton_tag_ancestorTruthStatus.push_back( ancestorTruthStatusAcc( *lep_itr ) ); }
          else   { m_lepton_tag_ancestorTruthStatus.push_back(0); }

	} else {
    	// fill lepton PROBE variables

	  m_lepton_probe_pt.push_back(  lep_itr->pt() );
    	  m_lepton_probe_eta.push_back(  lep_itr->eta() );
    	  m_lepton_probe_flavour.push_back( lep_flav );
	  m_lepton_probe_charge.push_back( lep_charge );

	  if ( isTrigMatchedLepAcc.isAvailable( *lep_itr ) )           { m_lepton_probe_isTrigMatched.push_back(  isTrigMatchedLepAcc( *lep_itr ) ); }
    	  else { m_lepton_probe_isTrigMatched .push_back(  -1 ); }
    	  if ( isTightAcc.isAvailable( *lep_itr ) )		    { m_lepton_probe_isTight.push_back(  isTightAcc( *lep_itr ) ); }
    	  else { m_lepton_probe_isTight.push_back(  -1 ); }
    	  if ( isMediumAcc.isAvailable( *lep_itr ) )		    { m_lepton_probe_isMedium.push_back(  isMediumAcc( *lep_itr ) ); }
    	  else { m_lepton_probe_isMedium.push_back(  -1 ); }
    	  if ( isTruthMatchedAcc.isAvailable( *lep_itr ) )	    { m_lepton_probe_isTruthMatched.push_back( isTruthMatchedAcc( *lep_itr ) ); }
    	  else {  m_lepton_probe_isTruthMatched.push_back( -1 ); }
    	  if ( isChFlipAcc.isAvailable( *lep_itr ) )		    { m_lepton_probe_isChFlip.push_back(  isChFlipAcc( *lep_itr ) ); }
    	  else { m_lepton_probe_isChFlip.push_back( -1 ); }
    	  if ( isBremAcc.isAvailable( *lep_itr ) )		    { m_lepton_probe_isBrem.push_back(  isBremAcc( *lep_itr ) ); }
    	  else { m_lepton_probe_isBrem.push_back( -1 ); }
   	  if ( truthPdgIdAcc.isAvailable( *lep_itr ) )	            { m_lepton_probe_truthPdgId.push_back(  truthPdgIdAcc( *lep_itr ) ); }
    	  else { m_lepton_probe_truthPdgId.push_back( -1 ); }
    	  if ( truthTypeAcc.isAvailable( *lep_itr ) )	            { m_lepton_probe_truthType.push_back(  truthTypeAcc( *lep_itr ) ); }
    	  else { m_lepton_probe_truthType.push_back( -1 ); }
    	  if ( truthOriginAcc.isAvailable( *lep_itr ) )	            { m_lepton_probe_truthOrigin.push_back(  truthOriginAcc( *lep_itr ) ); }
    	  else { m_lepton_probe_truthOrigin.push_back( -1 ); }
    	  if ( truthStatusAcc.isAvailable( *lep_itr ) )	            { m_lepton_probe_truthStatus.push_back(  truthStatusAcc( *lep_itr ) ); }
    	  else { m_lepton_probe_truthStatus.push_back( -1 ); }
          if ( ancestorTruthTypeAcc.isAvailable( *lep_itr ) )       { m_lepton_probe_ancestorTruthType.push_back( ancestorTruthTypeAcc( *lep_itr ) ); }
          else   { m_lepton_probe_ancestorTruthType.push_back(0); }
          if ( ancestorTruthPdgIdAcc.isAvailable( *lep_itr ) )      { m_lepton_probe_ancestorTruthPdgId.push_back( ancestorTruthPdgIdAcc( *lep_itr ) ); }
          else   { m_lepton_probe_ancestorTruthPdgId.push_back(0); }
          if ( ancestorTruthOriginAcc.isAvailable( *lep_itr ) )     { m_lepton_probe_ancestorTruthOrigin.push_back( ancestorTruthOriginAcc( *lep_itr ) ); }
          else   { m_lepton_probe_ancestorTruthOrigin.push_back(0); }
          if ( ancestorTruthStatusAcc.isAvailable( *lep_itr ) )     { m_lepton_probe_ancestorTruthStatus.push_back( ancestorTruthStatusAcc( *lep_itr ) ); }
          else   { m_lepton_probe_ancestorTruthStatus.push_back(0); }

    	}

      }
      else   { m_lepton_isTag.push_back(-1); }

      ++m_nlep;
  }

}
