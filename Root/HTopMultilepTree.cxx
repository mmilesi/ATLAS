#include "HTopMultilepAnalysis/HTopMultilepTree.h"

HTopMultilepTree :: HTopMultilepTree(xAOD::TEvent* event, TTree* tree, TFile* file, const float units, bool debug ) : 
  HelpTreeBase(event, tree, file, units, debug )
{
  Info("HTopMultilepTree", "Creating output TTree");
}

HTopMultilepTree :: ~HTopMultilepTree() 
{
}

void HTopMultilepTree::AddEventUser(const std::string detailStrUser)
{
  // event variables
  m_tree->Branch("isMC",              &m_is_mc, "isMC/I");  
  m_tree->Branch("ystar",             &m_ystar, "ystar/F");
  m_tree->Branch("categoryFlag",      &m_categoryFlag, "categoryFlag/i"); 
  m_tree->Branch("nMediumBjets",      &m_nBjetsMedium, "nMediumBjets/i");
  m_tree->Branch("isDilepton",        &m_isDilepton,  "isDilepton/I"); 
  m_tree->Branch("isTrilepton",       &m_isTrilepton, "isTrilepton/I"); 
  m_tree->Branch("isDileptonSS",      &m_isDileptonSS, "isDileptonSS/I"); 
  m_tree->Branch("isDileptonSS_elel", &m_isDileptonSS_elel, "isDileptonSS_elel/I");   
  m_tree->Branch("isDileptonSS_mumu", &m_isDileptonSS_mumu, "isDileptonSS_mumu/I");   
  m_tree->Branch("isDileptonSS_elmu", &m_isDileptonSS_elmu, "isDileptonSS_elmu/I");    
  m_tree->Branch("isSS01",            &m_isSS01, "isSS01/I"); 
  m_tree->Branch("isSS12",            &m_isSS12, "isSS12/I"); 
  m_tree->Branch("MMWeight",          &m_MMWeight);
  m_tree->Branch("FFWeight",          &m_FFWeight);
  m_tree->Branch("mll01",             &m_mll01, "mll01/F");
  m_tree->Branch("mll02",             &m_mll02, "mll02/F");
  m_tree->Branch("mll12",             &m_mll12, "mll12/F");
  m_tree->Branch("mlll012",           &m_mlll012, "mlll012/F");
  m_tree->Branch("isTT",              &m_isTT, "isTT/I"); 
  m_tree->Branch("isTL",              &m_isTL, "isTL/I"); 
  m_tree->Branch("isLT",              &m_isLT, "isLT/I"); 
  m_tree->Branch("isLL",              &m_isLL, "isLL/I");   
  m_tree->Branch("isNonTightEvent",   &m_isNonTightEvent,    "isNonTightEvent/I");
  m_tree->Branch("isProbeElEvent",    &m_isProbeElEvent,     "isProbeElEvent/I");
  m_tree->Branch("isProbeMuEvent",    &m_isProbeMuEvent,     "isProbeMuEvent/I");
}

void HTopMultilepTree::AddJetsUser(const std::string detailStrUser)
{
  // jet variables
  m_tree->Branch("jet_m",     &m_jet_m);    
  m_tree->Branch("jet_clean", &m_jet_clean);
}

void HTopMultilepTree::AddMuonsUser(const std::string detailStrUser)
{
  // muon variables  
  m_tree->Branch("muon_isCombined",                &m_muon_isCombined);
  m_tree->Branch("muon_isIsolated",                &m_muon_isIsolated);  
  m_tree->Branch("muon_isTight",                   &m_muon_isTight);
  m_tree->Branch("muon_isTruthMatchedMu",          &m_muon_isTruthMatched); 
  m_tree->Branch("muon_isTruthMatchedMuIso",       &m_muon_isTruthMatchedIso); 
  m_tree->Branch("muon_isTruthMatchedMuNonIso",	   &m_muon_isTruthMatchedNonIso); 
  m_tree->Branch("muon_isTruthMatchedMuSecondary", &m_muon_isTruthMatchedSecondary); 
  m_tree->Branch("muon_isTruthMatchedMuNoProdVtx", &m_muon_isTruthMatchedNoProdVtx);
  m_tree->Branch("muon_isTruthMatchedOther",	   &m_muon_isTruthMatchedOther);
  m_tree->Branch("muon_isTruthMatchedUnknown",	   &m_muon_isTruthMatchedUnknown); 
  m_tree->Branch("muon_isChFlip",	           &m_muon_isChFlip); 
  m_tree->Branch("muon_isBrem",	                   &m_muon_isBrem); 
  m_tree->Branch("muon_truthType",                 &m_muon_truthType);  
  m_tree->Branch("muon_truthOrigin",               &m_muon_truthOrigin);
  m_tree->Branch("muon_truthStatus",               &m_muon_truthStatus);  
  m_tree->Branch("muon_isTrigMatched",             &m_muon_isTrigMatched); 
  m_tree->Branch("muon_isOS",                      &m_muon_isOS);
  m_tree->Branch("muon_isClosestSS",               &m_muon_isClosestSS);
  m_tree->Branch("muon_isTag",                     &m_muon_isTag); 
  // muon TAG variables 
  m_tree->Branch("muon_tag_pt",	                         &m_muon_tag_pt			        );
  m_tree->Branch("muon_tag_eta",	                 &m_muon_tag_eta			);
  m_tree->Branch("muon_tag_isTight",	                 &m_muon_tag_isTight			); 
  m_tree->Branch("muon_tag_isTruthMatchedMu",            &m_muon_tag_isTruthMatched		); 
  m_tree->Branch("muon_tag_isTruthMatchedMuIso",         &m_muon_tag_isTruthMatchedIso        	); 
  m_tree->Branch("muon_tag_isTruthMatchedMuNonIso",      &m_muon_tag_isTruthMatchedNonIso	); 
  m_tree->Branch("muon_tag_isTruthMatchedMuSecondary",   &m_muon_tag_isTruthMatchedSecondary	); 
  m_tree->Branch("muon_tag_isTruthMatchedMuNoProdVtx",   &m_muon_tag_isTruthMatchedNoProdVtx	);
  m_tree->Branch("muon_tag_isTruthMatchedOther",         &m_muon_tag_isTruthMatchedOther	);
  m_tree->Branch("muon_tag_isTruthMatchedUnknown",       &m_muon_tag_isTruthMatchedUnknown	);  
  m_tree->Branch("muon_tag_isChFlip",	                 &m_muon_tag_isChFlip  		        );  
  m_tree->Branch("muon_tag_isBrem",	                 &m_muon_tag_isBrem			);  
  m_tree->Branch("muon_tag_truthType",                   &m_muon_tag_truthType 		        );  
  m_tree->Branch("muon_tag_truthOrigin",                 &m_muon_tag_truthOrigin		);  
  m_tree->Branch("muon_tag_truthStatus",                 &m_muon_tag_truthStatus		);  
  // muon PROBE variables 								
  m_tree->Branch("muon_probe_pt",	                 &m_muon_probe_pt                       );
  m_tree->Branch("muon_probe_eta",	                 &m_muon_probe_eta                      );
  m_tree->Branch("muon_probe_isTight",	                 &m_muon_probe_isTight                  );
  m_tree->Branch("muon_probe_isTruthMatchedMu",          &m_muon_probe_isTruthMatched           );
  m_tree->Branch("muon_probe_isTruthMatchedMuIso",       &m_muon_probe_isTruthMatchedIso        );
  m_tree->Branch("muon_probe_isTruthMatchedMuNonIso",	 &m_muon_probe_isTruthMatchedNonIso     );	  
  m_tree->Branch("muon_probe_isTruthMatchedMuSecondary", &m_muon_probe_isTruthMatchedSecondary  ); 
  m_tree->Branch("muon_probe_isTruthMatchedMuNoProdVtx", &m_muon_probe_isTruthMatchedNoProdVtx  ); 
  m_tree->Branch("muon_probe_isTruthMatchedOther",	 &m_muon_probe_isTruthMatchedOther      );
  m_tree->Branch("muon_probe_isTruthMatchedUnknown",	 &m_muon_probe_isTruthMatchedUnknown    );	  
  m_tree->Branch("muon_probe_isChFlip",	                 &m_muon_probe_isChFlip                 );
  m_tree->Branch("muon_probe_isBrem",	                 &m_muon_probe_isBrem                   );
  m_tree->Branch("muon_probe_truthType",                 &m_muon_probe_truthType                );
  m_tree->Branch("muon_probe_truthOrigin",               &m_muon_probe_truthOrigin              );
  m_tree->Branch("muon_probe_truthStatus",               &m_muon_probe_truthStatus              );    
}

void HTopMultilepTree::AddElectronsUser(const std::string detailStrUser)
{   		     
  // electron variables  
  m_tree->Branch("el_calo_eta",                  &m_electron_calo_eta);
  m_tree->Branch("el_crack",                     &m_electron_crack);
  m_tree->Branch("el_LHMedium",                  &m_electron_LHMedium);   
  m_tree->Branch("el_LHTight",                   &m_electron_LHTight);
  m_tree->Branch("el_LHVeryTight",               &m_electron_LHVeryTight);
  m_tree->Branch("el_mediumPP",                  &m_electron_MediumPP);   
  m_tree->Branch("el_tightPP",                   &m_electron_TightPP);  
  m_tree->Branch("el_isIsolated",                &m_electron_isIsolated); 
  m_tree->Branch("el_isTight",                   &m_electron_isTight); 
  m_tree->Branch("el_isTruthMatchedEl",          &m_electron_isTruthMatched); 
  m_tree->Branch("el_isTruthMatchedElIso",       &m_electron_isTruthMatchedIso); 
  m_tree->Branch("el_isTruthMatchedElNonIso",    &m_electron_isTruthMatchedNonIso); 
  m_tree->Branch("el_isTruthMatchedElSecondary", &m_electron_isTruthMatchedSecondary); 
  m_tree->Branch("el_isTruthMatchedElNoProdVtx", &m_electron_isTruthMatchedNoProdVtx);
  m_tree->Branch("el_isTruthMatchedOther",       &m_electron_isTruthMatchedOther); 
  m_tree->Branch("el_isTruthMatchedUnknown",     &m_electron_isTruthMatchedUnknown); 
  m_tree->Branch("el_isChFlip",                  &m_electron_isChFlip); 
  m_tree->Branch("el_isBrem",                    &m_electron_isBrem); 
  m_tree->Branch("el_truthType",	         &m_electron_truthType);  
  m_tree->Branch("el_truthOrigin",               &m_electron_truthOrigin);  
  m_tree->Branch("el_truthStatus",               &m_electron_truthStatus);  
  m_tree->Branch("el_isTrigMatched",             &m_electron_isTrigMatched);
  m_tree->Branch("el_isOS",	                 &m_electron_isOS);
  m_tree->Branch("el_isClosestSS",               &m_electron_isClosestSS);
  m_tree->Branch("el_isTag",                     &m_electron_isTag); 
  // electron TAG variables 
  m_tree->Branch("el_tag_pt",	                       &m_electron_tag_pt			  );
  m_tree->Branch("el_tag_eta",	                       &m_electron_tag_eta			  );
  m_tree->Branch("el_tag_isTight",	               &m_electron_tag_isTight 		          );
  m_tree->Branch("el_tag_isTruthMatchedEl",            &m_electron_tag_isTruthMatched  	          );
  m_tree->Branch("el_tag_isTruthMatchedElIso",         &m_electron_tag_isTruthMatchedIso	  );
  m_tree->Branch("el_tag_isTruthMatchedElNonIso",      &m_electron_tag_isTruthMatchedNonIso	  );
  m_tree->Branch("el_tag_isTruthMatchedElSecondary",   &m_electron_tag_isTruthMatchedSecondary    );
  m_tree->Branch("el_tag_isTruthMatchedElNoProdVtx",   &m_electron_tag_isTruthMatchedNoProdVtx    );
  m_tree->Branch("el_tag_isTruthMatchedOther",         &m_electron_tag_isTruthMatchedOther	  );
  m_tree->Branch("el_tag_isTruthMatchedUnknown",       &m_electron_tag_isTruthMatchedUnknown	  );
  m_tree->Branch("el_tag_isChFlip",	               &m_electron_tag_isChFlip		          );
  m_tree->Branch("el_tag_isBrem",	               &m_electron_tag_isBrem  		          );
  m_tree->Branch("el_tag_truthType",                   &m_electron_tag_truthType		  );
  m_tree->Branch("el_tag_truthOrigin",                 &m_electron_tag_truthOrigin		  );
  m_tree->Branch("el_tag_truthStatus",                 &m_electron_tag_truthStatus		  );
  // electron PROBE variables 								
  m_tree->Branch("el_probe_pt",	                       &m_electron_probe_pt                       );
  m_tree->Branch("el_probe_eta",	               &m_electron_probe_eta                      );
  m_tree->Branch("el_probe_isTight",	               &m_electron_probe_isTight                  );
  m_tree->Branch("el_probe_isTruthMatchedEl",          &m_electron_probe_isTruthMatched           );
  m_tree->Branch("el_probe_isTruthMatchedElIso",       &m_electron_probe_isTruthMatchedIso        );
  m_tree->Branch("el_probe_isTruthMatchedElNonIso",    &m_electron_probe_isTruthMatchedNonIso     );
  m_tree->Branch("el_probe_isTruthMatchedElSecondary", &m_electron_probe_isTruthMatchedSecondary  );
  m_tree->Branch("el_probe_isTruthMatchedElNoProdVtx", &m_electron_probe_isTruthMatchedNoProdVtx  );
  m_tree->Branch("el_probe_isTruthMatchedOther",       &m_electron_probe_isTruthMatchedOther      );
  m_tree->Branch("el_probe_isTruthMatchedUnknown",     &m_electron_probe_isTruthMatchedUnknown    );
  m_tree->Branch("el_probe_isChFlip",	               &m_electron_probe_isChFlip		  );
  m_tree->Branch("el_probe_isBrem",	               &m_electron_probe_isBrem		          );
  m_tree->Branch("el_probe_truthType",                 &m_electron_probe_truthType                );
  m_tree->Branch("el_probe_truthOrigin",               &m_electron_probe_truthOrigin              );
  m_tree->Branch("el_probe_truthStatus",               &m_electron_probe_truthStatus              );
}


void HTopMultilepTree::AddLeptons()
{   		     
  // always
  m_tree->Branch("nlep",   &m_nlep, "nlep/I");

  m_tree->Branch("lep_pt",          		   &m_lepton_pt);
  m_tree->Branch("lep_phi",         		   &m_lepton_phi);
  m_tree->Branch("lep_eta",         		   &m_lepton_eta);
  m_tree->Branch("lep_m",           		   &m_lepton_m);
  m_tree->Branch("lep_charge",      		   &m_lepton_charge);
  m_tree->Branch("lep_flavour",     		   &m_lepton_flavour);  
  m_tree->Branch("lep_isIsolated",  		   &m_lepton_isIsolated); 
  m_tree->Branch("lep_isTight",                    &m_lepton_isTight); 
  m_tree->Branch("lep_isTruthMatchedLep",          &m_lepton_isTruthMatched); 
  m_tree->Branch("lep_isTruthMatchedLepIso",       &m_lepton_isTruthMatchedIso); 
  m_tree->Branch("lep_isTruthMatchedLepNonIso",    &m_lepton_isTruthMatchedNonIso); 
  m_tree->Branch("lep_isTruthMatchedLepSecondary", &m_lepton_isTruthMatchedSecondary); 
  m_tree->Branch("lep_isTruthMatchedLepNoProdVtx", &m_lepton_isTruthMatchedNoProdVtx); 
  m_tree->Branch("lep_isTruthMatchedOther",        &m_lepton_isTruthMatchedOther); 
  m_tree->Branch("lep_isTruthMatchedUnknown",      &m_lepton_isTruthMatchedUnknown); 
  m_tree->Branch("lep_isChFlip",                   &m_lepton_isChFlip); 
  m_tree->Branch("lep_isBrem",                     &m_lepton_isBrem); 
  m_tree->Branch("lep_truthType",                  &m_lepton_truthType);  
  m_tree->Branch("lep_truthOrigin",                &m_lepton_truthOrigin);
  m_tree->Branch("lep_truthStatus",                &m_lepton_truthStatus);  
  m_tree->Branch("lep_isTrigMatched",              &m_lepton_isTrigMatched);
  m_tree->Branch("lep_isOS",	                   &m_lepton_isOS);
  m_tree->Branch("lep_isClosestSS",                &m_lepton_isClosestSS);
  m_tree->Branch("lep_isTag",                      &m_lepton_isTag); 
  // lepton TAG variables 
  m_tree->Branch("lep_tag_pt",	                        &m_lepton_tag_pt  			);
  m_tree->Branch("lep_tag_eta",	                        &m_lepton_tag_eta 			);
  m_tree->Branch("lep_tag_flavour",	                &m_lepton_tag_flavour 			);
  m_tree->Branch("lep_tag_charge",	                &m_lepton_tag_charge 			);  
  m_tree->Branch("lep_tag_isTight",	                &m_lepton_tag_isTight			); 
  m_tree->Branch("lep_tag_isTruthMatchedLep",           &m_lepton_tag_isTruthMatched	        ); 
  m_tree->Branch("lep_tag_isTruthMatchedLepIso",        &m_lepton_tag_isTruthMatchedIso         ); 
  m_tree->Branch("lep_tag_isTruthMatchedLepNonIso",     &m_lepton_tag_isTruthMatchedNonIso      ); 
  m_tree->Branch("lep_tag_isTruthMatchedLepSecondary",  &m_lepton_tag_isTruthMatchedSecondary   ); 
  m_tree->Branch("lep_tag_isTruthMatchedLepNoProdVtx",  &m_lepton_tag_isTruthMatchedNoProdVtx   );
  m_tree->Branch("lep_tag_isTruthMatchedOther",         &m_lepton_tag_isTruthMatchedOther 	);
  m_tree->Branch("lep_tag_isTruthMatchedUnknown",       &m_lepton_tag_isTruthMatchedUnknown	);  
  m_tree->Branch("lep_tag_isChFlip",	                &m_lepton_tag_isChFlip			);  
  m_tree->Branch("lep_tag_isBrem",	                &m_lepton_tag_isBrem			);  
  m_tree->Branch("lep_tag_truthType",                   &m_lepton_tag_truthType		        );  
  m_tree->Branch("lep_tag_truthOrigin",                 &m_lepton_tag_truthOrigin 		);  
  m_tree->Branch("lep_tag_truthStatus",                 &m_lepton_tag_truthStatus 		);  
  // lepton PROBE variables 								
  m_tree->Branch("lep_probe_pt",	                &m_lepton_probe_pt                      );
  m_tree->Branch("lep_probe_eta",	                &m_lepton_probe_eta                     );
  m_tree->Branch("lep_probe_flavour",	                &m_lepton_probe_flavour                 );
  m_tree->Branch("lep_probe_charge",	                &m_lepton_probe_charge                  );
  m_tree->Branch("lep_probe_isTight",	                &m_lepton_probe_isTight                 );
  m_tree->Branch("lep_probe_isTruthMatchedLep",         &m_lepton_probe_isTruthMatched 	        );
  m_tree->Branch("lep_probe_isTruthMatchedLepIso",      &m_lepton_probe_isTruthMatchedIso	);
  m_tree->Branch("lep_probe_isTruthMatchedLepNonIso",   &m_lepton_probe_isTruthMatchedNonIso	);	
  m_tree->Branch("lep_probe_isTruthMatchedLepSecondary",&m_lepton_probe_isTruthMatchedSecondary ); 
  m_tree->Branch("lep_probe_isTruthMatchedLepNoProdVtx",&m_lepton_probe_isTruthMatchedNoProdVtx ); 
  m_tree->Branch("lep_probe_isTruthMatchedOther",       &m_lepton_probe_isTruthMatchedOther     );
  m_tree->Branch("lep_probe_isTruthMatchedUnknown",     &m_lepton_probe_isTruthMatchedUnknown   );
  m_tree->Branch("lep_probe_isChFlip",	                &m_lepton_probe_isChFlip		);
  m_tree->Branch("lep_probe_isBrem",	                &m_lepton_probe_isBrem		        );
  m_tree->Branch("lep_probe_truthType",                 &m_lepton_probe_truthType               );
  m_tree->Branch("lep_probe_truthOrigin",               &m_lepton_probe_truthOrigin             );
  m_tree->Branch("lep_probe_truthStatus",               &m_lepton_probe_truthStatus             );
}

void HTopMultilepTree::AddTausUser(const std::string detailStrUser)
{   		     
  m_tree->Branch("tau_isBDTTight",	                &m_tau_isBDTTight );
}

void HTopMultilepTree::ClearEventUser() {
  m_MMWeight.clear();
  m_FFWeight.clear();
}
  
void HTopMultilepTree::ClearMuonsUser() {  
  // muon variables 
  m_muon_isCombined.clear(); 
  m_muon_isIsolated.clear();
  m_muon_isTight.clear();
  m_muon_isTruthMatched.clear(); 
  m_muon_isTruthMatchedIso.clear(); 
  m_muon_isTruthMatchedNonIso.clear(); 
  m_muon_isTruthMatchedSecondary.clear(); 
  m_muon_isTruthMatchedNoProdVtx.clear(); 
  m_muon_isTruthMatchedOther.clear(); 
  m_muon_isTruthMatchedUnknown.clear();
  m_muon_isChFlip.clear(); 
  m_muon_isBrem.clear();
  m_muon_truthType.clear();  
  m_muon_truthOrigin.clear();
  m_muon_truthStatus.clear();  
  m_muon_isTrigMatched.clear();
  m_muon_isOS.clear();
  m_muon_isClosestSS.clear();
  m_muon_isTag.clear(); 
  m_muon_tag_pt.clear();
  m_muon_tag_eta.clear();
  m_muon_tag_isTight.clear(); 
  m_muon_tag_isTruthMatched.clear(); 
  m_muon_tag_isTruthMatchedIso.clear(); 
  m_muon_tag_isTruthMatchedNonIso.clear(); 
  m_muon_tag_isTruthMatchedSecondary.clear(); 
  m_muon_tag_isTruthMatchedNoProdVtx.clear(); 
  m_muon_tag_isTruthMatchedOther.clear(); 
  m_muon_tag_isTruthMatchedUnknown.clear();
  m_muon_tag_isChFlip.clear(); 
  m_muon_tag_isBrem.clear();
  m_muon_tag_truthType.clear();  
  m_muon_tag_truthOrigin.clear(); 
  m_muon_tag_truthStatus.clear();  
  m_muon_probe_pt.clear();
  m_muon_probe_eta.clear();
  m_muon_probe_isTight.clear(); 
  m_muon_probe_isTruthMatched.clear(); 
  m_muon_probe_isTruthMatchedIso.clear(); 
  m_muon_probe_isTruthMatchedNonIso.clear(); 
  m_muon_probe_isTruthMatchedSecondary.clear(); 
  m_muon_probe_isTruthMatchedNoProdVtx.clear(); 
  m_muon_probe_isTruthMatchedOther.clear(); 
  m_muon_probe_isTruthMatchedUnknown.clear();
  m_muon_probe_isChFlip.clear(); 
  m_muon_probe_isBrem.clear();
  m_muon_probe_truthType.clear();  
  m_muon_probe_truthOrigin.clear(); 
  m_muon_probe_truthStatus.clear();  
}  

void HTopMultilepTree::ClearElectronsUser() {  
  // electron variables
  m_electron_calo_eta.clear();
  m_electron_crack.clear();
  m_electron_LHMedium.clear();   
  m_electron_LHTight.clear();
  m_electron_LHVeryTight.clear(); 
  m_electron_MediumPP.clear();  
  m_electron_TightPP.clear(); 
  m_electron_isIsolated.clear(); 
  m_electron_isTight.clear(); 
  m_electron_isTruthMatched.clear(); 
  m_electron_isTruthMatchedIso.clear(); 
  m_electron_isTruthMatchedNonIso.clear(); 
  m_electron_isTruthMatchedSecondary.clear(); 
  m_electron_isTruthMatchedNoProdVtx.clear();
  m_electron_isTruthMatchedOther.clear();
  m_electron_isTruthMatchedUnknown.clear();
  m_electron_isChFlip.clear(); 
  m_electron_isBrem.clear();
  m_electron_truthType.clear();      
  m_electron_truthOrigin.clear(); 
  m_electron_truthStatus.clear();      
  m_electron_isOS.clear();
  m_electron_isClosestSS.clear();
  m_electron_isTag.clear(); 
  m_electron_tag_pt.clear();
  m_electron_tag_eta.clear();
  m_electron_tag_isTight.clear(); 
  m_electron_tag_isTruthMatched.clear(); 
  m_electron_tag_isTruthMatchedIso.clear(); 
  m_electron_tag_isTruthMatchedNonIso.clear(); 
  m_electron_tag_isTruthMatchedSecondary.clear(); 
  m_electron_tag_isTruthMatchedNoProdVtx.clear(); 
  m_electron_tag_isTruthMatchedOther.clear(); 
  m_electron_tag_isTruthMatchedUnknown.clear();
  m_electron_tag_isChFlip.clear(); 
  m_electron_tag_isBrem.clear();
  m_electron_tag_truthType.clear();  
  m_electron_tag_truthOrigin.clear(); 
  m_electron_tag_truthStatus.clear();  
  m_electron_probe_pt.clear();
  m_electron_probe_eta.clear();
  m_electron_probe_isTight.clear(); 
  m_electron_probe_isTruthMatched.clear(); 
  m_electron_probe_isTruthMatchedIso.clear(); 
  m_electron_probe_isTruthMatchedNonIso.clear(); 
  m_electron_probe_isTruthMatchedSecondary.clear(); 
  m_electron_probe_isTruthMatchedNoProdVtx.clear(); 
  m_electron_probe_isTruthMatchedOther.clear(); 
  m_electron_probe_isTruthMatchedUnknown.clear();
  m_electron_probe_isChFlip.clear(); 
  m_electron_probe_isBrem.clear();
  m_electron_probe_truthType.clear();  
  m_electron_probe_truthOrigin.clear(); 
  m_electron_probe_truthStatus.clear();  
}


void HTopMultilepTree::ClearJetsUser() {
  // jet variables 
  m_jet_m.clear();
  m_jet_clean.clear();
}  


void HTopMultilepTree::ClearLeptons() {
  m_lepton_pt.clear();
  m_lepton_phi.clear();
  m_lepton_eta.clear();
  m_lepton_m.clear();
  m_lepton_charge.clear();
  m_lepton_flavour.clear();  
  m_lepton_isIsolated.clear();
  m_lepton_isTight.clear(); 
  m_lepton_isTruthMatched.clear(); 
  m_lepton_isTruthMatchedIso.clear(); 
  m_lepton_isTruthMatchedNonIso.clear(); 
  m_lepton_isTruthMatchedSecondary.clear(); 
  m_lepton_isTruthMatchedNoProdVtx.clear();
  m_lepton_isTruthMatchedOther.clear();
  m_lepton_isTruthMatchedUnknown.clear();
  m_lepton_isChFlip.clear(); 
  m_lepton_isBrem.clear();
  m_lepton_truthType.clear();      
  m_lepton_truthOrigin.clear(); 
  m_lepton_truthStatus.clear();      
  m_lepton_isOS.clear();
  m_lepton_isClosestSS.clear();
  m_lepton_isTag.clear(); 
  m_lepton_tag_pt.clear();
  m_lepton_tag_eta.clear();
  m_lepton_tag_flavour.clear();
  m_lepton_tag_charge.clear();
  m_lepton_tag_isTight.clear(); 
  m_lepton_tag_isTruthMatched.clear(); 
  m_lepton_tag_isTruthMatchedIso.clear(); 
  m_lepton_tag_isTruthMatchedNonIso.clear(); 
  m_lepton_tag_isTruthMatchedSecondary.clear(); 
  m_lepton_tag_isTruthMatchedNoProdVtx.clear(); 
  m_lepton_tag_isTruthMatchedOther.clear(); 
  m_lepton_tag_isTruthMatchedUnknown.clear();
  m_lepton_tag_isChFlip.clear(); 
  m_lepton_tag_isBrem.clear();
  m_lepton_tag_truthType.clear();  
  m_lepton_tag_truthOrigin.clear(); 
  m_lepton_tag_truthStatus.clear();  
  m_lepton_probe_pt.clear();
  m_lepton_probe_eta.clear();
  m_lepton_probe_flavour.clear();
  m_lepton_probe_charge.clear();
  m_lepton_probe_isTight.clear(); 
  m_lepton_probe_isTruthMatched.clear(); 
  m_lepton_probe_isTruthMatchedIso.clear(); 
  m_lepton_probe_isTruthMatchedNonIso.clear(); 
  m_lepton_probe_isTruthMatchedSecondary.clear(); 
  m_lepton_probe_isTruthMatchedNoProdVtx.clear(); 
  m_lepton_probe_isTruthMatchedOther.clear(); 
  m_lepton_probe_isTruthMatchedUnknown.clear();
  m_lepton_probe_isChFlip.clear(); 
  m_lepton_probe_isBrem.clear();
  m_lepton_probe_truthType.clear();  
  m_lepton_probe_truthOrigin.clear(); 
  m_lepton_probe_truthStatus.clear();   
}  

void HTopMultilepTree::ClearTausUser() {  
  m_tau_isBDTTight.clear();
}

void HTopMultilepTree::FillEventUser( const xAOD::EventInfo* eventInfo ) { 
  
  m_is_mc              =  ( eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) );
  m_ystar              =  ( eventInfo->isAvailable< float >( "ystar" ) )               ?  eventInfo->auxdecor< float >( "ystar" )                 :   999.0;
  m_nBjetsMedium       =  ( eventInfo->isAvailable< unsigned int >( "nBjetsMedium" ) ) ?  eventInfo->auxdecor< unsigned int >( "nBjetsMedium" )   :   -1;
  m_categoryFlag       =  ( eventInfo->isAvailable< unsigned int >( "categoryFlag" ) ) ?  eventInfo->auxdecor< unsigned int >( "categoryFlag" )   :   -1;
  m_isDilepton         =  ( eventInfo->isAvailable< char >( "isDilepton" ) )           ?  eventInfo->auxdecor< char >( "isDilepton" )             :   -1;
  m_isTrilepton        =  ( eventInfo->isAvailable< char >( "isTrilepton" ) )          ?  eventInfo->auxdecor< char >( "isTrilepton" )            :   -1;
  m_isDileptonSS       =  ( eventInfo->isAvailable< char >( "isDileptonSS" ) )         ?  eventInfo->auxdecor< char >( "isDileptonSS" )           :   -1;
  m_isDileptonSS_elel  =  ( eventInfo->isAvailable< char >( "isDileptonSS_elel" ) )    ?  eventInfo->auxdecor< char >( "isDileptonSS_elel" )      :   -1;
  m_isDileptonSS_elmu  =  ( eventInfo->isAvailable< char >( "isDileptonSS_elmu" ) )    ?  eventInfo->auxdecor< char >( "isDileptonSS_elmu" )      :   -1;
  m_isDileptonSS_mumu  =  ( eventInfo->isAvailable< char >( "isDileptonSS_mumu" ) )    ?  eventInfo->auxdecor< char >( "isDileptonSS_mumu" )      :   -1;
  m_isSS01   	       =  ( eventInfo->isAvailable< char >( "isSS01" ) )               ?  eventInfo->auxdecor< char >( "isSS01" )                 :   -1;
  m_isSS12   	       =  ( eventInfo->isAvailable< char >( "isSS12" ) )    	       ?  eventInfo->auxdecor< char >( "isSS12" )                 :   -1;
  m_MMWeight 	       =  ( eventInfo->isAvailable< std::vector<double> >( "MMWeight" ) )   ?  eventInfo->auxdecor< std::vector<double> >( "MMWeight" ) :  std::vector<double>( 5, 1.0 );
  m_FFWeight 	       =  ( eventInfo->isAvailable< std::vector<double> >( "FFWeight" ) )   ?  eventInfo->auxdecor< std::vector<double> >( "FFWeight" ) :  std::vector<double>( 5, 1.0 );
  m_mll01    	       =  ( eventInfo->isAvailable< float >( "mll01" ) )    	       ?  ( eventInfo->auxdecor< float >( "mll01" ) / m_units )   : -1.0;
  m_mll02    	       =  ( eventInfo->isAvailable< float >( "mll02" ) )    	       ?  ( eventInfo->auxdecor< float >( "mll02" ) / m_units )   : -1.0;
  m_mll12    	       =  ( eventInfo->isAvailable< float >( "mll12" ) )    	       ?  ( eventInfo->auxdecor< float >( "mll12" ) / m_units )   : -1.0;
  m_mlll012  	       =  ( eventInfo->isAvailable< float >( "mlll012" ) )  	       ?  ( eventInfo->auxdecor< float >( "mlll012" ) / m_units ) : -1.0;
  m_isTT   	       =  ( eventInfo->isAvailable< char >( "isTT" ) )                 ?  eventInfo->auxdecor< char >( "isTT" )                   :   -1;
  m_isTL   	       =  ( eventInfo->isAvailable< char >( "isTL" ) )                 ?  eventInfo->auxdecor< char >( "isTL" )                   :   -1;
  m_isLT   	       =  ( eventInfo->isAvailable< char >( "isLT" ) )                 ?  eventInfo->auxdecor< char >( "isLT" )                   :   -1;
  m_isLL   	       =  ( eventInfo->isAvailable< char >( "isLL" ) )                 ?  eventInfo->auxdecor< char >( "isLL" )                   :   -1;
  m_isNonTightEvent    =  ( eventInfo->isAvailable< char >( "isNonTightEvent" ) )      ?  eventInfo->auxdecor< char >( "isNonTightEvent" )        :  -1;
  m_isProbeElEvent     =  ( eventInfo->isAvailable< char >( "isProbeElEvent" ) )       ?  eventInfo->auxdecor< char >( "isProbeElEvent" )         :  -1;
  m_isProbeMuEvent     =  ( eventInfo->isAvailable< char >( "isProbeMuEvent" ) )       ?  eventInfo->auxdecor< char >( "isProbeMuEvent" )         :  -1;
}

void HTopMultilepTree::FillJetsUser( const xAOD::Jet* jet ) {
  
  static SG::AuxElement::Accessor< char > isCleanAcc("cleanJet");
  
  if ( isCleanAcc.isAvailable( *jet ) ) { m_jet_clean.push_back( isCleanAcc( *jet ) ); } 
  else { m_jet_clean.push_back(-1); }
  
  m_jet_m.push_back( jet->m() ); 
  
}

void HTopMultilepTree::FillMuonsUser( const xAOD::Muon* muon ) {
  
  int type = ( muon->muonType() == xAOD::Muon::Combined );
  m_muon_isCombined.push_back( type );

  static SG::AuxElement::Accessor< char > isIsolatedAcc("isIsolated");
  static SG::AuxElement::Accessor< char > isTightAcc("isTight");
  static SG::AuxElement::Accessor< char > isTrigMatchedAcc("isTrigMatched"); 
  static SG::AuxElement::Accessor< char > isTruthMatchedAcc("isTruthMatched");
  static SG::AuxElement::Accessor< char > isTruthMatchedIsoAcc("isTruthMatchedIso");
  static SG::AuxElement::Accessor< char > isTruthMatchedNonIsoAcc("isTruthMatchedNonIso");
  static SG::AuxElement::Accessor< char > isTruthMatchedSecondaryAcc("isTruthMatchedSecondary");
  static SG::AuxElement::Accessor< char > isTruthMatchedNoProdVtxAcc("isTruthMatchedNoProdVtx");   
  static SG::AuxElement::Accessor< char > isTruthMatchedOtherAcc("isTruthMatchedOther");   
  static SG::AuxElement::Accessor< char > isTruthMatchedUnknownAcc("isTruthMatchedUnknown");  
  static SG::AuxElement::Accessor< char > isChFlipAcc("isChFlip");
  static SG::AuxElement::Accessor< char > isBremAcc("isBrem");
  static SG::AuxElement::Accessor< int >  truthMatchTypeAcc("truthMatchType");    
  static SG::AuxElement::Accessor< int >  truthMatchOriginAcc("truthMatchOrigin");
  static SG::AuxElement::Accessor< int >  truthMatchStatusAcc("truthMatchStatus");    
  static SG::AuxElement::Accessor< char > isOSlepAcc("isOSlep");
  static SG::AuxElement::Accessor< char > isClosestSSlepAcc("isClosestSSlep");
  static SG::AuxElement::Accessor< char > isTagAcc("isTag");
 
  if (  isIsolatedAcc.isAvailable( *muon ) ) {  m_muon_isIsolated.push_back( isIsolatedAcc( *muon ) ); } 
  else {  m_muon_isIsolated.push_back(-1); }
  if (  isTightAcc.isAvailable( *muon ) ) {  m_muon_isTight.push_back( isTightAcc( *muon ) ); } 
  else {  m_muon_isTight.push_back(-1); }
  if ( isTrigMatchedAcc.isAvailable( *muon ) ) { m_muon_isTrigMatched.push_back( isTrigMatchedAcc( *muon ) ); }
  else   { m_muon_isTrigMatched.push_back(-1); }
  if ( isTruthMatchedAcc.isAvailable( *muon ) ) { m_muon_isTruthMatched.push_back( isTruthMatchedAcc( *muon ) ); }
  else   { m_muon_isTruthMatched.push_back(-1); }
  if ( isTruthMatchedIsoAcc.isAvailable( *muon ) ) { m_muon_isTruthMatchedIso.push_back( isTruthMatchedIsoAcc( *muon ) ); }
  else   { m_muon_isTruthMatchedIso.push_back(-1); }  
  if ( isTruthMatchedNonIsoAcc.isAvailable( *muon ) ) { m_muon_isTruthMatchedNonIso.push_back( isTruthMatchedNonIsoAcc( *muon ) ); }
  else   { m_muon_isTruthMatchedNonIso.push_back(-1); }  
  if ( isTruthMatchedSecondaryAcc.isAvailable( *muon ) ) { m_muon_isTruthMatchedSecondary.push_back( isTruthMatchedSecondaryAcc( *muon ) ); }
  else   { m_muon_isTruthMatchedSecondary.push_back(-1); }  
  if ( isTruthMatchedNoProdVtxAcc.isAvailable( *muon ) ) { m_muon_isTruthMatchedNoProdVtx.push_back( isTruthMatchedNoProdVtxAcc( *muon ) ); }
  else   { m_muon_isTruthMatchedNoProdVtx.push_back(-1); }    
  if ( isTruthMatchedOtherAcc.isAvailable( *muon ) ) { m_muon_isTruthMatchedOther.push_back( isTruthMatchedOtherAcc( *muon ) ); }
  else   { m_muon_isTruthMatchedOther.push_back(-1); }    
  if ( isTruthMatchedUnknownAcc.isAvailable( *muon ) ) { m_muon_isTruthMatchedUnknown.push_back( isTruthMatchedUnknownAcc( *muon ) ); }
  else   { m_muon_isTruthMatchedUnknown.push_back(-1); }
  if ( isChFlipAcc.isAvailable( *muon ) ) { m_muon_isChFlip.push_back( isChFlipAcc( *muon ) ); }
  else   { m_muon_isChFlip.push_back(-1); }  
  if ( isBremAcc.isAvailable( *muon ) ) { m_muon_isBrem.push_back( isBremAcc( *muon ) ); }
  else   { m_muon_isBrem.push_back(-1); }
  if ( truthMatchTypeAcc.isAvailable( *muon ) ) { m_muon_truthType.push_back( truthMatchTypeAcc( *muon ) ); }
  else   { m_muon_truthType.push_back(0); } 
  if ( truthMatchOriginAcc.isAvailable( *muon ) ) { m_muon_truthOrigin.push_back( truthMatchOriginAcc( *muon ) ); }
  else   { m_muon_truthOrigin.push_back(0); } 
  if ( truthMatchStatusAcc.isAvailable( *muon ) ) { m_muon_truthStatus.push_back( truthMatchStatusAcc( *muon ) ); }
  else   { m_muon_truthStatus.push_back(0); } 
  if (  isOSlepAcc.isAvailable( *muon ) ) { m_muon_isOS.push_back( isOSlepAcc( *muon ) ); }
  else  { m_muon_isOS.push_back(-1); }
  if (  isClosestSSlepAcc.isAvailable( *muon ) ) { m_muon_isClosestSS.push_back( isClosestSSlepAcc( *muon ) ); }
  else  { m_muon_isClosestSS.push_back(-1); }
  if ( isTagAcc.isAvailable( *muon ) ) { 
    
    m_muon_isTag.push_back( isTagAcc( *muon ) ); 
  
    // fill muon TAG variables
    if ( isTagAcc( *muon ) ){
      m_muon_tag_pt.push_back(  muon->pt() );
      m_muon_tag_eta.push_back(  muon->eta() );
      if ( isTightAcc.isAvailable( *muon ) )                 { m_muon_tag_isTight.push_back(  isTightAcc( *muon ) ); }
      else { m_muon_tag_isTight .push_back(  -1 ); }
      if ( isTruthMatchedAcc.isAvailable( *muon ) )          { m_muon_tag_isTruthMatched.push_back( isTruthMatchedAcc( *muon ) ); }
      else {  m_muon_tag_isTruthMatched .push_back( -1 ); }
      if ( isTruthMatchedIsoAcc.isAvailable( *muon ) )       { m_muon_tag_isTruthMatchedIso.push_back( isTruthMatchedIsoAcc( *muon ) ); }
      else {  m_muon_tag_isTruthMatchedIso .push_back( -1 ); }
      if ( isTruthMatchedNonIsoAcc.isAvailable( *muon ) )    { m_muon_tag_isTruthMatchedNonIso.push_back(  isTruthMatchedNonIsoAcc( *muon ) ); }
      else {  m_muon_tag_isTruthMatchedNonIso.push_back( -1 ); }
      if ( isTruthMatchedSecondaryAcc.isAvailable( *muon ) ) { m_muon_tag_isTruthMatchedSecondary.push_back(  isTruthMatchedSecondaryAcc( *muon ) ); }
      else {  m_muon_tag_isTruthMatchedSecondary.push_back( -1 ); }
      if ( isTruthMatchedNoProdVtxAcc.isAvailable( *muon ) ) { m_muon_tag_isTruthMatchedNoProdVtx.push_back(  isTruthMatchedNoProdVtxAcc( *muon ) ); }
      else {  m_muon_tag_isTruthMatchedNoProdVtx.push_back( -1 ); }
      if ( isTruthMatchedOtherAcc.isAvailable( *muon ) )     { m_muon_tag_isTruthMatchedOther.push_back(  isTruthMatchedOtherAcc( *muon ) ); }
      else {  m_muon_tag_isTruthMatchedOther.push_back( -1 ); }
      if ( isTruthMatchedUnknownAcc.isAvailable( *muon ) )   { m_muon_tag_isTruthMatchedUnknown.push_back(  isTruthMatchedUnknownAcc( *muon ) ); }
      else { m_muon_tag_isTruthMatchedUnknown.push_back( -1 ); }
      if ( isChFlipAcc.isAvailable( *muon ) )                { m_muon_tag_isChFlip.push_back(  isChFlipAcc( *muon ) ); }
      else { m_muon_tag_isChFlip.push_back( -1 ); }
      if ( isBremAcc.isAvailable( *muon ) )                  { m_muon_tag_isBrem.push_back(  isBremAcc( *muon ) ); }
      else { m_muon_tag_isBrem.push_back( -1 ); }
      if ( truthMatchTypeAcc.isAvailable( *muon ) )          { m_muon_tag_truthType.push_back(  truthMatchTypeAcc( *muon ) ); }
      else { m_muon_tag_truthType.push_back( -1 ); }
      if ( truthMatchOriginAcc.isAvailable( *muon ) )        { m_muon_tag_truthOrigin.push_back(  truthMatchOriginAcc( *muon ) ); }
      else { m_muon_tag_truthOrigin.push_back( -1 ); }
      if ( truthMatchStatusAcc.isAvailable( *muon ) )        { m_muon_tag_truthStatus.push_back(  truthMatchStatusAcc( *muon ) ); }
      else { m_muon_tag_truthStatus.push_back( -1 ); }
    } else {
    // fill muon PROBE variables
      m_muon_probe_pt.push_back(  muon->pt() );
      m_muon_probe_eta.push_back(  muon->eta() );
      if ( isTightAcc.isAvailable( *muon ) )                 { m_muon_probe_isTight.push_back(  isTightAcc( *muon ) ); }
      else { m_muon_probe_isTight.push_back(  -1 ); }
      if ( isTruthMatchedAcc.isAvailable( *muon ) )          { m_muon_probe_isTruthMatched.push_back( isTruthMatchedAcc( *muon ) ); }
      else {  m_muon_probe_isTruthMatched.push_back( -1 ); }
      if ( isTruthMatchedIsoAcc.isAvailable( *muon ) )       { m_muon_probe_isTruthMatchedIso.push_back( isTruthMatchedIsoAcc( *muon ) ); }
      else {  m_muon_probe_isTruthMatchedIso.push_back( -1 ); }
      if ( isTruthMatchedNonIsoAcc.isAvailable( *muon ) )    { m_muon_probe_isTruthMatchedNonIso.push_back(  isTruthMatchedNonIsoAcc( *muon ) ); }
      else {  m_muon_probe_isTruthMatchedNonIso.push_back( -1 ); }
      if ( isTruthMatchedSecondaryAcc.isAvailable( *muon ) ) { m_muon_probe_isTruthMatchedSecondary.push_back(  isTruthMatchedSecondaryAcc( *muon ) ); }
      else {  m_muon_probe_isTruthMatchedSecondary.push_back( -1 ); }
      if ( isTruthMatchedNoProdVtxAcc.isAvailable( *muon ) ) { m_muon_probe_isTruthMatchedNoProdVtx.push_back(  isTruthMatchedNoProdVtxAcc( *muon ) ); }
      else {  m_muon_probe_isTruthMatchedNoProdVtx.push_back( -1 ); }
      if ( isTruthMatchedOtherAcc.isAvailable( *muon ) )     { m_muon_probe_isTruthMatchedOther.push_back(  isTruthMatchedOtherAcc( *muon ) ); }
      else {  m_muon_probe_isTruthMatchedOther.push_back( -1 ); }
      if ( isTruthMatchedUnknownAcc.isAvailable( *muon ) )   { m_muon_probe_isTruthMatchedUnknown.push_back(  isTruthMatchedUnknownAcc( *muon ) ); }
      else { m_muon_probe_isTruthMatchedUnknown.push_back( -1 ); }
      if ( isChFlipAcc.isAvailable( *muon ) )                { m_muon_probe_isChFlip.push_back(  isChFlipAcc( *muon ) ); }
      else { m_muon_probe_isChFlip.push_back( -1 ); }
      if ( isBremAcc.isAvailable( *muon ) )                  { m_muon_probe_isBrem.push_back(  isBremAcc( *muon ) ); }
      else { m_muon_probe_isBrem.push_back( -1 ); }
      if ( truthMatchTypeAcc.isAvailable( *muon ) )          { m_muon_probe_truthType.push_back(  truthMatchTypeAcc( *muon ) ); }
      else { m_muon_probe_truthType.push_back( -1 ); }
      if ( truthMatchOriginAcc.isAvailable( *muon ) )        { m_muon_probe_truthOrigin.push_back(  truthMatchOriginAcc( *muon ) ); }
      else { m_muon_probe_truthOrigin.push_back( -1 ); }
      if ( truthMatchStatusAcc.isAvailable( *muon ) )        { m_muon_probe_truthStatus.push_back(  truthMatchStatusAcc( *muon ) ); }
      else { m_muon_probe_truthStatus.push_back( -1 ); }
    }
    
  }
  else   { m_muon_isTag.push_back(-1); }  

}

void HTopMultilepTree::FillElectronsUser( const xAOD::Electron* electron ) {
 
  static SG::AuxElement::Accessor< char > isLHMediumAcc("isLHMedium"); 
  static SG::AuxElement::Accessor< char > isLHTightAcc("isLHTight"); 
  static SG::AuxElement::Accessor< char > isLHVeryTightAcc("isLHVeryTight"); 
  static SG::AuxElement::Accessor< char > isCutBasedMediumPPAcc("isCutBasedMediumPP"); 
  static SG::AuxElement::Accessor< char > isCutBasedTightPPAcc("isCutBasedTightPP"); 
  static SG::AuxElement::Accessor< char > isIsolatedAcc("isIsolated");
  static SG::AuxElement::Accessor< char > isTightAcc("isTight");
  static SG::AuxElement::Accessor< char > isTrigMatchedAcc("isTrigMatched"); 
  static SG::AuxElement::Accessor< char > isTruthMatchedAcc("isTruthMatched");
  static SG::AuxElement::Accessor< char > isTruthMatchedIsoAcc("isTruthMatchedIso");
  static SG::AuxElement::Accessor< char > isTruthMatchedNonIsoAcc("isTruthMatchedNonIso");
  static SG::AuxElement::Accessor< char > isTruthMatchedSecondaryAcc("isTruthMatchedSecondary");
  static SG::AuxElement::Accessor< char > isTruthMatchedNoProdVtxAcc("isTruthMatchedNoProdVtx");   
  static SG::AuxElement::Accessor< char > isTruthMatchedOtherAcc("isTruthMatchedOther");  
  static SG::AuxElement::Accessor< char > isTruthMatchedUnknownAcc("isTruthMatchedUnknown");  
  static SG::AuxElement::Accessor< char > isChFlipAcc("isChFlip");
  static SG::AuxElement::Accessor< char > isBremAcc("isBrem");  
  static SG::AuxElement::Accessor< int >  truthMatchTypeAcc("truthMatchType");    
  static SG::AuxElement::Accessor< int >  truthMatchOriginAcc("truthMatchOrigin");   
  static SG::AuxElement::Accessor< int >  truthMatchStatusAcc("truthMatchStatus");    
  static SG::AuxElement::Accessor< char > isOSlepAcc("isOSlep");
  static SG::AuxElement::Accessor< char > isClosestSSlepAcc("isClosestSSlep");
  static SG::AuxElement::Accessor< char > isTagAcc("isTag");

  if ( electron->caloCluster() ) {  m_electron_calo_eta.push_back( electron->caloCluster()->eta() ); } 
  else {  m_electron_calo_eta.push_back(-999.); }

  if ( electron->caloCluster() ) { 
    if ( fabs( electron->caloCluster()->eta() ) > 1.37 && fabs( electron->caloCluster()->eta() ) < 1.52 ) {
      m_electron_crack.push_back(1);
    } else {
      m_electron_crack.push_back(0);
    }
  } else {   
    m_electron_crack.push_back(-1); 
  }

  if (  isLHMediumAcc.isAvailable( *electron ) ) {  m_electron_LHMedium.push_back( isLHMediumAcc( *electron ) ); } 
  else {  m_electron_LHMedium.push_back(-1); }
  if (  isLHTightAcc.isAvailable( *electron ) ) {  m_electron_LHTight.push_back( isLHTightAcc( *electron ) ); } 
  else {  m_electron_LHTight.push_back(-1); }
  if (  isLHVeryTightAcc.isAvailable( *electron ) ) {  m_electron_LHVeryTight.push_back( isLHVeryTightAcc( *electron ) ); } 
  else {  m_electron_LHVeryTight.push_back(-1); }
  if (  isCutBasedMediumPPAcc.isAvailable( *electron ) ) {  m_electron_MediumPP.push_back( isCutBasedMediumPPAcc( *electron ) ); } 
  else {  m_electron_MediumPP.push_back(-1); }
  if (  isCutBasedTightPPAcc.isAvailable( *electron ) ) {  m_electron_TightPP.push_back( isCutBasedTightPPAcc( *electron ) ); } 
  else {  m_electron_TightPP.push_back(-1); }
  if (  isIsolatedAcc.isAvailable( *electron ) ) {  m_electron_isIsolated.push_back( isIsolatedAcc( *electron ) ); } 
  else {  m_electron_isIsolated.push_back(-1); }
  if (  isTightAcc.isAvailable( *electron ) ) {  m_electron_isTight.push_back( isTightAcc( *electron ) ); } 
  else {  m_electron_isTight.push_back(-1); }
  if ( isTrigMatchedAcc.isAvailable( *electron ) ) { m_electron_isTrigMatched.push_back( isTrigMatchedAcc( *electron ) ); }
  else   { m_electron_isTrigMatched.push_back(-1); }
  if ( isTruthMatchedAcc.isAvailable( *electron ) ) { m_electron_isTruthMatched.push_back( isTruthMatchedAcc( *electron ) ); }
  else   { m_electron_isTruthMatched.push_back(-1); }
  if ( isTruthMatchedIsoAcc.isAvailable( *electron ) ) { /*std::cout << "branch el_isTruthMatchedIso: " << (int)isTruthMatchedIsoAcc( *electron ) << std::endl;*/ m_electron_isTruthMatchedIso.push_back( isTruthMatchedIsoAcc( *electron ) ); }
  else   { m_electron_isTruthMatchedIso.push_back(-1); }  
  if ( isTruthMatchedNonIsoAcc.isAvailable( *electron ) ) { /*std::cout << "branch el_isTruthMatchedNonIso: " << (int)isTruthMatchedNonIsoAcc( *electron ) << std::endl;*/ m_electron_isTruthMatchedNonIso.push_back( isTruthMatchedNonIsoAcc( *electron ) ); }
  else   { m_electron_isTruthMatchedNonIso.push_back(-1); }  
  if ( isTruthMatchedSecondaryAcc.isAvailable( *electron ) ) { /*std::cout << "branch el_isTruthMatchedSecondary: " << (int)isTruthMatchedSecondaryAcc( *electron ) << std::endl;*/ m_electron_isTruthMatchedSecondary.push_back( isTruthMatchedSecondaryAcc( *electron ) ); }
  else   { m_electron_isTruthMatchedSecondary.push_back(-1); } 
  if ( isTruthMatchedNoProdVtxAcc.isAvailable( *electron ) ) { m_electron_isTruthMatchedNoProdVtx.push_back( isTruthMatchedNoProdVtxAcc( *electron ) ); }
  else   { m_electron_isTruthMatchedNoProdVtx.push_back(-1); }    
  if ( isTruthMatchedOtherAcc.isAvailable( *electron ) ) { /*std::cout << "branch el_isTruthMatchedOther: " << (int)isTruthMatchedOtherAcc( *electron ) << std::endl;*/ m_electron_isTruthMatchedOther.push_back( isTruthMatchedOtherAcc( *electron ) ); }
  else   { m_electron_isTruthMatchedOther.push_back(-1); }    
  if ( isTruthMatchedUnknownAcc.isAvailable( *electron ) ) { m_electron_isTruthMatchedUnknown.push_back( isTruthMatchedUnknownAcc( *electron ) ); }
  else   { m_electron_isTruthMatchedUnknown.push_back(-1); }
  if ( isChFlipAcc.isAvailable( *electron ) ) { m_electron_isChFlip.push_back( isChFlipAcc( *electron ) ); }
  else   { m_electron_isChFlip.push_back(-1); }  
  if ( isBremAcc.isAvailable( *electron ) ) { m_electron_isBrem.push_back( isBremAcc( *electron ) ); }
  else   { m_electron_isBrem.push_back(-1); }
  if ( truthMatchTypeAcc.isAvailable( *electron ) ) { /*std::cout << "branch el_truthType: " << truthMatchTypeAcc( *electron ) << std::endl;*/ m_electron_truthType.push_back( truthMatchTypeAcc( *electron ) ); }
  else   { m_electron_truthType.push_back(0); } 
  if ( truthMatchOriginAcc.isAvailable( *electron ) ) { m_electron_truthOrigin.push_back( truthMatchOriginAcc( *electron ) ); }
  else   { m_electron_truthOrigin.push_back(0); } 
  if ( truthMatchStatusAcc.isAvailable( *electron ) ) { m_electron_truthStatus.push_back( truthMatchStatusAcc( *electron ) ); }
  else   { m_electron_truthStatus.push_back(0); } 
  if (  isOSlepAcc.isAvailable( *electron ) ) { m_electron_isOS.push_back( isOSlepAcc( *electron ) ); }
  else  { m_electron_isOS.push_back(-1); }
  if (  isClosestSSlepAcc.isAvailable( *electron ) ) { m_electron_isClosestSS.push_back( isClosestSSlepAcc( *electron ) ); }
  else  { m_electron_isClosestSS.push_back(-1); }
  if ( isTagAcc.isAvailable( *electron ) ) { 
    
    m_electron_isTag.push_back( isTagAcc( *electron ) ); 
  
    // fill electron TAG variables
    if ( isTagAcc( *electron ) ){
      m_electron_tag_pt.push_back(  electron->pt() );
      m_electron_tag_eta.push_back(  electron->eta() );
      if ( isTightAcc.isAvailable( *electron ) )                 { m_electron_tag_isTight.push_back(  isTightAcc( *electron ) ); }
      else { m_electron_tag_isTight .push_back(  -1 ); }
      if ( isTruthMatchedAcc.isAvailable( *electron ) )          { m_electron_tag_isTruthMatched.push_back( isTruthMatchedAcc( *electron ) ); }
      else {  m_electron_tag_isTruthMatched .push_back( -1 ); }
      if ( isTruthMatchedIsoAcc.isAvailable( *electron ) )       { m_electron_tag_isTruthMatchedIso.push_back( isTruthMatchedIsoAcc( *electron ) ); }
      else {  m_electron_tag_isTruthMatchedIso .push_back( -1 ); }
      if ( isTruthMatchedNonIsoAcc.isAvailable( *electron ) )    { m_electron_tag_isTruthMatchedNonIso.push_back(  isTruthMatchedNonIsoAcc( *electron ) ); }
      else {  m_electron_tag_isTruthMatchedNonIso.push_back( -1 ); }
      if ( isTruthMatchedSecondaryAcc.isAvailable( *electron ) ) { m_electron_tag_isTruthMatchedSecondary.push_back(  isTruthMatchedSecondaryAcc( *electron ) ); }
      else {  m_electron_tag_isTruthMatchedSecondary.push_back( -1 ); }
      if ( isTruthMatchedNoProdVtxAcc.isAvailable( *electron ) ) { m_electron_tag_isTruthMatchedNoProdVtx.push_back(  isTruthMatchedNoProdVtxAcc( *electron ) ); }
      else {  m_electron_tag_isTruthMatchedNoProdVtx.push_back( -1 ); }
      if ( isTruthMatchedOtherAcc.isAvailable( *electron ) )     { m_electron_tag_isTruthMatchedOther.push_back(  isTruthMatchedOtherAcc( *electron ) ); }
      else {  m_electron_tag_isTruthMatchedOther.push_back( -1 ); }
      if ( isTruthMatchedUnknownAcc.isAvailable( *electron ) )   { m_electron_tag_isTruthMatchedUnknown.push_back(  isTruthMatchedUnknownAcc( *electron ) ); }
      else { m_electron_tag_isTruthMatchedUnknown.push_back( -1 ); }
      if ( isChFlipAcc.isAvailable( *electron ) )                { m_electron_tag_isChFlip.push_back(  isChFlipAcc( *electron ) ); }
      else { m_electron_tag_isChFlip.push_back( -1 ); }
      if ( isBremAcc.isAvailable( *electron ) )                  { m_electron_tag_isBrem.push_back(  isBremAcc( *electron ) ); }
      else { m_electron_tag_isBrem.push_back( -1 ); }
      if ( truthMatchTypeAcc.isAvailable( *electron ) )          { m_electron_tag_truthType.push_back(  truthMatchTypeAcc( *electron ) ); }
      else { m_electron_tag_truthType.push_back( -1 ); }
      if ( truthMatchOriginAcc.isAvailable( *electron ) )        { m_electron_tag_truthOrigin.push_back(  truthMatchOriginAcc( *electron ) ); }
      else { m_electron_tag_truthOrigin.push_back( -1 ); }
      if ( truthMatchStatusAcc.isAvailable( *electron ) )        { m_electron_tag_truthStatus.push_back(  truthMatchStatusAcc( *electron ) ); }
      else { m_electron_tag_truthStatus.push_back( -1 ); }
    } else {
    // fill electron PROBE variables
      m_electron_probe_pt.push_back(  electron->pt() );
      m_electron_probe_eta.push_back(  electron->eta() );
      if ( isTightAcc.isAvailable( *electron ) )                 { m_electron_probe_isTight.push_back(  isTightAcc( *electron ) ); }
      else { m_electron_probe_isTight.push_back(  -1 ); }
      if ( isTruthMatchedAcc.isAvailable( *electron ) )          { m_electron_probe_isTruthMatched.push_back( isTruthMatchedAcc( *electron ) ); }
      else {  m_electron_probe_isTruthMatched.push_back( -1 ); }
      if ( isTruthMatchedIsoAcc.isAvailable( *electron ) )       { m_electron_probe_isTruthMatchedIso.push_back( isTruthMatchedIsoAcc( *electron ) ); }
      else {  m_electron_probe_isTruthMatchedIso.push_back( -1 ); }
      if ( isTruthMatchedNonIsoAcc.isAvailable( *electron ) )    { m_electron_probe_isTruthMatchedNonIso.push_back(  isTruthMatchedNonIsoAcc( *electron ) ); }
      else {  m_electron_probe_isTruthMatchedNonIso.push_back( -1 ); }
      if ( isTruthMatchedSecondaryAcc.isAvailable( *electron ) ) { m_electron_probe_isTruthMatchedSecondary.push_back(  isTruthMatchedSecondaryAcc( *electron ) ); }
      else {  m_electron_probe_isTruthMatchedSecondary.push_back( -1 ); }
      if ( isTruthMatchedNoProdVtxAcc.isAvailable( *electron ) ) { m_electron_probe_isTruthMatchedNoProdVtx.push_back(  isTruthMatchedNoProdVtxAcc( *electron ) ); }
      else {  m_electron_probe_isTruthMatchedNoProdVtx.push_back( -1 ); }
      if ( isTruthMatchedOtherAcc.isAvailable( *electron ) )     { m_electron_probe_isTruthMatchedOther.push_back(  isTruthMatchedOtherAcc( *electron ) ); }
      else {  m_electron_probe_isTruthMatchedOther.push_back( -1 ); }
      if ( isTruthMatchedUnknownAcc.isAvailable( *electron ) )   { m_electron_probe_isTruthMatchedUnknown.push_back(  isTruthMatchedUnknownAcc( *electron ) ); }
      else { m_electron_probe_isTruthMatchedUnknown.push_back( -1 ); }
      if ( isChFlipAcc.isAvailable( *electron ) )                { m_electron_probe_isChFlip.push_back(  isChFlipAcc( *electron ) ); }
      else { m_electron_probe_isChFlip.push_back( -1 ); }
      if ( isBremAcc.isAvailable( *electron ) )                  { m_electron_probe_isBrem.push_back(  isBremAcc( *electron ) ); }
      else { m_electron_probe_isBrem.push_back( -1 ); }
      if ( truthMatchTypeAcc.isAvailable( *electron ) )          { m_electron_probe_truthType.push_back(  truthMatchTypeAcc( *electron ) ); }
      else { m_electron_probe_truthType.push_back( -1 ); }
      if ( truthMatchOriginAcc.isAvailable( *electron ) )        { m_electron_probe_truthOrigin.push_back(  truthMatchOriginAcc( *electron ) ); }
      else { m_electron_probe_truthOrigin.push_back( -1 ); }
      if ( truthMatchStatusAcc.isAvailable( *electron ) )        { m_electron_probe_truthStatus.push_back(  truthMatchStatusAcc( *electron ) ); }
      else { m_electron_probe_truthStatus.push_back( -1 ); }
    }
    
  }
  else   { m_electron_isTag.push_back(-1); }  

}

void HTopMultilepTree::FillFatJetsUser( const xAOD::Jet* fatJet )
{
}

void HTopMultilepTree::FillTausUser( const xAOD::TauJet* tau ) {
  
  static SG::AuxElement::Accessor< char > isTauBDTTightAcc("isTauBDTTight");
  
  if ( isTauBDTTightAcc.isAvailable( *tau ) ) { m_tau_isBDTTight.push_back( isTauBDTTightAcc( *tau ) ); } 
  else { m_tau_isBDTTight.push_back(-1); }
  
}

void HTopMultilepTree::FillLeptons( const xAOD::IParticleContainer* leptons ) {

  this->ClearLeptons();

  m_nlep= 0;
  
  static SG::AuxElement::Accessor< char > isIsolatedAcc("isIsolated");
  static SG::AuxElement::Accessor< char > isTightAcc("isTight");
  static SG::AuxElement::Accessor< char > isTruthMatchedAcc("isTruthMatched");
  static SG::AuxElement::Accessor< char > isTruthMatchedIsoAcc("isTruthMatchedIso");
  static SG::AuxElement::Accessor< char > isTruthMatchedNonIsoAcc("isTruthMatchedNonIso");
  static SG::AuxElement::Accessor< char > isTruthMatchedSecondaryAcc("isTruthMatchedSecondary");
  static SG::AuxElement::Accessor< char > isTruthMatchedNoProdVtxAcc("isTruthMatchedNoProdVtx");   
  static SG::AuxElement::Accessor< char > isTruthMatchedOtherAcc("isTruthMatchedOther");  
  static SG::AuxElement::Accessor< char > isTruthMatchedUnknownAcc("isTruthMatchedUnknown");  
  static SG::AuxElement::Accessor< char > isChFlipAcc("isChFlip");
  static SG::AuxElement::Accessor< char > isBremAcc("isBrem");  
  static SG::AuxElement::Accessor< int >  truthMatchTypeAcc("truthMatchType");    
  static SG::AuxElement::Accessor< int >  truthMatchOriginAcc("truthMatchOrigin");   
  static SG::AuxElement::Accessor< int >  truthMatchStatusAcc("truthMatchStatus");    
  static SG::AuxElement::Accessor< char > isOSlepAcc("isOSlep");
  static SG::AuxElement::Accessor< char > isClosestSSlepAcc("isClosestSSlep");
  static SG::AuxElement::Accessor< char > isTagAcc("isTag");  
  
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
      
      if (  isIsolatedAcc.isAvailable( *lep_itr ) ) {  m_lepton_isIsolated.push_back( isIsolatedAcc( *lep_itr ) ); } 
      else {  m_lepton_isIsolated.push_back(-1); }
      if (  isTightAcc.isAvailable( *lep_itr ) ) {  m_lepton_isTight.push_back( isTightAcc( *lep_itr ) ); } 
      else {  m_lepton_isTight.push_back(-1); }
      if ( isTruthMatchedAcc.isAvailable( *lep_itr ) ) { m_lepton_isTruthMatched.push_back( isTruthMatchedAcc( *lep_itr ) ); }
      else   { m_lepton_isTruthMatched.push_back(-1); }
      if ( isTruthMatchedIsoAcc.isAvailable( *lep_itr ) ) { m_lepton_isTruthMatchedIso.push_back( isTruthMatchedIsoAcc( *lep_itr ) ); }
      else   { m_lepton_isTruthMatchedIso.push_back(-1); }  
      if ( isTruthMatchedNonIsoAcc.isAvailable( *lep_itr ) ) { m_lepton_isTruthMatchedNonIso.push_back( isTruthMatchedNonIsoAcc( *lep_itr ) ); }
      else   { m_lepton_isTruthMatchedNonIso.push_back(-1); }  
      if ( isTruthMatchedSecondaryAcc.isAvailable( *lep_itr ) ) { m_lepton_isTruthMatchedSecondary.push_back( isTruthMatchedSecondaryAcc( *lep_itr ) ); }
      else   { m_lepton_isTruthMatchedSecondary.push_back(-1); } 
      if ( isTruthMatchedNoProdVtxAcc.isAvailable( *lep_itr ) ) { m_lepton_isTruthMatchedNoProdVtx.push_back( isTruthMatchedNoProdVtxAcc( *lep_itr ) ); }
      else   { m_lepton_isTruthMatchedNoProdVtx.push_back(-1); }    
      if ( isTruthMatchedOtherAcc.isAvailable( *lep_itr ) ) { m_lepton_isTruthMatchedOther.push_back( isTruthMatchedOtherAcc( *lep_itr ) ); }
      else   { m_lepton_isTruthMatchedOther.push_back(-1); }    
      if ( isTruthMatchedUnknownAcc.isAvailable( *lep_itr ) ) { m_lepton_isTruthMatchedUnknown.push_back( isTruthMatchedUnknownAcc( *lep_itr ) ); }
      else   { m_lepton_isTruthMatchedUnknown.push_back(-1); }
      if ( isChFlipAcc.isAvailable( *lep_itr ) ) { m_lepton_isChFlip.push_back( isChFlipAcc( *lep_itr ) ); }
      else   { m_lepton_isChFlip.push_back(-1); }  
      if ( isBremAcc.isAvailable( *lep_itr ) ) { m_lepton_isBrem.push_back( isBremAcc( *lep_itr ) ); }
      else   { m_lepton_isBrem.push_back(-1); }
      if ( truthMatchTypeAcc.isAvailable( *lep_itr ) ) { m_lepton_truthType.push_back( truthMatchTypeAcc( *lep_itr ) ); }
      else   { m_lepton_truthType.push_back(0); } 
      if ( truthMatchOriginAcc.isAvailable( *lep_itr ) ) { m_lepton_truthOrigin.push_back( truthMatchOriginAcc( *lep_itr ) ); }
      else   { m_lepton_truthOrigin.push_back(0); } 
      if ( truthMatchStatusAcc.isAvailable( *lep_itr ) ) { m_lepton_truthStatus.push_back( truthMatchStatusAcc( *lep_itr ) ); }
      else   { m_lepton_truthStatus.push_back(0); } 
      if (  isOSlepAcc.isAvailable( *lep_itr ) ) { m_lepton_isOS.push_back( isOSlepAcc( *lep_itr ) ); }
      else  { m_lepton_isOS.push_back(-1); }
      if (  isClosestSSlepAcc.isAvailable( *lep_itr ) ) { m_lepton_isClosestSS.push_back( isClosestSSlepAcc( *lep_itr ) ); }
      else  { m_lepton_isClosestSS.push_back(-1); }
      if ( isTagAcc.isAvailable( *lep_itr ) ) { 
    	
    	m_lepton_isTag.push_back( isTagAcc( *lep_itr ) ); 
      
    	// fill lepton TAG variables
    	if ( isTagAcc( *lep_itr ) ){
    	  m_lepton_tag_pt.push_back(  lep_itr->pt() );
    	  m_lepton_tag_eta.push_back(  lep_itr->eta() );
    	  m_lepton_tag_flavour.push_back( lep_flav );
	  m_lepton_tag_charge.push_back( lep_charge );
    	  if ( isTightAcc.isAvailable( *lep_itr ) )		    { m_lepton_tag_isTight.push_back(  isTightAcc( *lep_itr ) ); }
    	  else { m_lepton_tag_isTight .push_back(  -1 ); }
    	  if ( isTruthMatchedAcc.isAvailable( *lep_itr ) )	    { m_lepton_tag_isTruthMatched.push_back( isTruthMatchedAcc( *lep_itr ) ); }
    	  else {  m_lepton_tag_isTruthMatched .push_back( -1 ); }
    	  if ( isTruthMatchedIsoAcc.isAvailable( *lep_itr ) )	    { m_lepton_tag_isTruthMatchedIso.push_back( isTruthMatchedIsoAcc( *lep_itr ) ); }
    	  else {  m_lepton_tag_isTruthMatchedIso .push_back( -1 ); }
    	  if ( isTruthMatchedNonIsoAcc.isAvailable( *lep_itr ) )    { m_lepton_tag_isTruthMatchedNonIso.push_back(  isTruthMatchedNonIsoAcc( *lep_itr ) ); }
    	  else {  m_lepton_tag_isTruthMatchedNonIso.push_back( -1 ); }
    	  if ( isTruthMatchedSecondaryAcc.isAvailable( *lep_itr ) ) { m_lepton_tag_isTruthMatchedSecondary.push_back(  isTruthMatchedSecondaryAcc( *lep_itr ) ); }
    	  else {  m_lepton_tag_isTruthMatchedSecondary.push_back( -1 ); }
    	  if ( isTruthMatchedNoProdVtxAcc.isAvailable( *lep_itr ) ) { m_lepton_tag_isTruthMatchedNoProdVtx.push_back(  isTruthMatchedNoProdVtxAcc( *lep_itr ) ); }
    	  else {  m_lepton_tag_isTruthMatchedNoProdVtx.push_back( -1 ); }
    	  if ( isTruthMatchedOtherAcc.isAvailable( *lep_itr ) )     { m_lepton_tag_isTruthMatchedOther.push_back(  isTruthMatchedOtherAcc( *lep_itr ) ); }
    	  else {  m_lepton_tag_isTruthMatchedOther.push_back( -1 ); }
    	  if ( isTruthMatchedUnknownAcc.isAvailable( *lep_itr ) )   { m_lepton_tag_isTruthMatchedUnknown.push_back(  isTruthMatchedUnknownAcc( *lep_itr ) ); }
    	  else { m_lepton_tag_isTruthMatchedUnknown.push_back( -1 ); }
    	  if ( isChFlipAcc.isAvailable( *lep_itr ) )		    { m_lepton_tag_isChFlip.push_back(  isChFlipAcc( *lep_itr ) ); }
    	  else { m_lepton_tag_isChFlip.push_back( -1 ); }
    	  if ( isBremAcc.isAvailable( *lep_itr ) )		    { m_lepton_tag_isBrem.push_back(  isBremAcc( *lep_itr ) ); }
    	  else { m_lepton_tag_isBrem.push_back( -1 ); }
    	  if ( truthMatchTypeAcc.isAvailable( *lep_itr ) )	    { m_lepton_tag_truthType.push_back(  truthMatchTypeAcc( *lep_itr ) ); }
    	  else { m_lepton_tag_truthType.push_back( -1 ); }
    	  if ( truthMatchOriginAcc.isAvailable( *lep_itr ) )	    { m_lepton_tag_truthOrigin.push_back(  truthMatchOriginAcc( *lep_itr ) ); }
    	  else { m_lepton_tag_truthOrigin.push_back( -1 ); }
    	  if ( truthMatchStatusAcc.isAvailable( *lep_itr ) )	    { m_lepton_tag_truthStatus.push_back(  truthMatchStatusAcc( *lep_itr ) ); }
    	  else { m_lepton_tag_truthStatus.push_back( -1 ); }
    	} else {
    	// fill lepton PROBE variables
    	  m_lepton_probe_pt.push_back(  lep_itr->pt() );
    	  m_lepton_probe_eta.push_back(  lep_itr->eta() );
    	  m_lepton_probe_flavour.push_back( lep_flav );
	  m_lepton_probe_charge.push_back( lep_charge );
    	  if ( isTightAcc.isAvailable( *lep_itr ) )		    { m_lepton_probe_isTight.push_back(  isTightAcc( *lep_itr ) ); }
    	  else { m_lepton_probe_isTight.push_back(  -1 ); }
    	  if ( isTruthMatchedAcc.isAvailable( *lep_itr ) )	    { m_lepton_probe_isTruthMatched.push_back( isTruthMatchedAcc( *lep_itr ) ); }
    	  else {  m_lepton_probe_isTruthMatched.push_back( -1 ); }
    	  if ( isTruthMatchedIsoAcc.isAvailable( *lep_itr ) )	    { m_lepton_probe_isTruthMatchedIso.push_back( isTruthMatchedIsoAcc( *lep_itr ) ); }
    	  else {  m_lepton_probe_isTruthMatchedIso.push_back( -1 ); }
    	  if ( isTruthMatchedNonIsoAcc.isAvailable( *lep_itr ) )    { m_lepton_probe_isTruthMatchedNonIso.push_back(  isTruthMatchedNonIsoAcc( *lep_itr ) ); }
    	  else {  m_lepton_probe_isTruthMatchedNonIso.push_back( -1 ); }
    	  if ( isTruthMatchedSecondaryAcc.isAvailable( *lep_itr ) ) { m_lepton_probe_isTruthMatchedSecondary.push_back(  isTruthMatchedSecondaryAcc( *lep_itr ) ); }
    	  else {  m_lepton_probe_isTruthMatchedSecondary.push_back( -1 ); }
    	  if ( isTruthMatchedNoProdVtxAcc.isAvailable( *lep_itr ) ) { m_lepton_probe_isTruthMatchedNoProdVtx.push_back(  isTruthMatchedNoProdVtxAcc( *lep_itr ) ); }
    	  else {  m_lepton_probe_isTruthMatchedNoProdVtx.push_back( -1 ); }
    	  if ( isTruthMatchedOtherAcc.isAvailable( *lep_itr ) )     { m_lepton_probe_isTruthMatchedOther.push_back(  isTruthMatchedOtherAcc( *lep_itr ) ); }
    	  else {  m_lepton_probe_isTruthMatchedOther.push_back( -1 ); }
    	  if ( isTruthMatchedUnknownAcc.isAvailable( *lep_itr ) )   { m_lepton_probe_isTruthMatchedUnknown.push_back(  isTruthMatchedUnknownAcc( *lep_itr ) ); }
    	  else { m_lepton_probe_isTruthMatchedUnknown.push_back( -1 ); }
    	  if ( isChFlipAcc.isAvailable( *lep_itr ) )		    { m_lepton_probe_isChFlip.push_back(  isChFlipAcc( *lep_itr ) ); }
    	  else { m_lepton_probe_isChFlip.push_back( -1 ); }
    	  if ( isBremAcc.isAvailable( *lep_itr ) )		    { m_lepton_probe_isBrem.push_back(  isBremAcc( *lep_itr ) ); }
    	  else { m_lepton_probe_isBrem.push_back( -1 ); }
    	  if ( truthMatchTypeAcc.isAvailable( *lep_itr ) )	    { m_lepton_probe_truthType.push_back(  truthMatchTypeAcc( *lep_itr ) ); }
    	  else { m_lepton_probe_truthType.push_back( -1 ); }
    	  if ( truthMatchOriginAcc.isAvailable( *lep_itr ) )	    { m_lepton_probe_truthOrigin.push_back(  truthMatchOriginAcc( *lep_itr ) ); }
    	  else { m_lepton_probe_truthOrigin.push_back( -1 ); }
    	  if ( truthMatchStatusAcc.isAvailable( *lep_itr ) )	    { m_lepton_probe_truthStatus.push_back(  truthMatchStatusAcc( *lep_itr ) ); }
    	  else { m_lepton_probe_truthStatus.push_back( -1 ); }
    	}
    	
      }
      else   { m_lepton_isTag.push_back(-1); }  

      ++m_nlep;
  }
  
}  
