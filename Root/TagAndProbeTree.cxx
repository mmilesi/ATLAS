#include "HTopMultilepAnalysis/TagAndProbeTree.h"

TagAndProbeTree :: TagAndProbeTree(xAOD::TEvent* event, TTree* tree, TFile* file, const float units) : 
  HelpTreeBase(event, tree, file, units)
{
  Info("TagAndProbeTree", "Creating output TTree");
}

TagAndProbeTree :: ~TagAndProbeTree() 
{
}

void TagAndProbeTree::AddEventUser(const std::string detailStrUser)
{
  // event variables
  m_tree->Branch("isMC",              &m_is_mc,              "isMC/I");  
  m_tree->Branch("nMediumBjets",      &m_nBjetsMedium,       "nMediumBjets/I");
  m_tree->Branch("isDileptonOS",      &m_isDileptonOS,       "isDileptonOS/I"); 
  m_tree->Branch("isDileptonSS",      &m_isDileptonSS,       "isDileptonSS/I"); 
  m_tree->Branch("isNonTightEvent",   &m_isNonTightEvent,    "isNonTightEvent/I");
  m_tree->Branch("isProbeElEvent",    &m_isProbeElEvent,     "isProbeElEvent/I");
  m_tree->Branch("isProbeMuEvent",    &m_isProbeMuEvent,     "isProbeMuEvent/I");
}

void TagAndProbeTree::AddJetsUser(const std::string detailStrUser)
{
  // jet variables
  m_tree->Branch("jet_m",       &m_jet_m);    
  m_tree->Branch("jet_clean",	&m_jet_clean);
}

void TagAndProbeTree::AddMuonsUser(const std::string detailStrUser)
{
  // muon variables  
  m_tree->Branch("muon_isTight",                   &m_muon_isTight); 
  m_tree->Branch("muon_isTrigMatched",             &m_muon_isTrigMatched); 
  m_tree->Branch("muon_isTag",                     &m_muon_isTag); 
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
  // muon TAG variables 
  m_tree->Branch("muon_tag_pt",	                   &m_muon_tag_pt);
  m_tree->Branch("muon_tag_eta",	           &m_muon_tag_eta);
  m_tree->Branch("muon_tag_isTight",	           &m_muon_tag_isTight); 
  m_tree->Branch("muon_tag_isChFlip",	           &m_muon_tag_isChFlip);     
  m_tree->Branch("muon_tag_isBrem",	           &m_muon_tag_isBrem);   
  m_tree->Branch("muon_tag_isTruthMatchedIso",	   &m_muon_tag_isTruthMatchedIso); 
  // muon PROBE variables 
  m_tree->Branch("muon_probe_pt",	           &m_muon_probe_pt);
  m_tree->Branch("muon_probe_eta",	           &m_muon_probe_eta);
  m_tree->Branch("muon_probe_isTight",	           &m_muon_probe_isTight); 
  m_tree->Branch("muon_probe_isChFlip",	           &m_muon_probe_isChFlip);	
  m_tree->Branch("muon_probe_isBrem",	           &m_muon_probe_isBrem); 
  m_tree->Branch("muon_probe_isTruthMatchedIso",   &m_muon_probe_isTruthMatchedIso); 
}

void TagAndProbeTree::AddElectronsUser(const std::string detailStrUser)
{   		     
  // electron variables  
  m_tree->Branch("el_calo_eta",                  &m_electron_calo_eta);
  m_tree->Branch("el_crack",                     &m_electron_crack);
  m_tree->Branch("el_isTight",                   &m_electron_isTight); 
  m_tree->Branch("el_isTrigMatched",             &m_electron_isTrigMatched);
  m_tree->Branch("el_isTag",                     &m_electron_isTag); 
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
  // electron TAG variables 
  m_tree->Branch("el_tag_pt",	                 &m_electron_tag_pt);
  m_tree->Branch("el_tag_eta",	                 &m_electron_tag_eta);
  m_tree->Branch("el_tag_isTight",	         &m_electron_tag_isTight); 
  m_tree->Branch("el_tag_isChFlip",	         &m_electron_tag_isChFlip);	
  m_tree->Branch("el_tag_isBrem",	         &m_electron_tag_isBrem);	
  m_tree->Branch("el_tag_isTruthMatchedIso",     &m_electron_tag_isTruthMatchedIso); 
  // electron PROBE variables 
  m_tree->Branch("el_probe_pt",	                 &m_electron_probe_pt);
  m_tree->Branch("el_probe_eta",	         &m_electron_probe_eta);
  m_tree->Branch("el_probe_isTight",	         &m_electron_probe_isTight); 
  m_tree->Branch("el_probe_isChFlip",	         &m_electron_probe_isChFlip);	  
  m_tree->Branch("el_probe_isBrem",	         &m_electron_probe_isBrem);	
  m_tree->Branch("el_probe_isTruthMatchedIso",   &m_electron_probe_isTruthMatchedIso); 

}

void TagAndProbeTree::ClearEventUser() {}
  
void TagAndProbeTree::ClearMuonsUser() {  
  // muon variables 
  m_muon_isTight.clear(); 
  m_muon_isTrigMatched.clear(); 
  m_muon_isTag.clear(); 
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
  m_muon_tag_pt.clear();
  m_muon_tag_eta.clear();
  m_muon_tag_isTight.clear(); 
  m_muon_tag_isChFlip.clear();     
  m_muon_tag_isBrem.clear();	
  m_muon_tag_isTruthMatchedIso.clear(); 
  m_muon_probe_pt.clear();
  m_muon_probe_eta.clear();
  m_muon_probe_isTight.clear(); 
  m_muon_probe_isChFlip.clear();     
  m_muon_probe_isBrem.clear(); 
  m_muon_probe_isTruthMatchedIso.clear(); 
}
  
void TagAndProbeTree::ClearElectronsUser() {  
  // electron variables
  m_electron_calo_eta.clear();
  m_electron_crack.clear();
  m_electron_isTight.clear(); 
  m_electron_isTrigMatched.clear(); 
  m_electron_isTag.clear(); 
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
  m_electron_tag_pt.clear();
  m_electron_tag_eta.clear();
  m_electron_tag_isTight.clear(); 
  m_electron_tag_isChFlip.clear();     
  m_electron_tag_isBrem.clear();	
  m_electron_tag_isTruthMatchedIso.clear(); 
  m_electron_probe_pt.clear();
  m_electron_probe_eta.clear();
  m_electron_probe_isTight.clear(); 
  m_electron_probe_isChFlip.clear();     
  m_electron_probe_isBrem.clear(); 
  m_electron_probe_isTruthMatchedIso.clear();   
}


void TagAndProbeTree::ClearJetsUser() {
  // jet variables 
  m_jet_m.clear();
  m_jet_clean.clear();
}  


void TagAndProbeTree::FillEventUser( const xAOD::EventInfo* eventInfo ) { 
  
  m_is_mc              =  ( eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) );
  m_nBjetsMedium       =  ( eventInfo->isAvailable< unsigned int >( "nBjetsMedium" ) ) ?  eventInfo->auxdecor< unsigned int >( "nBjetsMedium" )   :  -1;
  m_isDileptonOS       =  ( eventInfo->isAvailable< char >( "isDileptonOS" ) )         ?  eventInfo->auxdecor< char >( "isDileptonOS" )           :  -1;
  m_isDileptonSS       =  ( eventInfo->isAvailable< char >( "isDileptonSS" ) )         ?  eventInfo->auxdecor< char >( "isDileptonSS" )           :  -1;
  m_isNonTightEvent    =  ( eventInfo->isAvailable< char >( "isNonTightEvent" ) )      ?  eventInfo->auxdecor< char >( "isNonTightEvent" )        :  -1;
  m_isProbeElEvent     =  ( eventInfo->isAvailable< char >( "isProbeElEvent" ) )       ?  eventInfo->auxdecor< char >( "isProbeElEvent" )         :  -1;
  m_isProbeMuEvent     =  ( eventInfo->isAvailable< char >( "isProbeMuEvent" ) )       ?  eventInfo->auxdecor< char >( "isProbeMuEvent" )         :  -1;
  
}

void TagAndProbeTree::FillJetsUser( const xAOD::Jet* jet ) {
  
  static SG::AuxElement::Accessor< char > isCleanAcc("cleanJet");
  
  if( isCleanAcc.isAvailable( *jet ) ) { m_jet_clean.push_back( isCleanAcc( *jet ) ); } 
  else { m_jet_clean.push_back(-1); }
  
  m_jet_m.push_back( jet->m() ); 

}

void TagAndProbeTree::FillMuonsUser( const xAOD::Muon* muon ) {

  static SG::AuxElement::Accessor< char > isTightAcc("isTight");
  static SG::AuxElement::Accessor< char > isTrigMatchedAcc("isTrigMatched");
  static SG::AuxElement::Accessor< char > isTagAcc("isTag");
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
  
  if ( isTightAcc.isAvailable( *muon ) ) { m_muon_isTight.push_back( isTightAcc( *muon ) ); }
  else   { m_muon_isTight.push_back(-1); }
  if ( isTrigMatchedAcc.isAvailable( *muon ) ) { m_muon_isTrigMatched.push_back( isTrigMatchedAcc( *muon ) ); }
  else   { m_muon_isTrigMatched.push_back(-1); }
  if ( isTagAcc.isAvailable( *muon ) ) { 
    
    m_muon_isTag.push_back( isTagAcc( *muon ) ); 
  
    // fill muon TAG variables
    if ( isTagAcc( *muon ) ){
      m_muon_tag_pt.push_back( muon->pt() );
      m_muon_tag_eta.push_back( muon->eta() );
      if ( isTightAcc.isAvailable( *muon ) )           { m_muon_tag_isTight.push_back( isTightAcc( *muon ) ); }
      if ( isChFlipAcc.isAvailable( *muon ) )          { m_muon_tag_isChFlip.push_back( isChFlipAcc( *muon ) ); }     
      if ( isBremAcc.isAvailable( *muon ) )            { m_muon_tag_isBrem.push_back( isBremAcc( *muon ) ); }    
      if ( isTruthMatchedIsoAcc.isAvailable( *muon ) ) { m_muon_tag_isTruthMatchedIso.push_back( isTruthMatchedIsoAcc( *muon ) ); } 
    } else {
    // fill muon PROBE variables
      m_muon_probe_pt.push_back( muon->pt() );
      m_muon_probe_eta.push_back( muon->eta() );
      if ( isTightAcc.isAvailable( *muon ) )           { m_muon_probe_isTight.push_back( isTightAcc( *muon ) ); }
      if ( isChFlipAcc.isAvailable( *muon ) )          { m_muon_probe_isChFlip.push_back( isChFlipAcc( *muon ) ); }     
      if ( isBremAcc.isAvailable( *muon ) )            { m_muon_probe_isBrem.push_back( isBremAcc( *muon ) ); }    
      if ( isTruthMatchedIsoAcc.isAvailable( *muon ) ) { m_muon_probe_isTruthMatchedIso.push_back( isTruthMatchedIsoAcc( *muon ) ); } 
    }
    
  }
  else   { m_muon_isTag.push_back(-1); }  
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
}

void TagAndProbeTree::FillElectronsUser( const xAOD::Electron* electron ) {
  
  static SG::AuxElement::Accessor< char > isTightAcc("isTight");
  static SG::AuxElement::Accessor< char > isTrigMatchedAcc("isTrigMatched");
  static SG::AuxElement::Accessor< char > isTagAcc("isTag");  
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
  
  if ( isTightAcc.isAvailable( *electron ) ) { m_electron_isTight.push_back( isTightAcc( *electron ) ); }
  else   { m_electron_isTight.push_back(-1); }
  if ( isTrigMatchedAcc.isAvailable( *electron ) ) { m_electron_isTrigMatched.push_back( isTrigMatchedAcc( *electron ) ); }
  else   { m_electron_isTrigMatched.push_back(-1); }
  if ( isTagAcc.isAvailable( *electron ) ) { 
  
    m_electron_isTag.push_back( isTagAcc( *electron ) ); 
  
    // fill electron TAG variables
    if ( isTagAcc( *electron ) ){
      m_electron_tag_pt.push_back( electron->pt() );
      m_electron_tag_eta.push_back( electron->eta() );
      if ( isTightAcc.isAvailable( *electron ) )           { m_electron_tag_isTight.push_back( isTightAcc( *electron ) ); }
      if ( isChFlipAcc.isAvailable( *electron ) )          { m_electron_tag_isChFlip.push_back( isChFlipAcc( *electron ) ); }     
      if ( isBremAcc.isAvailable( *electron ) )            { m_electron_tag_isBrem.push_back( isBremAcc( *electron ) ); }    
      if ( isTruthMatchedIsoAcc.isAvailable( *electron ) ) { m_electron_tag_isTruthMatchedIso.push_back( isTruthMatchedIsoAcc( *electron ) ); } 
    } else {
    // fill electron PROBE variables
      m_electron_probe_pt.push_back( electron->pt() );
      m_electron_probe_eta.push_back( electron->eta() );
      if ( isTightAcc.isAvailable( *electron ) )           { m_electron_probe_isTight.push_back( isTightAcc( *electron ) ); }
      if ( isChFlipAcc.isAvailable( *electron ) )          { m_electron_probe_isChFlip.push_back( isChFlipAcc( *electron ) ); }     
      if ( isBremAcc.isAvailable( *electron ) )            { m_electron_probe_isBrem.push_back( isBremAcc( *electron ) ); }    
      if ( isTruthMatchedIsoAcc.isAvailable( *electron ) ) { m_electron_probe_isTruthMatchedIso.push_back( isTruthMatchedIsoAcc( *electron ) ); } 
    }
  
  }
  else   { m_electron_isTag.push_back(-1); }  
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
}

void TagAndProbeTree::FillFatJetsUser( const xAOD::Jet* fatJet )
{
}
