/*
 *
 * File     : Core/PostLoader.C
 * Authors  : KG <Kong.Guan.Tan@cern.ch>
 *
 * To be compiled within PyROOT with ACLiC.
 * Needed by PyROOT and ROOT after everything compiles
 *
 */

#include "CutFlow_Base.C"

#include<vector>
#include<map>
#include<string>

#ifdef __CINT__
#pragma link C++ class map<string,AnalysisFramework::CutFlows::EventTotal*>+;
#pragma link C++ struct pair<string,AnalysisFramework::CutFlows::EventTotal*>+;
#pragma link C++ class vector<AnalysisFramework::CutFlows::FlowItem*>+;
#else
template class std::map<string,AnalysisFramework::CutFlows::EventTotal*>;
template struct std::pair<string,AnalysisFramework::CutFlows::EventTotal*>;
template class std::vector<AnalysisFramework::CutFlows::FlowItem*>;
#endif
