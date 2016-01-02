/*
 *
 * File     : Core/Loader.C
 * Authors  : KG <Kong.Guan.Tan@cern.ch>
 *
 * To be compiled within PyROOT with ACLiC.
 * Needed by PyROOT and ROOT for subsequent compiles using ACLiC
 *
 */

#include<vector>
#include<string>

#ifdef __CINT__
#pragma link C++ class vector<vector<long> >+;
#pragma link C++ class vector<vector<short> >+;
#pragma link C++ class vector<vector<double> >+;
#pragma link C++ class vector<vector<float> >+;
#pragma link C++ class vector<vector<int> >+;
#pragma link C++ class vector<vector<unsigned int> >+;
#pragma link C++ class vector<vector<unsigned short> >+;
#pragma link C++ class vector<vector<unsigned long> >+;
#pragma link C++ class vector<vector<string> >+;
#pragma link C++ class pair<string,string>+;
#else
template class std::vector<vector<long> >;
template class std::vector<vector<short> >;
template class std::vector<vector<double> >;
template class std::vector<vector<float> >;
template class std::vector<vector<int> >;
template class std::vector<vector<unsigned int> >;
template class std::vector<vector<unsigned short> >;
template class std::vector<vector<unsigned long> >;
template class std::vector<vector<string> >;
template class std::pair<string,string>;
#endif
