/*
 *
 * File     : Core/CutFlowHist.C
 * Authors  : KG <Kong.Guan.Tan@cern.ch>
 *
 * To be compiled within PyROOT with ACLiC.
 * Manual ROOT Histogram is needed to store the CutFlow numbers because ROOT has a bug with their bin labels!
 *
 */

#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>

#include "TTree.h"
#include "TH1D.h"

namespace AnalysisFramework
{
namespace CutFlows
{

class CutFlowHist : public TH1D
{
public:
    CutFlowHist() : TH1D() { }
    CutFlowHist(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup) : TH1D(name, title, nbinsx, xlow, xup)
    {
        labels.reserve(nbinsx+2);
        bin_array.reserve(nbinsx+2);
//         labels.resize(nbinsx+2);
    }

    std::vector<std::string> labels;
    std::map<std::string, int> labels_map;
    std::vector<double> bin_array;

    void resizeBins(int newsize)
    {
        labels.resize(newsize);
        bin_array.resize(newsize);
    }

    virtual Int_t GetNbinsX()
    {
        return labels.size();
    }

    virtual Double_t GetBinContent(Int_t bin)
    {
        if ((int)(bin_array.size())< bin+1)
            resizeBins(std::max(bin+1, (int)(bin_array.size()*2)));
        return bin_array[bin];
    }

    std::string GetBinLabel(int bin)
    {
        if ((int)(labels.size())< bin+1)
            resizeBins(std::max(bin+1, (int)(labels.size()*2)));
        return labels[bin];
    }

    int GetBinNumber(std::string name)
    {
        std::map<std::string, int>::iterator it = labels_map.find(name);
        if (it == labels_map.end())
            return -1;
        return it->second;
    }

    void SetBinLabel(int bin, std::string name)
    {
        if ((int)(labels.size()) < bin+1)
            resizeBins(std::max(bin+1, (int)(labels.size()*2)));
        labels[bin] = name;
        labels_map[name] = bin;
    }

    virtual void SetBinContent(Int_t bin, double value)
    {
        if ((int)(bin_array.size())< bin+1)
            resizeBins(std::max(bin+1, (int)(bin_array.size()*2)));
        bin_array[bin] = value;
    }

    Long64_t Merge(TCollection *hlist)
    {
        // Uncomment to copy the labels from TAxis
        /*
        if (GetBinLabel(1).size() == 0)
        {
            TAxis *axis = GetXaxis();
            for (unsigned int i=1; i<=labels.size(); i++)
            {
                if (axis->GetBinLabel(i)[0] == '\0')
                    break;
                labels.at(i) = axis->GetBinLabel(i);
            }
        }*/

        if (hlist)
        {
            //TAxis *axis = GetXaxis();
            int unusedbin = 0;
            for (unsigned int i=1; i <= labels.size()+1; i++)
            {
                std::string label = GetBinLabel(i);
                //cout << "Bin: " << i << "Name: " << label << endl;
                if (label.size() == 0)
                {
                    unusedbin = i;
                    break;
                }
            }

            CutFlowHist *xh = 0;
            TIter nxh(hlist);
            while ( (xh = (CutFlowHist *) nxh()) )
            {
                //TAxis *extaxis = xh->GetXaxis();
                if (xh->labels.size() > labels.size())
                    resizeBins(xh->labels.size()*2);
                for (unsigned int i=1; i < xh->labels.size()+1; i++)
                {
                    std::string extlabel = xh->GetBinLabel(i);
                    if (extlabel.size() == 0) break;
                    int bin = GetBinNumber(extlabel);
                    if (bin<0)
                    {
                        SetBinLabel(unusedbin, extlabel);
                        bin = unusedbin++;
                    }
                    SetBinContent(bin, GetBinContent(bin) + xh->GetBinContent(i));
                }
            }
        }
        return (Long64_t)GetEntries();
    }
    ClassDef(CutFlowHist,2)
};

}
}
