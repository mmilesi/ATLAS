/*
 *
 * ROOT C++ part. Compiled within PyROOT with ACLiC.
 *
 */

#pragma once

#include <vector>
#include <string>
#include <algorithm>

namespace AnalysisFramework
{
namespace CutFlows
{

class WhiteBoard
{
public:
    WhiteBoard()
    {
        debugMode = true;

        isMC            = false;
        isEmbedding     = false;
        isEmbedding11   = false;
        isEmbedding12   = false;
        isEmbedding13   = false;
        isMC11a         = false;
        isMC11b         = false;
        isMC11c         = false;
        isMC12a         = false;
	isMC15          = false;
        isAFII          = false; // ATLFastII
        isZDY           = false; // Z DrellYan
        isWAlpgenPythia = false; // 

        isData    = false;
        isData11  = false;
        isData12  = false;
        isData15  = false;	
        isMuData  = false;
        isElData  = false;
        isTauData = false;

        // Data12 specific flags
        isPeriodC = false;
        isPeriodD = false;

        // Higgs samples specific
        isggF = false;
        higgsMass = 0;

        corrections = 0;
        systematics = 0;
        doDump = false;
        doRemoveUnselected = false;

        // MC Samples
        PeriodBtoD    = 180164;
        PeriodEtoH    = 183003;
        PeriodItoK1   = 185649; // mc11a only
        PeriodFuture  = 185761; // mc11a only
        PeriodIJK     = 186169; // mc11b only
        PeriodLM      = 189751; // mc11b only

        // Data (Start_ to _End)
        PeriodA_      = 177531; _PeriodA      = 177965;
        PeriodB_      = 177986; _PeriodB      = 178109;
        PeriodD_      = 179710; _PeriodD      = 180481;
        PeriodE_      = 180614; _PeriodE      = 180776;
        PeriodF_      = 182013; _PeriodF      = 182519;
        PeriodG_      = 182726; _PeriodG      = 183462;
        PeriodH_      = 183544; _PeriodH      = 184169;
        PeriodI_      = 185353; _PeriodI      = 186493;
        PeriodJ_      = 186516; _PeriodJ      = 186755;
        PeriodK1_     = 186873; _PeriodK1     = 186934; // needed for mc11a
        PeriodKRest_  = 186965; _PeriodKRest  = 187815; // needed for mc11a
        PeriodK_      = 186873; _PeriodK      = 187815;
        PeriodL_      = 188902; _PeriodL      = 190343;
        PeriodM_      = 190503; _PeriodM      = 191933;

        resetToDefaultValues();
    }

    bool debugMode;

    bool isMC;
    bool isEmbedding;
    bool isEmbedding11;
    bool isEmbedding12;
    bool isEmbedding13;
    bool isMC11a;
    bool isMC11b;
    bool isMC11c;
    bool isMC12a;
    bool isMC15;
    bool isAFII;
    bool isZDY;
    bool isWAlpgenPythia;

    bool isVBFFiltered;
    bool isggF;
    int higgsMass;

    bool isData;
    bool isData11;
    bool isData12;
    bool isData15;
    bool isElData;
    bool isMuData;
    bool isTauData;
    bool isPeriodC;
    bool isPeriodD;
   
    std::vector<std::string> *corrections;
    std::vector<std::string> *systematics;
    bool doDump;
    bool doRemoveUnselected;

    unsigned int PeriodBtoD, PeriodEtoH, PeriodItoK1, PeriodFuture, PeriodIJK, PeriodLM;
    unsigned int PeriodA_, _PeriodA, PeriodB_, _PeriodB, PeriodD_, _PeriodD, PeriodE_, _PeriodE, PeriodF_, _PeriodF,
                 PeriodG_, _PeriodG, PeriodH_, _PeriodH, PeriodI_, _PeriodI, PeriodJ_, _PeriodJ, PeriodK1_, _PeriodK1,
                 PeriodKRest_, _PeriodKRest, PeriodK_, _PeriodK, PeriodL_, _PeriodL, PeriodM_, _PeriodM;

    void resetToDefaultValues()
    {
    }

    bool doCorrection(std::string item)
    {
        return std::find(corrections->begin(), corrections->end(), item) != corrections->end();
    }

    bool doSystematic(std::string item)
    {
        return std::find(systematics->begin(), systematics->end(), item) != systematics->end();
    }

}; // End of WhiteBoard class

} // End of CutFlow namespace
} // End of AnalysisFramework namespace
