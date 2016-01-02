/*
 *
 * File     : Core/Loader.C
 * Authors  : KG <Kong.Guan.Tan@cern.ch>
 *
 * To be compiled within PyROOT with ACLiC.
 * Contains the basic and base classes for CutFlows
 *
 */

#pragma once

#include "../libs/Branches_generated.C"
#include "WhiteBoard.C"
#include <iostream>
#include <map>
#include <string>
#include <algorithm>

#include "TTree.h"

#define MSG(x) \
   std::cout << "\t@@@@@ " << this->name << "\tEventNumber: " << this->ao->evtinfo.EventNumber << "\t" << x << std::endl

#define MSG_(x) \
   std::cout << x

namespace AnalysisFramework
{
namespace CutFlows
{

struct EventTotal
{
    int passed;
    int failed;
    double passedW;
    double failedW;
    EventTotal(): passed(0), failed(0), passedW(0), failedW(0) {}
};

class FlowItem;
struct Passport
{
    double weight;
    bool previousPass;
    std::string streamlet;
    std::vector<std::string> streamletHistory;
    std::vector<FlowItem*> flowHistory;
    std::vector<std::string> systematicsList;
    short systematicsDirection;
    std::vector<TTree*> delayedNomTrees;
    std::vector<TTree*> delayedSysTrees;

    Passport()
    {
        weight = 1.0;
        previousPass = true;
        streamlet = "";
        systematicsDirection = 0;
    }

    void addSubStreamlet(std::string substreamlet)
    {
        if (streamlet.size()>0)
            streamlet += "&"+substreamlet;
        else
            streamlet += substreamlet;
        streamletHistory.push_back(substreamlet);
    }

    bool hasSubstreamlet(std::string substreamlet)
    {
        return std::find(streamletHistory.begin(), streamletHistory.end(), substreamlet) != streamletHistory.end();
    }

    bool checkSystematics(std::vector<std::string> &systematicsPartners)
    {
        if (systematicsList.size() == 0)
            return true;
        else if (systematicsPartners.size() == 0)
            return false;
        else
        {
            bool found = false;
            for (unsigned int i = 0; i < systematicsPartners.size(); i++)
            {
                if (std::find(systematicsList.begin(), systematicsList.end(), systematicsPartners.at(i)) != systematicsList.end())
                    found = true;
            }
            return found;
        }
    }
};

class FlowItem
{
public:
    std::string name;
    std::string substreamlet;
    WhiteBoard *wb;
    AnalysisFramework::Branches::AllObjects *ao;
    std::map<std::string, EventTotal*> totalEvents;
    std::vector<std::string> streamlets;
    TTree *treeNom;
    TTree *treeNomDelayed;

    bool verbose;
    bool running;
    Passport runningPassport;

    bool continueWhenFail;
    bool hasPrevious;
    bool hasNext;
    FlowItem* previous;
    FlowItem* next;

    FlowItem(std::string name_, std::string substreamlet_ = "")
    {
        name = name_;
        substreamlet = substreamlet_;
        treeNom = 0;
        treeNomDelayed = 0;
        wb = 0;
        ao = 0;

        verbose = false;
        running = false;
        continueWhenFail = true;
        hasPrevious = false;
        hasNext = false;
        previous = 0;
        next = 0;
    }

    virtual ~FlowItem()
    {
    }

    virtual void incrementTotalEvents(double weight, std::string streamlet, bool effectivePass)
    {
        EventTotal *totalEvents_ = 0;
        std::map<std::string, EventTotal*>::iterator it = totalEvents.find(streamlet);
        if (it == totalEvents.end())
        {
            totalEvents_ = new EventTotal();
            totalEvents[streamlet] = totalEvents_;
            streamlets.push_back(streamlet);
        }
        else
        {
            totalEvents_ = it->second;
        }

        if (effectivePass)
        {
            totalEvents_->passed++;
            totalEvents_->passedW += weight;;
        }
        else
        {
            totalEvents_->failed++;
            totalEvents_->failedW += weight;;
        }
    }

    int getFirstSelected(std::vector<bool> *selected, std::vector<int> *order = 0)
    {
        for (unsigned int s = 0; s < selected->size(); s++)
        {
            int r = s;
            if (order) if (order->size()>s) r = order->at(s);
            if (selected->at(r)) return r;
        }
        return -1;
    }

    int getSecondSelected(std::vector<bool> *selected, std::vector<int> *order = 0)
    {
        int first = -1;
        for (unsigned int s = 0; s < selected->size(); s++)
        {
            int r = s;
            if (order) if (order->size()>s) r = order->at(s);
            if (selected->at(r))
            {
                if (first == -1)
                    first = r;
                else
                    return r;
            }
        }
        return -1;
    }

    int getThirdSelected(std::vector<bool> *selected, std::vector<int> *order = 0)
    {
        int first = -1;
        int second = -1;
        for (unsigned int s = 0; s < selected->size(); s++)
        {
            int r = s;
            if (order) if (order->size()>s) r = order->at(s);
            if (selected->at(r))
            {
                if (first == -1)
                    first = r;
                else if (second == -1)
                    second = r;
                else
                    return r;
            }
        }
        return -1;
    }

    int countSelected(std::vector<bool> *selected)
    {
        int count = 0;
        for (unsigned int s = 0; s < selected->size(); s++)
        {
            if (selected->at(s)) count++;
        }
        return count;
    }

    template<class T>
    void orderBy(std::vector<T> *measureable, std::vector<int> *order, bool descending=true)
    {
        order->clear();
        order->reserve(measureable->size());
        for (unsigned int m = 0; m < measureable->size(); m++)
        {
            unsigned int o = 0;
            for (; o < m; o++)
            {
                if (descending && measureable->at(m) > measureable->at(order->at(o)))
                    break;
                if (!descending && measureable->at(m) < measureable->at(order->at(o)))
                    break;
            }
            if (o==m)
                order->push_back(m);
            else
                order->insert(order->begin()+o, m);
        }
        return;
    }

    bool hasFlow(std::vector<AnalysisFramework::CutFlows::FlowItem*> &flowVector, std::string flowname)
    {
        for (unsigned int f = 0; f < flowVector.size(); f++)
            if (flowVector.at(f)->name == flowname)
                return true;
        return false;
    }

    FlowItem* getFlow(std::vector<AnalysisFramework::CutFlows::FlowItem*> &flowVector, std::string flowname)
    {
        for (unsigned int f = 0; f < flowVector.size(); f++)
            if (flowVector.at(f)->name == flowname)
                return flowVector.at(f);
        return NULL;
    }

    virtual bool process(Passport passport = Passport())
    {
        // Make a copy of the passport
        runningPassport = passport;

        // Automatically add substreamlets if defined
        if (substreamlet.size() > 0)
            runningPassport.addSubStreamlet(substreamlet);
        runningPassport.flowHistory.push_back(this);

        // Execute the FlowItem
        running = true;
        bool pass = execute(runningPassport);
        running = false;

        // Increments the total events is if passes
        bool effectivePass = pass && runningPassport.previousPass;
        incrementTotalEvents(runningPassport.weight, runningPassport.streamlet, effectivePass);

        // Stops the streamlet if fails the flowitem and "continueWhenFail" option is set to false
        if (!pass && !continueWhenFail)
            return false;

        // Fill the nominal trees, and schedule delayed nominal trees for filling at the end of the streamlet
        if (wb->doDump && !runningPassport.systematicsDirection)
        {
            if (treeNom)
            {
                if (wb->doRemoveUnselected)
                    ao->removeUnselected();
                treeNom->Fill();
            }
            else if (treeNomDelayed)
            {
                runningPassport.delayedNomTrees.push_back(treeNomDelayed);
            }
        }

        // If has next item chained, run it, otherwise, fill the delayed trees
        bool success = false;
        if (hasNext)
        {
            runningPassport.previousPass = effectivePass;
            success = next->process(runningPassport);
        }
        else
        {
            success = pass;
            if (wb->doDump && !runningPassport.systematicsDirection && runningPassport.delayedNomTrees.size()>0)
            {
                if (wb->doRemoveUnselected)
                    ao->removeUnselected();
                for (unsigned int i=0; i<runningPassport.delayedNomTrees.size(); i++)
                    runningPassport.delayedNomTrees.at(i)->Fill();
            }
            else if (wb->doDump && runningPassport.systematicsDirection && runningPassport.delayedSysTrees.size()>0)
            {
                if (wb->doRemoveUnselected)
                    ao->removeUnselected();
                for (unsigned int i=0; i<runningPassport.delayedSysTrees.size(); i++)
                    runningPassport.delayedSysTrees.at(i)->Fill();
            }
        }
        return success;
    }

    virtual bool execute(Passport &)
    {
        // Function to be overloaded
        return true;
    }
}; // End of FlowItem class

class CorrectionItem : public FlowItem
{
public:
    CorrectionItem(std::string name_, std::string substreamlet_ = "") : FlowItem(name_, substreamlet_)
    {
    }
}; // End of CorrectionItem class

class CutItem : public FlowItem
{
public:
    CutItem(std::string name_, std::string substreamlet_ = "") : FlowItem(name_, substreamlet_)
    {
        continueWhenFail = false;
    }
}; // End of CutItem class


class ForkBase : public FlowItem
{
public:
    ForkBase(std::string name_, std::string substreamlet_ = "") : FlowItem(name_, substreamlet_)
    {
        keepModifiedValues = false;
    }

    bool keepModifiedValues;
    std::vector<std::string> nextNames;
    std::vector<FlowItem*> nextList;

    virtual void backupOriginalValues()
    {
    }

    virtual void restoreOriginalValues()
    {
    }

    virtual void saveValues()
    {
    }

    virtual void loadValues()
    {
    }

    template<class T>
    void copyVector(std::vector<T> *from, std::vector<T> *to)
    {
        *to = *from;
    }

    template<class T>
    void copyVectorVector(std::vector<std::vector<T> > *from, std::vector<std::vector<T> > *to)
    {
        *to = *from;
    }

    void addFlowItem(FlowItem *flow)
    {
        nextList.push_back(flow);
    }

    void addFlowItem(FlowItem *flow, std::string substreamlet_)
    {
        flow->substreamlet = substreamlet_;
        nextList.push_back(flow);
    }

    bool process(Passport passport = Passport())
    {
        // Make a copy of the passport
        runningPassport = passport;

        // Automatically add substreamlets if defined
        if (substreamlet.size() > 0)
            runningPassport.addSubStreamlet(substreamlet);
        runningPassport.flowHistory.push_back(this);

        // Execute the ForkItem (but not yet the streamlets associated to the fork)
        running = true;
        bool pass = execute(runningPassport);
        running = false;

        // Increments the total events is if passes
        bool effectivePass = pass && runningPassport.previousPass;
        incrementTotalEvents(runningPassport.weight, runningPassport.streamlet, effectivePass);

        // Stops the streamlet if fails the flowitem and "continueWhenFail" option is set to false
        if (!pass && !continueWhenFail)
            return false;

        // Fill the nominal trees, and schedule delayed nominal trees for filling at the end of the streamlet
        if (wb->doDump && !runningPassport.systematicsDirection)
        {
            if (treeNom)
            {
                if (wb->doRemoveUnselected)
                    ao->removeUnselected();
                treeNom->Fill();
            }
            else if (treeNomDelayed)
            {
                runningPassport.delayedNomTrees.push_back(treeNomDelayed);
            }
        }

        // If has next item chained, run it, otherwise, fill the delayed trees
        bool success = false;
        if (nextList.size() > 0)
        {
            runningPassport.previousPass = effectivePass;
            running = true;
            backupOriginalValues();
            running = false;
            for (unsigned int i = 0; i < nextList.size(); i++)
            {
                if (i > 0)
                {
                    running = true;
                    restoreOriginalValues();
                    running = false;
                }

                Passport tempPassport = runningPassport;
                if (nextList.at(i)->process(tempPassport))
                    success = true;
            }
        }
        else
        {
            success = pass;
            if (wb->doDump && !runningPassport.systematicsDirection && runningPassport.delayedNomTrees.size()>0)
            {
                if (wb->doRemoveUnselected)
                    ao->removeUnselected();
                for (unsigned int i=0; i<runningPassport.delayedNomTrees.size(); i++)
                    runningPassport.delayedNomTrees.at(i)->Fill();
            }
            else if (wb->doDump && runningPassport.systematicsDirection && runningPassport.delayedSysTrees.size()>0)
            {
                if (wb->doRemoveUnselected)
                    ao->removeUnselected();
                for (unsigned int i=0; i<runningPassport.delayedSysTrees.size(); i++)
                    runningPassport.delayedSysTrees.at(i)->Fill();
            }
        }
        return success;
    }
}; // End of ForkBase class

class ForkLoopBase : public FlowItem
{
public:
    ForkLoopBase(std::string name_, std::string substreamlet_ = "") : FlowItem(name_, substreamlet_)
    {
        keepModifiedValues = false;
    }

    bool keepModifiedValues;

    virtual void backupOriginalValues()
    {
    }

    virtual void restoreOriginalValues()
    {
    }

    virtual void saveValues()
    {
    }

    virtual void loadValues()
    {
    }

    virtual bool nextIteration(Passport &)
    {
        // Return true if the loop continues
        return false;
    }

    template<class T>
    void copyVector(std::vector<T> *from, std::vector<T> *to)
    {
        *to = *from;
    }

    template<class T>
    void copyVectorVector(std::vector<std::vector<T> > *from, std::vector<std::vector<T> > *to)
    {
        *to = *from;
    }

    bool process(Passport passport = Passport())
    {
        // Make a copy of the passport
        runningPassport = passport;

        // Automatically add substreamlets if defined
        if (substreamlet.size() > 0)
            runningPassport.addSubStreamlet(substreamlet);
        runningPassport.flowHistory.push_back(this);

        // Execute the ForkLoopItem and prepare/initialise the loops
        running = true;
        bool pass = execute(runningPassport);
        running = false;

        // Increments the total events is if passes
        bool effectivePass = pass && runningPassport.previousPass;
        incrementTotalEvents(runningPassport.weight, runningPassport.streamlet, effectivePass);

        // Stops the streamlet if fails the flowitem and "continueWhenFail" option is set to false
        if (!pass && !continueWhenFail)
            return false;

        // Fill the nominal trees, and schedule delayed nominal trees for filling at the end of the streamlet
        if (wb->doDump && !runningPassport.systematicsDirection)
        {
            if (treeNom)
            {
                if (wb->doRemoveUnselected)
                    ao->removeUnselected();
                treeNom->Fill();
            }
            else if (treeNomDelayed)
            {
                runningPassport.delayedNomTrees.push_back(treeNomDelayed);
            }
        }

        // If has next item chained, run it, otherwise, fill the delayed trees
        bool success = false;
        if (hasNext)
        {
            runningPassport.previousPass = effectivePass;
            Passport tempPassport = runningPassport;

            running = true;
            backupOriginalValues();
            running = false;
            while (nextIteration(tempPassport))
            {
                if (next->process(tempPassport))
                {
                    success = true;
                }
                tempPassport = runningPassport;
                running = true;
                restoreOriginalValues();
                running = false;
            }
        }
        else
        {
            success = pass;
            if (wb->doDump && !runningPassport.systematicsDirection && runningPassport.delayedNomTrees.size()>0)
            {
                if (wb->doRemoveUnselected)
                    ao->removeUnselected();
                for (unsigned int i=0; i<runningPassport.delayedNomTrees.size(); i++)
                    runningPassport.delayedNomTrees.at(i)->Fill();
            }
            else if (wb->doDump && runningPassport.systematicsDirection && runningPassport.delayedSysTrees.size()>0)
            {
                if (wb->doRemoveUnselected)
                    ao->removeUnselected();
                for (unsigned int i=0; i<runningPassport.delayedSysTrees.size(); i++)
                    runningPassport.delayedSysTrees.at(i)->Fill();
            }
        }
        return success;
    }
}; // End of ForkLoopBase class


class SystematicsBase : public FlowItem
{
public:
    SystematicsBase(std::string name_, std::string substreamlet_ = "") : FlowItem(name_, substreamlet_)
    {
        treeUp = 0;
        treeUpDelayed = 0;
        treeDown = 0;
        treeDownDelayed = 0;
        doMC = false;
        doData = false;
        doEmbedding = false;
        doNominalInternal = false;
        doNominalExternal = false;
        doNominalWhenOff = false;
        externalNominalFlow = 0;
    }

    TTree *treeUp;
    TTree *treeUpDelayed;
    TTree *treeDown;
    TTree *treeDownDelayed;
    std::vector<std::string> systematicsPartners;

    bool doMC;
    bool doData;
    bool doEmbedding;

    bool doNominalInternal;
    bool doNominalExternal;
    bool doNominalWhenOff;
    FlowItem* externalNominalFlow;

    virtual void backupOriginalValues()
    {
    }

    virtual void restoreOriginalValues()
    {
    }

    virtual void saveValues()
    {
    }

    virtual void loadValues()
    {
    }

    template<class T>
    void copyVector(std::vector<T> *from, std::vector<T> *to)
    {
        *to = *from;
    }

    template<class T>
    void copyVectorVector(std::vector<std::vector<T> > *from, std::vector<std::vector<T> > *to)
    {
        *to = *from;
    }

    bool processOne(Passport &passport, short systematicsDirection)
    {
        // Make a copy of the passport
        runningPassport = passport;

        // Add the substreamlets and passport flags if systematics
        if (systematicsDirection)
        {
            runningPassport.systematicsDirection = systematicsDirection;
            runningPassport.systematicsList.push_back(name);
            if (systematicsDirection == -1)
                runningPassport.addSubStreamlet(name + "Down");
            else if (systematicsDirection == 1)
                runningPassport.addSubStreamlet(name + "Up");
        }
        //else if (doNominalInternal || doNominalExternal)
        //    runningPassport.addSubStreamlet(name);

        // Execute the flow if systematics or internal nominal flag is set
        bool pass = true;
        if (systematicsDirection || doNominalInternal)
        {
            running = true;
            pass = executeSystematics(runningPassport, systematicsDirection);
            running = false;
        }
        else if (doNominalExternal && externalNominalFlow)
        {
            pass = externalNominalFlow->execute(runningPassport);
        }

        // Increments the total events is if passes
        bool effectivePass = pass && runningPassport.previousPass;
        if (systematicsDirection || doNominalInternal || doNominalExternal)
            incrementTotalEvents(runningPassport.weight, runningPassport.streamlet, effectivePass);

        // Stops the streamlet if fails the flowitem and "continueWhenFail" option is set to false
        if (!pass && !continueWhenFail)
            return false;

        // If has next item chained, run it
        bool success = false;
        if (hasNext)
        {
            runningPassport.previousPass = effectivePass;
            success = next->process(runningPassport);
        }
        else
        {
            success = pass;
        }
        return success;
    }

    bool process(Passport passport = Passport())
    {
        // Automatically add substreamlets if has one
        if (substreamlet.size() > 0)
            passport.addSubStreamlet(substreamlet);
        passport.flowHistory.push_back(this);

        // Figure out what needs to be run (nominal, up and/or down)
        bool doDown = false;
        bool doUp = false;
        bool doNominalThisTime = true;
        bool partnerUpstream = false;
        if (((doData && wb->isData) || (doMC && wb->isMC) || (doEmbedding && wb->isEmbedding)) && wb->systematics->size())
        {
            if (wb->doSystematic(name))
            {
                if (passport.checkSystematics(systematicsPartners))
                {
                    for (unsigned int f = 0; f < systematicsPartners.size(); f++)
                    {
                        if (hasFlow(passport.flowHistory, systematicsPartners.at(f)))
                        {
                            partnerUpstream = true;
                            break;
                        }
                    }
                    doNominalThisTime = false;
                    if (passport.systematicsDirection==-1)
                        doDown = true;
                    else if (passport.systematicsDirection==1)
                        doUp = true;
                    else
                    {
                        if (!partnerUpstream)
                        {
                            doUp = true;
                            doDown = true;
                        }
                        doNominalThisTime = true;
                    }
                }
            }
        }
        if (doUp + doDown + doNominalThisTime > 1)
        {
            running = true;
            runningPassport = passport;
            backupOriginalValues();
            running = false;
        }

        bool success = false;
        bool nominalSuccess = false;

        if (doDown)
        {
            Passport tempPassport = passport;
            if (wb->doDump && !partnerUpstream)
            {
                if (treeDownDelayed)
                    tempPassport.delayedSysTrees.push_back(treeDownDelayed);
            }
            if (processOne(tempPassport, -1))
            {
                success = true;
                if (wb->doDump && treeDown)
                {
                    if (wb->doRemoveUnselected)
                        ao->removeUnselected();
                    treeDown->Fill();
                }
            }
        }

        if (doUp)
        {
            running = true;
            if (doDown)
                restoreOriginalValues();
            running = false;
            Passport tempPassport = passport;
            if (wb->doDump && !partnerUpstream)
            {
                if (treeUpDelayed)
                    tempPassport.delayedSysTrees.push_back(treeUpDelayed);
            }
            if (processOne(tempPassport, 1))
            {
                success = true;
                if (wb->doDump && treeUp)
                {
                    if (wb->doRemoveUnselected)
                        ao->removeUnselected();
                    treeUp->Fill();
                }
            }
        }

        if (doNominalThisTime)
        {
            running = true;
            if (doDown || doUp)
                restoreOriginalValues();
            running = false;
            Passport tempPassport = passport;
            if (wb->doDump && !partnerUpstream)
            {
                if (treeNomDelayed)
                    tempPassport.delayedNomTrees.push_back(treeNomDelayed);
            }
            if (processOne(tempPassport, 0))
            {
                success = true;
                nominalSuccess = true;
                if (wb->doDump && treeNom)
                {
                    if (wb->doRemoveUnselected)
                        ao->removeUnselected();
                    treeNom->Fill();
                }
            }
        }

        if (partnerUpstream)
            return success;
        else
            return nominalSuccess;
    }

    virtual bool executeSystematics(Passport &, short &)
    {
        // Function to be overloaded
        return true;
    }
}; // End of SystematicsBase class


} // End of CutFlow namespace
} // End of AnalysisFramework namespace
