"""
File    : Core/CodeGenerator.py
Authors : KG <Kong.Guan.Tan@cern.ch>

Functions to generate .C files from the branches.py files and the base cutflow classes in Core/CutFlow_Base.C
TODO: List out files that are generated
"""

import os, sys
from Core import compileC, compiledir
filedir = os.path.dirname(os.path.abspath(__file__))

class PhysicsObjectBase:
    dumplist = None
    nodumplist = None
    modifiable = []
    unstable = []
    n = ''
    n_elems = []
    elems = []
    typePythonToD3PD = {
        int: 'Int_t',
        float: 'Double_t',
        bool: 'Bool_t',
    }
    typeD3PDToBranch = {
        'Char_t': '/B',
        'UChar_t': '/b',
        'Short_t': '/S',
        'UShort_t': '/s',
        'Int_t': '/I',
        'UInt_t': '/i',
        'Float_t': '/F',
        'Double_t': '/D',
        'Long_t': '/L',
        'ULong_t': '/l',
        'Bool_t': '/O',
    }

    @classmethod
    def elemsList(cls):
        bnames = []
        if cls.n:
            bnames.append('n')
        for el in cls.elems:
            bnames.append(el[0])
        for n_el in cls.n_elems:
            bnames.append(n_el[0])
        return bnames

    @classmethod
    def isNElems(cls, name):
        for n_el in cls.n_elems:
            if n_el[0] == name:
                return True
        return False

    @classmethod
    def findBranch(cls, name):
        if name == 'n' and cls.n:
            return cls.makeFullTuple( ('n', cls.n, 'Int_t') )
        for el in cls.elems:
            if el[0] == name:
                return cls.makeFullTuple(el)
        for n_el in cls.n_elems:
            if n_el[0] == name:
                if n_el[2].startswith('vector'):
                    t = 'vector<'+n_el[2]+' >'
                else:
                    t = 'vector<'+n_el[2]+'>'
                n_el = tuple(n_el[:2]) + (t,) + tuple(n_el[3:])
                return cls.makeFullTuple(n_el)

    @classmethod
    def makeFullTuple(cls, el):
        # Turn a (name, address, type/default) into a (name, address, type, default, branchtype)
        if len(el) == 3:
            t_or_d = el[2]
            t = cls.getD3PDType(t_or_d)
            d = cls.getBranchDefault(t_or_d)
        elif len(el) == 4:
            t = el[2]
            d = cls.getBranchDefault(el[3])
        elif len(el) == 5:
            return el
        else:
            print "ERROR: Invalid number of items in the tuple"
            exit(0)
        b = cls.getBranchType(t)
        return el[:2] + (t, d, b)

    @classmethod
    def getD3PDType(cls, value):
        t = type(value)
        if not t == str:
            value = cls.typePythonToD3PD[t]
        return value

    @classmethod
    def getBranchType(cls, value):
        t = type(value)
        if not t == str:
            value = cls.typePythonToD3PD[t]
        if not value in cls.typeD3PDToBranch:
            return ''
        return cls.typeD3PDToBranch[value]

    @classmethod
    def getBranchDefault(cls, value):
        t = type(value)
        if t in [int, float]:
            return value
        elif t is bool:
            return value.__str__().lower()
        elif value in ['Int_t', 'UInt_t', 'int', 'unsigned int']:
            return 0
        elif value in ['Float_t', 'Double_t', 'float', 'double']:
            return 0.
        elif value in ['Bool_t', 'bool']:
            return 'false'
        elif value.startswith('vector<vector<'):
            return 'CLEAR'
        elif value.startswith('vector<'):
            return cls.getBranchDefault( value[len('vector<'):-1].strip() )
        return None

AllObjects = None
obs = None

def setupBranches(modules, name=None, location=None, outputdump=[]):
    global AllObjects
    global obs

    if not outputdump:
        outputdump = []
    obs = {}

    # Populate physics objects, merging them if more than one provided
    for m in modules:
        for o in m.__dict__.values():
            if not hasattr(o, '__bases__'):
                continue
            if not PhysicsObjectBase in o.__bases__:
                continue

            if o.__name__ in obs:
                print "\tDuplicate objects found:", obs[o.__name__], "with", o, "... will merge contents"
                for attrname in o.__dict__:
                    if attrname.startswith('__') and attrname.endswith('__'):
                        continue
                    
                    storedattr = None
                    if hasattr(obs[o.__name__], attrname):
                        storedattr = getattr(obs[o.__name__], attrname)
                    currentattr = getattr(o, attrname)

                    if storedattr == currentattr:
                        continue
                    print "\t\tWARNING:", "'"+attrname+"'", "is inconsistent ... preferring the latter"

                    if storedattr is None:
                        print "\t\t\tElement", storedattr, "was replaced with", currentattr
                        setattr(obs[o.__name__], attrname, currentattr)
                    elif not storedattr and type(currentattr) is list:
                        print "\t\t\tElement", storedattr, "was replaced with a list of", len(currentattr), "items"
                        setattr(obs[o.__name__], attrname, currentattr)
                    elif type(storedattr) is list and currentattr is True:
                        print "\t\t\tElement with", len(storedattr), "items was replaced with", currentattr
                        setattr(obs[o.__name__], attrname, currentattr)
                    elif not type(storedattr) is list and not type(currentattr) is list:
                        print "\t\t\tElement", storedattr, "was replaced with", currentattr
                        setattr(obs[o.__name__], attrname, currentattr)
                    elif type(storedattr) is list and type(currentattr) is list:
                        for el in currentattr:
                            if type(el) is tuple:
                                found = False
                                for i in range(len(storedattr)):
                                    if el[0] == storedattr[i][0]:
                                        if not storedattr[i] == el:
                                            print "\t\t\tElement", storedattr[i], "was replaced with", el
                                            storedattr[i] = el
                                        found = True
                                        break
                                if not found:
                                    print "\t\t\tElement", el, "is added"
                                    storedattr.append(el)
                            else:
                                if not el in storedattr:
                                    print "\t\t\tElement '"+el+"'", "is added"
                                    storedattr.append(el)
                    else:
                        print "\t\tERROR: don't know how to handle inconsistencies in", "'"+attr+"'", "... report this problem to the devs!"
            else:
                obs[o.__name__] = o

    # Search for duplicates in D3PD addresses - this will cause memory issues!
    d3pdlist = {}
    d3pdduplicates = []
    for oname in obs:
        o = obs[oname]
        for attrname in o.__dict__:
            if attrname.startswith('__') and attrname.endswith('__'):
                continue
            currentattr = getattr(o, attrname)
            if type(currentattr) is list:
                for el in currentattr:
                    if not type(el) is tuple:
                        continue
                    if el[1] in d3pdlist:
                        d3pdduplicates.append(el[1])
                        d3pdlist[el[1]].append(oname + '.' + attrname + '.' + el[0])
                    else:
                        d3pdlist[el[1]] = [oname + '.' + attrname + '.' + el[0]]
    if d3pdduplicates:
        print "\tERROR: Duplicates found in D3PD branch addresses!"
        print "\tPlease resolve the following D3PD branches:"
        for d in d3pdduplicates:
            print "\t\t", d, "used in", ", ".join(d3pdlist[d])
        print "\tStopping program as this will cause memory issues with TChain.GetEntry!"
        sys.exit(1)

    # Make a list of branches and their types
    inputtree_branches = {}
    inputtree_branches_num = {}
    tree = None
    if name and location:
        from ROOT import TFile, AddressOf
        for filepath in location:
            treefile = TFile.Open(filepath)
            tree = treefile.Get(name)
            for leaf in tree.GetListOfLeaves():
                inputtree_branches[leaf.GetName()] = leaf.GetTypeName()
                inputtree_branches_num.setdefault(leaf.GetName(), 0)
                inputtree_branches_num[leaf.GetName()] += 1
            treefile.Close()
    numfiles = len(location)

    exist_in_inputtree = []
    if obs:
        error = False
        for name in sorted(obs.keys()):
            # Expand branches with wildcards
            elems_new = []
            n_elems_new = []
            for bname in obs[name].elemsList():
                bname, baddress, btype, bdefault, bbranch = obs[name].findBranch(bname)
                if baddress == obs[name].n:
                    continue
                if not baddress.endswith('*'):
                    if obs[name].isNElems(bname):
                        n_elems_new.append( (bname, baddress, btype[len('vector<'):-1].strip(), bdefault, bbranch) )
                    else:
                        elems_new.append( (bname, baddress, btype, bdefault, bbranch) )
                    continue

                if not tree:
                    print "@@@@@ WARNING: Need to load ntuple to declare wildcard branches! Skipping branch:", bd3pd
                    return False

                for l in inputtree_branches:
                    if l.startswith(baddress[:-1]):
                        bname2 = bname + l[len(baddress)-1:].replace('::', '__')
                        baddress2 = l
                        btype2 = inputtree_branches[l]
                        if l == obs[name].n:
                            continue
                        if obs[name].n:
                            n_elems_new.append( obs[name].makeFullTuple( (bname2, baddress2, btype2[len('vector<'):-1].strip()) ) )
                        else:
                            elems_new.append( obs[name].makeFullTuple( (bname2, baddress2, btype2) ) )
            obs[name].n_elems = n_elems_new
            obs[name].elems = elems_new

            # Find unstable branches
            for bname in obs[name].elemsList():
                bname, baddress, btype, bdefault, bbranch = obs[name].findBranch(bname)
                if baddress in inputtree_branches:
                    exist_in_inputtree.append(baddress)
                    if not inputtree_branches_num[baddress] == numfiles:
                        obs[name].unstable.append(bname)

            # Make dumplist
            if obs[name].dumplist is None and obs[name].nodumplist is None:
                obs[name].dumplist = True
            elif obs[name].dumplist is None:
                obs[name].dumplist = []
                for bname in obs[name].elemsList():
                    bname, baddress, btype, bdefault, bbranch = obs[name].findBranch(bname)
                    if not bname in obs[name].nodumplist:
                        obs[name].dumplist.append(bname)  
            elif obs[name].nodumplist is None:
                pass # already right configuration - do nothing
            else:
                print "@@@@@ ERROR: Only nodumplist OR dumplist may be set, not both at once. Object name:", name
                error = True
            if type(obs[name].dumplist) is list and obs[name].n and not 'n' in obs[name].dumplist:
                obs[name].dumplist.insert(0, 'n')
        if error:
            return False

        # Generate the code
        kw = {}
        kw['PhysicsObjects'] = ''
        kw['AllObjectsContent'] = ''
        kw['AllObjectsPrivateContent'] = ''
        kw['ResetDefaultValues'] = ''
        kw['ForceResetDefaultValues'] = ''
        kw['InitialiseVectors'] = ''
        kw['ResetSelected'] = ''
        kw['RemoveSelected'] = ''
        for name in sorted(obs.keys()):
            kw['StructName'] = name
            kw['StructContent'] = ''
            kw['InstanceName'] = name.lower()
            kw['CopyOrdered'] = ''
            kw['RemoveContentOrdered'] = ''
            kw['RemoveContent'] = ''
            kw['RemoveResize'] = ''
            missing_in_inputtree = []
            inconsistent_type = []

            if obs[name].n:
                kw['StructContent'] += oneContent('vector<bool>', 'Selected')
                kw['StructContent'] += oneContent('vector<int>', 'Order')
                kw['InitialiseVectors'] += oneInitVector(name.lower(), 'Selected', 'vector<bool>')
                kw['InitialiseVectors'] += oneInitVector(name.lower(), 'Order', 'vector<int>')
                kw['ResetSelected'] += oneSelected(name.lower())
            for bname in obs[name].elemsList():
                bname, baddress, btype, bdefault, bbranch = obs[name].findBranch(bname)
                kw['StructContent'] += oneContent(btype, bname)
                kw['InitialiseVectors'] += oneInitVector(name.lower(), bname, btype)
                if not baddress in exist_in_inputtree or bname in obs[name].unstable:
                    if btype.startswith('vector<'):
                        kw['ResetDefaultValues'] += oneDefault(name.lower(), bname, 'CLEAR')
                        if obs[name].n and not bdefault == 'CLEAR':
                            kw['ResetSelected'] += oneDefault(name.lower(), bname, bdefault, btype)
                    elif baddress in exist_in_inputtree and \
                            bname in obs[name].unstable and \
                            not (obs[name].modifiable is True or bname in obs[name].modifiable):
                        kw['ForceResetDefaultValues'] += oneDefault(name.lower(), bname, bdefault)
                    else:
                        kw['ResetDefaultValues'] += oneDefault(name.lower(), bname, bdefault)
                if not baddress in exist_in_inputtree:
                    missing_in_inputtree.append( (bname, baddress) )
                else:
                    btype2 = inputtree_branches[baddress]
                    if not btype2 == btype:
                        inconsistent_type.append( (bname, btype, btype2) )

                if obs[name].isNElems(bname) and (obs[name].dumplist is True or bname in obs[name].dumplist):
                    kw['AllObjectsPrivateContent'] += oneContent(btype, '__copy_'+name.lower()+'_'+bname)
                    kw['InitialiseVectors'] += oneInitVector2('__copy_'+name.lower()+'_'+bname, btype)
                    kw['CopyOrdered'] += oneCopy(name.lower()+'.'+bname, '__copy_'+name.lower()+'_'+bname, btype)
                    if not baddress in exist_in_inputtree or bname in obs[name].unstable:
                        kw['RemoveContentOrdered'] += oneRemoveContentOrderedCheck(name.lower(), bname)
                        kw['RemoveContent'] += oneRemoveContentCheck(name.lower(), bname)
                        kw['RemoveResize'] += oneRemoveResizeCheck(name.lower(), bname)
                    else:
                        kw['RemoveContentOrdered'] += oneRemoveContentOrdered(name.lower(), bname)
                        kw['RemoveContent'] += oneRemoveContent(name.lower(), bname)
                        kw['RemoveResize'] += oneRemoveResize(name.lower(), bname)

            kw['PhysicsObjects'] += template_struct % kw
            kw['AllObjectsContent'] += oneContent(name, name.lower())
            if obs[name].n and (outputdump is True or name in outputdump):
                kw['RemoveSelected'] += template_remove % kw
            if missing_in_inputtree:
                print "@@@@@ WARNING: The following branches are missing from input files in", name, ":"
                for btuple in missing_in_inputtree:
                    print "     %30s : %s" % btuple
            if inconsistent_type:
                print "@@@@@ WARNING: The following branches have inconsistent types in", name, ":"
                for btuple in inconsistent_type:
                    print "     %30s : Declared as: %-30s In input file: %-30s" % btuple

        full_str = template_full % kw
        filepath = filedir+'/../libs/Branches_generated.C'
        if not os.path.isfile(filepath) or not open(filepath, 'r').read() == full_str:
            open(filepath, 'w').write(full_str)
        compileC(filepath)
        from ROOT import AnalysisFramework
        AllObjects = AnalysisFramework.Branches.AllObjects
        return True
    return False

def setupForkItem(doRemoveUnselected = False, outputdump=[]):
    if obs:
        if not outputdump:
            outputdump = []
        kw = {}
        kw['DeclareFork'] = ''
        kw['InitVectors'] = ''
        kw['BackupContent'] = ''
        kw['RestoreContent'] = ''
        kw['SaveContent'] = ''
        kw['LoadContent'] = ''
        for name in sorted(obs.keys()):
            if obs[name].n and (obs[name].modifiable is True or 'Selected' in obs[name].modifiable):
                kw['DeclareFork'] += oneContent('vector<bool>', name.lower()+'_original_selected')
                kw['DeclareFork'] += oneContent('vector<bool>', name.lower()+'_saved_selected')
                kw['DeclareFork'] += oneContent('vector<int>', name.lower()+'_original_order')
                kw['DeclareFork'] += oneContent('vector<int>', name.lower()+'_saved_order')
                kw['InitVectors'] += oneInitVector2(name.lower()+'_original_selected', 'vector<bool>')
                kw['InitVectors'] += oneInitVector2(name.lower()+'_saved_selected', 'vector<bool>')
                kw['InitVectors'] += oneInitVector2(name.lower()+'_original_order', 'vector<int>')
                kw['InitVectors'] += oneInitVector2(name.lower()+'_saved_order', 'vector<int>')
                kw['BackupContent'] += oneCopy('ao->'+name.lower()+'.Selected', name.lower()+'_original_selected', 'vector<bool>')
                kw['RestoreContent'] += oneCopy(name.lower()+'_original_selected', 'ao->'+name.lower()+'.Selected', 'vector<bool>')
                kw['BackupContent'] += oneCopy('ao->'+name.lower()+'.Order', name.lower()+'_original_order', 'vector<int>')
                kw['RestoreContent'] += oneCopy(name.lower()+'_original_order', 'ao->'+name.lower()+'.Order', 'vector<int>')
                kw['SaveContent'] += oneCopy('ao->'+name.lower()+'.Selected', name.lower()+'_saved_selected', 'vector<bool>')
                kw['LoadContent'] += oneCopy(name.lower()+'_saved_selected', 'ao->'+name.lower()+'.Selected', 'vector<bool>')
                kw['SaveContent'] += oneCopy('ao->'+name.lower()+'.Order', name.lower()+'_saved_order', 'vector<int>')
                kw['LoadContent'] += oneCopy(name.lower()+'_saved_order', 'ao->'+name.lower()+'.Order', 'vector<int>')
            for bname in obs[name].elemsList():
                bname, baddress, btype, bdefault, bbranch = obs[name].findBranch(bname)
                if  obs[name].modifiable is True or \
                    bname in obs[name].modifiable or \
                    ( \
                        (obs[name].isNElems(bname) or baddress == obs[name].n) and \
                        doRemoveUnselected and \
                        (outputdump is True or name in outputdump) and \
                        (obs[name].dumplist is True or bname in obs[name].dumplist)
                    ):
                    kw['DeclareFork'] += oneContent(btype, name.lower()+'_original_'+bname)
                    kw['DeclareFork'] += oneContent(btype, name.lower()+'_saved_'+bname)
                    kw['InitVectors'] += oneInitVector2(name.lower()+'_original_'+bname, btype)
                    kw['InitVectors'] += oneInitVector2(name.lower()+'_saved_'+bname, btype)
                    kw['BackupContent'] += oneCopy('ao->'+name.lower()+'.'+bname, name.lower()+'_original_'+bname, btype)
                    kw['RestoreContent'] += oneCopy(name.lower()+'_original_'+bname, 'ao->'+name.lower()+'.'+bname, btype)
                    kw['SaveContent'] += oneCopy('ao->'+name.lower()+'.'+bname, name.lower()+'_saved_'+bname, btype)
                    kw['LoadContent'] += oneCopy(name.lower()+'_saved_'+bname, 'ao->'+name.lower()+'.'+bname, btype)
        full_str = template_fork % kw
        filepath = filedir+'/../libs/ForkItem_generated.C'
        if not os.path.isfile(filepath) or not open(filepath, 'r').read() == full_str:
            open(filepath, 'w').write(full_str)
        compileC(filepath)
        return True
    return False

###### Below this point is internal ######

template_full = """/*
 *
 * Do not modify - this file is automatically generated.
 *
 */

#pragma once

#include "TROOT.h"
#include "TTree.h"

#include<vector>

using namespace std;
namespace AnalysisFramework
{
namespace Branches
{
%(PhysicsObjects)s

class AllObjects
{
private:%(AllObjectsPrivateContent)s

public:%(AllObjectsContent)s

    AllObjects()
    {%(InitialiseVectors)s
    }

    void resetToDefaultValues()
    {%(ResetDefaultValues)s
    }

    void forceResetToDefaultValues()
    {%(ForceResetDefaultValues)s
    }

    void resetSelected()
    {%(ResetSelected)s
    }

    void removeUnselected()
    {%(RemoveSelected)s
    }
};

} // End namespace Branches
} // End namespace AnalysisFramework
"""


template_struct = """
struct %(StructName)s
{%(StructContent)s
};
"""

template_content = """
    %(Type)s %(Branch)s;"""

template_initvectors = """
        %(InstanceName)s.%(Branch)s = new %(Type)s();"""

template_initvectors2 = """
        %(Branch)s = new %(Type)s();"""

template_default = """
        %(InstanceName)s.%(Branch)s = %(Value)s;"""

template_defaultclear = """
        %(InstanceName)s.%(Branch)s->clear();"""

template_defaultvector = """
        if ( (Int_t)(%(InstanceName)s.%(Branch)s->size())!=%(InstanceName)s.n )
        {
            %(InstanceName)s.%(Branch)s->resize(%(InstanceName)s.n);
            for (int i = 0; i < %(InstanceName)s.n; i++)
                %(InstanceName)s.%(Branch)s->at(i) = %(Value)s;
        }"""

template_selected = """
        %(InstanceName)s.Order->clear();
        %(InstanceName)s.Selected->resize(%(InstanceName)s.n);
        for (int i = 0; i < %(InstanceName)s.n; i++)
            %(InstanceName)s.Selected->at(i) = true;"""

template_remove = """
        int %(InstanceName)sNewSize = 0;
        bool %(InstanceName)sOrder = %(InstanceName)s.Order->size();
        if (%(InstanceName)sOrder)
        {%(CopyOrdered)s
        }
        for (int i=0, size=%(InstanceName)s.Selected->size(); i < size; i++, %(InstanceName)sNewSize++)
        {
            int r = i;
            if (%(InstanceName)sOrder) r = %(InstanceName)s.Order->at(i);
            while (!%(InstanceName)s.Selected->at(r))
            {
                if (++i >= size) break;
                r = i;
                if (%(InstanceName)sOrder) r = %(InstanceName)s.Order->at(i);
            }
            if (i == %(InstanceName)sNewSize && !%(InstanceName)sOrder) continue;
            if (i >= size) break;
            if (%(InstanceName)sOrder)
            {%(RemoveContentOrdered)s
            }
            else
            {%(RemoveContent)s
            }
        }
        %(RemoveResize)s
        %(InstanceName)s.n = %(InstanceName)sNewSize;
"""

template_removecontent = """
                %(InstanceName)s.%(Branch)s->at(%(InstanceName)sNewSize) = %(InstanceName)s.%(Branch)s->at(r);"""
template_removecontent_check = """
                if ((Int_t)(%(InstanceName)s.%(Branch)s->size())==%(InstanceName)s.n) %(InstanceName)s.%(Branch)s->at(%(InstanceName)sNewSize) = %(InstanceName)s.%(Branch)s->at(r);"""

template_removecontent_ordered = """
                %(InstanceName)s.%(Branch)s->at(%(InstanceName)sNewSize) = __copy_%(InstanceName)s_%(Branch)s->at(r);"""
template_removecontent_ordered_check = """
                if ((Int_t)(%(InstanceName)s.%(Branch)s->size())==%(InstanceName)s.n) %(InstanceName)s.%(Branch)s->at(%(InstanceName)sNewSize) = __copy_%(InstanceName)s_%(Branch)s->at(r);"""

template_removeresize = """
        %(InstanceName)s.%(Branch)s->resize(%(InstanceName)sNewSize);"""
template_removeresize_check = """
        if ((Int_t)(%(InstanceName)s.%(Branch)s->size())==%(InstanceName)s.n) %(InstanceName)s.%(Branch)s->resize(%(InstanceName)sNewSize);"""

template_copyscalar = """
        %(To)s = %(From)s;"""

template_copyvector = """
        *(%(To)s) =  *(%(From)s);"""
        #copyVector(%(From)s, %(To)s);"""

template_copyvectorvector = """
        *(%(To)s) =  *(%(From)s);"""
        #copyVectorVector(%(From)s, %(To)s);"""

def oneContent(type_, branch_):
    if type_.startswith('vector'):
        branch_ = '*' + branch_
    kw = {'Type': type_, 'Branch': branch_}
    return template_content % (kw)

def oneInitVector(instance_, branch_, type_):
    if not type_.startswith('vector'):
        return ''
    kw = {'InstanceName': instance_, 'Branch': branch_, 'Type': type_}
    return template_initvectors % (kw)

def oneInitVector2(branch_, type_):
    if not type_.startswith('vector'):
        return ''
    kw = {'Branch': branch_, 'Type': type_}
    return template_initvectors2 % (kw)

def oneDefault(instance_, branch_, value_, type_=''):
    if value_ is None:
        return ''
    kw = {'InstanceName': instance_, 'Branch': branch_, 'Value': value_}
    if value_ == 'CLEAR':
        return template_defaultclear % (kw)
    elif type_.startswith('vector<'):
        return template_defaultvector % (kw)
    else:
        return template_default % (kw)

def oneSelected(instance_):
    kw = {'InstanceName': instance_}
    return template_selected % (kw)

def oneRemoveContentOrdered(instance_, branch_):
    kw = {'InstanceName': instance_, 'Branch': branch_}
    return template_removecontent_ordered % (kw)

def oneRemoveContentOrderedCheck(instance_, branch_):
    kw = {'InstanceName': instance_, 'Branch': branch_}
    return template_removecontent_ordered_check % (kw)

def oneRemoveContent(instance_, branch_):
    kw = {'InstanceName': instance_, 'Branch': branch_}
    return template_removecontent % (kw)

def oneRemoveContentCheck(instance_, branch_):
    kw = {'InstanceName': instance_, 'Branch': branch_}
    return template_removecontent_check % (kw)

def oneRemoveResize(instance_, branch_):
    kw = {'InstanceName': instance_, 'Branch': branch_}
    return template_removeresize % (kw)

def oneRemoveResizeCheck(instance_, branch_):
    kw = {'InstanceName': instance_, 'Branch': branch_}
    return template_removeresize_check % (kw)

def oneCopy(from_, to_, type_):
    if type_.startswith('vector<vector'):
        kw = {'From': from_, 'To': to_}
        return template_copyvectorvector % (kw)
    elif type_.startswith('vector'):
        kw = {'From': from_, 'To': to_}
        return template_copyvector % (kw)
    else:
        kw = {'From': from_, 'To': to_}
        return template_copyscalar % (kw)

template_fork = """/*
 *
 * Do not modify - this file is automatically generated.
 *
 */

#pragma once

#include "../Core/CutFlow_Base.C"
#include <iostream>

namespace AnalysisFramework
{
namespace CutFlows
{

class ForkItem : public ForkBase
{
public:
    ForkItem(std::string name_, std::string substreamlet_ = "") : ForkBase(name_, substreamlet_)
    {%(InitVectors)s
    }
    %(DeclareFork)s

    void backupOriginalValues()
    {%(BackupContent)s
    }

    void restoreOriginalValues()
    {%(RestoreContent)s
    }

    void saveValues()
    {%(SaveContent)s
    }

    void loadValues()
    {%(LoadContent)s
    }
}; // End of ForkItem class


class ForkLoopItem : public ForkLoopBase
{
public:
    ForkLoopItem(std::string name_, std::string substreamlet_ = "") : ForkLoopBase(name_, substreamlet_)
    {%(InitVectors)s
    }
    %(DeclareFork)s

    void backupOriginalValues()
    {%(BackupContent)s
    }

    void restoreOriginalValues()
    {%(RestoreContent)s
    }

    void saveValues()
    {%(SaveContent)s
    }

    void loadValues()
    {%(LoadContent)s
    }
}; // End of ForkLoopItem class


class SystematicsItem : public SystematicsBase
{
public:
    SystematicsItem(std::string name_, std::string substreamlet_ = "") : SystematicsBase(name_, substreamlet_)
    {%(InitVectors)s
    }
    %(DeclareFork)s

    void backupOriginalValues()
    {%(BackupContent)s
    }

    void restoreOriginalValues()
    {%(RestoreContent)s
    }

    void saveValues()
    {%(SaveContent)s
    }

    void loadValues()
    {%(LoadContent)s
    }
}; // End of SystematicsItem class


} // End namespace Branches
} // End namespace AnalysisFramework
"""
