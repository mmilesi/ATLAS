"""
File    : Core/DatasetManager.py
Authors : KG <Kong.Guan.Tan@cern.ch>

Functions to parse and check sample names in both the Files/*D3PD_datasets.txt and samples*.csv files
TODO: Show examples of usage?
"""

import os

tags = ['@', '#', '~', '=', '-']
comment = '//'
filedir = os.path.dirname(os.path.abspath(__file__))

class DatasetManager:
    def __init__(self, filelist = []):
        self.list_ = []
        for filename in filelist:
            current_tag = {}
            for t in tags:
                current_tag[t] = 'DEFAULT'
            for line in open(filename):
                line = line.strip()
                if line.startswith(comment):
                    continue
                for t in tags:
                    if line.startswith(t):
                        current_tag[t] = line.strip(t+' ')
                if line.startswith('mc') or line.startswith('data') or line.startswith('group') or line.startswith('user'):
                    line = line.strip('/')
                    self.list_.append( (line, current_tag.values()) )

    def getList(self, truncate=False, tags = [], filter_str=''):
        if type(tags) == str:
            tags = [tags]
        output = []
        for name_, tags_ in self.list_:
            triggered = False
            filtered = False
            if tags == []:
                triggered = True
            else:
                for item in tags:
                    if item in tags_:
                        triggered = True
                    else:
                        filtered = True
            if triggered and not filtered and filter_str in name_:
                if truncate:
                    name_ = self.truncate(name_)
                output.append(name_)
        return output

    def getListID(self, tags = [], filter_str=''):
        if type(tags) == str:
            tags = [tags]
        output = []
        for name_, tags_ in self.list_:
            triggered = False
            filtered = False
            if tags == []:
                triggered = True
            else:
                for item in tags:
                    if item in tags_:
                        triggered = True
                    else:
                        filtered = True
            if triggered and not filtered and filter_str in name_:
                name_ = self.truncate(name_)
                output.append(name_.split('.')[0])
        return output

    def getHiggsMass(self, pathname):
        pathname = self.truncate(pathname, remove_simrecotags=False)

        for name_, tags_ in self.list_:
            name_ = self.truncate(name_, remove_simrecotags=False)
            if pathname.startswith(name_):
                if not self.check(pathname, ['Signal']):
                    return None
                if "ggF" in tags_:
                    index = name_.find('_ggH')
                    if index < 0:
                        return None
                    return int(name_[index+4:].split('_')[0])
                elif "VBF" in tags_:
                    index = name_.find('_VBFH')
                    if index < 0:
                        return None
                    return int(name_[index+5:].split('_')[0])
                elif "WH" in tags_:
                    index = name_.find('WH')
                    if index < 0:
                        return None
                    return int(name_[index+2:].split('_')[0])
                elif "ZH" in tags_:
                    index = name_.find('ZH')
                    if index < 0:
                        return None
                    return int(name_[index+2:].split('_')[0])
                elif "MSSMggA" in tags_:
                    index = name_.find('MA')
                    if index < 0:
                        return None
                    return int(name_[index+2:].split('TB')[0])
                elif "MSSMbbA" in tags_:
                    index = name_.find('MA')
                    if index < 0:
                        return None
                    return int(name_[index+2:].split('TB')[0])
                elif "LFV" in tags_:
                    return 125
                else:
                    return None
        return None

    def getTags(self, pathname):
        pathname = self.truncate(pathname, remove_simrecotags=False)

        for name_, tags_ in self.list_:
            name_ = self.truncate(name_, remove_simrecotags=False)
            if pathname.startswith(name_):
                return tags_
        return []

    def getListSamples(self, samplesfile='Files/samples.csv', genericPath=False):
        samples = []
        h = {}
	
	filename = None
	if not genericPath:
	   filename = filedir+'/../'+samplesfile
	else:
	   filename = samplesfile
        
	for line in open(filename):
            tokens = line.strip().split(',')
            if not h:
                for i in range(len(tokens)):
                    h[tokens[i].strip('"')] = i
            else:
                entry = {}
                for header in h:
		    #print("token: {0}".format(tokens[h[header]]))
                    entry[header] = tokens[h[header]].strip('"')
                samples.append(entry)
        return samples

    def contains(self, pathname):
        pathname = self.truncate(pathname, remove_simrecotags=False)

        for name_, tags_ in self.list_:
            name_ = self.truncate(name_, remove_simrecotags=False)
            if pathname.startswith(name_):
                return True
        return False

    def check(self, pathname, tags=[]):
        if type(tags) == str:
            tags = [tags]
        pathname = self.truncate(pathname, remove_simrecotags=False)

        for name_, tags_ in self.list_:
            name_ = self.truncate(name_, remove_simrecotags=False)
            if pathname.startswith(name_):
                filtered = False
                if tags == []:
                    print "Need to provide list of tags to check!"
                    return False
                else:
                    for t in tags:
                        if not t in tags_:
                            filtered = True
                            break
                if not filtered:
                    return True
        return False

    def truncate(self, name, remove_simrecotags=False):
        def isnum(string):
            try:
                num = int(string)
                return True
            except:
                return False

        name = name.strip("/")
        name = name.split("/")[-1]
        removelist = [
            #'mc11_7TeV',
            #'mc12_8TeV',
            #'data11_7TeV',
            #'data12_8TeV',
            'group',
            'perf-tau',
            'phys-higgs',
            'LHSkim',
            'merge',
            'NTUP_TAU',
            'user',
            'tuna',
        ]
        tokens = name.split('.')
        newtokens = []
        for t in tokens:
            if t in removelist: continue
            found_e = found_s = found_r = False
            subtokens = t.split('_')
            for st in subtokens:
                if st.startswith('e') and isnum(st[1:]): found_e = True
                if st.startswith('s') and isnum(st[1:]): found_s = True
                if st.startswith('r') and isnum(st[1:]): found_r = True
            if found_e and found_s and found_r and remove_simrecotags:
                newsubtokens = [st for st in subtokens if st.startswith('p') and isnum(st[1:])]
                t = '_'.join(newsubtokens)
            newtokens.append(t)
        name = '.'.join(newtokens)
        if name.startswith('user.'):
            name = '.'.join(name.split('.')[2:])
        return name

    def printContents(self):
        for name_, tags_ in self.list_:
            for t in tags_:
                print "%-10s" % t,
            print "    ", name_

if __name__ == '__main__':
    d = DatasetManager(['../Files/MC_D3PD_datasets.txt', '../Files/Data_D3PD_datasets.txt', '../Files/Embedding_D3PD_datasets.txt'])
    d.printContents()
    pass
