#!/usr/bin/python

import pygtk
import gtk

class Selector:
    def __init__(self,tag,particle_ids,ranges):
        self.tag     = tag
        self.partids = particle_ids
        self.ranges  = ranges
        self.values  = {}
        self.on      = True

    def write(self,runfile):
        if not(self.on):
            return
        out = "   "+self.tag+" "
        if self.tag=="NJetFinder":
            expo = -1
            if self.values["Mode"]=="KT":
                expo = 1
            if self.values["Mode"]=="Cam-A":
                expo = 0
            if expo==-1 or expo==1 or expo==0:
                out = out+str(self.values["min_n"])+" "
                out = out+str(self.values["PT"])+" " +str(self.values["PT"])+" "
                out = out+str(self.values["R"])+" "+str(expo)+"\n"
            else:
                out = "%%% unknown mode, not implemented yet, will ignore it."
        if self.tag=="IsolationCut":
            out = out + str(" 22 ")
            out = out+str(self.values["R"])+" "
            out = out+str(self.values["expo"])+" "
            out = out+str(self.values["epsilon"])+"\n"
        if (self.tag=="Mass2" or self.tag=="DeltaR" or 
            self.tag=="Rapidity" or self.tag=="PT"):
            for part in self.partids:
                out = out+part+" "
            out = out+str(self.values["min"])+" "
            out = out+str(self.values["max"])+"\n"
        runfile.write(out)

    def getTag(self):
        return self.tag

    def getRanges(self):
        return self.ranges

    def getRange(self,keyword):
        return self.ranges[keyword]

    def getParticleIds(self):
        return self.partids

    def setValue(self,keyword,value):
        self.values[keyword] = value
        print "Set: ",keyword," = ",value," for ",self.tag,"."

    def getValue(self,keyword):
        return self.values[keyword]

    def switch(self):
        self.on = not(self.on)

    def isOn(self):
        return self.on



class SelectorBox:
    def __init__(self):
        self.initialiseDefaults(None,None)
        self.sellist = []
        pass

    def initialiseDefaults(self,collider,process,minjets=0):
        self.E_CMS     = 0.
        self.sellist   = []
        self.process   = process
        self.minjets   = minjets
        self.Leptons   = ["11","-11","13","-13","15","-15"]
        self.Neutrinos = ["12","-12","14","-14","16","-16"]
        print "In SelectorBox: minjets = ",self.minjets
        if collider!=None:
            self.E_CMS   = collider[3]+collider[5]
        if (self.process!=None):
            FS = self.constructFinalState()
            print "*** Final State: "
            for part in FS:
                print part," "
            self.fillSelectorList(FS)
            
    def write(self,runfile):
        runfile.write("%%% Selectors\n\n")
        runfile.write("(selector){\n")
        for sel in self.sellist:
            print "  write-out for ",sel.getTag()
            sel.write(runfile)
        runfile.write("}(selector)\n\n")

    def getSelectorList(self):
        return self.sellist

    def inSelectorList(self,tag,parts):
        for sel in self.sellist:
            if sel.getTag()==tag:
                if len(parts)==0:
                    return True
                for part in parts:
                    flag = True
                    if not(part in sel.getParticleIds()):
                        flag = False
                    if flag:
                        return True
        return False

    def fillSelectorList(self,FS):
        for part in FS:
            if (part!="93" and not(part in self.Neutrinos) and
                not(self.isDecayer(part))):
                if not(self.inSelectorList("PT",[part])):
                    pt = Selector("PT",[part],
                                  {"min":[1.,self.E_CMS,20.,5.],
                                   "max":[1.,self.E_CMS,self.E_CMS,5.]})
                    self.sellist.append(pt)
                if not(self.inSelectorList("Rapidity",[part])):
                    rap = Selector("Rapidity",[part],
                                   {"min":[-10.,10.,-5.,0.5],
                                    "max":[-10.,10., 5.,0.5]})
                    self.sellist.append(rap)
            if part=="93":
                if not(self.inSelectorList("NJetFinder",[])):
                    print part," --> init JF for ",FS.count("93")," jets"
                    njf = Selector("NJetFinder",[],
                                   {"PT":[0.,self.E_CMS,20.,10.],
                                    "R":[0.,2.0,0.6,0.1],
                                    "Mode":[]})
                    njf.setValue("min_n",FS.count("93"))
                    self.sellist.append(njf)
            if part=="22":
                if not(self.inSelectorList("IsolationCut",[22])):
                    print part," --> init FIC for ",FS.count("22")," photons"
                    fic = Selector("IsolationCut",[],
                                   {"R":[0.05,1.0,0.4,0.1],
                                    "expo":[0.0,2.0,1.0,0.1],
                                    "epsilon":[0.0,2.0,1.0,0.1]})
                    fic.setValue("min_n",FS.count("93"))
                    self.sellist.append(fic)
            if part in self.Leptons or part in self.Neutrinos:
                kfc  = -int(part)
                anti = str(kfc)
                if (FS.count(anti)>0 and
                    not(self.inSelectorList("NJetFinder",[part,anti]))):
                    m2 = Selector("Mass2",[part,anti],
                                  {"min":[1.,self.E_CMS,81.2,10.],
                                   "max":[1.,self.E_CMS,101.2,10.]})
                    self.sellist.append(m2)
            if part in self.Leptons:
                if (FS.count("93")>0 and
                    not(self.inSelectorList("DeltaR",[part,"93"]))):
                    dr = Selector("DeltaR",[part,"93"],
                                  {"min":[0.,2.0,0.4,0.1],
                                   "max":[0.,20.0,20.0,0.1]})
                    self.sellist.append(dr)
                kfc  = -int(part)
                if kfc<0: 
                    kfc = kfc-1
                else:
                    kfc = kfc+1
                anti = str(kfc)
                if (FS.count(anti)>0 and
                    not(self.inSelectorList("NJetFinder",[part,anti]))):
                    m2 = Selector("Mass2",[part,anti],
                                  {"min":[1.,self.E_CMS,70.4,10.],
                                   "max":[1.,self.E_CMS,90.4,10.]})
                    self.sellist.append(m2)
        print "Ended with ",len(self.sellist)," selectors."

    def constructFinalState(self):
        FS = []
        for line in self.process.declarationlines:
            print "* ",line
            for char in line:
                if char in ";":
                    line = line.replace(char,'')
            content = line.split(" ")
            take    = False 
            for part in content:
                if self.minjets>0 and part=="93{NJETS}":
                        part = "93"
                if take:
                    FS.append(part)
                else:
                    if part=="->":
                        take = True
        return FS

    def isDecayer(self,part):
        for char in part:
            if char=="[" or char=="{":
                return True
        return False
