#!/usr/bin/python

import pygtk
import gtk

class ProcessBox:
    def __init__(self):
        self.LCmap = {
            "Jets"                       : PI(["11 -11 -> 93 93 93{NJETS};"],
                                              ["MODEL = SM;"],
                                              6,4,0),
            }
        self.HCmap = {
            # Format for input:
            # [[<process_tag>,<decay_tags>],maxjets,nlojets,extra_FS_particles]
            # use -1,-1,-1 to steer alternatives, may have to use
            # the read-in facilities of Sherpa in a smarter way
            #
            "MinimumBias"                : PI(["MinimumBias"],[""],
                                              -1,-1,-1),
            "Jets"                       : PI(["93 93 -> 93 93 93{NJETS};"],
                                              ["MODEL = SM;"],
                                              4,2,1),
            "Gamma+Jets"                 : PI(["93 93 -> 22 93 93{NJETS};"],
                                              ["MODEL = SM;"],
                                              4,2,0),
            "(W^+ -> e^+ nu) + Jets"     : PI(["93 93 -> 24[a] 93{NJETS};",
                                               "24[a] -> -11 12;"],
                                              ["MODEL = SM;"],
                                              5,2,2),
            "(W^+ -> mu^+ nu) + Jets"    : PI(["93 93 -> 24[a] 93{NJETS};",
                                               "24[a] -> -13 14;"],
                                              ["MODEL = SM;"],
                                              5,2,2),
            "(W^+ -> tau^+ nu) + Jets"   : PI(["93 93 -> 24[a] 93{NJETS};",
                                               "24[a] -> -15 16;"],
                                              ["MODEL = SM;"],
                                              5,2,2),
            "(W^- -> e^- nu) + Jets"     : PI(["93 93 -> -24[a] 93{NJETS};",
                                               "-24[a] -> 11 -12;"],
                                              ["MODEL = SM;"],
                                              5,2,2),
            "(W^- -> mu^- nu) + Jets"    : PI(["93 93 -> -24[a] 93{NJETS};",
                                               "-24[a] -> 13 -14;"],
                                              ["MODEL = SM;"],
                                              5,2,2),
            "(W^- -> tau^- nu) + Jets"   : PI(["93 93 -> -24[a] 93{NJETS};",
                                               "-24[a] -> 15 -16;"],
                                              ["MODEL = SM;"],
                                              5,2,2),
            "(Z -> e^- e^+) + Jets"      : PI(["93 93 -> 23[a] 93{NJETS};",
                                               "23[a] -> 11 -11;"],
                                              ["MODEL = SM;"],
                                              5,2,2),
            "(Z -> mu^- mu^+) + Jets"    : PI(["93 93 -> 23[a] 93{NJETS};",
                                               "23[a] -> 13 -13;"],
                                              ["MODEL = SM;"],
                                              5,2,2),
            "(Z -> tau^- tau^+) + Jets"  : PI(["93 93 -> 23[a] 93{NJETS};",
                                               "23[a] -> 15 -15;"],
                                              ["MODEL = SM;"],
                                              5,2,2),
            "(Z -> nu nu) + Jets"        : PI(["93 93 -> 23[a] 93{NJETS};",
                                               "23[a] -> 11 -11;"],
                                              ["MODEL = SM;"],
                                              5,2,2),
            "(e^- e^+) + Jets"           : PI(["93 93 -> 11 -11 93{NJETS};"],
                                              ["MODEL = SM;"],
                                              5,2,2),
            "(mu^- mu^+) + Jets"         : PI(["93 93 -> 13 -13 93{NJETS};"],
                                              ["MODEL = SM;"],
                                              5,2,2),
            "(tau^- tau^+) + Jets"       : PI(["93 93 -> 15 -15 93{NJETS};"],
                                              ["MODEL = SM;"],
                                              5,2,2),
            }
        
        self.HC = gtk.TreeStore(str, bool)
        parent = self.HC.append(None,["Minimum Bias",False])
        self.HC.append(parent,["MinimumBias",True])

        parent = self.HC.append(None,["Hard QCD",False])
        self.HC.append(parent,["Jets",True])

        parent = self.HC.append(None,["V+jets",False])
        self.HC.append(parent,["Gamma+Jets",True])
        self.HC.append(parent,["(W^+ -> e^+ nu) + Jets",True])
        self.HC.append(parent,["(W^+ -> mu^+ nu) + Jets",True])
        self.HC.append(parent,["(W^+ -> tau^+ nu) + Jets",True])
        self.HC.append(parent,["(W^- -> e^- nu) + Jets",True])
        self.HC.append(parent,["(W^- -> mu^- nu) + Jets",True])
        self.HC.append(parent,["(W^- -> tau^- nu) + Jets",True])
        self.HC.append(parent,["(W^{+,-} -> mu^{+,-} nu) + Jets",True])
        self.HC.append(parent,["(Z -> e^- e^+) + Jets",True])
        self.HC.append(parent,["(Z -> mu^- mu^+) + Jets",True])
        self.HC.append(parent,["(Z -> tau^- tau^+) + Jets",True])
        self.HC.append(parent,["(Z -> nu nu) + Jets",True])
        self.HC.append(parent,["(e^- e^+) + Jets",True])
        self.HC.append(parent,["(mu^- mu^+) + Jets",True])
        self.HC.append(parent,["(tau^- tau^+) + Jets",True])

        parent = self.HC.append(None,["H+jets (gluon fusion)",False])
        self.HC.append(parent,["(H -> gamma gamma) + Jets",True])
        self.HC.append(parent,["(H -> tau tau) + Jets",True])
        self.HC.append(parent,["(H -> e^+ nu mu^- nu) + Jets",True])
        self.HC.append(parent,["(H -> e^+ e^- mu^+ mu^-) + Jets",True])
        self.HC.append(parent,["(H -> e^+ e^- nu nu) + Jets",True])
      
        self.default_process_tag = "Jets"
        self.default_process     = self.HCmap[self.default_process_tag]
        self.process_tag         = self.default_process_tag
        self.process             = self.HCmap[self.process_tag]
        self.processes           = self.HC
        self.totjets             = 0
        self.nlojets             = 0
        self.ckkw_param          = 10.

    def initialiseDefaults(self,collider):
        processes = None
        if ((collider[2]==2212 or collider[2]==-2212) and
            (collider[4]==2212 or collider[4]==-2212)):
            self.processes = self.HC
        if ((collider[2]==11 or collider[2]==-11) and
            (collider[4]==11 or collider[4]==-11)):
            self.processes = self.LC

        if self.processes!=None:
            for row in self.processes:
                if row[1]==self.default_process_tag:
                    self.default_process = row
                    self.process = self.default_process
        if self.process!=None:
            print "Initialised default process: ",self.default_process_tag
        
    def getNJets(self):
        print "getNJets = ",self.process.declarationlines
        if self.process!=None:
            print "yields ",self.process.totjets,self.process.nlojets
            return self.process.totjets,self.process.nlojets
        return None,None

    def getOptions(self):
        return self.processes

    def setProcess(self,process_tag):
        print "Set process = ",process_tag
        self.process_tag = process_tag
        self.process     = self.HCmap[self.process_tag]

    def isNLO(self):
        return (self.nlojets>0)

    def setJetMultis(self,minjets=0,totjets=0,nlojets=0):
        self.minjets = minjets
        self.totjets = totjets
        self.nlojets = nlojets
        if (self.totjets<self.nlojets):
            self.totjets=self.nlojets

    def setCKKW(self,ckkw_param):
        self.ckkw_param = ckkw_param

    def write(self,runfile):
        print "Write for ",self.process_tag
        runfile.write('\n')
        runfile.write('%%% Process setup: '+self.process_tag+'\n\n')
        runfile.write('(processes){\n')
        lineno = 0
        for line in self.process.declarationlines:
            if lineno==0:
                runfile.write('  Process '+line)
            else: 
                runfile.write('    Decay '+line)
            runfile.write('\n')
            lineno += 1
        if self.totjets>0:
            runfile.write('    CKKW sqr(%0.2f/E_CMS);\n' %(self.ckkw_param))
        if self.nlojets>=0:
            runfile.write('    NLO_QCD_PART BVIRS{')
            ljmin = self.process.fsparts
            ljmax = ljmin+self.nlojets
            for ljet in range (ljmin,ljmax):
                runfile.write('%i,' %(ljet))
            runfile.write('%i};\n' %(ljmax))
        lineno = 0
        runfile.write('\n')
        runfile.write('}(processes)\n\n')
        
        runfile.write('(model){\n')
        for line in self.process.modellines:
            runfile.write('  '+line)
            lineno += 1
            runfile.write('\n')
        runfile.write('}(model)\n')




class PI(object):
    def __init__(self,
                 declarationlines="",modellines="",
                 totjets=0,nlojets=0,fsparts=0):
        gobject.gobject.__init__(self)
        self.declarationlines = declarationlines
        self.modellines       = modellines
        self.totjets          = totjets
        self.nlojets          = nlojets
        self.fsparts          = fsparts

    def getList(self):
        return [self.declarationlines,self.modellines,
                self.totjets,self.nlojets,self.fsparts]

