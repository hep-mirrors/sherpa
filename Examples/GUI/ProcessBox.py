#!/usr/bin/python

import pygtk
import gtk

class ProcessBox:
    def __init__(self):
        self.LCmap = {}
        self.HCmap = {}
        self.fillMaps()
        
        self.NC = gtk.TreeStore(str, bool)
        self.NC.append(None,["None",False])

        self.LC = gtk.TreeStore(str, bool)
        self.fillLCStore()

        self.HC = gtk.TreeStore(str, bool)
        self.fillHCStore()
      
        self.type        = "None"
        self.process_tag = ""
        self.process     = None
        self.processes   = self.NC
        self.minjets     = 0
        self.totjets     = 0
        self.nlojets     = -1
        self.nlomax      = 3 
        self.ckkw_param  = 0.
        self.muFfactor   = 1.
        self.muRfactor   = 1.
        self.muQfactor   = 1.
        self.LOgen       = "Internal"
        self.loopgenerators = []
        for i in range(0,self.nlomax+1): 
            self.loopgenerators.append("None")

    def initialiseDefaults(self,collider):
        print "Proc::initialiseDefaults",collider[2],collider[4]
        self.processes = self.NC
        self.type      = "None"
        if ((collider[2]==2212 or collider[2]==-2212) and
            (collider[4]==2212 or collider[4]==-2212)):
            print "   select HC"
            self.processes = self.HC
            self.type      = "HC"
            return
        if ((collider[2]==11 or collider[2]==-11) and
            (collider[4]==11 or collider[4]==-11)):
            print "   select LC"
            self.processes = self.LC
            self.type      = "LC"
            return
        print"   select nothing"

    def getOptions(self):
        return self.processes

    def getProcess(self):
        return self.process

    def setProcess(self,process_tag):
        self.process     = None
        self.process_tag = process_tag
        if self.processes==self.HC:
            self.process = self.HCmap[process_tag]
        if self.processes==self.LC:
            self.process = self.LCmap[process_tag]
        if self.process==None:
            print "Set process = ",process_tag," --> not found."
            return
        print "Set process = ",process_tag," --> successful."
            
    def getType(self):
        return self.type

    def getNLOmax(self):
        return self.nlomax

    def getProcNJets(self):
        if self.process!=None:
            minj = self.process.minjets
            totj = self.process.totjets
            nloj = self.process.nlojets
            return minj,totj,nloj
        return None,None,None

    def setJetMultis(self,minjets=0,totjets=0,nlojets=0):
        if self.process!=None:
            if totjets>=self.process.minjets and totjets<=self.process.totjets:
                self.totjets = totjets
            if minjets>=self.process.minjets and minjets<=self.totjets:
                self.minjets = minjets
            if nlojets>=self.minjets and nlojets<=self.totjets:
                self.nlojets = nlojets

    def getMinJets(self):
        if self.process!=None:
            return self.minjets
        return 0

    def getNJets(self):
        if self.process!=None:
            return self.minjets,self.totjets,self.nlojets
        return None,None,None

    def setScaleFactors(self,muF,muR,muQ):
        self.muFfactor = muF
        self.muRfactor = muR
        self.muQfactor = muQ

    def getScaleFactors(self):
        return self.muFfactor,self.muRfactor,self.muQfactor

    def isNLO(self):
        return (self.process.nlojets>-1)

    def setCKKW(self,ckkw_param):
        self.ckkw_param = ckkw_param

    def setLOgen(self,gen):
        self.LOgen = gen

    def getLOGen(self):
        return self.LOgen

    def setLoopGen(self,n,loopgen):
        if n<len(self.loopgenerators):
            self.loopgenerators[n] = loopgen
            print "Set loopgen[",n,"] = ",loopgen

    def getLoopGens(self):
        print "Have ",len(self.loopgenerators)," LoopGens."
        return self.loopgenerators

    def write(self,runfile):
        if self.process==None:
            runfile.write("\n")
            runfile.write("%%% Process setup: none selected\n\n")
            return
        print "Write process stuff "
        runfile.write("\n")
        runfile.write("%%% Process setup: "+self.process_tag+"\n\n")
        runfile.write("(processes){\n")
        lineno = 0
        for line in self.process.declarationlines:
            if lineno==0:
                if self.minjets>0:
                    pos     = line.find("{NJETS}")
                    linel   = list(line)
                    print "**************************",pos,line
                    for i in range(0,self.minjets):
                        linel.insert(pos," 93")
                    line    = "".join(linel) 
                    print "**************************",pos,line
                runfile.write("  Process "+line)
            else: 
                runfile.write("    Decay "+line)
            runfile.write("\n")
            lineno += 1
        if self.totjets>0:
            runfile.write("    CKKW sqr(%0.2f/E_CMS);\n" %(self.ckkw_param))
        if self.nlojets>=0:
            runfile.write("    NLO_QCD_Mode = MC@NLO {")
            ljmin = self.process.fsparts
            ljmax = ljmin+self.nlojets
            for ljet in range (ljmin,ljmax):
                runfile.write("%i," %(ljet))
            runfile.write("%i};\n" %(ljmax))
            for n in range(0,len(self.loopgenerators)):
                gen = self.loopgenerators[n]
                if gen!="None":
                    runfile.write("    Loop_Generator = LOOPGEN%s {%s}\n" 
                                  %(str(n),str(n+self.process.fsparts)))
        runfile.write("}(processes)\n\n")
        
    def fillMaps(self):
        # Format for input:
        # [[<process_tag>,<decay_tags>],maxjets,nlojets,extra_FS_particles]
        # use -1,-1,-1 to steer alternatives, may have to use
        # the read-in facilities of Sherpa in a smarter way
        #
        self.LCmap["Jets"] = PI(
            ["11 -11 -> 93 93 93{NJETS};"],["MODEL = SM;"],6,4,2,2)

        self.HCmap["MinimumBias"] = PI(
            ["MinimumBias"],[""],-1,-1,-1)

        self.HCmap["Jets"] = PI(
            ["93 93 -> 93 93 93{NJETS};"],["MODEL = SM;"],4,2,2,2)

        self.HCmap["Gamma + Jets"] = PI(
            ["93 93 -> 22 93 93{NJETS};"],["MODEL = SM;"],4,2,1,2)

        self.HCmap["(W^+ -> e^+ nu) + Jets"] = PI(
            ["93 93 -> 24[a] 93{NJETS};","24[a] -> -11 12;"],
            ["MODEL = SM;"],5,2,0,2)
        self.HCmap["(W^+ -> mu^+ nu) + Jets"] = PI(
            ["93 93 -> 24[a] 93{NJETS};","24[a] -> -13 14;"],
            ["MODEL = SM;"],5,2,0,2)
        self.HCmap["(W^+ -> tau^+ nu) + Jets"] = PI(
            ["93 93 -> 24[a] 93{NJETS};","24[a] -> -15 16;"],
            ["MODEL = SM;"],5,2,0,2)
        self.HCmap["(W^- -> e^- nu) + Jets"] = PI(
            ["93 93 -> -24[a] 93{NJETS};","-24[a] -> 11 -12;"],
            ["MODEL = SM;"],5,2,0,2)
        self.HCmap["(W^- -> mu^- nu) + Jets"] = PI(
            ["93 93 -> -24[a] 93{NJETS};","-24[a] -> 13 -14;"],
            ["MODEL = SM;"],5,2,0,2)
        self.HCmap["(W^- -> tau^- nu) + Jets"] = PI(
            ["93 93 -> -24[a] 93{NJETS};","-24[a] -> 15 -16;"],
            ["MODEL = SM;"],5,2,0,2)

        self.HCmap["(Z -> e^- e^+) + Jets"] = PI(
            ["93 93 -> 23[a] 93{NJETS};","23[a] -> 11 -11;"],
            ["MODEL = SM;"],5,2,0,2)
        self.HCmap["(Z -> mu^- mu^+) + Jets"] = PI(
            ["93 93 -> 23[a] 93{NJETS};","23[a] -> 13 -13;"],
            ["MODEL = SM;"],5,2,0,2)
        self.HCmap["(Z -> tau^- tau^+) + Jets"] = PI(
            ["93 93 -> 23[a] 93{NJETS};","23[a] -> 15 -15;"],
            ["MODEL = SM;"],5,2,0,2)
        self.HCmap["(Z -> nu nu) + Jets"] = PI(
            ["93 93 -> 23[a] 93{NJETS};","23[a] -> 11 -11;"],
            ["MODEL = SM;"],5,2,0,2)

        self.HCmap["(e^- e^+) + Jets"] = PI(
            ["93 93 -> 11 -11 93{NJETS};"],["MODEL = SM;"],5,2,0,2)
        self.HCmap["(mu^- mu^+) + Jets"] = PI(
            ["93 93 -> 13 -13 93{NJETS};"],["MODEL = SM;"],5,2,0,2)
        self.HCmap["(tau^- tau^+) + Jets"] = PI(
            ["93 93 -> 15 -15 93{NJETS};"],["MODEL = SM;"],5,2,0,2)
        
    def fillLCStore(self):
        parent = self.LC.append(None,["Hard QCD",False])
        self.LC.append(parent,["Jets",True])

    def fillHCStore(self):
        parent = self.HC.append(None,["Minimum Bias",False])
        self.HC.append(parent,["MinimumBias",True])

        parent = self.HC.append(None,["Hard QCD",False])
        self.HC.append(parent,["Jets",True])

        parent = self.HC.append(None,["V + jets",False])
        self.HC.append(parent,["Gamma + Jets",True])
        daughter = self.HC.append(parent,["W + Jets",False])
        self.HC.append(daughter,["(W^+ -> e^+ nu) + Jets",True])
        self.HC.append(daughter,["(W^+ -> mu^+ nu) + Jets",True])
        self.HC.append(daughter,["(W^+ -> tau^+ nu) + Jets",True])
        self.HC.append(daughter,["(W^- -> e^- nu) + Jets",True])
        self.HC.append(daughter,["(W^- -> mu^- nu) + Jets",True])
        self.HC.append(daughter,["(W^- -> tau^- nu) + Jets",True])
        self.HC.append(daughter,["(W^{+,-} -> mu^{+,-} nu) + Jets",True])
        daughter = self.HC.append(parent,["Z + Jets",False])
        self.HC.append(daughter,["(Z -> e^- e^+) + Jets",True])
        self.HC.append(daughter,["(Z -> mu^- mu^+) + Jets",True])
        self.HC.append(daughter,["(Z -> tau^- tau^+) + Jets",True])
        self.HC.append(daughter,["(Z -> nu nu) + Jets",True])
        daughter = self.HC.append(parent,["Z/gamma^* + Jets",False])
        self.HC.append(daughter,["(e^- e^+) + Jets",True])
        self.HC.append(daughter,["(mu^- mu^+) + Jets",True])
        self.HC.append(daughter,["(tau^- tau^+) + Jets",True])

        parent = self.HC.append(None,["H+jets (gluon fusion)",False])
        self.HC.append(parent,["(H -> gamma gamma) + Jets",True])
        self.HC.append(parent,["(H -> tau tau) + Jets",True])
        self.HC.append(parent,["(H -> e^+ nu mu^- nu) + Jets",True])
        self.HC.append(parent,["(H -> e^+ e^- mu^+ mu^-) + Jets",True])
        self.HC.append(parent,["(H -> e^+ e^- nu nu) + Jets",True])





class PI():
    def __init__(self,
                 declarationlines="",modellines="",
                 totjets=0,nlojets=0,minjets=0,fsparts=0):
        self.declarationlines = declarationlines
        self.modellines       = modellines
        self.totjets          = totjets
        self.nlojets          = nlojets
        self.minjets          = minjets
        self.fsparts          = fsparts

    def getList(self):
        return [self.declarationlines,self.modellines,
                self.totjets,self.nlojets,self.fsparts]

    def getModelLines(self):
        return self.modellines
