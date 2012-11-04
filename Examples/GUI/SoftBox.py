#!/usr/bin/python

import pygtk
import gtk


class SoftBox:
    def __init__(self):
        self.initHadModels()
        self.initMPIModels()
        pass

    def initialiseDefaults(self,collider):
        self.hadmodel="Off"
        self.hadmode ="Normal"
        self.mpimodel="None"
        self.mpimode ="Normal"
        if ((collider[2]==11 or collider[2]==-11) and
            (collider[4]==11 or collider[4]==-11)):
            self.had = True
            self.mpi = False
        if ((collider[2]==2212 or collider[2]==-2212) and
            (collider[4]==2212 or collider[4]==-2212)):
            self.had = True
            self.mpi = True

    def getHadOn(self):
        return self.had

    def getHadModels(self):
        return self.hadmodels

    def getHadModel(self):
        return self.hadmodel

    def setHadModel(self,hadmodel):
        print "Set hadmodel = ",hadmodel
        self.hadmodel = hadmodel

    def getHadModes(self):
        return self.hadmodes

    def getHadMode(self):
        return self.hadmode

    def setHadMode(self,hadmode):
        print "Set hadmode = ",hadmode
        self.hadmode = hadmode

    def hasHadTunes(self,hadmodel):
        for had in self.hadmodels:
            if had[0]==hadmodel:
                return had[1]

    def getMPIOn(self):
        return self.mpi

    def getMPIModels(self):
        return self.mpimodels

    def getMPIModel(self):
        return self.mpimodel

    def setMPIModel(self,mpimodel):
        print "Set mpimodel = ",mpimodel
        self.mpimodel = mpimodel

    def getMPIModes(self):
        return self.mpimodes

    def getMPIMode(self):
        return self.mpimode

    def setMPIMode(self,mpimode):
        print "Set mpimode = ",mpimode
        self.mpimode = mpimode

    def hasMPITunes(self,mpimodel):
        for mpi in self.mpimodels:
            if mpi[0]==mpimodel:
                return mpi[1]


    def write(self,runfile):
        runfile.write("\n")
        runfile.write("%%% Soft setup \n\n")
        if self.hadmodel!="None" and self.had:
            runfile.write("FRAGMENTATION      = %s\n" %str(self.hadmodel))
            runfile.write("FRAGMENTATION_TUNE = %s\n" %str(self.hadmode))
        else:
            runfile.write("FRAGMENTATION      = Off\n" %str(self.hadmodel))
        runfile.write("MI_HANDLER    = %s\n" %str(self.mpimodel)) 
        if self.mpimodel!="None" and self.mpi:
            runfile.write("MI_TUNE       = %s\n" %str(self.mpimode))  
        runfile.write("\n")
        
    def initHadModels(self):
        self.had       = False
        self.hadmodels = []
        self.hadmodels.append(["None",False])
        self.hadmodels.append(["Ahadic",True])
        self.hadmodels.append(["Lund",False])
        self.hadmodel  = "None"
        self.hadmodes  = ["Low","Normal","High"]
        self.hadmode   = "Normal"

    def initMPIModels(self):
        self.mpi       = False
        self.mpimodels = []
        self.mpimodels.append(["None",False])
        self.mpimodels.append(["Amisic",True])
        self.mpimodel  = "None"
        self.mpimodes  = ["Low","Normal","High"]
        self.mpimode   = "Normal"
