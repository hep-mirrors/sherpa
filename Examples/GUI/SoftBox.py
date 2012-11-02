#!/usr/bin/python

import pygtk
import gtk


class SoftBox:
    def __init__(self):
        self.had       = False
        self.hadmodels = ["None","Ahadic","Lund"]
        self.hadmodel  = "None"
        self.hadmodes  = ["Low","Normal","High"]
        self.hadmode   = "Normal"
        self.mpi       = False
        self.mpimodels = ["None","Amisic"]
        self.mpimodel  = "None"
        self.mpimodes  = ["Low","Normal","High"]
        self.mpimode   = "Normal"
        pass

    def initialiseDefaults(self,collider):
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
        self.hadmodel = hadmodel

    def getHadModes(self):
        return self.hadmodes

    def getHadMode(self):
        return self.hadmode

    def setHadMode(self,hadmodel):
        self.hadmode = hadmode


    def getMPIOn(self):
        return self.mpi

    def getMPIModels(self):
        return self.mpimodels

    def getMPIModel(self):
        return self.mpimodel

    def setMPIModel(self,mpimodel):
        self.mpimodel = mpimodel

    def getMPIModes(self):
        return self.mpimodes

    def getMPIMode(self):
        return self.mpimode

    def setMPIMode(self,mpimodel):
        self.mpimode = mpimode


    def write(self,runfile,collider):
        runfile.write('\n')
        runfile.write('%%% Soft setup \n\n')
