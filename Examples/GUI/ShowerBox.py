#!/usr/bin/python

import pygtk
import gtk


class ShowerBox:
    def __init__(self):
        self.fs_pt2min = 1.
        self.is_pt2min = 1.
        self.fs_asfac  = 1.
        self.is_asfac  = 1.
        self.kinscheme = 0
        self.ison      = True

    def initialiseDefaults(self,nlo):
        print "Shower initdefaults(",nlo,")"
        if not nlo:
            self.fs_asfac  = 0.73
            self.is_asfac  = 1.27
            self.kinscheme = 0
        else:
            self.fs_asfac  = 1.
            self.is_asfac  = 1.
            self.kinscheme = 1
            
            
    def write(self,runfile):
        runfile.write("  %%% Shower settings\n\n")
        if not self.ison:
            runfile.write("  SHOWER_HANDLER = None")
        else:
            runfile.write("  CSS_FS_PT2MIN = %.2f;CSS_IS_PT2MIN = %.2f;\n"
                          %(self.fs_pt2min,self.is_pt2min))
            runfile.write("  CSS_FS_AS_FAC = %.2f;CSS_IS_AS_FAC = %.2f;\n"
                          %(self.fs_asfac,self.is_asfac))
            runfile.write("  CSS_KIN_SCHEME = "+str(self.kinscheme)+"\n")

    def getParams(self):
        return self.fs_pt2min,self.is_pt2min,self.fs_asfac,self.is_asfac

    def setFS_PT2Min(self,fspt2min):
        self.fs_pt2min = fspt2min

    def setIS_PT2Min(self,ispt2min):
        self.is_pt2min = ispt2min

    def setFS_AsFac(self,fsasfac):
        self.fs_asfac  = fsasfac

    def setIS_AsFac(self,isasfac):
        self.is_asfac  = isasfac

    def setKinScheme(self,kinscheme):
        self.kinscheme = kinscheme

    def getKinScheme(self):
        return self.kinscheme

    def switch(self):
        self.ison = not self.ison

    def isOn(self):
        return self.ison
