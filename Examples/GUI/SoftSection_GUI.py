#!/usr/bin/python

import pygtk
import gtk

import GUI_BaseClass as guibase
import ParameterBox
import SoftBox

class Softsection_gui(guibase.gui_object):
    def initialise(self):
        self.softbox = self.parameters.getSoftBox()
        self.initButtons()
        pass

    def getContent(self):
        table     = gtk.Table(2,1,False)
        table.attach(self.hadbox,0,1,0,1,gtk.FILL,gtk.FILL,0,0)
        table.attach(self.mpibox,0,1,1,2,gtk.FILL,gtk.FILL,0,0)
        table.show()
        return table

    def updateOptions(self):
        self.hadbox.set_sensitive(self.softbox.getHadOn())
        if self.hadbox.is_sensitive():
            hadmodel=self.softbox.getHadModel()
            if self.softbox.hasHadTunes(hadmodel):
                self.hadmodebox.set_sensitive(True)
            else:
                self.hadmodebox.set_sensitive(False)
        self.mpibox.set_sensitive(self.softbox.getMPIOn())
        if self.mpibox.is_sensitive():
            mpimodel=self.softbox.getMPIModel()
            if self.softbox.hasMPITunes(mpimodel):
                self.mpimodebox.set_sensitive(True)
            else:
                self.mpimodebox.set_sensitive(False)

    def buttonChanged(self,button,data):
        if (data[1]=="HadModel"):
            self.softbox.setHadModel(data[0])
        if (data[1]=="HadMode"):
            self.softbox.setHadMode(data[0])
        if (data[1]=="MPIModel"):
            self.softbox.setMPIModel(data[0])
        if (data[1]=="MPIMode"):
            self.softbox.setMPIMode(data[0])
        self.updateOptions()

    def initButtons(self):
        hadoptions = self.softbox.getHadModels()
        hadmodes   = self.softbox.getHadModes()
        haddefault = self.softbox.getHadModel()
        hadon      = self.softbox.getHadOn()
        self.hadbox = gtk.VBox(False,10)
        self.hadbox.set_border_width(10) 
        self.hadbox.pack_start(gtk.Label("Hadronisation"),False,False,2)
        self.hadmodelbox = gtk.HBox(False,10)
        self.hadmodelbox.set_border_width(10) 
        self.hadmodelbox.pack_start(gtk.Label("Model"),False,False,2)
        self.hadbox.pack_start(self.hadmodelbox,False,False,2)
        hadbutton = None
        for had in hadoptions:
            hadbutton = gtk.RadioButton(hadbutton,had[0])
            hadbutton.connect("toggled", self.buttonChanged, 
                              [had[0],"HadModel"])
            self.hadmodelbox.pack_start(hadbutton,False,False,2)
            if (had==haddefault):
                hadbutton.set_active(True)
        self.hadmodebox = gtk.HBox(False,10)
        self.hadmodebox.set_border_width(10) 
        self.hadmodebox.pack_start(gtk.Label("Activity"),False,False,2)
        self.hadbox.pack_start(self.hadmodebox,False,False,2)
        hadbutton = None
        for had in hadmodes:
            hadbutton = gtk.RadioButton(hadbutton,had)
            hadbutton.connect("toggled", self.buttonChanged, 
                              [had,"HadMode"])
            self.hadmodebox.pack_start(hadbutton,False,False,2)
            if (had=="Normal"):
                hadbutton.set_active(True)
        self.hadbox.set_sensitive(hadon)

        mpioptions = self.softbox.getMPIModels()
        mpimodes   = self.softbox.getMPIModes()
        mpidefault = self.softbox.getMPIModel()
        mpion      = self.softbox.getMPIOn()
        self.mpibox = gtk.VBox(False,10)
        self.mpibox.set_border_width(10) 
        self.mpibox.pack_start(gtk.Label("Multiple Interactions"),False,False,2)
        self.mpimodelbox = gtk.HBox(False,10)
        self.mpimodelbox.set_border_width(10) 
        self.mpimodelbox.pack_start(gtk.Label("Model"),False,False,2)
        self.mpibox.pack_start(self.mpimodelbox,False,False,2)
        mpibutton = None
        for mpi in mpioptions:
            mpibutton = gtk.RadioButton(mpibutton,mpi[0])
            mpibutton.connect("toggled", self.buttonChanged, 
                              [mpi[0],"MPIModel"])
            self.mpimodelbox.pack_start(mpibutton,False,False,2)
            if (mpi==mpidefault):
                mpibutton.set_active(True)
        self.mpimodebox = gtk.HBox(False,10)
        self.mpimodebox.set_border_width(10) 
        self.mpimodebox.pack_start(gtk.Label("Activity"),False,False,2)
        self.mpibox.pack_start(self.mpimodebox,False,False,2)
        mpibutton = None
        for mpi in mpimodes:
            mpibutton = gtk.RadioButton(mpibutton,mpi)
            mpibutton.connect("toggled", self.buttonChanged, 
                              [mpi,"MPIMode"])
            self.mpimodebox.pack_start(mpibutton,False,False,2)
            if (mpi=="Normal"):
                mpibutton.set_active(True)
        self.mpibox.set_sensitive(mpion)

