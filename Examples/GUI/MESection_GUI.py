#!/usr/bin/python

import pygtk
import gtk

import GUI_BaseClass as guibase
import ParameterBox
import ProcessBox

class MEsection_gui(guibase.gui_object):
    def initialise(self):
        self.procbox  = self.parameters.getProcBox()
        nbuttons      = self.procbox.getNLOmax()
        gens          = ["Internal","Amegic","Comix"]
        loopgens      = ["Internal","BlackHat","GoSam","OpenLoops"]
        self.initButtons(nbuttons,gens,loopgens)

    def initButtons(self,nbuttons,gens,loopgens):
        self.genbox = gtk.VBox(False,10)
        self.genbox.pack_start(gtk.Label("Tree-level generator"))
        self.genbox.set_border_width(10) 
        box    = gtk.HBox(False,10)
        box.set_border_width(10) 
        self.genbutton = None
        for gen in gens:
            self.genbutton = gtk.RadioButton(self.genbutton,gen)
            box.pack_start(self.genbutton,False,False,2)
            if (gen=="Comix"):
                self.genbutton.set_active(True)
        self.genbox.pack_start(box)

        self.loopgenbox   = gtk.VBox(False,10)
        self.loopgenfield = []
        self.loopgenbox.pack_start(gtk.Label("Loop generators"))
        self.loopgenbox.set_border_width(10) 
        self.loopgenbutton = []
        for i in range (0,nbuttons+1):
            box    = gtk.HBox(False,10)
            box.set_border_width(10) 
            button = gtk.RadioButton(None,"None")
            button.connect("toggled", self.loopgenChanged, ["None",i])
            box.pack_start(gtk.Label("%s extra jets:" %str(i)))
            box.pack_start(button,False,False,2)
            for gen in loopgens:
                button = gtk.RadioButton(button,gen)
                button.connect("toggled", self.loopgenChanged, [gen,i])
                box.pack_start(button,False,False,2)
            box.set_sensitive(False)
            self.loopgenfield.append(box)
            self.loopgenbutton.append(button)
            self.loopgenbox.pack_start(box)
    
    def loopgenChanged(self,button,data):
        print "Set loopgen[",data[1],"] = ",data[0]
        self.procbox.setLoopGen(data[1],data[0])

    def getContent(self):
        table     = gtk.Table(2,1,False)
        table.attach(self.genbox,0,1,0,1,gtk.FILL,gtk.FILL,0,0)
        table.attach(self.loopgenbox,0,1,1,2,gtk.FILL,gtk.FILL,0,0)
        table.show()
        return table

    def updateOptions(self):
        totjets,nlojets = self.procbox.getNJets()
        if totjets==None or nlojets==None:
            return
        print "MEs::updateOptions()",totjets,nlojets
        for field in self.loopgenfield:
            field.set_sensitive(False)
        for i in range(0,nlojets+1):
            self.loopgenfield[i].set_sensitive(True)

    def extractParameters(self):
        pass
