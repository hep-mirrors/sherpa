#!/usr/bin/python

import pygtk
import gtk

import GUI_BaseClass as guibase
import ParameterBox
import SelectorBox

class Selsection_gui(guibase.gui_object):
    def initialise(self):
        self.selbox   = self.parameters.getSelectorBox()
        self.label    = gtk.Label("No selectors for chosen process.")
        self.trows    = 8
        self.tcols    = 3
        self.boxtable = gtk.Table(self.tcols,self.trows,False)
        for col in range(0,self.tcols):
            self.boxtable.set_col_spacing(col,30)
        for row in range(0,self.trows):
            self.boxtable.set_row_spacing(row,30)
        self.resetTable()
        self.boxtable.show()

    def getContent(self):
        vbox          = gtk.VBox(False,0)
        vbox.set_border_width(10)
        hbox          = gtk.HBox(False,0)
        hbox.set_border_width(10)
        vbox.pack_start(hbox,False,False,2)
        hbox          = gtk.HBox(False,0)
        hbox.pack_start(self.boxtable,False,False,2)
        vbox.pack_start(hbox,False,False,2)
        return vbox

    def updateOptions(self,collider,process,minjets):
        self.selbox.initialiseDefaults(collider,process,minjets)
        row     = 0
        flag    = False
        print "SelectorGUI::updateOptions: reset table"
        print len(self.selbox.getSelectorList())," selectors to be added."
        self.resetTable()
        for sel in self.selbox.getSelectorList():
            print "Deal with ",sel.getTag()
            flag1 = self.addSelector(row,sel)
            flag  = flag or flag1
            row   = row+1
            self.boxtable.set_row_spacing(row,50)
        if not(flag):
            self.boxtable.attach(self.label,0,self.tcols-1,
                                 int(self.trows/2),int(self.trows/2)+1)


    def resetTable(self):
        children = self.boxtable.get_children()
        for child in children:
            self.boxtable.remove(child)

    def switchedOn(self,button,sel,pbox):
        print "will switch ",sel.getTag(),sel.isOn()
        sel.switch()
        print "switched ",sel.getTag(),sel.isOn()
        if sel.isOn():
            button.set_label("On")
            pbox.set_sensitive(True)
        else:
            button.set_label("Off")
            pbox.set_sensitive(False)

    def selectMode(self,button,data):
        print "select mode for ",data[0].getTag(),data[1]
        data[0].setValue("Mode",data[1])

    def paramChanged(self,button,data):
        print "change params for ",data[0].getTag(),", param=",data[1]
        data[0].setValue(data[1],button.get_value())

    def addSelector(self,row,sel):
        tag    = sel.getTag()
        print tag,len(sel.getRanges())
        if tag=="NJetFinder":
            return self.connectJetFinder(row,sel)
        if tag=="Mass2":
            return self.connectRangeSelector(row,sel)
        return False

    def buildBox(self,label,sel,row):
        box    = gtk.VBox(False,0)
        switch = gtk.Button()
        if sel.isOn():
            switch.set_label("On")
            box.set_sensitive(True)
        else:
            switch.set_label("Off")
            box.set_sensitive(False)
        switch.connect("clicked",self.switchedOn,sel,box)
        self.boxtable.attach(label,0,1,row,row+1)
        self.boxtable.attach(switch,1,2,row,row+1)
        label.show()
        switch.show()
        box.show()
        return box

    def connectRangeSelector(self,row,sel):
        label    = gtk.Label("Mass["+str(sel.getParticleIds()[0])+","+
                             str(sel.getParticleIds()[1])+"]")
        ranges   = sel.getRanges()
        mini     = ranges["min"]
        maxi     = ranges["max"]
        param   = gtk.HBox(False,0)
        print "adjustment for max: ",maxi[2]," in [",maxi[0],maxi[1],"]"
        sel.setValue("max",maxi[2])
        maxadj   = gtk.Adjustment(maxi[2],maxi[0],maxi[1],maxi[3])
        maxsel   = gtk.SpinButton(maxadj,1,2)
        maxsel.show()
        maxsel.connect("value_changed",self.paramChanged,[sel,"max"])
        param.pack_end(maxsel,False,False,2)
        param.pack_end(gtk.Label("Max:"),False,False,2)            
        print "adjustment for min: ",mini[2]," in [",mini[0],mini[1],"]"
        sel.setValue("min",mini[2])
        minadj   = gtk.Adjustment(mini[2],mini[0],mini[1],mini[3])
        minsel   = gtk.SpinButton(minadj,1,2)
        minsel.show()
        minsel.connect("value_changed",self.paramChanged,[sel,"min"])
        param.pack_end(minsel,False,False,2)
        param.pack_end(gtk.Label("Min:"),False,False,2)            
        param.show_all()

        params   = self.buildBox(label,sel,row)
        params.pack_start(param,False,False,2)
        self.boxtable.attach(params,2,3,row,row+1)
        return True

    def connectJetFinder(self,row,sel):
        label    = gtk.Label("Jet finder")

        param    = gtk.HBox(False,0)
        param.show()
        params   = self.buildBox(label,sel,row)
        params.pack_start(param,False,False,2)
        param.pack_start(gtk.Label("Algorithm:"),False,False,2)
        algo     = gtk.RadioButton(None,"KT")
        algo.connect("toggled",self.selectMode,[sel,"KT"])
        algo.show()
        param.pack_start(algo,False,False,2)
        algo     = gtk.RadioButton(algo,"Anti-KT")
        algo.connect("toggled",self.selectMode,[sel,"Anti-KT"])
        algo.show()
        param.pack_start(algo,False,False,2)
        algo     = gtk.RadioButton(algo,"SISCone")
        algo.connect("toggled",self.selectMode,[sel,"SISCone"])
        algo.show()
        param.pack_start(algo,False,False,2)

        ranges   = sel.getRanges()
        PT       = ranges["PT"]
        R        = ranges["R"]
        param   = gtk.HBox(False,0)
        param.show()
        params.pack_start(param,False,False,2)
        print "adjustment for R: ",R[2]," in [",R[0],R[1],"]"
        sel.setValue("R",R[2])
        Radj     = gtk.Adjustment(R[2],R[0],R[1],R[3])
        Rsize    = gtk.SpinButton(Radj,1,2)
        Rsize.show()
        Rsize.connect("value_changed",self.paramChanged,[sel,"R"])
        param.pack_end(Rsize,False,False,2)
        param.pack_end(gtk.Label("Cone size R:"),False,False,2)            
        print "adjustment for PT: ",PT[2]," in [",PT[0],PT[1],"]"
        sel.setValue("PT",PT[2])
        ptadj    = gtk.Adjustment(PT[2],PT[0],PT[1],PT[3])
        ptsize   = gtk.SpinButton(ptadj,1,0)
        ptsize.show()
        ptsize.connect("value_changed",self.paramChanged,[sel,"PT"])
        param.pack_end(ptsize,False,False,2)
        param.pack_end(gtk.Label("PT:"),False,False,2)            
        param.show_all()
        
        self.boxtable.attach(params,2,3,row,row+2)
        row = row+1
        return True
