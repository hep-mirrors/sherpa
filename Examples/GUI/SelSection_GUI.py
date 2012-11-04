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
            self.boxtable.set_col_spacing(col,20)
        for row in range(0,self.trows):
            self.boxtable.set_row_spacing(row,20)
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

    def updateOptions(self,collider,process):
        self.selbox.initialiseDefaults(collider,process)
        row     = 0
        flag    = False
        print "SelectorGUI::updateOptions: reset table"
        self.resetTable()
        for sel in self.selbox.getSelectorList():
            flag = flag or self.addSelector(row,sel)
            row = row+1
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

    def connectJetFinder(self,row,sel):
        label    = gtk.Label("Mass[",sel.getParticleIds()[0],",",
                             sel.getParticleIds()[1],"]")
        params   = gtk.VBox(False,0)
        switch   = gtk.Button()
        if sel.isOn():
            switch.set_label("On")
            params.set_sensitive(True)
        else:
            switch.set_label("Off")
            params.set_sensitive(False)
        switch.connect("clicked",self.switchedOn,sel,params)
        param    = gtk.HBox(False,0)
        param.show()
        params.pack_start(param,False,False,2)
        label.show()
        switch.show()
        params.show()
        self.boxtable.attach(label,0,1,row,row+1)
        self.boxtable.attach(switch,1,2,row,row+1)
        self.boxtable.attach(params,2,3,row,row+2)
        return True

    def connectJetFinder(self,row,sel):
        label    = gtk.Label("Jet finder")
        params   = gtk.VBox(False,0)
        switch   = gtk.Button()
        if sel.isOn():
            switch.set_label("On")
            params.set_sensitive(True)
        else:
            switch.set_label("Off")
            params.set_sensitive(False)
        switch.connect("clicked",self.switchedOn,sel,params)

        param    = gtk.HBox(False,0)
        param.show()
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
        Radj     = gtk.Adjustment(R[2],R[0],R[1],0.1)
        Rsize    = gtk.SpinButton(Radj,1,2)
        Rsize.show()
        Rsize.connect("value_changed",self.paramChanged,[sel,"R"])
        param.pack_end(Rsize,False,False,2)
        param.pack_end(gtk.Label("Cone size R:"),False,False,2)            
        print "adjustment for PT: ",PT[2]," in [",PT[0],PT[1],"]"
        sel.setValue("PT",PT[2])
        ptadj    = gtk.Adjustment(PT[2],PT[0],PT[1],5.0)
        ptsize   = gtk.SpinButton(ptadj,1,0)
        ptsize.show()
        ptsize.connect("value_changed",self.paramChanged,[sel,"PT"])
        param.pack_end(ptsize,False,False,2)
        param.pack_end(gtk.Label("PT:"),False,False,2)            
        param.show_all()
        
        label.show()
        switch.show()
        params.show()
        self.boxtable.attach(label,0,1,row,row+1)
        self.boxtable.attach(switch,1,2,row,row+1)
        self.boxtable.attach(params,2,3,row,row+2)
        row = row+1
        return True

        #self.sellist.append(["NJetFinder",True,
        #                    ["",[0.,self.E_CMS],[0.,2.0]]])
        #self.sellist.append(["PT 93",True,[0.,self.E_CMS]])
        #self.sellist.append(["Eta 93",True,[-10.,10.]])
