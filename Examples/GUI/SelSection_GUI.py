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
        self.trows    = 12
        self.tcols    = 3
        self.boxtable = gtk.Table(self.tcols,self.trows,False)
        for col in range(0,self.tcols):
            self.boxtable.set_col_spacing(col,30)
        for row in range(0,self.trows):
            self.boxtable.set_row_spacing(row,40)
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
            #self.boxtable.set_row_spacing(row,50)
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
        if tag=="IsolationCut":
            return self.connectIsolationCut(row,sel)
        if tag=="Mass2" or tag=="PT" or tag=="Rapidity" or tag=="DeltaR":
            return self.connectRangeSelector(row,sel)
        return False

    def buildBox(self,label,sel,row):
        box    = gtk.VBox(False,2)
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
        if sel.getTag()=="Mass2":
            label    = gtk.Label("Mass["+str(sel.getParticleIds()[0])+","+
                                 str(sel.getParticleIds()[1])+"]")
        if sel.getTag()=="DeltaR":
            label    = gtk.Label("DeltaR["+str(sel.getParticleIds()[0])+","+
                                 str(sel.getParticleIds()[1])+"]")
        if sel.getTag()=="PT":
            label    = gtk.Label("PT["+str(sel.getParticleIds()[0])+"]")
        if sel.getTag()=="Rapidity":
            label    = gtk.Label("Y["+str(sel.getParticleIds()[0])+"]")

        params   = self.buildBox(label,sel,row)
        ranges   = sel.getRanges()
        mini     = ranges["min"]
        maxi     = ranges["max"]

        box      = gtk.HBox(False,0)
        param    = gtk.HBox(False,0)
        print "adjustment for max: ",maxi[2]," in [",maxi[0],maxi[1],"]"
        sel.setValue("max",maxi[2])
        maxadj   = gtk.Adjustment(maxi[2],maxi[0],maxi[1],maxi[3])
        maxsel   = gtk.SpinButton(maxadj,1,2)
        maxsel.show()
        maxsel.connect("value_changed",self.paramChanged,[sel,"max"])
        param.pack_end(maxsel,False,False,2)
        param.pack_start(gtk.Label("Max:"),False,False,2)            
        param.show_all()
        box.pack_end(param)

        param    = gtk.HBox(False,0)
        print "adjustment for min: ",mini[2]," in [",mini[0],mini[1],"]"
        sel.setValue("min",mini[2])
        minadj   = gtk.Adjustment(mini[2],mini[0],mini[1],mini[3])
        minsel   = gtk.SpinButton(minadj,1,2)
        minsel.show()
        minsel.connect("value_changed",self.paramChanged,[sel,"min"])
        param.pack_end(minsel,False,False,2)
        param.pack_start(gtk.Label("Min:"),False,False,2)            
        param.show_all()
        box.pack_start(param)

        box.show_all()
        params.pack_start(box)

        self.boxtable.attach(params,2,3,row,row+1)

        return True

    def connectIsolationCut(self,row,sel):
        label    = gtk.Label("Frixione isolation")
        params   = self.buildBox(label,sel,row)
        param    = gtk.HBox(False,0)
        param.show()
        params.pack_start(param,False,False,2)

        ranges   = sel.getRanges()

        epsi     = ranges["epsilon"]
        print "adjustment for epsiislon: ",epsi[2]," in [",epsi[0],epsi[1],"]"
        sel.setValue("epsilon",epsi[2])
        epsiadj  = gtk.Adjustment(epsi[2],epsi[0],epsi[1],epsi[3])
        epsisize = gtk.SpinButton(epsiadj,1,2)
        epsisize.show()
        epsisize.connect("value_changed",self.paramChanged,[sel,"epsilon"])
        param.pack_end(epsisize,False,False,2)
        param.pack_end(gtk.Label("epsilon:"),False,False,2)            
        
        expo     = ranges["expo"]
        print "adjustment for expo: ",expo[2]," in [",expo[0],expo[1],"]"
        sel.setValue("expo",expo[2])
        expoadj  = gtk.Adjustment(expo[2],expo[0],expo[1],expo[3])
        exposize = gtk.SpinButton(expoadj,1,2)
        exposize.show()
        exposize.connect("value_changed",self.paramChanged,[sel,"expo"])
        param.pack_end(exposize,False,False,2)
        param.pack_end(gtk.Label("exponent:"),False,False,2)            

        R        = ranges["R"]
        print "adjustment for R: ",R[2]," in [",R[0],R[1],"]"
        sel.setValue("R",R[2])
        Radj     = gtk.Adjustment(R[2],R[0],R[1],R[3])
        Rsize    = gtk.SpinButton(Radj,1,2)
        Rsize.show()
        Rsize.connect("value_changed",self.paramChanged,[sel,"R"])
        param.pack_end(Rsize,False,False,2)
        param.pack_end(gtk.Label("Cone size R:"),False,False,2)            

        param.show_all()
        self.boxtable.attach(params,2,3,row,row+1)
        return True


    def connectJetFinder(self,row,sel):
        label    = gtk.Label("Jet finder")

        params   = self.buildBox(label,sel,row)

        param    = gtk.HBox(False,0)
        param.show()
        params.pack_start(param,False,False,2)
        param.pack_start(gtk.Label("Algorithm:"),False,False,2)
        algo     = None
        algo     = gtk.RadioButton(algo,"SISCone")
        algo.connect("toggled",self.selectMode,[sel,"SISCone"])
        algo.show()
        param.pack_end(algo,False,False,2)
        algo     = gtk.RadioButton(algo,"Anti-KT")
        algo.connect("toggled",self.selectMode,[sel,"Anti-KT"])
        algo.show()
        param.pack_end(algo,False,False,2)
        algo     = gtk.RadioButton(algo,"KT")
        algo.connect("toggled",self.selectMode,[sel,"KT"])
        algo.show()
        param.pack_end(algo,False,False,2)

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
