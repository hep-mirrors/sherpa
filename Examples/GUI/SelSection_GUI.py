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
        hbox          = gtk.HBox(True,0)
        hbox.set_border_width(10)
        vbox.pack_start(hbox,False,False,2)
        hbox          = gtk.HBox(True,0)
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
            print "Deal with ",sel.getTag(),row
            flag1 = self.addSelector(row,sel)
            flag  = flag or flag1
            row   = row+1
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
        label.set_size_request(120,30)
        label.show()
        hbox   = gtk.HBox(False,2)
        hbox.pack_start(label,False,False,2)
        hbox.set_size_request(120,30)
        hbox.show()

        switch = gtk.Button()
        switch.set_size_request(40,30)
        switch.show()
        sbox   = gtk.HBox(False,2)
        sbox.pack_start(switch,False,False,2)
        sbox.set_size_request(40,30)
        sbox.show()
        if sel.isOn():
            switch.set_label("On")
            box.set_sensitive(True)
        else:
            switch.set_label("Off")
            box.set_sensitive(False)
        switch.connect("clicked",self.switchedOn,sel,box)
        self.boxtable.attach(hbox,0,1,row,row+1)
        self.boxtable.attach(sbox,1,2,row,row+1)
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
        params.set_size_request(400,30)
        ranges   = sel.getRanges()
        mini     = ranges["min"]
        maxi     = ranges["max"]

        box      = gtk.HBox(True,0)
        box.set_size_request(400,30)

        param    = gtk.HBox(True,0)
        param.set_size_request(130,30)
        label = gtk.Label("Max:")
        label.set_size_request(40,30)
        param.pack_start(label,False,False,2)            
        sel.setValue("max",maxi[2])
        maxadj   = gtk.Adjustment(maxi[2],maxi[0],maxi[1],maxi[3])
        maxsel   = gtk.SpinButton(maxadj,1,2)
        maxsel.set_size_request(90,30)
        maxsel.show()
        maxsel.connect("value_changed",self.paramChanged,[sel,"max"])
        param.pack_start(maxsel,False,False,2)
        param.show_all()
        box.pack_end(param)

        param    = gtk.HBox(True,0)
        param.set_size_request(130,30)
        label = gtk.Label("Min:")
        label.set_size_request(40,30)
        param.pack_start(label,False,False,2)            
        sel.setValue("min",mini[2])
        minadj   = gtk.Adjustment(mini[2],mini[0],mini[1],mini[3])
        minsel   = gtk.SpinButton(minadj,1,2)
        minsel.set_size_request(90,30)
        minsel.show()
        minsel.connect("value_changed",self.paramChanged,[sel,"min"])
        param.pack_start(minsel,False,False,2)
        param.show_all()
        box.pack_start(param)

        box.show_all()
        params.pack_start(box)

        self.boxtable.attach(params,2,3,row,row+1)

        return True

    def connectIsolationCut(self,row,sel):
        label    = gtk.Label("Frixione iso")
        params   = self.buildBox(label,sel,row)
        param    = gtk.HBox(True,0)
        param.set_size_request(400,30)
        param.show()
        params.pack_start(param,False,False,2)

        ranges   = sel.getRanges()
        
        
        Rbox   = gtk.HBox(True,0)
        Rbox.set_size_request(120,30)
        label = gtk.Label("R:")
        label.set_size_request(40,30)
        Rbox.pack_start(label,False,False,2)            
        R        = ranges["R"]
        print "adjustment for R: ",R[2]," in [",R[0],R[1],"]"
        sel.setValue("R",R[2])
        Radj     = gtk.Adjustment(R[2],R[0],R[1],R[3])
        Rsize    = gtk.SpinButton(Radj,0,2)
        Rsize.set_size_request(90,30)
        Rsize.show()
        Rsize.connect("value_changed",self.paramChanged,[sel,"R"])
        Rbox.pack_start(Rsize,False,False,2)
        param.pack_end(Rbox,False,False,2)            

        expbox   = gtk.HBox(True,0)
        expbox.set_size_request(120,30)
        label = gtk.Label("expo:")
        label.set_size_request(40,30)
        expbox.pack_start(label,False,False,2)            
        expo     = ranges["expo"]
        print "adjustment for expo: ",expo[2]," in [",expo[0],expo[1],"]"
        sel.setValue("expo",expo[2])
        expoadj  = gtk.Adjustment(expo[2],expo[0],expo[1],expo[3])
        exposize = gtk.SpinButton(expoadj,0,2)
        exposize.set_size_request(80,30)
        exposize.show()
        exposize.connect("value_changed",self.paramChanged,[sel,"expo"])
        expbox.pack_start(exposize,False,False,2)
        param.pack_start(expbox,False,False,2)            


        epsbox   = gtk.HBox(True,0)
        epsbox.set_size_request(120,30)
        label = gtk.Label("eps:")
        label.set_size_request(40,30)
        epsbox.pack_start(label,False,False,2)            
        epsi     = ranges["epsilon"]
        print "adjustment for epsiislon: ",epsi[2]," in [",epsi[0],epsi[1],"]"
        sel.setValue("epsilon",epsi[2])
        epsiadj  = gtk.Adjustment(epsi[2],epsi[0],epsi[1],epsi[3])
        epsisize = gtk.SpinButton(epsiadj,0,2)
        epsisize.set_size_request(80,30)
        epsisize.show()
        epsisize.connect("value_changed",self.paramChanged,[sel,"epsilon"])
        epsbox.pack_start(epsisize,False,False,2)
        param.pack_start(epsbox,False,False,2)            

        param.show_all()
        self.boxtable.attach(params,2,3,row,row+1)
        return True


    def connectJetFinder(self,row,sel):
        label    = gtk.Label("Jet finder")

        params   = self.buildBox(label,sel,row)
        params.set_size_request(400,60)
        self.boxtable.attach(params,2,3,row,row+2)
        row = row+1

        param    = gtk.HBox(True,0)
        param.set_size_request(400,30)
        param.show()
        params.pack_start(param,False,False,2)
        label = gtk.Label("Algorithm:")
        label.set_size_request(80,30)
        label.show()
        sel.setValue("Mode","Anti-KT")
        param.pack_start(label,False,False,2)
        algo     = None
        algo     = gtk.RadioButton(algo,"Anti-KT")
        algo.set_size_request(80,30)
        algo.connect("toggled",self.selectMode,[sel,"Anti-KT"])
        algo.show()
        param.pack_end(algo,False,False,2)
        algo     = gtk.RadioButton(algo,"Cam-A")
        algo.set_size_request(80,30)
        algo.connect("toggled",self.selectMode,[sel,"Cam-A"])
        algo.show()
        param.pack_end(algo,False,False,2)
        algo     = gtk.RadioButton(algo,"KT")
        algo.set_size_request(50,30)
        algo.connect("toggled",self.selectMode,[sel,"KT"])
        algo.show()
        param.pack_end(algo,False,False,2)


        param   = gtk.HBox(True,0)
        param.set_size_request(400,30)
        param.show()
        params.pack_start(param,False,False,2)
        ranges   = sel.getRanges()
        PT       = ranges["PT"]
        R        = ranges["R"]

        Rbox   = gtk.HBox(True,0)
        Rbox.set_size_request(130,30)
        Rbox.show()
        label = gtk.Label("R:")
        label.set_size_request(40,30)
        Rbox.pack_start(label,False,False,2)            
        print "adjustment for R: ",R[2]," in [",R[0],R[1],"]"
        sel.setValue("R",R[2])
        Radj     = gtk.Adjustment(R[2],R[0],R[1],R[3])
        Rsize    = gtk.SpinButton(Radj,1,2)
        Rsize.set_size_request(90,30)
        Rsize.show()
        Rsize.connect("value_changed",self.paramChanged,[sel,"R"])
        Rbox.pack_start(Rsize,False,False,2)
        param.pack_end(Rbox,False,False,2)

        ptbox   = gtk.HBox(True,0)
        ptbox.set_size_request(130,30)
        ptbox.show()
        label = gtk.Label("PT:")
        label.set_size_request(40,30)
        ptbox.pack_start(label,False,False,2)            
        sel.setValue("PT",PT[2])
        ptadj    = gtk.Adjustment(PT[2],PT[0],PT[1],PT[3])
        ptsize   = gtk.SpinButton(ptadj,1,0)
        ptsize.set_size_request(90,30)
        ptsize.show()
        ptsize.connect("value_changed",self.paramChanged,[sel,"PT"])
        ptbox.pack_start(ptsize,False,False,2)
        param.pack_start(ptbox,False,False,10)
        param.show_all()
        return True
