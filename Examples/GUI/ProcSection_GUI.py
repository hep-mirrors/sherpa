#!/usr/bin/python

import pygtk
import gtk

import GUI_BaseClass as guibase
import ParameterBox
import ProcessBox

class Procsection_gui(guibase.gui_object):
    def initialise(self):
        self.procbox   = self.parameters.getProcBox()
        self.initTV(self.procbox.getOptions())
        self.initButtons()
        self.muF  = 1.
        self.muR  = 1.
        self.muQ  = 1.
        self.ckkw = 0.

    def getContent(self):
        table     = gtk.Table(1,2,False)
        table.attach(self.getProcSel(),0,1,0,1,gtk.FILL,gtk.FILL,0,0)
        table.attach(self.getOptSel(),1,2,0,1,gtk.FILL,gtk.FILL,0,0)
        table.set_col_spacing(0,10)
        table.show()
        return table

    def updateOptions(self):
        self.proc_view.set_model(self.procbox.getOptions())

    def procSelected(self,treesel):
        model,iterator = treesel.get_selected()
        if iterator==None:
            return
        tagname = model.get_value(iterator,0)
        print "   yields: ",tagname
        self.procbox.setProcess(tagname)
        LOjets,NLOjets = self.procbox.getNJets()
        self.minjetsbutton.get_adjustment().set_upper(LOjets)
        self.totjetsbutton.get_adjustment().set_upper(LOjets)
        self.nlojetsbutton.get_adjustment().set_upper(NLOjets)
        if self.procbox.getType()=="LC":
            self.ckkwbutton.get_adjustment().set_value(5.)
            self.ckkwbutton.get_adjustment().set_lower(2.)
            self.ckkwbutton.get_adjustment().set_upper(100.)
            self.ckkwbutton.get_adjustment().set_step_increment(2.)
        if self.procbox.getType()=="HC":
            self.ckkwbutton.get_adjustment().set_value(20.)
            self.ckkwbutton.get_adjustment().set_lower(10.)
            self.ckkwbutton.get_adjustment().set_upper(200.)
            self.ckkwbutton.get_adjustment().set_step_increment(5.)
        if self.procbox.getType()=="None":
            self.ckkwbutton.get_adjustment().set_value(0.)
            self.ckkwbutton.get_adjustment().set_lower(0.)
            self.ckkwbutton.get_adjustment().set_upper(0.)
            self.ckkwbutton.get_adjustment().set_step_increment(0.)
            self.ckkwbutton.set_sensitive(False)
        self.update("Proc")

    def ckkwChanged(self,button):
        self.ckkw = button.get_value()

    def scaleChanged(self,button,data):
        print " %s changed to %0.2f" %(data[0],data[1])
        if (data[0]=="muF"):
            self.muF = data[1]
        if (data[0]=="muR"):
            self.muR = data[1]
        if (data[0]=="muQ"):
            self.muQ = data[1]

    def buttonChanged(self,globalupdate,mode):
        self.minjets = self.minjetsbutton.get_value_as_int()
        print "new lower limits for jets:",self.minjets
        self.totjetsbutton.get_adjustment().set_lower(self.minjets)
        if (self.totjetsbutton.get_adjustment().get_value()<self.minjets):
            self.totjetsbutton.get_adjustment().set_value(self.minjets)
        self.totjets = self.totjetsbutton.get_value_as_int()
        if (self.nlojetsbutton.get_adjustment().get_value()>self.totjets):
            self.nlojetsbutton.get_adjustment().set_value(self.totjets)
        self.totjets = self.totjetsbutton.get_value_as_int()
        self.minjets = self.minjetsbutton.get_value_as_int()
        self.nlojets = self.nlojetsbutton.get_value_as_int()
        self.procbox.setJetMultis(self.minjets,self.totjets,self.nlojets)
        print "New multis in proc box:",self.procbox.getNJets()
        if (self.nlojets>-1):
            self.muqbox.set_sensitive(True)
        else:
            self.muqbox.set_sensitive(False)
        if self.totjets>self.minjets:
            self.ckkwbutton.set_sensitive(True)
        else:
            self.ckkwbutton.set_sensitive(False)            
        self.update("Proc")

    def extractParameters(self):
        self.procbox.setCKKW(self.ckkw)
        self.procbox.setScaleFactors(self.muF,self.muR,self.muQ)

    def initTV(self,options):
        self.proc_view = gtk.TreeView(options)
        tvcol          = gtk.TreeViewColumn('Process class')
        cell           = gtk.CellRendererText()
        cell.set_property("width",300)
        tvcol.pack_start(cell,True)
        tvcol.add_attribute(cell,'text',0)
        tvcol.set_sort_column_id(0)
        self.proc_view.set_search_column(0)
        self.proc_view.append_column(tvcol)
        self.proc_sel  = self.proc_view.get_selection()
        self.proc_sel.set_mode(gtk.SELECTION_SINGLE)
        self.proc_sel.connect("changed", self.procSelected)

    def initButtons(self):
        self.minjetsbutton = gtk.SpinButton()
        self.minjetsbutton.set_adjustment(gtk.Adjustment(0,0,0,1,1))
        self.minjetsbutton.set_numeric(True)
        self.minjetsbutton.connect("changed",self.buttonChanged,"Min")

        self.totjetsbutton = gtk.SpinButton()
        self.totjetsbutton.set_adjustment(gtk.Adjustment(0,0,0,1,1))
        self.totjetsbutton.set_numeric(True)
        self.totjetsbutton.connect("changed",self.buttonChanged,"Tot")

        self.nlojetsbutton = gtk.SpinButton()
        self.nlojetsbutton.set_adjustment(gtk.Adjustment(-1,-1,-1,1,1))
        self.nlojetsbutton.set_numeric(True)
        self.nlojetsbutton.connect("changed",self.buttonChanged,"NLO")

        self.ckkwbutton = gtk.SpinButton()
        self.ckkwbutton.set_adjustment(gtk.Adjustment(0.,0.,0.,5.,1))
        self.ckkwbutton.connect("changed",self.ckkwChanged)
        self.ckkwbutton.set_sensitive(False)

        self.mufbox = gtk.HBox(False,10)
        self.mufbox.set_border_width(10) 
        self.mufbox.pack_start(gtk.Label("mu_F2:"),False,False,2)
        self.mufbutton = gtk.RadioButton(None,"0.25")
        self.mufbutton.connect("toggled", self.scaleChanged, ["muF",0.25])
        self.mufbutton.show()
        self.mufbox.pack_start(self.mufbutton,False,False,2)
        self.mufbutton = gtk.RadioButton(self.mufbutton,"1.00")
        self.mufbutton.connect("toggled", self.scaleChanged, ["muF",1.00])
        self.mufbutton.set_active(True)
        self.mufbutton.show()
        self.mufbox.pack_start(self.mufbutton,False,False,2)
        self.mufbutton = gtk.RadioButton(self.mufbutton,"4.00")
        self.mufbutton.connect("toggled", self.scaleChanged, ["muF",4.00])
        self.mufbutton.show()
        self.mufbox.pack_start(self.mufbutton,False,False,2)

        self.murbox = gtk.HBox(False,10)
        self.murbox.set_border_width(10) 
        self.murbox.pack_start(gtk.Label("mu_R2:"),False,False,2)
        self.murbutton = gtk.RadioButton(None,"0.25")
        self.murbutton.connect("toggled", self.scaleChanged, ["muR",0.25])
        self.murbutton.show()
        self.murbox.pack_start(self.murbutton,False,False,2)
        self.murbutton = gtk.RadioButton(self.murbutton,"1.00")
        self.murbutton.connect("toggled", self.scaleChanged, ["muR",1.00])
        self.murbutton.set_active(True)
        self.murbutton.show()
        self.murbox.pack_start(self.murbutton,False,False,2)
        self.murbutton = gtk.RadioButton(self.murbutton,"4.00")
        self.murbutton.connect("toggled", self.scaleChanged, ["muR",4.00])
        self.murbutton.show()
        self.murbox.pack_start(self.murbutton,False,False,2)

        self.muqbox = gtk.HBox(False,10)
        self.muqbox.set_border_width(10) 
        self.muqbox.pack_start(gtk.Label("mu_Q2:"),False,False,2)
        self.muqbutton = gtk.RadioButton(None,"0.50")
        self.muqbutton.connect("toggled", self.scaleChanged, ["muQ",0.50])
        self.muqbutton.show()
        self.muqbox.pack_start(self.muqbutton,False,False,2)
        self.muqbutton = gtk.RadioButton(self.muqbutton,"1.00")
        self.muqbutton.connect("toggled", self.scaleChanged, ["muQ",1.00])
        self.muqbutton.set_active(True)
        self.muqbutton.show()
        self.muqbox.pack_start(self.muqbutton,False,False,2)
        self.muqbutton = gtk.RadioButton(self.muqbutton,"2.00")
        self.muqbutton.connect("toggled", self.scaleChanged, ["muQ",2.00])
        self.muqbox.set_sensitive(False)
        self.muqbutton.show()
        self.muqbox.pack_start(self.muqbutton,False,False,2)

    def getProcSel(self):
        procbox   = gtk.VBox(False,0)
        procbox.set_border_width(10)
        procbox.show()
        box = gtk.HBox(False,0)
        procbox.pack_start(box,False,False,2)
        box.pack_start(self.proc_view)
        return procbox

    def getOptSel(self):
        buttonbox = gtk.VBox(False,0)
        buttonbox.set_border_width(10)
        buttonbox.show()

        jetbox = gtk.VBox(False,0)
        jetbox.set_border_width(20)
        jetbox.show()
        jetbox.pack_start(gtk.Label("Number of additional jets"))
        minjets = gtk.HBox(False,0)
        minjets.set_border_width(10)
        minjets.show()
        minjets.pack_start(gtk.Label("Minimal:"),False,False,2)
        minjets.pack_end(self.minjetsbutton,False,False,2)
        totjets = gtk.HBox(False,0)
        totjets.set_border_width(10)
        totjets.show()
        totjets.pack_start(gtk.Label("Total:"),False,False,2)
        totjets.pack_end(self.totjetsbutton,False,False,2)
        NLOjets = gtk.HBox(False,0)
        NLOjets.set_border_width(10)
        NLOjets.show()
        NLOjets.pack_start(gtk.Label("NLO:"),False,False,2)
        NLOjets.pack_end(self.nlojetsbutton,False,False,2)
        jetbox.pack_start(minjets,False,False,2)
        jetbox.pack_start(totjets,False,False,2)
        jetbox.pack_start(NLOjets,False,False,2)

        ckkwbox = gtk.VBox(False,0)
        ckkwbox.set_border_width(20)
        ckkwbox.show()
        ckkwbox.pack_start(gtk.Label("CKKW merging parameter"))
        ckkw = gtk.HBox(False,0)
        ckkw.set_border_width(10)
        ckkw.show()
        ckkw.pack_start(gtk.Label("PT [GeV]:"),False,False,2)
        ckkw.pack_end(self.ckkwbutton,False,False,2)
        ckkwbox.pack_start(ckkw,False,False,2)

        scalebox = gtk.VBox(False,0)
        scalebox.set_border_width(20)
        scalebox.show()
        scalebox.pack_start(gtk.Label("Quadratic scale factors"))
        scalebox.pack_start(self.mufbox,False,False,2)
        scalebox.pack_start(self.murbox,False,False,2)
        scalebox.pack_start(self.muqbox,False,False,2)

        buttonbox.pack_start(jetbox)
        buttonbox.pack_start(ckkwbox)
        buttonbox.pack_start(scalebox)
        return buttonbox
