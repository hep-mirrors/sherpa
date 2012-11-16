#!/usr/bin/python

import pygtk
import gtk

import GUI_BaseClass as guibase
import ParameterBox
import BeamBox


class generalsettings_gui(guibase.gui_object):
    def initialise(self):
        self.params = self.parameters


    def getContent(self):
        vbox    = gtk.VBox(False,0)
        vbox.set_border_width(10)
        vbox.show()
        vbox.pack_start(self.makeEvtBox(),False,False,10)
        #vbox.pack_start(self.makeOutputBox(),False,False,10)
        return vbox

    def paramChanged(self,button,mode):
        if mode=="Events":
            self.params.setNEvents(button.get_value_as_int())
            return

    def selectMode(self,button,mode):
        if mode=="Weighted" or mode=="Unweighted":
            self.params.setEvtMode(mode)
            return

    def makeEvtBox(self):
        box     = gtk.VBox(False,0)
        box.set_border_width(10)
        ebox    = gtk.HBox(False,0)
        ebox.set_border_width(10)
        ebox.pack_start(gtk.Label("Number of events:"),False,False,10)
        button  = gtk.SpinButton(gtk.Adjustment(self.params.getNEvents(),
                                                0,1000000000,10000),1,0)
        button.connect("value_changed",self.paramChanged,"Events")
        ebox.pack_start(button,False,False,10)
        ebox.show_all()
        box.pack_start(ebox)

        mbox    = gtk.HBox(False,0)
        mbox.set_border_width(10)
        mbox.pack_start(gtk.Label("Event mode:"),False,False,10)
        button  = gtk.RadioButton(None,"Weighted")
        button.connect("toggled",self.selectMode,"Weighted")
        button.show()
        mbox.pack_start(button,False,False,10)
        button  = gtk.RadioButton(button,"Unweighted")
        button.connect("toggled",self.selectMode,"Unweighted")
        button.show()
        mbox.pack_start(button,False,False,10)
        mbox.show_all()
        box.pack_start(mbox)
        return box

    def makeOutputBox(self):
        pass
