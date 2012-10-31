#!/usr/bin/python

from gi.repository import Gtk
import os

import ParameterBox
import ProcessBox

class Procsection_gui_gtk():
    def __init__(self,ParameterBox):
        self.parameters       = ParameterBox
        self.procbox          = self.parameters.getProcBox()
        self.proc_sel         = None
        self.proc_jets        = 0
        self.proc_jets_at_NLO = 0
        self.minjets          = 0
        self.totjets          = 0
        self.nlojets          = 0
        self.ckkwbutton       = None
        self.mufbutton        = None
        self.murbutton        = None
        self.muqbutton        = None
        pass

    def selection(self,grid,globalupdate):
        renderer_text  = Gtk.CellRendererText()
        sel_field      = Gtk.Box(spacing=10) 
        label          = Gtk.Label()
        label.set_text('Process type')
        self.proc_sel  = Gtk.ComboBox.new_with_model(self.procbox.getOptions())
        self.proc_sel.connect("changed", self.proc_changed,globalupdate)
        self.proc_sel.pack_start(renderer_text, True)
        self.proc_sel.add_attribute(renderer_text, "text", 0)
        self.proc_sel.add_attribute(renderer_text, "sensitive", 1)
        sel_field.pack_start(label,False,False,2)
        sel_field.pack_start(self.proc_sel,False,False,20)


        njet_field     = Gtk.Box(spacing=10) 
        labelMinJets   = Gtk.Label()
        labelMinJets.set_text('Jet multiplicities: min njets ')
        labelMaxJets   = Gtk.Label()
        labelMaxJets.set_text('max njets ')
        labelNLOJets   = Gtk.Label()
        labelNLOJets.set_text('out of which NLO ')
        LOmax, NLOmax      = self.procbox.getNJets() 
        minjetlimits       = Gtk.Adjustment(0, 0, LOmax, 1, 1,-1)
        self.minjetsbutton = Gtk.SpinButton()
        self.minjetsbutton.set_adjustment(minjetlimits)
        self.minjetsbutton.set_numeric(True)
        totjetlimits       = Gtk.Adjustment(0, self.minjets, LOmax, 1, 1,-1)
        self.totjetsbutton = Gtk.SpinButton()
        self.totjetsbutton.set_adjustment(totjetlimits)
        self.totjetsbutton.set_numeric(True)
        NLOjetlimits       = Gtk.Adjustment(-1, -1, NLOmax, 1, 1,-1)
        self.nlojetsbutton = Gtk.SpinButton()
        self.nlojetsbutton.set_adjustment(NLOjetlimits)
        self.nlojetsbutton.set_numeric(True)
        self.totjetsbutton.connect("changed",self.button_changed)
        self.nlojetsbutton.connect("changed",self.button_changed)
        self.minjetsbutton.connect("changed",self.button_changed)
        njet_field.pack_start(labelMinJets,False,False,2)
        njet_field.pack_start(self.minjetsbutton,False,False,6)
        njet_field.pack_start(labelMaxJets,False,False,2)
        njet_field.pack_start(self.totjetsbutton,False,False,6)
        njet_field.pack_start(labelNLOJets,False,False,2)
        njet_field.pack_start(self.nlojetsbutton,False,False,6)
        
        ckkw_field      = Gtk.Box(spacing=10) 
        ckkw_label      = Gtk.Label()
        ckkw_label.set_text("CKKW merging value: (PT)")
        self.ckkwbutton = Gtk.SpinButton()
        ckkwlimits      = Gtk.Adjustment(20., 10., 500., 5., 1,-1)
        ckkw_field.pack_start(ckkw_label,False,False,2)
        ckkw_field.pack_start(self.ckkwbutton,False,False,10)

        scale_field      = Gtk.Box(spacing=10) 
        scf_label        = Gtk.Label()
        scf_label.set_text('Scale factor selection: mu_R ')
        scr_label        = Gtk.Label()
        scr_label.set_text('mu_F')
        scq_label        = Gtk.Label()
        scq_label.set_text('mu_Q')
        self.mufbutton = Gtk.RadioButton(None,"1")
        self.mufbutton.connect("toggled",self.callmuf,"1")
        self.mufbutton = Gtk.RadioButton(self.mufbutton,"0.5")
        self.mufbutton = Gtk.RadioButton(self.mufbutton,"2")
        #ckkw_field.pack_start(scf_label,False,False,2)
        #ckkw_field.pack_start(self.mufbutton,False,False,10)
        #ckkw_field.pack_start(scr_label,False,False,2)
        #ckkw_field.pack_start(self.murbutton,False,False,10)
        #ckkw_field.pack_start(scq_label,False,False,2)
        #ckkw_field.pack_start(self.muqbutton,False,False,10)


        grid.attach(sel_field,0,3,2,1)
        grid.attach(njet_field,0,4,3,1)
        #grid.attach(ckkw_field,0,5,2,1)


    def proc_changed(self,combo,globalupdate):
        print "In proc_changed"
        tree_iter = combo.get_active_iter()
        if tree_iter != None:
            model    = combo.get_model()
            proc_tag = model[tree_iter][0]
            print "Tag = ",proc_tag
            self.procbox.setProcess(proc_tag)
            LOjets,NLOjets = self.procbox.getNJets()
            self.minjetsbutton.get_adjustment().set_upper(LOjets)
            self.totjetsbutton.get_adjustment().set_upper(LOjets)
            self.nlojetsbutton.get_adjustment().set_upper(NLOjets)
            globalupdate("Procs")

    def button_changed(self,globalupdate):
        self.minjets = self.minjetsbutton.get_value_as_int()
        print "new lower limits for jets:",self.minjets
        self.totjetsbutton.get_adjustment().set_lower(self.minjets)
        #self.nlojetsbutton.get_adjustment().set_lower(self.minjets)
        if (self.totjetsbutton.get_adjustment().get_value()<self.minjets):
            self.totjetsbutton.get_adjustment().set_value(self.minjets)
        #if (self.nlojetsbutton.get_adjustment().get_value()<self.minjets):
        #    self.nlojetsbutton.get_adjustment().set_value(self.minjets)
        self.totjets = self.totjetsbutton.get_value_as_int()
        if (self.nlojetsbutton.get_adjustment().get_value()>self.totjets):
            self.nlojetsbutton.get_adjustment().set_value(self.totjets)
        self.nlojets = self.nlojetsbutton.get_value_as_int()
        self.procbox.setJetMultis(self.minjets,self.totjets,self.nlojets)

    def callmuf(slef):
        pass

    def updateOptions(self):
        pass

