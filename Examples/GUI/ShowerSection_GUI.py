#!/usr/bin/python

import pygtk
import gtk

import GUI_BaseClass as guibase
import ParameterBox
import ShowerBox

class Showersection_gui(guibase.gui_object):
    def initialise(self):
        self.showerbox  = self.parameters.getShowerBox()
        self.switchbox,self.parambox = self.initButtons()
        self.fillParamBox()
        pass

    def updateOptions(self,isNLO):
        print "####################### Shower updates options(",isNLO,")"
        self.showerbox.initialiseDefaults(isNLO)
        fs_pt2min,is_pt2min,fs_asfac,is_asfac = self.showerbox.getParams()
        kinmode = self.showerbox.getKinScheme()
        self.fsasfsel.get_adjustment().set_value(fs_asfac)
        self.isasfsel.get_adjustment().set_value(is_asfac)

    def getContent(self):
        table     = gtk.Table(2,2,False)
        table.set_row_spacing(0,30)
        table.attach(self.switchbox,0,1,1,2,gtk.FILL,gtk.FILL,0,0)
        table.attach(self.parambox,1,2,1,2,gtk.FILL,gtk.FILL,0,0)
        table.show()
        return table

    def switched(self,button,box):
        self.showerbox.switch()
        if self.showerbox.isOn():
            box.set_sensitive(True)
        else:
            box.set_sensitive(False)

    def paramChanged(self,button,data):
        if (data[0]=="FS_PT2Min"):
            self.showerbox.setFS_PT2Min(button.get_value())
        if (data[0]=="IS_PT2Min"):
            self.showerbox.setIS_PT2Min(button.get_value())
        if (data[0]=="FS_AsFac"):
            self.showerbox.setFS_AsFac(button.get_value())
        if (data[0]=="IS_AsFac"):
            self.showerbox.setIS_AsFac(button.get_value())
        if (data[0]=="KinScheme"):
            self.showerbox.setKinScheme(data[1])

    def fillParamBox(self):
        fs_pt2min,is_pt2min,fs_asfac,is_asfac = self.showerbox.getParams()
        kinmode = self.showerbox.getKinScheme()

        param    = gtk.HBox(True,0)
        param.set_size_request(540,30)
        label = gtk.Label("Final state:")
        label.set_size_request(150,30)
        param.pack_start(label,False,False,2)            
        label = gtk.Label("PT^2_Min = ")
        label.set_size_request(90,30)
        param.pack_start(label,False,False,2)            
        fspt2adj = gtk.Adjustment(fs_pt2min,0.5,4.0,0.1)
        self.fspt2sel = gtk.SpinButton(fspt2adj,1,2)
        self.fspt2sel.set_size_request(90,30)
        self.fspt2sel.show()
        self.fspt2sel.connect("value_changed",self.paramChanged,["FS_PT2Min",0])
        param.pack_start(self.fspt2sel,False,False,2)
        label = gtk.Label("As_Fac = ")
        label.set_size_request(90,30)
        param.pack_start(label,False,False,2)            
        fsasfadj = gtk.Adjustment(fs_asfac,0.5,4.0,0.1)
        self.fsasfsel = gtk.SpinButton(fsasfadj,1,2)
        self.fsasfsel.set_size_request(90,30)
        self.fsasfsel.show()
        self.fsasfsel.connect("value_changed",self.paramChanged,["FS_AsFac",0])
        param.pack_start(self.fsasfsel,False,False,2)
        param.show_all()
        self.parambox.pack_start(param)

        param    = gtk.HBox(True,0)
        param.set_size_request(540,30)
        label = gtk.Label("Initial state:")
        label.set_size_request(150,30)
        param.pack_start(label,False,False,2)            
        label = gtk.Label("PT^2_Min = ")
        label.set_size_request(90,30)
        param.pack_start(label,False,False,2)            
        ispt2adj = gtk.Adjustment(is_pt2min,0.5,4.0,0.1)
        self.ispt2sel = gtk.SpinButton(ispt2adj,1,2)
        self.ispt2sel.set_size_request(90,30)
        self.ispt2sel.show()
        self.ispt2sel.connect("value_changed",self.paramChanged,["IS_PT2Min",0])
        param.pack_start(self.ispt2sel,False,False,2)
        label = gtk.Label("As_Fac = ")
        label.set_size_request(90,30)
        param.pack_start(label,False,False,2)            
        isasfadj = gtk.Adjustment(is_asfac,0.5,4.0,0.1)
        self.isasfsel = gtk.SpinButton(isasfadj,1,2)
        self.isasfsel.set_size_request(90,30)
        self.isasfsel.show()
        self.isasfsel.connect("value_changed",self.paramChanged,["IS_AsFac",0])
        param.pack_start(self.isasfsel,False,False,2)
        param.show_all()
        self.parambox.pack_start(param)

        param    = gtk.HBox(True,0)
        param.set_size_request(540,30)
        label = gtk.Label("Kinematics:")
        label.set_size_request(120,30)
        param.pack_start(label,False,False,2)
        self.modebutton = gtk.RadioButton(None,"0")
        self.modebutton.set_size_request(40,30)
        self.modebutton.connect("toggled", self.paramChanged, 
                                ["KinScheme",0])
        if kinmode==0:
            self.modebutton.set_active(True)
        self.modebutton.show()
        param.pack_start(self.modebutton,False,False,2)
        self.modebutton = gtk.RadioButton(self.modebutton,"1")
        self.modebutton.set_size_request(40,30)
        self.modebutton.connect("toggled", self.paramChanged, 
                                ["KinScheme",1])
        if kinmode==1:
            self.modebutton.set_active(True)
        self.modebutton.show()
        param.pack_start(self.modebutton,False,False,2)
        label = gtk.Label("")
        label.set_size_request(350,30)
        param.pack_start(label,False,False,2)
        param.show_all()
        self.parambox.pack_start(param)

    def initButtons(self):
        box    = gtk.VBox(False,2)

        switch = gtk.Button()
        switch.set_size_request(40,30)
        switch.show()
        sbox   = gtk.HBox(False,2)
        sbox.pack_start(switch,False,False,2)
        sbox.set_size_request(40,30)
        sbox.show()
        if self.showerbox.isOn():
            switch.set_label("On")
            box.set_sensitive(True)
        else:
            switch.set_label("Off")
            box.set_sensitive(False)
        switch.connect("clicked",self.switched,box)
        box.show()
        return sbox,box

