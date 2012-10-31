#!/usr/bin/python

import pygtk
import gtk

import GUI_BaseClass as guibase
import ParameterBox
import BeamBox


class beamsection_gui(guibase.gui_object):
    def initialise(self):
        self.beams = self.parameters.getBeamBox()
        self.part1 = None
        self.part2 = None
        self.adj1  = gtk.Adjustment(4000.,0.,1000000.,1000.)
        self.adj2  = gtk.Adjustment(4000.,0.,1000000.,1000.)
        self.ene1  = None
        self.ene2  = None

    def getContent(self):
        colbox    = gtk.VBox(False,0)
        colbox.show()
        colbox.pack_start(self.makeColBox(),False,False,10)
        colbox.pack_start(self.makeUserBox(),False,False,10)
        return colbox

    def makeColBox(self):
        optbox    = gtk.HBox(False,spacing=10)
        optbox.set_border_width(10)
        optbox.show()
        label  = gtk.Label()
        label.set_text('Collider')
        opt    = gtk.OptionMenu()
        menu   = gtk.Menu()
        for col in self.beams.getColliders():
            print col[0],"|",col[1],"|"
            item = self.make_menu_item(col[1],self.col_select,col[0])
            menu.append(item)
        opt.set_menu(menu)
        optbox.pack_start(label,False,False,2)
        optbox.pack_start(opt,False,False,10)
        return optbox

    def makeUserBox(self):
        userbox    = gtk.HBox(False,0)
        userbox.set_border_width(10)
        beambox1   = gtk.HBox(False,0)
        beambox2   = gtk.HBox(False,0)
        label1     = gtk.Label("Beam 1")
        label2     = gtk.Label("Beam 2")
        opt1       = gtk.OptionMenu()
        self.part1 = gtk.Menu()
        for part in self.beams.getBeamParticles():
            print part[0],part[1]
            item = self.make_menu_item(part[1],self.part1_select,part[0])
            self.part1.append(item)
        self.part1.set_sensitive(False)
        opt1.set_menu(self.part1)
        opt2       = gtk.OptionMenu()
        self.part2 = gtk.Menu()
        for part in self.beams.getBeamParticles():
            print part[0],part[1]
            item = self.make_menu_item(part[1],self.part2_select,part[0])
            self.part2.append(item)
        self.part2.set_sensitive(False)
        opt2.set_menu(self.part2)
        self.ene1  = gtk.SpinButton(self.adj1,100.,2)
        self.ene1.set_sensitive(False)
        self.ene2  = gtk.SpinButton(self.adj2,100.,2)
        self.ene2.set_sensitive(False)
        beambox1.pack_start(label1,False,False,2)
        beambox1.pack_start(opt1,False,False,2)
        beambox1.pack_start(self.ene1,False,False,2)
        beambox2.pack_start(label2,False,False,2)
        beambox2.pack_start(opt2,False,False,2)
        beambox2.pack_start(self.ene2,False,False,2)
        userbox.pack_start(beambox1,False,False,2)
        userbox.pack_end(beambox2,False,False,2)
        return userbox

    def col_select(self,item,tag):
        print "tag = ",tag
        if tag==0:
            print "No collider selected."
        if tag==-1:
            print "User defined collider!"
            self.part1.set_sensitive(True)
            self.part2.set_sensitive(True)
            self.ene1.set_sensitive(True)
            self.ene2.set_sensitive(True)
            self.parameters.updateBeams(tag)
            self.beams.updateCollider(tag)
            self.update("Beams")
        else:
            print "Pre-defined collider!"
            self.part1.set_sensitive(False)
            self.part2.set_sensitive(False)
            self.ene1.set_sensitive(False)
            self.ene2.set_sensitive(False)
            self.parameters.updateBeams(tag)
            self.beams.printCollider()
            self.update("Beams")


    def part1_select(self,item,tag):
        print "tag = ",tag
        if tag==11 or tag==-11:
            print "Select lepton for beam 1"
            self.ene1.get_adjustment().set_lower(0.)
            self.ene1.get_adjustment().set_upper(100000.)
            self.ene1.get_adjustment().set_value(45.6)
            self.ene1.get_adjustment().set_step_increment(10.)
        if tag==2212 or tag==-2212:
            print "Select hadron for beam 1"
            self.ene1.get_adjustment().set_lower(0.)
            self.ene1.get_adjustment().set_upper(100000.)
            self.ene1.get_adjustment().set_value(4000.)
            self.ene1.get_adjustment().set_step_increment(1000.)
        self.beams.setBeam(1,tag,self.ene1.get_value())
        self.beams.printCollider()
        self.update("Beams")


    def part2_select(self,item,tag):
        print "tag = ",tag
        if tag==11 or tag==-11:
            print "Select lepton for beam 1"
            self.ene2.get_adjustment().set_lower(0.)
            self.ene2.get_adjustment().set_upper(100000.)
            self.ene2.get_adjustment().set_value(45.6)
            self.ene2.get_adjustment().set_step_increment(10.)
        if tag==2212 or tag==-2212:
            print "Select hadron for beam 1"
            self.ene2.get_adjustment().set_lower(0.)
            self.ene2.get_adjustment().set_upper(100000.)
            self.ene2.get_adjustment().set_value(4000.)
            self.ene2.get_adjustment().set_step_increment(1000.)
        self.beams.setBeam(2,tag,self.ene2.get_value())
        self.beams.printCollider()
        self.update("Beams")

    def extractValues(self):
        pass
