#!/usr/bin/python

import pygtk
import gtk

import GUI_BaseClass as guibase
import ParameterBox
import BeamBox


class beamsection_gui(guibase.gui_object):
    def initialise(self):
        self.beams = self.parameters.getBeamBox()
        self.part  = [gtk.Menu(),gtk.Menu()]
        self.part[0].set_sensitive(False)
        self.part[1].set_sensitive(False)
        self.adj   = [gtk.Adjustment(0.,0.,1000000.,1000.),
                      gtk.Adjustment(0.,0.,1000000.,1000.)]
        self.ene   = [None,None]

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
            menu.append(self.make_menu_item(col[1],self.selectCollider,col[0]))
        opt.set_menu(menu)
        optbox.pack_start(label,False,False,2)
        optbox.pack_start(opt,False,False,10)
        return optbox

    def makeUserBox(self):
        userbox    = gtk.HBox(False,0)
        userbox.set_border_width(10)
        beambox    = [gtk.HBox(False,0),gtk.HBox(False,0)]
        label      = [gtk.Label("Beam 1"),gtk.Label("Beam 2")]
        opt        = [gtk.OptionMenu(),gtk.OptionMenu()]
        for beam in range(0,2):
            for part in self.beams.getBeamParticles():
                print "Beam (",beam,") add: ",part[0],part[1]
                item = self.make_menu_item(part[1],self.selectPart,
                                           [part[0],beam])
                self.part[beam].append(item)
            opt[beam].set_menu(self.part[beam])
            self.ene[beam] = gtk.SpinButton(self.adj[beam],100.,2)
            self.ene[beam].set_sensitive(False)
            beambox[beam].pack_start(label[beam],False,False,2)
            beambox[beam].pack_start(opt[beam],False,False,2)
            beambox[beam].pack_start(self.ene[beam],False,False,2)
            userbox.pack_start(beambox[beam],False,False,2)
        return userbox

    def selectCollider(self,item,tag):
        print "tag = ",tag
        if tag==0:
            print "No collider selected."
        if tag==-1:
            print "User defined collider!"
            flag = True
        else:
            print "Pre-defined collider!"
            flag = False
        for beam in range(0,2):
            self.part[beam].set_sensitive(flag)
            self.ene[beam].set_sensitive(flag)
        self.beams.setCollider(tag)
        self.beams.printCollider()
        self.update("Beams")


    def selectPart(self,item,tags):
        print "BeamGUI::selectPart(tag = ",tags[0],tags[1],")"
        if tags[0]==11 or tags[0]==-11:
            print "Select lepton for beam 1"
            self.ene[tags[1]].get_adjustment().set_lower(0.)
            self.ene[tags[1]].get_adjustment().set_upper(100000.)
            self.ene[tags[1]].get_adjustment().set_value(45.6)
            self.ene[tags[1]].get_adjustment().set_step_increment(10.)
        if tags[0]==2212 or tags[0]==-2212:
            print "Select hadron for beam 1"
            self.ene[tags[1]].get_adjustment().set_lower(0.)
            self.ene[tags[1]].get_adjustment().set_upper(100000.)
            self.ene[tags[1]].get_adjustment().set_value(4000.)
            self.ene[tags[1]].get_adjustment().set_step_increment(1000.)
        self.beams.setBeam(tags[1],tags[0],self.ene[tags[1]].get_value())
        self.beams.printCollider()
        self.update("Beams")

    def extractParameters(self):
        print "Check this:"
        self.beams.printCollider()
        pass
