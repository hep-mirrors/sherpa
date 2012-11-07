#!/usr/bin/python

import pygtk
import gtk

import GUI_BaseClass as guibase
import ParameterBox
import BeamBox


class groupsection(guibase.gui_object):
    def __init__(self):
        print "Init group bit."
        self.type       = "gui_object" 
        self.id         = 666
        self.fname      = "About"
        self.flabel     = "About"
        self.update     = None

    def getContent(self):
        #vbox = gtk.VBox(False,0)
        #table = gtk.Table(2,2)
        #table.show()
        #vbox.pack_start(table)
        #entry = gtk.Entry()
        entry = gtk.VBox(False,0)
        self.addText(entry)
        tbox = gtk.ScrolledWindow()
        tbox.set_policy(gtk.POLICY_NEVER, gtk.POLICY_AUTOMATIC)
        tbox.add(entry)
        return tbox

    def getLabel(self):
        pixbuf = gtk.gdk.pixbuf_new_from_file_at_size("SherpaLogo.jpg",128,128)
        #pixbuf.scale_simple(10,10,gtk.gdk.INTERP_BILINEAR)
        image = gtk.Image()
        image.set_from_pixbuf(pixbuf)
        return image

    def getFrame(self,borderwidth=10,xsize=200,ysize=200,update=None):
        frame = gtk.Frame(self.flabel)
        frame.set_border_width(borderwidth)
        frame.set_size_request(xsize,ysize)
        frame.add(self.getContent())
        self.update = update
        return frame

    def addText(self,box):
        entry = gtk.TextView()
        entry.set_wrap_mode(gtk.WRAP_WORD)
        entry.set_justification(gtk.JUSTIFY_LEFT)
        entry.set_pixels_above_lines(8)
        entry.set_left_margin(40)
        entry.set_right_margin(40)
        textbuffer = entry.get_buffer()
        indent = "                     "
        ts =      "\nSHERPA 2.0.0\n\n"
        ts = ts + "Authors:    "
        ts = ts + "Stefan Hoeche,  Hendrik Hoeth,  Frank Krauss," 
        ts = ts + "Marek Schoenherr, \n"
        ts = ts + indent+"  Steffen Schumann,  Frank Siegert,  "
        ts = ts + "Jennifer Thompson,\n"
        ts = ts + indent+"  Jan Winter,  Korinna Zapp.\n\n"
        ts = ts + "Former Authors:    "
        ts = ts + "Timo Fischer,  Tanju Gleisberg,  Ralf Kuhn, \n"              
        ts = ts + indent+indent+"Thomas Laubrich,  Andreas Schaelicke.\n\n"

        ts = ts + "This program uses a lot of genuine and original research "
        ts = ts + "work by other people. Users are encouraged to refer to "
        ts = ts + "the various original publications.  We refer to the "
        ts = ts + "MCNet guidelines for event generator authors and users:\n"
        ts = ts + indent 
        ts = ts + "http://www.montecarlonet.org/index.php?p="
        ts = ts + "Publications/Guidelines\n\n"
                                                                             
        ts = ts + "Users are kindly asked to refer to the documentation "
        ts = ts + "in JHEP 02(2009)007.\n\n"
                                                                             
        ts = ts + "Please visit also our homepage\n\n"
        ts = ts + indent + indent + indent + "http://www.sherpa-mc.de\n\n"
        ts = ts + "for news, bugreports, updates and new releases." 

        textbuffer.set_text(ts)
        entry.show()
        box.show()
        box.pack_start(entry)
        return box
