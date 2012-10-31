#!/usr/bin/python

import pygtk
import gtk

class gui_object():
    def __init__(self,params,fname="gui_object",flabel="unspecified",gui_id=0):
        print "Init gui (",gui_id,") with name = ",fname,", label = ",flabel
        self.type       = "gui_object" 
        self.id         = gui_id
        self.fname      = fname
        self.flabel     = flabel
        self.parameters = params
        self.initialise()
        self.update     = None

    def initialise(self):
        pass

    def type(self):
        return self.type

    def getId(self):
        return self.id

    def getName(self):
        return str(self.fname)

    def getLabel(self):
        return self.flabel

    def getFrame(self,borderwidth=10,xsize=200,ysize=200,update=None):
        frame = gtk.Frame(self.flabel)
        frame.set_border_width(borderwidth)
        frame.set_size_request(xsize,ysize)
        frame.add(self.getContent())
        self.update = update
        return frame

    def getContent(self):
        return gtk.Label("needs to be filled")

    def make_menu_item(self, name, callback, data=None):
        item = gtk.MenuItem(name)
        item.connect("activate", callback, data)
        item.show()
        return item
