#!/usr/bin/python

import pygtk
import gtk

import GUI_BaseClass as guibase
import ParameterBox
import ProcessBox

class Selsection_gui(guibase.gui_object):
    def initialise(self):
        self.procbox   = self.parameters.getProcBox()


