#!/usr/bin/python

import pygtk
import gtk

import GUI_BaseClass as guibase
import ParameterBox
import PDFBox

class PDFsection_gui(guibase.gui_object):
    def initialise(self):
        self.pdfs = self.parameters.getPDFBox()

    def getContent(self):
        pdfbox    = gtk.VBox(False,0)
        pdfbox.show()
        pdfbox.pack_start(self.makePDFBox(),False,False,10)
        return pdfbox

    def makePDFBox(self):
        optbox    = gtk.HBox(False,spacing=10)
        optbox.set_border_width(10)
        optbox.show()
        label      = gtk.Label()
        label.set_text('PDFs')
        opt1       = gtk.OptionMenu()
        self.menu1 = gtk.Menu()
        opt1.set_menu(self.menu1)
        optbox.pack_start(label,False,False,2)
        optbox.pack_start(opt1,False,False,10)
        label2      = gtk.Label()
        label2.set_text('PDF for Beam 2')
        opt2       = gtk.OptionMenu()
        self.menu2 = gtk.Menu()
        opt2.set_menu(self.menu2)
        optbox.pack_end(opt2,False,False,10)
        optbox.pack_end(label2,False,False,10)
        return optbox

    def updateOptions(self):
        for i in self.menu1.get_children():
            self.menu1.remove(i)
        for i in self.menu2.get_children():
            self.menu2.remove(i)
        pdfs = self.pdfs.getOptions()
        if (len(pdfs)==2): 
            tag = 1
        else:
            tag = 0
        for pdf in pdfs[0]:
            print pdf[0],"|",pdf[1],"|"
            item = self.make_menu_item(pdf[1],self.pdf1_select,pdf[0])
            self.menu1.append(item)
        for pdf in pdfs[tag]:
            print pdf[0],"|",pdf[1],"|"
            item = self.make_menu_item(pdf[1],self.pdf2_select,pdf[0])
            self.menu2.append(item)
        if (tag==0):
            self.menu2.set_sensitive(False)
        else: 
            self.menu2.set_sensitive(True)

    def pdf1_select(self,item,tag):
        pass

    def pdf2_select(self,item,tag):
        pass
