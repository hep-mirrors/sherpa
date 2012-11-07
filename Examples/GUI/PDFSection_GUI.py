#!/usr/bin/python

import pygtk
import gtk

import GUI_BaseClass as guibase
import ParameterBox
import PDFBox

class PDFsection_gui(guibase.gui_object):
    def initialise(self):
        self.pdfs   = self.parameters.getPDFBox()
        self.opt    = [gtk.OptionMenu(),gtk.OptionMenu()]
        self.menu   = [gtk.Menu(),gtk.Menu()]
        for i in range(0,2):
            self.opt[i].set_menu(self.menu[i])
            self.menu[i].show()
            self.opt[i].show()
        self.label  = [gtk.Label('PDFs'),gtk.Label('PDF for Beam 2')]
        self.label[0].set_size_request(60,30)
        self.label[1].set_size_request(100,30)
        self.adj   = [gtk.Adjustment(0,0,0,1,0),
                      gtk.Adjustment(0,0,0,1,0)]
        self.set   = [gtk.SpinButton(self.adj[0],1,0),
                      gtk.SpinButton(self.adj[1],1,0)]
        self.set[0].set_size_request(50,30)
        self.set[1].set_size_request(50,30)
        self.setstr = [gtk.Label('Set'),gtk.Label('Set')]
        self.setstr[0].set_size_request(40,30)
        self.setstr[1].set_size_request(40,30)
        self.menu[0].set_sensitive(True)
        self.set[0].set_sensitive(True)
        self.menu[1].set_sensitive(False)
        self.set[1].set_sensitive(False)
        self.updateOptions()

    def getContent(self):
        pdfbox    = gtk.VBox(False,0)
        pdfbox.show()
        self.updateOptions()
        pdfbox.pack_start(self.makePDFBox(),False,False,10)
        return pdfbox

    def makePDFBox(self):
        optbox    = gtk.HBox(False,spacing=10)
        optbox.set_border_width(10)
        optbox.show()
        for i in range(0,2):
            if i==0:
                optbox.pack_start(self.label[i],False,False,2)
                optbox.pack_start(self.opt[i],False,False,10)
                optbox.pack_start(self.setstr[i],False,False,2)
                optbox.pack_start(self.set[i],False,False,10)
            else:
                optbox.pack_end(self.set[i],False,False,10)
                optbox.pack_end(self.setstr[i],False,False,2)
                optbox.pack_end(self.opt[i],False,False,10)
                optbox.pack_end(self.label[i],False,False,2)
        return optbox

    def updateOptions(self):
        pdfs = self.pdfs.getOptions()
        for beam in range(0,len(pdfs)):
            for i in self.menu[beam].get_children():
                self.menu[beam].remove(i)
            for pdf in pdfs[beam]:
                item = self.make_menu_item(pdf[1],
                                           self.selectPDF,[pdf[1],beam])
                self.menu[beam].append(item)
            self.adj[beam].set_value(pdf[2])
            self.adj[beam].set_lower(pdf[3])
            self.adj[beam].set_upper(pdf[4])
            if self.pdfs.getType(beam)=='lepton':
                self.setstr[beam] = gtk.Label('Ord')
            if self.pdfs.getType(beam)=='hadron':
                self.setstr[beam] = gtk.Label('Set')
        if (len(pdfs)==2): 
            self.menu[1].set_sensitive(True)
            self.set[1].set_sensitive(True)
        else:
            self.menu[1].set_sensitive(False)
            self.set[1].set_sensitive(False)


    def selectPDF(self,item,tags):
        print "Select PDF for beam ",tags[1],": ",tags[0]
        self.pdfs.setPDF(tags[1],tags[0])

    def extractParameters(self):
        for bid in range(0,2):
            if self.set[bid].is_sensitive():
                self.pdfs.setPDFset(bid,self.set[bid].get_value_as_int())
        self.pdfs.printStatus()
