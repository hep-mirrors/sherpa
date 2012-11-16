#!/usr/bin/python

import pygtk
import gtk


class PDFBox:
    def __init__(self):
        self.pdf_options = ["None","None"]
        self.pdf_type    = ["",""]
        self.pdf         = [None,None]
        self.pdf_tag     = ["",""]
        self.pdf_set     = [0,0]
        self.beam_id     = [0,0] 

        self.No_PDFs     = []
        self.Lepton_PDFs = []
        self.Hadron_PDFs = []
        self.fillPDFs()

        self.pdf_library = ""

    def initialiseDefaults(self,collider):
        self.beam_id[0] = collider[2]
        self.beam_id[1] = collider[4]
        print "PDF::Collider with beams: ",self.beam_id[0],self.beam_id[1]
        tag = None
        for bid in range (0,2):
            if (self.beam_id[bid]==0): 
                self.pdf_options[bid] = self.No_PDFs
                self.pdf_type[bid]    = "None"
            if (self.beam_id[bid]==11 or self.beam_id[bid]==-11): 
                self.pdf_options[bid] = self.Lepton_PDFs
                self.pdf_type[bid]    = "lepton"
            if (self.beam_id[bid]==2212 or self.beam_id[bid]==-2212): 
                self.pdf_options[bid] = self.Hadron_PDFs
                self.pdf_type[bid]    = "hadron"
            print "   bid = ",bid,self.pdf_type[bid],self.pdf_options[bid][0][0]
            self.pdf[bid]     = self.pdf_options[bid][0]
            self.pdf_tag[bid] = self.pdf_options[bid][0][1]
            self.pdf_set[bid] = self.pdf_options[bid][0][2]

        print "PDF::Initialised default pdfs for ",collider[2],collider[4]
        print "   ",self.pdf_type[0],": (",self.pdf_tag[0],self.pdf_set[0],")"
        print "   ",self.pdf_type[1],": (",self.pdf_tag[1],self.pdf_set[1],")"

    def getOptions(self):
        options = []
        options.append(self.pdf_options[0])
        if (self.pdf_type[0]!=self.pdf_type[1]):
            options.append(self.pdf_options[1])
        print "PDF::getOptions(",len(options),")"
        return options

    def setPDF(self,bid,pdf_tag,pdf_set=0):
        self.pdf_tag[bid] = pdf_tag
        self.pdf_set[bid] = pdf_set
        for pdf in self.pdf_options[bid]:
            if pdf[1]==pdf_tag:
                self.pdf[bid] = pdf
        if (self.pdf_type[0]==self.pdf_type[1]): 
            self.pdf_tag[1-bid] = self.pdf_tag[bid]
            self.pdf_set[1-bid] = self.pdf_set[bid]
            self.pdf[1-bid]     = self.pdf[bid]
        print "Tried to set PDFs from types",self.pdf_type[0],self.pdf_type[1] 
        print "   ==> flag = ",(self.pdf_type[0]==self.pdf_type[1]),bid
        print "Set PDF(",(bid),"):",self.pdf_tag[bid],self.pdf_set[bid]
        print "Set PDF(",(1-bid),"):",self.pdf_tag[1-bid],self.pdf_set[1-bid]

    def setPDFset(self,bid,pdf_set):
        self.pdf_set[bid] = pdf_set
        if (self.pdf_type[0]==self.pdf_type[1]): 
            self.pdf_set[1-bid] = pdf_set
        print "Set PDFset(",(bid),"):",self.pdf_tag[bid],self.pdf_set[bid]
        print "Set PDFset(",(1-bid),"):",self.pdf_tag[1-bid],self.pdf_set[1-bid]

    def getPDF(self,bid):
        return self.pdf_tag[bid]

    def getType(self,bid):
        return self.pdf_type[bid]

    def setLibrary(self,lib):
        if lib=="LHAPDF" or lib=="":
            self.pdf_library = lib

    def getLibrary(self):
        return self.pdf_library

    def printStatus(self):
        print "Settings in the PDF box:"
        print "Beam1: (",self.pdf_tag[0],self.pdf_set[0],")"
        print "Beam2: (",self.pdf_tag[1],self.pdf_set[1],")"

    def write(self,runfile,collider):
        runfile.write('\n')
        runfile.write('  %%% PDF setup \n\n')

        print "Check consistency for beams: ",collider[2],collider[4]
        if (collider[2]==collider[4] or
            collider[2]==-collider[4]):
            runfile.write("  PDF_SET = %20s; " 
                          % str(self.pdf_tag[0]))
            if self.pdf_set[0]!=0:
                runfile.write("  PDF_SET_VERSION = %3s; "
                              % str(self.pdf_set[0]))
        else:
            for beam in range(0,2):
                runfile.write("  PDF_SET_%i = %20s ; "
                              % (beam,str(self.pdf_tag[beam])))
                if self.pdf_set[beam]!=0:
                    runfile.write("PDF_SET_%i = %20s;\n" 
                                  % (beam,str(self.pdf_set[beam])))

        runfile.write('\n')

    def fillPDFs(self):
        self.No_PDFs.append([0,"No PDF",0,0,0])

        self.Lepton_PDFs.append([1,"PDFe",0,0,3])
        self.Lepton_PDFs.append([0,"None",0,0,0])

        self.Hadron_PDFs.append([0,"CT10",0,0,50])
        self.Hadron_PDFs.append([1,"NNPDF2.3",0,0,100])
        self.Hadron_PDFs.append([2,"MSTW10",0,0,50])
