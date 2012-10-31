#!/usr/bin/python

import pygtk
import gtk


class PDFBox:
    def __init__(self):
        self.pdf_options     = []
        self.pdf_type        = []
        self.pdf_tag         = []
        self.pdf_set         = []
        self.beam_id         = [] 
        self.pdf_options.append(None)
        self.pdf_options.append(None)
        self.pdf_type.append(None)
        self.pdf_type.append(None)
        self.pdf_tag.append(None)
        self.pdf_tag.append(None)
        self.pdf_set.append(None)
        self.pdf_set.append(None)
        self.beam_id.append(None)
        self.beam_id.append(None)

        self.Lepton_PDFs = []
        self.Lepton_PDFs.append([1,"PDFe",1,1,3])
        self.Lepton_PDFs.append([0,"None",0,0,0])

        self.Hadron_PDFs = []
        self.Hadron_PDFs.append([0,"CT10",0,0,50])
        self.Hadron_PDFs.append([1,"NNPDF2.3",0,0,100])
        self.Hadron_PDFs.append([2,"MSTW10",0,0,50])

        self.default_lepton_tag = "PDFe"
        self.default_hadron_tag = "CT10"
        self.pdf_options[0]     = None
        self.pdf_options[1]     = None
        self.pdf_type[0]        = None
        self.pdf_type[1]        = None

        self.pdf_tag_all        = None
        self.pdf_set_all        = None
        self.pdf_library        = ""

    def initialiseDefaults(self,collider):
        self.beam_id[0] = collider[2]
        self.beam_id[1] = collider[4]
        print "Collider with beams: ",self.beam_id[0],self.beam_id[1]
        tag = None
        for bid in range (0,2):
            if (self.beam_id[bid]==11 or self.beam_id[bid]==-11): 
                self.pdf_options[bid] = self.Lepton_PDFs
                self.pdf_type[bid]    = "lepton"
                tag                   = self.Lepton_PDFs[0][1]
            if (self.beam_id[bid]==2212 or self.beam_id[bid]==-2212): 
                self.pdf_options[bid] = self.Hadron_PDFs
                self.pdf_type[bid]    = "hadron"
                tag                   = self.Hadron_PDFs[0][1]
            if tag==None:
                pass
            else:
                for pdf in self.pdf_options[bid]:
                    if pdf[1]==tag:
                        self.pdf_tag[bid] = pdf[1]
                        self.pdf_set[bid] = pdf[2]

        print "Initialised default pdfs for ",collider[2],collider[4]
        print "   ",self.pdf_type[0],": (PDF = ",self.pdf_tag[0]," set = ",self.pdf_set[0],")"
        print "   ",self.pdf_type[1],": (PDF = ",self.pdf_tag[1]," set = ",self.pdf_set[1],")"

    def getOptions(self):
        options = []
        options.append(self.pdf_options[0])
        if (self.pdf_type[0]!=self.pdf_type[1]):
            options.append(self.pdf_options[1])
        return options

    def setPDF(self,bid,pdf_tag,pdf_set):
        self.pdf_tag[bid-1] = pdf_tag
        self.pdf_set[bid-1] = pdf_set
        print "Set PDF(",bid,") = ",pdf_tag,", set = ",pdf_set

    def getPDF(self,bid):
        return self.pdf_tag[bid-1],self.pdf_set[bid-1]

    def setLibrary(self,lib):
        if lib=="LHAPDF" or lib=="":
            self.pdf_library = lib

    def getLibrary(self):
        return self.pdf_library

    def printDefaults(self):
        print "Default settings in the PDF box:"
        print "Beam1: (",pdf_tag[0],pdf_set[0],")"
        print "Beam1: (",pdf_tag[1],pdf_set[1],")"

    def write(self,runfile,collider):
        runfile.write('\n')
        runfile.write('%%% PDF setup \n\n')
        pdf_tag1,pdf_set1 = self.getPDF(1)
        pdf_tag2,pdf_set2 = self.getPDF(2)

        print "Check consistency for beams: ",collider[2],collider[4]
        if (collider[2]==collider[4] or
            collider[2]==-collider[4]):
            print "Try to write PDFs: ",pdf_tag1,pdf_tag2
            if pdf_tag1!=pdf_tag2 or pdf_set1!=pdf_set2:
                print "Inconsistent PDFs, use PDF_SET_1 globally"
            runfile.write('PDF_SET = %20s ; ' % str(pdf_tag1))
        else:
            runfile.write('PDF_SET_1 = %20s ; ' % str(pdf_tag1))
            runfile.write('PDF_SET_2 = %20s ; ' % str(pdf_tag2))
        runfile.write('\n')
