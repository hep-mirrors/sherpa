#ifndef PDF_Handler_H
#define PDF_Handler_H

#include "PDF_Base.H"
#include "Flavour.H"

namespace ATOOLS { class Data_Reader; }

namespace PDF {
  class PDF_Handler {
    bool m_initlhapdf;
  public:
    PDF_Handler() : m_initlhapdf(0) {};
    PDF_Base * GetPDFLib(ATOOLS::Data_Reader *,ATOOLS::Flavour &,const int);
  };
  /*!
    \namespace PDF
  */
  /*!
    \file PDF_Handler.H
    contains the PDF_Handler.H
  */
  /*! 
    \class PDF_Handler
    \brief instantiates and returns PDFs according to an input file.

    This class has the only purpose of steering the selection of PDFs or structure
    functions once an input file is read in. It has only one method, GetPDFLib,
    that uses a Data_Read object to determine the specific parton density or
    structure function an user might want to instantiate and employ. It then returns
    a pointer to this object.
  */
  /*!
    \var bool PDF_Handler::m_initlhapdf
    Gives the status concerning LHAPDF to avoid double initialization. A similar thing has to
    be done for every Fortran PDF.
    \todo Add a flag for cteq6, Fortran version.
  */
  /*!
    \fn PDF_Handler::PDF_Handler()
    The constructor is basically empty, it just sets the "is already initialised"-flags.
  */
  /*!
    \fn PDF_Base * PDF_Handler::GetPDFLib(const ATOOLS::Data_Read *,ATOOLS::Flavour &,const int)
    
    This is the only method in the PDF_Handler. According to the file read in through
    the Data_Read object and according to the actual flavour, a suitable parton density or structure 
    function is instantiated and returned. So far only a structure function for electrons (that can 
    handle charged leptons) and various pdfs for protons are available.

    In case there is a mismatch of input and options, the program will terminate.

    \todo We will have to include also photon PDFs and maybe even pion PDFs.
          Also, we might want to use isospin arguments to map the proton PDFs to
	  neutron ones.
  */
}

#endif // PDF_Handler_H
