#ifndef PDF_Handler_H
#define PDF_Handler_H

#include "PDF_Base.H"
#include "Flavour.H"
#include "Data_Read.H"

namespace PDF {
  class PDF_Handler {
    bool m_initlhapdf;
  public:
    PDF_Handler() : m_initlhapdf(0) {};
    PDF_Base * GetPDFLib(AORGTOOLS::Data_Read *,APHYTOOLS::Flavour & flav,const int);
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
    \fn PDF_Base * PDF_Handler::GetPDFLib(const AORGTOOLS::Data_Read *,const int)
    
    This is the only method in the PDF_Handler. According to the file read in through
    the Data_Read object, a suitable parton density or structure function is instantiated and 
    returned. So far only a structure function for electrons (that can handle charged
    leptons) and protons is available.

    In case there is a mismatch of input and options, the program will terminate.

    \todo We will have to include also photon PDFs and maybe even pion PDFs.
          Also, we might want to use isospin argmunets to map the proton PDFs to
	  neutron ones.
  */
}

#endif // PDF_Handler_H
