#ifndef PDF_Handler_H
#define PDF_Handler_H

#include "PDF_Base.H"
#include "Flavour.H"

namespace PDF {

  //! provides a PDF Library for a given Flavour (Proton, Electron, ...)
  class PDF_Handler {
  public:
    PDF_Base * GetPDFLib(const APHYTOOLS::Flavour);
  };

}

#endif // PDF_Handler_H
