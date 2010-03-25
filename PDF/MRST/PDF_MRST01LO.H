#ifndef PDF_MRST_PDF_MRST01LO_H
#define PDF_MRST_PDF_MRST01LO_H

#include <iostream>
#include "PDF/Main/PDF_Base.H"

namespace PDF {

  class PDF_MRST01LO : public PDF_Base {
  private:

    std::string m_path;

    double p_xpdfv[2], p_xpdf[6];

    int  m_anti, m_mode;
    bool m_overscaled;

  public:

    PDF_MRST01LO(const ATOOLS::Flavour bunch); 

    PDF_Base * GetCopy();

    void   CalculateSpec(double,double);
    double GetXPDF(const ATOOLS::Flavour); 

  };// end of class PDF_MRST01LO

}  

#endif
