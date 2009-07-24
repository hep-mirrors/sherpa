#ifndef PDF_GRV_GRVph_Fortran_Interface_H
#define PDF_GRV_GRVph_Fortran_Interface_H

#include "PDF/Main/PDF_Base.H"
#include "ATOOLS/Phys/Flavour.H"

namespace PDF {
  class GRVph_Fortran_Interface : public PDF_Base {
  private:
    //std::string m_set,m_path;
    float   m_d,m_u,m_s,m_c,m_b,m_g;
  public:
    GRVph_Fortran_Interface(const ATOOLS::Flavour);
    PDF_Base * GetCopy();

    void   CalculateSpec(double,double);
    double GetXPDF(const ATOOLS::Flavour );

    double GetXMin()   { return 1.e-5; }
    double GetXMax()   { return 1.; }
    double GetQ2Min()  { return 0.25; }
    double GetQ2Max()  { return 1.e6; }
  };


  /*!
    \class GRVph_Fortran_Interface
    \brief An interface to the GRV photon distribution.
  */
  /*!
    \var std::string GRVph_Fortran_Interface::m_path
    The path from where the initial pdfs and grids are to be loaded in. This will be handed over to
    the corresponding Fortran subroutine as a char*.
  */
}
 
#endif
