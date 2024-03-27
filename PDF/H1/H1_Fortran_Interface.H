#ifndef PDF_H1_Fortran_Interface_H
#define PDF_H1_Fortran_Interface_H

#include "ATOOLS/Phys/Flavour.H"
#include "PDF/Main/Photon_PDF_Base.H"

namespace PDF {
  class H1_Fortran_Interface : public PDF_Base {
  private:
    double m_x, m_Q2;

  public:
    H1_Fortran_Interface(const std::string& _set);
    PDF_Base* GetCopy();

    void   CalculateSpec(const double&, const double&);
    double GetXPDF(const ATOOLS::Flavour& flavour);
    double GetXPDF(const kf_code& kf, bool anti);
  };

  /*!
    \class H1_Fortran_Interface
    \brief An interface to the H1 DPDF fit of parton distributions in the
    pomeron.
  */
}

#endif
