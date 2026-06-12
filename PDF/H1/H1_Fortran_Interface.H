#ifndef PDF_H1_Fortran_Interface_H
#define PDF_H1_Fortran_Interface_H

#include "ATOOLS/Phys/Flavour.H"
#include "PDF/Main/Photon_PDF_Base.H"

namespace PDF {
  class H1_Fortran_Interface : public PDF_Base {
  private:
    double m_x, m_Q2;

    static constexpr int NUM_FLAVORS = 13;
    std::array<double, NUM_FLAVORS> m_xpdf;
    std::array<bool, NUM_FLAVORS> m_flavor_calculated;

  public:
    H1_Fortran_Interface(const std::string& _set);
    PDF_Base* GetCopy();

    void   CalculateSpec(const double&, const double&);
    double GetXPDF(const ATOOLS::Flavour& flavour);
    double GetXPDF(const kf_code& kf, bool anti);
  };

  /*!
    \class H1_Fortran_Interface
    \brief Interface to H1 DPDF fit of parton distributions in the pomeron

    Reference: Eur.Phys.J.C 48 (2006) 715-748
  */
}

#endif
