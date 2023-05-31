#ifndef PDF_CJKL_CJKph_Fortran_Interface_H
#define PDF_CJKL_CJKph_Fortran_Interface_H

#include "ATOOLS/Phys/Flavour.H"
#include "PDF/Main/PDF_Base.H"

namespace PDF {
class CJKph_Fortran_Interface : public PDF_Base {
private:
  std::string m_path;
  double m_d, m_u, m_s, m_c, m_b, m_g;
  int m_nf, m_iset;

public:
  CJKph_Fortran_Interface(const ATOOLS::Flavour, const std::string _set);
  PDF_Base *GetCopy();

  void CalculateSpec(const double &, const double &);
  double GetXPDF(const ATOOLS::Flavour &);
  double GetXPDF(const kf_code&, bool);

    double GetXMin()   { return 1.e-5; }
    double GetXMax()   { return 1.; }
    double GetQ2Min()  { return 1.; }
    double GetQ2Max()  { return 2.e5; }
  };

  /*!
    \class CJKph_Fortran_Interface
    \brief An interface to the CJK photon distributions.
  */
  /*!
    \var std::string CJKph_Fortran_Interface::m_path
    The path from where the initial pdfs and grids are to be loaded in. This
    will be handed over to the corresponding Fortran subroutine as a char*.
  */
  }

#endif
