#ifndef PDF_CJKL_CJKph_Fortran_Interface_H
#define PDF_CJKL_CJKph_Fortran_Interface_H

#include "ATOOLS/Phys/Flavour.H"
#include "PDF/Main/Photon_PDF_Base.H"

namespace PDF {
class CJKph_Fortran_Interface : public Photon_PDF_Base {
private:
  std::string m_path;

public:
  CJKph_Fortran_Interface(const ATOOLS::Flavour, const std::string _set);
  PDF_Base *GetCopy();

  void CalculateSpec(const double &, const double &);
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
