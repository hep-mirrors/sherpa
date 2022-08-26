#ifndef PDF_GRS_GRSph_Fortran_Interface_H
#define PDF_GRS_GRSph_Fortran_Interface_H

#include "ATOOLS/Phys/Flavour.H"
#include "PDF/Main/Photon_PDF_Base.H"

namespace PDF {
class GRSph_Fortran_Interface : public Photon_PDF_Base {
private:
  std::string m_path;

public:
  GRSph_Fortran_Interface(const ATOOLS::Flavour, const std::string _set);
  PDF_Base *GetCopy();

  void CalculateSpec(const double &, const double &) override;
};

/*!
  \class GRSph_Fortran_Interface
  \brief An interface to the GRS photon distribution.
*/

} // namespace PDF

#endif
