#ifndef PDF_GRV_GRVph_Fortran_Interface_H
#define PDF_GRV_GRVph_Fortran_Interface_H

#include "ATOOLS/Phys/Flavour.H"
#include "PDF/Main/Photon_PDF_Base.H"

namespace PDF {
class GRVph_Fortran_Interface : public Photon_PDF_Base {
public:
  GRVph_Fortran_Interface(const ATOOLS::Flavour, const std::string);
  PDF_Base *GetCopy() override;

  void CalculateSpec(const double &, const double &) override;
};

/*!
  \class GRVph_Fortran_Interface
  \brief An interface to the GRV photon distribution.
*/
} // namespace PDF

#endif
