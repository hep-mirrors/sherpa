#ifndef PDF_SASG_SASGph_Fortran_Interface_H
#define PDF_SASG_SASGph_Fortran_Interface_H

#include "ATOOLS/Phys/Flavour.H"
#include "PDF/Main/Photon_PDF_Base.H"

namespace PDF {
class SASGph_Fortran_Interface : public Photon_PDF_Base {
private:
  // PDF parameters
  float P2 = 0;
  int IP2 = 0;

public:
  SASGph_Fortran_Interface(const ATOOLS::Flavour, const std::string _set);
  PDF_Base *GetCopy() override;

  void CalculateSpec(const double &, const double &) override;
};

/*!
 \class SASGph_Fortran_Interface
 \brief An interface to the SASG photon distribution.
 */
/*
 * See https://arxiv.org/pdf/hep-ph/9503384.pdf
 */
} // namespace PDF

#endif
