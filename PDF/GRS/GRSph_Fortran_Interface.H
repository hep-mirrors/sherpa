#ifndef PDF_GRS_GRSph_Fortran_Interface_H
#define PDF_GRS_GRSph_Fortran_Interface_H

#include "ATOOLS/Phys/Flavour.H"
#include "PDF/Main/PDF_Base.H"

namespace PDF {
class GRSph_Fortran_Interface : public PDF_Base {
private:
  int m_iset;
  std::string m_path;
  double m_d, m_u, m_s, m_g;
  int m_nf;

public:
  GRSph_Fortran_Interface(const ATOOLS::Flavour, const std::string _set);
  PDF_Base *GetCopy();

  void CalculateSpec(const double &, const double &);
  double GetXPDF(const ATOOLS::Flavour &);
  double GetXPDF(const kf_code &, bool);

  double GetXMin() { return 1.e-5; }
  double GetXMax() { return 1.; }
  double GetQ2Min() { return 0.4; }
  double GetQ2Max() { return 1.e6; }
};

/*!
  \class GRSph_Fortran_Interface
  \brief An interface to the GRS photon distribution.
*/

} // namespace PDF

#endif
