#ifndef PDF_SAL_SALPH_CPP_INTERFACE_H_
#define PDF_SAL_SALPH_CPP_INTERFACE_H_

#include "ATOOLS/Phys/Flavour.H"
#include "PDF/Main/Photon_PDF_Base.H"

namespace PDF {
class SALph_CPP_Interface : public Photon_PDF_Base {
public:
  SALph_CPP_Interface(const ATOOLS::Flavour);
  PDF_Base *GetCopy();

  void CalculateSpec(const double &, const double &) override;
};
} // namespace PDF

#endif /* PDF_SAL_SALPH_CPP_INTERFACE_H_ */
