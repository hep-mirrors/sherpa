#ifndef PDF_Electron_H
#define PDF_Electron_H

#include "PDF_Base.H"

namespace PDF {
  class PDF_Electron : public PDF_Base {
    // Electron or Positron
    APHYTOOLS::Flavour beam;
    std::vector<APHYTOOLS::Flavour> partons;    
    double mass;
    double alpha,beta;
    double xpdf;
  public:
    // initialisation
    PDF_Electron(int mode);
    // calculatation
    // returns x times the parton distribution corresponding to the given flavour
    void Calculate(const double x, const double Q2);

    // here only electron and positron is possible
    double GetXPDF(const APHYTOOLS::Flavour & fl) {
      if (fl==beam) return xpdf;
      return 0.;
    };
    PDF_Base * GetCopy() {
      if (beam==APHYTOOLS::Flavour(APHYTOOLS::kf::e)) return new PDF_Electron(0);
      return new PDF_Electron(1);
    }
    // return domain for x and Q2 of the PDF (guess)
    double GetXMin()                          { return 1.e-6; }
    double GetXMax()                          { return .999999; }
    double GetQ2Min()                         { return 0.25; }
    double GetQ2Max()                         { return 1.e14; }
    double Exponent()                         { return beta; }
    APHYTOOLS::Flavour Beam()                 { return beam; }
    std::vector<APHYTOOLS::Flavour> Partons() { return partons; }
  };
} // namespace PDF

#endif // PDF_Electron_H
