#ifndef PDF_Base_H
#define PDF_Base_H

#include "Flavour.H"

namespace PDF {
  class PDF_Base {
  public:
    // calculates distributions
    virtual void Calculate(const double x, const double Q)=0;
    // returns x times the parton distribution corresponding to the given flavour
    virtual double GetXPDF(const APHYTOOLS::Flavour &)=0;

    virtual PDF_Base                      * GetCopy()  = 0;
    virtual double                          GetXMin()  = 0;
    virtual double                          GetXMax()  = 0;
    virtual double                          GetQ2Min() = 0;
    virtual double                          GetQ2Max() = 0;
    virtual APHYTOOLS::Flavour              Beam()     = 0;
    virtual std::vector<APHYTOOLS::Flavour> Partons()  = 0;
    virtual double                          Exponent() { return 1.; }
  };
} // namespace PDF
 
#endif // PDF_Base_H
