#ifndef AddOns_OpenLoops_Color_Correlated_ME2_H
#define AddOns_OpenLoops_Color_Correlated_ME2_H

#include "PHASIC++/Process/Color_Correlated_ME2.H"

namespace OpenLoops {
  
  class Color_Correlated_ME2 : public PHASIC::Color_Correlated_ME2 {

  public:
    
    Color_Correlated_ME2(const PHASIC::External_ME_Args& args,
			 int ol_id, const AmplitudeType& type);

    ~Color_Correlated_ME2();

    /* Calculate full correlation matrxix */
    void Calc(const ATOOLS::Vec4D_Vector &p);

    /* GetValue(i,j) = <1,2,3,... | TiTj | ...,3,2,1 > */
    double GetValue(const size_t& i, const size_t& j) const;

    /* Squared uncorrelated amplitude */
    double GetBorn2() const;

  private:

    int m_ol_id;
    double m_symfac;
    AmplitudeType m_amptype;
    const AmplitudeType AmpType() const { return m_amptype; }

    /* Number of external legs */
    size_t m_dim;

    /* The fortran 2D arrays OpenLoops deals with are just 1D arrays
       with fancy index arithmetics. Therefore store color-correlation
       matrix in 1D array and do arithmetics in GetValue method */
    double* m_ccmatrix;

    /* When evaluating the color-correlation matrix, the squared born
       comes for free in OpenLoops. Store it here for convenience. */
    double m_born2;

  };

}

#endif
