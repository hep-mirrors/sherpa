#ifndef HELICITIES_Loops_Master_Integrals_H
#define HELICITIES_Loops_Master_Integrals_H

#include "ATOOLS/Math/MyComplex.H"
#include "HELICITIES/Loops/Divergence_Array.H"

namespace HELICITIES {

  /*! Convention on input:
      - \f$ p_i^2 \f$ labels the squared outer momenta
      - \f$ m_i^2 \f$ labels the masses of the internal lines
      - \f$ s_{ij}=(p_i + p_j)^2 \f$
      - \f$ \mu^2 \f$ labels the renormalisation scale
      Convention on output:
      - array with the following entries
        [0] - \f$ \sim \frac{1}{\epsilon_{UV}} \f$
        [1] - \f$ \sim \frac{1}{\epsilon_{IR}} \f$
        [2] - \f$ \sim \frac{1}{\epsilon_{IR}^2} \f$
        [3] - \f$ \sim 1 \f$
        [4] - \f$ \sim \epsilon \f$
      All terms \f$ \sim \epsilon \f$ are currently not calculated.
      This is suitable for 1-loop integrals, but nothing beyond.
   */
  //! A_0(m2;mu2)
  DivArrC
  Master_Tadpole(const Complex&,
                 const double&);
  //! B_0(s12;m02,m12;mu2)
  DivArrC
  Master_Bubble(const double&,
                const Complex&, const Complex&,
                const double&);
  //! C_0(p12,p22,s12;m02,m12,m22;mu2)
  DivArrC
  Master_Triangle(const double&,  const double&,  const double&,
                  const Complex&, const Complex&, const Complex&,
                  const double&);
  //! D_0(p12,p22,p32,p43,s12,s23;m12,m22,m32,m42;mu2)
  DivArrC
  Master_Box(const double&,  const double&,  const double&,  const double&,
             const double&,  const double&,
             const Complex&, const Complex&, const Complex&, const Complex&,
             const double&);


}

#endif
