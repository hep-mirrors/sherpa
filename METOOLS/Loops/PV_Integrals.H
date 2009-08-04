#ifndef METOOLS_Loops_PV_Integrals_H
#define METOOLS_Loops_PV_Integrals_H

#include "ATOOLS/Math/MyComplex.H"
#include "METOOLS/Loops/Divergence_Array.H"
#include "METOOLS/Loops/Master_Integrals.H"

namespace METOOLS {

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
        [5] - \f$ \sim \epsilon^2 \f$
      All terms \f$ \sim \epsilon \f$ are currently not calculated.
      This is suitable for 1-loop integrals, but nothing beyond.
   */
  /*! This file declares the scalar parts of higher order tensor
      integrals in the Passarino-Veltman reduction scheme.
   */


  //! A_munu = g_munu A_2
  //!   A_2(m2;mu2)
  DivArrC
  PV_Tadpole_2(const Complex&,
               double);


  //! B_mu = p_mu B_1
  //!   B_1(s12;m02,m12;mu2)
  DivArrC
  PV_Bubble_1(const double&,
              const Complex&, const Complex&,
              double);
  //! B_munu = p_mu p_nu B_21 + g_munu B_22
  //!   B_21(s12;m02,m12;mu2)
  DivArrC
  PV_Bubble_21(const double&,
               const Complex&, const Complex&,
               double);
  //!   B_22(s12;m02,m12;mu2)
  DivArrC
  PV_Bubble_22(const double&,
               const Complex&, const Complex&,
               double);
  //! B_munurho = p_mu p_nu p_rho B_31 + {g,p}_munurho B_32
  //!   B_31(s12;m02,m12;mu2)
  DivArrC
  PV_Bubble_31(const double&,
               const Complex&, const Complex&,
               double);
  //!   B_32(s12;m02,m12;mu2)
  DivArrC
  PV_Bubble_32(const double&,
               const Complex&, const Complex&,
               double);


  //! C_mu = p_mu C_10
  //!   C_10(p2,p2,0;m02,m12,m22;mu2)
  DivArrC
  PV_Triangle_10(const double&,  const double&,  const double&,
                 const Complex&, const Complex&, const Complex&,
                 double);
  //! C_mu = p_1,mu C_11 + p_2,nu C_12
  //!   C_11(p12,p22,s12;m02,m12,m22;mu2)
  DivArrC
  PV_Triangle_11(const double&,  const double&,  const double&,
                 const Complex&, const Complex&, const Complex&,
                 double);
  //!   C_12(p12,p22,s12;m02,m12,m22;mu2)
  DivArrC
  PV_Triangle_12(const double&,  const double&,  const double&,
                 const Complex&, const Complex&, const Complex&,
                 double);
  //! C_munu = p_mu p_nu C_20 + g_munu C_24
  //!   C_20(p12,p22,s12;m02,m12,m22;mu2)
  DivArrC
  PV_Triangle_20(const double&,  const double&,  const double&,
                 const Complex&, const Complex&, const Complex&,
                 double);
  //! C_munu = p_1,mu p_1.nu C_21 + p_2,mu p_2.nu C_22
  //!          + {p_1,p_2}_munu C_23 + g_munu C_24
  //!   C_21(p12,p22,s12;m02,m12,m22;mu2)
  DivArrC
  PV_Triangle_21(const double&,  const double&,  const double&,
                 const Complex&, const Complex&, const Complex&,
                 double);
  //!   C_22(p12,p22,s12;m02,m12,m22;mu2)
  DivArrC
  PV_Triangle_22(const double&,  const double&,  const double&,
                 const Complex&, const Complex&, const Complex&,
                 double);
  //!   C_23(p12,p22,s12;m02,m12,m22;mu2)
  DivArrC
  PV_Triangle_23(const double&,  const double&,  const double&,
                 const Complex&, const Complex&, const Complex&,
                 double);
  //!   C_24(p12,p22,s12;m02,m12,m22;mu2)
  DivArrC
  PV_Triangle_24(const double&,  const double&,  const double&,
                 const Complex&, const Complex&, const Complex&,
                 double);
  //! C_munurho = p_mu p_nu p_rho C_30 + {g,p}_munurho C_38
  //!   C_30(p12,p22,s12;m02,m12,m22;mu2)
  DivArrC
  PV_Triangle_30(const double&,  const double&,  const double&,
                 const Complex&, const Complex&, const Complex&,
                 double);
  //!   C_38(p12,p22,s12;m02,m12,m22;mu2)
  DivArrC
  PV_Triangle_38(const double&,  const double&,  const double&,
                 const Complex&, const Complex&, const Complex&,
                 double);
  //! C_munu = p_1,mu p_1,nu p_1,rho C_31 + {p_1,p_1,p_2}_munurho C_32
  //!           + {p_1,p_2,p_2}_munurho C_33 + p_2,mu p_2.nu p_1,rho C_34
  //!           + {g,p_1}_munurho C_35 + {g,p_2}_munurho C_36
  //!   C_31(p12,p22,s12;m02,m12,m22;mu2)
  DivArrC
  PV_Triangle_31(const double&,  const double&,  const double&,
                 const Complex&, const Complex&, const Complex&,
                 double);
  //!   C_32(p12,p22,s12;m02,m12,m22;mu2)
  DivArrC
  PV_Triangle_32(const double&,  const double&,  const double&,
                 const Complex&, const Complex&, const Complex&,
                 double);
  //!   C_33(p12,p22,s12;m02,m12,m22;mu2)
  DivArrC
  PV_Triangle_33(const double&,  const double&,  const double&,
                 const Complex&, const Complex&, const Complex&,
                 double);
  //!   C_34(p12,p22,s12;m02,m12,m22;mu2)
  DivArrC
  PV_Triangle_34(const double&,  const double&,  const double&,
                 const Complex&, const Complex&, const Complex&,
                 double);
  //!   C_35(p12,p22,s12;m02,m12,m22;mu2)
  DivArrC
  PV_Triangle_35(const double&,  const double&,  const double&,
                 const Complex&, const Complex&, const Complex&,
                 double);
  //!   C_36(p12,p22,s12;m02,m12,m22;mu2)
  DivArrC
  PV_Triangle_36(const double&,  const double&,  const double&,
                 const Complex&, const Complex&, const Complex&,
                 double);


}

#endif
