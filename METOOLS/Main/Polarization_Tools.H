#ifndef METOOLS_Main_Polarization_Tools_H
#define METOOLS_Main_Polarization_Tools_H

#ifndef SQRT_05
#define SQRT_05 0.70710678118654757
#endif

#include <map>
#include <iostream>
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Spinor.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Matrix.H"

namespace METOOLS
{
  typedef ATOOLS::Spinor<double> SpinorType;

  class Polarization_Vector : public std::vector<ATOOLS::Vec4C> {
    SpinorType m_kp;
    SpinorType m_km;
    ATOOLS::Vec4D m_k;

    void Init(ATOOLS::Vec4D p, ATOOLS::Vec4D ref_mom=ATOOLS::Vec4D(1.0,SQRT_05,SQRT_05,0.0));
    ATOOLS::Vec4C VT(const SpinorType &a,const SpinorType &b);
    ATOOLS::Vec4C EM(const ATOOLS::Vec4D &p);
    ATOOLS::Vec4C EP(const ATOOLS::Vec4D &p);
    ATOOLS::Vec4C EMM(const ATOOLS::Vec4D &p);
    ATOOLS::Vec4C EMP(const ATOOLS::Vec4D &p);
    ATOOLS::Vec4C EML(const ATOOLS::Vec4D &p);
  public:
    Polarization_Vector(ATOOLS::Vec4D p, double M2, bool anti=false, bool out=false);
    Polarization_Vector(ATOOLS::Vec4D p, bool anti=false, bool out=false);
    Polarization_Vector(ATOOLS::Vec4D p, ATOOLS::Vec4D ref_mom);
    void Test(ATOOLS::Vec4D p=ATOOLS::Vec4D(3.0399e+00, 1.7523e+00,-1.4883e-01, 3.1938e-01));

    vector<vector<Complex>> BasisTrafo(const Polarization_Vector& old_basis) const;
  };

  class Polarization_Tensor : public std::vector<ATOOLS::CMatrix> {
  public:
    Polarization_Tensor(ATOOLS::Vec4D p, double M2, bool anti=false, bool out=false);
    void Test(ATOOLS::Vec4D p=ATOOLS::Vec4D(3.0399e+00, 1.7523e+00,-1.4883e-01, 3.1938e-01));
  };

    /*!
      \file Polarization_Tools.H
      \brief Declares the class ATOOLS::Polarization_Vector and ATOOLS::Polarization_Tensor
      */
    /*!
      \class Polarization_Vectors
      \brief Calculation of polarization vectors used for matrix element calculation
      */
    /*!
    \fn Polarization_Vector::Init(ATOOLS::ATOOLS::Vec4D p, ATOOLS::Vec4D ref_mom=ATOOLS::Vec4D(1.0,SQRT_05,SQRT_05,0.0) )
    \brief Helper function of constructor; calculates complex circular polarization vectors of a massive spin-1 boson.

    Calculations according to arxiv: 9601359 eq. (25) or 9805445 eq. (3.19) (+ transformation back to vector
    representation);
    The ordering in the resulting Polarization_Vector is
        - \f$0 \to \epsilon^+\f$
        - \f$1 \to \epsilon^-\f$
        - \f$2 \to \epsilon^0\f$
        - \f$3 \to \epsilon^s\f$
        BUT: currently transverse polarizations seems to be exchanged
    Input: 
      - <var>p</var>: momentum of particle
      - <var>ref_mom</var>: arbitrary, light-like reference momentum, defines spin axis for massive particles
      .
    Output:
      - no output
      .
    */
  /*!
  \fn Polarization_Vector::BasisTrafo(const Polarization_Vector& old_basis) ()
  \brief Method determine transformation coefficients from one basis (e.g. spin basis given by reference vector in
         Weyl-van-der-Waerden formalism or reference frame for defining the polarisation) for defining polarisation
         vectors to another, if the new polarization vectors can be determined as a linear combination of the ones
         defined in the old basis;
  Input:
    - <var>old_basis</var>: Polarization vectors in old basis (in the most cases COMIX default polarization vectors)
    .
  Output:
    - vector<vector<Complex>>, which contains the transformation coefficients for each new polarization vector as
      one vector respectively, polarization ordering according to Polarization_Vector
   .
  */
}

#endif
