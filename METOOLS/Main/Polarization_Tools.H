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

    void Init(ATOOLS::Vec4D p);
    ATOOLS::Vec4C VT(const SpinorType &a,const SpinorType &b);
    ATOOLS::Vec4C EM(const ATOOLS::Vec4D &p);
    ATOOLS::Vec4C EP(const ATOOLS::Vec4D &p);
    ATOOLS::Vec4C EMM(const ATOOLS::Vec4D &p);
    ATOOLS::Vec4C EMP(const ATOOLS::Vec4D &p);
    ATOOLS::Vec4C EML(const ATOOLS::Vec4D &p);
  public:
    Polarization_Vector(ATOOLS::Vec4D p, double M2, bool anti=false, bool out=false);
    Polarization_Vector(ATOOLS::Vec4D p, bool anti=false, bool out=false);
    void Test(ATOOLS::Vec4D p=ATOOLS::Vec4D(3.0399e+00, 1.7523e+00,-1.4883e-01, 3.1938e-01));
  };

  class Polarization_Tensor : public std::vector<ATOOLS::CMatrix> {
  public:
    Polarization_Tensor(ATOOLS::Vec4D p, double M2, bool anti=false, bool out=false);
    void Test(ATOOLS::Vec4D p=ATOOLS::Vec4D(3.0399e+00, 1.7523e+00,-1.4883e-01, 3.1938e-01));
  };

    /*!
      \file Polarization_Tools.H
      \brief Declares the class ATOOLS::Polarization_Tools
      */
    /*!
      \class Polarization_Tools
      \brief A small collection of frequently used functions for ME's

      This class is completely static so that it does not need a constructor.
      Its methods can be called directly.
      */
    /*!
    \fn Polarization_Tools::ComplexBosonPolarizationVector(ATOOLS::ATOOLS::Vec4D p, int lambda, 
    double M2, ATOOLS::ATOOLS::Vec4D * eps )
    \brief Calculates complex linear polarization vectors of a massive spin-1 boson.

    Calculations according to the Diploma thesis by Tanju Gleisberg, p. 35
    Input: 
      - <var>p</var>: momentum of particle
      - <var>M2</var>: (offshell) mass
      - <var>lambda</var>
        - \f$0 \to \epsilon^+\f$
        - \f$1 \to \epsilon^-\f$
        - \f$2 \to \epsilon^0\f$
        - \f$3 \to \epsilon^s\f$
        .
      .
    Output:
      - real part of polarization vector as ATOOLS::Vec4D[0]
      - imaginary part of polarization vector as ATOOLS::Vec4D[1]
      .
    */
}

#endif
