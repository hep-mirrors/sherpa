#ifndef ATOOLS_Math_Tensor_H
#define ATOOLS_Math_Tensor_H

#include <iostream>
#include "ATOOLS/Math/Lorentz_Ten2.H"
#include "ATOOLS/Math/Lorentz_Ten3.H"
#include "ATOOLS/Math/Lorentz_Ten4.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Math/MathTools.H"

#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"


namespace ATOOLS {

  /*
   * declarations and specialisations for 4 dimensional 2nd, 3rd and 4th rank tensors of doubles
   */
  typedef Lorentz_Ten2<double>  Lorentz_Ten2D;
  typedef Lorentz_Ten3<double>  Lorentz_Ten3D;
  typedef Lorentz_Ten4<double>  Lorentz_Ten4D;

  bool IsEqual(const Lorentz_Ten2D& t1, const Lorentz_Ten2D& t2, const double crit=1.0e-12);
  bool IsEqual(const Lorentz_Ten3D& t1, const Lorentz_Ten3D& t2, const double crit=1.0e-12);
  bool IsEqual(const Lorentz_Ten4D& t1, const Lorentz_Ten4D& t2, const double crit=1.0e-12);

  inline bool operator==(const Lorentz_Ten2D& t1, const Lorentz_Ten2D& t2) {
    return IsEqual(t1,t2);
  }
  inline bool operator==(const Lorentz_Ten3D& t1, const Lorentz_Ten3D& t2) {
    return IsEqual(t1,t2);
  }
  inline bool operator==(const Lorentz_Ten4D& t1, const Lorentz_Ten4D& t2) {
    return IsEqual(t1,t2);
  }
  inline bool operator!=(const Lorentz_Ten2D& t1, const Lorentz_Ten2D& t2) {
    return !IsEqual(t1,t2);
  }
  inline bool operator!=(const Lorentz_Ten3D& t1, const Lorentz_Ten3D& t2) {
    return !IsEqual(t1,t2);
  }
  inline bool operator!=(const Lorentz_Ten4D& t1, const Lorentz_Ten4D& t2) {
    return !IsEqual(t1,t2);
  }

  // **
  // declarations and specialisations for 4 dimensional 2nd, 3rd and 4th rank tensors of Complex
  // **
  typedef Lorentz_Ten2<Complex> Lorentz_Ten2C;
  typedef Lorentz_Ten3<Complex> Lorentz_Ten3C;
  typedef Lorentz_Ten4<Complex> Lorentz_Ten4C;

  bool IsEqual(const Lorentz_Ten2C& t1, const Lorentz_Ten2C& t2, const double crit=1.0e-12);
  bool IsEqual(const Lorentz_Ten3C& t1, const Lorentz_Ten3C& t2, const double crit=1.0e-12);
  bool IsEqual(const Lorentz_Ten4C& t1, const Lorentz_Ten4C& t2, const double crit=1.0e-12);

  inline bool operator==(const Lorentz_Ten2C& t1, const Lorentz_Ten2C& t2) {
    return IsEqual(t1,t2);
  }
  inline bool operator==(const Lorentz_Ten3C& t1, const Lorentz_Ten3C& t2) {
    return IsEqual(t1,t2);
  }
  inline bool operator==(const Lorentz_Ten4C& t1, const Lorentz_Ten4C& t2) {
    return IsEqual(t1,t2);
  }
  inline bool operator!=(const Lorentz_Ten2C& t1, const Lorentz_Ten2C& t2) {
    return !IsEqual(t1,t2);
  }
  inline bool operator!=(const Lorentz_Ten3C& t1, const Lorentz_Ten3C& t2) {
    return !IsEqual(t1,t2);
  }
  inline bool operator!=(const Lorentz_Ten4C& t1, const Lorentz_Ten4C& t2) {
    return !IsEqual(t1,t2);
  }

  inline Lorentz_Ten2C conj(const Lorentz_Ten2C& t) {
      Complex x[4][4];
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          x[i][j] = std::conj(t.at(i,j));
      return Lorentz_Ten2C(x);
  }

  inline Lorentz_Ten3C conj(const Lorentz_Ten3C& t) {
      Complex x[4][4][4];
      Complex z = 0.;
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          for (unsigned short int k=0; k<4; ++k)
            x[i][j][k] = std::conj(t.at(i,j,k));
      return Lorentz_Ten3C(x);
  }

  inline Lorentz_Ten4C conj(const Lorentz_Ten4C& t) {
      Complex x[4][4][4][4];
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          for (unsigned short int k=0; k<4; ++k)
            for (unsigned short int l=0; l<4; ++l)
            x[i][j][k][l] = std::conj(t.at(i,j,k,l));
      return Lorentz_Ten4C(x);
  }


  // meaningful only for the transpostion of a pair of indizes (i<->j)
  inline Lorentz_Ten2C hermconj(const Lorentz_Ten2C& t) {
    return conj(t).Transpose();
  }

  inline Lorentz_Ten3C hermconj(const Lorentz_Ten3C& t, unsigned short int i, unsigned short int j) {
    return conj(t).Transpose(i,j);
  }

  inline Lorentz_Ten4C hermconj(const Lorentz_Ten4C& t, unsigned short int i, unsigned short int j) {
    return conj(t).Transpose(i,j);
  }





  /************************************************************************************/
  // define special tensors
  // all tensor defined with upper indices

  // metric tensor g^{mu,nu} = g_{mu,nu}, componentwise

  inline Lorentz_Ten2D MetricTensor() {
    return Lorentz_Ten2D( 1., 0., 0., 0.,
                          0.,-1., 0., 0.,
                          0., 0.,-1., 0.,
                          0., 0., 0.,-1.);
  }

  inline Lorentz_Ten4D EpsilonTensor() {
    return Lorentz_Ten4D(
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,-1., 0.,
        0., 0., 0., 0., 0., 0., 0.,-1., 0., 0., 0., 0., 0., 1., 0., 0.,
        0., 0., 0., 0., 0., 0., 1., 0., 0.,-1., 0., 0., 0., 0., 0., 0.,

        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,-1., 0., 0., 1., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,-1., 0., 0., 0.,
        0., 0.,-1., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,

        0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.,-1., 0., 0.,
        0., 0., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 1., 0., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,

        0., 0., 0., 0., 0., 0.,-1., 0., 0., 1., 0., 0., 0., 0., 0., 0.,
        0., 0., 1., 0., 0., 0., 0., 0.,-1., 0., 0., 0., 0., 0., 0., 0.,
        0.,-1., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.);
  }

  /***************************************************************************/
  // build tensor from vectors and other lower rank tensors

  template<typename Scal1, typename Scal2>                                // ok
  Lorentz_Ten2<PROMOTE(Scal1,Scal2)> 
  BuildTensor(const Vec4<Scal1>&, const Vec4<Scal2>&);

  template<typename Scal1, typename Scal2, typename Scal3>                // ok
  Lorentz_Ten3<PROMOTE(Scal1,PROMOTE(Scal2,Scal3))> 
  BuildTensor(const Vec4<Scal1>&, const Vec4<Scal2>&, const Vec4<Scal3>&);

  template<typename Scal1, typename Scal2>                                // ok
  Lorentz_Ten3<PROMOTE(Scal1,Scal2)> 
  BuildTensor(const Lorentz_Ten2<Scal1>&, const Vec4<Scal2>&);

  template<typename Scal1, typename Scal2, typename Scal3, typename Scal4>// ok
  Lorentz_Ten4<PROMOTE(Scal1,PROMOTE(Scal2,PROMOTE(Scal3,Scal4)))> 
  BuildTensor(const Vec4<Scal1>&, const Vec4<Scal2>&,
              const Vec4<Scal3>&, const Vec4<Scal4>&);

  template<typename Scal1, typename Scal2>                                // ok
  Lorentz_Ten4<PROMOTE(Scal1,Scal2)> 
  BuildTensor(const Lorentz_Ten2<Scal1>&, const Lorentz_Ten2<Scal2>&);

  template<typename Scal1, typename Scal2>                                // ok
  Lorentz_Ten4<PROMOTE(Scal1,Scal2)> 
  BuildTensor(const Lorentz_Ten3<Scal1>&, const Vec4<Scal2>&);


  /***************************************************************************/
  // tensor contractions

  // epsilon tensor contractions  

  template<typename Scalar>                                               // ok
  Lorentz_Ten3<Scalar>
  EpsilonTensorContraction(const Vec4<Scalar>&);

  template<typename Scal1, typename Scal2>                                // ok
  Lorentz_Ten2<PROMOTE(Scal1,Scal2)>
  EpsilonTensorContraction(const Vec4<Scal1>&, const Vec4<Scal2>&);

  template<typename Scal1, typename Scal2, typename Scal3>                // ok
  Vec4<PROMOTE(Scal1,PROMOTE(Scal2,Scal3))>
  EpsilonTensorContraction(const Vec4<Scal1>&, const Vec4<Scal2>&,
                           const Vec4<Scal3>&);

  template<typename Scal1, typename Scal2, typename Scal3, typename Scal4>// ok
  PROMOTE(Scal1,PROMOTE(Scal2,PROMOTE(Scal3,Scal4)))
  EpsilonTensorContraction(const Vec4<Scal1>&, const Vec4<Scal2>&,
                           const Vec4<Scal3>&, const Vec4<Scal4>&);

  template<typename Scalar>                                               // ok
  Lorentz_Ten2<Scalar>
  EpsilonTensorContraction(const Lorentz_Ten2<Scalar>&, int, int);

  template<typename Scal1, typename Scal2>                                // ok
  Vec4<PROMOTE(Scal1,Scal2)>
  EpsilonTensorContraction(const Lorentz_Ten2<Scal1>&, int, int,
                           const Vec4<Scal2>&, int);

  template<typename Scal1, typename Scal2>                                // ok
  PROMOTE(Scal1,Scal2)
  EpsilonTensorContraction(const Lorentz_Ten2<Scal1>&, int, int,
                           const Lorentz_Ten2<Scal2>&, int, int);

  template<typename Scalar>                                               // ok
  Vec4<Scalar>
  EpsilonTensorContraction(const Lorentz_Ten3<Scalar>&, int, int, int);

  template<typename Scal1, typename Scal2>                                // ok
  PROMOTE(Scal1,Scal2)
  EpsilonTensorContraction(const Lorentz_Ten3<Scal1>&, int, int, int,
                           const Vec4<Scal2>&, int);

  template<typename Scalar>                                               // ok
  Scalar
  EpsilonTensorContraction(const Lorentz_Ten4<Scalar>&, int, int, int, int);


  // others

  // tensor with itself
  template<typename Scalar>                                               // ok
  Scalar
  Contraction(const Lorentz_Ten2<Scalar>&);

  template<typename Scalar>                                               // ok
  Vec4<Scalar>
  Contraction(const Lorentz_Ten3<Scalar>&, int, int);

  template<typename Scalar>                                               // ok
  Lorentz_Ten2<Scalar>
  Contraction(const Lorentz_Ten4<Scalar>&, int, int);

  template<typename Scalar>                                               // ok
  Scalar
  Contraction(const Lorentz_Ten4<Scalar>&, int, int, int, int);


  // tensor with vector
  template<typename Scal1, typename Scal2>                                // ok
  Vec4<PROMOTE(Scal1,Scal2)>
  Contraction(const Lorentz_Ten2<Scal1>&, int, const Vec4<Scal2>&);

  template<typename Scal1, typename Scal2, typename Scal3>                // ok
  PROMOTE(PROMOTE(Scal1,Scal2),Scal3)
  Contraction(const Lorentz_Ten2<Scal1>&, const Vec4<Scal2>&,
              const Vec4<Scal3>&);

  template<typename Scal1, typename Scal2>                                // ok
  Lorentz_Ten2<PROMOTE(Scal1,Scal2)>
  Contraction(const Lorentz_Ten3<Scal1>&, int, const Vec4<Scal2>&);

  template<typename Scal1, typename Scal2>                                // ok
  Lorentz_Ten3<PROMOTE(Scal1,Scal2)>
  Contraction(const Lorentz_Ten4<Scal1>&, int, const Vec4<Scal2>&);


  // tensor with other tensor
  template<typename Scal1, typename Scal2>                                // ok
  PROMOTE(Scal1,Scal2)
  Contraction(const Lorentz_Ten2<Scal1>&, int, int,
              const Lorentz_Ten2<Scal2>&, int, int);

  template<typename Scal1, typename Scal2>                                // ok
  Lorentz_Ten2<PROMOTE(Scal1,Scal2)>
  Contraction(const Lorentz_Ten2<Scal1>&, int, const Lorentz_Ten2<Scal2>&, int);

  template<typename Scal1, typename Scal2>                                // ok
  Vec4<PROMOTE(Scal1,Scal2)>
  Contraction(const Lorentz_Ten3<Scal1>&, int, int,
              const Lorentz_Ten2<Scal2>&, int, int);

  template<typename Scal1, typename Scal2>                                // ok
  Lorentz_Ten3<PROMOTE(Scal1,Scal2)>
  Contraction(const Lorentz_Ten3<Scal1>&, int,
              const Lorentz_Ten2<Scal2>&, int);

  template<typename Scal1, typename Scal2>
  PROMOTE(Scal1,Scal2)
  Contraction(const Lorentz_Ten3<Scal1>&, int, int, int,
              const Lorentz_Ten3<Scal2>&, int, int, int);

  template<typename Scal1, typename Scal2>
  Lorentz_Ten2<PROMOTE(Scal1,Scal2)>
  Contraction(const Lorentz_Ten3<Scal1>&, int, int,
              const Lorentz_Ten3<Scal2>&, int, int);

  template<typename Scal1, typename Scal2>
  Lorentz_Ten4<PROMOTE(Scal1,Scal2)>
  Contraction(const Lorentz_Ten3<Scal1>&, int,
              const Lorentz_Ten3<Scal2>&, int);

  template<typename Scal1, typename Scal2>
  Lorentz_Ten2<PROMOTE(Scal1,Scal2)>
  Contraction(const Lorentz_Ten4<Scal1>&, int, int,
              const Lorentz_Ten2<Scal2>&, int, int);

  template<typename Scal1, typename Scal2>
  Lorentz_Ten4<PROMOTE(Scal1,Scal2)>
  Contraction(const Lorentz_Ten4<Scal1>&, int,
              const Lorentz_Ten2<Scal2>&, int);

  template<typename Scal1, typename Scal2>                                //(tbc)
  Vec4<PROMOTE(Scal1,Scal2)>
  Contraction(const Lorentz_Ten4<Scal1>&, int, int, int,
              const Lorentz_Ten3<Scal2>&, int, int, int);

  template<typename Scal1, typename Scal2>                                //(tbc)
  Lorentz_Ten3<PROMOTE(Scal1,Scal2)>
  Contraction(const Lorentz_Ten4<Scal1>&, int, int,
              const Lorentz_Ten3<Scal2>&, int, int);

  template<typename Scal1, typename Scal2>                                // tbc
  PROMOTE(Scal1,Scal2)
  Contraction(const Lorentz_Ten4<Scal1>&, int, int, int, int,
              const Lorentz_Ten4<Scal2>&, int, int, int, int);

  template<typename Scal1, typename Scal2>                                // tbc
  Lorentz_Ten2<PROMOTE(Scal1,Scal2)>
  Contraction(const Lorentz_Ten4<Scal1>&, int, int, int,
              const Lorentz_Ten4<Scal2>&, int, int, int);

  template<typename Scal1, typename Scal2>                                //(tbc)
  Lorentz_Ten4<PROMOTE(Scal1,Scal2)>
  Contraction(const Lorentz_Ten4<Scal1>&, int, int,
              const Lorentz_Ten4<Scal2>&, int, int);

}

//*********************************************************************************
//*********************************************************************************
// tensor algebra implementation
//*********************************************************************************

using namespace std;
using namespace ATOOLS;

// build tensors

#include "ATOOLS/Math/Tensor_Build.H"

// tensor contractions

#include "ATOOLS/Math/Tensor_Contractions_Epsilon.H"

// others

#include "ATOOLS/Math/Tensor_Contractions.H"

#endif
