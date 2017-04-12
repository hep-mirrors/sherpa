#ifndef ATOOLS_Math_Vector_H
#define ATOOLS_Math_Vector_H

#include <iostream>
#include <vector>
#include "ATOOLS/Math/Vec3.H"
#include "ATOOLS/Math/Vec4.H"
#include "ATOOLS/Math/MyComplex.H"

namespace ATOOLS {

  /*
   * declarations and specialisations for 3 and 4 dimensional vectors of doubles
   */
  typedef Vec3<double> Vec3D;
  typedef Vec4<double> Vec4D;
  typedef std::vector<Vec4D> Vec4D_Vector;

  template<> double Vec4D::CosPhi() const;
  template<> double Vec4D::SinPhi() const;
  template<> double Vec4D::Phi() const;
  template<> double Vec4D::CosTheta() const;
  template<> double Vec4D::SinTheta() const;
  template<> double Vec4D::Theta() const;
  template<> double Vec4D::Eta() const;
  template<> double Vec4D::CosTheta(const Vec4D& ref) const;
  template<> double Vec4D::Theta(const Vec4D& ref) const;
  template<> double Vec4D::Eta(const Vec4D& ref) const;
  template<> double Vec4D::CosDPhi(const Vec4D& ref) const;
  template<> double Vec4D::DPhi(const Vec4D& ref) const;
  template<> double Vec4D::DEta(const Vec4D& ref) const;
  template<> double Vec4D::DY(const Vec4D& ref) const;
  template<> double Vec4D::DR(const Vec4D& ref) const;

  std::istream& operator>>(std::istream& s,Vec4D& vec);
  std::istream& operator>>(std::istream& s,Vec3D& vec);

  bool IsEqual(const Vec4D& v1, const Vec4D& v2, const double crit=1.0e-12);
  inline bool operator==(const Vec4D& v1, const Vec4D& v2) {
    return IsEqual(v1,v2);
  }
  inline bool operator!=(const Vec4D& v1, const Vec4D& v2) {
    return !IsEqual(v1,v2);
  }
  bool IsEqual(const Vec3D& v1, const Vec3D& v2, const double crit=1.0e-12);
  inline bool operator==(const Vec3D& v1, const Vec3D& v2) {
    return IsEqual(v1,v2);
  }
  inline bool operator!=(const Vec3D& v1, const Vec3D& v2) {
    return !IsEqual(v1,v2);
  }


  /*
   * declarations and specialisations for 3 and 4 dimensional vectors of Complex
   */
  typedef Vec3<Complex> Vec3C;
  typedef Vec4<Complex> Vec4C;
  inline Vec4C conj(const Vec4C& v) {
    return Vec4C(std::conj(v[0]),std::conj(v[1]),
                 std::conj(v[2]),std::conj(v[3]));
  }

  // static members
  template<> const Vec4D Vec4D::XVEC;
  template<> const Vec4D Vec4D::YVEC;
  template<> const Vec4D Vec4D::ZVEC;
  template<> const Vec3D Vec3D::XVEC;
  template<> const Vec3D Vec3D::YVEC;
  template<> const Vec3D Vec3D::ZVEC;

}

/*!
 \file
 \brief   contains typedefs of Vec3D (3-dim. real vector), Vec3C (3-dim. complex
 vector), Vec4D (real Minkowski vector), Vec4C (complex Minkowski vector).

 Furthermore it contains specialisations of some member functions of the
 templated Vec3 and Vec4 classes.
*/

#endif
