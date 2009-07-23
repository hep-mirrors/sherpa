#ifndef ATOOLS_Math_Vec3_H
#define ATOOLS_Math_Vec3_H

#include "ATOOLS/Math/MathTools.H"
 
namespace ATOOLS {
  template<typename Scalar> class Vec4;
  template<typename Scalar> class Vec3;

  template<typename Scalar>
  class Vec3 {
    Scalar m_x[3];
  public:
    inline Vec3() {
      m_x[0]=m_x[1]=m_x[2]=Scalar(0.0);
    }
    inline Vec3(const Vec3<Scalar>& vec) {
      m_x[0]=vec.m_x[0];
      m_x[1]=vec.m_x[1];
      m_x[2]=vec.m_x[2];
    }
    inline Vec3(const Scalar& x, const Scalar& y, const Scalar& z) {
      m_x[0]=x;
      m_x[1]=y;
      m_x[2]=z;
    }
    inline Vec3(const Vec4<Scalar>& v) {
      m_x[0]=v[1];
      m_x[1]=v[2];
      m_x[2]=v[3];
    }

    bool Nan() const {
      for(short unsigned int i(0);i<3;++i)
        if (ATOOLS::IsNan<Scalar>(m_x[i])) return true;
      return false;
    }

    bool IsZero() const {
      for(short unsigned int i(0);i<3;++i)
        if (!ATOOLS::IsZero<Scalar>(m_x[i])) return false;
      return true;
    }

    inline const Scalar Abs() const { return sqrt(Sqr()); }
    inline const Scalar Sqr() const {
      return m_x[0]*m_x[0]+m_x[1]*m_x[1]+m_x[2]*m_x[2];
    }

    inline Scalar& operator[] (int i) { return m_x[--i]; }
    inline const Scalar operator[] (int i) const { return m_x[--i]; }

    template<typename Scalar2> inline PROMOTE(Scalar,Scalar2)
    operator*(const Vec3<Scalar2>& v2) const {
      return m_x[0]*v2[1]+m_x[1]*v2[2]+m_x[2]*v2[3];
    }

    template<typename Scalar2> inline Vec3<PROMOTE(Scalar,Scalar2)>
    operator+ (const Vec3<Scalar2>& v2) const {
      return Vec3<PROMOTE(Scalar,Scalar2)>
        (m_x[0]+v2[1],m_x[1]+v2[2],m_x[2]+v2[3]);
    }

    inline Vec3<Scalar>& operator+=(const Vec3<Scalar>& v) {
      m_x[0] += v[0];
      m_x[1] += v[1];
      m_x[2] += v[2];
      return *this;
    }

    template<typename Scalar2> inline Vec3<PROMOTE(Scalar,Scalar2)>
    operator- (const Vec3<Scalar2>& v2) const {
      return Vec3<PROMOTE(Scalar,Scalar2)>
        (m_x[0]-v2[1],m_x[1]-v2[2],m_x[2]-v2[3]);
    }

    template<typename Scalar2> inline Vec3<PROMOTE(Scalar,Scalar2)>
    operator* (const Scalar2& s) const {
      return Vec3<PROMOTE(Scalar,Scalar2)>
        (s*m_x[0],s*m_x[1],s*m_x[2]);
    }

    inline Vec3<Scalar>& operator*= (const Scalar& scal) {
      m_x[0] *= scal;
      m_x[1] *= scal;
      m_x[2] *= scal;
      return *this;
    }

    template<typename Scalar2> inline Vec3<PROMOTE(Scalar,Scalar2)>
    operator/ (const Scalar2& scal) const {
      return Vec3<PROMOTE(Scalar,Scalar2)>(m_x[0]/scal,m_x[1]/scal,m_x[2]/scal);
    }

    // standard vectors:
    const static Vec3<Scalar> XVEC;
    const static Vec3<Scalar> YVEC;
    const static Vec3<Scalar> ZVEC;
  };

  template<typename Scalar,typename Scal2> inline Vec3<PROMOTE(Scalar,Scal2)>
  operator* (const Scal2& s, const Vec3<Scalar>& v) {
    return v*s;
  }

  template<typename Scalar>
  inline Vec3<Scalar> cross(const Vec3<Scalar>& v1, const  Vec3<Scalar>& v2)
  {
    return Vec3<Scalar>(v1[2]*v2[3]-v1[3]*v2[2],
                        v1[3]*v2[1]-v1[1]*v2[3],
                        v1[1]*v2[2]-v1[2]*v2[1]);
  }

  template<typename Scalar>
  std::ostream& operator<<(std::ostream& s, const Vec3<Scalar>& vec) {
    return s<<'('<<vec[1]<<','<<vec[2]<<','<<vec[3]<<')';
  }

    
}

/*!
 \file
 \brief   contains class Vec3
*/

/*!
 \class Vec3<Scalar>
 \brief implementation of a 3 dimensional Euclidean vector and its algebraic
 operations

 This class can be used as a 3 dimensional vector with arbitrary scalar types.
 All necessary operations, e.g. addition, scalar product, cross product, 
 etc. are available.
 If you want to use this vector for a new scalar type, you might have to
 implement specific functions for this type in MathTools.H.
 "Fallback" types for operations with two vectors of different types, e. g.
 "complex and double become complex" have to be specified in MathTools.H as
 well, see the DECLARE_PROMOTE macro.
*/

/*!
 \fn Vec3::Vec3()
 \brief Standard Constructor
*/

/*!
 \fn Vec3::Vec3(double x, double y, double z){
 \brief Special Constructor taking 3 single components
*/

/*!
 \fn Vec3::Vec3(const Vec4D& v)
 \brief Special Constructor extracting the space part of a minkowski vector
*/

/*!
 \fn inline const double Vec3::Abs() const
 \brief returns \f$ \sqrt{x^2 + y^2 + z^2} \f$.
*/

/*!
 \fn inline const double Vec3::Sqr() const
 \brief returns \f$ x^2 + y^2 + z^2 \f$.
*/

/*!
 \fn   inline double& Vec3::operator[] (int i)
 \brief returns x,y,z for i=1,2,3 resp. (May be manipulated.)
*/

/*!
 \fn   inline const double Vec3::operator[] (int i) const
 \brief returns x,y,z for i=1,2,3 resp.
*/

#endif
