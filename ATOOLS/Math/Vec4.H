#ifndef ATOOLS_Math_Vec4_H
#define ATOOLS_Math_Vec4_H

#include "ATOOLS/Math/MathTools.H"

namespace ATOOLS {
  template<typename Scalar> class Vec4;
  template<typename Scalar> class Vec3;

  template<typename Scalar>
  class Vec4 {
    Scalar m_x[4];
    template <typename Scalar2> friend class Vec4;
  public:
    inline Vec4() {
      m_x[0]=m_x[1]=m_x[2]=m_x[3]=Scalar(0.0);
    }
    template<typename Scalar2>
    inline Vec4(const Vec4<Scalar2> & vec) {
      m_x[0] = vec[0]; m_x[1] = vec[1];
      m_x[2] = vec[2]; m_x[3] = vec[3];
    }
    inline Vec4(const Scalar& x0, const Scalar& x1,
                const Scalar& x2, const Scalar& x3) {
      m_x[0] = x0; m_x[1] = x1; m_x[2] = x2; m_x[3] = x3;
    }
    inline Vec4(const Scalar& E, const Vec3<Scalar>& xn) {
      m_x[0] = E;
      m_x[1] = xn[1];
      m_x[2] = xn[2];
      m_x[3] = xn[3];
    }

    inline Scalar& operator[] (int i) {
      return m_x[i];
    }

    inline const Scalar operator[] (int i) const {
      return m_x[i];
    }

    inline Vec4<Scalar>& operator+=(const Vec4<Scalar>& v) {
      m_x[0] += v[0];
      m_x[1] += v[1];
      m_x[2] += v[2];
      m_x[3] += v[3];
      return *this;
    }

    inline Vec4<Scalar>& operator-=(const Vec4<Scalar>& v) {
      m_x[0] -= v[0];
      m_x[1] -= v[1];
      m_x[2] -= v[2];
      m_x[3] -= v[3];
      return *this;
    }

    inline Vec4<Scalar>& operator*= (const Scalar& scal) {
      m_x[0] *= scal;
      m_x[1] *= scal;
      m_x[2] *= scal;
      m_x[3] *= scal;
      return *this;
    }

    inline Vec4<Scalar> operator-() const {
      return Vec4<Scalar>(-m_x[0],-m_x[1],-m_x[2],-m_x[3]);
    }

    template<typename Scalar2> inline PROMOTE(Scalar,Scalar2)
    operator*(const Vec4<Scalar2>& v2) const {
      return m_x[0]*v2[0]-m_x[1]*v2[1]-m_x[2]*v2[2]-m_x[3]*v2[3];
    }

    template<typename Scalar2> inline Vec4<PROMOTE(Scalar,Scalar2)>
    operator+ (const Vec4<Scalar2>& v2) const {
      return Vec4<PROMOTE(Scalar,Scalar2)>
        (m_x[0]+v2[0],m_x[1]+v2[1],m_x[2]+v2[2],m_x[3]+v2[3]);
    }

    template<typename Scalar2> inline Vec4<PROMOTE(Scalar,Scalar2)>
    operator- (const Vec4<Scalar2>& v2) const {
      return Vec4<PROMOTE(Scalar,Scalar2)>
        (m_x[0]-v2[0],m_x[1]-v2[1],m_x[2]-v2[2],m_x[3]-v2[3]);
    }

    template<typename Scalar2> inline Vec4<PROMOTE(Scalar,Scalar2)>
    operator* (const Scalar2& s) const {
      return Vec4<PROMOTE(Scalar,Scalar2)>
        (s*m_x[0],s*m_x[1],s*m_x[2],s*m_x[3]);
    }

    template<typename Scalar2> inline Vec4<PROMOTE(Scalar,Scalar2)>
    operator/ (const Scalar2& scal) const {
      return (*this)*(1./scal);
    }


    inline bool Nan() const {
      for(short unsigned int i(0);i<4;++i)
        if (ATOOLS::IsNan<Scalar>(m_x[i])) return true;
      return false;
    }
    inline bool IsZero() const {
      for(short unsigned int i(0);i<4;++i) 
        if (!ATOOLS::IsZero<Scalar>(m_x[i])) return false;
      return true;
    }


    Scalar Abs2() const { return m_x[0]*m_x[0]-PSpat2(); }      
    Scalar Abs() const { return sqrt(Abs2()); }
    double Mass() const {
      return sqrt(ATOOLS::Abs<Scalar>(Abs2()));
    }
    Scalar P() const { return PSpat(); }

    Scalar PPerp2() const { return m_x[1]*m_x[1]+m_x[2]*m_x[2]; }
    Scalar PPerp() const { return sqrt(PPerp2()); }
    Scalar MPerp2() const { return m_x[0]*m_x[0]-m_x[3]*m_x[3]; }
    Scalar MPerp() const { return sqrt(MPerp2()); }
    Scalar MPerp2(const Vec4<Scalar> &ref) const {
      return Abs2()+PPerp2(ref);
    }
    Scalar MPerp(const Vec4<Scalar> &ref) const {
      return sqrt(MPerp2(ref));
    }
    Scalar EPerp2() const {
      return m_x[0]*m_x[0]*PPerp2()/PSpat2();
    }
    Scalar EPerp() const { return sqrt(EPerp2()); }
    Scalar PPlus() const { return m_x[0]+m_x[3]; }
    Scalar PMinus() const { return m_x[0]-m_x[3]; }

    Scalar PSpat2() const {
      return m_x[1]*m_x[1]+m_x[2]*m_x[2]+m_x[3]*m_x[3];
    }
    Scalar PSpat() const { return sqrt(PSpat2()); }
    Scalar Y() const { return 0.5*log(PPlus()/PMinus()); }

    inline Vec4<Scalar> Perp() const { 
      return Vec4<Scalar>(0.,m_x[1],m_x[2],0.);
    }
    inline Vec4<Scalar> Long() const { 
      return Vec4<Scalar>(m_x[0],0.,0.,m_x[3]);
    }
    inline Vec4<Scalar> Plus() const {
      Scalar pplus=0.5*PPlus();
      return Vec4<Scalar>(pplus,0.0,0.0,pplus);
    }
    inline Vec4<Scalar> Minus() const {
      Scalar pminus=0.5*PMinus();
      return Vec4<Scalar>(pminus,0.0,0.0,-pminus);
    }
    
    Scalar PPerp2(const Vec4<Scalar>& ref) const {
      Scalar ppref = PPerp(ref);
      return ppref*ppref;
    }
    Scalar PPerp(const Vec4<Scalar>& ref) const {
      Vec3<Scalar> perp=1./ATOOLS::Max(Vec3<Scalar>(ref).Abs(),1.e-12)*
        Vec3<Scalar>(ref);
      perp=Vec3<Scalar>(*this)-perp*(perp*Vec3<Scalar>(*this));
      return perp.Abs();
    }

    // some functions which only make sense for doubles
    // and thus have to be specialised for the type
    Scalar CosPhi() const;
    Scalar SinPhi() const;
    Scalar Phi() const;
    Scalar CosTheta() const;
    Scalar SinTheta() const;
    Scalar Theta() const;
    Scalar Eta() const;
    Scalar CosTheta(const Vec4<Scalar>& ref) const;
    Scalar Theta(const Vec4<Scalar>& ref) const;
    Scalar Eta(const Vec4<Scalar>& ref) const;
    Scalar CosDPhi(const Vec4<Scalar>& ref) const;
    Scalar DPhi(const Vec4<Scalar>& ref) const;
    Scalar DEta(const Vec4<Scalar>& ref) const;
    Scalar DR(const Vec4<Scalar>& ref) const;

    // standard vectors:
    const static Vec4<Scalar> XVEC;
    const static Vec4<Scalar> YVEC;
    const static Vec4<Scalar> ZVEC;
  };

  template<typename Scalar,typename Scal2> inline Vec4<PROMOTE(Scalar,Scal2)>
  operator* (const Scal2& s, const Vec4<Scalar>& v) {
    return Vec4<PROMOTE(Scalar,Scal2)>(v*s);
  }

  template<typename S1,typename S2,typename S3> inline
  Vec4<PROMOTE(S1,PROMOTE(S2,S3))>
  cross(const Vec4<S1>&q, const Vec4<S2>&r, const Vec4<S3>&s)
  {
    PROMOTE(S1,PROMOTE(S2,S3)) r0s1mr1s0(r[0]*s[1]-r[1]*s[0]), 
      r0s2mr2s0(r[0]*s[2]-r[2]*s[0]), r0s3mr3s0(r[0]*s[3]-r[3]*s[0]),
      r1s2mr2s1(r[1]*s[2]-r[2]*s[1]), r1s3mr3s1(r[1]*s[3]-r[3]*s[1]),
      r2s3mr3s2(r[2]*s[3]-r[3]*s[2]);
    return Vec4<PROMOTE(S1,PROMOTE(S2,S3))>
      (-q[1]*r2s3mr3s2+q[2]*r1s3mr3s1-q[3]*r1s2mr2s1,
       -q[0]*r2s3mr3s2+q[2]*r0s3mr3s0-q[3]*r0s2mr2s0,
       +q[0]*r1s3mr3s1-q[1]*r0s3mr3s0+q[3]*r0s1mr1s0,
       -q[0]*r1s2mr2s1+q[1]*r0s2mr2s0-q[2]*r0s1mr1s0);
  }

  template<typename Scalar>
  std::ostream& operator<<(std::ostream& s, const Vec4<Scalar>& vec) {
    return s<<'('<<vec[0]<<','<<vec[1]<<','<<vec[2]<<','<<vec[3]<<')';
  }
}

/*!
 \file
 \brief   contains class Vec4
*/

/*!
 \class Vec4<Scalar>
 \brief implementation of a 4 dimensional Minkowski vector and its algebraic
 operations

 This class can be used as Minkowski vector with arbitrary scalar types.
 All necessary operations, e.g. addition, scalar product, cross product, 
 etc. are available.
 If you want to use this vector for a new scalar type, you might have to
 implement specific functions for this type in MathTools.H.
 "Fallback" types for operations with two vectors of different types, e. g.
 "complex and double become complex" have to be specified in MathTools.H as
 well, see the DECLARE_PROMOTE macro.
*/

/*!
 \fn Vec4::Vec4()
 \brief Standard Constructor
*/

/*!
 \fn Vec4::Vec4(Scalar x0,Scalar x1, Scalar x2, Scalar x3){
 \brief Special Constructor, templated in Scalar, taking 4 single components.
*/

/*!
 \fn Vec4::Vec4(const Scalar E, const Vec3D& v)
 \brief Special Constructor taking the space part and a Energie
*/

/*!
 \fn inline const Scalar Vec4::Abs2() const
 \brief returns \f$ x_0^2 - (x_1^2 + x_2^2 + x_3^2) \f$.
*/

/*!
 \fn   inline Scalar& Vec4::operator[] (int i)
 \brief returns \f$x_i\f$. (May be manipulated.)
*/

/*!
 \fn   inline const Scalar Vec4::operator[] (int i) const
 \brief returns \f$x_i\f$.
*/

#endif
