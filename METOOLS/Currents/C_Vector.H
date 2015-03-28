#ifndef METOOLS_Explicit_C_Vector_H
#define METOOLS_Explicit_C_Vector_H

#include "METOOLS/Explicit/C_Object.H"
#include "ATOOLS/Math/Vec4.H"
#include "ATOOLS/Org/STL_Tools.H"

#include <vector>
 
namespace METOOLS {

  template <class Scalar>
  class CVec4: public CObject {
  public:

    typedef std::complex<Scalar> SComplex;

  private:

    SComplex m_x[4];

    static double s_accu;

    static ATOOLS::AutoDelete_Vector<CVec4> s_objects;

    template <class _Scalar> friend std::ostream &
    operator<<(std::ostream &s,const CVec4<_Scalar> &vec);

  public:

    static CVec4 *New();
    static CVec4 *New(const CVec4 &s);
    static CVec4 *New(const Scalar &x0, const Scalar &x1, 
		      const Scalar &x2, const Scalar &x3,
		      const int c1=0,const int c2=0,
		      const size_t &hm=0,
		      const size_t &hp=0);

    CObject* Copy() const;

    void Delete();

    bool IsZero() const;

    inline CVec4()
    { 
      m_x[0]=m_x[1]=m_x[2]=m_x[3]=0.0; 
      m_c[0]=m_c[1]=0; 
      m_h=m_s=0; 
    }
    inline CVec4(const CVec4 &v)
    { 
      m_x[0]=v[0]; m_x[1]=v[1]; m_x[2]=v[2]; m_x[3]=v[3];
      m_c[0]=v(0); m_c[1]=v(1);
      m_h=v.m_h; m_s=v.m_s;
    }
    inline CVec4(const Scalar &x0, const Scalar &x1, 
		 const Scalar &x2, const Scalar &x3,
		 const int c1=0,const int c2=0,
		 const size_t &h=0,const size_t &s=0)
    { 
      m_x[0]=x0; m_x[1]=x1; m_x[2]=x2; m_x[3]=x3; 
      m_c[0]=c1; m_c[1]=c2; m_h=h; m_s=s;
    }
    inline CVec4(const int c1,const int c2,
		 const size_t &h=0,const size_t &s=0)
    { 
      m_x[0]=m_x[1]=m_x[2]=m_x[3]=0.0; 
      m_c[0]=c1; m_c[1]=c2; m_h=h; m_s=s;
    }
    inline CVec4(const SComplex &x0, const SComplex &x1, 
		 const SComplex &x2, const SComplex &x3,
		 const int c1=0,const int c2=0,
		 const size_t &h=0,const size_t &s=0)
    { 
      m_x[0]=x0; m_x[1]=x1; m_x[2]=x2; m_x[3]=x3; 
      m_c[0]=c1; m_c[1]=c2; m_h=h; m_s=s;
    }
    inline CVec4(const ATOOLS::Vec4<Scalar> &v,
		 const int c1=0,const int c2=0,
		 const size_t &h=0,const size_t &s=0)
    { 
      m_x[0]=v[0]; m_x[1]=v[1]; m_x[2]=v[2]; m_x[3]=v[3]; 
      m_c[0]=c1; m_c[1]=c2; m_h=h; m_s=s;
    }
    inline CVec4(const CVec4 &v,const Scalar &c)
    {
      m_x[0]=v[0]*c; m_x[1]=v[1]*c; m_x[2]=v[2]*c; m_x[3]=v[3]*c;
      m_c[0]=v(0); m_c[1]=v(1); m_h=v.m_h; m_s=v.m_s;
    }
    inline CVec4(const CVec4 &v,const SComplex &c)
    {
      m_x[0]=v[0]*c; m_x[1]=v[1]*c; m_x[2]=v[2]*c; m_x[3]=v[3]*c;
      m_c[0]=v(0); m_c[1]=v(1); m_h=v.m_h; m_s=v.m_s;
    }

    void Add(const CObject *c);
    void Divide(const double &d);
    void Multiply(const Complex &c);
    void Invert();

    inline SComplex &operator[](const int i) { return m_x[i]; }

    inline const SComplex &operator[](const int i) const { return m_x[i]; }

    inline CVec4 operator+(const CVec4 &v) const  
    { 
      return CVec4(m_x[0]+v[0],m_x[1]+v[1],m_x[2]+v[2],m_x[3]+v[3],
		   m_c[0],m_c[1],m_h,m_s); 
    }
    inline CVec4 operator-(const CVec4 &v) const
    { 
      return CVec4(m_x[0]-v[0],m_x[1]-v[1],m_x[2]-v[2],m_x[3]-v[3],
		   m_c[0],m_c[1],m_h,m_s); 
    }
    inline CVec4 operator-() const
    { 
      return CVec4(-m_x[0],-m_x[1],-m_x[2],-m_x[3],
		   m_c[0],m_c[1],m_h,m_s); 
    }

    inline CVec4& operator+=(const CVec4 &v) 
    {
      m_x[0]+=v[0]; m_x[1]+=v[1]; m_x[2]+=v[2]; m_x[3]+=v[3];
      return *this;
    }
    inline CVec4& operator-=(const CVec4 &v) 
    {
      m_x[0]-=v[0]; m_x[1]-=v[1]; m_x[2]-=v[2]; m_x[3]-=v[3];
      return *this;
    }
    inline CVec4& operator*=(const SComplex &c) 
    {
      m_x[0]*=c; m_x[1]*=c; m_x[2]*=c; m_x[3]*=c;
      return *this;
    }
  
    inline CVec4 Conj() const 
    {
      return CVec4(std::conj(m_x[0]),std::conj(m_x[1]),
		   std::conj(m_x[2]),std::conj(m_x[3]),
		   m_c[0],m_c[1],m_h,m_s);
    }
    inline SComplex Abs2() const 
    {
      return m_x[0]*m_x[0]-m_x[1]*m_x[1]-m_x[2]*m_x[2]-m_x[3]*m_x[3];
    }
    inline SComplex Abs() const 
    { 
      return sqrt(Abs2()); 
    }

    bool Nan() const;

    static void ResetAccu();

    inline static void   SetAccu(const double &accu) { s_accu=accu;   }
    inline static double Accu()                      { return s_accu; }

  };// end of class CVec4

  template <class Scalar> inline CVec4<Scalar> 
  operator*(const Scalar &c,const CVec4<Scalar> &v)
  { return CVec4<Scalar>(v,c); }
  template <class Scalar> inline CVec4<Scalar> 
  operator*(const CVec4<Scalar> &v,const Scalar &c)
  { return CVec4<Scalar>(v,c); }
  template <class Scalar> inline CVec4<Scalar> 
  operator*(const std::complex<Scalar> &c,const CVec4<Scalar> &v)
  { return CVec4<Scalar>(v,c); }
  template <class Scalar> inline CVec4<Scalar> 
  operator*(const CVec4<Scalar> &v,const std::complex<Scalar> &c)
  { return CVec4<Scalar>(v,c); }
  template <class Scalar> inline CVec4<Scalar> 
  operator/(const CVec4<Scalar> &v,const std::complex<Scalar> &c)
  { return CVec4<Scalar>(v,Scalar(1.0)/c); }

  template <class Scalar> inline std::complex<Scalar> 
  operator*(const CVec4<Scalar> &v1,const ATOOLS::Vec4<Scalar> &v2)
  { return v1[0]*v2[0]-v1[1]*v2[1]-v1[2]*v2[2]-v1[3]*v2[3]; }
  template <class Scalar> inline std::complex<Scalar>
  operator*(const ATOOLS::Vec4<Scalar> &v1,const CVec4<Scalar> &v2)
  { return v1[0]*v2[0]-v1[1]*v2[1]-v1[2]*v2[2]-v1[3]*v2[3]; }
  template <class Scalar> inline std::complex<Scalar>
  operator*(const CVec4<Scalar> &v1,const CVec4<Scalar> &v2)
  { return v1[0]*v2[0]-v1[1]*v2[1]-v1[2]*v2[2]-v1[3]*v2[3]; }

  template <class Scalar>
  std::ostream &operator<<(std::ostream &s,const CVec4<Scalar> &vec);

}// end of namespace ATOOLS

#define DCVec4D METOOLS::CVec4<double>
#define QCVec4D METOOLS::CVec4<long double>

#endif
