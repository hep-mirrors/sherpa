#ifndef METOOLS_Explicit_C_Tensor_H
#define METOOLS_Explicit_C_Tensor_H

#include "METOOLS/Currents/C_Vector.H"
#include "ATOOLS/Org/STL_Tools.H"
 
namespace METOOLS {

  template <class Scalar>
  class CAsT4: public CObject {
  public:

    typedef std::complex<Scalar> SComplex;

  private:

    SComplex m_x[6];

    static double s_accu;

    static ATOOLS::AutoDelete_Vector<CAsT4> s_objects;

    template <class _Scalar> friend std::ostream &
    operator<<(std::ostream &s,const CAsT4<_Scalar> &ten);

  public:

    static CAsT4 *New();
    static CAsT4 *New(const CAsT4 &s);

    CObject* Copy() const;

    void Delete();

    bool IsZero() const;

    inline CAsT4()
    { 
      m_x[0]=m_x[1]=m_x[2]=m_x[3]=m_x[4]=m_x[5]=0.0; 
      m_c[0]=m_c[1]=0; m_h=m_s=0;
    }
    inline CAsT4(const CAsT4 &t)
    { 
      m_x[0]=t[0]; m_x[1]=t[1]; m_x[2]=t[2]; 
      m_x[3]=t[3]; m_x[4]=t[4]; m_x[5]=t[5];
      m_c[0]=t(0); m_c[1]=t(1); m_h=t.m_h; m_s=t.m_s;
    }
    inline CAsT4(const Scalar &x0, const Scalar &x1, 
		 const Scalar &x2, const Scalar &x3,
		 const Scalar &x4, const Scalar &x5,
		 const int c1=0,const int c2=0,
		 const size_t &h=0,const size_t &s=0)
    { 
      m_x[0]=x0; m_x[1]=x1; m_x[2]=x2;
      m_x[3]=x3; m_x[4]=x4; m_x[5]=x5; 
      m_c[0]=c1; m_c[1]=c2; m_h=h; m_s=s;
    }
    inline CAsT4(const int c1,const int c2,
		 const size_t &h=0,const size_t &s=0)
    { 
      m_x[0]=m_x[1]=m_x[2]=m_x[3]=m_x[4]=m_x[5]=0.0; 
      m_c[0]=c1; m_c[1]=c2; m_h=h; m_s=s;
    }
    inline CAsT4(const SComplex &x0, const SComplex &x1, 
		 const SComplex &x2, const SComplex &x3,
		 const SComplex &x4, const SComplex &x5,
		 const int c1=-1,const int c2=-1,
		 const size_t &h=0,const size_t &s=0)
    { 
      m_x[0]=x0; m_x[1]=x1; m_x[2]=x2;
      m_x[3]=x3; m_x[4]=x4; m_x[5]=x5; 
      m_c[0]=c1; m_c[1]=c2; m_h=h; m_s=s;
    }

    CAsT4(const CVec4<Scalar> &v1,const CVec4<Scalar> &v2);
    CAsT4(const CAsT4 &v,const SComplex &c);
    CAsT4(const CAsT4 &v,const Scalar &c);

    void Add(const CObject *c);
    void Divide(const double &d);
    void Multiply(const Complex &c);
    void Invert();

    inline SComplex &operator[](const int i) { return m_x[i]; }

    inline const SComplex &operator[](const int i) const { return m_x[i]; }

    inline CAsT4 operator+(const CAsT4 &v) const  
    { 
      return CAsT4(m_x[0]+v[0],m_x[1]+v[1],m_x[2]+v[2],
		   m_x[3]+v[3],m_x[4]+v[4],m_x[5]+v[5],
		   m_c[0],m_c[1],m_h,m_s); 
    }
    inline CAsT4 operator-(const CAsT4 &v) const
    { 
      return CAsT4(m_x[0]-v[0],m_x[1]-v[1],m_x[2]-v[2],
		   m_x[3]-v[3],m_x[4]-v[4],m_x[5]-v[5],
		   m_c[0],m_c[1],m_h,m_s); 
    }
    inline CAsT4 operator-() const
    { 
      return CAsT4(-m_x[0],-m_x[1],-m_x[2],-m_x[3],-m_x[4],-m_x[5],
		   m_c[0],m_c[1],m_h,m_s); 
    }

    inline CAsT4& operator+=(const CAsT4 &v) 
    {
      m_x[0]+=v[0]; m_x[1]+=v[1]; m_x[2]+=v[2];
      m_x[3]+=v[3]; m_x[4]+=v[4]; m_x[5]+=v[5];
      return *this;
    }
    inline CAsT4& operator-=(const CAsT4 &v) 
    {
      m_x[0]-=v[0]; m_x[1]-=v[1]; m_x[2]-=v[2];
      m_x[3]-=v[3]; m_x[4]-=v[4]; m_x[5]-=v[5];
      return *this;
    }

    CAsT4& operator*=(const SComplex &c);
  
    inline CAsT4 Conj() const 
    {
      return CAsT4(std::conj(m_x[0]),std::conj(m_x[1]),std::conj(m_x[2]),
		   std::conj(m_x[3]),std::conj(m_x[4]),std::conj(m_x[5]),
		   m_c[0],m_c[1],m_h,m_s);
    }

    bool Nan() const;

    static void ResetAccu();

    inline static void   SetAccu(const double &accu) { s_accu=accu;   }
    inline static double Accu()                      { return s_accu; }

  };// end of class CAsT4

  template <class Scalar> inline CAsT4<Scalar> 
  operator*(const Scalar &c,const CAsT4<Scalar> &t)
  { return CAsT4<Scalar>(t,c); }
  template <class Scalar> inline CAsT4<Scalar> 
  operator*(const CAsT4<Scalar> &t,const Scalar &c)
  { return CAsT4<Scalar>(t,c); }
  template <class Scalar> inline CAsT4<Scalar> 
  operator*(const std::complex<Scalar> &c,const CAsT4<Scalar> &t)
  { return CAsT4<Scalar>(t,c); }
  template <class Scalar> inline CAsT4<Scalar> 
  operator*(const CAsT4<Scalar> &t,const std::complex<Scalar> &c)
  { return CAsT4<Scalar>(t,c); }
  template <class Scalar> inline CAsT4<Scalar> 
  operator/(const CAsT4<Scalar> &t,const std::complex<Scalar> &c)
  { return CAsT4<Scalar>(t,1.0/c); }
  template <class Scalar> CVec4<Scalar> 
  operator*(const CVec4<Scalar> &v,const CAsT4<Scalar> &t);
  template <class Scalar> inline CVec4<Scalar>
  operator*(const CAsT4<Scalar> &t,const CVec4<Scalar> &v)
  { return -(v*t); }

  template <class Scalar> 
  std::ostream &operator<<(std::ostream &s,const CAsT4<Scalar> &ten);

}// end of namespace ATOOLS

#define DCAsT4D METOOLS::CAsT4<double>
#define QCAsT4D METOOLS::CAsT4<long double>

#endif
