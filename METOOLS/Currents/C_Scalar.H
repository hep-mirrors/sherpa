#ifndef METOOLS_Explicit_C_Scalar_H
#define METOOLS_Explicit_C_Scalar_H

#include "METOOLS/Explicit/C_Object.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/STL_Tools.H"

#include <vector>
#include <iostream>
 
namespace METOOLS {

  template <class Scalar>
  class CScalar: public CObject {
  public:

    typedef std::complex<Scalar> SComplex;

  private:

    SComplex m_x;

    static double s_accu;

    static ATOOLS::AutoDelete_Vector<CScalar> s_objects;

    template <class _Scalar> friend std::ostream &
    operator<<(std::ostream &,const CScalar<_Scalar> &);

  public:

    static CScalar *New();
    static CScalar *New(const CScalar &s);

    CObject* Copy() const;

    void Delete();

    bool IsZero() const;

    inline CScalar(): m_x(Scalar(0.0))
    { 
      m_h=m_s=0; m_c[0]=m_c[1]=0;
    }
    inline CScalar(const CScalar &s): m_x(s.m_x)
    { 
      m_h=s.m_h; m_s=s.m_s; m_c[0]=s.m_c[0]; m_c[1]=s.m_c[1];
    }
    inline CScalar(const Scalar &x,
		   const int cr=0,const int ca=0,
		   const size_t &h=0,const size_t &s=0): m_x(x)
    { 
      m_h=h; m_s=s; m_c[0]=cr; m_c[1]=ca;
    }
    inline CScalar(const SComplex &x,
		   const int cr=0,const int ca=0,
		   const size_t &h=0,const size_t &s=0): m_x(x)
    { 
      m_h=h; m_s=s; m_c[0]=cr; m_c[1]=ca;
    }
    inline CScalar(const CScalar &s,const SComplex &c): m_x(s.m_x*c)
    {
      m_h=s.m_h; m_s=s.m_s; m_c[0]=s.m_c[0]; m_c[1]=s.m_c[1];
    }

    void Add(const CObject *c);
    void Divide(const double &d);
    void Invert();

    inline SComplex &operator[](const int i) { return m_x;    }

    inline const SComplex &operator[](const int i) const { return m_x;    }

    inline CScalar operator+(const CScalar &s) const  
    { 
      return CScalar(m_x+s.m_x,m_c[0],m_c[1],m_h,m_s); 
    }
    inline CScalar operator-(const CScalar &s) const
    { 
      return CScalar(m_x-s.m_x,m_c[0],m_c[1],m_h,m_s); 
    }
    inline CScalar operator-() const
    { 
      return CScalar(-m_x,m_c[0],m_c[1],m_h,m_s); 
    }

    inline CScalar& operator+=(const CScalar &s) 
    {
      m_x+=s.m_x;
      return *this;
    }
    inline CScalar& operator-=(const CScalar &s) 
    {
      m_x-=s.m_x;
      return *this;
    }
    inline CScalar& operator*=(const SComplex &c) 
    {
      m_x*=c;
      return *this;
    }
  
    inline CScalar Conj() const 
    {
      return CScalar(std::conj(m_x),m_c[0],m_c[1],m_h,m_s);
    }
    inline SComplex Abs2() const 
    {
      return m_x*m_x;
    }
    inline SComplex Abs() const 
    { 
      return sqrt(Abs2()); 
    }

    inline bool Nan() const { return ATOOLS::IsNan<SComplex>(m_x); }

    static void ResetAccu();

    inline static void   SetAccu(const Scalar &accu) { s_accu=accu;   }
    inline static Scalar Accu()                      { return s_accu; }

  };// end of class CScalar

  template <class Scalar> inline std::complex<Scalar>
  operator*(const Scalar &d,const CScalar<Scalar> &s)
  { return CScalar<Scalar>(s,d); }
  template <class Scalar> inline CScalar<Scalar> 
  operator*(const CScalar<Scalar> &s,const Scalar &d)
  { return CScalar<Scalar>(s,d); }
  template <class Scalar> inline CScalar<Scalar>
  operator*(const std::complex<Scalar> &c,const CScalar<Scalar> &s)
  { return CScalar<Scalar>(s,c); }
  template <class Scalar> inline CScalar<Scalar> 
  operator*(const CScalar<Scalar> &s,const std::complex<Scalar> &c)
  { return CScalar<Scalar>(s,c); }
  template <class Scalar> inline CScalar<Scalar>
  operator/(const CScalar<Scalar> &s,const std::complex<Scalar> &c)
  { return CScalar<Scalar>(s,1.0/c); }

  template <class Scalar> inline std::complex<Scalar> 
  operator*(const CScalar<Scalar> &s1,const CScalar<Scalar> &s2) 
  { return s1[0]*s2[0]; }

  template <class Scalar>
  std::ostream &operator<<(std::ostream &str,const CScalar<Scalar> &s);

}// end of namespace ATOOLS

#define DCScalar METOOLS::CScalar<double>
#define QCScalar METOOLS::CScalar<long double>

#endif
