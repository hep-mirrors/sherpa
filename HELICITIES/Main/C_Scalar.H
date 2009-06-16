#ifndef HELICITIES_Main_C_Scalar_H
#define HELICITIES_Main_C_Scalar_H

#include "ATOOLS/Math/MathTools.H"

#include <iostream>
 
namespace ATOOLS {

  template <class Scalar>
  class CScalar {
  public:

    typedef std::complex<Scalar> SComplex;

  private:

    SComplex m_x;
    int      m_c;

    long unsigned int m_h[2];

    static double s_accu;

    template <class _Scalar> friend std::ostream &
    operator<<(std::ostream &str,const CScalar<_Scalar> &s);

  public:

    inline CScalar(): m_x(Scalar(0.0)), m_c(0)
    { 
      m_h[0]=m_h[1]=0; 
    }
    inline CScalar(const CScalar &s): m_x(s.m_x), m_c(0)
    { 
      m_h[0]=s.m_h[0]; m_h[1]=s.m_h[1];
    }
    inline CScalar(const Scalar &x,
		   const long unsigned int &hm=0,
		   const long unsigned int &hp=0): m_x(x), m_c(0)
    { 
      m_h[0]=hm; m_h[1]=hp;
    }
    inline CScalar(const SComplex &x,
		   const long unsigned int &hm=0,
		   const long unsigned int &hp=0): m_x(x), m_c(0)
    { 
      m_h[0]=hm; m_h[1]=hp;
    }
    inline CScalar(const CScalar &s,const SComplex &c): m_x(s.m_x*c), m_c(0)
    {
      m_h[0]=s.m_h[0]; m_h[1]=s.m_h[1];
    }

    inline SComplex &operator[](const int i) { return m_x; }
    inline int      &operator()(const int i) { return m_c; }

    inline const SComplex &operator[](const int i) const { return m_x; }
    inline const int      &operator()(const int i) const { return m_c; }

    inline void SetH(const long unsigned int &h,const int i) { m_h[i]=h; }

    inline long unsigned int H(const int i) const { return m_h[i]; }

    inline CScalar operator+(const CScalar &s) const  
    { 
      return CScalar(m_x+s.m_x,m_h[0],m_h[1]); 
    }
    inline CScalar operator-(const CScalar &s) const
    { 
      return CScalar(m_x-s.m_x,m_h[0],m_h[1]); 
    }
    inline CScalar operator-() const
    { 
      return CScalar(-m_x,m_h[0],m_h[1]); 
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
      return CScalar(std::conj(m_x),m_h[0],m_h[1]);
    }
    inline SComplex Abs2() const 
    {
      return m_x*m_x;
    }
    inline SComplex Abs() const 
    { 
      return sqrt(Abs2()); 
    }

    inline bool Nan() const    { return ATOOLS::IsNan<SComplex>(m_x);  }
    inline bool IsZero() const { return ATOOLS::IsZero<SComplex>(m_x); }

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

#define DCScalar ATOOLS::CScalar<double>
#define QCScalar ATOOLS::CScalar<long double>

#endif
