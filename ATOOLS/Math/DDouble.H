#ifndef ATOOLS_Math_DDouble_H
#define ATOOLS_Math_DDouble_H

#include <cmath>
#include <complex>
#include <iostream>
#include <limits>

namespace ATOOLS {

  /*!
    \class DDouble

    A value is carried as an exact sum m_h + m_l of two doubles with
    |m_l| <= ulp(m_h)/2, so the pair represents roughly twice the significand
    of a double at about 1e-32 relative accuracy. The arithmetic is the
    classical Dekker/Knuth error-free-transformation set: every operation
    computes the double-precision result together with its exact rounding
    error and folds the error back in.

    This exists because Comix's accuracy on the YFS real emission is limited by
    cancellation inside the current recursion, and `long double` is not a
    portable answer to that: it is 80-bit on x86 but a plain double on arm64,
    so a PMODE built on it gives different physics on different machines. This
    type is identical everywhere and needs no compiler or library support
  */
  class DDouble {
  private:

    double m_h, m_l;

    // Error-free transformations. two_sum is Knuth's; quick_two_sum needs
    // |a| >= |b|
    static inline void QuickTwoSum(double a,double b,double &s,double &e)
    { s=a+b; e=b-(s-a); }
    static inline void TwoSum(double a,double b,double &s,double &e)
    { s=a+b; const double bb(s-a); e=(a-(s-bb))+(b-bb); }
    static inline void TwoProd(double a,double b,double &p,double &e)
    { p=a*b; e=std::fma(a,b,-p); }

  public:

    inline DDouble(): m_h(0.0), m_l(0.0) {}
    inline DDouble(const double &h): m_h(h), m_l(0.0) {}
    inline DDouble(const int &i): m_h(i), m_l(0.0) {}
    inline DDouble(const double &h,const double &l): m_h(h), m_l(l) {}

    inline double High() const { return m_h; }
    inline double Low() const  { return m_l; }

    inline explicit operator double() const { return m_h+m_l; }
    inline explicit operator long double() const
    { return (long double)m_h+(long double)m_l; }

    inline DDouble operator-() const { return DDouble(-m_h,-m_l); }

    friend inline DDouble operator+(const DDouble &a,const DDouble &b)
    {
      double s,e,s2,e2;
      TwoSum(a.m_h,b.m_h,s,e);
      TwoSum(a.m_l,b.m_l,s2,e2);
      e+=s2;
      QuickTwoSum(s,e,s,e);
      e+=e2;
      QuickTwoSum(s,e,s,e);
      return DDouble(s,e);
    }

    friend inline DDouble operator-(const DDouble &a,const DDouble &b)
    { return a+(-b); }

    friend inline DDouble operator*(const DDouble &a,const DDouble &b)
    {
      double p,e;
      TwoProd(a.m_h,b.m_h,p,e);
      e+=a.m_h*b.m_l+a.m_l*b.m_h;
      QuickTwoSum(p,e,p,e);
      return DDouble(p,e);
    }

    friend inline DDouble operator/(const DDouble &a,const DDouble &b)
    {
      const double q1(a.m_h/b.m_h);
      DDouble r(a-b*DDouble(q1));
      const double q2(r.m_h/b.m_h);
      r=r-b*DDouble(q2);
      const double q3(r.m_h/b.m_h);
      double s,e;
      QuickTwoSum(q1,q2,s,e);
      e+=q3;
      QuickTwoSum(s,e,s,e);
      return DDouble(s,e);
    }

    inline DDouble &operator+=(const DDouble &b) { return *this=*this+b; }
    inline DDouble &operator-=(const DDouble &b) { return *this=*this-b; }
    inline DDouble &operator*=(const DDouble &b) { return *this=*this*b; }
    inline DDouble &operator/=(const DDouble &b) { return *this=*this/b; }

    friend inline bool operator==(const DDouble &a,const DDouble &b)
    { return a.m_h==b.m_h && a.m_l==b.m_l; }
    friend inline bool operator!=(const DDouble &a,const DDouble &b)
    { return !(a==b); }
    friend inline bool operator<(const DDouble &a,const DDouble &b)
    { return a.m_h<b.m_h || (a.m_h==b.m_h && a.m_l<b.m_l); }
    friend inline bool operator>(const DDouble &a,const DDouble &b)
    { return b<a; }
    friend inline bool operator<=(const DDouble &a,const DDouble &b)
    { return !(b<a); }
    friend inline bool operator>=(const DDouble &a,const DDouble &b)
    { return !(a<b); }

  };// end of class DDouble

  // Declaring abs/sqrt/fabs in ATOOLS would otherwise hide the global ones
  // from every unqualified call inside the namespace: Vec3<double>::Abs()
  // calling sqrt(double) would convert its argument to DDouble and hand back
  // a DDouble. Pulling the standard overloads in keeps the set complete.
  using std::abs;
  using std::fabs;
  using std::sqrt;
  using std::isnan;
  using std::isinf;
  using std::isfinite;

  inline DDouble abs(const DDouble &a)
  { return a<DDouble(0.0)?-a:a; }
  inline DDouble fabs(const DDouble &a) { return abs(a); }
  inline DDouble dabs(const DDouble &a) { return abs(a); }

  inline bool isnan(const DDouble &a)
  { return std::isnan(a.High())||std::isnan(a.Low()); }
  inline bool isinf(const DDouble &a)
  { return std::isinf(a.High())||std::isinf(a.Low()); }
  inline bool isfinite(const DDouble &a)
  { return std::isfinite(a.High())&&std::isfinite(a.Low()); }

  inline DDouble sqrt(const DDouble &a)
  {
    if (a.High()==0.0 && a.Low()==0.0) return DDouble(0.0);
    if (a.High()<0.0)
      return DDouble(std::numeric_limits<double>::quiet_NaN());
    // one Newton step on 1/sqrt(a) in double, corrected once in DDouble:
    // x is accurate to ~1e-16, the correction lifts it to ~1e-32.
    const double x(1.0/std::sqrt(a.High())), ax(a.High()*x);
    const DDouble d(a-DDouble(ax)*DDouble(ax));
    return DDouble(ax)+DDouble(d.High()*x*0.5);
  }

  inline std::ostream &operator<<(std::ostream &s,const DDouble &a)
  { return s<<(double)a; }



  inline DDouble abs(const std::complex<DDouble> &z)
  {
    const DDouble ax(abs(z.real())), ay(abs(z.imag()));
    const DDouble m(ax>ay?ax:ay);
    if (m==DDouble(0.0)) return m;
    const DDouble rx(ax/m), ry(ay/m);
    return m*sqrt(rx*rx+ry*ry);
  }

  inline std::complex<DDouble> sqrt(const std::complex<DDouble> &z)
  {
    const DDouble zero(0.0), half(0.5);
    if (z.real()==zero && z.imag()==zero) return std::complex<DDouble>(zero,zero);
    // principal branch via the stable half-angle form
    const DDouble r(abs(z));
    DDouble re(sqrt((r+z.real())*half));
    DDouble im(sqrt((r-z.real())*half));
    if (z.imag()<zero) im=-im;
    return std::complex<DDouble>(re,im);
  }

  typedef ATOOLS::DDouble DD;

}// end of namespace ATOOLS

#endif
