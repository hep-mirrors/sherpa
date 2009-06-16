#ifndef ATOOLS_Phys_Spinor_H
#define ATOOLS_Phys_Spinor_H

#include "ATOOLS/Math/Vec4.H"

#include <complex>
#include <vector>

namespace ATOOLS {

  template <class Scalar>
  class Spinor {
  public:

    typedef std::complex<Scalar> SComplex;

  protected:
    
    static double s_accu;
    static unsigned int s_r1, s_r2, s_r3;

    int m_r;

    SComplex m_u1, m_u2;

    template <class _Scalar> friend std::ostream &operator<<
      (std::ostream &ostr,const Spinor<_Scalar> &s); 

  public:

    // constructor
    inline Spinor(const int &r=1): 
      m_r(r), m_u1(0.0), m_u2(0.0) {}
    inline Spinor(const int &r,const SComplex &u1,const SComplex &u2):
      m_r(r), m_u1(u1), m_u2(u2) {}
    inline Spinor(const int &r,const Vec4<Scalar> &p):
      m_r(r) { Construct(p); }
    
    // member functions
    void Construct(const Vec4<Scalar> &p);

    SComplex operator*(const Spinor &s) const;

    Spinor operator*(const Scalar &d) const;
    Spinor operator*(const SComplex &c) const;
    Spinor operator/(const Scalar &d) const;
    Spinor operator/(const SComplex &c) const;

    Spinor operator*=(const Scalar &d); 
    Spinor operator*=(const SComplex &c); 
    Spinor operator/=(const Scalar &d);
    Spinor operator/=(const SComplex &c); 

    Spinor operator+(const Spinor &s) const;
    Spinor operator-(const Spinor &s) const;

    Spinor operator+=(const Spinor &s); 
    Spinor operator-=(const Spinor &s);

    bool operator==(const Spinor &s) const;

    static void SetGauge(const int gauge);

    // inline functions
    inline SComplex &operator[](const size_t &i) 
    { return i==0?m_u1:m_u2; }
    inline SComplex operator()(const size_t &i) const 
    { return i==0?m_u1:m_u2; }

    inline SComplex U1() const { return m_u1; }
    inline SComplex U2() const { return m_u2; }

    inline int R() const { return m_r; }

    inline Scalar PPlus(const Vec4<Scalar> &p) const  
    { return p[0]+p[s_r3]; }
    inline Scalar PMinus(const Vec4<Scalar> &p) const 
    { return p[0]-p[s_r3]; }

    inline SComplex PT(const Vec4<Scalar> &p) const  
    { return SComplex(p[s_r1],p[s_r2]); }
    inline SComplex PTC(const Vec4<Scalar> &p) const  
    { return SComplex(p[s_r1],-p[s_r2]); }

    inline Spinor operator-() const
    { return Spinor(m_r,-m_u1,-m_u2); }

    inline static void SetAccuracy(const double &accu) 
    { s_accu=accu; }
    inline static void ResetAccuracy() 
    { s_accu=1.0e-12; }

    inline static double Accuracy() { return s_accu; }

  };// end of class Spinor

  template <class Scalar> std::ostream &operator<<
    (std::ostream &ostr,const Spinor<Scalar> &s); 

}// end of namespace ATOOLS

#define DWSpinor ATOOLS::Spinor<double>
#define QWSpinor ATOOLS::Spinor<long double>

#endif
