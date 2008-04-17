#ifndef Spinor_H
#define Spinor_H

#include "MyComplex.H"
#include "Vector.H"

#include <vector>

namespace ATOOLS {

  class Spinor {
  protected:
    
    static double s_accu;
    static unsigned int s_r1, s_r2, s_r3;

    int m_r;

    Complex m_u1, m_u2;

    friend std::ostream &operator<<(std::ostream &ostr,const Spinor &s); 

  public:

    // constructor
    inline Spinor(const int &r=1): 
      m_r(r), m_u1(0.0), m_u2(0.0) {}
    inline Spinor(const int &r,const Complex &u1,const Complex &u2):
      m_r(r), m_u1(u1), m_u2(u2) {}
    inline Spinor(const int &r,const ATOOLS::Vec4D &p):
      m_r(r) { Construct(p); }
    
    // member functions
    void Construct(const ATOOLS::Vec4D &p);

    Complex operator*(const Spinor &s) const;

    Spinor operator*(const double &d) const;
    Spinor operator*(const Complex &c) const;
    Spinor operator/(const double &d) const;
    Spinor operator/(const Complex &c) const;

    Spinor operator*=(const double &d); 
    Spinor operator*=(const Complex &c); 
    Spinor operator/=(const double &d);
    Spinor operator/=(const Complex &c); 

    Spinor operator+(const Spinor &s) const;
    Spinor operator-(const Spinor &s) const;

    Spinor operator+=(const Spinor &s); 
    Spinor operator-=(const Spinor &s);

    bool operator==(const Spinor &s) const;

    static void SetGauge(const int gauge);

    // inline functions
    inline Complex &operator[](const size_t &i) 
    { return i==0?m_u1:m_u2; }
    inline Complex operator()(const size_t &i) const 
    { return i==0?m_u1:m_u2; }

    inline Complex U1() const { return m_u1; }
    inline Complex U2() const { return m_u2; }

    inline int R() const { return m_r; }

    inline double PPlus(const ATOOLS::Vec4D &p) const  
    { return p[0]+p[s_r3]; }
    inline double PMinus(const ATOOLS::Vec4D &p) const 
    { return p[0]-p[s_r3]; }

    inline Complex PT(const ATOOLS::Vec4D &p) const  
    { return Complex(p[s_r1],p[s_r2]); }
    inline Complex PTC(const ATOOLS::Vec4D &p) const  
    { return Complex(p[s_r1],-p[s_r2]); }

    inline Spinor operator-() const
    { return Spinor(m_r,-m_u1,-m_u2); }

    inline static void SetAccuracy(const double &accu) 
    { s_accu=accu; }
    inline static void ResetAccuracy() 
    { s_accu=1.0e-12; }

    inline static double Accuracy() { return s_accu; }

  };// end of class Spinor

  std::ostream &operator<<(std::ostream &ostr,const Spinor &s); 
  std::ostream &operator<<(std::ostream &ostr,
			   const std::vector<Spinor> &v); 

}// end of namespace ATOOLS

#endif
