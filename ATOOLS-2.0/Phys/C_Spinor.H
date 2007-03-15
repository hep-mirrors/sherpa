#ifndef CSpinor_H
#define CSpinor_H

#include "MyComplex.H"
#include "Vector.H"

#include <vector>

namespace ATOOLS {

  class CSpinor {
  protected:
    
    static double s_accu;

    int m_r, m_c;

    Complex m_u[4];

    friend std::ostream &operator<<(std::ostream &ostr,const CSpinor &s); 

  public:

    // constructor
    inline CSpinor(const int &r=1,const int c=0): 
      m_r(r), m_c(c) { m_u[0]=m_u[1]=m_u[2]=m_u[3]=0.0; }
    inline CSpinor(const int &r,const Complex &u1,const Complex &u2,
		   const Complex &u3,const Complex &u4,const int c=0):
      m_r(r), m_c(c) { m_u[0]=u1; m_u[1]=u2; m_u[2]=u3; m_u[3]=u4; }
    inline CSpinor(const int &r,const int &h,const Vec4D &p,const int c=0):
      m_r(r), m_c(c) { Construct(h,p); }
    
    // member functions
    void Construct(const int h,const Vec4D &p);

    Complex operator*(const CSpinor &s) const;

    CSpinor operator*(const double &d) const;
    CSpinor operator*(const Complex &c) const;
    CSpinor operator/(const double &d) const;
    CSpinor operator/(const Complex &c) const;

    CSpinor operator*=(const double &d); 
    CSpinor operator*=(const Complex &c); 
    CSpinor operator/=(const double &d);
    CSpinor operator/=(const Complex &c); 

    CSpinor operator+(const CSpinor &s) const;
    CSpinor operator-(const CSpinor &s) const;

    CSpinor operator+=(const CSpinor &s); 
    CSpinor operator-=(const CSpinor &s);

    bool operator==(const CSpinor &s) const;

    // inline functions
    inline Complex &operator[](const size_t &i) 
    { return m_u[i]; }
    inline const Complex &operator[](const size_t &i) const 
    { return m_u[i]; }
    inline int &operator()() 
    { return m_c; }
    inline const int &operator()() const 
    { return m_c; }

    inline int R() const { return m_r; }

    inline double PPlus(const Vec4D &p) const  
    { return p[0]+p[3]; }
    inline double PMinus(const Vec4D &p) const 
    { return p[0]-p[3]; }

    inline Complex PT(const Vec4D &p) const  
    { return Complex(p[1],p[2]);  }

    inline CSpinor operator-() const
    { return CSpinor(m_r,-m_u[0],-m_u[1],-m_u[2],-m_u[3],m_c); }

    inline CSpinor Bar() const
    { return CSpinor(-m_r,std::conj(m_u[0]),std::conj(m_u[1]),
		     -std::conj(m_u[2]),-std::conj(m_u[3]),m_c); }
    
    inline static void SetAccuracy(const double &accu) 
    { s_accu=accu; }
    inline static void ResetAccuracy() 
    { s_accu=1.0e-12; }

    inline static double Accuracy() { return s_accu; }

  };// end of class CSpinor

  std::ostream &operator<<(std::ostream &ostr,const CSpinor &s); 
  std::ostream &operator<<(std::ostream &ostr,
			   const std::vector<CSpinor> &v); 

}// end of namespace BCF

#endif
