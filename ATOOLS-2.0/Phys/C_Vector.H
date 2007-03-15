#ifndef C_Vector_H
#define C_Vector_H

#include "Vector.H"
 
namespace ATOOLS {

  class CVec4D {
  private:

    Complex m_x[4];
    int     m_c[2];

    static double s_accu;

    friend std::ostream &operator<<(std::ostream &s,const CVec4D &vec);

  public:

    inline CVec4D()
    { 
      m_x[0]=m_x[1]=m_x[2]=m_x[3]=0.0; 
      m_c[0]=m_c[1]=-1; 
    }
    inline CVec4D(const CVec4D &v)
    { 
      m_x[0]=v[0]; m_x[1]=v[1]; m_x[2]=v[2]; m_x[3]=v[3];
      m_c[0]=v(0); m_c[1]=v(1);
    }
    inline CVec4D(const double &x0, const double &x1, 
		  const double &x2, const double &x3,
		  const int c1=-1,const int c2=-1)
    { 
      m_x[0]=x0; m_x[1]=x1; m_x[2]=x2; m_x[3]=x3; 
      m_c[0]=c1; m_c[1]=c2; 
    }
    inline CVec4D(const int c1,const int c2)
    { 
      m_x[0]=m_x[1]=m_x[2]=m_x[3]=0.0; 
      m_c[0]=c1; m_c[1]=c2; 
    }
    inline CVec4D(const Complex &x0, const Complex &x1, 
		  const Complex &x2, const Complex &x3,
		  const int c1=-1,const int c2=-1)
    { 
      m_x[0]=x0; m_x[1]=x1; m_x[2]=x2; m_x[3]=x3; 
      m_c[0]=c1; m_c[1]=c2; 
    }
    inline CVec4D(const Vec4D &v,const int c1=-1,const int c2=-1)
    { 
      m_x[0]=v[0]; m_x[1]=v[1]; m_x[2]=v[2]; m_x[3]=v[3]; 
      m_c[0]=c1; m_c[1]=c2; 
    }
    inline CVec4D(const CVec4D &v,const Complex &c)
    {
      m_x[0]=v[0]*c; m_x[1]=v[1]*c; m_x[2]=v[2]*c; m_x[3]=v[3]*c;
      m_c[0]=v(0); m_c[1]=v(1);
    }

    inline Complex &operator[](const int i) { return m_x[i]; }
    inline int     &operator()(const int i) { return m_c[i]; }

    inline const Complex &operator[](const int i) const { return m_x[i]; }
    inline const int     &operator()(const int i) const { return m_c[i]; }

    inline CVec4D operator+(const CVec4D &v) const  
    { 
      return CVec4D(m_x[0]+v[0],m_x[1]+v[1],
		    m_x[2]+v[2],m_x[3]+v[3],m_c[0],m_c[1]); 
    }
    inline CVec4D operator-(const CVec4D &v) const
    { 
      return CVec4D(m_x[0]-v[0],m_x[1]-v[1],
		    m_x[2]-v[2],m_x[3]-v[3],m_c[0],m_c[1]); 
    }
    inline CVec4D operator-() const
    { 
      return CVec4D(-m_x[0],-m_x[1],-m_x[2],-m_x[3],m_c[0],m_c[1]); 
    }

    inline CVec4D& operator+=(const CVec4D &v) 
    {
      m_x[0]+=v[0]; m_x[1]+=v[1]; m_x[2]+=v[2]; m_x[3]+=v[3];
      return *this;
    }
    inline CVec4D& operator-=(const CVec4D &v) 
    {
      m_x[0]-=v[0]; m_x[1]-=v[1]; m_x[2]-=v[2]; m_x[3]-=v[3];
      return *this;
    }
    inline CVec4D& operator*=(const Complex &c) 
    {
      m_x[0]*=c; m_x[1]*=c; m_x[2]*=c; m_x[3]*=c;
      return *this;
    }
  
    inline CVec4D Conj() const 
    {
      return CVec4D(std::conj(m_x[0]),std::conj(m_x[1]),
		    std::conj(m_x[2]),std::conj(m_x[3]),m_c[0],m_c[1]);
    }
    inline Complex Abs2() const 
    {
      return m_x[0]*m_x[0]-m_x[1]*m_x[1]-m_x[2]*m_x[2]-m_x[3]*m_x[3];
    }
    inline Complex Abs() const 
    { 
      return sqrt(Abs2()); 
    }

    bool Nan() const;
    bool IsZero() const;

    static void ResetAccu();

    inline static void   SetAccu(const double &accu) { s_accu=accu;   }
    inline static double Accu()                      { return s_accu; }

  };// end of class CVec4D

  inline CVec4D operator*(const Complex &c,const CVec4D &v)
  { return CVec4D(v,c); }
  inline CVec4D operator*(const CVec4D &v,const Complex &c)
  { return CVec4D(v,c); }
  inline CVec4D operator/(const CVec4D &v,const Complex &c)
  { return CVec4D(v,1.0/c); }

  inline Complex operator*(const CVec4D &v1,const CVec4D &v2) 
  { return v1[0]*v2[0]-v1[1]*v2[1]-v1[2]*v2[2]-v1[3]*v2[3]; }

}// end of namespace ATOOLS

#endif
