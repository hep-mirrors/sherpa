#ifndef HELICITIES_Main_C_Tensor_H
#define HELICITIES_Main_C_Tensor_H

#include "HELICITIES/Main/C_Vector.H"
 
namespace ATOOLS {

  template <class Scalar>
  class CAsT4 {
  public:

    typedef std::complex<Scalar> SComplex;

  private:

    SComplex m_x[6];
    int      m_c[2];

    long unsigned int m_h[2];

    static double s_accu;

    template <class _Scalar> friend std::ostream &
    operator<<(std::ostream &s,const CAsT4<_Scalar> &ten);

  public:

    inline CAsT4()
    { 
      m_x[0]=m_x[1]=m_x[2]=m_x[3]=m_x[4]=m_x[5]=0.0; 
      m_c[0]=m_c[1]=-1;
      m_h[0]=m_h[1]=0;
    }
    inline CAsT4(const CAsT4 &t)
    { 
      m_x[0]=t[0]; m_x[1]=t[1]; m_x[2]=t[2]; 
      m_x[3]=t[3]; m_x[4]=t[4]; m_x[5]=t[5];
      m_c[0]=t(0); m_c[1]=t(1);
      m_h[0]=t.m_h[0]; m_h[1]=t.m_h[1];
    }
    inline CAsT4(const Scalar &x0, const Scalar &x1, 
		  const Scalar &x2, const Scalar &x3,
		  const Scalar &x4, const Scalar &x5,
		  const int c1=-1,const int c2=-1,
		  const long unsigned int &hm=0,
		  const long unsigned int &hp=0)
    { 
      m_x[0]=x0; m_x[1]=x1; m_x[2]=x2;
      m_x[3]=x3; m_x[4]=x4; m_x[5]=x5; 
      m_c[0]=c1; m_c[1]=c2;
      m_h[0]=hm; m_h[1]=hp;
    }
    inline CAsT4(const int c1,const int c2,
		  const long unsigned int &hm=0,
		  const long unsigned int &hp=0)
    { 
      m_x[0]=m_x[1]=m_x[2]=m_x[3]=m_x[4]=m_x[5]=0.0; 
      m_c[0]=c1; m_c[1]=c2; 
      m_h[0]=hm; m_h[1]=hp;
    }
    inline CAsT4(const SComplex &x0, const SComplex &x1, 
		  const SComplex &x2, const SComplex &x3,
		  const SComplex &x4, const SComplex &x5,
		  const int c1=-1,const int c2=-1,
		  const long unsigned int &hm=0,
		  const long unsigned int &hp=0)
    { 
      m_x[0]=x0; m_x[1]=x1; m_x[2]=x2;
      m_x[3]=x3; m_x[4]=x4; m_x[5]=x5; 
      m_c[0]=c1; m_c[1]=c2; 
      m_h[0]=hm; m_h[1]=hp;
    }
    inline CAsT4(const CVec4<Scalar> &v1,const CVec4<Scalar> &v2)
    { 
      m_x[0]=v1[0]*v2[1]-v1[1]*v2[0];// 0,1 
      m_x[1]=v1[1]*v2[2]-v1[2]*v2[1];// 1,2 
      m_x[2]=v1[2]*v2[3]-v1[3]*v2[2];// 2,3 
      m_x[3]=v1[0]*v2[2]-v1[2]*v2[0];// 0,2 
      m_x[4]=v1[1]*v2[3]-v1[3]*v2[1];// 1,3 
      m_x[5]=v1[0]*v2[3]-v1[3]*v2[0];// 0,3 
      m_c[0]=v2(0); m_c[1]=v1(1); 
      m_h[0]=v1.H(0)+v2.H(0); m_h[1]=v1.H(1)+v2.H(1);
    }
    inline CAsT4(const CAsT4 &v,const SComplex &c)
    {
      m_x[0]=v[0]*c; m_x[1]=v[1]*c; m_x[2]=v[2]*c;
      m_x[3]=v[3]*c; m_x[4]=v[4]*c; m_x[5]=v[5]*c;
      m_c[0]=v(0); m_c[1]=v(1); m_h[0]=v.m_h[0]; m_h[1]=v.m_h[1];
    }
    inline CAsT4(const CAsT4 &v,const Scalar &c)
    {
      m_x[0]=v[0]*c; m_x[1]=v[1]*c; m_x[2]=v[2]*c;
      m_x[3]=v[3]*c; m_x[4]=v[4]*c; m_x[5]=v[5]*c;
      m_c[0]=v(0); m_c[1]=v(1); m_h[0]=v.m_h[0]; m_h[1]=v.m_h[1];
    }

    inline SComplex &operator[](const int i) { return m_x[i]; }
    inline int     &operator()(const int i) { return m_c[i]; }

    inline const SComplex &operator[](const int i) const { return m_x[i]; }
    inline const int     &operator()(const int i) const { return m_c[i]; }

    inline void SetH(const long unsigned int &h,const int i) { m_h[i]=h; }

    inline long unsigned int H(const int i) const { return m_h[i]; }

    inline CAsT4 operator+(const CAsT4 &v) const  
    { 
      return CAsT4(m_x[0]+v[0],m_x[1]+v[1],m_x[2]+v[2],
		    m_x[3]+v[3],m_x[4]+v[4],m_x[5]+v[5],
		    m_c[0],m_c[1],m_h[0],m_h[1]); 
    }
    inline CAsT4 operator-(const CAsT4 &v) const
    { 
      return CAsT4(m_x[0]-v[0],m_x[1]-v[1],m_x[2]-v[2],
		    m_x[3]-v[3],m_x[4]-v[4],m_x[5]-v[5],
		    m_c[0],m_c[1],m_h[0],m_h[1]); 
    }
    inline CAsT4 operator-() const
    { 
      return CAsT4(-m_x[0],-m_x[1],-m_x[2],-m_x[3],-m_x[4],-m_x[5],
		    m_c[0],m_c[1],m_h[0],m_h[1]); 
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
    inline CAsT4& operator*=(const SComplex &c) 
    {
      m_x[0]*=c; m_x[1]*=c; m_x[2]*=c; m_x[3]*=c; m_x[4]*=c; m_x[5]*=c;
      return *this;
    }
  
    inline CAsT4 Conj() const 
    {
      return CAsT4(std::conj(m_x[0]),std::conj(m_x[1]),std::conj(m_x[2]),
		    std::conj(m_x[3]),std::conj(m_x[4]),std::conj(m_x[5]),
		    m_c[0],m_c[1],m_h[0],m_h[1]);
    }

    bool Nan() const;
    bool IsZero() const;

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
  template <class Scalar> inline CVec4<Scalar> 
  operator*(const CVec4<Scalar> &v,const CAsT4<Scalar> &t)
  { 
    CVec4<Scalar> j;
    j[0]=-t[0]/*-0,1*/*v[1]-t[3]/*-0,2*/*v[2]-t[5]/*-0,3*/*v[3];
    j[1]=-t[0]/*-0,1*/*v[0]-t[1]/*-1,2*/*v[2]-t[4]/*-1,3*/*v[3];
    j[2]=-t[3]/*-0,2*/*v[0]+t[1]/* 1,2*/*v[1]-t[2]/*-2,3*/*v[3];
    j[3]=-t[5]/*-0,3*/*v[0]+t[4]/* 1,3*/*v[1]+t[2]/* 2,3*/*v[2];
    return j;
  }
  template <class Scalar> inline CVec4<Scalar>
  operator*(const CAsT4<Scalar> &t,const CVec4<Scalar> &v)
  { return -(v*t); }

  template <class Scalar> 
  std::ostream &operator<<(std::ostream &s,const CAsT4<Scalar> &ten);

}// end of namespace ATOOLS

#define DCAsT4D ATOOLS::CAsT4<double>
#define QCAsT4D ATOOLS::CAsT4<long double>

#endif
