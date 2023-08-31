#ifndef METOOLS_Explicit_C_RaritaSchwinger_H
#define METOOLS_Explicit_C_RaritaSchwinger_H

#include "METOOLS/Explicit/C_Object.H"
#include "ATOOLS/Math/Vec4.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Math/Vector.H"

#include <vector>
 //TODO: Zum schluss nochmal durchgehen durch Vector und Spinor und schauen, ob wir alle nötigen Funktionen haben...
 // TODO: Bar() und CConj() Funktionen implementieren analog zu Spin 1/2
 // TODO: Majorana-3/2 Teilchen?

namespace METOOLS {

  template <class Scalar>
  class CRaritaSchwinger: public CObject {
  public:

    typedef std::complex<Scalar> SComplex;

  private:

    SComplex m_x[16];

    // TODO: was sind m_b, m_on bei C_Spinor; m_s, m_c bei CObject? Ist m_h wirklich die Helicität, m_r, ob positive
    //  oder negative Energie
    // positive or negative energy vector-spinor
    // bar or non-bar vector-spinor
    int m_r, m_b;

    static double s_accu;

    static ATOOLS::AutoDelete_Vector<CRaritaSchwinger> s_objects;

    template <class _Scalar> friend std::ostream &
    operator<<(std::ostream &s,const CRaritaSchwinger<_Scalar> &vec);

  public:

    static CRaritaSchwinger *New();
    static CRaritaSchwinger *New(const CRaritaSchwinger &s);
    // TODO: notwendige Parameter RS anpassen!!!
    static CRaritaSchwinger *New(const Scalar &x0, const Scalar &x1,
		      const Scalar &x2, const Scalar &x3,
		      const int c1=0,const int c2=0,
		      const size_t &hm=0,
		      const size_t &hp=0);

    CObject* Copy() const;

    void Delete();

    bool IsZero() const;

    // TODO: Sinnvolle inlines ausdenken entsprechend der für RS nötigen Parameter!!!
    // TODO: Parameter anpassen, wenn Bedeutung geklärt!!! Bislang m_c, m_h, m_s, m_r, m_b da nötig (m_r, m_h, m_b) oder
    //       schon in CObject gesetzt
   inline CRaritaSchwinger(const int r=1, const int b=1, const int cr=0,const int ca=0,
                           const size_t &h=0,const size_t &s=0)
    { 
      for (auto entry : m_x){
        entry = 0.0;
      }
      m_c[0]=cr; m_c[1]=ca;
      m_h=h; m_s=s; m_r=r; m_b=b;
    }
    inline CRaritaSchwinger(const CRaritaSchwinger &v)
    {
      for (size_t i(0); i<16; ++i){
        m_x[i] = v[i];
      }
      m_c[0]=v(0); m_c[1]=v(1);
      m_h=v.m_h; m_s=v.m_s; m_r=v.m_r; m_b=v.m_b;
    }
    /*inline CRaritaSchwinger(const Scalar &x0, const Scalar &x1,
		 const Scalar &x2, const Scalar &x3,
		 const int c1=0,const int c2=0,
		 const size_t &h=0,const size_t &s=0)
    { 
      m_x[0]=x0; m_x[1]=x1; m_x[2]=x2; m_x[3]=x3; 
      m_c[0]=c1; m_c[1]=c2; m_h=h; m_s=s;
    }
    inline CRaritaSchwinger(const int c1,const int c2,
		 const size_t &h=0,const size_t &s=0)
    { 
      m_x[0]=m_x[1]=m_x[2]=m_x[3]=0.0; 
      m_c[0]=c1; m_c[1]=c2; m_h=h; m_s=s;
    }
    inline CRaritaSchwinger(const SComplex &x0, const SComplex &x1,
		 const SComplex &x2, const SComplex &x3,
		 const int c1=0,const int c2=0,
		 const size_t &h=0,const size_t &s=0)
    { 
      m_x[0]=x0; m_x[1]=x1; m_x[2]=x2; m_x[3]=x3; 
      m_c[0]=c1; m_c[1]=c2; m_h=h; m_s=s;
    }
    inline CRaritaSchwinger(const ATOOLS::Vec4<Scalar> &v,
		 const int c1=0,const int c2=0,
		 const size_t &h=0,const size_t &s=0)
    { 
      m_x[0]=v[0]; m_x[1]=v[1]; m_x[2]=v[2]; m_x[3]=v[3]; 
      m_c[0]=c1; m_c[1]=c2; m_h=h; m_s=s;
    }
    inline CRaritaSchwinger(const CRaritaSchwinger &v,const Scalar &c)
    {
      m_x[0]=v[0]*c; m_x[1]=v[1]*c; m_x[2]=v[2]*c; m_x[3]=v[3]*c;
      m_c[0]=v(0); m_c[1]=v(1); m_h=v.m_h; m_s=v.m_s;
    }
    inline CRaritaSchwinger(const CRaritaSchwinger &v,const SComplex &c)
    {
      m_x[0]=v[0]*c; m_x[1]=v[1]*c; m_x[2]=v[2]*c; m_x[3]=v[3]*c;
      m_c[0]=v(0); m_c[1]=v(1); m_h=v.m_h; m_s=v.m_s;
    }*/

// TODO: Funktionen füllen!!!
/*void Add(const CObject *c);
    void Divide(const double &d);
    void Multiply(const Complex &c);
    void Invert();*/

// TODO: Spinorwellenfunktionsklasse hat im Vergleich zur Vektorklasse noch ein Construct, SetOn, CConj -  brauchen wir
//       das hier???
// TODO: Welche Operationen brauchen wir für Rarita-Schwinger-Teilchen? Gern auch mit anderen Wellenfunktionen C_...
//       vergleichen!!!

   inline SComplex &operator[](const int i) { return m_x[i]; }
   inline const SComplex &operator[](const int i) const { return m_x[i]; }

/*     inline CVec4 operator+(const CVec4 &v) const
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
    }*/

    inline int R() const { return m_r; }
    inline int B() const { return m_b; }
    inline void SetR(const int &r) { m_r=r; }
    inline void SetB(const int &b) { m_b=b; }

    bool Nan() const;

    static void ResetAccu();

    inline static void   SetAccu(const double &accu) { s_accu=accu;   }
    inline static double Accu()                      { return s_accu; }

    bool Test_Properties(const ATOOLS::Vec4D &p);

  };
// end of class CVec4
// TODO: Brauchen wir Templates und wenn ja, welche?
  /*template <class Scalar> inline CVec4<Scalar>
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
  { return v1[0]*v2[0]-v1[1]*v2[1]-v1[2]*v2[2]-v1[3]*v2[3]; }*/

  template <class Scalar>
  std::ostream &operator<<(std::ostream &s,const CRaritaSchwinger<Scalar> &vec);
}// end of namespace ATOOLS

#define DCRaritaSchwinger METOOLS::CRaritaSchwinger<double>
#define QCRaritaSchwinger METOOLS::CRaritaSchwinger<long double>

#endif
