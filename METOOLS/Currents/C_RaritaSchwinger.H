#ifndef METOOLS_Explicit_C_RaritaSchwinger_H
#define METOOLS_Explicit_C_RaritaSchwinger_H

#include "METOOLS/Explicit/C_Object.H"
#include "ATOOLS/Math/Vec4.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Matrix.H"
#include "ATOOLS/Org/Exception.H"
#include "METOOLS/Currents/C_Spinor.H"

#define ZERO std::complex<Scalar>(0.0,0.0)

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
    // positive or negative energy vector-spinor (particle or anti-particle)
    // bar or non-bar vector-spinor
    // massless (m_on=1,2) or massive spin-3/2 particle (m_on=3): helper variable to left physically ZERO components
    // also ZERO numerically
    int m_r, m_b, m_on;

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
                           const size_t &h=0,const size_t &s=0, const int &on=3)
    { 
      for (auto entry : m_x){
        entry = 0.0;
      }
      m_c[0]=cr; m_c[1]=ca;
      m_h=h; m_s=s; m_r=r; m_b=b; m_on=on;
    }
    inline CRaritaSchwinger(const CRaritaSchwinger &v)
    {
      for (size_t i(0); i<16; ++i){
        m_x[i] = v[i];
      }
      m_c[0]=v(0); m_c[1]=v(1);
      m_h=v.m_h; m_s=v.m_s; m_r=v.m_r; m_b=v.m_b; m_on=v.m_on;
    }
    inline CRaritaSchwinger(const int &r, const int &b, const SComplex &x0, const SComplex &x1, const SComplex &x2,
                            const SComplex &x3, const SComplex &x4, const SComplex &x5, const SComplex &x6,
                            const SComplex &x7, const SComplex &x8, const SComplex &x9, const SComplex &x10,
                            const SComplex &x11, const SComplex &x12, const SComplex &x13, const SComplex &x14,
                            const SComplex &x15, const int c1=0,const int c2=0, const size_t &h=0,const size_t &s=0,
                            const int &on=3)
    {
      m_x[0]=x0; m_x[1]=x1; m_x[2]=x2; m_x[3]=x3; m_x[4]=x4; m_x[5]=x5; m_x[6]=x6; m_x[7]=x7; m_x[8]=x8; m_x[9]=x9;
      m_x[10]=x10; m_x[11]=x11; m_x[12]=x12; m_x[13]=x13; m_x[14]=x14; m_x[15]=x15;
      m_r=r, m_b=b, m_on =on;
      m_c[0]=c1; m_c[1]=c2; m_h=h; m_s=s;
    }

    /*inline CRaritaSchwinger(const Scalar &x0, const Scalar &x1,
		 const Scalar &x2, const Scalar &x3,
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
    void Add(const CObject *c);
    void Divide(const double &d);
    void Multiply(const Complex &c);
    void Invert();

// TODO: Spinorwellenfunktionsklasse hat im Vergleich zur Vektorklasse noch ein Construct, SetOn, CConj -  brauchen wir
//       das hier???
// TODO: Welche Operationen brauchen wir für Rarita-Schwinger-Teilchen? Gern auch mit anderen Wellenfunktionen C_...
//       vergleichen!!!

   inline SComplex &operator[](const int i) { return m_x[i]; }
   inline const SComplex &operator[](const int i) const { return m_x[i]; }
   ATOOLS::TCMatrix<Scalar> operator*(const CRaritaSchwinger<Scalar> &rs) const;

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
    inline int On() const { return m_on; }
    inline void SetR(const int &r) { m_r=r; }
    inline void SetB(const int &b) { m_b=b; }
    bool SetOn();

    bool Nan() const;

    //TODO: Das ist jetzt Spinbasisabhängig, wir könnten das Ganze alternativ mithilfe der SpinFuncs implementieren,
    // dann könnte man einfacher die Eichung wechseln, NACHTEILE: recht viele (triviale) Multiplikationen nötig + evtl.
    // sind dann die Komponenten nicht mehr exakt die gleichen wegen einer zusätzlichen numerischen Rechnung?
    inline CRaritaSchwinger<Scalar> Bar(){
      return CRaritaSchwinger<Scalar>(1, -1, std::conj(m_x[2]), std::conj(m_x[3]),
                                      std::conj(m_x[0]), std::conj(m_x[1]), std::conj(m_x[6]), std::conj(m_x[7]),
                                      std::conj(m_x[4]), std::conj(m_x[5]), std::conj(m_x[10]), std::conj(m_x[11]),
                                      std::conj(m_x[8]), std::conj(m_x[9]), std::conj(m_x[14]), std::conj(m_x[15]),
                                      std::conj(m_x[12]), std::conj(m_x[13]), m_c[0], m_c[1], m_h, m_s, m_on); }

    static void ResetAccu();

    inline static void   SetAccu(const double &accu) { s_accu=accu;   }
    inline static double Accu()                      { return s_accu; }

    bool Test_Properties(const ATOOLS::Vec4D &p, int r);

  };


  // end of class CRaritaSchwinger

  // TODO: Sollten wir noch irgendwelche Templates/Operatoren o.ä. implementieren, die wir bislang noch nicht gebraucht
  //       haben?
  template <class Scalar> inline CRaritaSchwinger<Scalar>
  operator*(const Scalar &c,const CRaritaSchwinger<Scalar> &v)
  { return CRaritaSchwinger<Scalar>(v,c); }
  template <class Scalar> inline CRaritaSchwinger<Scalar>
  operator*(const CRaritaSchwinger<Scalar> &v,const Scalar &c)
  { return CRaritaSchwinger<Scalar>(v,c); }
  template <class Scalar> inline CRaritaSchwinger<Scalar>
  operator*(const std::complex<Scalar> &c,const CRaritaSchwinger<Scalar> &v)
  { return CRaritaSchwinger<Scalar>(v,c); }
  template <class Scalar> inline CRaritaSchwinger<Scalar>
  operator*(const CRaritaSchwinger<Scalar> &v,const std::complex<Scalar> &c)
  { return CRaritaSchwinger<Scalar>(v,c); }
  template <class Scalar> inline CRaritaSchwinger<Scalar>
  operator/(const CRaritaSchwinger<Scalar> &v,const std::complex<Scalar> &c)
  { return CRaritaSchwinger<Scalar>(v,Scalar(1.0)/c); }
  /*
  template <class Scalar> inline std::complex<Scalar> 
  operator*(const CVec4<Scalar> &v1,const ATOOLS::Vec4<Scalar> &v2)
  { return v1[0]*v2[0]-v1[1]*v2[1]-v1[2]*v2[2]-v1[3]*v2[3]; }
  template <class Scalar> inline std::complex<Scalar>
  operator*(const ATOOLS::Vec4<Scalar> &v1,const CVec4<Scalar> &v2)
  { return v1[0]*v2[0]-v1[1]*v2[1]-v1[2]*v2[2]-v1[3]*v2[3]; }
  template <class Scalar> inline std::complex<Scalar>
  operator*(const CVec4<Scalar> &v1,const CVec4<Scalar> &v2)
  { return v1[0]*v2[0]-v1[1]*v2[1]-v1[2]*v2[2]-v1[3]*v2[3]; } */

  template <class Scalar>
  std::ostream &operator<<(std::ostream &s,const CRaritaSchwinger<Scalar> &vec);
}// end of namespace ATOOLS

#define DCRaritaSchwinger METOOLS::CRaritaSchwinger<double>
#define QCRaritaSchwinger METOOLS::CRaritaSchwinger<long double>

#endif
