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
 // TODO: Majorana-3/2 Teilchen?

namespace METOOLS {

  template <class Scalar>
  class CRaritaSchwinger: public CObject {
  public:

    typedef std::complex<Scalar> SComplex;

  private:

    SComplex m_x[16];

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
    static CRaritaSchwinger *New(const int r, const int b, const Scalar &x0, const Scalar &x1, const Scalar &x2,
                                 const Scalar &x3, const Scalar &x4, const Scalar &x5, const Scalar &x6,
                                 const Scalar &x7, const Scalar &x8, const Scalar &x9, const Scalar &x10,
                                 const Scalar &x11, const Scalar &x12, const Scalar &x13, const Scalar &x14,
                                 const Scalar &x15, const int c1=0,const int c2=0, const size_t &h=0,
                                 const size_t &s=0, const int on=3);

    CObject* Copy() const;

    void Delete();

    bool IsZero() const;

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

    inline CRaritaSchwinger(const CRaritaSchwinger &rs, const Scalar &c)
    {
      if (rs.m_on==1){
        m_x[0]= rs[0] * c; m_x[1]= rs[1] * c; m_x[2]= ZERO; m_x[3]= ZERO; m_x[4]= rs[4] * c; m_x[5]= rs[5] * c;
        m_x[6]= ZERO; m_x[7]= ZERO; m_x[8]= rs[8] * c; m_x[9]= rs[9] * c; m_x[10]= ZERO; m_x[11]= ZERO;
        m_x[12]= rs[12] * c; m_x[13]= rs[13] * c; m_x[14]= ZERO; m_x[15]= ZERO;
      }
      else if (rs.m_on==2){
        m_x[0]= ZERO; m_x[1]= ZERO; m_x[2]= rs[2] * c; m_x[3]= rs[3] * c; m_x[4]= ZERO; m_x[5]= ZERO; m_x[6]= rs[6] * c;
        m_x[7]= rs[7] * c; m_x[8]= ZERO; m_x[9]= ZERO; m_x[10]= rs[10] * c; m_x[11]= rs[11] * c; m_x[12]= ZERO;
        m_x[13]= ZERO; m_x[14]= rs[14] * c; m_x[15]= rs[15] * c;
      }
      else{
        m_x[0]= rs[0] * c; m_x[1]= rs[1] * c; m_x[2]= rs[2] * c; m_x[3]= rs[3] * c; m_x[4]= rs[4] * c;
        m_x[5]= rs[5] * c; m_x[6]= rs[6] * c; m_x[7]= rs[7] * c; m_x[8]= rs[8] * c; m_x[9]= rs[9] * c;
        m_x[10]= rs[10] * c; m_x[11]= rs[11] * c; m_x[12]= rs[12] * c; m_x[13]= rs[13] * c; m_x[14]= rs[14] * c;
        m_x[15]= rs[15] * c;
      }
      m_c[0]=rs(0); m_c[1]=rs(1); m_h=rs.m_h; m_s=rs.m_s; m_on=rs.m_on; m_r=rs.m_r; m_b=rs.m_b;
    }

    inline CRaritaSchwinger(const CRaritaSchwinger &rs, const SComplex &c)
    {
      if (m_on==1){
        m_x[0]= rs[0] * c; m_x[1]= rs[1] * c; m_x[2]= ZERO; m_x[3]= ZERO; m_x[4]= rs[4] * c; m_x[5]= rs[5] * c;
        m_x[6]= ZERO; m_x[7]= ZERO; m_x[8]= rs[8] * c; m_x[9]= rs[9] * c; m_x[10]= ZERO; m_x[11]= ZERO;
        m_x[12]= rs[12] * c; m_x[13]= rs[13] * c; m_x[14]= ZERO; m_x[15]= ZERO;
      }
      else if (m_on==2){
        m_x[0]= ZERO; m_x[1]= ZERO; m_x[2]= rs[2] * c; m_x[3]= rs[3] * c; m_x[4]= ZERO; m_x[5]= ZERO; m_x[6]= rs[6] * c;
        m_x[7]= rs[7] * c; m_x[8]= ZERO; m_x[9]= ZERO; m_x[10]= rs[10] * c; m_x[11]= rs[11] * c; m_x[12]= ZERO;
        m_x[13]= ZERO; m_x[14]= rs[14] * c; m_x[15]= rs[15] * c;
      }
      else{
        m_x[0]= rs[0] * c; m_x[1]= rs[1] * c; m_x[2]= rs[2] * c; m_x[3]= rs[3] * c; m_x[4]= rs[4] * c;
        m_x[5]= rs[5] * c; m_x[6]= rs[6] * c; m_x[7]= rs[7] * c; m_x[8]= rs[8] * c; m_x[9]= rs[9] * c;
        m_x[10]= rs[10] * c; m_x[11]= rs[11] * c; m_x[12]= rs[12] * c; m_x[13]= rs[13] * c; m_x[14]= rs[14] * c;
        m_x[15]= rs[15] * c;
      }
      m_c[0]=rs(0); m_c[1]=rs(1); m_h=rs.m_h; m_s=rs.m_s; m_on=rs.m_on; m_r=rs.m_r; m_b=rs.m_b;
    }

    void Add(const CObject *c);
    void Divide(const double &d);
    void Multiply(const Complex &c);
    void Invert();

   inline SComplex &operator[](const int i) { return m_x[i]; }
   inline const SComplex &operator[](const int i) const { return m_x[i]; }
   ATOOLS::TCMatrix<Scalar> Contract4Index(const CRaritaSchwinger<Scalar> &rs) const;
   ATOOLS::TCMatrix<Scalar> ContractSpinorIndex(const CRaritaSchwinger<Scalar> &rs) const;

   SComplex operator*(const CRaritaSchwinger<Scalar> &rs) const;
   inline CRaritaSchwinger operator+(const CRaritaSchwinger &rs) const
    {
      if (rs.m_b!=m_b || rs.m_r!=m_r) THROW(fatal_error, "Summing particle and anti-particle wave functions / bar "
                                                           "and non-bar wave functions is not allowed!")
      return CRaritaSchwinger(m_r, m_b, m_x[0]+rs[0], m_x[1]+rs[1], m_x[2]+rs[2], m_x[3]+rs[3], m_x[4]+rs[4],
                              m_x[5]+rs[5], m_x[6]+rs[6], m_x[7]+rs[7], m_x[8]+rs[8], m_x[9]+rs[9], m_x[10]+rs[10],
                              m_x[11]+rs[11], m_x[12]+rs[12], m_x[13]+rs[13], m_x[14]+rs[14], m_x[15]+rs[15],
                              m_c[0], m_c[1], m_h, m_s, m_on|rs.m_on);
    }
    inline CRaritaSchwinger operator-(const CRaritaSchwinger &rs) const
    {
      if (rs.m_b!=m_b || rs.m_r!=m_r) THROW(fatal_error, "Summing particle and anti-particle wave functions / bar "
                                                           "and non-bar wave functions is not allowed!")
      return CRaritaSchwinger(m_r, m_b, m_x[0]-rs[0], m_x[1]-rs[1], m_x[2]-rs[2], m_x[3]-rs[3], m_x[4]-rs[4],
                              m_x[5]-rs[5], m_x[6]-rs[6], m_x[7]-rs[7], m_x[8]-rs[8], m_x[9]-rs[9], m_x[10]-rs[10],
                              m_x[11]-rs[11], m_x[12]-rs[12], m_x[13]-rs[13], m_x[14]-rs[14], m_x[15]-rs[15],
                              m_c[0], m_c[1], m_h, m_s, m_on|rs.m_on);
    }
    inline CRaritaSchwinger operator-() const
    {
      return CRaritaSchwinger(m_r, m_b, -m_x[0], -m_x[1], -m_x[2], -m_x[3], -m_x[4], -m_x[5], -m_x[6], -m_x[7],
                              -m_x[8], -m_x[9], -m_x[10],-m_x[11], -m_x[12], -m_x[13], -m_x[14], -m_x[15],
                              m_c[0], m_c[1], m_h, m_s);
    }

    inline CRaritaSchwinger& operator+=(const CRaritaSchwinger &rs)
    {
      if (rs.m_b!=m_b || rs.m_r!=m_r) THROW(fatal_error, "Summing particle and anti-particle wave functions / bar "
                                                           "and non-bar wave functions is not allowed!")
      for (size_t i(0); i<16; ++i) m_x[i] += rs[i];
      m_on|=rs.m_on;
      return *this;
    }
    inline CRaritaSchwinger& operator-=(const CRaritaSchwinger &rs)
    {
      if (rs.m_b!=m_b || rs.m_r!=m_r) THROW(fatal_error, "Summing particle and anti-particle wave functions / bar "
                                                           "and non-bar wave functions is not allowed!")
      for (size_t i(0); i<16; ++i) m_x[i] -= rs[i];
      m_on|=rs.m_on;
      return *this;
    }
    /*inline CVec4& operator*=(const SComplex &c)
    {
      m_x[0]*=c; m_x[1]*=c; m_x[2]*=c; m_x[3]*=c;
      return *this;
    }*/
  // TODO: Stimmt das, dass sich m_r und m_b hier so ändern? Was soll CConj eigentlich machen? Bei Spinoren werden auch
  //       noch Komponenten paarweise vertauscht!!!
    inline CRaritaSchwinger CConj() const
    {
      return CRaritaSchwinger(-m_r, -m_b, std::conj(m_x[0]),std::conj(m_x[1]), std::conj(m_x[2]),std::conj(m_x[3]),
                              std::conj(m_x[4]),std::conj(m_x[5]), std::conj(m_x[6]),std::conj(m_x[7]),
                              std::conj(m_x[8]),std::conj(m_x[9]), std::conj(m_x[10]),std::conj(m_x[10]),
                              std::conj(m_x[11]),std::conj(m_x[12]), std::conj(m_x[13]),std::conj(m_x[14]),
		                          m_c[0], m_c[1], m_h, m_s, m_on);
    }

    inline int R() const { return m_r; }
    inline int B() const { return m_b; }
    inline int On() const { return m_on; }
    inline void SetR(const int &r) { m_r=r; }
    inline void SetB(const int &b) { m_b=b; }
    bool SetOn();

    bool Nan() const;

    inline CRaritaSchwinger<Scalar> Bar() const {
      return CRaritaSchwinger<Scalar>(m_r, -m_b, std::conj(m_x[2]), std::conj(m_x[3]),
                                      std::conj(m_x[0]), std::conj(m_x[1]), std::conj(m_x[6]), std::conj(m_x[7]),
                                      std::conj(m_x[4]), std::conj(m_x[5]), std::conj(m_x[10]), std::conj(m_x[11]),
                                      std::conj(m_x[8]), std::conj(m_x[9]), std::conj(m_x[14]), std::conj(m_x[15]),
                                      std::conj(m_x[12]), std::conj(m_x[13]), m_c[0], m_c[1], m_h, m_s,
                                      (m_on&1)<<1|(m_on&2)>>1); }

    static void ResetAccu();

    inline static void   SetAccu(const double &accu) { s_accu=accu;   }
    inline static double Accu()                      { return s_accu; }

    bool Test_Properties(const ATOOLS::Vec4D &p, int r, int b, int dir);

  };

  // end of class CRaritaSchwinger

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

  template <class Scalar>
  std::ostream &operator<<(std::ostream &s,const CRaritaSchwinger<Scalar> &vec);
}// end of namespace ATOOLS

#define DCRaritaSchwinger METOOLS::CRaritaSchwinger<double>
#define QCRaritaSchwinger METOOLS::CRaritaSchwinger<long double>

#endif
