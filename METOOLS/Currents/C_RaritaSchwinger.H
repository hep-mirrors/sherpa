#ifndef METOOLS_Explicit_C_RaritaSchwinger_H
#define METOOLS_Explicit_C_RaritaSchwinger_H

#include "METOOLS/Explicit/C_Object.H"
#include "ATOOLS/Math/Vec4.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Matrix.H"
#include "ATOOLS/Org/Exception.H"
#include "METOOLS/Currents/C_Spinor.H"
#include <vector>

#define ZERO std::complex<Scalar>(0.0,0.0)

 // TODO: Majorana-3/2 Teilchen?

namespace METOOLS {

  template <class Scalar>
  class CRaritaSchwinger: public CObject {
  public:

    typedef std::complex<Scalar> SComplex;

  private:

    SComplex m_x[16];

    int m_r, m_b, m_on;

    static double s_accu;

    static ATOOLS::AutoDelete_Vector<CRaritaSchwinger> s_objects;

    template <class _Scalar> friend std::ostream &
    operator<<(std::ostream &s,const CRaritaSchwinger<_Scalar> &vec);

  public:

    static CRaritaSchwinger *New();
    static CRaritaSchwinger *New(const CRaritaSchwinger &s);
    static CRaritaSchwinger *New(const int r, const int b, const Scalar &x0, const Scalar &x1, const Scalar &x2,
                                 const Scalar &x3, const Scalar &x4, const Scalar &x5, const Scalar &x6,
                                 const Scalar &x7, const Scalar &x8, const Scalar &x9, const Scalar &x10,
                                 const Scalar &x11, const Scalar &x12, const Scalar &x13, const Scalar &x14,
                                 const Scalar &x15, const int c1=0,const int c2=0, const size_t &h=0,
                                 const size_t &s=0, const int on=3);

    CObject* Copy() const;

    void Delete();

    bool IsZero() const;

   inline CRaritaSchwinger(const int r=1, const int b=1, const int cr=0,const int ca=0,
                           const size_t &h=0,const size_t &s=0, const int &on=3)
    {
      for (auto entry : m_x){
        entry = Scalar(0.0);
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

   // operator return the component i of the RS wave function
   inline SComplex &operator[](const size_t &i) { return m_x[i]; }
   inline const SComplex &operator[](const size_t &i) const { return m_x[i]; }

   ATOOLS::TCMatrix<Scalar> Contract4Index(const CRaritaSchwinger<Scalar> &rs) const;
   ATOOLS::TCMatrix<Scalar> ContractSpinorIndex(const CRaritaSchwinger<Scalar> &rs) const;

   // operator multiplies two RS wave functions component-wise and returns the result as a new RS wave function object
   SComplex operator*(const CRaritaSchwinger<Scalar> &rs) const;
   // operator adds two RS wave functions together and returns the result as a new RS wave function object
   inline CRaritaSchwinger operator+(const CRaritaSchwinger &rs) const
    {
      if (rs.m_b!=m_b || rs.m_r!=m_r) THROW(fatal_error, "Summing particle and anti-particle wave functions / bar "
                                                           "and non-bar wave functions is not allowed!")
      return CRaritaSchwinger(m_r, m_b, m_x[0]+rs[0], m_x[1]+rs[1], m_x[2]+rs[2], m_x[3]+rs[3], m_x[4]+rs[4],
                              m_x[5]+rs[5], m_x[6]+rs[6], m_x[7]+rs[7], m_x[8]+rs[8], m_x[9]+rs[9], m_x[10]+rs[10],
                              m_x[11]+rs[11], m_x[12]+rs[12], m_x[13]+rs[13], m_x[14]+rs[14], m_x[15]+rs[15],
                              m_c[0], m_c[1], m_h, m_s, m_on|rs.m_on);
    }
    // operator subtracts two RS wave functions from each other and returns the result as a new RS wave function object
    inline CRaritaSchwinger operator-(const CRaritaSchwinger &rs) const
    {
      if (rs.m_b!=m_b || rs.m_r!=m_r) THROW(fatal_error, "Summing particle and anti-particle wave functions / bar "
                                                           "and non-bar wave functions is not allowed!")
      return CRaritaSchwinger(m_r, m_b, m_x[0]-rs[0], m_x[1]-rs[1], m_x[2]-rs[2], m_x[3]-rs[3], m_x[4]-rs[4],
                              m_x[5]-rs[5], m_x[6]-rs[6], m_x[7]-rs[7], m_x[8]-rs[8], m_x[9]-rs[9], m_x[10]-rs[10],
                              m_x[11]-rs[11], m_x[12]-rs[12], m_x[13]-rs[13], m_x[14]-rs[14], m_x[15]-rs[15],
                              m_c[0], m_c[1], m_h, m_s, m_on|rs.m_on);
    }
    // operator multiplies the (*this) with a factor (-1)
    inline CRaritaSchwinger operator-() const
    {
      return CRaritaSchwinger(m_r, m_b, -m_x[0], -m_x[1], -m_x[2], -m_x[3], -m_x[4], -m_x[5], -m_x[6], -m_x[7],
                              -m_x[8], -m_x[9], -m_x[10],-m_x[11], -m_x[12], -m_x[13], -m_x[14], -m_x[15],
                              m_c[0], m_c[1], m_h, m_s);
    }
    // operator adds a RS wave function to (*this)
    inline CRaritaSchwinger& operator+=(const CRaritaSchwinger &rs)
    {
      if (rs.m_b!=m_b || rs.m_r!=m_r) THROW(fatal_error, "Summing particle and anti-particle wave functions / bar "
                                                           "and non-bar wave functions is not allowed!")
      for (size_t i(0); i<16; ++i) m_x[i] += rs[i];
      m_on|=rs.m_on;
      return *this;
    }
    // operator subtracts a RS wave function from (*this)
    inline CRaritaSchwinger& operator-=(const CRaritaSchwinger &rs)
    {
      if (rs.m_b!=m_b || rs.m_r!=m_r) THROW(fatal_error, "Summing particle and anti-particle wave functions / bar "
                                                           "and non-bar wave functions is not allowed!")
      for (size_t i(0); i<16; ++i) m_x[i] -= rs[i];
      m_on|=rs.m_on;
      return *this;
    }
    inline CRaritaSchwinger& operator*=(const SComplex &c)
    {
      if (m_on&1) {
        m_x[0]*=c;
        m_x[1]*=c;
        m_x[4]*=c;
        m_x[5]*=c;
        m_x[8]*=c;
        m_x[9]*=c;
        m_x[12]*=c;
        m_x[13]*=c;
      }
      if (m_on&2) {
        m_x[2]*=c;
        m_x[3]*=c;
        m_x[6]*=c;
        m_x[7]*=c;
        m_x[10]*=c;
        m_x[11]*=c;
        m_x[14]*=c;
        m_x[15]*=c;
      }
      return *this;
    }
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
    // Getter of particle/anti-particle parameter (m_r), adjoint wave function parameter (m_b) and m_on
    inline int R() const { return m_r; }
    inline int B() const { return m_b; }
    inline int On() const { return m_on; }
    // Setter for particle/anti-particle parameter (m_r), adjoint wave function parameter (m_b) and m_on
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

    bool Test_Properties(const ATOOLS::Vec4D &p, int dir);

  };

  // end of class CRaritaSchwinger
  // operators for multiplying a scalar or a complex number to a Rarita-Schwinger wavefunction from left or right
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

  // Output operator, enabling std::cout for Rarita-Schwinger wave functions
  template <class Scalar>
  std::ostream &operator<<(std::ostream &s,const CRaritaSchwinger<Scalar> &vec);
}// end of namespace ATOOLS

#define DCRaritaSchwinger METOOLS::CRaritaSchwinger<double>
#define QCRaritaSchwinger METOOLS::CRaritaSchwinger<long double>

#endif

/*!
   \class CRaritaSchwinger
   \brief Class whose objects are wave functions of spin 3/2 particles

     This class describes the wave functions of spin 3/2 particles and corresponding operations which can be applied
     on it. It inherits from the general CObject class. The wave function is based on the Rarita-Schwinger
     representation describing spin 3/2 particles as an object with a 4-Spinor index and a Lorentz-index resulting in
     16 components for the wave function. The components are calculated in RaSc.C.
   */
/*!
\var SComplex m_x
\brief C++ array that stores the 16 components of the Rarita-Schwinger wave function
*/
/*!
\var int m_r
\brief Integer characterizes whether the current object is the wave function of a particle (m_r=1) or of an anti-
       particle (-1). Depending on that value, particle or antiparticle spinors were used to construct the wave
       function.
*/
/*!
\var int m_b
\brief Integer denotes whether spinors or adjoint spinors were used to calculate the wave function.
*/
/*!
\var int m_on
\brief Integer characterizes the type of the wave function:
       m_on = 1: massless right-handed particles
       m_on = 2: massless left-handed particles
       m_on = 3: massive particles
       For the first two half of the 16 components are exactly ZERO, respectively which can also be ensured numerically
       in later steps of the calcuation by checking this parameter.
*/
/*!
\var static double s_acc
\brief Contains the accuracy with which the wave function can be calculated, double and long double precision are
       currently supported.
*/
/*!
\var static ATOOLS::AutoDelete_Vector<CRaritaSchwinger> s_objects
\brief
*/

/*!
 \fn CRaritaSchwinger::New()
\brief New method for Rarita-Schwinger wave function, generate pointer to an empty Rarita-Schwinger wave function
Input:
  - no input
    .
Output:
  - CRaritaSchwinger *: Pointer to an empty Rarita-Schwinger wave function
    .
*/
/*!
 \fn CRaritaSchwinger::New(const CRaritaSchwinger &s)
\brief New method for Rarita-Schwinger wave function, generate pointer to an Rarita-Schwinger wave function identical
       to a wave function given as input
Input:
  - <var>s</var>: Rarita-Schwinger wave function
    .
Output:
  - CRaritaSchwinger *: Pointer to a Rarita-Schwinger wave function identical to the given input
   .
*/
/*!
\fn CRaritaSchwinger::New(const int r, const int b, const Scalar &x0, const Scalar &x1, const Scalar &x2,
                           const Scalar &x3, const Scalar &x4, const Scalar &x5, const Scalar &x6,
                           const Scalar &x7, const Scalar &x8, const Scalar &x9, const Scalar &x10,
                           const Scalar &x11, const Scalar &x12, const Scalar &x13, const Scalar &x14,
                           const Scalar &x15, const int c1=0,const int c2=0, const size_t &h=0,
                           const size_t &s=0, const int on=3)
\brief New method for Rarita-Schwinger wave function, generate pointer to a Rarita-Schwinger wave function with the
 parameter values given as input
  Input:
- <var>r</var>: desired value for m_r (particle / antiparticle)
- <var>b</var>: desired value for m_b (spinor / adjoint spinor)
- <var>x0</var>-<var>x15</var>: desired values for the 16 components of the Rarita-Schwinger wave function
- optional <var>c1</var>: color factor, desired value for m_c[0] of CObject (default : 0)
- optional <var>c2</var>: color factor, desired value for m_c[1] of CObject (default : 0)
- optional <var>h</var>: helicity number, desired value for m_h of CObject (default : 0)
- optional <var>s</var>: subtraction, desired value for m_s of CObject (default : 0)
- optional <var>on</var>: massless right-(m_on=1)/left(m_on=2)-handed or massive (m_on=3) particle, desired value for
                          m_on (default m_on=3)

.
Output:
- CRaritaSchwinger *: Pointer to a Rarita-Schwinger wave function identical to a Rarita-Schwinger wave function with the
 parameter values given as input
     .
*/
/*!
 \fn CRaritaSchwinger::Copy() const
\brief New method for Rarita-Schwinger wave function, generate pointer to an Rarita-Schwinger wave function identical
       to a wave function given as input
Input:
  - <var>s</var>: Rarita-Schwinger wave function
    .
Output:
  - CRaritaSchwinger *: Pointer to a Rarita-Schwinger wave function identical to the given input
    .
*/
/*!
 \fn CRaritaSchwinger::Delete()
\brief Delete method for Rarita-Schwinger wave function
Input:
  - no input
    .
Output:
  - no input
    .
*/
/*!
 \fn CRaritaSchwinger::IsZero() const
\brief Method tests whether the Rarita-Schwinger wave function only contains zeros
Input:
  - no input
    .
Output:
  - bool: whether all components are zero or not
   .
*/
/*!
 \fn CRaritaSchwinger::CRaritaSchwinger(const int r=1, const int b=1, const int cr=0,const int ca=0,
                                        const size_t &h=0,const size_t &s=0, const int &on=3)
\brief Constructor of a Rarita-Schwinger wave function filled with zeros
Input:
  - optional <var>r</var>: desired value for m_r (particle / antiparticle)
  - optional <var>b</var>: desired value for m_b (spinor / adjoint spinor)
  - optional <var>c1</var>: color factor, desired value for m_c[0] of CObject (default : 0)
  - optional <var>c2</var>: color factor, desired value for m_c[1] of CObject (default : 0)
  - optional <var>h</var>: helicity number, desired value for m_h of CObject (default : 0)
  - optional <var>s</var>: subtraction, desired value for m_s of CObject (default : 0)
  - optional <var>on</var>: massless right-(m_on=1)/left(m_on=2)-handed or massive (m_on=3) particle, desired value for
                            m_on (default m_on=3)
    .
*/
/*!
 \fn CRaritaSchwinger::CRaritaSchwinger(const CRaritaSchwinger &v)
\brief Constructor of a Rarita-Schwinger wave function identical to the provided one
Input:
  - <var>v</var>: CRaritaSchwinger which should be used to construct a new Rarita-Schwinger wave function
    .
*/
/*!
 \fn CRaritaSchwinger::CRaritaSchwinger(const int &r, const int &b, const SComplex &x0, const SComplex &x1, const SComplex &x2,
                            const SComplex &x3, const SComplex &x4, const SComplex &x5, const SComplex &x6,
                            const SComplex &x7, const SComplex &x8, const SComplex &x9, const SComplex &x10,
                            const SComplex &x11, const SComplex &x12, const SComplex &x13, const SComplex &x14,
                            const SComplex &x15, const int c1=0,const int c2=0, const size_t &h=0,const size_t &s=0,
                            const int &on=3)
\brief Constructor of a Rarita-Schwinger wave function, all parameters can be provided via input values
Input:
  - <var>r</var>: desired value for m_r (particle / antiparticle)
  - <var>b</var>: desired value for m_b (spinor / adjoint spinor)
  - <var>x0</var>-<var>x15</var>: desired values for the 16 components of the Rarita-Schwinger wave function
  - optional <var>c1</var>: color factor, desired value for m_c[0] of CObject (default : 0)
  - optional <var>c2</var>: color factor, desired value for m_c[1] of CObject (default : 0)
  - optional <var>h</var>: helicity number, desired value for m_h of CObject (default : 0)
  - optional <var>s</var>: subtraction, desired value for m_s of CObject (default : 0)
  - optional <var>on</var>: massless right-(m_on=1)/left(m_on=2)-handed or massive (m_on=3) particle, desired value for
                            m_on (default m_on=3)
    .
*/
/*!
 \fn CRaritaSchwinger::CRaritaSchwinger(const CRaritaSchwinger &rs, const Scalar &c)
\brief Constructor of a Rarita-Schwinger wave function identical to the provided one multiplied by a Scalar
Input:
  - <var>v</var>: CRaritaSchwinger which should be used to construct a new Rarita-Schwinger wave function
  - <var>c</var>: Scalar to multiply with the input Rarita-Schwinger wave function to receive the desired one
    .
*/
/*!
 \fn CRaritaSchwinger::CRaritaSchwinger(const CRaritaSchwinger &rs, const SComplex &c)
\brief Constructor of a Rarita-Schwinger wave function identical to the provided one multiplied by a complex number
Input:
  - <var>v</var>: CRaritaSchwinger which should be used to construct a new Rarita-Schwinger wave function
  - <var>c</var>: Complex number to multiply with the input Rarita-Schwinger wave function to receive the desired one
    .
*/
/*!
 \fn CRaritaSchwinger::Add(const CObject *c)
\brief Method adds one Rarita-Schwinger wave function to (*this)
Input:
  - <var>c</var>: pointer to the Rarita-Schwinger wave function which should be added to (*this)
    .
Output:
  - no output
    .
*/
/*!
 \fn CRaritaSchwinger::Divide(const double &d)
\brief Method divides (*this) by a factor
Input:
  - <var>d</var>: factor, (*this) should be divided with
    .
Output:
  - no output
    .
*/
/*!
 \fn CRaritaSchwinger::Multiply(const Complex &d)
\brief Method multiplies (*this) by a factor
Input:
  - <var>d</var>: factor, (*this) should be multiplied with
    .
Output:
  - no output
    .
*/
/*!
 \fn CRaritaSchwinger::Invert()
\brief Method multiplies (*this) by (-1)
Input:
  - no input
    .
Output:
  - no output
    .
*/
/*!
 \fn CRaritaSchwinger::Contract4Index(const CRaritaSchwinger<Scalar> &rs) const
\brief Method sums over the Lorentz indices (*this) and the input one
Input:
  - <var>rs</var>: Rarita-Schwinger wave function, (*this) should be contracted with
    .
Output:
  - ATOOLS::TCMatrix<Scalar>: 4x4 matrix, containing the result of the contraction, indices of the matrix are the
                              remaining spinor indices of the contracted wave functions
    .
*/
/*!
 \fn CRaritaSchwinger::ContractSpinorIndex(const CRaritaSchwinger<Scalar> &rs) const
\brief Method sums over the spinor indices of the (*this) and the input Rarita-Schwinger wave function
Input:
  - <var>rs</var>: Rarita-Schwinger wave function, (*this) should be contracted with
    .
Output:
  - ATOOLS::TCMatrix<Scalar>: 4x4 matrix, containing the result of the contraction, indices of the matrix are the
                              remaining Lorentz indices of the contracted wave functions
    .
*/
/*!
 \fn CRaritaSchwinger::Nan() const
\brief Method checks whether at least one entries of (*this) is nan
Input:
  - no input
    .
Output:
  - bool: true, if at least one entry of (*this) is none, false otherwise
  .
*/
/*!
\fn CRaritaSchwinger::Bar() const
\brief Method  calculates the adjoint Rarita-Schwinger wave functions for (*this)
Input:
  - no input
    .
Output:
  - CRaritaSchwinger<Scalar>: adjoint Rarita-Schwinger wave function of (*this)
    .
*/
/*!
\fn CRaritaSchwinger::ResetAccu()
\brief Resets numerical precision of (*this) to the default value
Input:
  - no input
    .
Output:
  - no output
    .
*/
/*!
\fn CRaritaSchwinger::Test_Properties(const ATOOLS::Vec4D &p, int r, int b, int dir)
\brief Tests the properties of (*this), i.e. irreducibility conditions (gamma_mu*rs=0 and p_mu*rs=0) and the Dirac
       equation
Input:
  - <var>p</var>: momentum of the particle corresponding to (*this)
  - <var>dir</var>: direction parameter, i.e. whether the particles corresponding to (*this) is originally incoming
    (dir>0) or outgoing (dir<0)
    .
Output:
  - bool: true, if all tests are passed, false otherwise
    .
*/
