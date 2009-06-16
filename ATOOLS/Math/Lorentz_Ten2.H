#ifndef ATOOLS_Math_Lorentz_Ten2_H
#define ATOOLS_Math_Lorentz_Ten2_H

#include "ATOOLS/Math/MathTools.H"

using namespace ATOOLS;

namespace ATOOLS {
  template<typename Scalar> class Vec4;

  template<typename Scalar>
  class Lorentz_Ten2 {
    Scalar m_x[4][4];
  public:
    inline Lorentz_Ten2() {
      m_x[0][0]=m_x[1][0]=m_x[2][0]=m_x[3][0]=
      m_x[0][1]=m_x[1][1]=m_x[2][1]=m_x[3][1]=
      m_x[0][2]=m_x[1][2]=m_x[2][2]=m_x[3][2]=
      m_x[0][3]=m_x[1][3]=m_x[2][3]=m_x[3][3]=Scalar(0.0);
    }
    template<typename Scalar2>
    inline Lorentz_Ten2(const Lorentz_Ten2<Scalar2>& ten) {
      m_x[0][0]=ten.at(0,0); m_x[1][0]=ten.at(1,0); m_x[2][0]=ten.at(2,0); m_x[3][0]=ten.at(3,0);
      m_x[0][1]=ten.at(0,1); m_x[1][1]=ten.at(1,1); m_x[2][1]=ten.at(2,1); m_x[3][1]=ten.at(3,1);
      m_x[0][2]=ten.at(0,2); m_x[1][2]=ten.at(1,2); m_x[2][2]=ten.at(2,2); m_x[3][2]=ten.at(3,2);
      m_x[0][3]=ten.at(0,3); m_x[1][3]=ten.at(1,3); m_x[2][3]=ten.at(2,3); m_x[3][3]=ten.at(3,3);
    }
    inline Lorentz_Ten2(const Scalar& x00, const Scalar& x10, const Scalar& x20, const Scalar& x30,
                const Scalar& x01, const Scalar& x11, const Scalar& x21, const Scalar& x31,
                const Scalar& x02, const Scalar& x12, const Scalar& x22, const Scalar& x32,
                const Scalar& x03, const Scalar& x13, const Scalar& x23, const Scalar& x33) {
      m_x[0][0] = x00; m_x[1][0] = x10; m_x[2][0] = x20; m_x[3][0] = x30;
      m_x[0][1] = x01; m_x[1][1] = x11; m_x[2][1] = x21; m_x[3][1] = x31;
      m_x[0][2] = x02; m_x[1][2] = x12; m_x[2][2] = x22; m_x[3][2] = x32;
      m_x[0][3] = x03; m_x[1][3] = x13; m_x[2][3] = x23; m_x[3][3] = x33;
    }
    inline Lorentz_Ten2(const Scalar ten[4][4]) {
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          m_x[i][j] = ten[i][j];
    }

    // element extraction
    inline const Scalar at(unsigned short int i, unsigned short int j) const {
      return m_x[i][j];
    }

    inline Scalar at(unsigned short int i, unsigned short int j) {
      return m_x[i][j];
    }

    // sign flip
    inline Lorentz_Ten2<Scalar> operator-() const {
      return Lorentz_Ten2<Scalar>(-m_x[0][0],-m_x[1][0],-m_x[2][0],-m_x[3][0],
                          -m_x[0][1],-m_x[1][1],-m_x[2][1],-m_x[3][1],
                          -m_x[0][2],-m_x[1][2],-m_x[2][2],-m_x[3][2],
                          -m_x[0][3],-m_x[1][3],-m_x[2][3],-m_x[3][3]);
    }

    // transpose
    inline Lorentz_Ten2<Scalar> Transpose() const {
      return Lorentz_Ten2<Scalar>(m_x[0][0],m_x[0][1],m_x[0][2],m_x[0][3],
                                  m_x[1][0],m_x[1][1],m_x[1][2],m_x[1][3],
                                  m_x[2][0],m_x[2][1],m_x[2][2],m_x[2][3],
                                  m_x[3][0],m_x[3][1],m_x[3][2],m_x[3][3]);
    }

    // addition/subtraction operators
    inline Lorentz_Ten2<Scalar>& operator+=(const Lorentz_Ten2<Scalar>& t) {
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          m_x[i][j] += t.at(i,j);
      return *this;
    }

    inline Lorentz_Ten2<Scalar>& operator-=(const Lorentz_Ten2<Scalar>& t) {
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          m_x[i][j] -= t.at(i,j);
      return *this;
    }

    template<typename Scalar2> inline Lorentz_Ten2<PROMOTE(Scalar,Scalar2)>
    operator+ (const Lorentz_Ten2<Scalar2>& t) const {
      return Lorentz_Ten2<PROMOTE(Scalar,Scalar2)>
                (m_x[0][0]+t.at(0,0),m_x[1][0]+t.at(1,0),
                                      m_x[2][0]+t.at(2,0),m_x[3][0]+t.at(3,0),
                 m_x[0][1]+t.at(0,1),m_x[1][1]+t.at(1,1),
                                      m_x[2][1]+t.at(2,1),m_x[3][1]+t.at(3,1),
                 m_x[0][2]+t.at(0,2),m_x[1][2]+t.at(1,2),
                                      m_x[2][2]+t.at(2,2),m_x[3][2]+t.at(3,2),
                 m_x[0][3]+t.at(0,3),m_x[1][3]+t.at(1,3),
                                      m_x[2][3]+t.at(2,3),m_x[3][3]+t.at(3,3));
    }

    template<typename Scalar2> inline Lorentz_Ten2<PROMOTE(Scalar,Scalar2)>
    operator- (const Lorentz_Ten2<Scalar2>& t2) const {
      return Lorentz_Ten2<PROMOTE(Scalar,Scalar2)>
                (m_x[0][0]-t2.at(0,0),m_x[1][0]-t2.at(1,0),
                                      m_x[2][0]-t2.at(2,0),m_x[3][0]-t2.at(3,0),
                 m_x[0][1]-t2.at(0,1),m_x[1][1]-t2.at(1,1),
                                      m_x[2][1]-t2.at(2,1),m_x[3][1]-t2.at(3,1),
                 m_x[0][2]-t2.at(0,2),m_x[1][2]-t2.at(1,2),
                                      m_x[2][2]-t2.at(2,2),m_x[3][2]-t2.at(3,2),
                 m_x[0][3]-t2.at(0,3),m_x[1][3]-t2.at(1,3),
                                      m_x[2][3]-t2.at(2,3),m_x[3][3]-t2.at(3,3));
    }

    // multiplication operators
    inline Lorentz_Ten2<Scalar>& operator*= (const Scalar& scal) {
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          m_x[i][j] *= scal;
      return *this;
    }

    template<typename Scalar2> inline Lorentz_Ten2<PROMOTE(Scalar,Scalar2)>
    operator* (const Scalar2 scal) const {
      return Lorentz_Ten2<PROMOTE(Scalar,Scalar2)>
        (scal*m_x[0][0],scal*m_x[1][0],scal*m_x[2][0],scal*m_x[3][0],
         scal*m_x[0][1],scal*m_x[1][1],scal*m_x[2][1],scal*m_x[3][1],
         scal*m_x[0][2],scal*m_x[1][2],scal*m_x[2][2],scal*m_x[3][2],
         scal*m_x[0][3],scal*m_x[1][3],scal*m_x[2][3],scal*m_x[3][3]);
    }

    template<typename Scalar2> inline Lorentz_Ten2<PROMOTE(Scalar,Scalar2)>
    operator/ (const Scalar2 scal) const {
      return (*this)*(1./scal);
    }

    inline bool Nan() const {
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          if (ATOOLS::IsNan<Scalar>(m_x[i][j])) return true;
      return false;
    }

    inline bool IsZero() const {
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          if (!ATOOLS::IsZero<Scalar>(m_x[i][j])) return false;
      return true;
    }
  };

  template<typename Scalar,typename Scal2> inline Lorentz_Ten2<PROMOTE(Scalar,Scal2)>
  operator* (const Scal2 s, const Lorentz_Ten2<Scalar>& t) {
    return Lorentz_Ten2<PROMOTE(Scalar,Scal2)>(t*s);
  }

  template<typename Scalar>
  std::ostream& operator<<(std::ostream& s, const Lorentz_Ten2<Scalar>& ten) {
    return s<<'('<<ten.at(0,0)<<','<<ten.at(1,0)<<','<<ten.at(2,0)<<','<<ten.at(3,0)<<' '<<std::endl
            <<' '<<ten.at(0,1)<<','<<ten.at(1,1)<<','<<ten.at(2,1)<<','<<ten.at(3,1)<<' '<<std::endl
            <<' '<<ten.at(0,2)<<','<<ten.at(1,2)<<','<<ten.at(2,2)<<','<<ten.at(3,2)<<' '<<std::endl
            <<' '<<ten.at(0,3)<<','<<ten.at(1,3)<<','<<ten.at(2,3)<<','<<ten.at(3,3)<<')'<<std::endl;
  }
}

/*!
 \file
 \brief   contains class Lorentz_Ten2
*/

/*!
 \class Lorentz_Ten2<Scalar>
 \brief implementation of a 4 dimensional Minkowski 2nd rank tensor and its algebraic
 operations

 This class can be used as Minkowski 2nd rank tensor with arbitrary scalar types.
 All necessary operations, e.g. addition, scalar multiplication, contractions
 with other 2nd rank tensor or vectors, etc. are available.
 If you want to use this 2nd rank tensor for a new scalar type, you might have to
 implement specific functions for this type in MathTools.H.
 "Fallback" types for operations with two 2nd rank tensors of different types, e. g.
 "complex and double become complex" have to be specified in MathTools.H as
 well, see the DECLARE_PROMOTE macro.
*/

/*!
 \fn Lorentz_Ten2::Lorentz_Ten2()
 \brief Standard Constructor
*/

/*!
 \fn Lorentz_Ten2::Lorentz_Ten2(Scalar x00, Scalar x10, Scalar x20, Scalar x30, Scalar x01, Scalar x11, Scalar x21, Scalar x31, Scalar x02, Scalar x12, Scalar x22, Scalar x32, Scalar x03, Scalar x13, Scalar x23, Scalar x33){
 \brief Special Constructor, templated in Scalar, taking 16 single components.
*/

/*!
 \fn   inline Scalar* Lorentz_Ten2::operator[] (int i)
 \brief returns array \f$x_{ij}\f$, \f$j\in[0,3]\f$. (May be manipulated.)
*/

/*!
 \fn   inline const Scalar* Lorentz_Ten2::operator[] (int i) const
 \brief returns array \f$x_{ij}\f$, \f$j\in[0,3]\f$.
*/

#endif
