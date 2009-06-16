#ifndef ATOOLS_Math_Lorentz_Ten3_H
#define ATOOLS_Math_Lorentz_Ten3_H

#include "ATOOLS/Math/MathTools.H"

using namespace ATOOLS;

namespace ATOOLS {

  template<typename Scalar> class Lorentz_Ten2;

  template<typename Scalar>
  class Lorentz_Ten3 {
    Scalar m_x[4][4][4];
  public:
    inline Lorentz_Ten3() {
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          for (unsigned short int k=0; k<4; ++k)
            m_x[i][j][k] = Scalar(0.0);
    }
    template<typename Scalar2>
    inline Lorentz_Ten3(const Lorentz_Ten3<Scalar2> & ten) {
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          for (unsigned short int k=0; k<4; ++k)
            m_x[i][j][k] = ten.at(i,j,k);
    }
    inline Lorentz_Ten3(
          const Scalar& x000, const Scalar& x100, const Scalar& x200, const Scalar& x300,
          const Scalar& x010, const Scalar& x110, const Scalar& x210, const Scalar& x310,
          const Scalar& x020, const Scalar& x120, const Scalar& x220, const Scalar& x320,
          const Scalar& x030, const Scalar& x130, const Scalar& x230, const Scalar& x330,

          const Scalar& x001, const Scalar& x101, const Scalar& x201, const Scalar& x301,
          const Scalar& x011, const Scalar& x111, const Scalar& x211, const Scalar& x311,
          const Scalar& x021, const Scalar& x121, const Scalar& x221, const Scalar& x321,
          const Scalar& x031, const Scalar& x131, const Scalar& x231, const Scalar& x331,

          const Scalar& x002, const Scalar& x102, const Scalar& x202, const Scalar& x302,
          const Scalar& x012, const Scalar& x112, const Scalar& x212, const Scalar& x312,
          const Scalar& x022, const Scalar& x122, const Scalar& x222, const Scalar& x322,
          const Scalar& x032, const Scalar& x132, const Scalar& x232, const Scalar& x332,

          const Scalar& x003, const Scalar& x103, const Scalar& x203, const Scalar& x303,
          const Scalar& x013, const Scalar& x113, const Scalar& x213, const Scalar& x313,
          const Scalar& x023, const Scalar& x123, const Scalar& x223, const Scalar& x323,
          const Scalar& x033, const Scalar& x133, const Scalar& x233, const Scalar& x333
            ) {
      m_x[0][0][0]=x000; m_x[1][0][0]=x100; m_x[2][0][0]=x200; m_x[3][0][0]=x300;
      m_x[0][1][0]=x010; m_x[1][1][0]=x110; m_x[2][1][0]=x210; m_x[3][1][0]=x310;
      m_x[0][2][0]=x020; m_x[1][2][0]=x120; m_x[2][2][0]=x220; m_x[3][2][0]=x320;
      m_x[0][3][0]=x030; m_x[1][3][0]=x130; m_x[2][3][0]=x230; m_x[3][3][0]=x330;

      m_x[0][0][1]=x001; m_x[1][0][1]=x101; m_x[2][0][1]=x201; m_x[3][0][1]=x301;
      m_x[0][1][1]=x011; m_x[1][1][1]=x111; m_x[2][1][1]=x211; m_x[3][1][1]=x311;
      m_x[0][2][1]=x021; m_x[1][2][1]=x121; m_x[2][2][1]=x221; m_x[3][2][1]=x321;
      m_x[0][3][1]=x031; m_x[1][3][1]=x131; m_x[2][3][1]=x231; m_x[3][3][1]=x331;

      m_x[0][0][2]=x002; m_x[1][0][2]=x102; m_x[2][0][2]=x202; m_x[3][0][2]=x302;
      m_x[0][1][2]=x012; m_x[1][1][2]=x112; m_x[2][1][2]=x212; m_x[3][1][2]=x312;
      m_x[0][2][2]=x022; m_x[1][2][2]=x122; m_x[2][2][2]=x222; m_x[3][2][2]=x322;
      m_x[0][3][2]=x032; m_x[1][3][2]=x132; m_x[2][3][2]=x232; m_x[3][3][2]=x332;

      m_x[0][0][3]=x003; m_x[1][0][3]=x103; m_x[2][0][3]=x203; m_x[3][0][3]=x303;
      m_x[0][1][3]=x013; m_x[1][1][3]=x113; m_x[2][1][3]=x213; m_x[3][1][3]=x313;
      m_x[0][2][3]=x023; m_x[1][2][3]=x123; m_x[2][2][3]=x223; m_x[3][2][3]=x323;
      m_x[0][3][3]=x033; m_x[1][3][3]=x133; m_x[2][3][3]=x233; m_x[3][3][3]=x333;
    }

    inline Lorentz_Ten3(const Scalar ten[4][4][4]) {
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          for (unsigned short int k=0; k<4; ++k)
            m_x[i][j][k] = ten[i][j][k];
    }

    // element extraction
    inline const Scalar at(unsigned short int i, unsigned short int j,
                           unsigned short int k) const {
      return m_x[i][j][k];
    }

    inline Scalar at(unsigned short int i, unsigned short int j,
                           unsigned short int k) {
      return m_x[i][j][k];
    }

    // sign flip
    inline Lorentz_Ten3<Scalar> operator-() const {
      Scalar x[4][4][4];
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          for (unsigned short int k=0; k<4; ++k)
            x[i][j][k] = -m_x[i][j][k];
      return Lorentz_Ten3<Scalar>(x);
    }

//     inline Lorentz_Ten3<Scalar> operator-() const {
//       for (unsigned short int i=0; i<4; ++i)
//         for (unsigned short int j=0; i<4; ++i)
//           for (unsigned short int k=0; i<4; ++i)
//             m_x[i][j][k] *= -1.;
//       return *this;
//     }

    // transpose pair of indizes
    inline Lorentz_Ten3<Scalar> Transpose(unsigned short int a, unsigned short int b) const {
      Scalar x[4][4][4];
      if      (((a==1) && (b==2)) || ((a==2) && (b==1)))
        for (unsigned short int i=0; i<4; ++i)
          for (unsigned short int j=0; j<4; ++j)
            for (unsigned short int k=0; k<4; ++k)
                x[i][j][k] = m_x[j][i][k];
      else if (((a==1) && (b==3)) || ((a==3) && (b==1)))
        for (unsigned short int i=0; i<4; ++i)
          for (unsigned short int j=0; j<4; ++j)
            for (unsigned short int k=0; k<4; ++k)
                x[i][j][k] = m_x[k][j][i];
      else if (((a==2) && (b==3)) || ((a==3) && (b==2)))
        for (unsigned short int i=0; i<4; ++i)
          for (unsigned short int j=0; j<4; ++j)
            for (unsigned short int k=0; k<4; ++k)
                x[i][j][k] = m_x[i][k][j];
      return Lorentz_Ten3<Scalar>(x);
    }

    // addition/subtraction operators
    inline Lorentz_Ten3<Scalar>& operator+=(const Lorentz_Ten3<Scalar>& t) {
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          for (unsigned short int k=0; k<4; ++k)
            m_x[i][j][k] += t.at(i,j,k);
      return *this;
    }

    inline Lorentz_Ten3<Scalar>& operator-=(const Lorentz_Ten3<Scalar>& t) {
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          for (unsigned short int k=0; k<4; ++k)
            m_x[i][j][k] -= t.at(i,j,k);
      return *this;
    }

    template<typename Scalar2> 
    inline Lorentz_Ten3<PROMOTE(Scalar,Scalar2)>
    operator+ (const Lorentz_Ten3<Scalar2>& ten) const {
      PROMOTE(Scalar,Scalar2) x[4][4][4];
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          for (unsigned short int k=0; k<4; ++k)
            x[i][j][k] = m_x[i][j][k]+ten.at(i,j,k);
      return Lorentz_Ten3<PROMOTE(Scalar,Scalar2)>(x);
    }

    template<typename Scalar2> 
    inline Lorentz_Ten3<PROMOTE(Scalar,Scalar2)>
    operator- (const Lorentz_Ten3<Scalar2>& ten) const {
      PROMOTE(Scalar,Scalar2) x[4][4][4];
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          for (unsigned short int k=0; k<4; ++k)
            x[i][j][k] = m_x[i][j][k]-ten.at(i,j,k);
      return Lorentz_Ten3<PROMOTE(Scalar,Scalar2)>(x);
    }

    // multiplication operators
    inline Lorentz_Ten3<Scalar>& operator*= (const Scalar& scal) {
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          for (unsigned short int k=0; k<4; ++k)
            m_x[i][j][k] *= scal;
      return *this;
    }

    template<typename Scalar2> 
    inline Lorentz_Ten3<PROMOTE(Scalar,Scalar2)> operator* (const Scalar2 scal) const {
      PROMOTE(Scalar,Scalar2) x[4][4][4];
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          for (unsigned short int k=0; k<4; ++k)
            x[i][j][k] = scal*m_x[i][j][k];
      return Lorentz_Ten3<PROMOTE(Scalar,Scalar2)>(x);
    }

    template<typename Scalar2> 
    inline Lorentz_Ten3<PROMOTE(Scalar,Scalar2)> operator/ (const Scalar2 scal) const {
      return (*this)*(1./scal);
    }

    inline bool Nan() const {
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          for (unsigned short int k=0; k<4; ++k)
            if (ATOOLS::IsNan<Scalar>(m_x[i][j][k])) return true;
      return false;
    }

    inline bool IsZero() const {
      for (unsigned short int i=0; i<4; ++i)
        for (unsigned short int j=0; j<4; ++j)
          for (unsigned short int k=0; k<4; ++k)
            if (!ATOOLS::IsZero<Scalar>(m_x[i][j][k])) return false;
      return true;
    }
  };

  template<typename Scalar,typename Scal2> 
  inline Lorentz_Ten3<PROMOTE(Scalar,Scal2)> 
  operator* (const Scal2 s, const Lorentz_Ten3<Scalar>& t) {
    return Lorentz_Ten3<PROMOTE(Scalar,Scal2)>(t*s);
  }

  template<typename Scalar>
  std::ostream& operator<<(std::ostream& s, const Lorentz_Ten3<Scalar>& ten) {
//     return s<<"not implemented for tensors of third rank ... "
//             <<"mainly due to insufficient depth of the screen ... "
//             <<"please project to second rank and display then ..."<<std::endl;
    return s<<'('<<ten.at(0,0,0)<<','<<ten.at(1,0,0)<<','<<ten.at(2,0,0)<<','<<ten.at(3,0,0)<<' '<<std::endl
            <<' '<<ten.at(0,1,0)<<','<<ten.at(1,1,0)<<','<<ten.at(2,1,0)<<','<<ten.at(3,1,0)<<' '<<std::endl
            <<' '<<ten.at(0,2,0)<<','<<ten.at(1,2,0)<<','<<ten.at(2,2,0)<<','<<ten.at(3,2,0)<<' '<<std::endl
            <<' '<<ten.at(0,3,0)<<','<<ten.at(1,3,0)<<','<<ten.at(2,3,0)<<','<<ten.at(3,3,0)<<')'<<std::endl
            <<'('<<ten.at(0,0,1)<<','<<ten.at(1,0,1)<<','<<ten.at(2,0,1)<<','<<ten.at(3,0,1)<<' '<<std::endl
            <<' '<<ten.at(0,1,1)<<','<<ten.at(1,1,1)<<','<<ten.at(2,1,1)<<','<<ten.at(3,1,1)<<' '<<std::endl
            <<' '<<ten.at(0,2,1)<<','<<ten.at(1,2,1)<<','<<ten.at(2,2,1)<<','<<ten.at(3,2,1)<<' '<<std::endl
            <<' '<<ten.at(0,3,1)<<','<<ten.at(1,3,1)<<','<<ten.at(2,3,1)<<','<<ten.at(3,3,1)<<')'<<std::endl
            <<'('<<ten.at(0,0,2)<<','<<ten.at(1,0,2)<<','<<ten.at(2,0,2)<<','<<ten.at(3,0,2)<<' '<<std::endl
            <<' '<<ten.at(0,1,2)<<','<<ten.at(1,1,2)<<','<<ten.at(2,1,2)<<','<<ten.at(3,1,2)<<' '<<std::endl
            <<' '<<ten.at(0,2,2)<<','<<ten.at(1,2,2)<<','<<ten.at(2,2,2)<<','<<ten.at(3,2,2)<<' '<<std::endl
            <<' '<<ten.at(0,3,2)<<','<<ten.at(1,3,2)<<','<<ten.at(2,3,2)<<','<<ten.at(3,3,2)<<')'<<std::endl
            <<'('<<ten.at(0,0,3)<<','<<ten.at(1,0,3)<<','<<ten.at(2,0,3)<<','<<ten.at(3,0,3)<<' '<<std::endl
            <<' '<<ten.at(0,1,3)<<','<<ten.at(1,1,3)<<','<<ten.at(2,1,3)<<','<<ten.at(3,1,3)<<' '<<std::endl
            <<' '<<ten.at(0,2,3)<<','<<ten.at(1,2,3)<<','<<ten.at(2,2,3)<<','<<ten.at(3,2,3)<<' '<<std::endl
            <<' '<<ten.at(0,3,3)<<','<<ten.at(1,3,3)<<','<<ten.at(2,3,3)<<','<<ten.at(3,3,3)<<')'<<std::endl;
  }
}

/*!
 \file
 \brief   contains class Lorentz_Ten3
*/

/*!
 \class Lorentz_Ten3<Scalar>
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
 \fn Lorentz_Ten3::Lorentz_Ten3()
 \brief Standard Constructor
*/

/*!
 \fn Lorentz_Ten3::Lorentz_Ten3(Scalar x00, Scalar x10, Scalar x20, Scalar x30, Scalar x01, Scalar x11, Scalar x21, Scalar x31, Scalar x02, Scalar x12, Scalar x22, Scalar x32, Scalar x03, Scalar x13, Scalar x23, Scalar x33){
 \brief Special Constructor, templated in Scalar, taking 16 single components.
*/

/*!
 \fn   inline Scalar* Lorentz_Ten3::operator[] (int i)
 \brief returns array \f$x_{ij}\f$, \f$j\in[0,3]\f$. (May be manipulated.)
*/

/*!
 \fn   inline const Scalar* Lorentz_Ten3::operator[] (int i) const
 \brief returns array \f$x_{ij}\f$, \f$j\in[0,3]\f$.
*/

#endif
