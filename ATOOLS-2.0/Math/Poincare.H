#ifndef Poincare_H
#define Poincare_H

#include "Vector.H"

namespace ATOOLS {

  class Poincare {
  private:

    int    m_status;               // 0 no operation, 1 boost, 2 rot, 3 boost+rot
    double m_mat[4][4];
    Vec4D  m_beta;
    double m_rsq;

  public:

    Poincare();                    // standard constructor Unity
    Poincare(Vec4D v1, Vec4D v2);  // rotations constructor
    Poincare(Vec4D v);             // boost constuctor 
 
    void Boost(Vec4D& v);          // boosts vectors in CMS (of vector given at initialization)
    void BoostBack(Vec4D& v);      // boost back to LAB Frame

    void Rotate(Vec4D& v);         // rotate
    void RotateBack(Vec4D& v);     // rotate back

    bool CheckBoost();             // test boost 
    bool CheckRotation();          // test rotation

    void Invert();  

    inline Vec4D operator*(const Vec4D& vin) {
      Vec4D v(vin);
      if ((m_status)&(1)) Boost(v);
      if ((m_status)&(2)) Rotate(v);
      return v;
    }

  };

  /*!
    \file
    \brief  contains the class Poincare
  */

  /*!
    \class Poincare
    \brief   provides Lorentz transformation (Boosts) and rotation for Vec4D
    
    This class can be used to transform 4-vectors (Vec4D) from one frame 
    to another (Lorentz transformation) or to rotate the 3-vector part.
  */

  /*!
    \fn Poincare::Poincare()
    \brief Standard constructor (Unity transformation)
  */

  /*!
    \fn   Poincare::Poincare(Vec4D v1, Vec4D v2)
    \brief  Special constructor initializes a rotation

    This constructor initializes a rotation matrix 
    from the space parts of \em v1 and \em v2 (length of all vectors is preserved).
  */

  /*!
    \fn Poincare::Poincare(Vec4D v)
    \brief Special constuctor initiaizes a boost

    This constructor initializes a boost matrix that would
    transform \em v in its center of mass frame.
  */
 
  /*!
    \fn void Poincare::Boost(Vec4D& v)
    \brief boosts a vector in a CMS (of vector given at initialization)
  */

  /*!
    \fn void Poincare::BoostBack(Vec4D& v)
    \brief applies the invers boost transformation (i.e. back to LAB Frame)
  */

  /*!
    \fn void Poincare::Rotate(Vec4D& v)
    \brief rotates a vector (as given by a vector pair during initalization)
  */

  /*!
    \fn void Poincare::RotateBack(Vec4D& v)
    \brief applies the invers rotation
  */

  /*!
    \fn Vec4D Poincare::operator*(const Vec4D& vin)
    \brief returns a boosted or rotated vector (depends on initialization)

    With this operation a boost or rotation transformation can be
    performed by multiplying a Poincare object with a 4-vector, meaning
    \f[
       p'_\mu = \Lambda_\mu^\nu p_\nu
    \f]

    A transformation of two vectors in their common CMS frame looks as follows:
    \verbatim

      Vec4D  a(45.,0.,0.,45.);
      Vec4D  b(21.,0.,0.,-21.);
      Poincare lambda(a+b);
      a = lambda*a;
      b = lambda*b;
    \endverbatim
  */

}

#endif
