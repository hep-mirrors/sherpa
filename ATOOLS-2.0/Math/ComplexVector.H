#ifndef ComplexVector_H
#define ComplexVector_H

#include <iostream>
#include "MathTools.H"
#include "Vector.H"
 
namespace ATOOLS {

  class ComplexVec4D {
    Complex m_x[4];
  public:
    ComplexVec4D() {
      m_x[0]=m_x[1]=m_x[2]=m_x[3]= Complex(0., 0.);
    };
    ComplexVec4D(const Complex x0, const Complex x1, const Complex x2, const Complex x3);
    ComplexVec4D(const ATOOLS::Vec4D real, const ATOOLS::Vec4D imag);
    ComplexVec4D(const ATOOLS::Vec4D real);
    Complex& operator[](int i);
    const Complex& operator[] (int i) const;
    ComplexVec4D& operator+= (const ComplexVec4D& v);
    ComplexVec4D& operator-= (const ComplexVec4D& v);
    ComplexVec4D& operator*= (const Complex& c);
    ComplexVec4D& operator*= (const double& c);
    ComplexVec4D& operator/= (const double& d);
    
    ComplexVec4D Conjugate();
  };
  
  Complex operator* (const ComplexVec4D& c1, const ComplexVec4D& c2);
  ComplexVec4D operator+ (const ComplexVec4D& c1, const ComplexVec4D& c2);
  ComplexVec4D operator- (const ComplexVec4D& c1, const ComplexVec4D& c2);
  
  ComplexVec4D operator* (const Complex scal, const ComplexVec4D& c1);
  ComplexVec4D operator* (const ComplexVec4D& c1, const Complex scal);
  
  ComplexVec4D operator* (const double scal, const ComplexVec4D& c1);
  ComplexVec4D operator* (const ComplexVec4D& c1, const double scal);

  Complex operator* (const ComplexVec4D& c1, const Vec4D& c2);
  Complex operator* (const Vec4D& c1, const ComplexVec4D& c2);
  ComplexVec4D operator+ (const ComplexVec4D& c1, const Vec4D& c2);
  ComplexVec4D operator+ (const Vec4D& c1, const ComplexVec4D& c2);
  ComplexVec4D operator- (const ComplexVec4D& c1, const Vec4D& c2);
  ComplexVec4D operator- (const Vec4D& c1, const ComplexVec4D& c2);

  ComplexVec4D cross(const ComplexVec4D& a, const  Vec4D& b, const  Vec4D& c);
}

std::ostream& operator<<(std::ostream& s, const ATOOLS::ComplexVec4D& c);

#endif
