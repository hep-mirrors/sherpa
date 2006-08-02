#include "ComplexVector.H"

using namespace ATOOLS;

ComplexVec4D::ComplexVec4D(const Complex x0, const Complex x1, const Complex x2, const Complex x3) {
  m_x[0] = x0; 
  m_x[1] = x1; 
  m_x[2] = x2; 
  m_x[3] = x3; 
}


ComplexVec4D::ComplexVec4D(const ATOOLS::Vec4D real, const ATOOLS::Vec4D imag) {
  m_x[0] = Complex(real[0], imag[0]) ;
  m_x[1] = Complex(real[1], imag[1]) ;
  m_x[2] = Complex(real[2], imag[2]) ;
  m_x[3] = Complex(real[3], imag[3]) ;
}


ComplexVec4D::ComplexVec4D(const ATOOLS::Vec4D real) {
  m_x[0] = Complex(real[0], 0.0) ;
  m_x[1] = Complex(real[1], 0.0) ;
  m_x[2] = Complex(real[2], 0.0) ;
  m_x[3] = Complex(real[3], 0.0) ;
}


Complex& ComplexVec4D::operator[] (int i) 
{
#ifdef CHECK
  if(i<0 || i>3) {
    cerr<<"ComplexVec4D: out of bound.\n";
    return m_x[0];
  }
#endif
  return m_x[i];
}


const Complex& ComplexVec4D::operator[] (int i) const
{
#ifdef CHECK
  if(i<0 || i>3) {
    cerr<<"ComplexVec4D: out of bound.\n";
    return m_x[0];
  }
#endif
    return m_x[i];
}


ComplexVec4D ComplexVec4D::Conjugate()
{
  return ComplexVec4D( conj(m_x[0]), conj(m_x[1]), conj(m_x[2]), conj(m_x[3]) );
}


ComplexVec4D& ComplexVec4D::operator+= (const ComplexVec4D& v)
{
  m_x[0] += v[0];
  m_x[1] += v[1];
  m_x[2] += v[2];
  m_x[3] += v[3];
  return *this;
}

ComplexVec4D ATOOLS::operator+ (const ComplexVec4D& c1, const ComplexVec4D& c2)
{
  return ComplexVec4D(c1[0]+c2[0], c1[1]+c2[1], c1[2]+c2[2], c1[3]+c2[3]);
}

ComplexVec4D ATOOLS::operator+ (const ComplexVec4D& c1, const Vec4D& c2)
{
  return ComplexVec4D(c1[0]+c2[0], c1[1]+c2[1], c1[2]+c2[2], c1[3]+c2[3]);
}

ComplexVec4D ATOOLS::operator+ (const Vec4D& c1, const ComplexVec4D& c2)
{
  return c2+c1;
}


ComplexVec4D& ComplexVec4D::operator-= (const ComplexVec4D& v)
{
  m_x[0] -= v[0];
  m_x[1] -= v[1];
  m_x[2] -= v[2];
  m_x[3] -= v[3];
  return *this;
}

ComplexVec4D ATOOLS::operator- (const ComplexVec4D& c1, const ComplexVec4D& c2)
{
  return ComplexVec4D(c1[0]-c2[0], c1[1]-c2[1], c1[2]-c2[2], c1[3]-c2[3]);
}

ComplexVec4D ATOOLS::operator- (const ComplexVec4D& c1, const Vec4D& c2)
{
  return ComplexVec4D(c1[0]-c2[0], c1[1]-c2[1], c1[2]-c2[2], c1[3]-c2[3]);
}

ComplexVec4D ATOOLS::operator- (const Vec4D& c1, const ComplexVec4D& c2)
{
  return ComplexVec4D(c1[0]-c2[0], c1[1]-c2[1], c1[2]-c2[2], c1[3]-c2[3]);
}


ComplexVec4D& ComplexVec4D::operator*= (const Complex& c)
{
  m_x[0] *= c;
  m_x[1] *= c;
  m_x[2] *= c;
  m_x[3] *= c;
  return *this;
}

ComplexVec4D& ComplexVec4D::operator*= (const double& c)
{
  m_x[0] *= c;
  m_x[1] *= c;
  m_x[2] *= c;
  m_x[3] *= c;
  return *this;
}

Complex ATOOLS::operator* (const ComplexVec4D& c1, const ComplexVec4D& c2) 
{
  return c1[0]*c2[0]-c1[1]*c2[1]-c1[2]*c2[2]-c1[3]*c2[3];
}

ComplexVec4D ATOOLS::operator* (const Complex scal, const ComplexVec4D& c1)
{
  return ComplexVec4D(scal*c1[0], scal*c1[1], scal*c1[2], scal*c1[3]);
}

ComplexVec4D ATOOLS::operator* (const ComplexVec4D& c1, const Complex scal)
{
  return scal*c1;
}

ComplexVec4D ATOOLS::operator* (const double scal, const ComplexVec4D& c1)
{
  return ComplexVec4D(scal*c1[0], scal*c1[1], scal*c1[2], scal*c1[3]);
}

ComplexVec4D ATOOLS::operator* (const ComplexVec4D& c1, const double scal)
{
  return scal*c1;
}

Complex ATOOLS::operator* (const ComplexVec4D& c1, const Vec4D& c2)
{
  return c1[0]*c2[0]-c1[1]*c2[1]-c1[2]*c2[2]-c1[3]*c2[3];
}

Complex ATOOLS::operator* (const Vec4D& c1, const ComplexVec4D& c2)
{
  return c2*c1;
}


ComplexVec4D& ComplexVec4D::operator/= (const double& d)
{
  m_x[0] /= d;
  m_x[1] /= d;
  m_x[2] /= d;
  m_x[3] /= d;
  return *this;
}


ComplexVec4D ATOOLS::cross(const ComplexVec4D& v1, const Vec4D& v2, const Vec4D& v3)
{
  Complex x0 = -v1[1]*v2[2]*v3[3] -v1[2]*v2[3]*v3[1] -v1[3]*v2[1]*v3[2] +v1[1]*v2[3]*v3[2] +v1[3]*v2[2]*v3[1] +v1[2]*v2[1]*v3[3];
  Complex x1 = -v1[0]*v2[2]*v3[3] -v1[2]*v2[3]*v3[0] -v1[3]*v2[0]*v3[2] +v1[0]*v2[3]*v3[2] +v1[3]*v2[2]*v3[0] +v1[2]*v2[0]*v3[3];
  Complex x2 =  v1[0]*v2[1]*v3[3] +v1[1]*v2[3]*v3[0] +v1[3]*v2[0]*v3[1] -v1[0]*v2[3]*v3[1] -v1[3]*v2[1]*v3[0] -v1[1]*v2[0]*v3[3];
  Complex x3 = -v1[0]*v2[1]*v3[2] -v1[1]*v2[2]*v3[0] -v1[2]*v2[0]*v3[1] +v1[0]*v2[2]*v3[1] +v1[2]*v2[1]*v3[0] +v1[1]*v2[0]*v3[2];
  return ComplexVec4D(x0,x1,x2,x3);
}


std::ostream& operator<< (std::ostream& s, const ATOOLS::ComplexVec4D& c) {
  s<<"[\t"<<c[0];
  for(int i=1; i<4; i++) {
    s<<",\t"<<c[i];
  }
  s<<"\t]";
  return s;
}
