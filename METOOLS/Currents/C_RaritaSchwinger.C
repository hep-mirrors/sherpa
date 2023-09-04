#include "METOOLS/Currents/C_RaritaSchwinger.H"

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Exception.H"
#include "METOOLS/Main/SpinFuncs.H"
#include "METOOLS/Currents/C_Vector.H"

using namespace METOOLS;

template <class Scalar>
double CRaritaSchwinger<Scalar>::s_accu(1.0e-12);

// allows the output of the Rarita-Schwinger vector-spinor with std::cout
template <class Scalar> std::ostream &
METOOLS::operator<<(std::ostream &s,const CRaritaSchwinger<Scalar> &rs)
{
  //encoding m_h
  std::string helicity;
  if (rs.H()==2) helicity = "++";
  else if (rs.H()==1) helicity = "+";
  else if (rs.H()==-1) helicity = "-";
  else if (rs.H()==-2) helicity = "--";
  else THROW(fatal_error, "The value of the helicity of the Rarita-Schwinger particle is not permitted.")

  return s<<'{'<<(rs.B()>0?(rs.R()>0?"Ubar$"+helicity:"$Vbar$"+helicity):
  (rs.R()>0?"U"+helicity:"V"+helicity)) <<","<<rs.S()<<";"<<rs(0)<<","<<rs(1)<<'|'
	  <<rs[0]<<','<<rs[1]<<','<<rs[2]<<','<<rs[3]<<','<<rs[4]<<','<<rs[5]<<','<<rs[6]<<','<<rs[7]<<','<<rs[8]<<','
    <<rs[9]<<','<<rs[10]<<','<<rs[11]<<','<<rs[12]<<','<<rs[13]<<','<<rs[14]<<','<<rs[15]<<'}';
}

template <class Scalar>
bool CRaritaSchwinger<Scalar>::Nan() const
{
  for(auto & i : m_x) {
    if (ATOOLS::IsNan(i)) return true;
  }
  return false;
}

template <class Scalar>
void CRaritaSchwinger<Scalar>::ResetAccu()
{ 
  s_accu=Accu(); 
}

/*template <class Scalar>
void CRaritaSchwinger<Scalar>::Add(const CObject *c)
{
  const CVec4 *v(static_cast<const CVec4*>(c));
  m_x[0]+=v->m_x[0]; 
  m_x[1]+=v->m_x[1];
  m_x[2]+=v->m_x[2]; 
  m_x[3]+=v->m_x[3];
}

template <class Scalar>
void CVec4<Scalar>::Divide(const double &d)
{
  m_x[0]/=Scalar(d);
  m_x[1]/=Scalar(d); 
  m_x[2]/=Scalar(d); 
  m_x[3]/=Scalar(d);
}

template <class Scalar>
void CVec4<Scalar>::Multiply(const Complex &c)
{
  m_x[0]*=SComplex(c);
  m_x[1]*=SComplex(c); 
  m_x[2]*=SComplex(c); 
  m_x[3]*=SComplex(c);
}

template <class Scalar>
void CVec4<Scalar>::Invert()
{
  m_x[0]=-m_x[0]; 
  m_x[1]=-m_x[1]; 
  m_x[2]=-m_x[2]; 
  m_x[3]=-m_x[3]; 
}*/

template <class Scalar>
bool CRaritaSchwinger<Scalar>::IsZero() const
{
  for (size_t i(0); i<16; ++i){
    if (m_x[i]!=Scalar(0.0)) return false;
  }
  return true;
}

// TODO: Was machen diese Funktionen? Bislang analog zu Spinor und Vector -  stimmt das?
template <class Scalar>
typename ATOOLS::AutoDelete_Vector<CRaritaSchwinger<Scalar> >
CRaritaSchwinger<Scalar>::s_objects;

template <class Scalar>
CRaritaSchwinger<Scalar> *CRaritaSchwinger<Scalar>::New()
{
#ifndef USING__Threading
  if (s_objects.empty())
#endif
    return new CRaritaSchwinger();
  CRaritaSchwinger *v(s_objects.back());
  s_objects.pop_back();
  return v;
}

template <class Scalar>
CRaritaSchwinger<Scalar> *CRaritaSchwinger<Scalar>::New(const CRaritaSchwinger &s)
{
#ifndef USING__Threading
  if (s_objects.empty())
#endif
    return new CRaritaSchwinger(s);
  CRaritaSchwinger *v(s_objects.back());
  s_objects.pop_back();
  *v=s;
  return v;
}

// TODO: dieses new anpassen, wenn wir verstehen, was die einzelnen Parameter bedeuten...
/*template <class Scalar>
CVec4<Scalar> *CVec4<Scalar>::New
(const Scalar &x0, const Scalar &x1, 
 const Scalar &x2, const Scalar &x3,
 const int c1,const int c2,
 const size_t &h,const size_t &s)
{
#ifndef USING__Threading
  if (s_objects.empty())
#endif
    return new CVec4(x0,x1,x2,x3,c1,c2,h,s);
  CVec4 *v(s_objects.back());
  s_objects.pop_back();
  v->m_x[0]=x0;
  v->m_x[1]=x1;
  v->m_x[2]=x2;
  v->m_x[3]=x3; 
  v->m_c[0]=c1;
  v->m_c[1]=c2;
  v->m_h=h;
  v->m_s=s;
  return v;
}*/

template <class Scalar>
CObject *CRaritaSchwinger<Scalar>::Copy() const
{
  return CRaritaSchwinger::New(*this);
}

template <class Scalar>
void CRaritaSchwinger<Scalar>::Delete()
{
#ifndef USING__Threading
  s_objects.push_back(this);
#else
  delete this;
#endif
}

template <class Scalar> bool CRaritaSchwinger<Scalar>::SetOn()
{
  m_on=0;
  if (m_x[0]!=ZERO || m_x[1]!=ZERO || m_x[2]!=ZERO || m_x[3]!=ZERO || m_x[4]!=ZERO || m_x[5]!=ZERO || m_x[6]!=ZERO ||
      m_x[7]!=ZERO) m_on|=1;
  if (m_x[8]!=ZERO || m_x[9]!=ZERO || m_x[10]!=ZERO || m_x[11]!=ZERO || m_x[12]!=ZERO || m_x[13]!=ZERO ||
      m_x[14]!=ZERO || m_x[15]!=ZERO) m_on|=2;
  return m_on&3;
}

// TODO: Springt er hier wirklich zum richtigen Operator?
template<class Scalar>
bool CRaritaSchwinger<Scalar>::Test_Properties(const ATOOLS::Vec4D &p, int r) {
  // Dirac equation (gamma_mu * p_mu -m)^A_B RS^B, nu -> auch +m? für V statt U?
  METOOLS::Gamma gammavec = Gamma();
  ATOOLS::CMatrix intermediate = (gammavec * p + Complex(r) * p.Abs2() * ATOOLS::CMatrix(4, 1));
  std::vector<Scalar> result1(16, ZERO);
  for(int i(0); i<4; ++i) {
    result1[0+2*i] = intermediate[0][0] * (*this)[2*i] + intermediate[0][1] * (*this)[1+2*i] + intermediate[0][2] * (*this)[8+2*i] + intermediate[0][3] * (*this)[9+2*i];
    result1[1+2*i] = intermediate[1][0] * (*this)[2*i] + intermediate[1][1] * (*this)[1+2*i] + intermediate[1][2] * (*this)[8+2*i] + intermediate[1][3] * (*this)[9+2*i];
    result1[8+2*i] = intermediate[2][0] * (*this)[2*i] + intermediate[2][1] * (*this)[1+2*i] + intermediate[2][2] * (*this)[8+2*i] + intermediate[2][3] * (*this)[9+2*i];
    result1[9+2*i] = intermediate[3][0] * (*this)[2*i] + intermediate[3][1] * (*this)[1+2*i] + intermediate[3][2] * (*this)[8+2*i] + intermediate[3][3] * (*this)[9+2*i];
  }
  /* Implemented for filling when the four Dirac spinors are above each other in the spin-3/2 wave function
 * for (int i(0); i<4; ++i){
      for (size_t j(0); j<4; ++j){
        if (rs.On()!=1){
          comps[i] += gamma[i][j] * rs[j];
          comps[i+4] += gamma[i][j] * rs[j+4];
        }
        if (rs.On()!=2){
          comps[i+8] += gamma[i][j] * rs[j+8];
          comps[i+12] += gamma[i][j] * rs[j+12];
        }
      }
    }*/
  for (size_t j(0); j<16; ++j){
    if (result1[j]>s_accu) return false;
  }

  // gamma_mu^A_B times RS^B,mu = 0
  // ÜBERPRÜFEN!!!
  //METOOLS::CVec4<Scalar> result2 = CVec4<Scalar>();
  std::vector<Scalar> result2(4, ZERO);
  for (int i(0); i<4; ++i){
    result2[0] += gammavec[i][0][0] * (*this)[2*i] + gammavec[i][0][1] * (*this)[1+2*i] + gammavec[i][0][2] * (*this)[8+2*i] + gammavec[i][0][3] * (*this)[9+2*i];
    result2[1] += gammavec[i][1][0] * (*this)[2*i] + gammavec[i][1][1] * (*this)[1+2*i] + gammavec[i][1][2] * (*this)[8+2*i] + gammavec[i][1][3] * (*this)[9+2*i];
    result2[2] += gammavec[i][2][0] * (*this)[2*i] + gammavec[i][2][1] * (*this)[1+2*i] + gammavec[i][2][2] * (*this)[8+2*i] + gammavec[i][2][3] * (*this)[9+2*i];
    result2[3] += gammavec[i][3][0] * (*this)[2*i] + gammavec[i][3][1] * (*this)[1+2*i] + gammavec[i][3][2] * (*this)[8+2*i] + gammavec[i][3][3] * (*this)[9+2*i];
    // result2 += gammavec[i] * METOOLS::CVec4<Scalar>((*this)[i+i*4], (*this)[i+i*4+1], (*this)[i+i*4+2], (*this)[i+i*4+3]);
  }
  if (!(result2[0]<s_accu && result2[1]<s_accu && result2[2]<s_accu && result2[3]<s_accu)) return false;
  
  // p_mu times RS = 0 ? aus partielle Ableitung_mu RS=0?
  std::vector<Scalar> result3(4, ZERO);
  for (int i(0); i<2; ++i){
    result3[i] = p[0] * (*this)[i] + p[1] * (*this)[2+i] + p[2] * (*this)[4+i] + p[3] * (*this)[6+i];
    result3[i+2] = p[0] * (*this)[i+8] + p[1] * (*this)[10+i] + p[2] * (*this)[12+i] + p[3] * (*this)[14+i];
  }
  if (!(result3[0]<s_accu && result3[1]<s_accu && result3[2]<s_accu && result3[3]<s_accu)) return false;

  // normalizations???

  // completeness -> in Stromklasse testen!!!
  return true;
}

// TODO: kompliert nicht :-(
namespace METOOLS {

  //template class DCRaritaSchwinger;
  //template std::ostream &operator<<(std::ostream &ostr,const DCRaritaSchwinger &s);

  //template class QCRaritaSchwinger;
  //template std::ostream &operator<<(std::ostream &ostr,const QCRaritaSchwinger &s);

}
