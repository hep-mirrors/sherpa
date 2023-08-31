#include "METOOLS/Currents/C_RaritaSchwinger.H"

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Exception.H"
#include "METOOLS/Main/SpinFuncs.H"

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

template<class Scalar>
bool CRaritaSchwinger<Scalar>::Test_Properties(const ATOOLS::Vec4D &p) {
  bool passed = true;
  // Dirac equation (gamma_mu * p_mu -m)^A_B RS^B, mu -> auch +m? für V statt U?
  METOOLS::Gamma gammavec = Gamma();
  ATOOLS::Vec4D<METOOLS::CMatrix> intermediate = gammavec * p - p.Abs2();
  // gamma_mu^A_B times RS^B,mu = 0



  // p_mu times RS = 0 ? aus partielle Ableitung_mu RS=0?

  // normalizations???

  // completeness -> in Stromklasse testen!!!
  return passed;
}

// TODO: kompliert nicht :-(
namespace METOOLS {

  //template class DCRaritaSchwinger;
  //template std::ostream &operator<<(std::ostream &ostr,const DCRaritaSchwinger &s);

  //template class QCRaritaSchwinger;
  //template std::ostream &operator<<(std::ostream &ostr,const QCRaritaSchwinger &s);

}
