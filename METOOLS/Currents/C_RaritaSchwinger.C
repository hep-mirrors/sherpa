#include "METOOLS/Currents/C_RaritaSchwinger.H"

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Exception.H"
#include "METOOLS/Main/SpinFuncs.H"
#include "METOOLS/Currents/C_Vector.H"

using namespace METOOLS;

template <class Scalar>
double CRaritaSchwinger<Scalar>::s_accu(1.0e-12);

// allows the output of the Rarita-Schwinger vector-spinor with std::cout
// TODO: Wie sind die Hs nun eigentlich zu setzen? Nur positiv oder ganzzahlige Spinorwerte?
template <class Scalar> std::ostream &
METOOLS::operator<<(std::ostream &s,const CRaritaSchwinger<Scalar> &rs)
{
   return s<<'{'<<(rs.B()<0?(rs.R()>0?"Ubar(":"Vbar("):
  (rs.R()>0?"U(":"V(")) << rs.On() << "),"<< rs.H() << "," << rs.S()<<";"<<rs(0)<<","<<rs(1)<<'|'
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

//TODO: Bessere Lösung als THROW?
template<class Scalar>
void CRaritaSchwinger<Scalar>::Add(const CObject *c) {
  const CRaritaSchwinger *rs(static_cast<const CRaritaSchwinger*>(c));
  if (rs->m_b!=m_b || rs->m_r!=m_r) THROW(fatal_error, "Summing particle and anti-particle wave functions / bar "
                                                       "and non-bar wave functions is not allowed!")
  m_on|=rs->m_on;
  if (m_on&1) {
    m_x[0]+=rs->m_x[0]; m_x[1]+=rs->m_x[1]; m_x[4]+=rs->m_x[4]; m_x[5]+=rs->m_x[5]; m_x[8]+=rs->m_x[8];
    m_x[9]+=rs->m_x[9]; m_x[12]+=rs->m_x[12]; m_x[13]+=rs->m_x[13];
  }
  if (m_on&2) {
    m_x[2]+=rs->m_x[2]; m_x[3]+=rs->m_x[3]; m_x[6]+=rs->m_x[6]; m_x[7]+=rs->m_x[7]; m_x[10]+=rs->m_x[10];
    m_x[11]+=rs->m_x[11]; m_x[14]+=rs->m_x[14]; m_x[15]+=rs->m_x[15];
  }
}
template<class Scalar>
void CRaritaSchwinger<Scalar>::Multiply(const Complex &c) {
  if (m_on&1) {
    m_x[0]*=SComplex(c) ; m_x[1]*=SComplex(c); m_x[4]*=SComplex(c); m_x[5]*=SComplex(c); m_x[8]*=SComplex(c);
    m_x[9]*=SComplex(c); m_x[12]*=SComplex(c); m_x[13]*=SComplex(c);
  }
  if (m_on&2) {
    m_x[2]*=SComplex(c); m_x[3]*=SComplex(c); m_x[6]*=SComplex(c); m_x[7]*=SComplex(c); m_x[10]*=SComplex(c);
    m_x[11]*=SComplex(c); m_x[14]*=SComplex(c); m_x[15]*=SComplex(c);
  }
}

template<class Scalar>
void CRaritaSchwinger<Scalar>::Divide(const double &d) {
  if (d==0.0) THROW(fatal_error, "Zero Division");
  if (m_on&1) {
    m_x[0]/=SComplex(d) ; m_x[1]/=SComplex(d); m_x[4]/=SComplex(d); m_x[5]/=SComplex(d); m_x[8]/=SComplex(d);
    m_x[9]/=SComplex(d); m_x[12]/=SComplex(d); m_x[13]/=SComplex(d);
  }
  if (m_on&2) {
    m_x[2]/=SComplex(d); m_x[3]/=SComplex(d); m_x[6]/=SComplex(d); m_x[7]/=SComplex(d); m_x[10]/=SComplex(d);
    m_x[11]/=SComplex(d); m_x[14]/=SComplex(d); m_x[15]/=SComplex(d);
  }
}
template<class Scalar>
void CRaritaSchwinger<Scalar>::Invert() {
  for (size_t i(0); i<16; ++i) m_x[i] = -m_x[i];
}

template<class Scalar>
ATOOLS::TCMatrix<Scalar> CRaritaSchwinger<Scalar>::Contract4Index(const CRaritaSchwinger<Scalar> &rs) const {
  SComplex** intermediate = new SComplex*[4];
  for (int i=0;i<4;++i) intermediate[i] = new SComplex[4];
  for (size_t i(0); i<4; ++i){
    for (size_t j(0); j<4; ++j){
      intermediate[i][j] = (*this)[i] * rs[j] - (*this)[i+4] * rs[j+4] - (*this)[i+8] * rs[j+8] - (*this)[i+12] * rs[j+12];
    }
  }
  return ATOOLS::TCMatrix<Scalar>(intermediate, 4);
}

template<class Scalar>
ATOOLS::TCMatrix<Scalar> CRaritaSchwinger<Scalar>::ContractSpinorIndex(const CRaritaSchwinger<Scalar> &rs) const {
  SComplex** intermediate = new SComplex*[4];
  for (int i=0;i<4;++i) intermediate[i] = new SComplex[4];
  for (size_t i(0); i<4; ++i){
    for (size_t j(0); j<4; ++j){
      intermediate[i][j] = (*this)[4*i] * rs[4*j] + (*this)[4*i+1] * rs[4*j+1] + (*this)[4*i+2] * rs[4*j+2] +
        (*this)[4*i+3] * rs[4*j+3];
    }
  }
  return ATOOLS::TCMatrix<Scalar>(intermediate, 4);
}
  // TODO: Abprüfung von Eigenschaften wie bar oder Teilchen/Antiteilchen in Operatoren!
  // TODO: Bessere Lösung als THROW?
template <class Scalar> std::complex<Scalar>
CRaritaSchwinger<Scalar>::operator*(const CRaritaSchwinger<Scalar> &rs) const
{
  //if (rs.m_b==m_b) return (*this)*rs.CConj();
  if (rs.m_r!=m_r) THROW(fatal_error, "Multiplying particle and anti-particle wave functions is not allowed!")
  std::complex<Scalar> result(0.0, 0.0);
  std::complex<Scalar> sign(1);
  for (size_t i(0); i<16; ++i){
    if (i>=4) sign=-1;
    result+=sign*m_x[i]*rs.m_x[i];
  }
  return result;
}

template <class Scalar>
bool CRaritaSchwinger<Scalar>::IsZero() const
{
  for (size_t i(0); i<16; ++i){
    if (m_x[i]!=Scalar(0.0)) return false;
  }
  return true;
}

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
  if (m_x[0]!=ZERO || m_x[1]!=ZERO || m_x[4]!=ZERO || m_x[5]!=ZERO || m_x[8]!=ZERO || m_x[9]!=ZERO || m_x[12]!=ZERO ||
      m_x[13]!=ZERO) m_on|=1;
  if (m_x[2]!=ZERO || m_x[3]!=ZERO || m_x[6]!=ZERO || m_x[7]!=ZERO || m_x[10]!=ZERO || m_x[11]!=ZERO ||
      m_x[14]!=ZERO || m_x[15]!=ZERO) m_on|=2;
  return m_on&3;
}

template<class Scalar>
bool CRaritaSchwinger<Scalar>::Test_Properties(const ATOOLS::Vec4D &p, int r, int b) {
  // Dirac equation (gamma_mu * p_mu -m)^A_B RS^B, nu -> auch +m? für V statt U?
  std::cout<<METHOD<<": Testing Dirac equation..."<<std::endl;
  METOOLS::Gamma gammavec = Gamma<Scalar>();
  bool testresult(true);
  bool single_testresult(true);
  // SComplex(-r) for distinguishing between particle and anti-particle Dirac equation
  ATOOLS::TCMatrix<Scalar> intermediate = ATOOLS::TCMatrix<Scalar>(gammavec * p);
  if (sqrt(p.Abs2())>1e-6)
    intermediate += ATOOLS::TCMatrix<Scalar>(SComplex(-r) * SComplex(sqrt(p.Abs2())) * ATOOLS::TCMatrix<Scalar>(4, true));
  std::vector<SComplex> result1(16);
  // TODO: Sollte man die exakten ZEROS auch hier explizit implementieren?
  for (int i(0); i<4; ++i){
    for (size_t j(0); j<4; ++j){
        if (b>0){
          result1[i] += intermediate[i][j] * (*this)[j];
          result1[i+4] += intermediate[i][j] * (*this)[j+4];
          result1[i+8] += intermediate[i][j] * (*this)[j+8];
          result1[i+12] += intermediate[i][j] * (*this)[j+12];
        }
        else{
          result1[j] += (*this)[i] * intermediate[i][j];
          result1[j+4] += (*this)[i+4] * intermediate[i][j];
          result1[j+8] += (*this)[i+8] * intermediate[i][j];
          result1[j+12] += (*this)[i+12] * intermediate[i][j];
        }
    }
  }

  for (size_t j(0); j<16; ++j){
    if (std::abs(result1[j].real())>s_accu || std::abs(result1[j].imag())>s_accu) {
      msg_Out()<<"Component " << j << " of resulting Rarita-Schwinger wave function is " << result1[j] << " instead of zero!"
      << std::endl;
      testresult =false;
      single_testresult=false;
    }
  }
  if (single_testresult) msg_Out()<< "passed" << std::endl;

  // gamma_mu^A_B times RS^B,mu = 0
  std::cout<<METHOD<<": Testing gamma_mu Psi^mu = 0..."<<std::endl;
  single_testresult = true;
  std::vector<SComplex> result2(4);
  for (size_t i(0); i<4; ++i){
    for (size_t j(0); j<4; ++j){
      if (b>0) result2[i] += gammavec[0][i][j] * (*this)[j] - gammavec[1][i][j] * (*this)[j+4] -
                             gammavec[2][i][j] * (*this)[j+8] - gammavec[3][i][j] * (*this)[j+12];
      else{
        for (size_t k(0); k<4; ++k){
          result2[k] += (*this)[i] * gammavec[0][i][j] * gammavec[0][j][k] -
                        (*this)[i+4] * gammavec[1][i][j] * gammavec[0][j][k] -
                        (*this)[i+8] * gammavec[2][i][j] * gammavec[0][j][k] -
                        (*this)[i+12] * gammavec[3][i][j] * gammavec[0][j][k];
        }
      }
    }
  }

  if (!(std::abs(result2[0].real())<s_accu && std::abs(result2[1].real())<s_accu && std::abs(result2[2].real())<s_accu && std::abs(result2[3].real())<s_accu)){
    msg_Out() << "gamma_mu Psi^mu is " << result2[0] << result2[1] << result2[2] << result2[3] << " not zero!" << std::endl;
    testresult = false; single_testresult = false;
  }
  if (!(std::abs(result2[0].imag())<s_accu && std::abs(result2[1].imag())<s_accu && std::abs(result2[2].imag())<s_accu && std::abs(result2[3].imag())<s_accu)){
    msg_Out() << "gamma_mu Psi^mu is " << result2[0] << result2[1] << result2[2] << result2[3] << " not zero!" << std::endl;
    testresult = false; single_testresult = false;
  }
  if (single_testresult) msg_Out()<< "passed" << std::endl;

  // p_mu times RS = 0
  std::cout<<METHOD<<": Testing p_mu Psi^mu = 0..."<<std::endl;
  single_testresult = true;
  std::vector<SComplex> result3(4);
  for (int i(0); i<4; ++i){
    result3[i] = SComplex(p[0]) * (*this)[i] + SComplex(-p[1]) * (*this)[4+i] + SComplex(-p[2]) * (*this)[8+i] + SComplex(-p[3]) * (*this)[12+i];
  }
  if (!(std::abs(result3[0].real())<s_accu && std::abs(result3[1].real())<s_accu && std::abs(result3[2].real())<s_accu && std::abs(result3[3].real())<s_accu)){
    msg_Out() << "p_mu Psi^mu is " << result3[0] << result3[1] << result3[2] << result3[3] << " not zero!" << std::endl;
    testresult = false; single_testresult = false;
  }
  if (!(std::abs(result3[0].imag())<s_accu && std::abs(result3[1].imag())<s_accu && std::abs(result3[2].imag())<s_accu && std::abs(result3[3].imag())<s_accu)) {
    msg_Out() << "p_mu Psi^mu is " << result3[0] << result3[1] << result3[2] << result3[3] << " not zero!" << std::endl;
    testresult = false; single_testresult = false;
  }
  if (single_testresult) msg_Out()<< "passed" << std::endl;
  return testresult;
}

namespace METOOLS {

  template class DCRaritaSchwinger;
  template std::ostream &operator<<(std::ostream &ostr,const DCRaritaSchwinger &s);

  template class QCRaritaSchwinger;
  template std::ostream &operator<<(std::ostream &ostr,const QCRaritaSchwinger &s);

}
