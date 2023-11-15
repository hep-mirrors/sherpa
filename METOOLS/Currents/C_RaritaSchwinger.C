#include "METOOLS/Currents/C_RaritaSchwinger.H"

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Exception.H"
#include "METOOLS/Main/SpinFuncs.H"
#include "METOOLS/Currents/C_Vector.H"

using namespace METOOLS;

template <class Scalar>
double CRaritaSchwinger<Scalar>::s_accu(1.0e-12);

// allows the output of the Rarita-Schwinger vector-spinor with std::cout
// TODO: Outputfunktion testen
// TODO: Wie sind die Hs nun eigentlich zu setzen? Nur positiv oder ganzzahlige Spinorwerte?
template <class Scalar> std::ostream &
METOOLS::operator<<(std::ostream &s,const CRaritaSchwinger<Scalar> &rs)
{
  //encoding m_h
  std::string helicity;
  if (rs.H()==3) helicity = "++";
  else if (rs.H()==1) helicity = "+";
  else if (rs.H()==-1) helicity = "-";
  else if (rs.H()==-3) helicity = "--";
  // TODO: Wieder reinnehmen, wenn geklärt ist, wie H() nun aussehen soll!!!
  //else THROW(fatal_error, "The value of the helicity of the Rarita-Schwinger particle is not permitted.")

  return s<<'{'<<(rs.B()<0?(rs.R()>0?"Ubar"+helicity:"Vbar"+helicity):
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

//TODO: Bessere Lösung als THROW? Was genau bedeutet das auskommentierte? -> Abgeschrieben von Spinor
// TODO: FUnktioniert das auch für m_on =3? Sind dann beide ifs true?
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

// TODO: WANN LÖSCHE ICH DAS ZEUG HIER WIEDER?
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

  //TODO: Bessere Lösung als THROW? Was genau bedeutet das auskommentierte? -> Abgeschrieben von Spinor
template <class Scalar> std::complex<Scalar>
CRaritaSchwinger<Scalar>::operator*(const CRaritaSchwinger<Scalar> &rs) const
{
  //if (rs.m_b==m_b) return (*this)*rs.CConj();
  if (rs.m_b!=m_b || rs.m_r!=m_r) THROW(fatal_error, "Multiplying particle and anti-particle wave functions / bar "
                                                       "and non-bar wave functions is not allowed!")
  std::complex<Scalar> result(0.0, 0.0);
  std::complex<Scalar> sign(1);
  for (size_t i(0); i<16; ++i){
    if (i>=4) sign=-1;
    result+=sign*m_x[i]*rs.m_x[i];
  }
  return result;
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

// Für Novaes Füllung
/*template <class Scalar> bool CRaritaSchwinger<Scalar>::SetOn()
{
  m_on=0;
  if (m_x[0]!=ZERO || m_x[1]!=ZERO || m_x[2]!=ZERO || m_x[3]!=ZERO || m_x[4]!=ZERO || m_x[5]!=ZERO || m_x[6]!=ZERO ||
      m_x[7]!=ZERO) m_on|=1;
  if (m_x[8]!=ZERO || m_x[9]!=ZERO || m_x[10]!=ZERO || m_x[11]!=ZERO || m_x[12]!=ZERO || m_x[13]!=ZERO ||
      m_x[14]!=ZERO || m_x[15]!=ZERO) m_on|=2;
  return m_on&3;
}*/

template <class Scalar> bool CRaritaSchwinger<Scalar>::SetOn()
{
  m_on=0;
  if (m_x[0]!=ZERO || m_x[1]!=ZERO || m_x[4]!=ZERO || m_x[5]!=ZERO || m_x[8]!=ZERO || m_x[9]!=ZERO || m_x[12]!=ZERO ||
      m_x[13]!=ZERO) m_on|=1;
  if (m_x[2]!=ZERO || m_x[3]!=ZERO || m_x[6]!=ZERO || m_x[7]!=ZERO || m_x[10]!=ZERO || m_x[11]!=ZERO ||
      m_x[14]!=ZERO || m_x[15]!=ZERO) m_on|=2;
  return m_on&3;
}

// TODO: Springt er hier wirklich zum richtigen Operator? UNBEDINGT BEI ALLEN OPERATOREN TESTEN IM DEBUG-MODE
template<class Scalar>
bool CRaritaSchwinger<Scalar>::Test_Properties(const ATOOLS::Vec4D &p, int r, int b) {
  std::vector<SComplex> test_spinor(4, SComplex());
  test_spinor[0] = sqrt(p[0]-p.P())*(1/(sqrt(2*p.P()*(p.P()+p[3]))))*(p.P()+p[3]);
  test_spinor[1] = SComplex(sqrt(p[0]-p.P())*(1/(sqrt(2*p.P()*(p.P()+p[3])))))*(SComplex(p[1])+SComplex(0,p[2]));
  test_spinor[2] = sqrt(p[0]+p.P())*(1/(sqrt(2*p.P()*(p.P()+p[3]))))*(p.P()+p[3]);
  test_spinor[3] = SComplex(sqrt(p[0]+p.P())*(1/(sqrt(2*p.P()*(p.P()+p[3])))))*(SComplex(p[1])+SComplex(0,p[2]));
  std::vector<SComplex> test_polvec(4, SComplex());
  test_polvec[1] = SComplex(1/sqrt(2))*(SComplex(-(1/(p.P()*sqrt(p[1]*p[1]+p[2]*p[2])))*(p[1]*p[3]))+SComplex(0,1)*SComplex((1/(sqrt(p[1]*p[1]+p[2]*p[2])))*p[2]));
  test_polvec[2] = SComplex(1/sqrt(2))*(SComplex(-(1/(p.P()*sqrt(p[1]*p[1]+p[2]*p[2])))*(p[2]*p[3]))-SComplex(0,1)*SComplex((1/(sqrt(p[1]*p[1]+p[2]*p[2])))*p[1]));
  test_polvec[3] = SComplex(1/sqrt(2))*(SComplex((1/(p.P()*sqrt(p[1]*p[1]+p[2]*p[2])))*(sqrt(p[1]*p[1]+p[2]*p[2])*sqrt(p[1]*p[1]+p[2]*p[2]))));

  CRaritaSchwinger<Scalar> test_u(CRaritaSchwinger(1, -1, test_spinor[0]*test_polvec[0], test_spinor[1]*test_polvec[0],
                                                   test_spinor[2]*test_polvec[0], test_spinor[3]*test_polvec[0],
                                                   test_spinor[0]*test_polvec[1], test_spinor[1]*test_polvec[1],
                                                   test_spinor[2]*test_polvec[1], test_spinor[3]*test_polvec[1],
                                                   test_spinor[0]*test_polvec[2], test_spinor[1]*test_polvec[2],
                                                   test_spinor[2]*test_polvec[2], test_spinor[2]*test_polvec[3],
                                                   test_spinor[3]*test_polvec[0], test_spinor[3]*test_polvec[1],
                                                   test_spinor[3]*test_polvec[2], test_spinor[3]*test_polvec[3]));
  std::cout << "Test-u" << std::endl;
  // gamma_mu^A_B times RS^B,mu = 0
  std::cout<<METHOD<<": Testing gamma_mu Psi^mu = 0..."<<std::endl;
  METOOLS::Gamma gammavec1 = Gamma<Scalar>();
  std::vector<SComplex> result_test(4);
  for (size_t i(0); i<4; ++i){
    for (size_t j(0); j<4; ++j){
      if (b>0) result_test[i] += gammavec1[0][i][j] * test_u[j] - gammavec1[1][i][j] * test_u[j+4] -
                             gammavec1[2][i][j] * test_u[j+8] - gammavec1[3][i][j] * test_u[j+12];
      else{
        for (size_t k(0); k<4; ++k){
          result_test[k] += test_u[i] * gammavec1[0][i][j] * gammavec1[0][j][k] -
                        test_u[i+4] * gammavec1[1][i][j] * gammavec1[0][j][k] -
                        test_u[i+8] * gammavec1[2][i][j] * gammavec1[0][j][k] -
                        test_u[i+12] * gammavec1[3][i][j] * gammavec1[0][j][k];
        }
      }
    }
  }
  /*for (int i(0); i<4; ++i){
    result2[0] += gammavec[i][0][0] * (*this)[2*i] + gammavec[i][0][1] * (*this)[1+2*i] + gammavec[i][0][2] * (*this)[8+2*i] + gammavec[i][0][3] * (*this)[9+2*i];
    result2[1] += gammavec[i][1][0] * (*this)[2*i] + gammavec[i][1][1] * (*this)[1+2*i] + gammavec[i][1][2] * (*this)[8+2*i] + gammavec[i][1][3] * (*this)[9+2*i];
    result2[2] += gammavec[i][2][0] * (*this)[2*i] + gammavec[i][2][1] * (*this)[1+2*i] + gammavec[i][2][2] * (*this)[8+2*i] + gammavec[i][2][3] * (*this)[9+2*i];
    result2[3] += gammavec[i][3][0] * (*this)[2*i] + gammavec[i][3][1] * (*this)[1+2*i] + gammavec[i][3][2] * (*this)[8+2*i] + gammavec[i][3][3] * (*this)[9+2*i];
    // result2 += gammavec[i] * METOOLS::CVec4<Scalar>((*this)[i+i*4], (*this)[i+i*4+1], (*this)[i+i*4+2], (*this)[i+i*4+3]);
  }*/
  if (!(std::abs(result_test[0].real())<s_accu && std::abs(result_test[1].real())<s_accu && std::abs(result_test[2].real())<s_accu && std::abs(result_test[3].real())<s_accu)){
    msg_Out() << "gamma_mu Psi^mu is " << result_test[0] << result_test[1] << result_test[2] << result_test[3] << " not zero!" << std::endl;
  }
  if (!(std::abs(result_test[0].imag())<s_accu && std::abs(result_test[1].imag())<s_accu && std::abs(result_test[2].imag())<s_accu && std::abs(result_test[3].imag())<s_accu)){
    msg_Out() << "gamma_mu Psi^mu is " << result_test[0] << result_test[1] << result_test[2] << result_test[3] << " not zero!" << std::endl;
  }



  // Dirac equation (gamma_mu * p_mu -m)^A_B RS^B, nu -> auch +m? für V statt U?
  std::cout<<METHOD<<": Testing Dirac equation..."<<std::endl;
  METOOLS::Gamma gammavec = Gamma<Scalar>();
  bool testresult(true);
  // SComplex(-r) for distinguishing between particle and anti-particle Dirac equation
  ATOOLS::TCMatrix<Scalar> intermediate = ATOOLS::TCMatrix<Scalar>(gammavec * p);
  if (sqrt(p.Abs2())>1e-6)
    intermediate += ATOOLS::TCMatrix<Scalar>(SComplex(-r) * SComplex(sqrt(p.Abs2())) * ATOOLS::TCMatrix<Scalar>(4, true));
  std::vector<SComplex> result1(16);
  // TODO: Sollte man die exakten ZEROS auch hier explizit implementieren (zwei Komponenten jedes Dirac spinors und damit
  //       die Hälfte der Komponenten der RaSC sind ja leer im masselosen Fall!)
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
  /*for(int i(0); i<4; ++i) {
    result1[0+2*i] = intermediate[0][0] * (*this)[2*i] + intermediate[0][1] * (*this)[1+2*i] + intermediate[0][2] * (*this)[8+2*i] + intermediate[0][3] * (*this)[9+2*i];
    result1[1+2*i] = intermediate[1][0] * (*this)[2*i] + intermediate[1][1] * (*this)[1+2*i] + intermediate[1][2] * (*this)[8+2*i] + intermediate[1][3] * (*this)[9+2*i];
    result1[8+2*i] = intermediate[2][0] * (*this)[2*i] + intermediate[2][1] * (*this)[1+2*i] + intermediate[2][2] * (*this)[8+2*i] + intermediate[2][3] * (*this)[9+2*i];
    result1[9+2*i] = intermediate[3][0] * (*this)[2*i] + intermediate[3][1] * (*this)[1+2*i] + intermediate[3][2] * (*this)[8+2*i] + intermediate[3][3] * (*this)[9+2*i];
  }*/
   /*Implemented for filling when the four Dirac spinors are above each other in the spin-3/2 wave function
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
    if (std::abs(result1[j].real())>s_accu || std::abs(result1[j].imag())>s_accu) {
      msg_Out()<<"Component " << j << " of resulting Rarita-Schwinger wave function is " << result1[j] << " instead of zero!"
      << std::endl;
      //return false;
      testresult =false;
    }
  }
  if (testresult) msg_Out()<< "passed" << std::endl;

  // gamma_mu^A_B times RS^B,mu = 0
  std::cout<<METHOD<<": Testing gamma_mu Psi^mu = 0..."<<std::endl;
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
  /*for (int i(0); i<4; ++i){
    result2[0] += gammavec[i][0][0] * (*this)[2*i] + gammavec[i][0][1] * (*this)[1+2*i] + gammavec[i][0][2] * (*this)[8+2*i] + gammavec[i][0][3] * (*this)[9+2*i];
    result2[1] += gammavec[i][1][0] * (*this)[2*i] + gammavec[i][1][1] * (*this)[1+2*i] + gammavec[i][1][2] * (*this)[8+2*i] + gammavec[i][1][3] * (*this)[9+2*i];
    result2[2] += gammavec[i][2][0] * (*this)[2*i] + gammavec[i][2][1] * (*this)[1+2*i] + gammavec[i][2][2] * (*this)[8+2*i] + gammavec[i][2][3] * (*this)[9+2*i];
    result2[3] += gammavec[i][3][0] * (*this)[2*i] + gammavec[i][3][1] * (*this)[1+2*i] + gammavec[i][3][2] * (*this)[8+2*i] + gammavec[i][3][3] * (*this)[9+2*i];
    // result2 += gammavec[i] * METOOLS::CVec4<Scalar>((*this)[i+i*4], (*this)[i+i*4+1], (*this)[i+i*4+2], (*this)[i+i*4+3]);
  }*/
  if (!(std::abs(result2[0].real())<s_accu && std::abs(result2[1].real())<s_accu && std::abs(result2[2].real())<s_accu && std::abs(result2[3].real())<s_accu)){
    msg_Out() << "gamma_mu Psi^mu is " << result2[0] << result2[1] << result2[2] << result2[3] << " not zero!" << std::endl;
    testresult = false;
  }
  if (!(std::abs(result2[0].imag())<s_accu && std::abs(result2[1].imag())<s_accu && std::abs(result2[2].imag())<s_accu && std::abs(result2[3].imag())<s_accu)){
    msg_Out() << "gamma_mu Psi^mu is " << result2[0] << result2[1] << result2[2] << result2[3] << " not zero!" << std::endl;
    testresult = false;
  }
  if (testresult) msg_Out()<< "passed" << std::endl;

  // p_mu times RS = 0 ? aus partielle Ableitung_mu RS=0?
  std::cout<<METHOD<<": Testing p_mu Psi^mu = 0..."<<std::endl;
  std::vector<SComplex> result3(4);
  for (int i(0); i<4; ++i){
    result3[i] = SComplex(p[0]) * (*this)[i] + SComplex(-p[1]) * (*this)[4+i] + SComplex(-p[2]) * (*this)[8+i] + SComplex(-p[3]) * (*this)[12+i];
  }
  /*for (int i(0); i<2; ++i){
    result3[i] = SComplex(p[0]) * (*this)[i] + SComplex(p[1]) * (*this)[2+i] + SComplex(p[2]) * (*this)[4+i] + SComplex(p[3]) * (*this)[6+i];
    result3[i+2] = SComplex(p[0]) * (*this)[i+8] + SComplex(p[1]) * (*this)[10+i] + SComplex(p[2]) * (*this)[12+i] + SComplex (p[3]) * (*this)[14+i];
  }*/
  if (!(std::abs(result3[0].real())<s_accu && std::abs(result3[1].real())<s_accu && std::abs(result3[2].real())<s_accu && std::abs(result3[3].real())<s_accu)){
    msg_Out() << "p_mu Psi^mu is " << result3[0] << result3[1] << result3[2] << result3[3] << " not zero!" << std::endl;
    testresult = false;
  } ;
  if (!(std::abs(result3[0].imag())<s_accu && std::abs(result3[1].imag())<s_accu && std::abs(result3[2].imag())<s_accu && std::abs(result3[3].imag())<s_accu)) {
    msg_Out() << "p_mu Psi^mu is " << result3[0] << result3[1] << result3[2] << result3[3] << " not zero!" << std::endl;
    testresult = false;
  }
  if (testresult) msg_Out()<< "passed" << std::endl;
  return testresult;
}

namespace METOOLS {

  template class DCRaritaSchwinger;
  template std::ostream &operator<<(std::ostream &ostr,const DCRaritaSchwinger &s);

  template class QCRaritaSchwinger;
  template std::ostream &operator<<(std::ostream &ostr,const QCRaritaSchwinger &s);

}
