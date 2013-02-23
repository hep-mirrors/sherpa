#include "METOOLS/Explicit/Lorentz_Calculator.H"

#include "METOOLS/Explicit/Vertex.H"
#include "METOOLS/Currents/C_Scalar.H"
#include "METOOLS/Currents/C_Vector.H"
#include "MODEL/Interaction_Models/Vertex.H"

namespace METOOLS {

  template <typename SType>
  class SSV_Calculator: public Lorentz_Calculator {
  public:

    typedef std::complex<SType> SComplex;

    typedef CVec4<SType> CVec4Type;
    typedef std::vector<CVec4Type*> CVec4Type_Vector;

    typedef CScalar<SType> CScalarType;
    typedef std::vector<CScalarType*> CScalarType_Vector;

  private:

    SComplex m_cpl;

    int m_dir, m_swap;
    
  public:
    
    SSV_Calculator(const Vertex_Key &key);

    std::string Label() const;
    
    void Evaluate();

  };// end of class SSV_Calculator

}// end of namespace METOOLS

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"

#define ZERO SComplex(0.0,0.0)

using namespace METOOLS;
using namespace ATOOLS;

template <typename SType>
SSV_Calculator<SType>::SSV_Calculator(const Vertex_Key &key):
  Lorentz_Calculator(key),  
  m_dir(key.p_b->Flav().IsScalar()?
	(key.p_a->Flav().IsScalar()?0:2):1) 
{
  if (key.p_mv->Lorentz[key.m_n]->
      Type().find("SSV")==std::string::npos)
    THROW(not_implemented,"Invalid Lorentz structure");
  int a(key.p_mv->Lorentz[key.m_n]->ParticleArg(0));
  int b(key.p_mv->Lorentz[key.m_n]->ParticleArg(1));
  m_swap=0;
  if (m_dir==1) m_swap=a>b;
  else m_swap=a<b;
  m_cpl=SComplex(p_v->Coupling(0)*p_cc->Coupling());
}
 
template <typename SType>
void SSV_Calculator<SType>::Evaluate()
{
  p_v->SetZero();
  if (p_v->JA()->Zero()||p_v->JB()->Zero()) return;
#ifdef DEBUG__BG
  msg_Debugging()<<*p_v->JA()<<"(+)"<<*p_v->JB()<<" SSV("<<m_dir
		 <<"): m_cpl = "<<m_cpl<<"\n";
  msg_Indent();
#endif
  size_t i(0);
  for (typename CObject_Matrix::const_iterator 
	 jait(p_v->JA()->J().begin());jait!=p_v->JA()->J().end();++jait) {
    for (typename CObject_Matrix::const_iterator 
	   jbit(p_v->JB()->J().begin());jbit!=p_v->JB()->J().end();++jbit,++i) {
  switch(m_dir) {
  case 0: {
    Vec4D p(p_v->JB()->P()-p_v->JA()->P());
    const CScalarType_Vector *ca(jait->Get<CScalarType>());
    const CScalarType_Vector *cb(jbit->Get<CScalarType>());
    for (typename CScalarType_Vector::const_iterator 
	   ait(ca->begin());ait!=ca->end();++ait)
      for (typename CScalarType_Vector::const_iterator 
	     bit(cb->begin());bit!=cb->end();++bit)
	if (p_cc->Evaluate(*ait,*bit)) {
#ifdef DEBUG__BG
	  msg_Debugging()<<"    "<<**ait<<"\n";
	  msg_Debugging()<<"    "<<**bit<<"\n";
#endif
	  CVec4Type *j(CVec4Type::New(p[0],p[1],p[2],p[3]));
	  *j*=(**ait)[0]*(**bit)[0]*SComplex(m_swap?-m_cpl:m_cpl);
	  j->SetH(p_v->H(i));
	  j->SetS((*ait)->S()|(*bit)->S());
	  p_cc->AddJ(j);
	  p_v->SetZero(false);
	}
    break;
  }
  case 1: {
    CVec4Type p(2.0*p_v->JA()->P()+p_v->JB()->P());
    const CScalarType_Vector *ca(jait->Get<CScalarType>());
    const CVec4Type_Vector *cb(jbit->Get<CVec4Type>());
    for (typename CScalarType_Vector::const_iterator 
	   ait(ca->begin());ait!=ca->end();++ait)
      for (typename CVec4Type_Vector::const_iterator 
	     bit(cb->begin());bit!=cb->end();++bit)
	if (p_cc->Evaluate(*ait,*bit)) {
#ifdef DEBUG__BG
	  msg_Debugging()<<"    "<<**ait<<"\n";
	  msg_Debugging()<<"    "<<**bit<<"\n";
#endif
	  CScalarType *j(CScalarType::New(**ait));
	  *j*=SComplex(m_swap?-m_cpl:m_cpl)*(**bit*p);
	  j->SetH(p_v->H(i));
	  j->SetS((*ait)->S()|(*bit)->S());
	  p_cc->AddJ(j);
	  p_v->SetZero(false);
	}
    break;
  }
  case 2: {
    CVec4Type p(-2.0*p_v->JB()->P()-p_v->JA()->P());
    const CVec4Type_Vector *ca(jait->Get<CVec4Type>());
    const CScalarType_Vector *cb(jbit->Get<CScalarType>());
    for (typename CVec4Type_Vector::const_iterator 
	   ait(ca->begin());ait!=ca->end();++ait)
      for (typename CScalarType_Vector::const_iterator 
	     bit(cb->begin());bit!=cb->end();++bit)
	if (p_cc->Evaluate(*ait,*bit)) {
#ifdef DEBUG__BG
	  msg_Debugging()<<"    "<<**bit<<"\n";
	  msg_Debugging()<<"    "<<**ait<<"\n";
#endif
	  CScalarType *j(CScalarType::New(**bit));
	  *j*=SComplex(m_swap?-m_cpl:m_cpl)*(**ait*p);
	  j->SetH(p_v->H(i));
	  j->SetS((*ait)->S()|(*bit)->S());
	  p_cc->AddJ(j);
	  p_v->SetZero(false);
	}
    break;
  }
  default:
    THROW(fatal_error,"Internal error");
  }
    }
  }
}
 
template <typename SType>
std::string SSV_Calculator<SType>::Label() const
{
  return "SSV["+ToString(m_cpl)+"]";
}

namespace METOOLS {

  template class SSV_Calculator<double>;

}

DECLARE_GETTER(SSV_Calculator<double>,"DSSV",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,SSV_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new SSV_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    SSV_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"SSV vertex"; }
