#include "METOOLS/Explicit/Lorentz_Calculator.H"

#include "METOOLS/Explicit/Vertex.H"
#include "METOOLS/Currents/C_Scalar.H"
#include "METOOLS/Currents/C_Spinor.H"

namespace METOOLS {

  template <typename SType>
  class SFF_Calculator: public Lorentz_Calculator {
  public:

    typedef std::complex<SType> SComplex;

    typedef CSpinor<SType> CSpinorType;
    typedef std::vector<CSpinorType*> CSpinorType_Vector;

    typedef CScalar<SType> CScalarType;
    typedef std::vector<CScalarType*> CScalarType_Vector;

    CScalar<SType> LorentzLeft(CSpinorType a,CSpinorType b);
    CScalar<SType> LorentzRight(CSpinorType a,CSpinorType b);

    CSpinor<SType> LorentzLeft(const CSpinorType &a,const CScalarType &b);
    CSpinor<SType> LorentzRight(const CSpinorType &a,const CScalarType &b);

  private:

    SComplex m_cpll, m_cplr;

    int m_dir;
    
  public:
    
    SFF_Calculator(const Vertex_Key &key);

    std::string Label() const;
    
    void Evaluate();

  };// end of class SFF_Calculator

}// end of namespace METOOLS

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"

#define ZERO SComplex(0.0,0.0)

using namespace METOOLS;
using namespace ATOOLS;

template <typename SType>
SFF_Calculator<SType>::SFF_Calculator(const Vertex_Key &key):
  Lorentz_Calculator(key),  
  m_dir(key.p_b->Flav().IsFermion()?
	(key.p_a->Flav().IsFermion()?0:2):1) 
{
  m_cpll=SComplex(-p_v->Coupling(0)*p_cc->Coupling());
  m_cplr=SComplex(-p_v->Coupling(1)*p_cc->Coupling());
}

 
template <typename SType>
void SFF_Calculator<SType>::Evaluate()
{
  p_v->SetZero();
  if (p_v->JA()->Zero()||p_v->JB()->Zero()) return;
#ifdef DEBUG__BG
  msg_Debugging()<<*p_v->JA()<<"(+)"<<*p_v->JB()<<" SFF("<<m_dir
		 <<"): m_cpll = "<<m_cpll<<", m_cplr = "<<m_cplr<<"\n";
  msg_Indent();
#endif
  size_t i(0);
  for (typename CObject_Matrix::const_iterator 
	 jait(p_v->JA()->J().begin());jait!=p_v->JA()->J().end();++jait) {
    for (typename CObject_Matrix::const_iterator 
	   jbit(p_v->JB()->J().begin());jbit!=p_v->JB()->J().end();++jbit,++i) {
  switch(m_dir) {
  case 0: {
    const CSpinorType_Vector *ca(jait->Get<CSpinorType>());
    const CSpinorType_Vector *cb(jbit->Get<CSpinorType>());
    for (typename CSpinorType_Vector::const_iterator 
	   ait(ca->begin());ait!=ca->end();++ait)
      for (typename CSpinorType_Vector::const_iterator 
	     bit(cb->begin());bit!=cb->end();++bit)
	if (p_cc->Evaluate(*ait,*bit)) {
#ifdef DEBUG__BG
	  msg_Debugging()<<"    "<<**ait<<"\n";
	  msg_Debugging()<<"    "<<**bit<<"\n";
#endif
	  CScalarType *j(CScalarType::New
			 (LorentzLeft(**ait,**bit)+
			  LorentzRight(**ait,**bit)));
	  j->SetH(p_v->H(i));
	  p_cc->AddJ(j);
	  p_v->SetZero(false);
	}
    break;
  }
  case 1: {
    const CSpinorType_Vector *ca(jait->Get<CSpinorType>());
    const CScalarType_Vector *cb(jbit->Get<CScalarType>());
    for (typename CSpinorType_Vector::const_iterator 
	   ait(ca->begin());ait!=ca->end();++ait)
      for (typename CScalarType_Vector::const_iterator 
	     bit(cb->begin());bit!=cb->end();++bit)
	if (p_cc->Evaluate(*ait,*bit)) {
#ifdef DEBUG__BG
	  msg_Debugging()<<"    "<<**ait<<"\n";
	  msg_Debugging()<<"    "<<**bit<<"\n";
#endif
	  CSpinorType *j(CSpinorType::New
			 (LorentzLeft(**ait,**bit)+
			  LorentzRight(**ait,**bit)));
	  j->SetH(p_v->H(i));
	  p_cc->AddJ(j);
	  p_v->SetZero(false);
	}
    break;
  }
  case 2: {
    const CScalarType_Vector *ca(jait->Get<CScalarType>());
    const CSpinorType_Vector *cb(jbit->Get<CSpinorType>());
    for (typename CScalarType_Vector::const_iterator 
	   ait(ca->begin());ait!=ca->end();++ait)
      for (typename CSpinorType_Vector::const_iterator 
	     bit(cb->begin());bit!=cb->end();++bit)
	if (p_cc->Evaluate(*ait,*bit)) {
#ifdef DEBUG__BG
	  msg_Debugging()<<"    "<<**bit<<"\n";
	  msg_Debugging()<<"    "<<**ait<<"\n";
#endif
	  CSpinorType *j(CSpinorType::New
			 (LorentzLeft(**bit,**ait)+
			  LorentzRight(**bit,**ait)));
	  j->SetH(p_v->H(i));
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
 
template <typename SType> CScalar<SType>
SFF_Calculator<SType>::LorentzLeft(CSpinorType a,CSpinorType b)
{
  if (a.B()>0) std::swap<CSpinorType>(a,b);
#ifdef DEBUG__BG
  msg_Debugging()<<"<> L "<<a<<"\n";
  msg_Debugging()<<"     "<<b<<"\n";
#endif
  return CScalarType((a[2]*b[2]+a[3]*b[3])*m_cpll,0,a.S()|b.S());
}

template <typename SType> CScalar<SType>
SFF_Calculator<SType>::LorentzRight(CSpinorType a,CSpinorType b)
{
  if (a.B()>0) std::swap<CSpinorType>(a,b);
#ifdef DEBUG__BG
  msg_Debugging()<<"<> R "<<a<<"\n";
  msg_Debugging()<<"     "<<b<<"\n";
#endif
  return CScalarType((a[0]*b[0]+a[1]*b[1])*m_cplr,0,a.S()|b.S());
}

template <typename SType> CSpinor<SType>
SFF_Calculator<SType>::LorentzLeft(const CSpinorType &a,
				   const CScalarType &b)
{
  switch (a.B()) {
  case -1: {
#ifdef DEBUG__BG
    msg_Debugging()<<"<| R "<<a<<"\n";
    msg_Debugging()<<"     "<<b<<"\n";
#endif
    CSpinorType j(a.R(),a.B(),0,0,a.H()|b.H(),a.S()|b.S(),1);
    j[0]=a[0]*b[0]*m_cplr;
    j[1]=a[1]*b[0]*m_cplr;
    j[3]=j[2]=SComplex(0.0,0.0);
    return j;
  }
  case 1: {
#ifdef DEBUG__BG
    msg_Debugging()<<"|> L "<<a<<"\n";
    msg_Debugging()<<"     "<<b<<"\n";
#endif
    CSpinorType j(a.R(),a.B(),0,0,a.H()|b.H(),a.S()|b.S(),2);
    j[1]=j[0]=SComplex(0.0,0.0);
    j[2]=a[2]*b[0]*m_cpll;
    j[3]=a[3]*b[0]*m_cpll;
    return j;
  }
  default:
    THROW(fatal_error,"Internal error");
  }
  return CSpinorType();
}

template <typename SType> CSpinor<SType>
SFF_Calculator<SType>::LorentzRight(const CSpinorType &a,
				    const CScalarType &b)
{
  switch (a.B()) {
  case -1: {
#ifdef DEBUG__BG
    msg_Debugging()<<"<| L "<<a<<"\n";
    msg_Debugging()<<"     "<<b<<"\n";
#endif
    CSpinorType j(a.R(),a.B(),0,0,a.H()|b.H(),a.S()|b.S(),2);
    j[2]=a[2]*b[0]*m_cpll;
    j[3]=a[3]*b[0]*m_cpll;
    j[1]=j[0]=SComplex(0.0,0.0);
    return j;
  }
  case 1: {
#ifdef DEBUG__BG
    msg_Debugging()<<"|> R "<<a<<"\n";
    msg_Debugging()<<"     "<<b<<"\n";
#endif
    CSpinorType j(a.R(),a.B(),0,0,a.H()|b.H(),a.S()|b.S(),1);
    j[3]=j[2]=SComplex(0.0,0.0);
    j[0]=a[0]*b[0]*m_cplr;
    j[1]=a[1]*b[0]*m_cplr;
    return j;
  }
  default:
    THROW(fatal_error,"Internal error");
  }
  return CSpinorType();
}

template <typename SType>
std::string SFF_Calculator<SType>::Label() const
{
  return "SFF["+ToString(m_cpll)+","+ToString(m_cplr)+"]";
}

namespace METOOLS {

  template class SFF_Calculator<double>;

}

DECLARE_GETTER(SFF_Calculator<double>,"DFFS",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,SFF_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new SFF_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    SFF_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"SFF vertex"; }
