#include "METOOLS/Explicit/Lorentz_Calculator.H"

#include "METOOLS/Explicit/Vertex.H"
#include "METOOLS/Currents/C_Scalar.H"
#include "METOOLS/Currents/C_Spinor.H"

namespace METOOLS {

  template <typename SType>
  class TAUPI_Calculator: public Lorentz_Calculator {
  public:

    typedef std::complex<SType> SComplex;

    typedef CSpinor<SType> CSpinorType;
    typedef std::vector<CSpinorType*> CSpinorType_Vector;

    typedef CScalar<SType> CScalarType;
    typedef std::vector<CScalarType*> CScalarType_Vector;

    typedef CVec4<SType> CVec4Type;
    typedef std::vector<CVec4Type*> CVec4Type_Vector;

  private:

    SComplex m_cpll, m_cplr;

    int m_dir, m_cl, m_cr;
    
    CSpinorType LorentzLeft(const CSpinorType &a, const CVec4Type &b);
    CSpinorType LorentzRight(const CSpinorType &a,const CVec4Type &b);

  public:
    
    TAUPI_Calculator(const Vertex_Key &key);

    std::string Label() const;
    
    void Evaluate();

    inline SComplex PPlus(const CVec4<SType> &p) const
    { return p[0]+p[ATOOLS::Spinor<SType>::R3()]; }
    inline SComplex PMinus(const CVec4<SType> &p) const
    { return p[0]-p[ATOOLS::Spinor<SType>::R3()]; }

    inline SComplex PT(const CVec4<SType> &p) const
    { return p[ATOOLS::Spinor<SType>::R1()]+
	SComplex(0.0,1.0)*p[ATOOLS::Spinor<SType>::R2()]; }
    inline SComplex PTC(const CVec4<SType> &p) const
    { return p[ATOOLS::Spinor<SType>::R1()]-
	SComplex(0.0,1.0)*p[ATOOLS::Spinor<SType>::R2()]; }

  };// end of class TAUPI_Calculator

}// end of namespace METOOLS

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"

#define ZERO SComplex(0.0,0.0)

using namespace METOOLS;
using namespace ATOOLS;

template <typename SType>
TAUPI_Calculator<SType>::TAUPI_Calculator(const Vertex_Key &key):
  Lorentz_Calculator(key),  
  m_dir(key.p_b->Flav().IsFermion()?
	(key.p_a->Flav().IsFermion()?0:2):1)
{
  if (m_dir!=0 && key.p_c->Flav().IsAnti()) {
    m_cpll=SComplex(-p_v->Coupling(0)*p_cc->Coupling());
    m_cplr=SComplex(-p_v->Coupling(1)*p_cc->Coupling());
  }
  else {
    m_cpll=SComplex(p_v->Coupling(1)*p_cc->Coupling());
    m_cplr=SComplex(p_v->Coupling(0)*p_cc->Coupling());
  }
}
 
template <typename SType>
void TAUPI_Calculator<SType>::Evaluate()
{
  p_v->SetZero();
  if (p_v->JA()->Zero()||p_v->JB()->Zero()) return;
  size_t i(0);
  for (typename CObject_Matrix::const_iterator 
	 jait(p_v->JA()->J().begin());jait!=p_v->JA()->J().end();++jait) {
    for (typename CObject_Matrix::const_iterator 
	   jbit(p_v->JB()->J().begin());jbit!=p_v->JB()->J().end();++jbit,++i) {
  switch(m_dir) {
  case 0: {
    THROW(not_implemented, "Case 0.")
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
          CSpinorType *j(CSpinorType::New(CSpinorType((**ait).R(),(**ait).B(),(**ait)(0),(**ait)(1),0,0,0)));
          *j+=LorentzLeft((**ait),CVec4Type(p_v->JB()->P()))*m_cpll;
          *j+=LorentzRight((**ait),CVec4Type(p_v->JB()->P()))*m_cplr;
          *j*=(**bit)[0];
	  j->SetH(p_v->H(i));
	  p_cc->AddJ(j);
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
          CSpinorType *j(CSpinorType::New(CSpinorType((**bit).R(),(**bit).B(),(**bit)(0),(**bit)(1),0,0,0)));
          *j+=LorentzLeft((**bit),CVec4Type(p_v->JA()->P()))*m_cpll;
          *j+=LorentzRight((**bit),CVec4Type(p_v->JA()->P()))*m_cplr;
          *j*=(**ait)[0];
	  j->SetH(p_v->H(i));
	  p_cc->AddJ(j);
	}
    break;
  }
  default:
    THROW(fatal_error,"Internal error");
  }
    }
  }
}

template <typename SType> CSpinor<SType>
TAUPI_Calculator<SType>::LorentzRight(const CSpinorType &a,
				   const CVec4Type &b)
{
  switch (a.B()) {
  case -1: {
#ifdef DEBUG__BG
    msg_Debugging()<<"<|g L "<<a<<"\n";
    msg_Debugging()<<"      "<<b<<"\n";
#endif
    CSpinorType j(a.R(),a.B(),0,0,0,a.S()|b.S(),1);
    SComplex jp(PPlus(b)), jm(PMinus(b)), jt(PT(b)), jtc(PTC(b));
    j[0]=a[2]*jp+a[3]*jt;
    j[1]=a[2]*jtc+a[3]*jm;
    j[3]=j[2]=SComplex(0.0,0.0);
    return j;
  }
  case 1: {
#ifdef DEBUG__BG
    msg_Debugging()<<"g|> L "<<a<<"\n";
    msg_Debugging()<<"      "<<b<<"\n";
#endif
    CSpinorType j(a.R(),a.B(),0,0,0,a.S()|b.S(),2);
    SComplex jp(PPlus(b)), jm(PMinus(b)), jt(PT(b)), jtc(PTC(b));
    j[1]=j[0]=SComplex(0.0,0.0);
    j[2]=a[0]*jp+a[1]*jtc;
    j[3]=a[0]*jt+a[1]*jm;
    return j;
  }
  default:
    THROW(fatal_error,"Internal error");
  }
  return CSpinorType();
}

template <typename SType> CSpinor<SType>
TAUPI_Calculator<SType>::LorentzLeft(const CSpinorType &a,
                                      const CVec4Type &b)
{
  switch (a.B()) {
  case -1: {
#ifdef DEBUG__BG
    msg_Debugging()<<"<|g R "<<a<<"\n";
    msg_Debugging()<<"      "<<b<<"\n";
#endif
    CSpinorType j(a.R(),a.B(),0,0,0,a.S()|b.S(),2);
    SComplex jp(PPlus(b)), jm(PMinus(b)), jt(PT(b)), jtc(PTC(b));
    j[1]=j[0]=SComplex(0.0,0.0);
    j[2]=a[0]*jm-a[1]*jt;
    j[3]=-a[0]*jtc+a[1]*jp;
    return j;
  }
  case 1: {
#ifdef DEBUG__BG
    msg_Debugging()<<"g|> R "<<a<<"\n";
    msg_Debugging()<<"      "<<b<<"\n";
#endif
    CSpinorType j(a.R(),a.B(),0,0,0,a.S()|b.S(),1);
    SComplex jp(PPlus(b)), jm(PMinus(b)), jt(PT(b)), jtc(PTC(b));
    j[0]=a[2]*jm-a[3]*jtc;
    j[1]=-a[2]*jt+a[3]*jp;
    j[3]=j[2]=SComplex(0.0,0.0);
    return j;
  }
  default:
    THROW(fatal_error,"Internal error");
  }
  return CSpinorType();
}


 
template <typename SType>
std::string TAUPI_Calculator<SType>::Label() const
{
  return "TAUPI["+ToString(m_cpll)+",,"+ToString(m_cplr)+"]";
}

namespace METOOLS {

  template class TAUPI_Calculator<double>;

}

DECLARE_GETTER(TAUPI_Calculator<double>,"DTAUPI",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,TAUPI_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new TAUPI_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    TAUPI_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"TAUPI vertex"; }
