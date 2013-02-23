#include "METOOLS/Explicit/Lorentz_Calculator.H"

#include "METOOLS/Explicit/Vertex.H"
#include "METOOLS/Currents/C_Scalar.H"
#include "ATOOLS/Org/MyStrStream.H"

namespace METOOLS {

  template <typename SType>
  class SSSS_Calculator: public Lorentz_Calculator {
  public:

    typedef std::complex<SType> SComplex;

    typedef CScalar<SType> CScalarType;
    typedef std::vector<CScalarType*> CScalarType_Vector;

  private:

    SComplex m_cpl;

  public:
    
    SSSS_Calculator(const Vertex_Key &key);
    
    std::string Label() const;

    void Evaluate();

  };// end of class SSSS_Calculator

}// end of namespace METOOLS

#include "ATOOLS/Org/Message.H"

using namespace METOOLS;
using namespace ATOOLS;

template <typename SType>
SSSS_Calculator<SType>::SSSS_Calculator(const Vertex_Key &key): 
  Lorentz_Calculator(key) 
{
  m_cpl=SComplex(p_v->Coupling(0)*p_cc->Coupling());
}

template <typename SType>
std::string SSSS_Calculator<SType>::Label() const
{
  return "SSSS["+ToString(m_cpl)+"]";
}

template <typename SType>
void SSSS_Calculator<SType>::Evaluate()
{
  p_v->SetZero();
  if (p_v->JA()->Zero()||p_v->JB()->Zero()||p_v->JE()->Zero()) return;
#ifdef DEBUG__BG
  msg_Debugging()<<*p_v->JA()<<"(+)"<<*p_v->JB()<<"(+)"<<*p_v->JE()<<" SSSS\n";
  msg_Indent();
#endif
  size_t i(0);
  const CObject_Matrix &cca(p_v->JA()->J()),
    &ccb(p_v->JB()->J()), &cce(p_v->JE()->J());
  for (typename CObject_Matrix::const_iterator 
	 jait(cca.begin());jait!=cca.end();++jait) {
    const CScalarType_Vector *ca(jait->Get<CScalarType>());
    for (typename CObject_Matrix::const_iterator 
	   jbit(ccb.begin());jbit!=ccb.end();++jbit) {
      const CScalarType_Vector *cb(jbit->Get<CScalarType>());
      for (typename CObject_Matrix::const_iterator 
	     jeit(cce.begin());jeit!=cce.end();++jeit,++i) {
	const CScalarType_Vector *ce(jeit->Get<CScalarType>());
	for (typename CScalarType_Vector::const_iterator 
	       eit(ce->begin());eit!=ce->end();++eit)
	  for (typename CScalarType_Vector::const_iterator 
		 bit(cb->begin());bit!=cb->end();++bit)
	    for (typename CScalarType_Vector::const_iterator 
		   ait(ca->begin());ait!=ca->end();++ait)
	      if (p_cc->Evaluate(*ait,*bit,*eit)) {
#ifdef DEBUG__BG
		msg_Debugging()<<"    "<<**ait<<"\n";
		msg_Debugging()<<"    "<<**bit<<"\n";
		msg_Debugging()<<"    "<<**eit<<"\n";
#endif
		CScalarType *j(CScalarType::New(**ait***bit***eit));
		*j*=SComplex(m_cpl);
		j->SetH(p_v->H(i));
		j->SetS((*ait)->S()|(*bit)->S()|(*eit)->S());
		p_cc->AddJ(j);
		p_v->SetZero(false);
	      }
      }
    }
  }
}

namespace METOOLS {

  template class SSSS_Calculator<double>;

}

DECLARE_GETTER(SSSS_Calculator<double>,"DSSSS",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,SSSS_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new SSSS_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    SSSS_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"SSSS vertex"; }
