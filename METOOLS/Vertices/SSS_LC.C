#include "METOOLS/Explicit/Lorentz_Calculator.H"

#include "METOOLS/Explicit/Vertex.H"
#include "METOOLS/Currents/C_Scalar.H"
#include "ATOOLS/Org/MyStrStream.H"

namespace METOOLS {

  template <typename SType>
  class SSS_Calculator: public Lorentz_Calculator {
  public:

    typedef std::complex<SType> SComplex;

    typedef CScalar<SType> CScalarType;
    typedef std::vector<CScalarType*> CScalarType_Vector;

  private:

    SComplex m_cpl;

  public:
    
    SSS_Calculator(const Vertex_Key &key);
    
    std::string Label() const;

    void Evaluate();

  };// end of class SSS_Calculator

}// end of namespace METOOLS

#include "ATOOLS/Org/Message.H"

using namespace METOOLS;
using namespace ATOOLS;

template <typename SType>
SSS_Calculator<SType>::SSS_Calculator(const Vertex_Key &key): 
  Lorentz_Calculator(key) 
{
  m_cpl=SComplex(p_v->Coupling(0)*p_cc->Coupling());
}

template <typename SType>
std::string SSS_Calculator<SType>::Label() const
{
  return "SSS["+ToString(m_cpl)+"]";
}

template <typename SType>
void SSS_Calculator<SType>::Evaluate()
{
  p_v->SetZero();
  if (p_v->JA()->Zero()||p_v->JB()->Zero()) return;
#ifdef DEBUG__BG
  msg_Debugging()<<*p_v->JA()<<"(+)"<<*p_v->JB()<<" SSS\n";
  msg_Indent();
#endif
  size_t i(0);
  const CObject_Matrix &cca(p_v->JA()->J()), &ccb(p_v->JB()->J());
  for (typename CObject_Matrix::const_iterator 
	 jait(cca.begin());jait!=cca.end();++jait) {
    const CScalarType_Vector *ca(jait->Get<CScalarType>());
    for (typename CObject_Matrix::const_iterator 
	   jbit(ccb.begin());jbit!=ccb.end();++jbit,++i) {
      const CScalarType_Vector *cb(jbit->Get<CScalarType>());
      for (typename CScalarType_Vector::const_iterator 
	     bit(cb->begin());bit!=cb->end();++bit)
	for (typename CScalarType_Vector::const_iterator 
	       ait(ca->begin());ait!=ca->end();++ait)
	  if (p_cc->Evaluate(*ait,*bit)) {
#ifdef DEBUG__BG
	    msg_Debugging()<<"    "<<**ait<<"\n";
	    msg_Debugging()<<"    "<<**bit<<"\n";
#endif
	    CScalarType *j(CScalarType::New(**ait***bit));
	    *j*=SComplex(m_cpl);
	    j->SetH(p_v->H(i));
	    j->SetS((*ait)->S()|(*bit)->S());
	    p_cc->AddJ(j);
	    p_v->SetZero(false);
	  }
    }
  }
}

namespace METOOLS {

  template class SSS_Calculator<double>;
  template class SSS_Calculator<long double>;

}

DECLARE_TEMPLATE_GETTER(SSS_Getter,"SSS",Lorentz_Calculator,Vertex_Key);

template <typename SType,char STag> Lorentz_Calculator *
SSS_Getter<SType,STag>::operator()(const Vertex_Key &key) const
{
  return new SSS_Calculator<SType>(key);
}

template <typename SType,char STag>
void SSS_Getter<SType,STag>::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"SSS vertex";
}

template class SSS_Getter<double,'D'>;
template class SSS_Getter<long double,'Q'>;
