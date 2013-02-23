#include "METOOLS/Explicit/Lorentz_Calculator.H"

#include "METOOLS/Explicit/Vertex.H"
#include "METOOLS/Currents/C_Tensor.H"

namespace METOOLS {

  template <typename SType>
  class VVP_Calculator: public Lorentz_Calculator {
  public:

    typedef std::complex<SType> SComplex;

    typedef CAsT4<SType> CAsT4Type;
    typedef std::vector<CAsT4Type*> CAsT4Type_Vector;

    typedef CVec4<SType> CVec4Type;
    typedef std::vector<CVec4Type*> CVec4Type_Vector;

  protected:

    SComplex m_cpl;

    int m_dir;
    
  public:
    
    VVP_Calculator(const Vertex_Key &key);

    std::string Label() const;
    
    void Evaluate();

  };// end of class VVP_Calculator

}// end of namespace METOOLS

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"

#define ZERO SComplex(0.0,0.0)

using namespace METOOLS;
using namespace ATOOLS;

template <typename SType>
VVP_Calculator<SType>::VVP_Calculator(const Vertex_Key &key):
  Lorentz_Calculator(key),  
  m_dir(key.p_b->Flav().IsVector()?
	(key.p_a->Flav().IsVector()?0:2):1)
{
  m_cpl=SComplex(p_v->Coupling(0)*p_cc->Coupling());
}
 
template <typename SType>
void VVP_Calculator<SType>::Evaluate()
{
  p_v->SetZero();
  if (p_v->JA()->Zero()||p_v->JB()->Zero()) return;
#ifdef DEBUG__BG
  msg_Debugging()<<*p_v->JA()<<"(+)"<<*p_v->JB()<<" VVP("<<m_dir<<")\n";
  msg_Indent();
#endif
  size_t i(0);
  for (typename CObject_Matrix::const_iterator 
	 jait(p_v->JA()->J().begin());jait!=p_v->JA()->J().end();++jait) {
    for (typename CObject_Matrix::const_iterator 
	   jbit(p_v->JB()->J().begin());jbit!=p_v->JB()->J().end();++jbit,++i) {
  switch(m_dir) {
  case 0: {
    const CVec4Type_Vector *ca(jait->Get<CVec4Type>());
    const CVec4Type_Vector *cb(jbit->Get<CVec4Type>());
    for (typename CVec4Type_Vector::const_iterator 
	   ait(ca->begin());ait!=ca->end();++ait)
      for (typename CVec4Type_Vector::const_iterator 
	     bit(cb->begin());bit!=cb->end();++bit)
	if (p_cc->Evaluate(*ait,*bit)) {
#ifdef DEBUG__BG
	  msg_Debugging()<<"    "<<**ait<<"\n";
	  msg_Debugging()<<"    "<<**bit<<"\n";
#endif
	  CAsT4Type *j(CAsT4Type::New(CAsT4Type(**ait,**bit)*m_cpl));
	  j->SetH(p_v->H(i));
	  j->SetS((*ait)->S()|(*bit)->S());
	  p_cc->AddJ(j);
	  p_v->SetZero(false);
	}
    break;
  }
  case 1: {
    const CVec4Type_Vector *ca(jait->Get<CVec4Type>());
    const CAsT4Type_Vector *cb(jbit->Get<CAsT4Type>());
    for (typename CVec4Type_Vector::const_iterator 
	   ait(ca->begin());ait!=ca->end();++ait)
      for (typename CAsT4Type_Vector::const_iterator 
	     bit(cb->begin());bit!=cb->end();++bit)
	if (p_cc->Evaluate(*ait,*bit)) {
#ifdef DEBUG__BG
	  msg_Debugging()<<"    "<<**ait<<"\n";
	  msg_Debugging()<<"    "<<**bit<<"\n";
#endif
	  CVec4Type *j(CVec4Type::New((**ait***bit)*m_cpl));
	  j->SetH(p_v->H(i));
	  j->SetS((*ait)->S()|(*bit)->S());
	  p_cc->AddJ(j);
	  p_v->SetZero(false);
	}
    break;
  }
  case 2: {
    const CAsT4Type_Vector *ca(jait->Get<CAsT4Type>());
    const CVec4Type_Vector *cb(jbit->Get<CVec4Type>());
    for (typename CAsT4Type_Vector::const_iterator 
	   ait(ca->begin());ait!=ca->end();++ait)
      for (typename CVec4Type_Vector::const_iterator 
	     bit(cb->begin());bit!=cb->end();++bit)
	if (p_cc->Evaluate(*ait,*bit)) {
#ifdef DEBUG__BG
	  msg_Debugging()<<"    "<<**bit<<"\n";
	  msg_Debugging()<<"    "<<**ait<<"\n";
#endif
	  CVec4Type *j(CVec4Type::New((**bit***ait)*m_cpl));
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
std::string VVP_Calculator<SType>::Label() const
{
  return "VVP["+ToString(m_cpl)+"]";
}

namespace METOOLS {

  template class VVP_Calculator<double>;

}

DECLARE_GETTER(VVP_Calculator<double>,"DGaugeP4",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,VVP_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new VVP_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    VVP_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"VVP vertex"; }
