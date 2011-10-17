#include "METOOLS/Explicit/Lorentz_Calculator.H"

#include "METOOLS/Explicit/Vertex.H"
#include "METOOLS/Currents/C_Scalar.H"
#include "METOOLS/Currents/C_Vector.H"
#include "ATOOLS/Org/MyStrStream.H"

namespace METOOLS {

  template <typename SType>
  class SSVV_Calculator: public Lorentz_Calculator {
  public:

    typedef std::complex<SType> SComplex;

    typedef CScalar<SType> CScalarType;
    typedef std::vector<CScalarType*> CScalarType_Vector;
    typedef CVec4<SType> CVec4Type;
    typedef std::vector<CVec4Type*> CVec4Type_Vector;

  private:

    SComplex m_cpl;

    int m_dir, m_n[3];

  public:
    
    SSVV_Calculator(const Vertex_Key &key);
    
    std::string Label() const;

    void Evaluate();

  };// end of class SSVV_Calculator

}// end of namespace METOOLS

#include "MODEL/Interaction_Models/Single_Vertex.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

using namespace METOOLS;
using namespace ATOOLS;

template <typename SType>
SSVV_Calculator<SType>::SSVV_Calculator(const Vertex_Key &key): 
  Lorentz_Calculator(key) 
{
  for (size_t i(0);i<3;++i) m_n[i]=-1;
  if (key.p_c->Flav().IntSpin()==0) {
    m_dir=0;
    for (size_t i(0);i<3;++i)
      if (key.J(i)->Flav().IntSpin()==0) {
	if (m_n[0]>=0) THROW(fatal_error,"Invalid flavours");
	m_n[0]=i;
      }
      else {
	if (m_n[1]<0) {
	  m_n[1]=i;
	}
	else {
	  if (m_n[2]>=0) THROW(fatal_error,"Invalid flavours");
	  m_n[2]=i;
	}
      }
  }
  else {
    m_dir=1;
    for (size_t i(0);i<3;++i) {
      DEBUG_VAR(i<<" "<<*key.J(i));
      if (key.J(i)->Flav().IntSpin()>0) {
	if (m_n[2]>=0) THROW(fatal_error,"Invalid flavours");
	m_n[2]=i;
      }
      else {
	if (m_n[0]<0) {
	  m_n[0]=i;
	}
	else {
	  if (m_n[1]>=0) THROW(fatal_error,"Invalid flavours");
	  m_n[1]=i;
	}
      }}
  }
  m_cpl=SComplex(p_v->Coupling(0)*p_cc->Coupling());
}

template <typename SType>
std::string SSVV_Calculator<SType>::Label() const
{
  return "SSVV["+ToString(m_cpl)+"]";
}

template <typename SType>
void SSVV_Calculator<SType>::Evaluate()
{
  p_v->SetZero();
  if (p_v->JA()->Zero()||p_v->JB()->Zero()||p_v->JE()->Zero()) return;
#ifdef DEBUG__BG
  msg_Debugging()<<*p_v->J(m_n[0])<<"(+)"<<*p_v->J(m_n[1])
		 <<"(+)"<<*p_v->J(m_n[2])<<" SSVV("<<m_dir<<")\n";
  msg_Indent();
#endif
  size_t i(0);
  typename CObject_Matrix::const_iterator jit[3];
  const CObject_Matrix &cca(p_v->JA()->J()),
    &ccb(p_v->JB()->J()), &cce(p_v->JE()->J());
  for (jit[0]=cca.begin();jit[0]!=cca.end();++jit[0]) {
    for (jit[1]=ccb.begin();jit[1]!=ccb.end();++jit[1]) {
      for (jit[2]=cce.begin();jit[2]!=cce.end();++jit[2],++i) {
  if (m_dir==0) {
    const CVec4Type_Vector *ca(jit[m_n[1]]->Get<CVec4Type>());
    const CVec4Type_Vector *cb(jit[m_n[2]]->Get<CVec4Type>());
    const CScalarType_Vector *ce(jit[m_n[0]]->Get<CScalarType>());
    for (typename CVec4Type_Vector::const_iterator 
	   ait(ca->begin());ait!=ca->end();++ait)
      for (typename CVec4Type_Vector::const_iterator 
	     bit(cb->begin());bit!=cb->end();++bit)
	for (typename CScalarType_Vector::const_iterator 
	       eit(ce->begin());eit!=ce->end();++eit)
	  if (p_cc->Evaluate(*ait,*bit,*eit)) {
#ifdef DEBUG__BG
	    msg_Debugging()<<"  a "<<**ait<<"\n";
	    msg_Debugging()<<"  b "<<**bit<<"\n";
	    msg_Debugging()<<"  e "<<**eit<<"\n";
#endif
	    CScalarType *j(CScalarType::New((**ait***bit)***eit));
	    *j*=SComplex(m_cpl);
	    j->SetH(p_v->H(i));
	    j->SetS((*ait)->S()|(*bit)->S()|(*eit)->S());
	    p_cc->AddJ(j);
	    p_v->SetZero(false);
	  }
  }
  else {
    const CScalarType_Vector *ca(jit[m_n[0]]->Get<CScalarType>());
    const CScalarType_Vector *cb(jit[m_n[1]]->Get<CScalarType>());
    const CVec4Type_Vector *ce(jit[m_n[2]]->Get<CVec4Type>());
    for (typename CScalarType_Vector::const_iterator 
	   ait(ca->begin());ait!=ca->end();++ait)
      for (typename CScalarType_Vector::const_iterator 
	     bit(cb->begin());bit!=cb->end();++bit)
	for (typename CVec4Type_Vector::const_iterator 
	       eit(ce->begin());eit!=ce->end();++eit)
	  if (p_cc->Evaluate(*ait,*bit,*eit)) {
#ifdef DEBUG__BG
	    msg_Debugging()<<"  a "<<**ait<<"\n";
	    msg_Debugging()<<"  b "<<**bit<<"\n";
	    msg_Debugging()<<"  e "<<**eit<<"\n";
#endif
	    CVec4Type *j(CVec4Type::New((**ait***bit)***eit));
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
}

namespace METOOLS {

  template class SSVV_Calculator<double>;
  template class SSVV_Calculator<long double>;

}

DECLARE_TEMPLATE_GETTER(SSVV_Getter,"VVSS",Lorentz_Calculator,Vertex_Key);

template <typename SType,char STag> Lorentz_Calculator *
SSVV_Getter<SType,STag>::operator()(const Vertex_Key &key) const
{
  return new SSVV_Calculator<SType>(key);
}

template <typename SType,char STag>
void SSVV_Getter<SType,STag>::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"SSVV vertex";
}

template class SSVV_Getter<double,'D'>;
template class SSVV_Getter<long double,'Q'>;

