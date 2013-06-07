#include "METOOLS/Explicit/Lorentz_Calculator.H"

#include "METOOLS/Explicit/Vertex.H"
#include "METOOLS/Currents/C_Vector.H"
#include "ATOOLS/Org/MyStrStream.H"

namespace METOOLS {

  template <typename SType>
  class VVVV_Calculator: public Lorentz_Calculator {
  public:

    typedef std::complex<SType> SComplex;

    typedef CVec4<SType> CVec4Type;
    typedef std::vector<CVec4Type*> CVec4Type_Vector;

  private:

    SComplex m_cpl;

    size_t m_n[3];

  public:
    
    VVVV_Calculator(const Vertex_Key &key);
    
    std::string Label() const;

    void Evaluate();

  };// end of class VVVV_Calculator

}// end of namespace METOOLS

#include "MODEL/Interaction_Models/Single_Vertex.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

using namespace METOOLS;
using namespace ATOOLS;

template <typename SType>
VVVV_Calculator<SType>::VVVV_Calculator(const Vertex_Key &key): 
  Lorentz_Calculator(key) 
{
  if (key.p_mv->Lorentz[key.m_n]->
      Type().find("Gluon4")!=std::string::npos) {
    for (size_t i(1);i<4;++i)
      m_n[i-1]=key.p_mv->Lorentz[key.m_n]->ParticleArg(i)-1;
    std::swap<size_t>(m_n[1],m_n[2]);
  }
  else {
    for (size_t i(0);i<3;++i) m_n[i]=i;
    if (key.p_a->Flav()==key.p_c->Flav().Bar())
      std::swap<size_t>(m_n[0],m_n[2]);
    else if (key.p_b->Flav()==key.p_c->Flav().Bar())
      std::swap<size_t>(m_n[1],m_n[2]);
  }
  m_cpl=SComplex(p_v->Coupling(0)*p_cc->Coupling());
}

template <typename SType>
std::string VVVV_Calculator<SType>::Label() const
{
  return "VVVV["+ToString(m_cpl)+"]";
}

template <typename SType>
void VVVV_Calculator<SType>::Evaluate()
{
  p_v->SetZero();
  if (p_v->JA()->Zero()||p_v->JB()->Zero()||p_v->JE()->Zero()) return;
#ifdef DEBUG__BG
  msg_Debugging()<<*p_v->J(m_n[0])<<"(+)"<<*p_v->J(m_n[1])
		 <<"(+)"<<*p_v->J(m_n[2])<<" VVVV\n";
  msg_Indent();
#endif
  size_t i(0);
  typename CObject_Matrix::const_iterator jit[3];
  const CObject_Matrix &cca(p_v->JA()->J()),
    &ccb(p_v->JB()->J()), &cce(p_v->JE()->J());
  for (jit[0]=cca.begin();jit[0]!=cca.end();++jit[0]) {
    for (jit[1]=ccb.begin();jit[1]!=ccb.end();++jit[1]) {
      for (jit[2]=cce.begin();jit[2]!=cce.end();++jit[2],++i) {
	const CVec4Type_Vector *ca(jit[m_n[0]]->template Get<CVec4Type>());
	const CVec4Type_Vector *cb(jit[m_n[1]]->template Get<CVec4Type>());
	const CVec4Type_Vector *ce(jit[m_n[2]]->template Get<CVec4Type>());
	for (typename CVec4Type_Vector::const_iterator 
	       eit(ce->begin());eit!=ce->end();++eit)
	  for (typename CVec4Type_Vector::const_iterator 
		 bit(cb->begin());bit!=cb->end();++bit)
	    for (typename CVec4Type_Vector::const_iterator 
		   ait(ca->begin());ait!=ca->end();++ait)
	      if (p_cc->Evaluate(*ait,*bit,*eit)) {
#ifdef DEBUG__BG
		msg_Debugging()<<"  a "<<**ait<<"\n";
		msg_Debugging()<<"  b "<<**bit<<"\n";
		msg_Debugging()<<"  e "<<**eit<<"\n";
#endif
		CVec4Type *j(CVec4Type::New
			     (SComplex(2.0)*(**ait***bit)***eit
			      -(**eit***bit)***ait-(**eit***ait)***bit));
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

  template class VVVV_Calculator<double>;
  template class VVVV_Calculator<long double>;

}

DECLARE_TEMPLATE_GETTER(Gauge4_Getter,"Gauge4",Lorentz_Calculator,Vertex_Key);

template <typename SType,char STag> Lorentz_Calculator *
Gauge4_Getter<SType,STag>::operator()(const Vertex_Key &key) const
{
  return new VVVV_Calculator<SType>(key);
}

template <typename SType,char STag>
void Gauge4_Getter<SType,STag>::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"VVVV vertex";
}

template class Gauge4_Getter<double,'D'>;
template class Gauge4_Getter<long double,'Q'>;

DECLARE_TEMPLATE_GETTER(Gluon4_Getter,"Gluon4",Lorentz_Calculator,Vertex_Key);

template <typename SType,char STag> Lorentz_Calculator *
Gluon4_Getter<SType,STag>::operator()(const Vertex_Key &key) const
{
  return new VVVV_Calculator<SType>(key);
}

template <typename SType,char STag>
void Gluon4_Getter<SType,STag>::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"VVVV vertex";
}

template class Gluon4_Getter<double,'D'>;
template class Gluon4_Getter<long double,'Q'>;

