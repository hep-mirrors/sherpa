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

    int m_n[3], m_mode;


    inline CVec4Type *Lorentz
    (const CVec4Type &a,const CVec4Type &e,const CVec4Type&b)
    {
      if (m_mode==1) return CVec4Type::New((e*b)*a-(e*a)*b);
      return CVec4Type::New(SComplex(2.0)*(a*b)*e-(e*b)*a-(e*a)*b);
    }

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
  m_mode=0;
  if (key.p_mv->Lorentz[key.m_n]->
      Type().find("Gluon4")!=std::string::npos) {
    m_mode=1;
    for (size_t i(1);i<4;++i)
      m_n[i-1]=key.p_mv->Lorentz[key.m_n]->ParticleArg(i)-1;
  }
  else {
    int n[4];
    for (size_t i(0);i<4;++i)
      n[i]=key.p_mv->Lorentz[key.m_n]->ParticleArg(i)-1;
    if (n[0]<0) { m_n[1]=n[1]; m_n[0]=n[2]; m_n[2]=n[3]; }
    if (n[1]<0) { m_n[1]=n[0]; m_n[0]=n[2]; m_n[2]=n[3]; }
    if (n[2]<0) { m_n[1]=n[3]; m_n[0]=n[0]; m_n[2]=n[1]; }
    if (n[3]<0) { m_n[1]=n[2]; m_n[0]=n[0]; m_n[2]=n[1]; }
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
		 <<"(+)"<<*p_v->J(m_n[2])<<" VVVV("<<m_mode<<")\n";
  msg_Indent();
#endif
  size_t i(0);
  const CObject_Matrix &cca(p_v->JA()->J()),
    &ccb(p_v->JB()->J()), &cce(p_v->JE()->J());
  for (typename CObject_Matrix::const_iterator 
	 jait(cca.begin());jait!=cca.end();++jait) {
    for (typename CObject_Matrix::const_iterator 
	   jbit(ccb.begin());jbit!=ccb.end();++jbit) {
      for (typename CObject_Matrix::const_iterator 
	     jeit(cce.begin());jeit!=cce.end();++jeit,++i) {
	typename CObject_Vector::const_iterator cit[3];
	for (cit[2]=jeit->begin();cit[2]!=jeit->end();++cit[2])
	  for (cit[1]=jbit->begin();cit[1]!=jbit->end();++cit[1])
	    for (cit[0]=jait->begin();cit[0]!=jait->end();++cit[0])
	      if (p_cc->Evaluate(*cit[0],*cit[1],*cit[2])) {
		const CVec4Type *ait((CVec4Type*)*cit[m_n[0]]); 
		const CVec4Type *bit((CVec4Type*)*cit[m_n[1]]); 
		const CVec4Type *eit((CVec4Type*)*cit[m_n[2]]); 
#ifdef DEBUG__BG
		msg_Debugging()<<"  a "<<*ait<<"\n";
		msg_Debugging()<<"  b "<<*bit<<"\n";
		msg_Debugging()<<"  e "<<*eit<<"\n";
#endif
		CVec4Type *j(Lorentz(*ait,*bit,*eit));
		*j*=SComplex(m_cpl);
		j->SetH(p_v->H(i));
		j->SetS(ait->S()|bit->S()|eit->S());
		p_cc->AddJ(j);
		p_v->SetZero(false);
	      }
      }
    }
  }
}

namespace METOOLS {

  template class VVVV_Calculator<double>;

  template <typename SType>
  class Gluon4_Calculator: public VVVV_Calculator<SType> {};

  template class Gluon4_Calculator<double>;

}

DECLARE_GETTER(VVVV_Calculator<double>,"DGauge4",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,VVVV_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new VVVV_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    VVVV_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"VVVV vertex"; }

DECLARE_GETTER(Gluon4_Calculator<double>,"DGluon4",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,Gluon4_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new VVVV_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    Gluon4_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"VVVV vertex"; }
