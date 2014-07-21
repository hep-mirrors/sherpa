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

    int m_dir, m_v[2], m_s[2];

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
  m_cpl=SComplex(p_v->Coupling(0)*p_cc->Coupling());
  m_v[0]=key.p_mv->Lorentz[key.m_n]->ParticleArg(0)-1;
  m_v[1]=key.p_mv->Lorentz[key.m_n]->ParticleArg(1)-1;
  if (m_v[0]<0) std::swap<int>(m_v[0],m_v[1]);
  m_s[0]=m_s[1]=3;
  for (int i(-1);i<3;++i)
    if (i!=m_v[0] && i!=m_v[1]) m_s[m_s[0]<3?1:0]=i;
  if (m_s[0]<0) std::swap<int>(m_s[0],m_s[1]);
  if (m_v[0]<0 || m_s[0]<0) THROW(fatal_error,"Invalid call");
  m_dir=(m_s[0]>=0)&&(m_s[1]>=0);
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
  msg_Debugging()<<*p_v->J(0)<<"(+)"<<*p_v->J(1)
		 <<"(+)"<<*p_v->J(2)<<" SSVV("<<m_dir<<")\n";
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
		if (m_dir==0) {
		  const CVec4Type *ait((CVec4Type*)*cit[m_v[0]]); 
		  const CVec4Type *bit((CVec4Type*)*cit[m_v[1]]); 
		  const CScalarType *eit((CScalarType*)*cit[m_s[0]]); 
#ifdef DEBUG__BG
		  msg_Debugging()<<"  a "<<*ait<<"\n";
		  msg_Debugging()<<"  b "<<*bit<<"\n";
		  msg_Debugging()<<"  e "<<*eit<<"\n";
#endif
		  CScalarType *j(CScalarType::New((*ait**bit)**eit));
		  *j*=SComplex(m_cpl);
		  j->SetH(p_v->H(i));
		  j->SetS(ait->S()|bit->S()|eit->S());
		  p_cc->AddJ(j);
		}
		else {
		  const CScalarType *ait((CScalarType*)*cit[m_s[0]]); 
		  const CScalarType *bit((CScalarType*)*cit[m_s[1]]); 
		  const CVec4Type *eit((CVec4Type*)*cit[m_v[0]]); 
#ifdef DEBUG__BG
		  msg_Debugging()<<"  a "<<*ait<<"\n";
		  msg_Debugging()<<"  b "<<*bit<<"\n";
		  msg_Debugging()<<"  e "<<*eit<<"\n";
#endif
		  CVec4Type *j(CVec4Type::New((*ait**bit)**eit));
		  *j*=SComplex(m_cpl);
		  j->SetH(p_v->H(i));
		  j->SetS(ait->S()|bit->S()|eit->S());
		  p_cc->AddJ(j);
		}
		p_v->SetZero(false);
	      }
      }
    }
  }
}

namespace METOOLS {

  template class SSVV_Calculator<double>;

}

DECLARE_GETTER(SSVV_Calculator<double>,"DVVSS",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,SSVV_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new SSVV_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    SSVV_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"SSVV vertex"; }
