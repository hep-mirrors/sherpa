#include "METOOLS/Explicit/Lorentz_Calculator.H"

#include "METOOLS/Explicit/Vertex.H"
#include "METOOLS/Explicit/Dipole_Kinematics.H"
#include "METOOLS/Explicit/Dipole_Color.H"
#include "METOOLS/Explicit/Dipole_Terms.H"
#include "METOOLS/Currents/C_Vector.H"
#include "METOOLS/Currents/C_Scalar.H"
#include "ATOOLS/Phys/Spinor.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"

namespace METOOLS {

  template <typename SType>
  class Box_Calculator: public Lorentz_Calculator {
  public:

    typedef std::complex<SType> SComplex;

    typedef ATOOLS::Spinor<SType> SpinorType;

    typedef CVec4<SType> CVec4Type;
    typedef std::vector<CVec4Type*> CVec4Type_Vector;

    typedef CScalar<SType> CScalarType;
    typedef std::vector<CScalarType*> CScalarType_Vector;

  private:

    SComplex m_cpl;

    int m_mode, m_n[3];

    CVec4Type Lorentz(const CVec4Type &a,const ATOOLS::Vec4D &pa,
		      const CVec4Type &b,const ATOOLS::Vec4D &pb,
		      const CScalarType &e,const ATOOLS::Vec4D &pe);

  public:
    
    Box_Calculator(const Vertex_Key &key);
    
    std::string Label() const;

    void Evaluate();

    void SetGauge(const ATOOLS::Vec4D &k);

  };// end of class Box_Calculator

}// end of namespace METOOLS

#include "MODEL/Interaction_Models/Single_Vertex.H" 
#include "ATOOLS/Org/Message.H"

using namespace METOOLS;
using namespace ATOOLS;

template <typename SType>
Box_Calculator<SType>::Box_Calculator(const Vertex_Key &key): 
  Lorentz_Calculator(key), m_mode(0)
{
  m_cpl=SComplex(p_v->Coupling(0)*p_cc->Coupling());
  for (size_t i(0);i<3;++i)
    m_n[i]=key.p_mv->Lorentz[key.m_n]->ParticleArg(i)-1;
  if (m_n[0]>0 && m_n[1]>0 && m_n[2]>0) {
    THROW(not_implemented,"Implement me!");
  }
  else {
    if (m_n[1]<0 || m_n[2]<0) THROW(fatal_error,"Invalid rotation");
    if (m_n[1]==0) m_n[0]=(m_n[2]==1)?2:1;
    if (m_n[1]==1) m_n[0]=(m_n[2]==0)?2:0;
    if (m_n[1]==2) m_n[0]=(m_n[2]==1)?0:1;
  }
}

template <typename SType>
std::string Box_Calculator<SType>::Label() const
{
  return "Box["+ToString(m_cpl)+"]";
}

template <typename SType>
void Box_Calculator<SType>::Evaluate()
{
  p_v->SetZero();
  if (p_v->JA()->Zero()||p_v->JB()->Zero()||p_v->JE()->Zero()) return;
  if (m_mode==1) {
    THROW(not_implemented,"Implement me!");
  }
  else {
#ifdef DEBUG__BG
    msg_Debugging()<<*p_v->J(m_n[0])<<"(+)"<<*p_v->J(m_n[1])
     		   <<"(+)"<<*p_v->J(m_n[2])<<" Box("<<m_mode
		   <<"), m_cpl = "<<m_cpl<<"\n";
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
		  const CScalarType *eit((CScalarType*)*cit[m_n[0]]); 
		  const CVec4Type *ait((CVec4Type*)*cit[m_n[1]]); 
		  const CVec4Type *bit((CVec4Type*)*cit[m_n[2]]); 
#ifdef DEBUG__BG
		  msg_Debugging()<<"  a "<<*ait<<"\n";
		  msg_Debugging()<<"  b "<<*bit<<"\n";
		  msg_Debugging()<<"  e "<<*eit<<"\n";
#endif
		  CVec4Type *j(CVec4Type::New
			       (Lorentz(*ait,p_v->J(m_n[1])->P(),
					*bit,p_v->J(m_n[2])->P(),
					*eit,p_v->J(m_n[0])->P())
				*SComplex(-m_cpl)));
		  j->SetH(p_v->H(i));
		  j->SetS(ait->S()|bit->S()|eit->S());
		  p_cc->AddJ(j);
		  p_v->SetZero(false);
		}
	}
      }
    }
  }
}

template <typename SType> CVec4<SType>
Box_Calculator<SType>::Lorentz(const CVec4Type &a,const Vec4D &pa,
			       const CVec4Type &b,const Vec4D &pb,
			       const CScalarType &e,const ATOOLS::Vec4D &pe)
{
  return e[0]*((a*b)*CVec4Type(pa-pb)
	       +(a*ATOOLS::Vec4<SType>(pb+pb+pa+pe))*b
	       -(b*ATOOLS::Vec4<SType>(pa+pa+pb+pe))*a);
}

namespace METOOLS {

  template class Box_Calculator<double>;

}

DECLARE_GETTER(Box_Calculator<double>,"DBox",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,Box_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new Box_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    Box_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"Box vertex"; }
