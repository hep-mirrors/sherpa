#include "METOOLS/Explicit/Lorentz_Calculator.H"
#include "METOOLS/Currents/C_Vector.H"
#include "METOOLS/Currents/C_Scalar.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Phys/Spinor.H"
#include "MODEL/SMH/gggH.h"

using namespace ATOOLS;

namespace METOOLS {

  template <typename SType>
  class HGGG_Calculator: public Lorentz_Calculator {
  private:

    int m_n[3];
    double m_mh;
    std::vector<double> m_mq;

  public:
    
    typedef CVec4<SType> CVec4Type;
    typedef CScalar<SType> CScalarType;

    HGGG_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key)
    {
      m_mh=Flavour(kf_h0).Mass();
      for (size_t i(0);i<6;++i)
	if (Flavour((kf_code)i).Yuk())
	  m_mq.push_back(Flavour((kf_code)i).Yuk());
      if (p_v->V()->id.back()==3) { m_n[0]=0; m_n[1]=1; m_n[2]=2; }
      if (p_v->V()->id.back()==2) { m_n[0]=1; m_n[1]=0; m_n[2]=2; }
      if (p_v->V()->id.back()==1) { m_n[0]=2; m_n[1]=0; m_n[2]=1; }
    }

    std::string Label() const { return "HGGG"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      // arXiv:1712.06549 [hep-ph], Eqs. (2.9) & (2.10)
      if (p_v->V()->id.back()==0) {
	const CVec4Type &a(*jj[0]->template Get<CVec4Type>());
	const CVec4Type &b(*jj[1]->template Get<CVec4Type>());
	const CVec4Type &c(*jj[2]->template Get<CVec4Type>()); 
	Vec4D pa(p_v->J(0)->P()), pb(p_v->J(1)->P()), pc(p_v->J(2)->P());
	double sab((pa+pb).Abs2()), sbc((pb+pc).Abs2()), sac((pa+pc).Abs2());
	Complex F[4]={0.,0.,0.,0.};
#ifdef USING__HEFT
	F[0]=-1.-sac/sab-(sab+sac)/sbc;
	F[1]=-1.-sbc/sac-(sac+sbc)/sab;
	F[2]=-1.-sab/sbc-(sab+sbc)/sac;
	F[3]=2./sab+2./sac+2./sbc;
	for (int i(0);i<4;++i) F[i]*=2./3.;
#else
	for (size_t i(0);i<m_mq.size();++i)
	  Fggg(m_mq[i],m_mh,sab,sac,sbc,F);
#endif
	CScalarType *j(CScalarType::New
		       ((a*b)*(c*pb)*F[0]+
			(b*c)*(a*pc)*F[1]+
			(c*a)*(b*pa)*F[2]+
			(a*pc)*(b*pa)*(c*pb)*F[3]));
	j->SetS(a.S()|b.S()|c.S());
	return j;
      }
      const CVec4Type &a(*jj[m_n[1]]->template Get<CVec4Type>());
      const CVec4Type &b(*jj[m_n[2]]->template Get<CVec4Type>());
      const CScalarType &e(*jj[m_n[0]]->template Get<CScalarType>()); 
      Vec4D pa(p_v->J(m_n[1])->P()), pb(p_v->J(m_n[2])->P());
      Vec4D pe(p_v->J(m_n[0])->P()), pc(-pa-pb-pe);
      double sab((pa+pb).Abs2()), sbc((pb+pc).Abs2()), sac((pa+pc).Abs2());
      Complex F[4]={0.,0.,0.,0.};
#ifdef USING__HEFT
      F[0]=-1.-sac/sab-(sab+sac)/sbc;
      F[1]=-1.-sbc/sac-(sac+sbc)/sab;
      F[2]=-1.-sab/sbc-(sab+sbc)/sac;
      F[3]=2./sab+2./sac+2./sbc;
      for (int i(0);i<4;++i) F[i]*=2./3.;
#else
      for (size_t i(0);i<m_mq.size();++i)
	Fggg(m_mq[i],m_mh,sab,sac,sbc,F);
#endif
      CVec4Type *j(CVec4Type::New
		   (e[0]*
		    ((a*b)*CVec4Type(pb)*F[0]
		     +(b*ATOOLS::Vec4<SType>(pa))*a*F[1]
		     +(a*ATOOLS::Vec4<SType>(pc))*b*F[2]
		     +(a*ATOOLS::Vec4<SType>(pc))
		     *(b*ATOOLS::Vec4<SType>(pa))*CVec4Type(pb)*F[3])));
      j->SetS(a.S()|b.S()|e.S());
      return j;
    }

  };// end of class HGGG_Calculator

  template class HGGG_Calculator<double>;

}// end of namespace METOOLS

using namespace METOOLS;

DECLARE_GETTER(HGGG_Calculator<double>,"DHGGG",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,HGGG_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new HGGG_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    HGGG_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"HGGG vertex"; }
