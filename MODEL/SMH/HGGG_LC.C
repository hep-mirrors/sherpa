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
	for (size_t i(0);i<m_mq.size();++i)
	  Fggg(m_mq[i],m_mh,sab,sac,sbc,F);
	CScalarType *j(CScalarType::New
		       ((a*b)*(c*pb)*F[0]+
			(b*c)*(a*pc)*F[1]+
			(c*a)*(b*pa)*F[2]+
			(a*pc)*(b*pa)*(c*pb)*F[3]));
	j->SetS(a.S()|b.S()|c.S());
	return j;
      }
      if (0) {// HEFT in Comix notation
	const CVec4Type &a(*jj[m_n[1]]->template Get<CVec4Type>());
	const CVec4Type &b(*jj[m_n[2]]->template Get<CVec4Type>());
	const CScalarType &e(*jj[m_n[0]]->template Get<CScalarType>());
	CVec4Type pa(p_v->J(m_n[1])->P()), pb(p_v->J(m_n[2])->P());
	CVec4Type pe(p_v->J(m_n[0])->P()), pc(-pa-pb-pe);
	CVec4Type pab(pa+pb), pae(pa+pe), pbe(pb+pe);
	Complex sab(pab.Abs2()), sae(pae.Abs2()), sbe(pbe.Abs2());
	CVec4Type ab(((a*b)*(pa-pb)+(a*(pb+pab))*b-(b*(pab+pa))*a)/sab);
	CVec4Type ae((e[0]*(a*(pa*pae)-(a*pae)*pa))/-sae);
	CVec4Type be((e[0]*(b*(pb*pbe)-(b*pbe)*pb))/-sbe);
	CVec4Type *abe(CVec4Type::New
		       ((ae*b)*(pae-pb)+(ae*(pb-pc))*b+(b*(pc-pae))*ae
			+(a*be)*(pa-pbe)+(a*(pbe-pc))*be+(be*(pc-pa))*a
			+e[0]*(ab*(pab*pc)-(ab*pc)*pab)
			+e[0]*((a*b)*(pa-pb)+(a*(pb-pc))*b+(b*(pc-pa))*a)));
	abe->SetS(a.S()|b.S()|e.S());
	return abe;
      }
      if (0) {// HEFT in form factor notation
	const CVec4Type &a(*jj[m_n[1]]->template Get<CVec4Type>());
	const CVec4Type &b(*jj[m_n[2]]->template Get<CVec4Type>());
	const CScalarType &e(*jj[m_n[0]]->template Get<CScalarType>());
	Vec4D pa(p_v->J(m_n[1])->P()), pb(p_v->J(m_n[2])->P());
	Vec4D pe(p_v->J(m_n[0])->P()), pc(-pa-pb-pe);
	Vec4D pab(pa+pb), pbc(pb+pc), pac(pa+pc);
	Complex sab(pab.Abs2()), sac(pac.Abs2()), sbc(pbc.Abs2());
	Complex F1(-1.-sac/sab-(sab+sac)/sbc);
	Complex F2(-1.-sbc/sac-(sac+sbc)/sab);
	Complex F3(-1.-sab/sbc-(sab+sbc)/sac);
	Complex F4(2./sab+2./sac+2./sbc);
	Spinor<double> kp1(1,pa), km1(-1,pa);
	Spinor<double> kp2(1,pb), km2(-1,pb);
	Spinor<double> kp3(1,pc), km3(-1,pc);
	Complex fac=Complex(0.,1.)*p_v->V()->Coupling(0)/sqrt(2.);
	Complex Appp(fac*sqr(m_mh)/(sqrt(2.)*(kp1*kp2)*(kp2*kp3)*(kp3*kp1))*
		      sab*sbc/sqr(m_mh)*(F1+sac/sbc*F2+sac/sab*F3+sac/2.*F4));
	Complex Appm1(fac*pow(km1*km2,3)/(sqrt(2.)*(km1*km3)*(km2*km3)*sqr(m_mh))*
		      sqr(m_mh)*sbc/sab*(F1+sac/2.*F4));
	Complex Appm2(fac*pow(km2*km3,3)/(sqrt(2.)*(km2*km1)*(km3*km1)*sqr(m_mh))*
		      sqr(m_mh)*sac/sbc*(F3+sab/2.*F4));
	Complex Appm3(fac*pow(km3*km1,3)/(sqrt(2.)*(km3*km2)*(km1*km2)*sqr(m_mh))*
		      sqr(m_mh)*sab/sac*(F2+sbc/2.*F4));
	DEBUG_VAR(Appp<<" "<<std::abs(Appp));
	DEBUG_VAR(Appm1<<" "<<std::abs(Appm1));
	DEBUG_VAR(Appm2<<" "<<std::abs(Appm2));
	DEBUG_VAR(Appm3<<" "<<std::abs(Appm3));
	CVec4Type *j(CVec4Type::New
		     (e[0]*((a*b)*CVec4Type(pb)*F1
			    +(b*ATOOLS::Vec4<SType>(pa))*a*F2
			    +(a*ATOOLS::Vec4<SType>(pc))*b*F3
			    +(a*ATOOLS::Vec4<SType>(pc))
			    *(b*ATOOLS::Vec4<SType>(pa))*CVec4Type(pb)*F4)));
	j->SetS(a.S()|b.S()|e.S());
	return j;
      }
      const CVec4Type &a(*jj[m_n[1]]->template Get<CVec4Type>());
      const CVec4Type &b(*jj[m_n[2]]->template Get<CVec4Type>());
      const CScalarType &e(*jj[m_n[0]]->template Get<CScalarType>()); 
      Vec4D pa(p_v->J(m_n[1])->P()), pb(p_v->J(m_n[2])->P());
      Vec4D pe(p_v->J(m_n[0])->P()), pc(-pa-pb-pe);
      double sab((pa+pb).Abs2()), sbc((pb+pc).Abs2()), sac((pa+pc).Abs2());
      Complex F[4]={0.,0.,0.,0.};
      for (size_t i(0);i<m_mq.size();++i) {
	Fggg(m_mq[i],m_mh,sab,sac,sbc,F);
	/*
	Complex F1(F[0]), F2(F[1]), F3(F[2]),F4(F[3]);
	// std::cout<<"Oppp(s,t,u) "<<sab*sbc/sqr(m_mh)*(F1+sac/sbc*F2+sac/sab*F3+sac/2.*F4)
	// 	 <<" "<<sab*sbc/sqr(m_mh)*(F[0]+sac/sbc*F[1]+sac/sab*F[2]+sac/2.*F[3])/sqr(m_mq[i])<<std::endl;
	// std::cout<<"Oppm(s,t,u) "<<sqr(m_mh)*sbc/sab*(F1+sac/2.*F4)
	// 	 <<" "<<sqr(m_mh)*sbc/sab*(F[0]+sac/2.*F[3])/sqr(m_mq[i])<<std::endl;
	// std::cout<<"Oppm(t,u,s) "<<sqr(m_mh)*sab/sac*(F2+sbc/2.*F4)
	// 	 <<" "<<sqr(m_mh)*sab/sac*(F[1]+sbc/2.*F[3])/sqr(m_mq[i])<<std::endl;
	// std::cout<<"Oppm(u,s,t) "<<sqr(m_mh)*sac/sbc*(F3+sab/2.*F4)
	// 	 <<" "<<sqr(m_mh)*sac/sbc*(F[2]+sab/2.*F[3])/sqr(m_mq[i])<<std::endl;
	Spinor<double> kp1(1,pa), km1(-1,pa);
	Spinor<double> kp2(1,pb), km2(-1,pb);
	Spinor<double> kp3(1,pc), km3(-1,pc);
	Complex fac=Complex(0.,1.)*p_v->V()->Coupling(0)/sqrt(2.);
	Complex Appp(fac*sqr(m_mh)/(sqrt(2.)*(kp1*kp2)*(kp2*kp3)*(kp3*kp1))*
		      sab*sbc/sqr(m_mh)*(F1+sac/sbc*F2+sac/sab*F3+sac/2.*F4));
	Complex Appm1(fac*pow(km1*km2,3)/(sqrt(2.)*(km1*km3)*(km2*km3)*sqr(m_mh))*
		      sqr(m_mh)*sbc/sab*(F1+sac/2.*F4));
	Complex Appm2(fac*pow(km2*km3,3)/(sqrt(2.)*(km2*km1)*(km3*km1)*sqr(m_mh))*
		      sqr(m_mh)*sac/sbc*(F3+sab/2.*F4));
	Complex Appm3(fac*pow(km3*km1,3)/(sqrt(2.)*(km3*km2)*(km1*km2)*sqr(m_mh))*
		      sqr(m_mh)*sab/sac*(F2+sbc/2.*F4));
	DEBUG_VAR(Appp<<" "<<std::abs(Appp));
	DEBUG_VAR(Appm1<<" "<<std::abs(Appm1));
	DEBUG_VAR(Appm2<<" "<<std::abs(Appm2));
	DEBUG_VAR(Appm3<<" "<<std::abs(Appm3));
	*/
      }
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
