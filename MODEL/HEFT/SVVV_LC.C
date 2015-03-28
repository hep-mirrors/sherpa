#include "METOOLS/Explicit/Lorentz_Calculator.H"
#include "METOOLS/Currents/C_Vector.H"
#include "METOOLS/Currents/C_Scalar.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;

namespace METOOLS {

  template <typename SType>
  class SVVV_Calculator: public Lorentz_Calculator {
  private:

    int m_n[3];

  public:
    
    typedef CVec4<SType> CVec4Type;
    typedef CScalar<SType> CScalarType;

    SVVV_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key)
    {
      if (p_v->V()->id.back()==3) { m_n[0]=0; m_n[1]=1; m_n[2]=2; }
      if (p_v->V()->id.back()==2) { m_n[0]=1; m_n[1]=0; m_n[2]=2; }
      if (p_v->V()->id.back()==1) { m_n[0]=2; m_n[1]=0; m_n[2]=1; }
    }

    std::string Label() const { return "HVVV"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      if (p_v->V()->id.back()==0) {
	const CVec4Type &a(*jj[0]->Get<CVec4Type>());
	const CVec4Type &b(*jj[1]->Get<CVec4Type>());
	const CVec4Type &c(*jj[2]->Get<CVec4Type>()); 
	Vec4D pa(p_v->J(0)->P()), pb(p_v->J(1)->P()), pc(p_v->J(2)->P());
	CScalarType *j(CScalarType::New
		       ((a*b)*(c*(pa-pb))+
			(b*c)*(a*(pb-pc))+
			(c*a)*(b*(pc-pa))));
	j->SetS(a.S()|b.S()|c.S());
	return j;
      }
      const CVec4Type &a(*jj[m_n[1]]->Get<CVec4Type>());
      const CVec4Type &b(*jj[m_n[2]]->Get<CVec4Type>());
      const CScalarType &e(*jj[m_n[0]]->Get<CScalarType>()); 
      Vec4D pa(p_v->J(m_n[1])->P()), pb(p_v->J(m_n[2])->P());
      Vec4D pe(p_v->J(m_n[0])->P());
      CVec4Type *j(CVec4Type::New
		   (e[0]*((a*b)*CVec4Type(pa-pb)
			  +(a*ATOOLS::Vec4<SType>(pb+pb+pa+pe))*b
			  -(b*ATOOLS::Vec4<SType>(pa+pa+pb+pe))*a)));
      j->SetS(a.S()|b.S()|e.S());
      return j;
    }

  };// end of class SVVV_Calculator

  template class SVVV_Calculator<double>;

}// end of namespace METOOLS

using namespace METOOLS;

DECLARE_GETTER(SVVV_Calculator<double>,"DHVVV",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,SVVV_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new SVVV_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    SVVV_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"HVVV vertex"; }
