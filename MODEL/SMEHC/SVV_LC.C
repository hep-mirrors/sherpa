#include "METOOLS/Explicit/Lorentz_Calculator.H"
#include "METOOLS/Currents/C_Vector.H"
#include "METOOLS/Currents/C_Scalar.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Single_Vertex.H"

using namespace ATOOLS;

namespace METOOLS {

  template <typename SType>
  class HVV_Calculator: public Lorentz_Calculator {
  public:
    
    typedef CVec4<SType> CVec4Type;
    typedef CScalar<SType> CScalarType;

    HVV_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "HVV"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      if (p_v->V()->id.back()==0) {
	const CVec4Type &a(*jj[0]->Get<CVec4Type>());
	const CVec4Type &b(*jj[1]->Get<CVec4Type>());
	CVec4Type pa(p_v->J(0)->P()), pb(p_v->J(1)->P());
	CScalarType *j(CScalarType::New(-(a*b)*(pa*pb)+(a*pb)*(b*pa)));
	j->SetS(a.S()|b.S());
	return j;
      }
      if (p_v->V()->id.back()==2) {
	const CScalarType &a(*jj[0]->Get<CScalarType>());
	const CVec4Type &b(*jj[1]->Get<CVec4Type>());
	CVec4Type pa(p_v->J(0)->P()), pb(p_v->J(1)->P());
	CVec4Type s(b*(pb*(pa+pb))-(b*(pa+pb))*pb);
	CVec4Type *j(CVec4Type::New(s*a[0]));
	j->SetS(a.S()|b.S());
	return j;
      }
      if (p_v->V()->id.back()==1) {
	const CScalarType &a(*jj[1]->Get<CScalarType>());
	const CVec4Type &b(*jj[0]->Get<CVec4Type>());
	CVec4Type pa(p_v->J(1)->P()), pb(p_v->J(0)->P());
	CVec4Type s(b*(pb*(pa+pb))-(b*(pa+pb))*pb);
	CVec4Type *j(CVec4Type::New(s*a[0]));
	j->SetS(a.S()|b.S());
	return j;
      }
      return NULL;
    }

  };// end of class HVV_Calculator

  template class HVV_Calculator<double>;

}// end of namespace METOOLS

using namespace METOOLS;

DECLARE_GETTER(HVV_Calculator<double>,"DHVV",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,HVV_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new HVV_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    HVV_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"HVV vertex"; }
