#include "METOOLS/Explicit/Lorentz_Calculator.H"
#include "METOOLS/Currents/C_Vector.H"
#include "METOOLS/Currents/C_Scalar.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Single_Vertex.H"

using namespace ATOOLS;

namespace METOOLS {

  template <typename SType>
  class VVS_Calculator: public Lorentz_Calculator {
  public:
    
    typedef CVec4<SType> CVec4Type;
    typedef CScalar<SType> CScalarType;

    VVS_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "VVS"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      if (p_v->V()->id.back()==2) {
	const CVec4Type &a(*jj[0]->Get<CVec4Type>());
	const CVec4Type &b(*jj[1]->Get<CVec4Type>());
	CScalarType *j(CScalarType::New(a*b));
	j->SetS(a.S()|b.S());
	return j;
      }
      if (p_v->V()->id.back()==1) {
	const CScalarType &a(*jj[0]->Get<CScalarType>());
	const CVec4Type &b(*jj[1]->Get<CVec4Type>());
	CVec4Type *j(CVec4Type::New(b*a[0]));
	j->SetS(a.S()|b.S());
	return j;
      }
      if (p_v->V()->id.back()==0) {
	const CScalarType &a(*jj[1]->Get<CScalarType>());
	const CVec4Type &b(*jj[0]->Get<CVec4Type>());
	CVec4Type *j(CVec4Type::New(b*a[0]));
	j->SetS(a.S()|b.S());
	return j;
      }
      return NULL;
    }

  };// end of class VVS_Calculator

  template class VVS_Calculator<double>;

}// end of namespace METOOLS

using namespace METOOLS;

DECLARE_GETTER(VVS_Calculator<double>,"DVVS",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,VVS_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new VVS_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    VVS_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"VVS vertex"; }
