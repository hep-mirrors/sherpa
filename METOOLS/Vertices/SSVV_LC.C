#include "METOOLS/Explicit/Lorentz_Calculator.H"
#include "METOOLS/Currents/C_Scalar.H"
#include "METOOLS/Currents/C_Spinor.H"
#include "METOOLS/Currents/C_Vector.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

namespace METOOLS {

  template <typename SType>
  class VVSS_Calculator: public Lorentz_Calculator {
  public:
    
    typedef CVec4<SType> CVec4Type;
    typedef CScalar<SType> CScalarType;

    VVSS_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "VVSS"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      if (p_v->V()->id.back()==0){
	const CVec4Type &a(*jj[0]->Get<CVec4Type>());
	const CScalarType &e(*jj[1]->Get<CScalarType>());
	const CScalarType &b(*jj[2]->Get<CScalarType>());
	CVec4Type *j(CVec4Type::New(a*(e[0]*b[0])));
	j->SetS(a.S()|e.S()|b.S());
	return j;
      }
      if (p_v->V()->id.back()==1){
	const CVec4Type &a(*jj[2]->Get<CVec4Type>());
	const CScalarType &e(*jj[0]->Get<CScalarType>());
	const CScalarType &b(*jj[1]->Get<CScalarType>());
	CVec4Type *j(CVec4Type::New(a*(e[0]*b[0])));
	j->SetS(a.S()|e.S()|b.S());
	return j;
      }
      if (p_v->V()->id.back()==2){
	const CVec4Type &a(*jj[1]->Get<CVec4Type>());
	const CVec4Type &e(*jj[2]->Get<CVec4Type>());
	const CScalarType &b(*jj[0]->Get<CScalarType>());
	CScalarType *j(CScalarType::New((a*e)*b[0]));
	j->SetS(a.S()|e.S()|b.S());
	return j;
      }
      if (p_v->V()->id.back()==3){
	const CVec4Type &a(*jj[0]->Get<CVec4Type>());
	const CVec4Type &e(*jj[1]->Get<CVec4Type>());
	const CScalarType &b(*jj[2]->Get<CScalarType>());
	CScalarType *j(CScalarType::New((a*e)*b[0]));
	j->SetS(a.S()|e.S()|b.S());
	return j;
      }

    }

  };// end of class VVSS_Calculator

  template class VVSS_Calculator<double>;

}// end of namespace METOOLS

using namespace METOOLS;

DECLARE_GETTER(VVSS_Calculator<double>,"DVVSS",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,VVSS_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new VVSS_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    VVSS_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"VVSS vertex"; }
