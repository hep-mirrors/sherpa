#include "METOOLS/Explicit/Lorentz_Calculator.H"
#include "METOOLS/Currents/C_Vector.H"
#include "METOOLS/Currents/C_Pseudo.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Single_Vertex.H"

using namespace ATOOLS;

namespace METOOLS {

  template <typename SType>
  class VVP_Calculator: public Lorentz_Calculator {
  public:
    
    typedef std::complex<SType> SComplex;

    typedef CVec4<SType> CVec4Type;
    typedef CAsT4<SType> CAsT4Type;

    VVP_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "VVP"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      if (p_v->V()->id.back()==0){
	const CVec4Type &a(*jj[0]->Get<CVec4Type>());
	const CAsT4Type &b(*jj[1]->Get<CAsT4Type>());
	CVec4Type *j(CVec4Type::New(a*b));
	j->SetS(a.S()|b.S());
	return j;
      }
      if (p_v->V()->id.back()==1){
	const CAsT4Type &a(*jj[0]->Get<CAsT4Type>());
	const CVec4Type &b(*jj[1]->Get<CVec4Type>());
	CVec4Type *j(CVec4Type::New(-(b*a)));
	j->SetS(a.S()|b.S());
	return j;
      }
      const CVec4Type &a(*jj[0]->Get<CVec4Type>());
      const CVec4Type &b(*jj[1]->Get<CVec4Type>());
      CAsT4Type *j(CAsT4Type::New(CAsT4Type(a,b)));
      j->SetS(a.S()|b.S());
      return j;
    }

  };// end of class VVP_Calculator

  template class VVP_Calculator<double>;

}// end of namespace METOOLS

using namespace METOOLS;

DECLARE_GETTER(VVP_Calculator<double>,"DVVP",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,VVP_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new VVP_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    VVP_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"VVP vertex"; }
