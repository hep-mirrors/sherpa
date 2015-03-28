#include "METOOLS/Explicit/Lorentz_Calculator.H"
#include "METOOLS/Currents/C_Vector.H"
#include "METOOLS/Currents/C_Pseudo.H"
#include "METOOLS/Currents/C_Scalar.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;

namespace METOOLS {

  template <typename SType>
  class SVVP_Calculator: public Lorentz_Calculator {
  public:
    
    typedef CVec4<SType> CVec4Type;
    typedef CAsT4<SType> CAsT4Type;
    typedef CScalar<SType> CScalarType;

    SVVP_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "HVVP"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      if (p_v->V()->id.back()!=3) THROW(not_implemented,"Help!");
      const CScalarType &e(*jj[0]->Get<CScalarType>()); 
      const CVec4Type &a(*jj[1]->Get<CVec4Type>());
      const CVec4Type &b(*jj[2]->Get<CVec4Type>());
      CAsT4Type *j(CAsT4Type::New(CAsT4Type(a,b)*e[0]));
      j->SetS(a.S()|b.S()|e.S());
      return j;
    }

  };// end of class SVVP_Calculator

  template class SVVP_Calculator<double>;

}// end of namespace METOOLS

using namespace METOOLS;

DECLARE_GETTER(SVVP_Calculator<double>,"DHVVP",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,SVVP_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new SVVP_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    SVVP_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"HVVP vertex"; }
