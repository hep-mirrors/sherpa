#include "METOOLS/Explicit/Lorentz_Calculator.H"
#include "METOOLS/Currents/C_Vector.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Single_Vertex.H"

using namespace ATOOLS;

namespace METOOLS {

  template <typename SType>
  class VVV_Calculator: public Lorentz_Calculator {
  public:
    
    typedef CVec4<SType> CVec4Type;

    VVV_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "VVV"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      const CVec4Type &ja(*jj[0]->Get<CVec4Type>());
      const CVec4Type &jb(*jj[1]->Get<CVec4Type>());
      const ATOOLS::Vec4D &pa(p_v->J(0)->P()), &pb(p_v->J(1)->P());
      CVec4Type *j(CVec4Type::New((ja*jb)*CVec4Type(pa-pb)
				  +(ja*Vec4<SType>(pb+pb+pa))*jb
				  -(jb*Vec4<SType>(pa+pa+pb))*ja));
      j->SetS(ja.S()|jb.S());
      return j;
    }

  };// end of class VVV_Calculator

  template class VVV_Calculator<double>;

}// end of namespace METOOLS

using namespace METOOLS;

DECLARE_GETTER(VVV_Calculator<double>,"DVVV",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,VVV_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new VVV_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    VVV_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"VVV vertex"; }
