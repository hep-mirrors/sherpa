#include "METOOLS/Explicit/Lorentz_Calculator.H"
#include "METOOLS/Currents/C_Scalar.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Single_Vertex.H"

using namespace ATOOLS;

namespace METOOLS {

  template <typename SType>
  class GGG_Calculator: public Lorentz_Calculator {
  public:
    
    typedef CScalar<SType> CScalarType;

    GGG_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "GGG"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      const CScalarType &a(*jj[0]->Get<CScalarType>());
      const CScalarType &b(*jj[1]->Get<CScalarType>());
      CScalarType *j(CScalarType::New(a*b));
      j->SetS(a.S()|b.S());
      return j;
    }

  };// end of class GGG_Calculator

  template class GGG_Calculator<double>;

}// end of namespace METOOLS

using namespace METOOLS;

DECLARE_GETTER(GGG_Calculator<double>,"DGGG",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,GGG_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new GGG_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    GGG_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"GGG vertex"; }
