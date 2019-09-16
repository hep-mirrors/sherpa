#include "METOOLS/Explicit/Lorentz_Calculator.H"
#include "METOOLS/Currents/C_Scalar.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Single_Vertex.H"

using namespace ATOOLS;

namespace METOOLS {

  template <typename SType>
  class GGGG_Calculator: public Lorentz_Calculator {
  public:
    
    typedef CScalar<SType> CScalarType;

    GGGG_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "GGGG"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      const CScalarType &a(*jj[0]->Get<CScalarType>());
      const CScalarType &e(*jj[1]->Get<CScalarType>());
      const CScalarType &b(*jj[2]->Get<CScalarType>());
      CScalarType *j(CScalarType::New(a*e*b));
      j->SetS(a.S()|e.S()|b.S());
      return j;
    }

  };// end of class GGGG_Calculator

  template class GGGG_Calculator<double>;

}// end of namespace METOOLS

using namespace METOOLS;

DECLARE_GETTER(GGGG_Calculator<double>,"DGGGG",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,GGGG_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new GGGG_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    GGGG_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"GGGG vertex"; }
