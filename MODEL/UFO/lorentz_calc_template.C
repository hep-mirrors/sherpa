#include "METOOLS/Explicit/Lorentz_Calculator.H"
#include "METOOLS/Currents/C_Scalar.H"
#include "METOOLS/Currents/C_Spinor.H"
#include "METOOLS/Currents/C_Vector.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

typedef std::complex<double> complex;

namespace METOOLS {

  template <typename SType>
  class ${vertex_name}_Calculator: public Lorentz_Calculator {
  public:
    
    ${vertex_name}_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "${vertex_name}"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
${implementation}
    }

  };// end of class ${vertex_name}_Calculator

  template class ${vertex_name}_Calculator<double>;

}// end of namespace METOOLS

using namespace METOOLS;

DECLARE_GETTER(${vertex_name}_Calculator<double>,"D${vertex_name}",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,${vertex_name}_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new ${vertex_name}_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    ${vertex_name}_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"${vertex_name} vertex"; }
