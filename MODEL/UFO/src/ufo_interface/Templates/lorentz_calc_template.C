#include "METOOLS/Explicit/Lorentz_Calculator.H"
#include "METOOLS/Currents/C_Scalar.H"
#include "METOOLS/Currents/C_Spinor.H"
#include "METOOLS/Currents/C_Vector.H"
#include "METOOLS/Explicit/Vertex.H"
#include "METOOLS/Explicit/Form_Factor.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/MyComplex.H"


namespace METOOLS {

  template <typename SType>
  class ${vertex_name}_Calculator: public Lorentz_Calculator {
  public:

    typedef std::complex<SType> SComplex;

    const SComplex I = SComplex(0.0,1.0);
    
    ${vertex_name}_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {
${form_factor_impl}
    }

    std::string Label() const { return "${vertex_name}"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
${implementation}
      THROW(fatal_error, "Internal error in Lorentz calculator");
      return NULL;
    }

${form_factor_decl}
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
