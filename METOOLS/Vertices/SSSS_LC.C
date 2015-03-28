#include "METOOLS/Explicit/Lorentz_Calculator.H"
#include "METOOLS/Currents/C_Scalar.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Single_Vertex.H"

using namespace ATOOLS;

namespace METOOLS {

  template <typename SType>
  class SSSS_Calculator: public Lorentz_Calculator {
  public:
    
    typedef CScalar<SType> CScalarType;

    SSSS_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "SSSS"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      const CScalarType &a(*jj[0]->Get<CScalarType>());
      const CScalarType &e(*jj[1]->Get<CScalarType>());
      const CScalarType &b(*jj[2]->Get<CScalarType>());
      CScalarType *j(CScalarType::New(a*e*b));
      j->SetS(a.S()|e.S()|b.S());
      return j;
    }

  };// end of class SSSS_Calculator

  template class SSSS_Calculator<double>;

}// end of namespace METOOLS

using namespace METOOLS;

DECLARE_GETTER(SSSS_Calculator<double>,"DSSSS",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,SSSS_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new SSSS_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    SSSS_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"SSSS vertex"; }
