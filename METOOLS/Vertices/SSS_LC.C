#include "METOOLS/Explicit/Lorentz_Calculator.H"
#include "METOOLS/Currents/C_Scalar.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Single_Vertex.H"

using namespace ATOOLS;

namespace METOOLS {

  template <typename SType>
  class SSS_Calculator: public Lorentz_Calculator {
  public:
    
    typedef CScalar<SType> CScalarType;

    SSS_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "SSS"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      const CScalarType &a(*jj[0]->Get<CScalarType>());
      const CScalarType &b(*jj[1]->Get<CScalarType>());
      CScalarType *j(CScalarType::New(a*b));
      j->SetS(a.S()|b.S());
      return j;
    }

  };// end of class SSS_Calculator

  template class SSS_Calculator<double>;

}// end of namespace METOOLS

using namespace METOOLS;

DECLARE_GETTER(SSS_Calculator<double>,"DSSS",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,SSS_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new SSS_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    SSS_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"SSS vertex"; }

// ---- QPREC_BEGIN: long-double instantiation ----
DECLARE_GETTER(SSS_Calculator<long double>,"QSSS",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,SSS_Calculator<long double> >::
operator()(const Vertex_Key &key) const
{ return new SSS_Calculator<long double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    SSS_Calculator<long double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"SSS vertex"; }
// ---- QPREC_END ----
// ---- XPREC_BEGIN: double-double instantiation ----
DECLARE_GETTER(SSS_Calculator<ATOOLS::DDouble>,"XSSS",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,SSS_Calculator<ATOOLS::DDouble> >::
operator()(const Vertex_Key &key) const
{ return new SSS_Calculator<ATOOLS::DDouble>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    SSS_Calculator<ATOOLS::DDouble> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"SSS vertex"; }
// ---- XPREC_END ----
