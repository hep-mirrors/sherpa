#include "METOOLS/Explicit/Lorentz_Calculator.H"
#include "METOOLS/Currents/C_Scalar.H"
#include "METOOLS/Currents/C_Vector.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

namespace METOOLS {

  template <typename SType>
  class HVVVVA_Calculator: public Lorentz_Calculator {
  public:

    typedef CVec4<SType> CVec4Type;
    typedef CScalar<SType> CScalarType;
    
    HVVVVA_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "HVVVVA"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      if (p_v->V()->id.back()!=4) THROW(not_implemented,"Help!");
      const CVec4Type &a(*jj[0]->Get<CVec4Type>());
      const CVec4Type &e(*jj[2]->Get<CVec4Type>());
      const CVec4Type &b(*jj[1]->Get<CVec4Type>());
      const CVec4Type &d(*jj[3]->Get<CVec4Type>());
      CScalarType *j(CScalarType::New((e*b)*(a*d)-(e*a)*(b*d)));
      j->SetS(a.S()|e.S()|b.S()|d.S());
      return j;
    }

  };// end of class HVVVVA_Calculator

  template class HVVVVA_Calculator<double>;

  template <typename SType>
  class HVVVVB_Calculator: public Lorentz_Calculator {
  public:
    
    typedef CVec4<SType> CVec4Type;
    typedef CScalar<SType> CScalarType;

    HVVVVB_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "HVVVVB"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      if (p_v->V()->id.back()!=4) THROW(not_implemented,"Help!");
      const CVec4Type &a(*jj[0]->Get<CVec4Type>());
      const CVec4Type &e(*jj[1]->Get<CVec4Type>());
      const CVec4Type &b(*jj[2]->Get<CVec4Type>());
      const CVec4Type &d(*jj[3]->Get<CVec4Type>());
      CScalarType *j(CScalarType::New((e*b)*(a*d)-(e*a)*(b*d)));
      j->SetS(a.S()|e.S()|b.S()|d.S());
      return j;
    }

  };// end of class HVVVVB_Calculator

  template class HVVVVB_Calculator<double>;

  template <typename SType>
  class HVVVVC_Calculator: public Lorentz_Calculator {
  public:
    
    typedef CVec4<SType> CVec4Type;
    typedef CScalar<SType> CScalarType;

    HVVVVC_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "HVVVVC"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      if (p_v->V()->id.back()!=4) THROW(not_implemented,"Help!");
      const CVec4Type &a(*jj[1]->Get<CVec4Type>());
      const CVec4Type &e(*jj[0]->Get<CVec4Type>());
      const CVec4Type &b(*jj[2]->Get<CVec4Type>());
      const CVec4Type &d(*jj[3]->Get<CVec4Type>());
      CScalarType *j(CScalarType::New((e*b)*(a*d)-(e*a)*(b*d)));
      j->SetS(a.S()|e.S()|b.S()|d.S());
      return j;
    }

  };// end of class HVVVVC_Calculator

  template class HVVVVC_Calculator<double>;

}// end of namespace METOOLS

using namespace METOOLS;

DECLARE_GETTER(HVVVVA_Calculator<double>,"DHVVVVA",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,HVVVVA_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new HVVVVA_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    HVVVVA_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"HVVVVA vertex"; }

DECLARE_GETTER(HVVVVB_Calculator<double>,"DHVVVVB",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,HVVVVB_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new HVVVVB_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    HVVVVB_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"HVVVVB vertex"; }

DECLARE_GETTER(HVVVVC_Calculator<double>,"DHVVVVC",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,HVVVVC_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new HVVVVC_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    HVVVVC_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"HVVVVC vertex"; }
