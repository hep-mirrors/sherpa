#include "METOOLS/Explicit/Lorentz_Calculator.H"
#include "METOOLS/Currents/C_Vector.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Single_Vertex.H"

using namespace ATOOLS;

namespace METOOLS {

  template <typename SType>
  class VVVVA_Calculator: public Lorentz_Calculator {
  public:
    
    typedef CVec4<SType> CVec4Type;

    VVVVA_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "VVVVA"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      const CVec4Type &a(*jj[0]->Get<CVec4Type>());
      const CVec4Type &e(*jj[2]->Get<CVec4Type>());
      const CVec4Type &b(*jj[1]->Get<CVec4Type>());
      CVec4Type *j(CVec4Type::New((e*b)*a-(e*a)*b));
      j->SetS(a.S()|e.S()|b.S());
      return j;
    }

  };// end of class VVVVA_Calculator

  template class VVVVA_Calculator<double>;

  template <typename SType>
  class VVVVB_Calculator: public Lorentz_Calculator {
  public:
    
    typedef std::complex<SType> SComplex;

    typedef CVec4<SType> CVec4Type;

    VVVVB_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "VVVVB"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      const CVec4Type &a(*jj[0]->Get<CVec4Type>());
      const CVec4Type &e(*jj[1]->Get<CVec4Type>());
      const CVec4Type &b(*jj[2]->Get<CVec4Type>());
      CVec4Type *j(CVec4Type::New((e*b)*a-(e*a)*b));
      j->SetS(a.S()|e.S()|b.S());
      return j;
    }

  };// end of class VVVVB_Calculator

  template class VVVVB_Calculator<double>;

  template <typename SType>
  class VVVVC_Calculator: public Lorentz_Calculator {
  public:
    
    typedef std::complex<SType> SComplex;

    typedef CVec4<SType> CVec4Type;

    VVVVC_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "VVVVC"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      const CVec4Type &a(*jj[1]->Get<CVec4Type>());
      const CVec4Type &e(*jj[0]->Get<CVec4Type>());
      const CVec4Type &b(*jj[2]->Get<CVec4Type>());
      CVec4Type *j(CVec4Type::New((e*b)*a-(e*a)*b));
      j->SetS(a.S()|e.S()|b.S());
      return j;
    }

  };// end of class VVVVC_Calculator

  template class VVVVC_Calculator<double>;

  template <typename SType>
  class VVVV_Calculator: public Lorentz_Calculator {
  public:
    
    typedef std::complex<SType> SComplex;

    typedef CVec4<SType> CVec4Type;

    VVVV_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "VVVV"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      if (p_v->V()->id.back()%2) {
	const CVec4Type &a(*jj[0]->Get<CVec4Type>());
	const CVec4Type &e(*jj[2]->Get<CVec4Type>());
	const CVec4Type &b(*jj[1]->Get<CVec4Type>());
	CVec4Type *j(CVec4Type::New((e*b)*a+(e*a)*b-SComplex(2.0)*(a*b)*e));
	j->SetS(a.S()|e.S()|b.S());
	return j;
      }
      const CVec4Type &a(*jj[1]->Get<CVec4Type>());
      const CVec4Type &e(*jj[0]->Get<CVec4Type>());
      const CVec4Type &b(*jj[2]->Get<CVec4Type>());
      CVec4Type *j(CVec4Type::New((e*b)*a+(e*a)*b-SComplex(2.0)*(a*b)*e));
      j->SetS(a.S()|e.S()|b.S());
      return j;
    }

  };// end of class VVVV_Calculator

  template class VVVV_Calculator<double>;

}// end of namespace METOOLS

using namespace METOOLS;

DECLARE_GETTER(VVVVA_Calculator<double>,"DVVVVA",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,VVVVA_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new VVVVA_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    VVVVA_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"VVVVA vertex"; }

DECLARE_GETTER(VVVVB_Calculator<double>,"DVVVVB",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,VVVVB_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new VVVVB_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    VVVVB_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"VVVVB vertex"; }

DECLARE_GETTER(VVVVC_Calculator<double>,"DVVVVC",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,VVVVC_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new VVVVC_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    VVVVC_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"VVVVC vertex"; }

DECLARE_GETTER(VVVV_Calculator<double>,"DVVVV",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,VVVV_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new VVVV_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    VVVV_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"VVVV vertex"; }
