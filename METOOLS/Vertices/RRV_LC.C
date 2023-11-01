#include "METOOLS/Explicit/Lorentz_Calculator.H"
#include "METOOLS/Currents/C_Vector.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "METOOLS/Currents/C_RaritaSchwinger.H"
#include "METOOLS/Main/SpinFuncs.H"

using namespace ATOOLS;

namespace METOOLS {

  template <typename SType>
  class RRV_Calculator: public Lorentz_Calculator {
  public:
    
    typedef CRaritaSchwinger<SType> CRaritaSchwingerType;
    typedef CVec4<SType> CVec4Type;

    RRV_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "RRV"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      // TODO: Bar of a necessary?
      // TODO: Welche Helizitäten für den Vierervektor?
      CRaritaSchwingerType *a(jj[1]->Get<CRaritaSchwingerType>());
      CRaritaSchwingerType *b(jj[0]->Get<CRaritaSchwingerType>());
      Gamma<SType> gamma = Gamma<SType>();
      CVec4Type *j(CVec4Type::New(0.0,0.0,0.0,0.0,0,0,0,a->S()|b->S()));
      for (size_t i(0); i<4; ++i){
        for (size_t k(0); k<16; ++k){
          for (size_t l(0); l<16; ++l){
            std::vector<std::complex<SType>> intermediate(16);
            (*j)[i] += (*a)[l] * std::complex<SType>(gamma[i][k%4][l%4]) * (*b)[l];
          }
        }
      }
      return j;
    }
  };// end of class RRV_Calculator

  template class RRV_Calculator<double>;

}// end of namespace METOOLS

using namespace METOOLS;

DECLARE_GETTER(RRV_Calculator<double>,"DRRV",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,RRV_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new RRV_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    RRV_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"RRV vertex"; }
