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
      // TODO: Helizitäten des Vierervektors setzen!!!
      CRaritaSchwingerType *a(jj[1]->Get<CRaritaSchwingerType>());
      CRaritaSchwingerType *b(jj[0]->Get<CRaritaSchwingerType>());
      Gamma<SType> gamma = Gamma<SType>();
      CVec4Type *j(CVec4Type::New(0.0,0.0,0.0,0.0,0,0,0,a->S()|b->S()));
      // Explicit expressions for psinbar gamma_mu psi'
      (*j)[0] = (*a)[0]*(*b)[2]+(*a)[1]*(*b)[3]+(*a)[2]*(*b)[0]+(*a)[3]*(*b)[1]-(*a)[4]*(*b)[6]
                -(*a)[5]*(*b)[7]-(*a)[6]*(*b)[4]-(*a)[7]*(*b)[5]-(*a)[8]*(*b)[10]-(*a)[9]*(*b)[11]-(*a)[10]*(*b)[8]
                -(*a)[11]*(*b)[9]-(*a)[12]*(*b)[14]-(*a)[13]*(*b)[15]-(*a)[14]*(*b)[12]-(*a)[15]*(*b)[13];
      (*j)[1] = (*a)[0]*(*b)[3]+(*a)[1]*(*b)[2]-(*a)[2]*(*b)[1]-(*a)[3]*(*b)[0]-(*a)[4]*(*b)[7]
                -(*a)[5]*(*b)[6]+(*a)[6]*(*b)[5]+(*a)[7]*(*b)[4]-(*a)[8]*(*b)[11]-(*a)[9]*(*b)[10]+(*a)[10]*(*b)[9]
                +(*a)[11]*(*b)[8]-(*a)[12]*(*b)[15]-(*a)[13]*(*b)[14]+(*a)[14]*(*b)[13]+(*a)[15]*(*b)[12];
      (*j)[2] = std::complex<SType>(0,1)*(std::complex<SType>(-1,0)*(*a)[0]*(*b)[3]+(*a)[1]*(*b)[2]+(*a)[2]*(*b)[1]-(*a)[3]*(*b)[0]+(*a)[4]*(*b)[7]
                                          -(*a)[5]*(*b)[6]-(*a)[6]*(*b)[5]+(*a)[7]*(*b)[4]+(*a)[8]*(*b)[11]-(*a)[9]*(*b)[10]-(*a)[10]*(*b)[9]
                                          +(*a)[11]*(*b)[8]+(*a)[12]*(*b)[15]-(*a)[13]*(*b)[14]-(*a)[14]*(*b)[13]+(*a)[15]*(*b)[12]);
      (*j)[3] =(*a)[0]*(*b)[2]-(*a)[1]*(*b)[3]-(*a)[2]*(*b)[0]+(*a)[3]*(*b)[1]-(*a)[4]*(*b)[6]
               +(*a)[5]*(*b)[7]+(*a)[6]*(*b)[4]-(*a)[7]*(*b)[5]-(*a)[8]*(*b)[10]+(*a)[9]*(*b)[11]+(*a)[10]*(*b)[8]
               -(*a)[11]*(*b)[9]-(*a)[12]*(*b)[14]+(*a)[13]*(*b)[15]+(*a)[14]*(*b)[12]-(*a)[15]*(*b)[13];
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
