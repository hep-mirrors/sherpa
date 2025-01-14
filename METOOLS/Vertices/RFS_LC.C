#include "METOOLS/Explicit/Lorentz_Calculator.H"
#include "METOOLS/Currents/C_Scalar.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "METOOLS/Currents/C_RaritaSchwinger.H"
#include "METOOLS/Main/SpinFuncs.H"

using namespace ATOOLS;

namespace METOOLS {

  template <typename SType>
  class RFS_Calculator: public Lorentz_Calculator {
  public:
    
    typedef CRaritaSchwinger<SType> CRaritaSchwingerType;
    typedef CSpinor<SType> CSpinorType;
    typedef CScalar<SType> CScalarType;

    RFS_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "RFS"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      if (p_v->V()->id.back()==2) {
        CRaritaSchwingerType *j0(jj[1]->Get<CRaritaSchwingerType>());
        CSpinorType *j1(jj[0]->Get<CSpinorType>());
        DEBUG_VAR(*j0);
        DEBUG_VAR(*j1);
        const ATOOLS::Vec4<SType> &q(p_v->J(1)->P() + p_v->J(0)->P()); // Scalar four momentum
        DEBUG_VAR(q);
        std::complex<SType> result(0);
        for (int i(0); i < 16; ++i) {
          std::complex<SType> factor = (i/4>0)?-1:1;
          result += factor * (*j1)[i % 4] * q[i / 4] * (*j0)[i];
        }
        CScalarType *j(CScalarType::New(result));
        DEBUG_VAR(*j);
        return j;
      }
      else if (p_v->V()->id.back()==0){
        CRaritaSchwingerType *j0(jj[0]->Get<CRaritaSchwingerType>());
        CScalarType *j1(jj[1]->Get<CScalarType>());
        DEBUG_VAR(*j0);
        DEBUG_VAR(*j1);
        const ATOOLS::Vec4<SType> &q(p_v->J(1)->P()); // Scalar four momentum
        DEBUG_VAR(q);
        std::complex<SType> result(0);
        CSpinorType *j(CSpinorType::New(j0->R(),j0->B(),0,0,j0->H()|j1->H(),j0->S()|j1->S(),1));
        for (int i(0); i < 16; ++i) {
          std::complex<SType> factor = (i/4>0)?-1:1;
          (*j)[i%4] += factor * (*j0)[i] * q[i / 4] * (*j1)[0];
        }
        DEBUG_VAR(*j);
        return j;
      }
      else
        THROW(not_implemented, "Vertex rotation in RFS not implemented!")
    }
  };// end of class RFS_Calculator

  template class RFS_Calculator<double>;

}// end of namespace METOOLS

using namespace METOOLS;

DECLARE_GETTER(RFS_Calculator<double>,"DRFS",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,RFS_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new RFS_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    RFS_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"RFS vertex"; }
