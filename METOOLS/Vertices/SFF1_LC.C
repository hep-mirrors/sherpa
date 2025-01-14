#include "METOOLS/Explicit/Lorentz_Calculator.H"
#include "METOOLS/Currents/C_Scalar.H"
#include "METOOLS/Currents/C_Spinor.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "METOOLS/Main/SpinFuncs.H"

using namespace ATOOLS;

namespace METOOLS {

  template <typename SType>
  class SFF1_Calculator: public Lorentz_Calculator {
  public:

    typedef CSpinor<SType> CSpinorType;
    typedef CScalar<SType> CScalarType;

    SFF1_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "SSF1"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      if (p_v->V()->id.back()==2) {
        const CSpinor <SType> & j0 = (*jj[0]->Get< CSpinor <SType> >()) ;
        const CSpinor <SType> & j1 = (*jj[1]->Get< CSpinor <SType> >()) ;

        // TODO: q hier richtig berechnet?
        const ATOOLS::Vec4<SType> &q(-p_v->J(1)->P() - p_v->J(0)->P()); // Scalar four momentum
        DEBUG_VAR(q);
        DEBUG_VAR(j0);
        DEBUG_VAR(j1);
        METOOLS::Gamma gammavec = Gamma<SType>();
        ATOOLS::TCMatrix<SType> gammamugamma5q = (gammavec[0] * q[0] + std::complex<SType>(-1., 0) * gammavec[1] * q[1]
                                                  + std::complex<SType>(-1., 0) * gammavec[2] * q[2] +
                                                  std::complex<SType>(-1., 0) * gammavec[3] * q[3]) * gammavec.Gamma5();
        std::complex<SType> result(0);
        for (int i(0); i < 4; ++i) {
          for (int k(0); k < 4; ++k) {
            if (j1.B()==-1 and j0.B()==1)
            result += j1[i] * gammamugamma5q[i][k] * j0[k];
            else if (j1.B()==1 and j0.B()==-1){
              result += j0[i] * gammamugamma5q[i][k] * j1[k];
            }
            else
              THROW(inconsistent_option, "Both spinors adjoint or non-adjoint")
          }
        }
        CScalarType *j(CScalarType::New(result));
        DEBUG_VAR(*j);
        return j;
      }
      else if (p_v->V()->id.back()==0) {
        const CSpinorType j0(*(jj[0]->Get< CSpinor <SType> >()));
        const CScalarType j1(*(jj[1 - p_v->V()->id.back()]->Get<CScalarType>()));
        const ATOOLS::Vec4<SType> &q(p_v->J(1)->P()); // Scalar four momentum
        DEBUG_VAR(q);
        DEBUG_VAR(j0);
        DEBUG_VAR(j1);
        METOOLS::Gamma gammavec = Gamma<SType>();
        ATOOLS::TCMatrix<SType> gammamugamma5q = (gammavec[0] * q[0] + std::complex<SType>(-1., 0) * gammavec[1] * q[1]
                                                  + std::complex<SType>(-1., 0) * gammavec[2] * q[2] +
                                                  std::complex<SType>(-1., 0) * gammavec[3] * q[3]) * gammavec.Gamma5();
        CSpinorType *j(CSpinorType::New(j0.R(), j0.B(), 0, 0, j0.H() | j1.H(), j0.S() | j1.S(), 3));
        switch (j0.B()){
          case -1: {
            for (int i(0); i < 4; ++i) {
              for (int k(0); k < 4; ++k) {
                (*j)[i] += j0[k] * gammamugamma5q[k][i] * j1[0];
              }
            }
            DEBUG_VAR(*j);
            return j;
          }
          case 1:{
            for (int i(0); i < 4; ++i) {
              for (int k(0); k < 4; ++k) {
                (*j)[i] += gammamugamma5q[i][k] * j0[k] * j1[0];
              }
            }
            DEBUG_VAR(*j);
            return j;
          }
        }
      }
      else
        THROW(not_implemented, "SSF1_LC vertex rotation not implemented!")
    }
  };// end of class SFF1_Calculator

  template class SFF1_Calculator<double>;

}// end of namespace METOOLS

using namespace METOOLS;

DECLARE_GETTER(SFF1_Calculator<double>,"DSFF1",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,SFF1_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new SFF1_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    SFF1_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"SFF1 vertex"; }
