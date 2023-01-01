#include "METOOLS/Explicit/Lorentz_Calculator.H"
#include "METOOLS/Currents/C_Scalar.H"
#include "METOOLS/Currents/C_Spinor.H"
#include "METOOLS/Currents/C_Vector.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/MyComplex.H"


namespace METOOLS {

  template <typename SType>
  class VV3_Calculator: public Lorentz_Calculator {
  public:

    typedef std::complex<SType> SComplex;

    const SComplex I = SComplex(0.0,1.0);
    
    VV3_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "VV3"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {

// if outgoing index is 0
if (p_v->V()->id.back()==0){
const CVec4 <SType> & j1 = *(jj[0]->Get< CVec4 <SType> >());
const SComplex & j10 = j1[0];
const SComplex & j11 = j1[ATOOLS::Spinor<SType>::R1()];
const SComplex & j12 = j1[ATOOLS::Spinor<SType>::R2()];
const SComplex & j13 = j1[ATOOLS::Spinor<SType>::R3()];
const ATOOLS::Vec4D & p1 = p_v->J(0)->P();
const double& p10 = p1[0];
const double& p11 = p1[ATOOLS::Spinor<SType>::R1()];
const double& p12 = p1[ATOOLS::Spinor<SType>::R2()];
const double& p13 = p1[ATOOLS::Spinor<SType>::R3()];
ATOOLS::Vec4D p0 = -p1;
const double& p00 = p0[0];
const double& p01 = p0[ATOOLS::Spinor<SType>::R1()];
const double& p02 = p0[ATOOLS::Spinor<SType>::R2()];
const double& p03 = p0[ATOOLS::Spinor<SType>::R3()];
CVec4<SType>* j0 = NULL;
j0 = CVec4<SType>::New();
(*j0)[0] = 1.0*j10*pow(p10, 2) - 1.0*j10*pow(p11, 2) - 1.0*j10*pow(p12, 2) - 1.0*j10*pow(p13, 2);
(*j0)[ATOOLS::Spinor<SType>::R1()] = 1.0*j11*pow(p10, 2) - 1.0*j11*pow(p11, 2) - 1.0*j11*pow(p12, 2) - 1.0*j11*pow(p13, 2);
(*j0)[ATOOLS::Spinor<SType>::R2()] = 1.0*j12*pow(p10, 2) - 1.0*j12*pow(p11, 2) - 1.0*j12*pow(p12, 2) - 1.0*j12*pow(p13, 2);
(*j0)[ATOOLS::Spinor<SType>::R3()] = 1.0*j13*pow(p10, 2) - 1.0*j13*pow(p11, 2) - 1.0*j13*pow(p12, 2) - 1.0*j13*pow(p13, 2);
j0->SetS(j1.S());
return j0;

}

// if outgoing index is 1
if (p_v->V()->id.back()==1){
const CVec4 <SType> & j0 = *(jj[0]->Get< CVec4 <SType> >());
const SComplex & j00 = j0[0];
const SComplex & j01 = j0[ATOOLS::Spinor<SType>::R1()];
const SComplex & j02 = j0[ATOOLS::Spinor<SType>::R2()];
const SComplex & j03 = j0[ATOOLS::Spinor<SType>::R3()];
const ATOOLS::Vec4D & p0 = p_v->J(0)->P();
const double& p00 = p0[0];
const double& p01 = p0[ATOOLS::Spinor<SType>::R1()];
const double& p02 = p0[ATOOLS::Spinor<SType>::R2()];
const double& p03 = p0[ATOOLS::Spinor<SType>::R3()];
ATOOLS::Vec4D p1 = -p0;
const double& p10 = p1[0];
const double& p11 = p1[ATOOLS::Spinor<SType>::R1()];
const double& p12 = p1[ATOOLS::Spinor<SType>::R2()];
const double& p13 = p1[ATOOLS::Spinor<SType>::R3()];
CVec4<SType>* j1 = NULL;
j1 = CVec4<SType>::New();
(*j1)[0] = 1.0*j00*pow(p10, 2) - 1.0*j00*pow(p11, 2) - 1.0*j00*pow(p12, 2) - 1.0*j00*pow(p13, 2);
(*j1)[ATOOLS::Spinor<SType>::R1()] = 1.0*j01*pow(p10, 2) - 1.0*j01*pow(p11, 2) - 1.0*j01*pow(p12, 2) - 1.0*j01*pow(p13, 2);
(*j1)[ATOOLS::Spinor<SType>::R2()] = 1.0*j02*pow(p10, 2) - 1.0*j02*pow(p11, 2) - 1.0*j02*pow(p12, 2) - 1.0*j02*pow(p13, 2);
(*j1)[ATOOLS::Spinor<SType>::R3()] = 1.0*j03*pow(p10, 2) - 1.0*j03*pow(p11, 2) - 1.0*j03*pow(p12, 2) - 1.0*j03*pow(p13, 2);
j1->SetS(j0.S());
return j1;

}

      THROW(fatal_error, "Internal error in Lorentz calculator");
      return NULL;
    }

  };// end of class VV3_Calculator

  template class VV3_Calculator<double>;

}// end of namespace METOOLS


using namespace METOOLS;

DECLARE_GETTER(VV3_Calculator<double>,"DVV3",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,VV3_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new VV3_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    VV3_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"VV3 vertex"; }
