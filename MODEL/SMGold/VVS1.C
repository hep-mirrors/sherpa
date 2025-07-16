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
  class VVS1_Calculator: public Lorentz_Calculator {
  public:

    typedef std::complex<SType> SComplex;

    const SComplex I = SComplex(0.0,1.0);
    
    VVS1_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "VVS1"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {

// if outgoing index is 0
if (p_v->V()->id.back()==0){
const CVec4 <SType> & j1 = *(jj[0]->Get< CVec4 <SType> >());
const SComplex & j10 = j1[0];
const SComplex & j11 = j1[ATOOLS::Spinor<SType>::R1()];
const SComplex & j12 = j1[ATOOLS::Spinor<SType>::R2()];
const SComplex & j13 = j1[ATOOLS::Spinor<SType>::R3()];
const CScalar <SType> & j2 = *(jj[1]->Get< CScalar <SType> >());
const SComplex & j20 = j2[0];
CVec4<SType>* j0 = NULL;
j0 = CVec4<SType>::New();
(*j0)[0] = 1.0*j10*j20;
(*j0)[ATOOLS::Spinor<SType>::R1()] = 1.0*j11*j20;
(*j0)[ATOOLS::Spinor<SType>::R2()] = 1.0*j12*j20;
(*j0)[ATOOLS::Spinor<SType>::R3()] = 1.0*j13*j20;
j0->SetS(j1.S()|j2.S());
return j0;

}

// if outgoing index is 1
if (p_v->V()->id.back()==1){
const CVec4 <SType> & j0 = *(jj[1]->Get< CVec4 <SType> >());
const SComplex & j00 = j0[0];
const SComplex & j01 = j0[ATOOLS::Spinor<SType>::R1()];
const SComplex & j02 = j0[ATOOLS::Spinor<SType>::R2()];
const SComplex & j03 = j0[ATOOLS::Spinor<SType>::R3()];
const CScalar <SType> & j2 = *(jj[0]->Get< CScalar <SType> >());
const SComplex & j20 = j2[0];
CVec4<SType>* j1 = NULL;
j1 = CVec4<SType>::New();
(*j1)[0] = 1.0*j00*j20;
(*j1)[ATOOLS::Spinor<SType>::R1()] = 1.0*j01*j20;
(*j1)[ATOOLS::Spinor<SType>::R2()] = 1.0*j02*j20;
(*j1)[ATOOLS::Spinor<SType>::R3()] = 1.0*j03*j20;
j1->SetS(j0.S()|j2.S());
return j1;

}

// if outgoing index is 2
if (p_v->V()->id.back()==2){
const CVec4 <SType> & j0 = *(jj[0]->Get< CVec4 <SType> >());
const SComplex & j00 = j0[0];
const SComplex & j01 = j0[ATOOLS::Spinor<SType>::R1()];
const SComplex & j02 = j0[ATOOLS::Spinor<SType>::R2()];
const SComplex & j03 = j0[ATOOLS::Spinor<SType>::R3()];
const CVec4 <SType> & j1 = *(jj[1]->Get< CVec4 <SType> >());
const SComplex & j10 = j1[0];
const SComplex & j11 = j1[ATOOLS::Spinor<SType>::R1()];
const SComplex & j12 = j1[ATOOLS::Spinor<SType>::R2()];
const SComplex & j13 = j1[ATOOLS::Spinor<SType>::R3()];
CScalar<SType>* j2 = NULL;
j2 = CScalar<SType>::New(1.0*j00*j10 - 1.0*j01*j11 - 1.0*j02*j12 - 1.0*j03*j13);
j2->SetS(j0.S()|j1.S());
return j2;

}

      THROW(fatal_error, "Internal error in Lorentz calculator");
      return NULL;
    }

  };// end of class VVS1_Calculator

  template class VVS1_Calculator<double>;

}// end of namespace METOOLS


using namespace METOOLS;

DECLARE_GETTER(VVS1_Calculator<double>,"DVVS1",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,VVS1_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new VVS1_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    VVS1_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"VVS1 vertex"; }
