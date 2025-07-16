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
  class VVVV4_Calculator: public Lorentz_Calculator {
  public:

    typedef std::complex<SType> SComplex;

    const SComplex I = SComplex(0.0,1.0);
    
    VVVV4_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "VVVV4"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {

// if outgoing index is 0
if (p_v->V()->id.back()==0){
const CVec4 <SType> & j1 = *(jj[0]->Get< CVec4 <SType> >());
const SComplex & j10 = j1[0];
const SComplex & j11 = j1[ATOOLS::Spinor<SType>::R1()];
const SComplex & j12 = j1[ATOOLS::Spinor<SType>::R2()];
const SComplex & j13 = j1[ATOOLS::Spinor<SType>::R3()];
const CVec4 <SType> & j2 = *(jj[1]->Get< CVec4 <SType> >());
const SComplex & j20 = j2[0];
const SComplex & j21 = j2[ATOOLS::Spinor<SType>::R1()];
const SComplex & j22 = j2[ATOOLS::Spinor<SType>::R2()];
const SComplex & j23 = j2[ATOOLS::Spinor<SType>::R3()];
const CVec4 <SType> & j3 = *(jj[2]->Get< CVec4 <SType> >());
const SComplex & j30 = j3[0];
const SComplex & j31 = j3[ATOOLS::Spinor<SType>::R1()];
const SComplex & j32 = j3[ATOOLS::Spinor<SType>::R2()];
const SComplex & j33 = j3[ATOOLS::Spinor<SType>::R3()];
CVec4<SType>* j0 = NULL;
j0 = CVec4<SType>::New();
(*j0)[0] = 1.0*j10*j20*j30 - 1.0*j10*j21*j31 - 1.0*j10*j22*j32 - 1.0*j10*j23*j33;
(*j0)[ATOOLS::Spinor<SType>::R1()] = 1.0*j11*j20*j30 - 1.0*j11*j21*j31 - 1.0*j11*j22*j32 - 1.0*j11*j23*j33;
(*j0)[ATOOLS::Spinor<SType>::R2()] = 1.0*j12*j20*j30 - 1.0*j12*j21*j31 - 1.0*j12*j22*j32 - 1.0*j12*j23*j33;
(*j0)[ATOOLS::Spinor<SType>::R3()] = 1.0*j13*j20*j30 - 1.0*j13*j21*j31 - 1.0*j13*j22*j32 - 1.0*j13*j23*j33;
j0->SetS(j1.S()|j2.S()|j3.S());
return j0;

}

// if outgoing index is 1
if (p_v->V()->id.back()==1){
const CVec4 <SType> & j0 = *(jj[2]->Get< CVec4 <SType> >());
const SComplex & j00 = j0[0];
const SComplex & j01 = j0[ATOOLS::Spinor<SType>::R1()];
const SComplex & j02 = j0[ATOOLS::Spinor<SType>::R2()];
const SComplex & j03 = j0[ATOOLS::Spinor<SType>::R3()];
const CVec4 <SType> & j2 = *(jj[0]->Get< CVec4 <SType> >());
const SComplex & j20 = j2[0];
const SComplex & j21 = j2[ATOOLS::Spinor<SType>::R1()];
const SComplex & j22 = j2[ATOOLS::Spinor<SType>::R2()];
const SComplex & j23 = j2[ATOOLS::Spinor<SType>::R3()];
const CVec4 <SType> & j3 = *(jj[1]->Get< CVec4 <SType> >());
const SComplex & j30 = j3[0];
const SComplex & j31 = j3[ATOOLS::Spinor<SType>::R1()];
const SComplex & j32 = j3[ATOOLS::Spinor<SType>::R2()];
const SComplex & j33 = j3[ATOOLS::Spinor<SType>::R3()];
CVec4<SType>* j1 = NULL;
j1 = CVec4<SType>::New();
(*j1)[0] = 1.0*j00*j20*j30 - 1.0*j00*j21*j31 - 1.0*j00*j22*j32 - 1.0*j00*j23*j33;
(*j1)[ATOOLS::Spinor<SType>::R1()] = 1.0*j01*j20*j30 - 1.0*j01*j21*j31 - 1.0*j01*j22*j32 - 1.0*j01*j23*j33;
(*j1)[ATOOLS::Spinor<SType>::R2()] = 1.0*j02*j20*j30 - 1.0*j02*j21*j31 - 1.0*j02*j22*j32 - 1.0*j02*j23*j33;
(*j1)[ATOOLS::Spinor<SType>::R3()] = 1.0*j03*j20*j30 - 1.0*j03*j21*j31 - 1.0*j03*j22*j32 - 1.0*j03*j23*j33;
j1->SetS(j0.S()|j2.S()|j3.S());
return j1;

}

// if outgoing index is 2
if (p_v->V()->id.back()==2){
const CVec4 <SType> & j0 = *(jj[1]->Get< CVec4 <SType> >());
const SComplex & j00 = j0[0];
const SComplex & j01 = j0[ATOOLS::Spinor<SType>::R1()];
const SComplex & j02 = j0[ATOOLS::Spinor<SType>::R2()];
const SComplex & j03 = j0[ATOOLS::Spinor<SType>::R3()];
const CVec4 <SType> & j1 = *(jj[2]->Get< CVec4 <SType> >());
const SComplex & j10 = j1[0];
const SComplex & j11 = j1[ATOOLS::Spinor<SType>::R1()];
const SComplex & j12 = j1[ATOOLS::Spinor<SType>::R2()];
const SComplex & j13 = j1[ATOOLS::Spinor<SType>::R3()];
const CVec4 <SType> & j3 = *(jj[0]->Get< CVec4 <SType> >());
const SComplex & j30 = j3[0];
const SComplex & j31 = j3[ATOOLS::Spinor<SType>::R1()];
const SComplex & j32 = j3[ATOOLS::Spinor<SType>::R2()];
const SComplex & j33 = j3[ATOOLS::Spinor<SType>::R3()];
CVec4<SType>* j2 = NULL;
j2 = CVec4<SType>::New();
(*j2)[0] = 1.0*j00*j10*j30 - 1.0*j01*j11*j30 - 1.0*j02*j12*j30 - 1.0*j03*j13*j30;
(*j2)[ATOOLS::Spinor<SType>::R1()] = 1.0*j00*j10*j31 - 1.0*j01*j11*j31 - 1.0*j02*j12*j31 - 1.0*j03*j13*j31;
(*j2)[ATOOLS::Spinor<SType>::R2()] = 1.0*j00*j10*j32 - 1.0*j01*j11*j32 - 1.0*j02*j12*j32 - 1.0*j03*j13*j32;
(*j2)[ATOOLS::Spinor<SType>::R3()] = 1.0*j00*j10*j33 - 1.0*j01*j11*j33 - 1.0*j02*j12*j33 - 1.0*j03*j13*j33;
j2->SetS(j0.S()|j1.S()|j3.S());
return j2;

}

// if outgoing index is 3
if (p_v->V()->id.back()==3){
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
const CVec4 <SType> & j2 = *(jj[2]->Get< CVec4 <SType> >());
const SComplex & j20 = j2[0];
const SComplex & j21 = j2[ATOOLS::Spinor<SType>::R1()];
const SComplex & j22 = j2[ATOOLS::Spinor<SType>::R2()];
const SComplex & j23 = j2[ATOOLS::Spinor<SType>::R3()];
CVec4<SType>* j3 = NULL;
j3 = CVec4<SType>::New();
(*j3)[0] = 1.0*j00*j10*j20 - 1.0*j01*j11*j20 - 1.0*j02*j12*j20 - 1.0*j03*j13*j20;
(*j3)[ATOOLS::Spinor<SType>::R1()] = 1.0*j00*j10*j21 - 1.0*j01*j11*j21 - 1.0*j02*j12*j21 - 1.0*j03*j13*j21;
(*j3)[ATOOLS::Spinor<SType>::R2()] = 1.0*j00*j10*j22 - 1.0*j01*j11*j22 - 1.0*j02*j12*j22 - 1.0*j03*j13*j22;
(*j3)[ATOOLS::Spinor<SType>::R3()] = 1.0*j00*j10*j23 - 1.0*j01*j11*j23 - 1.0*j02*j12*j23 - 1.0*j03*j13*j23;
j3->SetS(j0.S()|j1.S()|j2.S());
return j3;

}

      THROW(fatal_error, "Internal error in Lorentz calculator");
      return NULL;
    }

  };// end of class VVVV4_Calculator

  template class VVVV4_Calculator<double>;

}// end of namespace METOOLS


using namespace METOOLS;

DECLARE_GETTER(VVVV4_Calculator<double>,"DVVVV4",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,VVVV4_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new VVVV4_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    VVVV4_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"VVVV4 vertex"; }
