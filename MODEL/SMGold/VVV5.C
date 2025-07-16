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
  class VVV5_Calculator: public Lorentz_Calculator {
  public:

    typedef std::complex<SType> SComplex;

    const SComplex I = SComplex(0.0,1.0);
    
    VVV5_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "VVV5"; }

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
const CVec4 <SType> & j2 = *(jj[1]->Get< CVec4 <SType> >());
const SComplex & j20 = j2[0];
const SComplex & j21 = j2[ATOOLS::Spinor<SType>::R1()];
const SComplex & j22 = j2[ATOOLS::Spinor<SType>::R2()];
const SComplex & j23 = j2[ATOOLS::Spinor<SType>::R3()];
const ATOOLS::Vec4D & p2 = p_v->J(1)->P();
const double& p20 = p2[0];
const double& p21 = p2[ATOOLS::Spinor<SType>::R1()];
const double& p22 = p2[ATOOLS::Spinor<SType>::R2()];
const double& p23 = p2[ATOOLS::Spinor<SType>::R3()];
ATOOLS::Vec4D p0 = -p1-p2;
const double& p00 = p0[0];
const double& p01 = p0[ATOOLS::Spinor<SType>::R1()];
const double& p02 = p0[ATOOLS::Spinor<SType>::R2()];
const double& p03 = p0[ATOOLS::Spinor<SType>::R3()];
CVec4<SType>* j0 = NULL;
j0 = CVec4<SType>::New();
(*j0)[0] = 1.0*j10*j20*p10 - 1.0*j11*j20*p11 - 1.0*j12*j20*p12 - 1.0*j13*j20*p13;
(*j0)[ATOOLS::Spinor<SType>::R1()] = 1.0*j10*j21*p10 - 1.0*j11*j21*p11 - 1.0*j12*j21*p12 - 1.0*j13*j21*p13;
(*j0)[ATOOLS::Spinor<SType>::R2()] = 1.0*j10*j22*p10 - 1.0*j11*j22*p11 - 1.0*j12*j22*p12 - 1.0*j13*j22*p13;
(*j0)[ATOOLS::Spinor<SType>::R3()] = 1.0*j10*j23*p10 - 1.0*j11*j23*p11 - 1.0*j12*j23*p12 - 1.0*j13*j23*p13;
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
const ATOOLS::Vec4D & p0 = p_v->J(1)->P();
const double& p00 = p0[0];
const double& p01 = p0[ATOOLS::Spinor<SType>::R1()];
const double& p02 = p0[ATOOLS::Spinor<SType>::R2()];
const double& p03 = p0[ATOOLS::Spinor<SType>::R3()];
const CVec4 <SType> & j2 = *(jj[0]->Get< CVec4 <SType> >());
const SComplex & j20 = j2[0];
const SComplex & j21 = j2[ATOOLS::Spinor<SType>::R1()];
const SComplex & j22 = j2[ATOOLS::Spinor<SType>::R2()];
const SComplex & j23 = j2[ATOOLS::Spinor<SType>::R3()];
const ATOOLS::Vec4D & p2 = p_v->J(0)->P();
const double& p20 = p2[0];
const double& p21 = p2[ATOOLS::Spinor<SType>::R1()];
const double& p22 = p2[ATOOLS::Spinor<SType>::R2()];
const double& p23 = p2[ATOOLS::Spinor<SType>::R3()];
ATOOLS::Vec4D p1 = -p0-p2;
const double& p10 = p1[0];
const double& p11 = p1[ATOOLS::Spinor<SType>::R1()];
const double& p12 = p1[ATOOLS::Spinor<SType>::R2()];
const double& p13 = p1[ATOOLS::Spinor<SType>::R3()];
CVec4<SType>* j1 = NULL;
j1 = CVec4<SType>::New();
(*j1)[0] = 1.0*j00*j20*p10 - 1.0*j01*j21*p10 - 1.0*j02*j22*p10 - 1.0*j03*j23*p10;
(*j1)[ATOOLS::Spinor<SType>::R1()] = 1.0*j00*j20*p11 - 1.0*j01*j21*p11 - 1.0*j02*j22*p11 - 1.0*j03*j23*p11;
(*j1)[ATOOLS::Spinor<SType>::R2()] = 1.0*j00*j20*p12 - 1.0*j01*j21*p12 - 1.0*j02*j22*p12 - 1.0*j03*j23*p12;
(*j1)[ATOOLS::Spinor<SType>::R3()] = 1.0*j00*j20*p13 - 1.0*j01*j21*p13 - 1.0*j02*j22*p13 - 1.0*j03*j23*p13;
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
const ATOOLS::Vec4D & p0 = p_v->J(0)->P();
const double& p00 = p0[0];
const double& p01 = p0[ATOOLS::Spinor<SType>::R1()];
const double& p02 = p0[ATOOLS::Spinor<SType>::R2()];
const double& p03 = p0[ATOOLS::Spinor<SType>::R3()];
const CVec4 <SType> & j1 = *(jj[1]->Get< CVec4 <SType> >());
const SComplex & j10 = j1[0];
const SComplex & j11 = j1[ATOOLS::Spinor<SType>::R1()];
const SComplex & j12 = j1[ATOOLS::Spinor<SType>::R2()];
const SComplex & j13 = j1[ATOOLS::Spinor<SType>::R3()];
const ATOOLS::Vec4D & p1 = p_v->J(1)->P();
const double& p10 = p1[0];
const double& p11 = p1[ATOOLS::Spinor<SType>::R1()];
const double& p12 = p1[ATOOLS::Spinor<SType>::R2()];
const double& p13 = p1[ATOOLS::Spinor<SType>::R3()];
ATOOLS::Vec4D p2 = -p0-p1;
const double& p20 = p2[0];
const double& p21 = p2[ATOOLS::Spinor<SType>::R1()];
const double& p22 = p2[ATOOLS::Spinor<SType>::R2()];
const double& p23 = p2[ATOOLS::Spinor<SType>::R3()];
CVec4<SType>* j2 = NULL;
j2 = CVec4<SType>::New();
(*j2)[0] = 1.0*j00*j10*p10 - 1.0*j00*j11*p11 - 1.0*j00*j12*p12 - 1.0*j00*j13*p13;
(*j2)[ATOOLS::Spinor<SType>::R1()] = 1.0*j01*j10*p10 - 1.0*j01*j11*p11 - 1.0*j01*j12*p12 - 1.0*j01*j13*p13;
(*j2)[ATOOLS::Spinor<SType>::R2()] = 1.0*j02*j10*p10 - 1.0*j02*j11*p11 - 1.0*j02*j12*p12 - 1.0*j02*j13*p13;
(*j2)[ATOOLS::Spinor<SType>::R3()] = 1.0*j03*j10*p10 - 1.0*j03*j11*p11 - 1.0*j03*j12*p12 - 1.0*j03*j13*p13;
j2->SetS(j0.S()|j1.S());
return j2;

}

      THROW(fatal_error, "Internal error in Lorentz calculator");
      return NULL;
    }

  };// end of class VVV5_Calculator

  template class VVV5_Calculator<double>;

}// end of namespace METOOLS


using namespace METOOLS;

DECLARE_GETTER(VVV5_Calculator<double>,"DVVV5",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,VVV5_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new VVV5_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    VVV5_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"VVV5 vertex"; }
