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
  class VVV3_Calculator: public Lorentz_Calculator {
  public:

    typedef std::complex<SType> SComplex;

    const SComplex I = SComplex(0.0,1.0);
    
    VVV3_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "VVV3"; }

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
(*j0)[0] = 1.0*j10*j20*p20 - 1.0*j10*j21*p21 - 1.0*j10*j22*p22 - 1.0*j10*j23*p23;
(*j0)[ATOOLS::Spinor<SType>::R1()] = 1.0*j11*j20*p20 - 1.0*j11*j21*p21 - 1.0*j11*j22*p22 - 1.0*j11*j23*p23;
(*j0)[ATOOLS::Spinor<SType>::R2()] = 1.0*j12*j20*p20 - 1.0*j12*j21*p21 - 1.0*j12*j22*p22 - 1.0*j12*j23*p23;
(*j0)[ATOOLS::Spinor<SType>::R3()] = 1.0*j13*j20*p20 - 1.0*j13*j21*p21 - 1.0*j13*j22*p22 - 1.0*j13*j23*p23;
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
(*j1)[0] = 1.0*j00*j20*p20 - 1.0*j00*j21*p21 - 1.0*j00*j22*p22 - 1.0*j00*j23*p23;
(*j1)[ATOOLS::Spinor<SType>::R1()] = 1.0*j01*j20*p20 - 1.0*j01*j21*p21 - 1.0*j01*j22*p22 - 1.0*j01*j23*p23;
(*j1)[ATOOLS::Spinor<SType>::R2()] = 1.0*j02*j20*p20 - 1.0*j02*j21*p21 - 1.0*j02*j22*p22 - 1.0*j02*j23*p23;
(*j1)[ATOOLS::Spinor<SType>::R3()] = 1.0*j03*j20*p20 - 1.0*j03*j21*p21 - 1.0*j03*j22*p22 - 1.0*j03*j23*p23;
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
(*j2)[0] = 1.0*j00*j10*p20 - 1.0*j01*j11*p20 - 1.0*j02*j12*p20 - 1.0*j03*j13*p20;
(*j2)[ATOOLS::Spinor<SType>::R1()] = 1.0*j00*j10*p21 - 1.0*j01*j11*p21 - 1.0*j02*j12*p21 - 1.0*j03*j13*p21;
(*j2)[ATOOLS::Spinor<SType>::R2()] = 1.0*j00*j10*p22 - 1.0*j01*j11*p22 - 1.0*j02*j12*p22 - 1.0*j03*j13*p22;
(*j2)[ATOOLS::Spinor<SType>::R3()] = 1.0*j00*j10*p23 - 1.0*j01*j11*p23 - 1.0*j02*j12*p23 - 1.0*j03*j13*p23;
j2->SetS(j0.S()|j1.S());
return j2;

}

      THROW(fatal_error, "Internal error in Lorentz calculator");
      return NULL;
    }

  };// end of class VVV3_Calculator

  template class VVV3_Calculator<double>;

}// end of namespace METOOLS


using namespace METOOLS;

DECLARE_GETTER(VVV3_Calculator<double>,"DVVV3",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,VVV3_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new VVV3_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    VVV3_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"VVV3 vertex"; }
