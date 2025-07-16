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
  class FFV3_Calculator: public Lorentz_Calculator {
  public:

    typedef std::complex<SType> SComplex;

    const SComplex I = SComplex(0.0,1.0);
    
    FFV3_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "FFV3"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {

// if outgoing index is 0
if (p_v->V()->id.back()==0){
const CSpinor <SType> & j1 = ((jj[0]->Get< CSpinor <SType> >())->B() == -1) ? (*(jj[0]->Get< CSpinor <SType> >())) : (*(jj[0]->Get< CSpinor <SType> >())).CConj() ;
const SComplex & j10 = j1[0];
const SComplex & j11 = j1[1];
const SComplex & j12 = j1[2];
const SComplex & j13 = j1[3];
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
CSpinor<SType>* j0 = NULL;
switch(+(j1.On()<<(2))){
case 4:
j0 = CSpinor<SType>::New(m_r[0],-1,0,0,0,0,2);
(*j0)[0] = 0;
(*j0)[1] = 0;
(*j0)[2] = 1.0*j10*j20 - 1.0*j10*j23 - 1.0*j11*j21 - 1.0*I*j11*j22;
(*j0)[3] = -1.0*j10*j21 + 1.0*I*j10*j22 + 1.0*j11*j20 + 1.0*j11*j23;
j0->SetS(j1.S()|j2.S());
return j0;
case 8:
return j0;
case 12:
j0 = CSpinor<SType>::New(m_r[0],-1,0,0,0,0,2);
(*j0)[0] = 0;
(*j0)[1] = 0;
(*j0)[2] = 1.0*j10*j20 - 1.0*j10*j23 - 1.0*j11*j21 - 1.0*I*j11*j22;
(*j0)[3] = -1.0*j10*j21 + 1.0*I*j10*j22 + 1.0*j11*j20 + 1.0*j11*j23;
j0->SetS(j1.S()|j2.S());
return j0;
default:
 THROW(fatal_error, "Massless spinor optimization error in Lorentz calculator");
}

}

// if outgoing index is 1
if (p_v->V()->id.back()==1){
const CSpinor <SType> & j0 = ((jj[1]->Get< CSpinor <SType> >())->B() == 1) ? (*(jj[1]->Get< CSpinor <SType> >())) : (*(jj[1]->Get< CSpinor <SType> >())).CConj() ;
const SComplex & j00 = j0[0];
const SComplex & j01 = j0[1];
const SComplex & j02 = j0[2];
const SComplex & j03 = j0[3];
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
CSpinor<SType>* j1 = NULL;
switch(+(j0.On()<<(0))){
case 1:
return j1;
case 2:
j1 = CSpinor<SType>::New(m_r[1],1,0,0,0,0,1);
(*j1)[0] = 1.0*j02*j20 - 1.0*j02*j23 - 1.0*j03*j21 + 1.0*I*j03*j22;
(*j1)[1] = -1.0*j02*j21 - 1.0*I*j02*j22 + 1.0*j03*j20 + 1.0*j03*j23;
(*j1)[2] = 0;
(*j1)[3] = 0;
j1->SetS(j0.S()|j2.S());
return j1;
case 3:
j1 = CSpinor<SType>::New(m_r[1],1,0,0,0,0,1);
(*j1)[0] = 1.0*j02*j20 - 1.0*j02*j23 - 1.0*j03*j21 + 1.0*I*j03*j22;
(*j1)[1] = -1.0*j02*j21 - 1.0*I*j02*j22 + 1.0*j03*j20 + 1.0*j03*j23;
(*j1)[2] = 0;
(*j1)[3] = 0;
j1->SetS(j0.S()|j2.S());
return j1;
default:
 THROW(fatal_error, "Massless spinor optimization error in Lorentz calculator");
}

}

// if outgoing index is 2
if (p_v->V()->id.back()==2){
const CSpinor <SType> & j0 = ((jj[0]->Get< CSpinor <SType> >())->B() == 1) ? (*(jj[0]->Get< CSpinor <SType> >())) : (*(jj[0]->Get< CSpinor <SType> >())).CConj() ;
const SComplex & j00 = j0[0];
const SComplex & j01 = j0[1];
const SComplex & j02 = j0[2];
const SComplex & j03 = j0[3];
const ATOOLS::Vec4D & p0 = p_v->J(0)->P();
const double& p00 = p0[0];
const double& p01 = p0[ATOOLS::Spinor<SType>::R1()];
const double& p02 = p0[ATOOLS::Spinor<SType>::R2()];
const double& p03 = p0[ATOOLS::Spinor<SType>::R3()];
const CSpinor <SType> & j1 = ((jj[1]->Get< CSpinor <SType> >())->B() == -1) ? (*(jj[1]->Get< CSpinor <SType> >())) : (*(jj[1]->Get< CSpinor <SType> >())).CConj() ;
const SComplex & j10 = j1[0];
const SComplex & j11 = j1[1];
const SComplex & j12 = j1[2];
const SComplex & j13 = j1[3];
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
switch(+(j0.On()<<(0))+(j1.On()<<(2))){
case 5:
return j2;
case 6:
j2 = CVec4<SType>::New();
(*j2)[0] = 1.0*j02*j10 + 1.0*j03*j11;
(*j2)[ATOOLS::Spinor<SType>::R1()] = 1.0*j02*j11 + 1.0*j03*j10;
(*j2)[ATOOLS::Spinor<SType>::R2()] = 1.0*I*j02*j11 - 1.0*I*j03*j10;
(*j2)[ATOOLS::Spinor<SType>::R3()] = 1.0*j02*j10 - 1.0*j03*j11;
j2->SetS(j0.S()|j1.S());
return j2;
case 7:
j2 = CVec4<SType>::New();
(*j2)[0] = 1.0*j02*j10 + 1.0*j03*j11;
(*j2)[ATOOLS::Spinor<SType>::R1()] = 1.0*j02*j11 + 1.0*j03*j10;
(*j2)[ATOOLS::Spinor<SType>::R2()] = 1.0*I*j02*j11 - 1.0*I*j03*j10;
(*j2)[ATOOLS::Spinor<SType>::R3()] = 1.0*j02*j10 - 1.0*j03*j11;
j2->SetS(j0.S()|j1.S());
return j2;
case 9:
return j2;
case 10:
return j2;
case 11:
return j2;
case 13:
return j2;
case 14:
j2 = CVec4<SType>::New();
(*j2)[0] = 1.0*j02*j10 + 1.0*j03*j11;
(*j2)[ATOOLS::Spinor<SType>::R1()] = 1.0*j02*j11 + 1.0*j03*j10;
(*j2)[ATOOLS::Spinor<SType>::R2()] = 1.0*I*j02*j11 - 1.0*I*j03*j10;
(*j2)[ATOOLS::Spinor<SType>::R3()] = 1.0*j02*j10 - 1.0*j03*j11;
j2->SetS(j0.S()|j1.S());
return j2;
case 15:
j2 = CVec4<SType>::New();
(*j2)[0] = 1.0*j02*j10 + 1.0*j03*j11;
(*j2)[ATOOLS::Spinor<SType>::R1()] = 1.0*j02*j11 + 1.0*j03*j10;
(*j2)[ATOOLS::Spinor<SType>::R2()] = 1.0*I*j02*j11 - 1.0*I*j03*j10;
(*j2)[ATOOLS::Spinor<SType>::R3()] = 1.0*j02*j10 - 1.0*j03*j11;
j2->SetS(j0.S()|j1.S());
return j2;
default:
 THROW(fatal_error, "Massless spinor optimization error in Lorentz calculator");
}

}

      THROW(fatal_error, "Internal error in Lorentz calculator");
      return NULL;
    }

  };// end of class FFV3_Calculator

  template class FFV3_Calculator<double>;

}// end of namespace METOOLS


using namespace METOOLS;

DECLARE_GETTER(FFV3_Calculator<double>,"DFFV3",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,FFV3_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new FFV3_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    FFV3_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"FFV3 vertex"; }
