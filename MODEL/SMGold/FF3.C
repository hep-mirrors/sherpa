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
  class FF3_Calculator: public Lorentz_Calculator {
  public:

    typedef std::complex<SType> SComplex;

    const SComplex I = SComplex(0.0,1.0);
    
    FF3_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "FF3"; }

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
ATOOLS::Vec4D p0 = -p1;
const double& p00 = p0[0];
const double& p01 = p0[ATOOLS::Spinor<SType>::R1()];
const double& p02 = p0[ATOOLS::Spinor<SType>::R2()];
const double& p03 = p0[ATOOLS::Spinor<SType>::R3()];
CSpinor<SType>* j0 = NULL;
switch(+(j1.On()<<(2))){
case 4:
return j0;
case 8:
j0 = CSpinor<SType>::New(m_r[0],-1,0,0,0,0,1);
(*j0)[0] = 1.0*j12*p00 + 1.0*j12*p03 + 1.0*j13*p01 + 1.0*I*j13*p02;
(*j0)[1] = 1.0*j12*p01 - 1.0*I*j12*p02 + 1.0*j13*p00 - 1.0*j13*p03;
(*j0)[2] = 0;
(*j0)[3] = 0;
j0->SetS(j1.S());
return j0;
case 12:
j0 = CSpinor<SType>::New(m_r[0],-1,0,0,0,0,1);
(*j0)[0] = 1.0*j12*p00 + 1.0*j12*p03 + 1.0*j13*p01 + 1.0*I*j13*p02;
(*j0)[1] = 1.0*j12*p01 - 1.0*I*j12*p02 + 1.0*j13*p00 - 1.0*j13*p03;
(*j0)[2] = 0;
(*j0)[3] = 0;
j0->SetS(j1.S());
return j0;
default:
 THROW(fatal_error, "Massless spinor optimization error in Lorentz calculator");
}

}

// if outgoing index is 1
if (p_v->V()->id.back()==1){
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
ATOOLS::Vec4D p1 = -p0;
const double& p10 = p1[0];
const double& p11 = p1[ATOOLS::Spinor<SType>::R1()];
const double& p12 = p1[ATOOLS::Spinor<SType>::R2()];
const double& p13 = p1[ATOOLS::Spinor<SType>::R3()];
CSpinor<SType>* j1 = NULL;
switch(+(j0.On()<<(0))){
case 1:
j1 = CSpinor<SType>::New(m_r[1],1,0,0,0,0,2);
(*j1)[0] = 0;
(*j1)[1] = 0;
(*j1)[2] = 1.0*j00*p00 + 1.0*j00*p03 + 1.0*j01*p01 - 1.0*I*j01*p02;
(*j1)[3] = 1.0*j00*p01 + 1.0*I*j00*p02 + 1.0*j01*p00 - 1.0*j01*p03;
j1->SetS(j0.S());
return j1;
case 2:
return j1;
case 3:
j1 = CSpinor<SType>::New(m_r[1],1,0,0,0,0,2);
(*j1)[0] = 0;
(*j1)[1] = 0;
(*j1)[2] = 1.0*j00*p00 + 1.0*j00*p03 + 1.0*j01*p01 - 1.0*I*j01*p02;
(*j1)[3] = 1.0*j00*p01 + 1.0*I*j00*p02 + 1.0*j01*p00 - 1.0*j01*p03;
j1->SetS(j0.S());
return j1;
default:
 THROW(fatal_error, "Massless spinor optimization error in Lorentz calculator");
}

}

      THROW(fatal_error, "Internal error in Lorentz calculator");
      return NULL;
    }

  };// end of class FF3_Calculator

  template class FF3_Calculator<double>;

}// end of namespace METOOLS


using namespace METOOLS;

DECLARE_GETTER(FF3_Calculator<double>,"DFF3",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,FF3_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new FF3_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    FF3_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"FF3 vertex"; }
