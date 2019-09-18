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
  class FFS1_Calculator: public Lorentz_Calculator {
  public:

    typedef std::complex<SType> SComplex;

    const SComplex I = SComplex(0.0,1.0);
    
    FFS1_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "FFS1"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {

      // if outgoing index is 0
      if (p_v->V()->id.back()==0){
	const CSpinor <SType> & j1 = ((jj[0]->Get< CSpinor <SType> >())->B() == -1) ? (*(jj[0]->Get< CSpinor <SType> >())) : (*(jj[0]->Get< CSpinor <SType> >())).CConj() ;
	const SComplex & j10 = j1[0];
	const SComplex & j11 = j1[1];
	const SComplex & j12 = j1[2];
	const SComplex & j13 = j1[3];
	const CScalar <SType> & j2 = *(jj[1]->Get< CScalar <SType> >());
	const SComplex & j20 = j2[0];
	CSpinor<SType>* j0 = NULL;
	switch(+(j1.On()<<(2))){
	case 4:
	  j0 = CSpinor<SType>::New(m_r[0],-1,0,0,0,0,1);
	  (*j0)[0] = 1.0*j10*j20;
	  (*j0)[1] = 1.0*j11*j20;
	  (*j0)[2] = 0;
	  (*j0)[3] = 0;
	  j0->SetS(j1.S()|j2.S());
	  return j0;
	case 8:
	  return j0;
	case 12:
	  j0 = CSpinor<SType>::New(m_r[0],-1,0,0,0,0,1);
	  (*j0)[0] = 1.0*j10*j20;
	  (*j0)[1] = 1.0*j11*j20;
	  (*j0)[2] = 0;
	  (*j0)[3] = 0;
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
	const CScalar <SType> & j2 = *(jj[0]->Get< CScalar <SType> >());
	const SComplex & j20 = j2[0];
	CSpinor<SType>* j1 = NULL;
	switch(+(j0.On()<<(0))){
	case 1:
	  j1 = CSpinor<SType>::New(m_r[1],1,0,0,0,0,1);
	  (*j1)[0] = 1.0*j00*j20;
	  (*j1)[1] = 1.0*j01*j20;
	  (*j1)[2] = 0;
	  (*j1)[3] = 0;
	  j1->SetS(j0.S()|j2.S());
	  return j1;
	case 2:
	  return j1;
	case 3:
	  j1 = CSpinor<SType>::New(m_r[1],1,0,0,0,0,1);
	  (*j1)[0] = 1.0*j00*j20;
	  (*j1)[1] = 1.0*j01*j20;
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
	const CSpinor <SType> & j1 = ((jj[1]->Get< CSpinor <SType> >())->B() == -1) ? (*(jj[1]->Get< CSpinor <SType> >())) : (*(jj[1]->Get< CSpinor <SType> >())).CConj() ;
	const SComplex & j10 = j1[0];
	const SComplex & j11 = j1[1];
	const SComplex & j12 = j1[2];
	const SComplex & j13 = j1[3];
	CScalar<SType>* j2 = NULL;
	switch(+(j0.On()<<(0))+(j1.On()<<(2))){
	case 5:
	  j2 = CScalar<SType>::New(1.0*j00*j10 + 1.0*j01*j11);
	  j2->SetS(j0.S()|j1.S());
	  return j2;
	case 6:
	  return j2;
	case 7:
	  j2 = CScalar<SType>::New(1.0*j00*j10 + 1.0*j01*j11);
	  j2->SetS(j0.S()|j1.S());
	  return j2;
	case 9:
	  return j2;
	case 10:
	  return j2;
	case 11:
	  return j2;
	case 13:
	  j2 = CScalar<SType>::New(1.0*j00*j10 + 1.0*j01*j11);
	  j2->SetS(j0.S()|j1.S());
	  return j2;
	case 14:
	  return j2;
	case 15:
	  j2 = CScalar<SType>::New(1.0*j00*j10 + 1.0*j01*j11);
	  j2->SetS(j0.S()|j1.S());
	  return j2;
	default:
	  THROW(fatal_error, "Massless spinor optimization error in Lorentz calculator");
	}

      }

      THROW(fatal_error, "Internal error in Lorentz calculator");
      return NULL;
    }

  };// end of class FFS1_Calculator

  template class FFS1_Calculator<double>;

}// end of namespace METOOLS


using namespace METOOLS;

DECLARE_GETTER(FFS1_Calculator<double>,"DFFS1",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,FFS1_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new FFS1_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    FFS1_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"FFS1 vertex"; }


namespace METOOLS {

  template <typename SType>
  class FFS3_Calculator: public Lorentz_Calculator {
  public:

    typedef std::complex<SType> SComplex;

    const SComplex I = SComplex(0.0,1.0);
    
    FFS3_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "FFS3"; }

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
	const CScalar <SType> & j2 = *(jj[1]->Get< CScalar <SType> >());
	const SComplex & j20 = j2[0];
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
	  return j0;
	case 8:
	  j0 = CSpinor<SType>::New(m_r[0],-1,0,0,0,0,2);
	  (*j0)[0] = 0;
	  (*j0)[1] = 0;
	  (*j0)[2] = 1.0*j12*j20;
	  (*j0)[3] = 1.0*j13*j20;
	  j0->SetS(j1.S()|j2.S());
	  return j0;
	case 12:
	  j0 = CSpinor<SType>::New(m_r[0],-1,0,0,0,0,2);
	  (*j0)[0] = 0;
	  (*j0)[1] = 0;
	  (*j0)[2] = 1.0*j12*j20;
	  (*j0)[3] = 1.0*j13*j20;
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
	const CScalar <SType> & j2 = *(jj[0]->Get< CScalar <SType> >());
	const SComplex & j20 = j2[0];
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
	  j1 = CSpinor<SType>::New(m_r[1],1,0,0,0,0,2);
	  (*j1)[0] = 0;
	  (*j1)[1] = 0;
	  (*j1)[2] = 1.0*j02*j20;
	  (*j1)[3] = 1.0*j03*j20;
	  j1->SetS(j0.S()|j2.S());
	  return j1;
	case 3:
	  j1 = CSpinor<SType>::New(m_r[1],1,0,0,0,0,2);
	  (*j1)[0] = 0;
	  (*j1)[1] = 0;
	  (*j1)[2] = 1.0*j02*j20;
	  (*j1)[3] = 1.0*j03*j20;
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
	CScalar<SType>* j2 = NULL;
	switch(+(j0.On()<<(0))+(j1.On()<<(2))){
	case 5:
	  return j2;
	case 6:
	  return j2;
	case 7:
	  return j2;
	case 9:
	  return j2;
	case 10:
	  j2 = CScalar<SType>::New(1.0*j02*j12 + 1.0*j03*j13);
	  j2->SetS(j0.S()|j1.S());
	  return j2;
	case 11:
	  j2 = CScalar<SType>::New(1.0*j02*j12 + 1.0*j03*j13);
	  j2->SetS(j0.S()|j1.S());
	  return j2;
	case 13:
	  return j2;
	case 14:
	  j2 = CScalar<SType>::New(1.0*j02*j12 + 1.0*j03*j13);
	  j2->SetS(j0.S()|j1.S());
	  return j2;
	case 15:
	  j2 = CScalar<SType>::New(1.0*j02*j12 + 1.0*j03*j13);
	  j2->SetS(j0.S()|j1.S());
	  return j2;
	default:
	  THROW(fatal_error, "Massless spinor optimization error in Lorentz calculator");
	}

      }

      THROW(fatal_error, "Internal error in Lorentz calculator");
      return NULL;
    }

  };// end of class FFS3_Calculator

  template class FFS3_Calculator<double>;

}// end of namespace METOOLS


using namespace METOOLS;

DECLARE_GETTER(FFS3_Calculator<double>,"DFFS3",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,FFS3_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new FFS3_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    FFS3_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"FFS3 vertex"; }
