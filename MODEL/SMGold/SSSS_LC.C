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
  class SSSS1_Calculator: public Lorentz_Calculator {
  public:

    typedef std::complex<SType> SComplex;

    const SComplex I = SComplex(0.0,1.0);
    
    SSSS1_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "SSSS1"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {

      // if outgoing index is 0
      if (p_v->V()->id.back()==0){
	const CScalar <SType> & j1 = *(jj[0]->Get< CScalar <SType> >());
	const SComplex & j10 = j1[0];
	const CScalar <SType> & j2 = *(jj[1]->Get< CScalar <SType> >());
	const SComplex & j20 = j2[0];
	const CScalar <SType> & j3 = *(jj[2]->Get< CScalar <SType> >());
	const SComplex & j30 = j3[0];
	CScalar<SType>* j0 = NULL;
	j0 = CScalar<SType>::New(j10*j20*j30);
	j0->SetS(j1.S()|j2.S()|j3.S());
	return j0;

      }

      // if outgoing index is 1
      if (p_v->V()->id.back()==1){
	const CScalar <SType> & j0 = *(jj[2]->Get< CScalar <SType> >());
	const SComplex & j00 = j0[0];
	const CScalar <SType> & j2 = *(jj[0]->Get< CScalar <SType> >());
	const SComplex & j20 = j2[0];
	const CScalar <SType> & j3 = *(jj[1]->Get< CScalar <SType> >());
	const SComplex & j30 = j3[0];
	CScalar<SType>* j1 = NULL;
	j1 = CScalar<SType>::New(j00*j20*j30);
	j1->SetS(j0.S()|j2.S()|j3.S());
	return j1;

      }

      // if outgoing index is 2
      if (p_v->V()->id.back()==2){
	const CScalar <SType> & j0 = *(jj[1]->Get< CScalar <SType> >());
	const SComplex & j00 = j0[0];
	const CScalar <SType> & j1 = *(jj[2]->Get< CScalar <SType> >());
	const SComplex & j10 = j1[0];
	const CScalar <SType> & j3 = *(jj[0]->Get< CScalar <SType> >());
	const SComplex & j30 = j3[0];
	CScalar<SType>* j2 = NULL;
	j2 = CScalar<SType>::New(j00*j10*j30);
	j2->SetS(j0.S()|j1.S()|j3.S());
	return j2;

      }

      // if outgoing index is 3
      if (p_v->V()->id.back()==3){
	const CScalar <SType> & j0 = *(jj[0]->Get< CScalar <SType> >());
	const SComplex & j00 = j0[0];
	const CScalar <SType> & j1 = *(jj[1]->Get< CScalar <SType> >());
	const SComplex & j10 = j1[0];
	const CScalar <SType> & j2 = *(jj[2]->Get< CScalar <SType> >());
	const SComplex & j20 = j2[0];
	CScalar<SType>* j3 = NULL;
	j3 = CScalar<SType>::New(j00*j10*j20);
	j3->SetS(j0.S()|j1.S()|j2.S());
	return j3;

      }

      THROW(fatal_error, "Internal error in Lorentz calculator");
      return NULL;
    }

  };// end of class SSSS1_Calculator

  template class SSSS1_Calculator<double>;

}// end of namespace METOOLS


using namespace METOOLS;

DECLARE_GETTER(SSSS1_Calculator<double>,"DSSSS1",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,SSSS1_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new SSSS1_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    SSSS1_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"SSSS1 vertex"; }
