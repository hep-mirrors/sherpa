#include "METOOLS/Explicit/Lorentz_Calculator.H"
#include "METOOLS/Currents/C_Spinor.H"
#include "METOOLS/Currents/C_Scalar.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;

namespace METOOLS {

  template <typename SType>
  class FFS_Worker {
  public:
    
    typedef std::complex<SType> SComplex;

    typedef CSpinor<SType> CSpinorType;
    typedef CScalar<SType> CScalarType;

    inline SComplex PPlus(const CVec4<SType> &p) const  
    { return p[0]+p[ATOOLS::Spinor<SType>::R3()]; }
    inline SComplex PMinus(const CVec4<SType> &p) const 
    { return p[0]-p[ATOOLS::Spinor<SType>::R3()]; }

    inline SComplex PT(const CVec4<SType> &p) const  
    { return p[ATOOLS::Spinor<SType>::R1()]+
	SComplex(0.0,1.0)*p[ATOOLS::Spinor<SType>::R2()]; }
    inline SComplex PTC(const CVec4<SType> &p) const  
    { return p[ATOOLS::Spinor<SType>::R1()]-
	SComplex(0.0,1.0)*p[ATOOLS::Spinor<SType>::R2()]; }

    CScalarType LorentzLeft(const CSpinorType &a,const CSpinorType &b)
    {
#ifdef DEBUG__BG
      msg_Debugging()<<"<> L "<<a<<"\n";
      msg_Debugging()<<"     "<<b<<"\n";
#endif
      return CScalarType(a[2]*b[2]+a[3]*b[3],0,a.S()|b.S());
    }

    CScalarType LorentzRight(const CSpinorType &a,const CSpinorType &b)
    {
#ifdef DEBUG__BG
      msg_Debugging()<<"<> R "<<a<<"\n";
      msg_Debugging()<<"     "<<b<<"\n";
#endif
      return CScalarType(a[0]*b[0]+a[1]*b[1],0,a.S()|b.S());
    }

    CSpinorType LorentzLeft(const CSpinorType &a,const CScalarType &b)
    {
      switch (a.B()) {
      case -1: {
	CSpinorType j(a.R(),a.B(),0,0,a.H()|b.H(),a.S()|b.S(),1);
#ifdef DEBUG__BG
	msg_Debugging()<<"<| R "<<a<<"\n";
	msg_Debugging()<<"     "<<b<<"\n";
#endif
	j[0]=a[0]*b[0];
	j[1]=a[1]*b[0];
	j[3]=j[2]=SComplex(0.0,0.0);
	return j;
      }
      case 1: {
	CSpinorType j(a.R(),a.B(),0,0,a.H()|b.H(),a.S()|b.S(),2);
#ifdef DEBUG__BG
	msg_Debugging()<<"|> L "<<a<<"\n";
	msg_Debugging()<<"     "<<b<<"\n";
#endif
	j[1]=j[0]=SComplex(0.0,0.0);
	j[2]=a[2]*b[0];
	j[3]=a[3]*b[0];
	return j;
      }
      }
      return CSpinorType();
    }

    CSpinorType LorentzRight(const CSpinorType &a,const CScalarType &b)
    {
      switch (a.B()) {
      case -1: {
	CSpinorType j(a.R(),a.B(),0,0,a.H()|b.H(),a.S()|b.S(),2);
#ifdef DEBUG__BG
	msg_Debugging()<<"<| L "<<a<<"\n";
	msg_Debugging()<<"     "<<b<<"\n";
#endif
	j[2]=a[2]*b[0];
	j[3]=a[3]*b[0];
	j[1]=j[0]=SComplex(0.0,0.0);
	return j;
      }
      case 1: {
	CSpinorType j(a.R(),a.B(),0,0,a.H()|b.H(),a.S()|b.S(),1);
#ifdef DEBUG__BG
	msg_Debugging()<<"|> R "<<a<<"\n";
	msg_Debugging()<<"     "<<b<<"\n";
#endif
	j[3]=j[2]=SComplex(0.0,0.0);
	j[0]=a[0]*b[0];
	j[1]=a[1]*b[0];
	return j;
      }
      }
      return CSpinorType();
    }

  };// end of class FFS_Worker

  template class FFS_Worker<double>;

  template <typename SType>
  class FFS_Calculator: public Lorentz_Calculator, 
			public FFS_Worker<SType> {
  public:
    
    typedef CSpinor<SType> CSpinorType;
    typedef CScalar<SType> CScalarType;

    FFS_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "FFS"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      if (p_v->V()->id.back()==2) {
	CSpinorType *a(jj[1]->Get<CSpinorType>());
	CSpinorType *b(jj[0]->Get<CSpinorType>());
	return CScalarType::New
	  (this->LorentzLeft(*a,*b)+this->LorentzRight(*a,*b));
      }
      const CSpinorType &a(*jj[p_v->V()->id.back()]->Get<CSpinorType>());
      const CScalarType &b(*jj[1-p_v->V()->id.back()]->Get<CScalarType>());
      return CSpinorType::New
	(this->LorentzLeft(a,b)+this->LorentzRight(a,b));
    }

  };// end of class FFS_Calculator

  template class FFS_Calculator<double>;

  template <typename SType>
  class FFSL_Calculator: public Lorentz_Calculator, 
			 public FFS_Worker<SType> {
  public:

    typedef CSpinor<SType> CSpinorType;
    typedef CScalar<SType> CScalarType;
    
    FFSL_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "FFSL"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      if (p_v->V()->id.back()==2) {
	CSpinorType *a(jj[1]->Get<CSpinorType>());
	CSpinorType *b(jj[0]->Get<CSpinorType>());
	return CScalarType::New(this->LorentzLeft(*a,*b));
      }
      const CSpinorType &a(*jj[p_v->V()->id.back()]->Get<CSpinorType>());
      const CScalarType &b(*jj[1-p_v->V()->id.back()]->Get<CScalarType>());
      return CSpinorType::New(this->LorentzLeft(a,b));
    }

  };// end of class FFSL_Calculator

  template class FFSL_Calculator<double>;

  template <typename SType>
  class FFSR_Calculator: public Lorentz_Calculator, 
			 public FFS_Worker<SType> {
  public:
    
    typedef CSpinor<SType> CSpinorType;
    typedef CScalar<SType> CScalarType;

    FFSR_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "FFSR"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      if (p_v->V()->id.back()==2) {
	CSpinorType *a(jj[1]->Get<CSpinorType>());
	CSpinorType *b(jj[0]->Get<CSpinorType>());
	return CScalarType::New(this->LorentzRight(*a,*b));
      }
      const CSpinorType &a(*jj[p_v->V()->id.back()]->Get<CSpinorType>());
      const CScalarType &b(*jj[1-p_v->V()->id.back()]->Get<CScalarType>());
      return CSpinorType::New(this->LorentzRight(a,b));
    }

  };// end of class FFSR_Calculator

  template class FFSR_Calculator<double>;

}// end of namespace METOOLS

using namespace METOOLS;

DECLARE_GETTER(FFS_Calculator<double>,"DFFS",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,FFS_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new FFS_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    FFS_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"FFS vertex"; }

DECLARE_GETTER(FFSL_Calculator<double>,"DFFSL",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,FFSL_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new FFSL_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    FFSL_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"FFSL vertex"; }

DECLARE_GETTER(FFSR_Calculator<double>,"DFFSR",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,FFSR_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new FFSR_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    FFSR_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"FFSR vertex"; }
