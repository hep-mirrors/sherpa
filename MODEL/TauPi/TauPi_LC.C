#include "METOOLS/Explicit/Lorentz_Calculator.H"
#include "METOOLS/Currents/C_Spinor.H"
#include "METOOLS/Currents/C_Vector.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;

namespace METOOLS {

  template <typename SType>
  class TauPi_Worker {
  public:
    
    typedef std::complex<SType> SComplex;

    typedef CSpinor<SType> CSpinorType;
    typedef CVec4<SType> CVec4Type;

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

    inline bool CalcLeft(const CSpinorType &a,
			 const CSpinorType &b) 
    { return a.B()<0 ? a.On()&2 && b.On()&1 : a.On()&1 && b.On()&2; }
    inline bool CalcRight(const CSpinorType &a,
			  const CSpinorType &b) 
    { return a.B()<0 ? a.On()&1 && b.On()&2 : a.On()&2 && b.On()&1; }

    inline bool CalcLeft(const CSpinorType &a) 
    { return a.B()<0 ? a.On()&2 : a.On()&1; }
    inline bool CalcRight(const CSpinorType &a) 
    { return a.B()<0 ? a.On()&1 : a.On()&2; }

    CVec4<SType> *LorentzLeft(const CSpinorType &a,const CSpinorType &b)
    {
#ifdef DEBUG__BG
      msg_Debugging()<<"<> L "<<a<<"\n";
      msg_Debugging()<<"     "<<b<<"\n";
#endif
      SComplex j01(a[3]*b[1]), j02(a[2]*b[0]);
      SComplex j11(-a[2]*b[1]), j12(-a[3]*b[0]), j112(j11-j12);
      CVec4Type *j(CVec4Type::New(0.0,0.0,0.0,0.0,0,0,0,a.S()|b.S()));
      (*j)[0]=(j01+j02);
      (*j)[Spinor<SType>::R3()]=(j01-j02);
      (*j)[Spinor<SType>::R1()]=(j11+j12);
      (*j)[Spinor<SType>::R2()]=SComplex(j112.imag(),-j112.real());
      return j;
    }

    CVec4<SType> *LorentzRight(const CSpinorType &a,const CSpinorType &b)
    {
#ifdef DEBUG__BG
      msg_Debugging()<<"<> R "<<a<<"\n";
      msg_Debugging()<<"     "<<b<<"\n";
#endif
      SComplex j01(a[0]*b[2]), j02(a[1]*b[3]);
      SComplex j11(a[0]*b[3]), j12(a[1]*b[2]), j112(j11-j12);
      CVec4Type *j(CVec4Type::New(0.0,0.0,0.0,0.0,0,0,0,a.S()|b.S()));
      (*j)[0]=(j01+j02);
      (*j)[Spinor<SType>::R3()]=(j01-j02);
      (*j)[Spinor<SType>::R1()]=(j11+j12);
      (*j)[Spinor<SType>::R2()]=SComplex(j112.imag(),-j112.real());
      return j;
    }

    CVec4<SType> *LorentzLeftRight(const CSpinorType &a,const CSpinorType &b)
    {
#ifdef DEBUG__BG
      msg_Debugging()<<"<> LR "<<a<<"\n";
      msg_Debugging()<<"      "<<b<<"\n";
#endif
      SComplex l01(a[3]*b[1]), l02(a[2]*b[0]);
      SComplex l11(-a[2]*b[1]), l12(-a[3]*b[0]), l112(l11-l12);
      SComplex r01(a[0]*b[2]), r02(a[1]*b[3]);
      SComplex r11(a[0]*b[3]), r12(a[1]*b[2]), r112(r11-r12);
      CVec4Type *j(CVec4Type::New(0.0,0.0,0.0,0.0,0,0,0,a.S()|b.S()));
      (*j)[0]=(l01+l02)+(r01+r02);
      (*j)[Spinor<SType>::R3()]=(l01-l02)+(r01-r02);
      (*j)[Spinor<SType>::R1()]=(l11+l12)+(r11+r12);
      (*j)[Spinor<SType>::R2()]=
	SComplex(l112.imag(),-l112.real())+
	SComplex(r112.imag(),-r112.real());
      return j;
    }

    CSpinor<SType> *LorentzLeft(const CSpinorType &a,const CVec4Type &b)
    {
      switch (a.B()) {
      case -1: {
#ifdef DEBUG__BG
	msg_Debugging()<<"<|g L "<<a<<"\n";
	msg_Debugging()<<"      "<<b<<"\n";
#endif
	CSpinorType *j(CSpinorType::New(a.R(),a.B(),0,0,0,a.S()|b.S(),1));
	SComplex jp(PPlus(b)), jm(PMinus(b)), jt(PT(b)), jtc(PTC(b));
	(*j)[0]=(a[2]*jp+a[3]*jt);
	(*j)[1]=(a[2]*jtc+a[3]*jm);
	(*j)[3]=(*j)[2]=SComplex(0.0,0.0);
	return j;
      }
      case 1: {
#ifdef DEBUG__BG
	msg_Debugging()<<"g|> L "<<a<<"\n";
	msg_Debugging()<<"      "<<b<<"\n";
#endif
	CSpinorType *j(CSpinorType::New(a.R(),a.B(),0,0,0,a.S()|b.S(),2));
	SComplex jp(PPlus(b)), jm(PMinus(b)), jt(PT(b)), jtc(PTC(b));
	(*j)[1]=(*j)[0]=SComplex(0.0,0.0);
	(*j)[2]=(a[0]*jp+a[1]*jtc);
	(*j)[3]=(a[0]*jt+a[1]*jm);
	return j;
      }
      }
      return NULL;
    }

    CSpinor<SType> *LorentzRight(const CSpinorType &a,const CVec4Type &b)
    {
      switch (a.B()) {
      case -1: {
#ifdef DEBUG__BG
	msg_Debugging()<<"<|g R "<<a<<"\n";
	msg_Debugging()<<"      "<<b<<"\n";
#endif
	CSpinorType *j(CSpinorType::New(a.R(),a.B(),0,0,0,a.S()|b.S(),2));
	SComplex jp(PPlus(b)), jm(PMinus(b)), jt(PT(b)), jtc(PTC(b));
	(*j)[1]=(*j)[0]=SComplex(0.0,0.0);
	(*j)[2]=(a[0]*jm-a[1]*jt);
	(*j)[3]=(-a[0]*jtc+a[1]*jp);
	return j;
      }
      case 1: {
#ifdef DEBUG__BG
	msg_Debugging()<<"g|> R "<<a<<"\n";
	msg_Debugging()<<"      "<<b<<"\n";
#endif
	CSpinorType *j(CSpinorType::New(a.R(),a.B(),0,0,0,a.S()|b.S(),1));
	SComplex jp(PPlus(b)), jm(PMinus(b)), jt(PT(b)), jtc(PTC(b));
	(*j)[0]=(a[2]*jm-a[3]*jtc);
	(*j)[1]=(-a[2]*jt+a[3]*jp);
	(*j)[3]=(*j)[2]=SComplex(0.0,0.0);
	return j;
      }
      }
      return NULL;
    }
    
    CSpinor<SType> *LorentzLeftRight(const CSpinorType &a,const CVec4Type &b)
    {
      switch (a.B()) {
      case -1: {
#ifdef DEBUG__BG
	msg_Debugging()<<"<|g LR "<<a<<"\n";
	msg_Debugging()<<"       "<<b<<"\n";
#endif
	CSpinorType *j(CSpinorType::New(a.R(),a.B(),0,0,0,a.S()|b.S(),3));
	SComplex jp(PPlus(b)), jm(PMinus(b)), jt(PT(b)), jtc(PTC(b));
	(*j)[0]=(a[2]*jp+a[3]*jt);
	(*j)[1]=(a[2]*jtc+a[3]*jm);
	(*j)[2]=(a[0]*jm-a[1]*jt);
	(*j)[3]=(-a[0]*jtc+a[1]*jp);
	return j;
      }
      case 1: {
#ifdef DEBUG__BG
	msg_Debugging()<<"g|> LR "<<a<<"\n";
	msg_Debugging()<<"       "<<b<<"\n";
#endif
	CSpinorType *j(CSpinorType::New(a.R(),a.B(),0,0,0,a.S()|b.S(),3));
	SComplex jp(PPlus(b)), jm(PMinus(b)), jt(PT(b)), jtc(PTC(b));
	(*j)[0]=(a[2]*jm-a[3]*jtc);
	(*j)[1]=(-a[2]*jt+a[3]*jp);
	(*j)[2]=(a[0]*jp+a[1]*jtc);
	(*j)[3]=(a[0]*jt+a[1]*jm);
	return j;
      }
      }
      return NULL;
    }

  };// end of class TauPi_Worker

  template class TauPi_Worker<double>;

  template <typename SType>
  class TauPi_Calculator: public Lorentz_Calculator, 
			public TauPi_Worker<SType> {
  public:
    
    typedef CSpinor<SType> CSpinorType;
    typedef CVec4<SType> CVec4Type;

    TauPi_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "TauPi"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      if (p_v->V()->id.back()==2) THROW(not_implemented,"Implement me!");
      const CSpinorType &a(*jj[p_v->V()->id.back()]->Get<CSpinorType>());
      CVec4Type b(p_v->J(1-p_v->V()->id.back())->P());
      bool cl(this->CalcLeft(a)), cr(this->CalcRight(a));
      if (!(cl || cr)) return NULL;
      CSpinorType *j(NULL);
      if (cl && cr) j=this->LorentzLeftRight(a,b);
      else if (cl) j=this->LorentzLeft(a,b);
      else if (cr) j=this->LorentzRight(a,b);
      return j;
    }

  };// end of class TauPi_Calculator

  template class TauPi_Calculator<double>;

  template <typename SType>
  class TauPiL_Calculator: public Lorentz_Calculator, 
			 public TauPi_Worker<SType> {
  public:

    typedef CSpinor<SType> CSpinorType;
    typedef CVec4<SType> CVec4Type;
    
    TauPiL_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "TauPiL"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      if (p_v->V()->id.back()==2) THROW(not_implemented,"Implement me!");
      const CSpinorType &a(*jj[p_v->V()->id.back()]->Get<CSpinorType>());
      CVec4Type b(p_v->J(1-p_v->V()->id.back())->P());
      if (!this->CalcLeft(a)) return NULL;
      return this->LorentzLeft(a,b);
    }

  };// end of class TauPiL_Calculator

  template class TauPiL_Calculator<double>;

  template <typename SType>
  class TauPiR_Calculator: public Lorentz_Calculator, 
			 public TauPi_Worker<SType> {
  public:
    
    typedef CSpinor<SType> CSpinorType;
    typedef CVec4<SType> CVec4Type;

    TauPiR_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "TauPiR"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      if (p_v->V()->id.back()==2) THROW(not_implemented,"Implement me!");
      const CSpinorType &a(*jj[p_v->V()->id.back()]->Get<CSpinorType>());
      CVec4Type b(p_v->J(1-p_v->V()->id.back())->P());
      if (!this->CalcRight(a)) return NULL;
      return this->LorentzRight(a,b);
    }

  };// end of class TauPiR_Calculator

  template class TauPiR_Calculator<double>;

}// end of namespace METOOLS

using namespace METOOLS;

DECLARE_GETTER(TauPi_Calculator<double>,"DTauPi",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,TauPi_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new TauPi_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    TauPi_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"TauPi vertex"; }

DECLARE_GETTER(TauPiL_Calculator<double>,"DTauPiL",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,TauPiL_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new TauPiL_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    TauPiL_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"TauPiL vertex"; }

DECLARE_GETTER(TauPiR_Calculator<double>,"DTauPiR",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,TauPiR_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new TauPiR_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    TauPiR_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"TauPiR vertex"; }
