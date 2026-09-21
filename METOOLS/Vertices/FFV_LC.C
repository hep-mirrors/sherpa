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
  class FFV_Worker {
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

  };// end of class FFV_Worker

  template class FFV_Worker<double>;

  template <typename SType>
  class FFV_Calculator: public Lorentz_Calculator, 
			public FFV_Worker<SType> {
  public:
    
    typedef CSpinor<SType> CSpinorType;
    typedef CVec4<SType> CVec4Type;

    FFV_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "FFV"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      if (p_v->V()->id.back()==2) {
	CSpinorType *a(jj[1]->Get<CSpinorType>());
	CSpinorType *b(jj[0]->Get<CSpinorType>());
	bool cl(this->CalcLeft(*a,*b)), cr(this->CalcRight(*a,*b));
	if (!(cl || cr)) return NULL;
	CVec4Type *j(NULL);
	if (cl && cr) j=this->LorentzLeftRight(*a,*b);
	else if (cl) j=this->LorentzLeft(*a,*b);
	else if (cr) j=this->LorentzRight(*a,*b);
	return j;
      }
      const CSpinorType &a(*jj[p_v->V()->id.back()]->Get<CSpinorType>());
      const CVec4Type &b(*jj[1-p_v->V()->id.back()]->Get<CVec4Type>());
      bool cl(this->CalcLeft(a)), cr(this->CalcRight(a));
      if (!(cl || cr)) return NULL;
      CSpinorType *j(NULL);
      if (cl && cr) j=this->LorentzLeftRight(a,b);
      else if (cl) j=this->LorentzLeft(a,b);
      else if (cr) j=this->LorentzRight(a,b);
      return j;
    }

  };// end of class FFV_Calculator

  template class FFV_Calculator<double>;

  template <typename SType>
  class FFVL_Calculator: public Lorentz_Calculator, 
			 public FFV_Worker<SType> {
  public:

    typedef CSpinor<SType> CSpinorType;
    typedef CVec4<SType> CVec4Type;
    
    FFVL_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "FFVL"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      if (p_v->V()->id.back()==2) {
	CSpinorType *a(jj[1]->Get<CSpinorType>());
	CSpinorType *b(jj[0]->Get<CSpinorType>());
	if (!this->CalcLeft(*a,*b)) return NULL;
	return this->LorentzLeft(*a,*b);
      }
      const CSpinorType &a(*jj[p_v->V()->id.back()]->Get<CSpinorType>());
      const CVec4Type &b(*jj[1-p_v->V()->id.back()]->Get<CVec4Type>());
      if (!this->CalcLeft(a)) return NULL;
      return this->LorentzLeft(a,b);
    }

  };// end of class FFVL_Calculator

  template class FFVL_Calculator<double>;

  template <typename SType>
  class FFVR_Calculator: public Lorentz_Calculator, 
			 public FFV_Worker<SType> {
  public:
    
    typedef CSpinor<SType> CSpinorType;
    typedef CVec4<SType> CVec4Type;

    FFVR_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "FFVR"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      if (p_v->V()->id.back()==2) {
	CSpinorType *a(jj[1]->Get<CSpinorType>());
	CSpinorType *b(jj[0]->Get<CSpinorType>());
	if (!this->CalcRight(*a,*b)) return NULL;
	return this->LorentzRight(*a,*b);
      }
      const CSpinorType &a(*jj[p_v->V()->id.back()]->Get<CSpinorType>());
      const CVec4Type &b(*jj[1-p_v->V()->id.back()]->Get<CVec4Type>());
      if (!this->CalcRight(a)) return NULL;
      return this->LorentzRight(a,b);
    }

  };// end of class FFVR_Calculator

  template class FFVR_Calculator<double>;

}// end of namespace METOOLS

using namespace METOOLS;

DECLARE_GETTER(FFV_Calculator<double>,"DFFV",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,FFV_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new FFV_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    FFV_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"FFV vertex"; }

DECLARE_GETTER(FFVL_Calculator<double>,"DFFVL",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,FFVL_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new FFVL_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    FFVL_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"FFVL vertex"; }

DECLARE_GETTER(FFVR_Calculator<double>,"DFFVR",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,FFVR_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new FFVR_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    FFVR_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"FFVR vertex"; }
