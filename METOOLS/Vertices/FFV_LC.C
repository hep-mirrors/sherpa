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

    CVec4<SType> *LorentzMagnetic
    (const CSpinorType &a,const CSpinorType &b,const Vec4D pab)
    {
#ifdef DEBUG__BG
      msg_Debugging()<<"<> sig "<<a<<"\n";
      msg_Debugging()<<"       "<<b<<"\n";
      msg_Debugging()<<"       "<<pab<<"\n";
#endif
      const SComplex I(0,1);
      SComplex p0(pab[0]), p3(-pab[Spinor<SType>::R3()]);
      SComplex p1(-pab[Spinor<SType>::R1()]), p2(-pab[Spinor<SType>::R2()]);
      SComplex l12(a[0]*b[0]-a[1]*b[1]+a[2]*b[2]-a[3]*b[3]);
      SComplex l13(I*(a[0]*b[1]-a[1]*b[0]+a[2]*b[3]-a[3]*b[2]));
      SComplex l23(a[0]*b[1]+a[1]*b[0]+a[2]*b[3]+a[3]*b[2]);
      SComplex l01(I*(-a[0]*b[1]-a[1]*b[0]+a[2]*b[3]+a[3]*b[2]));
      SComplex l02(-a[0]*b[1]+a[1]*b[0]+a[2]*b[3]-a[3]*b[2]);
      SComplex l03(I*(-a[0]*b[0]+a[1]*b[1]+a[2]*b[2]-a[3]*b[3]));
      CVec4Type *j(CVec4Type::New(0.0,0.0,0.0,0.0,0,0,0,a.S()|b.S()));
      (*j)[0]=l01*p1+l02*p2+l03*p3;
      (*j)[Spinor<SType>::R1()]=-l01*p0+l12*p2+l13*p3;
      (*j)[Spinor<SType>::R2()]=-l02*p0-l12*p1+l23*p3;
      (*j)[Spinor<SType>::R3()]=-l03*p0-l13*p1-l23*p2;
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

    CSpinor<SType> *LorentzMagnetic
    (const CSpinorType &a,const CVec4Type &b,const Vec4D &pb)
    {
      switch (a.B()) {
      case -1: {
#ifdef DEBUG__BG
	msg_Debugging()<<"<|sig  "<<a<<"\n";
	msg_Debugging()<<"       "<<b<<"\n";
	msg_Debugging()<<"       "<<pb<<"\n";
#endif
	CSpinorType *j(CSpinorType::New(a.R(),a.B(),0,0,0,a.S()|b.S(),3));
	SComplex jp(PPlus(b)), jm(PMinus(b)), jt(PT(b)), jtc(PTC(b));
	SComplex j0(a[2]*jp+a[3]*jt);
	SComplex j1(a[2]*jtc+a[3]*jm);
	SComplex j2(a[0]*jm-a[1]*jt);
	SComplex j3(-a[0]*jtc+a[1]*jp);
	SComplex pe(pb*b), I(0,1);
	SComplex pp(PPlus(pb)), pm(PMinus(pb)), pt(PT(pb)), ptc(PTC(pb));
	(*j)[0]=I*(j2*pp+j3*pt-pe*a[0]);
	(*j)[1]=I*(j2*ptc+j3*pm-pe*a[1]);
	(*j)[2]=I*(j0*pm-j1*pt-pe*a[2]);
	(*j)[3]=I*(-j0*ptc+j1*pp-pe*a[3]);
#ifdef CHECK__Sigma
	{
        const SComplex I(0,1);
	CSpinorType *k(CSpinorType::New(a.R(),a.B(),0,0,0,a.S()|b.S(),3));
	SComplex ep01(-I*(b[0]*pb[Spinor<SType>::R1()]-b[Spinor<SType>::R1()]*pb[0]));
	SComplex ep02(-(b[0]*pb[Spinor<SType>::R2()]-b[Spinor<SType>::R2()]*pb[0]));
	SComplex ep03(-I*(b[0]*pb[Spinor<SType>::R3()]-b[Spinor<SType>::R3()]*pb[0]));
	SComplex ep12(b[Spinor<SType>::R1()]*pb[Spinor<SType>::R2()]
		      -b[Spinor<SType>::R2()]*pb[Spinor<SType>::R1()]);
	SComplex ep13(I*(b[Spinor<SType>::R1()]*pb[Spinor<SType>::R3()]
			 -b[Spinor<SType>::R3()]*pb[Spinor<SType>::R1()]));
	SComplex ep23(b[Spinor<SType>::R2()]*pb[Spinor<SType>::R3()]
		      -b[Spinor<SType>::R3()]*pb[Spinor<SType>::R2()]);
	(*k)[0]=a[0]*(ep12-ep03)+a[1]*(-ep13+ep23-ep01+ep02);
	(*k)[1]=a[1]*(-ep12+ep03)+a[0]*(ep13+ep23-ep01-ep02);
	(*k)[2]=a[2]*(ep12+ep03)+a[3]*(-ep13+ep23+ep01-ep02);
	(*k)[3]=a[3]*(-ep12-ep03)+a[2]*(ep13+ep23+ep01+ep02);
	DEBUG_VAR(*j);
	DEBUG_VAR(*k);
	}
#endif
	return j;
      }
      case 1: {
#ifdef DEBUG__BG
	msg_Debugging()<<"sig|>  "<<a<<"\n";
	msg_Debugging()<<"       "<<b<<"\n";
	msg_Debugging()<<"       "<<pb<<"\n";
#endif
	CSpinorType *j(CSpinorType::New(a.R(),a.B(),0,0,0,a.S()|b.S(),3));
	SComplex pp(PPlus(pb)), pm(PMinus(pb)), pt(PT(pb)), ptc(PTC(pb));
	SComplex j0(a[2]*pm-a[3]*ptc);
	SComplex j1(-a[2]*pt+a[3]*pp);
	SComplex j2(a[0]*pp+a[1]*ptc);
	SComplex j3(a[0]*pt+a[1]*pm);
	SComplex pe(pb*b), I(0,1);
	SComplex jp(PPlus(b)), jm(PMinus(b)), jt(PT(b)), jtc(PTC(b));
	(*j)[0]=I*(j2*jm-j3*jtc-pe*a[0]);
	(*j)[1]=I*(-j2*jt+j3*jp-pe*a[1]);
	(*j)[2]=I*(j0*jp+j1*jtc-pe*a[2]);
	(*j)[3]=I*(j0*jt+j1*jm-pe*a[3]);
#ifdef CHECK__Sigma
	{
        const SComplex I(0,1);
	CSpinorType *k(CSpinorType::New(a.R(),a.B(),0,0,0,a.S()|b.S(),3));
	SComplex ep01(-I*(b[0]*pb[Spinor<SType>::R1()]-b[Spinor<SType>::R1()]*pb[0]));
	SComplex ep02(-(b[0]*pb[Spinor<SType>::R2()]-b[Spinor<SType>::R2()]*pb[0]));
	SComplex ep03(-I*(b[0]*pb[Spinor<SType>::R3()]-b[Spinor<SType>::R3()]*pb[0]));
	SComplex ep12(b[Spinor<SType>::R1()]*pb[Spinor<SType>::R2()]
		      -b[Spinor<SType>::R2()]*pb[Spinor<SType>::R1()]);
	SComplex ep13(I*(b[Spinor<SType>::R1()]*pb[Spinor<SType>::R3()]
			 -b[Spinor<SType>::R3()]*pb[Spinor<SType>::R1()]));
	SComplex ep23(b[Spinor<SType>::R2()]*pb[Spinor<SType>::R3()]
		      -b[Spinor<SType>::R3()]*pb[Spinor<SType>::R2()]);
	(*k)[0]=a[0]*(ep12-ep03)+a[1]*(ep13+ep23-ep01-ep02);
	(*k)[1]=a[1]*(-ep12+ep03)+a[0]*(-ep13+ep23-ep01+ep02);
	(*k)[2]=a[2]*(ep12+ep03)+a[3]*(ep13+ep23+ep01+ep02);
	(*k)[3]=a[3]*(-ep12-ep03)+a[2]*(-ep13+ep23+ep01-ep02);
	DEBUG_VAR(*j);
	DEBUG_VAR(*k);
	}
#endif
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

  template <typename SType>
  class FFVM_Calculator: public Lorentz_Calculator,
			 public FFV_Worker<SType> {
  public:

    typedef CSpinor<SType> CSpinorType;
    typedef CVec4<SType> CVec4Type;

    FFVM_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "FFVM"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      if (p_v->V()->id.back()==2) {
	CSpinorType *a(jj[1]->Get<CSpinorType>());
	CSpinorType *b(jj[0]->Get<CSpinorType>());
	Vec4D pab(p_v->J(0)->P()+p_v->J(1)->P());
	return this->LorentzMagnetic(*a,*b,pab);
      }
      const CSpinorType &a(*jj[p_v->V()->id.back()]->Get<CSpinorType>());
      const CVec4Type &b(*jj[1-p_v->V()->id.back()]->Get<CVec4Type>());
      Vec4D pb(p_v->J(1)->P());
      return this->LorentzMagnetic(a,b,pb);
    }

  };// end of class FFVM_Calculator

  template class FFVM_Calculator<double>;

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

DECLARE_GETTER(FFVM_Calculator<double>,"DFFVM",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,FFVM_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new FFVM_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    FFVM_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"FFVM vertex"; }
