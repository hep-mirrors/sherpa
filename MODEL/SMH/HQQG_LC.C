#include "METOOLS/Explicit/Lorentz_Calculator.H"
#include "METOOLS/Currents/C_Spinor.H"
#include "METOOLS/Currents/C_Vector.H"
#include "METOOLS/Currents/C_Scalar.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "MODEL/SMH/qqgH.h"

using namespace ATOOLS;

namespace METOOLS {

  template <typename SType>
  class HQQG_Worker {
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

  };// end of class HQQG_Worker

  template class HQQG_Worker<double>;

  template <typename SType>
  class HQQG_Calculator: public Lorentz_Calculator, 
			public HQQG_Worker<SType> {
  private:

    int m_n[3];
    double m_mh;
    std::vector<double> m_mq;

  public:
    
    typedef CSpinor<SType> CSpinorType;
    typedef CVec4<SType> CVec4Type;
    typedef CScalar<SType> CScalarType;

    HQQG_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key)
    {
      m_mh=Flavour(kf_h0).Mass();
      for (size_t i(0);i<6;++i)
	if (Flavour((kf_code)i).Yuk())
	  m_mq.push_back(Flavour((kf_code)i).Yuk());
      if (p_v->V()->id.back()==3) { m_n[0]=0; m_n[1]=2;/*q*/ m_n[2]=1;/*Q*/ }
      if (p_v->V()->id.back()==2) { m_n[0]=1; m_n[1]=2;/*Q*/ m_n[2]=0;/*g*/ }
      if (p_v->V()->id.back()==1) { m_n[0]=2; m_n[1]=0;/*q*/ m_n[2]=1;/*g*/ }
    }

    std::string Label() const { return "HQQG"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      if (p_v->V()->id.back()==0) THROW(fatal_error,"Invalid call");
      if (p_v->V()->id.back()==3) {
	const CSpinorType &a(*jj[m_n[1]]->Get<CSpinorType>());
	const CSpinorType &b(*jj[m_n[2]]->Get<CSpinorType>());
	const CScalarType &e(*jj[m_n[0]]->template Get<CScalarType>());
	bool cl(this->CalcLeft(a,b)), cr(this->CalcRight(a,b));
	if (!(cl || cr)) return NULL;
	Vec4D pa(p_v->J(m_n[1])->P()), pb(p_v->J(m_n[2])->P());
	Vec4D pe(p_v->J(m_n[0])->P()), pc(-pa-pb-pe);
	Complex F[2]={0.,0.};
#ifdef USING__HEFT
	double sab((pa+pb).Abs2());
	F[0]=Complex(0.,2./3.)*1./sab;
	F[1]=Complex(0.,2./3.)*1./sab;
#else
	double sac((pa+pc).Abs2()), sbc((pb+pc).Abs2());
	for (size_t i(0);i<m_mq.size();++i) {
	  F[0]+=F1qqg(m_mq[i],m_mh,sac,sbc);
	  F[1]+=F2qqg(m_mq[i],m_mh,sac,sbc);
	}
#endif
	CVec4Type *j1(NULL);
	if (cl && cr) j1=this->LorentzLeftRight(a,b);
	else if (cl) j1=this->LorentzLeft(a,b);
	else if (cr) j1=this->LorentzRight(a,b);
	CVec4Type j2(e[0]*(*j1*pc)*(CVec4Type(pa)*F[0]+CVec4Type(pb)*F[1]));
	j1->Multiply(e[0]*(-(pc*pa)*F[0]-(pc*pb)*F[1]));
	j1->Add(&j2);
	return j1;
      }
      const CSpinorType &a(*jj[m_n[1]]->Get<CSpinorType>());
      const CVec4Type &c(*jj[m_n[2]]->Get<CVec4Type>());
      const CScalarType &e(*jj[m_n[0]]->template Get<CScalarType>());
      bool cl(this->CalcLeft(a)), cr(this->CalcRight(a));
      if (!(cl || cr)) return NULL;
      Vec4D pa(p_v->J(m_n[1])->P()), pc(p_v->J(m_n[2])->P());
      Vec4D pe(p_v->J(m_n[0])->P()), pb(-pa-pc-pe);
      CSpinorType *j1(NULL), *j2(NULL);
      if (cl && cr) {
	j1=this->LorentzLeftRight(a,CVec4Type(pc));
	j2=this->LorentzLeftRight(a,c);
      }
      else if (cl) {
	j1=this->LorentzLeft(a,CVec4Type(pc));
	j2=this->LorentzLeft(a,c);
      }
      else if (cr) {
	j1=this->LorentzRight(a,CVec4Type(pc));
	j2=this->LorentzRight(a,c);
      }
      Complex F[2]={0.,0.};
#ifdef USING__HEFT
      double sab((pa+pb).Abs2());
      F[0]=Complex(0.,2./3.)*1./sab;
      F[1]=Complex(0.,2./3.)*1./sab;
#else
      double sac((pa+pc).Abs2()), sbc((pb+pc).Abs2());
      for (size_t i(0);i<m_mq.size();++i) {
	F[0]+=F1qqg(m_mq[i],m_mh,sac,sbc);
	F[1]+=F2qqg(m_mq[i],m_mh,sac,sbc);
      }
#endif
      j1->Multiply(e[0]*((c*pa)*F[0]+(c*pb)*F[1]));
      j2->Multiply(e[0]*(-(pc*pa)*F[0]-(pc*pb)*F[1]));
      j1->Add(j2);
      j2->Delete();
      return j1;
    }

  };// end of class HQQG_Calculator

}// end of namespace METOOLS

using namespace METOOLS;

DECLARE_GETTER(HQQG_Calculator<double>,"DHQQG",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,HQQG_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new HQQG_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    HQQG_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"HQQG vertex"; }
