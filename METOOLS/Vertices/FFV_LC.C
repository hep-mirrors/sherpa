#include "METOOLS/Explicit/Lorentz_Calculator.H"

#include "METOOLS/Explicit/Vertex.H"
#include "METOOLS/Explicit/Dipole_Kinematics.H"
#include "METOOLS/Explicit/Dipole_Color.H"
#include "METOOLS/Currents/C_Spinor.H"

namespace METOOLS {

  template <typename SType>
  class FFV_Calculator: public Lorentz_Calculator {
  public:

    typedef std::complex<SType> SComplex;

    typedef CSpinor<SType> CSpinorType;
    typedef std::vector<CSpinorType*> CSpinorType_Vector;
    typedef std::vector<CSpinorType_Vector> CSpinorType_Matrix;

    typedef CVec4<SType> CVec4Type;
    typedef std::vector<CVec4Type*> CVec4Type_Vector;

  private:

    SComplex m_cpll, m_cplr;

    int m_dir, m_cl, m_cr;

    void Evaluate(const CSpinorType &a,const CSpinorType &b,const size_t &h);
    void Evaluate(const CSpinorType &a,const CVec4Type &b,const size_t &h);

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

    CVec4Type *LorentzLeft(CSpinorType a,CSpinorType b);
    CVec4Type *LorentzRight(CSpinorType a,CSpinorType b);
    CVec4Type *LorentzLeftRight(CSpinorType a,CSpinorType b);

    CSpinorType *LorentzLeft(const CSpinorType &a,const CVec4Type &b);
    CSpinorType *LorentzRight(const CSpinorType &a,const CVec4Type &b);
    CSpinorType *LorentzLeftRight(const CSpinorType &a,const CVec4Type &b);

    CVec4Type *GetPol(const ATOOLS::Vec4D &p,
		      const ATOOLS::Vec4D &q,const int mode);
    CSpinorType *GetPol(const ATOOLS::Vec4D &q,
			const double &q2,const int mode);

  public:
    
    FFV_Calculator(const Vertex_Key &key);

    std::string Label() const;
    
    void Evaluate();

    void ConstructFFSDipole();
    void ConstructFVSDipole();

    void ConstructFVIDipole();

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

  };// end of class FFV_Calculator

}// end of namespace METOOLS

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"

#define ZERO SComplex(0.0,0.0)

using namespace METOOLS;
using namespace ATOOLS;

template <typename SType>
FFV_Calculator<SType>::FFV_Calculator(const Vertex_Key &key): 
  Lorentz_Calculator(key),  
  m_dir(key.FlB().IsFermion()?
	(key.FlA().IsFermion()?0:2):1) 
{
  bool nuc(p_v->Info() && p_v->Info()->Mode()&2);
  if (m_dir!=0 && key.p_c->Flav().IsAnti()) {
    m_cpll=SComplex(-(nuc?1.0:p_v->Coupling(0))*p_cc->Coupling());
    m_cplr=SComplex(-(nuc?1.0:p_v->Coupling(1))*p_cc->Coupling());
  }
  else {
    m_cpll=SComplex((nuc?1.0:p_v->Coupling(1))*p_cc->Coupling());
    m_cplr=SComplex((nuc?1.0:p_v->Coupling(0))*p_cc->Coupling());
  }
  m_cl=m_cpll!=SComplex(0.0,0.0);
  m_cr=m_cplr!=SComplex(0.0,0.0);
}

template <typename SType> CVec4<SType> *
FFV_Calculator<SType>::GetPol
(const Vec4D &p,const Vec4D &q,const int mode)
{
  static double sqrttwo(sqrt(2.0));
  Vec3D qxp(cross(Vec3D(q),Vec3D(p))), qxpxp(p*q[0]-q*p[0]);
  CVec4Type e1(Vec4D(0.0,qxp/(sqrttwo*qxp.Abs())));
  CVec4Type e2(Vec4D(0.0,qxpxp/(sqrttwo*qxpxp.Abs())));
  return CVec4Type::New(e1+SComplex(0.0,mode?1.0:-1.0)*e2); 
}

template <typename SType>
void FFV_Calculator<SType>::ConstructFFSDipole()
{
  const CSpinorType_Vector *ca(p_v->JA()->J().front().Get<CSpinorType>());
  const CSpinorType_Vector *cb(p_v->JB()->J().front().Get<CSpinorType>());
  const CObject_Vector *cc(&p_v->Kin()->JK()->J().front());
  if (!p_cc->Evaluate(ca->front(),cb->front(),cc->front())) return;
  Vec4D p(p_v->JC()->P()), q;
  double A(0.0), B(0.0), t(0.0);
  if (p_v->Kin()->Type()==0) {
    double zi(p_v->Kin()->Z()), zj(1.0-zi), rv(1.0);
    Vec4D pi(p_v->Kin()->PI()), pj(p_v->Kin()->PJ());
    double pij2((pi+pj).Abs2()), mt2(0.0), zim(zi), zjm(zj);
    if (p_v->Kin()->Massive()) {
      double mi(p_v->Kin()->JI()->Flav().Mass());
      double mj(p_v->Kin()->JJ()->Flav().Mass());
      double mk(p_v->Kin()->JK()->Flav().Mass());
      double y(p_v->Kin()->Y()), Q2(p_v->Kin()->Q2());
      double mi2(mi*mi), mj2(mj*mj), mk2(mk*mk), eps(Q2-mi2-mj2-mk2);
      double viji(sqrt(sqr(eps*y)-sqr(2.0*mi*mj))/(eps*y+2.0*mi2));
      rv=sqrt(sqr(2.0*mk2+eps*(1.0-y))-4.0*mk2*Q2)/(eps*(1.0-y));
      double zc(0.5*(2.0*mi2+eps*y)/(mi2+mj2+eps*y));
      double zm(zc*(1.0-viji*rv)), zp(zc*(1.0+viji*rv));
      mt2=2.0*p_v->Info()->Kappa()*(zp*zm-mi2/pij2);
      zim-=0.5*(1.0-rv);
      zjm-=0.5*(1.0-rv);
    }
    q=zim*pi-zjm*pj;
    A=1.0-mt2;
    B=-4.0*q.Abs2()/pij2;
    t=rv*pij2;
  }
  else if (p_v->Kin()->Type()==2) {
    double zi(p_v->Kin()->Z()), zj(1.0-zi);
    Vec4D pi(p_v->Kin()->PI()), pj(p_v->Kin()->PJ());
    double pij2((pi+pj).Abs2());
    q=zi*pi-zj*pj;
    A=1.0;
    B=-4.0*q.Abs2()/pij2;
    t=pij2*(1.0-p_v->Kin()->Y());
  }
  else if (p_v->Kin()->Type()==1) {
    double x(p_v->Kin()->Z()), ui(p_v->Kin()->Y());
    Vec4D pi(p_v->Kin()->PJ()), pk(p_v->Kin()->PK());
    A=x;
    B=-4.0*(1.0-x)/x;
    q=pi/ui-pk/(1.0-ui);
    t=-2.0*(pi*p_v->Kin()->PI())*x;
  }
  else {
    double x(p_v->Kin()->Z()), vi(p_v->Kin()->Y());
    Vec4D pi(p_v->Kin()->PJ()), pk(-p_v->Kin()->PK());
    A=x;
    B=-4.0*(1.0-x)/x;
    q=pi-vi*pk;
    t=-2.0*(pi*p_v->Kin()->PI())*x;
  }
  double At(A-B/2.0);
  p_v->Kin()->SetPhase(1.0/(2.0*A/B-1.0),0);
  p_v->Kin()->SetPhase(1.0/(2.0*A/B-1.0),1);
#ifdef CHECK__pols
  DEBUG_VAR(p_v->Kin()->Type()<<" "<<A<<" "<<B<<" "<<q<<" "<<q.Abs2());
  CVec4<SType> *pols[2]={GetPol(p,q,0),GetPol(p,q,1)};
  CVec4<SType> e1(SType(sqrt(0.5))*(*pols[0]+*pols[1]));
  CVec4<SType> e2(SComplex(0.0,sqrt(0.5))*(*pols[0]-*pols[1]));
  SComplex mat[4][4], ref[4][4];
  Vec4D n(p[0],Vec3D(-p)), qt(p*q[0]-q*p[0]);
  for (size_t i(0);i<4;++i)
    for (size_t j(0);j<4;++j) {
      mat[i][j]=0.0;
      for (size_t k(0);k<2;++k) {
	mat[i][j]+=SType(At)*
	  ((*pols[k])[i]*std::conj((*pols[k])[j])+
	   SType(p_v->Kin()->Phase(0))*
	   (*pols[k])[i]*std::conj((*pols[1-k])[j]));
 	ref[i][j]=SType(A)*e1[i]*e1[j]+SType(A-B)*e2[i]*e2[j];
      }
      SComplex sb(B*qt[i]*qt[j]/qt.Abs2());
      if (i==j) sb-=i==0?A:-A;
      sb+=A*(p[i]*n[j]+p[j]*n[i])/(p*n);
      DEBUG_VAR(i<<" "<<j<<" "<<mat[i][j]<<" "<<ref[i][j]<<" "<<sb
		<<" "<<(mat[i][j]/sb-SComplex(1.0,0.0)));
      if (!IsEqual(mat[i][j].real(),sb.real(),1.0e-6)) {
	PRINT_VAR(i<<" "<<j<<" pol error "<<mat[i][j]<<" "<<sb
		  <<" "<<(mat[i][j]/sb-SComplex(1.0,0.0)));
      }
    }
#endif
  for (size_t cp(0);cp<2;++cp) {
    CVec4Type *j(GetPol(p,q,cp));
    *j*=m_cpll;
    j->SetH(cp+1);
    static_cast<Dipole_Color*>(p_cc)->AddJI(j,0);
    static_cast<Dipole_Color*>(p_cc)->AddJI(j,1);
    *j*=2.0/t*At;
    p_cc->AddJ(j);
    p_v->SetZero(false);
  }
#ifdef DEBUG__BG
  p_v->JC()->Print();
#endif
}

template <typename SType> CSpinor<SType> *
FFV_Calculator<SType>::GetPol(const Vec4D &q,const double &q2,
			      const int mode)
{
  int dir(p_v->Kin()->JI()->Direction());
  if (mode==0) {
    CSpinorType j(p_v->JC()->Flav().IsAnti()^(dir>0)?
		  CSpinorType(-1,-dir,1,q,0,0,0,0,q2):
		  CSpinorType(1,dir,1,q,0,0,0,0,q2));
#ifdef DEBUG__BG
    msg_Debugging()<<METHOD<<"(): "<<(dir>0?'I':'O')<<"+ "<<j<<" "
		   <<(dir>0?p_v->JC()->Flav().Bar():p_v->JC()->Flav())
		   <<", q2 = "<<q2<<" ("<<q.Mass()<<")\n";
#endif
    return CSpinorType::New(j);
  }
  else if (mode==1) {
    CSpinorType j(p_v->JC()->Flav().IsAnti()^(dir>0)?
		  CSpinorType(-1,-dir,-1,q,0,0,0,0,q2):
		  CSpinorType(1,dir,-1,q,0,0,0,0,q2));
#ifdef DEBUG__BG
    msg_Debugging()<<METHOD<<"(): "<<(dir>0?'I':'O')<<"- "<<j<<" "
		   <<(dir>0?p_v->JC()->Flav().Bar():p_v->JC()->Flav())
		   <<", q2 = "<<q2<<" ("<<q.Mass()<<")\n";
#endif
    return CSpinorType::New(j);
  }
  return NULL;
}

template <typename SType>
void FFV_Calculator<SType>::ConstructFVSDipole()
{
  const CSpinorType_Vector *ca(NULL);
  const CVec4Type_Vector *cb(NULL);
  if (m_dir==1) {
    ca=p_v->JA()->J().front().Get<CSpinorType>();
    cb=p_v->JB()->J().front().Get<CVec4Type>();
  }
  else {
    ca=p_v->JB()->J().front().Get<CSpinorType>();
    cb=p_v->JA()->J().front().Get<CVec4Type>();
  }
  const CObject_Vector *cc(&p_v->Kin()->JK()->J().front());
  if (m_dir==1 && !p_cc->Evaluate
      (ca->front(),cb->front(),cc->front())) return;
  if (m_dir==2 && !p_cc->Evaluate
      (cb->front(),ca->front(),cc->front())) return;
  bool iisf(p_v->Kin()->JI()->Flav().IsFermion());
  double A(0.0), t(0.0);
  if (p_v->Kin()->Type()==0) {
    double zi(p_v->Kin()->Z()), y(p_v->Kin()->Y());
    double pipj(p_v->Kin()->PI()*p_v->Kin()->PJ());
    double rv(1.0), mt2(0.0);
    if (p_v->Kin()->Massive()) {
      double mij2(sqr(p_v->Kin()->JIJT()->Flav().Mass()));
      double mk2(sqr(p_v->Kin()->JKT()->Flav().Mass()));
      double pij2(2.0*pipj+mij2), Q2(p_v->Kin()->Q2());
      rv=(Q2-pij2-mk2)/(Q2-mij2-mk2)*
	sqrt((sqr(Q2-mij2-mk2)-4.0*mij2*mk2)/
	     (sqr(Q2-pij2-mk2)-4.0*pij2*mk2));
      mt2=mij2/pipj;
    }
    if (iisf) A=2.0/(1.0-zi*(1.0-y))-rv*(1.0+zi+mt2);
    else A=2.0/(1.0-(1.0-zi)*(1.0-y))-rv*(2.0-zi+mt2);
    t=2.0*pipj;
  }
  else if (p_v->Kin()->Type()==2) {
    double zi(p_v->Kin()->Z()), y(p_v->Kin()->Y());
    double pipj(p_v->Kin()->PI()*p_v->Kin()->PJ());
    double mt2(p_v->Kin()->JIJT()->Flav().Mass());
    if (mt2) mt2*=mt2/pipj;
    if (iisf) A=2.0/(1.0-zi+y)-(1.0+zi+mt2);
    else A=2.0/(1.0-(1.0-zi)+y)-(2.0-zi+mt2);
    t=2.0*pipj*(1.0-y);
  }
  else if (p_v->Kin()->Type()==1) {
    double x(p_v->Kin()->Z()), ui(p_v->Kin()->Y());
    if (iisf) A=2.0/(1.0-x+ui)-(1.0+x);
    else A=1.0-2.0*x*(1.0-x);
    t=-2.0*(p_v->Kin()->PI()*p_v->Kin()->PJ())*x;
  }
  else {
    double x(p_v->Kin()->Z());
    if (iisf) A=2.0/(1.0-x)-(1.0+x);
    else A=1.0-2.0*x*(1.0-x);
    t=-2.0*(p_v->Kin()->PI()*p_v->Kin()->PJ())*x;
  }
  for (size_t cp(0);cp<2;++cp) {
    CSpinorType *j(GetPol(p_v->JC()->P(),sqr(p_v->JC()->Mass()),cp));
    *j*=m_cpll;
    j->SetH(cp+1);
    static_cast<Dipole_Color*>(p_cc)->AddJI(j,0);
    *j*=2.0*A/t;
    p_cc->AddJ(j);
    p_v->SetZero(false);
  }
#ifdef DEBUG__BG
  p_v->JC()->Print();
#endif
}

template <typename SType>
void FFV_Calculator<SType>::ConstructFVIDipole()
{
  Current *cj(m_dir==1?p_v->JA():p_v->JB());
  p_v->Kin()->JIJT()->SetP(cj->P());
  p_v->Kin()->JKT()->SetP(p_v->Kin()->JK()->P());
  const CSpinorType_Matrix *c(cj->J().Get<CSpinorType>());
  const CObject_Vector *cc(&p_v->Kin()->JK()->J().front());
  if (!p_cc->Evaluate(c->front().front(),cc->front())) return;
  double d(p_v->Info()->DRMode()?0.5:0.0);
  double a(p_v->Info()->AMax()), loga(log(a));
  p_v->Kin()->SetRes(1.0,2);
  p_v->Kin()->SetRes(3.0/2.0,1);
  p_v->Kin()->SetRes(5.0-sqr(M_PI)/2.0-sqr(loga)+1.5*(a-1-loga)-d,0);
  if (p_v->Kin()->JK()->Flav().IsGluon()) {
    double nf(Flavour(kf_quark).Size()/2);
    double Vqq(-16.0/9.0-2.0/3.0*(a-1.0-loga)-d);
    double Vgg(50.0/9.0-0.5*sqr(M_PI)-sqr(loga)+11.0/6.0*(a-1.0-loga)-d);
    p_v->Kin()->AddRes(1.0,2);
    p_v->Kin()->AddRes(11.0/6.0-2.0/3.0*0.5/3.0*nf,1);
    p_v->Kin()->AddRes(Vgg+0.5/3.0*Vqq*nf,0);
  }
  else {
    p_v->Kin()->AddRes(1.0,2);
    p_v->Kin()->AddRes(3.0/2.0,1);
    p_v->Kin()->AddRes(5.0-sqr(M_PI)/2.0-sqr(loga)+1.5*(a-1-loga)-d,0);
  }
  for (size_t i(0);i<c->size();++i) {
    CSpinorType *j((CSpinorType*)(*c)[i].front()->Copy());
    *j*=m_cpll*std::conj(m_cpll);
    static_cast<Dipole_Color*>(p_cc)->AddJI(j,0);
    static_cast<Dipole_Color*>(p_cc)->AddJI(j,1);
#ifndef DEBUG__BGS_AMAP
    j->Delete();
#else
    *j*=SComplex(1.0)/(m_cpll*std::conj(m_cpll));
    p_v->AddJ(j);
#endif
    p_v->SetZero(false);
  }
#ifdef DEBUG__BG
  p_v->JC()->Print();
#endif
}

template <typename SType>
void FFV_Calculator<SType>::Evaluate()
{
  p_v->SetZero();
  if (p_v->Info()) {
    if (p_v->Info()->Mode()==1) {
      if (p_v->JA()->Zero()||p_v->JB()->Zero()) return;
      switch (m_dir) {
      case 0: ConstructFFSDipole(); return;
      case 1: ConstructFVSDipole(); return;
      case 2: ConstructFVSDipole(); return;
      default: THROW(fatal_error,"Internal error");
      }
    }
    else {
      switch (m_dir) {
      case 1: ConstructFVIDipole(); return;
      case 2: ConstructFVIDipole(); return;
      default: THROW(fatal_error,"Internal error");
      }
    }
  }
  if (p_v->JA()->Zero()||p_v->JB()->Zero()) return;
#ifdef DEBUG__BG
  msg_Debugging()<<*p_v->JA()<<"(+)"<<*p_v->JB()<<" FFV("<<m_dir
		 <<"): m_cpll = "<<m_cpll<<", m_cplr = "<<m_cplr<<"\n";
  msg_Indent();
#endif
  const CObject_Matrix &cca(p_v->JA()->J()), &ccb(p_v->JB()->J());
  switch (m_dir) {
  case 0: {
    size_t i(0);
    for (CObject_Matrix::const_iterator 
	   jait(cca.begin());jait!=cca.end();++jait) {
      const CSpinorType_Vector *ca(jait->Get<CSpinorType>());
      if (ca->empty()) i+=ccb.size(); else
      for (CObject_Matrix::const_iterator 
	     jbit(ccb.begin());jbit!=ccb.end();++jbit,++i) {
	const CSpinorType_Vector *cb(jbit->Get<CSpinorType>());
	for (typename CSpinorType_Vector::const_iterator 
	       bit(cb->begin());bit!=cb->end();++bit)
	  for (typename CSpinorType_Vector::const_iterator 
		 ait(ca->begin());ait!=ca->end();++ait)
	    if (p_cc->Evaluate(*ait,*bit)) Evaluate(**ait,**bit,p_v->H(i));
      }
    }
    break;
  }
  case 1: {
    size_t i(0);
    for (CObject_Matrix::const_iterator 
	   jait(cca.begin());jait!=cca.end();++jait) {
      const CSpinorType_Vector *ca(jait->Get<CSpinorType>());
      if (ca->empty()) i+=ccb.size(); else
      for (CObject_Matrix::const_iterator 
	     jbit(ccb.begin());jbit!=ccb.end();++jbit,++i) {
	const CVec4Type_Vector *cb(jbit->Get<CVec4Type>());
	for (typename CVec4Type_Vector::const_iterator 
	     bit(cb->begin());bit!=cb->end();++bit)
	  for (typename CSpinorType_Vector::const_iterator 
		 ait(ca->begin());ait!=ca->end();++ait)
	    if (p_cc->Evaluate(*ait,*bit)) Evaluate(**ait,**bit,p_v->H(i));
      }
    }
    break;
  }
  case 2: {
    size_t i(0);
    for (CObject_Matrix::const_iterator 
	   jait(cca.begin());jait!=cca.end();++jait) {
      const CVec4Type_Vector *ca(jait->Get<CVec4Type>());
      if (ca->empty()) i+=ccb.size(); else
      for (CObject_Matrix::const_iterator 
	     jbit(ccb.begin());jbit!=ccb.end();++jbit,++i) {
	const CSpinorType_Vector *cb(jbit->Get<CSpinorType>());
	for (typename CSpinorType_Vector::const_iterator 
	       bit(cb->begin());bit!=cb->end();++bit)
	  for (typename CVec4Type_Vector::const_iterator 
		 ait(ca->begin());ait!=ca->end();++ait)
	    if (p_cc->Evaluate(*ait,*bit)) Evaluate(**bit,**ait,p_v->H(i));
      }
    }
    break;
  }
  default:
    THROW(fatal_error,"Internal error");
  }
}

template <typename SType> void FFV_Calculator<SType>::Evaluate
(const CSpinorType &a,const CSpinorType &b,const size_t &h)
{
  bool cl(m_cl&&CalcLeft(a,b)), cr(m_cr&&CalcRight(a,b));
  if (!(cl || cr)) return;
  CVec4Type *j(NULL);
  if (cl && cr) j=LorentzLeftRight(a,b);
  else if (cl) j=LorentzLeft(a,b);
  else if (cr) j=LorentzRight(a,b);
  j->SetH(h);
  p_cc->AddJ(j);
  p_v->SetZero(false);
}

template <typename SType> void FFV_Calculator<SType>::Evaluate
(const CSpinorType &a,const CVec4Type &b,const size_t &h)
{
  bool cl(m_cl&&CalcLeft(a)), cr(m_cr&&CalcRight(a));
  if (!(cl || cr)) return;
  CSpinorType *j(NULL);
  if (cl && cr) j=LorentzLeftRight(a,b);
  else if (cl) j=LorentzLeft(a,b);
  else if (cr) j=LorentzRight(a,b);
  j->SetH(h);
  p_cc->AddJ(j);
  p_v->SetZero(false);
}

template <typename SType> CVec4<SType> *
FFV_Calculator<SType>::LorentzLeft(CSpinorType a,CSpinorType b)
{
  if (a.B()>0) std::swap<CSpinorType>(a,b);
#ifdef DEBUG__BG
  msg_Debugging()<<"<> L "<<a<<"\n";
  msg_Debugging()<<"     "<<b<<"\n";
#endif
  SComplex j01(a[3]*b[1]), j02(a[2]*b[0]);
  SComplex j11(-a[2]*b[1]), j12(-a[3]*b[0]), j112(j11-j12);
  CVec4Type *j(CVec4Type::New(0.0,0.0,0.0,0.0,0,0,0,a.S()|b.S()));
  (*j)[0]=(j01+j02)*m_cpll;
  (*j)[Spinor<SType>::R3()]=(j01-j02)*m_cpll;
  (*j)[Spinor<SType>::R1()]=(j11+j12)*m_cpll;
  (*j)[Spinor<SType>::R2()]=SComplex(j112.imag(),-j112.real())*m_cpll;
  return j;
}

template <typename SType> CVec4<SType> *
FFV_Calculator<SType>::LorentzRight(CSpinorType a,CSpinorType b)
{
  if (a.B()>0) std::swap<CSpinorType>(a,b);
#ifdef DEBUG__BG
  msg_Debugging()<<"<> R "<<a<<"\n";
  msg_Debugging()<<"     "<<b<<"\n";
#endif
  SComplex j01(a[0]*b[2]), j02(a[1]*b[3]);
  SComplex j11(a[0]*b[3]), j12(a[1]*b[2]), j112(j11-j12);
  CVec4Type *j(CVec4Type::New(0.0,0.0,0.0,0.0,0,0,0,a.S()|b.S()));
  (*j)[0]=(j01+j02)*m_cplr;
  (*j)[Spinor<SType>::R3()]=(j01-j02)*m_cplr;
  (*j)[Spinor<SType>::R1()]=(j11+j12)*m_cplr;
  (*j)[Spinor<SType>::R2()]=SComplex(j112.imag(),-j112.real())*m_cplr;
  return j;
}

template <typename SType> CVec4<SType> *
FFV_Calculator<SType>::LorentzLeftRight(CSpinorType a,CSpinorType b)
{
  if (a.B()>0) std::swap<CSpinorType>(a,b);
#ifdef DEBUG__BG
  msg_Debugging()<<"<> LR "<<a<<"\n";
  msg_Debugging()<<"      "<<b<<"\n";
#endif
  SComplex l01(a[3]*b[1]), l02(a[2]*b[0]);
  SComplex l11(-a[2]*b[1]), l12(-a[3]*b[0]), l112(l11-l12);
  SComplex r01(a[0]*b[2]), r02(a[1]*b[3]);
  SComplex r11(a[0]*b[3]), r12(a[1]*b[2]), r112(r11-r12);
  CVec4Type *j(CVec4Type::New(0.0,0.0,0.0,0.0,0,0,0,a.S()|b.S()));
  (*j)[0]=(l01+l02)*m_cpll+(r01+r02)*m_cplr;
  (*j)[Spinor<SType>::R3()]=(l01-l02)*m_cpll+(r01-r02)*m_cplr;
  (*j)[Spinor<SType>::R1()]=(l11+l12)*m_cpll+(r11+r12)*m_cplr;
  (*j)[Spinor<SType>::R2()]=
    SComplex(l112.imag(),-l112.real())*m_cpll+
    SComplex(r112.imag(),-r112.real())*m_cplr;
  return j;
}

template <typename SType> CSpinor<SType> *
FFV_Calculator<SType>::LorentzLeft(const CSpinorType &a,
				   const CVec4Type &b)
{
  switch (a.B()) {
  case -1: {
#ifdef DEBUG__BG
    msg_Debugging()<<"<|g L "<<a<<"\n";
    msg_Debugging()<<"      "<<b<<"\n";
#endif
    CSpinorType *j(CSpinorType::New(a.R(),a.B(),0,0,0,a.S()|b.S(),1));
    SComplex jp(PPlus(b)), jm(PMinus(b)), jt(PT(b)), jtc(PTC(b));
    (*j)[0]=(a[2]*jp+a[3]*jt)*m_cpll;
    (*j)[1]=(a[2]*jtc+a[3]*jm)*m_cpll;
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
    (*j)[2]=(a[0]*jp+a[1]*jtc)*m_cpll;
    (*j)[3]=(a[0]*jt+a[1]*jm)*m_cpll;
    return j;
  }
  default:
    THROW(fatal_error,"Internal error");
  }
  return NULL;
}

template <typename SType> CSpinor<SType> *
FFV_Calculator<SType>::LorentzRight(const CSpinorType &a,
				    const CVec4Type &b)
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
    (*j)[2]=(a[0]*jm-a[1]*jt)*m_cplr;
    (*j)[3]=(-a[0]*jtc+a[1]*jp)*m_cplr;
    return j;
  }
  case 1: {
#ifdef DEBUG__BG
    msg_Debugging()<<"g|> R "<<a<<"\n";
    msg_Debugging()<<"      "<<b<<"\n";
#endif
    CSpinorType *j(CSpinorType::New(a.R(),a.B(),0,0,0,a.S()|b.S(),1));
    SComplex jp(PPlus(b)), jm(PMinus(b)), jt(PT(b)), jtc(PTC(b));
    (*j)[0]=(a[2]*jm-a[3]*jtc)*m_cplr;
    (*j)[1]=(-a[2]*jt+a[3]*jp)*m_cplr;
    (*j)[3]=(*j)[2]=SComplex(0.0,0.0);
    return j;
  }
  default:
    THROW(fatal_error,"Internal error");
  }
  return NULL;
}

template <typename SType> CSpinor<SType> *
FFV_Calculator<SType>::LorentzLeftRight(const CSpinorType &a,
					const CVec4Type &b)
{
  switch (a.B()) {
  case -1: {
#ifdef DEBUG__BG
    msg_Debugging()<<"<|g LR "<<a<<"\n";
    msg_Debugging()<<"       "<<b<<"\n";
#endif
    CSpinorType *j(CSpinorType::New(a.R(),a.B(),0,0,0,a.S()|b.S(),3));
    SComplex jp(PPlus(b)), jm(PMinus(b)), jt(PT(b)), jtc(PTC(b));
    (*j)[0]=(a[2]*jp+a[3]*jt)*m_cpll;
    (*j)[1]=(a[2]*jtc+a[3]*jm)*m_cpll;
    (*j)[2]=(a[0]*jm-a[1]*jt)*m_cplr;
    (*j)[3]=(-a[0]*jtc+a[1]*jp)*m_cplr;
    return j;
  }
  case 1: {
#ifdef DEBUG__BG
    msg_Debugging()<<"g|> LR "<<a<<"\n";
    msg_Debugging()<<"       "<<b<<"\n";
#endif
    CSpinorType *j(CSpinorType::New(a.R(),a.B(),0,0,0,a.S()|b.S(),3));
    SComplex jp(PPlus(b)), jm(PMinus(b)), jt(PT(b)), jtc(PTC(b));
    (*j)[0]=(a[2]*jm-a[3]*jtc)*m_cplr;
    (*j)[1]=(-a[2]*jt+a[3]*jp)*m_cplr;
    (*j)[2]=(a[0]*jp+a[1]*jtc)*m_cpll;
    (*j)[3]=(a[0]*jt+a[1]*jm)*m_cpll;
    return j;
  }
  default:
    THROW(fatal_error,"Internal error");
  }
  return NULL;
}

template <typename SType>
std::string FFV_Calculator<SType>::Label() const
{
  return "FFV["+ToString(m_cpll)+",,"+ToString(m_cplr)+"]";
}

namespace METOOLS {

  template class FFV_Calculator<double>;
  template class FFV_Calculator<long double>;

}

DECLARE_TEMPLATE_GETTER(FFV_Getter,"Gamma",Lorentz_Calculator,Vertex_Key);

template <typename SType,char STag> Lorentz_Calculator *
FFV_Getter<SType,STag>::operator()(const Vertex_Key &key) const
{
  return new FFV_Calculator<SType>(key);
}

template <typename SType,char STag>
void FFV_Getter<SType,STag>::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"FFV vertex";
}

template class FFV_Getter<double,'D'>;
template class FFV_Getter<long double,'Q'>;
