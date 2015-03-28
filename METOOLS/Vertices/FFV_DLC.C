#include "METOOLS/Explicit/Lorentz_Calculator.H"

#include "METOOLS/Explicit/Vertex.H"
#include "METOOLS/Explicit/Dipole_Kinematics.H"
#include "METOOLS/Explicit/Dipole_Color.H"
#include "METOOLS/Explicit/Dipole_Terms.H"
#include "METOOLS/Currents/C_Spinor.H"
#include "MODEL/Main/Single_Vertex.H"

namespace METOOLS {

  template <typename SType>
  class FFV_DCalculator: public Lorentz_Calculator {
  public:

    typedef std::complex<SType> SComplex;

    typedef CSpinor<SType> CSpinorType;
    typedef std::vector<CSpinorType*> CSpinorType_Vector;
    typedef std::vector<CSpinorType_Vector> CSpinorType_Matrix;

    typedef CVec4<SType> CVec4Type;
    typedef std::vector<CVec4Type*> CVec4Type_Vector;

  private:

    Color_Calculator *p_cc;

    SComplex m_cpll, m_cplr;

    int m_dir, m_cl, m_cr;

    double m_mi, m_mi2, m_mj, m_mj2, m_mk, m_mk2, m_mij, m_mij2;

    CVec4Type *GetPol(const ATOOLS::Vec4D &p,
		      const ATOOLS::Vec4D &q,const int mode);
    CSpinorType *GetPol(const ATOOLS::Vec4D &q,
			const double &q2,const int mode);

  public:
    
    FFV_DCalculator(const Vertex_Key &key);

    std::string Label() const { return "XFFV"; }

    CObject *Evaluate(const CObject_Vector &j) { return NULL; }
    
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

  };// end of class FFV_DCalculator

}// end of namespace METOOLS

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"

#define ZERO SComplex(0.0,0.0)

using namespace METOOLS;
using namespace ATOOLS;

template <typename SType>
FFV_DCalculator<SType>::FFV_DCalculator(const Vertex_Key &key): 
  Lorentz_Calculator(key), p_cc(key.p_cc),
  m_dir(key.Fl(1).IsFermion()?
	(key.Fl(0).IsFermion()?0:2):1),
  m_mi(-1.0), m_mi2(-1.0), m_mj(-1.0), m_mj2(-1.0),
  m_mij(-1.0), m_mij2(-1.0), m_mk(-1.0), m_mk2(-1.0)
{
  if (p_v->Kin()) {
    if (p_v->Kin()->JI())
      m_mi2=sqr(m_mi=p_v->Kin()->JI()->Flav().Mass());
    if (p_v->Kin()->JJ())
      m_mj2=sqr(m_mj=p_v->Kin()->JJ()->Flav().Mass());
    m_mk2=sqr(m_mk=p_v->Kin()->JK()->Flav().Mass());
    m_mij2=sqr(m_mij=p_v->Kin()->JIJT()->Flav().Mass());
  }
  bool nuc(p_v->Info() && p_v->Info()->Mode()&2);
  m_cpll=SComplex((nuc?1.0:p_v->V()->cpl.front().Value())*p_cc->Coupling());
  m_cplr=SComplex((nuc?1.0:p_v->V()->cpl.front().Value())*p_cc->Coupling());
  m_cl=m_cpll!=SComplex(0.0,0.0);
  m_cr=m_cplr!=SComplex(0.0,0.0);
}

template <typename SType> CVec4<SType> *
FFV_DCalculator<SType>::GetPol
(const Vec4D &p,const Vec4D &q,const int mode)
{
  static double sqrttwo(sqrt(2.0));
  Vec3D qxp(cross(Vec3D(q),Vec3D(p))), qxpxp(p*q[0]-q*p[0]);
  CVec4Type e1(Vec4D(0.0,qxp/(sqrttwo*qxp.Abs())));
  CVec4Type e2(Vec4D(0.0,qxpxp/(sqrttwo*qxpxp.Abs())));
  return CVec4Type::New(e1+SComplex(0.0,mode?1.0:-1.0)*e2); 
}

template <typename SType>
void FFV_DCalculator<SType>::ConstructFFSDipole()
{
  const CSpinorType_Vector *ca(p_v->J(0)->J().front().Get<CSpinorType>());
  const CSpinorType_Vector *cb(p_v->J(1)->J().front().Get<CSpinorType>());
  CObject_Vector cj(3);
  cj[0]=ca->front();
  cj[1]=cb->front();
  cj[2]=p_v->Kin()->JK()->J().front().front();
  if (!p_cc->Evaluate(cj)) return;
  Vec4D p(p_v->JC()->P()), q;
  double A(0.0), B(0.0), t(0.0);
  if (p_v->Kin()->Type()==0) {
    double zi(p_v->Kin()->Z()), zj(1.0-zi), rv(1.0);
    Vec4D pi(p_v->Kin()->PI()), pj(p_v->Kin()->PJ());
    double pij2((pi+pj).Abs2()), mt2(0.0), zim(zi), zjm(zj);
    p_v->Kin()->SetA(zim*zjm);
    if (p_v->Kin()->Massive()) {
      double y(p_v->Kin()->Y()), Q2(p_v->Kin()->Q2());
      double s(Q2-m_mi2-m_mj2-m_mk2);
      double viji(sqrt(sqr(s*y)-sqr(2.0*m_mi*m_mj))/(s*y+2.0*m_mi2));
      rv=sqrt(sqr(2.0*m_mk2+s*(1.0-y))-4.0*m_mk2*Q2)/(s*(1.0-y));
      double zc(0.5*(2.0*m_mi2+s*y)/(m_mi2+m_mj2+s*y));
      double zm(zc*(1.0-viji*rv)), zp(zc*(1.0+viji*rv));
      mt2=2.0*p_v->Info()->Kappa()*(zp*zm-m_mi2/pij2);
      zim-=0.5*(1.0-rv);
      zjm-=0.5*(1.0-rv);
      p_v->Kin()->SetA(zim*zjm-zp*zm);
    }
    q=zim*pi-zjm*pj;
    A=1.0-mt2;
    B=-4.0*q.Abs2()/pij2;
    t=rv*pij2;
    p_v->Kin()->SetA(A-2.0*p_v->Kin()->A());
  }
  else if (p_v->Kin()->Type()==2) {
    double zi(p_v->Kin()->Z()), zj(1.0-zi);
    Vec4D pi(p_v->Kin()->PI()), pj(p_v->Kin()->PJ());
    double pij2((pi+pj).Abs2());
    q=zi*pi-zj*pj;
    A=1.0;
    B=-4.0*q.Abs2()/pij2;
    t=pij2*(1.0-p_v->Kin()->Y());
    p_v->Kin()->SetA((1.0-zi)*zi);
    if (p_v->Kin()->Massive()) {
      double Q2((pi+pj+p_v->Kin()->PK()).Abs2());
      double mui2(m_mi2/Q2), y(p_v->Kin()->Y());
      double eps(sqrt(sqr(y-2.0*mui2)-4.0*mui2*mui2)/y);
      p_v->Kin()->SetA((0.5*(1.0+eps)-zi)*(zi-0.5*(1.0-eps)));
    }
    p_v->Kin()->SetA(A-2.0*p_v->Kin()->A());
  }
  else if (p_v->Kin()->Type()==1) {
    double x(p_v->Kin()->Z()), ui(p_v->Kin()->Y());
    Vec4D pi(p_v->Kin()->PJ()), pk(p_v->Kin()->PK());
    q=pi/ui-pk/(1.0-ui);
    A=x;
    B=2.0*(1.0-x)/x*ui*(1.0-ui)*q.Abs2()/(pi*pk);
    t=-2.0*(pi*p_v->Kin()->PI())*x;
    p_v->Kin()->SetA((1.0-x)/x);
    if (p_v->Kin()->Massive()) {
      double Q2(2.0*(p_v->Kin()->JKT()->P()*p_v->JC()->P()));
      p_v->Kin()->SetA((1.0-x)/x-pk.Abs2()/Q2*ui/(1.0-ui));
    }
    p_v->Kin()->SetA(A+2.0*p_v->Kin()->A());
  }
  else {
    double x(p_v->Kin()->Z()), vi(p_v->Kin()->Y());
    Vec4D pi(p_v->Kin()->PJ()), pk(-p_v->Kin()->PK());
    A=x;
    B=-4.0*(1.0-x)/x;
    q=pi-vi*pk;
    t=-2.0*(pi*p_v->Kin()->PI())*x;
    p_v->Kin()->SetA(A+2.0*(1.0-x)/x);
  }
  p_v->Kin()->CheckKT2Min(); 
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
FFV_DCalculator<SType>::GetPol(const Vec4D &q,const double &q2,
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
void FFV_DCalculator<SType>::ConstructFVSDipole()
{
  const CSpinorType_Vector *ca(NULL);
  const CVec4Type_Vector *cb(NULL);
  if (m_dir==1) {
    ca=p_v->J(0)->J().front().Get<CSpinorType>();
    cb=p_v->J(1)->J().front().Get<CVec4Type>();
  }
  else {
    ca=p_v->J(1)->J().front().Get<CSpinorType>();
    cb=p_v->J(0)->J().front().Get<CVec4Type>();
  }
  CObject_Vector cj(3);
  cj[0]=m_dir==1?(CObject*)ca->front():(CObject*)cb->front();
  cj[1]=m_dir==1?(CObject*)cb->front():(CObject*)ca->front();
  cj[2]=p_v->Kin()->JK()->J().front().front();
  if (!p_cc->Evaluate(cj)) return;
  bool iisf(p_v->Kin()->JI()->Flav().IsFermion());
  double A(0.0), t(0.0);
  if (p_v->Kin()->Type()==0) {
    double zi(p_v->Kin()->Z()), y(p_v->Kin()->Y());
    double pipj(p_v->Kin()->PI()*p_v->Kin()->PJ());
    double rv(1.0), mt2(0.0);
    if (p_v->Kin()->Massive()) {
      double pij2(2.0*pipj+m_mij2), Q2(p_v->Kin()->Q2());
      rv=(Q2-pij2-m_mk2)/(Q2-m_mij2-m_mk2)*
	sqrt((sqr(Q2-m_mij2-m_mk2)-4.0*m_mij2*m_mk2)/
	     (sqr(Q2-pij2-m_mk2)-4.0*pij2*m_mk2));
      mt2=m_mij2/pipj;
    }
    if (iisf) A=2.0/(1.0-zi*(1.0-y))-rv*(1.0+zi+mt2);
    else A=2.0/(1.0-(1.0-zi)*(1.0-y))-rv*(2.0-zi+mt2);
    t=2.0*pipj;
    p_v->Kin()->SetA(A);
  }
  else if (p_v->Kin()->Type()==2) {
    double zi(p_v->Kin()->Z()), y(p_v->Kin()->Y());
    double pipj(p_v->Kin()->PI()*p_v->Kin()->PJ());
    double mt2(m_mij?m_mij2/pipj:0.0);
    if (iisf) A=2.0/(1.0-zi+y)-(1.0+zi+mt2);
    else A=2.0/(1.0-(1.0-zi)+y)-(2.0-zi+mt2);
    t=2.0*pipj*(1.0-y);
    p_v->Kin()->SetA(A);
  }
  else if (p_v->Kin()->Type()==1) {
    double x(p_v->Kin()->Z()), ui(p_v->Kin()->Y());
    if (iisf) A=2.0/(1.0-x+ui)-(1.0+x);
    else A=1.0-2.0*x*(1.0-x);
    t=-2.0*(p_v->Kin()->PI()*p_v->Kin()->PJ())*x;
    p_v->Kin()->SetA(A);
  }
  else {
    double x(p_v->Kin()->Z());
    if (iisf) A=2.0/(1.0-x)-(1.0+x);
    else A=1.0-2.0*x*(1.0-x);
    t=-2.0*(p_v->Kin()->PI()*p_v->Kin()->PJ())*x;
    p_v->Kin()->SetA(A);
  }
  p_v->Kin()->CheckKT2Min(); 
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
void FFV_DCalculator<SType>::ConstructFVIDipole()
{
  Current *cj(p_v->J(0));
  p_v->Kin()->JIJT()->SetP(cj->P());
  p_v->Kin()->JKT()->SetP(p_v->Kin()->JK()->P());
  const CSpinorType_Matrix *c(cj->J().Get<CSpinorType>());
  CObject_Vector cc(2);
  cc[0]=c->front().front();
  cc[1]=p_v->Kin()->JK()->J().front().front();
  if (!p_cc->Evaluate(cc)) return;
  double d(p_v->Info()->DRMode()?0.5:0.0);
  I_Args ia(p_v->Kin()->JIJT()->P(),
	    p_v->Kin()->JKT()->P(),m_mij,m_mk);
  NLO_Value iv(FFQQ(ia,p_v->Info()));
  p_v->Kin()->SetRes(iv.m_e2,2);
  p_v->Kin()->SetRes(iv.m_e1,1);
  p_v->Kin()->SetRes(iv.m_f-d,0);
  ia.Swap();
  if (!p_v->Kin()->JK()->Flav().IsGluon()) {
    iv=FFQQ(ia,p_v->Info());
  }
  else {
    d=p_v->Info()->DRMode()?1.0/6.0:0.0;
    double nf(Flavour(kf_quark).Size()/2);
    iv=0.5/3.0*nf*FFGQ(ia,p_v->Info(),0.0);
    for (size_t i(nf+1);i<=p_v->Info()->Nf();++i)
      iv+=0.5/3.0*FFGQ(ia,p_v->Info(),Flavour(i).Mass());
    iv+=FFGG(ia,p_v->Info());
  }
  p_v->Kin()->AddRes(iv.m_e2,2);
  p_v->Kin()->AddRes(iv.m_e1,1);
  p_v->Kin()->AddRes(iv.m_f-d,0);
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
void FFV_DCalculator<SType>::Evaluate()
{
  p_v->SetZero();
  if (p_v->Info()==NULL) THROW(fatal_error,"Invalid call");
  if (p_v->Info()->Mode()==1) {
    if (p_v->J(0)->Zero()||p_v->J(1)->Zero()) return;
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

namespace METOOLS {

  template class FFV_DCalculator<double>;

}

DECLARE_GETTER(FFV_DCalculator<double>,"DXFFV",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,FFV_DCalculator<double> >::
operator()(const Vertex_Key &key) const
{ return new FFV_DCalculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    FFV_DCalculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"FFV dipole vertex"; }
