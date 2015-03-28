#include "METOOLS/Explicit/Lorentz_Calculator.H"

#include "METOOLS/Explicit/Vertex.H"
#include "METOOLS/Explicit/Dipole_Kinematics.H"
#include "METOOLS/Explicit/Dipole_Color.H"
#include "METOOLS/Explicit/Dipole_Terms.H"
#include "METOOLS/Currents/C_Vector.H"
#include "ATOOLS/Phys/Spinor.H"

namespace METOOLS {

  template <typename SType>
  class VVV_SCalculator: public Lorentz_Calculator {
  public:

    typedef std::complex<SType> SComplex;

    typedef ATOOLS::Spinor<SType> SpinorType;

    typedef CVec4<SType> CVec4Type;
    typedef std::vector<CVec4Type*> CVec4Type_Vector;
    typedef std::vector<CVec4Type_Vector> CVec4Type_Matrix;

    CVec4Type *GetPol(const ATOOLS::Vec4D &p,
		      const ATOOLS::Vec4D &q,const int mode);

  private:

    Color_Calculator *p_cc;

    SComplex m_cpl;

    double m_mk, m_mk2;

    void ConstructSDipole();
    void ConstructIDipole();

  public:
    
    VVV_SCalculator(const Vertex_Key &key);
    
    std::string Label() const { return "XVVV"; }

    CObject *Evaluate(const CObject_Vector &j) { return NULL; }

    void Evaluate();

    void SetGauge(const ATOOLS::Vec4D &k);

  };// end of class VVV_SCalculator

}// end of namespace METOOLS

#include "MODEL/Main/Single_Vertex.H" 
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

using namespace METOOLS;
using namespace ATOOLS;

template <typename SType>
VVV_SCalculator<SType>::VVV_SCalculator(const Vertex_Key &key): 
  Lorentz_Calculator(key), p_cc(key.p_cc), m_mk(-1.0), m_mk2(-1.0)
{
  if (p_v->Kin()) m_mk2=sqr(m_mk=p_v->Kin()->JK()->Flav().Mass());
  bool nuc(p_v->Info() && p_v->Info()->Mode()&2);
  m_cpl=SComplex(-(nuc?1.0:p_v->V()->cpl.front().Value())*p_cc->Coupling());
}

template <typename SType> CVec4<SType> *
VVV_SCalculator<SType>::GetPol
(const Vec4D &p,const Vec4D &q,const int mode)
{
  static double sqrttwo(sqrt(2.0));
  Vec3D qxp(cross(Vec3D(q),Vec3D(p))), qxpxp(p*q[0]-q*p[0]);
  CVec4Type e1(Vec4D(0.0,qxp/(sqrttwo*qxp.Abs())));
  CVec4Type e2(Vec4D(0.0,qxpxp/(sqrttwo*qxpxp.Abs())));
  return CVec4Type::New(e1+SComplex(0.0,mode?1.0:-1.0)*e2); 
}

template <typename SType>
void VVV_SCalculator<SType>::ConstructSDipole()
{
  const CVec4Type_Vector *ca(p_v->J(0)->J().front().Get<CVec4Type>());
  const CVec4Type_Vector *cb(p_v->J(1)->J().front().Get<CVec4Type>());
  CObject_Vector cc(3);
  cc[0]=ca->front();
  cc[1]=cb->front();
  cc[2]=p_v->Kin()->JK()->J().front().front();
  if (!p_cc->Evaluate(cc)) return;
  Vec4D p(p_v->JC()->P()), q;
  double Ai(0.0), Aj(0.0), B(0.0), t(0.0);
  if (p_v->Kin()->Type()==0) {
    double zi(p_v->Kin()->Z()), zj(1.0-zi), y(p_v->Kin()->Y());
    Vec4D pi(p_v->Kin()->PI()), pj(p_v->Kin()->PJ());
    double pipj(pi*pj), rv(1.0), sl(1.0), zim(zi), zjm(zj);
    p_v->Kin()->SetA(zim*zjm);
    if (p_v->Kin()->Massive()) {
      double y(p_v->Kin()->Y()), Q2(p_v->Kin()->Q2()), s(Q2-m_mk2);
      rv=sqrt(sqr(2.0*m_mk2+s*(1.0-y))-4.0*m_mk2*Q2)/(s*(1.0-y));
      double zm(0.5*(1.0-rv)), zp(0.5*(1.0+rv));
      sl=(1.0-0.5*p_v->Info()->Kappa()*zp*zm)/rv;
      zim-=0.5*(1.0-rv);
      zjm-=0.5*(1.0-rv);
      p_v->Kin()->SetA((zim*zjm-zp*zm)/rv);
    }
    Ai=2.0*(1.0/(1.0-zi*(1.0-y))-sl);
    Aj=2.0*(1.0/(1.0-zj*(1.0-y))-sl);
    if (p_v->Kin()->Swap()) std::swap<double>(Ai,Aj);
    q=zim*pi-zjm*pj;
    B=q.Abs2()/pipj/rv;
    t=2.0*pipj;
    p_v->Kin()->SetA(Ai+Aj+2.0*p_v->Kin()->A());
  }
  else if (p_v->Kin()->Type()==2) {
    double zi(p_v->Kin()->Z()), zj(1.0-zi), y(p_v->Kin()->Y());
    Vec4D pi(p_v->Kin()->PI()), pj(p_v->Kin()->PJ());
    Ai=2.0*(zi-y)/(1.0-zi+y);
    Aj=2.0*(zj-y)/(1.0-zj+y);
    if (p_v->Kin()->Swap()) std::swap<double>(Ai,Aj);
    B=-2.0*zi*zj;
    q=zi*pi-zj*pj;
    t=2.0*(pi*pj)*(1.0-y);
    p_v->Kin()->SetA(Ai+Aj+2.0*zi*zj);
  }
  else if (p_v->Kin()->Type()==1) {
    double x(p_v->Kin()->Z()), ui(p_v->Kin()->Y());
    Vec4D pi(p_v->Kin()->PJ()), pk(p_v->Kin()->PK());
    Ai=2.0*(x-ui)/(1.0-x+ui);
    Aj=2.0*x*(1.0-x);
    if (p_v->Kin()->Swap()) std::swap<double>(Ai,Aj);
    q=pi/ui-pk/(1.0-ui);
    B=(1.0-x)/x*ui*(1.0-ui)*q.Abs2()/(pi*pk);
    t=-2.0*(pi*p_v->Kin()->PI())*x;
    p_v->Kin()->SetA((1.0-x)/x);
    if (p_v->Kin()->Massive()) {
      double Q2(2.0*(p_v->Kin()->JKT()->P()*p_v->JC()->P()));
      p_v->Kin()->SetA((1.0-x)/x-pk.Abs2()/Q2*ui/(1.0-ui));
    }
    p_v->Kin()->SetA(Ai+Aj+2.0*p_v->Kin()->A());
  }
  else {
    double x(p_v->Kin()->Z()), vi(p_v->Kin()->Y());
    Vec4D pi(p_v->Kin()->PJ()), pk(-p_v->Kin()->PK());
    Ai=2.0*x/(1.0-x);
    Aj=2.0*x*(1.0-x);
    if (p_v->Kin()->Swap()) std::swap<double>(Ai,Aj);
    B=-2.0*(1.0-x)/x;
    q=pi-vi*pk;
    t=-2.0*(pi*p_v->Kin()->PI())*x;
    p_v->Kin()->SetA(Ai+Aj+2.0*(1.0-x)/x);
  }
  p_v->Kin()->CheckKT2Min(); 
  double Ait(Ai-B/2.0), Ajt(Aj-B/2.0);
  p_v->Kin()->SetPhase(1.0/(2.0*Ai/B-1.0),0);
  p_v->Kin()->SetPhase(1.0/(2.0*Aj/B-1.0),1);
#ifdef CHECK__pols
  DEBUG_VAR(p_v->Kin()->Type()<<" "<<Ai<<" "<<Aj<<" "<<B<<" "<<q);
  CVec4<SType> *pols[2]={GetPol(p,q,0),GetPol(p,q,1)};
  SComplex mat[4][4], tam[4][4];
  Vec4D n(p[0],Vec3D(-p)), qt(p*q[0]-q*p[0]);
  for (size_t i(0);i<4;++i)
    for (size_t j(0);j<4;++j) {
      mat[i][j]=tam[i][j]=0.0;
      for (size_t k(0);k<2;++k) {
	mat[i][j]+=SType(Ait)*
	  ((*pols[k])[i]*std::conj((*pols[k])[j])+
	   SType(p_v->Kin()->Phase(0))*
	   (*pols[k])[i]*std::conj((*pols[1-k])[j]));
	tam[i][j]+=SType(Ajt)*
	  ((*pols[k])[i]*std::conj((*pols[k])[j])+
	   SType(p_v->Kin()->Phase(1))*
	   (*pols[k])[i]*std::conj((*pols[1-k])[j]));
      }
      SComplex sb(B*qt[i]*qt[j]/qt.Abs2());
      SComplex bs(B*qt[i]*qt[j]/qt.Abs2());
      if (i==j) {
	sb-=i==0?Ai:-Ai;
	bs-=i==0?Aj:-Aj;
      }
      sb+=Ai*(p[i]*n[j]+p[j]*n[i])/(p*n);
      bs+=Aj*(p[i]*n[j]+p[j]*n[i])/(p*n);
      DEBUG_VAR(i<<" "<<j<<" "<<mat[i][j]<<" "<<sb
		<<" "<<(mat[i][j]/sb-SComplex(1.0,0.0)));
      DEBUG_VAR(i<<" "<<j<<" "<<tam[i][j]<<" "<<bs
		<<" "<<(tam[i][j]/bs-SComplex(1.0,0.0)));
      if (!IsEqual(mat[i][j].real(),sb.real(),1.0e-6)) {
	PRINT_VAR(i<<" "<<j<<" pol error "<<mat[i][j]<<" "<<sb
		  <<" "<<(mat[i][j]/sb-SComplex(1.0,0.0)));
      }
      if (!IsEqual(tam[i][j].real(),bs.real(),1.0e-6)) {
	PRINT_VAR(i<<" "<<j<<" pol error "<<tam[i][j]<<" "<<bs
		  <<" "<<(tam[i][j]/bs-SComplex(1.0,0.0)));
      }
    }
#endif
  for (size_t cp(0);cp<2;++cp) {
    CVec4Type *j(GetPol(p,q,cp));
    *j*=m_cpl;
    j->SetH(cp+1);
    if (Ait) {
    *j*=Ait;
    static_cast<Dipole_Color*>(p_cc)->AddJI(j,0);
    *j*=1.0/Ait;
    }
    if (Ajt) {
    *j*=Ajt;
    static_cast<Dipole_Color*>(p_cc)->AddJI(j,1);
    *j*=1.0/Ajt;
    }
    *j*=2.0/t;
    p_cc->AddJ(j);
    p_v->SetZero(false);
  }
#ifdef DEBUG__BG
  p_v->JC()->Print();
#endif
}

template <typename SType>
void VVV_SCalculator<SType>::ConstructIDipole()
{
  Current *cj(p_v->J(0));
  p_v->Kin()->JIJT()->SetP(cj->P());
  p_v->Kin()->JKT()->SetP(p_v->Kin()->JK()->P());
  const CVec4Type_Matrix *c(cj->J().Get<CVec4Type>());
  CObject_Vector cc(2);
  cc[0]=c->front().front();
  cc[1]=p_v->Kin()->JK()->J().front().front();
  if (!p_cc->Evaluate(cc)) return;
  I_Args ia(p_v->Kin()->JIJT()->P(),
	    p_v->Kin()->JKT()->P(),0.0,m_mk);
  double nf(Flavour(kf_quark).Size()/2);
  double d(p_v->Info()->DRMode()?1.0/6.0:0.0);
  NLO_Value iv(0.5/3.0*nf*FFGQ(ia,p_v->Info(),0.0));
  for (size_t i(nf+1);i<=p_v->Info()->Nf();++i)
    iv+=0.5/3.0*FFGQ(ia,p_v->Info(),Flavour(i).Mass());
  iv+=FFGG(ia,p_v->Info());
  p_v->Kin()->SetRes(iv.m_e2,2);
  p_v->Kin()->SetRes(iv.m_e1,1);
  p_v->Kin()->SetRes(iv.m_f-d,0);
  ia.Swap();
  if (p_v->Kin()->JK()->Flav().IsGluon()) {
    double nf(Flavour(kf_quark).Size()/2);
    iv=0.5/3.0*nf*FFGQ(ia,p_v->Info(),0.0);
    for (size_t i(nf+1);i<=p_v->Info()->Nf();++i)
      iv+=0.5/3.0*FFGQ(ia,p_v->Info(),Flavour(i).Mass());
    iv+=FFGG(ia,p_v->Info());
  }
  else {
    d=p_v->Info()->DRMode()?0.5:0.0;
    iv=FFQQ(ia,p_v->Info());
  }
  p_v->Kin()->AddRes(iv.m_e2,2);
  p_v->Kin()->AddRes(iv.m_e1,1);
  p_v->Kin()->AddRes(iv.m_f-d,0);
  for (size_t i(0);i<c->size();++i) {
    CVec4Type *j((CVec4Type*)(*c)[i].front()->Copy());
    *j*=m_cpl*std::conj(m_cpl);
    static_cast<Dipole_Color*>(p_cc)->AddJI(j,0);
    static_cast<Dipole_Color*>(p_cc)->AddJI(j,1);
#ifndef DEBUG__BGS_AMAP
    j->Delete();
#else
    *j*=SComplex(1.0)/(m_cpl*std::conj(m_cpl));
    p_v->AddJ(j);
#endif
    p_v->SetZero(false);
  }
#ifdef DEBUG__BG
  p_v->JC()->Print();
#endif
}

template <typename SType>
void VVV_SCalculator<SType>::Evaluate()
{
  p_v->SetZero();
  if (p_v->Kin()==NULL) THROW(fatal_error,"Invalid call");
  if (p_v->Info()->Mode()==1) {
    if (p_v->J(0)->Zero()||p_v->J(1)->Zero()) return;
    ConstructSDipole();
    return;
  }
  else {
    ConstructIDipole();
    return;
  }
}

namespace METOOLS {

  template class VVV_SCalculator<double>;

}

DECLARE_GETTER(VVV_SCalculator<double>,"DXVVV",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,VVV_SCalculator<double> >::
operator()(const Vertex_Key &key) const
{ return new VVV_SCalculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    VVV_SCalculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"VVV subtraction vertex"; }
