#include "METOOLS/Explicit/Lorentz_Calculator.H"

#include "METOOLS/Explicit/Vertex.H"
#include "METOOLS/Explicit/Dipole_Kinematics.H"
#include "METOOLS/Explicit/Dipole_Color.H"
#include "METOOLS/Currents/C_Vector.H"
#include "ATOOLS/Phys/Spinor.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"

namespace METOOLS {

  template <typename SType>
  class VVV_Calculator: public Lorentz_Calculator {
  public:

    typedef std::complex<SType> SComplex;

    typedef ATOOLS::Spinor<SType> SpinorType;

    typedef CVec4<SType> CVec4Type;
    typedef std::vector<CVec4Type*> CVec4Type_Vector;
    typedef std::vector<CVec4Type_Vector> CVec4Type_Matrix;

    CVec4Type *GetPol(const ATOOLS::Vec4D &p,
		      const ATOOLS::Vec4D &q,const int mode);

  private:

    SComplex m_cpl;

    CVec4Type Lorentz(const CVec4Type &a,const CVec4Type &b);

    void ConstructSDipole();
    void ConstructIDipole();

  public:
    
    VVV_Calculator(const Vertex_Key &key);
    
    std::string Label() const;

    void Evaluate();

    void SetGauge(const ATOOLS::Vec4D &k);

  };// end of class VVV_Calculator

}// end of namespace METOOLS

#include "ATOOLS/Org/Message.H"

using namespace METOOLS;
using namespace ATOOLS;

template <typename SType>
VVV_Calculator<SType>::VVV_Calculator(const Vertex_Key &key): 
  Lorentz_Calculator(key) 
{
  bool nuc(p_v->Info() && p_v->Info()->Mode()&2);
  m_cpl=SComplex((nuc?1.0:p_v->Coupling(0))*p_cc->Coupling());
}

template <typename SType>
std::string VVV_Calculator<SType>::Label() const
{
  return "VVV["+ToString(m_cpl)+"]";
}

template <typename SType> CVec4<SType> *
VVV_Calculator<SType>::GetPol
(const Vec4D &p,const Vec4D &q,const int mode)
{
  static double sqrttwo(sqrt(2.0));
  Vec3D qxp(cross(Vec3D(q),Vec3D(p))), qxpxp(p*q[0]-q*p[0]);
  CVec4Type e1(Vec4D(0.0,qxp/(sqrttwo*qxp.Abs())));
  CVec4Type e2(Vec4D(0.0,qxpxp/(sqrttwo*qxpxp.Abs())));
  return CVec4Type::New(e1+SComplex(0.0,mode?1.0:-1.0)*e2); 
}

template <typename SType>
void VVV_Calculator<SType>::ConstructSDipole()
{
  const CVec4Type_Vector *ca(p_v->JA()->J().front().Get<CVec4Type>());
  const CVec4Type_Vector *cb(p_v->JB()->J().front().Get<CVec4Type>());
  const CObject_Vector *cc(&p_v->Kin()->JK()->J().front());
  if (!p_cc->Evaluate(ca->front(),cb->front(),cc->front())) return;
  Vec4D p(p_v->JC()->P()), q;
  double Ai(0.0), Aj(0.0), B(0.0), t(0.0);
  if (p_v->Kin()->Type()==0) {
    double zi(p_v->Kin()->Z()), zj(1.0-zi), y(p_v->Kin()->Y());
    Vec4D pi(p_v->Kin()->PI()), pj(p_v->Kin()->PJ());
    Ai=2.0*zi*(1.0-y)/(1.0-zi*(1.0-y));
    Aj=2.0*zj*(1.0-y)/(1.0-zj*(1.0-y));
    if (p_v->Kin()->Swap()) std::swap<double>(Ai,Aj);
    B=-2.0*zi*zj;
    q=zi*pi-zj*pj;
    t=2.0*(pi*pj);
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
  }
  else if (p_v->Kin()->Type()==1) {
    double x(p_v->Kin()->Z()), ui(p_v->Kin()->Y());
    Vec4D pi(p_v->Kin()->PJ()), pk(p_v->Kin()->PK());
    Ai=2.0*(x-ui)/(1.0-x+ui);
    Aj=2.0*x*(1.0-x);
    if (p_v->Kin()->Swap()) std::swap<double>(Ai,Aj);
    B=-2.0*(1.0-x)/x;
    q=pi/ui-pk/(1.0-ui);
    t=-2.0*(pi*p_v->Kin()->PI())*x;
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
  }
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
    *j*=Ait;
    static_cast<Dipole_Color*>(p_cc)->AddJI(j,0);
    *j*=Ajt/Ait;
    static_cast<Dipole_Color*>(p_cc)->AddJI(j,1);
    *j*=2.0/t/Ajt;
    p_cc->AddJ(j);
    p_v->SetZero(false);
  }
#ifdef DEBUG__BG
  p_v->JC()->Print();
#endif
}

template <typename SType>
void VVV_Calculator<SType>::ConstructIDipole()
{
  Current *cj(p_v->JA());
  p_v->Kin()->JIJT()->SetP(cj->P());
  p_v->Kin()->JKT()->SetP(p_v->Kin()->JK()->P());
  const CVec4Type_Matrix *c(cj->J().Get<CVec4Type>());
  const CObject_Vector *cc(&p_v->Kin()->JK()->J().front());
  if (!p_cc->Evaluate(c->front().front(),cc->front())) return;
  double nf(Flavour(kf_quark).Size()/2);
  double d(p_v->Info()->DRMode()?1.0/6.0:0.0);
  double a(p_v->Info()->AMax()), loga(log(a));
  double Vqq(-16.0/9.0-2.0/3.0*(a-1.0-loga)-d);
  double Vgg(50.0/9.0-0.5*sqr(M_PI)-sqr(loga)+11.0/6.0*(a-1.0-loga)-d);
  p_v->Kin()->SetRes(1.0,2);
  p_v->Kin()->SetRes(11.0/6.0-2.0/3.0*0.5/3.0*nf,1);
  p_v->Kin()->SetRes(Vgg+0.5/3.0*Vqq*nf,0);
  if (p_v->Kin()->JK()->Flav().IsGluon()) {
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
    CVec4Type *j((CVec4Type*)(*c)[i].front()->Copy());
    *j*=m_cpl;
    static_cast<Dipole_Color*>(p_cc)->AddJI(j,0);
    static_cast<Dipole_Color*>(p_cc)->AddJI(j,1);
    p_v->AddJ(j);
    p_v->SetZero(false);
  }
#ifdef DEBUG__BG
  p_v->JC()->Print();
#endif
}

template <typename SType>
void VVV_Calculator<SType>::Evaluate()
{
  p_v->SetZero();
  if (p_v->Kin()) {
    if (p_v->Info()->Mode()==1) {
      if (p_v->JA()->Zero()||p_v->JB()->Zero()) return;
      ConstructSDipole();
      return;
    }
    else {
      ConstructIDipole();
      return;
    }
  }
  if (p_v->JA()->Zero()||p_v->JB()->Zero()) return;
#ifdef DEBUG__BG
  msg_Debugging()<<*p_v->JA()<<"(+)"<<*p_v->JB()<<" VVV\n";
  msg_Indent();
#endif
  size_t i(0);
  const CObject_Matrix &cca(p_v->JA()->J()), &ccb(p_v->JB()->J());
  for (typename CObject_Matrix::const_iterator 
	 jait(cca.begin());jait!=cca.end();++jait) {
    const CVec4Type_Vector *ca(jait->Get<CVec4Type>());
    if (ca->empty()) i+=ccb.size(); else
    for (typename CObject_Matrix::const_iterator 
	   jbit(ccb.begin());jbit!=ccb.end();++jbit,++i) {
      const CVec4Type_Vector *cb(jbit->Get<CVec4Type>());
      for (typename CVec4Type_Vector::const_iterator 
	     bit(cb->begin());bit!=cb->end();++bit)
	for (typename CVec4Type_Vector::const_iterator 
	       ait(ca->begin());ait!=ca->end();++ait)
	  if (p_cc->Evaluate(*ait,*bit)) {
#ifdef DEBUG__BG
	    msg_Debugging()<<"  a "<<**ait<<"\n";
	    msg_Debugging()<<"  b "<<**bit<<"\n";
#endif
	    CVec4Type *j(CVec4Type::New(Lorentz(**ait,**bit)*m_cpl));
	    j->SetH(p_v->H(i));
	    j->SetS((*ait)->S()|(*bit)->S());
	    p_cc->AddJ(j);
	    p_v->SetZero(false);
	  }
    }
  }
}

template <typename SType> CVec4<SType>
VVV_Calculator<SType>::Lorentz(const CVec4Type &a,
			       const CVec4Type &b)
{
  CVec4Type j((a*b)*CVec4Type(p_v->JA()->P()-p_v->JB()->P())
	      +(a*ATOOLS::Vec4<SType>(p_v->JB()->P()+p_v->JB()->P()+p_v->JA()->P()))*b
	      -(b*ATOOLS::Vec4<SType>(p_v->JA()->P()+p_v->JA()->P()+p_v->JB()->P()))*a);
  return j;
}

namespace METOOLS {

  template class VVV_Calculator<double>;
  template class VVV_Calculator<long double>;

}

DECLARE_TEMPLATE_GETTER(VVV_Getter,"Gauge3",Lorentz_Calculator,Vertex_Key);

template <typename SType,char STag> Lorentz_Calculator *
VVV_Getter<SType,STag>::operator()(const Vertex_Key &key) const
{
  return new VVV_Calculator<SType>(key);
}

template <typename SType,char STag>
void VVV_Getter<SType,STag>::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"VVV vertex";
}

template class VVV_Getter<double,'D'>;
template class VVV_Getter<long double,'Q'>;
