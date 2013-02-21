#include "MCATNLO/Showers/Splitting_Function_Base.H"

#define COMPILE__Getter_Function
#define PARAMETER_TYPE MCATNLO::SF_Key
#define OBJECT_TYPE MCATNLO::SF_Lorentz
#define SORT_CRITERION std::less<std::string>
#include "ATOOLS/Org/Getter_Function.C"

template class Getter_Function
<MCATNLO::SF_Coupling,MCATNLO::SF_Key,SORT_CRITERION>;

template class Getter_Function
<void,MCATNLO::SFC_Filler_Key,SORT_CRITERION>;

#include "MCATNLO/Tools/Parton.H"
#include "MODEL/Interaction_Models/Lorentz_Function.H"
#include "MODEL/Interaction_Models/Color_Function.H"
#include "MODEL/Interaction_Models/Single_Vertex.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "PDF/Main/PDF_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "MCATNLO/Showers/Shower.H"

using namespace MCATNLO;
using namespace MODEL;
using namespace ATOOLS;

double SF_Lorentz::s_pdfcut=1.0e-6;
double SF_Lorentz::s_kappa=2.0/3.0;

SF_Lorentz::SF_Lorentz(const SF_Key &key):
  p_ms(key.p_ms), p_cf(key.p_cf), m_col(0)
{
  m_flavs[0]=key.p_v->in[0];
  if (key.m_mode==0) {
    m_flavs[1]=key.p_v->in[1];
    m_flavs[2]=key.p_v->in[2];
  }
  else {
    m_flavs[1]=key.p_v->in[2];
    m_flavs[2]=key.p_v->in[1];
  }
}

SF_Lorentz::~SF_Lorentz() {}

double SF_Lorentz::Lambda
(const double &a,const double &b,const double &c) const
{
  return a*a+b*b+c*c-2.*(a*b+a*c+b*c);
}

double SF_Lorentz::AsymmetryFactor
(const double z,const double y,const double Q2)
{
  return 1.0;
}

double SF_Coupling::s_qfac=1.0;

SF_Coupling::SF_Coupling(const SF_Key &key):
  p_lf(NULL), m_type(key.m_type),
  m_cplfac(1.0), m_kfmode(key.m_kfmode) 
{
}

SF_Coupling::~SF_Coupling() {}

double SF_Coupling::CplFac(const double &scale) const
{
  return m_cplfac;
}

void SF_Coupling::ColorPoint(Parton *const p) const
{
}

double SF_Coupling::ColorWeight(const Color_Info &ci) const
{
  return 1.0;
}

Splitting_Function_Base::Splitting_Function_Base():
  p_lf(NULL), p_cf(NULL), m_type(cstp::none),
  m_on(1), m_qcd(-1)
{
}

Splitting_Function_Base::Splitting_Function_Base(const SF_Key &key):
  p_lf(NULL), p_cf(NULL), m_type(key.m_type),
  m_symf(1.0), m_polfac(1.0), m_lpdf(1.0),
  m_on(1), m_qcd(-1)
{
  SF_Key ckey(key);
  ckey.p_cf=p_cf = SFC_Getter::GetObject(ckey.ID(0),ckey);
  if (p_cf==NULL) {
    ckey.p_cf=p_cf = SFC_Getter::GetObject(ckey.ID(1),ckey);
    if (p_cf==NULL) {
      m_on=-1;
      return;
    }
  }
  p_lf = SFL_Getter::GetObject(ckey.p_v->Lorentz[0]->Type(),ckey);
  if (p_lf==NULL) {
    m_on=-1;
    return;
  }
  p_cf->SetLF(p_lf);
  p_lf->SetSF(this);
  m_qcd=p_lf->FlA().Strong()&&p_lf->FlB().Strong()&&p_lf->FlC().Strong();
  m_on=PureQCD();// so far only qcd evolution
  if (!m_on && (ckey.m_ewmode&1) &&
      (p_lf->FlA().IsPhoton() || p_lf->FlB().IsPhoton() ||
       p_lf->FlC().IsPhoton())) m_on=true;
  Data_Reader read(" ",";","#","=");
  bool massive_splittings=read.GetValue<int>("MCATNLO_MASSIVE_SPLITTINGS",1);
  if (!massive_splittings &&
      (p_lf->FlA().IsMassive() || p_lf->FlB().IsMassive() ||
       p_lf->FlC().IsMassive())) {
    m_on=false;
  }
  if (key.p_v->in[1].Mass()>10.0 &&
      key.p_v->in[2].Mass()>10.0) m_on=0;
  if (key.p_v->in[1]==key.p_v->in[2] &&
      (key.m_type==cstp::FF || key.m_type==cstp::FI)) m_symf=2.0;
  m_polfac=key.p_v->in[0].IntSpin()+1;
  if (key.p_v->in[0].IntSpin()==2 && IsZero(key.p_v->in[0].Mass())) m_polfac=2.0;
  msg_Debugging()<<"Init("<<m_on<<") "<<p_lf->FlA()<<"->"
		 <<p_lf->FlB()<<","<<p_lf->FlC()
		 <<" => ("<<Demangle(typeid(*p_lf).name()).substr(8)
		 <<","<<Demangle(typeid(*p_cf).name()).substr(8)
		 <<"), sf="<<m_symf<<", polfac="<<m_polfac
		 <<", col="<<p_lf->Col();
}

Splitting_Function_Base::~Splitting_Function_Base()
{
  if (p_lf) delete p_lf;
  if (p_cf) delete p_cf;
}

double Splitting_Function_Base::MEPSWeight
(const double &z,const double &y,const double &eta,
 const double &scale,const double &Q2,const Cluster_Amplitude *const sub) const
{
  double ma2(p_lf->MS()->Mass2(p_lf->FlA())), mb2(p_lf->MS()->Mass2(p_lf->FlB()));
  double mk2(p_lf->MS()->Mass2(p_lf->FlSpec())), mc2(p_lf->MS()->Mass2(p_lf->FlC()));
  switch (m_type) {
  case cstp::FF:
    return (8.0*M_PI)/(Q2*y)/p_lf->JFF(y,mb2/Q2,mc2/Q2,mk2/Q2,ma2/Q2);
  case cstp::FI:
    return (8.0*M_PI)/((Q2+mb2+mc2)*y)/p_lf->JFI(y,eta,scale,sub);
  case cstp::IF:
    return (8.0*M_PI)/((Q2+mk2)*y)/p_lf->JIF(z,y,eta,scale,sub);
  case cstp::II:
    return (8.0*M_PI)/(Q2*y)/p_lf->JII(z,y,eta,scale,sub);
  case cstp::none: break;
  }
  return 0.0;
}

double Splitting_Function_Base::operator()
  (const double z,const double y,const double eta,
   const double scale,const double Q2,const Color_Info &ci,
   Cluster_Amplitude *const sub)
{
  double sf((*p_lf)(z,y,eta,scale,Q2,sub));
  if (sf/p_lf->AsymmetryFactor(z,y,Q2)<0.0) return 0.0;
  if (sf<0.0 && sub==NULL) sf=-sf; 
  if (IsBad(sf)) {
    PRINT_INFO("Invalid weight in MC@NLO "+
               Demangle(std::string(typeid(*p_lf).name()).substr(9))+"|"+
               Demangle(std::string(typeid(*p_cf).name()).substr(9)));
    return 0.0;
  }
  return sf*(ci.m_new<0?1.0:p_cf->ColorWeight(ci))/m_symf/m_polfac;
}

double Splitting_Function_Base::AsymmetryFactor
(const double z,const double y,const double Q2)
{
  return p_lf->AsymmetryFactor(z,y,Q2);
}

double Splitting_Function_Base::OverIntegrated
(const double zmin,const double zmax,const double scale,const double xbj)
{
  double lastint = p_lf->OverIntegrated(zmin,zmax,scale,xbj)/m_symf/m_polfac;
  if (!(IsBad(lastint)||lastint<0.0)) {
    m_lastint+=lastint;
  }
  else {
    msg_Error()<<METHOD<<"(): Integral is "<<lastint<<" in ("<<m_type<<") "
	       <<p_lf->FlA()<<"->"<<p_lf->FlB()<<p_lf->FlC()<<std::endl;
    return 0.0;
  }
  return lastint;
}

double Splitting_Function_Base::Overestimated(const double z,const double y)
{
  return p_lf->OverEstimated(z,y)/m_symf/m_polfac;
}

double Splitting_Function_Base::Z()
{
  return p_lf->Z();
}
        
double Splitting_Function_Base::RejectionWeight
(const double z,const double y,const double eta,
 const double scale,const double Q2) 
{
  double res = operator()(z,y,eta,scale,Q2)/Overestimated(z,y);
#ifdef CHECK_rejection_weight
  if (res>1.0) {
    msg_Error()<<METHOD<<"(): Weight is "<<res<<" in ("<<m_type<<") "
	       <<p_lf->FlA()<<"->"<<p_lf->FlB()<<p_lf->FlC()
	       <<" at z = "<<z<<", y = "<<y<<", x = "
	       <<eta<<", Q = "<<sqrt(Q2)<<std::endl;
  }
#endif
  return res;
}

void Splitting_Function_Base::ColorPoint(Parton *const p) const
{
  p_cf->ColorPoint(p);
}
        
Parton *Splitting_Function_Base::SetSpec(Parton *const spec)
{
  SetFlavourSpec(spec->GetFlavour());
  return p_spec=spec;
}

Parton *Splitting_Function_Base::SelectSpec()
{
  if (m_specs.empty()) return NULL;
  double disc=ran->Get()*m_specs.size();
  return SetSpec(m_specs[Min(m_specs.size()-1,(size_t)disc)]);
}

void Splitting_Function_Base::ClearSpecs()
{
  m_specs.clear();
}

double Splitting_Function_Base::GetXPDF
(const double &scale,const double &x,const ATOOLS::Flavour &a,
 const int beam,const int mode)
{
  if (p_pdf[beam]==NULL) return 1.0;
  if (!p_pdf[beam]->Contains(a)) {
    if (a.Strong() || a.Mass()<10.0) return 0.0;
    return 1.0;
  }
  if (mode==1) return m_lpdf==-1.0?0.0:p_pdf[beam]->GetXPDF(a);
  if (IsNan(scale) || IsNan(x)) return 0.0;
  double Q2(scale*p_cf->CplFac(scale));
  if (Q2<p_lf->MS()->Mass2(a) || x<p_pdf[beam]->XMin() ||
      x>p_pdf[beam]->XMax()*p_pdf[beam]->RescaleFactor() ||
      Q2<p_pdf[beam]->Q2Min() || Q2>p_pdf[beam]->Q2Max())
    return m_lpdf=-1.0;
  p_pdf[beam]->Calculate(x,Q2);
  return m_lpdf=p_pdf[beam]->GetXPDF(a);
}

double SF_Lorentz::JFF(const double &y,const double &mui2,
		       const double &muj2,const double &muk2,
		       const double &muij2) const
{ 
  return (1.-y)*sqr(1.0-mui2-muj2-muk2)/sqrt(Lambda(1.0,muij2,muk2));
}

double SF_Lorentz::JFI(const double &y,const double &eta,
		       const double &scale,const Cluster_Amplitude *const sub) const
{ 
  if (sub) return 1.0;
  double fresh = p_sf->GetXPDF(scale,eta/(1.0-y),m_flspec,m_beam);
  double old = p_sf->GetXPDF(scale,eta,m_flspec,m_beam);
  if (fresh<0.0 || old<0.0 || IsZero(old,s_pdfcut) || IsZero(fresh,s_pdfcut)) return 0.; 
  return (1.0-y) * fresh/old;
}

double SF_Lorentz::JIF(const double &z,const double &y,const double &eta,
		       const double &scale,const Cluster_Amplitude *const sub) const
{ 
  if (sub) return 1.0/z;
  double fresh = p_sf->GetXPDF(scale,eta/z,m_flavs[0],m_beam);
  double old = p_sf->GetXPDF(scale,eta,m_flavs[1],m_beam);
  if (fresh<0.0 || old<0.0 || IsZero(old,s_pdfcut) || IsZero(fresh,s_pdfcut)) return 0.; 
  return fresh/old;
}

double SF_Lorentz::JII(const double &z,const double &y,const double &eta,
		       const double &scale,const Cluster_Amplitude *const sub) const
{ 
  if (sub) return 1.0/z;
  double fresh = p_sf->GetXPDF(scale,eta/z,m_flavs[0],m_beam);
  double old = p_sf->GetXPDF(scale,eta,m_flavs[1],m_beam);
  if (fresh<0.0 || old<0.0 || IsZero(old,s_pdfcut) || IsZero(fresh,s_pdfcut)) return 0.; 
  return fresh/old;
}

void Splitting_Function_Base::ResetLastInt()
{
  m_lastint=0.0;
}

double Splitting_Function_Base::Phi(double z) const
{
  return 2.*M_PI*ATOOLS::ran->Get();
}

const Flavour & Splitting_Function_Base::GetFlavourA() const
{
  return p_lf->FlA();
}

const Flavour & Splitting_Function_Base::GetFlavourB() const
{
  return p_lf->FlB();
}

const Flavour & Splitting_Function_Base::GetFlavourC() const
{
  return p_lf->FlC();
}

const Flavour & Splitting_Function_Base::GetFlavourSpec() const
{
  return p_lf->FlSpec();
}

int Splitting_Function_Base::GetCol() const
{
  return p_lf->Col();
}

bool Splitting_Function_Base::PureQCD() const
{ 
  if (m_qcd<0) THROW(fatal_error,"Invalid request");
  return m_qcd;
}

std::string SF_Key::ID(const int mode) const
{
  if ((m_mode==1)^(mode==1))
    return "{"+ToString(p_v->in[0])+"}{"
      +ToString(p_v->in[2])+"}{"+ToString(p_v->in[1])+"}";
  return "{"+ToString(p_v->in[0])+"}{"
    +ToString(p_v->in[1])+"}{"+ToString(p_v->in[2])+"}";
}

namespace MCATNLO {

  std::ostream &operator<<(std::ostream &str,const SF_Key &k)
  {
    if (k.m_mode==0) 
      return str<<k.m_type<<" "<<k.p_v->in[0]<<"->"<<k.p_v->in[1]<<","<<k.p_v->in[2];
    return str<<k.m_type<<" "<<k.p_v->in[0]<<"->"<<k.p_v->in[2]<<","<<k.p_v->in[1];
  }

  std::ostream& operator<<(std::ostream& str, const Splitting_Function_Base &base) {
    str<<"  "<<base.GetFlavourA()<<" -> "<<base.GetFlavourB()<<" + "<<base.GetFlavourC()
       <<" : "<<base.m_lastint<<std::endl;
    return str;
  }

}

