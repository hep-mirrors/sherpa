#include "CSSHOWER++/Showers/Splitting_Function_Base.H"

#define COMPILE__Getter_Function
#define PARAMETER_TYPE CSSHOWER::SF_Key
#define OBJECT_TYPE CSSHOWER::SF_Lorentz
#define SORT_CRITERION std::less<std::string>
#include "ATOOLS/Org/Getter_Function.C"

template class Getter_Function
<CSSHOWER::SF_Coupling,CSSHOWER::SF_Key,SORT_CRITERION>;

template class Getter_Function
<void,const MODEL::Model_Base*,SORT_CRITERION>;

#include "MODEL/Interaction_Models/Lorentz_Function.H"
#include "MODEL/Interaction_Models/Color_Function.H"
#include "MODEL/Interaction_Models/Single_Vertex.H"
#include "PDF/Main/PDF_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Shell_Tools.H"

using namespace CSSHOWER;
using namespace MODEL;
using namespace ATOOLS;

double SF_Lorentz::s_pdfcut=1.0e-6;

SF_Lorentz::SF_Lorentz(const SF_Key &key):
  p_ms(key.p_ms), p_cf(key.p_cf) 
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
(const double &a,const double &b,const double &c)
{
  return a*a+b*b+c*c-2.*(a*b+a*c+b*c);
}

SF_Coupling::~SF_Coupling() {}

double SF_Coupling::CplFac(const double &scale) const
{
  return m_cplfac;
}

Splitting_Function_Base::Splitting_Function_Base():
  p_lf(NULL), p_cf(NULL), m_type(cstp::none),
  m_on(1), m_bwon(0), m_qcd(-1)
{
}

Splitting_Function_Base::Splitting_Function_Base(const SF_Key &key):
  p_lf(NULL), p_cf(NULL), m_type(key.m_type),
  m_symf(1.0), m_polfac(1.0), m_on(1), m_bwon(0), m_qcd(-1)
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
  m_bwon=((key.m_type==cstp::FF || key.m_type==cstp::FI)
	  && key.p_v->in[0].Width()) ||
    ((key.m_type==cstp::IF || key.m_type==cstp::II)
     && key.p_v->in[1].Width());
  if (key.p_v->in[1].Mass()>10.0 &&
      key.p_v->in[2].Mass()>10.0) m_on=0;
  if (key.p_v->in[1]==key.p_v->in[2] &&
      (key.m_type==cstp::FF || key.m_type==cstp::FI)) m_symf=2.0;
  m_polfac=key.p_v->in[0].IntSpin()+1;
  if (key.p_v->in[0].IntSpin()==2 && IsZero(key.p_v->in[0].Mass())) m_polfac=2.0;
  msg_Debugging()<<"Init("<<m_on<<") "<<key
		 <<" => ("<<Demangle(typeid(*p_lf).name()).substr(10)
		 <<","<<Demangle(typeid(*p_cf).name()).substr(10)
		 <<"), sf="<<m_symf<<", polfac="<<m_polfac;
}

Splitting_Function_Base::~Splitting_Function_Base()
{
  if (p_lf) delete p_lf;
  if (p_cf) delete p_cf;
}

double Splitting_Function_Base::BWFactor
(const double &z,const double &y,const double &eta,
 const double &scale,const double &Q2,int mode) const
{
  if (!m_bwon) return 1.0;
  switch (m_type) {
  case cstp::FF: {
    double mij2 = sqr(p_lf->MS()->Mass(p_lf->FlA()));
    double mi2 = sqr(p_lf->MS()->Mass(p_lf->FlB()));
    double mj2 = sqr(p_lf->MS()->Mass(p_lf->FlC()));
    double mk2 = sqr(p_lf->MS()->Mass(p_lf->FlSpec()));
    return scale/
      sqrt(sqr(y*(Q2-mk2)+(1.0-y)*(mi2+mj2)-mij2)
	   +mij2*sqr(p_lf->FlA().Width()));
  }
  case cstp::FI: {
    double mij2 = sqr(p_lf->MS()->Mass(p_lf->FlA()));
    double mi2 = sqr(p_lf->MS()->Mass(p_lf->FlB()));
    double mj2 = sqr(p_lf->MS()->Mass(p_lf->FlC()));
    double ma2 = sqr(p_lf->MS()->Mass(p_lf->FlSpec()));
    return scale/
      sqrt(sqr((y*(Q2+ma2)+(mi2+mj2))/(1.0-y)-mij2)
	   +mij2*sqr(p_lf->FlA().Width()));
  }
  case cstp::IF: {
    double mai2 = sqr(p_lf->MS()->Mass(p_lf->FlB()));
    double ma2 = sqr(p_lf->MS()->Mass(p_lf->FlA()));
    double mi2 = sqr(p_lf->MS()->Mass(p_lf->FlC()));
    double mk2 = sqr(p_lf->MS()->Mass(p_lf->FlSpec()));
    return scale/
      sqrt(sqr((-y*(Q2+mk2)+(z-y)*(mi2+ma2))/z-mai2)
	   +mai2*sqr(p_lf->FlB().Width()));
  }
  case cstp::II: {
    double mai2 = sqr(p_lf->MS()->Mass(p_lf->FlB()));
    double ma2 = sqr(p_lf->MS()->Mass(p_lf->FlA()));
    double mi2 = sqr(p_lf->MS()->Mass(p_lf->FlC()));
    double mb2 = sqr(p_lf->MS()->Mass(p_lf->FlSpec()));
    return scale/
      sqrt(sqr((-y*(Q2-mb2)+(z+y)*(mi2+ma2))/z-mai2)
	   +mai2*sqr(p_lf->FlB().Width()));
  }
  case cstp::none: break;
  }
  return 0.0;
}

double Splitting_Function_Base::operator()
  (const double z,const double y,const double eta,
   const double scale,const double Q2,int mode)
{
  return (*p_lf)(z,y,eta,scale,Q2,mode)
    *BWFactor(z,y,eta,scale,Q2,mode)/m_symf/m_polfac;
}

double Splitting_Function_Base::OverIntegrated
(const double zmin,const double zmax,const double scale,const double xbj)
{
  double lastint = p_lf->OverIntegrated(zmin,zmax,scale,xbj)/m_symf/m_polfac;
  m_lastint+=lastint;
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
  return operator()(z,y,eta,scale,Q2)/Overestimated(z,y);
}

Parton *Splitting_Function_Base::SelectSpec() const
{
  if (m_specs.empty()) return NULL;
  double disc=ran.Get()*m_specs.size();
  return m_specs[Min(m_specs.size()-1,(size_t)disc)];
}

void Splitting_Function_Base::ClearSpecs()
{
  m_specs.clear();
}

double Splitting_Function_Base::GetXPDF
(const double &scale,const double &x,const ATOOLS::Flavour &a,
 const int beam,const int mode)
{
  if (p_pdf[beam]==NULL || !p_pdf[beam]->Contains(a)) return 0.0;
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

double SF_Lorentz::JFF(const double &y) const
{ 
  return (1.-y);
}

double SF_Lorentz::JFI(const double &y,const double &eta,
		       const double &scale) const
{ 
  double fresh = p_sf->GetXPDF(scale,eta/(1.0-y),m_flspec,m_beam);
  double old = p_sf->GetXPDF(scale,eta,m_flspec,m_beam);
  if (fresh<0.0 || old<0.0 || IsZero(old,s_pdfcut) || IsZero(fresh,s_pdfcut)) return 0.; 
  return (1.0-y) * fresh/old;
}

double SF_Lorentz::JIF(const double &z,const double &y,
		       const double &eta,const double &scale) const
{ 
  double fresh = p_sf->GetXPDF(scale,eta/z,m_flavs[0],m_beam);
  double old = p_sf->GetXPDF(scale,eta,m_flavs[1],m_beam);
  if (fresh<0.0 || old<0.0 || IsZero(old,s_pdfcut) || IsZero(fresh,s_pdfcut)) return 0.; 
  return fresh/old;
}

double SF_Lorentz::JII(const double &z,const double &y,
		       const double &eta,const double &scale) const
{ 
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
  return 2.*M_PI*ATOOLS::ran.Get();
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

namespace CSSHOWER {

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

