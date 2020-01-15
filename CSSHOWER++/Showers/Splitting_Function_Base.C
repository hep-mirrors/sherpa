#include "CSSHOWER++/Showers/Splitting_Function_Base.H"

#include "CSSHOWER++/Tools/Parton.H"
#include "MODEL/Main/Color_Function.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "PDF/Main/PDF_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "CSSHOWER++/Showers/Shower.H"

#include <algorithm>

using namespace CSSHOWER;
using namespace MODEL;
using namespace ATOOLS;

Splitting_Function_Base::Splitting_Function_Base():
  p_lf(NULL), p_cf(NULL), m_type(cstp::none), m_mth(0.0),
  m_on(1), m_qcd(-1)
{
  Settings& s = Settings::GetMainSettings();
  m_evol = s["CSS_EVOLUTION_SCHEME"].Get<int>();
}

void Splitting_Function_Base::SetEFac(Shower *const shower)
{
  std::string id="S{"+ToString(p_lf->FlA())+"}{"
    +ToString(p_lf->FlB())+"}{"+ToString(p_lf->FlC())+"}";
  m_efac=shower->EFac(id);
  if (m_efac!=1.0) msg_Info()<<"Enhance "<<id<<" with "<<m_efac<<"\n";
  if (m_efac==0.0) m_on=0;
}

SF_Lorentz* Splitting_Function_Base::InitLorentzCalc(const MODEL::Single_Vertex& vertex,
						     const SF_Key& sf_key)
{
  std::vector<int> spins;
  for(const auto& fl : vertex.in) 
    spins.push_back(fl.IntSpin());
  std::sort(spins.begin(),spins.end());

  SF_Lorentz* lf(NULL);
  if(spins[0]==0 && spins[1]==0 && spins[2]==2)
    lf = SFL_Getter::GetObject("SSV",sf_key);
  else if(spins[0]==0 && spins[1]==2 && spins[2]==2)
    lf = SFL_Getter::GetObject("HVV",sf_key);
  else if(spins[0]==1 && spins[1]==1 && spins[2]==2)
    lf = SFL_Getter::GetObject("FFV1",sf_key);
  else if(spins[0]==2 && spins[1]==2 && spins[2]==2)
    lf = SFL_Getter::GetObject("VVV",sf_key);

  return lf;
}

Splitting_Function_Base::Splitting_Function_Base(const SF_Key &key):
  p_lf(NULL), p_cf(NULL), m_type(key.m_type),
  m_symf(1.0), m_polfac(1.0), m_lpdf(1.0), m_efac(1.0), m_mth(0.0),
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
  p_lf = InitLorentzCalc(*ckey.p_v, ckey);
  if (p_lf==NULL) {
    m_on=-1;
    return;
  }
  p_cf->SetLF(p_lf);
  p_lf->SetSF(this);
  m_qcd=p_lf->FlA().Strong()&&p_lf->FlB().Strong()&&p_lf->FlC().Strong();
  m_on=PureQCD()&&(ckey.m_qcdmode&1);// so far only qcd evolution
  if (!m_on && (ckey.m_ewmode&1) &&
      (p_lf->FlA().IsPhoton() || p_lf->FlB().IsPhoton() ||
       p_lf->FlC().IsPhoton())) m_on=true;
  if (key.p_v->in[1].Mass()>10.0 &&
      key.p_v->in[2].Mass()>10.0) m_on=0;
  if (key.p_v->in[1]==key.p_v->in[2] &&
      (key.m_type==cstp::FF || key.m_type==cstp::FI)) m_symf=2.0;
  m_polfac=key.p_v->in[0].IntSpin()+1;
  if (key.p_v->in[0].IntSpin()==2 && IsZero(key.p_v->in[0].Mass())) m_polfac=2.0;
  msg_Debugging()<<"Init("<<m_on<<") "<<p_lf->FlA()<<"->"
		 <<p_lf->FlB()<<","<<p_lf->FlC()
		 <<" => ("<<Demangle(typeid(*p_lf).name()).substr(10)
		 <<","<<Demangle(typeid(*p_cf).name()).substr(10)
		 <<"), sf="<<m_symf<<", polfac="<<m_polfac
		 <<", col="<<p_lf->Col();
  Settings& s = Settings::GetMainSettings();
  m_evol = s["CSS_EVOLUTION_SCHEME"].Get<int>();
}

Splitting_Function_Base::~Splitting_Function_Base()
{
  if (p_lf) delete p_lf;
  if (p_cf) delete p_cf;
}

double Splitting_Function_Base::MEPSWeight
(const double &z,const double &y,const double &eta,
 const double &scale,const double &Q2) const
{
  double ma2(p_lf->MS()->Mass2(p_lf->FlA())), mb2(p_lf->MS()->Mass2(p_lf->FlB()));
  double mk2(p_lf->MS()->Mass2(p_lf->FlSpec())), mc2(p_lf->MS()->Mass2(p_lf->FlC()));
  switch (m_type) {
  case cstp::FF:
    return (8.0*M_PI)/(Q2*y)/p_lf->JFF(y,mb2/Q2,mc2/Q2,mk2/Q2,ma2/Q2);
  case cstp::FI:
    return (8.0*M_PI)/((Q2+mb2+mc2)*y)/p_lf->JFI(y,eta,scale);
  case cstp::IF:
    return (8.0*M_PI)/((Q2+mk2)*y)/p_lf->JIF(z,y,eta,scale);
  case cstp::II:
    return (8.0*M_PI)/(Q2*y)/p_lf->JII(z,y,eta,scale);
  case cstp::none: break;
  }
  return 0.0;
}

double Splitting_Function_Base::operator()
  (const double z,const double y,const double eta,
   const double _scale,const double Q2)
{
  double scale(_scale);
  if (scale>0.0) scale=p_lf->Scale(z,y,scale,Q2);
  m_lastscale = scale;
  double sf((*p_lf)(z,y,eta,scale,Q2)/m_symf/m_polfac);
  if (IsBad(sf)) {
    PRINT_INFO("Invalid weight in CSS "+
               Demangle(std::string(typeid(*p_lf).name()).substr(12))+"|"+
               Demangle(std::string(typeid(*p_cf).name()).substr(11)));
    return 0.0;
  }
  return Max(0.0,sf);
}

double Splitting_Function_Base::OverIntegrated
(const double zmin,const double zmax,const double scale,const double xbj)
{
  if (m_mth && (m_type==cstp::FF || m_type==cstp::FI)) {
    if (p_lf->FlA().Mass(true)<m_mth &&
	sqr(p_lf->FlA().Mass(true))>scale) return 0.0;
    if (p_lf->FlB().Mass(true)<m_mth &&
	p_lf->FlC().Mass(true)<m_mth &&
	sqr(p_lf->FlB().Mass(true)+
	    p_lf->FlC().Mass(true))>scale) return 0.0;
  }
  double lastint = p_lf->OverIntegrated(zmin,zmax,scale,xbj)/m_symf/m_polfac;
  if (!(IsBad(lastint)||lastint<0.0)) {
    if (m_efac!=1.0) lastint*=m_efac;
    m_lastint+=lastint;
  }
  else {
    msg_Error()<<METHOD<<"(): Integral is "<<lastint<<" in ("<<m_type<<") "
	       <<p_lf->FlA()<<"->"<<p_lf->FlB()<<p_lf->FlC()<<std::endl;
    return 0.0;
  }
  return lastint;
}

double Splitting_Function_Base::EFac() const
{
  if (m_efac!=1.0) return m_efac; 
  return 1.0;
}

double Splitting_Function_Base::Overestimated(const double z,const double y,const double phi)
{
  return p_lf->OverEstimated(z,y,phi)/m_symf/m_polfac;
}

double Splitting_Function_Base::Z()
{
  return p_lf->Z();
}
        
double Splitting_Function_Base::Jacobian(const double &z, const double &vi, const double &phi,
       const double &alpha, const double &paipb, const double &Q2) const
{
  const double mw2 = sqr(Flavour(24).Mass());
  const double Qprime2 = Q2-mw2;
  const double A = alpha*vi*Qprime2/2.;
  const double B = paipb*(vi*(z-1-mw2/Qprime2)+1-z);
  const double kt2 = vi*Qprime2/(2.*paipb)*(A+B-2.*sqrt(A*B)*cos(phi));

  const double jacobian = 1./(kt2/vi + vi*Qprime2/(2.*paipb)*(alpha*Qprime2/2.
               +paipb*(z-1-mw2/Qprime2)-2.*cos(phi)/(2.*sqrt(A*B))
               *(alpha*Qprime2/2.*B + A*paipb*(z-1-mw2/Qprime2))));
  return dabs(jacobian*kt2/vi);
}

double Splitting_Function_Base::RejectionWeight
(const double z,const double y,const double eta,
 const double _scale,const double Q2,const double phi, const double vi,
 const double alpha, const double paipb)
{
  double scale(_scale);
  if (scale>0.0) scale=p_lf->Scale(z,y,scale,Q2);
  m_lastacceptwgt = operator()(z,y,eta,scale,Q2)/Overestimated(z,y,phi);
  if(p_lf->m_dipole_case==EXTAMP::IDa && m_evol==1)
    m_lastacceptwgt *= Jacobian(z,vi,phi,alpha,paipb,Q2);
#ifdef CHECK_rejection_weight
  if (m_lastacceptwgt > 1.0) {
    msg_Error()<<METHOD<<"(): Weight is "<<m_lastacceptwgt<<" in ("<<m_type<<") "
	       <<p_lf->FlA()<<"->"<<p_lf->FlB()<<p_lf->FlC()
	       <<" at z = "<<z<<", y = "<<y<<", x = "
	       <<eta<<", Q = "<<sqrt(Q2)<<std::endl;
  }
#endif
  return m_lastacceptwgt;
}

void Splitting_Function_Base::SetMS(const ATOOLS::Mass_Selector *const ms)
{
  p_lf->SetMS(ms);
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
  double Q2(scale);
  if (Q2<p_lf->MS()->Mass2(a) || x<p_pdf[beam]->XMin() ||
      x>p_pdf[beam]->XMax()*p_pdf[beam]->RescaleFactor() ||
      Q2<p_pdf[beam]->Q2Min() || Q2>p_pdf[beam]->Q2Max())
    return m_lpdf=-1.0;
  p_pdf[beam]->Calculate(x,Q2);
  return m_lpdf=p_pdf[beam]->GetXPDF(a);
}

bool Splitting_Function_Base::CheckPDF
(const double &x,const ATOOLS::Flavour &a,const int beam)
{
  if (p_pdf[beam]==NULL) return true;
  return x<=p_pdf[beam]->XMax()*p_pdf[beam]->RescaleFactor();
}

void Splitting_Function_Base::ResetLastInt()
{
  m_lastint=0.0;
}

double Splitting_Function_Base::Phi() const
{
  switch(p_lf->m_dipole_case){
    case EXTAMP::CS:
      return 2.*M_PI*ATOOLS::ran->Get();
    case EXTAMP::IDa:
      return 2.*atan(sqrt(1.-sqr(m_K))/(1+m_K)*tan(M_PI*ran->Get()));
  }
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
    return "{"+ToString(p_v->in[0].Bar())+"}{"
      +ToString(p_v->in[2])+"}{"+ToString(p_v->in[1])+"}";
  return "{"+ToString(p_v->in[0].Bar())+"}{"
    +ToString(p_v->in[1])+"}{"+ToString(p_v->in[2])+"}";
}

namespace CSSHOWER {

  std::ostream &operator<<(std::ostream &str,const SF_Key &k)
  {
    if (k.m_mode==0) 
      return str<<k.m_type<<" "<<k.p_v->in[0].Bar()<<"->"<<k.p_v->in[1]<<","<<k.p_v->in[2];
    return str<<k.m_type<<" "<<k.p_v->in[0].Bar()<<"->"<<k.p_v->in[2]<<","<<k.p_v->in[1];
  }

  std::ostream& operator<<(std::ostream& str, const Splitting_Function_Base &base) {
    str<<"  "<<base.GetFlavourA()<<" -> "<<base.GetFlavourB()<<" + "<<base.GetFlavourC()
       <<" : "<<base.m_lastint<<std::endl;
    return str;
  }

}

