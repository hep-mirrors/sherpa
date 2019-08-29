#include "MCATNLO/Showers/Splitting_Function_Base.H"

#include "MCATNLO/Tools/Parton.H"
#include "MODEL/Main/Color_Function.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "PDF/Main/PDF_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "MCATNLO/Showers/Shower.H"

using namespace MCATNLO;
using namespace MODEL;
using namespace ATOOLS;

Splitting_Function_Base::Splitting_Function_Base():
  p_lf(NULL), p_cf(NULL), m_type(cstp::none), m_mth(0.0),
  m_on(1), m_qcd(-1)
{
}

Splitting_Function_Base::Splitting_Function_Base(const SF_Key &key):
  p_lf(NULL), p_cf(NULL), m_type(key.m_type),
  m_symf(1.0), m_polfac(1.0), m_lpdf(1.0), m_mth(0.0),
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
  p_lf = SFL_Getter::GetObject(ckey.p_v->Lorentz[0],ckey);
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
  // exclude all massive partons to split, but also Q->Qg
  Settings& s = Settings::GetMainSettings();
  bool massive_splittings =
    s["MCATNLO_MASSIVE_SPLITTINGS"].SetDefault(1).Get<int>();
  if (!massive_splittings &&
      (p_lf->FlA().IsMassive() || p_lf->FlB().IsMassive() ||
       p_lf->FlC().IsMassive())) {
    m_on=false;
  }
  // do not add a new massive FS (not implemented in as ME dipoles in SRC)
  // IS a->bc (b to hard process): forbid c to be massive
  //                               if b massless, this forbids Q->gQ
  //                               if b is massive, this forbids g->QQ
  //                                                and allows Q->Qg
  // FS a->bc (a from hard process): forbid b massive && b==cbar
  //                                 if a is massless, this forbids g->QQ
  //                                 if a is massive, this allows Q->Qg, Q->gQ
  bool splitintomassive =
    s["MCATNLO_SPLIT_INTO_MASSIVE"].SetDefault(0).Get<int>();
  if (!splitintomassive &&
      (((key.m_type==cstp::IF || key.m_type==cstp::II) &&
        p_lf->FlC().IsMassive()) ||
       ((key.m_type==cstp::FF || key.m_type==cstp::FI) &&
        (p_lf->FlB().IsMassive() && p_lf->FlC()==p_lf->FlB().Bar())))) {
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
   const double _scale,const double Q2,const Color_Info &ci,
   Cluster_Amplitude *const sub)
{
  double scale(_scale);
  if (scale>0.0) scale=p_lf->Scale(z,y,scale,Q2);
  m_lastscale = scale;
  double sf((*p_lf)(z,y,eta,scale,Q2,sub));
  if (p_lf->m_dipole_case == EXTAMP::IDa) sf = dabs(sf);
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
(const double zmin,const double zmax,const double scale,const double xbj, const double phimin)
{
  if (m_mth && (m_type==cstp::FF || m_type==cstp::FI)) {
    if (p_lf->FlA().Mass(true)<m_mth &&
  	sqr(p_lf->FlA().Mass(true))>scale) return 0.0;
    if (p_lf->FlB().Mass(true)<m_mth &&
  	p_lf->FlC().Mass(true)<m_mth &&
  	sqr(p_lf->FlB().Mass(true)+
  	    p_lf->FlC().Mass(true))>scale) return 0.0;
  }
  double lastint = p_lf->OverIntegrated(zmin,zmax,scale,xbj,phimin)/m_symf/m_polfac;
  if (!(IsBad(lastint)||lastint<0.0)) {
    m_lastint+=lastint;
    m_lastints.push_back(m_lastint);
  }
  else {
    msg_Error()<<METHOD<<"(): Integral is "<<lastint<<" in ("<<m_type<<") "
	       <<p_lf->FlA()<<"->"<<p_lf->FlB()<<p_lf->FlC()<<std::endl;
    return 0.0;
  }
  return lastint;
}

double Splitting_Function_Base::Overestimated(const double z,const double y, const double phi)
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
  if(p_lf->m_dipole_case==EXTAMP::IDa) m_lastacceptwgt *= Jacobian(z,vi,phi,alpha,paipb,Q2);
  // TODO: there are m_lastacceptwgt>1
#ifdef CHECK_rejection_weight
  if (m_lastacceptwgt>1.0) {
    msg_Error()<<METHOD<<"(): Weight is "<<m_lastacceptwgt<<" in ("<<m_type<<") "
	       <<p_lf->FlA()<<"->"<<p_lf->FlB()<<p_lf->FlC()
	       <<" at z = "<<z<<", y = "<<y<<", x = "
	       <<eta<<", Q = "<<sqrt(Q2)<<std::endl;
  }
#endif
  return m_lastacceptwgt;
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
  double disc=ran->Get()*m_lastints.back();
  for (size_t i(0);i<m_lastints.size();++i)
    if (m_lastints[i]>=disc) return SetSpec(m_specs[i]);
  THROW(fatal_error, "Spectators ints are not ordered as expected.");
}

void Splitting_Function_Base::ClearSpecs()
{
  m_specs.clear();
  m_lastints.clear();
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

void Splitting_Function_Base::ResetLastInt()
{
  m_lastint=0.0;
}

double Splitting_Function_Base::Phi(const double &phimin) const
{
  switch(p_lf->m_dipole_case){
    case EXTAMP::CS:
      return 2.*M_PI*ATOOLS::ran->Get();
    case EXTAMP::IDa:
      return 2.*atan(tan(phimin/2.) / (1-2.*ATOOLS::ran->Get()));
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

namespace MCATNLO {

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

