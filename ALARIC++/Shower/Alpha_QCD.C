#include "ALARIC++/Shower/Alpha_QCD.H"

#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "ALARIC++/Shower/Shower.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Scoped_Settings.H"

#include <algorithm>

#define ZETA3 1.202056903159594

using namespace ALARIC;
using namespace ATOOLS;

Alpha_QCD::Alpha_QCD(const Kernel_Key &key):
  Gauge(key), p_cpl(p_sk->PS()->AlphaS()), m_override(0)
{
  Settings& s = Settings::GetMainSettings();
  m_lc=s["SHOWER"]["CMODE"].Get<unsigned int>();
  m_Nc=s["SHOWER"]["NCOL"].Get<unsigned int>();
  m_CF=(m_Nc*m_Nc-1.)/(2.0*m_Nc);
  m_CA=m_Nc;
  m_TR=1.0/2.0;
}

void Alpha_QCD::SetLimits()
{
  Shower *ps(p_sk->PS());
  m_fac=ps->CplFac((m_type&1)?1:0);
  double scale(ps->TMin((m_type&1)?1:0));
  double scl(CplFac(scale)*scale*p_sk->PS()->MuR2Factor());
  m_max=(*p_cpl)(Max(p_cpl->CutQ2(),scl));
}

double Alpha_QCD::B0(const double &nf) const
{
  return 11.0/6.0*m_CA-2.0/3.0*m_TR*nf;
}

double Alpha_QCD::B1(const double &nf) const
{
  return 17.0/6.0*sqr(m_CA)-(5.0/3.0*m_CA+m_CF)*m_TR*nf;
}

double Alpha_QCD::G2(const double &nf) const
{
  return m_CA*(67./18.-sqr(M_PI)/6.)-10./9.*m_TR*nf;
}

double Alpha_QCD::G3(const double &nf) const
{
  double G3=sqr(m_CA)*(245./6.-134./27.*sqr(M_PI)+11./45.*pow(M_PI,4)+22./3.*ZETA3)
    +m_CA*m_TR*nf*(-418./27.+40./27.*sqr(M_PI)-56./3.*ZETA3)
    +m_CF*m_TR*nf*(-55./3.+16.*ZETA3)-16./27.*sqr(m_TR*nf);
  return G3/4.0;
}

double Alpha_QCD::K(const Splitting &s) const
{
  if (!(s.m_kfac&1)) return 0.0;
  m_override=1;
  double asf=Coupling(s)/(2.0*M_PI), nf=Nf(s);
  m_override=0;
  if (!(s.m_kfac&4)) return asf*G2(nf);
  return asf*G2(nf)+sqr(asf)*G3(nf);
}

double Alpha_QCD::KMax(const Splitting &s) const
{
  if (!(s.m_kfac&1)) return 0.0;
  double asf=CplMax(s)/(2.0*M_PI);
  if (!(s.m_kfac&4)) return asf*G2(3.);
  return asf*G2(3.)+sqr(asf)*G3(3.);
}

double Alpha_QCD::Nf(const Splitting &s) const
{
  return p_cpl->Nf(Scale(s));
}

double Alpha_QCD::CplFac(const double &scale) const
{
  return m_fac;
}

double Alpha_QCD::TrueScale(const Splitting &s) const
{
  double scale(Scale(s));
  return CplFac(scale)*scale;
}

double Alpha_QCD::Coupling(const Splitting &s) const
{
  if (m_override==0) {
    if (s.m_clu&1) return 1.0;
    if (s.m_clu&2) return (*p_cpl)(s.m_t1);
  }
  double murf(p_sk->PS()->MuR2Factor());
  double scl(TrueScale(s)*murf);
  if (scl<murf*p_cpl->CutQ2()) return 0.0;
  double cpl=(*p_cpl)(scl);
  return cpl;
}

double Alpha_QCD::RenCT(const Splitting &s) const
{
  if (m_override==0) {
    if (s.m_clu&3) return 0.0;
  }
  if (!(s.m_kfac&8)) return 0.0;
  double scale(Scale(s)), murf(p_sk->PS()->MuR2Factor());
  double scl(TrueScale(s)*murf);
  if (scl<murf*p_cpl->CutQ2()) return 0.0;
  double cpl=(*p_cpl)(scl);
  double ct=0, ocpl=cpl;
#ifdef DEBUG__AlphaS
  msg_Debugging()<<"t="<<scale<<", \\mu_R^2="<<scl<<std::endl;
  msg_Debugging()<<"as(t)="<<cpl<<std::endl;
#endif
  if (!IsEqual(scl,scale)) {
#ifdef DEBUG__AlphaS
    msg_Debugging()<<"as(\\mu_R^2)="<<cpl<<std::endl;
#endif
    std::vector<double> ths(p_cpl->Thresholds(scale,scl));
    ths.push_back((scl>scale)?scl:scale);
    ths.insert(ths.begin(),(scl>scale)?scale:scl);
    if (scale<scl) std::reverse(ths.begin(),ths.end());
#ifdef DEBUG__AlphaS
    msg_Debugging()<<"thresholds: "<<ths<<std::endl;
#endif
    double fac(1.);
    switch (p_sk->PS()->ScaleVarScheme()) {
    case 1:
    case 5:
    case 9:
      // local counterterms and redefinition at threshold
      for (size_t i(1);i<ths.size();++i) {
        ct=cpl/M_PI*p_cpl->Beta0((ths[i]+ths[i-1])/2.0)*log(ths[i]/ths[i-1]);
	cpl*=1.0-ct;
      }
      ct=1.0-cpl/ocpl;
      break;
    case 2:
    case 6:
    case 10:
      // replace as(t) -> as(t)*[1-sum as/2pi*beta(nf)*log(th[i]/th[i-1])]
      for (size_t i(1);i<ths.size();++i)
        ct+=cpl/M_PI*p_cpl->Beta0((ths[i]+ths[i-1])/2.0)*log(ths[i]/ths[i-1]);
      fac=1.0-ct;
      break;
    default:
      fac=1.;
      break;
    }
#ifdef DEBUG__AlphaS
    msg_Debugging()<<"ct="<<ct<<std::endl;
#endif
    cpl*=fac;
#ifdef DEBUG__AlphaS
    msg_Debugging()<<"as(\\mu_R^2)*(1-ct)="<<cpl
		   <<" vs. "<<(*p_cpl)(scl/murf)<<std::endl;
#endif
    ct=-ct;
  }
  return ct;
}

double Alpha_QCD::CplMax(const Splitting &s) const
{
  return m_max;
}

double Alpha_QCD::Solve(const double &as) const
{
  double t0=p_sk->PS()->TMin(m_type&1);
  t0=Max(p_cpl->CutQ2(),CplFac(t0)*t0);
  double mur2=p_cpl->WDBSolve(as,t0,sqr(rpa->gen.Ecms()));
  msg_Debugging()<<"\\alpha_s("<<sqrt(mur2)<<") = "
		 <<(*p_cpl)(mur2)<<" / "<<as<<"\n";
  return mur2;
}
