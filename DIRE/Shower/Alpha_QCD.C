#include "DIRE/Shower/Alpha_QCD.H"

#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "DIRE/Shower/Shower.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"

#include <algorithm>

#define ZETA3 1.202056903159594

using namespace DIRE;
using namespace ATOOLS;

Alpha_QCD::Alpha_QCD(const Kernel_Key &key):
  Gauge(key)
{
  p_cpl=(MODEL::Running_AlphaS*)
    p_sk->PS()->Model()->GetScalarFunction("alpha_S");
  m_lc=key.p_rd->GetValue<unsigned int>("CSS_CMODE",1);
  m_Nc=key.p_rd->GetValue<unsigned int>("CSS_NCOL",3);
  m_CF=(m_Nc*m_Nc-1.)/(2.0*m_Nc);
  m_CA=m_Nc;
  m_TR=1.0/2.0;
}

void Alpha_QCD::SetLimits()
{
  Shower *ps(p_sk->PS());
  m_fac=(m_type&1)?ps->CplFac(1):ps->CplFac(0);
  double scale=(m_type&1)?ps->TMin(1):ps->TMin(0);
  double scl(Max(p_cpl->CutQ2(),CplFac(scale)*scale*p_sk->PS()->MuR2Factor()));
  double ct(0.);
  if (p_sk->PS()->MuR2Factor()>1.) // only for f>1 cpl gets larger
    ct=-(*p_cpl)(scl)/(2.*M_PI)*B0(0)*log(p_sk->PS()->MuR2Factor());
  m_max=(*p_cpl)(scl)*(1.-ct);
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
  double asf=Coupling(s)/(2.0*M_PI), nf=Nf(s);
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

double Alpha_QCD::Coupling(const Splitting &s) const
{
  DEBUG_FUNC("");
  if (s.m_clu&1) return 1.0;
  double scale(Scale(s)), mur2f(p_sk->PS()->MuR2Factor());
  double t(CplFac(scale)*scale), scl(CplFac(scale)*scale*mur2f);
  double cpl=(*p_cpl)(scl);
  msg_Debugging()<<"t="<<t<<", \\mu_R^2="<<scl<<std::endl;
  msg_Debugging()<<"as(t)="<<(*p_cpl)(t)<<std::endl;
  if (!IsEqual(scl,t)) {
    msg_Debugging()<<"as(\\mu_R^2)="<<cpl<<std::endl;
    std::vector<double> ths(p_cpl->Thresholds(t,scl));
    ths.push_back((scl>t)?scl:t);
    ths.insert(ths.begin(),(scl>t)?t:scl);
    if (t<scl) std::reverse(ths.begin(),ths.end());
    msg_Debugging()<<"thresholds: "<<ths<<std::endl;
    double fac(1.),ct(0.);
    switch (p_sk->PS()->ScaleVariationScheme()) {
    case 1:
      // replace as(t) -> as(t)*prod[1-as/2pi*beta(nf)*log(th[i]/th[i-1])]
      for (size_t i(1);i<ths.size();++i) {
        ct=cpl/(2.*M_PI)*B0(p_cpl->Nf((ths[i]+ths[i-1])/2.0))*log(ths[i]/ths[i-1]);
        fac*=1.0-ct;
      }
      break;
    case 2:
      // replace as(t) -> as(t)*[1-sum as/2pi*beta(nf)*log(th[i]/th[i-1])]
      for (size_t i(1);i<ths.size();++i)
        ct+=cpl/(2.*M_PI)*B0(p_cpl->Nf((ths[i]+ths[i-1])/2.0))*log(ths[i]/ths[i-1]);
      fac=1.-ct;
      break;
    default:
      fac=1.;
      break;
    }
    msg_Debugging()<<"ct="<<ct<<std::endl;
    if (fac<0.) {
      msg_Tracking()<<METHOD<<"(): Renormalisation term too large. Remove."
                    <<std::endl;
      fac=1.;
    }
    cpl*=fac;
    msg_Debugging()<<"as(\\mu_R^2)*(1-ct)="<<cpl<<std::endl;
  }
  if (cpl>m_max) {
    msg_Error()<<METHOD<<"(): Value exceeds maximum at t = "
               <<sqrt(t)<<" -> \\mu_R = "<<sqrt(scl)
               <<", qmin = "<<sqrt(p_cpl->CutQ2())<<std::endl;
  }
  return cpl;
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
