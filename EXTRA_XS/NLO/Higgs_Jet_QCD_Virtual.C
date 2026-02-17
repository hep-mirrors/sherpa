#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace PHASIC;
using namespace ATOOLS;

namespace EXTRAXS {

  const double NC(3.0);
  const double CA(3.0);
  const double CF(0.5*(NC*NC-1.0)/NC);
  const double TR(0.5);

  class GGHG_QCD_Virtual : public Virtual_ME2_Base { //gg->hg virtual
    bool   m_flip, m_con;
    double m_b0, m_nlf;
  public:
    GGHG_QCD_Virtual(const Process_Info& pi, const Flavour_Vector& flavs,
                     bool flip, bool con) :
      Virtual_ME2_Base(pi, flavs), m_flip(flip), m_con(con), m_b0(0.), m_nlf(0.)
    {
      Flavour lq(kf_quark);
      double nlf(0.);
      for (size_t i(0); i<lq.Size(); ++i) if (!lq[i].IsMassive()) nlf++;
      m_nlf=nlf/2.;
      m_b0 = (NC*11.-2.*m_nlf)/6.;
    }

    ~GGHG_QCD_Virtual() {}

    void Calc(const ATOOLS::Vec4D_Vector& momenta);

    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom) {
      return 4.*M_PI;
    }
  };

  class QQHG_QCD_Virtual : public Virtual_ME2_Base { //qq->hg virtual
    bool   m_flip, m_con;
    double m_b0, m_nlf;
  public:
    QQHG_QCD_Virtual(const Process_Info& pi, const Flavour_Vector& flavs,
                     bool flip, bool con) :
      Virtual_ME2_Base(pi, flavs), m_flip(flip), m_con(con), m_b0(0.), m_nlf(0.)
    {
      Flavour lq(kf_quark);
      double nlf(0.);
      for (size_t i(0); i<lq.Size(); ++i) if (!lq[i].IsMassive()) nlf++;
      m_nlf=nlf/2.;
      m_b0 = (NC*11.-2.*m_nlf)/6.;
    }

    ~QQHG_QCD_Virtual() {}

    void Calc(const ATOOLS::Vec4D_Vector& momenta);

    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom) {
      return 4.*M_PI;
    }
  };

  class GQHQ_QCD_Virtual : public Virtual_ME2_Base { //gq->hq virtual
    bool   m_flip, m_con;
    double m_b0, m_nlf;
  public:
    GQHQ_QCD_Virtual(const Process_Info& pi, const Flavour_Vector& flavs,
                     bool flip, bool con) :
      Virtual_ME2_Base(pi, flavs), m_flip(flip), m_con(con), m_b0(0.), m_nlf(0.)
    {
      Flavour lq(kf_quark);
      double nlf(0.);
      for (size_t i(0); i<lq.Size(); ++i) if (!lq[i].IsMassive()) nlf++;
      m_nlf=nlf/2.;
      m_b0 = (NC*11.-2.*m_nlf)/6.;
    }

    ~GQHQ_QCD_Virtual() {}

    void Calc(const ATOOLS::Vec4D_Vector& momenta);

    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom) {
      return 4.*M_PI;
    }
  };

  class HGQQ_QCD_Virtual : public Virtual_ME2_Base{ // h->gqq virtual
    bool m_flip, m_con;
    double m_b0, m_nlf;

  public:
    HGQQ_QCD_Virtual(const Process_Info &pi, const Flavour_Vector &flavs,
                     bool flip, bool con) : Virtual_ME2_Base(pi, flavs), m_flip(flip), m_con(con), m_b0(0.), m_nlf(0.)
    {
      msg_Out() << "HGQQ_QCD_Virtual(): Using h->gqq virtual matrix element." << std::endl;
      Flavour lq(kf_quark);
      double nlf(0.);
      for (size_t i(0); i < lq.Size(); ++i)
        if (!lq[i].IsMassive())
          nlf++;
      m_nlf = nlf / 2.;
      m_b0 = (NC * 11. - 2. * m_nlf) / 6.;
    }

    ~HGQQ_QCD_Virtual() {}

    void Calc(const ATOOLS::Vec4D_Vector &momenta);

    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector &mom)
    {
      return 4. * M_PI;
    }
  };

  class HGGG_QCD_Virtual : public Virtual_ME2_Base{ // h->ggg virtual
    bool m_flip, m_con;
    double m_b0, m_nlf;

  public:
    HGGG_QCD_Virtual(const Process_Info &pi, const Flavour_Vector &flavs,
                     bool flip, bool con) : Virtual_ME2_Base(pi, flavs), m_flip(flip), m_con(con), m_b0(0.), m_nlf(0.)
    {
      msg_Out() << "HGGG_QCD_Virtual(): Using h->ggg virtual matrix element." << std::endl;
      Flavour lq(kf_quark);
      double nlf(0.);
      for (size_t i(0); i < lq.Size(); ++i)
        if (!lq[i].IsMassive())
          nlf++;
      m_nlf = nlf / 2.;
      m_b0 = (NC * 11. - 2. * m_nlf) / 6.;
    }

    ~HGGG_QCD_Virtual() {}

    void Calc(const ATOOLS::Vec4D_Vector &momenta);

    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector &mom)
    {
      return 4. * M_PI;
    }
  };
}

using namespace EXTRAXS;

void GGHG_QCD_Virtual::Calc(const Vec4D_Vector& momenta) {
  double s((momenta[0]+momenta[1]).Abs2());
  double t((momenta[1]-momenta.back()).Abs2());
  double u((momenta[0]-momenta.back()).Abs2());
  // check
  msg_Out() << "GGHG_QCD_Virtual::Calc(): Calculating gg->hg virtual matrix element for: " << "\n";
  msg_Out() << "   s = " << s << "\n";
  msg_Out() << "   t = " << t << "\n";
  msg_Out() << "   u = " << u << "\n";
  msg_Out() << "   mom[0]     = " << momenta[0] << "\n";
  msg_Out() << "   mom[1]     = " << momenta[1] << "\n";
  msg_Out() << "   mom.back() = " << momenta.back() << "\n";
  msg_Out() << "   mom[0].mass     = " << momenta[0].Mass() << "\n";
  msg_Out() << "   mom[1].mass     = " << momenta[1].Mass() << "\n";
  msg_Out() << "   mom.back().mass = " << momenta.back().Mass() << "\n";

  double mh2(s+t+u);
  double logg((sqr(mh2*mh2)+sqr(s*s)+sqr(t*t)+sqr(u*u))/(s*t*u));
  double lnm(log(m_mur2/mh2));
  double lns(log(s/mh2)), lnt(log(-t/mh2)), lnu(log(-u/mh2));
  double ln2t(log((mh2-t)/mh2)), ln2u(log((mh2-u)/mh2));
  double Li2s(DiLog((s-mh2)/s)), Li2t(DiLog(t/mh2)), Li2u(DiLog(u/mh2));
  // 1/epsIR
  m_res.IR()=NC*(lns+lnt+lnu-3.*lnm)
             -3.*m_b0;
  // 1/epsIR2
  m_res.IR2()=-3.*NC;
  // finite -> slight difference to MCFM
  m_res.Finite()=NC*(2.*(Li2t+Li2u+Li2s)
                     +lnm*(lns+lnt+lnu)-lns*lnt-lns*lnu-lnt*lnu
                     +0.5*(sqr(lns)-sqr(lnt)-sqr(lnu))-1.5*sqr(lnm)
                     +2.*(lnu*ln2u+lnt*ln2t)+4./3.*sqr(M_PI))
                 +(NC-m_nlf)/3.*mh2*(1.+mh2/s+mh2/t+mh2/u)/logg
                 +(m_con?0.:11.);
}

void QQHG_QCD_Virtual::Calc(const Vec4D_Vector& momenta) {
  double s((momenta[0]+momenta[1]).Abs2());
  double t((momenta[1]-momenta.back()).Abs2());
  double u((momenta[0]-momenta.back()).Abs2());
  double mh2(s+t+u);
  double loqq(NC*CF/s*(sqr(t)+sqr(u)));
  double lnm(log(m_mur2/mh2));
  double lns(log(s/mh2)), lnt(log(-t/mh2)), lnu(log(-u/mh2));
  double ln2t(log((mh2-t)/mh2)), ln2u(log((mh2-u)/mh2));
  double Li2s(DiLog((s-mh2)/s)), Li2t(DiLog(t/mh2)), Li2u(DiLog(u/mh2));
  // 1/epsIR
  m_res.IR()=-2./3.*m_nlf
             +NC*(13./6.-2.*lnm+lnt+lnu)
             +1./NC*(1.5-lns+lnm)
             -3.*m_b0;
  // 1/epsIR2
  m_res.IR2()=-2.*NC+1./NC;
  // finite
  m_res.Finite()=m_nlf*(-10./9.+2./3.*lns-2./3.*lnm)
                 +NC*(40./9.+Li2t+Li2u+2.*Li2s-13./6.*(lns-lnm)
                      +(lnm-lns)*(lnt+lnu)
                      +sqr(lns)-sqr(lnm)-0.5*sqr(lnt)-0.5*sqr(lnu)
                      +lnt*ln2t+lnu*ln2u)
                 +1./NC*(4.-Li2t-Li2u-1.5*(lns-lnm)+0.5*sqr(lns-lnm)
                         +lnt*lnu-lnt*ln2t-lnu*ln2u)
                 -4./3.*sqr(M_PI)/NC
                 -0.25*(NC*sqr(NC)-1./NC)*(t+u)/loqq
                 +(m_con?0.:11.);
}

void GQHQ_QCD_Virtual::Calc(const Vec4D_Vector& momenta) {
  double s((momenta[0]+momenta[1]).Abs2());
  double t((momenta[1]-momenta.back()).Abs2());
  double u((momenta[0]-momenta.back()).Abs2());
  if (m_flip) std::swap(t,u);
  double mh2(s+t+u);
  double logq(-NC*CF/t*(sqr(s)+sqr(u)));
  double lnm(log(m_mur2/mh2));
  double lns(log(s/mh2)), lnt(log(-t/mh2)), lnu(log(-u/mh2));
  double ln2t(log((mh2-t)/mh2)), ln2u(log((mh2-u)/mh2));
  double Li2s(DiLog((s-mh2)/s)), Li2t(DiLog(t/mh2)), Li2u(DiLog(u/mh2));
  // 1/epsIR
  m_res.IR()=-2./3.*m_nlf
             +NC*(13./6.+lns-2.*lnm+lnu)
             +1./NC*(3./2.+lnm-lnt)
             -3.*m_b0;
  // 1/epsIR2
  m_res.IR2()=-2.*NC+1./NC;
  // finite
  m_res.Finite()=m_nlf*(-10./9.-2./3.*lnm+2./3.*lnt)
                 +NC*(40./9.+Li2u+2.*Li2t+Li2s
                            +lns*lnm-lns*lnt-13./6.*(lnt-lnm)
                            +lnm*lnu-sqr(lnm)-lnt*lnu-0.5*sqr(lnu)
                            +2.*lnt*ln2t+lnu*ln2u)
                 +1./NC*(4.-Li2u-Li2s+lns*lnu+0.5*sqr(lnt)-0.5*sqr(lns)
                             -lnm*lnt+0.5*sqr(lnm)-lnu*ln2u-1.5*(lnt-lnm))
                 +4./3.*sqr(M_PI)*NC
                 +0.25*(NC*sqr(NC)-1./NC)*(u+s)/logq
                 +(m_con?0.:11.);
}

void HGQQ_QCD_Virtual::Calc(const Vec4D_Vector &momenta) // h->gqq virtual calc
{ 
  double s((momenta[0] + momenta[1]).Abs2());
  double t((momenta[1] - momenta.back()).Abs2());
  double u((momenta[0] - momenta.back()).Abs2());
  // check
  msg_Out() << "HGQQ_QCD_Virtual::Calc(): Calculating h->gqq virtual matrix element for: " << "\n";
  msg_Out() << "   s = " << s << "\n";
  msg_Out() << "   t = " << t << "\n";
  msg_Out() << "   u = " << u << "\n";
  msg_Out() << "   mom[0]     = " << momenta[0] << "\n";
  msg_Out() << "   mom[1]     = " << momenta[1] << "\n";
  msg_Out() << "   mom.back() = " << momenta.back() << "\n";
  msg_Out() << "   mom[0].mass     = " << momenta[0].Mass() << "\n";
  msg_Out() << "   mom[1].mass     = " << momenta[1].Mass() << "\n";
  msg_Out() << "   mom.back().mass = " << momenta.back().Mass() << "\n";
  double mh2(s + t + u);
  double lnm(log(m_mur2 / mh2));
  double lns(log(s / mh2)), lnt(log(-t / mh2)), lnu(log(-u / mh2));
  double ln2t(log((mh2 - t) / mh2)), ln2u(log((mh2 - u) / mh2));
  double Li2s(DiLog((s - mh2) / s)), Li2t(DiLog(t / mh2)), Li2u(DiLog(u / mh2));

  // qq->hg kinematics
  double loqq(NC * CF / s * (sqr(t) + sqr(u)));
  // 1/epsIR
  m_res.IR() = -2. / 3. * m_nlf + NC * (13. / 6. - 2. * lnm + lnt + lnu) + 1. / NC * (1.5 - lns + lnm) - 3. * m_b0;
  // finite
  m_res.Finite() = m_nlf * (-10. / 9. - 2. / 3. * lnm + 2. / 3. * lns) + NC * (40. / 9. + Li2u + 2. * Li2s + Li2t - 13. / 6. * (lns - lnm) + (lnm - lns) * (lnt + lnu) - sqr(lnm) - 0.5 * sqr(lnt) - 0.5 * sqr(lnu) + sqr(lns) + lnt * ln2t + lnu * ln2u) 
                      + 1. / NC * (4. - Li2t - Li2u - 1.5 * (lns - lnm) + 0.5 * sqr(lns - lnm) + lnt * lnu - lnt * ln2t - lnu * ln2u) - 4. / 3. * sqr(M_PI) / NC - 0.25 * (NC * sqr(NC) - 1. / NC) * (t + u) / loqq + (m_con ? 0. : 11.);

  // gq->hq kinematics (swaps lns->lnt/lnu and lnt->lns)
  if (m_flip)
    std::swap(t, u);
  double logq(-NC * CF / t * (sqr(s) + sqr(u)));
  // 1/epsIR
  m_res.IR() = -2. / 3. * m_nlf + NC * (13. / 6. - 2. * lnm + lns + lnu) + 1. / NC * (1.5 - lnt + lnm) - 3. * m_b0;
  
  // both kinematics
  // 1/epsIR2
  m_res.IR2() = -2. * NC + 1. / NC;
  // finite
  m_res.Finite() = m_nlf * (-10. / 9. - 2. / 3. * lnm + 2. / 3. * lnt) + NC * (40. / 9. + Li2u + 2. * Li2t + Li2s - 13. / 6. * (lnt - lnm) + (lnm - lnt) * (lns + lnu) - sqr(lnm) - 0.5 * sqr(lnu) + 2. * lnt * ln2t + lnu * ln2u) + 
                      +1. / NC * (4. - Li2u - Li2s + lns * lnu + 0.5 * sqr(lnt) - 0.5 * sqr(lns) - lnm * lnt + 0.5 * sqr(lnm) - lnu * ln2u - 1.5 * (lnt - lnm)) + 4. / 3. * sqr(M_PI) * NC + 0.25 * (NC * sqr(NC) - 1. / NC) * (u + s) / logq + (m_con ? 0. : 11.);
}

void HGGG_QCD_Virtual::Calc(const Vec4D_Vector &momenta)
{ // h->ggg virtual calc
  double s((momenta[0] + momenta[1]).Abs2());
  double t((momenta[1] - momenta.back()).Abs2());
  double u((momenta[0] - momenta.back()).Abs2());
  // check
  msg_Out() << "HGGG_QCD_Virtual::Calc(): Calculating h->ggg virtual matrix element for: " << "\n";
  msg_Out() << "   s = " << s << "\n";
  msg_Out() << "   t = " << t << "\n";
  msg_Out() << "   u = " << u << "\n";
  msg_Out() << "   mom[0]     = " << momenta[0] << "\n";
  msg_Out() << "   mom[1]     = " << momenta[1] << "\n";
  msg_Out() << "   mom.back() = " << momenta.back() << "\n";
  msg_Out() << "   mom[0].mass     = " << momenta[0].Mass() << "\n";
  msg_Out() << "   mom[1].mass     = " << momenta[1].Mass() << "\n";
  msg_Out() << "   mom.back().mass = " << momenta.back().Mass() << "\n";
  double mh2(s + t + u);
  double logg((sqr(mh2 * mh2) + sqr(s * s) + sqr(t * t) + sqr(u * u)) / (s * t * u));
  double lnm(log(m_mur2 / mh2));
  double lns(log(s / mh2)), lnt(log(-t / mh2)), lnu(log(-u / mh2));
  double ln2t(log((mh2 - t) / mh2)), ln2u(log((mh2 - u) / mh2));
  double Li2s(DiLog((s - mh2) / s)), Li2t(DiLog(t / mh2)), Li2u(DiLog(u / mh2));
  // 1/epsIR
  m_res.IR() = NC * (lns + lnt + lnu - 3. * lnm) - 3. * m_b0;
  // 1/epsIR2
  m_res.IR2() = -3. * NC;
  // finite -> slight difference to MCFM
  m_res.Finite() = NC * (2. * (Li2t + Li2u + Li2s) + lnm * (lns + lnt + lnu) - lns * lnt - lns * lnu - lnt * lnu + 0.5 * (sqr(lns) - sqr(lnt) - sqr(lnu)) - 1.5 * sqr(lnm) + 2. * (lnu * ln2u + lnt * ln2t) + 4. / 3. * sqr(M_PI)) + (NC - m_nlf) / 3. * mh2 * (1. + mh2 / s + mh2 / t + mh2 / u) / logg + (m_con ? 0. : 11.);
}

DECLARE_VIRTUALME2_GETTER(EXTRAXS::GGHG_QCD_Virtual,"Higgs_Jet_QCD_Virtual")
Virtual_ME2_Base *ATOOLS::Getter
<PHASIC::Virtual_ME2_Base,PHASIC::Process_Info,EXTRAXS::GGHG_QCD_Virtual>::
operator()(const Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  if (pi.m_loopgenerator!="Internal") return NULL;
  if (pi.m_fi.m_nlotype&nlo_type::loop) {
    if (pi.m_fi.m_nlocpl[1]!=0.) return NULL;
    if (pi.m_fi.m_asscontribs!=asscontrib::none) {
      msg_Error()<<"Higgs_Jet_QCD_Virtual(): Error: cannot provide requested "
                 <<"associated contributions "<<pi.m_fi.m_asscontribs<<std::endl;
      return NULL;
    }
    Settings& s = Settings::GetMainSettings();
    int con = s["HNNLO_KF_MODE"].Get<int>();
    Flavour_Vector fl=pi.ExtractFlavours();
    if (fl[0].IsGluon() && fl[1].IsGluon() && 
        pi.m_fi.m_ps.size()==2 && 
        pi.m_fi.m_ps[0].m_fl.Kfcode()==kf_h0 &&
        pi.m_fi.m_ps[1].m_fl.IsGluon()) {
      for (size_t i=2; i<fl.size()-1; ++i) {
        if (fl[i].Strong()) return NULL;
      }
      return new GGHG_QCD_Virtual(pi, fl, false, con);
    }
    if (fl[0].IsGluon() && fl[1].IsQuark() && 
        pi.m_fi.m_ps.size()==2 && 
        pi.m_fi.m_ps[0].m_fl.Kfcode()==kf_h0 &&
        pi.m_fi.m_ps[1].m_fl==fl[1]) {
      for (size_t i=2; i<fl.size()-1; ++i) {
        if (fl[i].Strong()) return NULL;
      }
      return new GQHQ_QCD_Virtual(pi, fl, false, con);
    }
    if (fl[1].IsGluon() && fl[0].IsQuark() && 
        pi.m_fi.m_ps.size()==2 && 
        pi.m_fi.m_ps[0].m_fl.Kfcode()==kf_h0 &&
        pi.m_fi.m_ps[1].m_fl==fl[0]) {
      for (size_t i=2; i<fl.size()-1; ++i) {
        if (fl[i].Strong()) return NULL;
      }
      return new GQHQ_QCD_Virtual(pi, fl, true, con);
    }
    if (fl[0].IsQuark() && fl[1]==fl[0].Bar() && 
        pi.m_fi.m_ps.size()==2 && 
        pi.m_fi.m_ps[0].m_fl.Kfcode()==kf_h0 &&
        pi.m_fi.m_ps[1].m_fl.IsGluon()) {
      for (size_t i=2; i<fl.size()-1; ++i) {
        if (fl[i].Strong()) return NULL;
      }
      return new QQHG_QCD_Virtual(pi, fl, false, con);
    }
    if (fl[0].Kfcode() == kf_e)
    {
      msg_Out() << "in HGQQ unchecked getter. you have been warned \n";
      
      return new HGQQ_QCD_Virtual(pi, fl, false, con);
    }
    if (fl[0].Kfcode() == kf_e)
    {
      msg_Out() << "in HGGG unchecked getter. you have been warned \n";
      
      return new HGGG_QCD_Virtual(pi, fl, false, con);
    }
  }
  return NULL;
}

