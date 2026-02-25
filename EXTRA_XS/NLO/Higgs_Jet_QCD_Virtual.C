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
      msg_Out() << "GGHG_QCD_Virtual(): Using gg->hg virtual matrix element." << std::endl;

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
      msg_Out() << "QQHG_QCD_Virtual(): Using qq->hg virtual matrix element." << std::endl;

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
      msg_Out() << "GQHQ_QCD_Virtual(): Using gq->hq virtual matrix element." << std::endl;

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
      msg_Tracking() << "HGQQ_QCD_Virtual(): Using h->gqq virtual matrix element." << "\n";
      msg_Debugging() << "flavs[0] = " << flavs[0] << ", flavs[1] = " << flavs[1] << ", flavs[2] = " << flavs[2] << ", flavs[3] = " << flavs[3] << ", flavs[4] = " << flavs[4] << "\n";

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
      msg_Tracking() << "HGGG_QCD_Virtual(): Using h->ggg virtual matrix element." << "\n";
      msg_Debugging() << "flavs[0] = " << flavs[0] << ", flavs[1] = " << flavs[1] << ", flavs[2] = " << flavs[2] << ", flavs[3] = " << flavs[3] << ", flavs[4] = " << flavs[4] << "\n";
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

void HGQQ_QCD_Virtual::Calc(const Vec4D_Vector &momenta) // ee->h->gqq virtual calc
{ 
  // msg_Out() << "HGQQ_QCD_Virtual::Calc(): Calculating h->gqq virtual matrix element for: " << "\n";
  // Vec4D ph = (momenta[0] + momenta[1]); // higgs momentum
  // Vec4D pg = momenta[2];                // gluon momentum 
  // Vec4D pq = momenta[3];                // quark momentum
  // Vec4D pqbar = momenta[4];             // antiquark momentum
  // Vec4D pgtilde = -pg;                  // so that 2-2 kinematics are h(g~)->qq instead of h->gqq

  // double s((ph + pgtilde).Abs2()); 
  // double t((ph - pq).Abs2());
  // double u((ph - pqbar).Abs2());
  // double mh2(s+t+u);

  // if (s != ((pq + pqbar).Abs2()))
  // {
  //   msg_Error() << "HGQQ_QCD_Virtual::Calc(): Kinematics Error: s must be equal to (pq+pqbar)^2, but is s=" << s << " and (pq+pqbar)^2=" << (pq + pqbar).Abs2() << "\n";
  //   msg_Error() << "   mh2=s+t+u =" << mh2 << "\n";
  //   throw ATOOLS::Exception("Invalid kinematics");
  // }
  // if (t != ((pgtilde - pqbar).Abs2()))
  // {
  //   msg_Error() << "HGQQ_QCD_Virtual::Calc(): Kinematics Error: t must be equal to (pgtilde-pqbar)^2, but is t=" << t << " and (pgtilde-pqbar)^2=" << (pgtilde - pqbar).Abs2() << "\n";
  //   throw ATOOLS::Exception("Invalid kinematics");
  // }
  // if (u != ((pgtilde - pq).Abs2()))
  // {
  //   msg_Error() << "HGQQ_QCD_Virtual::Calc(): Kinematics Error: u must be equal to (pgtilde-pq)^2, but is u=" << u << " and (pgtilde-pq)^2=" << (pgtilde - pq).Abs2() << "\n";
  //   throw ATOOLS::Exception("Invalid kinematics");
  // }

  msg_Out() << "HGQQ_QCD_Virtual::Calc(): Calculating h->gqq virtual matrix element for: " << "\n";
  Vec4D p0 = -momenta[3];                // "incoming" quark
  Vec4D p1 = -momenta[4];                // "incoming" antiquark
  Vec4D p2 = -(momenta[0]+momenta[1]);   // "outgoing" higgs
  Vec4D p3 = momenta[2];                 // outgoing gluon
  
  double s((p0 + p1).Abs2());
  double t((p1 - p3).Abs2());
  double u((p0 - p3).Abs2());
  double mh2(s + t + u);
  msg_Tracking() << "HGQQ_QCD_Virtual::Calc(): Calculated kinematic invariants: s=" << s << ", t=" << t << ", u=" << u << ", mh2=" << mh2 << "\n";

  if (s != ((p2+p3).Abs2()))
  {
    msg_Error() << "HGQQ_QCD_Virtual::Calc(): Kinematics Error: s must be equal to (p2+p3)^2, but is s=" << s << " and (p2+p3)^2=" << (p2+p3).Abs2() << "\n";
    throw ATOOLS::Exception("Invalid kinematics");
  }
  if (t != ((p0-p2).Abs2()))
  {
    msg_Error() << "HGQQ_QCD_Virtual::Calc(): Kinematics Error: t must be equal to (p0-p2)^2, but is t=" << t << " and (p0-p2)^2=" << (p0-p2).Abs2() << "\n";
    throw ATOOLS::Exception("Invalid kinematics");
  }
  if (u != ((p1-p2).Abs2()))
  {
    msg_Error() << "HGQQ_QCD_Virtual::Calc(): Kinematics Error: u must be equal to (p1-p2)^2, but is u=" << u << " and (p1-p2)^2=" << (p1-p2).Abs2() << "\n";
    throw ATOOLS::Exception("Invalid kinematics");
  }

  // check
  // msg_Out() << "   mom[0]     = " << momenta[0] << "\n";
  // msg_Out() << "   mom[1]     = " << momenta[1] << "\n";
  // msg_Out() << "   mom[2]     = " << momenta[2] << "\n";
  // msg_Out() << "   mom[3]     = " << momenta[3] << "\n";
  // msg_Out() << "   mom[4]     = " << momenta[4] << "\n";
  // msg_Out() << "   mom[0].mass     = " << momenta[0].Mass() << "\n";
  // msg_Out() << "   mom[1].mass     = " << momenta[1].Mass() << "\n";
  // msg_Out() << "   momH.mass       = " << (momenta[0]+momenta[1]).Mass() << "\n";
  // msg_Out() << "   mom[2].mass     = " << momenta[2].Mass() << "\n";
  // msg_Out() << "   momGtilde.mass  = " << pgtilde.Mass() << "\n";
  // msg_Out() << "   mom[3].mass     = " << momenta[3].Mass() << "\n";
  // msg_Out() << "   mom[4].mass     = " << momenta[4].Mass() << "\n";

  // try again? 
  // Vec4D pk = momenta[2];
  // Vec4D pij = momenta[3] + momenta[4]; 

  // double s((momenta[0]+momenta[1]).Abs2());
  // double t((momenta[0]-pk).Abs2());
  // double u((momenta[1]-pk).Abs2());
  // double mh2( s );
  // msg_Tracking() << "HGQQ_QCD_Virtual::Calc(): Calculated kinematic invariants: s=" << s << ", t=" << t << ", u=" << u << ", mh2=" << mh2 << "\n";
  // if (s != ((pk+pij).Abs2()))
  // {
  //   msg_Error() << "HGQQ_QCD_Virtual::Calc(): Kinematics Error: s must be equal to (pq+pqbar)^2, but is s=" << s << " and (pq+pqbar)^2=" << (pk + pij).Abs2() << "\n";
  //   msg_Error() << "   mh2=s+t+u =" << mh2 << "\n";
  //   throw ATOOLS::Exception("Invalid kinematics");
  // }
  // if (t != ((momenta[1]-pij).Abs2()))
  // {
  //   msg_Error() << "HGQQ_QCD_Virtual::Calc(): Kinematics Error: t must be equal to (momenta[1]-pij)^2, but is t=" << t << " and (momenta[1]-pij)^2=" << (momenta[1]-pij).Abs2() << "\n";
  //   throw ATOOLS::Exception("Invalid kinematics");
  // }
  // if (u != ((momenta[0]-pij).Abs2()))
  // {
  //   msg_Error() << "HGQQ_QCD_Virtual::Calc(): Kinematics Error: u must be equal to (momenta[0]-pij)^2, but is u=" << u << " and (momenta[0]-pij)^2=" << (momenta[0]-pij).Abs2() << "\n";
  //   throw ATOOLS::Exception("Invalid kinematics");
  // }

  if (m_mur2 <= 0.|| mh2 <= 0.) {
    msg_Error() << "HGQQ_QCD_Virtual::Calc(): Log Error: renormalization scale squared must be positive, but is " << m_mur2 << "\n";
    msg_Error() << "   mh2 = " << mh2 << "\n";
    throw ATOOLS::Exception("Invalid renormalization scale");
  }
  double lnm(log(m_mur2 / mh2));
  if (s <= 0. || t >= 0. || u >= 0.)
  {
    msg_Error() << "HGQQ_QCD_Virtual::Calc(): Log Error: ! kinematic invariants must satisfy s>0, t<0, u<0, but are s=" << s << ", t=" << t << ", u=" << u << "\n";
    // throw ATOOLS::Exception("Invalid kinematics"); //breaks here
  }
  double lns(log(s / mh2)), lnt(log(-t / mh2)), lnu(log(-u / mh2));
  if (mh2 - t <= 0. || mh2 - u <= 0.) {
    msg_Error() << "HGQQ_QCD_Virtual::Calc(): Log Error: kinematic invariants must satisfy mh2-t>0, mh2-u>0, but are mh2-t=" << mh2 - t << ", mh2-u=" << mh2 - u << "\n";
    throw ATOOLS::Exception("Invalid kinematics");
  }
  double ln2t(log((mh2 - t) / mh2)), ln2u(log((mh2 - u) / mh2));
  double Li2s(DiLog((s - mh2) / s)), Li2t(DiLog(t / mh2)), Li2u(DiLog(u / mh2));

  // qq<->hg kinematics
  double loqq(NC * CF / s * (sqr(t) + sqr(u)));
  msg_Tracking() << "HGQQ_QCD_Virtual::Calc(): Calculated leading order term: loqq = " << loqq << "\n";
  // 1/epsIR
  m_res.IR() = -2. / 3. * m_nlf + NC * (13. / 6. - 2. * lnm + lnt + lnu) + 1. / NC * (1.5 - lns + lnm) - 3. * m_b0;
  // 1/epsIR2
  m_res.IR2() = -2. * NC + 1. / NC;
  // finite
  m_res.Finite() = m_nlf * (-10. / 9. - 2. / 3. * lnm + 2. / 3. * lnt) + NC * (40. / 9. + Li2u + 2. * Li2t + Li2s - 13. / 6. * (lnt - lnm) + (lnm - lnt) * (lns + lnu) - sqr(lnm) - 0.5 * sqr(lnu) + 2. * lnt * ln2t + lnu * ln2u) + 
                      +1. / NC * (4. - Li2u - Li2s + lns * lnu + 0.5 * sqr(lnt) - 0.5 * sqr(lns) - lnm * lnt + 0.5 * sqr(lnm) - lnu * ln2u - 1.5 * (lnt - lnm)) + 4. / 3. * sqr(M_PI) * NC + 0.25 * (NC * sqr(NC) - 1. / NC) * (u + s) / loqq + (m_con ? 0. : 11.);
  msg_Tracking() << "HGQQ_QCD_Virtual::Calc(): Calculated virtual correction: IR2=" << m_res.IR2() << ", IR=" << m_res.IR() << ", Finite=" << m_res.Finite() << "\n";
  throw ATOOLS::Exception("Break here for debugging");
}

void HGGG_QCD_Virtual::Calc(const Vec4D_Vector &momenta)
{ // ee->h->ggg virtual calc
  Vec4D ph = momenta[0] + momenta[1]; // higgs momentum
  Vec4D pgtilde = -momenta[2];        // flipped gluon momentum so that 2-2 kinematics are h(g~)->gg instead of h->ggg
  double s((ph + pgtilde).Abs2());
  double t((ph - momenta[3]).Abs2());
  double u((ph - momenta[4]).Abs2());
  // check
  // msg_Out() << "HGGG_QCD_Virtual::Calc(): Calculating h->ggg virtual matrix element for: " << "\n";
  // msg_Out() << "   mom[0]     = " << momenta[0] << "\n";
  // msg_Out() << "   mom[1]     = " << momenta[1] << "\n";
  // msg_Out() << "   mom[2]     = " << momenta[2] << "\n";
  // msg_Out() << "   mom[3]     = " << momenta[3] << "\n";
  // msg_Out() << "   mom[4]     = " << momenta[4] << "\n";
  // msg_Out() << "   mom[0].mass     = " << momenta[0].Mass() << "\n";
  // msg_Out() << "   mom[1].mass     = " << momenta[1].Mass() << "\n";
  // msg_Out() << "   momH.mass       = " << (momenta[0] + momenta[1]).Mass() << "\n";
  // msg_Out() << "   mom[2].mass     = " << momenta[2].Mass() << "\n";
  // msg_Out() << "   momGtilde.mass  = " << pgtilde.Mass() << "\n";
  // msg_Out() << "   mom[3].mass     = " << momenta[3].Mass() << "\n";
  // msg_Out() << "   mom[4].mass     = " << momenta[4].Mass() << "\n";

  double mh2(s + t + u);
  if (m_mur2 <= 0. || mh2 <= 0.)
  {
    msg_Error() << "HGGG_QCD_Virtual::Calc(): Log Error: renormalization scale squared must be positive, but is " << m_mur2 << "\n";
    msg_Error() << "   mh2 = " << mh2 << "\n";
    throw ATOOLS::Exception("Invalid renormalization scale");
  }
  double lnm(log(m_mur2 / mh2));
  if (s <= 0. || t >= 0. || u >= 0.)
  {
    msg_Error() << "HGGG_QCD_Virtual::Calc(): Log Error: kinematic invariants must satisfy s>0, t<0, u<0, but are s=" << s << ", t=" << t << ", u=" << u << "\n";
    throw ATOOLS::Exception("Invalid kinematics");
  }
  double lns(log(s / mh2)), lnt(log(-t / mh2)), lnu(log(-u / mh2));
  double ln2t(log((mh2 - t) / mh2)), ln2u(log((mh2 - u) / mh2));
  double Li2s(DiLog((s - mh2) / s)), Li2t(DiLog(t / mh2)), Li2u(DiLog(u / mh2));
  // gg<->hg kinematics
  double logg((sqr(mh2 * mh2) + sqr(s * s) + sqr(t * t) + sqr(u * u)) / (s * t * u));
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
    if (fl[0].Kfcode() == kf_e && fl[1].Kfcode() == kf_e // ee->h->gqq
                                                         // && pi.m_fi.m_fl.Kfcode() == kf_h0   // explicit Higgs intermediate state
        && pi.m_fi.m_ps.size() == 3 
        && pi.m_fi.m_ps[1].m_fl.IsQuark() && pi.m_fi.m_ps[2].m_fl == pi.m_fi.m_ps[1].m_fl.Bar() // check quark/gluon ordering
        && pi.m_fi.m_ps[0].m_fl.IsGluon())
    {
      msg_Tracking() << "\n In HGQQ unchecked getter. you have been warned \n";
      msg_Debugging() << " fl[0] = " << fl[0] << ", fl[1] = " << fl[1] << ", fl[2] = " << fl[2] << ", fl[3] = " << fl[3] << ", fl[4] = " << fl[4] << "\n";
      if (pi.m_maxcpl[0] == 4 && pi.m_maxcpl[1] == 0 && pi.m_maxcpl[2] == 2 &&
          pi.m_mincpl[0] == 4 && pi.m_mincpl[1] == 0 && pi.m_mincpl[2] == 2)
      {
        for (size_t i = 0; i < fl.size() - 3; ++i)
        {
          if (fl[i].Strong())
          {
            msg_Error() << "HGQQ_QCD_Virtual getter: Warning: strong particle in initial state: " << fl[i] << "\n";
            return NULL;
          }
        }
        for (size_t i = 2; i < fl.size()-2; ++i)
        { // only strong stuff in final state
          if (!fl[i].Strong()) {
            msg_Error() << "HGQQ_QCD_Virtual getter: Warning: strong particles in wrong order in final state: " << fl[i] << "\n";
            return NULL;}
        }
        return new HGQQ_QCD_Virtual(pi, fl, false, con);
      }
    }
    if (fl[0].Kfcode() == kf_e && fl[1].Kfcode() == kf_e // ee->h->ggg
        // && pi.m_fi.m_fl.Kfcode() == kf_h0   // explicit Higgs intermediate state
        && pi.m_fi.m_ps.size() == 3 &&
        pi.m_fi.m_ps[0].m_fl.IsGluon() && pi.m_fi.m_ps[1].m_fl.IsGluon() && pi.m_fi.m_ps[2].m_fl.IsGluon())
    {
      msg_Tracking() << "\n In HGGG unchecked getter. you have been warned \n";
      msg_Debugging() << " fl[0] = " << fl[0] << ", fl[1] = " << fl[1] << ", fl[2] = " << fl[2] << ", fl[3] = " << fl[3] << ", fl[4] = " << fl[4] << "\n";
      if (pi.m_maxcpl[0] == 4 && pi.m_maxcpl[1] == 0 && pi.m_maxcpl[2] == 2 &&
          pi.m_mincpl[0] == 4 && pi.m_mincpl[1] == 0 && pi.m_mincpl[2] == 2)
      {
        for (size_t i = 0; i < fl.size()-3; ++i) 
        { 
          if (fl[i].Strong())
          {
            msg_Error() << "HGGG_QCD_Virtual getter: Warning: strong particle in initial state: " << fl[i] << "\n";
            return NULL;
          }
        }
        for (size_t i = 2; i < fl.size(); ++i) 
        { // only strong stuff in final state
          if (!fl[i].Strong())
          {
            msg_Error() << "HGGG_QCD_Virtual getter: Warning: non-strong particle in final state: " << fl[i] << "\n";
            return NULL;
          }
        }
        return new HGGG_QCD_Virtual(pi, fl, false, con);
      }
    }
  }
  return NULL;
}

