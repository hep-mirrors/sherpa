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
  const Complex I = Complex(0., 1.);

  class GGHG_QCD_Virtual : public Virtual_ME2_Base
  { // gg->hg virtual
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
      // =(3*11-2*m_nlf)/6 
      // ?= (11.*CA -4.*TR*NF) = (11*3-4*0.5*NF)/6
      // yes if NF=m_nlf
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
      msg_Out() << "HGQQ_QCD_Virtual(): Using h->gqq virtual matrix element." << "\n";
      msg_Out() << "flavs[0] = " << flavs[0] << ", flavs[1] = " << flavs[1] << ", flavs[2] = " << flavs[2] << ", flavs[3] = " << flavs[3] << ", flavs[4] = " << flavs[4] << "\n";

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
      msg_Out() << "HGGG_QCD_Virtual(): Using h->ggg virtual matrix element." << "\n";
      msg_Out() << "flavs[0] = " << flavs[0] << ", flavs[1] = " << flavs[1] << ", flavs[2] = " << flavs[2] << ", flavs[3] = " << flavs[3] << ", flavs[4] = " << flavs[4] << "\n";
      
      Flavour lq(kf_quark);
      double nlf(0.);
      for (size_t i(0); i < lq.Size(); ++i)
        if (!lq[i].IsMassive())
          nlf++;
      m_nlf = nlf / 2.;
      m_b0 = (11.*CA -4.*TR*m_nlf)/6.;
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
                     +2.*(lnu*ln2u+lnt*ln2t)
                     
                     +4./3.*sqr(M_PI))

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
  // finite: qg->hq 
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

void HGQQ_QCD_Virtual::Calc(const Vec4D_Vector &momenta)
{ // ee->h->gqq virtual calc
  msg_Tracking() << METHOD << " ==============================================" << "\n";

  Vec4D ph = momenta[0] + momenta[1]; // higgs momentum
  // Dalitz variables
  // quark     : mom 3 or 4 : index 1
  // antiquark : mom 3 or 4 : index 2
  // gluon     : mom 2      : index 3
  double s12((momenta[3] + momenta[4]).Abs2());
  double s13((momenta[2] + momenta[3]).Abs2());
  double s23((momenta[2] + momenta[4]).Abs2()); 

  double mh2(s12 + s13 + s23); // should be mH^2 if quarks massless and higgs onshell
  msg_Tracking() << METHOD << ": ph.Abs2 =" << ph.Abs2() << ", s123 = mh^2 =" << mh2 << "\n";
  //  Normalized Dalitz variables
  double x(s12 / mh2);
  double y(s13 / mh2);
  double z(s23 / mh2);

  msg_Tracking() << METHOD << ": Calculated kinematic invariants: x=" << x << ", y=" << y << ", z=" << z << "\n";

  double lnm(log(m_mur2 / mh2));
  double lnxyz(log(x * y * z));
  double lnx(log(x)), lny(log(y)), lnz(log(z));
  double ln1x(log(y + z)), ln1y(log(x + z)), ln1z(log(x + y));
  double Li2x(DiLog(y + z)), Li2y(DiLog(x + z)), Li2z(DiLog(x + y));
  double sum4(pow(x, 4) + pow(y, 4) + pow(z, 4) + 1);
  double IR2(1/NC -2*NC);

  // 1/epsIR
  m_res.IR() = - 3.*m_b0 - (2./3.)*m_nlf +NC*(-2*lnm +lny+lnz +13./6.) +(1./NC)*(lnm -lnx +3./2.) ;
  // 1/epsIR2
  m_res.IR2() = IR2;
  // finite
  m_res.Finite() = (1./2.)*IR2*lnm*lnm + lnm*( (m_b0/2.) -(lnx/NC) +NC*(lny+lnz) -3.*CF) 
                    + m_nlf*( -lnx/3. -(lny+lnz)/4. +5./9.)
                    + NC*Li2x + CF*(Li2y+Li2z)
                    + sqr(M_PI)*(-1./2.)*IR2
                    + sqr(lnx)/(2*NC) -(sqr(lny)+sqr(lnz))*(NC/2.)
                    -(3./(2.*NC))*lnx +(5.*NC/2.)*(lny+lnz)
                    +(1./2.)*( NC*lnx*(lny+lnz) -(1./NC)*lny*lnz)
                    +(NC+1./NC)*(x/(y*z))*(pow(y,3)+pow(z,3))/((pow(y,2)+pow(z,2)))
                    -NC*20./9. -(1./NC)*2./9. ;

  msg_Tracking() << METHOD << ": Calculated virtual correction: IR2=" << m_res.IR2() << ", IR=" << m_res.IR() << ", Finite=" << m_res.Finite() << "\n";
}

void HGGG_QCD_Virtual::Calc(const Vec4D_Vector &momenta)
{ // ee->h->ggg virtual calc
  msg_Tracking() << METHOD << " ==============================================" << "\n";

  Vec4D ph = momenta[0] + momenta[1]; // higgs momentum
  // Dalitz variables
  double s12((momenta[2]+momenta[3]).Abs2());
  double s13((momenta[2]+momenta[4]).Abs2());
  double s23((momenta[3]+momenta[4]).Abs2());
  double mh2(s12+s13+s23); // should be mH^2 if higgs onshell
  msg_Tracking() << METHOD << ": ph.Abs2 =" << ph.Abs2() << ", s123 = mh^2 =" << mh2 << "\n";
  //  Normalized Dalitz variables
  double x(s12/mh2);
  double y(s13/mh2);
  double z(s23/mh2);

  msg_Tracking() << METHOD << ": Calculated kinematic invariants: x=" << x << ", y=" << y << ", z=" << z << "\n";
  
  double lnm(log(m_mur2 / mh2));
  double lnxyz(log(x*y*z));
  double lnx(log(x)), lny(log(y)), lnz(log(z));
  double ln1x(log(y+z)), ln1y(log(x+z)), ln1z(log(x+y));
  double Li2x(DiLog(y+z)), Li2y(DiLog(x+z)), Li2z(DiLog(x+y));
  double sum4(pow(x,4)+pow(y,4)+pow(z,4)+1);
  double IR2(-3.*NC);

  // 1/epsIR
  m_res.IR()  = NC * (lnxyz - lnm) - 3. * m_b0;
  // 1/epsIR2
  m_res.IR2() = IR2; 
  // finite 
  m_res.Finite() = (1./2.)*(IR2)*lnm*lnm +lnm*(lnxyz*NC + m_b0/2.) 
                    +(1./(6.*sum4))*( (NC-m_nlf)*( x*y*z*z +z*(x+y)*(x*y+1) +x*y ))
                    +2*sqr(M_PI)*NC
                    +(1./2.)*m_b0*lnxyz
                    -NC*(Li2x+Li2y+Li2z)
                    -(NC/2.)*(sqr(lnx)+sqr(lny)+sqr(lnz) +lnx*lny +lnx*lnz +lny*lnz);

  msg_Tracking() << METHOD << ": Calculated virtual correction: IR2=" << m_res.IR2() << ", IR=" << m_res.IR() << ", Finite=" << m_res.Finite() << "\n";
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
        { // gluons then quarks 
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
        { // only strong final state
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

