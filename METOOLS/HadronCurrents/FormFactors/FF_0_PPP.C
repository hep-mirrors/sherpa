#include "METOOLS/HadronCurrents/FormFactors/FF_0_PPP.H"
#include "METOOLS/HadronCurrents/FormFactors/A1_Decays.H"
#include "METOOLS/HadronCurrents/Tools.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

namespace {
  long int SignedKfcode(const Flavour & flav) {
    const long int kfcode = static_cast<long int>(flav.Kfcode());
    return flav.IsAnti() ? -kfcode : kfcode;
  }

  bool Match(const Flavour_Vector & flavs,const vector<int> & indices,
             const long int kf0,const long int kf1,const long int kf2) {
    return SignedKfcode(flavs[indices[0]])==kf0 &&
           SignedKfcode(flavs[indices[1]])==kf1 &&
           SignedKfcode(flavs[indices[2]])==kf2;
  }

  FF_0_PPP_mode ModeFromFlavours(const Flavour_Vector & flavs,
                                 const vector<int> & indices) {
    if (indices.size()!=3) return FF_0_PPP_mode::unknown;

    const long int pi0  =  static_cast<long int>(kf_pi);
    const long int pip  =  static_cast<long int>(kf_pi_plus);
    const long int pim  = -static_cast<long int>(kf_pi_plus);
    const long int K0   =  static_cast<long int>(kf_K);
    const long int K0b  = -static_cast<long int>(kf_K);
    const long int Kp   =  static_cast<long int>(kf_K_plus);
    const long int Km   = -static_cast<long int>(kf_K_plus);
    const long int KS   =  static_cast<long int>(kf_K_S);
    const long int KL   =  static_cast<long int>(kf_K_L);

    // The pion currents use the identical particles in slots q1,q2 and the
    // odd-charge pion in q3; the KS-to-METOOLS basis conversion relies on it.
    if (flavs[indices[0]].Kfcode()==kf_pi &&
        flavs[indices[1]].Kfcode()==kf_pi &&
        flavs[indices[2]].Kfcode()==kf_pi_plus) {
      return FF_0_PPP_mode::piP_pi0_pi0;
    }
    if (flavs[indices[0]].Kfcode()==kf_pi_plus &&
        flavs[indices[1]].Kfcode()==kf_pi_plus &&
        flavs[indices[2]].Kfcode()==kf_pi_plus &&
        SignedKfcode(flavs[indices[0]])==SignedKfcode(flavs[indices[1]]) &&
        SignedKfcode(flavs[indices[0]])==-SignedKfcode(flavs[indices[2]])) {
      return FF_0_PPP_mode::piM_piP_piP;
    }

    // Finkemeier-Mirkes 1995 table ordering: q1 q2 q3 == a b c.
    if (Match(flavs,indices,Km,pim,Kp))  return FF_0_PPP_mode::K_pi_K;
    if (Match(flavs,indices,K0,pim,K0b)) return FF_0_PPP_mode::K0_pi_K0b;
    if (Match(flavs,indices,KS,pim,KS))  return FF_0_PPP_mode::KS_pi_KS;
    if (Match(flavs,indices,KS,pim,KL))  return FF_0_PPP_mode::KS_pi_KL;
    if (Match(flavs,indices,KL,pim,KL))  return FF_0_PPP_mode::KL_pi_KL;
    if (Match(flavs,indices,Km,pi0,K0))  return FF_0_PPP_mode::K_pi0_K0;
    if (Match(flavs,indices,pi0,pi0,Km)) return FF_0_PPP_mode::pi0_pi0_K;
    if (Match(flavs,indices,Km,pim,pip)) return FF_0_PPP_mode::K_pi_pi;
    if (Match(flavs,indices,pim,K0b,pi0)) return FF_0_PPP_mode::pi_K0b_pi0;

    return FF_0_PPP_mode::unknown;
  }
}


void FF_0_PPP_Base::FixMode() {
  m_mode = ModeFromFlavours(m_flavs,m_pi);
}

Complex FF_0_PPP_Base::operator()(const ATOOLS::Vec4D_Vector& moms) {
  Vec4D p1    = moms[m_pi[0]], p2 = moms[m_pi[1]], p3 = moms[m_pi[2]];
  Vec4D q     = p1+p2+p3;
  double s123 = q.Abs2();
  double s12  = (p1+p2).Abs2();
  double s13  = (p1+p3).Abs2();
  double s23  = (p2+p3).Abs2();

  switch (m_ffmodel) {
  case ff_model::none:
    return Complex(1.,0.);

  case ff_model::KS:
    return m_norm * FF_KS(s123,s12,s13,s23);

  case ff_model::RChiPT:
    return m_norm * FF_RChiPT(s123,s12,s13,s23);

  case ff_model::unknown:
  default:
    break;
  }

  return Complex(0.,0.);
}

Complex FF_0_PPP_Base::
FF_KS(const double & s123,const double & s12,const double & s13,
      const double & s23) {
  msg_Error()<<"Error in "<<METHOD<<": KS form factor not implemented for "
	     <<m_flavs[m_pi[0]]<<"+"<<m_flavs[m_pi[1]]<<"+"<<m_flavs[m_pi[2]]
	     <<" form factor.\n"
	     <<"   Will exit the run.\n";
  THROW(critical_error,"No form factor for 3 pseudoscalars.")
}

Complex FF_0_PPP_Base::
FF_RChiPT(const double & s123,const double & s12,const double & s13,
          const double & s23) {
  msg_Error()<<"Error in "<<METHOD<<": RChiPT not available for "
	     <<m_flavs[m_pi[0]]<<"+"<<m_flavs[m_pi[1]]<<"+"<<m_flavs[m_pi[2]]
	     <<" form factor.\n"
	     <<"   Will exit the run.\n";
  THROW(critical_error,"No form factor for 3 pseudoscalars.")
}


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//
// Three-pseudoscalar tau currents:
// - pi pi pi channels use a Kuhn-Santamaria-style a1 -> rho pi current.
// - strange PPP channels use the Finkemeier-Mirkes 1995 form factors.
// The RChiPT hook exists in the base interface but is not implemented here.
//
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

class F1_0_PiPlusPiZeroPiZero : public FF_0_PPP_Base {
protected:
  bool    m_isF2;
  double  m_fpi;
  double  m_mpi2, m_FV, m_GV, m_FA;
  double  m_l0, m_lp, m_lpp;
  double  m_mrho, m_Grho, m_mrhop, m_Grhop, m_ma1;
  double  m_beta_rhop, m_r3pi;
  Complex m_alpha, m_gamma, m_delta;

  Total_Width_Base  * p_a1_ks_width;
  Summed_Propagator * p_a1s, * p_rhos;

  void    FixParameters(const FF_Parameters & params);
  void    Construct();
  Complex FF_KS(const double & s123,const double & s12,
                const double & s13,const double & s23);
  Complex FF_RChiPT(const double & s123,const double & s12,
                    const double & s13,const double & s23);
  double  GammaRho(const double & s,const double & mass,
                   const double & width) const;
  double  GammaA1(const double & s) const;
  Complex RhoPropagator(const double & s) const;
  Complex F1RChiPT(const double & q2,const double & s1,
                   const double & s2) const;

public:
  F1_0_PiPlusPiZeroPiZero(const FF_Parameters & params);
  ~F1_0_PiPlusPiZeroPiZero();
};

F1_0_PiPlusPiZeroPiZero::F1_0_PiPlusPiZeroPiZero(const FF_Parameters & params) :
  FF_0_PPP_Base(params),
  p_a1_ks_width(NULL),
  p_a1s(NULL),
  p_rhos(NULL),
  m_isF2(false),
  m_fpi((*params.p_model)("fpi",0.1307)/sqrt(2.)),
  m_mpi2(sqr((2.*Flavour(kf_pi_plus).HadMass()+Flavour(kf_pi).HadMass())/3.))
{
  if (params.m_name=="F2_0_PPP") m_isF2 = true;
  FixParameters(params);
  Construct();
}

F1_0_PiPlusPiZeroPiZero::~F1_0_PiPlusPiZeroPiZero() {
  if (p_a1s)         { delete p_a1s;         p_a1s         = NULL; }
  if (p_rhos)        { delete p_rhos;        p_rhos        = NULL; }
  if (p_a1_ks_width) { delete p_a1_ks_width; p_a1_ks_width = NULL; }
}

void F1_0_PiPlusPiZeroPiZero::FixParameters(const FF_Parameters & params) {
  if (m_mode==FF_0_PPP_mode::piP_pi0_pi0 ||
      m_mode==FF_0_PPP_mode::piM_piP_piP) {
    if (m_ffmodel==ff_model::RChiPT) {
      m_norm = (*params.p_model)("Vud", Tools::Vud)/m_fpi;
    }
    else {
      m_norm = Complex(0., -((2.*sqrt(2)*(*params.p_model)("Vud", Tools::Vud)) /
			     (3.*m_fpi)));
    }

    if (m_ffmodel==ff_model::KS) {
      m_gamma = Complex((*params.p_model)("KS3Pi_gamma_rho1450", -0.14500),
                        (*params.p_model)("KS3Pi_gamma_rho1450_im", 0.0000));
      m_delta = Complex((*params.p_model)("KS3Pi_delta_rho1700",   0.00000),
                        (*params.p_model)("KS3Pi_delta_rho1700_im", 0.0000));
      m_alpha = Complex((*params.p_model)("KS3Pi_alpha_rho_omega",  0.00185),
                        (*params.p_model)("KS3Pi_alpha_rho_omega_im", 0.0000));
    }
    else if (m_ffmodel==ff_model::RChiPT) {
      m_FV = (*params.p_model)("RChiPT_FV", sqrt(2.)*m_fpi);
      m_GV = (*params.p_model)("RChiPT_GV", sqr(m_fpi)/m_FV);
      m_FA = (*params.p_model)("RChiPT_FA", m_fpi);
      m_lp = (*params.p_model)("RChiPT_lambdap",
                               sqr(m_fpi)/(2.*sqrt(2.)*m_FA*m_GV));
      m_lpp = (*params.p_model)("RChiPT_lambdapp",
                                -(1.-2.*sqr(m_GV)/sqr(m_fpi))*m_lp);
      m_l0 = (*params.p_model)("RChiPT_lambda0", 0.25*(m_lp+m_lpp));

      m_mrho = (*params.p_model)("RChiPT_m_rho", 0.77554);
      m_Grho = (*params.p_model)("RChiPT_Gamma_rho", 0.1491);
      m_mrhop = (*params.p_model)("RChiPT_m_rhop", 1.465);
      m_Grhop = (*params.p_model)("RChiPT_Gamma_rhop", 0.4);
      m_ma1 = (*params.p_model)("RChiPT_m_a1", 1.12);
      m_beta_rhop = (*params.p_model)("RChiPT_beta_rhop", -0.145);
      m_r3pi = (m_mode==FF_0_PPP_mode::piM_piP_piP) ? 1. : -1.;
    }
  }
}

void F1_0_PiPlusPiZeroPiZero::Construct() {
  if (m_ffmodel==ff_model::RChiPT) {
    p_a1_ks_width = new KS_A1_1260_plus_Width();
  }
  if (m_ffmodel==ff_model::KS) {
    p_a1_ks_width = new KS_A1_1260_plus_Width();
    if (m_mode==FF_0_PPP_mode::piP_pi0_pi0) {
      Propagator_Base * a11260 =
        new BreitWigner(p_a1_ks_width);
      Propagator_Base * rho770 =
        new BreitWigner(LineShapes->Get(Flavour(kf_rho_770_plus)));
      Propagator_Base * rho1450 =
        new BreitWigner(LineShapes->Get(Flavour(kf_rho_1450_plus)));
      Propagator_Base * rho1700 =
        new BreitWigner(LineShapes->Get(Flavour(kf_rho_1700_plus)));

      p_a1s = new Summed_Propagator();
      p_a1s->Add(a11260, Complex(1.,0.));

      p_rhos = new Summed_Propagator();
      p_rhos->Add(rho770,  Complex(1.,0.));
      p_rhos->Add(rho1450, m_gamma);
      p_rhos->Add(rho1700, m_delta);
    }
    else if (m_mode==FF_0_PPP_mode::piM_piP_piP) {
      Propagator_Base * a11260 =
        new BreitWigner(p_a1_ks_width);
      Propagator_Base * rho770 =
        new BreitWigner(LineShapes->Get(Flavour(kf_rho_770)));
      Propagator_Base * rho770_1 =
        new BreitWigner(LineShapes->Get(Flavour(kf_rho_770)));
      Propagator_Base * omega782 =
        new BreitWigner(LineShapes->Get(Flavour(kf_omega_782)));
      Propagator_Base * rho1450 =
        new BreitWigner(LineShapes->Get(Flavour(kf_rho_1450)));
      Propagator_Base * rho1700 =
        new BreitWigner(LineShapes->Get(Flavour(kf_rho_1700)));

      p_a1s = new Summed_Propagator();
      p_a1s->Add(a11260, Complex(1.,0.));

      Summed_Propagator * rho = new Summed_Propagator();
      rho->Add(rho770, Complex(1.,0.));

      Multiplied_Propagator * inter = new Multiplied_Propagator();
      inter->Add(rho770_1, Complex(1.,0.));
      inter->Add(omega782, Complex(1.,0.));

      rho->Add(inter, m_alpha);

      p_rhos = new Summed_Propagator();
      p_rhos->Add(rho,     Complex(1.,0.));
      p_rhos->Add(rho1450, m_gamma);
      p_rhos->Add(rho1700, m_delta);
    }
  }
}

double F1_0_PiPlusPiZeroPiZero::
GammaRho(const double & s,const double & mass,const double & width) const {
  return Tools::OffShellMassWidth(s,sqr(mass),width,m_mpi2);
}

double F1_0_PiPlusPiZeroPiZero::GammaA1(const double & s) const {
  if (p_a1_ks_width==NULL || s<=0.) return m_ma1*Flavour(kf_a_1_1260_plus).Width();
  return sqrt(s)*(*p_a1_ks_width)(s);
}

Complex F1_0_PiPlusPiZeroPiZero::RhoPropagator(const double & s) const {
  Complex rho = 1./Complex(s-sqr(m_mrho),-GammaRho(s,m_mrho,m_Grho));
  Complex rhop = 1./Complex(s-sqr(m_mrhop),-GammaRho(s,m_mrhop,m_Grhop));
  return (rho+m_beta_rhop*rhop)/(1.+m_beta_rhop);
}

Complex F1_0_PiPlusPiZeroPiZero::
F1RChiPT(const double & q2,const double & s1,const double & s2) const {
  const double s3(q2-s1-s2+3.*m_mpi2);
  const Complex rho1(RhoPropagator(s1)), rho2(RhoPropagator(s2));
  const double gvfv(2.*m_GV/m_FV-1.);

  const Complex fchi = -2.*sqrt(2.)/3.;
  const Complex fR = sqrt(2.)*m_FV*m_GV/(3.*sqr(m_fpi)) *
    (3.*s1*rho1 - gvfv*((2.*q2-2.*s1-s3)*rho1 + (s3-s1)*rho2));

  const Complex h1 = -m_l0*m_mpi2/q2 + m_lp*s1/q2 + m_lpp;
  const Complex h2 = -m_l0*m_mpi2/q2 + m_lp*s2/q2 + m_lpp;
  const Complex fa = q2/Complex(q2-sqr(m_ma1),-GammaA1(q2));
  const Complex fRR = 4.*m_FA*m_GV/(3.*sqr(m_fpi)) * fa *
    (-(m_lp+m_lpp)*3.*s1*rho1 +
     h1*(2.*q2+s1-s3)*rho1 + h2*(s3-s1)*rho2);

  return m_r3pi*(fchi+fR+fRR);
}

Complex F1_0_PiPlusPiZeroPiZero::
FF_KS(const double & s123,const double & s12,const double & s13,
      const double & s23) {
  if (p_a1s==NULL || p_rhos==NULL) return Complex(0.,0.);

  const Complex a1 = (*p_a1s)(s123);
  // KS writes the axial current in the basis
  // (p1-p3)_T F1_KS + (p2-p3)_T F2_KS.  VA_0_PiPiPi contracts
  // form factors with (p2-p1)_T and (p3-p1)_T instead, giving
  // F_current1 = F2_KS and F_current2 = -F1_KS - F2_KS.
  const Complex F1_KS = a1 * (*p_rhos)(s13);
  const Complex F2_KS = a1 * (*p_rhos)(s23);

  if (m_isF2) return -F1_KS-F2_KS;
  return F2_KS;
}

Complex F1_0_PiPlusPiZeroPiZero::
FF_RChiPT(const double & s123,const double & s12,const double & s13,
          const double & s23) {
  const double s1(s23), s2(s13);
  const Complex F1 = F1RChiPT(s123,s1,s2);
  const Complex F2 = F1RChiPT(s123,s2,s1);
  // Sherpa's basis vectors: $v_1 = (p_2 - p_1)_T$ and $v_2 = (p_3 - p_1)_T$.
  // Paper's basis vectors: $u_1 = (p_2 - p_3)_T$ and $u_2 = (p_3 - p_1)_T$.
  return m_isF2 ? -F1-F2 : F1;
}


//////////////////////////////////////////////////////////////////////////////////
//
// Finkemeier-Mirkes, hep-ph/9503474, Table I axial form factors.
// Their basis is (q1-q3)_T F1 + (q2-q3)_T F2, while VA_0_PiPiPi uses
// (q2-q1)_T and (q3-q1)_T.  Hence F_current1 = F_paper2 and
// F_current2 = -F_paper1 - F_paper2.
//
//////////////////////////////////////////////////////////////////////////////////

class F1_0_FM95 : public FF_0_PPP_Base {
  bool    m_isF2;
  double  m_fpi, m_mpi2, m_mK2;
  double  m_beta_rho, m_beta_Kstar, m_xi_K1;
  double  m_mrho, m_Grho, m_mrhop, m_Grhop;
  double  m_mKstar, m_GKstar, m_mKstarp, m_GKstarp;
  double  m_mK1_1400, m_GK1_1400, m_mK1_1270, m_GK1_1270;

  Total_Width_Base * p_a1_ks_width;
  Propagator_Base  * p_a1;

  void    FixParameters(const FF_Parameters & params);
  void    Construct();
  bool    UsesA1() const;
  double  Acoeff() const;
  Complex BW_Running(const double & s,const double & mass,const double & width,
                     const double & m12,const double & m22) const;
  Complex BW_Fixed(const double & s,const double & mass,const double & width) const;
  Complex A1(const double & s) const;
  Complex K1_a(const double & s) const;
  Complex K1_b(const double & s) const;
  Complex T_rho_1(const double & s) const;
  Complex T_Kstar_1(const double & s) const;
  Complex G1(const double & Q2,const double & s1,const double & s2,
             const double & s3) const;
  Complex G2(const double & Q2,const double & s1,const double & s2,
             const double & s3) const;
  Complex FF_KS(const double & s123,const double & s12,
                const double & s13,const double & s23);

public:
  F1_0_FM95(const FF_Parameters & params);
  ~F1_0_FM95();
};

F1_0_FM95::F1_0_FM95(const FF_Parameters & params) :
  FF_0_PPP_Base(params),
  m_isF2(false),
  m_fpi((*params.p_model)("fpi",0.1307)/sqrt(2.)),
  m_mpi2(sqr(Flavour(kf_pi_plus).HadMass())),
  m_mK2(sqr(Flavour(kf_K_plus).HadMass())),
  p_a1_ks_width(NULL),
  p_a1(NULL)
{
  if (params.m_name=="F2_0_PPP") m_isF2 = true;
  m_norm = Complex(1.,0.);
  FixParameters(params);
  Construct();
}

F1_0_FM95::~F1_0_FM95() {
  if (p_a1)          { delete p_a1;          p_a1          = NULL; }
  if (p_a1_ks_width) { delete p_a1_ks_width; p_a1_ks_width = NULL; }
}

void F1_0_FM95::FixParameters(const FF_Parameters & params) {
  m_beta_rho   = (*params.p_model)("FM95_beta_rho",   -0.145);
  m_beta_Kstar = (*params.p_model)("FM95_beta_Kstar", -0.135);
  m_xi_K1      = (*params.p_model)("FM95_xi_K1",       0.330);

  m_mrho  = (*params.p_model)("FM95_m_rho",       0.773);
  m_Grho  = (*params.p_model)("FM95_Gamma_rho",   0.145);
  m_mrhop = (*params.p_model)("FM95_m_rhop",      1.370);
  m_Grhop = (*params.p_model)("FM95_Gamma_rhop",  0.510);

  m_mKstar  = (*params.p_model)("FM95_m_Kstar",      0.892);
  m_GKstar  = (*params.p_model)("FM95_Gamma_Kstar",  0.050);
  m_mKstarp = (*params.p_model)("FM95_m_Kstarp",     1.412);
  m_GKstarp = (*params.p_model)("FM95_Gamma_Kstarp", 0.227);

  m_mK1_1400 = (*params.p_model)("FM95_m_K1_1400",     1.463);
  m_GK1_1400 = (*params.p_model)("FM95_Gamma_K1_1400", 0.300);
  m_mK1_1270 = (*params.p_model)("FM95_m_K1_1270",     1.254);
  m_GK1_1270 = (*params.p_model)("FM95_Gamma_K1_1270", 0.260);
}

void F1_0_FM95::Construct() {
  if (UsesA1()) {
    p_a1_ks_width = new KS_A1_1260_plus_Width();
    p_a1          = new BreitWigner(p_a1_ks_width);
  }
}

bool F1_0_FM95::UsesA1() const {
  return m_mode==FF_0_PPP_mode::K_pi_K      ||
         m_mode==FF_0_PPP_mode::K0_pi_K0b   ||
         m_mode==FF_0_PPP_mode::KS_pi_KS    ||
         m_mode==FF_0_PPP_mode::KS_pi_KL    ||
         m_mode==FF_0_PPP_mode::KL_pi_KL    ||
         m_mode==FF_0_PPP_mode::K_pi0_K0;
}

double F1_0_FM95::Acoeff() const {
  switch (m_mode) {
  case FF_0_PPP_mode::K_pi_K:
  case FF_0_PPP_mode::K0_pi_K0b:
    return -(*p_model)("Vud", Tools::Vud)/2.;
  case FF_0_PPP_mode::KS_pi_KS:
  case FF_0_PPP_mode::KS_pi_KL:
    return -(*p_model)("Vud", Tools::Vud)/4.;
  case FF_0_PPP_mode::KL_pi_KL:
    return  (*p_model)("Vud", Tools::Vud)/4.;
  case FF_0_PPP_mode::K_pi0_K0:
    return  3.*(*p_model)("Vud", Tools::Vud)/(2.*sqrt(2.));
  case FF_0_PPP_mode::pi0_pi0_K:
    return  (*p_model)("Vus", Tools::Vus)/4.;
  case FF_0_PPP_mode::K_pi_pi:
    return -(*p_model)("Vus", Tools::Vus)/2.;
  case FF_0_PPP_mode::pi_K0b_pi0:
    return  3.*(*p_model)("Vus", Tools::Vus)/(2.*sqrt(2.));
  default:
    return 0.;
  }
}

Complex F1_0_FM95::BW_Running(const double & s,const double & mass,
                              const double & width,const double & m12,
                              const double & m22) const {
  const double mass2 = sqr(mass);
  return Tools::BreitWigner(s,mass2,
                            Tools::OffShellMassWidth(s,mass2,width,m12,m22));
}

Complex F1_0_FM95::BW_Fixed(const double & s,const double & mass,
                            const double & width) const {
  return Tools::BreitWignerFix(s,sqr(mass),mass*width);
}

Complex F1_0_FM95::A1(const double & s) const {
  return p_a1 ? (*p_a1)(s) : Complex(0.,0.);
}

Complex F1_0_FM95::K1_a(const double & s) const {
  return (BW_Fixed(s,m_mK1_1400,m_GK1_1400) +
          m_xi_K1*BW_Fixed(s,m_mK1_1270,m_GK1_1270))/(1.+m_xi_K1);
}

Complex F1_0_FM95::K1_b(const double & s) const {
  return BW_Fixed(s,m_mK1_1270,m_GK1_1270);
}

Complex F1_0_FM95::T_rho_1(const double & s) const {
  return (BW_Running(s,m_mrho,m_Grho,m_mpi2,m_mpi2) +
          m_beta_rho*BW_Running(s,m_mrhop,m_Grhop,m_mpi2,m_mpi2)) /
         (1.+m_beta_rho);
}

Complex F1_0_FM95::T_Kstar_1(const double & s) const {
  return (BW_Running(s,m_mKstar,m_GKstar,m_mpi2,m_mK2) +
          m_beta_Kstar*BW_Running(s,m_mKstarp,m_GKstarp,m_mpi2,m_mK2)) /
         (1.+m_beta_Kstar);
}

Complex F1_0_FM95::G1(const double & Q2,const double & s1,
                      const double & s2,const double & s3) const {
  switch (m_mode) {
  case FF_0_PPP_mode::K_pi_K:
  case FF_0_PPP_mode::K0_pi_K0b:
    return A1(Q2)*T_rho_1(s2);
  case FF_0_PPP_mode::KS_pi_KS:
  case FF_0_PPP_mode::KL_pi_KL:
    return A1(Q2)*T_Kstar_1(s3);
  case FF_0_PPP_mode::KS_pi_KL:
    return A1(Q2)*(2.*T_rho_1(s2)+T_Kstar_1(s3));
  case FF_0_PPP_mode::K_pi0_K0:
    return A1(Q2)*(2./3.*T_rho_1(s2)+1./3.*T_Kstar_1(s3));
  case FF_0_PPP_mode::pi0_pi0_K:
    return K1_a(Q2)*T_Kstar_1(s2);
  case FF_0_PPP_mode::K_pi_pi:
    return K1_a(Q2)*T_Kstar_1(s2);
  case FF_0_PPP_mode::pi_K0b_pi0:
    return 2./3.*K1_b(Q2)*T_rho_1(s2)+
           1./3.*K1_a(Q2)*T_Kstar_1(s3);
  default:
    return Complex(0.,0.);
  }
}

Complex F1_0_FM95::G2(const double & Q2,const double & s1,
                      const double & s2,const double & s3) const {
  switch (m_mode) {
  case FF_0_PPP_mode::K_pi_K:
  case FF_0_PPP_mode::K0_pi_K0b:
    return A1(Q2)*T_Kstar_1(s1);
  case FF_0_PPP_mode::KS_pi_KS:
  case FF_0_PPP_mode::KL_pi_KL:
    return -A1(Q2)*(T_Kstar_1(s1)+T_Kstar_1(s3));
  case FF_0_PPP_mode::KS_pi_KL:
    return A1(Q2)*(T_Kstar_1(s1)-T_Kstar_1(s3));
  case FF_0_PPP_mode::K_pi0_K0:
    return 1./3.*A1(Q2)*(T_Kstar_1(s1)-T_Kstar_1(s3));
  case FF_0_PPP_mode::pi0_pi0_K:
    return K1_a(Q2)*T_Kstar_1(s1);
  case FF_0_PPP_mode::K_pi_pi:
    return K1_b(Q2)*T_rho_1(s1);
  case FF_0_PPP_mode::pi_K0b_pi0:
    return 1./3.*K1_a(Q2)*(T_Kstar_1(s1)-T_Kstar_1(s3));
  default:
    return Complex(0.,0.);
  }
}

Complex F1_0_FM95::
FF_KS(const double & s123,const double & s12,const double & s13,
      const double & s23) {
  const double s1 = s23, s2 = s13, s3 = s12;
  const Complex coeff = 2.*sqrt(2.)*Acoeff()/(3.*m_fpi);
  const Complex F1_paper = coeff*G1(s123,s1,s2,s3);
  const Complex F2_paper = coeff*G2(s123,s1,s2,s3);

  if (m_isF2) return -F1_paper-F2_paper;
  return F2_paper;
}


//////////////////////////////////////////////////////////////////////////////////
//
// Form factors for
//   * pi pi pi vanish due to identical masses
// - Kuehn-Santamaria
//   * K pi: Z.Phys.C 72 (1996) 619 (scalar form factor)
//
// Todo: implement all scalar form factors here!
//
//////////////////////////////////////////////////////////////////////////////////

class F3_0_PiPlusPiZeroPiZero : public FF_0_PPP_Base {
  double  m_fpi;
  double  m_mpi2, m_FV, m_GV;
  double  m_mrho, m_Grho, m_mrhop, m_Grhop;
  double  m_beta_rhop, m_r3pi, m_kappa;
  Complex m_gamma, m_delta;

  void    FixParameters(const FF_Parameters & params);
  void    Construct();
  double  GammaRho(const double & s,const double & mass,
                   const double & width) const;
  Complex RhoPropagator(const double & s) const;
  Complex Alpha2(const double & q2,const double & s1,
                 const double & s2) const;

  Complex FF_KS(const double & s123,const double & s12,
                const double & s13,const double & s23) {
    return Complex(0.,0.);
  }

  Complex FF_RChiPT(const double & s123,const double & s12,
                    const double & s13,const double & s23);

public:
  F3_0_PiPlusPiZeroPiZero(const FF_Parameters & params) :
    FF_0_PPP_Base(params),
    m_fpi((*params.p_model)("fpi",0.1307)/sqrt(2.)),
    m_mpi2(sqr((2.*Flavour(kf_pi_plus).HadMass()+Flavour(kf_pi).HadMass())/3.))
  {
    FixParameters(params);
    Construct();
  }
};

void F3_0_PiPlusPiZeroPiZero::FixParameters(const FF_Parameters & params) {
  if (m_ffmodel==ff_model::RChiPT &&
      (m_mode==FF_0_PPP_mode::piP_pi0_pi0 ||
       m_mode==FF_0_PPP_mode::piM_piP_piP)) {
    m_norm = (*params.p_model)("Vud", Tools::Vud)/m_fpi;
    m_FV = (*params.p_model)("RChiPT_FV", sqrt(2.)*m_fpi);
    m_GV = (*params.p_model)("RChiPT_GV", sqr(m_fpi)/m_FV);
    m_mrho = (*params.p_model)("RChiPT_m_rho", 0.77554);
    m_Grho = (*params.p_model)("RChiPT_Gamma_rho", 0.1491);
    m_mrhop = (*params.p_model)("RChiPT_m_rhop", 1.465);
    m_Grhop = (*params.p_model)("RChiPT_Gamma_rhop", 0.4);
    m_beta_rhop = (*params.p_model)("RChiPT_beta_rhop", -0.145);
    m_r3pi = (m_mode==FF_0_PPP_mode::piM_piP_piP) ? 1. : -1.;
    m_kappa = (m_mode==FF_0_PPP_mode::piM_piP_piP) ? 1. : 0.5;
  }
}

void F3_0_PiPlusPiZeroPiZero::Construct() {}

double F3_0_PiPlusPiZeroPiZero::
GammaRho(const double & s,const double & mass,const double & width) const {
  return Tools::OffShellMassWidth(s,sqr(mass),width,m_mpi2);
}

Complex F3_0_PiPlusPiZeroPiZero::RhoPropagator(const double & s) const {
  Complex rho = 1./Complex(s-sqr(m_mrho),-GammaRho(s,m_mrho,m_Grho));
  Complex rhop = 1./Complex(s-sqr(m_mrhop),-GammaRho(s,m_mrhop,m_Grhop));
  return (rho+m_beta_rhop*rhop)/(1.+m_beta_rhop);
}

Complex F3_0_PiPlusPiZeroPiZero::
Alpha2(const double & q2,const double & s1,const double & s2) const {
  const double s3(q2-s1-s2+3.*m_mpi2);
  return 3.*m_GV/m_FV * s1/q2 * m_mpi2/(q2-m_mpi2) * (s3-s2) * RhoPropagator(s1);
}

Complex F3_0_PiPlusPiZeroPiZero::
FF_RChiPT(const double & s123,const double & s12,const double & s13,
          const double & s23) {
  const double s1(s23), s2(s13);
  const double s3(s12);
  const Complex f4chi = 2.*sqrt(2.)/3. * m_mpi2 *
    (3.*(s3-m_mpi2)-s123*(1.+2.*m_kappa*m_r3pi)) / (2.*s123*(s123-m_mpi2));
  const Complex f4R = -sqrt(2.)*m_FV*m_GV/(3.*sqr(m_fpi)) *
    (Alpha2(s123,s2,s1)+Alpha2(s123,s1,s2));
  return m_r3pi*(f4chi+f4R);
}


//////////////////////////////////////////////////////////////////////////////////
//
// Form factors for
//   * pi pi pi vanish due to identical masses
// - Kuehn-Santamaria
//   * K pi: Z.Phys.C 72 (1996) 619 (scalar form factor)
//
// Todo: implement all scalar form factors here!
//
//////////////////////////////////////////////////////////////////////////////////

class FS_0_PiPlusPiZeroPiZero : public FF_0_PPP_Base {
  double  m_fpi;
  Complex m_gamma, m_delta;

  void    FixParameters(const FF_Parameters & params);
  void    Construct();

  Complex FF_KS(const double & s123,const double & s12,
                const double & s13,const double & s23) {
    return Complex(0.,0.);
  }

  Complex FF_RChiPT(const double & s123,const double & s12,
                    const double & s13,const double & s23) {
    return Complex(0.,0.);
  }

public:
  FS_0_PiPlusPiZeroPiZero(const FF_Parameters & params) :
    FF_0_PPP_Base(params),
    m_fpi((*params.p_model)("fpi",0.1307)/sqrt(2.))
  {
    FixParameters(params);
    Construct();
  }
};

void FS_0_PiPlusPiZeroPiZero::FixParameters(const FF_Parameters & params) {}

void FS_0_PiPlusPiZeroPiZero::Construct() {}


//////////////////////////////////////////////////////////////////////////////////
//
// Finkemeier-Mirkes, hep-ph/9503474, Table II vector form factor.
// In the paper the vector-current term is i epsilon(q1,q2,q3) F3.  The
// METOOLS current multiplies FS directly by epsilon(q1,q2,q3), so this class
// returns i*F3_paper.
//
//////////////////////////////////////////////////////////////////////////////////

class FS_0_FM95 : public FF_0_PPP_Base {
  double  m_fpi, m_mpi2, m_mK2;
  double  m_beta_rho, m_beta_Kstar, m_lambda, m_mu, m_epsilon;
  double  m_mrho1, m_Grho1, m_mrhop1, m_Grhop1;
  double  m_mKstar1, m_GKstar1, m_mKstarp1, m_GKstarp1;
  double  m_mrho2, m_Grho2, m_mrhop2, m_Grhop2, m_mrhopp2, m_Grhopp2;
  double  m_mKstar2, m_GKstar2, m_mKstarp2, m_GKstarp2, m_mKstarpp2, m_GKstarpp2;
  double  m_momega, m_Gomega, m_mphi, m_Gphi;

  void    FixParameters(const FF_Parameters & params);
  double  Acoeff() const;
  Complex BW_Running(const double & s,const double & mass,const double & width,
                     const double & m12,const double & m22) const;
  Complex BW_Fixed(const double & s,const double & mass,const double & width) const;
  Complex T_rho_1(const double & s) const;
  Complex T_Kstar_1(const double & s) const;
  Complex T_rho_2(const double & s) const;
  Complex T_Kstar_2(const double & s) const;
  Complex T_omega(const double & s) const;
  Complex G3(const double & Q2,const double & s1,const double & s2,
             const double & s3) const;
  Complex FF_KS(const double & s123,const double & s12,
                const double & s13,const double & s23);

public:
  FS_0_FM95(const FF_Parameters & params);
};

FS_0_FM95::FS_0_FM95(const FF_Parameters & params) :
  FF_0_PPP_Base(params),
  m_fpi((*params.p_model)("fpi",0.1307)/sqrt(2.)),
  m_mpi2(sqr(Flavour(kf_pi_plus).HadMass())),
  m_mK2(sqr(Flavour(kf_K_plus).HadMass()))
{
  m_norm = Complex(1.,0.);
  FixParameters(params);
}

void FS_0_FM95::FixParameters(const FF_Parameters & params) {
  m_beta_rho   = (*params.p_model)("FM95_beta_rho",   -0.145);
  m_beta_Kstar = (*params.p_model)("FM95_beta_Kstar", -0.135);
  m_lambda     = (*params.p_model)("FM95_lambda_vector", -0.250);
  m_mu         = (*params.p_model)("FM95_mu_vector",     -0.038);
  m_epsilon    = (*params.p_model)("FM95_epsilon_phi",    0.050);

  m_mrho1  = (*params.p_model)("FM95_m_rho",       0.773);
  m_Grho1  = (*params.p_model)("FM95_Gamma_rho",   0.145);
  m_mrhop1 = (*params.p_model)("FM95_m_rhop",      1.370);
  m_Grhop1 = (*params.p_model)("FM95_Gamma_rhop",  0.510);

  m_mKstar1  = (*params.p_model)("FM95_m_Kstar",      0.892);
  m_GKstar1  = (*params.p_model)("FM95_Gamma_Kstar",  0.050);
  m_mKstarp1 = (*params.p_model)("FM95_m_Kstarp",     1.412);
  m_GKstarp1 = (*params.p_model)("FM95_Gamma_Kstarp", 0.227);

  m_mrho2   = (*params.p_model)("FM95_m_rho_vector",       0.773);
  m_Grho2   = (*params.p_model)("FM95_Gamma_rho_vector",   0.145);
  m_mrhop2  = (*params.p_model)("FM95_m_rhop_vector",      1.500);
  m_Grhop2  = (*params.p_model)("FM95_Gamma_rhop_vector",  0.220);
  m_mrhopp2 = (*params.p_model)("FM95_m_rhopp_vector",     1.750);
  m_Grhopp2 = (*params.p_model)("FM95_Gamma_rhopp_vector", 0.120);

  m_mKstar2   = (*params.p_model)("FM95_m_Kstar_vector",       0.892);
  m_GKstar2   = (*params.p_model)("FM95_Gamma_Kstar_vector",   0.050);
  m_mKstarp2  = (*params.p_model)("FM95_m_Kstarp_vector",      1.412);
  m_GKstarp2  = (*params.p_model)("FM95_Gamma_Kstarp_vector",  0.227);
  m_mKstarpp2 = (*params.p_model)("FM95_m_Kstarpp_vector",     1.714);
  m_GKstarpp2 = (*params.p_model)("FM95_Gamma_Kstarpp_vector", 0.323);

  m_momega = (*params.p_model)("FM95_m_omega",     0.782);
  m_Gomega = (*params.p_model)("FM95_Gamma_omega", 0.00843);
  m_mphi   = (*params.p_model)("FM95_m_phi",       1.020);
  m_Gphi   = (*params.p_model)("FM95_Gamma_phi",   0.00443);
}

double FS_0_FM95::Acoeff() const {
  switch (m_mode) {
  case FF_0_PPP_mode::K_pi_K:
    return -(*p_model)("Vud", Tools::Vud);
  case FF_0_PPP_mode::K0_pi_K0b:
    return  (*p_model)("Vud", Tools::Vud);
  case FF_0_PPP_mode::KS_pi_KS:
    return -(*p_model)("Vud", Tools::Vud)/2.;
  case FF_0_PPP_mode::KL_pi_KL:
    return  (*p_model)("Vud", Tools::Vud)/2.;
  case FF_0_PPP_mode::KS_pi_KL:
    return  (*p_model)("Vud", Tools::Vud)/2.;
  case FF_0_PPP_mode::K_pi0_K0:
    return -(*p_model)("Vud", Tools::Vud)/sqrt(2.);
  case FF_0_PPP_mode::pi0_pi0_K:
  case FF_0_PPP_mode::K_pi_pi:
    return  (*p_model)("Vus", Tools::Vus);
  case FF_0_PPP_mode::pi_K0b_pi0:
    return  sqrt(2.)*(*p_model)("Vus", Tools::Vus);
  default:
    return 0.;
  }
}

Complex FS_0_FM95::BW_Running(const double & s,const double & mass,
                              const double & width,const double & m12,
                              const double & m22) const {
  const double mass2 = sqr(mass);
  return Tools::BreitWigner(s,mass2,
                            Tools::OffShellMassWidth(s,mass2,width,m12,m22));
}

Complex FS_0_FM95::BW_Fixed(const double & s,const double & mass,
                            const double & width) const {
  return Tools::BreitWignerFix(s,sqr(mass),mass*width);
}

Complex FS_0_FM95::T_rho_1(const double & s) const {
  return (BW_Running(s,m_mrho1,m_Grho1,m_mpi2,m_mpi2) +
          m_beta_rho*BW_Running(s,m_mrhop1,m_Grhop1,m_mpi2,m_mpi2)) /
         (1.+m_beta_rho);
}

Complex FS_0_FM95::T_Kstar_1(const double & s) const {
  return (BW_Running(s,m_mKstar1,m_GKstar1,m_mpi2,m_mK2) +
          m_beta_Kstar*BW_Running(s,m_mKstarp1,m_GKstarp1,m_mpi2,m_mK2)) /
         (1.+m_beta_Kstar);
}

Complex FS_0_FM95::T_rho_2(const double & s) const {
  return (BW_Running(s,m_mrho2,m_Grho2,m_mpi2,m_mpi2) +
          m_lambda*BW_Running(s,m_mrhop2,m_Grhop2,m_mpi2,m_mpi2) +
          m_mu*BW_Running(s,m_mrhopp2,m_Grhopp2,m_mpi2,m_mpi2)) /
         (1.+m_lambda+m_mu);
}

Complex FS_0_FM95::T_Kstar_2(const double & s) const {
  return (BW_Running(s,m_mKstar2,m_GKstar2,m_mpi2,m_mK2) +
          m_lambda*BW_Running(s,m_mKstarp2,m_GKstarp2,m_mpi2,m_mK2) +
          m_mu*BW_Running(s,m_mKstarpp2,m_GKstarpp2,m_mpi2,m_mK2)) /
         (1.+m_lambda+m_mu);
}

Complex FS_0_FM95::T_omega(const double & s) const {
  return (BW_Fixed(s,m_momega,m_Gomega)+
          m_epsilon*BW_Fixed(s,m_mphi,m_Gphi))/(1.+m_epsilon);
}

Complex FS_0_FM95::G3(const double & Q2,const double & s1,
                      const double & s2,const double & s3) const {
  const double root2 = sqrt(2.);

  switch (m_mode) {
  case FF_0_PPP_mode::K_pi_K:
  case FF_0_PPP_mode::K0_pi_K0b:
    return T_rho_2(Q2)*(root2-1.)*(root2*T_omega(s2)+T_Kstar_1(s1));
  case FF_0_PPP_mode::KS_pi_KS:
  case FF_0_PPP_mode::KL_pi_KL:
    return T_rho_2(Q2)*(root2-1.)*(T_Kstar_1(s1)-T_Kstar_1(s3));
  case FF_0_PPP_mode::KS_pi_KL:
    return T_rho_2(Q2)*(root2-1.)*
           (2.*root2*T_omega(s2)+T_Kstar_1(s1)+T_Kstar_1(s3));
  case FF_0_PPP_mode::K_pi0_K0:
    return T_rho_2(Q2)*(root2-1.)*(T_Kstar_1(s3)-T_Kstar_1(s1));
  case FF_0_PPP_mode::pi0_pi0_K:
    return 1./4.*T_Kstar_2(Q2)*(T_Kstar_1(s1)-T_Kstar_1(s2));
  case FF_0_PPP_mode::K_pi_pi:
    return 1./2.*T_Kstar_2(Q2)*(T_rho_1(s1)+T_Kstar_1(s2));
  case FF_0_PPP_mode::pi_K0b_pi0:
    return 1./4.*T_Kstar_2(Q2)*
           (2.*T_rho_1(s2)+T_Kstar_1(s1)+T_Kstar_1(s3));
  default:
    return Complex(0.,0.);
  }
}

Complex FS_0_FM95::
FF_KS(const double & s123,const double & s12,const double & s13,
      const double & s23) {
  const double s1 = s23, s2 = s13, s3 = s12;
  const Complex coeff = Acoeff()/(2.*sqrt(2.)*sqr(M_PI)*pow(m_fpi,3.));
  return Complex(0.,1.)*coeff*G3(s123,s1,s2,s3);
}


//////////////////////////////////////////////////////////////////////////////////
//
// RChiPT KKpi form factors (arXiv:1203.3955)
//
//////////////////////////////////////////////////////////////////////////////////

class F1_0_RChiPT_KKpi : public FF_0_PPP_Base {
protected:
  bool    m_isF2;
  double  m_fpi;
  double  m_mpi, m_mK, m_meta;
  double  m_mpi2, m_mK2;
  double  m_FV, m_GV, m_FA;
  double  m_lp, m_lpp, m_l0;
  double  m_mrho, m_Grho, m_mrhop, m_Grhop, m_ma1;
  double  m_mKstar, m_GKstar;
  double  m_beta_rhop;

  Total_Width_Base * p_a1_ks_width;

  void    FixParameters(const FF_Parameters & params);
  void    Construct();
  Complex FF_KS(const double & s123,const double & s12,
                const double & s13,const double & s23);
  Complex FF_RChiPT(const double & s123,const double & s12,
                    const double & s13,const double & s23);
  Complex EvaluateF1F2(const double & s123,const double & s1,
                       const double & s2,const double & s3,bool get_F2);
  double  GammaRho(const double & s,const double & mass,
                   const double & width) const;
  double  GammaKstar(const double & s) const;
  double  GammaA1(const double & s) const;
  Complex RhoPropagator(const double & s) const;
  Complex KstarPropagator(const double & s) const;
  Complex A1Propagator(const double & s) const;

  double  KLambda(double a,double b,double c) const;
  Complex H_func(const Complex & x,const Complex & y) const;

  Complex A_R(const double & q2,const double & x,const double & y,
              const double & m12,const double & m22,const double & m32) const;
  Complex B_R(const double & x,const double & y,const double & m12,
              const double & m22) const;
  Complex A_RR(const double & q2,const double & x,const double & y,
               const double & m12,const double & m22,const double & m32) const;
  Complex B_RR(const double & q2,const double & x,const double & y,const double & z,
               const double & m12,const double & m22,const double & m32) const;

public:
  F1_0_RChiPT_KKpi(const FF_Parameters & params);
  ~F1_0_RChiPT_KKpi();
};

F1_0_RChiPT_KKpi::F1_0_RChiPT_KKpi(const FF_Parameters & params) :
  FF_0_PPP_Base(params),
  p_a1_ks_width(NULL),
  m_isF2(false),
  m_fpi((*params.p_model)("fpi",0.1307)/sqrt(2.)),
  m_mpi((Flavour(kf_pi_plus).HadMass()+Flavour(kf_pi).HadMass())/2.),
  m_mK((Flavour(kf_K_plus).HadMass()+Flavour(kf_K).HadMass())/2.),
  m_meta(Flavour(kf_eta).HadMass()),
  m_mpi2(sqr(m_mpi)),
  m_mK2(sqr(m_mK))
{
  if (params.m_name=="F2_0_PPP") m_isF2 = true;
  FixParameters(params);
  Construct();
}

F1_0_RChiPT_KKpi::~F1_0_RChiPT_KKpi() {
  if (p_a1_ks_width) { delete p_a1_ks_width; p_a1_ks_width = NULL; }
}

void F1_0_RChiPT_KKpi::FixParameters(const FF_Parameters & params) {
  m_FV = (*params.p_model)("RChiPT_FV", sqrt(2.)*m_fpi);
  m_GV = (*params.p_model)("RChiPT_GV", sqr(m_fpi)/m_FV);
  m_FA = (*params.p_model)("RChiPT_FA", m_fpi);
  m_lp = (*params.p_model)("RChiPT_lambdap",
                           sqr(m_fpi)/(2.*sqrt(2.)*m_FA*m_GV));
  m_lpp = (*params.p_model)("RChiPT_lambdapp",
                            -(1.-2.*sqr(m_GV)/sqr(m_fpi))*m_lp);
  m_l0 = (*params.p_model)("RChiPT_lambda0", 0.25*(m_lp+m_lpp));

  m_mrho = (*params.p_model)("RChiPT_m_rho", 0.77554);
  m_Grho = (*params.p_model)("RChiPT_Gamma_rho", 0.1491);
  m_mrhop = (*params.p_model)("RChiPT_m_rhop", 1.465);
  m_Grhop = (*params.p_model)("RChiPT_Gamma_rhop", 0.4);
  m_ma1 = (*params.p_model)("RChiPT_m_a1", 1.12);
  m_beta_rhop = (*params.p_model)("RChiPT_beta_rhop", -0.145);
  m_mKstar = (*params.p_model)("RChiPT_m_Kstar", 0.892);
  m_GKstar = (*params.p_model)("RChiPT_Gamma_Kstar", 0.050);

  m_norm = (*params.p_model)("Vud", Tools::Vud)/m_fpi;
}

void F1_0_RChiPT_KKpi::Construct() {
  p_a1_ks_width = new KS_A1_1260_plus_Width();
}

double F1_0_RChiPT_KKpi::KLambda(double a,double b,double c) const {
  return sqr(a-b-c) - 4.*b*c;
}

Complex F1_0_RChiPT_KKpi::H_func(const Complex & x,const Complex & y) const {
  return -m_l0*y + m_lp*x + m_lpp;
}

double F1_0_RChiPT_KKpi::
GammaRho(const double & s,const double & mass,const double & width) const {
  return Tools::OffShellMassWidth(s,sqr(mass),width,m_mpi2);
}

double F1_0_RChiPT_KKpi::GammaKstar(const double & s) const {
  if (s<=sqr(m_mK+m_mpi)) return 0.;
  double lam_Kpi = KLambda(s, m_mK2, m_mpi2);
  double sum = pow(lam_Kpi, 1.5);
  if (s>sqr(m_mK+m_meta)) {
    double lam_Keta = KLambda(s, m_mK2, sqr(m_meta));
    sum += pow(lam_Keta, 1.5);
  }
  sum /= (s*s);

  double lam_Kpi_peak = KLambda(sqr(m_mKstar), m_mK2, m_mpi2);
  double sum_peak = pow(lam_Kpi_peak, 1.5);
  if (sqr(m_mKstar)>sqr(m_mK+m_meta)) {
    double lam_Keta_peak = KLambda(sqr(m_mKstar), m_mK2, sqr(m_meta));
    sum_peak += pow(lam_Keta_peak, 1.5);
  }
  sum_peak /= pow(m_mKstar, 4);

  return m_GKstar * sum / sum_peak;
}

double F1_0_RChiPT_KKpi::GammaA1(const double & s) const {
  if (p_a1_ks_width==NULL || s<=0.) return m_ma1*Flavour(kf_a_1_1260_plus).Width();
  return sqrt(s)*(*p_a1_ks_width)(s);
}

Complex F1_0_RChiPT_KKpi::RhoPropagator(const double & s) const {
  return -1./Complex(s-sqr(m_mrho),-m_mrho*GammaRho(s,m_mrho,m_Grho));
}

Complex F1_0_RChiPT_KKpi::KstarPropagator(const double & s) const {
  return -1./Complex(s-sqr(m_mKstar),-m_mKstar*GammaKstar(s));
}

Complex F1_0_RChiPT_KKpi::A1Propagator(const double & s) const {
  return -1./Complex(s-sqr(m_ma1),-GammaA1(s));
}

Complex F1_0_RChiPT_KKpi::
A_R(const double & q2,const double & x,const double & y,
    const double & m12,const double & m22,const double & m32) const {
  return 3.*x + m12 - m32 + (1.-2.*m_GV/m_FV)*(2.*q2-2.*x-y+m32-m22);
}

Complex F1_0_RChiPT_KKpi::
B_R(const double & x,const double & y,const double & m12,
    const double & m22) const {
  return 2.*(m22-m12) + (1.-2.*m_GV/m_FV)*(y-x+m12-m22);
}

Complex F1_0_RChiPT_KKpi::
A_RR(const double & q2,const double & x,const double & y,
     const double & m12,const double & m22,const double & m32) const {
  return (m_lp+m_lpp)*(-3.*x+m32-m12) + (2.*q2+x-y+m12-m22)*H_func(x/q2,m22/q2);
}

Complex F1_0_RChiPT_KKpi::
B_RR(const double & q2,const double & x,const double & y,const double & z,
     const double & m12,const double & m22,const double & m32) const {
  return 2.*(m_lp+m_lpp)*(m12-m22) + (y-x+m22-m12)*H_func(z/q2,m32/q2);
}

Complex F1_0_RChiPT_KKpi::
FF_KS(const double & s123,const double & s12,const double & s13,
      const double & s23) {
  return Complex(0.,0.);
}

Complex F1_0_RChiPT_KKpi::
EvaluateF1F2(const double & s123,const double & s1,const double & s2,
             const double & s3,bool get_F2) {
  const Complex f1_chi = -sqrt(2.)/3.;
  const Complex rho_denom = -RhoPropagator(s2);
  const Complex kstar_denom = -KstarPropagator(s1);
  const Complex f1_R = -sqrt(2.)/6.*m_FV*m_GV/sqr(m_fpi) *
    (B_R(s1,s3,m_mK2,m_mK2)*rho_denom +
     A_R(s123,s1,s3,m_mK2,m_mK2,m_mpi2)*kstar_denom);

  const Complex a1_denom = -A1Propagator(s123);
  const Complex f1_RR = 2./3.*m_FA*m_GV/sqr(m_fpi)*s123*a1_denom *
    (B_RR(s123,s1,s3,s2,m_mK2,m_mK2,m_mpi2)*rho_denom +
     A_RR(s123,s1,s3,m_mK2,m_mK2,m_mpi2)*kstar_denom);

  const Complex F1_paper = f1_chi + f1_R + f1_RR;

  const Complex f2_chi = f1_chi;
  const Complex f2_R = -sqrt(2.)/6.*m_FV*m_GV/sqr(m_fpi) *
    (A_R(s123,s2,s3,m_mK2,m_mpi2,m_mK2)*rho_denom +
     B_R(s2,s3,m_mK2,m_mpi2)*kstar_denom);
  const Complex f2_RR = 2./3.*m_FA*m_GV/sqr(m_fpi)*s123*a1_denom *
    (A_RR(s123,s2,s3,m_mK2,m_mpi2,m_mK2)*rho_denom +
     B_RR(s123,s2,s3,s1,m_mK2,m_mpi2,m_mK2)*kstar_denom);

  const Complex F2_paper = f2_chi + f2_R + f2_RR;

  if (get_F2) return -F1_paper-F2_paper;
  return F1_paper;
}

Complex F1_0_RChiPT_KKpi::
FF_RChiPT(const double & s123,const double & s12,const double & s13,
          const double & s23) {
  if (m_mode==FF_0_PPP_mode::K_pi_K || m_mode==FF_0_PPP_mode::K0_pi_K0b) {
    const double s1(s23), s2(s13), s3(s12);
    return EvaluateF1F2(s123,s1,s2,s3,m_isF2);
  }
  else if (m_mode==FF_0_PPP_mode::KS_pi_KS) {
    const double s1(s23), s2(s13), s3(s12);
    Complex f1_orig = EvaluateF1F2(s123,s1,s2,s3,false);
    Complex f2_orig = EvaluateF1F2(s123,s1,s2,s3,true);
    Complex f1_swap = EvaluateF1F2(s123,s3,s2,s1,false);
    Complex f2_swap = EvaluateF1F2(s123,s3,s2,s1,true);
    if (m_isF2) return 0.5*(f2_orig-f2_swap);
    return 0.5*(f1_orig+f1_swap+f2_swap);
  }
  else if (m_mode==FF_0_PPP_mode::KS_pi_KL) {
    const double s1(s23), s2(s13), s3(s12);
    Complex f1_orig = EvaluateF1F2(s123,s1,s2,s3,false);
    Complex f2_orig = EvaluateF1F2(s123,s1,s2,s3,true);
    Complex f1_swap = EvaluateF1F2(s123,s3,s2,s1,false);
    Complex f2_swap = EvaluateF1F2(s123,s3,s2,s1,true);
    if (m_isF2) return 0.5*(f2_orig+f2_swap);
    return 0.5*(f1_orig-f1_swap-f2_swap);
  }
  else if (m_mode==FF_0_PPP_mode::KL_pi_KL) {
    const double s1(s23), s2(s13), s3(s12);
    Complex f1_orig = EvaluateF1F2(s123,s1,s2,s3,false);
    Complex f2_orig = EvaluateF1F2(s123,s1,s2,s3,true);
    Complex f1_swap = EvaluateF1F2(s123,s3,s2,s1,false);
    Complex f2_swap = EvaluateF1F2(s123,s3,s2,s1,true);
    if (m_isF2) return -0.5*(f2_orig+f2_swap);
    return -0.5*(f1_orig-f1_swap-f2_swap);
  }
  else if (m_mode==FF_0_PPP_mode::K_pi0_K0) {
    const double s1(s23), s2(s13), s3(s12);

    const Complex f2_chi = -1.;
    const Complex rho_denom = -RhoPropagator(s2);
    const Complex kstar3_denom = -KstarPropagator(s3);
    const Complex kstar1_denom = -KstarPropagator(s1);

    const Complex f2_R = -1./6.*m_FV*m_GV/sqr(m_fpi) *
      (B_R(s2,s1,m_mK2,m_mpi2)*kstar3_denom +
       2.*A_R(s123,s2,s1,m_mK2,m_mpi2,m_mK2)*rho_denom +
       A_R(s123,s1,s2,m_mpi2,m_mK2,m_mK2)*kstar1_denom);

    const Complex a1_denom = -A1Propagator(s123);
    const Complex f2_RR = sqrt(2.)/3.*m_FA*m_GV/sqr(m_fpi)*s123*a1_denom *
      (B_RR(s123,s2,s1,s3,m_mK2,m_mpi2,m_mK2)*kstar3_denom +
       2.*A_RR(s123,s2,s1,m_mK2,m_mpi2,m_mK2)*rho_denom +
       A_RR(s123,s1,s2,m_mpi2,m_mK2,m_mK2)*kstar1_denom);

    const Complex F2_paper = f2_chi + f2_R + f2_RR;

    const Complex f3_chi = 0.;
    const Complex f3_R = -1./6.*m_FV*m_GV/sqr(m_fpi) *
      (A_R(s123,s3,s1,m_mK2,m_mK2,m_mpi2)*kstar3_denom +
       2.*B_R(s3,s1,m_mK2,m_mK2)*rho_denom -
       A_R(s123,s1,s3,m_mK2,m_mK2,m_mpi2)*kstar1_denom);
    const Complex f3_RR = sqrt(2.)/3.*m_FA*m_GV/sqr(m_fpi)*s123*a1_denom *
      (A_RR(s123,s3,s1,m_mK2,m_mK2,m_mpi2)*kstar3_denom +
       2.*B_RR(s123,s3,s1,s2,m_mK2,m_mK2,m_mpi2)*rho_denom -
       A_RR(s123,s1,s3,m_mK2,m_mK2,m_mpi2)*kstar1_denom);

    const Complex F3_paper = f3_chi + f3_R + f3_RR;

    if (m_isF2) return F2_paper;
    return F3_paper;
  }
  return Complex(0.,0.);
}

class FS_0_RChiPT_KKpi : public FF_0_PPP_Base {
protected:
  double  m_fpi;
  double  m_mpi, m_mK, m_meta;
  double  m_mpi2, m_mK2;
  double  m_FV, m_GV;
  double  m_mrho, m_Grho;
  double  m_mKstar, m_GKstar;
  double  m_momega, m_Gomega, m_mphi, m_Gphi;

  double  m_c125, m_c1256, m_c1235, m_c4, m_g4, m_g5, m_g2, m_g13, m_d3, m_d123;
  double  m_sin_thetaV, m_cos_thetaV;

  void    FixParameters(const FF_Parameters & params);
  void    Construct();
  Complex FF_KS(const double & s123,const double & s12,
                const double & s13,const double & s23);
  Complex FF_RChiPT(const double & s123,const double & s12,
                    const double & s13,const double & s23);
  Complex EvaluateFS(const double & s123,const double & s1,
                     const double & s2,const double & s3);
  double  GammaRho(const double & s,const double & mass,
                   const double & width) const;
  double  GammaKstar(const double & s) const;
  Complex RhoPropagator(const double & s) const;
  Complex KstarPropagator(const double & s) const;

  double  KLambda(double a,double b,double c) const;

  Complex C_R(const double & q2,const double & x,const double & m12,
              const double & m22,const double & m32) const;
  Complex C_RR(const double & q2,const double & x,const double & m2) const;
  Complex D_R(const double & q2,const double & x,const double & y) const;

public:
  FS_0_RChiPT_KKpi(const FF_Parameters & params);
  ~FS_0_RChiPT_KKpi();
};

FS_0_RChiPT_KKpi::FS_0_RChiPT_KKpi(const FF_Parameters & params) :
  FF_0_PPP_Base(params),
  m_fpi((*params.p_model)("fpi",0.1307)/sqrt(2.)),
  m_mpi((Flavour(kf_pi_plus).HadMass()+Flavour(kf_pi).HadMass())/2.),
  m_mK((Flavour(kf_K_plus).HadMass()+Flavour(kf_K).HadMass())/2.),
  m_meta(Flavour(kf_eta).HadMass()),
  m_mpi2(sqr(m_mpi)),
  m_mK2(sqr(m_mK))
{
  FixParameters(params);
  Construct();
}

FS_0_RChiPT_KKpi::~FS_0_RChiPT_KKpi() {}

void FS_0_RChiPT_KKpi::FixParameters(const FF_Parameters & params) {
  m_FV = (*params.p_model)("RChiPT_FV", sqrt(2.)*m_fpi);
  m_GV = (*params.p_model)("RChiPT_GV", sqr(m_fpi)/m_FV);

  m_mrho = (*params.p_model)("RChiPT_m_rho", 0.77554);
  m_Grho = (*params.p_model)("RChiPT_Gamma_rho", 0.1491);
  m_mKstar = (*params.p_model)("RChiPT_m_Kstar", 0.892);
  m_GKstar = (*params.p_model)("RChiPT_Gamma_Kstar", 0.050);

  m_momega = 0.78266;
  m_Gomega = 0.00849;
  m_mphi = 1.01946;
  m_Gphi = 0.00425;

  double theta_V = atan(1./sqrt(2.));
  m_sin_thetaV = sin(theta_V);
  m_cos_thetaV = cos(theta_V);

  m_c125 = 0.0;
  m_c4 = -0.07;
  m_g4 = -0.72;
  m_g5 = 0.84;
  m_d123 = 0.05;

  double root2 = sqrt(2.);
  m_g2 = m_mrho/(192.*sqr(M_PI)*root2*m_FV);
  m_g13 = -2.*m_mrho/(192.*sqr(M_PI)*root2*m_FV);
  m_c1256 = -3.*m_FV*m_mrho/(96.*sqr(M_PI)*root2*sqr(m_fpi));
  m_d3 = -sqr(m_mrho)/(64.*sqr(M_PI)*sqr(m_fpi));
  m_c1235 = 0.0;

  m_norm = (*params.p_model)("Vud", Tools::Vud)/m_fpi;
}

void FS_0_RChiPT_KKpi::Construct() {}

double FS_0_RChiPT_KKpi::KLambda(double a,double b,double c) const {
  return sqr(a-b-c) - 4.*b*c;
}

double FS_0_RChiPT_KKpi::
GammaRho(const double & s,const double & mass,const double & width) const {
  return Tools::OffShellMassWidth(s,sqr(mass),width,m_mpi2);
}

double FS_0_RChiPT_KKpi::GammaKstar(const double & s) const {
  if (s<=sqr(m_mK+m_mpi)) return 0.;
  double lam_Kpi = KLambda(s, m_mK2, m_mpi2);
  double sum = pow(lam_Kpi, 1.5);
  if (s>sqr(m_mK+m_meta)) {
    double lam_Keta = KLambda(s, m_mK2, sqr(m_meta));
    sum += pow(lam_Keta, 1.5);
  }
  sum /= (s*s);

  double lam_Kpi_peak = KLambda(sqr(m_mKstar), m_mK2, m_mpi2);
  double sum_peak = pow(lam_Kpi_peak, 1.5);
  if (sqr(m_mKstar)>sqr(m_mK+m_meta)) {
    double lam_Keta_peak = KLambda(sqr(m_mKstar), m_mK2, sqr(m_meta));
    sum_peak += pow(lam_Keta_peak, 1.5);
  }
  sum_peak /= pow(m_mKstar, 4);

  return m_GKstar * sum / sum_peak;
}

Complex FS_0_RChiPT_KKpi::RhoPropagator(const double & s) const {
  return -1./Complex(s-sqr(m_mrho),-m_mrho*GammaRho(s,m_mrho,m_Grho));
}

Complex FS_0_RChiPT_KKpi::KstarPropagator(const double & s) const {
  return -1./Complex(s-sqr(m_mKstar),-m_mKstar*GammaKstar(s));
}

Complex FS_0_RChiPT_KKpi::
C_R(const double & q2,const double & x,const double & m12,
    const double & m22,const double & m32) const {
  return m_c125*q2 - m_c1256*x + m_c1235*m32 + 8.*m_c4*(m12-m22);
}

Complex FS_0_RChiPT_KKpi::C_RR(const double & q2,const double & x,
                               const double & m2) const {
  return m_d3*(q2+x) + m_d123*m2;
}

Complex FS_0_RChiPT_KKpi::
D_R(const double & q2,const double & x,const double & y) const {
  return (m_g13+2.*m_g2)*(x+y) - 2.*m_g2*(q2+m_mK2) - m_g13*(3.*m_mK2+m_mpi2) +
    2.*m_g4*(m_mK2+m_mpi2) + 2.*m_g5*m_mK2;
}

Complex FS_0_RChiPT_KKpi::
FF_KS(const double & s123,const double & s12,const double & s13,
      const double & s23) {
  return Complex(0.,0.);
}

Complex FS_0_RChiPT_KKpi::
EvaluateFS(const double & s123,const double & s1,const double & s2,
           const double & s3) {
  const Complex omega_denom = 1./Complex(sqr(m_momega)-s2,-m_momega*m_Gomega);
  const Complex phi_denom = 1./Complex(sqr(m_mphi)-s2,-m_mphi*m_Gphi);

  const Complex omega_phi_mix = sqr(m_sin_thetaV) *
    (1. + sqrt(2.)/m_sin_thetaV*m_cos_thetaV)*omega_denom +
    sqr(m_cos_thetaV)*(1. - sqrt(2.)/m_cos_thetaV*m_sin_thetaV)*phi_denom;

  const Complex f5_chi = sqrt(2.);
  const Complex kstar_denom = -KstarPropagator(s1);
  const Complex rho_denom = -RhoPropagator(s123);

  const Complex f5_R = 16.*sqr(M_PI)*m_GV/m_mrho *
    (C_R(s123,s2,m_mK2,m_mK2,m_mpi2)*omega_phi_mix +
     C_R(s123,s1,m_mK2,m_mpi2,m_mK2)*kstar_denom -
     2.*m_FV/m_GV*D_R(s123,s2,s1)*rho_denom);

  const Complex f5_RR = -16.*sqrt(2.)*sqr(M_PI)*m_FV*m_GV*rho_denom *
    (C_RR(s123,s1,m_mK2)*kstar_denom + C_RR(s123,s2,m_mpi2)*omega_phi_mix);

  const Complex F5_paper = f5_chi + f5_R + f5_RR;

  return Complex(0.,-1./(4.*sqr(M_PI)*sqr(m_fpi))) * F5_paper;
}

Complex FS_0_RChiPT_KKpi::
FF_RChiPT(const double & s123,const double & s12,const double & s13,
          const double & s23) {
  const double s1(s23), s2(s13), s3(s12);
  if (m_mode==FF_0_PPP_mode::K_pi_K || m_mode==FF_0_PPP_mode::K0_pi_K0b) {
    return EvaluateFS(s123,s1,s2,s3);
  }
  else if (m_mode==FF_0_PPP_mode::KS_pi_KS) {
    Complex fs_orig = EvaluateFS(s123,s1,s2,s3);
    Complex fs_swap = EvaluateFS(s123,s3,s2,s1);
    return 0.5*(fs_orig-fs_swap);
  }
  else if (m_mode==FF_0_PPP_mode::KS_pi_KL) {
    Complex fs_orig = EvaluateFS(s123,s1,s2,s3);
    Complex fs_swap = EvaluateFS(s123,s3,s2,s1);
    return 0.5*(fs_orig+fs_swap);
  }
  else if (m_mode==FF_0_PPP_mode::KL_pi_KL) {
    Complex fs_orig = EvaluateFS(s123,s1,s2,s3);
    Complex fs_swap = EvaluateFS(s123,s3,s2,s1);
    return -0.5*(fs_orig+fs_swap);
  }
  else if (m_mode==FF_0_PPP_mode::K_pi0_K0) {
    const Complex f5_chi = 0.;
    const Complex kstar3_denom = -KstarPropagator(s3);
    const Complex kstar1_denom = -KstarPropagator(s1);

    const Complex f5_R = 8.*sqrt(2.)*sqr(M_PI)*m_GV/m_mrho *
      (C_R(s123,s3,m_mK2,m_mpi2,m_mK2)*kstar3_denom -
       C_R(s123,s1,m_mK2,m_mpi2,m_mK2)*kstar1_denom);

    const Complex f5_RR = 16.*sqr(M_PI)*m_FV*m_GV/(s123-sqr(m_mrho)) *
      (-C_RR(s123,s3,m_mK2)*kstar3_denom + C_RR(s123,s1,m_mK2)*kstar1_denom);

    const Complex F5_paper = f5_chi + f5_R + f5_RR;

    return Complex(0.,1./(4.*sqr(M_PI)*sqr(m_fpi))) * F5_paper;
  }
  return Complex(0.,0.);
}


DECLARE_FF_GETTER(FF_0_PPP_Base,"FF_0_PPP")

FormFactor_Base * ATOOLS::Getter<FormFactor_Base,FF_Parameters,
				 FF_0_PPP_Base>::
operator()(const METOOLS::FF_Parameters &params) const
{
  size_t Nmesons = 0;
  for (size_t i=0;i<params.m_pi.size();i++) {
    if (params.m_flavs[params.m_pi[i]].IsMeson()) Nmesons++;
  }

  if (Nmesons!=3) return NULL;

  const FF_0_PPP_mode mode = ModeFromFlavours(params.m_flavs,params.m_pi);

  if (mode==FF_0_PPP_mode::piP_pi0_pi0 ||
      mode==FF_0_PPP_mode::piM_piP_piP) {
    if (params.m_name=="F1_0_PPP") return new F1_0_PiPlusPiZeroPiZero(params);
    if (params.m_name=="F2_0_PPP") return new F1_0_PiPlusPiZeroPiZero(params);
    if (params.m_name=="F3_0_PPP") return new F3_0_PiPlusPiZeroPiZero(params);
    if (params.m_name=="FS_0_PPP") return new FS_0_PiPlusPiZeroPiZero(params);
  }
  else if (mode!=FF_0_PPP_mode::unknown) {
    if (params.m_ffmodel==ff_model::RChiPT &&
        (mode==FF_0_PPP_mode::K_pi_K || mode==FF_0_PPP_mode::K0_pi_K0b ||
         mode==FF_0_PPP_mode::K_pi0_K0 || mode==FF_0_PPP_mode::KS_pi_KS ||
         mode==FF_0_PPP_mode::KS_pi_KL || mode==FF_0_PPP_mode::KL_pi_KL)) {
      if (params.m_name=="F1_0_PPP") return new F1_0_RChiPT_KKpi(params);
      if (params.m_name=="F2_0_PPP") return new F1_0_RChiPT_KKpi(params);
      if (params.m_name=="F3_0_PPP") return new F3_0_PiPlusPiZeroPiZero(params);
      if (params.m_name=="FS_0_PPP") return new FS_0_RChiPT_KKpi(params);
    }
    else {
      if (params.m_name=="F1_0_PPP") return new F1_0_FM95(params);
      if (params.m_name=="F2_0_PPP") return new F1_0_FM95(params);
      if (params.m_name=="F3_0_PPP") return new F3_0_PiPlusPiZeroPiZero(params);
      if (params.m_name=="FS_0_PPP") return new FS_0_FM95(params);
    }
  }

  return NULL;
}
