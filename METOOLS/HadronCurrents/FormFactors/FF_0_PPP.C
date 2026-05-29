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

    if (flavs[indices[0]].Kfcode()==kf_pi &&
        flavs[indices[1]].Kfcode()==kf_pi &&
        flavs[indices[2]].Kfcode()==kf_pi_plus) {
      return FF_0_PPP_mode::piP_pi0_pi0;
    }
    if (flavs[indices[0]].Kfcode()==kf_pi_plus &&
        flavs[indices[1]].Kfcode()==kf_pi_plus &&
        flavs[indices[2]].Kfcode()==kf_pi_plus) {
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

    // case ff_model::RChiPT:
    //   return m_norm * FF_RChiPT(s123,s12,s13,s23);

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
// {pi pi pi}^+, {K pi pi}^+, {K K pi}^+  form factors.
// Todos:
//  - add more theory - i.e. Gounaris-Sakurai forms and variations of RChiPT
//  - mirror parameters to run card/decay yaml files
//
// Form factors from:
// - KS
//     * pi pi pi from
// - none, 0: no form factor
//
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

class F1_0_PiPlusPiZeroPiZero : public FF_0_PPP_Base {
protected:
  bool    m_isF2;
  double  m_fpi;
  Complex m_alpha, m_gamma, m_delta;

  Total_Width_Base  * p_a1_ks_width;
  Summed_Propagator * p_a1s, * p_rhos;

  void    FixParameters(const FF_Parameters & params);
  void    Construct();
  Complex FF_KS(const double & s123,const double & s12,
                const double & s13,const double & s23);
  Complex FF_RChiPT(const double & s123,const double & s12,
                    const double & s13,const double & s23);

  /*
    Complex A(const double & m2,const double & s,const double & mu2);
    double  Gamma_V(const double & s);
    double  Gamma_Vp(const double & s);
    double  Gamma_Vpp(const double & s);
  */

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
  m_fpi((*params.p_model)("fpi",0.1307)/sqrt(2.))
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
    m_norm = Complex(0., -((2.*sqrt(2)*(*params.p_model)("Vud", Tools::Vud)) /
			   (3.*m_fpi)));

    if (m_ffmodel==ff_model::KS) {
      m_gamma = Complex(-0.14500,0.0000);
      m_delta = Complex( 0.00000,0.0000);
      m_alpha = Complex( 0.00185,0.0000);
    }
  }
}

void F1_0_PiPlusPiZeroPiZero::Construct() {
  if (m_ffmodel==ff_model::KS) {
    p_a1_ks_width = new KS_A1_1260_plus_Width();
    if (m_mode==FF_0_PPP_mode::piP_pi0_pi0) {
      Propagator_Base * a11260 =
        new BreitWigner(p_a1_ks_width);
      Propagator_Base * rho770 =
        new BreitWigner(LineShapes->Get(Flavour(kf_rho_770_plus)));
      Propagator_Base * rho1450 =
        new BreitWigner(LineShapes->Get(Flavour(kf_rho_1450_plus)));

      p_a1s = new Summed_Propagator();
      p_a1s->Add(a11260, Complex(1.,0.));

      p_rhos = new Summed_Propagator();
      p_rhos->Add(rho770,  Complex(1.,0.));
      p_rhos->Add(rho1450, m_gamma);
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
        new BreitWigner(LineShapes->Get(Flavour(kf_rho_1450_plus)));

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
    }
  }
}

Complex F1_0_PiPlusPiZeroPiZero::
FF_KS(const double & s123,const double & s12,const double & s13,
      const double & s23) {
  if (p_a1s==NULL || p_rhos==NULL) return Complex(0.,0.);

  if (m_isF2) {
    return m_norm * (*p_a1s)(s123) * (*p_rhos)(s12);
  }

  return m_norm * (*p_a1s)(s123) * (*p_rhos)(s13);
}

Complex F1_0_PiPlusPiZeroPiZero::
FF_RChiPT(const double & s123,const double & s12,const double & s13,
          const double & s23) {
  return Complex(0.,0.);
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

  m_mK1_1400 = (*params.p_model)("FM95_m_K1_1400",     1.402);
  m_GK1_1400 = (*params.p_model)("FM95_Gamma_K1_1400", 0.174);
  m_mK1_1270 = (*params.p_model)("FM95_m_K1_1270",     1.270);
  m_GK1_1270 = (*params.p_model)("FM95_Gamma_K1_1270", 0.090);
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
  F3_0_PiPlusPiZeroPiZero(const FF_Parameters & params) :
    FF_0_PPP_Base(params),
    m_fpi((*params.p_model)("fpi",0.1307)/sqrt(2.))
  {
    FixParameters(params);
    Construct();
  }
};

void F3_0_PiPlusPiZeroPiZero::FixParameters(const FF_Parameters & params) {}

void F3_0_PiPlusPiZeroPiZero::Construct() {}


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
    if (params.m_name=="F1_0_PPP") return new F1_0_FM95(params);
    if (params.m_name=="F2_0_PPP") return new F1_0_FM95(params);
    if (params.m_name=="F3_0_PPP") return new F3_0_PiPlusPiZeroPiZero(params);
    if (params.m_name=="FS_0_PPP") return new FS_0_FM95(params);
  }

  return NULL;
}
