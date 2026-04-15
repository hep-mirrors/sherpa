#include "METOOLS/HadronCurrents/FormFactors/FF_0_PP.H"
#include "METOOLS/HadronCurrents/Tools.H"
#include "ATOOLS/Org/Message.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;


void FF_0_PP_Base::FixMode() {
  if (m_flavs[m_pi[0]].Kfcode()==kf_pi &&
      m_flavs[m_pi[1]].Kfcode()==kf_pi_plus)     m_mode = FF_0_PP_mode::pipi_plus;
  else if (m_flavs[m_pi[0]].Kfcode()==kf_K &&
	    m_flavs[m_pi[1]].Kfcode()==kf_K_plus) m_mode = FF_0_PP_mode::KK_plus;
  else if (m_flavs[m_pi[0]].Kfcode()==kf_pi_plus &&
	    (m_flavs[m_pi[1]].Kfcode()==kf_K   ||
	     m_flavs[m_pi[1]].Kfcode()==kf_K_S ||
	     m_flavs[m_pi[1]].Kfcode()==kf_K_L))     m_mode = FF_0_PP_mode::Kpi_plus;
  else if (m_flavs[m_pi[0]].Kfcode()==kf_pi &&
	    m_flavs[m_pi[1]].Kfcode()==kf_K_plus)    m_mode = FF_0_PP_mode::Kpi_plus;
//  msg_Out() << "FF_0_PP_Base::FixMode selected mode = " << int(m_mode) << std::endl;
}

Complex FF_0_PP_Base::operator()(const ATOOLS::Vec4D_Vector& moms) {
  double Q2 = (moms[m_pi[0]]+moms[m_pi[1]]).Abs2();
  switch (m_ffmodel) {
  case ff_model::none:    return Complex(1.,0.);
  case ff_model::KS:      return ( p_props!=NULL ?
				   m_norm * (*p_props)(Q2) : Complex(0.,0.) );
  case ff_model::RChiPT:  return m_norm * FF_RChiPT(Q2);
  case ff_model::unknown:
  default:
    break;
  }
  return Complex(0.,0.);
}

Complex FF_0_PP_Base::FF_RChiPT(const double & Q2) {
  msg_Error()<<"Error in "<<METHOD<<": RChiPT not available for "
	     <<m_flavs[m_pi[0]]<<"+"<<m_flavs[m_pi[1]]<<" form factor.\n"
	     <<"   Will exit the run.\n";
  return Complex(1.,0.);
}


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//
// pi^0 pi^+ form factors.
// Todos:
//  - add more theory - i.e. Gounaris-Sakurai forms and variations of RChiPT
//  - mirror parameters to run card/decay yaml files
//
// Form factors from:
// - KS, 1: (Kuehn-Santamaria model):
//   * pi pi (original version): Z.Phys.C 48 (1990) 445-452
//     (https://doi.org/10.1007/BF01572024)
//   * pi pi (TODO: KS with Gounaris-Sakurai form factors):
//     from Belle measurement & analysis:
//     (https://arxiv.org/pdf/0805.3773) 
//   * K pi: Z.Phys.C 69 (1996) 243 (vector form factor)
//     (https://doi.org/10.1007/s002880050024)
// - RChT, 2: Resonance Chiral Perturbation Theory
//   * pi pi/KK fit: Eur.Phys.J.C 79 (2019) 5, 436
//     (https://doi.org/10.1140/epjc/s10052-019-6943-9)
//           - We could also have some alternative dispersive approach, 
//             maybe as form factor model 3, based on RchT
//             (https://arxiv.org/pdf/1112.0962)
//
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

class Fplus_0_PiZeroPiPlus : public FF_0_PP_Base {
  double  m_fpi, m_fK, m_mpi2, m_mK2, m_meta2;
  double  m_mV, m_mV2, m_mVp, m_mVp2, m_mVpp, m_mVpp2, m_GVp, m_GVpp;
  Complex m_gamma, m_delta;
  
  void    FixParameters(const FF_Parameters & params);
  void    Construct();
  Complex FF_RChiPT(const double & s);
  Complex A(const double & m2,const double & s,const double & mu2);
  double  Gamma_V(const double & s);
  double  Gamma_Vp(const double & s);
  double  Gamma_Vpp(const double & s);
  Complex Jbar_PQ(const double & mP2, const double & mQ2, const double & s);
  double  Jbar_prime0(const double & mP2, const double & mQ2);
  Complex H_tilde(const double & s, const double & mP2, const double & mQ2, const double & mu2);
  double Gamma_Kstar(const double & s);
  double Gamma_Kstarp(const double & s);
public :
  Fplus_0_PiZeroPiPlus(const FF_Parameters & params);
};

Fplus_0_PiZeroPiPlus::Fplus_0_PiZeroPiPlus(const FF_Parameters & params)  :
  FF_0_PP_Base(params),
  m_fpi((*params.p_model)("fpi",0.1307)/sqrt(2.)),
  m_fK((*params.p_model)("fK",  0.1562)/sqrt(2.)),
  m_mpi2(sqr(Flavour(kf_pi_plus).HadMass())),
  m_mK2(sqr(Flavour(kf_K_plus).HadMass())),
  m_meta2(sqr(Flavour(kf_eta).HadMass()))
{
  FixParameters(params);
  Construct();
}

//Fplus_0_PiZeroPiPlus::Fplus_0_PiZeroPiPlus(const FF_Parameters & params)  :
//  FF_0_PP_Base(params),
//  m_fpi((*params.p_model)("fpi",0.1307)/sqrt(2.)),
//  m_mpi2(sqr(Flavour(kf_pi_plus).HadMass())),
//  m_mK2(sqr(Flavour(kf_K_plus).HadMass()))
//{
//  msg_Out() << "Fplus ctor: start" << std::endl;
//  FixParameters(params);
//  msg_Out() << "Fplus ctor: after FixParameters" << std::endl;
//  Construct();
//  msg_Out() << "Fplus ctor: after Construct" << std::endl;
//}

void Fplus_0_PiZeroPiPlus::FixParameters(const FF_Parameters & params)  {
  if (m_mode==FF_0_PP_mode::pipi_plus) {
    m_norm    = (*params.p_model)("Vud", Tools::Vud);
    if (m_ffmodel==ff_model::KS) {
      m_gamma = Complex(-0.167,0.000);
      m_delta = Complex(0.050,0.000);
    }
    else if (m_ffmodel==ff_model::RChiPT) {
      m_mV    = Flavour(kf_rho_770_plus).HadMass();  m_mV2   = sqr(m_mV);
      m_mVp   = 1.438; m_mVp2  = sqr(m_mVp);  m_GVp  = 0.535;
      m_mVpp  = 1.754; m_mVpp2 = sqr(m_mVpp); m_GVpp = 0.412;
      m_gamma = Complex(-0.15*cos(-0.36),-0.15*sin(-0.36));
      m_delta = Complex( 0.12*cos(-0.02), 0.12*sin(-0.02));
      //m_mVp   = Flavour(kf_rho_1450_plus).HadMass(); m_mVp2  = sqr(m_mVp);
      //m_mVpp  = Flavour(kf_rho_1700_plus).HadMass(); m_mVpp2 = sqr(m_mVpp); 
      //m_GVp   = Flavour(kf_rho_1450_plus).Width();
      //m_GVpp  = Flavour(kf_rho_1700_plus).Width();
    }
  }
  else if (m_mode==FF_0_PP_mode::KK_plus) {
    m_norm    = (*params.p_model)("Vud", Tools::Vud)/sqrt(2.);
    if (m_ffmodel==ff_model::KS) {
      m_gamma = Complex(-0.167,0.000);
      m_delta = Complex(0.050,0.000);
    }
    else if (m_ffmodel==ff_model::RChiPT) {
      m_mV    = Flavour(kf_rho_770_plus).HadMass();  m_mV2   = sqr(m_mV);
      m_mVp   = 1.467; m_mVp2  = sqr(m_mVp);  m_GVp  = 0.415;
      m_mVpp  = 1.754; m_mVpp2 = sqr(m_mVpp); m_GVpp = 0.412;
      m_gamma = Complex( 0.09*cos(-1.88),-0.15*sin(-1.88));
      m_delta = Complex( 0.00*cos( 0.00), 0.00*sin( 0.00));
    }
  }
  else if (m_mode==FF_0_PP_mode::Kpi_plus) {
    m_norm    = (*params.p_model)("Vus", Tools::Vus)/sqrt(2.);
    if (m_ffmodel==ff_model::KS) {
      m_gamma = Complex(-0.135,0.000);
    }
    else if (m_ffmodel==ff_model::RChiPT) {
      m_mV    = Flavour(kf_K_star_892_plus).HadMass();  m_mV2   = sqr(m_mV);
      m_mVp   = Flavour(kf_K_star_1410_plus).HadMass(); m_mVp2  = sqr(m_mVp);
      m_GVp   = Flavour(kf_K_star_1410_plus).Width();
      m_gamma = Complex((*params.p_model)("gamma_Kpi", 0.0), 0.0);
    }
  }
}

void Fplus_0_PiZeroPiPlus::Construct() {
  if (m_ffmodel==ff_model::KS) {
    if (m_mode==FF_0_PP_mode::pipi_plus ||
	m_mode==FF_0_PP_mode::KK_plus) {
      Propagator_Base * rho770  =
	new BreitWigner(LineShapes->Get(Flavour(kf_rho_770_plus)));
      Propagator_Base * rho1450 =
	new BreitWigner(LineShapes->Get(Flavour(kf_rho_1450_plus)));
      Propagator_Base * rho1700 =
	new BreitWigner(LineShapes->Get(Flavour(kf_rho_1700_plus)));
      p_props = new Summed_Propagator();
      p_props->Add(rho770,  Complex(1., 0.));
      p_props->Add(rho1450, m_gamma);
      p_props->Add(rho1700, m_delta);
    }
    else if (m_mode==FF_0_PP_mode::Kpi_plus) {
      Propagator_Base * Kstar892  =
	new BreitWigner(LineShapes->Get(Flavour(kf_K_star_892_plus)));
      Propagator_Base * Kstar1410 =
	new BreitWigner(LineShapes->Get(Flavour(kf_K_star_1410_plus)));
      p_props = new Summed_Propagator();
      p_props->Add(Kstar892,  Complex(1., 0.));
      p_props->Add(Kstar1410, m_gamma);
    }
  }
}

//void Fplus_0_PiZeroPiPlus::Construct() {
////  msg_Out() << "Construct(): m_ffmodel=" << int(m_ffmodel)
////            << " m_mode=" << int(m_mode) << std::endl;
//
//  if (m_ffmodel == ff_model::KS &&
//      m_mode == FF_0_PP_mode::Kpi_plus) {
//
////    msg_Out() << "A" << std::endl;
//
//    auto ls892 = LineShapes->Get(Flavour(kf_K_star_892_plus));
////    msg_Out() << "B ls892=" << ls892 << std::endl;
//
//    Propagator_Base * Kstar892 = new BreitWigner(ls892);
////    msg_Out() << "C Kstar892=" << Kstar892 << std::endl;
//
//    p_props = new Summed_Propagator();
////    msg_Out() << "D p_props=" << p_props << std::endl;
//
//    p_props->Add(Kstar892, Complex(1.,0.));
////    msg_Out() << "E added 892" << std::endl;
//  }
//}

Complex Fplus_0_PiZeroPiPlus::FF_RChiPT(const double & s) {
  //msg_Out() << "FF_RChiPT entered with s = " << s << std::endl;
  if (m_mode == FF_0_PP_mode::pipi_plus || m_mode == FF_0_PP_mode::KK_plus) {
    Complex V   = ( (m_mV2 - s*(m_gamma+m_delta)) *
                    Complex(m_mV2-s, m_mV*Gamma_V(s)) /
                    (sqr(m_mV2-s)+sqr(m_mV*Gamma_V(s))) *
                    exp(-s/(96*sqr(M_PI*m_fpi)) *
                        (A(m_mpi2,s,m_mV2)+1./2.*A(m_mK2,s,m_mV2)).real()) );
    Complex Vp  = ( s*m_gamma *
                    Complex(m_mVp2-s, m_mVp*Gamma_Vp(s)) /
                    (sqr(m_mVp2-s)+sqr(m_mVp*Gamma_Vp(s))) *
                    exp(-s*Gamma_Vp(m_mVp2) /
                        (M_PI*pow(m_mVp*sqrt(1.-4.*m_mpi2/m_mVp2),3./2.)) *
                        A(m_mpi2,s,m_mV2).real()) );
    Complex Vpp = ( s*m_delta *
                    Complex(m_mVpp2-s, m_mVpp*Gamma_Vpp(s)) /
                    (sqr(m_mVpp2-s)+sqr(m_mVpp*Gamma_Vpp(s))) *
                    exp(-s*Gamma_Vpp(m_mVpp2) /
                        (M_PI*pow(m_mVpp*sqrt(1.-4.*m_mpi2/m_mVpp2),3./2.)) *
                        A(m_mpi2,s,m_mV2).real()) );
    return V + Vp + Vpp;
  }
  if (m_mode == FF_0_PP_mode::Kpi_plus) {
    Complex V  = (m_mV2  + s*m_gamma) /
                 Complex(m_mV2  - s, -m_mV  * Gamma_Kstar(s));
    Complex Vp = (s*m_gamma) /
                 Complex(m_mVp2 - s, -m_mVp * Gamma_Kstarp(s));
    Complex exp_factor = exp(1.5 * (H_tilde(s, m_mK2, m_mpi2, m_mV2) +
                                H_tilde(s, m_mK2, m_meta2, m_mV2)).real());
    return (V - Vp) * exp_factor;
  }
  return Complex(0.0,0.0);
}

Complex Fplus_0_PiZeroPiPlus::
A(const double & m2,const double & s,const double & mu2) {
  Complex sigma = csqrt(1.-4.*m2/s);
  return log(m2/mu2)+8.*m2/s-5./3.+pow(sigma,3.)*log((sigma+1.)/(sigma-1.));
}

double Fplus_0_PiZeroPiPlus::Gamma_V(const double & s) {
  if (s<4.*m_mpi2) return 0.;
  double pref = m_mV*s/(96.*M_PI*sqr(m_fpi));
  double arg  = (pow(1.-4.*m_mpi2/s,3./2.) +
		 ((s>4.*m_mK2) ? 1./2. * pow(1.-4.*m_mK2/s,3./2.) : 0.));
  return pref*arg;
}

double Fplus_0_PiZeroPiPlus::Gamma_Vp(const double & s) {
  if (s<4.*m_mpi2) return 0.;
  return (m_GVp * s/m_mVp2 *
	  pow((s-4.*m_mpi2)/(m_mVp2-4.*m_mpi2),3./2.));
}

double Fplus_0_PiZeroPiPlus::Gamma_Vpp(const double & s) {
  if (s<4.*m_mpi2) return 0.;
  return (m_GVpp * s/m_mVpp2 *
	  pow((s-4.*m_mpi2)/(m_mVpp2-4.*m_mpi2),3./2.));
}
namespace {
  inline double Lambda(const double x, const double y, const double z) {
    return x*x + y*y + z*z - 2.0*(x*y + x*z + y*z);
  }

  inline double MomentumPQ(const double s, const double mP2, const double mQ2) {
    if (s <= 0.0) return 0.0;
    const double lam = Lambda(s, mP2, mQ2);
    if (lam <= 0.0) return 0.0;
    return sqrt(lam) / (2.0 * sqrt(s));
  }
}

double Fplus_0_PiZeroPiPlus::Gamma_Kstar(const double & s) {
  const double sth = sqr(sqrt(m_mK2) + sqrt(m_mpi2));
  if (s <= sth) return 0.0;

  const double p_s   = MomentumPQ(s,     m_mK2, m_mpi2);
  const double p_mV2 = MomentumPQ(m_mV2, m_mK2, m_mpi2);
  if (p_mV2 <= 0.0) return 0.0;

  const double G0 = Flavour(kf_K_star_892_plus).Width();

  return G0 * (m_mV2 / s) * pow(p_s / p_mV2, 3);
}

double Fplus_0_PiZeroPiPlus::Gamma_Kstarp(const double & s) {
  const double sth = sqr(sqrt(m_mK2) + sqrt(m_mpi2));
  if (s <= sth) return 0.0;

  const double p_s    = MomentumPQ(s,      m_mK2, m_mpi2);
  const double p_mVp2 = MomentumPQ(m_mVp2, m_mK2, m_mpi2);
  if (p_mVp2 <= 0.0) return 0.0;

  return m_GVp * (m_mVp2 / s) * pow(p_s / p_mVp2, 3);
}

Complex Fplus_0_PiZeroPiPlus::Jbar_PQ(const double & mP2,
                                      const double & mQ2,
                                      const double & s)
{
  const double Delta = mP2 - mQ2;
  if (std::abs(Delta) < 1e-10) return Complex(0., 0.);
  const Complex nu = std::sqrt(
    Complex((s - sqr(std::sqrt(mP2) + std::sqrt(mQ2))) *
            (s - sqr(std::sqrt(mP2) - std::sqrt(mQ2))), 0.0)
            );

  const Complex logterm =
    std::log((s + nu + Delta)*(s + nu - Delta) / ((s - nu + Delta)*(s - nu - Delta)));

  return (1.0 / (32.0 * sqr(M_PI))) *
          ( 2.0
           + (Delta / s - (mP2 + mQ2) / Delta) * log(mQ2 / mP2)
           - (nu / s) * logterm );
}

double Fplus_0_PiZeroPiPlus::Jbar_prime0(const double & mP2, const double & mQ2)
{
  // GL (A.10): Jbar'(0) for unequal masses
  const double Delta  = mP2 - mQ2;
  const double Sigma  = mP2 + mQ2;
  if (std::abs(Delta) < 1e-10) {
    // equal mass limit from GL (A.11): Jbar'(0) = 1/(96*pi^2*M^2)
    return 1./(96.*sqr(M_PI)*mP2);
  }
  return (1./(32.*sqr(M_PI))) *
         (Sigma/sqr(Delta) + 2.*mP2*mQ2/pow(Delta,3.) * log(mQ2/mP2));
}

Complex Fplus_0_PiZeroPiPlus::H_tilde(const double & s,
                                      const double & mP2,
                                      const double & mQ2,
                                      const double & mu2)
{
  const double Delta   = mP2 - mQ2;
  const double Sigma   = mP2 + mQ2;
  const Complex Jbar   = Jbar_PQ(mP2, mQ2, s);
  const Complex Jdbar  = Jbar - s*Jbar_prime0(mP2, mQ2);  // GL (8.10)

  // scale k at mu = M_Kstar, GL (8.9)
  const double k   = (std::abs(Delta) < 1e-10) ?
      (1./(32.*sqr(M_PI))) * (log(mP2/mu2) + 1.) :
      (1./(32.*sqr(M_PI))) * (mP2*log(mP2/mu2) - mQ2*log(mQ2/mu2)) / Delta;

  // sM^r - L using GL (8.9)
  const Complex sMr_minus_L =
      ((s - 2.*Sigma)/12.) * Jbar
      + (sqr(Delta)/(3.*s)) * Jdbar
      - (sqr(Delta)/(4.*s)) * Jbar
      - s*k/6.
      + s/(288.*sqr(M_PI));

  return sMr_minus_L / (m_fK * m_fpi);
}


//////////////////////////////////////////////////////////////////////////////////
//
// Form factors for
//   * pi pi & K K vanish due to identical masses
// - Kuehn-Santamaria
//   * K pi: Z.Phys.C 72 (1996) 619 (scalar form factor)
//     (https://doi.org/10.1007/s002880050284), equation references below refer to this paper
//
// Todo: implement all scalar form factors here!
//
//////////////////////////////////////////////////////////////////////////////////

class Fzero_0_PiZeroPiPlus : public FF_0_PP_Base {
  double  m_mpi2, m_mK2;
  double  m_mV2, m_mVp2, m_mS2;
  Complex m_beta, m_cS;
  Propagator_Base * p_BW_V;    // K*(892)    -- vector off-shell contribution
  Propagator_Base * p_BW_Vp;   // K*(1410)   -- vector off-shell contribution
  Propagator_Base * p_BW_S;    // K*_0(1430) -- direct scalar contribution

  void    FixParameters(const FF_Parameters & params);
  void    Construct();
  Complex FF_RChiPT(const double & Q2) { return Complex(0.,0.); }
  Complex operator()(const ATOOLS::Vec4D_Vector& moms);
public :
  Fzero_0_PiZeroPiPlus(const FF_Parameters & params) :
    FF_0_PP_Base(params),
    m_mpi2(sqr(Flavour(kf_pi_plus).HadMass())),
    m_mK2(sqr(Flavour(kf_K_plus).HadMass())),
    p_BW_V(NULL), p_BW_Vp(NULL), p_BW_S(NULL)
  {
    FixParameters(params);
    Construct();
  }
  ~Fzero_0_PiZeroPiPlus() {
    if (p_BW_V)  delete p_BW_V;
    if (p_BW_Vp) delete p_BW_Vp;
    if (p_BW_S)  delete p_BW_S;
  }
};

void Fzero_0_PiZeroPiPlus::FixParameters(const FF_Parameters & params) {
  if (m_mode==FF_0_PP_mode::Kpi_plus) {
    m_norm  = (*params.p_model)("Vus", Tools::Vus)/sqrt(2.);
    if (m_ffmodel==ff_model::KS) {
      m_mV2   = sqr(Flavour(kf_K_star_892_plus).HadMass());
      m_mVp2  = sqr(Flavour(kf_K_star_1410_plus).HadMass());
      m_mS2   = sqr(Flavour(kf_K_0_star_1430_plus).HadMass());
      // beta_{K*}: relative K*(1410) contribution (same as in the vector FF)
      m_beta  = Complex((*params.p_model)("beta_Kpi", -0.135), 0.);
      // c_S: scalar resonance strength, fixed by matching to ChPT (eq. 27)
      m_cS    = Complex((*params.p_model)("cS_Kpi",    1.7  ), 0.);
    }
  }
}

void Fzero_0_PiZeroPiPlus::Construct() {
  if (m_ffmodel==ff_model::KS) {
    if (m_mode==FF_0_PP_mode::Kpi_plus) {
      p_BW_V  = new BreitWigner(LineShapes->Get(Flavour(kf_K_star_892_plus)));
      p_BW_Vp = new BreitWigner(LineShapes->Get(Flavour(kf_K_star_1410_plus)));
      p_BW_S  = new BreitWigner(LineShapes->Get(Flavour(kf_K_0_star_1430_plus)));
    }
  }
}

Complex Fzero_0_PiZeroPiPlus::operator()(const Vec4D_Vector& moms) {
  // Implements the scalar Kpi form factor from Finkemeier & Mirkes,
  // Z.Phys.C 72 (1996) 619, eq. (18).
  if (m_ffmodel != ff_model::KS || m_mode != FF_0_PP_mode::Kpi_plus)
    return Complex(0., 0.);
  double Q2 = (moms[m_pi[0]]+moms[m_pi[1]]).Abs2();
  if (Q2 <= 0.) return Complex(0., 0.);

  double  deltam2 = m_mK2 - m_mpi2;
  Complex beta    = m_beta;   // beta_{K*} = -0.135
  Complex cS      = m_cS;     // c_S = 1.7

  // Vector resonance off-shell scalar projection (eq. 18, lines 1-2):
  // (1/Q^2) * 1/(1+beta) * [ (mV^2-Q^2)/mV^2 * BW_V + beta*(mVp^2-Q^2)/mVp^2 * BW_Vp ]
  Complex BW_V     = (*p_BW_V)(Q2);
  Complex BW_Vp    = (*p_BW_Vp)(Q2);
  Complex vec_part = ( (m_mV2  - Q2)/m_mV2  * BW_V  +
		       beta * (m_mVp2 - Q2)/m_mVp2 * BW_Vp )
		     / (Q2 * (1. + beta));

  // Direct scalar resonance contribution (eq. 18, last line):
  // c_S / m_{K*_0}^2 * BW_{K*_0}
  Complex BW_S        = (*p_BW_S)(Q2);
  Complex scalar_part = cS / m_mS2 * BW_S;

  Complex result = m_norm * deltam2 * (vec_part + scalar_part);
  //msg_Out() << "DEBUG Fzero_0_PiZeroPiPlus: Q2=" << Q2
  //	    << "  vec=" << vec_part << "  scalar=" << scalar_part
  //	    << "  F0=" << result << "\n";
  return result;
}


DECLARE_FF_GETTER(FF_0_PP_Base,"FF_0_PP")

FormFactor_Base * ATOOLS::Getter<FormFactor_Base,FF_Parameters,
				 FF_0_PP_Base>:: 
operator()(const METOOLS::FF_Parameters &params) const
{
  msg_Out()<<METHOD<<":\n";
  size_t Nmesons = 0;
  for (size_t i=0;i<params.m_pi.size();i++) {
    if (params.m_flavs[params.m_pi[i]].IsMeson()) Nmesons++;
  }
  if (Nmesons!=2) return NULL;
  // Below a first round of decays/currents for which we have both
  // Kuehn-Santamaria and RChiPT parametrizations
  if (//   pi^+ pi^0
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_pi &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_pi_plus)    ||
      //   K^+ K^0 (note - this could be K_S or K_L - just to play it safe
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_K &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_K_plus)     ||
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_K_S &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_K_plus)     ||
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_K_L &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_K_plus)     ||
       //   pi^+ K^0 (note - this could be K_S or K_L - just to play it safe
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_K &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_pi_plus)    ||
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_pi_plus &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_K)          ||
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_K_S &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_pi_plus)    ||
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_pi_plus &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_K_S)        ||
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_K_L &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_pi_plus)    ||
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_pi_plus &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_K_L)        ||
       //   K^+ pi^0
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_K_plus &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_pi)         ||
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_pi &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_K_plus)) {
      if (params.m_name=="F+_0_PP") {
//        msg_Out() << "Getter: about to build Fplus" << std::endl;
        return new Fplus_0_PiZeroPiPlus(params);
      }
      if (params.m_name=="F0_0_PP") {
//        msg_Out() << "Getter: about to build Fzero" << std::endl;
        return new Fzero_0_PiZeroPiPlus(params);
      }
    //if (params.m_name=="F+_0_PP") return new Fplus_0_PiZeroPiPlus(params);
    //if (params.m_name=="F0_0_PP") return new Fzero_0_PiZeroPiPlus(params);
  }
  return NULL;
}


/*


class Fplus_PiPlusPiZero : public FF_0_PP_Base {
  void Construct();
public:
  Fplus_PiPlusPiZero(const FF_Parameters & params) : FF_0_PP_Base(params) {
    Construct();
  }
};
  

class Fzero_PiPlusPiZero : public FF_0_PP_Base {
public:
  Fzero_PiPlusPiZero(const FF_Parameters & params) : FF_0_PP_Base(params) {
    Construct();
  }
};

void FF_0_PP_Base::FixMode() {
  if (m_flavs[m_pi[0]].Kfcode()==kf_pi &&
      m_flavs[m_pi[1]].Kfcode()==kf_pi_plus)
    m_PSmode = PSmode::pipi_plus;
  else if ( (m_flavs[m_pi[0]].Kfcode()==kf_K_L ||
	     m_flavs[m_pi[0]].Kfcode()==kf_K_S ||
	     m_flavs[m_pi[0]].Kfcode()==kf_K) &&
	    m_flavs[m_pi[1]].Kfcode()==kf_K_plus)
    m_PSmode = PSmode::KK_plus;
  else if ( (m_flavs[m_pi[0]].Kfcode()==kf_pi_plus && 
	     (m_flavs[m_pi[1]].Kfcode()==kf_K_S ||
	      m_flavs[m_pi[1]].Kfcode()==kf_K_L ||
	      m_flavs[m_pi[1]].Kfcode()==kf_K) ) ||
	    (m_flavs[m_pi[1]].Kfcode()==kf_pi_plus && 
	     (m_flavs[m_pi[0]].Kfcode()==kf_K_S ||
	      m_flavs[m_pi[0]].Kfcode()==kf_K_L ||
	      m_flavs[m_pi[0]].Kfcode()==kf_K) ) ||
	    (m_flavs[m_pi[0]].Kfcode()==kf_pi && 
	     m_flavs[m_pi[1]].Kfcode()==kf_K_plus) ||
	    (m_flavs[m_pi[1]].Kfcode()==kf_pi && 
	     m_flavs[m_pi[0]].Kfcode()==kf_K_plus) )
    m_PSmode = PSmode::Kpi_plus;
  else if ( m_flavs[m_pi[0]].Kfcode()==kf_pi_plus &&
	    m_flavs[m_pi[1]].Kfcode()==kf_eta )
    m_PSmode = PSmode::etapi_plus;
  else if ( m_flavs[m_pi[0]].Kfcode()==kf_pi_plus &&
	    m_flavs[m_pi[1]].Kfcode()==kf_eta_prime_958)
    m_PSmode = PSmode::etaprimepi_plus;
  else if ( m_flavs[m_pi[0]].Kfcode()==kf_eta &&
	    m_flavs[m_pi[1]].Kfcode()==kf_K_plus)
    m_PSmode = PSmode::Keta_plus;
  else if ( m_flavs[m_pi[0]].Kfcode()==kf_eta_prime_958 &&
	    m_flavs[m_pi[1]].Kfcode()==kf_K_plus)
    m_PSmode = PSmode::Ketaprime_plus;
  if (m_PSmode==PSmode::unknown) {
    msg_Out()<<"Weird flavours: "
	     <<m_flavs[m_pi[0]]<<" + "<<m_flavs[m_pi[1]]<<"\n";
    THROW(fatal_error,"Current called for illegal flavour combination.");
  }
}

void FF_0_PP_Base::FixNorm() {
  // global pre-factor: 1/sqrt(2) for pi_0 wave-function, V_ud for the
  // quark-level coupling producing a rho (or rho-resonance), 1/sqrt(2)
  // for the overall normalisation.
  double iso = 1.;
  switch (int(m_PSmode)) {
  case int(PSmode::KK_plus):
    iso = 1./sqrt(2.);
    break;
  case int(PSmode::Kpi_plus):
    if      (m_flavs[m_pi[0]].Kfcode()==kf_pi_plus) iso = 1./2.;
    else if (m_flavs[m_pi[0]].Kfcode()==kf_pi)      iso = 1./sqrt(2.);
    break;
  case int(PSmode::pipi_plus):
    iso = 1.;
    break;
  default: break;
  }
  double CKM = 1.;
  if (m_PSmode==PSmode::pipi_plus  || m_PSmode==PSmode::KK_plus ||
      m_PSmode==PSmode::etapi_plus || m_PSmode==PSmode::etaprimepi_plus)
    CKM = (*p_model)("Vud", Tools::Vud);
  else if (m_PSmode==PSmode::Kpi_plus || m_PSmode==PSmode::Keta_plus ||
	   m_PSmode==PSmode::Ketaprime_plus) 
    CKM = (*p_model)("Vus", Tools::Vus);
  m_norm = iso * CKM;
  if (m_norm<=0.) THROW(fatal_error,"Current with zero norm.");
}

*/

/*
class F0_0_PiPlusPiZero : public FF_0_PP_Base {
  void Construct() {};
public:
  F0_0_PiPlusPiZero(const FF_Parameters & params) : 
    FF_0_PP_Base(params) { Construct(); }
  ~F0_0_PiPlusPiZero() {}
};
*/

///////////////////////////////////////////////////////////////////////////
//
// Getters
//
///////////////////////////////////////////////////////////////////////////

