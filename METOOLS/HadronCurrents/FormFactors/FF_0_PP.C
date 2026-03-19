#include "METOOLS/HadronCurrents/FormFactors/FF_0_PP.H"
#include "METOOLS/HadronCurrents/Tools.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;


void FF_0_PP_Base::FixMode() {
  if (m_flavs[m_pi[0]].Kfcode()==kf_pi &&
      m_flavs[m_pi[1]].Kfcode()==kf_pi_plus)     m_mode = FF_0_PP_mode::pipi_plus;
  else if (m_flavs[m_pi[0]].Kfcode()==kf_K &&
	   m_flavs[m_pi[1]].Kfcode()==kf_K_plus) m_mode = FF_0_PP_mode::KK_plus;
}

Complex FF_0_PP_Base::operator()(const ATOOLS::Vec4D_Vector& moms) {
  double Q2 = (moms[m_pi[0]]+moms[m_pi[1]]).Abs2();
  switch (m_ffmodel) {
  case ff_model::none:    return Complex(1.,0.);
  case ff_model::KS:      return m_norm * (*p_props)(Q2);
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
  double  m_fpi, m_mpi2, m_mK2;
  double  m_mV, m_mV2, m_mVp, m_mVp2, m_mVpp, m_mVpp2, m_GVp, m_GVpp;
  Complex m_gamma, m_delta;
  
  void    FixParameters(const FF_Parameters & params);
  void    Construct();
  Complex FF_RChiPT(const double & s);
  Complex A(const double & m2,const double & s,const double & mu2);
  double  Gamma_V(const double & s);
  double  Gamma_Vp(const double & s);
  double  Gamma_Vpp(const double & s);
public :
  Fplus_0_PiZeroPiPlus(const FF_Parameters & params);
};

Fplus_0_PiZeroPiPlus::Fplus_0_PiZeroPiPlus(const FF_Parameters & params)  :
  FF_0_PP_Base(params),
  m_fpi((*params.p_model)("fpi",0.1307)/sqrt(2.)),
  m_mpi2(sqr(Flavour(kf_pi_plus).HadMass())),
  m_mK2(sqr(Flavour(kf_K_plus).HadMass()))
{
  FixParameters(params);
  Construct();
}

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

Complex Fplus_0_PiZeroPiPlus::FF_RChiPT(const double & s) {
  Complex V   = ( (m_mV2 - s*(m_gamma+m_delta))                       *
		  Complex(m_mV2-s, m_mV*Gamma_V(s))/
		  (sqr(m_mV2-s)+sqr(m_mV*Gamma_V(s)))             *
		  exp(-s/(96*sqr(M_PI*m_fpi)) *
		      (A(m_mpi2,s,m_mV2)+1./2.*A(m_mK2,s,m_mV2)).real()) );
  Complex Vp  = ( s*m_gamma                                             *
		  Complex(m_mVp2-s, m_mVp*Gamma_Vp(s))/
		  (sqr(m_mVp2-s)+sqr(m_mVp*Gamma_Vp(s)))           *
		  exp(-s*Gamma_Vp(m_mVp2)/
		      (M_PI*pow(m_mVp*
				sqrt(1.-4.*m_mpi2/m_mVp2),3./2.)) *
		      A(m_mpi2,s,m_mV2).real())                           );
  Complex Vpp = ( s*m_delta                                             *
		  Complex(m_mVpp2-s, m_mVpp*Gamma_Vpp(s))/
		  (sqr(m_mVpp2-s)+sqr(m_mVpp*Gamma_Vpp(s)))        *
		  exp(-s*Gamma_Vpp(m_mVpp2)/
		      (M_PI*pow(m_mVpp*
				sqrt(1.-4.*m_mpi2/m_mVpp2),3./2.)) *
		      A(m_mpi2,s,m_mV2).real())                          );
  return V + Vp + Vpp;
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

//////////////////////////////////////////////////////////////////////////////////
//
// Form factors for
//   * pi pi & K K vanish due to identical masses
// - Kuehn-Santamaria
//   * K pi: Z.Phys.C 72 (1996) 619 (scalar form factor)
//     (https://doi.org/10.1007/s002880050284)
//
// Todo: implement all scalar form factors here!
//
//////////////////////////////////////////////////////////////////////////////////

class Fzero_0_PiZeroPiPlus : public FF_0_PP_Base {
  double  m_fpi, m_mpi2, m_mK2;
  double  m_mV, m_mV2, m_mVp, m_mVp2, m_mVpp, m_mVpp2, m_GVp, m_GVpp;
  Complex m_gamma, m_delta;
  
  void    FixParameters(const FF_Parameters & params);
  void    Construct();
  Complex FF_RChiPT(const double & Q2) { return complex(0.,0.); }
public :
  Fzero_0_PiZeroPiPlus(const FF_Parameters & params)  :
    FF_0_PP_Base(params),
    m_fpi((*params.p_model)("fpi",0.1307)/sqrt(2.)),
    m_mpi2(sqr(Flavour(kf_pi_plus).HadMass())),
    m_mK2(sqr(Flavour(kf_K_plus).HadMass()))
 {
    FixParameters(params);
    Construct();
  }
};

void Fzero_0_PiZeroPiPlus::FixParameters(const FF_Parameters & params)  {}

void Fzero_0_PiZeroPiPlus::Construct()  {}


DECLARE_FF_GETTER(FF_0_PP_Base,"FF_0_PP")

FormFactor_Base * ATOOLS::Getter<FormFactor_Base,FF_Parameters,
				 FF_0_PP_Base>:: 
operator()(const METOOLS::FF_Parameters &params) const
{
  msg_Out()<<METHOD<<":\n";
  size_t Nmesons = 0;
  for (size_t i=0;i<params.m_flavs.size();i++) {
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
    if (params.m_name=="F+_0_PP") return new Fplus_0_PiZeroPiPlus(params);
    if (params.m_name=="F0_0_PP") return new Fzero_0_PiZeroPiPlus(params);
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

