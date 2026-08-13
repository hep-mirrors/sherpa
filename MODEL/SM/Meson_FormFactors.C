#include "METOOLS/Explicit/Form_Factor.H"
#include "METOOLS/Explicit/Current.H"
#include "METOOLS/Explicit/Vertex.H"
#include "METOOLS/FormFactors/Propagator.H"
#include "METOOLS/FormFactors/Line_Shapes.H"

#include "MODEL/Main/Single_Vertex.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Settings.H"
#include <set>

namespace METOOLS {
  enum FF_0_PP_mode {
    pipi_plus       = 1,
    KK_plus         = 2,
    Kpi_plus        = 11,
    pipi_0          = 101,
    KK_0            = 102,
    pipipi_0        = 201,
    V_pi            = 301,
    unknown         = 999
  };

  class FFVMD: public Form_Factor {
  private:
    METOOLS::Summed_Propagator m_props;
    FF_0_PP_mode m_mode;
    int          m_pos;
    void FixMode(const Vertex_Key &key);
    void Construct();
    void ConstructPionFormFactor();
    void ConstructKaonFormFactor();
    void ConstructThreePionFormFactor();
    void ConstructVectorPionFormFactor();
  public:
    FFVMD(const Vertex_Key &key);
    Complex FF();
  };// end of class FFPoint
}// end of namespace METOOLS


using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

/////////////////////////////////////////////////////////////////////
//
// Gounaris-Sakurai parameters of the pion form factor.
//
// The amplitudes and phases below come out of a fit that determined
// the masses and widths listed here at the same time, so the two sets
// only make sense together.  In particular they must not be taken from
// the particle table of the low-energy model (rho(770): m = 0.77,
// Gamma = 0.1507, i.e. the PDG values): that shifts the rho peak by
// about 5 MeV and undershoots the BESIII pi+ pi- cross section by up
// to 10% above 0.75 GeV.
//
/////////////////////////////////////////////////////////////////////
namespace {
  // Any subset of these can be tuned from the run card, keeping the defaults
  // for everything left out, e.g.
  //
  //   PION_FORM_FACTOR:
  //     rho(770):   {Mass: 0.7755, Width: 0.1494}
  //     omega(782): {Amplitude: 0.00205, Phase: 0.287}
  //
  // Phases are in radians.  Note that the defaults belong together as one
  // fit: retuning a mass or width in isolation will generally need the
  // amplitudes and phases refitted along with it.
  class GS_Parameters {
    const std::string m_block, m_tag;
    const double m_m0, m_Gamma0, m_c0, m_phi0;  // fitted defaults, never change
    double       m_m,  m_Gamma,  m_c,  m_phi;   // values in use
  public:
    GS_Parameters(const std::string & block,const std::string & tag,
		  const double & m,const double & Gamma,
		  const double & c,const double & phi) :
      m_block(block), m_tag(tag), m_m0(m), m_Gamma0(Gamma), m_c0(c), m_phi0(phi),
      m_m(m), m_Gamma(Gamma), m_c(c), m_phi(phi) {}

    // Always registers the fitted value as the default, so this stays
    // idempotent however often the form factor is constructed.
    void Read() {
      ATOOLS::Scoped_Settings s
	{ ATOOLS::Settings::GetMainSettings()[m_block][m_tag] };
      m_m     = s["Mass"     ].SetDefault(m_m0    ).Get<double>();
      m_Gamma = s["Width"    ].SetDefault(m_Gamma0).Get<double>();
      m_c     = s["Amplitude"].SetDefault(m_c0    ).Get<double>();
      m_phi   = s["Phase"    ].SetDefault(m_phi0  ).Get<double>();
      msg_Debugging()<<"  "<<m_tag<<": m = "<<m_m<<", Gamma = "<<m_Gamma
		     <<", c = "<<m_c<<", phi = "<<m_phi<<"\n";
    }
    inline const double & Mass()  const { return m_m;     }
    inline const double & Width() const { return m_Gamma; }
    inline Complex Weight() const { return m_c*Complex(cos(m_phi),sin(m_phi)); }
  };
  //                            tag            m       Gamma       c        phi
  GS_Parameters s_rho      { "PION_FORM_FACTOR", "rho(770)"  , 0.77456, 0.14832, 1.0    , 0.     };
  GS_Parameters s_rho1450  { "PION_FORM_FACTOR", "rho(1450)" , 1.4859 , 0.37360, 0.14104, 3.7797 };
  GS_Parameters s_rho1700  { "PION_FORM_FACTOR", "rho(1700)" , 1.8668 , 0.30334, 0.0614 , 1.429  };
  GS_Parameters s_rho2150  { "PION_FORM_FACTOR", "rho(2150)" , 2.2645 , 0.11327, 0.0047 , 0.921  };
  GS_Parameters s_omega782 { "PION_FORM_FACTOR", "omega(782)", 0.78248, 0.00855, 0.00158, 0.075  };
  GS_Parameters s_phi1020  { "PION_FORM_FACTOR", "phi(1020)" , 1.01947, 0.00425, 0.00045, 2.888  };

  /////////////////////////////////////////////////////////////////////
  // Masses and widths are the PDG ones.  The phi amplitude is OZI
  // suppressed and is tuned, not fitted: 0.042 with a relative phase of pi
  // reproduces the measured ratio of the phi and omega peaks.  Together with
  // g_gammarhopi in LowEnergy_Model.C this puts both peaks within 10% of
  // BESIII_2019_I1773081.  It is a two-point tune and nothing more - see
  // doc/examples/Soft_QCD/LowEnergy/NOTES-3pi.md.
  /////////////////////////////////////////////////////////////////////
  GS_Parameters s_v_omega782
    { "THREE_PION_FORM_FACTOR", "omega(782)", 0.78266 , 0.00868 , 1.0, 0.   };
  GS_Parameters s_v_phi1020
    { "THREE_PION_FORM_FACTOR", "phi(1020)" , 1.019461, 0.004249, 0.042, M_PI };

  /////////////////////////////////////////////////////////////////////
  // Kaon form factor, hep-ph/0409080 eq. 64.  K+ = u sbar, so unlike the
  // pion form factor the photon sees THREE separate, additive SU(3)
  // channels with fixed relative weights 1/2 (isovector, rho family) and
  // 1/6, 1/3 (isoscalar, omega and phi families) - not one isospin current
  // with a small isospin-violating admixture, the way rho-omega/rho-phi
  // mixing works for pi+ pi-.
  //
  // The coefficients in the reference are real, with no explicit phase
  // factor; a negative one is encoded here as amplitude > 0 with
  // phase = pi, so the existing GS_Parameters Weight() = c*exp(i*phi)
  // machinery can be reused unchanged.
  //
  /////////////////////////////////////////////////////////////////////
  //                             tag              m         Gamma      c      phi
  GS_Parameters s_K_rho       { "KAON_FORM_FACTOR", "rho(770)"   , 0.77456 , 0.14832 , 1.195, 0.     };
  GS_Parameters s_K_rho1450   { "KAON_FORM_FACTOR", "rho(1450)"  , 1.4859  , 0.37360 , 0.112, M_PI   };
  GS_Parameters s_K_rho1700   { "KAON_FORM_FACTOR", "rho(1700)"  , 1.8668  , 0.30334 , 0.083, M_PI   };
  // c_rho(2150) = 1 - c_rho - c_rho(1450) - c_rho(1700) = 0 for the
  // defaults above: the sum rule that enforces F(0) = 1 leaves this term
  // with no contribution by default, but it stays tunable.
  GS_Parameters s_K_rho2150   { "KAON_FORM_FACTOR", "rho(2150)"  , 2.2645  , 0.11327 , 0.0  , 0.     };
  GS_Parameters s_K_omega     { "KAON_FORM_FACTOR", "omega(782)" , 0.78248 , 0.00855 , 1.195, 0.     };
  GS_Parameters s_K_omega1420 { "KAON_FORM_FACTOR", "omega(1420)", 1.410   , 0.290   , 0.112, M_PI   };
  GS_Parameters s_K_omega1650 { "KAON_FORM_FACTOR", "omega(1650)", 1.67    , 0.315   , 0.083, M_PI   };
  GS_Parameters s_K_phi       { "KAON_FORM_FACTOR", "phi(1020)"  , 1.01947 , 0.00425 , 1.018, 0.     };
  GS_Parameters s_K_phi1680   { "KAON_FORM_FACTOR", "phi(1680)"  , 1.680   , 0.150   , 0.018, M_PI   };
}

/////////////////////////////////////////////////////////////////////
//
// VMD form factors, basically assuming some (interfering)
// vector-meson propagators.
//
/////////////////////////////////////////////////////////////////////
FFVMD::FFVMD(const Vertex_Key &key):
  Form_Factor("VMD",key),
  m_props(), m_mode(FF_0_PP_mode::unknown),
  m_pos(-1) {
  // Find the position of the (elementary) particle:
  // Normally this would be the photon, but I'll assume it is anything
  // that is not a hadron.
  // Also: extract the flavours in the vertex to look up parameters
  // in the model look-up table.
  std::set<kf_code> kfs;
  msg_Out()<<METHOD<<"("<<key.m_j.size()<<"): "<<key.m_j[0]->Flav()<<", "
	   <<key.m_j[1]->Flav()<<"\n";
  msg_Out()<<key.p_mv->in[0]<<", "<<key.p_mv->in[1]<<", "<<key.p_mv->in[2]<<"\n";
  for (size_t i(0);i<key.m_j.size();++i) {
    if (!key.m_j[i]->Flav().IsHadron()) m_pos=i;
    else kfs.insert(key.m_j[i]->Flav().Kfcode());
  }
  msg_Out()<<METHOD<<"("<<key.m_j.size()<<"): "<<key.m_j[0]->Flav()<<", "
	   <<key.m_j[1]->Flav()<<"\n"; //<<", "<<key.m_j[2]->Flav()
  FixMode(key);
  Construct();
  msg_Out()<<METHOD<<":\n";
}

void FFVMD::FixMode(const Vertex_Key &key) {
  // The three-pion vertex has to be identified from the vertex legs rather
  // than from key.m_j: m_j only holds the incoming currents, so a three-point
  // vertex appears there with two of its three legs, and which two depends on
  // the recursion step.
  std::multiset<kf_code> kfs;
  for (size_t i(0);i<key.p_mv->in.size();++i)
    if (key.p_mv->in[i].IsHadron()) kfs.insert(key.p_mv->in[i].Kfcode());
  if (kfs==std::multiset<kf_code>{kf_pi,kf_pi_plus,kf_pi_plus}) {
    m_mode = FF_0_PP_mode::pipipi_0;
    return;
  }
  // gamma rho pi: the anomalous vertex of the three-pion chain.  This is the
  // only vertex there that sees the photon virtuality, so it carries the
  // whole s dependence.
  if (kfs==std::multiset<kf_code>{kf_rho_770,kf_pi} ||
      kfs==std::multiset<kf_code>{kf_rho_770_plus,kf_pi_plus}) {
    m_mode = FF_0_PP_mode::V_pi;
    return;
  }
  if      (kfs==std::multiset<kf_code>{kf_pi_plus,kf_pi_plus}) m_mode = FF_0_PP_mode::pipi_0;
  else if (kfs==std::multiset<kf_code>{kf_K_plus ,kf_K_plus }) m_mode = FF_0_PP_mode::KK_0;
}

void FFVMD::Construct() {
  switch (int(m_mode)) {
  case int(FF_0_PP_mode::pipi_0):
    ConstructPionFormFactor();
    break;
  case int(FF_0_PP_mode::pipipi_0):
    ConstructThreePionFormFactor();
    break;
  case int(FF_0_PP_mode::V_pi):
    ConstructVectorPionFormFactor();
    break;
  case int(FF_0_PP_mode::KK_0):
    ConstructKaonFormFactor();
    break;
  case int(FF_0_PP_mode::unknown):
  default:
    msg_Out()<<METHOD<<" yields no form factor.\n";
    break;
  }
}

void FFVMD::ConstructPionFormFactor() {
  msg_Debugging()<<METHOD<<": Gounaris-Sakurai parameters\n";
  for (GS_Parameters * p : { &s_rho, &s_rho1450, &s_rho1700, &s_rho2150,
			     &s_omega782, &s_phi1020 }) p->Read();
  // rho-omega and rho-phi mixing.  Both mixing terms carry an explicit
  // factor s/m^2 and therefore vanish at s = 0, so this sum is already
  // normalised - dividing it by (1 + c_omega + c_phi) would spoil F(0) = 1.
  METOOLS::Summed_Propagator * rhofac = new METOOLS::Summed_Propagator(false);
  rhofac->Add(new METOOLS::Unity(),Complex(1.,0.));
  rhofac->Add(new METOOLS::WeightedBreitWigner(LineShapes->Get(Flavour(kf_omega_782)),
					       resonance_type::fixed,
					       s_omega782.Mass(),s_omega782.Width()),
	      s_omega782.Weight());
  // No lineshape for the phi(1020) yet, so it runs off a fixed width.
  rhofac->Add(new METOOLS::WeightedBreitWigner(NULL,resonance_type::fixed,
					       s_phi1020.Mass(),s_phi1020.Width()),
	      s_phi1020.Weight());
  METOOLS::Multiplied_Propagator * rho = new METOOLS::Multiplied_Propagator();
  rho->Add(new METOOLS::GounarisSakurai(LineShapes->Get(Flavour(kf_rho_770)),
					resonance_type::GS,
					s_rho.Mass(),s_rho.Width()));
  rho->Add(rhofac);
  m_props.Add(rho,s_rho.Weight());
  m_props.Add(new METOOLS::GounarisSakurai(LineShapes->Get(Flavour(kf_rho_1450)),
					   resonance_type::GS,
					   s_rho1450.Mass(),s_rho1450.Width()),
	      s_rho1450.Weight());
  m_props.Add(new METOOLS::GounarisSakurai(LineShapes->Get(Flavour(kf_rho_1700)),
					   resonance_type::GS,
					   s_rho1700.Mass(),s_rho1700.Width()),
	      s_rho1700.Weight());
  // No lineshape for the rho(2150) either - the GS width is analytic anyway.
  m_props.Add(new METOOLS::GounarisSakurai(NULL,resonance_type::GS,
					   s_rho2150.Mass(),s_rho2150.Width()),
	      s_rho2150.Weight());
}


void FFVMD::ConstructKaonFormFactor() {
  msg_Debugging()<<METHOD<<": Gounaris-Sakurai parameters\n";
  for (GS_Parameters * p : { &s_K_rho, &s_K_rho1450, &s_K_rho1700, &s_K_rho2150,
			     &s_K_omega, &s_K_omega1420, &s_K_omega1650,
			     &s_K_phi, &s_K_phi1680 }) p->Read();
  //////////////////////////////////////////////////////////////////////
  // Three independent, additive SU(3) channels (see the comment on the
  // s_K_* parameters above) - not a single normalised isospin current, so
  // this is built with normalise=false and added to m_props as a single
  // term of weight 1.  m_props always renormalises by its own weight sum,
  // but with only one term of weight 1 that division is a no-op; the
  // F(0) = 1 normalisation here comes entirely from the sum rule
  // c_rho(2150) = 1 - c_rho - c_rho(1450) - c_rho(1700) baked into the
  // coefficients themselves.
  //////////////////////////////////////////////////////////////////////
  METOOLS::Summed_Propagator * kaon = new METOOLS::Summed_Propagator(false);
  // rho family: the same resonances, and the same pi pi running width, as
  // the pion channel - rho decays via pi pi regardless of which form
  // factor is asking about it.
  kaon->Add(new METOOLS::GounarisSakurai(LineShapes->Get(Flavour(kf_rho_770)),
					 resonance_type::GS,
					 s_K_rho.Mass(),s_K_rho.Width()),
	    0.5*s_K_rho.Weight());
  kaon->Add(new METOOLS::GounarisSakurai(LineShapes->Get(Flavour(kf_rho_1450)),
					 resonance_type::GS,
					 s_K_rho1450.Mass(),s_K_rho1450.Width()),
	    0.5*s_K_rho1450.Weight());
  kaon->Add(new METOOLS::GounarisSakurai(LineShapes->Get(Flavour(kf_rho_1700)),
					 resonance_type::GS,
					 s_K_rho1700.Mass(),s_K_rho1700.Width()),
	    0.5*s_K_rho1700.Weight());
  kaon->Add(new METOOLS::GounarisSakurai(NULL,resonance_type::GS,
					 s_K_rho2150.Mass(),s_K_rho2150.Width()),
	    0.5*s_K_rho2150.Weight());
  // omega family: fixed widths.  None of these decay predominantly to two
  // pseudoscalars, so a running width makes no more sense here than it
  // does for omega in the pion channel.
  kaon->Add(new METOOLS::FixedBreitWigner(NULL,s_K_omega.Mass(),s_K_omega.Width()),
	    s_K_omega.Weight()/6.);
  kaon->Add(new METOOLS::FixedBreitWigner(NULL,s_K_omega1420.Mass(),s_K_omega1420.Width()),
	    s_K_omega1420.Weight()/6.);
  kaon->Add(new METOOLS::FixedBreitWigner(NULL,s_K_omega1650.Mass(),s_K_omega1650.Width()),
	    s_K_omega1650.Weight()/6.);
  // phi family: also fixed width here.  phi(1020) sits only 33 MeV above
  // the K+K- threshold itself, where a running width from the K Kbar loop
  // (rather than the pi pi one GounarisSakurai assumes) would matter most -
  // hep-ph/0409080 does include that.  Left as the obvious next
  // refinement; the fixed-width approximation is the same simplification
  // already made for phi in the pion channel.
  kaon->Add(new METOOLS::FixedBreitWigner(NULL,s_K_phi.Mass(),s_K_phi.Width()),
	    s_K_phi.Weight()/3.);
  kaon->Add(new METOOLS::FixedBreitWigner(NULL,s_K_phi1680.Mass(),s_K_phi1680.Width()),
	    s_K_phi1680.Weight()/3.);
  m_props.Add(kaon,Complex(1.,0.));
}



void FFVMD::ConstructThreePionFormFactor() {
  msg_Debugging()<<METHOD<<": Gounaris-Sakurai parameters\n";
  for (GS_Parameters * p : { &s_rho, &s_rho1450, &s_rho1700 }) p->Read();
  //////////////////////////////////////////////////////////////////////
  // FIRST STAB - the rho-family propagator sum in the two-pion
  // sub-channel, i.e. the prop_ij(s_ij) factors of the rho-pi picture that
  // METOOLS::V1minus_PPP_Arg::dg() already uses for omega -> 3 pi.  At this
  // vertex the off-shell leg is the charged pion going to the photon, so its
  // virtuality is the invariant mass squared of the other two pions, which
  // is the argument a rho propagator wants.
  //
  // This is NOT the gamma* -> 3 pi form factor. the
  // omega/phi peak right needs the SSS vertex replaced by a proper
  // gamma-pi-pi-pi current (compare VVP_LC / TAUPI_LC), not a form factor.
  //////////////////////////////////////////////////////////////////////
  m_props.Add(new METOOLS::GounarisSakurai(LineShapes->Get(Flavour(kf_rho_770)),
					   resonance_type::GS,
					   s_rho.Mass(),s_rho.Width()),
	      s_rho.Weight());
  m_props.Add(new METOOLS::GounarisSakurai(LineShapes->Get(Flavour(kf_rho_1450)),
					   resonance_type::GS,
					   s_rho1450.Mass(),s_rho1450.Width()),
	      s_rho1450.Weight());
  m_props.Add(new METOOLS::GounarisSakurai(LineShapes->Get(Flavour(kf_rho_1700)),
					   resonance_type::GS,
					   s_rho1700.Mass(),s_rho1700.Width()),
	      s_rho1700.Weight());
}

void FFVMD::ConstructVectorPionFormFactor() {
  msg_Debugging()<<METHOD<<": three-pion isoscalar parameters\n";
  for (GS_Parameters * p : { &s_v_omega782, &s_v_phi1020 }) p->Read();
  // gamma -> omega/phi -> rho pi.  Fixed widths: neither resonance is narrow
  // enough here for the shape of a running width to matter much, and there
  // is no phi lineshape to run one off anyway.  Normalised to F(0) = 1 by
  // Summed_Propagator, which divides by the sum of the weights.
  m_props.Add(new METOOLS::FixedBreitWigner(NULL,s_v_omega782.Mass(),
					    s_v_omega782.Width()),
	      s_v_omega782.Weight());
  m_props.Add(new METOOLS::FixedBreitWigner(NULL,s_v_phi1020.Mass(),
					    s_v_phi1020.Width()),
	      s_v_phi1020.Weight());
}

Complex FFVMD::FF() {
  Current *j = m_pos<0?p_v->JC():p_v->J(m_pos);
  // there was a minus sign before in Q2.
  double Q2  = j->P().Abs2();
  return m_props(Q2);
}



DECLARE_GETTER(FFVMD,"VMD",Form_Factor,Vertex_Key);
Form_Factor *Getter<Form_Factor,Vertex_Key,FFVMD>::
operator()(const Vertex_Key &args) const {
  msg_Out()<<METHOD<<"(size = "<<args.m_j.size()<<")\n";
  return new FFVMD(args);
}

void ATOOLS::Getter<Form_Factor,Vertex_Key,FFVMD>::
PrintInfo(std::ostream &str,const size_t width) const { str<<"VMD"; }

