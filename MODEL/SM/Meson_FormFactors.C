#include "METOOLS/Explicit/Form_Factor.H"
#include "METOOLS/Explicit/Current.H"
#include "METOOLS/Explicit/Vertex.H"
#include "METOOLS/FormFactors/Propagator.H"
#include "METOOLS/FormFactors/Line_Shapes.H"

#include "MODEL/Main/Single_Vertex.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include <set>

namespace METOOLS {
  enum FF_0_PP_mode {
    pipi_plus       = 1,
    KK_plus         = 2,
    Kpi_plus        = 11,
    pipi_0          = 101,
    KK_0            = 102,
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
  struct GS_Parameters { double m, Gamma, c, phi; };
  //                                    m       Gamma       c        phi
  const GS_Parameters s_rho      = { 0.77456, 0.14832, 1.0    , 0.     };
  const GS_Parameters s_rho1450  = { 1.4859 , 0.37360, 0.14104, 3.7797 };
  const GS_Parameters s_rho1700  = { 1.8668 , 0.30334, 0.0614 , 1.429  };
  const GS_Parameters s_rho2150  = { 2.2645 , 0.11327, 0.0047 , 0.921  };
  const GS_Parameters s_omega782 = { 0.78248, 0.00855, 0.00158, 0.075  };
  const GS_Parameters s_phi1020  = { 1.01947, 0.00425, 0.00045, 2.888  };

  inline Complex Weight(const GS_Parameters & p) {
    return p.c*Complex(cos(p.phi),sin(p.phi));
  }
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
  if (key.m_j[0]->Flav()==Flavour(kf_pi_plus) &&
      key.m_j[1]->Flav()==Flavour(kf_pi_plus).Bar() )     m_mode = FF_0_PP_mode::pipi_0;
  else if (key.m_j[0]->Flav()==Flavour(kf_K_plus) &&
	   key.m_j[1]->Flav()==Flavour(kf_K_plus).Bar() ) m_mode = FF_0_PP_mode::KK_0;
}

void FFVMD::Construct() {
  switch (int(m_mode)) {
  case int(FF_0_PP_mode::pipi_0):
    ConstructPionFormFactor();
    break;
  case int(FF_0_PP_mode::KK_0):
    // Recognised, but not parametrised yet.  Falling through to a silent
    // F = 1 here would look exactly like a working form factor.
    THROW(not_implemented,"No VMD form factor for K+ K- yet.");
  case int(FF_0_PP_mode::unknown):
  default:
    msg_Out()<<METHOD<<" yields no form factor.\n";
    break;
  }
}

void FFVMD::ConstructPionFormFactor() {
  // rho-omega and rho-phi mixing.  Both mixing terms carry an explicit
  // factor s/m^2 and therefore vanish at s = 0, so this sum is already
  // normalised - dividing it by (1 + c_omega + c_phi) would spoil F(0) = 1.
  METOOLS::Summed_Propagator * rhofac = new METOOLS::Summed_Propagator(false);
  rhofac->Add(new METOOLS::Unity(),Complex(1.,0.));
  rhofac->Add(new METOOLS::WeightedBreitWigner(LineShapes->Get(Flavour(kf_omega_782)),
					       resonance_type::fixed,
					       s_omega782.m,s_omega782.Gamma),
	      Weight(s_omega782));
  // No lineshape for the phi(1020) yet, so it runs off a fixed width.
  rhofac->Add(new METOOLS::WeightedBreitWigner(NULL,resonance_type::fixed,
					       s_phi1020.m,s_phi1020.Gamma),
	      Weight(s_phi1020));
  METOOLS::Multiplied_Propagator * rho = new METOOLS::Multiplied_Propagator();
  rho->Add(new METOOLS::GounarisSakurai(LineShapes->Get(Flavour(kf_rho_770)),
					resonance_type::GS,
					s_rho.m,s_rho.Gamma));
  rho->Add(rhofac);
  m_props.Add(rho,Weight(s_rho));
  m_props.Add(new METOOLS::GounarisSakurai(LineShapes->Get(Flavour(kf_rho_1450)),
					   resonance_type::GS,
					   s_rho1450.m,s_rho1450.Gamma),
	      Weight(s_rho1450));
  m_props.Add(new METOOLS::GounarisSakurai(LineShapes->Get(Flavour(kf_rho_1700)),
					   resonance_type::GS,
					   s_rho1700.m,s_rho1700.Gamma),
	      Weight(s_rho1700));
  // No lineshape for the rho(2150) either - the GS width is analytic anyway.
  m_props.Add(new METOOLS::GounarisSakurai(NULL,resonance_type::GS,
					   s_rho2150.m,s_rho2150.Gamma),
	      Weight(s_rho2150));
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

