#include "METOOLS/HadronCurrents/FormFactors/FF_0_PP.H"
#include "METOOLS/HadronCurrents/FormFactors/Line_Shapes.H"
#include "METOOLS/HadronCurrents/Tools.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;


FF_0_PP_Base::FF_0_PP_Base(const FF_Parameters & params) :
  FormFactor_Base(params),
  m_PSmode(PSmode::unknown), m_norm(0.), p_props(NULL) {
  FixMode();
  FixNorm();
  Construct();
}

FF_0_PP_Base::~FF_0_PP_Base() {
  if (p_props) delete p_props;
}


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

Complex FF_0_PP_Base::operator()(const ATOOLS::Vec4D_Vector& moms) {
  if (p_props) {
    double Q2 = (moms[m_pi[0]]+moms[m_pi[1]]).Abs2();
    msg_Out()<<METHOD<<"("<<m_flavs[m_pi[0]]<<"+"<<m_flavs[m_pi[1]]<<", "
	     <<"model = "<<int(m_ffmodel)<<"):\n"
	     <<"* moms = "<<moms[m_pi[0]]<<"+"<<moms[m_pi[1]]<<" --> Q^2 = "<<Q2<<"\n";
    switch (m_ffmodel) {
    case ff_model::none:
      return Complex(1.,0.);
    case ff_model::KS:
      msg_Out()<<"trying to evaluate propagator structure: ["<<p_props<<"].\n";
      return (*p_props)(Q2);
    case ff_model::RChiPT:
    case ff_model::unknown:
    default:
      break;
    }
  }
  return Complex(0.,0.);
}

///////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////
Fplus_0_PP::Fplus_0_PP(const FF_Parameters & params) :
  FF_0_PP_Base(params) {
  msg_Out()<<METHOD<<"("<<m_name<<"): "
	   <<"["<<m_flavs[m_pi[0]]<<" + "<<m_flavs[m_pi[1]]<<"] "
	   <<"--> mode = "<<int(m_PSmode)<<", norm = "<<m_norm<<"\n";
  Construct();
}

void Fplus_0_PP::Construct() {
  msg_Out()<<METHOD<<"("<<int(m_PSmode)<<")\n";
  if (m_PSmode==PSmode::pipi_plus) {
    Propagator_Base * rho770  =
      new BreitWigner(LineShapes->Get(Flavour(kf_rho_770_plus)));
    Propagator_Base * rho1450 =
      new BreitWigner(LineShapes->Get(Flavour(kf_rho_1450_plus)));
    Propagator_Base * rho1700 =
      new BreitWigner(LineShapes->Get(Flavour(kf_rho_1700_plus)));
    p_props = new Summed_Propagator();
    p_props->Add(rho770,  Complex(  1.000, 0.000));
    p_props->Add(rho1450, Complex( -0.103, 0.000));
    p_props->Add(rho1700, Complex( -0.037, 0.000));
  }
}

///////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////
F0_0_PP::F0_0_PP(const FF_Parameters & params) :
  FF_0_PP_Base(params) {
  msg_Out()<<METHOD<<"("<<m_name<<"): "
	   <<"["<<m_flavs[m_pi[0]]<<" + "<<m_flavs[m_pi[1]]<<"] "
	   <<"--> mode = "<<int(m_PSmode)<<", norm = "<<m_norm<<"\n";
}

DEFINE_FF_GETTER(Fplus_0_PP,"F+_0_PP")
DEFINE_FF_GETTER(F0_0_PP,"F0_0_PP")
