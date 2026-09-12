#include "METOOLS/HadronCurrents/FormFactors/FF_0_Isoscalar3Pi.H"
#include "METOOLS/HadronCurrents/FormFactors/Line_Shapes.H"
#include "METOOLS/HadronCurrents/Tools.H"
#include "ATOOLS/Phys/Flavour.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

FF_0_Isoscalar3Pi::FF_0_Isoscalar3Pi(const FF_Parameters & params) :
  FormFactor_Base(params),
  p_r770_0(NULL),  p_r770_c(NULL),
  p_r1450_0(NULL), p_r1450_c(NULL),
  p_r1700_0(NULL), p_r1700_c(NULL),
  m_i_plus(-1), m_i_minus(-1), m_i_zero(-1)
{
  FixParameters(params);
  Construct();
}

FF_0_Isoscalar3Pi::~FF_0_Isoscalar3Pi() {
  while (!m_iso.empty()) { delete m_iso.back(); m_iso.pop_back(); }
  if (p_r770_0)  { delete p_r770_0;  p_r770_0  = NULL; }
  if (p_r770_c)  { delete p_r770_c;  p_r770_c  = NULL; }
  if (p_r1450_0) { delete p_r1450_0; p_r1450_0 = NULL; }
  if (p_r1450_c) { delete p_r1450_c; p_r1450_c = NULL; }
  if (p_r1700_0) { delete p_r1700_0; p_r1700_0 = NULL; }
  if (p_r1700_c) { delete p_r1700_c; p_r1700_c = NULL; }
}

void FF_0_Isoscalar3Pi::FixParameters(const FF_Parameters & params) {
  // Identify the pions by charge rather than by position, so the caller may
  // hand the indices over in any order.
  for (size_t i=0;i<m_pi.size();i++) {
    const Flavour & fl = m_flavs[m_pi[i]];
    if      (fl.Kfcode()==kf_pi_plus && fl.IntCharge()>0) m_i_plus  = m_pi[i];
    else if (fl.Kfcode()==kf_pi_plus && fl.IntCharge()<0) m_i_minus = m_pi[i];
    else if (fl.Kfcode()==kf_pi)                          m_i_zero  = m_pi[i];
  }
  if (m_i_plus<0 || m_i_minus<0 || m_i_zero<0)
    THROW(fatal_error,
	  "FF_0_Isoscalar3Pi needs exactly one pi+, one pi- and one pi0.");

  // A--F refitted against Belle 2024 with the line shapes below.  The
  // published values are -0.77, -1.12, -0.59 for C, D, F; they were fitted
  // against that paper's own inline propagators and do not transfer.
  const double c[6] = { 18.20, -0.87, -0.5785, -1.2062, -0.72, -0.3947 };
  const string tags[6] = { "A", "B", "C", "D", "E", "F" };
  for (size_t i(0);i<6;++i) {
    m_c[i] = c[i];
    if (p_model) m_c[i] = (*p_model)("EE3Pi_"+tags[i],m_c[i]);
  }

  // Isoscalar poles.  Defaults: the particle-database mass with the width the
  // model wants; omega' and omega'' additionally take hep-ph/0512180's own
  // masses, which differ from the PDG states by ~40 MeV.
  const long int kf[4] = { kf_omega_782, kf_phi_1020,
			   kf_omega_1420, kf_omega_1600 };
  const double defM[4] = { -1., -1., 1.375, 1.631 };   // <0 -> particle database
  const double defG[4] = { 0.00868, 0.004249, 0.250, 0.245 };
  const string nm[4] = { "omega", "phi", "omegaP", "omegaPP" };
  for (size_t i(0);i<4;++i) {
    m_M[i] = defM[i]>0. ? defM[i] : Flavour(kf[i]).HadMass();
    m_G[i] = defG[i];
    if (p_model) {
      m_M[i] = (*p_model)("EE3Pi_M_"+nm[i],m_M[i]);
      m_G[i] = (*p_model)("EE3Pi_G_"+nm[i],m_G[i]);
    }
  }
}

void FF_0_Isoscalar3Pi::Construct() {
  // The four isoscalars, as CONSTANT-width Breit-Wigners at model-level poles.
  // See the header: their registered widths run, which is right near the pole
  // and badly wrong far above it, and the paper's omega'/omega'' are its own
  // fit parameters rather than the PDG states.
  const long int kf[4] = { kf_omega_782, kf_phi_1020,
			   kf_omega_1420, kf_omega_1600 };
  for (size_t i(0);i<4;++i) {
    m_iso.push_back(new FixedBreitWigner(m_M[i],m_G[i]));
    msg_Out()<<"FF_0_Isoscalar3Pi: "<<Flavour(kf[i])
	     <<" -> constant-width pole M = "<<m_M[i]
	     <<", Gamma = "<<m_G[i]<<"\n";
  }

  // The rho tower FROM THE REGISTRY, per flavour rather than as a bundled
  // channel: Eq. (10) pairs rho(1450) with the phi and rho(1700) with the
  // omega'' specifically, so the three recurrences must stay separable.
  //
  // On this branch the registry hands out widths rather than line shapes, so
  // the shape is chosen here at the call site.  rho(770) uses
  // GounarisSakuraiM, NOT BreitWigner's resonance_type::GS: the two build the
  // denominator differently (M*Gamma against sqrt(s)*Gamma) and the couplings
  // above were fitted against this one.  See Propagator.H.
  const double mpi = Flavour(kf_pi_plus).HadMass();
  p_r770_0  = new GounarisSakuraiM(Flavour(kf_rho_770),mpi);
  p_r770_c  = new GounarisSakuraiM(Flavour(kf_rho_770_plus),mpi);
  p_r1450_0 = new BreitWigner(LineShapes->Get(Flavour(kf_rho_1450)));
  p_r1450_c = new BreitWigner(LineShapes->Get(Flavour(kf_rho_1450_plus)));
  p_r1700_0 = new BreitWigner(LineShapes->Get(Flavour(kf_rho_1700)));
  p_r1700_c = new BreitWigner(LineShapes->Get(Flavour(kf_rho_1700_plus)));

  msg_Out()<<"FF_0_Isoscalar3Pi: hep-ph/0512180 Eq.(10) over registry rho "
	   <<"line shapes; A-F = ";
  for (size_t i(0);i<6;++i) msg_Out()<<m_c[i]<<(i<5?", ":"\n");
}

Complex FF_0_Isoscalar3Pi::operator()(const Vec4D_Vector & moms) {
  const Vec4D & pp=moms[m_i_plus], & pm=moms[m_i_minus], & p0=moms[m_i_zero];
  const double q2 =(pp+pm+p0).Abs2();
  const double spm=(pp+pm).Abs2(), sp0=(pp+p0).Abs2(), sm0=(pm+p0).Abs2();

  // All-plus over the three pairings; see the header on why that is forced.
  const Complex h770  = (*p_r770_0 )(spm)+(*p_r770_c )(sp0)+(*p_r770_c )(sm0);
  const Complex h1450 = (*p_r1450_0)(spm)+(*p_r1450_c)(sp0)+(*p_r1450_c)(sm0);
  const Complex h1700 = (*p_r1700_0)(spm)+(*p_r1700_c)(sp0)+(*p_r1700_c)(sm0);

  const Complex w782 =(*m_iso[0])(q2), wphi =(*m_iso[1])(q2);
  const Complex w1420=(*m_iso[2])(q2), w1650=(*m_iso[3])(q2);

  // Eq. (10): a sum of products.  E pairs rho(1450) with the phi alone and F
  // pairs rho(1700) with the omega'' alone -- that selectivity is the point.
  return ( m_c[0]*w782 + m_c[1]*wphi + m_c[2]*w1420 + m_c[3]*w1650 ) * h770
       +   m_c[4]*wphi  * h1450
       +   m_c[5]*w1650 * h1700;
}

Vec4C FF_0_Isoscalar3Pi::Current(const Vec4D_Vector & moms) {
  const Vec4D & pp=moms[m_i_plus], & pm=moms[m_i_minus], & p0=moms[m_i_zero];
  return Vec4C(cross(pp,pm,p0))*(*this)(moms);
}
