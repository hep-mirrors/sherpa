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
  m_i_plus(-1), m_i_minus(-1), m_i_zero(-1), p_amp(NULL)
{
  FixParameters(params);
  Construct();
}

FF_0_Isoscalar3Pi::~FF_0_Isoscalar3Pi() {
  if (p_amp) { delete p_amp; p_amp = NULL; }
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

  BuildAmplitude();

  msg_Out()<<"FF_0_Isoscalar3Pi: hep-ph/0512180 Eq.(10) over registry rho "
	   <<"line shapes; A-F = ";
  for (size_t i(0);i<6;++i) msg_Out()<<m_c[i]<<(i<5?", ":"\n");
}

void FF_0_Isoscalar3Pi::BuildAmplitude() {
  // sub_0 = the pi+ pi- pairing, sub_1 = pi+ pi0, sub_2 = pi- pi0.  The
  // neutral rho takes the first, the charged one the other two -- which is
  // exactly the assignment a single-invariant sum could not express.
  using ia = invariant_arg;
  Sum_Term * h770  = new Sum_Term();
  h770 ->Add(p_r770_0, ia::sub_0); h770 ->Add(p_r770_c, ia::sub_1);
  h770 ->Add(p_r770_c, ia::sub_2);
  Sum_Term * h1450 = new Sum_Term();
  h1450->Add(p_r1450_0,ia::sub_0); h1450->Add(p_r1450_c,ia::sub_1);
  h1450->Add(p_r1450_c,ia::sub_2);
  Sum_Term * h1700 = new Sum_Term();
  h1700->Add(p_r1700_0,ia::sub_0); h1700->Add(p_r1700_c,ia::sub_1);
  h1700->Add(p_r1700_c,ia::sub_2);

  // A..D share rho(770); E takes rho(1450) with the phi alone and F takes
  // rho(1700) with the omega'' alone.  That selectivity is the whole point.
  Sum_Term * iso = new Sum_Term();
  for (size_t i(0);i<4;++i) iso->Add(m_iso[i],ia::total,m_c[i]);

  Product_Term * t1 = new Product_Term(); t1->Add(iso);   t1->Add(h770);
  Product_Term * t2 = new Product_Term();
  t2->Add(m_iso[1],ia::total,m_c[4]); t2->Add(h1450);
  Product_Term * t3 = new Product_Term();
  t3->Add(m_iso[3],ia::total,m_c[5]); t3->Add(h1700);

  p_amp = new Sum_Term();
  p_amp->Add(t1); p_amp->Add(t2); p_amp->Add(t3);
}

Complex FF_0_Isoscalar3Pi::operator()(const Vec4D_Vector & moms) {
  const Vec4D & pp=moms[m_i_plus], & pm=moms[m_i_minus], & p0=moms[m_i_zero];
  const Invariants k((pp+pm+p0).Abs2(),
                     { (pp+pm).Abs2(), (pp+p0).Abs2(), (pm+p0).Abs2() });
  return (*p_amp)(k);
}

Vec4C FF_0_Isoscalar3Pi::Current(const Vec4D_Vector & moms) {
  const Vec4D & pp=moms[m_i_plus], & pm=moms[m_i_minus], & p0=moms[m_i_zero];
  return Vec4C(cross(pp,pm,p0))*(*this)(moms);
}
