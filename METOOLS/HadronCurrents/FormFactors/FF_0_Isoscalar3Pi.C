#include "METOOLS/HadronCurrents/FormFactors/FF_0_Isoscalar3Pi.H"
#include "METOOLS/HadronCurrents/FormFactors/Line_Shapes.H"
#include "METOOLS/HadronCurrents/FormFactors/KT3Pi_Tables.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "METOOLS/HadronCurrents/Tools.H"
#include "ATOOLS/Phys/Flavour.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

FF_0_Isoscalar3Pi::FF_0_Isoscalar3Pi(const FF_Parameters & params) :
  FormFactor_Base(params),
  p_r770_0(NULL),  p_r770_c(NULL),
  p_r1450_0(NULL), p_r1450_c(NULL),
  p_r1700_0(NULL), p_r1700_c(NULL), p_amp_nc(NULL),
  m_i_plus(-1), m_i_minus(-1), m_i_zero(-1), p_amp(NULL)
{
  FixParameters(params);
  Construct();
}

FF_0_Isoscalar3Pi::~FF_0_Isoscalar3Pi() {
  if (p_amp)    { delete p_amp;    p_amp    = NULL; }
  if (p_amp_nc) { delete p_amp_nc; p_amp_nc = NULL; }
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

  // Which parameterisation.  100 is Eq. (10), and the other 100s fall back to
  // it, as the family-fallback policy in FormFactor_Base.H says; 300 is
  // Khuri-Treiman, see the header.  Nothing else is implemented, and falling
  // back to a constant form factor would give a cross section that looks
  // plausible and means nothing -- so anything else throws.  Decided first,
  // because the default couplings depend on it.
  m_kt = (GetModel()==Model::KT);
  if (!m_kt && GetModel()!=Model::None && GetModel()!=Model::Point && GetModel()!=Model::KS)
    THROW(fatal_error,"FF_0_Isoscalar3Pi implements FormFactor 100 "
	  "(hep-ph/0512180 Eq. (10)) and 300 (Khuri-Treiman), not "+
	  ToString(int(m_ffmodel))+".");

  // A--F refitted against Belle 2024 with the line shapes below.  The
  // published values are -0.77, -1.12, -0.59 for C, D, F; they were fitted
  // against that paper's own inline propagators and do not transfer.
  //
  // FITTED BELOW 2.5 GeV ONLY, and that restriction is deliberate.  Eq. (10)
  // cannot describe 1.6-2.5 GeV and the 2.5-3.5 GeV tail at the same time: its
  // cross section falls by a factor 2.7 between 2.5 and 3.5 GeV where the data
  // falls by 5.5, so one F has to serve both and no value does.  F carries
  // 50-80% of the cross section above 1.6 GeV, so a fit that includes the tail
  // lowers F to rescue it and drags 1.6-2.5 down to 0.64-0.80 of the data.
  // Fitting below 2.5 GeV instead gives 0.98 there (chi2/N = 0.77 over
  // 1.05-2.5) and leaves the tail a factor ~2.6 high.  That is the right trade
  // for a generator: the region below 2.5 GeV carries ~20x the cross section,
  // and Belle's tail errors are 50-100% per point.  ABOVE 2.5 GeV THIS MODEL
  // IS AN EXTRAPOLATION AND KNOWN TO BE HIGH.
  //
  // E is poorly constrained -- it contributes 0.1-0.5% of the cross section
  // anywhere and fits to -0.20 +- 0.20.  Holding it at the paper's -0.72
  // instead costs almost nothing (chi2/N 0.92 against 0.83) and returns D and
  // F to within 1% of published, which is a useful consistency check but not a
  // reason to freeze a parameter the data does not support.
  //
  // KHURI-TREIMAN MODE: its own A and B.  The dispersive basis is normalised
  // differently from the rho(770) sum, so Eq. (10)'s couplings put the cross
  // section a factor ~2 high at the omega and right at the phi.  Refitted
  // against the same Belle 2024 data over 0.740-1.050 GeV, with Belle's full
  // statistical and systematic covariance: chi2/ndf 374/122 = 3.07, from 242
  // with the Eq. (10) values.  (Diagonal errors only, as the Eq. (10) fit
  // used, give A = 12.93 and the same B.)  C and D keep their Eq. (10) values
  // -- omega' and omega'' sit far above the window, where their tails are
  // degenerate (correlation -0.999) and nothing separates them.  E and F are
  // unused.  Four things the fit showed:
  //  * BaBar 2021 lies 6-12% below Belle across the window (the known Belle II
  //    / BaBar tension in this channel); fitted to BaBar instead, A and B come
  //    out 4% and 8% smaller.  The choice of data sets the normalisation at
  //    that level.
  //  * Belle's cross section is dressed -- vacuum polarisation included --
  //    and this model uses alpha(0), so A and B absorb |1 - Pi|^-2, a few %.
  //  * The worst region is the omega peak, and the cause is the pole mass: the
  //    particle database has 781.94 MeV, 0.7 below PDG.  The PDG mass alone
  //    takes chi2 from 374 to 261, and moves A by 0.4% and B by 0.2%.
  //  * The fit is to Sherpa's own cross section, not to the Rivet histogram.
  //    BELLE_2024_I2775022 vetoes events with an FSR photon, 2-4% of them here
  //    and ~10% by 3.5 GeV, so its histogram sits below the cross section.
  const double c_eq10[6] = { 18.20, -0.87,  -0.5068, -1.1262, -0.1954, -0.6368 };
  const double c_kt[6]   = { 13.14, -0.895, -0.5068, -1.1262, -0.1954, -0.6368 };
  const double * c = m_kt ? c_kt : c_eq10;
  const string tags[6] = { "A", "B", "C", "D", "E", "F" };
  for (size_t i(0);i<6;++i) {
    m_c[i] = c[i];
    if (p_model) m_c[i] = (*p_model)("EE3Pi_"+tags[i],m_c[i]);
  }

  // Contact term, off unless asked for.  See the header for the convention,
  // the measured values, and the sign trap in the phase.
  m_contact_a     = 0.;
  m_contact_phase = -2.91;
  m_contact_norm  = 7.52;
  m_contact_preserve = true;
  if (p_model) {
    m_contact_a     = (*p_model)("EE3Pi_contact_a",    m_contact_a);
    m_contact_phase = (*p_model)("EE3Pi_contact_phase",m_contact_phase);
    m_contact_norm  = (*p_model)("EE3Pi_contact_norm", m_contact_norm);
    m_contact_preserve =
      ((*p_model)("EE3Pi_contact_preserve_norm",1.)>0.5);
  }
  // The pi is CMD-2's Eq.(42) sign convention; see the header.
  m_contact_c = m_contact_norm * m_contact_a *
                std::polar(1.,m_contact_phase+M_PI);

  m_kt_warned = false;
  // r = 0.03298 m_pi^-2 was fixed on the BESIII omega -> 3 pi Dalitz
  // parameter alpha.
  m_kt_r = 0.03298;
  if (p_model) m_kt_r = (*p_model)("EE3Pi_KT_r",m_kt_r);
  if (m_kt && m_contact_a!=0.)
    THROW(fatal_error,"FF_0_Isoscalar3Pi: the Khuri-Treiman mode has no "
	  "contact term -- its second subtraction supplies the same freedom. "
	  "Set Contact_a to 0.");

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
  //
  // The Khuri-Treiman mode has no rho line shapes at all, so none are made.
  if (!m_kt) {
    const double mpi = Flavour(kf_pi_plus).HadMass();
    p_r770_0  = new GounarisSakuraiM(Flavour(kf_rho_770),mpi);
    p_r770_c  = new GounarisSakuraiM(Flavour(kf_rho_770_plus),mpi);
    p_r1450_0 = new BreitWigner(LineShapes->Get(Flavour(kf_rho_1450)));
    p_r1450_c = new BreitWigner(LineShapes->Get(Flavour(kf_rho_1450_plus)));
    p_r1700_0 = new BreitWigner(LineShapes->Get(Flavour(kf_rho_1700)));
    p_r1700_c = new BreitWigner(LineShapes->Get(Flavour(kf_rho_1700_plus)));
  }

  BuildAmplitude();

  if (m_kt) {
    msg_Out()<<"FF_0_Isoscalar3Pi: Khuri-Treiman (FormFactor 300), r = "
	     <<m_kt_r<<" m_pi^-2, tabulated for sqrt(q^2) = "
	     <<std::round(kKT3Pi_M[0]*kKT3Pi_mpi*1.e3)/1.e3<<" - "
	     <<std::round(kKT3Pi_M[kKT3Pi_NM-1]*kKT3Pi_mpi*1.e3)/1.e3
	     <<" GeV; A-D = ";
    for (size_t i(0);i<4;++i) msg_Out()<<m_c[i]<<(i<3?", ":"; E, F unused\n");
    return;
  }

  if (m_contact_a!=0.)
    msg_Out()<<"FF_0_Isoscalar3Pi: contact term ON, a = "<<m_contact_a
	     <<", phase = "<<m_contact_phase<<" (CMD-2 convention), norm = "
	     <<m_contact_norm<<"\n";
  msg_Out()<<"FF_0_Isoscalar3Pi: hep-ph/0512180 Eq.(10) over registry rho "
	   <<"line shapes; A-F = ";
  for (size_t i(0);i<6;++i) msg_Out()<<m_c[i]<<(i<5?", ":"\n");
}

Sum_Term * FF_0_Isoscalar3Pi::MakeAmplitude(const bool & with_contact) {
  // sub_0 = the pi+ pi- pairing, sub_1 = pi+ pi0, sub_2 = pi- pi0.  The
  // neutral rho takes the first, the charged one the other two -- which is
  // exactly the assignment a single-invariant sum could not express.
  //
  // Built TWICE, with and without the contact constant, so the shape-
  // preserving normalisation has an exact "off" amplitude to integrate
  // against.  The propagators are borrowed by both trees, never owned, so the
  // second tree costs a handful of small term objects and nothing else.
  using ia = invariant_arg;
  if (m_kt) {
    // The isoscalar bracket times ONE dispersive function summed over the
    // three pairings.  The same function serves all three: the solver works
    // in the isospin limit, with one pion mass, so the neutral/charged split
    // the rho line shapes carry has no counterpart here.  with_contact is moot
    // -- FixParameters has already refused a contact term in this mode.
    Sum_Term * kt = new Sum_Term();
    kt->Add(new KT3Pi_Term(ia::sub_0,m_kt_r));
    kt->Add(new KT3Pi_Term(ia::sub_1,m_kt_r));
    kt->Add(new KT3Pi_Term(ia::sub_2,m_kt_r));
    Sum_Term * iso = new Sum_Term();
    for (size_t i(0);i<4;++i) iso->Add(m_iso[i],ia::total,m_c[i]);
    Product_Term * t = new Product_Term(); t->Add(iso); t->Add(kt);
    Sum_Term * amp = new Sum_Term();
    amp->Add(t);
    return amp;
  }
  Sum_Term * h770  = new Sum_Term();
  h770 ->Add(p_r770_0, ia::sub_0); h770 ->Add(p_r770_c, ia::sub_1);
  h770 ->Add(p_r770_c, ia::sub_2);
  // The contact term belongs INSIDE this sum, not added to the amplitude at
  // the end: it interferes with the rho-pi mechanism, and it then inherits the
  // isoscalar q^2 dependence of the A--D bracket, so at the phi it tracks
  // BW_phi exactly as CMD-2's overall normalisation does.
  if (with_contact && m_contact_a!=0.)
    h770->Add(new Constant_Term(m_contact_c));
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

  Sum_Term * amp = new Sum_Term();
  amp->Add(t1); amp->Add(t2); amp->Add(t3);
  return amp;
}

void FF_0_Isoscalar3Pi::BuildAmplitude() {
  p_amp    = MakeAmplitude(true);
  p_amp_nc = MakeAmplitude(false);
}

double FF_0_Isoscalar3Pi::DalitzIntegral(const double & q2,
                                         const bool & with_contact) {
  // Integrate |eps|^2 |F|^2 over the Dalitz plot at fixed q^2, on a grid in
  // (s_+-, s_+0).  No momenta are built: in the three-pion rest frame
  //     E_i    = (q2 + m^2 - s_jk) / (2 sqrt(q2)),   |p_i|^2 = E_i^2 - m^2,
  //     p_+ . p_- = E_+ E_- - (s_+- - 2 m^2)/2,
  // so the epsilon weight |p_+ x p_-|^2 = |p_+|^2 |p_-|^2 - (p_+ . p_-)^2
  // follows from the invariants alone, and its sign marks the boundary.
  const double m2 = sqr(Flavour(kf_pi_plus).HadMass()), rq = sqrt(q2);
  const double lo = 4.*m2, hi = sqr(rq-sqrt(m2));
  if (hi<=lo) return 0.;
  const size_t N = 120;
  const double d = (hi-lo)/double(N);
  double sum = 0.;
  for (size_t i(0);i<N;++i) {
    const double s1 = lo + (double(i)+0.5)*d;
    for (size_t j(0);j<N;++j) {
      const double s2 = lo + (double(j)+0.5)*d;
      const double s3 = q2 + 3.*m2 - s1 - s2;
      if (s3<=lo || s3>=hi) continue;
      const double Ep = (q2+m2-s3)/(2.*rq), Em = (q2+m2-s2)/(2.*rq);
      const double pp2 = Ep*Ep-m2, pm2 = Em*Em-m2;
      if (pp2<=0. || pm2<=0.) continue;
      const double dot = Ep*Em - 0.5*(s1-2.*m2);
      const double eps2 = pp2*pm2 - dot*dot;      // |p_+ x p_-|^2
      if (eps2<=0.) continue;                     // outside the Dalitz plot
      const Invariants k(q2,{s1,s2,s3});
      const Complex F = with_contact ? (*p_amp)(k) : (*p_amp_nc)(k);
      sum += eps2 * std::norm(F) * d * d;
    }
  }
  return sum;
}

double FF_0_Isoscalar3Pi::ContactNorm(const double & q2) {
  if (m_contact_a==0. || !m_contact_preserve) return 1.;
  const double key = std::floor(q2*1e6+0.5)*1e-6;
  std::map<double,double>::const_iterator it = m_norm_cache.find(key);
  if (it!=m_norm_cache.end()) return it->second;
  const double on = DalitzIntegral(q2,true), off = DalitzIntegral(q2,false);
  const double f = (on>0. && off>0.) ? sqrt(off/on) : 1.;
  m_norm_cache[key] = f;
  return f;
}

Complex FF_0_Isoscalar3Pi::operator()(const Vec4D_Vector & moms) {
  if (GetModel()==Model::None) return Complex(0.,0.);
  if (GetModel()==Model::Point) return Complex(1.,0.);
  const Vec4D & pp=moms[m_i_plus], & pm=moms[m_i_minus], & p0=moms[m_i_zero];
  const double q2 = (pp+pm+p0).Abs2();
  if (m_kt) {
    // Above the window the dropped rho(1450)/rho(1700) terms matter, so the
    // answer would be wrong, not merely imprecise: refuse.  Below it the
    // tables are missing but the model is fine, so the amplitude is zero and
    // the run carries on.  See the header.
    const double M(sqrt(q2)/kKT3Pi_mpi);
    if (M>kKT3Pi_M[kKT3Pi_NM-1]*(1.+1.e-9))
      THROW(fatal_error,"FF_0_Isoscalar3Pi: the Khuri-Treiman mode was asked "
	    "for sqrt(q^2) = "+ToString(sqrt(q2))+" GeV, above its tabulated "
	    "window; it drops the rho(1450) and rho(1700) terms, which are not "
	    "negligible there.");
    if (M<kKT3Pi_M[0]*(1.-1.e-9) && !m_kt_warned) {
      msg_Error()<<"FF_0_Isoscalar3Pi: the Khuri-Treiman mode was asked for "
		 <<"sqrt(q^2) = "<<sqrt(q2)<<" GeV, below its tabulated window, "
		 <<"where the amplitude is zero.  Reported once.\n";
      m_kt_warned = true;
    }
  }
  const Invariants k(q2,{ (pp+pm).Abs2(), (pp+p0).Abs2(), (pm+p0).Abs2() });
  // ContactNorm is 1 unless the contact term is on AND shape preservation is
  // asked for, so this is a multiply by a cached double in every other case.
  return (*p_amp)(k) * ContactNorm(q2);
}

Vec4C FF_0_Isoscalar3Pi::Current(const Vec4D_Vector & moms) {
  const Vec4D & pp=moms[m_i_plus], & pm=moms[m_i_minus], & p0=moms[m_i_zero];
  return Vec4C(cross(pp,pm,p0))*(*this)(moms);
}
