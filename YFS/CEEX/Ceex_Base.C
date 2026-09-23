#include "YFS/CEEX/Ceex_Base.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "METOOLS/Main/Spin_Structure.H"
#include "PHASIC++/Process/Process_Base.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "EXTAMP/External_ME_Interface.H"
#include "PHASIC++/Process/External_ME_Args.H"

using namespace YFS;

int Amplitude::s_nlegs = 4;

std::vector<int> Amplitude::s_bit;

void Amplitude::SetLegs(const ATOOLS::Flavour_Vector &flavs)
{
  s_bit.assign(flavs.size(), -1);
  int nb(0);
  for (size_t i(0); i < flavs.size(); ++i) {
    // the same rule METOOLS uses: a massless vector has two states, anything
    // else has 2s+1, so a scalar has one and is not packed at all
    const int ns(flavs[i].IsVector() && !flavs[i].IsMassive()
                 ? 2 : flavs[i].IntSpin() + 1);
    if (ns > 1) s_bit[i] = nb++;
  }
  SetLegs(nb);
}

void Amplitude::SetLegs(int n) {
  if (n < 1 || n > s_maxlegs) {
    msg_Error()<<METHOD<<"(): "<<n<<" legs requested, capacity is "<<s_maxlegs
               <<". Raise Amplitude::s_maxlegs."<<std::endl;
    return;
  }
  s_nlegs = n;
}

Amplitude::Amplitude() {
  // Only the active extent: zeroing the whole capacity would cost 4 KB per
  // construction, and one Amplitude is built per partition.
  const int n(NHel());
  for (int i = 0; i < n; ++i) m_A[i] = Complex(0, 0);
}


Ceex_Base::Ceex_Base(const Flavour_Vector &flavs)
{
  // The amplitude container's extent. 2 -> 2 gives 4 legs and 16 helicity
  // entries, which is what the fixed m_A[2][2][2][2] used to hold.
  // Counts the legs that carry two helicity states, not the legs: a scalar
  // in the final state is not packed, which is what makes the container size
  // agree with Comix's Spin_Amplitudes.
  Amplitude::SetLegs(flavs);

  RegisterDefaults();
  Scoped_Settings s{ Settings::GetMainSettings()["CEEX"] };
  Settings& ss = Settings::GetMainSettings();
  m_onlyz = s["ONLYZ"].Get<int>();
  m_onlyg = s["ONLYG"].Get<int>();
  m_checkxs = s["CHECK_XS"].Get<int>();
  m_comixreal = s["COMIX_REAL"].Get<int>();
  m_comixflip = s["COMIX_REAL_FLIP"].Get<int>();
  m_comixnorm = s["COMIX_REAL_NORM"].Get<double>();
  m_comixphoflip = s["COMIX_REAL_PHOTON_FLIP"].Get<int>();
  m_perphoton    = s["COMIX_REAL_PER_PHOTON"].Get<int>();
  m_comixborn    = s["COMIX_BORN"].Get<int>();
  m_vpon         = s["VIRT_PARTITION_CHECK"].Get<int>() != 0;
  /*
    Any path that asks Comix for the PARTITION Born hands it CEEX's own
    arguments: full-energy beam spinors with the physical pair, and the
    propagator moved to the partition scale. Those momenta do not conserve -
    the photons carry the difference, measured at over 100 GeV in energy - and
    that is not a defect, it is what the partition sum is: the Born as a
    function of a scale decoupled from the kinematics.

    Comix itself is content with that; its recursion carries the propagator as
    a separate factor. What is not content is Amplitude::SetMomenta, which
    checks conservation only under DEBUG__BG and then calls
    ProjectWideMomenta unconditionally. Handed a non-conserving set that
    projection returns NaN, silently - no error, no warning, just a matrix
    element that is not a number.

    So refuse the combination up front rather than let a run produce NaN and
    be debugged twice. With MOMENTUM_PROJECTION off the same call reproduces
    the hand-coded partition Born to 1.7e-8 on every partition.
  */
  if ((s["PARTITION_BORN_CHECK"].Get<int>() != 0 || m_comixborn || m_perphoton)
      && ss["COMIX"]["MOMENTUM_PROJECTION"].SetDefault(true).Get<bool>())
    THROW(fatal_error,
          "CEEX asks Comix for the partition Born, whose arguments do not"
          " conserve momentum by construction. ProjectWideMomenta turns that"
          " into a silent NaN. Set COMIX: {MOMENTUM_PROJECTION: 0}, and note"
          " that it is also Comix's collinear-stability fix, so weigh what"
          " else in the run depends on it.");
  string widthscheme = ss["WIDTH_SCHEME"].Get<string>();
  m_fixedwidth = (widthscheme == "Fixed" || widthscheme == "CMS");
  m_flavs = flavs;
  /*
    CEEX past 2 -> 2 is no longer blocked by the Born. BornAmplitude() used to
    index k[0..3] blindly, which on more legs built the amplitude out of the
    first four - for H l+ l- that is (e-, e+, H, mu-), a scalar in a fermion
    slot. It now takes the outgoing FERMION pair, and the Comix Born can be
    substituted wholesale (CEEX: COMIX_BORN), which for H l+ l- agrees with the
    hand-coded one to 1e-5 at every partition.

    What that leaves is a PER-PROCESS question rather than a structural one.
    The hand-coded spinor algebra is the four-fermion T/U structure, so it is
    right exactly when the Born collapses to that - true for H l+ l-, where
    the Z propagator numerator and the ZZH vertex contract to the 2 -> 2
    current-current form, and NOT something to assume for the next final
    state. CEEX's own virtual is still 2 -> 2 (see CeexOwnVirtual), and the
    Comix REAL covers one photon.

    So DEV_MULTILEG stays a development switch: it says "this process has not
    been shown to reduce to the structure CEEX hand-codes", not "the number is
    meaningless".
  */
  static const int devmultileg(s["DEV_MULTILEG"].SetDefault(0).Get<int>());
  if (flavs.size() != 4) {
    if (!devmultileg)
      THROW(fatal_error, "CEEX is only for 2->2");
    msg_Error()<<METHOD<<"(): CEEX: DEV_MULTILEG is set and this process has "
               <<flavs.size()<<" legs. CEEX's hand-coded amplitude is the"
               <<" four-fermion spinor structure; whether this process reduces"
               <<" to it has to be shown per process, not assumed. Cross-check"
               <<" against CEEX: COMIX_BORN, which is derived rather than"
               <<" hand-coded, before believing the CEEX column. Note also"
               <<" that CEEX's own virtual is 2 -> 2 (use YFS:"
               <<" CEEX_Virtual: external) and the Comix real covers one"
               <<" photon."<<std::endl;
  }
  /*
    The outgoing FERMION pair, by inspection rather than by position. Neutral
    fermions count - e+e- -> nu nubar is a legitimate CEEX process - so the
    test is on being a fermion, not on carrying charge. At 2 -> 2 this
    returns 2 and 3 and nothing downstream changes.
  */
  { size_t n(0);
    for (size_t i(2); i < flavs.size() && n < 2; ++i)
      if (flavs[i].IsFermion()) { (n == 0 ? m_if1 : m_if2) = i; ++n; }
    if (n < 2)
      msg_Error()<<METHOD<<"(): no outgoing fermion pair found among "
                 <<flavs.size()<<" legs; the electroweak couplings will be "
                 <<"those of legs 2 and 3, which is almost certainly wrong."
                 <<std::endl;
  }

  /*
    CEEX's own virtual is a 2 -> 2 object - the vertex and box functions in
    Ceex_Virtual.C are analytic expressions for e+e- -> f fbar. Used at any
    other multiplicity it does not fail, it returns a number: measured on
    e+e- -> H mu+ mu- it gave <rho1/rho0 - 1> = +22, a 2200% "correction",
    while every other column of the same run was correct to a few percent.
    That is the failure this refuses.
  */
  if (m_flavs.size() != 4 && m_useceex && m_ceexvirtsrc == ceexvirt::ceex)
    THROW(fatal_error,
          "CEEX's own virtual is 2 -> 2 only, and this process has "
          + ATOOLS::ToString(m_flavs.size()) + " legs. Set YFS: CEEX_Virtual:"
          " external to take the helicity-summed virtual from the loop"
          " provider instead (which needs V in NLO_Part).");

  if (flavs[m_if1].IsNeutrino() && flavs[m_if2].IsNeutrino()) {
    m_onlyz = true;
  }

  m_Q1Q2I = flavs[0].Charge() * flavs[1].Charge();
  m_QIQF  = flavs[0].Charge() * flavs[m_if1].Charge();
  // Bhabha: the final pair is the beam pair, so a t-channel gamma/Z is exchanged
  // between the two fermion lines on top of the s-channel annihilation.
  m_bhabha = (flavs[0] == flavs[m_if1]);
  m_Q1Q2F = flavs[m_if1].Charge() * flavs[m_if2].Charge();
  m_MZ = Flavour(kf_Z).Mass();
  m_gZ = Flavour(kf_Z).Width();
  double mw = Flavour(kf_Wplus).Mass();
  double MH = Flavour(kf_h0).Mass();
  double  GH  = Flavour(kf_h0).Width();
  double  GW  = Flavour(kf_Wplus).Width();
  double  GZ  = Flavour(kf_Z).Width();
  m_I   = Complex(0., 1.);

  double F_L = 0.;
  double F_R = 0.;

  m_sin2tw = MODEL::s_model->ComplexConstant("csin2_thetaW");
  if (Settings::GetMainSettings()["CEEX"]["REAL_SIN2THETAW"]
      .SetDefault(0).Get<int>() != 0)
    m_sin2tw = Complex(m_sin2tw.real(), 0.);
  m_sW = m_sin2tw;
  m_e = sqrt(4.*M_PI * m_alpha);
  m_cW = 1. - m_sW;
  m_norm = sqrt(16. * m_sW * (1. - m_sW));
  m_qe       = m_flavs[0].Charge();
  m_qf       = m_flavs[m_if1].Charge();
  m_Q1Q2I = m_flavs[0].Charge() * m_flavs[1].Charge();
  m_ae       = 2.*m_flavs[0].IsoWeak();
  m_af       = 2.*m_flavs[m_if1].IsoWeak();
  // Keep 2*T3 and 4*Q*sw^2 separately: the electroweak kappa factors multiply
  // only the sin^2 piece (GPS_EWFFact), so the two cannot be pre-combined.
  m_t3e2 = m_ae;
  m_t3f2 = m_af;
  m_qesw = 4.*m_qe * m_sin2tw;
  m_qfsw = 4.*m_qf * m_sin2tw;
  m_ve       = (m_ae - 4.*m_qe * m_sin2tw) / m_norm;
  m_vf       = (m_af - 4.*m_qf * m_sin2tw) / m_norm;
  m_ae /= m_norm;
  m_af /= m_norm;
  m_weak = s["WEAK"].Get<int>();
  m_mass_I = flavs[0].Mass();
  m_mass_F = flavs[m_if1].Mass();
  // full EW couplings
  m_I_L = -m_I * sqrt(4 * M_PI * m_alpha) / (2.*m_sW * m_cW) * (2.*flavs[0].IsoWeak()
          - 2.*flavs[0].Charge() * m_sW * m_sW);

  m_I_R = -m_I * sqrt(4 * M_PI * m_alpha) / (2.*m_sW * m_cW) * (-2.*flavs[0].Charge() * m_sW * m_sW);

  m_F_L = -m_I * sqrt(4 * M_PI * m_alpha) / (2.*m_sW * m_cW) * (2.*flavs[m_if1].IsoWeak()
          - 2.*flavs[m_if1].Charge() * m_sW * m_sW);

  m_F_R = -m_I * sqrt(4 * M_PI * m_alpha) / (2.*m_sW * m_cW) * (-2.*flavs[m_if1].Charge() * m_sW * m_sW);
  m_cL = -m_I * sqrt(4 * M_PI * m_alpha) / (2.*m_sW * m_cW) * (2.*m_flavs[1].IsoWeak()
         - 2.*m_flavs[1].Charge() * m_sW * m_sW) / m_norm;
  m_cR = -m_I * sqrt(4 * M_PI * m_alpha) / (2.*m_sW * m_cW) * (-2.*m_flavs[1].Charge() * m_sW * m_sW) / m_norm;
  m_zeta = {1, 1, 0, 0};
  m_eta = {0, 0, 1, 0};
  m_b  = {0.0,  0.8723e0, -0.7683e0, 0.3348e0};
}


void Ceex_Base::RegisterDefaults()
{
  Scoped_Settings s{ Settings::GetMainSettings()["CEEX"] };
  s["ONLYZ"].SetDefault(0);
  s["ONLYG"].SetDefault(0);
  s["CHECK_XS"].SetDefault(0);
  s["WEAK"].SetDefault(1);
  // Take the O(alpha) real (beta_1) from Comix's helicity amplitudes instead
  // of the hand-coded spinor algebra. Off by default.
  s["COMIX_REAL"].SetDefault(0);
  /*
    Both of these used to be fitted numbers (26 and 0.5). They are now DERIVED
    at the Born, by Ceex_Base::CalibrateComixMap, which is why the defaults
    are sentinels rather than values:

      COMIX_REAL_FLIP  < 0  derive the fermion bits from the Born and take
                            the photon bit from COMIX_REAL_PHOTON_FLIP
      COMIX_REAL_NORM <= 0  derive from sum|A_comix|^2 / sum|e^2 A_hand|^2

    A positive value overrides the derivation, which is how a disagreement
    with the calibration gets investigated rather than papered over.
  */
  s["COMIX_REAL_FLIP"].SetDefault(-1);
  s["COMIX_REAL_PHOTON_FLIP"].SetDefault(1);
  s["COMIX_REAL_NORM"].SetDefault(-1.);
  /*
    Let the Comix amplitude be used above one photon. OFF: the exact n-photon
    amplitude is not infrared subtracted and carries soft content that exp(Y)
    and the crude S-factors already hold, so it is a different matching scheme
    rather than an extension of this one. Kept reachable because the machinery
    is in place and the comparison is worth being able to make.
  */
  s["COMIX_REAL_MULTIPHOTON"].SetDefault(0);
  /*
    Take the one-photon real from Comix once per PHOTON, at every
    multiplicity, instead of once per event with every photon attached. This
    is the structure CEEX actually has - beta_1 is a sum over photons - and it
    needs only the one-photon process, which always exists.
  */
  s["COMIX_REAL_PER_PHOTON"].SetDefault(0);
  // How many events the soft probe walks down the lambda ladder. Each rung
  // costs a Comix evaluation and perturbs the random sequence, so it is a
  // diagnostic budget, not something to leave large.
  s["SOFT_PROBE_EVENTS"].SetDefault(10);
  // @@@ PARTBORN: the partition Born against the Born at the configuration
  // that realises the partition's scale.
  s["PARTITION_BORN_CHECK"].SetDefault(0);
  /*
    Take the partition Born from Comix instead of the hand-coded spinor
    algebra. That algebra is the 2 -> 2 specific part of CEEX, so this is the
    step that makes the partition sum process independent. Requires
    COMIX: MOMENTUM_PROJECTION: 0 - see the check in the constructor.
  */
  s["COMIX_BORN"].SetDefault(0);
  // @@@ BORNALIGN: is the Comix -> CEEX Born alignment scale independent?
  s["BORN_ALIGN_CHECK"].SetDefault(0);
  // @@@ BETA1: the Comix hard remainder against the hand-coded one, vs E_gamma
  s["BETA1_CHECK"].SetDefault(0);
  // Hand Comix CEEX's own (non-conserving) beta_1 arguments rather than a
  // mapped conserving configuration. Needs COMIX: MOMENTUM_PROJECTION: 0.
  s["COMIX_REAL_CEEX_ARGS"].SetDefault(0);
  // Write the |beta_1|/beta_0 vs E_gamma scan (the real-validation figure).
  s["BETA1_SCAN"].SetDefault(0);
  // @@@ VIRTPART: does the virtual factor V(h) depend on the partition?
  s["VIRT_PARTITION_CHECK"].SetDefault(0);
  /*
    Diagnostics. All off by default, each costing one branch on a cached
    static once the run is going. They live here rather than in the
    environment so that a run is fully specified by its YAML card: an
    environment variable does not appear in the card, is not echoed in the
    run summary and is not carried to a batch node with the job, so a
    diagnostic run could not be reproduced from what was kept of it.
  */
  s["BORNNORM_CHECK"].SetDefault(0);   // @@@ BORNANG, Born normalisation
  s["PIN_PHOTON_HEL"].SetDefault(0);   // pin photon helicity (+1/-1) for sums
  s["DUMP_XMIN"].SetDefault(0.0);      // min x_gamma for the point dump
  s["DUMP_NPHOT"].SetDefault(1);       // photon multiplicity to dump at
  s["COMIX_CHECK"].SetDefault(0);      // @@@ CEEXCX, hand-coded vs Comix
  s["GOLDEN"].SetDefault(0);           // @@@ CEEXGOLD, regression stream
  s["WEIGHT_PROBE"].SetDefault(0);     // @@@ CEEXWT, what makes a heavy event
}



void Ceex_Base::Init(const Vec4D_Vector &p)
{
  m_momenta = p;
  m_cms = Poincare(m_momenta[0] + m_momenta[1]);
  Poincare cms = m_cms;
  Poincare Rot = Poincare(Vec4D(0., 0., 0., 1.));
  for (size_t i(0); i < p.size(); ++i) {
    cms.Boost(m_momenta[i]);
    // cms.Boost(m_bornmomenta[i]);
  }
  m_crude = 2.0 / (4.0 * M_PI);
  for (size_t i(0); i < m_isrphotons.size(); ++i) {
    m_crude /= pow(2 * M_PI, 3);
  }
  if (m_momenta.size() >= 6) {
    m_sp = (m_momenta[4] + m_momenta[5]).Abs2();
    m_sQ = m_sp;
  } else {
    m_sp = (m_momenta[2] + m_momenta[3]).Abs2();
    m_sQ = m_sp;
  }
  m_T = 0;
}



void Ceex_Base::MakeProp()
{
  if (m_fixedwidth) {
    m_propZ =  1. / Complex(m_sp - sqr(m_MZ), m_gZ * m_MZ);
  }
  else {
    m_propZ =   1. / Complex(m_sp - sqr(m_MZ), m_gZ * m_sp / m_MZ);
  }
  m_propG = 1. / Complex(m_sp,0);
  if (m_onlyz)  m_propG = 0;
  if (m_onlyg)  m_propZ = 0;
  m_prop =  m_propZ + m_propG;
}


void Ceex_Base::MakePropT(const Vec4D_Vector &p)
{
  if (!m_bhabha) {
    m_propGt = m_propZt = Complex(0., 0.);
    return;
  }
  m_tinv = (p[0] - p[2]).Abs2();
  m_propZt = 1. / Complex(m_tinv - sqr(m_MZ), m_gZ * m_MZ);
  m_propGt = 1. / Complex(m_tinv, 0.);
  if (m_onlyz) m_propGt = 0.;
  if (m_onlyg) m_propZt = 0.;
}


void Ceex_Base::MakeEWFF(double svar, double costhd)
{
  m_kapE = m_kapF = m_kapEF = Complex(1., 0.);
  m_rhoEW = m_gamVPi = Complex(1., 0.);
  m_vvcor = Complex(1., 0.);
  if (!m_weak) {           // tree couplings == KKMC with KeyElw = 0
    m_ve = (m_t3e2 - m_qesw) / m_norm;
    m_vf = (m_t3f2 - m_qfsw) / m_norm;
    return;
  }
  // --- weak form factors go here; not yet implemented ---
  m_ve = (m_t3e2 - m_qesw * m_kapE) / m_norm;
  m_vf = (m_t3f2 - m_qfsw * m_kapF) / m_norm;
  // Angle-dependent double-vector correction; kapEF carries the box content.
  const Complex vvcef((m_t3e2*m_t3f2
                       - m_qesw*m_t3f2*m_kapE
                       - m_qfsw*m_t3e2*m_kapF
                       + m_qesw*m_qfsw*m_kapEF) / (m_norm*m_norm));
  m_vvcor = (std::abs(m_ve*m_vf) > 0.) ? vvcef/(m_ve*m_vf) : Complex(1., 0.);
}


Complex Ceex_Base::CouplingZ(double  j, int mode) {
  if (m_onlyg) return 0.;
  Complex zcpl;
  if (mode == 1) {
    zcpl = m_ve * m_vf * m_vvcor - dcmplx(j) * m_ae * m_vf + dcmplx(j) * m_ve * m_af - m_ae * m_af;
  }
  else if (mode == 0) {
    zcpl = m_ve * m_vf * m_vvcor - dcmplx(j) * m_ae * m_vf - dcmplx(j) * m_ve * m_af + m_af * m_ae;
  }
  else msg_Error() << METHOD << "\n wrong mode\n";

  if (zcpl == 0.) {
    msg_Error() << "Z coupling is Zero!\n";
  }
  return zcpl;
}



Complex Ceex_Base::CouplingG() {
  m_gcpl = Complex(m_QIQF, 0);
  return m_gcpl;
}


void Ceex_Base::BuildCeexMomenta()
{
  /*
    m_momenta is {beam, beam, Born final state, lab final state}, MakeCEEX
    appending the lab set. CEEX wants the beams and the LAB legs, those being
    the ones that balance against the photons; the Born set is the fallback
    when MakeCEEX has not appended yet.

    Generalised from the fixed pair to nf = m_flavs.size() - 2 final legs. At
    nf = 2 the offset is 4 or 2 exactly as before.
  */
  const size_t nf(m_flavs.size() >= 2 ? m_flavs.size() - 2 : 0);
  m_pceex.clear();
  m_pceex.push_back(m_momenta[0]);
  m_pceex.push_back(m_momenta[1]);
  const bool havelab(m_momenta.size() >= 2 + 2*nf);
  const size_t off(havelab ? 2 + nf : 2);
  for (size_t i(0); i < nf && off + i < m_momenta.size(); ++i)
    m_pceex.push_back(m_momenta[off + i]);
}


void Ceex_Base::ZerAmplit() {
  // Every helicity entry the container actually holds, not the 16 a 2 -> 2
  // final state happens to need.
  const int nh(Amplitude::NHel());
  for (int f = 0; f < nh; ++f) {
    m_AmpExpo0.m_A[f]    = Complex(0., 0.);
    m_AmpExpo1.m_A[f]    = Complex(0., 0.);
    m_AmpBornVirt.m_A[f] = Complex(0., 0.);
    m_AmpBornReal.m_A[f] = Complex(0., 0.);
    m_snapBorn.m_A[f]    = Complex(0., 0.);
    m_snapVirt.m_A[f]    = Complex(0., 0.);
    m_snapReal.m_A[f]    = Complex(0., 0.);
  }
}



void Ceex_Base::MakeRho() {
  double sum0(0.), sum1(0.);
  const int nh(Amplitude::NHel());
  for (int f = 0; f < nh; ++f) {
    sum0 += std::real(m_AmpExpo0.m_A[f] * conj(m_AmpExpo0.m_A[f]));
    sum1 += std::real(m_AmpExpo1.m_A[f] * conj(m_AmpExpo1.m_A[f]));
  }
  // Average over the four initial-state helicity configurations.
  m_result0 = sum0 / 4.;
  m_result  = sum1 / 4.;
  double sumbv(0.), sumbr(0.);
  for (int j1 = 0; j1 <= 1; ++j1)
    for (int j2 = 0; j2 <= 1; ++j2)
      for (int j3 = 0; j3 <= 1; ++j3)
        for (int j4 = 0; j4 <= 1; ++j4) {
          sumbv += std::real(m_AmpBornVirt.m_A[Idx(j1,j2,j3,j4)]
                             * conj(m_AmpBornVirt.m_A[Idx(j1,j2,j3,j4)]));
          sumbr += std::real(m_AmpBornReal.m_A[Idx(j1,j2,j3,j4)]
                             * conj(m_AmpBornReal.m_A[Idx(j1,j2,j3,j4)]));
        }
  m_resultbv = sumbv / 4.;
  m_resultbr = sumbr / 4.;
  m_rho0sum += m_result0;
  m_rho1sum += m_result;
  m_rhobvsum += m_resultbv;
  m_rhobrsum += m_resultbr;
  m_rhocrudsum += m_rhocrud;
  ++m_rhon;
}


void Ceex_Base::Reset() {
  m_result = 0;
}


double Ceex_Base::Xi(const Vec4D p, const Vec4D q) {
  return sqrt((m_zeta * p) / (q * m_zeta));
}

double Ceex_Base::RealFactorPhoton(size_t j) const
{
  if (j >= m_realphot.size() || m_result0 <= 0.) return 0.;
  double sum(0.);
  for (int a = 0; a <= 1; ++a)
    for (int b = 0; b <= 1; ++b)
      for (int c = 0; c <= 1; ++c)
        for (int d = 0; d <= 1; ++d) {
          const Complex z(m_AmpExpo0.m_A[Idx(a,b,c,d)] + m_realphot[j].m_A[Idx(a,b,c,d)]);
          sum += std::real(z * conj(z));
        }
  return sum/4./m_result0 - 1.;
}
