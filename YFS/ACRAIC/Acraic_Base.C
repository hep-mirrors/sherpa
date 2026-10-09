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
  m_onlyz = s["ONLYZ"].Get<bool>();
  m_onlyg = s["ONLYG"].Get<bool>();
  m_checkxs = s["CHECK_XS"].Get<int>();
  m_comixreal = s["COMIX_REAL"].Get<bool>();
  m_comixflip = s["COMIX_REAL_FLIP"].Get<int>();
  m_comixnorm = s["COMIX_REAL_NORM"].Get<double>();
  m_comixphoflip = s["COMIX_REAL_PHOTON_FLIP"].Get<int>();
  m_perphoton    = s["COMIX_REAL_PER_PHOTON"].Get<bool>();
  m_comixborn    = s["COMIX_BORN"].Get<bool>();
  m_wstages      = s["W_STAGES"].Get<tristate::code>();
  m_weikonal     = s["W_EIKONAL"].Get<weikonal::code>();
  m_crudegen     = s["CRUDE_FROM_GENERATOR"].Get<crudegen::code>();
  m_vpon         = s["VIRT_PARTITION_CHECK"].Get<int>() != 0;
  m_order        = s["ORDER"].Get<int>();
  if (m_order != 1 && m_order != 2) {
    msg_Error()<<METHOD<<"(): CEEX: ORDER = "<<m_order<<" is not 1 or 2; "
               <<"using 1."<<std::endl;
    m_order = 1;
  }
  if (m_order == 2)
    msg_Info()<<"CEEX: ORDER 2 - tree-level double real beta_2 from Comix "
              <<"two-photon amplitudes, pairs with x = 2E/sqrt(s) > "
              <<s["BETA2_XCUT"].Get<double>()<<"."<<std::endl;
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

    This used to be REFUSED up front, with the user told to set
    COMIX: MOMENTUM_PROJECTION: 0 by hand. Set it here instead: a run that asks
    for the Comix Born should not also have to know what that implies for
    Comix. With the projection off the same call reproduces the hand-coded
    partition Born to 1.7e-8 on every partition.

    It is still a GLOBAL switch, and it is also Comix's collinear-stability
    fix, so this degrades every other Comix evaluation in the run. Scoping it
    to the CEEX calls is the right answer and needs a hook PHASIC++ owns and
    COMIX fulfils: Real_Correction holds a Process_Base, so a COMIX-side
    override is not reached for the real, and the real is where it matters most
    (measured on H l+ l-: +101% unscoped against +0.29%). Until that exists,
    say so loudly rather than let the setting be forgotten.
  */
  /*
    COMIX: MOMENTUM_PROJECTION used to be switched off GLOBALLY here, because
    CEEX hands Comix legs that do not balance. That also removed Comix's
    collinear-stability fix from every other evaluation in the run, and for
    e+e- -> mu mu tau tau the fixed-order real came out with an 87% error
    from outliers. The projection is now suspended only inside CEEX's own
    calls (COMIX::Amplitude::Scoped_No_Projection on every CEEX entry point
    of Single_Process), so the setting is left alone.
  */
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
  static const bool devmultileg(s["DEV_MULTILEG"].SetDefault(false).Get<bool>());
  /*
    The 2 -> 2 special path (hand-coded Born, its Comix alignment and map,
    CEEX's own virtual, the fermion-pair m_if1/m_if2) is a description of
    e+e- -> f fbar. It used to be selected by the leg count alone, so
    e+e- -> gamma gamma - four legs, no outgoing fermion - took m_if1/m_if2
    as the BEAMS, calibrated the Comix map against a fermion-pair Born that
    does not exist for it, failed ("Born calibration failed"), and dropped
    every beta_1. The selection is the final state's structure, not its size.
  */
  m_ffbar = (flavs.size() == 4 && flavs[2].IsFermion()
             && flavs[3] == flavs[2].Bar());
  if (!m_ffbar) {
    /*
      The throw is about the BORN, so it applies only when the hand-coded Born
      is the one selected. That used to be unconditional, which was right while
      the hand-coded Born was the default - it indexes a four-fermion T/U
      structure and on more legs quietly builds the wrong amplitude. Now that
      COMIX_BORN is on by default the Born is derived per process and general,
      so the binding reason is gone and only the warning below is owed.

      DEV_MULTILEG remains the override for the hand-coded case.
    */
    if (!m_comixborn && !devmultileg)
      THROW(fatal_error, "CEEX past 2 -> 2 needs the Comix Born. Either leave"
            " CEEX: COMIX_BORN at its default, or set CEEX: DEV_MULTILEG to"
            " use the hand-coded 2 -> 2 spinor structure anyway.");
    msg_Error()<<METHOD<<"(): this process has "<<flavs.size()<<" legs."
               <<" The Born is taken from Comix and is general, but CEEX's own"
               <<" virtual is 2 -> 2 (use YFS: CEEX_Virtual: external) and the"
               <<" Comix real covers one photon, above which the hand-coded"
               <<" beta_1 runs. Whether this process reduces to CEEX's"
               <<" four-fermion structure has to be shown, not assumed."
               <<std::endl;
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
    if (n < 2) {
      // no fermion pair (gamma gamma): legs 2 and 3. Only the hand-coded
      // 2 -> 2 couplings read these, and that path is not taken here.
      m_if1 = 2; m_if2 = 3;
      msg_Debugging()<<METHOD<<"(): no outgoing fermion pair among "
                     <<flavs.size()<<" legs; m_if1/m_if2 = 2/3, general "
                     <<"Comix branch."<<std::endl;
    }
  }

  /*
    CEEX's own virtual is a 2 -> 2 object - the vertex and box functions in
    Ceex_Virtual.C are analytic expressions for e+e- -> f fbar. Used at any
    other multiplicity it does not fail, it returns a number: measured on
    e+e- -> H mu+ mu- it gave <rho1/rho0 - 1> = +22, a 2200% "correction",
    while every other column of the same run was correct to a few percent.
    That is the failure this refuses.
  */
  if (!m_ffbar && m_useceex && m_ceexvirtsrc == ceexvirt::ceex)
    THROW(fatal_error,
          "CEEX's own virtual is e+e- -> f fbar only, and this process has "
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
      .SetDefault(false).Get<bool>())
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
  m_weak = s["WEAK"].Get<bool>();
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
  s["ONLYZ"].SetDefault(false);
  s["ONLYG"].SetDefault(false);
  s["CHECK_XS"].SetDefault(0);
  // false: tree-level couplings in CEEX's hand-coded Born; the weak virtual
  // comes from the loop provider (CEEX_Virtual: auto/external)
  s["WEAK"].SetDefault(false);
  s["COMIX_REAL"].SetDefault(true);
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
  s["COMIX_REAL_MULTIPHOTON"].SetDefault(false);
  /*
    Take the one-photon real from Comix once per PHOTON, at every
    multiplicity, instead of once per event with every photon attached. This
    is the structure CEEX actually has - beta_1 is a sum over photons - and it
    needs only the one-photon process, which always exists.
  */
  s["COMIX_REAL_PER_PHOTON"].SetDefault(false);
  s["BETA1_PARTITION_CHECK"].SetDefault(0);  // @@@ B1PART
  s["BETA1_CLOSURE"].SetDefault(0);          // @@@ B1CLOSE
  s["CLOSURE_1PHOT"].SetDefault(0);          // @@@ CLOS1
  s["BETA1_XCUT"].SetDefault(1e-3);          // beta_1 soft threshold
  s["BETA1_TRACE"].SetDefault(0);            // @@@ B1TRACE on the dump event
  /*
    Heavy-event trace: > 0 runs the beta_1 trace on EVERY event into a buffer
    that YFS_Handler prints only when rho_1/rho_crude exceeds this value,
    together with the event's legs and photons (@@@ HEAVY). Costly (one extra
    Born and M_1 per photon per partition); diagnostics only.
  */
  s["TRACE_FACTOR_ABOVE"].SetDefault(0.);
  /*
    x_gamma = 2E/sqrt(s) below which a photon is not enumerated in the
    partition sum (see m_fixedstage). 0 enumerates every photon (KKMC).
  */
  s["SOFT_PARTITION_CUT"].SetDefault(1e-3);
  /*
    MOMENTUM_REPAIR (true): project the final legs and all photons onto exact
    momentum conservation against the beams and onto their mass shells before
    anything is evaluated on them; false = use the event's momenta as handed
    over. See Ceex_Base::RepairMomentumBalance.
  */
  s["MOMENTUM_REPAIR"].SetDefault(true);
  /*
    Legs for beta_1's M_1 (name or old integer):
      physical (0, default): the physical legs with the s-channel propagator
        momenta shifted by the other photons per side (KKMC's construction);
      balanced (1): a rebuilt balanced point (PartitionLegs).
    See ComixInfraredSubtracted_1_0.
  */
  s["BETA1_LEGS"].SetDefault(beta1legs::physical);
  /*
    Space-like exchange lines (Bhabha's t-channel boson, the t/u-channel
    electrons of e+e- -> gamma gamma, a t-channel neutrino: any current with
    one initial leg and part of the final state, detected by leg content in
    COMIX::Amplitude::SetPropShifts) of the partition Born and of M_1 at the
    partition's REDUCED invariant rather than at the unreduced one Comix's
    root-0 recursion leaves them at. on (1), off (0), auto (-1, default): on
    for 2 -> 2 only (ExchangeLineShiftsOn). See
    Ceex_Base::AddExchangeLineShifts.
  */
  s["TCHANNEL_SHIFT"].SetDefault(tristate::automatic);
  /*
    W stages for e+e- -> W+W- -> 4f (NOTES-w-stages-2026-09-27.md). The
    charged final legs of a WW-type final state are split into two DECAY
    stages, one per W, each with the W as an incoming leg and its charged
    daughter as the outgoing one, and the W's join the beams as outgoing legs
    of the production stage. Each stage is charge neutral; a photon on a
    decay stage shifts only that W's propagator (COMIX::Amplitude::
    SetPropShifts with the W's daughter mask), a photon on the production
    stage the s-channel line. Without it the two charged leptons of such a
    state form one flat stage whose photons shift BOTH W lines at once
    (SetPropShifts's "partly contained, root side" rule), which dragged an
    off-shell W onto its pole: cc_em_mup at 161 GeV, CEEX factors up to
    4777 with TCHANNEL_SHIFT off, rho_0/rho_crude 267 against the 2^n bound.
    Tristate, name or old integer:
      off (0): every process as before, bit for bit.
      on (1, default since 2026-09-28): whenever DipoleSet::FindWW recognises
        the final state; a process without W's is unchanged bit for bit.
        cc_mum_taup CEEX/NLO 0.989 -> 1.000 +- 0.02, cc_em_mup 0.93 -> 0.975
        +- 0.04 (with CRUDE_FROM_GENERATOR), NOTES-w-stages-2026-09-27.md 8.2.
      auto (-1): on when in addition both W's are within YFS:
        CLUSTERING_THRESHOLD widths of the pole (the pole scheme's own
        window).
    Needs the handler to hand over the W groups (YFS_Handler::CEEXStageGroups).
  */
  s["W_STAGES"].SetDefault(tristate::on);
  /*
    The W momentum a decay- or production-stage eikonal uses (see
    StageLegMomentum; name or old integer). daughters (0, default): the W at
    its daughters, a per-event quantity, so
    the soft-factor table is computed once; the soft limit of Comix's M_1
    is then missed by O(K_decay/M_W) on the W term when other hard photons
    sit on that decay stage. partition (1): the pole momentum follows the
    partition
    (daughters + that partition's other decay photons), which is what
    Comix's shifted propagator carries; the resonance stages' soft factors
    are recomputed per partition. Cost negligible against the Comix calls.
  */
  s["W_EIKONAL"].SetDefault(weikonal::daughters);
  /*
    The crude the CEEX weight divides by, built on the GENERATOR's stages
    (the flat radiating-dipole groups) rather than on CEEX's own; see
    Ceex_Base::CrudeFromGenerator. With CEEX's stages equal to the
    generator's dipoles - every process without W stages - it is the same
    number as the per-partition crude, which is the gate it has to pass
    (compare prints both per event as @@@ CRUDEGEN). Name or old integer:
    off (0); on (1), use it; compare (2), compare only; auto (-1, default),
    use it exactly when W stages are active - the case where CEEX's stages
    and the generator's dipoles differ - and leave every other process
    untouched.
  */
  s["CRUDE_FROM_GENERATOR"].SetDefault(crudegen::automatic);
  /*
    beta_0(X_wp) on the REAL phase-space point whose invariant is X_wp^2 -
    beams at X_wp, radiating pair at X_wp minus the spectators, per partition
    - with natural propagators and no pseudo-flux, instead of KKMC's physical
    spinors with the pole pinned to m_sp times svarY/svarQ. See
    InfraredSubtractedME_0_0 for the two forms and what the earlier,
    single-point version of this switch got wrong.
  */
  s["BORN_AT_SPRIME"].SetDefault(false);
  /*
    Born processes with a space-like exchange line (a current holding one
    initial leg and part of the final state: the t/u-channel electrons of
    e+e- -> gamma gamma, Bhabha's t-channel boson, a t-channel neutrino).
    For those the physical-spinor partition Born with only its poles moved
    is not the reduced-point Born times a flux, as it is for an s-channel
    Born: the numerators do not scale with the poles. Tristate, name or old
    integer:
      on (1): beta_0 is the Born at the partition's REAL reduced point
         (the BORN_AT_SPRIME form, BornLegsAt(X_wp)), the crude that Born
         times s/X_wp^2 (the generator's density), and every photon's M_1
         in beta_1 is evaluated on the balanced point with the other
         photons taken out of the beams (the BETA1_LEGS: 1 form,
         PartitionLegs) - genuine amplitudes throughout, so the soft
         cancellations between beta_0 and beta_1 hold between like objects.
      off (0): the physical-spinor forms, as for s-channel Borns.
      auto (-1, default): on exactly when the Born has such a line, detected once
         per Born process (Ceex_Base::BornHasExchangeLine), AND there is a
         single radiating stage (no charged final state, e.g. gamma gamma).
         Every s-channel process keeps the old numbers bit for bit. With a
         final-state stage (Bhabha) the one-photon closure fails with 1:
         +10.6% on the Bhabha CEEX column; see Ceex_Partitions.C.
    Measured on e+e- -> gamma gamma at the Z pole (pT > 1 GeV, 20k events,
    2026-09-26): with 0 the CEEX factor at one photon deviates from the exact
    |M_1|^2/density by factors 0.04-1.7 at wide angle and the column is
    6829 +- 68% pb from multi-photon events with weights up to 3e5 x the
    crude; with 1 the one-photon factor equals the exact one to 1e-4 at every
    x and the column is 158.4 +- 1.1% pb, largest single-event share 0.6%.
  */
  s["TCHANNEL_REDUCED_BORN"].SetDefault(tristate::automatic);
  /*
    CEEX: TCHANNEL_MULTIPHOTON - beta_1 beyond one photon when beta_0 is the
    reduced t-channel Born (TCHANNEL_REDUCED_BORN active: gamma gamma).
    Nothing changes at one photon, nor for any process without the reduced
    Born (every s-channel Born, Bhabha). Name or old integer:
      partition_legs (0, default): M_1 on PartitionLegs, subtraction with
        the eikonal on those legs.
      one_photon_point (1): M_1 at the photon's one-photon point (OnePhotonScaledLegs, YFS.NLO's
        REAL_MAP 2 point carried into the generator's frame), as the ratio
        M_1/s(point) times the physical eikonal, so the subtraction is
        exactly the Born term's s_phys B_0.
      factorised (2): one_photon_point, and the partition's beta_1 terms
        combined in factorised form
        along the Born helicity vector (AddFactorisedRemainder): identical
        at O(alpha^1), with the factorised beta_2 and higher added.
    Why (e+e- -> gamma gamma, Z pole, 2026-09-28): with a hard wide-angle ISR
    photon and soft companions, PartitionLegs tilts the companions' beams
    by the hard photon's pT - their eikonal there was 1e-5..18 times the
    physical one - and in events whose generator Born sits on the reduced
    frame's t-channel pole (m_born 1e3-1e7) that left |A_1| 5-18 times
    |A_0| where the exact ME is 1e-4 of it: ACRAIC 1.65 x YFS.NLO in the
    photon-tagged region. NOTES-aa-multiisr-2026-09-28.md.
  */
  s["TCHANNEL_MULTIPHOTON"].SetDefault(tchmultiphoton::partition_legs);
  s["BETA1_BORNLEGS"].SetDefault(1);         // 1 = reduced, 0 = physical
  /*
    The pseudo-flux svarY/svarQ on beta_0 (name or old integer):
    rho0_and_rho1 (0) = in rho_0 and rho_1 (KKMC's beta_0 without KKMC's
    compensation), neither (1) = in neither, rho0_only (2) = in rho_0 only.
    KKMC's O(alpha^1) amplitude is flux-free (its (1-CKine) terms cancel the
    flux per final-state photon up to 2k_i.k_j/Q^2) while its RhoExp0 keeps
    it, so rho0_only is KKMC's own rho_1/rho_0: -0.16525 against KKMC's
    -0.16557 on the seed-11 n=2 point, and the Z-pole cross section within
    1.1% of YFS.NLO (rho0_and_rho1: -2.8%, neither: +24%). Default rho0_only.
  */
  s["NO_PSEUDOFLUX"].SetDefault(pseudoflux::rho0_only);
  s["SOFT_NORM_CHECK"].SetDefault(0);        // @@@ SOFTNORM
  s["SOFT_LIMIT_TEST"].SetDefault(0);        // @@@ SOFTLIM
  s["BETA1_COMPARE"].SetDefault(0);          // @@@ B1CMP
  s["SFAC_CALIB"].SetDefault(0);             // @@@ SFACCAL
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
  /*
    ON by default: the Comix Born is the DERIVED object and the hand-coded
    2 -> 2 spinor algebra is the fallback. At 2 -> 2 it reproduces the
    hand-coded CEEX weight to round-off; past 2 -> 2 the hand-coded Born is not
    the right amplitude at all. COMIX_BORN: false selects the hand-coded one.
  */
  s["COMIX_BORN"].SetDefault(true);
  // @@@ BORNALIGN: is the Comix -> CEEX Born alignment scale independent?
  s["BORN_ALIGN_CHECK"].SetDefault(0);
  // @@@ BETA1: the Comix hard remainder against the hand-coded one, vs E_gamma
  s["BETA1_CHECK"].SetDefault(0);
  // Hand Comix CEEX's own (non-conserving) beta_1 arguments rather than a
  // mapped conserving configuration. Needs COMIX: MOMENTUM_PROJECTION: 0.
  s["COMIX_REAL_CEEX_ARGS"].SetDefault(false);
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
  s["DUMP_XMIN_EACH"].SetDefault(0.0); // min x_gamma of EVERY dumped photon
  s["COMIX_CHECK"].SetDefault(0);      // @@@ CEEXCX, hand-coded vs Comix
  s["GOLDEN"].SetDefault(0);           // @@@ CEEXGOLD, regression stream
  s["WEIGHT_PROBE"].SetDefault(0);     // @@@ CEEXWT, what makes a heavy event
  /*
    The perturbative order of the CEEX amplitude. 1 (default): beta_0 +
    beta_1, today's weight bit for bit. 2: adds the tree-level double real
    beta_2 (hep-ph/0006359 eq. 95, KKMC's GPS_HiiPlus/HffPlus/HifPlus) built
    from Comix's two-photon amplitude per partition, Ceex_Beta2.C. It needs
    the two-photon real provider (YFS: NLO_Part containing W, which builds
    the e+e- -> X + 2 photons process); without it beta_2 is zero and says
    so once. Neither the real-virtual nor the two-loop virtual is included.
  */
  s["ORDER"].SetDefault(1);
  /*
    beta_2 only for pairs with BOTH photons above x = 2E/sqrt(s) in the CEEX
    frame (the beam CMS). KKMC's vcut2 = xpar(42) = 0.05 is the same variable
    (E/E_beam). Below the cut beta_2 = 0, the honest value, as for
    BETA1_XCUT. Only photons ENUMERATED in the partition sum can pair, so
    the effective cut is max(BETA2_XCUT, SOFT_PARTITION_CUT).
  */
  s["BETA2_XCUT"].SetDefault(0.05);
  /*
    External virtual at ORDER 2: false (default) keeps the O(alpha^1)
    composition sum_h |A_2 + (v/2) A_0|^2; true puts (1+v/2) on beta_0 AND
    beta_1,
    sum_h |A_2 + (v/2) A_1|^2, which is KKMC's O(alpha^2) structure
    ((1+d_I)(1+d_F) on r in HiniPlus/HfinPlus). The difference is the
    O(alpha^2) real-virtual v x beta_1 interference; see
    NOTES-ceex-order2-2026-09-26.md. Ignored at ORDER 1 (CEEX: REAL_VIRTUAL
    factorisable is the same composition at ORDER 1).
  */
  s["ORDER2_VIRTUAL_ON_BETA1"].SetDefault(false);
  /*
    CEEX: REAL_VIRTUAL (default off; name or old integer) - the O(alpha^2)
    real-virtual on the
    one-photon residuals, so that ACRAIC carries what YFS.NLO carries with
    YFS: VIRTUAL_COMBINE product and YFS: RV_MODE remainder
    (NOTES-yfsnlo-realvirtual-2026-09-28.md, ACRAIC section). External/auto
    virtual only.
    off (0): unchanged, sum_h |A_1 + (v/2) A_0|^2 (at ORDER 2 as
       ORDER2_VIRTUAL_ON_BETA1 says).
    factorisable (1): the factorisable part, sum_h |(1 + v/2) A_1|^2 (A_2 at ORDER 2): v/2
       on every beta_1 as well as on beta_0 - ORDER2_VIRTUAL_ON_BETA1 true, now
       also at ORDER 1. KKMC CEEX2's (1+d_I)(1+d_F) on beta_1, and YFS.NLO's
       (1 + v) x real.
    averaged (2): factorisable plus the non-factorisable, helicity-averaged
       remainder: photon j's M_1 (in every
       partition that carries it) gets (v_{n+1,j} - v_B)/2 on top, the
       helicity-averaged loop over tree of the (n+1)-body point minus the
       Born's, both IR subtracted on their own legs. The numbers are YFS.NLO's
       (NLO_Base, YFS: RV_MODE remainder: same loop call, photons matched by
       lab momentum), so ACRAIC needs no loop call of its own; it needs
       NLO_Part with E and YFS: RV_MODE remainder. In |A|^2 this adds sum_j dv_j Re<A_1, M_1j>
       (RealVirtualRemainderRho), at one photon dv |M_1|^2 = YFS.NLO's
       rho dv on the same event.
  */
  s["REAL_VIRTUAL"].SetDefault(ceexrv::off);
  s["BETA2_CLOSURE"].SetDefault(0);    // @@@ B2CLOS, n = 2 closure
  s["BETA2_SOFT_TEST"].SetDefault(0);  // @@@ B2SOFT, soft limits (N events)
  s["BETA2_TRACE"].SetDefault(0);      // @@@ B2TRACE, per pair and partition
  s["KKMC_FLUX_EMULATION"].SetDefault(0); // diagnostic, see Ceex_Base.H
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

  /*
    A zero here is physics, not an error: for a neutrino pair v_f = a_f, so
    the coupling (v_e + a_e)(v_f - a_f) of one helicity vanishes identically.
    This layer is the hand-coded 2 -> 2 Born; with the Comix Born it only
    feeds the 2 -> 2 alignment and diagnostics, and nothing at all beyond
    2 -> 2 (unit alignment), so it is reported at debugging level only.
  */
  if (zcpl == 0.) msg_Debugging() << METHOD << ": Z coupling is zero (j=" << j
                                  << ", mode=" << mode << ")\n";
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
    m_AmpBeta2.m_A[f]    = Complex(0., 0.);
    m_AmpFluxAll.m_A[f]  = Complex(0., 0.);
    m_AmpFluxSoft.m_A[f] = Complex(0., 0.);
    m_AmpExpo2.m_A[f]    = Complex(0., 0.);
  }
}



void Ceex_Base::MakeRho() {
  double sum0(0.), sum1(0.), sum01(0.);
  const int nh(Amplitude::NHel());
  for (int f = 0; f < nh; ++f) {
    sum0 += std::real(m_AmpExpo0.m_A[f] * conj(m_AmpExpo0.m_A[f]));
    sum1 += std::real(m_AmpExpo1.m_A[f] * conj(m_AmpExpo1.m_A[f]));
    sum01 += std::real(conj(m_AmpExpo0.m_A[f]) * m_AmpExpo1.m_A[f]);
  }
  // Average over the four initial-state helicity configurations.
  m_result0  = sum0 / 4.;
  m_result   = sum1 / 4.;
  m_result01 = sum01 / 4.;
  m_result1    = m_result;
  if (!m_ifi_coherent) {   // CEEX: IFI 0, the incoherent partition sum
    m_result0  = m_inc00 / 4.;
    m_result   = m_inc11 / 4.;
    m_result01 = m_inc01 / 4.;
    m_result1  = m_result;
  }
  m_result01o1 = m_result01;
  m_resultV    = m_result0;
  /*
    YFS: ME_PROBE - one-photon events: rho_1 (the coherent assembly), the
    map-independent sum_f |M_1|^2 of the Comix table for the drawn photon
    helicity (m_comixM1, the mapped copy, is a permutation of that plane) and
    the crude, next to the photon's x and the event's yfs photon count; the
    fixed-order side prints r/(S~ B) for the same photon (@@@ SUB8).
  */
  { static const bool mp(ATOOLS::Settings::GetMainSettings()["YFS"]
                         ["ME_PROBE"].SetDefault(0).Get<int>() != 0);
    if (mp && NPhot() == 1) {
      double m1(0.);
      for (int f = 0; f < nh; ++f) m1 += std::norm(m_comixM1.m_A[f]);
      std::ostringstream o;
      o<<std::setprecision(10)<<"@@@ CEEXMP x="<<std::setprecision(6)
       <<2.*m_allphotons[0][0]/sqrt(m_s)<<std::setprecision(10)
       <<" rho1="<<m_result<<" rho0="<<m_result0<<" M1sq="<<m1/4.
       <<" rhocr="<<m_rhocrud;
      // the crude's final-stage pieces for this photon: |s_F|^2 on the
      // physical legs, the pseudo-flux (q+k)^2/q^2 of the F assignment, and
      // |s_F|^2 on the pre-emission legs (the generator's)
      for (int g(0); g < m_nstages && g < (int)m_Sfac.size(); ++g) {
        if (g == m_initstage || m_Sfac[g].empty()) continue;
        Vec4D q; Complex spre(0., 0.);
        bool ok(m_prefsr.size() == m_flavs.size());
        for (size_t l(0); l < m_stagelegs[g].size(); ++l) {
          const int lg(m_stagelegs[g][l].leg);
          if (lg < (int)m_pceex.size()) q += m_pceex[lg];
          if (ok && lg < (int)m_prefsr.size())
            spre += m_stagelegs[g][l].w * SfactorLeg(m_prefsr[lg], m_allphotons[0], m_PhoHel[0]);
          else ok = false;
        }
        const double q2(q.Abs2()), pf(q2 > 0. ? (q + m_allphotons[0]).Abs2()/q2 : 0.);
        o<<" g="<<g<<" sFpost2="<<std::norm(m_Sfac[g][0])<<" pflux="<<pf
         <<" sFpre2="<<(ok ? std::norm(spre) : -1.);
      }
      if (m_initstage < (int)m_Sfac.size() && !m_Sfac[m_initstage].empty())
        o<<" sI2="<<std::norm(m_Sfac[m_initstage][0]);
      o<<"\n";
      std::cerr<<o.str();
    } }
  /*
    ORDER 2: A_2 = A_1 + the beta_2 sum, and the weight (m_result) becomes
    rho_2. m_result01/m_resultV are what an external virtual composes with
    (YFS_Handler): v Re<A_v, A_2> + (v/2)^2 |A_v|^2 with A_v = A_0, or A_1
    under ORDER2_VIRTUAL_ON_BETA1. At ORDER 1 none of this runs and every
    number above is the one it always was.
  */
  // CEEX: REAL_VIRTUAL >= 1 at ORDER 1: v/2 on A_1 as a whole
  if (m_order == 1 && RealVirtualMode() != ceexrv::off) {
    m_result01 = m_result1;
    m_resultV  = m_result1;
  }
  if (m_order == 2) {
    static const bool vb1(ATOOLS::Settings::GetMainSettings()["CEEX"]
                          ["ORDER2_VIRTUAL_ON_BETA1"].Get<bool>()
                          || RealVirtualMode() != ceexrv::off);
    double s2(0.), s12(0.), s02(0.);
    for (int f = 0; f < nh; ++f) {
      m_AmpExpo2.m_A[f] = m_AmpExpo1.m_A[f] + m_AmpBeta2.m_A[f];
      s2  += std::norm(m_AmpExpo2.m_A[f]);
      s12 += std::real(conj(m_AmpExpo1.m_A[f]) * m_AmpExpo2.m_A[f]);
      s02 += std::real(conj(m_AmpExpo0.m_A[f]) * m_AmpExpo2.m_A[f]);
    }
    m_result2  = s2 / 4.;
    static const int kkflux(ATOOLS::Settings::GetMainSettings()["CEEX"]
                            ["KKMC_FLUX_EMULATION"].Get<int>());
    if (kkflux) {
      // diagnostic: rho_1 and rho_2 as KKMC builds them (see Ceex_Base.H)
      double k1(0.), k2(0.);
      for (int f = 0; f < nh; ++f) {
        k1 += std::norm(m_AmpExpo1.m_A[f] + m_AmpFluxAll.m_A[f]);
        k2 += std::norm(m_AmpExpo2.m_A[f] + m_AmpFluxSoft.m_A[f]);
      }
      m_result1 = k1 / 4.;
      m_result2 = k2 / 4.;
    }
    m_result12 = s12 / 4.;
    m_result02 = s02 / 4.;
    m_result   = m_result2;
    m_result01 = vb1 ? m_result12 : m_result02;
    m_resultV  = vb1 ? m_result1 : m_result0;
  }
  double sumbv(0.), sumbr(0.);
  for (int f = 0; f < nh; ++f) {
    sumbv += std::norm(m_AmpBornVirt.m_A[f]);
    sumbr += std::norm(m_AmpBornReal.m_A[f]);
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

ceexrv::code Ceex_Base::RealVirtualMode()
{
  static const ceexrv::code m(ATOOLS::Settings::GetMainSettings()["CEEX"]
    ["REAL_VIRTUAL"].SetDefault(ceexrv::off).Get<ceexrv::code>());
  return m;
}

bool Ceex_Base::PhotonHasM1(size_t j) const
{
  if (j >= m_realphotM1.size()) return false;
  for (int f = 0; f < Amplitude::NHel(); ++f)
    if (std::norm(m_realphotM1[j].m_A[f]) > 0.) return true;
  return false;
}

double Ceex_Base::RealVirtualRemainderRho(const std::vector<double> &dv) const
{
  const Amplitude &A(m_order == 2 ? m_AmpExpo2 : m_AmpExpo1);
  double sum(0.);
  const int nh(Amplitude::NHel());
  for (size_t j(0); j < dv.size() && j < m_realphotM1.size(); ++j) {
    if (dv[j] == 0.) continue;
    double re(0.);
    for (int f = 0; f < nh; ++f)
      re += std::real(conj(A.m_A[f]) * m_realphotM1[j].m_A[f]);
    sum += dv[j]*re;
  }
  return sum/4.;
}

double Ceex_Base::RealFactorPhoton(size_t j) const
{
  if (j >= m_realphot.size() || m_result0 <= 0.) return 0.;
  double sum(0.);
  const int nh(Amplitude::NHel());
  for (int f = 0; f < nh; ++f)
    sum += std::norm(m_AmpExpo0.m_A[f] + m_realphot[j].m_A[f]);
  return sum/4./m_result0 - 1.;
}
