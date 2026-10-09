#include "YFS/ACRAIC/Acraic_Base.H"
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

Acraic_Base::Acraic_Base(const Flavour_Vector &flavs)
{
  // The amplitude container's extent. 2 -> 2 gives 4 legs and 16 helicity
  // entries, which is what the fixed m_A[2][2][2][2] used to hold.
  // Counts the legs that carry two helicity states, not the legs: a scalar
  // in the final state is not packed, which is what makes the container size
  // agree with Comix's Spin_Amplitudes.
  Amplitude::SetLegs(flavs);

  RegisterDefaults();
  Scoped_Settings s{ Settings::GetMainSettings()["ACRAIC"] };
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
  m_resexpmode   = s["RESONANCE_EXPONENTIATION"].Get<acraicresexp::code>();
  if (m_resexpmode == acraicresexp::automatic &&
      (s["BETA1_LEGS"].Get<beta1legs::code>() == beta1legs::balanced ||
       s["BORN_AT_SPRIME"].Get<bool>())) {
    msg_Info()<<"ACRAIC: RESONANCE_EXPONENTIATION auto - none: BETA1_LEGS "
              <<"balanced or BORN_AT_SPRIME carry no propagator shifts."<<std::endl;
    m_resexpmode = acraicresexp::none;
  }
  if (ResonanceExponentiation()) {
    // charged resonance stages (ACRAIC: W_STAGES, active only for W+W- -> 4f)
    // are refused by the first event that activates them
    // (BuildResonanceExponents)
    if (s["BETA1_LEGS"].Get<beta1legs::code>() == beta1legs::balanced)
      THROW(fatal_error, "ACRAIC: RESONANCE_EXPONENTIATION: greco needs "
            "ACRAIC: BETA1_LEGS physical; the balanced legs evaluate M_1 "
            "without propagator shifts, which carry the factor.");
    if (s["BORN_AT_SPRIME"].Get<bool>())
      THROW(fatal_error, "ACRAIC: RESONANCE_EXPONENTIATION: greco needs "
            "ACRAIC: BORN_AT_SPRIME false; the real-point Born carries no "
            "propagator shifts, which carry the factor.");
    msg_Info()<<"ACRAIC: RESONANCE_EXPONENTIATION greco - the resonant "
              <<"initial-final log of each final stage's resonance is "
              <<"exponentiated into its propagator in every amplitude and "
              <<"taken out of the virtual."<<std::endl;
  }
  if (m_order != 1 && m_order != 2) {
    msg_Error()<<METHOD<<"(): ACRAIC: ORDER = "<<m_order<<" is not 1 or 2; "
               <<"using 1."<<std::endl;
    m_order = 1;
  }
  if (m_order == 2)
    msg_Info()<<"ACRAIC: ORDER 2 - tree-level double real beta_2 from Comix "
              <<"two-photon amplitudes, pairs with x = 2E/sqrt(s) > "
              <<s["BETA2_XCUT"].Get<double>()<<"."<<std::endl;
  // ACRAIC: TREE_BETA_MAX: -1 (default) follows ORDER, so nothing changes
  // unless it is set; above 2 it needs ORDER 2 (A_2 is where beta_m lives)
  m_treebetamax = s["TREE_BETA_MAX"].Get<int>();
  if (m_treebetamax < 0) m_treebetamax = m_order;
  if (m_treebetamax > 1 && m_order != 2)
    THROW(fatal_error, "ACRAIC: TREE_BETA_MAX > 1 needs ACRAIC: ORDER 2: the tree "
          "residuals beta_m, m >= 2, enter only the ORDER 2 amplitude A_2.");
  if (m_order == 2 && m_treebetamax != 2)
    msg_Info()<<"ACRAIC: TREE_BETA_MAX "<<m_treebetamax<<" - tree residuals "
              <<"beta_m for 2 <= m <= "<<m_treebetamax<<" from Comix m-photon "
              <<"amplitudes, photons with x > BETA2_XCUT."<<std::endl;
  m_treeprobe = s["TREE_BETA_PROBE"].Get<int>() > 0;
   /*
    COMIX: MOMENTUM_PROJECTION used to be switched off GLOBALLY here, because
    ACRAIC hands Comix legs that do not balance. That also removed Comix's
    collinear-stability fix from every other evaluation in the run
  */
  string widthscheme = ss["WIDTH_SCHEME"].Get<string>();
  m_fixedwidth = (widthscheme == "Fixed" || widthscheme == "CMS");
  m_flavs = flavs;

  static const bool devmultileg(s["DEV_MULTILEG"].SetDefault(false).Get<bool>());
  /*
    The 2 -> 2 special path (hand-coded Born, its Comix alignment and map,
    ACRAIC's own virtual, the fermion-pair m_if1/m_if2) is a description of
    e+e- -> f fbar.
  */
  m_ffbar = (flavs.size() == 4 && flavs[2].IsFermion()
             && flavs[3] == flavs[2].Bar());
  if (!m_ffbar) {
    /*
      The throw is about the BORN, so it applies only when the hand-coded Born
      is the one selected..

      DEV_MULTILEG remains the override for the hand-coded case.
    */
    if (!m_comixborn && !devmultileg)
      THROW(fatal_error, "ACRAIC past 2 -> 2 needs the Comix Born. Either leave"
            " ACRAIC: COMIX_BORN at its default, or set ACRAIC: DEV_MULTILEG to"
            " use the hand-coded 2 -> 2 spinor structure anyway.");
    msg_Error()<<METHOD<<"(): this process has "<<flavs.size()<<" legs."
               <<" The Born is taken from Comix and is general, but ACRAIC's own"
               <<" virtual is 2 -> 2 (use YFS: ACRAIC_Virtual: external) and the"
               <<" Comix real covers one photon, above which the hand-coded"
               <<" beta_1 runs. Whether this process reduces to ACRAIC's"
               <<" four-fermion structure has to be shown, not assumed."
               <<std::endl;
  }

  FindOutgoingFermionPair(flavs);

 
  if (!m_ffbar && m_useacraic && m_acraicvirtsrc == acraicvirt::acraic)
    THROW(fatal_error,
          "ACRAIC's own virtual is e+e- -> f fbar only, and this process has "
          + ATOOLS::ToString(m_flavs.size()) + " legs. Set YFS: ACRAIC_Virtual:"
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
  if (Settings::GetMainSettings()["ACRAIC"]["REAL_SIN2THETAW"]
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


void Acraic_Base::RegisterDefaults()
{
  Scoped_Settings s{ Settings::GetMainSettings()["ACRAIC"] };
  s["ONLYZ"].SetDefault(false);
  s["ONLYG"].SetDefault(false);
  s["CHECK_XS"].SetDefault(0);
  // false: tree-level couplings in ACRAIC's hand-coded Born; the weak virtual
  // comes from the loop provider (ACRAIC_Virtual: auto/external)
  s["WEAK"].SetDefault(false);
  s["COMIX_REAL"].SetDefault(true);
  /*
      COMIX_REAL_FLIP  < 0  derive the fermion bits from the Born and take
                            the photon bit from COMIX_REAL_PHOTON_FLIP
      COMIX_REAL_NORM <= 0  derive from sum|A_comix|^2 / sum|e^2 A_hand|^2

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
    is the structure ACRAIC actually has - beta_1 is a sum over photons - and it
    needs only the one-photon process, which always exists.
  */
  s["COMIX_REAL_PER_PHOTON"].SetDefault(false);
  s["BETA1_PARTITION_CHECK"].SetDefault(0);  // @@@ B1PART
  s["BETA1_CLOSURE"].SetDefault(0);          // @@@ B1CLOSE
  s["CLOSURE_1PHOT"].SetDefault(0);          // @@@ CLOS1
  s["BETA1_XCUT"].SetDefault(1e-3);          // beta_1 soft threshold
  /*
    NEXT_TO_SOFT_FLUX (default true): a photon that gets no beta_1 in any
    partition (y below BETA1_XCUT even with every other photon initial-state)
    carries the next-to-soft flux of its stage.
  */
  s["NEXT_TO_SOFT_FLUX"].SetDefault(true);
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
    partition sum but fixed to the stage with the largest eikonal 
  */
  s["SOFT_PARTITION_CUT"].SetDefault(0.);
  /*
    MOMENTUM_REPAIR (true): project the final legs and all photons onto exact
    momentum conservation against the beams and onto their mass shells before
    anything is evaluated on them; false = use the event's momenta as handed
    over. See Acraic_Base::RepairMomentumBalance.
  */
  s["MOMENTUM_REPAIR"].SetDefault(true);
  /*
    Legs for beta_1's M_1 (name or old integer):
      physical (0, default): the physical legs with the s-channel propagator
        momenta shifted by the other photons per side
      balanced (1): a rebuilt balanced point (PartitionLegs).
    See ComixInfraredSubtracted_1_0.
  */
  s["BETA1_LEGS"].SetDefault(beta1legs::physical);
  /*
    Space-like exchange lines (Bhabha's t-channel boson, the t/u-channel
    electrons of e+e- -> gamma gamma, a t-channel neutrino: any current with
    one initial leg and part of the final state, detected by leg content in
    COMIX::Amplitude::SetPropShifts) of the partition Born and of M_1 at the
    partition's REDUCED invariant rather than at the unreduced one
    on (1), off (0), auto (-1, default): on
    for 2 -> 2 only (ExchangeLineShiftsOn). See
    Acraic_Base::AddExchangeLineShifts.
  */
  s["TCHANNEL_SHIFT"].SetDefault(tristate::automatic);
  /*
    W stages for e+e- -> W+W- -> 4f. The charged final legs of a 
    WW-type final state are split into two decay
    stages, one per W, each with the W as an incoming leg and its charged
    daughter as the outgoing one, and the W's join the beams as outgoing legs
    of the production stage. Each stage is charge neutral; a photon on a
    decay stage shifts only that W's propagator (COMIX::Amplitude::
    SetPropShifts with the W's daughter mask), a photon on the production
    stage the s-channel line. Without it the two charged leptons of such a
    state form one flat stage whose photons shift BOTH W lines at once

      off (0): every process as before, bit for bit.
      on (1): whenever DipoleSet::FindWW recognises
        the final state; a process without W's is unchanged bit for bit.
       
      auto (-1): on when in addition both W's are within YFS:
        CLUSTERING_THRESHOLD widths of the pole (the pole scheme's own
        window).
    Needs the handler to hand over the W groups (YFS_Handler::ACRAICStageGroups).
  */
  s["W_STAGES"].SetDefault(tristate::on);
  /*
    The W momentum a decay- or production-stage eikonal uses. 
    daughters (0, default): the W at
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
    The crude the ACRAIC weight divides by, built on the GENERATOR's stages
    rather than on ACRAIC's own; see
    Acraic_Base::CrudeFromGenerator.
    off (0); on (1), use it; compare (2), compare only; auto (-1, default),
    use it exactly when W stages are active - the case where ACRAIC's stages
    and the generator's dipoles differ - and leave every other process
    untouched.
  */
  s["CRUDE_FROM_GENERATOR"].SetDefault(crudegen::automatic);
  /*
    beta_0(X_wp) on the REAL phase-space point whose invariant is X_wp^2 -
    beams at X_wp, radiating pair at X_wp minus the spectators, per partition
  */
  s["BORN_AT_SPRIME"].SetDefault(false);
  /*
    Born processes with a space-like exchange line.
    For those the physical-spinor partition Born with only its poles moved
    is not the reduced-point Born times a flux, as it is for an s-channel
    Born: the numerators do not scale with the poles. Tristate, name or old
    integer:
      on (1): beta_0 is the Born at the partition's real reduced point
         (the BORN_AT_SPRIME form, BornLegsAt(X_wp)), the crude that Born
         times s/X_wp^2 (the generator's density), and every photon's M_1
         in beta_1 is evaluated on the balanced point with the other
         photons taken out of the beams (the BETA1_LEGS: 1 form,
         PartitionLegs) - genuine amplitudes throughout, so the soft
         cancellations between beta_0 and beta_1 hold between like objects.
      off (0): the physical-spinor forms, as for s-channel Borns.
      auto (-1, default): on exactly when the Born has such a line, detected once
         per Born process (Acraic_Base::BornHasExchangeLine), and there is a
         single radiating stage (no charged final state, e.g. gamma gamma).
         Every s-channel process keeps the old numbers bit for bit. With a
         final-state stage
  */
  s["TCHANNEL_REDUCED_BORN"].SetDefault(tristate::automatic);
  /*
    ACRAIC: TCHANNEL_MULTIPHOTON - beta_1 beyond one photon when beta_0 is the
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
  */
  s["TCHANNEL_MULTIPHOTON"].SetDefault(tchmultiphoton::partition_legs);
  s["BETA1_BORNLEGS"].SetDefault(1);         // 1 = reduced, 0 = physical

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
    algebra. That algebra is the 2 -> 2 specific part of ACRAIC, so this is the
    step that makes the partition sum process independent. Requires
    COMIX: MOMENTUM_PROJECTION: 0 - see the check in the constructor.
  */
  /*
    ON by default: the Comix Born is the DERIVED object and the hand-coded
    2 -> 2 spinor algebra is the fallback. At 2 -> 2 it reproduces the
    hand-coded ACRAIC weight to round-off; past 2 -> 2 the hand-coded Born is not
    the right amplitude at all. COMIX_BORN: false selects the hand-coded one.
  */
  s["COMIX_BORN"].SetDefault(true);
  // @@@ BORNALIGN: is the Comix -> ACRAIC Born alignment scale independent?
  s["BORN_ALIGN_CHECK"].SetDefault(0);
  // @@@ BETA1: the Comix hard remainder against the hand-coded one, vs E_gamma
  s["BETA1_CHECK"].SetDefault(0);
  // Hand Comix ACRAIC's own (non-conserving) beta_1 arguments rather than a
  // mapped conserving configuration. Needs COMIX: MOMENTUM_PROJECTION: 0.
  s["COMIX_REAL_ACRAIC_ARGS"].SetDefault(false);
  // Write the |beta_1|/beta_0 vs E_gamma scan (the real-validation figure).
  s["BETA1_SCAN"].SetDefault(0);
  // @@@ VIRTPART: does the virtual factor V(h) depend on the partition?
  s["VIRT_PARTITION_CHECK"].SetDefault(0);
  /*
    Diagnostics. All off by default,
  */
  s["BORNNORM_CHECK"].SetDefault(0);   // @@@ BORNANG, Born normalisation
  s["PIN_PHOTON_HEL"].SetDefault(0);   // pin photon helicity (+1/-1) for sums
  s["DUMP_XMIN"].SetDefault(0.0);      // min x_gamma for the point dump
  s["DUMP_NPHOT"].SetDefault(1);       // photon multiplicity to dump at
  s["DUMP_XMIN_EACH"].SetDefault(0.0); // min x_gamma of EVERY dumped photon
  s["COMIX_CHECK"].SetDefault(0);      // @@@ ACRAICCX, hand-coded vs Comix
  s["GOLDEN"].SetDefault(0);           // @@@ ACRAICGOLD, regression stream
  s["WEIGHT_PROBE"].SetDefault(0);     // @@@ ACRAICWT, what makes a heavy event
  /*
    The perturbative order of the ACRAIC amplitude. 1 (default): beta_0 +
    beta_1
  */
  s["ORDER"].SetDefault(1);

  s["BETA2_XCUT"].SetDefault(0.05);
  /*
    ACRAIC: TREE_BETA_MAX (integer; -1, the default, = ORDER):
  */
  s["TREE_BETA_MAX"].SetDefault(-1);
  /*
    ACRAIC: TREE_BETA_PROBE N: for the first N events (ORDER 2) print the
    truncation ladder |A_N|^2, N = 0..TREE_BETA_MAX, against the exact tree
    |M_n|^2 from Comix at the event's point (@@@ BNCLOS). Closure at
    N = n needs every photon in every set: BETA1_XCUT, BETA2_XCUT and
    SOFT_PARTITION_CUT 0, TREE_BETA_MAX >= n.
  */
  s["TREE_BETA_PROBE"].SetDefault(0);
  /*
    External virtual at ORDER 2: false (default) keeps the O(alpha^1)
    composition sum_h |A_2 + (v/2) A_0|^2; true puts (1+v/2) on beta_0 AND
    beta_1,
    sum_h |A_2 + (v/2) A_1|^2, 
  */
  s["ORDER2_VIRTUAL_ON_BETA1"].SetDefault(false);
  /*
    ACRAIC: REAL_VIRTUAL (default off; name or old integer) - the O(alpha^2)
    real-virtual on the
    one-photon residuals, so that ACRAIC carries what YFS.NLO carries with
    YFS: VIRTUAL_COMBINE product and YFS: RV_MODE remainder
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
  s["REAL_VIRTUAL"].SetDefault(acraicrv::automatic);
  /*
    ACRAIC: RESONANCE_EXPONENTIATION (none | greco | auto; default auto) - the
    resonant initial-final log of every final stage's resonance taken out
    of the provider's virtual and exponentiated into that resonance's
    propagator in every amplitude of the partition sum (Acraic_Resonance.C).
  */
  s["RESONANCE_EXPONENTIATION"].SetDefault(acraicresexp::automatic);
  s["RESONANCE_EXPONENTIATION_PROBE"].SetDefault(0);
  s["BETA2_CLOSURE"].SetDefault(0);    // @@@ B2CLOS, n = 2 closure
  s["BETA2_SOFT_TEST"].SetDefault(0);  // @@@ B2SOFT, soft limits (N events)
  s["BETA2_TRACE"].SetDefault(0);      // @@@ B2TRACE, per pair and partition
  s["KKMC_FLUX_EMULATION"].SetDefault(0); // diagnostic, see Acraic_Base.H
}



void Acraic_Base::Init(const Vec4D_Vector &p)
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



void Acraic_Base::MakeProp(bool resonance)
{
  if (m_fixedwidth) {
    m_propZ =  1. / Complex(m_sp - sqr(m_MZ), m_gZ * m_MZ);
  }
  else {
    m_propZ =   1. / Complex(m_sp - sqr(m_MZ), m_gZ * m_sp / m_MZ);
  }
  // ACRAIC: RESONANCE_EXPONENTIATION on the hand-coded 2 -> 2 Born, where it
  // is the Born (ACRAIC: COMIX_BORN false): KKMC's PropZet * exp(IntReson)
  if (resonance && !m_comixborn && m_resexp.size() == 1) {
    const Complex mb2(sqr(m_MZ), -m_MZ * m_gZ);
    m_propZ *= exp(m_resexp[0].a * log((mb2 - m_sp) / mb2));
  }
  m_propG = 1. / Complex(m_sp,0);
  if (m_onlyz)  m_propG = 0;
  if (m_onlyg)  m_propZ = 0;
  m_prop =  m_propZ + m_propG;
}


void Acraic_Base::MakePropT(const Vec4D_Vector &p)
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


void Acraic_Base::MakeEWFF(double svar, double costhd)
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


Complex Acraic_Base::CouplingZ(double  j, int mode) {
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



Complex Acraic_Base::CouplingG() {
  m_gcpl = Complex(m_QIQF, 0);
  return m_gcpl;
}


/*
  The four momenta the ACRAIC spin amplitudes are built from:
  {p_a, p_b, p_c, p_d} with p_c,p_d the PHYSICAL outgoing fermions.

  m_momenta is {beam, beam, Born pair, lab pair} - MakeACRAIC appends the
  lab final state after the four Born momenta - and only the lab pair
  satisfies p_a + p_b = p_c + p_d + sum(k). The Born pair sits at full
  beam energy, so using it puts the final-state spinors at a different
  kinematic point from the propagator, which eq.(one-photon) evaluates at
  X = P - sum(k).
*/
void Acraic_Base::BuildAcraicMomenta()
{
  /*
    m_momenta is {beam, beam, Born final state, lab final state}, MakeACRAIC
    appending the lab set. ACRAIC wants the beams and the LAB legs, those being
    the ones that balance against the photons; the Born set is the fallback
    when MakeACRAIC has not appended yet.

    Generalised from the fixed pair to nf = m_flavs.size() - 2 final legs. At
    nf = 2 the offset is 4 or 2 exactly as before.
  */
  const size_t nf(m_flavs.size() >= 2 ? m_flavs.size() - 2 : 0);
  m_pacraic.clear();
  m_pacraic.push_back(m_momenta[0]);
  m_pacraic.push_back(m_momenta[1]);
  const bool havelab(m_momenta.size() >= 2 + 2*nf);
  const size_t off(havelab ? 2 + nf : 2);
  for (size_t i(0); i < nf && off + i < m_momenta.size(); ++i)
    m_pacraic.push_back(m_momenta[off + i]);
}


void Acraic_Base::ZerAmplit() {
  // Every helicity entry the container actually holds, not the 16 a 2 -> 2
  // final state happens to need.
  const int nh(Amplitude::NHel());
  for (int f = 0; f < nh; ++f) {
    m_AmpExpo0.m_A[f]    = Complex(0., 0.);
    m_AmpExpo1.m_A[f]    = Complex(0., 0.);
    m_AmpBornVirt.m_A[f] = Complex(0., 0.);
    m_AmpBornOnly.m_A[f] = Complex(0., 0.);
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



void Acraic_Base::MakeRho() {
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
  if (!m_ifi_coherent) {   // ACRAIC: IFI 0, the incoherent partition sum
    m_result0  = m_inc00 / 4.;
    m_result   = m_inc11 / 4.;
    m_result01 = m_inc01 / 4.;
    m_result1  = m_result;
  }
  m_result01o1 = m_result01;
  m_resultV    = m_result0;
  m_result0x   = m_result01;
  m_result0v   = m_result0;
  m_vonbeta1   = false;
  ProbeMakeRho(nh);
  /*
    ORDER 2: A_2 = A_1 + the beta_2 sum, and the weight (m_result) becomes
    rho_2. m_result01/m_resultV are what an external virtual composes with
    (YFS_Handler): v Re<A_v, A_2> + (v/2)^2 |A_v|^2 with A_v = A_0, or A_1
    under ORDER2_VIRTUAL_ON_BETA1. At ORDER 1 none of this runs and every
    number above is the one it always was.
  */
  // ACRAIC: REAL_VIRTUAL >= 1 at ORDER 1: v/2 on A_1 as a whole
  if (m_order == 1 && RealVirtualMode() != acraicrv::off) {
    m_result01 = m_result1;
    m_resultV  = m_result1;
    m_result0v = m_result01o1;
    m_vonbeta1 = true;
  }
  if (m_order == 2) {
    static const bool vb1(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                          ["ORDER2_VIRTUAL_ON_BETA1"].Get<bool>()
                          || RealVirtualMode() != acraicrv::off);
    double s2(0.), s12(0.), s02(0.);
    for (int f = 0; f < nh; ++f) {
      m_AmpExpo2.m_A[f] = m_AmpExpo1.m_A[f] + m_AmpBeta2.m_A[f];
      s2  += std::norm(m_AmpExpo2.m_A[f]);
      s12 += std::real(conj(m_AmpExpo1.m_A[f]) * m_AmpExpo2.m_A[f]);
      s02 += std::real(conj(m_AmpExpo0.m_A[f]) * m_AmpExpo2.m_A[f]);
    }
    m_result2  = s2 / 4.;
    static const int kkflux(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                            ["KKMC_FLUX_EMULATION"].Get<int>());
    if (kkflux) {
      // diagnostic: rho_1 and rho_2 as KKMC builds them (see Acraic_Base.H)
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
    if (!m_ifi_coherent) {   // YFS: IFI false, the incoherent partition sum
      m_result2  = m_inc22 / 4.;
      m_result12 = m_inc12 / 4.;
      m_result02 = m_inc02 / 4.;
    }
    m_result   = m_result2;
    m_result01 = vb1 ? m_result12 : m_result02;
    m_resultV  = vb1 ? m_result1 : m_result0;
    m_result0x = m_result02;
    m_result0v = vb1 ? m_result01o1 : m_result0;
    m_vonbeta1 = vb1;
  }
  double sumbv(0.), sumbr(0.);
  for (int f = 0; f < nh; ++f) {
    sumbv += std::norm(m_AmpBornVirt.m_A[f]);
    sumbr += std::norm(m_AmpBornReal.m_A[f]);
  }
  m_resultbv = sumbv / 4.;
  m_resultbr = sumbr / 4.;
  double sumbo(0.);
  for (int f = 0; f < nh; ++f) sumbo += std::norm(m_AmpBornOnly.m_A[f]);
  m_resultbo = sumbo / 4.;
  // YFS: IFI false: VirtualFactor()'s two sums incoherent, as m_result0
  if (!m_ifi) { m_resultbv = m_incbv / 4.; m_resultbo = m_incbo / 4.; }
  m_rho0sum += m_result0;
  m_rho1sum += m_result;
  m_rhobvsum += m_resultbv;
  m_rhobrsum += m_resultbr;
  m_rhocrudsum += m_rhocrud;
  ++m_rhon;
}


void Acraic_Base::Reset() {
  m_result = 0;
}


double Acraic_Base::Xi(const Vec4D p, const Vec4D q) {
  return sqrt((m_zeta * p) / (q * m_zeta));
}

acraicrv::code Acraic_Base::RealVirtualMode()
{
  static const acraicrv::code set(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
    ["REAL_VIRTUAL"].SetDefault(acraicrv::automatic).Get<acraicrv::code>());
  static const acraicrv::code m(set != acraicrv::automatic ? set :
    ATOOLS::Settings::GetMainSettings()["ACRAIC"]["ORDER"].SetDefault(1).Get<int>() == 2
    ? acraicrv::averaged : acraicrv::factorisable);
  return m;
}

bool Acraic_Base::PhotonHasM1(size_t j) const
{
  if (j >= m_realphotM1.size()) return false;
  for (int f = 0; f < Amplitude::NHel(); ++f)
    if (std::norm(m_realphotM1[j].m_A[f]) > 0.) return true;
  return false;
}

double Acraic_Base::RealVirtualRemainderRho(const std::vector<double> &dv) const
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

double Acraic_Base::RealFactorPhoton(size_t j) const
{
  if (j >= m_realphot.size() || m_result0 <= 0.) return 0.;
  double sum(0.);
  const int nh(Amplitude::NHel());
  for (int f = 0; f < nh; ++f)
    sum += std::norm(m_AmpExpo0.m_A[f] + m_realphot[j].m_A[f]);
  return sum/4./m_result0 - 1.;
}

/*
  Where the outgoing fermion pair sits. At 2 -> 2 this is legs 2 and 3 and
  the whole coupling layer said so literally. With anything else in the
  final state it is not: for e+e- -> H l+ l- leg 2 is the HIGGS, so
  m_qf = flavs[2].Charge() made the final-state Z coupling that of a
  neutral scalar and CouplingZ reported "Z coupling is Zero!".
*/

/*
  The outgoing FERMION pair, by inspection rather than by position. Neutral
  fermions count - e+e- -> nu nubar is a legitimate ACRAIC process - so the
  test is on being a fermion, not on carrying charge. At 2 -> 2 this
  returns 2 and 3 and nothing downstream changes.
*/
void Acraic_Base::FindOutgoingFermionPair(const ATOOLS::Flavour_Vector &flavs)
{
  size_t n(0);
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
