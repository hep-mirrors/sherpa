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
  // The amplitude container's extent: 2 -> 2 gives 4 legs, 16 helicity entries.
  Amplitude::SetLegs(flavs);

  RegisterDefaults();
  Scoped_Settings s{ Settings::GetMainSettings()["ACRAIC"] };
  Settings& ss = Settings::GetMainSettings();
  m_onlyz = s["ONLYZ"].Get<bool>();
  m_onlyg = s["ONLYG"].Get<bool>();
  m_checkxs = s["CHECK_XS"].Get<int>();
  m_comixflip = s["COMIX_REAL_FLIP"].Get<int>();
  m_comixnorm = s["COMIX_REAL_NORM"].Get<double>();
  rpa->gen.AddCitation(1, "ACRAIC builds on the coherent exclusive exponentiation of"
                       " \\cite{Jadach:2000ir} as implemented in KKMC \\cite{Jadach:1999vf}.");
  m_tnp          = s["TNP_STREAMS"].Get<bool>();
  if (m_tnp)
    rpa->gen.AddCitation(1, "The theory-nuisance streams of ACRAIC: TNP_STREAMS follow"
                         " \\cite{Tackmann:2024kci}.");
  m_tnpifi       = m_tnp && s["TNP_IFI"].Get<bool>();
  m_tnpifichk    = s["TNP_IFI_CHECK"].Get<int>();
  m_comixborn    = s["COMIX_BORN"].Get<bool>();
  m_wstages      = s["W_STAGES"].Get<tristate::code>();
  m_weikonal     = s["W_EIKONAL"].Get<weikonal::code>();
  m_crudegen     = s["CRUDE_FROM_GENERATOR"].Get<crudegen::code>();
  m_pfluxmassive = s["PSEUDOFLUX_MASSIVE"].Get<bool>();
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
  m_treebetamax = s["TREE_BETA_MAX"].Get<int>();
  if (m_treebetamax < 0) m_treebetamax = m_order;
  if (m_treebetamax > 1 && m_order != 2)
    THROW(fatal_error, "ACRAIC: TREE_BETA_MAX > 1 needs ACRAIC: ORDER 2: the tree "
          "residuals beta_m, m >= 2, enter only the ORDER 2 amplitude A_2.");
  m_loopbetamax = s["LOOP_BETA_MAX"].Get<int>();
  m_loopbetaxcut = s["LOOP_BETA_XCUT"].Get<double>();
  // a pair below BETA2_XCUT has no tree beta_2 to carry the residual
  m_loopbetaxcut = Max(m_loopbetaxcut, s["BETA2_XCUT"].Get<double>());
  if (m_loopbetamax < 1) m_loopbetamax = 1;
  if (m_loopbetamax > 1) {
    if (m_order != 2)
      THROW(fatal_error, "ACRAIC: LOOP_BETA_MAX > 1 needs ACRAIC: ORDER 2.");
    if (m_loopbetamax > m_treebetamax)
      THROW(fatal_error, "ACRAIC: LOOP_BETA_MAX > TREE_BETA_MAX: a one-loop "
            "residual beta_m^[1] needs the tree beta_m of the same set.");
    if (RealVirtualMode() != acraicrv::helicity)
      THROW(fatal_error, "ACRAIC: LOOP_BETA_MAX > 1 needs ACRAIC: REAL_VIRTUAL helicity "
            "(the per-photon one-loop remainders and the Born's per-helicity "
            "factors enter the two-photon residual).");
    if (m_loopbetamax > 3)
      THROW(fatal_error, "ACRAIC: LOOP_BETA_MAX > 3: the four-photon one-loop "
            "amplitude is out of reach (NOTES-acraic-oneloop-nphoton.md).");
    msg_Info()<<"ACRAIC: LOOP_BETA_MAX "<<m_loopbetamax<<" - one-loop "
              <<m_loopbetamax<<"-photon providers built (x > "<<m_loopbetaxcut
              <<"). The two-photon one-loop residual enters the amplitude."<<std::endl;
  }
  if (m_order == 2 && m_treebetamax != 2)
    msg_Info()<<"ACRAIC: TREE_BETA_MAX "<<m_treebetamax<<" - tree residuals "
              <<"beta_m for 2 <= m <= "<<m_treebetamax<<" from Comix m-photon "
              <<"amplitudes, photons with x > BETA2_XCUT."<<std::endl;
  m_treeprobe = s["TREE_BETA_PROBE"].Get<int>() > 0;
  string widthscheme = ss["WIDTH_SCHEME"].Get<string>();
  m_fixedwidth = (widthscheme == "Fixed" || widthscheme == "CMS");
  m_flavs = flavs;

  static const bool devmultileg(s["DEV_MULTILEG"].SetDefault(false).Get<bool>());
  /* The 2 -> 2 special path (hand-coded Born, own virtual, m_if1/m_if2) describes
     e+e- -> f fbar. */
  m_ffbar = (flavs.size() == 4 && flavs[2].IsFermion()
             && flavs[3] == flavs[2].Bar());
  if (!m_ffbar) {
    /* The throw concerns the Born: it applies only to the hand-coded one.
       DEV_MULTILEG overrides it. */
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

  m_QIQF  = flavs[0].Charge() * flavs[m_if1].Charge();
  m_bhabha = (flavs[0] == flavs[m_if1]);
  m_MZ = Flavour(kf_Z).Mass();
  m_gZ = Flavour(kf_Z).Width();
  m_I   = Complex(0., 1.);

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
  m_ae       = 2.*m_flavs[0].IsoWeak();
  m_af       = 2.*m_flavs[m_if1].IsoWeak();
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
  s["WEAK"].SetDefault(false);
  /*
    COMIX_REAL_FLIP < 0: derive the fermion bits from the Born, set the photon bit.
    COMIX_REAL_NORM <= 0: derive from sum|A_comix|^2 / sum|e^2 A_hand|^2.
  */
  s["COMIX_REAL_FLIP"].SetDefault(-1);
  s["COMIX_REAL_NORM"].SetDefault(-1.);
  /* TNP_STREAMS: per-event theory-nuisance-parameter weights. */
  s["TNP_STREAMS"].SetDefault(false);
  // ACRAIC: TNP_IFI: IF template T_IFI = Odd[w_k] - Odd[w_0], 2 -> 2 only
  s["TNP_IFI"].SetDefault(true);
  s["TNP_IFI_CHECK"].SetDefault(0);
  s["BETA1_XCUT"].SetDefault(1e-3);          // beta_1 soft threshold
  /* Photons with x_gamma = 2E/sqrt(s) below this are fixed to the stage with the
     largest eikonal instead of being enumerated. */
  s["SOFT_PARTITION_CUT"].SetDefault(0.);
  /*
    Legs for beta_1's M_1 (ComixInfraredSubtracted_1_0):
      physical (0, default): physical legs, s-channel propagator momenta shifted
        by the other photons per side.
      balanced (1): a rebuilt balanced point (PartitionLegs).
  */
  s["BETA1_LEGS"].SetDefault(beta1legs::physical);
  /*
    Evaluate space-like exchange lines (t-channel bosons or fermions) of the
    partition Born and M_1 at the partition's reduced invariant (see
    AddExchangeLineShifts). on (1), off (0), auto (-1, default): on for 2 -> 2.
  */
  s["TCHANNEL_SHIFT"].SetDefault(tristate::automatic);
  /*
    W stages for e+e- -> W+W- -> 4f: each W decay is its own charge-neutral stage
    (W incoming, charged daughter outgoing) and the W's are outgoing legs of the
    production stage, so a decay-stage photon shifts only its W propagator.
      off (0): no W stages.  on (1): when DipoleSet::FindWW recognises the state.
      auto (-1): on when both W's are within CLUSTERING_THRESHOLD widths of the pole.
  */
  s["W_STAGES"].SetDefault(tristate::on);
  /*
    W momentum in a decay- or production-stage eikonal:
      daughters (0, default): W at its daughters (per event, soft-factor table
        computed once; misses Comix's soft limit by O(K_decay/M_W)).
      partition (1): the pole follows the partition's decay photons, as Comix's
        shifted propagator does; soft factors recomputed per partition.
  */
  s["W_EIKONAL"].SetDefault(weikonal::daughters);
  /*
    Crude built on the generator's stages instead of ACRAIC's (see
    CrudeFromGenerator): off (0), on (1), compare only (2), auto (-1, default):
    on exactly when W stages are active.
  */
  s["CRUDE_FROM_GENERATOR"].SetDefault(crudegen::automatic);
  /* beta_0(X_wp) on the real phase-space point of invariant X_wp^2 (beams at
     X_wp, radiating pair at X_wp minus the spectators), per partition. */
  s["BORN_AT_SPRIME"].SetDefault(false);
  /*
    Born with a space-like exchange line: the physical-spinor partition Born with
    only its poles moved is not the reduced-point Born times a flux, since the
    numerators do not scale with the poles.
      on (1): beta_0 is the Born at the real reduced point (BornLegsAt(X_wp)),
        the crude is that Born times s/X_wp^2, and M_1 is evaluated on the
        balanced point (BETA1_LEGS 1), so soft cancellations hold between like
        objects.
      off (0): the physical-spinor forms, as for s-channel Borns.
      auto (-1, default): on when the Born has such a line (BornHasExchangeLine)
        and there is a single radiating stage.
  */
  s["TCHANNEL_REDUCED_BORN"].SetDefault(tristate::automatic);
  /*
    TCHANNEL_MULTIPHOTON: beta_1 beyond one photon with the reduced t-channel
    Born (TCHANNEL_REDUCED_BORN active).
      partition_legs (0, default): M_1 on PartitionLegs, eikonal on those legs.
      one_photon_point (1): M_1 at the one-photon point (OnePhotonScaledLegs) as
        M_1/s(point) times the physical eikonal.
      factorised (2): one_photon_point with the beta_1 terms combined along the
        Born helicity vector (AddFactorisedRemainder).
  */
  s["TCHANNEL_MULTIPHOTON"].SetDefault(tchmultiphoton::partition_legs);
  /* BETA1_DRESSED_FSR: beta_1 of a photon on a two-leg final stage with
     companions is taken at the pair dressed with them, with its own Born
     subtraction (ComixInfraredSubtracted_1_0). */
  s["BETA1_DRESSED_FSR"].SetDefault(false);

  s["NO_PSEUDOFLUX"].SetDefault(pseudoflux::rho0_only);
  /*
    PSEUDOFLUX_MASSIVE: the final stage's pseudo-flux is s_pre/s_Q for massless
    legs; for massive legs it is (s_pre/s_Q) beta(s_pre)/beta(s_Q). Pairs with
    FSR_SHELL_JACOBIAN.
  */
  s["PSEUDOFLUX_MASSIVE"].SetDefault(false);
  /*
    CRUDE_BORN_LEGS: Born point of the crude. generator (default) keeps the
    pair direction in the pair rest frame (GeneratorBornAt); reduced uses the
    reduced system's frame (BornLegsAt), which differs by the recoil against a
    hard final-state photon.
  */
  s["CRUDE_BORN_LEGS"].SetDefault(crudebornlegs::generator);
  /* Partition Born from Comix instead of the hand-coded 2 -> 2 spinors.
     Requires COMIX: MOMENTUM_PROJECTION: 0 (checked in the constructor). */
  /* Default on: the hand-coded 2 -> 2 spinor Born is the fallback
     (COMIX_BORN false). */
  s["COMIX_BORN"].SetDefault(true);
  /*
    ISR_VERTEX_AT_S: put the external virtual's initial-state vertex at s as
    ACRAIC's own virtual and KKMC's CEEX do, adding 2 Re[delta_I(s) - delta_I(s')]
    to the virtual evaluated at the post-ISR s'.
  */
  s["ISR_VERTEX_AT_S"].SetDefault(true);
  /* Diagnostics, all off by default. */
  s["DUMP_XMIN"].SetDefault(0.0);      // min x_gamma for the point dump
  s["DUMP_NPHOT"].SetDefault(1);       // photon multiplicity to dump at
  /* Perturbative order of the amplitude: 1 (default) beta_0 + beta_1, 2 adds beta_2. */
  s["ORDER"].SetDefault(1);

  s["BETA2_XCUT"].SetDefault(0.05);
  /* TREE_BETA_MAX (-1 = ORDER): highest tree residual beta_m. */
  s["TREE_BETA_MAX"].SetDefault(-1);
  /*
    TREE_BETA_PROBE N: for the first N events (ORDER 2) print the truncation
    ladder |A_N|^2 against the exact Comix tree |M_n|^2 (@@@ BNCLOS). Closure
    needs BETA1_XCUT, BETA2_XCUT, SOFT_PARTITION_CUT 0 and TREE_BETA_MAX >= n.
  */
  s["TREE_BETA_PROBE"].SetDefault(0);
  /*
    LOOP_BETA_MAX (default 1): highest photon multiplicity with a one-loop
    amplitude. 2 builds the two-photon provider (needs ORDER 2, <= TREE_BETA_MAX)
    and enters A_2 as beta_2^[1] (Acraic_LoopN.C).
    LOOP_BETA_XCUT (-1 = BETA2_XCUT, never below it): softer photons get no
    one-loop residual.
  */
  s["LOOP_BETA_MAX"].SetDefault(1);
  s["LOOP_BETA_XCUT"].SetDefault(-1.);
  /* LOOP_LADDER_PROBE N: print the one-loop truncation ladder (@@@ LLCLOS)
     for the first N events (LoopLadderProbe). */
  s["LOOP_LADDER_PROBE"].SetDefault(0);
  /*
    External virtual at ORDER 2: false (default) sum_h |A_2 + (v/2) A_0|^2;
    true sum_h |A_2 + (v/2) A_1|^2.
  */
  s["ORDER2_VIRTUAL_ON_BETA1"].SetDefault(false);
  /*
    REAL_VIRTUAL: the O(alpha^2) real-virtual on the one-photon residuals.
      off (0): sum_h |A_1 + (v/2) A_0|^2 (at ORDER 2 per ORDER2_VIRTUAL_ON_BETA1).
      factorisable (1): sum_h |(1 + v/2) A_1|^2, v/2 on every beta_1 as on beta_0
        (as CEEX2's (1+d_I)(1+d_F)).
      averaged (2): factorisable plus the helicity-averaged non-factorisable
        remainder (v_{n+1,j} - v_B)/2 on photon j's M_1, taken from YFS.NLO
        (needs NLO_Part with E and YFS: RV_MODE remainder); adds
        sum_j dv_j Re<A_1, M_1j> (RealVirtualRemainderRho).
      isr_nll (5): factorisable plus, on initial-state photons, the LL+NLL
        single-emission one-loop term V_2 of KKMC's MakeVini (Phys. Rev. D65 (2002) 073030),
        dv_j = 2 Re V_2(r_1, r_2), r_i = p_i.k/p_1.p_2. No loop call. Final-state
        photons keep the factorisable factor.
      hybrid (6): isr_nll on initial-state photons, averaged on final-state ones
        (needs what averaged needs). HYBRID_ISR_RMIN r_c >= 0: V_2 only where
        min(p_i.k)/p_1.p_2 < r_c.
  */
  s["REAL_VIRTUAL"].SetDefault(acraicrv::automatic);
  s["HYBRID_ISR_RMIN"].SetDefault(-1.);
  /* RESONANCE_EXPONENTIATION (none | greco | auto): exponentiate each final
     resonance's initial-final log into its propagator in every amplitude, and
     remove it from the virtual (Acraic_Resonance.C). */
  s["RESONANCE_EXPONENTIATION"].SetDefault(acraicresexp::automatic);
  s["TNP_IFI_VIRTUAL"].SetDefault(tnpifivirtual::automatic);
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

  /* A zero is physics (neutrino pair: v_f = a_f); this hand-coded layer is
     diagnostic only with the Comix Born. */
  if (zcpl == 0.) msg_Debugging() << METHOD << ": Z coupling is zero (j=" << j
                                  << ", mode=" << mode << ")\n";
  return zcpl;
}



Complex Acraic_Base::CouplingG() {
  m_gcpl = Complex(m_QIQF, 0);
  return m_gcpl;
}


/*
  The momenta the ACRAIC spin amplitudes are built from: the beams and the
  physical (lab) outgoing legs, the only set with p_a + p_b = sum p_f + sum(k).
*/
void Acraic_Base::BuildAcraicMomenta()
{
  /* m_momenta is {beam, beam, Born final state, lab final state}; use the lab
     set, falling back to the Born set if it has not been appended. */
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
    for (int c(0); c < tnp_nclass; ++c) m_tnpamp[c].m_A[f] = Complex(0., 0.);
    m_tnpflip0.m_A[f] = m_tnpflip1.m_A[f] = Complex(0., 0.);
    m_tnpabs[f] = 0.;
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
