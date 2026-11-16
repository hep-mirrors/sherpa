#ifndef YFS_Main_YFS_Options_H
#define YFS_Main_YFS_Options_H

/*!
  \file YFS_Options.H

  Named values for the multi-valued YFS: and ACRAIC: run-card switches. A card
  may write the name (case-insensitive) or the integer. Pattern: struct
  <option> { enum code {...}; } with operator>> and operator<< for Settings.
*/

#include <iosfwd>
#include <string>
#include <vector>

namespace YFS {

  //! One spelling of a named option value (a name that is an integer is an alias).
  struct Option_Name {
    const char *name;
    int         value;
  };

  //! The value of `tag` (name or integer) in `names`; fatal error naming `option` otherwise.
  int ParseOptionName(const std::string &tag,
                      const std::vector<Option_Name> &names,
                      const std::string &option);
  //! The first name of `value` in `names` (the integer if it has none).
  std::string OptionName(int value, const std::vector<Option_Name> &names);

  //! The -1/0/1 switches: auto (decided per process or event), off, on.
  struct tristate {
    enum code { automatic = -1, off = 0, on = 1 };
  };

  //! YFS: REAL_MAP - the (n+1)-body point for initial-state photons.
  struct realmap {
    enum code { rest_frame = 0, beam_axis = 1, scaled = 2, invariant = 3 };
  };
  //! YFS: REAL_FSR_MAP - the (n+1)-body point for final-state photons.
  struct realfsrmap {
    enum code { rest_frame = 0, rescale_all = 1, dipole = 2, pre_emission = 3 };
  };
  //! YFS: REAL_SUB_EIK - which eikonal beta_1 subtracts and divides by.
  struct realsubeik {
    enum code {
      coherent_point      = 0,
      crude               = 1,
      coherent_born       = 2,
      crude_fsr           = 3,
      post_emission_fsr   = 4,
      coherent_born_fsr   = 5,
      assignment_born_fsr = 6,
      multichannel        = 7,
      multichannel_prefsr = 8
    };
  };
  //! YFS: SUB8_FSR_LEGS - the FSR channel density of REAL_SUB_EIK multichannel_prefsr.
  struct sub8fsrlegs {
    enum code { generation = 0, post_emission = 1, event_density = 2,
                point_density = 3 };
  };
  //! YFS: SUB8_FSRSUB - the FSR subtraction of REAL_SUB_EIK multichannel_prefsr.
  struct sub8fsrsub {
    enum code { generation = 0, point = 1, blend = 4 };
  };
  //! YFS: REAL_FSR_FLUX - the flux a final-state photon's real carries.
  struct realfsrflux {
    enum code { event_flux = 0, no_flux = 1, flux_squared = 2, on_crude = 3,
                own_dipole = 4 };
  };
  //! YFS: REAL_COMBINE - how the photons' real corrections are combined.
  struct realcombine {
    enum code { automatic = -1, sum = 0, product = 1 };
  };
  //! YFS: VIRTUAL_COMBINE - where the virtual enters the squared-level weight.
  struct virtualcombine {
    enum code { sum = 0, product = 1 };
  };
  //! YFS: REAL_BORN_PHOTON_MULTICHANNEL - the Born-photon multichannel weight.
  struct bornphotonmc {
    enum code { off = 0, on = 1, report = 2 };
  };
  //! YFS: BORN_PHOTON_CHANNELS (see NLO_Base::BornPhotonChannelSum)
  struct bornphotonchannels {
    enum code { isr = 0, stage = 1 };
  };
  //! YFS: RV_MODE - how the squared-level real-virtual is built (NLO_Base_Internal.H).
  struct rvmode {
    enum code { legacy = 0, remainder = 1, remainder_own_legs = 2 };
  };
  //! YFS: EW_HIGHER_ORDER - EW corrections beyond the provider's NLO, added to
  //! the Born-level virtual. none: nothing; griffin: GRIFFIN NNLO - NLO (2->2 s-channel).
  struct ewhigherorder {
    enum code { none = 0, griffin = 1 };
  };
  //! YFS: COULOMB_HIGHER_ORDER - W-pair Coulomb correction beyond one loop
  //! (flat scheme), minus its O(alpha) term: none, second (O(alpha^2), eq. (7)
  //! of hep-ph/9507422) or all (non-relativistic all orders, eq. (2)).
  struct coulombho {
    enum code { none = 0, second = 1, all = 2 };
  };
  //! YFS: BETA0_2 and ACRAIC: BETA0_2 - the photonic two-loop beta_0^(2)
  //! (Beta0_2.H). none: off; ll: (1/2) delta_f^2 [+ delta_f w]; llnll: plus the
  //! per-leg NLL; exp: exp(delta_f)(1 + NLL) - 1 - delta_f [+ w(...)]; two_loop:
  //! llnll plus the constant and boxes.
  struct beta02 {
    enum code { none = 0, ll = 1, llnll = 2, exp = 3, two_loop = 4 };
  };
  //! YFS: BETA0_2_PAIRS - the charged pairs making up delta_f: all, or
  //! same_side (both incoming or both outgoing).
  struct beta02pairs {
    enum code { all = 0, same_side = 1 };
  };
  //! YFS: IF_VIRTUAL - the virtual B of an initial-final (one in, one out) dipole, in
  //! the exponent (IFForFac) and in the virtual subtraction alike.
  //! tbvirt_exact_ir (default): KKMC's TBvirt with the exact massive IR coefficient,
  //! exact poles for any masses and TBvirt's finite part (small remainder).
  //! exact: the massive one-in-one-out B (BVirtGeneral / BVirtGeneralEps), any masses;
  //! its finite part leaves a large remainder when a leg is light.
  //! tbvirt: KKMC's TBvirt, whose IR coefficient assumes one leg with m^2 << |t|.
  //! s_channel: the s-channel form in the subtraction, TBvirt in the exponent
  //! (the old, inconsistent IF_VIRTUAL_TCHANNEL false).
  struct ifvirtual {
    enum code { exact = 0, tbvirt = 1, s_channel = 2, tbvirt_exact_ir = 3 };
  };
  //! YFS: BETA0_2_HIGHER - BETA0_2 two_loop beyond O(alpha^2). none: delta^(2)
  //! stops at O(alpha^2). exponentiated: adds the higher orders of the
  //! factorisable exponent, exp(delta_f)(1 + NLL) - 1 - delta_f minus its
  //! O(alpha^2) part, i.e. the LL delta_f^3/6 and beyond (Beta0_2.H).
  struct beta02higher {
    enum code { none = 0, exponentiated = 1 };
  };
  //! ACRAIC: BETA0_2_ON - which ACRAIC amplitude delta^(2) multiplies. born:
  //! A_0, the Born-level terms only. beta1: A_1, also the one-hard-photon beta_1
  //! terms, i.e. the factorised one-photon two-loop delta^(2) beta_1 (O(alpha^3),
  //! exact in the soft and collinear limits). Born-weighted path only.
  struct beta02on {
    enum code { born = 0, beta1 = 1 };
  };
  //! YFS: BETA0_2_ALPHA - coupling of delta_f: alpha0 = alpha(0), loop = the provider's alpha.
  struct beta02alpha {
    enum code { alpha0 = 0, loop = 1 };
  };
  //! YFS: DYSON_VP - Dyson-resummed photon vacuum polarisation (Dyson_VP.H).
  //! none: off; resummed: beyond the provider's first order (full-EW provider);
  //! complete: the whole (1 - Delta r)/(1 - Pi) - 1 (provider without the photon
  //! self-energy).
  struct dysonvp {
    enum code { none = 0, resummed = 1, complete = 2 };
  };
  //! YFS: M1_SQUARE - |m_1|^2 of YFS.NLO's beta_0^(2). averaged: (v/2)^2 from
  //! the helicity-summed v; helicity: Born-weighted sum_h |delta(h)|^2 (averaged
  //! where helicity factors are missing).
  struct m1square {
    enum code { averaged = 0, helicity = 1 };
  };
  //! ACRAIC: EW_HIGHER_ORDER_HELICITY - EW_HIGHER_ORDER in the helicity virtual.
  //! none: not in delta(h); averaged: v_HO/2 on every helicity; helicity:
  //! GRIFFIN NNLO - NLO per helicity, scaled so its Born-weighted 2 Re is v_HO.
  struct ewhohel {
    enum code { none = 0, averaged = 1, helicity = 2 };
  };
  //! YFS: RR_MODE - how the squared-level double real is built (NLO_Base_Internal.H).
  struct rrmode {
    enum code { legacy = 0, exact = 1, singles_product = 2 };
  };
  //! YFS: RV_LOOP_FRAME - the frame the (n+1)-point loop is evaluated in.
  struct rvloopframe {
    enum code { lab = 0, canonical = 1, tilted = 2, beam_axis = 3 };
  };
  //! YFS: VIRTUAL_SUB_LEGS - the dipoles of the Born virtual's IR subtraction:
  //! own = the loop's own point (post-ISR beams), event = the event's dipoles
  //! (poles no longer cancel after initial-state radiation, diagnostic only),
  //! own_if_exponent (default) = own, plus the initial-final exponent's leg difference
  //! Y_IF(own legs) - Y_IF(event legs) in v. The exponent's IF dipoles take the
  //! beams before ISR, the loop the post-ISR beams; with the difference in v an
  //! IF_VIRTUAL change cancels between exponent and v in every event.
  struct virtsublegs {
    enum code { event = 0, own = 1, own_if_exponent = 2 };
  };
  //! YFS: VIRTUAL_MASSIFICATION - the external fermions the loop provider sees massless.
  struct virtmassif {
    //! light_fermions: e, mu, u, d, s. light_fermions_charm: plus c (W-pair
    //! energies and above). Massification's 4 m^2/s guard applies.
    enum code { off = 0, electrons = 1, leptons = 2, light_fermions = 3,
                light_fermions_charm = 4 };
  };
  //! YFS: RV_PHOTON_CT - the emitted photon's charge counterterm in the RV.
  struct rvphotonct {
    enum code { none = 0, calibrated = 1, analytic = 2 };
  };
  //! YFS: Flux_Mode - which photon the real's initial-state flux is taken at.
  struct fluxmode {
    enum code { event = 0, mapped = 1, average = 2 };
  };

  /*!
    YFS: RR_CONVENTIONS - the pieces of the legacy double real (RR_MODE legacy)
    taken from the single real of each photon. A list of names, or an integer
    mask (1, 2, 4, 8 in the order below); `legacy` (0) takes none.
  */
  struct rrconventions {
    bool single_flux        = false;  //!< 1: the single real's flux
    bool single_subtraction = false;  //!< 2: its subtraction eikonal
    bool single_denominator = false;  //!< 4: its denominator
    bool soft_limit_legs    = false;  //!< 8: eikonals on the other photon's legs
    //! none of the above: the old assembly
    bool Legacy() const
    { return !single_flux && !single_subtraction && !single_denominator
             && !soft_limit_legs; }
  };
  rrconventions ReadRRConventions(const std::vector<std::string> &tags);

  //! ACRAIC: CRUDE_FROM_GENERATOR - the crude on the generator's own stages.
  struct crudegen {
    enum code { automatic = -1, off = 0, on = 1, compare = 2 };
  };
  //! ACRAIC: TCHANNEL_MULTIPHOTON - beta_1 beyond one photon with the reduced
  //! t-channel Born.
  struct tchmultiphoton {
    enum code { partition_legs = 0, one_photon_point = 1, factorised = 2 };
  };
  //! ACRAIC: NO_PSEUDOFLUX - where the pseudo-flux svarY/svarQ on beta_0 sits.
  struct pseudoflux {
    enum code { rho0_and_rho1 = 0, neither = 1, rho0_only = 2 };
  };
  //! ACRAIC: REAL_VIRTUAL - the O(alpha^2) real-virtual: off, factorisable,
  //! averaged (helicity-averaged remainder), helicity (remainder per helicity
  //! from Recola), isr_nll (factorisable plus analytic LL+NLL initial-state
  //! remainder), hybrid (isr_nll for ISR photons, averaged for FSR photons).
  struct acraicrv {
    enum code { off = 0, factorisable = 1, averaged = 2, automatic = 3,
                helicity = 4, isr_nll = 5, hybrid = 6 };
  };
  //! ACRAIC: CRUDE_BORN_LEGS - the Born point of the crude: reduced, generator
  //! (pre-emission pair direction in its rest frame) or stage (pair rebuilt
  //! along its physical direction in the reduced system's frame).
  struct crudebornlegs {
    enum code { reduced = 0, generator = 1, stage = 2 };
  };
  //! ACRAIC: BETA1_LEGS - the legs beta_1's M_1 is evaluated on.
  struct beta1legs {
    enum code { physical = 0, balanced = 1 };
  };
  //! ACRAIC: W_EIKONAL - the W momentum a W-stage eikonal uses.
  struct weikonal {
    enum code { daughters = 0, partition = 1 };
  };
  //! ACRAIC: RESONANCE_EXPONENTIATION - the resonant initial-final log of a
  //! stage's s-channel resonance: left in the provider's virtual (none), or
  //! exponentiated into the resonance propagator in every partition amplitude
  //! (greco). See Acraic_Base::ResonanceExponents.
  struct acraicresexp {
    enum code { none = 0, greco = 1, automatic = 2 };
  };

  //! ACRAIC: TNP_IFI_VIRTUAL - the initial-final virtual of the flipped weight:
  //! ACRAIC's own boxes (e+e- -> f fbar), the mirrored loop, both by process
  //! (auto), or none (virtual even, no boxes)
  struct tnpifivirtual {
    enum code { automatic = 0, own_boxes = 1, mirrored = 2, none = 3 };
  };

  std::istream &operator>>(std::istream &str, tristate::code &c);
  std::ostream &operator<<(std::ostream &str, const tristate::code &c);
  std::istream &operator>>(std::istream &str, realmap::code &c);
  std::ostream &operator<<(std::ostream &str, const realmap::code &c);
  std::istream &operator>>(std::istream &str, realfsrmap::code &c);
  std::ostream &operator<<(std::ostream &str, const realfsrmap::code &c);
  std::istream &operator>>(std::istream &str, realsubeik::code &c);
  std::ostream &operator<<(std::ostream &str, const realsubeik::code &c);
  std::istream &operator>>(std::istream &str, sub8fsrlegs::code &c);
  std::ostream &operator<<(std::ostream &str, const sub8fsrlegs::code &c);
  std::istream &operator>>(std::istream &str, sub8fsrsub::code &c);
  std::ostream &operator<<(std::ostream &str, const sub8fsrsub::code &c);
  std::istream &operator>>(std::istream &str, realfsrflux::code &c);
  std::ostream &operator<<(std::ostream &str, const realfsrflux::code &c);
  std::istream &operator>>(std::istream &str, realcombine::code &c);
  std::ostream &operator<<(std::ostream &str, const realcombine::code &c);
  std::istream &operator>>(std::istream &str, virtualcombine::code &c);
  std::ostream &operator<<(std::ostream &str, const virtualcombine::code &c);
  std::istream &operator>>(std::istream &str, bornphotonmc::code &c);
  std::ostream &operator<<(std::ostream &str, const bornphotonmc::code &c);
  std::istream &operator>>(std::istream &str, bornphotonchannels::code &c);
  std::ostream &operator<<(std::ostream &str, const bornphotonchannels::code &c);
  std::istream &operator>>(std::istream &str, ewhigherorder::code &c);
  std::ostream &operator<<(std::ostream &str, const ewhigherorder::code &c);
  std::istream &operator>>(std::istream &str, coulombho::code &c);
  std::ostream &operator<<(std::ostream &str, const coulombho::code &c);
  std::istream &operator>>(std::istream &str, beta02::code &c);
  std::ostream &operator<<(std::ostream &str, const beta02::code &c);
  std::istream &operator>>(std::istream &str, beta02pairs::code &c);
  std::ostream &operator<<(std::ostream &str, const beta02pairs::code &c);
  std::istream &operator>>(std::istream &str, ifvirtual::code &c);
  std::ostream &operator<<(std::ostream &str, const ifvirtual::code &c);
  std::istream &operator>>(std::istream &str, beta02on::code &c);
  std::ostream &operator<<(std::ostream &str, const beta02on::code &c);
  std::istream &operator>>(std::istream &str, beta02higher::code &c);
  std::ostream &operator<<(std::ostream &str, const beta02higher::code &c);
  std::istream &operator>>(std::istream &str, beta02alpha::code &c);
  std::ostream &operator<<(std::ostream &str, const beta02alpha::code &c);
  std::istream &operator>>(std::istream &str, dysonvp::code &c);
  std::ostream &operator<<(std::ostream &str, const dysonvp::code &c);
  std::istream &operator>>(std::istream &str, m1square::code &c);
  std::ostream &operator<<(std::ostream &str, const m1square::code &c);
  std::istream &operator>>(std::istream &str, ewhohel::code &c);
  std::ostream &operator<<(std::ostream &str, const ewhohel::code &c);
  std::istream &operator>>(std::istream &str, rvmode::code &c);
  std::ostream &operator<<(std::ostream &str, const rvmode::code &c);
  std::istream &operator>>(std::istream &str, rrmode::code &c);
  std::ostream &operator<<(std::ostream &str, const rrmode::code &c);
  std::istream &operator>>(std::istream &str, rvloopframe::code &c);
  std::ostream &operator<<(std::ostream &str, const rvloopframe::code &c);
  std::istream &operator>>(std::istream &str, virtsublegs::code &c);
  std::ostream &operator<<(std::ostream &str, const virtsublegs::code &c);
  std::istream &operator>>(std::istream &str, virtmassif::code &c);
  std::ostream &operator<<(std::ostream &str, const virtmassif::code &c);
  std::istream &operator>>(std::istream &str, rvphotonct::code &c);
  std::ostream &operator<<(std::ostream &str, const rvphotonct::code &c);
  std::istream &operator>>(std::istream &str, fluxmode::code &c);
  std::ostream &operator<<(std::ostream &str, const fluxmode::code &c);
  std::istream &operator>>(std::istream &str, crudegen::code &c);
  std::ostream &operator<<(std::ostream &str, const crudegen::code &c);
  std::istream &operator>>(std::istream &str, tchmultiphoton::code &c);
  std::ostream &operator<<(std::ostream &str, const tchmultiphoton::code &c);
  std::istream &operator>>(std::istream &str, pseudoflux::code &c);
  std::ostream &operator<<(std::ostream &str, const pseudoflux::code &c);
  std::istream &operator>>(std::istream &str, acraicrv::code &c);
  std::ostream &operator<<(std::ostream &str, const acraicrv::code &c);
  std::istream &operator>>(std::istream &str, crudebornlegs::code &c);
  std::ostream &operator<<(std::ostream &str, const crudebornlegs::code &c);
  std::istream &operator>>(std::istream &str, beta1legs::code &c);
  std::ostream &operator<<(std::ostream &str, const beta1legs::code &c);
  std::istream &operator>>(std::istream &str, weikonal::code &c);
  std::ostream &operator<<(std::ostream &str, const weikonal::code &c);
  std::istream &operator>>(std::istream &str, acraicresexp::code &c);
  std::ostream &operator<<(std::ostream &str, const acraicresexp::code &c);
  std::istream &operator>>(std::istream &str, tnpifivirtual::code &c);
  std::ostream &operator<<(std::ostream &str, const tnpifivirtual::code &c);

}

#endif
