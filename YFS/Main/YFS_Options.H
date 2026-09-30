#ifndef YFS_Main_YFS_Options_H
#define YFS_Main_YFS_Options_H

/*
  Named values for the multi-valued YFS: and CEEX: run-card switches.

  Every option here used to be a bare integer (REAL_SUB_EIK: 8,
  SUB8_FSR_LEGS: 3, ...). Each value now has a name, and a card may write
  either the name (case-insensitive) or the old integer, so old cards and the
  Ares scripts keep working unchanged. What each value does is documented
  where the option's default is set; the enum below only lists the names.

  Pattern: struct <option> { enum code {...}; }, as for ceexvirt, wwscheme
  and yfsmode in YFS_Base.H, with operator>> and operator<< so that the
  Settings can read and print the value.
*/

#include <iosfwd>
#include <string>
#include <vector>

namespace YFS {

  //! One spelling of a named option value. A name that is an integer
  //! ("4") is an alias for an old integer whose behaviour another value has.
  struct Option_Name {
    const char *name;
    int         value;
  };

  //! The value named `tag` - a name, case-insensitive, or the old integer -
  //! in `names`; throws a fatal error naming `option` otherwise.
  int ParseOptionName(const std::string &tag,
                      const std::vector<Option_Name> &names,
                      const std::string &option);
  //! The first name of `value` in `names` (the integer if it has none).
  std::string OptionName(int value, const std::vector<Option_Name> &names);

  //! The -1/0/1 switches: auto (decided per process or event), off, on.
  struct tristate {
    enum code { automatic = -1, off = 0, on = 1 };
  };

  // ------------------------------------------------------------ YFS: (NLO)

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
  //! YFS: RV_MODE - how the squared-level real-virtual is built.
  struct rvmode {
    enum code { legacy = 0, remainder = 1 };
  };
  //! YFS: RR_MODE - how the squared-level double real is built.
  struct rrmode {
    enum code { legacy = 0, exact = 1 };
  };
  //! YFS: RV_LOOP_FRAME - the frame the (n+1)-point loop is evaluated in.
  struct rvloopframe {
    enum code { lab = 0, canonical = 1, tilted = 2, beam_axis = 3 };
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
    YFS: RR_CONVENTIONS - the pieces of the legacy double real (RR_MODE
    legacy) taken from the single real's own calculation of each photon.
    Read from a list of names, or from the old integer mask (1, 2, 4, 8 in
    the order below); `legacy` (0) takes none of them.
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
  std::string   RRConventionsNames(const rrconventions &c);

  // ------------------------------------------------------------ CEEX:

  //! CEEX: CRUDE_FROM_GENERATOR - the crude on the generator's own stages.
  struct crudegen {
    enum code { automatic = -1, off = 0, on = 1, compare = 2 };
  };
  //! CEEX: TCHANNEL_MULTIPHOTON - beta_1 beyond one photon with the reduced
  //! t-channel Born.
  struct tchmultiphoton {
    enum code { partition_legs = 0, one_photon_point = 1, factorised = 2 };
  };
  //! CEEX: NO_PSEUDOFLUX - where the pseudo-flux svarY/svarQ on beta_0 sits.
  struct pseudoflux {
    enum code { rho0_and_rho1 = 0, neither = 1, rho0_only = 2 };
  };
  //! CEEX: REAL_VIRTUAL - the O(alpha^2) real-virtual in ACRAIC. The names
  //! are those of the ap-port-yfs-nlo-helrv work, which adds helicity = 3.
  struct ceexrv {
    enum code { off = 0, factorisable = 1, averaged = 2 };
  };
  //! CEEX: BETA1_LEGS - the legs beta_1's M_1 is evaluated on.
  struct beta1legs {
    enum code { physical = 0, balanced = 1 };
  };
  //! CEEX: W_EIKONAL - the W momentum a W-stage eikonal uses.
  struct weikonal {
    enum code { daughters = 0, partition = 1 };
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
  std::istream &operator>>(std::istream &str, rvmode::code &c);
  std::ostream &operator<<(std::ostream &str, const rvmode::code &c);
  std::istream &operator>>(std::istream &str, rrmode::code &c);
  std::ostream &operator<<(std::ostream &str, const rrmode::code &c);
  std::istream &operator>>(std::istream &str, rvloopframe::code &c);
  std::ostream &operator<<(std::ostream &str, const rvloopframe::code &c);
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
  std::istream &operator>>(std::istream &str, ceexrv::code &c);
  std::ostream &operator<<(std::ostream &str, const ceexrv::code &c);
  std::istream &operator>>(std::istream &str, beta1legs::code &c);
  std::ostream &operator<<(std::ostream &str, const beta1legs::code &c);
  std::istream &operator>>(std::istream &str, weikonal::code &c);
  std::ostream &operator<<(std::ostream &str, const weikonal::code &c);

}

#endif
