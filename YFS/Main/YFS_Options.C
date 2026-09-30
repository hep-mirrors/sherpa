#include "YFS/Main/YFS_Options.H"

#include "ATOOLS/Org/Exception.H"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <istream>
#include <ostream>

using namespace YFS;

namespace {

  std::string Lower(std::string s)
  {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return s;
  }

  bool ParseInteger(const std::string &tag, int &value)
  {
    if (tag.empty()) return false;
    char *end(nullptr);
    const long v(std::strtol(tag.c_str(), &end, 10));
    if (end == tag.c_str() || *end != '\0') return false;
    value = (int)v;
    return true;
  }

  std::string NameList(const std::vector<Option_Name> &names)
  {
    std::string list;
    for (const Option_Name &n : names) {
      int dummy;
      if (ParseInteger(n.name, dummy)) continue;   // integer aliases
      if (!list.empty()) list += ", ";
      list += std::string(n.name) + " (" + std::to_string(n.value) + ")";
    }
    return list;
  }

  template <typename Code>
  std::istream &ReadOption(std::istream &str, Code &code,
                           const std::vector<Option_Name> &names,
                           const char *option)
  {
    std::string tag;
    str >> tag;
    code = static_cast<Code>(ParseOptionName(tag, names, option));
    return str;
  }

  template <typename Code>
  std::ostream &WriteOption(std::ostream &str, const Code &code,
                            const std::vector<Option_Name> &names)
  {
    return str << OptionName(static_cast<int>(code), names);
  }

  // The names of every option, in the enum's order; the first name of a
  // value is the one printed. Integer names are aliases for old integers.

  const std::vector<Option_Name> s_tristate{
    {"auto", tristate::automatic}, {"off", tristate::off}, {"on", tristate::on},
    {"false", tristate::off}, {"true", tristate::on}};

  const std::vector<Option_Name> s_realmap{
    {"rest_frame", realmap::rest_frame}, {"beam_axis", realmap::beam_axis},
    {"scaled", realmap::scaled}, {"invariant", realmap::invariant}};

  // 4 was a separate value once; the code has treated it as 1 since the
  // per-dipole map (2) became the default, and old cards still carry it.
  const std::vector<Option_Name> s_realfsrmap{
    {"rest_frame", realfsrmap::rest_frame},
    {"rescale_all", realfsrmap::rescale_all},
    {"dipole", realfsrmap::dipole},
    {"pre_emission", realfsrmap::pre_emission},
    {"4", realfsrmap::rescale_all}};

  const std::vector<Option_Name> s_realsubeik{
    {"coherent_point", realsubeik::coherent_point},
    {"crude", realsubeik::crude},
    {"coherent_born", realsubeik::coherent_born},
    {"crude_fsr", realsubeik::crude_fsr},
    {"post_emission_fsr", realsubeik::post_emission_fsr},
    {"coherent_born_fsr", realsubeik::coherent_born_fsr},
    {"assignment_born_fsr", realsubeik::assignment_born_fsr},
    {"multichannel", realsubeik::multichannel},
    {"multichannel_prefsr", realsubeik::multichannel_prefsr}};

  const std::vector<Option_Name> s_sub8fsrlegs{
    {"generation", sub8fsrlegs::generation},
    {"post_emission", sub8fsrlegs::post_emission},
    {"event_density", sub8fsrlegs::event_density},
    {"point_density", sub8fsrlegs::point_density}};

  const std::vector<Option_Name> s_sub8fsrsub{
    {"generation", sub8fsrsub::generation}, {"point", sub8fsrsub::point},
    {"blend", sub8fsrsub::blend}};

  const std::vector<Option_Name> s_realfsrflux{
    {"event_flux", realfsrflux::event_flux}, {"no_flux", realfsrflux::no_flux},
    {"flux_squared", realfsrflux::flux_squared},
    {"on_crude", realfsrflux::on_crude}, {"own_dipole", realfsrflux::own_dipole}};

  const std::vector<Option_Name> s_realcombine{
    {"auto", realcombine::automatic}, {"sum", realcombine::sum},
    {"product", realcombine::product}};

  const std::vector<Option_Name> s_virtualcombine{
    {"sum", virtualcombine::sum}, {"product", virtualcombine::product}};

  const std::vector<Option_Name> s_bornphotonmc{
    {"off", bornphotonmc::off}, {"on", bornphotonmc::on},
    {"report", bornphotonmc::report}};

  const std::vector<Option_Name> s_rvmode{
    {"legacy", rvmode::legacy}, {"remainder", rvmode::remainder}};

  const std::vector<Option_Name> s_rrmode{
    {"legacy", rrmode::legacy}, {"exact", rrmode::exact}};

  const std::vector<Option_Name> s_rvloopframe{
    {"lab", rvloopframe::lab}, {"canonical", rvloopframe::canonical},
    {"tilted", rvloopframe::tilted}, {"beam_axis", rvloopframe::beam_axis}};

  const std::vector<Option_Name> s_rvphotonct{
    {"none", rvphotonct::none}, {"calibrated", rvphotonct::calibrated},
    {"analytic", rvphotonct::analytic}};

  const std::vector<Option_Name> s_fluxmode{
    {"event", fluxmode::event}, {"mapped", fluxmode::mapped},
    {"average", fluxmode::average}};

  const std::vector<Option_Name> s_crudegen{
    {"auto", crudegen::automatic}, {"off", crudegen::off},
    {"on", crudegen::on}, {"compare", crudegen::compare}};

  const std::vector<Option_Name> s_tchmultiphoton{
    {"partition_legs", tchmultiphoton::partition_legs},
    {"one_photon_point", tchmultiphoton::one_photon_point},
    {"factorised", tchmultiphoton::factorised}};

  const std::vector<Option_Name> s_pseudoflux{
    {"rho0_and_rho1", pseudoflux::rho0_and_rho1},
    {"neither", pseudoflux::neither}, {"rho0_only", pseudoflux::rho0_only}};

  const std::vector<Option_Name> s_ceexrv{
    {"off", ceexrv::off}, {"factorisable", ceexrv::factorisable},
    {"averaged", ceexrv::averaged}};

  const std::vector<Option_Name> s_beta1legs{
    {"physical", beta1legs::physical}, {"balanced", beta1legs::balanced}};

  const std::vector<Option_Name> s_weikonal{
    {"daughters", weikonal::daughters}, {"partition", weikonal::partition}};

  // RR_CONVENTIONS: the old integer mask's values of the four pieces
  const std::vector<Option_Name> s_rrconventions{
    {"single_flux", 1}, {"single_subtraction", 2}, {"single_denominator", 4},
    {"soft_limit_legs", 8}};

}

int YFS::ParseOptionName(const std::string &tag,
                         const std::vector<Option_Name> &names,
                         const std::string &option)
{
  const std::string low(Lower(tag));
  for (const Option_Name &n : names)
    if (low == Lower(n.name)) return n.value;
  int value(0);
  if (ParseInteger(tag, value))
    for (const Option_Name &n : names)
      if (n.value == value) return value;
  THROW(fatal_error, "Unknown value '" + tag + "' for " + option
                     + "; use one of " + NameList(names) + ".");
}

std::string YFS::OptionName(int value, const std::vector<Option_Name> &names)
{
  for (const Option_Name &n : names) {
    int dummy;
    if (n.value == value && !ParseInteger(n.name, dummy)) return n.name;
  }
  return std::to_string(value);
}

rrconventions YFS::ReadRRConventions(const std::vector<std::string> &tags)
{
  rrconventions c;
  for (const std::string &tag : tags) {
    int mask(0);
    if (Lower(tag) == "legacy") continue;
    if (!ParseInteger(tag, mask))
      mask = ParseOptionName(tag, s_rrconventions, "YFS: RR_CONVENTIONS");
    if (mask < 0 || mask > 15)
      THROW(fatal_error, "YFS: RR_CONVENTIONS " + tag + " is not a mask of "
                         + NameList(s_rrconventions) + ".");
    // the old integer mask, decoded once here: bit value -> named piece
    auto has = [mask](int bit) { return (mask & bit) != 0; };
    c.single_flux        = c.single_flux        || has(1);
    c.single_subtraction = c.single_subtraction || has(2);
    c.single_denominator = c.single_denominator || has(4);
    c.soft_limit_legs    = c.soft_limit_legs    || has(8);
  }
  return c;
}

std::string YFS::RRConventionsNames(const rrconventions &c)
{
  if (c.Legacy()) return "legacy";
  std::string s;
  auto add = [&s](bool on, const char *name) {
    if (!on) return;
    if (!s.empty()) s += ", ";
    s += name;
  };
  add(c.single_flux, "single_flux");
  add(c.single_subtraction, "single_subtraction");
  add(c.single_denominator, "single_denominator");
  add(c.soft_limit_legs, "soft_limit_legs");
  return s;
}

std::istream &YFS::operator>>(std::istream &s, tristate::code &c)
{ return ReadOption(s, c, s_tristate, "a YFS/CEEX auto/off/on switch"); }
std::ostream &YFS::operator<<(std::ostream &s, const tristate::code &c)
{ return WriteOption(s, c, s_tristate); }

std::istream &YFS::operator>>(std::istream &s, realmap::code &c)
{ return ReadOption(s, c, s_realmap, "YFS: REAL_MAP"); }
std::ostream &YFS::operator<<(std::ostream &s, const realmap::code &c)
{ return WriteOption(s, c, s_realmap); }

std::istream &YFS::operator>>(std::istream &s, realfsrmap::code &c)
{ return ReadOption(s, c, s_realfsrmap, "YFS: REAL_FSR_MAP"); }
std::ostream &YFS::operator<<(std::ostream &s, const realfsrmap::code &c)
{ return WriteOption(s, c, s_realfsrmap); }

std::istream &YFS::operator>>(std::istream &s, realsubeik::code &c)
{ return ReadOption(s, c, s_realsubeik, "YFS: REAL_SUB_EIK"); }
std::ostream &YFS::operator<<(std::ostream &s, const realsubeik::code &c)
{ return WriteOption(s, c, s_realsubeik); }

std::istream &YFS::operator>>(std::istream &s, sub8fsrlegs::code &c)
{ return ReadOption(s, c, s_sub8fsrlegs, "YFS: SUB8_FSR_LEGS"); }
std::ostream &YFS::operator<<(std::ostream &s, const sub8fsrlegs::code &c)
{ return WriteOption(s, c, s_sub8fsrlegs); }

std::istream &YFS::operator>>(std::istream &s, sub8fsrsub::code &c)
{ return ReadOption(s, c, s_sub8fsrsub, "YFS: SUB8_FSRSUB"); }
std::ostream &YFS::operator<<(std::ostream &s, const sub8fsrsub::code &c)
{ return WriteOption(s, c, s_sub8fsrsub); }

std::istream &YFS::operator>>(std::istream &s, realfsrflux::code &c)
{ return ReadOption(s, c, s_realfsrflux, "YFS: REAL_FSR_FLUX"); }
std::ostream &YFS::operator<<(std::ostream &s, const realfsrflux::code &c)
{ return WriteOption(s, c, s_realfsrflux); }

std::istream &YFS::operator>>(std::istream &s, realcombine::code &c)
{ return ReadOption(s, c, s_realcombine, "YFS: REAL_COMBINE"); }
std::ostream &YFS::operator<<(std::ostream &s, const realcombine::code &c)
{ return WriteOption(s, c, s_realcombine); }

std::istream &YFS::operator>>(std::istream &s, virtualcombine::code &c)
{ return ReadOption(s, c, s_virtualcombine, "YFS: VIRTUAL_COMBINE"); }
std::ostream &YFS::operator<<(std::ostream &s, const virtualcombine::code &c)
{ return WriteOption(s, c, s_virtualcombine); }

std::istream &YFS::operator>>(std::istream &s, bornphotonmc::code &c)
{ return ReadOption(s, c, s_bornphotonmc, "YFS: REAL_BORN_PHOTON_MULTICHANNEL"); }
std::ostream &YFS::operator<<(std::ostream &s, const bornphotonmc::code &c)
{ return WriteOption(s, c, s_bornphotonmc); }

std::istream &YFS::operator>>(std::istream &s, rvmode::code &c)
{ return ReadOption(s, c, s_rvmode, "YFS: RV_MODE"); }
std::ostream &YFS::operator<<(std::ostream &s, const rvmode::code &c)
{ return WriteOption(s, c, s_rvmode); }

std::istream &YFS::operator>>(std::istream &s, rrmode::code &c)
{ return ReadOption(s, c, s_rrmode, "YFS: RR_MODE"); }
std::ostream &YFS::operator<<(std::ostream &s, const rrmode::code &c)
{ return WriteOption(s, c, s_rrmode); }

std::istream &YFS::operator>>(std::istream &s, rvloopframe::code &c)
{ return ReadOption(s, c, s_rvloopframe, "YFS: RV_LOOP_FRAME"); }
std::ostream &YFS::operator<<(std::ostream &s, const rvloopframe::code &c)
{ return WriteOption(s, c, s_rvloopframe); }

std::istream &YFS::operator>>(std::istream &s, rvphotonct::code &c)
{ return ReadOption(s, c, s_rvphotonct, "YFS: RV_PHOTON_CT"); }
std::ostream &YFS::operator<<(std::ostream &s, const rvphotonct::code &c)
{ return WriteOption(s, c, s_rvphotonct); }

std::istream &YFS::operator>>(std::istream &s, fluxmode::code &c)
{ return ReadOption(s, c, s_fluxmode, "YFS: Flux_Mode"); }
std::ostream &YFS::operator<<(std::ostream &s, const fluxmode::code &c)
{ return WriteOption(s, c, s_fluxmode); }

std::istream &YFS::operator>>(std::istream &s, crudegen::code &c)
{ return ReadOption(s, c, s_crudegen, "CEEX: CRUDE_FROM_GENERATOR"); }
std::ostream &YFS::operator<<(std::ostream &s, const crudegen::code &c)
{ return WriteOption(s, c, s_crudegen); }

std::istream &YFS::operator>>(std::istream &s, tchmultiphoton::code &c)
{ return ReadOption(s, c, s_tchmultiphoton, "CEEX: TCHANNEL_MULTIPHOTON"); }
std::ostream &YFS::operator<<(std::ostream &s, const tchmultiphoton::code &c)
{ return WriteOption(s, c, s_tchmultiphoton); }

std::istream &YFS::operator>>(std::istream &s, pseudoflux::code &c)
{ return ReadOption(s, c, s_pseudoflux, "CEEX: NO_PSEUDOFLUX"); }
std::ostream &YFS::operator<<(std::ostream &s, const pseudoflux::code &c)
{ return WriteOption(s, c, s_pseudoflux); }

std::istream &YFS::operator>>(std::istream &s, ceexrv::code &c)
{ return ReadOption(s, c, s_ceexrv, "CEEX: REAL_VIRTUAL"); }
std::ostream &YFS::operator<<(std::ostream &s, const ceexrv::code &c)
{ return WriteOption(s, c, s_ceexrv); }

std::istream &YFS::operator>>(std::istream &s, beta1legs::code &c)
{ return ReadOption(s, c, s_beta1legs, "CEEX: BETA1_LEGS"); }
std::ostream &YFS::operator<<(std::ostream &s, const beta1legs::code &c)
{ return WriteOption(s, c, s_beta1legs); }

std::istream &YFS::operator>>(std::istream &s, weikonal::code &c)
{ return ReadOption(s, c, s_weikonal, "CEEX: W_EIKONAL"); }
std::ostream &YFS::operator<<(std::ostream &s, const weikonal::code &c)
{ return WriteOption(s, c, s_weikonal); }
