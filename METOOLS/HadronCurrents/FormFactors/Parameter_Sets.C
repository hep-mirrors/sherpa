#include "METOOLS/HadronCurrents/FormFactors/Parameter_Sets.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/YAML/yaml-cpp/yaml.h"
#include <algorithm>
#include <set>
#include <dirent.h>

using namespace METOOLS;
using namespace ATOOLS;

namespace {

  typedef std::vector<long> Canon;
  struct Parameter_Set { std::string file; std::map<std::string,double> values; };
  typedef std::map<std::pair<Canon,int>,Parameter_Set> Parameter_Set_Map;

  Canon CanonKey(const std::string & key) {
    Canon out;
    std::string item;
    for (char c: key+",") {
      if (c==',') { if (!item.empty()) out.push_back(std::stol(item)); item.clear(); }
      else if (c!=' ') item += c;
    }
    std::sort(out.begin(),out.end());
    return out;
  }

  const Parameter_Set_Map & ParameterSets(const std::string & dir) {
    static std::map<std::string,Parameter_Set_Map> cache;
    auto cached = cache.find(dir);
    if (cached!=cache.end()) return cached->second;
    Parameter_Set_Map & out = cache[dir];
    DIR * d = opendir(dir.c_str());
    if (d==NULL)
      THROW(fatal_error,"Parameter_Sets_Path '"+dir+"' cannot be opened.");
    std::vector<std::string> names;
    while (dirent * e = readdir(d)) {
      const std::string n = e->d_name;
      if (n.size()>5 && n.compare(n.size()-5,5,".yaml")==0) names.push_back(n);
    }
    closedir(d);
    std::sort(names.begin(),names.end());
    for (const auto & name: names) {
      SHERPA_YAML::Node node = SHERPA_YAML::LoadFile(dir+"/"+name);
      if (!node["sherpa"] || !node["applies_to"]) {
        msg_Info()<<"### parameter set "<<name<<" has no generated sherpa block "
                  <<"- skipped; run export_sherpa_blocks.py\n";
        continue;
      }
      Parameter_Set set;
      set.file = name;
      // A set may RECORD more channels than its generated block implements -
      // FM95 covers ten but only its three K pi pi currents have a local path,
      // because the K K pi ones need a width the paper never prints.  When the
      // block says so, only those channels select it; the rest fall back to the
      // registry rather than selecting a set that would be half-applied.
      std::set<Canon> covered;
      const auto only = node["sherpa"]["sherpa_applies_to"];
      if (only) for (const auto & c: only) covered.insert(CanonKey(c.as<std::string>()));
      for (const auto & kv: node["sherpa"]) {
        const std::string key = kv.first.as<std::string>();
        if (key=="sherpa_applies_to") continue;
        set.values[key] = kv.second.as<double>();
      }
      for (const auto & a: node["applies_to"]) {
        const std::string ch = a["channel"].as<std::string>();
        const int model = a["model"].as<int>();
        const auto key = std::make_pair(CanonKey(ch),model);
        if (!covered.empty() && !covered.count(key.first)) continue;
        if (out.count(key))
          THROW(fatal_error,"("+ch+", "+std::to_string(model)+") is claimed by both "
                +out[key].file+" and "+name+" - a (channel, model) pair must "
                "select exactly one set.");
        out[key] = set;
      }
    }
    return out;
  }

  //! Flavour codes of the given entries, sorted, plus the charge conjugate.
  void Keys(const std::vector<long> & codes,Canon & key,Canon & conj) {
    key.clear(); conj.clear();
    for (const long & c: codes) { key.push_back(c); conj.push_back(-c); }
    std::sort(key.begin(),key.end());
    std::sort(conj.begin(),conj.end());
  }

}

void METOOLS::ApplyParameterSet(GeneralModel * model,
                                const ATOOLS::Flavour_Vector & flavs,
                                const std::vector<int> & hadrons,
                                const int & ffmodel)
{
  if (model==NULL) return;
  auto hds = Settings::GetMainSettings()["HADRON_DECAYS"];
  const std::string policy =
    hds["Resonance_Parameters"].SetDefault("global").Get<std::string>();
  if (policy=="global") return;
  if (policy!="local" && policy!="local_else_global")
    THROW(fatal_error,"HADRON_DECAYS:Resonance_Parameters must be global, "
          "local or local_else_global, not '"+policy+"'.");
  const std::string dir =
    hds["Parameter_Sets_Path"].SetDefault("").Get<std::string>();
  if (dir.empty())
    THROW(fatal_error,"Resonance_Parameters = "+policy+" needs "
          "HADRON_DECAYS:Parameter_Sets_Path.");

  // Primary key: every flavour but the first, which for a decay is the list of
  // products.  Secondary: the hadrons this form factor actually serves, which
  // is what a production process has to be keyed on.
  std::vector<long> products, hadcodes;
  std::string label;
  for (size_t i=1; i<flavs.size(); ++i) {
    const long code = (long int)flavs[i];
    products.push_back(code);
    label += (i>1 ? "," : "")+std::to_string(code);
  }
  for (const int & i: hadrons)
    if (i>=0 && size_t(i)<flavs.size()) hadcodes.push_back((long int)flavs[i]);

  Canon pk, pc, hk, hc;
  Keys(products,pk,pc);
  Keys(hadcodes,hk,hc);

  const Parameter_Set_Map & sets = ParameterSets(dir);
  auto it = sets.end();
  for (const Canon * k: { &pk, &pc, &hk, &hc }) {
    it = sets.find(std::make_pair(*k,ffmodel));
    if (it!=sets.end()) break;
  }

  // One form factor per slot means one current can land here several times
  // (F1, F2, F3, FS all share a channel); report each pair once.
  static std::set<std::pair<Canon,int> > reported;
  const bool first = reported.insert(std::make_pair(pk,ffmodel)).second;

  if (it==sets.end()) {
    if (policy=="local")
      THROW(fatal_error,"Resonance_Parameters = local, but no parameter set "
            "exists for ("+label+", FORM_FACTOR = "+std::to_string(ffmodel)+").");
    if (first)
      msg_Out()<<"### parameter set for ("<<label<<", FORM_FACTOR = "<<ffmodel
               <<"): none - registry [global]\n";
    return;
  }
  // Merge every key the set exports.  Different models consume different
  // names: the KS/GS tower reads PS_*, while RChL2012 reads the same
  // parameter names it already uses (Mrho_rchl, gammaMag_pipi_RChL2012, ...)
  // and a RCHL_* marker.  The form factor decides what to do with them.
  size_t n = 0;
  for (const auto & kv: it->second.values) { (*model)[kv.first] = kv.second; ++n; }
  if (first)
    msg_Out()<<"### parameter set for ("<<label<<", FORM_FACTOR = "<<ffmodel
             <<"): "<<it->second.file<<" ["<<n<<" values, local]\n";
}
