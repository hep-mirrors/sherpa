#include "PHASIC++/Selectors/Fastjet_Selector_Base.H"

#include "ATOOLS/Org/Run_Parameter.H"

using namespace PHASIC;
using namespace ATOOLS;

Fastjet_Selector_Base::Fastjet_Selector_Base(const std::string& name,
                                             Process_Base* const proc,
                                             Scoped_Settings s):
  Selector_Base(name, proc),
  m_eekt(0), p_jdef(0)
{
  // parameter/mode settings
  const auto algo = s["Algorithm"]          .SetDefault("")  .Get<std::string>();
  const auto reco = s["RecombinationScheme"].SetDefault("E") .Get<std::string>();
  m_delta_r       = s["DR"]                 .SetDefault(0.4) .Get<double>();
  m_f             = s["f"]                  .SetDefault(0.75).Get<double>();

  // min/max settings
  m_nj    = s["N"]     .SetDefault("None").UseZeroReplacements()     .Get<size_t>();
  m_ptmin = s["PTMin"] .SetDefault("None").UseZeroReplacements()     .Get<double>();
  m_etmin = s["ETMin"] .SetDefault("None").UseZeroReplacements()     .Get<double>();
  m_eta   = s["EtaMax"].SetDefault("None").UseMaxDoubleReplacements().Get<double>();
  m_y     = s["YMax"]  .SetDefault("None").UseMaxDoubleReplacements().Get<double>();

  fjcore::RecombinationScheme recom;
  if      (reco=="E")     recom=fjcore::E_scheme;
  else if (reco=="pt")    recom=fjcore::pt_scheme;
  else if (reco=="pt2")   recom=fjcore::pt2_scheme;
  else if (reco=="Et")    recom=fjcore::Et_scheme;
  else if (reco=="Et2")   recom=fjcore::Et2_scheme;
  else if (reco=="BIpt")  recom=fjcore::BIpt_scheme;
  else if (reco=="BIpt2") recom=fjcore::BIpt2_scheme;
  else THROW(fatal_error, "Unknown recombination scheme \"" + reco + "\".");

  bool ee(rpa->gen.Beam1().IsLepton() && rpa->gen.Beam2().IsLepton());

  fjcore::JetAlgorithm ja(fjcore::kt_algorithm);
  if (algo=="cambridge") ja = fjcore::cambridge_algorithm;
  if (algo=="antikt")    ja = fjcore::antikt_algorithm;

  if (ee) {
    p_jdef=new fjcore::JetDefinition(fjcore::ee_kt_algorithm);
    m_eekt=1;
  }
  else p_jdef=new fjcore::JetDefinition(ja,m_delta_r);

  m_smin = Max(sqr(m_ptmin),sqr(m_etmin));
}

Fastjet_Selector_Base::~Fastjet_Selector_Base()
{
  delete p_jdef;
}

void Fastjet_Selector_Base::PrintCommonInfoLines(std::ostream& str, size_t width)
{
  str<<width<<"  Algorithm: kt (default)|antikt|cambridge|siscone   # hadron colliders\n"
     <<width<<"  Algorithm: eekt (default)|jade|eecambridge|siscone # lepton colliders\n"
     <<width<<"  N: number of jets\n"
     <<width<<"  # optional settings:\n"
     <<width<<"  PTMin: minimum jet pT\n"
     <<width<<"  ETMin: minimum jet eta\n"
     <<width<<"  DR: jet distance parameter\n"
     <<width<<"  f: Siscone f parameter (default: 0.75)\n"
     <<width<<"  EtaMax: maximum jet eta (default: None)\n"
     <<width<<"  YMax: maximum jet rapidity (default: None)\n";
}
