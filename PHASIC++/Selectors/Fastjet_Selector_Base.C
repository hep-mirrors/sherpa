#include "PHASIC++/Selectors/Fastjet_Selector_Base.H"

#ifdef USING__FASTJET

#include "ATOOLS/Org/Run_Parameter.H"

using namespace PHASIC;
using namespace ATOOLS;

Fastjet_Selector_Base::Fastjet_Selector_Base(const std::string& name,
                                             Process_Base* const proc,
                                             Scoped_Settings s):
  Selector_Base(name, proc),
  m_eekt(0), p_jdef(0),
  p_siscplug(nullptr), p_eecamplug(nullptr), p_jadeplug(nullptr)
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

  fastjet::RecombinationScheme recom;
  if      (reco=="E")     recom=fastjet::E_scheme;
  else if (reco=="pt")    recom=fastjet::pt_scheme;
  else if (reco=="pt2")   recom=fastjet::pt2_scheme;
  else if (reco=="Et")    recom=fastjet::Et_scheme;
  else if (reco=="Et2")   recom=fastjet::Et2_scheme;
  else if (reco=="BIpt")  recom=fastjet::BIpt_scheme;
  else if (reco=="BIpt2") recom=fastjet::BIpt2_scheme;
  else THROW(fatal_error, "Unknown recombination scheme \"" + reco + "\".");

  bool ee(rpa->gen.Beam1().IsLepton() && rpa->gen.Beam2().IsLepton());

  fastjet::JetAlgorithm ja(fastjet::kt_algorithm);
  if (algo=="cambridge") ja = fastjet::cambridge_algorithm;
  if (algo=="antikt")    ja = fastjet::antikt_algorithm;
  if (algo=="siscone")   p_siscplug = new fastjet::SISConePlugin(m_delta_r, m_f);
  if (ee) {
    if (algo=="eecambridge") p_eecamplug = new fastjet::EECambridgePlugin(m_delta_r);
    if (algo=="jade")        p_jadeplug  = new fastjet::JadePlugin();
  }

  if (p_siscplug) p_jdef=new fastjet::JetDefinition(p_siscplug);
  else if (p_eecamplug) p_jdef=new fastjet::JetDefinition(p_eecamplug);
  else if (p_jadeplug) p_jdef=new fastjet::JetDefinition(p_jadeplug);
  else if (ee) {
    p_jdef=new fastjet::JetDefinition(fastjet::ee_kt_algorithm);
    m_eekt=1;
  }
  else p_jdef=new fastjet::JetDefinition(ja,m_delta_r,recom);

  m_smin = Max(sqr(m_ptmin),sqr(m_etmin));
}

Fastjet_Selector_Base::~Fastjet_Selector_Base()
{
  delete p_jdef;
  if (p_siscplug) delete p_siscplug;
  if (p_eecamplug) delete p_eecamplug;
  if (p_jadeplug) delete p_jadeplug;
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

#endif
