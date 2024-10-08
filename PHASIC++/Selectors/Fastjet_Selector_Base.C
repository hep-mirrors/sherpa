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
  rpa->gen.AddCitation(1, "FastJet is published under \\cite{Cacciari:2011ma}.");

  // parameter/mode settings
  const auto algo = s["Algorithm"]          .SetDefault("")  .Get<std::string>();
  const auto reco = s["RecombinationScheme"].SetDefault("E") .Get<std::string>();
  m_delta_r       = s["DR"]                 .SetDefault(0.4) .Get<double>();

  // min/max settings
  m_nj    = s["N"]     .SetDefault("None").UseZeroReplacements()     .Get<size_t>();
  m_ptmin = s["PTMin"] .SetDefault("None").UseZeroReplacements()     .Get<double>();
  m_etmin = s["ETMin"] .SetDefault("None").UseZeroReplacements()     .Get<double>();
  m_eta   = s["EtaMax"].SetDefault("None").UseMaxDoubleReplacements().Get<double>();
  m_y     = s["YMax"]  .SetDefault("None").UseMaxDoubleReplacements().Get<double>();

  fjcore::RecombinationScheme recom = ATOOLS::GetRecombinationScheme(reco);

  bool ee(rpa->gen.Bunch(0).IsLepton() && rpa->gen.Bunch(1).IsLepton());

  fjcore::JetAlgorithm ja(fjcore::kt_algorithm);
  if (algo=="cambridge") ja = fjcore::cambridge_algorithm;
  if (algo=="antikt")    ja = fjcore::antikt_algorithm;

  if (ee) {
    p_jdef=new fjcore::JetDefinition(fjcore::ee_kt_algorithm);
    m_eekt=1;
  } else
    p_jdef = new fjcore::JetDefinition(ja, m_delta_r, recom);

  m_smin = Max(sqr(m_nj*m_ptmin),sqr(m_nj*m_etmin));
}

Fastjet_Selector_Base::~Fastjet_Selector_Base()
{
  delete p_jdef;
}

void Fastjet_Selector_Base::PrintCommonInfoLines(std::ostream& str, size_t width)
{
  str << width
      << "  Algorithm: kt (default)|antikt|cambridge   # hadron colliders\n"
      << width << "  Algorithm: eekt (default) # lepton colliders\n"
      << width << "  N: number of jets\n"
      << width << "  # optional settings:\n"
      << width << "  PTMin: minimum jet pT\n"
      << width << "  ETMin: minimum jet eta\n"
      << width << "  DR: jet distance parameter\n"
      << width << "  f: Siscone f parameter (default: 0.75)\n"
      << width << "  EtaMax: maximum jet eta (default: None)\n"
      << width << "  YMax: maximum jet rapidity (default: None)\n";
}
