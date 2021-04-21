#include "ATOOLS/Phys/Particle_Dresser.H"
#include "PHASIC++/Selectors/Selector.H"
#include "ATOOLS/Org/MyStrStream.H"

namespace PHASIC {
  class Process_Base;

  class DressedParticleSelector : public Selector_Base {
    ATOOLS::Particle_Dresser *               p_dresser;
  public:
    DressedParticleSelector(const Selector_Key &key);

    ~DressedParticleSelector();


    bool   Trigger(ATOOLS::Selector_List &);

    void   BuildCuts(Cut_Data *);
  };
}

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"


using namespace PHASIC;
using namespace ATOOLS;

DressedParticleSelector::DressedParticleSelector(const Selector_Key &key) :
  Selector_Base("DressedParticleSelector",key.p_proc), p_dresser(NULL)
{
  DEBUG_FUNC("");
  auto s = key.m_settings["DressedParticleSelector"];
  const auto algoparams = s["DressingAlgorithm"]
    .SetDefault<std::string>({})
    .GetVector<std::string>();
  if (algoparams.size() != 2)
    THROW(fatal_error, "DressingAlgorithm needs two values");
  const auto algo = algoparams[0];
  const auto dR = ToType<double>(algoparams[1]);
  p_dresser = new Particle_Dresser(algo,dR);
  if (!p_dresser) THROW(fatal_error, "Invalid dressing algorithm");
  const auto flavradiusparams = s["FlavourDependentRadius"]
    .SetDefault<std::string>({})
    .GetVector<std::string>();
  if (!flavradiusparams.empty()) {
    if (flavradiusparams.size() != 2)
      THROW(fatal_error, "FlavourDependentRadius needs two values");
    const auto kf = ToType<kf_code>(flavradiusparams[0]);
    const auto fDR = ToType<double>(flavradiusparams[1]);
    p_dresser->SetFlavourDependentCone(kf, fDR);
  }
  p_dresser->CompleteConeLists();
  ReadInSubSelectors(key);
  if (msg_LevelIsDebugging()) {
    msg_Out()<<"Additional Selectors:\n";
    for (size_t i(0);i<m_sels.size();++i)
      msg_Out()<<"  "<<m_sels[i]->Name()<<std::endl;
  }
}

DressedParticleSelector::~DressedParticleSelector() {
  if (p_dresser) delete p_dresser;
  while (m_sels.size()>0) {
    delete *m_sels.begin();
    m_sels.erase(m_sels.begin());
  }
}

bool DressedParticleSelector::Trigger(Selector_List &sl)
{
  DEBUG_FUNC((p_proc?p_proc->Flavours():Flavour_Vector()));
  p_dresser->Dress(sl);
  for (size_t k=0;k<m_sels.size();++k) {
    if (!m_sels[k]->Trigger(sl)) {
      msg_Debugging()<<"Point discarded"<<std::endl;
      m_sel_log->Hit(true);
      return false;
    }
  }
  msg_Debugging()<<"Point passed"<<std::endl;
  m_sel_log->Hit(false);
  return true;
}

void DressedParticleSelector::BuildCuts(Cut_Data * cuts)
{
  for (size_t i(0);i<m_sels.size();++i) m_sels[i]->BuildCuts(cuts);
}

DECLARE_GETTER(DressedParticleSelector, "DressedParticleSelector",
               Selector_Base, Selector_Key);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,
                              DressedParticleSelector>::operator()
(const Selector_Key &key) const
{
  DressedParticleSelector *dpsel(new DressedParticleSelector(key));
  return dpsel;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,DressedParticleSelector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  std::string w(width+4,' ');
  str<<"DressedParticleSelector:\n"
     <<w<<"  DressingAlgorithm: [<Cone|Recombination>, <dR>, <exp>]  # exp is optional\n"
     <<w<<"  # optional settings:\n"
     <<w<<"  FlavourDependentCone: [<kf>, <dR>]\n"
     <<w<<"  Subselectors: [ ... ]";
}
