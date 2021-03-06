#ifndef PHASIC_Selectors_MinSelector_h
#define PHASIC_Selectors_MinSelector_h

#include "ATOOLS/Phys/Particle_List.H"
#include "PHASIC++/Selectors/Selector.H"
#include "ATOOLS/Math/Poincare.H"

namespace PHASIC {
  class MinSelector : public Selector_Base {
  public:
    MinSelector(const Selector_Key &key);

    ~MinSelector();

    bool   Trigger(ATOOLS::Selector_List &);

    void   BuildCuts(Cut_Data *);
  };
}

#endif

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"


using namespace PHASIC;
using namespace ATOOLS;

/*---------------------------------------------------------------------

  General form - flavours etc are unknown, will operate on a Particle_List.

  --------------------------------------------------------------------- */

MinSelector::MinSelector(const Selector_Key &key) : 
  Selector_Base("MinSelector",key.p_proc)
{
  DEBUG_FUNC("");
  ReadInSubSelectors(key);
}


MinSelector::~MinSelector() {
  while (m_sels.size()>0) {
    delete *m_sels.begin();
    m_sels.erase(m_sels.begin());
  }
}

bool MinSelector::Trigger(Selector_List &sl)
{
  for (size_t k=0;k<m_sels.size();++k) {
    if (m_sels[k]->Trigger(sl)) {
      m_sel_log->Hit(0);
      return 1;
    }
  }
  m_sel_log->Hit(1);
  return 0;
}

void MinSelector::BuildCuts(Cut_Data * cuts) 
{
  // maybe be smarter here and take minimum of cuts build in subselectors
  return;
}

DECLARE_GETTER(MinSelector,"MinSelector",Selector_Base,Selector_Key);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,MinSelector>::
operator()(const Selector_Key &key) const
{
  MinSelector *msel(new MinSelector(key));
  return msel;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,MinSelector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  std::string w(width+4,' ');
  str<<"MinSelector:\n"
     <<w<<"  Subselectors: [...]";
}
