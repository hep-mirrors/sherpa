#ifndef PHASIC_Selectors_MinSelector_h
#define PHASIC_Selectors_MinSelector_h

#include "ATOOLS/Phys/Particle_List.H"
#include "PHASIC++/Selectors/Selector.H"
#include "ATOOLS/Math/Poincare.H"

namespace PHASIC {
  class MinSelector : public Selector_Base {
    std::vector<Selector_Base *> m_sels; 
  public:
    MinSelector(const Selector_Key &key);

    ~MinSelector();


    bool   NoJetTrigger(const ATOOLS::Vec4D_Vector &);
    bool   Trigger(const ATOOLS::Vec4D_Vector &);
    bool   JetTrigger(const ATOOLS::Vec4D_Vector &,
		      ATOOLS::NLO_subevtlist *const);

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
#include "ATOOLS/Org/Data_Reader.H"


using namespace PHASIC;
using namespace ATOOLS;

/*---------------------------------------------------------------------

  General form - flavours etc are unknown, will operate on a Particle_List.

  --------------------------------------------------------------------- */

MinSelector::MinSelector(const Selector_Key &key) : 
  Selector_Base("MinSelector")
{
  for (size_t k=0;k<key.size();++k) {
    Selector_Key subkey(key.p_proc,key.p_read);
    subkey.push_back(std::vector<std::string>(key[k].size()-1));
    for (size_t j=1;j<key[k].size();++j) subkey.back()[j-1]=key[k][j];

    subkey.m_key=key[k][0];
    Selector_Base *sel(Selector_Getter::GetObject(key[k][0],subkey));
    if (sel!=NULL) m_sels.push_back(sel);
  }
  m_sel_log    = new Selector_Log(m_name);
}


MinSelector::~MinSelector() {
  while (m_sels.size()>0) {
    delete *m_sels.begin();
    m_sels.erase(m_sels.begin());
  }
}


bool MinSelector::Trigger(const Vec4D_Vector &p)
{
  for (size_t k=0;k<m_sels.size();++k) {
    if (m_sels[k]->Trigger(p)) {
      m_sel_log->Hit(0);
      return 1;
    }
  }
  m_sel_log->Hit(1);
  return 0;
}

bool MinSelector::JetTrigger(const Vec4D_Vector &p,NLO_subevtlist *const subs)
{
  for (size_t k=0;k<m_sels.size();++k) {
    if (m_sels[k]->JetTrigger(p,subs)) {
      m_sel_log->Hit(0);
      return 1;
    }
  }
  m_sel_log->Hit(1);
  return 0;
}

bool MinSelector::NoJetTrigger(const Vec4D_Vector & mom)
{
  for (size_t k=0;k<m_sels.size();++k) {
    if (m_sels[k]->NoJetTrigger(mom)) {
      m_sel_log->Hit(0);
      return 1;
    }
  }
  m_sel_log->Hit(1);
  return 0;
}

void MinSelector::BuildCuts(Cut_Data * cuts) 
{
  return;
}

DECLARE_ND_GETTER(MinSelector,"MinSelector",Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,MinSelector>::
operator()(const Selector_Key &key) const
{
  MinSelector *msel(new MinSelector(key));
  return msel;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,MinSelector>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"MinSelector {\n"
     <<"                          Selector 1\n"
     <<"                          Selector 2\n"
     <<"                          ...\n"
     <<"                        }"; 
}
