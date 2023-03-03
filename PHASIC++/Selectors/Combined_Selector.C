#include "PHASIC++/Selectors/Combined_Selector.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Phys/Selector_List.H"

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

Combined_Selector::Combined_Selector(Process_Base *const proc):
  Selector_Base("Combined_Selector",proc), m_count(0), m_on(1), m_res(0)
{
}

Combined_Selector::~Combined_Selector()
{
  while (m_sels.size()>0) {
    delete *m_sels.begin();
    m_sels.erase(m_sels.begin());
  }
}

bool Combined_Selector::Initialize(const Selector_Key &key)
{
  for (auto s : key.GetSelectors()) {
    std::string name;
    if (s.IsList()) {
      name = s.SetDefault<std::string>({}).GetVector<std::string>()[0];
    } else {
      const auto keys = s.GetKeys();
      if (keys.size() != 1) {
        THROW(fatal_error, std::string("Each selector mapping must have ")
            + "exactly one key-value pair, where the key gives the selector "
            + "type.");
      }
      name = keys[0];
    }
    Selector_Key subkey;
    subkey.m_settings = s;
    subkey.p_proc = p_proc;
    auto* sel = Selector_Getter::GetObject(name, subkey);
    if (sel) {
      m_sels.push_back(sel);
      msg_Debugging() << "new Selector_Base(\"" << name << "\")\n";
    } else {
      msg_Out()<<endl;
      THROW(fatal_error, "Did not find selector \"" + name + "\".");
    }
  }
  return true;
}

bool Combined_Selector::Trigger(const Vec4D_Vector &p,
                                const Flavour *fl, size_t n)
{
  DEBUG_FUNC(p.size()<<" momenta, "<<n<<" flavours");
  Selector_List sl(n?Selector_List(fl,n,p,m_nin):
                     Selector_List(p_proc->Flavours(),p,m_nin));
  return Trigger(sl);
}

bool Combined_Selector::Trigger(Selector_List& sl)
{
  DEBUG_FUNC("");
  // BEWARE: Selector_List will be modified
  m_res=1;
  if (!m_on) return m_res;
  for (size_t i=0; i<m_sels.size(); ++i) {
    msg_Debugging()<<m_sels[i]->Name()<<std::endl;
    if (!(m_sels[i]->Trigger(sl))) {
      msg_Debugging()<<"Point discarded"<<std::endl;
      m_res=0;
      return m_res;
    }
  }
  msg_Debugging()<<"Point passed"<<std::endl;
  return m_res;
}

bool Combined_Selector::RSTrigger(NLO_subevtlist *const subs)
{
  int pass(0);
  for (size_t n(0);n<subs->size();++n) {
    p_sub=(*subs)[n];
    Vec4D_Vector mom(p_sub->p_mom,&p_sub->p_mom[p_sub->m_n]);
    for (size_t i(0);i<m_nin;++i)
      if (mom[i][0]<0.0) mom[i]=-mom[i];
    Selector_List sl=Selector_List
      (p_sub->p_fl,p_sub->m_n,mom,m_nin);
    p_sub->m_trig=Trigger(sl);
    if (p_sub->m_trig) pass=1;
    p_sub=NULL;
  }
  m_rsres=pass;
  return pass;
}

bool Combined_Selector::Pass() const
{
  for (size_t i=0;i<m_sels.size();++i)
    if (!m_sels[i]->Pass()) return false;
  return true;
}

void Combined_Selector::BuildCuts(Cut_Data * cuts)
{
  for (size_t i=0; i<m_sels.size(); ++i) m_sels[i]->BuildCuts(cuts);

  for (size_t i=0; i<m_osc.size(); ++i) cuts->Setscut(m_osc[i].first,m_osc[i].second);
  cuts->Complete();
  for (size_t i=0; i<m_osc.size(); ++i) cuts->Setscut(m_osc[i].first,m_osc[i].second);
}

void Combined_Selector::AddOnshellCondition(size_t s,double d)
{
  m_osc.push_back(std::pair<size_t,double>(s,d));
}

void Combined_Selector::Output()
{
  msg_Debugging()<<"========================================="<<std::endl
                 <<"Efficiency of the Selector : "<<m_name<<std::endl;
  for (size_t i=0; i<m_sels.size(); ++i) m_sels[i]->Output();
  msg_Debugging()<<"========================================="<<std::endl;
}

Selector_Base * Combined_Selector::GetSelector(const std::string &name) const
{
  for (size_t i=0; i<m_sels.size(); ++i) 
    if (m_sels[i]->Name()==name) return m_sels[i];
  return 0;
  
}

void Combined_Selector::ListSelectors() const
{
  msg_Info()<<"Selectors:"<<std::endl;
  for (size_t i=0; i<m_sels.size(); ++i)
    msg_Info()<<m_sels[i]->Name()<<std::endl;
}

std::vector<Weights_Map> Combined_Selector::CombinedResults() const
{
  std::vector<Weights_Map> res = {Weights_Map{}};
  for (auto& sel : m_sels) {
    std::vector<Weights_Map> other = sel->Results();
    if (other.size() == 1) {
      for (auto& weights : res) {
        weights *= other[0];
      }
    } else if (res.size() == 1) {
      Weights_Map currentweights = res[0];
      res = other;
      for (auto& weights : res) {
        weights *= currentweights;
      }
    } else {
      assert(res.size() == other.size());
      for (int i {0}; i < res.size(); ++i) {
        res[i] *= other[i];
      }
    }
  }
  return res;
}
