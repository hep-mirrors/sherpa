#ifndef PHASIC_Selectors_DressedParticleSelector_h
#define PHASIC_Selectors_DressedParticleSelector_h

#include "ATOOLS/Phys/Particle_Dresser.H"
#include "PHASIC++/Selectors/Selector.H"
#include "ATOOLS/Org/MyStrStream.H"

namespace PHASIC {
  class Process_Base;

  class DressedParticleSelector : public Selector_Base {
    ATOOLS::Particle_Dresser *               p_dresser;
    std::vector<ATOOLS::Particle_Dresser * > m_dressers;
  public:
    DressedParticleSelector(const Selector_Key &key);

    ~DressedParticleSelector();


    bool   Trigger(const ATOOLS::Vec4D_Vector &,
                   ATOOLS::NLO_subevt *const=NULL);

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

DressedParticleSelector::DressedParticleSelector(const Selector_Key &key) :
  Selector_Base("DressedParticleSelector",key.p_proc), p_dresser(NULL)
{
  DEBUG_FUNC("");
  NLO_subevtlist * subs(p_proc->GetSubevtList());
  if (subs) m_dressers.resize(subs->size(),NULL);
  double dR(0.),exp(1.);
  for (size_t k=0;k<key.size();++k) {
    if (key[k].size()>2 && key[k][0]=="DressingAlgorithm") {
      if (p_dresser) THROW(fatal_error,"Too many dressing algorithms.");
      std::string algo(ToType<std::string>(key[k][1]));
      dR=ToType<double>(key.p_read->Interpreter()->Interprete(key[k][2]));
      if (key[k].size()>3)
        exp=ToType<double>(key.p_read->Interpreter()->Interprete(key[k][3]));
      p_dresser = new Particle_Dresser(p_fl,m_nin,m_nout,algo,dR,exp);
      for (size_t i(0);i<m_dressers.size();++i) {
        NLO_subevt * sub((*subs)[i]);
        m_dressers[i] = new Particle_Dresser(sub->p_fl,m_nin,sub->m_n-m_nin,
                                            algo,dR,exp);
      }
    }
    else if (key[k].size()>2 && key[k][0]=="FlavourDependentCone") {
      if (!p_dresser) THROW(fatal_error,"No dressing algorithm set yet.");
      kf_code kf=ToType<kf_code>(key.p_read->Interpreter()->Interprete(key[k][1]));
      dR=ToType<double>(key.p_read->Interpreter()->Interprete(key[k][2]));
      p_dresser->SetFlavourDependentCone(kf,dR);
      for (size_t i(0);i<m_dressers.size();++i)
        m_dressers[i]->SetFlavourDependentCone(kf,dR);
    }
    else {
      if (!p_dresser) THROW(fatal_error,"No dressing algorithm defined.");
      p_dresser->CompleteConeLists();
      for (size_t i(0);i<m_dressers.size();++i)
        m_dressers[i]->CompleteConeLists();
      ReadInSubSelectors(key,k);
      break;
    }
  }
  if (msg_LevelIsDebugging()) {
    msg_Out()<<"Additional Selectors:\n";
    for (size_t i(0);i<m_sels.size();++i)
      msg_Out()<<"  "<<m_sels[i]->Name()<<std::endl;
  }
}

DressedParticleSelector::~DressedParticleSelector() {
  if (p_dresser) delete p_dresser;
  for (size_t i(0);i<m_dressers.size();++i) delete m_dressers[i];
  m_dressers.resize(0);
  while (m_sels.size()>0) {
    delete *m_sels.begin();
    m_sels.erase(m_sels.begin());
  }
}

bool DressedParticleSelector::Trigger(const Vec4D_Vector &p,
                                      NLO_subevt *const sub)
{
  DEBUG_FUNC((p_proc?p_proc->Flavours():Flavour_Vector()));
  Particle_Dresser * dresser(sub?m_dressers[sub->m_idx]:p_dresser);
  Vec4D_Vector pd(dresser->Dress(p));
  if (msg_LevelIsIODebugging()) {
    for (size_t i(0);i<p.size();++i)
      msg_Debugging()<<i<<": "<<p[i]<<" -> "<<pd[i]<<std::endl;
  }
  for (size_t k=0;k<m_sels.size();++k) {
    if (!m_sels[k]->Trigger(pd,sub)) {
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

DECLARE_ND_GETTER(DressedParticleSelector,"DressedParticleSelector",
                  Selector_Base,Selector_Key,true);

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
  str<<"DressedParticleSelector {\n"
     <<w<<"  DressingAlgorithm <Cone|Recombination> <dR> [<exp>]\n"
     <<w<<"  [FlavourDependentCone <kf> <dR>]\n"
     <<w<<"  Selector 1\n"
     <<w<<"  Selector 2\n"
     <<w<<"  ...\n"
     <<w<<"}";
}
