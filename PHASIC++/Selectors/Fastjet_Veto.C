#include "PHASIC++/Selectors/Fastjet_Selector_Base.H"

namespace PHASIC {
  class Fastjet_Veto : public Fastjet_Selector_Base {
    int m_nb, m_nb2;

  public:
    Fastjet_Veto(Process_Base* const proc, ATOOLS::Scoped_Settings s,
                 int nb, int nb2);
    bool Trigger(ATOOLS::Selector_List&);
    void BuildCuts(Cut_Data* cuts) {
      cuts->smin= ATOOLS::Max(cuts->smin,m_smin);
    }
  };
}

#include "PHASIC++/Main/Process_Integrator.H"

using namespace PHASIC;
using namespace ATOOLS;


/*---------------------------------------------------------------------

  General form - flavours etc are unknown, will operate on a Particle_List.

  --------------------------------------------------------------------- */

Fastjet_Veto::Fastjet_Veto(Process_Base* const proc, Scoped_Settings s,
                           int nb, int nb2):
  Fastjet_Selector_Base("Fastjetfinder", proc, s),
  m_nb(nb), m_nb2(nb2)
{
}


bool Fastjet_Veto::Trigger(Selector_List &sl)
{
  DEBUG_FUNC((p_proc?p_proc->Flavours():Flavour_Vector()));
  std::vector<fjcore::PseudoJet> input,jets;
  for (size_t i(m_nin);i<sl.size();++i) {
    if (ToBeClustered(sl[i].Flavour(), (m_nb>0 || m_nb2>0))) {
      input.push_back(MakePseudoJet(sl[i].Flavour(), sl[i].Momentum()));
    }
  }

  fjcore::ClusterSequence cs(input,*p_jdef);
  jets=fjcore::sorted_by_pt(cs.inclusive_jets());
  msg_Debugging()<<"njets(ini)="<<jets.size()<<std::endl;

  if (m_eekt) {
    int n(0);
    for (size_t i(0);i<input.size();++i)
      if (cs.exclusive_dmerge_max(i)>sqr(m_ptmin)) ++n;
    return (1-m_sel_log->Hit(1-(n>=m_nj)));
  }

  int n(0), nb(0), nb2(0);
  for (size_t i(0);i<jets.size();++i) {
    Vec4D pj(jets[i].E(),jets[i].px(),jets[i].py(),jets[i].pz());
    msg_Debugging()<<"Jet "<<i<<": pT="<<pj.PPerp()<<", |eta|="<<dabs(pj.Eta())
                   <<", |y|="<<dabs(pj.Y())<<std::endl;
    if (pj.PPerp() > m_ptmin
        && pj.EPerp() > m_etmin
        && dabs(pj.Eta()) < m_eta
        && dabs(pj.Y()) < m_y) {
      n++;
      if (BTag(jets[i], 1)) nb++;
      if (BTag(jets[i], 2)) nb2++;
    }
  }
  msg_Debugging()<<"njets(fin)="<<n<<std::endl;

  bool trigger(true);
  if (n<m_nj)   trigger=false;
  if (nb<m_nb) trigger=false;
  if (nb2<m_nb2) trigger=false;
  trigger=!trigger;

  if (!trigger) {
    msg_Debugging()<<"Point discarded by jet veto"<<std::endl;
  } else {
    msg_Debugging()<<"Point passed"<<std::endl;
  }
  return (1-m_sel_log->Hit(1-trigger));
}


DECLARE_GETTER(Fastjet_Veto,"FastjetVeto",Selector_Base,Selector_Key);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,Fastjet_Veto>::
operator()(const Selector_Key &key) const
{
  auto s = key.m_settings["FastjetVeto"];

  // b tagging
  const auto nb  = s["Nb"] .SetDefault(-1).Get<int>();
  const auto nb2 = s["Nb2"].SetDefault(-1).Get<int>();

  return new Fastjet_Veto(key.p_proc,s,nb,nb2);
}

void ATOOLS::Getter<Selector_Base,Selector_Key,Fastjet_Veto>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"FastjetVeto:\n";
  Fastjet_Selector_Base::PrintCommonInfoLines(str, width);
  str<<width<<"  Nb: number of jets with b quarks\n"
     <<width<<"  Nb2: number of jets with non-vanishing b content";
}
