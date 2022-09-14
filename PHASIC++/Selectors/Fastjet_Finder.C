#include "PHASIC++/Selectors/Fastjet_Selector_Base.H"

namespace PHASIC {
  class Fastjet_Finder : public Fastjet_Selector_Base {
    int m_nb, m_nb2;

  public:
    Fastjet_Finder(Process_Base* const proc, ATOOLS::Scoped_Settings s,
                   int nb, int nb2);
    bool Trigger(ATOOLS::Selector_List&);
    void BuildCuts(Cut_Data* cuts) {
      cuts->smin= ATOOLS::Max(cuts->smin,m_smin);
    }
  };
}

#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace PHASIC;
using namespace ATOOLS;


/*---------------------------------------------------------------------

  General form - flavours etc are unknown, will operate on a Particle_List.

  --------------------------------------------------------------------- */

Fastjet_Finder::Fastjet_Finder(Process_Base* const proc, Scoped_Settings s,
                               int nb, int nb2):
  Fastjet_Selector_Base("FastjetFinder", proc, s),
  m_nb(nb), m_nb2(nb2)
{}


bool Fastjet_Finder::Trigger(Selector_List &sl)
{
  if (mpi->Rank()) fjcore::ClusterSequence::set_fastjet_banner_stream(NULL);
  bool ee(rpa->gen.Beam1().IsLepton() && rpa->gen.Beam2().IsLepton());

  DEBUG_FUNC((p_proc?p_proc->Flavours():Flavour_Vector()));

  std::vector<ATOOLS::Vec4D> p(sl.size());
  for (size_t i(0);i<p.size();++i) p[i]=sl[i].Momentum();
  if (sl[0].Flavour().IsLepton()&&rpa->gen.Beam2().IsHadron()) {
    msg_Debugging()<<METHOD<<"(): Boost to Breit frame {\n";
    Vec4D pp(rpa->gen.PBeam(1)), qq(p[0]-p[2]);
    Poincare cms(pp+qq);
    double Q2(-qq.Abs2()), x(Q2/(2.0*pp*qq)), E(sqrt(Q2)/(2.0*x));
    Vec4D P(sqrt(E*E+pp.Abs2()),0.0,0.0,-E);
    Vec4D q(0.0,0.0,0.0,2.0*x*E);
    cms.Boost(pp);
    cms.Boost(qq);
    Poincare zrot(pp,-Vec4D::ZVEC);
    zrot.Rotate(pp);
    zrot.Rotate(qq);
    Poincare breit(P+q);
    breit.BoostBack(pp);
    breit.BoostBack(qq);
    if (!IsEqual(pp,P,1.0e-3) || !IsEqual(qq,q,1.0e-3))
      msg_Error()<<METHOD<<"(): Boost error."<<std::endl;
    for (int i(0);i<p.size();++i) {
      msg_Debugging()<<"  "<<i<<": "<<p[i];
      cms.Boost(p[i]);
      zrot.Rotate(p[i]);
      breit.BoostBack(p[i]);
      msg_Debugging()<<" -> "<<p[i]<<"\n";
    }
    msg_Debugging()<<"}\n";
  }

  std::vector<fjcore::PseudoJet> input,jets;
  for (size_t i(m_nin);i<sl.size();++i) {
    if (ToBeClustered(sl[i].Flavour(), (m_nb>0 || m_nb2>0))) {
      input.push_back(MakePseudoJet(sl[i].Flavour(), p[i]));
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
  if (nb<m_nb)   trigger=false;
  if (nb2<m_nb2) trigger=false;

  if (!trigger) {
    msg_Debugging()<<"Point discarded by jet finder"<<std::endl;
  } else {
    msg_Debugging()<<"Point passed"<<std::endl;
  }

  return (1-m_sel_log->Hit(1-trigger));
}


DECLARE_GETTER(Fastjet_Finder,"FastjetFinder",Selector_Base,Selector_Key);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,Fastjet_Finder>::
operator()(const Selector_Key &key) const
{
  auto s = key.m_settings["FastjetFinder"];

  // b tagging
  const auto nb  = s["Nb"] .SetDefault(-1).Get<int>();
  const auto nb2 = s["Nb2"].SetDefault(-1).Get<int>();

  return new Fastjet_Finder(key.p_proc, s, nb, nb2);
}

void ATOOLS::Getter<Selector_Base,Selector_Key,Fastjet_Finder>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"FastjetFinder:\n";
  Fastjet_Selector_Base::PrintCommonInfoLines(str, width);
  str<<width<<"  Nb: number of jets with b quarks\n"
     <<width<<"  Nb2: number of jets with non-vanishing b content";
}

