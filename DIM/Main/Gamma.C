#include "DIM/Main/Gamma.H"

#include "DIM/Main/MCatNLO.H"
#include "DIM/Shower/Kernel.H"
#include "PHASIC++/Process/MCatNLO_Process.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PDF/Main/Jet_Criterion.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"

using namespace DIM;
using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;

#define DEBUG__Trial_Weight
// #define DEBUG__Differential

Gamma::Gamma(MCatNLO *const dire,Shower *const shower):
  p_dire(dire), p_shower(shower)
{
}

Weight_Value Gamma::Differential
(Cluster_Amplitude *const ampl,const nlo_type::code type,
 const std::string add) const
{
#ifndef DEBUG__Differential
  int olv(msg->Level());
  msg->SetLevel(2);
#endif
  NLOTypeStringProcessMap_Map *procs
    (ampl->Procs<NLOTypeStringProcessMap_Map>());
  Process_Base::SortFlavours(ampl);
  std::string pname(Process_Base::GenerateName(ampl));
  StringProcess_Map::const_iterator pit((*(*procs)[type]).find(pname+add));
  if (pit==(*(*procs)[type]).end()) return Weight_Value();
  Weight_Value meps(pit->second);
  bool kon(pit->second->KFactorSetter(true)->On());
  pit->second->KFactorSetter(true)->SetOn(false);
  meps.m_b = meps.m_me = static_cast<double>(
      pit->second->Differential(*ampl, Weight_Type::nominal, 1 | 2 | 4));
  pit->second->KFactorSetter(true)->SetOn(kon);
  meps.m_me*=pit->second->SymFac();
  meps.m_muf2[0]=ampl->MuF2(0);
  meps.m_muf2[1]=ampl->MuF2(1);
  meps.m_mur2=ampl->MuR2();
  meps.m_muq2=ampl->MuQ2();
#ifndef DEBUG__Differential
  msg->SetLevel(olv);
#endif
  return meps;
}

Weight_Map Gamma::CalculateWeight(Cluster_Amplitude *const ampl)
{
  Splitting s(p_shower->LastSplitting());
  s.m_clu=2;
  s.m_kfac=0;
  s.m_t1=ampl->MuR2();
  Cluster_Leg *li(ampl->IdLeg(1<<(s.p_c->Id()-1)));
  Cluster_Leg *lk(ampl->IdLeg(1<<(s.p_s->Id()-1)));
  Cluster_Leg *lj(ampl->IdLeg(ampl->IdNew()));
#ifdef DEBUG__Trial_Weight
  DEBUG_FUNC(ID(li->Id())<<","<<ID(lj->Id())<<"<->"<<ID(lk->Id()));
  msg_Debugging()<<*ampl<<"\n";
#endif
  Cluster_Amplitude *bampl(p_dire->GetBornAmplitude());
  const Kernel *cdip(s.p_sk);
#ifdef DEBUG__Trial_Weight
  msg_Debugging()<<"B config -> "<<*bampl<<" -> "<<s<<" ( "
		 <<cdip->LF()->Flav(0)<<" -> "<<cdip->LF()->Flav(1)
		 <<" "<<cdip->LF()->Flav(2)<<" )\n";
#endif
  Weight_Value meps(Differential(bampl));
  if (meps.p_proc==NULL) return Weight_Map();
  meps.p_sf=cdip;
  meps.m_me/=cdip->LF()->AsymmetryFactor(s);
#ifdef DEBUG__Trial_Weight
  double me=meps.m_me;
#endif
  meps.m_me*=cdip->Value(s)*cdip->LF()->MEPSWeight(s);
  if (meps.m_me==0.0) {
#ifdef DEBUG__Trial_Weight
    msg_Debugging()<<"zero matrix element\n";
#endif
    return Weight_Map();
  }
#ifdef DEBUG__Trial_Weight
  msg_Debugging()<<"add ( x = "<<s.m_x<<", y = "<<s. m_y<<", z = "
		 <<s.m_z<<", kt = "<<sqrt(s.m_t)<<" ) {\n  "<<*li
		 <<"\n  "<<*lj<<"\n  "<<*lk<<"\n} -> w = "
		 <<me<<" * "<<meps.m_me/me<<" -> "<<meps.m_me
		 <<" ( S = "<<cdip->LF()->AsymmetryFactor(s)<<" )\n";
#endif
  Weight_Map ws;
  ws[Weight_Key(li->Id()|lj->Id(),lk->Id())]=meps;
  return ws;
}

MC_Weight Gamma::TrialWeight(Cluster_Amplitude *const ampl)
{
  DEBUG_FUNC("");
  p_ms=ampl->MS();
  p_shower->SetMS(p_ms);
  Weight_Map ws(CalculateWeight(ampl));
  if (ws.empty()) return MC_Weight(0.0,1.0,1.0);
  const Splitting &s(p_shower->LastSplitting());
  size_t idij(0), idk(0);
  double wgt(0.0);
  Weight_Value wact;
  Weight_Map::const_iterator ait;
#ifdef DEBUG__Trial_Weight
  msg_Debugging()<<"Accumulate weights {\n";
#endif
  for (Weight_Map::const_iterator
	 wit(ws.begin());wit!=ws.end();++wit) {
#ifdef DEBUG__Trial_Weight
    msg_Debugging()<<"  "<<wit->first<<" -> "<<wit->second;
#endif
    wgt+=wit->second.m_me;
    if ((wit->first.m_ij==((1<<(s.p_c->Id()-1))|ampl->IdNew())) &&
	(wit->first.m_k==(1<<(s.p_s->Id()-1)))) {
      ait=wit;
      wact=ait->second;
      idij=wit->first.m_ij;
      idk=wit->first.m_k;
#ifdef DEBUG__Trial_Weight
      msg_Debugging()<<" <- active";
#endif
    }
#ifdef DEBUG__Trial_Weight
    msg_Debugging()<<"\n";
#endif
  }
#ifdef DEBUG__Trial_Weight
  msg_Debugging()<<"} -> w = "<<wgt<<"\n";
#endif
  if (!wact.p_sf || wact.m_me==-1.0)
    THROW(fatal_error,"No active splitting weight");
  ampl->SetMuF2(0,wact.m_muf2[0]);
  ampl->SetMuF2(1,wact.m_muf2[1]);
  ampl->SetMuR2(wact.m_mur2);
  ampl->SetKT2(ampl->MuQ2());
  int i(-1), j(-1), k(-1);
  for (size_t l(0);l<ampl->Legs().size();++l)
    if (ampl->Leg(l)->Id()&idk) k=l;
    else if (ampl->Leg(l)->Id()&idij) {
      if (i<0) i=l;
      else j=l;
    }
  std::string nadd("__QCD(S)_RS");
  nadd+=ToString(i)+"_"+ToString(j)+"_"+ToString(k);
  double rme(Differential(ampl,nlo_type::rsub,nadd).m_me);
  msg_Debugging()<<"me / ecss = "<<rme<<" / "<<wact.m_me
		 <<" = "<<rme/wact.m_me<<"\n";
  double h(wact.m_me), g(p_shower->OEF()*rme);
  g*=Max(1.0,h/dabs(rme));
  if (IsEqual(rme,h,1.0e-6) || rme==0.0) g=h;
  return MC_Weight(rme,g,h);
}

bool Gamma::Reject()
{
  DEBUG_FUNC("");
  if (p_dire->PSMode()) {
    m_weight=1.0;
    return false;
  }
  Cluster_Amplitude *rampl=p_dire->GetRealEmissionAmplitude(1);
  MC_Weight wgt(TrialWeight(rampl));
  rampl->Delete();
  if (wgt.MC()>ran->Get()) {
    m_weight=wgt.Accept();
    msg_Debugging()<<"w = "<<wgt.MC()<<" "<<wgt<<" -> accept\n";
    return false;
  }
  m_weight=wgt.Reject();
  msg_Debugging()<<"w = "<<wgt.MC()<<" "<<wgt<<" -> reject\n";
  return true;
}

namespace DIM {

  std::ostream &operator<<(std::ostream &str,const Weight_Key &k)
  {
    return str<<"["<<ATOOLS::ID(k.m_ij)<<","<<ATOOLS::ID(k.m_k)<<"]";
  }

  std::ostream &operator<<(std::ostream &str,const Weight_Value &w)
  {
    return str<<w.m_me<<"  "<<w.p_proc->Name()<<" [ "
	      <<w.p_sf->LF()->Flav(0)<<" -> "<<w.p_sf->LF()->Flav(1)
	      <<" "<<w.p_sf->LF()->Flav(2)<<" ] ( \\mu_F = {"
	      <<sqrt(w.m_muf2[0])<<","<<sqrt(w.m_muf2[1])
	      <<"}, \\mu_R = "<<sqrt(w.m_mur2)<<" ) ";
  }

}
