#include "ATOOLS/Org/Message.H"
#include "AddOns/GoSam/GoSam_Born.H"
#include "AddOns/GoSam/GoSam_Interface.H"

#include <algorithm>

using namespace GoSam;
using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

GoSam_Born::GoSam_Born(const PHASIC::External_ME_Args& args,
                       bool active, int proc_id) :
Tree_ME2_Base(args),
  m_pn(args.m_inflavs.size()+args.m_outflavs.size()),
  m_nin(args.m_inflavs.size()), m_active(active),
  m_proc_id(proc_id)
{
  Flavour_Vector flavs = args.m_inflavs;
  for(Flavour fl: args.m_outflavs) flavs.push_back(fl);
  m_Norm=Flavour::FSSymmetryFactor(args.m_outflavs)*Flavour::ISSymmetryFactor(args.m_inflavs);
  m_oqcd = args.m_orders[0];
  m_oqed = args.m_orders[1];
  if(!m_active) return;

  m_g2qed = 4.0*M_PI*MODEL::s_model->ScalarConstant("alpha_QED");
  m_g2qcd = 4.0*M_PI*MODEL::s_model->ScalarConstant("alpha_S");

  m_momenta = vector<double>(m_pn*5);
  for (size_t i=0;i<m_pn;i++) m_momenta[4+i*5]=flavs[i].Mass();
  for (size_t i=0;i<3;i++) m_result[i]=0.;
  m_result[3]=1.;

#ifdef VALIDATE_GOSAM_INTERFACE
  auto s = Settings::GetMainSettings()["GOSAM"];
  m_errorOnFailedCheck = s["ErrorOnFailedCheck"].SetDefault(false).Get<bool>();
  m_printAllChecks = s["PrintAllChecks"].SetDefault(false).Get<bool>();
#endif
}

double GoSam_Born::Calc(const Vec4D_Vector& pp) {
  DEBUG_FUNC(m_proc_id);
  if (!m_active) return 0;
  if (m_proc_id<0) return 1;

  // TODO: Do we need CMS boost?
  Vec4D_Vector momenta(pp);
  Poincare cms(momenta[0]+momenta[1]);
  for (size_t i(0);i<momenta.size();++i) cms.Boost(momenta[i]);

  for (size_t i=0;i<m_pn;i++) {
    m_momenta[0+i*5]=momenta[i][0];
    m_momenta[1+i*5]=momenta[i][1];
    m_momenta[2+i*5]=momenta[i][2];
    m_momenta[3+i*5]=momenta[i][3];
  }
  double scale = p_aqcd->Scale();
  double acc;
  int ierr = 1;
  GoSam::GoSam_Interface::SetParameter("alphaS", AlphaQCD(), ierr);
  GoSam::GoSam_Interface::SetParameter("alpha", AlphaQED(), ierr);
  GoSam::OLP_EvalSubProcess2(&m_proc_id, m_momenta.data(), &scale, m_result.data(), &acc);
  // born, either tree or loop induced, one of them should always be zero
  double born = m_Norm*(m_result[2]+m_result[3]);
#ifdef VALIDATE_GOSAM_INTERFACE
  if(p_check) {
    MODEL::Coupling_Map cpls;
    cpls.insert({"Alpha_QCD",p_aqcd});
    cpls.insert({"Alpha_QED",p_aqed});
    p_check->SetCouplings(cpls);
    double checkBorn = p_check->Calc(pp);
    if(!IsEqual(born,checkBorn,1e-2)) {
      msg_Out()<<"Checking GoSam ("<<born<<") vs. Internal ("<<checkBorn<<") failed.\n";
      if(m_errorOnFailedCheck) THROW(fatal_error,"Failed check is not acceptable.");
    }
    else if(m_printAllChecks) {
      msg_Out()<<"Check GoSam ("<<born<<") vs. Internal ("<<checkBorn<<") succesful.\n";
    }
  }
#endif
  return born;
}

Flavour GoSam_Born::CombineQCD(const Flavour &a, const Flavour &b)
{
  if (a.IsGluon() && b.IsGluon()) return Flavour(kf_gluon);
  if (a.IsQuark() && b.IsQuark() && a.Bar()==b) return Flavour(kf_gluon);
  if (a.IsGluon() && b.IsQuark()) return b;
  if (b.IsGluon() && a.IsQuark()) return a;
  return Flavour(kf_none);
}

void GoSam_Born::FillCombinations
(set<pair<size_t,size_t> > &combs,map<size_t,Flavour_Vector> &fls)
{
  DEBUG_FUNC(m_flavs);
  // legs in the all-outgoing convention used by the cluster amplitudes
  vector<pair<size_t,Flavour> > legs;
  for (size_t i(0);i<m_flavs.size();++i)
    legs.push_back(make_pair(size_t(1)<<i,i<m_nin?m_flavs[i].Bar():m_flavs[i]));
  set<vector<size_t> > visited;
  FillCombinations(legs,combs,fls,visited);
  CompleteCombinations(combs,fls);
  msg_Debugging()<<combs.size()<<" combinations, "
                 <<fls.size()<<" combined flavours\n";
}

void GoSam_Born::FillCombinations
(const vector<pair<size_t,Flavour> > &legs,set<pair<size_t,size_t> > &combs,
 map<size_t,Flavour_Vector> &fls,set<vector<size_t> > &visited) const
{
  vector<size_t> key;
  for (const auto &l: legs) key.push_back(l.first);
  sort(key.begin(),key.end());
  // insert into std::set checks if key already existed
  if (!visited.insert(key).second) return;
  const size_t inmask((size_t(1)<<m_nin)-1);
  for (size_t i(0);i<legs.size();++i) {
    for (size_t j(i+1);j<legs.size();++j) {
      const size_t ida(legs[i].first), idb(legs[j].first);
      if ((ida&inmask) && (idb&inmask)) continue; // II doesnt cluster together
      const Flavour &fa(legs[i].second), &fb(legs[j].second);
      if (!fa.Strong() || !fb.Strong()) continue; // no EW cluserings
      const Flavour fab(CombineQCD(fa,fb));
      if (fab.Kfcode()==kf_none) continue;
      combs.insert(make_pair(ida,idb));
      Flavour_Vector &cf(fls[ida|idb]);
      if (find(cf.begin(),cf.end(),fab)==cf.end()) cf.push_back(fab);
      msg_Debugging()<<"comb "<<ID(ida)<<" & "<<ID(idb)<<" -> "<<fab<<"\n";
      // recursion, combinations for list with combined flavs
      vector<pair<size_t,Flavour> > next;
      for (size_t k(0);k<legs.size();++k)
        if (k!=i && k!=j) next.push_back(legs[k]);
      next.push_back(make_pair(ida|idb,fab));
      FillCombinations(next,combs,fls,visited);
    }
  }
}

DECLARE_TREEME2_GETTER(GoSam_Born,
		       "GoSam_Born")

Tree_ME2_Base* ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			      PHASIC::External_ME_Args,
			      GoSam::GoSam_Born>::
operator()(const External_ME_Args &args) const
{
  if (!args.m_source.empty() && args.m_source != "GoSam") return nullptr;
  GoSam_Born * ret = new GoSam_Born(args, true, GoSam_Interface::RegisterProcess(args,true));
#ifdef VALIDATE_GOSAM_INTERFACE
  External_ME_Args argscheck(args.m_inflavs,args.m_outflavs,args.m_orders,"Internal");
  ret->p_check = Tree_ME2_Base::GetME2(argscheck);
  if(ret->p_check)
    msg_Out()<<"Initialized Internal ME to validate GoSam interface.\n";
  else
    msg_Out()<<"Failed to initialize Internal ME to validate GoSam interface.\n";
#endif
  if(ret) {
    msg_Debugging()<<"Return GoSam Born.\n";
  }
  return ret;
}

int GoSam_Born::OrderQCD(const int &id) const
{ return m_oqcd; }

int GoSam_Born::OrderEW(const int &id) const
{ return m_oqed; }
