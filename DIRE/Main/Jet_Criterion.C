#include "PDF/Main/Jet_Criterion.H"
#include "PDF/Main/Shower_Base.H"
#include "PDF/Main/Cluster_Definitions_Base.H"
#include "PHASIC++/Selectors/Jet_Finder.H"
#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Org/Message.H"

using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;

namespace DIRE {

  class Jet_Criterion: public PDF::Jet_Criterion {
  private:

    struct Q2_Value {
      double m_q2;
      ATOOLS::Flavour m_fl;
      int m_i, m_j, m_k;
      inline Q2_Value(const double &q2,
		      const ATOOLS::Flavour &fl,
		      int i, int j, int k):
	m_q2(q2), m_fl(fl), m_i(i), m_j(j), m_k(k) {}
      inline bool operator<(const Q2_Value &v) const
      { return m_q2<v.m_q2; }
    };

    PDF::Cluster_Definitions_Base *p_clus;

  public:

    Jet_Criterion(const JetCriterion_Key &args):
      p_clus(args.p_shower->GetClusterDefinitions()) {}

    double Qij2(const Vec4D &pi,const Vec4D &pj,const Vec4D &pk,
		const Flavour &fi,const Flavour &fj) const
    {
      // arXiv:2002.11114 [hep-ph]
      const double beta(0.5);
      double t1(2.0*(pi*pj)*(pj*pk)/(pi*pk));
      double t2(2.0*(pj*pi)*(pi*pk)/(pj*pk));
      double xi1(dabs((pi*pj)/(pk*pj)));
      double xi2(dabs((pj*pi)/(pk*pi)));
      t1*=pow(Max(xi1,1.0/xi1),-beta/2.0);
      t2*=pow(Max(xi2,1.0/xi2),-beta/2.0);
      if (pi[0]<0.0) return dabs(t1);
      if (pj[0]<0.0) return dabs(t2);
      return Min(t1,t2);
    }

    double Value(Cluster_Amplitude *ampl,int mode)
    {
      DEBUG_FUNC("mode = "<<mode);
      msg_Debugging()<<*ampl<<"\n";
      NLO_subevtlist *subs(NULL);
      if (mode) subs=ampl->Proc<PHASIC::Process_Base>()->GetRSSubevtList();
      size_t noem(0), nospec(0);
      for (size_t i(0);i<ampl->Decays().size();++i) {
	noem|=ampl->Decays()[i]->m_id;
	if (!ampl->Decays()[i]->m_fl.Strong())
	  nospec|=ampl->Decays()[i]->m_id;
      }
      msg_Debugging()<<"noem = "<<ID(noem)<<", nospec = "<<ID(nospec)<<"\n";
      std::set<Q2_Value> q2list;
      for (size_t i(0);i<ampl->Legs().size();++i) {
	Cluster_Leg *li(ampl->Leg(i));
	if (li->Id()&noem) continue;
	Flavour fi(i<ampl->NIn()?li->Flav().Bar():li->Flav());
	for (size_t j(Max(i+1,ampl->NIn()));j<ampl->Legs().size();++j) {
	  Cluster_Leg *lj(ampl->Leg(j));
	  if (lj->Id()&noem) continue;
	  Flavour fj(j<ampl->NIn()?lj->Flav().Bar():lj->Flav());
	  for (size_t k(0);k<ampl->Legs().size();++k) {
	    if (k==i || k==j) continue;
	    if (subs) {
	      bool found(false);
	      for (size_t l(0);l<subs->size()-1;++l) {
		NLO_subevt *sub((*subs)[l]);
		if (k==sub->m_k && ((i==sub->m_i && j==sub->m_j) ||
				    (i==sub->m_j && j==sub->m_i))) {
		  found=true;
		  break;
		}
	      }
	      if (!found) continue;
	    }
	    Cluster_Leg *lk(ampl->Leg(k));
	    if (lk->Id()&nospec) continue;
	    Flavour fk(k<ampl->NIn()?lk->Flav().Bar():lk->Flav());
	    if (lk->Flav().Strong() &&
		li->Flav().Strong() && lj->Flav().Strong() &&
		(li->Flav().IsGluon() || lj->Flav().IsGluon() ||
		 li->Flav()==lj->Flav().Bar())) {
	      double q2ijk(Qij2(li->Mom(),lj->Mom(),lk->Mom(),
				li->Flav(),lj->Flav()));
	      msg_Debugging()<<"Q_{"<<ID(li->Id())<<ID(lj->Id())
			     <<","<<ID(lk->Id())<<"} = "<<sqrt(q2ijk)<<"\n";
	      if (q2ijk<0.0) continue;
	      Flavour mofl=Flavour(kf_gluon);
	      if (li->Flav().IsGluon()) mofl=lj->Flav();
	      if (lj->Flav().IsGluon()) mofl=li->Flav();
	      q2list.insert(Q2_Value(q2ijk,mofl,i,j,k));
	    }
	    else {
	      msg_IODebugging()<<"No kernel for "<<fi<<" "<<fj<<" <-> "<<fk<<"\n";
	    }
	  }
	}
      }
      if (mode==0) {
	double q2min(std::numeric_limits<double>::max());
	if (q2list.size()) q2min=q2list.begin()->m_q2;
	msg_Debugging()<<"--- "<<sqrt(q2min)<<" ---\n";
	return q2min;
      }
      while (q2list.size()) { 
	Flavour mofl(q2list.begin()->m_fl);
	size_t imin(q2list.begin()->m_i);
	size_t jmin(q2list.begin()->m_j);
	size_t kmin(q2list.begin()->m_k);
	q2list.erase(q2list.begin());
	Cluster_Param cp=p_clus->Cluster
	  (Cluster_Config(ampl,imin,jmin,kmin,mofl,ampl->MS(),NULL,1));
	if (cp.m_pijt==Vec4D())
	  cp=p_clus->Cluster(Cluster_Config(ampl,jmin,imin,kmin,mofl,ampl->MS(),NULL,1));
	if (cp.m_pijt==Vec4D()) continue;
	Cluster_Amplitude *bampl(Cluster_Amplitude::New());
	bampl->SetProc(ampl->Proc<void>());
	bampl->SetNIn(ampl->NIn());
	bampl->SetJF(ampl->JF<void>());
	for (int i(0), j(0);i<ampl->Legs().size();++i) {
	  if (i==jmin) continue;
	  if (i==imin) {
	    bampl->CreateLeg(cp.m_pijt,mofl,ampl->Leg(i)->Col());
	    bampl->Legs().back()->SetId(ampl->Leg(imin)->Id()|ampl->Leg(jmin)->Id());
	    bampl->Legs().back()->SetK(ampl->Leg(kmin)->Id());	
	  }
	  else {
	    bampl->CreateLeg(i==kmin?cp.m_pkt:cp.m_lam*ampl->Leg(i)->Mom(),
			     ampl->Leg(i)->Flav(),ampl->Leg(i)->Col());
	  }
	  ++j;
	}
	double res=Value(bampl,0);
	bampl->Delete();
	if (res==std::numeric_limits<double>::max()) continue;
	return res;
      }
      msg_Debugging()<<METHOD<<"(): Combine failed. Use R configuration."<<std::endl;
      return Value(ampl,0);
    }

  };// end of class DIRE_Jet_Criterion

}// end of namespace DIRE

DECLARE_GETTER(DIRE::Jet_Criterion,"Dire",Jet_Criterion,JetCriterion_Key);

Jet_Criterion *Getter<Jet_Criterion,JetCriterion_Key,DIRE::Jet_Criterion>::
operator()(const JetCriterion_Key &args) const
{
  return new DIRE::Jet_Criterion(args);
}

void Getter<Jet_Criterion,JetCriterion_Key,DIRE::Jet_Criterion>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"The DIRE jet criterion";
}
