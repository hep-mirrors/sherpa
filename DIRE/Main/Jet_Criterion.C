#include "PDF/Main/Jet_Criterion.H"
#include "PDF/Main/Shower_Base.H"
#include "PDF/Main/Cluster_Definitions_Base.H"
#include "PHASIC++/Selectors/Jet_Finder.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Org/Message.H"

using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;

namespace DIRE {

  class Jet_Criterion: public PDF::Jet_Criterion {
  private:

    PDF::Cluster_Definitions_Base *p_clus;

  public:

    Jet_Criterion(const JetCriterion_Key &args):
      p_clus(args.p_shower->GetClusterDefinitions()) {}

    double Qij2(const Vec4D &pi,const Vec4D &pj,const Vec4D &pk,
		const Flavour &fi,const Flavour &fj) const
    {
      Vec4D npi(pi), npj(pj);
      if (npi[0]<0.0) npi=-pi-pj;
      if (npj[0]<0.0) npj=-pj-pi;
      double pipj(dabs(npi*npj)), pipk(dabs(npi*pk)), pjpk(dabs(npj*pk));
      double Cij(pipk/(pipj+pjpk)), Cji(pjpk/(pipj+pipk));
      return 2.0*dabs(pi*pj)/(Cij+Cji);
    }

    bool Jets(Cluster_Amplitude *ampl,int mode)
    {
      DEBUG_FUNC("mode = "<<mode);
      msg_Debugging()<<*ampl<<"\n";
      Jet_Finder *jf(ampl->JF<Jet_Finder>());
      double q2cut(jf->Ycut()*sqr(rpa->gen.Ecms()));
      size_t noem(0), nospec(0);
      for (size_t i(0);i<ampl->Decays().size();++i) {
	noem|=ampl->Decays()[i]->m_id;
	if (!ampl->Decays()[i]->m_fl.Strong())
	  nospec|=ampl->Decays()[i]->m_id;
      }
      msg_Debugging()<<"noem = "<<ID(noem)<<", nospec = "<<ID(nospec)<<"\n";
      double q2min(std::numeric_limits<double>::max());
      size_t imin(0), jmin(0), kmin(0);
      Flavour mofl;
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
	    Cluster_Leg *lk(ampl->Leg(k));
	    if (lk->Id()&nospec) continue;
	    Flavour fk(k<ampl->NIn()?lk->Flav().Bar():lk->Flav());
	    if (lk->Flav().Strong() &&
		li->Flav().Strong() && lj->Flav().Strong()) {
	      double q2ijk(Qij2(li->Mom(),lj->Mom(),lk->Mom(),
				li->Flav(),lj->Flav()));
	      msg_Debugging()<<"Q_{"<<ID(li->Id())<<ID(lj->Id())
			     <<","<<ID(lk->Id())<<"} = "<<sqrt(q2ijk)<<"\n";
	      if (q2ijk<0.0) continue;
	      if (mode==0) {
		if (q2ijk<q2cut) return false;
	      }
	      else {
		if (q2ijk<q2min) {
		  q2min=q2ijk;
		  mofl=Flavour(kf_gluon);
		  if (li->Flav().IsGluon()) mofl=lj->Flav();
		  if (lj->Flav().IsGluon()) mofl=li->Flav();
		  imin=i;
		  jmin=j;
		  kmin=k;
		}
	      }
	    }
	    else {
	      msg_Debugging()<<"No kernel for "<<fi<<" "<<fj<<" <-> "<<fk<<"\n";
	    }
	  }
	}
      }
      if (mode!=0 && imin!=jmin) {
	Vec4D_Vector p=p_clus->Combine(*ampl,imin,jmin,kmin,mofl,ampl->MS(),1);
	if (p.empty()) {
	  msg_Error()<<METHOD<<"(): Combine failed. Use R configuration."<<std::endl;
	  return Jets(ampl,0);
	}
	Cluster_Amplitude *bampl(Cluster_Amplitude::New());
	bampl->SetProc(ampl->Proc<void>());
	bampl->SetNIn(ampl->NIn());
	bampl->SetJF(ampl->JF<void>());
	for (int i(0), j(0);i<ampl->Legs().size();++i) {
	  if (i==jmin) continue;
	  if (i==imin) {
	    bampl->CreateLeg(p[j],mofl,ampl->Leg(i)->Col());
	    bampl->Legs().back()->SetId(ampl->Leg(imin)->Id()|ampl->Leg(jmin)->Id());
	    bampl->Legs().back()->SetK(ampl->Leg(kmin)->Id());	
	  }
	  else {
	    bampl->CreateLeg(p[j],ampl->Leg(i)->Flav(),ampl->Leg(i)->Col());
	  }
	  ++j;
	}
	bool res=Jets(bampl,0);
	bampl->Delete();
	return res;
      }
      msg_Debugging()<<"--- Jet veto ---\n";
      return true;
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
