#include "PHASIC++/Scales/Core_Scale_Setter.H"

#include "PHASIC++/Process/Single_Process.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

namespace PHASIC {

  class E_Gamma_Core_Scale: public Core_Scale_Setter {
  public:

    E_Gamma_Core_Scale(const Core_Scale_Arguments &args):
      Core_Scale_Setter(args) {}

    PDF::Cluster_Param Calculate(ATOOLS::Cluster_Amplitude *const ampl);

    ATOOLS::Cluster_Amplitude *Cluster
    (ATOOLS::Cluster_Amplitude *const ampl) const;

  };// end of class Scale_Setter_Base

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

PDF::Cluster_Param E_Gamma_Core_Scale::Calculate(Cluster_Amplitude *const ampl)
{
  DEBUG_FUNC("");
  msg_Debugging()<<*ampl<<"\n";
  if (ampl->Legs().size()==3 && ampl->NIn()==2) {
    double kt2cmin(ampl->Leg(2)->Mom().Abs2());
    return PDF::Cluster_Param(NULL,kt2cmin,kt2cmin,kt2cmin,-1);
  }
  Cluster_Amplitude *campl(Cluster(ampl->Copy()));
  if (campl->Legs().size()!=ampl->Legs().size())
    msg_Debugging()<<*campl<<"\n";
  // Q^2 is the momentum transfer on the lepton line, i.e. the invariant mass
  // of the sum of all (incoming and outgoing) lepton momenta. Identifying the
  // leptons by flavour rather than by a fixed leg index is robust against the
  // leg ordering changing during clustering and works for both the direct and
  // the resolved channel. H_{T,had} is the scalar sum of the final-state
  // hadronic transverse masses.
  Vec4D plep(0.,0.,0.,0.);
  bool haslep(false);
  double ht(0.0);
  for (size_t i(0);i<campl->Legs().size();++i) {
    if (campl->Leg(i)->Flav().IsLepton()) {
      plep+=campl->Leg(i)->Mom();
      haslep=true;
    }
    else if (i>=campl->NIn())
      ht+=sqrt(dabs(campl->Leg(i)->Mom().MPerp2()));
  }
  const double Q2(haslep?dabs(plep.Abs2()):0.0);
  // 2->2 DIS-like core: the scale is the pure momentum transfer Q^2.
  // For higher multiplicities the hadronic H_{T,had}^2 is added.
  const double mu2((campl->Legs().size()==4 && haslep) ?
		   Q2 : (Q2+sqr(ht))/4.0);
  campl->Delete();
  msg_Debugging()<<"\\mu_f = \\mu_r = \\mu_q = "<<sqrt(mu2)<<"\n";
  return PDF::Cluster_Param(NULL,mu2,mu2,mu2,-1);
}

Cluster_Amplitude *E_Gamma_Core_Scale::Cluster
(Cluster_Amplitude *const ampl) const
{
  struct Combination { size_t i, j; Flavour fl;
    inline Combination(const size_t &_i=0,const size_t &_j=0,
		       const Flavour &_fl=kf_none):
      i(_i), j(_j), fl(_fl) {} };// end of struct
  if (ampl->Legs().size()==ampl->NIn()+2) return ampl;
  Single_Process *proc(ampl->Proc<Single_Process>());
  std::map<double,Combination,std::less<double> > tij;
  for (size_t i(ampl->NIn());i<ampl->Legs().size();++i) {
    Cluster_Leg *li(ampl->Leg(i));
    for (size_t j(i+1);j<ampl->Legs().size();++j) {
      Cluster_Leg *lj(ampl->Leg(j));
      if (proc->Combinable(li->Id(),lj->Id())) {
	Flavour_Vector fls(proc->CombinedFlavour(li->Id()|lj->Id()));
	for (size_t k(0);k<fls.size();++k) {
	  double t((li->Mom()+lj->Mom()).Abs2());
	  double p(sqr(t-sqr(fls[k].Mass()))+
		   sqr(fls[k].Mass()*fls[k].Width()));
	  msg_Debugging()<<"check "<<ID(li->Id())<<"&"<<ID(lj->Id())
			 <<"["<<fls[k]<<"] -> m = "<<sqrt(dabs(t))
			 <<", 1/p = "<<sqrt(p)<<"\n";
	  tij[p]=Combination(i,j,fls[k]);
	}
      }
    }
  }
  for (std::map<double,Combination,std::less<double> >::
	 const_iterator it(tij.begin());it!=tij.end();++it) {
    Cluster_Leg *li(ampl->Leg(it->second.i));
    Cluster_Leg *lj(ampl->Leg(it->second.j));
    bool dec(false);
    for (size_t l(0);l<ampl->Decays().size();++l)
      if (ampl->Decays()[l]->m_id==(li->Id()|lj->Id())) {
	dec=true;
	break;
      }
    if ((!li->Flav().Strong() && !lj->Flav().Strong() &&
	 !it->second.fl.Strong()) || dec) {
      msg_Debugging()<<"combine "<<ID(li->Id())<<"&"<<ID(lj->Id())
		     <<"->"<<it->second.fl<<" ("<<dec<<")\n";
      li->SetFlav(it->second.fl);
      li->SetMom(li->Mom()+lj->Mom());
      li->SetId(li->Id()|lj->Id());
      lj->Delete();
      for (ClusterLeg_Vector::iterator lit(ampl->Legs().begin());
	   lit!=ampl->Legs().end();++lit)
	if (*lit==lj) {
	  ampl->Legs().erase(lit);
	  break;
	}
      return Cluster(ampl);
    }
  }
  return ampl;
}

DECLARE_ND_GETTER(E_Gamma_Core_Scale,"E_Gamma",
		  Core_Scale_Setter,Core_Scale_Arguments,true);

Core_Scale_Setter *ATOOLS::Getter
<Core_Scale_Setter,Core_Scale_Arguments,E_Gamma_Core_Scale>::
operator()(const Core_Scale_Arguments &args) const
{
  return new E_Gamma_Core_Scale(args);
}

void ATOOLS::Getter<Core_Scale_Setter,Core_Scale_Arguments,
		    E_Gamma_Core_Scale>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"E_Gamma core scale";
}
