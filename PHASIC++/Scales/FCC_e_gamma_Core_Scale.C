#include "PHASIC++/Scales/Core_Scale_Setter.H"

#include "PHASIC++/Process/Single_Process.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

namespace PHASIC {

  class FCC_e_gamma_Core_Scale: public Core_Scale_Setter {
  public:

    FCC_e_gamma_Core_Scale(const Core_Scale_Arguments &args):
      Core_Scale_Setter(args) {}

    PDF::Cluster_Param Calculate(ATOOLS::Cluster_Amplitude *const ampl);

    ATOOLS::Cluster_Amplitude *Cluster
    (ATOOLS::Cluster_Amplitude *const ampl) const;

  };// end of class Scale_Setter_Base

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

PDF::Cluster_Param FCC_e_gamma_Core_Scale::Calculate(Cluster_Amplitude *const ampl)
{
  DEBUG_FUNC("");
  msg_Debugging()<<*ampl<<"\n";
  if (ampl->Legs().size()==3 && ampl->NIn()==2) {
    double kt2cmin(ampl->Leg(2)->Mom().Abs2());
    return PDF::Cluster_Param(NULL,kt2cmin,kt2cmin,kt2cmin,-1);
  }
  double muf2(0.0), mur2(0.0), muq2(0.0);
  Cluster_Amplitude *campl(Cluster(ampl->Copy()));
  if (campl->Legs().size()!=ampl->Legs().size())
    msg_Debugging()<<*campl<<"\n";
  if (campl->Legs().size()!=4) {
    msg_Debugging()<<"more than 4 legs, use sqrt(Q^2+H_T'^2)/2 as scale"<<std::endl;
    double q=0.0;
    for (size_t i(0);i<campl->Legs().size();++i)
      if (!campl->Leg(i)->Flav().IsLepton()) q+=sqrt(dabs(campl->Leg(i)->Mom().MPerp2()));
    q = q*q; //< this now is H_{T,hadr}^2
    if (campl->Leg(0)->Flav().IsLepton() && campl->Leg(2)->Flav().IsLepton())
      q += dabs((campl->Leg(0)->Mom() + campl->Leg(2)->Mom()).Abs2());
    campl->Delete();
    return PDF::Cluster_Param(NULL,q/4.0,q/4.0,q/4.0,-1);
  }
  Flavour_Vector fl; fl.resize(4);
  fl[0]=campl->Leg(0)->Flav();
  fl[1]=campl->Leg(1)->Flav();
  fl[2]=campl->Leg(2)->Flav();
  fl[3]=campl->Leg(3)->Flav();
  if (fl[0].IsLepton() && fl[2].IsLepton()) { // DIS-like
    muq2=muf2=mur2=dabs((campl->Leg(0)->Mom()+campl->Leg(2)->Mom()).Abs2());
  } else {
    muq2=muf2=mur2=0.;
  }
  campl->Delete();
  msg_Debugging()<<"\\mu_f = "<<sqrt(muf2)<<"\n"
		 <<"\\mu_r = "<<sqrt(mur2)<<"\n"
		 <<"\\mu_q = "<<sqrt(muq2)<<"\n";
  return PDF::Cluster_Param(NULL,muq2,muf2,mur2,-1);
}

Cluster_Amplitude *FCC_e_gamma_Core_Scale::Cluster
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

DECLARE_ND_GETTER(FCC_e_gamma_Core_Scale,"FCC_e_gamma",
		  Core_Scale_Setter,Core_Scale_Arguments,true);

Core_Scale_Setter *ATOOLS::Getter
<Core_Scale_Setter,Core_Scale_Arguments,FCC_e_gamma_Core_Scale>::
operator()(const Core_Scale_Arguments &args) const
{
  return new FCC_e_gamma_Core_Scale(args);
}

void ATOOLS::Getter<Core_Scale_Setter,Core_Scale_Arguments,
		    FCC_e_gamma_Core_Scale>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"FCC_e_gamma core scale";
}
