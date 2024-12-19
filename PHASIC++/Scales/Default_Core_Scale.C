#include "PHASIC++/Scales/Core_Scale_Setter.H"

#include "PHASIC++/Process/Single_Process.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

namespace PHASIC {

  class Default_Core_Scale: public Core_Scale_Setter {
  public:

    Default_Core_Scale(const Core_Scale_Arguments &args):
      Core_Scale_Setter(args) {}

    PDF::Cluster_Param Calculate(ATOOLS::Cluster_Amplitude *const ampl);

    ATOOLS::Cluster_Amplitude *Cluster
    (ATOOLS::Cluster_Amplitude *const ampl) const;

  };// end of class Scale_Setter_Base

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

PDF::Cluster_Param Default_Core_Scale::Calculate(Cluster_Amplitude *const ampl)
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
    msg_Debugging()<<"more than 4 legs, use HT'/2 as scale"<<std::endl;
    double q=0.0;
    Vec4D ewsum;
    for (size_t i(0);i<campl->Legs().size();++i)
      if (!campl->Leg(i)->Flav().Strong()) ewsum+=campl->Leg(i)->Mom();
      else q+=sqrt(dabs(campl->Leg(i)->Mom().MPerp2()));
    q+=sqrt(dabs(ewsum.MPerp2()));
    campl->Delete();
    return PDF::Cluster_Param(NULL,q*q/4.0,q*q/4.0,q*q/4.0,-1);
  }
  Flavour_Vector fl; fl.resize(4);
  fl[0]=campl->Leg(0)->Flav();
  fl[1]=campl->Leg(1)->Flav();
  fl[2]=campl->Leg(2)->Flav();
  fl[3]=campl->Leg(3)->Flav();
  if (fl[0].Strong() && fl[1].Strong()) {// hh collision
    if (fl[2].Strong() && fl[3].Strong()) {
      msg_Debugging()<<"pure QCD like\n";
      muq2=muf2=mur2=(campl->Leg(2)->Mom().PPerp2()+
		      campl->Leg(3)->Mom().PPerp2())/8.;
    }
    else if (!fl[2].Strong() && !fl[3].Strong()) {
      msg_Debugging()<<"DY like\n";
      muq2=muf2=mur2=(campl->Leg(2)->Mom()+
		      campl->Leg(3)->Mom()).Abs2();
    }
    else if (fl[2].Strong() && !fl[3].Strong()) {
      msg_Debugging()<<"jV like\n";
      muq2=muf2=Max(campl->Leg(3)->Mom().Abs2(),
		    campl->Leg(2)->Mom().PPerp2());
      mur2=campl->Leg(2)->Mom().PPerp2()/4.;
      if (fl[3].Kfcode()==25) {
	msg_Debugging()<<"H special\n";
	mur2=pow(mur2*pow(fl[3].Mass(),4.),1./3.); 
      }
    }
    else if (!fl[2].Strong() && fl[3].Strong()) {
      msg_Debugging()<<"Vj like\n";
      muq2=muf2=Max(campl->Leg(2)->Mom().Abs2(),
		    campl->Leg(3)->Mom().PPerp2());
      mur2=campl->Leg(3)->Mom().PPerp2()/4.;
      if (fl[2].Kfcode()==25) {
	msg_Debugging()<<"H special\n";
	mur2=pow(mur2*pow(fl[2].Mass(),4.),1./3.); 
      }
    }
    else THROW(fatal_error,"Internal error.");
  }
  else if (!fl[0].Strong() && !fl[1].Strong()) {// ll collision
    if (fl[2].Strong() && fl[3].Strong()) {
      msg_Debugging()<<"ll->jets like\n";
    } else {
      msg_Debugging()<<"ll->unknown, Mandelstam s will be used as the scale\n";
    }
    muq2=muf2=mur2=(campl->Leg(0)->Mom()+campl->Leg(1)->Mom()).Abs2();
  }
  else {
    if (!fl[0].Strong() && !fl[2].Strong()) {
      msg_Debugging()<<"DIS like\n";
      muq2=muf2=mur2=dabs((campl->Leg(0)->Mom()+campl->Leg(2)->Mom()).Abs2());
    } else {
      msg_Debugging()<<"QCD Compton like, i.e. q+gamma -> q+gluon\n";
      muq2=muf2=mur2=dabs(sqrt(campl->Leg(2)->Mom().MPerp2()*
			       campl->Leg(3)->Mom().MPerp2()));
    }
  }
  campl->Delete();
  msg_Debugging()<<"\\mu_f = "<<sqrt(muf2)<<"\n"
		 <<"\\mu_r = "<<sqrt(mur2)<<"\n"
		 <<"\\mu_q = "<<sqrt(muq2)<<"\n";
  return PDF::Cluster_Param(NULL,muq2,muf2,mur2,-1);
}

Cluster_Amplitude *Default_Core_Scale::Cluster
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

DECLARE_ND_GETTER(Default_Core_Scale,"Default",
		  Core_Scale_Setter,Core_Scale_Arguments,true);

Core_Scale_Setter *ATOOLS::Getter
<Core_Scale_Setter,Core_Scale_Arguments,Default_Core_Scale>::
operator()(const Core_Scale_Arguments &args) const
{
  return new Default_Core_Scale(args);
}

void ATOOLS::Getter<Core_Scale_Setter,Core_Scale_Arguments,
		    Default_Core_Scale>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"default core scale"; 
}
