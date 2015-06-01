#include "PHASIC++/Scales/Core_Scale_Setter.H"

#include "PHASIC++/Process/Single_Process.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

namespace PHASIC {

  class Default_Core_Scale: public Core_Scale_Setter {
  public:

    Default_Core_Scale(const Core_Scale_Arguments &args):
      Core_Scale_Setter(args) {}

    PDF::CParam Calculate(ATOOLS::Cluster_Amplitude *const ampl);

    ATOOLS::Cluster_Amplitude *Cluster
    (ATOOLS::Cluster_Amplitude *const ampl) const;

  };// end of class Scale_Setter_Base

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

PDF::CParam Default_Core_Scale::Calculate(Cluster_Amplitude *const ampl)
{
  DEBUG_FUNC("");
  msg_Debugging()<<*ampl<<"\n";
  if (ampl->Legs().size()==3 && ampl->NIn()==2) {
    double kt2cmin(ampl->Leg(2)->Mom().Abs2());
    return PDF::CParam(kt2cmin,kt2cmin,0.0,kt2cmin,-1);
  }
  double muf2(0.0), mur2(0.0), muq2(0.0);
  Cluster_Amplitude *campl(Cluster(ampl->Copy()));
  if (campl->Legs().size()!=ampl->Legs().size())
    msg_Debugging()<<*campl<<"\n";
  if (campl->Legs().size()!=4) {
    double q=0.0;
    Vec4D ewsum;
    for (size_t i(0);i<campl->Legs().size();++i)
      if (!campl->Leg(i)->Flav().Strong()) ewsum+=campl->Leg(i)->Mom();
      else q+=campl->Leg(i)->Mom().MPerp();
    q+=ewsum.MPerp();
    campl->Delete();
    return PDF::CParam(q*q/4.0,q*q/4.0,0.0,q*q/4.0,-1);
  }
  Flavour_Vector fl; fl.resize(4);
  fl[0]=campl->Leg(0)->Flav();
  fl[1]=campl->Leg(1)->Flav();
  fl[2]=campl->Leg(2)->Flav();
  fl[3]=campl->Leg(3)->Flav();
  if (fl[0].Strong() && fl[1].Strong()) {// hh collision
    if (fl[2].Strong() && fl[3].Strong()) {
      msg_Debugging()<<"pure QCD like\n";
      double s(2.0*campl->Leg(0)->Mom()*campl->Leg(1)->Mom());
      double t(2.0*campl->Leg(0)->Mom()*campl->Leg(2)->Mom());
      double u(2.0*campl->Leg(0)->Mom()*campl->Leg(3)->Mom());
      muq2=muf2=mur2=-1.0/(1.0/s+1.0/t+1.0/u)/4.0;
    }
    else if (!fl[2].Strong() && !fl[3].Strong()) {
      msg_Debugging()<<"DY like\n";
      muq2=muf2=mur2=(campl->Leg(0)->Mom()+campl->Leg(1)->Mom()).Abs2();
    }
    else if (fl[2].Strong() && !fl[3].Strong()) {
      msg_Debugging()<<"jV like\n";
      muq2=muf2=mur2=campl->Leg(3)->Mom().MPerp2()/4.0;
    }
    else if (!fl[2].Strong() && fl[3].Strong()) {
      msg_Debugging()<<"Vj like\n";
      muq2=muf2=mur2=campl->Leg(2)->Mom().MPerp2()/4.0;
    }
    else THROW(fatal_error,"Internal error.");
  }
  else if (!fl[0].Strong() && !fl[1].Strong()) {// ll collision
    if (fl[2].Strong() && fl[3].Strong()) {
      msg_Debugging()<<"jets like\n";
      muq2=muf2=mur2=(campl->Leg(0)->Mom()+campl->Leg(1)->Mom()).Abs2();
    }
  }
  else {
    if (!fl[0].Strong() && !fl[2].Strong()) {
      msg_Debugging()<<"DIS like\n";
      muq2=muf2=mur2=dabs((campl->Leg(fl[0].Strong()?1:0)->Mom()+
			   campl->Leg(fl[2].Strong()?3:2)->Mom()).Abs2());
    }
    else {
      msg_Debugging()<<"QCD Compton like\n";
      muq2=muf2=mur2=dabs(sqrt(campl->Leg(2)->Mom().MPerp2()*
			       campl->Leg(3)->Mom().MPerp2()));
    }
  }
  campl->Delete();
  msg_Debugging()<<"\\mu_f = "<<sqrt(muf2)<<"\n"
		 <<"\\mu_r = "<<sqrt(mur2)<<"\n"
		 <<"\\mu_q = "<<sqrt(muq2)<<"\n";
  return PDF::CParam(muf2,muq2,0.0,mur2,-1);
}

Cluster_Amplitude *Default_Core_Scale::Cluster
(Cluster_Amplitude *const ampl) const
{
  if (ampl->Legs().size()==ampl->NIn()+2) return ampl;
  Single_Process *proc(ampl->Proc<Single_Process>());
  for (size_t i(ampl->NIn());i<ampl->Legs().size();++i) {
    Cluster_Leg *li(ampl->Leg(i));
    for (size_t j(i+1);j<ampl->Legs().size();++j) {
      Cluster_Leg *lj(ampl->Leg(j));
      if (proc->Combinable(li->Id(),lj->Id())) {
	Flavour_Vector fls(proc->CombinedFlavour(li->Id()|lj->Id()));
	for (size_t k(0);k<fls.size();++k) {
	  bool dec(false);
	  for (size_t l(0);l<ampl->Decays().size();++l)
	    if (ampl->Decays()[l]->m_id==(li->Id()|lj->Id())) {
	      dec=true;
	      break;
	    }
	  if ((!li->Flav().Strong() && !lj->Flav().Strong() &&
	       !fls[k].Strong()) || dec) {
	    msg_Debugging()<<"combine "<<ID(li->Id())<<"&"<<ID(lj->Id())
			   <<"->"<<fls[k]<<" ("<<dec<<")\n";
	    li->SetFlav(fls[k]);
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
      }
    }
  }
  return ampl;
}

DECLARE_ND_GETTER(Default_Core_Scale,"DEFAULT",
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
