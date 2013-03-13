#include "PHASIC++/Scales/Core_Scale_Setter.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Process/Process_Base.H"

#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

namespace PHASIC {

  class WBF_Core_Scale: public PHASIC::Core_Scale_Setter {
  public:

    WBF_Core_Scale(const PHASIC::Core_Scale_Arguments &args):
      Core_Scale_Setter(args) {}

    PDF::CParam Calculate(ATOOLS::Cluster_Amplitude *const ampl);

  };// end of class WBF_Core_Scale

}// end of namespace PHASIC

using namespace PHASIC;
using namespace PHASIC;
using namespace ATOOLS;

PDF::CParam WBF_Core_Scale::Calculate(Cluster_Amplitude *const ampl)
{
  double muf2(-1.), mur2(-1.), q2(-1.);
  size_t nlegs(ampl->Legs().size());
  if (nlegs<5) {
    msg_Error()<<"ERROR in "<<METHOD<<": "
	       <<"Amplitude with "<<nlegs<<" legs only.\n"
	       <<"   Will return all scales = -1 and hope for the best.\n"; 
    return PDF::CParam(muf2,q2,0.0,mur2,-1);
  }
  //msg_Out()<<"==========================================================\n"
  //	   <<METHOD<<"("<<nlegs<<" legs)\n";
  std::list<Cluster_Leg *> stronglegs;
  for (size_t i(0);i<ampl->Legs().size();i++) {
    Cluster_Leg * leg(ampl->Leg(i));
    if (leg->Flav().Strong()) stronglegs.push_back(leg);
  }
  Vec4D maxvec;
  for (std::list<Cluster_Leg *>::iterator leg=stronglegs.begin();
       leg!=stronglegs.end();leg++) {
    if (ID((*leg)->Id())[0]<2) {
      //msg_Out()<<"Add: ";
      maxvec += (*leg)->Mom();
    }
    //else msg_Out()<<"     ";
    //msg_Out()<<(**leg)<<"  ("<<sqrt((*leg)->Mom().Abs2())<<")\n";
  }
  double propscale(sqr(maxvec.Abs2())), t12, t34, test;
  bool wbfkin(false), wbftest;
  std::vector<size_t> links;
  Single_Process *proc(p_proc->Get<Single_Process>());
  for (std::list<Cluster_Leg *>::iterator leg1=stronglegs.begin();
       leg1!=stronglegs.end();leg1++) {
    for (std::list<Cluster_Leg *>::iterator leg2=leg1;
	 leg2!=stronglegs.end();leg2++) {
      if ((*leg2)==(*leg1) ||
	  !proc->Combinable((*leg1)->Id(),(*leg2)->Id())) continue;
      t12 = ((*leg1)->Mom()+(*leg2)->Mom()).Abs2();
      for (std::list<Cluster_Leg *>::iterator leg3=leg1;
	   leg3!=stronglegs.end();leg3++) {
	if ((*leg3)==(*leg1)||(*leg3)==(*leg2)) continue;
	for (std::list<Cluster_Leg *>::iterator leg4=leg3;
	     leg4!=stronglegs.end();leg4++) {
	  if ((*leg4)==(*leg1) || (*leg4)==(*leg2) || (*leg4)==(*leg3)||
	      !proc->Combinable((*leg3)->Id(),(*leg4)->Id())) continue;
	  t34 = ((*leg3)->Mom()+(*leg4)->Mom()).Abs2();
	  test    = 1.;
	  if ((*leg1)->Flav()==(*leg2)->Flav().Bar()) 
	    test *= dabs(t12-sqr(91.2));
	  else 
	    test *= dabs(t12-sqr(80.4));
	  if ((*leg3)->Flav()==(*leg4)->Flav().Bar()) 
	    test *= dabs(t34-sqr(91.2));
	  else 
	    test *= dabs(t34-sqr(80.4));
	  if (t12<0. && t34<0.) {
	    wbftest = true;
	    msg_Debugging()<<"   --> tt ["<<ID((*leg1)->Id())[0]<<", "
			   <<ID((*leg2)->Id())[0]<<"] "
			   <<" and ["<<ID((*leg3)->Id())[0]<<", "
			   <<ID((*leg4)->Id())[0]<<"]: "<<sqrt(test)<<"\n";
	  }
	  else if (t12>0. && t34>0.) {
	    wbftest = false;
	    msg_Debugging()<<"   --> ss ["<<ID((*leg1)->Id())[0]<<", "
			   <<ID((*leg2)->Id())[0]<<"] "
			   <<" and ["<<ID((*leg3)->Id())[0]<<", "
			   <<ID((*leg4)->Id())[0]<<"]: "<<sqrt(test)<<"\n";
	  }
	  else {
	    msg_Out()<<"   --> st ["<<ID((*leg1)->Id())[0]<<", "
		     <<ID((*leg2)->Id())[0]<<"] "
		     <<" and ["<<ID((*leg3)->Id())[0]<<", "
		     <<ID((*leg4)->Id())[0]<<"]: "<<t12<<"/"<<t34<<"\n";
	  }
	  if (test<propscale) {
	    propscale  = dabs(t12*t34);
	    wbfkin = wbftest;
	    links.clear();
	    links.push_back((*leg1)->Id());
	    links.push_back((*leg2)->Id());
	    links.push_back((*leg3)->Id());
	    links.push_back((*leg4)->Id());
	  }
	}
      }
    }
  }

  if (wbfkin) {
    muf2 = mur2 = sqrt(dabs(t12*t34));
    q2   = Max(dabs(t12),dabs(t34));
  }
  else {
    muf2 = t12;
    mur2 = sqrt(dabs(t12*t34));
    q2   = Max(dabs(t12),dabs(t34));
  }

  //msg_Out()<<"Winner for ["<<ID(links[0])[0]<<", "<<ID(links[1])[0]<<"] "
  //	   <<" and ["<<ID(links[2])[0]<<", "<<ID(links[3])[0]<<"] : "
  //	   <<"muF = "<<sqrt(muf2)<<", "<<"muR = "<<sqrt(mur2)<<", "
  //	   <<"Q = "<<sqrt(q2)<<", wbf = "<<wbfkin<<".\n";
  return PDF::CParam(muf2,q2,0.0,mur2,-1);
}

DECLARE_ND_GETTER(WBF_Core_Scale,"WBF",
		  Core_Scale_Setter,Core_Scale_Arguments,true);

Core_Scale_Setter *ATOOLS::Getter
<Core_Scale_Setter,Core_Scale_Arguments,WBF_Core_Scale>::
operator()(const Core_Scale_Arguments &args) const
{
  return new WBF_Core_Scale(args);
}

void ATOOLS::Getter<Core_Scale_Setter,Core_Scale_Arguments,
		    WBF_Core_Scale>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"WBF core scale"; 
}
