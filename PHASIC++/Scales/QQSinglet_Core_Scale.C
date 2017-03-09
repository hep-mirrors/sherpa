#include "PHASIC++/Scales/Core_Scale_Setter.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Process/Process_Base.H"

#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

namespace PHASIC {

  class QQSinglet_Core_Scale: public PHASIC::Core_Scale_Setter {
  public:

    QQSinglet_Core_Scale(const PHASIC::Core_Scale_Arguments &args):
      Core_Scale_Setter(args) {}

    PDF::CParam Calculate(ATOOLS::Cluster_Amplitude *const ampl);

  };// end of class QQSinglet_Core_Scale

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

PDF::CParam QQSinglet_Core_Scale::Calculate(Cluster_Amplitude *const ampl)
{
  double muf2(-1.), mur2(-1.), q2(-1.);
  size_t nlegs(ampl->Legs().size());
  // need at least 2 incoming and 3 outgoing legs
  if (nlegs<5) {
    msg_Error()<<"ERROR in "<<METHOD<<": "
	       <<"Amplitude with "<<nlegs<<" legs only.\n"
	       <<"   Will return all scales = -1 and hope for the best.\n"; 
    return PDF::CParam(muf2,q2,0.0,mur2,-1);
  }
  std::list<Cluster_Leg *> stronglegs, weaklegs, lightlegs;
  std::vector<Cluster_Leg *> heavylegs;
  Vec4D weakvec(0.,0.,0.,0.);
  // find weak particles and sum momentum
  for (size_t i(0);i<ampl->Legs().size();i++) {
    Cluster_Leg * leg(ampl->Leg(i));
    if (leg->Flav().Strong()) stronglegs.push_back(leg);
    else {
      weaklegs.push_back(leg);
      weakvec += leg->Mom();
    }
  }
  double pperpqcd(0);
  // find strong particles and fing pT heavy strong particles
  for (std::list<Cluster_Leg *>::iterator leg=stronglegs.begin();
       leg!=stronglegs.end();leg++){
    if ((*leg)->Flav().Mass()){
      heavylegs.push_back((*leg));
      pperpqcd+=(*leg)->Mom().PPerp();
    }
    else lightlegs.push_back((*leg));
  }
  // check have 2 heavy quarks
  if (heavylegs.size()!=2) {
    msg_Error()<<"ERROR in "<<METHOD<<": "
	       <<"Amplitude with "<<heavylegs.size()<<" heavy quarks only.\n"
	       <<"   Will return all scales = -1 and hope for the best.\n"; 
    return PDF::CParam(muf2,q2,0.0,mur2,-1);
  }
  //list of top transverse momenta
  double top_perps[2];
  for (int i=0;i<2;i++){
    top_perps[i]=heavylegs[i]->Mom().PPerp();
  }
  // define scales
  muf2 = mur2 = pow(pperpqcd,2);
  q2   = Max(pow(top_perps[0],2),pow(top_perps[1],2));
  return PDF::CParam(muf2,q2,0.0,mur2,-1);
}

DECLARE_ND_GETTER(QQSinglet_Core_Scale,"QQSinglet",
		  Core_Scale_Setter,Core_Scale_Arguments,true);

Core_Scale_Setter *ATOOLS::Getter
<Core_Scale_Setter,Core_Scale_Arguments,QQSinglet_Core_Scale>::
operator()(const Core_Scale_Arguments &args) const
{
  return new QQSinglet_Core_Scale(args);
}

void ATOOLS::Getter<Core_Scale_Setter,Core_Scale_Arguments,
		    QQSinglet_Core_Scale>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"QQSinglet core scale"; 
}

