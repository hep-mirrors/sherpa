#include "Transitions.H"
#include "Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


All_Single_Transitions::All_Single_Transitions(All_Hadron_Multiplets * multis) :
  m_stmode(stm::massXwaves), 
  p_multiplets(multis),
  p_transitions(new Single_Transition_Map)
{
  Hadron_WF_Map * allwaves = p_multiplets->GetWaveFunctions();
  WFcomponent * waves;
  Single_Transition_Miter   stiter;
  Single_Transition_List  * stl;
  Flavour hadron;
  for (Hadron_WF_Miter wf=allwaves->begin();wf!=allwaves->end();wf++) {
    hadron = wf->first;
    waves  = wf->second->GetWaves();
    for (WFcompiter singlewave=waves->begin();singlewave!=waves->end();singlewave++) {
      stiter=p_transitions->find((*(singlewave->first)));
      if (stiter!=p_transitions->end()) {
	stiter->second->insert(hadron);
      }
      else {
	stl = new Single_Transition_List;
	stl->insert(hadron);
	(*p_transitions)[(*(singlewave->first))] = stl;
      }
    }
  }
}

All_Single_Transitions::~All_Single_Transitions()
{
  if (p_transitions) {
    for (Single_Transition_Miter stiter=p_transitions->begin();
	 stiter!=p_transitions->end();stiter++) {
      stiter->second->clear();
      delete stiter->second;
    }
    p_transitions->clear();
    p_transitions = NULL;
  }
}

bool All_Single_Transitions::MustTransit(Cluster * cluster,Flavour & hadron,const double offset)
{
  FlavPair fpair;
  fpair.first  = cluster->GetFlav(1);
  fpair.second = cluster->GetFlav(2);
  Single_Transition_Miter stiter = p_transitions->find(fpair);
  if (stiter==p_transitions->end()) {
    msg.Error()<<"Potential error in  All_Single_Transitions::MustTransit :"<<endl
	       <<"   Did not find any entry for "<<fpair.first<<"/"<<fpair.second<<";"<<endl
	       <<"   will continue and hope for the best."<<endl;
  }
  double mass                  = cluster->Mass(0);
  Single_Transition_List * stl = stiter->second;
  double massdiscrim           = stl->begin()->Mass()+offset; 
  //cout<<"         Must transit "<<mass<<" "<<(*stl->begin())<<" "<<offset;
  if (mass<massdiscrim) {
    //cout<<" -> IN "<<endl;
    switch (m_stmode) {
    case (stm::massXwaves) :
      MassTimesWavefunction(cluster,stl,hadron);
      break;
    case (stm::simplemass) :
    default :
      SimpleMassCriterion(cluster,stl,hadron);
    }
    return true;
  }
  //cout<<" -> OUT "<<endl;
  hadron = Flavour(kf::none);
  return false;
}

bool All_Single_Transitions::NextLightest(Cluster * cluster,Flavour & had)
{
  FlavPair fpair;
  fpair.first  = cluster->GetFlav(1);
  fpair.second = cluster->GetFlav(2);
  Single_Transition_Miter stiter = p_transitions->find(fpair);
  if (stiter==p_transitions->end()) {
    msg.Error()<<"Potential error in  All_Single_Transitions::NextLightest :"<<endl
	       <<"   Did not find any entry for "<<fpair.first<<"/"<<fpair.second<<";"<<endl
	       <<"   will continue and hope for the best."<<endl;
  }
  Single_Transition_List * stl  = stiter->second;
  Single_Transition_Siter siter;
  for (siter=stl->begin();siter!=stl->end();siter++) {
    if ((*siter)==had) break;
  }
  if (siter==stl->end()) {
    msg.Error()<<"Potential error in  All_Single_Transitions::NextLightest :"<<endl
	       <<"   Did not find any entry for "<<had<<" in hadron list;"<<(*stl->begin())<<"; "<<stl->size()<<endl
	       <<"   will abort."<<endl;
    abort();
  }
  siter++;
  if (siter!=stl->end()) { 
    had=(*siter); 
    return true; 
  }
  return false;
}


void All_Single_Transitions::SimpleMassCriterion(Cluster * cluster,Single_Transition_List * stl,Flavour & hadron)
{
  double mass  = cluster->Mass(0), massdisc = 1.e6;
  for (Single_Transition_Siter siter=stl->begin();siter!=stl->end();siter++) {
    if (dabs(mass-(*siter).Mass())<massdisc) {
      hadron   = (*siter);
      massdisc = dabs(mass-(*siter).Mass());
    }
  }
}

void All_Single_Transitions::MassTimesWavefunction(Cluster * cluster,Single_Transition_List * stl,Flavour & hadron)
{
  double mass2  = sqr(cluster->Mass(0)), massdisc = 0.;
  Flavour testhad;
  Hadron_Wave_Function * waves;
  double weight, test;
  for (Single_Transition_Siter siter=stl->begin();siter!=stl->end();siter++) {
    testhad = (*siter);
    waves   = p_multiplets->GetWaveFunction(testhad);
    if (waves!=NULL) {
      weight  = sqr(waves->Weight(cluster->GetFlav(1),cluster->GetFlav(2)));
      test    = weight/dabs(mass2-sqr(testhad.Mass()));
      if (test>massdisc) {
	hadron   = testhad;
	massdisc = test;
      }
    }
  }
}

void All_Single_Transitions::PrintSingleTransitions()
{
  for (Single_Transition_Miter stiter=p_transitions->begin();
       stiter!=p_transitions->end();stiter++) {
    msg.Out()<<"("<<stiter->first.first<<","<<stiter->first.second<<") : ";
    for (Single_Transition_Siter sit=stiter->second->begin();
	 sit!=stiter->second->end();sit++) msg.Out()<<(*sit)<<" ";
    msg.Out()<<endl;
  }
  msg.Out()<<"-------- END OF ALL_SINGLE_TRANSITIONS -----"<<endl;  
}
