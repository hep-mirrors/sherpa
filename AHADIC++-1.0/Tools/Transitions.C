#include "Transitions.H"
#include "Hadronisation_Parameters.H"
#include "Random.H"
#include "Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


All_Single_Transitions::All_Single_Transitions(All_Hadron_Multiplets * multis) :
  m_stmode(stm::masswidthXwaves),
  p_multiplets(multis),
  p_transitions(new Single_Transition_Map),
  m_offset1(hadpars.Get(string("Offset_C->H"))),
  m_offset2(hadpars.Get(string("Offset_C->HH")))
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

bool All_Single_Transitions::MustTransit(Cluster * cluster,Flavour & hadron,
					 const double offset,const bool lighter)
{
  double   mass = cluster->Mass(0);
  FlavPair fpair;
  fpair.first  = cluster->GetFlav(1);
  fpair.second = cluster->GetFlav(2);
  hadron = Flavour(kf::none);
  Single_Transition_Miter stiter = p_transitions->find(fpair);

  if (stiter==p_transitions->end()) {
    msg.Error()<<"Potential error in  All_Single_Transitions::MustTransit :"<<endl
	       <<"   Did not find any entry for "<<fpair.first<<"/"<<fpair.second
	       <<", mass = "<<cluster->Mass()<<";"<<endl<<(*cluster->GetPrev())
	       <<"   will continue and hope for the best."<<endl;
    return false;
  }
  Single_Transition_List * stl = stiter->second;
  double massdiscrim           = max(stl->begin()->Mass()+max(offset,m_offset1),
				     hadpars.GetConstituents()->Mass(fpair.first)+
				     hadpars.GetConstituents()->Mass(fpair.second)+
				     2.*hadpars.GetConstituents()->MinMass()+m_offset2);

  //   if (cluster->GetFlav(1)==Flavour(kf::c) || cluster->GetFlav(1)==Flavour(kf::b) ||
  //       cluster->GetFlav(2)==Flavour(kf::c).Bar() || cluster->GetFlav(2)==Flavour(kf::b).Bar()) {
  //     cout<<METHOD<<"  "<<(*(stl->begin()))<<":"<<stl->begin()->Mass()<<"+"<<offset
  // 	<<" vs. "<<hadpars.GetConstituents()->Mass(fpair.first)<<"+"
  // 	<<hadpars.GetConstituents()->Mass(fpair.second)<<"+"
  // 	<<2.*hadpars.GetConstituents()->MinMass()<<"+"<<m_offset2
  // 	<<" ---> "<<massdiscrim<<endl;
  //   }
  //if ((fpair.first.IsDiQuark() || fpair.second.IsDiQuark()) &&
  //   (fpair.first.Kfcode()==4 || fpair.second.Kfcode()==4))
  // cout<<METHOD<<" for {"<<fpair.first<<","<<fpair.second<<"} : "
  //	<<(*stl->begin())<<" "<<mass<<"  "<<massdiscrim<<endl;
  bool success(true);
  if (mass<massdiscrim) {
    do {
      success=true;
      switch (m_stmode) {
      case (stm::masswidthXwaves) :
	MWTimesWavefunction(cluster,stl,hadron,lighter);
	break;
      case (stm::masswidth) :
	MWCriterion(cluster,stl,hadron,lighter);
	break;
      case (stm::massXwaves) :
	MassTimesWavefunction(cluster,stl,hadron,lighter);
	break;
      case (stm::simplemass) :
      default :
	SimpleMassCriterion(cluster,stl,hadron,lighter);
      }
      if (lighter) { if (hadron.Mass()>mass) success=false; }
    } while (!success);
    return true;
  }
  return false;
}

Return_Value::code
All_Single_Transitions::NextLightest(Cluster * cluster,Flavour & had)
{
  FlavPair fpair;
  fpair.first  = cluster->GetFlav(1);
  fpair.second = cluster->GetFlav(2);
  Single_Transition_Miter stiter = p_transitions->find(fpair);
  if (stiter==p_transitions->end())  {
    msg.Error()<<"Potential error in  "<<METHOD<<" :"<<endl
	       <<"   Did not find any entry for "<<fpair.first<<"/"<<fpair.second
	       <<" in transition list;"<<endl
	       <<"   will continue and hope for the best."<<endl;
    rvalue.IncError(METHOD);
    return Return_Value::Error;
  }

  Single_Transition_List * stl  = stiter->second;
  if (had==Flavour(kf::none)) {
    had   = (*stl->begin());
    return Return_Value::Success;
  }

  Single_Transition_Siter siter;
  for (siter=stl->begin();siter!=stl->end();siter++) if ((*siter)==had) break;
  if (siter==stl->end()) {
    msg.Error()<<"Potential error in  "<<METHOD<<" :"<<endl
	       <<"   Did not find any entry for "<<had<<" in hadron list;"
	       <<(*stl->begin())<<"; "<<stl->size()<<endl
	       <<"   will continue and hope for the best."<<endl;
    rvalue.IncError(METHOD);
    return Return_Value::Error;
  }
  siter++;
  if (siter!=stl->end()) { 
    had=(*siter); 
    return Return_Value::Success;
  }
  rvalue.IncWarning(METHOD);
  return Return_Value::Warning;
}


void All_Single_Transitions::SimpleMassCriterion(Cluster * cluster,
						 Single_Transition_List * stl,
						 Flavour & hadron, bool lighter)
{
  double mass(cluster->Mass(0)), massdisc(1.e6);
  Single_Transition_Siter start(stl->begin());
  if (hadron!=Flavour(kf::none)) {
    for (;start!=stl->end();start++) { if ((*start)==hadron) break; }
    if (lighter && mass<hadron.Mass() && start!=stl->end()) hadron = (*start);
    start++; 
  }
  for (Single_Transition_Siter siter=start;siter!=stl->end();siter++) {
    if (dabs(mass-(*siter).Mass())<massdisc) {
      hadron   = (*siter);
      massdisc = dabs(mass-(*siter).Mass());
    }
  }
}

void All_Single_Transitions::MWCriterion(Cluster * cluster,
					 Single_Transition_List * stl,
					 Flavour & hadron, bool lighter)
{
  double mass(cluster->Mass(0)), crit(0.), minwidth(1.e6), disc(1.e100);
  Single_Transition_Siter start(stl->begin());
  if (hadron!=Flavour(kf::none)) {
    for (;start!=stl->end();start++) { if ((*start)==hadron) break; }
    if (lighter &&  mass<hadron.Mass() && start!=stl->end()) hadron = (*start);
    start++;
  }
  for (Single_Transition_Siter siter=start;siter!=stl->end();siter++) {
    if (!(*siter).IsStable() && ((*siter).Width()<minwidth)) 
      minwidth = (*siter).Width();
  }
  if (minwidth==1.e6) {
    SimpleMassCriterion(cluster,stl,hadron);
    return;
  }

  for (Single_Transition_Siter siter=start;siter!=stl->end();siter++) {
    if ((*siter).IsStable()) 
      crit = sqr(sqr(mass)-sqr((*siter).Mass()))/((*siter).Mass()*minwidth);
    else 
      crit = sqr(sqr(mass)-sqr((*siter).Mass()))/((*siter).Mass()*(*siter).Width());
    if (crit<disc) {
      hadron = (*siter);
      disc   = crit;
    }
  }
}

void All_Single_Transitions::MassTimesWavefunction(Cluster * cluster,
						   Single_Transition_List * stl,
						   Flavour & hadron, bool lighter)
{
  double mass2(sqr(cluster->Mass(0))), disc(1.e100), test;
  Hadron_Wave_Function * waves;
  Single_Transition_Siter start(stl->begin());
  if (hadron!=Flavour(kf::none)) {
    for (;start!=stl->end();start++) { if ((*start)==hadron) break; }
    if (lighter && sqrt(mass2)<hadron.Mass() && start!=stl->end()) hadron = (*start);
    start++;
  }
  for (Single_Transition_Siter siter=start;siter!=stl->end();siter++) {
    waves   = p_multiplets->GetWaveFunction((*siter));
    if (waves!=NULL) {
      test  = dabs(mass2-sqr((*siter).Mass()))/
	sqr(waves->WaveWeight(cluster->GetFlav(1),cluster->GetFlav(2)));
      if (test<disc) {
	hadron = (*siter);
	disc   = test;
      }
    }
  }
}

void All_Single_Transitions::MWTimesWavefunction(Cluster * cluster,
						 Single_Transition_List * stl,
						 Flavour & hadron, bool lighter)
{
  double mass(cluster->Mass(0)), crit(0.), minwidth(1.e6), disc(1.e100);
  Hadron_Wave_Function * waves;
  double WFweight,PSweight;
  Single_Transition_Siter start(stl->begin());
  if (hadron!=Flavour(kf::none)) {
    for (;start!=stl->end();start++) { if ((*start)==hadron) break; }
    if (lighter && mass<hadron.Mass() && start!=stl->end()) hadron = (*start);
    start++;
  }
  for (Single_Transition_Siter siter=start;siter!=stl->end();siter++) {
    if (!(*siter).IsStable() && ((*siter).Width()<minwidth)) 
      minwidth = (*siter).Width();
  }
  if (minwidth==1.e6) {
    MassTimesWavefunction(cluster,stl,hadron);
    return;
  }
  for (Single_Transition_Siter siter=start;siter!=stl->end();siter++) {
    waves   = p_multiplets->GetWaveFunction((*siter));
    if (waves!=NULL) {
      WFweight  = 1./sqr(waves->WaveWeight(cluster->GetFlav(1),cluster->GetFlav(2)));
      if ((*siter).IsStable()) 
	PSweight = sqr(sqr(mass)-sqr((*siter).Mass()))/((*siter).Mass()*1.e-6);
      else 
	PSweight = sqr(sqr(mass)-sqr((*siter).Mass()))/((*siter).Mass()*(*siter).Width());
      if (WFweight*PSweight<disc) {
	hadron = (*siter);
	disc   = WFweight*PSweight;
      }
    }
  }
}

void All_Single_Transitions::PrintSingleTransitions()
{
  Hadron_Wave_Function * wave;
  double wt,mpletwt;
  map<Flavour,double> checkit;
  for (Single_Transition_Miter stiter=p_transitions->begin();
       stiter!=p_transitions->end();stiter++) {
    msg.Debugging()<<"("<<stiter->first.first<<","<<stiter->first.second<<") : "<<endl;
    for (Single_Transition_Siter sit=stiter->second->begin();
	 sit!=stiter->second->end();sit++) {
      wave    = p_multiplets->GetWaveFunction((*sit));
      wt      = sqr(wave->WaveWeight(stiter->first.first,stiter->first.second));
      mpletwt = wave->MultipletWeight();
      msg.Debugging()<<"   "<<(*sit)<<" = "<<wt<<" * "<<mpletwt<<endl;
      if (checkit.find(stiter->first.first)==checkit.end()) 
	checkit[stiter->first.first] = wt;
      else checkit[stiter->first.first] += wt;
      if (checkit.find(stiter->first.second)==checkit.end()) 
       	checkit[stiter->first.second] = wt;
      else checkit[stiter->first.second] += wt;
    }
  }
  msg.Out()<<"In total (summed weights per hadron):"<<endl;
  for (map<Flavour,double>::iterator it=checkit.begin();it!=checkit.end();it++) 
    msg.Out()<<"     -> "<<it->first<<" : "<<it->second<<endl;
  msg.Out()<<"-------- END OF ALL_SINGLE_TRANSITIONS -----"<<endl;  
}



Flavour All_Single_Transitions::GetSU3Pseudoscalar(double mass) {
  Hadron_Multiplet * pseudos = p_multiplets->GetSU3Pseudoscalars();
  double help(100.);
  if (mass>0.) help=mass;

  FlavourSet all;
  FlSetIter  fliter;
  Flavour    flav;
  for (fliter=pseudos->GetElements()->begin();fliter!=pseudos->GetElements()->end();fliter++) {
    flav = (*fliter);
    if (flav.Mass()<help) all.insert(flav);
  }
  if (all.size()==0) all.insert(Flavour(kf::photon));

  fliter = all.begin();
  for (int i=0;i<int(.9999999*all.size()*ran.Get());i++) fliter++;
  return (*fliter);
}

Flavour All_Single_Transitions::GetSU3NeutralPseudoscalar(double mass) {
  Hadron_Multiplet * pseudos = p_multiplets->GetSU3Pseudoscalars();
  double help(100.);
  if (mass>0.) help=mass;

  FlavourSet neutrals;
  FlSetIter  fliter;
  Flavour    flav;
  for (fliter=pseudos->GetElements()->begin();fliter!=pseudos->GetElements()->end();fliter++) {
    flav = (*fliter);
    if (flav.Charge()==0 && flav.Mass()<help) neutrals.insert(flav);
  }
  if (neutrals.size()==0) neutrals.insert(Flavour(kf::photon));

  fliter = neutrals.begin();
  for (int i=0;i<int(.9999999*neutrals.size()*ran.Get());i++) fliter++;
  return (*fliter);
}

