#include "Transitions.H"
#include "Hadronisation_Parameters.H"
#include "Random.H"
#include "Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


All_Single_Transitions::All_Single_Transitions(All_Hadron_Multiplets * multis) :
  // m_stmode(stm::masswidth), 
  m_stmode(stm::masswidthXwaves),
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

bool All_Single_Transitions::MustTransit(Cluster * cluster,Flavour & hadron,
					 const double offset,const bool lighter)
{
  FlavPair fpair;
  fpair.first  = cluster->GetFlav(1);
  fpair.second = cluster->GetFlav(2);
  hadron = Flavour(kf::none);
  Single_Transition_Miter stiter = p_transitions->find(fpair);
  if (stiter==p_transitions->end()) {
    msg.Error()<<"Potential error in  All_Single_Transitions::MustTransit :"<<endl
	       <<"   Did not find any entry for "<<fpair.first<<"/"<<fpair.second<<";"<<endl
	       <<"   will continue and hope for the best."<<endl;
    return false;
  }
  double mass                  = cluster->Mass(0);
  Single_Transition_List * stl = stiter->second;
  double massdiscrim           = stl->begin()->Mass()+offset; 

  if ((fpair.first.IsDiQuark() || fpair.second.IsDiQuark()) &&
      (fpair.first.Kfcode()==4 || fpair.second.Kfcode()==4))
    cout<<METHOD<<" for {"<<fpair.first<<","<<fpair.second<<"} : "
	<<(*stl->begin())<<" "<<mass<<"  "<<massdiscrim<<endl;
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
    Flavour flav1 = cluster->GetFlav(1), flav2 = cluster->GetFlav(2);
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
  double weight;
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
      weight  = 1./sqr(waves->WaveWeight(cluster->GetFlav(1),cluster->GetFlav(2)));
      if ((*siter).IsStable()) 
	weight *= sqr(sqr(mass)-sqr((*siter).Mass()))/((*siter).Mass()*1.e-6);
      else 
	weight *= sqr(sqr(mass)-sqr((*siter).Mass()))/((*siter).Mass()*(*siter).Width());
      if (weight<disc) {
	hadron = (*siter);
	disc   = weight;
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

All_Double_Transitions::All_Double_Transitions(All_Hadron_Multiplets * multis) :
  m_dtmode(dtm::waves_PS), 
  p_multiplets(multis),
  p_transitions(new Double_Transition_Map),
  p_alloweds(&(hadpars.GetConstituents()->CCMap))
{
  Flavour   had1,had2,cchelp;
  FlavPair  flpair,hadpair;
  FlavCCMap_Iterator cc;
  double    wt;
  WFcomponent   * waves1, * waves2;
  Hadron_WF_Map * allwaves = p_multiplets->GetWaveFunctions();

  Double_Transition_Miter   dtiter;
  Double_Transition_List  * dtl;

  for (Hadron_WF_Miter wf1=allwaves->begin();wf1!=allwaves->end();wf1++) {
    had1          = wf1->first;
    waves1        = wf1->second->GetWaves();
    hadpair.first = had1;
    for (Hadron_WF_Miter wf2=allwaves->begin();wf2!=allwaves->end();wf2++) {
      had2           = wf2->first;
      waves2         = wf2->second->GetWaves();
      wt             = wf1->second->MultipletWeight()*wf2->second->MultipletWeight();
      hadpair.second = had2;
      for (WFcompiter swv1=waves1->begin();swv1!=waves1->end();swv1++) {
	flpair.first = swv1->first->first;
	for (WFcompiter swv2=waves2->begin();swv2!=waves2->end();swv2++) {
	  flpair.second = swv2->first->second;
	  if (swv1->first->second!=swv2->first->first.Bar() ||
	      swv1->first->second==Flavour(kf::gluon) ||
	      swv1->first->second==Flavour(kf::c) ||
	      swv1->first->second==Flavour(kf::b)) continue;
	  cchelp = swv2->first->first;
	  if (cchelp.IsDiQuark()) cchelp = cchelp.Bar();
	  cc     = p_alloweds->find(cchelp);
	  if (cc==p_alloweds->end() ||
	      cc->second->TotWeight()<1.e-6) continue;
	  wt    *= cc->second->TotWeight();
	  dtiter = p_transitions->find(flpair);
	  if (dtiter!=p_transitions->end()) {
	    (*dtiter->second)[hadpair] += wt*sqr(swv1->second*swv2->second);
	  }
	  else {
	    if (wt*sqr(swv1->second*swv2->second)>0.) {
	      dtl                      = new Double_Transition_List;
	      (*dtl)[hadpair]          = wt*sqr(swv1->second*swv2->second);
	      (*p_transitions)[flpair] = dtl;
	    }
	  }
	}
      }
    }
  }
}

All_Double_Transitions::~All_Double_Transitions()
{
  if (p_transitions) {
    for (Double_Transition_Miter dtiter=p_transitions->begin();
	 dtiter!=p_transitions->end();dtiter++) {
      dtiter->second->clear();
      delete dtiter->second;
    }
    p_transitions->clear();
    p_transitions = NULL;
  }
}  

bool All_Double_Transitions::MustTransit(Cluster * cluster,Flavour & dec1,Flavour & dec2,
					 const double offset, bool lighter)
{
  dec1 = dec2 = Flavour(kf::none); 

  FlavPair flpair;
  flpair.first  = cluster->GetFlav(1);
  flpair.second = cluster->GetFlav(2);

  double wt(0.), cmass2(cluster->Mass2()),mmax(0.),mmin(sqrt(cmass2)), ps,m1,m2;
  Double_Transition_Miter dtliter = p_transitions->find(flpair);
  if (dtliter==p_transitions->end()) {
    msg.Error()<<"ERROR in "<<METHOD<<" : "<<endl
	       <<"   No transition table found for "<<flpair.first<<"/"<<flpair.second<<endl
	       <<"   Return 'false' and hope for the best."<<std::endl;
    return false;
  }
  for (Double_Transition_Siter decit=dtliter->second->begin();
       decit!=dtliter->second->end();decit++) {
    m1  = decit->first.first.Mass();
    m2  = decit->first.second.Mass();
    ps  = (cmass2-sqr(m1+m2)>0 ? sqrt((cmass2-sqr(m1+m2))*(cmass2-sqr(m1-m2))):0.);
    if (ps>0.) {
      wt += ps*decit->second;
      if (m1+m2>mmax) { 
	mmax=m1+m2;
	dec1 = decit->first.first;
	dec2 = decit->first.second;
      }
    }
    if (m1+m2<mmin) mmin=m1+m2;
  }
  if (cmass2>sqr(mmax+offset) || wt==0.) return false;
  cout<<"       "<<METHOD<<" : "<<sqrt(cmass2)<<" "<<mmax<<" + "<<offset<<" = "<<wt
      <<" for "<<dec1<<" & "<<dec2<<" from "<<flpair.first<<" & "<<flpair.second<<endl;
  dec1 = dec2 = Flavour(kf::none); 
  wt *= ran.Get();
  for (Double_Transition_Siter decit=dtliter->second->begin();
       decit!=dtliter->second->end();decit++) {
    m1  = decit->first.first.Mass();
    m2  = decit->first.second.Mass();
    ps  = (cmass2-sqr(m1+m2)>0 ? sqrt((cmass2-sqr(m1+m2))*(cmass2-sqr(m1-m2))):0.);
    if (ps>0.) wt -= ps*decit->second;
    if (wt<0.) {
      dec1 = decit->first.first;
      dec2 = decit->first.second;
      cout<<"       "<<METHOD<<" select cluster ("<<flpair.first<<","<<flpair.second<<"),"
	  <<" m = "<<sqrt(cmass2)<<" -> "<<dec1<<" & "<<dec2<<endl;
      return true;
    }
  }
  return false;
}

bool All_Double_Transitions::IsoDecay(Cluster * cluster,Flavour & dec1,Flavour & dec2)
{
  FlavPair flpair;
  flpair.first  = cluster->GetFlav(1);
  flpair.second = cluster->GetFlav(2);
  Double_Transition_Miter dtliter = p_transitions->find(flpair);
  if (dtliter!=p_transitions->end()) {
    double PS=0., cmass2=cluster->Mass2(), ps,m1,m2;
    for (Double_Transition_Siter decit=dtliter->second->begin();
	 decit!=dtliter->second->end();decit++) {
      m1  = decit->first.first.Mass();
      m2  = decit->first.second.Mass();
      ps  = (cmass2-sqr(m1+m2))*(cmass2-sqr(m1-m2));
      if (ps>0.) ps = sqrt(ps);
            else ps = 0.;
      PS += ps*decit->second;
    }
    PS *= ran.Get();
    for (Double_Transition_Siter decit=dtliter->second->begin();
	 decit!=dtliter->second->end();decit++) {
      m1  = decit->first.first.Mass();
      m2  = decit->first.second.Mass();
      ps  = (cmass2-sqr(m1+m2))*(cmass2-sqr(m1-m2));
      if (ps>0.) ps = sqrt(ps);
            else ps = 0.;
      PS -= ps*decit->second;
      if (PS<0.) { 
	dec1 = decit->first.first;
	dec2 = decit->first.second;
	return true;
      }
    }
  }
  else PrintDoubleTransitions();
  dec1 = dec2 = Flavour(kf::none);
  return false;
}

void All_Double_Transitions::PrintDoubleTransitions() 
{
  map<Flavour,double> checkit;
  for (Double_Transition_Miter dtiter=p_transitions->begin();
       dtiter!=p_transitions->end();dtiter++) {
    msg.Out()<<"Transitions for <"
	     <<dtiter->first.first<<", "<<dtiter->first.second<<"> : "<<endl;
    for (Double_Transition_Siter dtit=dtiter->second->begin();
	 dtit!=dtiter->second->end();dtit++) {
      msg.Out()<<"   -> {"<<dtit->first.first<<", "<<dtit->first.second<<" }"
	       <<" with "<<dtit->second<<endl;
      if (checkit.find(dtit->first.first)==checkit.end()) 
	checkit[dtit->first.first] = dtit->second;
      else checkit[dtit->first.first] += dtit->second;
      if (checkit.find(dtit->first.second)==checkit.end()) 
	checkit[dtit->first.second] = dtit->second;
      else checkit[dtit->first.second] += dtit->second;
    }
  }
  msg.Out()<<"In total (summed weights per hadron):"<<endl;
  for (map<Flavour,double>::iterator it=checkit.begin();it!=checkit.end();it++) 
    msg.Out()<<"     -> "<<it->first<<" : "<<it->second<<endl;
  msg.Out()<<"-------- END OF ALL_DOUBLE_TRANSITIONS -----"<<endl;
}
