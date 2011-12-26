#include "AHADIC++/Tools/Transitions.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


Single_Transitions::Single_Transitions() :
  m_singletsuppression2(sqr(hadpars->Get("Singlet_Suppression"))),
  p_transitions(new Single_Transition_Map)
{
  Constituents * constituents = hadpars->GetConstituents();

  double mass=100.;
  for (FlavCCMap_Iterator flit=constituents->CCMap.begin();
       flit!=constituents->CCMap.end();flit++) {
    if (flit->second->Mass()<mass && flit->first!=Flavour(kf_gluon)) {
      mass                   = flit->second->Mass();
      m_lightest_constituent = flit->first;
    }
  }    

  Hadron_WF_Map           * allwaves = hadpars->GetMultiplets()->GetWaveFunctions();
  WFcomponent             * waves;
  Single_Transition_Miter   stiter;
  Single_Transition_List  * stl;
  Flavour                   hadron;
  double                    wt, wtprod, costh, sinth;
  int                       flnum,spin,lp;
  for (Hadron_WF_Miter wf=allwaves->begin();wf!=allwaves->end();wf++) {
    hadron = wf->first;
    wt     = wf->second->MultipletWeight();
    if (wt>0.) {
      waves  = wf->second->GetWaves();
      for (WFcompiter singlewave=waves->begin();singlewave!=waves->end();singlewave++) {
	stiter = p_transitions->find((*(singlewave->first)));
	wtprod = wt * sqr(singlewave->second);
	if (!hadron.IsBaryon() && 
	    singlewave->first->first.Kfcode()==singlewave->first->second.Kfcode() &&
	    singlewave->first->first.Kfcode()<4) {
	  flnum  = int(hadron.Kfcode()/100)-10*int(hadron.Kfcode()/1000);
	  lp     = int(hadron.Kfcode()/10000);
	  spin   = int(hadron.Kfcode())-10*int(hadron.Kfcode()/10);
	  hadpars->GetMultiplets()->LookUpAngles(lp,spin,costh,sinth);
	  if (flnum==2) wtprod *= sqrt(m_singletsuppression2*costh*costh+sinth*sinth);
	  if (flnum==3) wtprod *= sqrt(m_singletsuppression2*sinth*sinth+costh*costh);
	  //if (flnum==2 || flnum==3) {
	  // msg_Out()<<"Singlet suppression = "<<sqrt(m_singletsuppression2)<<" for "
	  //	     <<hadron<<" with "<<costh<<"/"<<sinth<<" --> "<<wtprod
	  //	     <<" ("<<singlewave->first->first.Kfcode()<<")\n";
	  //}
	}
	if (wtprod>0. && stiter!=p_transitions->end()) {
	  (*stiter->second)[hadron] += wtprod;
	}
	else {
	  if (wtprod>0.) {
	    stl = new Single_Transition_List;
	    (*stl)[hadron] = wtprod;
	    (*p_transitions)[(*(singlewave->first))] = stl;
	  }
	}
      }
    }
  }
  //PrintSingleTransitions();
}

Single_Transitions::~Single_Transitions()
{
  if (p_transitions) {
    for (Single_Transition_Miter stiter=p_transitions->begin();
	 stiter!=p_transitions->end();stiter++) {
      delete stiter->second;
    }
    delete p_transitions;
  }
}

Flavour Single_Transitions::GetLightestTransition(const Flavour_Pair & fpair) {
  Flavour had = Flavour(kf_none);
  Single_Transition_Miter stiter = p_transitions->find(fpair);
  if (stiter==p_transitions->end())  return had;
  Single_Transition_List * stl  = stiter->second;
  if (stl->empty()) return had;
  return (--stl->end())->first;
}

Flavour Single_Transitions::GetHeaviestTransition(const Flavour_Pair & fpair) {
  Flavour had = Flavour(kf_none);
  Single_Transition_Miter stiter = p_transitions->find(fpair);
  if (stiter!=p_transitions->end()) had =  stiter->second->begin()->first;
  return had;
}

double Single_Transitions::GetLightestMass(const Flavour_Pair & fpair) {
  Flavour had = GetLightestTransition(fpair);
  if (had==Flavour(kf_none)) return -1.;
  return had.HadMass();
}

double Single_Transitions::GetHeaviestMass(const Flavour_Pair & fpair) {
  Flavour had = GetHeaviestTransition(fpair);
  if (had==Flavour(kf_none)) return -1.;
  return had.HadMass();
}

void Single_Transitions::PrintSingleTransitions()
{
  Hadron_Wave_Function * wave;
  double wt,mpletwt;
  map<Flavour,double> checkit;
  for (Single_Transition_Miter stiter=p_transitions->begin();
       stiter!=p_transitions->end();stiter++) {
    msg_Out()<<"("<<stiter->first.first<<","<<stiter->first.second<<") : "<<endl;
    for (Single_Transition_Siter sit=stiter->second->begin();
	 sit!=stiter->second->end();sit++) {
      wave    = hadpars->GetMultiplets()->GetWaveFunction(sit->first);
      wt      = sqr(wave->WaveWeight(stiter->first.first,stiter->first.second));
      mpletwt = wave->MultipletWeight();
      if (mpletwt<=0.) continue;
      msg_Out()<<"   "<<sit->first<<" ("<<sit->first.HadMass()<<" ) = "
	       <<wt<<" * "<<mpletwt<<" = "<<sit->second<<endl;
      if (checkit.find(stiter->first.first)==checkit.end()) 
	checkit[stiter->first.first] = wt;
      else checkit[stiter->first.first] += wt;
      if (checkit.find(stiter->first.second)==checkit.end()) 
	checkit[stiter->first.second] = wt;
      else checkit[stiter->first.second] += wt;
    }
  }
  msg_Out()<<"In total (summed weights per hadron):"<<endl;
  for (map<Flavour,double>::iterator it=checkit.begin();it!=checkit.end();it++) 
    msg_Out()<<"     -> "<<it->first<<" : "<<it->second<<endl;
  msg_Out()<<"-------- END OF ALL_SINGLE_TRANSITIONS -----"<<endl;  
}







std::ostream & AHADIC::operator<<(std::ostream & s, const Double_Transition_List & dtl) {
  for (Double_Transition_List::const_iterator diter=dtl.begin();diter!=dtl.end();diter++)
    s<<"  {"<<diter->first.first<<", "<<diter->first.second<<"} = "<<diter->second<<std::endl;
  return s;
}



Double_Transitions::Double_Transitions() :
  p_transitions(new Double_Transition_Map)
{
  FlavCCMap constituents          = hadpars->GetConstituents()->CCMap;
  Single_Transition_Map * singles = hadpars->GetSingleTransitions()->GetTransitions();
  Single_Transition_List * hads1, * hads2;
  Flavour trip,anti,popped;
  Flavour_Pair flpair,hadpair,hadpair1,wf1,wf2;
  double popwt,weight,norm;
  Double_Transition_List * dtl;
  Double_Transition_Miter dtiter;
  for (FlavCCMap_Iterator iter1=constituents.begin();
       iter1!=constituents.end();iter1++) {
    trip = iter1->first;
    if (trip.IsDiQuark()) trip = trip.Bar();
    flpair.first = wf1.first = trip;
    for (FlavCCMap_Iterator iter2=constituents.begin();
	 iter2!=constituents.end();iter2++) {
      anti = iter2->first;
      if (anti.IsQuark()) anti = anti.Bar();
      flpair.second = wf2.second = anti;
      //if (trip==Flavour(kf_u) && anti==Flavour(kf_u).Bar()) 
      //std::cout<<"  "<<trip<<" / "<<anti<<std::endl;

      for (FlavCCMap_Iterator iter3=constituents.begin();
	   iter3!=constituents.end();iter3++) {
	norm  = 0.;
	popwt = iter3->second->TotWeight();
	if (popwt==0.) continue;
	popped = iter3->first;
	if (popped.IsQuark()) popped = popped.Bar();
	wf1.second = popped;
	wf2.first  = popped.Bar();
	if (singles->find(wf1)!=singles->end() &&
	    singles->find(wf2)!=singles->end()) {
	  hads1 = (*singles)[wf1];
	  hads2 = (*singles)[wf2];
	  for (Single_Transition_Siter haditer1=hads1->begin();
	       haditer1!=hads1->end();haditer1++) {
	    for (Single_Transition_Siter haditer2=hads2->begin();
		 haditer2!=hads2->end();haditer2++) {
	      norm += haditer1->second * haditer2->second;
	    }
	  }
	  for (Single_Transition_Siter haditer1=hads1->begin();
	       haditer1!=hads1->end();haditer1++) {
	    for (Single_Transition_Siter haditer2=hads2->begin();
		 haditer2!=hads2->end();haditer2++) {
	      weight = haditer1->second * haditer2->second * popwt;
	      //if (trip==Flavour(kf_u) && anti==Flavour(kf_u).Bar()) 
	      //	std::cout<<"---add : "<<weight/popwt<<" (pop = "<<popwt<<") "
	      //		 <<"("<<(haditer1->second*haditer2->second)<<") for "
	      //		 <<haditer1->first<<"("<<haditer1->second<<") * "
	      //		 <<haditer2->first<<"("<<haditer2->second<<")."<<std::endl;
	      hadpair.first  = haditer1->first;
	      hadpair.second = haditer2->first;
	      if (weight>0.) {
		dtiter = p_transitions->find(flpair);
		if (dtiter!=p_transitions->end()) {
		  dtl   = dtiter->second;
		  if (dtl->find(hadpair)==dtl->end()) {
		    //if (trip==Flavour(kf_u) && anti==Flavour(kf_u).Bar())
		    // std::cout<<"new entry: "<<hadpair.first<<" / "<<hadpair.second<<"."<<std::endl;
		    (*dtl)[hadpair] = weight;
		  }
		  else {
		    //if (trip==Flavour(kf_u) && anti==Flavour(kf_u).Bar())
		    // std::cout<<"mod entry."<<std::endl;
		    (*dtl)[hadpair] += weight;
		  }
		  //if (trip==Flavour(kf_u) && anti==Flavour(kf_u).Bar() && weight>0.) {
		  // std::cout<<"  insert "<<hadpair.first<<" / "<<hadpair.second
		  //	     <<" --> wt = "<<weight
		  //	     <<" into --> "<<(*dtl)[hadpair];
		  // Double_Transition_List::iterator didi = dtl->find(hadpair);
		  // if (didi==dtl->end()) std::cout<<" not found.";
		  // else std::cout<<" "<<didi->second;
		  // hadpair1.first  = hadpair.second;
		  // hadpair1.second = hadpair.first;
		  // if (dtl->find(hadpair1)!=dtl->end())
		  //   std::cout<<" wt1 = "<<(*dtl)[hadpair1];
		  // std::cout<<std::endl<<(*dtl)<<std::endl;
		  //}
		}
		else {
		  dtl                      = new Double_Transition_List;
		  (*dtl)[hadpair]          = weight;
		  (*p_transitions)[flpair] = dtl;
		  //if (trip==Flavour(kf_u) && anti==Flavour(kf_u).Bar() && weight>0.) 
		  // std::cout<<"  new list "<<hadpair.first<<" / "<<hadpair.second
		  //	     <<", wt = "<<weight<<"."<<std::endl;
		}
	      }
	    }
	  }
	}
      }
    }
  }

  map<Flavour,double> checkit;
  for (Double_Transition_Miter dtiter=p_transitions->begin();
       dtiter!=p_transitions->end();dtiter++) {
    //if (dtiter->first.first==Flavour(kf_u) && 
    //	dtiter->first.second==Flavour(kf_u).Bar()) 
    // std::cout<<"--- Transitions for <"
    //	       <<dtiter->first.first<<", "<<dtiter->first.second<<"> : "<<endl;
    for (Double_Transition_Siter dtit=dtiter->second->begin();
	 dtit!=dtiter->second->end();dtit++) {
      //if (dtiter->first.first==Flavour(kf_u) && 
      //	  dtiter->first.second==Flavour(kf_u).Bar()) 
      //	std::cout<<"   -> {"<<dtit->first.first<<", "<<dtit->first.second<<" }"
      //		 <<" with "<<dtit->second<<endl;
      if (checkit.find(dtit->first.first)==checkit.end()) 
	checkit[dtit->first.first] = dtit->second;
      else checkit[dtit->first.first] += dtit->second;
      if (checkit.find(dtit->first.second)==checkit.end()) 
	checkit[dtit->first.second] = dtit->second;
      else checkit[dtit->first.second] += dtit->second;
    }
    //break;
  }
}

Double_Transitions::~Double_Transitions() {
  if (p_transitions) {
    while (!p_transitions->empty()) {
      delete (p_transitions->begin()->second); 
      p_transitions->erase(p_transitions->begin());
    }
    delete p_transitions;
  }
}

Flavour_Pair Double_Transitions::GetLightestTransition(const Flavour_Pair & fpair) {
  Flavour_Pair pair;
  pair.first = pair.second = Flavour(kf_none);
  Double_Transition_Miter dtiter = p_transitions->find(fpair);
  if (dtiter==p_transitions->end())  return pair;
  Double_Transition_List * dtl  = dtiter->second;
  if (dtl->empty()) return pair;
  return (--dtl->end())->first;
}

Flavour_Pair Double_Transitions::GetHeaviestTransition(const Flavour_Pair & fpair) {
  Flavour_Pair pair;
  pair.first = pair.second = Flavour(kf_none);
  Double_Transition_Miter dtiter = p_transitions->find(fpair);
  if (dtiter!=p_transitions->end()) pair = dtiter->second->begin()->first;
  return pair;
}

double Double_Transitions::GetLightestMass(const Flavour_Pair & fpair) {
  Flavour_Pair pair = GetLightestTransition(fpair);
  if (pair.first==Flavour(kf_none) || pair.second==Flavour(kf_none)) return -1.;
  return pair.first.HadMass()+pair.second.HadMass();
}

double Double_Transitions::GetHeaviestMass(const Flavour_Pair & fpair) {
  Flavour_Pair pair = GetHeaviestTransition(fpair);
  if (pair.first==Flavour(kf_none) || pair.second==Flavour(kf_none)) return -1.;
  return pair.first.HadMass()+pair.second.HadMass();
}

void Double_Transitions::PrintDoubleTransitions() 
{
  map<Flavour,double> checkit;
  for (Double_Transition_Miter dtiter=p_transitions->begin();
       dtiter!=p_transitions->end();dtiter++) {
    msg_Out()<<"Transitions for <"
	     <<dtiter->first.first<<", "<<dtiter->first.second<<"> : "<<endl;
    for (Double_Transition_Siter dtit=dtiter->second->begin();
	 dtit!=dtiter->second->end();dtit++) {
      msg_Out()<<"   -> {"<<dtit->first.first<<", "<<dtit->first.second<<" }"
	       <<" with "<<dtit->second<<endl;
      if (checkit.find(dtit->first.first)==checkit.end()) 
	checkit[dtit->first.first] = dtit->second;
      else checkit[dtit->first.first] += dtit->second;
      if (checkit.find(dtit->first.second)==checkit.end()) 
	checkit[dtit->first.second] = dtit->second;
      else checkit[dtit->first.second] += dtit->second;
    }
  }
  msg_Out()<<"In total (summed weights per hadron):"<<endl;
  for (map<Flavour,double>::iterator it=checkit.begin();it!=checkit.end();it++) {
    msg_Out()<<"     -> "<<it->first<<" : "<<it->second<<endl;
    if (it->first.Bar()!=it->first) {
      if (!IsEqual(it->second,checkit[it->first.Bar()])) 
	msg_Out()<<"   ----->>>>> look: "<<it->first
		 <<" wt/wtbar = "<<it->second<<"/"<<checkit[it->first.Bar()]<<"."<<std::endl;
    }
  }
  msg_Out()<<"-------- END OF ALL_DOUBLE_TRANSITIONS -----"<<endl;
}
