#include "AHADIC++/Tools/Transitions.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


Single_Transitions::Single_Transitions() :
  p_transitions(new Single_Transition_Map)
{
  Constituents * constituents = hadpars.GetConstituents();

  double mass=100.;
  for (FlavCCMap_Iterator flit=constituents->CCMap.begin();
       flit!=constituents->CCMap.end();flit++) {
    if (flit->second->Mass()<mass && flit->first!=Flavour(kf_gluon)) {
      mass                   = flit->second->Mass();
      m_lightest_constituent = flit->first;
    }
  }    

  Hadron_WF_Map           * allwaves = hadpars.GetMultiplets()->GetWaveFunctions();
  WFcomponent             * waves;
  Single_Transition_Miter   stiter;
  Single_Transition_List  * stl;
  Flavour                   hadron;
  double                    wt, wtprod;
  for (Hadron_WF_Miter wf=allwaves->begin();wf!=allwaves->end();wf++) {
    hadron = wf->first;
    waves  = wf->second->GetWaves();
    wt     = wf->second->MultipletWeight();
    //std::cout<<"Try to include new hadron : "<<hadron<<" --> "<<wt<<"."<<std::endl;
    for (WFcompiter singlewave=waves->begin();singlewave!=waves->end();singlewave++) {
      stiter = p_transitions->find((*(singlewave->first)));
      wtprod = wt * sqr(singlewave->second);
      if (wtprod>0. && stiter!=p_transitions->end()) {
	(*stiter->second)[hadron] += wtprod;
      }
      else {
	//std::cout<<"New Single Transition List for "
	//	 <<singlewave->first->first<<" "<<singlewave->first->second<<" --> "
	//	 <<singlewave->second<<" in "<<hadron<<std::endl;
	if (wtprod>0.) {
	  stl = new Single_Transition_List;
	  (*stl)[hadron] = wtprod;
	  (*p_transitions)[(*(singlewave->first))] = stl;
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

bool Single_Transitions::MustDesintegrate(Cluster * cluster,Flavour & had1,Flavour & had2)
{
  if (!(cluster->GetTrip()->m_flav.IsDiQuark() && 
	cluster->GetAnti()->m_flav.IsDiQuark()))   return false;

  double  mass = cluster->Mass();
  Flavour_Pair fpair1,fpair2;

  fpair1.first = cluster->GetTrip()->m_flav; fpair1.second = m_lightest_constituent.Bar();
  fpair2.first = m_lightest_constituent; fpair2.second = cluster->GetAnti()->m_flav;

  had1 = GetLightestTransition(fpair1);
  had2 = GetLightestTransition(fpair2);
  if (had1.HadMass()+had2.HadMass()<mass) return false;

  int kfc1 = int(cluster->GetTrip()->m_flav.Kfcode());
  int kfc2 = int(cluster->GetAnti()->m_flav.Kfcode());
  
  fpair1.second  = Flavour(int(kfc1/1000)).Bar(); 
  fpair2.second  = Flavour(int((kfc1-int(kfc1/1000)*1000)/100)).Bar();
  if (ran.Get()>0.5) {
    fpair1.first = Flavour(int(kfc2/1000)); 
    fpair2.first = Flavour(int((kfc2-int(kfc2/1000)*1000)/100));
  }
  else {
    fpair2.first = Flavour(int(kfc2/1000)); 
    fpair1.first = Flavour(int((kfc2-int(kfc2/1000)*1000)/100));
  }
  had1 = GetLightestTransition(fpair1);
  had2 = GetLightestTransition(fpair2);
  if (had1.HadMass()+had2.HadMass()<mass) return true;
  return false;
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
      wave    = hadpars.GetMultiplets()->GetWaveFunction(sit->first);
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










Double_Transitions::Double_Transitions() :
  p_transitions(new Double_Transition_Map)
{
  FlavCCMap constituents          = hadpars.GetConstituents()->CCMap;
  Single_Transition_Map * singles = hadpars.GetSingleTransitions()->GetTransitions();
  Single_Transition_List * hads1, * hads2;
  Flavour trip,anti,popped;
  Flavour_Pair flpair,hadpair,wf1,wf2;
  double popwt,weight,norm;
  bool   found;
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
      //std::cout<<"  "<<trip<<" / "<<anti<<std::endl;

      for (FlavCCMap_Iterator iter3=constituents.begin();
	   iter3!=constituents.end();iter3++) {
	norm = 0.;
	popwt = iter3->second->TotWeight();
	//std::cout<<"  pop "<<iter3->first<<" ( wt = "<<popwt<<") "
	//	 <<" between "<<trip<<" and "<<anti<<"."<<std::endl;
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
	      weight = haditer1->second * haditer2->second/norm * popwt;
	      //std::cout<<"   add : "<<weight/popwt<<" for "
	      //	       <<haditer1->first<<"("<<haditer1->second<<") * "
	      //	       <<haditer2->first<<"("<<haditer2->second<<")."<<std::endl;
	      hadpair.first  = haditer1->first;
	      hadpair.second = haditer2->first;
	      if (weight>0.) {
		dtiter = p_transitions->find(flpair);
		if (dtiter!=p_transitions->end()) {
		  dtl   = dtiter->second;
		  found = false;
		  for (Double_Transition_Siter dit=dtl->begin();
		       dit!=dtl->end();dit++) {
		    if (hadpair.first==dit->first.first &&
			hadpair.second==dit->first.second) {
		      dit->second += weight;
		      found = true;
		      break;
		    }
		  }
		  if (!found) (*dtl)[hadpair] = weight;
		}
		else {
		  dtl                      = new Double_Transition_List;
		  (*dtl)[hadpair]          = weight;
		  (*p_transitions)[flpair] = dtl;
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
    //msg_Out()<<"Transitions for <"
    //	     <<dtiter->first.first<<", "<<dtiter->first.second<<"> : "<<endl;
    for (Double_Transition_Siter dtit=dtiter->second->begin();
	 dtit!=dtiter->second->end();dtit++) {
      //msg_Out()<<"   -> {"<<dtit->first.first<<", "<<dtit->first.second<<" }"
      //	       <<" with "<<dtit->second<<endl;
      if (checkit.find(dtit->first.first)==checkit.end()) 
	checkit[dtit->first.first] = dtit->second;
      else checkit[dtit->first.first] += dtit->second;
      if (checkit.find(dtit->first.second)==checkit.end()) 
	checkit[dtit->first.second] = dtit->second;
      else checkit[dtit->first.second] += dtit->second;
    }
    break;
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
  for (map<Flavour,double>::iterator it=checkit.begin();it!=checkit.end();it++) 
    msg_Out()<<"     -> "<<it->first<<" : "<<it->second<<endl;
  msg_Out()<<"-------- END OF ALL_DOUBLE_TRANSITIONS -----"<<endl;
}
