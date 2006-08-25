#include "Hadron_Part.H"
#include "Hadronisation_Parameters.H"
#include "Poincare.H"
#include "Message.H"
#include "Random.H"


using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


Hadron_Part::Hadron_Part() :
  m_dtmode(dtm::waves_PS_pop), 
  p_transitions(new Double_Transition_Map),
  p_popper(hadpars.GetPopper())
{ 
  Hadron_WF_Map * allwaves = hadpars.GetMultiplets()->GetWaveFunctions();
  FlavCCMap     * alloweds = (&(hadpars.GetConstituents()->CCMap));
  FlavCCMap_Iterator cc;

  Flavour         had1, had2, cchelp;
  FlavPair        flpair, hadpair;
  double          wt;
  WFcomponent   * waves1, * waves2;

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
	  if (swv1->first->second!=swv2->first->first.Bar()) continue;
	  cchelp = swv2->first->first;
	  if (cchelp.IsDiQuark()) cchelp = cchelp.Bar();
	  cc     = alloweds->find(cchelp);
	  if (cc==alloweds->end() ||
	      cc->second->TotWeight()<1.e-6) {
	    continue;
	  }
	  if (m_dtmode==dtm::waves_PS_pop ||
	      m_dtmode==dtm::waves_PS) 
	    wt  *= cc->second->TotWeight();
	  if (m_dtmode==dtm::waves_PS_pop ||
	      m_dtmode==dtm::PS_pop) 
	    wt  *= p_popper->PopWeight(cchelp);
 
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

Hadron_Part::~Hadron_Part()
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

bool Hadron_Part::FixHHDecay(Cluster * cluster,
			     ATOOLS::Flavour & had1,ATOOLS::Flavour & had2)
{
  double M       = cluster->Mass(), M2 = M*M;
  double m12     = sqr(had1.PSMass()), m22 = sqr(had2.PSMass());
  double ptmax   = sqrt(sqr(M2-m12-m22)-4.*m12*m22)/(2.*M); 
  double pt      = p_popper->SelectPT(ptmax);

  cluster->BoostInCMSAndRotateOnZ();
  double E1      = (M2+m12-m22)/(2.*M);
  double pl1     = sqrt(sqr(E1)-sqr(pt)-m12);
  double cosphi  = cos(2.*M_PI*ran.Get()), sinphi = sqrt(1.-cosphi*cosphi);
  Vec4D  p1      = Vec4D(E1,pt*cosphi,pt*sinphi,pl1);
  Vec4D  p2      = cluster->Momentum()-p1;

  Particle * part;
  Cluster * clus;
  clus = new Cluster();
  clus->SetMomentum(0,p1);
  clus->SetPrev(cluster);
  cluster->SetLeft(clus);

  clus = new Cluster();
  clus->SetMomentum(0,p2);
  clus->SetPrev(cluster);
  cluster->SetRight(clus);

  cluster->RotateAndBoostBack();

  part = new Particle(-1,had1,cluster->GetLeft()->Momentum());
  part->SetNumber();
  part->SetInfo('P');
  part->SetStatus(part_status::active);
  part->SetFinalMass(had1.PSMass());
  control::s_AHAparticles++;
  cluster->GetLeft()->SetSelf(part);
  
  part = new Particle(-1,had2,cluster->GetRight()->Momentum());
  part->SetNumber();
  part->SetInfo('P');
  part->SetStatus(part_status::active);
  part->SetFinalMass(had2.PSMass());
  control::s_AHAparticles++;
  cluster->GetRight()->SetSelf(part);
  
  //   cout<<METHOD<<" check this decay "
  //       <<cluster->Number()<<"["<<cluster->Mass()<<"] --> "<<had1<<" & "<<had2<<endl
  //       <<(*cluster)<<endl;
  
  return true;
}


bool Hadron_Part::MustTransit(Cluster * cluster,Flavour & dec1,Flavour & dec2,
			      const double offset)
{
  dec1 = dec2 = Flavour(kf::none); 

  FlavPair flpair;
  flpair.first  = cluster->GetFlav(1);
  flpair.second = cluster->GetFlav(2);

  double wt(0.), MC(cluster->Mass()), MC2(MC*MC),mmax(0.),mmin(MC), ps,m1,m2;
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
    if (m1+m2>mmax) mmax=m1+m2;
    if (m1+m2<mmin) mmin=m1+m2;
    if (m1+m2>MC-offset) continue;
    ps  = sqrt((MC2-sqr(m1+m2))*(MC2-sqr(m1-m2)));
    wt += ps*decit->second;
  }
  if (MC>mmax+offset || wt==0.) {
    return false;
  }
  dec1 = dec2 = Flavour(kf::none); 
  wt *= ran.Get();
  for (Double_Transition_Siter decit=dtliter->second->begin();
       decit!=dtliter->second->end();decit++) {
    m1  = decit->first.first.Mass();
    m2  = decit->first.second.Mass();
    if (m1+m2>MC-offset) continue;
    ps  = sqrt((MC2-sqr(m1+m2))*(MC2-sqr(m1-m2)));
    wt -= ps*decit->second;
    if (wt<0.) {
      dec1 = decit->first.first;
      dec2 = decit->first.second;
      return true;
    }
  }
  return false;
}

void Hadron_Part::PrintDoubleTransitions() 
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

