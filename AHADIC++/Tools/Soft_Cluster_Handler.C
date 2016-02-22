#include "AHADIC++/Tools/Soft_Cluster_Handler.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Phys/Momenta_Stretcher.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Return_Value.H"
#include <sys/stat.h>

using namespace AHADIC;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Soft_Cluster_Handler::Soft_Cluster_Handler(bool ana) :
  p_singletransitions(hadpars->GetSingleTransitions()), 
  p_doubletransitions(hadpars->GetDoubleTransitions()),
  p_as((Strong_Coupling*)s_model->GetScalarFunction("strong_cpl")),
  m_transoffset(hadpars->Get(std::string("Offset_C->H"))), 
  m_decayoffset(hadpars->Get(std::string("Offset_C->HH"))), 
  m_minmass(2.1*hadpars->GetConstituents()->MinMass()),
  m_kappa(hadpars->Get(std::string("MassExponent_C->H"))), 
  m_lambda(hadpars->Get(std::string("WidthExponent_C->H"))), 
  m_chi(hadpars->Get(std::string("MassExponent_C->HH"))), 
  m_pt2max(sqr(hadpars->Get(string("ptmax")))),
  m_pt2maxfac(sqr(hadpars->Get(std::string("ptmax_factor")))),
  m_pt02(hadpars->Get(std::string("pt02"))), 
  m_transitions(0), m_dtransitions(0), m_decays(0), 
  m_forceddecays(0), m_lists(0), m_update(0),
  m_ana(ana), m_out(false)
{
  if (m_ana) {
    m_histograms[string("PT_HH")]  = new Histogram(0,0.,10.,100);
    m_histograms[string("PT2_HH")] = new Histogram(0,0.,100.,2000);
    m_histograms[string("Method_TreatCL")] = new Histogram(0,-0.5,1.5,2);
    m_histograms[string("Method_CLFT")] = new Histogram(0,-0.5,1.5,2);
    m_histograms[string("Method_CheckCL")] = new Histogram(0,-0.5,4.5,5);
    m_histograms[string("DecayWeight")] = new Histogram(0, -1.5,2.5,31);
    m_histograms[string("TransWeight")] = new Histogram(0, -1.5,2.5,31);
    //m_histograms[string("MassTransition")]       = new Histogram(0,0.,8.,100);
    //m_histograms[string("HadronMassTransition")] = new Histogram(0,0.,8.,100);
  }
}

Soft_Cluster_Handler::~Soft_Cluster_Handler() 
{
  if (m_ana) {
    Histogram * histo;
    for (map<string,Histogram *>::iterator hit=m_histograms.begin();
         hit!=m_histograms.end();hit++) {
      histo = hit->second;
      delete histo;
    }
    m_histograms.clear();
  }
}

/// Master function
bool Soft_Cluster_Handler::TreatClusterList(Cluster_List * clin, Blob * blob)
{
  // if (m_out) 
    //msg_Out()<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  // 	     <<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
  // 	     <<"++++++ "<<METHOD<<"("<<clin->size()<<" clusters):\n";
  // checks for transitions to hadrons and attaches them, if neccessary
  if (!CheckListForTreatment(clin) && m_out) { // HS: FIXME this currently never evaluates to true --- BUG?
    // msg_Out()<<"++++++ No hadrons produced.  Just continue.\n"
    // 	     <<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    // 	     <<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    //
    if (m_ana) m_histograms["Method_TreatCL"]->Insert(0);
    return false; //
  }
  if (m_ana) m_histograms["Method_TreatCL"]->Insert(1);
  // if (m_out) 
  //   msg_Out()<<"++++++ Hadrons produced, will attach to blob.\n";
  return AttachHadronsToBlob(clin,blob); // That method always returns true
}

bool Soft_Cluster_Handler::CheckListForTreatment(Cluster_List * clin) {
  // Iterate over all clusters and check number of hadrons from
  // possible decays, administered by Checkcluster.
  // Overall result 0 means no hadrons replacing clusters.
  int count(0);
  for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
    Cluster * cluster(*cit);
    if (cluster==NULL || !cluster->Active()) continue;
    count += CheckCluster(cluster);
  }
  if (count==0) {
    if (m_ana) m_histograms["Method_CLFT"]->Insert(0);
    return false;
  }
  if (m_ana) m_histograms["Method_CLFT"]->Insert(1);
  return true;
}

int Soft_Cluster_Handler::CheckCluster(Cluster * cluster) {
  // if (m_out) 
  //   msg_Out()<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++"
  // 	     <<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
  // 	     <<"++++++ "<<METHOD<<" for:\n"<<(*cluster);
  cluster->clear();
  // Init target hadron flavours
  Flavour haddec1(Flavour(kf_none)), haddec2(Flavour(kf_none));
  Flavour hadtrans(Flavour(kf_none));
  // Get weights. Methods also overwrite Falvours initialised above
  double decayweight(DecayWeight(cluster,haddec1,haddec2));
  double transweight(TransformWeight(cluster,hadtrans));
  if (m_ana) m_histograms["DecayWeight"]->Insert(decayweight);
  if (m_ana) m_histograms["TransWeight"]->Insert(transweight);
  // if (m_out) 
  //   msg_Out()<<"++++++ "<<METHOD<<"["<<cluster->Mass()<<" "
  // 	     <<"("<<cluster->GetTrip()->m_flav<<" + "
  // 	     <<cluster->GetAnti()->m_flav<<"] --> "
  // 	     <<"(dec = "<<decayweight<<", "
  // 	     <<"trans = "<<transweight<<").\n";
  
  if (decayweight == 0 && transweight ==0) {
    // no decay and no transition: both weight equal 0.
    // if (m_out) msg_Out()<<"++++++ no decay.\n"
    cluster->clear();
    if (m_ana) m_histograms["Method_CheckCL"]->Insert(0);
    return 0;
  }

  if (decayweight>0.) {
    if (transweight>0.) {
      // Competition between decay and transition
      double totweight(decayweight+transweight);
      if (totweight*ran->Get()*0.9999999 < decayweight) {    
	// Decay wins
	cluster->push_back(haddec1);
	cluster->push_back(haddec2);
	m_decays      += 1;
	// if (m_out) msg_Out()<<"++++++ decays to "<<haddec1<<" + "<<haddec2<<".\n"
        if (m_ana) m_histograms["Method_CheckCL"]->Insert(1);
      }
      else {
	// Transition wins
	cluster->push_back(hadtrans);
	cluster->push_back(Flavour(kf_photon));
	// if (m_out) msg_Out()<<"++++++ decays to "<<hadtrans<<" + photon.\n"
	m_transitions += 1;
        if (m_ana) m_histograms["Method_CheckCL"]->Insert(2);
      }
    }
    else if (transweight<=0.) {
      // regular decay, and no simple transition open - so no competition
      cluster->push_back(haddec1);
      cluster->push_back(haddec2);
      // if (m_out) msg_Out()<<"++++++ decays to "<<haddec1<<" + "<<haddec2<<".\n"
      m_decays      += 1;
      if (m_ana) m_histograms["Method_CheckCL"]->Insert(3);
    }
  }
  else if (decayweight<0.) {
    // transition - if neccessary enforced.
    if (transweight<=0.) transweight = TransformWeight(cluster,hadtrans,true);
    cluster->push_back(hadtrans);
    cluster->push_back(Flavour(kf_photon));
    // if (m_out)  msg_Out()<<"++++++ decays to "<<hadtrans<<" + photon.\n"
    m_transitions += 1;
    if (m_ana) m_histograms["Method_CheckCL"]->Insert(4);
  }
  if (cluster->size()!=2) {
    msg_Error() << "Two-body cluster decay failed, retrying event\n";
    throw Return_Value::Retry_Event;
  }
  return 2;
}

bool Soft_Cluster_Handler::
AttachHadronsToBlob(Cluster_List * clin,Blob * blob)
{
  // Iterate over cluster list, find those that decayed
  Cluster_Iterator cit(clin->begin());
  Cluster  * cluster;
  while (cit!=clin->end()) {
    cluster = (*cit);
    if (cluster->size()==2) {
      // Fix kinematics, remove cluster from list
      FixHHDecay(cluster,blob,(*cluster)[0],(*cluster)[1]);
      delete cluster->GetTrip();
      delete cluster->GetAnti();
      delete cluster;
      cit = clin->erase(cit);
    }
    else cit++;
  }
  // if (m_out) 
  //   msg_Out()<<"++++++ "<<METHOD<<" was successful:"
  // 	     <<blob->CheckMomentumConservation()<<"\n"
  return true;
}

double Soft_Cluster_Handler::
TransformWeight(Cluster * cluster,Flavour & hadron,const bool & enforce)
{
  Flavour_Pair fpair;
  fpair.first  = cluster->GetTrip()->m_flav;
  fpair.second = cluster->GetAnti()->m_flav;

  if (fpair.first.IsDiQuark() && fpair.second.IsDiQuark()) return 0.;
  if (p_singletransitions->GetTransitions()->find(fpair)==
      p_singletransitions->GetTransitions()->end()) {
    msg_Error()<<"Error in "<<METHOD<<" for cluster\n"
	       <<(*cluster)
	       <<"   illegal flavour combination.\n"
	       <<"   Will return 0 and hope for the best.\n";
    return 0.;
  }

  double MC(cluster->Mass());
  double critM(p_singletransitions->GetLightestMass(fpair)*(1.-m_transoffset)+
	       p_singletransitions->GetHeaviestMass(fpair)*m_transoffset);
  if (!enforce && MC>critM) {
    hadron = Flavour(kf_none);
    return 0.;
  }

  Single_Transition_Miter stiter = 
    p_singletransitions->GetTransitions()->find(fpair);
  if (stiter==p_singletransitions->GetTransitions()->end()) {
    hadron = Flavour(kf_none);
    return -1.;
  }
  Single_Transition_List * stl(stiter->second);

  double wt(0.),totweight(0.);
  for (Single_Transition_Siter siter=stl->begin();siter!=stl->end();siter++) {
    if (siter->first.Mass()>MC) continue;
    wt  = TransformKin(MC,siter->first,enforce);
    wt *= siter->second;
    totweight += wt;
  }

  double disc(totweight * 0.9999999999*ran->Get());
  for (Single_Transition_Siter siter=stl->begin();siter!=stl->end();siter++) {
    if (siter->first.Mass()>MC) continue;
    wt  = TransformKin(MC,siter->first,enforce);
    wt *= siter->second;
    disc -= wt;
    if (disc<=0.) {
      hadron = siter->first;
      break;
    }
  }
  return totweight/(16.*M_PI*MC)/137.;
}

double Soft_Cluster_Handler::
TransformKin(const double MC,const Flavour & flav,const bool & enforce) {
  double mass2(sqr(flav.HadMass()));
  double width2(sqr(Max(flav.Width(),1.e-8)));
  return
    pow(sqr(mass2)/(sqr(MC*MC-mass2) + mass2*width2),m_kappa) * 
    pow(mass2*width2/(sqr(MC*MC-mass2) + mass2*width2),1.+m_lambda);
}



double Soft_Cluster_Handler::
DecayWeight(Cluster * cluster,Flavour & had1,Flavour & had2)
{
  Flavour_Pair flpair;
  flpair.first  = cluster->GetTrip()->m_flav;
  flpair.second = cluster->GetAnti()->m_flav;
  Double_Transition_Miter dtliter = 
    p_doubletransitions->GetTransitions()->find(flpair);
  if (dtliter==p_doubletransitions->GetTransitions()->end()) {
    msg_Error()<<"Error in "<<METHOD<<" for cluster\n"
	       <<(*cluster)
	       <<"   illegal flavour combination.\n"
	       <<"   Will return 0 and hope for the best.\n";
    return 0.;
  }
  double MC(cluster->Mass()),MC2(MC*MC);
  if (MC<p_doubletransitions->GetLightestMass(flpair)) {
    if (flpair.first.IsDiQuark() && flpair.second.IsDiQuark()) {
      // cluster consists of two diquarks, but is too light for regular decay
      // must therefore annihilate the beast.
      if (Annihilation(cluster,had1,had2)) {
	//msg_Out()<<"   --> "<<METHOD<<" cluster too light, annihilate.\n";
	return 1.;
      }
      else {
	msg_Error()<<"ERROR in "<<METHOD<<":\n"
		   <<"   Found cluster that MUST annihilate, but couldn't.\n"
		   <<"   Will return -1 and hope for the best.\n";
	abort();
	return -1;
      }
    }
    else {
      //msg_Out()<<"   --> "<<METHOD<<" cluster too light, return -1.\n";
      had1 = had2 = Flavour(kf_none);
      return -1.;
    }
  }
  double critM(p_doubletransitions->GetLightestMass(flpair)*(1.-m_decayoffset)+
	       p_doubletransitions->GetHeaviestMass(flpair)*m_decayoffset);
  if (MC>critM &&
      MC>cluster->GetTrip()->m_mass+cluster->GetAnti()->m_mass+m_minmass) {
    had1 = had2 = Flavour(kf_none);
    return 0.;
  }
  double totweight(0.),m1,m2,wt(1.),wfweight(0.),wfmax(0.);
  double tm(cluster->GetTrip()->m_flav.HadMass());
  double am(cluster->GetTrip()->m_flav.HadMass());
  Flavour max1, max2;
  for (Double_Transition_Siter decit=dtliter->second->begin();
       decit!=dtliter->second->end();decit++) {
    m1  = decit->first.first.HadMass();
    m2  = decit->first.second.HadMass();
    if (m1+m2<MC) {
      wt = (sqrt((MC2-sqr(m1+m2))*(MC2-sqr(m1-m2))) * 
	    pow(sqr(m1+m2)/MC2,m_chi) * decit->second);
      if (wfweight>wfmax) {
	max1  = decit->first.first;
	max2  = decit->first.second;
	wfmax = wfweight;
      }
      totweight += wt;
    }
  }
  if (totweight<=0.) return 0.;
  had1 = had2 = Flavour(kf_none); 
  double disc(totweight * 0.9999999999*ran->Get());
  for (Double_Transition_Siter decit=dtliter->second->begin();
       decit!=dtliter->second->end();decit++) {
    m1    = decit->first.first.HadMass();
    m2    = decit->first.second.HadMass();
    if (m1+m2<MC) {
      wt = (sqrt((MC2-sqr(m1+m2))*(MC2-sqr(m1-m2))) * 
	    pow(sqr(m1+m2)/MC2,m_chi) * decit->second);
      disc -= wt;
      if (disc<0.) {
	had1 = decit->first.first;
	had2 = decit->first.second;
	break;
      }
    }
  }
  return totweight/(16.*M_PI*MC*MC*MC);
}


bool Soft_Cluster_Handler::
Annihilation(Cluster * cluster,Flavour & had1,Flavour & had2) {
  int kfc1(int(cluster->GetTrip()->m_flav.Kfcode())); 
  int kfc2(int(cluster->GetAnti()->m_flav.Kfcode())); 
  kf_code kfc11(kfc1/1000),kfc12((kfc1-kfc11*1000)/100);
  kf_code kfc21(kfc2/1000),kfc22((kfc2-kfc21*1000)/100);
  Flavour fl1(kfc11), fl2(kfc12), fl3(kfc21), fl4(kfc22);
  //msg_Out()<<METHOD<<"("<<fl1<<", "<<fl2<<") & ("<<fl3<<", "<<fl4<<").\n";
  fl1 = fl1.Bar();
  fl2 = fl2.Bar();
  bool order(ran->Get()>0.5?true:false);
  Proto_Particle *pp1(new Proto_Particle(fl1,cluster->GetTrip()->m_mom/2.,'l'));
  Proto_Particle *pp2(new Proto_Particle(fl2,cluster->GetAnti()->m_mom/2.,'l'));
  Proto_Particle *pp3(new Proto_Particle(fl3,(order?
					      cluster->GetTrip()->m_mom/2.:
					      cluster->GetAnti()->m_mom/2.),
					 'l'));
  Proto_Particle *pp4(new Proto_Particle(fl4,(order?
					      cluster->GetAnti()->m_mom/2.:
					      cluster->GetTrip()->m_mom/2.),
					 'l'));
  Cluster cluster1((order?pp3:pp4),pp1), cluster2((order?pp4:pp3),pp2);
  double mass(cluster->Mass());
  //msg_Out()<<cluster1<<cluster2;
  double wt1(TransformWeight(&cluster1,had1,true));
  double wt2(TransformWeight(&cluster2,had2,true));
  //msg_Out()<<"  --> "<<had1<<" + "<<had2<<".\n";
  if (had1.Mass()+had2.Mass()>mass) return false;
  return true;
}


/// Safety mechanism for Cluster dacay into two hadrons
void Soft_Cluster_Handler::FixHHDecay(Cluster * cluster,Blob * blob,
				      const Flavour had1,const Flavour had2,
				      const bool & constrained)
{
  double M   = cluster->Mass(), M2 = M*M;
  // Target hadron masses
  double m1  = sqr(had1.HadMass()), m2 = sqr(had2.HadMass());

  cluster->BoostInCMSAndRotateOnZ();

  double E1((M2+m1-m2)/(2.*M)), pl2(sqr(E1)-m1); ///FIXME comment longitud momentum**2
  double masscor(m_pt02/Max(m_pt02,m1) * m_pt02/Max(m_pt02,m2)); /// FIXME comment mass correction
  /// The higher the hadron mass, the smaller the masscor,
  // -> the heavier the hadrons, the more likely they are blabla original direction
  // -> makes b-frag sufficiently hard TODO check do b-hadrons follow the initial quark?
  //
  double stheta, pt2; /// FIXME: initialize to 0, pl2 ?
  /// Black magic
  do { 
    stheta = 1.-2.*ran->Get(); // TODO chage/test cos
    pt2    = pl2*sqr(stheta);
  } while (pt2>m_pt2max*m_pt2maxfac*masscor || 
	   sqr((*p_as)(pt2,false)/p_as->MaxValue())<ran->Get()); // get rid of 2nd cond?
  double pt     = sqrt(pt2);
  int sign      = cluster->GetTrip()->m_mom[3]<0?-1:1;
  double pl1    = sign*sqrt(sqr(E1)-sqr(pt)-m1);
  double cosphi = cos(2.*M_PI*ran->Get()), sinphi = sqrt(1.-cosphi*cosphi);
  Vec4D  p1     = Vec4D(E1,pt*cosphi,pt*sinphi,pl1);
  Vec4D  p2     = cluster->Momentum()-p1;

  if (p1[0]<0. || p2[0]<0.) {
      throw Return_Value::Retry_Event;
  }

  // if (cluster->GetTrip()->m_flav==Flavour(kf_b) ||
  //     cluster->GetAnti()->m_flav==Flavour(kf_b).Bar()) {
  //   msg_Out()<<"\n\n\n"
  //            <<"======================================================\n"
  // 	     <<METHOD<<" for M = "<<M<<", sign = "<<sign<<", pt = "<<pt<<":\n"
  // 	     <<(*cluster)
  // 	     <<"   "<<cluster->Momentum()<<" --> \n"
  // 	     <<"   "<<p1<<" ("<<had1<<") + "<<p2<<" ("<<had2<<")\n."
  // 	     <<"======================================================\n";
  // }
  cluster->RotateAndBoostBack(p1);
  cluster->RotateAndBoostBack(p2);
  cluster->RotateAndBoostBack();

  Particle * left(new Particle(-1,had1,p1));
  left->SetNumber();
  left->SetInfo('P');
  left->SetFinalMass(had1.HadMass());
  Particle * right(new Particle(-1,had2,p2));
  right->SetNumber();
  right->SetInfo('P');
  right->SetFinalMass(had2.HadMass());
  control::s_AHAparticles+=2; // For test if all is deleted in the end



  if (blob!=NULL) {
    blob->AddToOutParticles(left);
    blob->AddToOutParticles(right);
  }
  //if (cluster->GetTrip()->m_info=='B' || cluster->GetAnti()->m_info=='B') {
  //  msg_Out()<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
  // 	     <<"Cluster decay (pt + "<<pt<<") for cluster \n"<<(*cluster)<<"\n"
  // 	     <<"++> "<<left->Momentum()<<" + "<<right->Momentum()<<" for "
  // 	     <<left->Flav()<<" + "<<right->Flav()<<"\n"
  // 	     <<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  //}
  if (m_ana) {
    m_histograms["PT_HH"]->Insert(pt);
    m_histograms["PT2_HH"]->Insert(pt*pt);
  }
}

Vec4D Soft_Cluster_Handler::SumMomentum(Cluster_List * clin) {
  Cluster_Iterator cit;
  Vec4D listmom(0.,0.,0.,0.);
  for (cit=clin->begin();cit!=clin->end();cit++) listmom += (*cit)->Momentum(); 
  return listmom;
}

