#include "AHADIC++/Formation/Gluon_Decayer.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"

using namespace AHADIC;
using namespace ATOOLS;


Gluon_Decayer::Gluon_Decayer(Dipole_Splitter * splitter,bool ana) :
  p_splitter(splitter), 
  m_pt2max(sqr(hadpars.Get(std::string("ptmax")))), 
  m_analyse(ana),m_tot(0),m_s(0),m_u(0),m_d(0)
{ 
  double norm(0.);
  for (FlavCCMap_Iterator fdit=hadpars.GetConstituents()->CCMap.begin();
       fdit!=hadpars.GetConstituents()->CCMap.end();fdit++) {
    if (hadpars.GetConstituents()->TotWeight(fdit->first)>norm)
      norm = hadpars.GetConstituents()->TotWeight(fdit->first);
  }
  DecaySpecs * decspec;
  for (FlavCCMap_Iterator fdit=hadpars.GetConstituents()->CCMap.begin();
       fdit!=hadpars.GetConstituents()->CCMap.end();fdit++) {
    if (!fdit->first.IsAnti()) {
      decspec = new DecaySpecs;
      decspec->popweight = hadpars.GetConstituents()->TotWeight(fdit->first)/norm;
      decspec->massmin   = hadpars.GetConstituents()->Mass(fdit->first);
      m_options.insert(std::make_pair(fdit->first,decspec));
      msg_Debugging()<<"Insert option : g->"<<fdit->first<<" "<<fdit->first.Bar()<<std::endl;
    }
  }
  if (m_analyse) {
    m_histograms[std::string("PT_Gluon")]          = new Histogram(0,0.,2.,100);
    m_histograms[std::string("PT_Rescue")]         = new Histogram(0,0.,2.,100);
    m_histograms[std::string("Flavour_Gluon")]     = new Histogram(0,0.,15.,15);
    m_histograms[std::string("Flavour_Rescue")]    = new Histogram(0,0.,15.,15);
    m_histograms[std::string("DecayedDipoleMass")] = new Histogram(0,0.,15.,30);
    m_histograms[std::string("MergedMassBefore")]  = new Histogram(0,0.,15.,30);
    m_histograms[std::string("MergedMassAfter")]   = new Histogram(0,0.,30.,60);
  }

  p_splitter->SetOptions(&m_options);

  msg_Tracking()<<"------------- END OF GLUON_DECAYER --------------"<<std::endl;
  if (m_options.empty()) {
    msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
	       <<"   No decay channels found for gluons, will abort the run."<<std::endl
	       <<"   Please contact the Sherpa group for further assistance."<<std::endl;
    exit(0);
  }
}

Gluon_Decayer::~Gluon_Decayer() {
  if (m_analyse) {
    Histogram * histo;
    std::string name;
    for (std::map<std::string,Histogram *>::iterator hit=m_histograms.begin();
	 hit!=m_histograms.end();hit++) {
      histo = hit->second;
      name  = std::string("Fragmentation_Analysis/")+hit->first+std::string(".dat");
      histo->Output(name);
      delete histo;
    }
    m_histograms.clear();
  }
  //msg_Out()<<"Flavour production in "<<m_tot<<" gluon splittings: "
  //	   <<"s_rate = "<<double(m_s)/double(m_tot)<<", "
  //	   <<"u_rate = "<<double(m_u)/double(m_tot)<<", "
  //	   <<"d_rate = "<<double(m_d)/double(m_tot)<<"."<<std::endl;

  for (FDIter fdit=m_options.begin();fdit!=m_options.end();++fdit) {
    if (fdit->second!=NULL) { delete fdit->second; fdit->second=NULL; }
  }
  m_options.clear();
}

bool Gluon_Decayer::DecayList(Proto_Particle_List * plin)
{
  if (plin==NULL || plin->empty()) return true;
  if (!m_dipoles.empty()) {
    while (!m_dipoles.empty()) {
      m_dipoles.pop_front();
    }
  }

  msg_Tracking()<<"------------------------------------------------------------"<<std::endl
		<<"------------------------------------------------------------"<<std::endl
		<<"   "<<METHOD<<" : incoming particle list."<<std::endl<<(*plin);

#ifdef AHAmomcheck
  Vec4D checkbef(0.,0.,0.,0.), checkaft(0.,0.,0.,0.);
  for (PPL_Iterator pit=plin->begin();pit!=plin->end();pit++) checkbef += (*pit)->m_mom;
#endif

  if (!FillDipoleList(plin)) return false;
  if (!DecayDipoles()) return false;
  UpdatePPList(plin);

#ifdef AHAmomcheck
  for (PPL_Iterator pit=plin->begin();pit!=plin->end();pit++) checkaft += (*pit)->m_mom;
  if (dabs((checkbef-checkaft).Abs2()/checkbef.Abs2())>1.e-6) {
    msg_Error()<<METHOD<<" yields momentum violation : "<<std::endl
	       <<checkbef<<" - "<<checkaft<<" --> "<<(checkbef-checkaft).Abs2()<<std::endl;
  }
  else msg_Tracking()<<METHOD<<" conserves momentum."<<std::endl;
#endif

  msg_Tracking()<<"   "<<METHOD<<" : outgoing particle list."<<std::endl<<(*plin);
  if (true) {
    for (PPL_Iterator pit=plin->begin();pit!=plin->end();pit++) {
      if (((*pit)->m_mom.Abs2()-sqr((*pit)->m_flav.HadMass()))/
	  ((*pit)->m_mom.Abs2()+sqr((*pit)->m_flav.HadMass())) > 1.e-6) { 
	msg_Tracking()<<"   Potential problem: "<<(**pit)
		      <<"   mass = "<<(*pit)->m_flav.HadMass()<<" vs. "
		      <<sqrt((*pit)->m_mom.Abs2())<<"."<<std::endl;
      }
    }
  }
  msg_Tracking()<<"------------------------------------------------------------"
		<<"------------------------------------------------------------"
		<<std::endl<<std::endl;
  
  return true;
}

bool Gluon_Decayer::FillDipoleList(Proto_Particle_List * plin)
{
  PPL_Iterator pit(plin->begin()), pit1(plin->begin());
  pit1++;
  Proto_Particle * begin(*pit);
  Dipole * dip;
  do {
    (*pit)->m_kt2max  = (*pit)->m_info=='B'?(*pit)->m_mom.PPerp2():m_pt2max;
    (*pit1)->m_kt2max = (*pit1)->m_info=='B'?(*pit1)->m_mom.PPerp2():m_pt2max;
    dip = new Dipole(*pit,*pit1);
    m_dipoles.push_back(dip);
    pit = pit1;
    pit1++;
    PrintDipoleList();
  } while (pit1!=plin->end());
  if ((*pit)->m_flav.IsGluon()) {
    if (begin->m_flav.IsGluon()) {
      dip = new Dipole(*pit,begin);
      m_dipoles.push_back(dip);
    }
    else {
      msg_Error()<<"ERROR in "<<METHOD<<":"<<std::endl
		 <<"    Last flavour in list = "<<(*pit)->m_flav
		 <<" but first flavour = "<<begin->m_flav<<"."<<std::endl
		 <<"   Don't know what to do, try new event."<<std::endl;
      return false;
    }
  }
  
  msg_Debugging()<<"------------------------------------------------------------"<<std::endl
		 <<"--------------------------------------------------"<<std::endl
		 <<"   "<<METHOD<<": yields dipole list: "<<std::endl;
  PrintDipoleList();

  return true;
}

void Gluon_Decayer::UpdatePPList(Proto_Particle_List * plin)
{
  if (plin==NULL || plin->empty()) return;
  plin->clear();
  for (DipIter dip=m_dipoles.begin();dip!=m_dipoles.end();dip++) {
    plin->push_back((*dip)->Triplet());
    plin->push_back((*dip)->AntiTriplet());
    delete (*dip);
  }
  m_dipoles.clear();
}

bool Gluon_Decayer::DecayDipoles() {
  msg_Debugging()<<"------------------------------------------------------------"<<std::endl
		 <<"------------------------------------------------------------"<<std::endl
		 <<"   "<<METHOD<<" for "<<m_dipoles.size()<<" dipoles:"<<std::endl;
  if (msg->LevelIsDebugging()) {
    for (DipIter diter=m_dipoles.begin();diter!=m_dipoles.end();diter++) {
      msg_Out()<<(*diter)<<" : ";
      (*diter)->Output();
    }
  }
  DipIter dipiter;
  do {
    dipiter = SelectDipole(); 
    if (dipiter==m_dipoles.end()) {
      msg_Debugging()<<METHOD<<" : all dipoles done!"<<std::endl;
      return true;
    }
#ifdef AHAmomcheck
    Vec4D checkbef(0.,0.,0.,0.);
    for (DipIter diter=m_dipoles.begin();diter!=m_dipoles.end();diter++) {
      if ((*diter)->AntiTriplet()->m_flav!=Flavour(kf_gluon)) 
	checkbef += (*diter)->Momentum();
      else checkbef += (*diter)->Triplet()->m_mom;
    }
#endif
    if (msg->LevelIsDebugging()) {
      msg_Out()<<"--- "<<METHOD<<" splits the following dipole :"<<(*dipiter)<<std::endl;
      (*dipiter)->Output();
    }
    if (!p_splitter->SplitDipole((*dipiter),true,false)) {
      switch (Rescue(dipiter)) {
      case -1:
	msg_Debugging()<<"--- Rescue failed."<<std::endl;
	return false;
      case 0:  
	dipiter=m_dipoles.begin(); continue; 
	break;
      case 1:
	if (m_analyse) {
	  Histogram* histo((m_histograms.find(std::string("PT_Rescue")))->second);
	  histo->Insert(sqrt(p_splitter->PT2()));
	}
      default:
	break;
      }
    }
    else {
      if (m_analyse) {
	Histogram* histo((m_histograms.find(std::string("PT_Gluon")))->second);
	histo->Insert(sqrt(p_splitter->PT2()));
	Histogram* histo2((m_histograms.find(std::string("DecayedDipoleMass")))->second);
	histo2->Insert(sqrt(p_splitter->Mass2()));
      }
      AfterSplit(dipiter);
    }
#ifdef AHAmomcheck
    SplitIt(dipiter,checkbef);
#else
    SplitIt(dipiter);
#endif
    if (msg->LevelIsDebugging()) {
      msg_Out()<<"--- "<<METHOD<<" (after splitting) for "<<m_dipoles.size()<<" dipoles:"<<std::endl;
      for (DipIter diter=m_dipoles.begin();diter!=m_dipoles.end();diter++) (*diter)->Output();
    }
  } while (dipiter!=m_dipoles.end());

  return true;
}

DipIter Gluon_Decayer::SelectDipole() {
  double smax(-1.),smin(1.e20);
  DipIter winner=m_dipoles.end();
  for (DipIter dip=m_dipoles.begin();dip!=m_dipoles.end();dip++) {
    if ((*dip)->MustDecay() && (*dip)->MassBar2()>smax) {
      smax   = (*dip)->MassBar2();
    }
    if ((*dip)->MustDecay() && (*dip)->MassBar2()<smin) {
      smin   = (*dip)->MassBar2();
      winner = dip;
    }
  }
  return winner;
}

int Gluon_Decayer::Rescue(DipIter & dip) {
  if (msg->LevelIsDebugging()) {
    msg_Out()<<"--- "<<METHOD<<" ---"<<std::endl;
    (*dip)->Output();
  }
  DipIter partner=dip, dip1,dip2;
  if ((*dip)->Triplet()->m_flav.IsGluon() &&
      (*dip)->AntiTriplet()->m_flav.IsGluon()) {
    if ((*dip)!=(*m_dipoles.rbegin())) {
      partner++;
      if (p_splitter->SplitDipole((*partner),true,false)) {
	AfterSplit(partner);
	dip = partner;
	return 1;
      }
    }
    if (dip!=m_dipoles.begin()) {
      partner = dip;
      partner--;
      if (p_splitter->SplitDipole((*partner),true,false)) {
	AfterSplit(partner);
	dip = partner;
	return 1;
      }
    }
    if ((*dip)==(*m_dipoles.rbegin())) { 
      dip1 = dip2 = dip; dip1--; 
    }
    else if (dip==m_dipoles.begin()) { 
      dip1 = dip2 = dip; dip2++; 
    }
    else { 
      dip1 = dip2 = dip; 
      if (ran.Get()>0.5) dip2++;
      else dip1--;
    }
  }
  else if (!(*dip)->Triplet()->m_flav.IsGluon() &&
	   (*dip)->AntiTriplet()->m_flav.IsGluon()) {
    partner++;
    if (p_splitter->SplitDipole((*partner),true,false)) {
      AfterSplit(partner);
      dip = partner;
      return 1;
    }
    dip1 = dip; dip2 = partner;
  }
  else if (!(*dip)->AntiTriplet()->m_flav.IsGluon() &&
	   (*dip)->Triplet()->m_flav.IsGluon()) {
    partner--;
    if (p_splitter->SplitDipole((*partner),true,false)) {
      AfterSplit(partner);
      dip = partner;
      return 1;
    }
    dip1 = partner; dip2 = dip;
  }
  if (m_dipoles.size()==2 &&
      (*dip1)->Triplet()->m_flav.IsGluon() &&
      (*dip1)->AntiTriplet()->m_flav.IsGluon() &&
      (*dip1)->Triplet() == (*dip2)->AntiTriplet() &&
      (*dip1)->AntiTriplet() == (*dip2)->Triplet()) {
    msg_Debugging()<<"--- Warning in "<<METHOD<<" : "<<std::endl
		   <<"   A gluon-gluon singlet with low mass : "<<sqrt((*dip1)->Mass2())<<","<<std::endl
		   <<"   Do not know how to handle this, trigger new event."<<std::endl;
    return -1;
  }
  if (!MergeDipoles(dip1,dip2)) return -1; 
  return 0;
}

bool Gluon_Decayer::MergeDipoles(DipIter & dip1,DipIter & dip2) {
  msg_Debugging()<<"--- "<<METHOD<<" for : "<<(*dip1)<<" + "<<(*dip2)<<std::endl
		 <<sqrt((*dip1)->Mass2())<<" vs. "<<sqrt((*dip2)->Mass2())<<std::endl;
  if (m_analyse) {
    Histogram* histo((m_histograms.find(std::string("MergedMassBefore")))->second);
    histo->Insert(sqrt((*dip1)->Mass2()));
    histo->Insert(sqrt((*dip2)->Mass2()));
  }
  Dipole save1(new Proto_Particle((*(*dip1)->Triplet())),
	       new Proto_Particle((*(*dip1)->AntiTriplet())));
  Dipole save2(new Proto_Particle((*(*dip2)->Triplet())),
	       new Proto_Particle((*(*dip2)->AntiTriplet())));
#ifdef memchecker
  std::cout<<"### New Proto_Particle ("
	   <<save1.Triplet()->m_flav<<"/"<<save1.Triplet()<<") from "<<METHOD<<"."<<std::endl;
  std::cout<<"### New Proto_Particle ("
	   <<save1.AntiTriplet()->m_flav<<"/"<<save1.AntiTriplet()<<") from "<<METHOD<<"."<<std::endl;
  std::cout<<"### New Proto_Particle ("
	   <<save2.Triplet()->m_flav<<"/"<<save2.Triplet()<<") from "<<METHOD<<"."<<std::endl;
  std::cout<<"### New Proto_Particle ("
	   <<save2.AntiTriplet()->m_flav<<"/"<<save2.AntiTriplet()<<") from "<<METHOD<<"."<<std::endl;
#endif


  Vec4D   Q(0.,0.,0.,0.),pi,pj,pk;
  Q += pi = (*dip1)->Triplet()->m_mom;
  Q += pj = (*dip2)->Triplet()->m_mom;
  Q += pk = (*dip2)->AntiTriplet()->m_mom;
  double Q2   = Q.Abs2();
  double pij2 = (pi+pj).Abs2(), pjk2 = (pj+pk).Abs2();
  double mij2 = sqr(hadpars.GetConstituents()->Mass((*dip1)->Triplet()->m_flav)), mi2(mij2);
  double mk2  = sqr(hadpars.GetConstituents()->Mass((*dip2)->AntiTriplet()->m_flav)), mjk2(mk2);
  double aij  = (sqr(Q2-mij2-mk2)-4.*mij2*mk2), bij = (sqr(Q2-pij2-mk2)-4.*pij2*mk2);
  double ajk  = (sqr(Q2-mjk2-mi2)-4.*mjk2*mi2), bjk = (sqr(Q2-pjk2-mi2)-4.*pjk2*mi2);
  if (aij/bij<0. && ajk/bjk<0.) {
    msg_Error()<<"Error in "<<METHOD<<"."<<std::endl
	       <<"   Cannot merge dipoles, kinematics does not work out."<<std::endl
	       <<"   Try new event."<<std::endl;
    return false;
  }
  if (ajk/bjk<0. || (pij2>pjk2 && aij/bij>0.)) {
    Vec4D  pkt  = sqrt(aij/bij) * (pk - (Q*pk)/Q2*Q) + (Q2 + mk2-mij2)/(2.*Q2)*Q;  
    Vec4D  pijt = Q-pkt;
    msg_Debugging()<<"--- "<<METHOD<<" merge :"<<std::endl
		   <<"   pi = "<<pi<<", "<<pi.Abs2()<<" + pj = "<<pj<<", "<<pj.Abs2()<<std::endl
		   <<"   + pk = "<<pk<<", "<<pk.Abs2()<<" --> Q = "<<Q<<", "<<Q.Abs2()<<std::endl
		   <<"     into "<<pijt<<"("<<pijt.Abs2()<<") + "<<pkt<<"("<<pkt.Abs2()<<") "
		   <<"from "<<aij<<"/"<<bij
		   <<" and mij = "<<sqrt(mij2)<<" and mk = "<<sqrt(mk2)<<std::endl;
    
    (*dip1)->Triplet()->m_mom      = pijt;
    (*dip2)->AntiTriplet()->m_mom  = pkt;
  }
  else {
    Vec4D  pit  = sqrt(ajk/bjk) * (pi - (Q*pi)/Q2*Q) + (Q2 + mi2-mjk2)/(2.*Q2)*Q;  
    //std::cout<<"Check 0 = "<<((pi - (Q*pi)/Q2*Q) * Q * (Q2 + mi2-mjk2)/(2.*Q2))<<"."<<std::endl
    //	     <<"    -bjk/(4Q^2) = "<<(pi - (Q*pi)/Q2*Q).Abs2()<<" = "<<(-bjk/(4.*Q2))<<"."<<std::endl
    //	     <<"    mi^2 = "<<pi.Abs2()<<" = "<<mi2<<", mk^2 = "<<pk.Abs2()<<" = "<<mk2<<"."<<std::endl;
    Vec4D  pjkt = Q-pit;
    msg_Debugging()<<"--- "<<METHOD<<" merge :"<<std::endl
		   <<"   pi = "<<pi<<", "<<pi.Abs2()<<" + pj = "<<pj<<", "<<pj.Abs2()<<std::endl
		   <<"   + pk = "<<pk<<", "<<pk.Abs2()<<" --> Q = "<<Q<<", "<<Q.Abs2()<<std::endl
		   <<"     into "<<pit<<"("<<pit.Abs2()<<") + "<<pjkt<<"("<<pjkt.Abs2()<<") "
		   <<"from "<<ajk<<"/"<<bjk
		   <<" and mi = "<<sqrt(mi2)<<" and mjk = "<<sqrt(mjk2)<<std::endl;
    
    
    (*dip1)->Triplet()->m_mom      = pit;
    (*dip2)->AntiTriplet()->m_mom  = pjkt;
  }
  (*dip1)->SetAntiTriplet((*dip2)->AntiTriplet());

  m_dipoles.erase(dip2);
  for (DipIter dipiter=m_dipoles.begin();dipiter!=m_dipoles.end();dipiter++) {
    (*dipiter)->Update();
  }

  if (msg->LevelIsDebugging()) (*dip1)->Output();
  if (m_analyse) {
    Histogram* histo((m_histograms.find(std::string("MergedMassAfter")))->second);
    histo->Insert(sqrt((*dip1)->Mass2()));
  }
#ifdef AHAmomcheck
  Vec4D Qafter = (*dip1)->Momentum();
  if (dabs((Q-Qafter).Abs2())>1.e-12) {
    msg_Error()<<METHOD<<" yields momentum violation : "<<std::endl
	       <<Q<<" - "<<Qafter<<" --> "<<(Q-Qafter).Abs2()<<std::endl;
    save1.Output();
    save2.Output();
  }
  else msg_Debugging()<<METHOD<<" conserves momentum."<<std::endl;
#endif
  delete save1.Triplet();
  delete save1.AntiTriplet();
  delete save2.Triplet();
  delete save2.AntiTriplet();
  return true;
}

void Gluon_Decayer::AfterSplit(DipIter dip) {
  DipIter partner(dip);
  if ((*dip)->IsSwitched()) {
    (*dip)->SetTriplet(NULL);
    if (dip!=m_dipoles.begin()) {
      partner--;
      (*partner)->SetAntiTriplet(NULL);
    }
    else (*m_dipoles.rbegin())->SetAntiTriplet(NULL);
    return;
  }
  else {
    (*dip)->SetAntiTriplet(NULL);
    if ((*dip)!=(*m_dipoles.rbegin())) partner++;
    else partner=m_dipoles.begin();
    (*partner)->SetTriplet(NULL);
    return;
  }
}

void Gluon_Decayer::SplitIt(DipIter dipiter,Vec4D checkbef) {
  //msg_Out()<<"--- "<<METHOD<<"(in): two new proto_particles:"<<std::endl;
  Proto_Particle * new1, * new2;
  p_splitter->GetNewParticles(new1,new2);
  //msg_Out()<<" new1: "<<(*new1)<<std::endl<<" new2: "<<(*new2)<<std::endl;

  if (new2->m_flav==Flavour(kf_d)) m_d++;
  else if (new2->m_flav==Flavour(kf_u)) m_u++;
  else if (new2->m_flav==Flavour(kf_s)) m_s++;
  m_tot++; 

  DipIter partner;
  Dipole * dip((*dipiter));
#ifdef AHAmomcheck
  Vec4D checkaft(0.,0.,0.,0.);
#endif
  if (m_dipoles.begin()==dipiter || (*m_dipoles.begin())->Triplet()==NULL ||
      (!(*m_dipoles.begin())->Triplet()->m_flav.IsQuark() &&
       !(*m_dipoles.begin())->Triplet()->m_flav.IsDiQuark())) {
    msg_Debugging()<<"--- "<<METHOD<<" g-g-g-g "<<dip->IsSwitched()<<std::endl;
    while (m_dipoles.begin()!=dipiter) {
      m_dipoles.push_back((*m_dipoles.begin()));
      m_dipoles.pop_front();
    }
    if (dip->IsSwitched()) {
      dip->SetTriplet(new2);
      partner = m_dipoles.end();
      partner--;
      (*partner)->SetAntiTriplet(new1);
    }
    else {
      m_dipoles.push_back((*m_dipoles.begin()));
      m_dipoles.pop_front();
      dip->SetAntiTriplet(new1);
      partner = m_dipoles.begin();
      (*partner)->SetTriplet(new2);
    }
  }
  else {
    partner = dipiter;
    if (msg->LevelIsDebugging()) {
      msg_Out()<<"--- "<<METHOD<<" q-g-g-q "<<dip->IsSwitched()<<std::endl
	       <<"--- "<<METHOD<<" check this: test = begin "
	       <<(partner==m_dipoles.begin())<<std::endl;
      (*partner)->Output();
      (*(dip)).Output();
    }

    if (dip->IsSwitched()) {
      partner--;
      if (msg->LevelIsDebugging()) {
	msg_Out()<<"--- "<<METHOD<<" check for neighbour after (switched): "<<std::endl;
	(*partner)->Output();
      }
      dip->SetTriplet(new2);
      (*partner)->SetAntiTriplet(new1);
    }
    else {
      partner++;
      if (msg->LevelIsDebugging()) {
	msg_Out()<<"--- "<<METHOD<<" check for neighbour after (original): "<<std::endl;
	(*partner)->Output();
      }
      dip->SetAntiTriplet(new1);
      (*partner)->SetTriplet(new2);      
    }
    if (msg->LevelIsDebugging()) {
      msg_Out()<<"--- "<<METHOD<<" check for neighbour after setting: "<<std::endl;
      (*partner)->Output();
      (*(dip)).Output();
    }
  }

  for (DipIter diter=m_dipoles.begin();diter!=m_dipoles.end();diter++) {
    (*diter)->Update();
    if (msg->LevelIsDebugging()) (*diter)->Output();
#ifdef AHAmomcheck
    if ((*diter)->AntiTriplet()->m_flav!=Flavour(kf_gluon)) checkaft += (*diter)->Momentum();
    else checkaft += (*diter)->Triplet()->m_mom;
#endif
  }
#ifdef AHAmomcheck
  if (dabs((checkbef-checkaft).Abs2())>1.e-12) {
    msg_Error()<<METHOD<<" yields momentum violation : "<<std::endl
	       <<checkbef<<" - "<<checkaft<<" --> "<<(checkbef-checkaft).Abs2()<<std::endl;
    for (DipIter diter=m_dipoles.begin();diter!=m_dipoles.end();diter++)
      (*diter)->Output();
  }
  else msg_Debugging()<<METHOD<<" conserves momentum."<<std::endl;
#endif
  msg_Debugging()<<"--- "<<METHOD<<"(out)."<<std::endl;
}

void Gluon_Decayer::PrintDipoleList()
{
  for (DipIter dip=m_dipoles.begin();dip!=m_dipoles.end();dip++) {
    msg_Debugging()<<"Dipole("<<((*dip)->MustDecay())<<", "
		   <<sqrt((*dip)->Mass2())<<") : "<<std::endl
		   <<"  "<<(*dip)->Triplet()->m_flav<<"("<<(*dip)->Triplet()->m_mom<<"), "
		   <<" "<<hadpars.GetConstituents()->Mass((*dip)->Triplet()->m_flav)
		   <<" vs. "<<sqrt(Max((*dip)->Triplet()->m_mom.Abs2(),0.0))
		   <<" --> kt^2_max = "<<(*dip)->Triplet()->m_kt2max<<";"<<std::endl
		   <<"  "<<(*dip)->AntiTriplet()->m_flav<<"("<<(*dip)->AntiTriplet()->m_mom<<"),"
		   <<" "<<hadpars.GetConstituents()->Mass((*dip)->AntiTriplet()->m_flav)
		   <<" vs. "<<sqrt(Max((*dip)->AntiTriplet()->m_mom.Abs2(),0.0))
		   <<" --> kt^2_max = "<<(*dip)->AntiTriplet()->m_kt2max<<"."<<std::endl;
  }
}


