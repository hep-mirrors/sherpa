#include "Gluon_Decayer.H"
#include "Hadronisation_Parameters.H"
#include "Message.H"

using namespace AHADIC;
using namespace ATOOLS;


Gluon_Decayer::Gluon_Decayer(Dipole_Splitter * splitter) :
  p_splitter(splitter)
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
      decspec->massmin   = 2.*hadpars.GetConstituents()->Mass(fdit->first);
      m_options.insert(std::make_pair(fdit->first,decspec));
      msg_Tracking()<<"Insert option : g->"<<fdit->first<<" "<<fdit->first.Bar()<<std::endl;
    }
  }
  p_splitter->SetOptions(&m_options);

  msg_Tracking()<<"------------- END OF GLUON_DECAYER --------------"<<std::endl;
  if (m_options.empty()) {
    msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
	       <<"   No decay channels found for gluons, will abort."<<std::endl;
    abort();
  }
}

Gluon_Decayer::~Gluon_Decayer() {
  for (FDIter fdit=m_options.begin();fdit!=m_options.end();++fdit) {
    if (fdit->second!=NULL) { delete fdit->second; fdit->second=NULL; }
  }
  m_options.clear();
}

bool Gluon_Decayer::DecayList(Proto_Particle_List * plin)
{
  if (plin==NULL || plin->empty()) return true;

  if (msg->LevelIsDebugging()) {
    msg_Out()<<"------------------------------------------------------------"<<std::endl
	     <<"   "<<METHOD<<" : incoming particle list."<<std::endl<<(*plin)<<std::endl;
  }

#ifdef AHAmomcheck
  Vec4D checkbef(0.,0.,0.,0.), checkaft(0.,0.,0.,0.);
  for (PPL_Iterator pit=plin->begin();pit!=plin->end();pit++) checkbef += (*pit)->m_mom;
#endif

  FillDipoleList(plin);
  DecayDipoles();
  UpdatePPList(plin);

#ifdef AHAmomcheck
  for (PPL_Iterator pit=plin->begin();pit!=plin->end();pit++) checkaft += (*pit)->m_mom;
  if (dabs((checkbef-checkaft).Abs2())>1.e-12) {
    msg_Out()<<METHOD<<" yields momentum violation : "<<std::endl
  	     <<checkbef<<" - "<<checkaft<<" --> "<<(checkbef-checkaft).Abs2()<<std::endl;
  }
  else msg_Out()<<METHOD<<" conserves momentum."<<std::endl;
#endif

  if (msg->LevelIsDebugging()) {
    msg_Out()<<"   "<<METHOD<<" : outgoing particle list."<<std::endl<<(*plin)
	     <<"------------------------------------------------------------"<<std::endl<<std::endl;
  }

  return true;
}

bool Gluon_Decayer::FillDipoleList(Proto_Particle_List * plin)
{
  Vec4D intot(0.,0.,0.,0.);
  PPL_Iterator pit(plin->begin()), pit1(plin->begin());
  pit1++;
  Proto_Particle * begin(*pit);
  Dipole * dip;
  do {
    intot += (*pit)->m_mom;
    dip = new Dipole(*pit,*pit1);
    m_dipoles.push_back(dip);
    pit = pit1;
    pit1++;
  } while (pit1!=plin->end());
  intot += (*pit)->m_mom;
  if ((*pit)->m_flav.IsGluon()) {
    if (begin->m_flav.IsGluon()) {
      dip = new Dipole(*pit,begin);
      m_dipoles.push_back(dip);
    }
    else {
      msg_Error()<<"ERROR in "<<METHOD<<":"<<std::endl
		 <<"    Last flavour in list = "<<(*pit)->m_flav
		 <<" but first flavour = "<<begin->m_flav<<"."<<std::endl
		 <<"   Don't know what to do, abort."<<std::endl;
      abort();
    }
  }
  return true;
}

void Gluon_Decayer::UpdatePPList(Proto_Particle_List * plin)
{
  if (!plin || plin->empty()) return;
  plin->clear();
  for (DipIter dip=m_dipoles.begin();dip!=m_dipoles.end();dip++) {
    plin->push_back((*dip)->Triplet());
    plin->push_back((*dip)->AntiTriplet());
    delete (*dip);
  }
  m_dipoles.clear();
}

bool Gluon_Decayer::DecayDipoles() {
  DipIter dip;
  do { 
    dip = SelectDipole(); 
    if (dip==m_dipoles.end()) return true;
#ifdef AHAmomcheck
    msg_Out()<<"~~~~~~~~~~~~~~ "<<METHOD<<"("<<m_dipoles.size()<<") ~~~~~~~~~~~~~~~~~~"<<std::endl;
    Vec4D checkbef(0.,0.,0.,0.);
    for (DipIter dipiter=m_dipoles.begin();dipiter!=m_dipoles.end();dipiter++) {
      if ((*dipiter)->AntiTriplet()->m_flav!=Flavour(kf::gluon)) checkbef += (*dipiter)->Momentum();
      else checkbef += (*dipiter)->Triplet()->m_mom;
    }
    //for (DipIter dipiter=m_dipoles.begin();dipiter!=m_dipoles.end();dipiter++)
    //  (*dipiter)->Output();
#endif
    //std::cout<<METHOD<<" splits the following dipole:"<<std::endl;
    //(*dip)->Output();
    if (!p_splitter->SplitDipole((*dip))) {
      if (!Rescue(dip)) { 
	//	msg_Out()<<"............... Rescue failed ..................."<<std::endl;
	dip=m_dipoles.begin(); continue; 
      }
    }
#ifdef AHAmomcheck
    SplitIt((*dip),checkbef);
#else
    SplitIt((*dip));
#endif
  } while (dip!=m_dipoles.end());
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

bool Gluon_Decayer::Rescue(DipIter & dip) {
  DipIter partner(dip);
  if ((*dip)->Triplet()->m_flav.IsGluon() &&
      (*dip)->AntiTriplet()->m_flav.IsGluon()) {
    partner++;
    if (p_splitter->SplitDipole((*partner))) {
      dip = partner;
      return true;
    }
    else {
      partner = dip;
      partner--;
      if (p_splitter->SplitDipole((*partner))) {
	dip = partner;
	return true;
      }
    }
    if (MergeDipoles(partner,dip)) return false;
  }
  else if (!(*dip)->Triplet()->m_flav.IsGluon() &&
	   (*dip)->AntiTriplet()->m_flav.IsGluon()) {
    partner++;
    if (p_splitter->SplitDipole((*partner))) {
      dip = partner;
      return true;
    }
    if (MergeDipoles(dip,partner)) return false;
  }
  else if (!(*dip)->AntiTriplet()->m_flav.IsGluon() &&
	   (*dip)->Triplet()->m_flav.IsGluon()) {
    partner--;
    if (p_splitter->SplitDipole((*partner))) {
      dip = partner;
      return true;
    }
  }  
  if (MergeDipoles(partner,dip)) return false;
}

bool Gluon_Decayer::MergeDipoles(DipIter & dip1,DipIter & dip2) {
  Vec4D   Q(0.,0.,0.,0.),pi,pj,pk;
  if (msg->LevelIsDebugging()) {
    msg_Out()<<METHOD<<" : "<<std::endl
	     <<sqrt((*dip1)->Mass2())<<" vs. "<<sqrt((*dip2)->Mass2())<<std::endl;
  }
  Q += pi = (*dip1)->Triplet()->m_mom;
  Q += pj = (*dip2)->Triplet()->m_mom;
  Q += pk = (*dip2)->AntiTriplet()->m_mom;
  double Q2   = Q.Abs2();
  double mij2 = sqr(hadpars.GetConstituents()->Mass((*dip1)->Triplet()->m_flav));
  double mk2  = sqr(hadpars.GetConstituents()->Mass((*dip2)->AntiTriplet()->m_flav));
  double pij2 = (pi+pj).Abs2();
  double aij  = (sqr(Q2-mij2-mk2)-4.*mij2*mk2), bij = (sqr(Q2-pij2-mk2)-4.*pij2*mk2);
  Vec4D  pkt  = sqrt(aij/bij) * (pk - (Q*pk)/Q2*Q) + (Q2 + mk2-mij2)/(2.*Q2)*Q;  
  Vec4D  pijt = Q-pkt;
  if (msg->LevelIsDebugging()) {
    msg_Out()<<"    merge "<<pi<<", "<<pi.Abs2()<<" + "<<pj<<", "<<pj.Abs2()<<std::endl
	     <<"        + "<<pk<<", "<<pk.Abs2()<<" = "<<Q<<", "<<Q.Abs2()<<std::endl
	     <<"     into "<<pijt<<" + "<<pkt<<" from "<<aij<<"/"<<bij
	     <<" and mij = "<<sqrt(mij2)<<" and mk = "<<sqrt(mk2)<<std::endl;
  }
  (*dip1)->Triplet()->m_mom      = pijt;
  (*dip2)->AntiTriplet()->m_mom  = pkt;

  if (msg->LevelIsDebugging()) {
    msg_Out()<<"      check this (before): "
	     <<(*dip1)->Triplet()<<"/"<<(*dip1)->AntiTriplet()<<" and "
	     <<(*dip2)->Triplet()<<"/"<<(*dip2)->AntiTriplet()<<"."<<std::endl;
  }

  (*dip1)->SetAntiTriplet((*dip2)->AntiTriplet());
  delete (*dip2)->Triplet();
  delete (*dip2);

  m_dipoles.erase(dip2);
  for (DipIter dipiter=m_dipoles.begin();dipiter!=m_dipoles.end();dipiter++) {
    (*dipiter)->Update();
  }

  if (msg->LevelIsDebugging()) {
    (*dip1)->Output();
  }
#ifdef AHAmomcheck
  Vec4D Qafter = (*dip1)->Momentum();
  if (dabs((Q-Qafter).Abs2())>1.e-12) {
    msg_Out()<<METHOD<<" yields momentum violation : "<<std::endl
  	     <<Q<<" - "<<Qafter<<" --> "<<(Q-Qafter).Abs2()<<std::endl;    
  }
  else msg_Out()<<METHOD<<" conserves momentum."<<std::endl;
#endif
  return true;
}

void Gluon_Decayer::SplitIt(Dipole * dip,Vec4D checkbef) {
  Proto_Particle * new1, * new2;
  p_splitter->GetNewParticles(new1,new2);
  DipIter neighbour;
#ifdef AHAmomcheck
  Vec4D checkaft(0.,0.,0.,0.);
#endif
  if ((*m_dipoles.begin())->Triplet()->m_flav.IsGluon()) {
    //std::cout<<METHOD<<" g-g-g-g "<<dip->IsSwitched()<<std::endl;
    while ((*m_dipoles.begin())!=dip) {
      m_dipoles.push_back((*m_dipoles.begin()));
      m_dipoles.pop_front();
    }
    if (dip->IsSwitched()) {
      dip->SetTriplet(new2);
      neighbour = m_dipoles.end();
      neighbour--;
      (*neighbour)->SetAntiTriplet(new1);
    }
    else {
      m_dipoles.push_back((*m_dipoles.begin()));
      m_dipoles.pop_front();
      dip->SetAntiTriplet(new1);
      neighbour = m_dipoles.begin();
      (*neighbour)->SetTriplet(new2);
    }
  }
  else {
    //std::cout<<METHOD<<" q-g-g-q "<<dip->IsSwitched()<<std::endl;
    DipIter test;
    for (test=m_dipoles.begin();test!=m_dipoles.end();test++) {
      if ((*test)==dip) break;
    }
    //std::cout<<"    check this: "<<(*test)<<"/"<<dip<<" test = begin "
    //	     <<(test==m_dipoles.begin())<<std::endl;
    neighbour = test;
    if (dip->IsSwitched()) {
      dip->SetTriplet(new2);
      neighbour--;
      (*neighbour)->SetAntiTriplet(new1);
    }
    else {
      dip->SetAntiTriplet(new1);
      neighbour++;
      (*neighbour)->SetTriplet(new2);      
    }
  }  

  for (DipIter dipiter=m_dipoles.begin();dipiter!=m_dipoles.end();dipiter++) {
    (*dipiter)->Update();
#ifdef AHAmomcheck
    if ((*dipiter)->AntiTriplet()->m_flav!=Flavour(kf::gluon)) checkaft += (*dipiter)->Momentum();
    else checkaft += (*dipiter)->Triplet()->m_mom;
#endif
  }
#ifdef AHAmomcheck
  if (dabs((checkbef-checkaft).Abs2())>1.e-12) {
    msg_Out()<<METHOD<<" yields momentum violation : "<<std::endl
  	     <<checkbef<<" - "<<checkaft<<" --> "<<(checkbef-checkaft).Abs2()<<std::endl;
    for (DipIter dipiter=m_dipoles.begin();dipiter!=m_dipoles.end();dipiter++)
      (*dipiter)->Output();
  }
  else msg_Out()<<METHOD<<" conserves momentum."<<std::endl;
#endif
}

void Gluon_Decayer::PrintDipoleList()
{
  for (DipIter dip=m_dipoles.begin();dip!=m_dipoles.end();dip++) {
    msg_Out()<<"Dipole("<<((*dip)->MustDecay())<<", "<<sqrt((*dip)->Mass2())<<") : "<<std::endl
	     <<"  "<<(*dip)->Triplet()->m_flav<<"("<<(*dip)->Triplet()->m_mom<<"), "
	     <<" "<<hadpars.GetConstituents()->Mass((*dip)->Triplet()->m_flav)
	     <<" vs. "<<sqrt(Max((*dip)->Triplet()->m_mom.Abs2(),0.0))<<";"<<std::endl
	     <<"  "<<(*dip)->AntiTriplet()->m_flav<<"("<<(*dip)->AntiTriplet()->m_mom<<"),"
	     <<" "<<hadpars.GetConstituents()->Mass((*dip)->AntiTriplet()->m_flav)
	     <<" vs. "<<sqrt(Max((*dip)->AntiTriplet()->m_mom.Abs2(),0.0))<<"."<<std::endl;
  }
}


