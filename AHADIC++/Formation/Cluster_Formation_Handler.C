#include <cassert>
#include "AHADIC++/Formation/Cluster_Formation_Handler.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"

namespace AHADIC {
  bool triplet(Proto_Particle * pp) {
    return ((pp->m_flav.IsQuark() && !pp->m_flav.IsAnti()) ||
	    (pp->m_flav.IsDiQuark() && pp->m_flav.IsAnti()) );
  }
  
  bool antitriplet(Proto_Particle * pp) {
    return ((pp->m_flav.IsQuark() && pp->m_flav.IsAnti()) ||
	    (pp->m_flav.IsDiQuark() && !pp->m_flav.IsAnti()) );
  }
}

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Cluster_Formation_Handler::Cluster_Formation_Handler(Cluster_List* clulist,
						     bool ana) :
  m_single_cr(true), m_double_cr(false),
  p_gludecayer(new Gluon_Decayer(hadpars.GetSplitter(),ana)), 
  p_cformer(new Cluster_Former()),
  p_recons(new Colour_Reconnections(0,0,1.)), 
  p_softclusters(hadpars.GetSoftClusterHandler()),
  p_clulist(clulist), m_analyse(ana)
{ 
  if (m_analyse) {
    m_histograms[string("Cluster_Mass_Formation")]     = new Histogram(0,0.,100.,200);
    m_histograms[string("Cluster_Mass_Reconnections")] = new Histogram(0,0.,100.,200);
    m_histograms[string("Cluster_Mass_Transformed")]   = new Histogram(0,0.,100.,200);
    m_histograms[string("Cluster_Number_Formation")]   = new Histogram(0,0.,20.,20);
    m_histograms[string("Cluster_Number_Transformed")] = new Histogram(0,0.,20.,20);
    m_histograms[string("XB_bquark_before_shift")]     = new Histogram(0,0.,1.,100);
    m_histograms[string("XB_bquark_after_shift")]      = new Histogram(0,0.,1.,100);
    m_histograms[string("XB_bquark_in_cluster")]       = new Histogram(0,0.,1.,100);
  }
}


Cluster_Formation_Handler::~Cluster_Formation_Handler() {

  if (m_analyse) {
    Histogram * histo;
    string name;
    for (map<string,Histogram *>::iterator hit=m_histograms.begin();
	 hit!=m_histograms.end();hit++) {
      histo = hit->second;
      name  = string("Fragmentation_Analysis/")+hit->first+std::string(".dat");
      histo->Output(name);
      delete histo;
    }
    m_histograms.clear();
  }

  Reset();
  if (p_gludecayer) delete p_gludecayer;
  if (p_cformer)    delete p_cformer;
  if (p_recons)     delete p_recons;

}

int Cluster_Formation_Handler::FormClusters(Blob * blob) {

  if (blob==NULL) return 1;
  assert(m_partlists.empty() && m_clulists.empty());

  msg_Debugging()<<std::endl<<std::endl<<std::endl<<std::endl
		 <<"=================================================================="<<std::endl
		 <<"=================================================================="<<std::endl
		 <<"=================================================================="<<std::endl
		 <<"In "<<METHOD<<": hadronize "<<blob->NInP()<<" partons."<<std::endl
		 <<(*blob)<<std::endl;
  if (!ExtractSinglets(blob))      { return -1; }
  Histogram   * histoxb;
  LPPL_Iterator pplit;
  PPL_Iterator  pit;
  for(pplit=m_partlists.begin(); pplit!=m_partlists.end(); ++pplit) {
    for(pit=(*pplit)->begin(); pit!=(*pplit)->end(); ++pit) {
      if (m_analyse && (*pit)->m_flav.Kfcode()==5) {
	histoxb = (m_histograms.find(string("XB_bquark_before_shift")))->second;
	histoxb->Insert((*pit)->m_mom[0]/45.625);
      }
    }
  }
  if (!ShiftOnMassShells())        { return -1; }
  for(pplit=m_partlists.begin(); pplit!=m_partlists.end(); ++pplit) {
    for(pit=(*pplit)->begin(); pit!=(*pplit)->end(); ++pit) {
      if (m_analyse && (*pit)->m_flav.Kfcode()==5) { 
	histoxb = (m_histograms.find(string("XB_bquark_after_shift")))->second;
	histoxb->Insert((*pit)->m_mom[0]/45.625);
      }
    }
  }
  if (!FormOriginalClusters())     { return -1; }
  if (!ApplyColourReconnections()) { return 0; }
  if (!MergeClusterListsIntoOne()) { return 0; }
  for(Cluster_Iterator cit=p_clulist->begin(); cit!=p_clulist->end(); cit++) {
    if (m_analyse && (*cit)->GetTrip()->m_flav.Kfcode()==5) { 
      histoxb = (m_histograms.find(string("XB_bquark_in_cluster")))->second;
      histoxb->Insert((*cit)->GetTrip()->m_mom[0]/45.625);
      msg_Debugging()<<"Clusters (b)    : "<<(*cit)->GetTrip()->m_mom<<"."<<std::endl;
    }
    if (m_analyse && (*cit)->GetAnti()->m_flav.Kfcode()==5) {
      histoxb = (m_histograms.find(string("XB_bquark_in_cluster")))->second;
      histoxb->Insert((*cit)->GetAnti()->m_mom[0]/45.625);
      msg_Debugging()<<"Clusters (bbar) : "<<(*cit)->GetAnti()->m_mom<<"."<<std::endl;
    }
  }
  if (!ClustersToHadrons(blob))    { return -1; }

  return 1;
}


void Cluster_Formation_Handler::Reset() {
  if(!m_partlists.empty()) {
    while (!m_partlists.empty()) {
      delete m_partlists.front();
      m_partlists.pop_front();
    }
    m_partlists.clear();
  }
  if(!m_clulists.empty()) {
    for(size_t j=0; j<m_clulists.size(); ++j) {
      assert(m_clulists[j]);
      Cluster_List& clist=*m_clulists[j];
      while(!clist.empty()) {
#ifdef memchecker
	msg_Debugging()<<"@@@ Delete cluster "<<clist.back()<<" in "<<METHOD<<"."<<std::endl;
#endif
        clist.pop_back();
      }
      delete &clist;
    }
    m_clulists.clear();
  }
}


bool Cluster_Formation_Handler::ExtractSinglets(Blob * blob)
{
  Proto_Particle_List * pli(NULL);
  bool       construct(false);
  int        col1, col2;
  Particle * part;
  for (int i=0;i<blob->NInP();i++) {
    part = blob->InParticle(i); 
    if ((part->Status()!=part_status::active && part->Status()!=part_status::fragmented) || 
	(part->GetFlow(1)==0 && part->GetFlow(2)==0)) continue;
    if (construct) {
      if (part->GetFlow(2)==col1) {
	Proto_Particle * copy = new Proto_Particle(part->Flav(),part->Momentum(),'L');
#ifdef memchecker
	std::cout<<"### New Proto_Particle ("
		 <<copy<<"/"<<part->Flav()<<") from "<<METHOD<<"."<<std::endl;
#endif
	pli->push_back(copy);
	col1 = part->GetFlow(1);
	if (col1==col2) construct = false;
      }
      else {
	msg_Tracking()<<"ERROR in "<<METHOD<<":"<<std::endl
		      <<"   Assumed everything okay with blob."
		      <<std::endl<<(*blob)<<std::endl
		      <<"   will try new event."<<std::endl;
	Reset();
	return false;
      }
    }
    else {
      col1 = part->GetFlow(1);
      col2 = part->GetFlow(2);
      pli  = new Proto_Particle_List;
      Proto_Particle * copy = new Proto_Particle(part->Flav(),part->Momentum(),'L');
#ifdef memchecker
      std::cout<<"### New Proto_Particle ("
	       <<copy<<"/"<<part->Flav()<<") from "<<METHOD<<"."<<std::endl;
#endif
      pli->push_back(copy);
      m_partlists.push_back(pli);
      construct = true;
    }
  }
  return true;
}

bool Cluster_Formation_Handler::ShiftOnMassShells() {
  ListOfPPLs shiftables, nonshiftables;
  LPPL_Iterator pplit;
  PPL_Iterator  pit;
  for(pplit=m_partlists.begin(); pplit!=m_partlists.end(); ++pplit) {
    Vec4D  mom(0.,0.,0.,0.);
    double mass(0.);
    for(pit=(*pplit)->begin(); pit!=(*pplit)->end(); ++pit) {
      mom  += (*pit)->m_mom;
      mass += hadpars.GetConstituents()->Mass((*pit)->m_flav);
    }
    
    if(mom.Abs2()>sqr(mass)) shiftables.push_back(new Proto_Particle_List(**pplit));
    else nonshiftables.push_back(new Proto_Particle_List(**pplit));
  }
  Proto_Particle_List * pplin;
  while(!nonshiftables.empty()) {
    bool takefromshift(false);
    if(nonshiftables.size()==1) {
      if(shiftables.empty()) {
	Reset();
	delete nonshiftables.front();
	return false;
      }
      pplin=SelectFromList(&shiftables);
      takefromshift=true;
    }
    else pplin=new Proto_Particle_List;
    while(!nonshiftables.empty()) {
      pplin->splice(pplin->end(),*nonshiftables.front());
      nonshiftables.pop_front();
    }
    Vec4D  mom(0.,0.,0.,0.);
    double mass(0.);
    for(pit=pplin->begin(); pit!=pplin->end(); ++pit) {
      mom  += (*pit)->m_mom;
      mass += hadpars.GetConstituents()->Mass((*pit)->m_flav);
    }
    if(mom.Abs2()<sqr(mass)) {
      if(takefromshift) {
 	shiftables.remove(pplin);
	nonshiftables.push_back(pplin);
      }
      else nonshiftables.push_back(pplin);
    }
    else { if(!takefromshift) shiftables.push_back(pplin);}
  }

  /*
    for(pplit=m_partlists.begin();pplit!=m_partlists.end();++pplit)
    std::cout<<**pplit<<std::endl;
    PRINT_INFO("--------------");
    for(pplit=shiftables.begin();pplit!=shiftables.end();++pplit)
    std::cout<<**pplit<<std::endl;
  */

  assert(nonshiftables.empty());

  while(!shiftables.empty()) {
    Proto_Particle_List * pplist=shiftables.front();
    if(!ShiftList(pplist)) {
      while(!shiftables.empty()) {
	delete shiftables.front();
	shiftables.pop_front();
      }
      Reset();
      return false;
    }
    shiftables.pop_front();
  }

  return true;
}

Proto_Particle_List * Cluster_Formation_Handler::SelectFromList(ListOfPPLs * lppl,
								Proto_Particle_List * ppl)
{
  double maxmass(0.0);
  Proto_Particle_List * winner(NULL);
  for (LPPL_Iterator pplit=lppl->begin();pplit!=lppl->end();pplit++) {
    if (ppl!=NULL && (*pplit)==ppl) continue;
    Vec4D mom(0.,0.,0.,0.);
    for (PPL_Iterator pit=(*pplit)->begin();pit!=(*pplit)->end();pit++) {
      mom += (*pit)->m_mom;
    }
    if (mom.Abs2()>maxmass) {
      winner  = (*pplit);
      maxmass = mom.Abs2();
    }
  }
  return winner;
}

bool Cluster_Formation_Handler::ShiftList(Proto_Particle_List * pl)
{
  size_t number(pl->size());
  if (number<2) return true;
  std::vector<Vec4D>  momenta(number);
  std::vector<double> masses(number);
  int k(0);
  Flavour flav;
  PPL_Iterator pit;

#ifdef AHAmomcheck
  Vec4D checkbef(0.,0.,0.,0.), checkaft(0.,0.,0.,0.);
#endif
  for (pit=pl->begin();pit!=pl->end();++pit,++k) {
    flav       = (*pit)->m_flav;
    momenta[k] = (*pit)->m_mom;
#ifdef AHAmomcheck
    checkbef  += momenta[k];
#endif
    masses[k]  = hadpars.GetConstituents()->Mass(flav);
  }
  if (!hadpars.AdjustMomenta(number,&momenta.front(),&masses.front()))  {
    return false;
  }
  k = 0;
  for (pit=pl->begin();pit!=pl->end();++pit,++k) {
    (*pit)->m_mom = momenta[k]; 
#ifdef AHAmomcheck
    checkaft += (*pit)->m_mom;
#endif
  }

#ifdef AHAmomcheck
  if (dabs((checkbef-checkaft).Abs2())>1.e-12) {
    msg_Error()<<METHOD<<" yields momentum violation : "<<std::endl
	       <<checkbef<<" - "<<checkaft<<" --> "<<(checkbef-checkaft).Abs2()<<std::endl;
  }
  else msg_Tracking()<<METHOD<<" conserves momentum."<<std::endl;
#endif

  return true;
}

bool Cluster_Formation_Handler::FormOriginalClusters()
{
  Cluster_List* clist=NULL;
  LPPL_Iterator pplit;

  while (!m_partlists.empty()) {
    pplit=m_partlists.begin();
    msg_Tracking()<<"======= "<<METHOD<<" for :"<<std::endl
		  <<(**pplit)<<std::endl;
    if(p_gludecayer->DecayList(*pplit)) {
      clist = new Cluster_List;
      p_cformer->ConstructClusters(*pplit,clist);
      m_clulists.push_back(clist);
      pplit=m_partlists.erase(pplit);
    }
    else {
      msg_Error()<<"WARNING in "<<METHOD<<":"<<std::endl
		 <<"   Could not form a suitable list after gluon decays from :"<<std::endl
		 <<(**pplit)
		 <<"   Try a new event."<<std::endl;
      Reset();
      exit(1);
      return false;
    }
  }

  if(m_analyse) {
    for(size_t j=0; j<m_clulists.size(); ++j) {
      clist=m_clulists[j];
      Histogram* histomass=
	(m_histograms.find(string("Cluster_Mass_Formation")))->second;
      Histogram* histonumb=
	(m_histograms.find(string("Cluster_Number_Formation")))->second;
      histonumb->Insert(clist->size());
      for(Cluster_Iterator cit=clist->begin(); cit!=clist->end(); cit++) {
	histomass->Insert((*cit)->Mass());
      }
    }
  }

  return true;
}


bool Cluster_Formation_Handler::ApplyColourReconnections()
{
  std::vector<Cluster_List *>::iterator clit1, clit2;
  if (m_single_cr) {
    for (clit1=m_clulists.begin();clit1!=m_clulists.end();clit1++) p_recons->Singlet_CR((*clit1));
  }
  if (m_double_cr && m_clulists.size()>1) {
    clit1 = m_clulists.begin(); 
    do {
      clit2 = clit1; clit2++;
      do {
	p_recons->Two_Singlet_CR((*clit1),(*clit2));
	clit2++;
      } while (clit2!=m_clulists.end());
      clit1++;
    } while (clit1!=m_clulists.end());
  }

  Histogram * histomass;
  if (m_analyse) {
    histomass = (m_histograms.find(string("Cluster_Mass_Reconnections")))->second;
    for (clit1=m_clulists.begin();clit1!=m_clulists.end();clit1++) {
      for (Cluster_Iterator cit=(*clit1)->begin();cit!=(*clit1)->end();cit++) {
	histomass->Insert((*cit)->Mass());
      }
    }
  }
  return true;
}

bool Cluster_Formation_Handler::MergeClusterListsIntoOne() {
  assert(p_clulist->empty());
  for(size_t j=0; j<m_clulists.size(); ++j)
    p_clulist->splice(p_clulist->end(),*m_clulists[j]);
  for(size_t j=0; j<m_clulists.size(); ++j)
    delete m_clulists[j];
  m_clulists.clear();

  return true;

}


bool Cluster_Formation_Handler::ClustersToHadrons(Blob * blob)
{
#ifdef AHAmomcheck
  Vec4D checkbef(0.,0.,0.,0.);
  for (Cluster_Iterator cit=p_clulist->begin();cit!=p_clulist->end();cit++) {
    checkbef += (*cit)->Momentum();
  }
#endif

  if (!p_softclusters->TreatClusterList(p_clulist,blob)) {
    msg_Debugging()<<"Error in "<<METHOD<<" : "<<std::endl
	       <<"   Did not find a kinematically allowed solution for the cluster list."<<std::endl
	       <<"   Will trigger a new event."<<std::endl;
    return false;
  }

#ifdef AHAmomcheck
  Vec4D checkaft(0.,0.,0.,0.);
  for (Cluster_Iterator cit=p_clulist->begin();cit!=p_clulist->end();cit++) {
    checkaft += (*cit)->Momentum();
  }
  for (short int i=0;i<blob->NOutP();i++) {
    checkaft += blob->OutParticle(i)->Momentum();
  }
  double Q2(dabs((checkbef-checkaft).Abs2()));
  if (Q2>1.e-12 || IsNan(Q2)) {
    msg_Error()<<METHOD<<" yields a momentum violation : "<<std::endl
	       <<"   "<<checkbef<<" - "<<checkaft<<" --> "
	       <<(checkbef-checkaft).Abs2()<<"."<<std::endl;
  }
  else msg_Tracking()<<METHOD<<" conserves momentum:"
		     <<checkbef<<" --> "<<checkaft<<"."<<std::endl;
#endif

  Histogram * histomass, * histonumb;
  if (m_analyse) {
    histomass = (m_histograms.find(string("Cluster_Mass_Transformed")))->second;
    histonumb = (m_histograms.find(string("Cluster_Number_Transformed")))->second;
    int numb  = p_clulist->size();
    for (Cluster_Iterator cit=p_clulist->begin();cit!=p_clulist->end();cit++) {
      histomass->Insert((*cit)->Mass());
    }
    histonumb->Insert(numb);
  }
  return true;
}
