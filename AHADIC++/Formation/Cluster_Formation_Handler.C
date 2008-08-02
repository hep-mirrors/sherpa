#include <cassert>
#include "Cluster_Formation_Handler.H"
#include "Hadronisation_Parameters.H"

namespace AHADIC {
  bool triplet(SP(Proto_Particle) pp) {
    return ((pp->m_flav.IsQuark() && !pp->m_flav.IsAnti()) ||
	    (pp->m_flav.IsDiQuark() && pp->m_flav.IsAnti()) );
  }
  
  bool antitriplet(SP(Proto_Particle) pp) {
    return ((pp->m_flav.IsQuark() && pp->m_flav.IsAnti()) ||
	    (pp->m_flav.IsDiQuark() && !pp->m_flav.IsAnti()) );
  }
}

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Cluster_Formation_Handler::Cluster_Formation_Handler(Cluster_List* clulist,
						     Cluster_List* prilist,
						     bool ana) :
  m_single_cr(true), m_double_cr(false),
  p_gludecayer(new Gluon_Decayer(hadpars.GetSplitter(),ana)), 
  p_cformer(new Cluster_Former()),
  p_recons(new Colour_Reconnections(0,0,1.)), 
  p_softclusters(hadpars.GetSoftClusterHandler()),
  p_clulist(clulist), p_primaries(prilist),
  m_analyse(ana)
{ 
  if (m_analyse) {
    m_histograms[string("Cluster_Mass_Formation")]     = new Histogram(0,0.,100.,200);
    m_histograms[string("Cluster_Mass_Reconnections")] = new Histogram(0,0.,100.,200);
    m_histograms[string("Cluster_Mass_Transformed")]   = new Histogram(0,0.,100.,200);
    m_histograms[string("Cluster_Number_Formation")]   = new Histogram(0,0.,20.,20);
    m_histograms[string("Cluster_Number_Transformed")] = new Histogram(0,0.,20.,20);
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

  if (!ExtractSinglets(blob))      return -1;
  if (!ShiftOnMassShells())        return -1;
  if (!FormOriginalClusters())     return -1;
  if (!ApplyColourReconnections()) return 0;
  if (!MergeClusterListsIntoOne()) return 0;
  if (!ClustersToHadrons(blob))    return -1;

  return 1;

}


void Cluster_Formation_Handler::Reset() {
  if(!m_partlists.empty()) {
    m_partlists.clear();
  }
  if(!m_clulists.empty()) {
    for(size_t j=0; j<m_clulists.size(); ++j) {
      assert(m_clulists[j]);
      Cluster_List& clist=*m_clulists[j];
      while(!clist.empty()) {
#ifdef memchecker
	std::cout<<"@@@ Delete cluster "<<clist.back()<<" in "<<METHOD<<"."<<std::endl;
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
  SP(Proto_Particle_List) pli(NULL);
  bool       construct(false);
  int        col1, col2;
  Particle * part;
  for (int i=0;i<blob->NInP();i++) {
    part = blob->InParticle(i); 
    if ((part->Status()!=part_status::active && part->Status()!=part_status::fragmented) || 
	(part->GetFlow(1)==0 && part->GetFlow(2)==0)) continue;
    if (construct) {
      if (part->GetFlow(2)==col1) {
	SP(Proto_Particle) copy = new Proto_Particle(part->Flav(),part->Momentum(),'L');
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
      SP(Proto_Particle) copy = new Proto_Particle(part->Flav(),part->Momentum(),'L');
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
    SP(Proto_Particle_List) copy=new Proto_Particle_List(**pplit);
    if(mom.Abs2()>sqr(mass)) shiftables.push_back(copy);
    else nonshiftables.push_back(copy);
  }
  SP(Proto_Particle_List) pplin;
  while(!nonshiftables.empty()) {
    bool takefromshift(false);
    if(nonshiftables.size()==1) {
      if(shiftables.empty()) {
	Reset();
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
    SP(Proto_Particle_List) pplist=shiftables.front();
    if(!ShiftList(pplist)) {
      shiftables.pop_front();
      while(!shiftables.empty()) {
	shiftables.pop_front();}
      Reset();
      return false;
    }
    shiftables.pop_front();
  }

  return true;
}

SP(Proto_Particle_List) Cluster_Formation_Handler::SelectFromList(ListOfPPLs * lppl,
								  SP(Proto_Particle_List) ppl)
{
  double maxmass(0.0);
  SP(Proto_Particle_List) winner(NULL);
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

bool Cluster_Formation_Handler::ShiftList(SP(Proto_Particle_List) pl)
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
    msg_Out()<<METHOD<<" yields momentum violation : "<<std::endl
	     <<checkbef<<" - "<<checkaft<<" --> "<<(checkbef-checkaft).Abs2()<<std::endl;
  }
  else msg_Debugging()<<METHOD<<" conserves momentum."<<std::endl;
#endif

  return true;
}

bool Cluster_Formation_Handler::FormOriginalClusters()
{
  Cluster_List* clist=NULL;
  LPPL_Iterator pplit;

  while (!m_partlists.empty()) {
    pplit=m_partlists.begin();
    if(p_gludecayer->DecayList(*pplit)) {
      clist = new Cluster_List;
      p_cformer->ConstructClusters(*pplit,clist);
      m_clulists.push_back(clist);
      pplit=m_partlists.erase(pplit);
    }
    else {
      msg_Tracking()<<"WARNING in "<<METHOD<<":"<<std::endl
		    <<"   Could not form a suitable list after gluon decays from :"<<std::endl
		    <<(**pplit)
		    <<"   Try a new event."<<std::endl;
      Reset();
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

  assert(p_clulist->empty() && p_primaries->empty());

  for(size_t j=0; j<m_clulists.size(); ++j)
    p_clulist->splice(p_clulist->end(),*m_clulists[j]);
  for(size_t j=0; j<m_clulists.size(); ++j)
    delete m_clulists[j];
  m_clulists.clear();

  //Additional book-keeping to have the memory under control.
  for(Cluster_Iterator cit=p_clulist->begin(); cit!=p_clulist->end(); cit++)
    p_primaries->push_back(*cit);

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
    msg_Tracking()<<"Error in "<<METHOD<<" : "<<std::endl
		  <<"   Did not find a kinematically allowed solution for the cluster list."<<std::endl
		  <<"   Will trigger a new event."<<std::endl;
    return false;
  }
  if (msg->LevelIsDebugging()) {
    msg_Out()<<"         remaining cluster list with size "<<p_clulist->size()<<std::endl
	     <<"======================================================="<<std::endl;
  }

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
