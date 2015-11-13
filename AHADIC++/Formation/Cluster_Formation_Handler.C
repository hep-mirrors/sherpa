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
  m_single_cr(hadpars->Get(string("colour_reconnections"))==1), 
  m_double_cr(false),
  p_gludecayer(new Gluon_Decayer(ana)), 
  p_cformer(new Cluster_Former()),
  p_recons(new Colour_Reconnections(2,1,hadpars->Get(string("pt02")))), 
  p_softclusters(hadpars->GetSoftClusterHandler()),
  p_clulist(clulist), m_analyse(true)
{ 
  if (m_analyse) {
    m_histograms[string("Cluster_Mass_Formation")]     = new Histogram(0,0.,100.,200);
    m_histograms[string("Cluster_Mass_Reconnections")] = new Histogram(0,0.,100.,200);
    m_histograms[string("Cluster_Mass_Transformed")]   = new Histogram(0,0.,100.,200);
    m_histograms[string("Cluster_Number_Formation")]   = new Histogram(0,0.,20.,20);
    m_histograms[string("Cluster_Number_Transformed")] = new Histogram(0,0.,20.,20);
    m_histograms[string("Forward_Number")] = new Histogram(0,0.,20.,20);
    m_histograms[string("Central_Number")] = new Histogram(0,0.,20.,20);
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
  //msg_Out()<<"##################################################\n"
  //  	   <<"##################################################\n"
  //  	   <<"##################################################\n"
  // 	   <<(*blob)<<"\n";
  if (!m_partlists.empty() || !m_clulists.empty()) {
    // This might introduce a tiny (!) memory leak, but should fix a
    // double free error when retrying the event.
    m_partlists.clear();
    m_clulists.clear();
  }
   Vec4D cms(0.,0.,0.,0.);
   for (size_t i=0;i<blob->NInP();i++) cms += blob->InParticle(i)->Momentum();
   // msg_Out()<<"\n\n\n\n"
   // 	    <<"====================================================\n"
   // 	    <<"====================================================\n"
   // 	    <<"====================================================\n"
   // 	    <<"In "<<METHOD<<": hadronize "<<blob->NInP()<<" partons "
   // 	    <<"with E = "<<sqrt(cms.Abs2())<<".\n"<<(*blob)<<"\n";
  if (!ExtractSinglets(blob))      { Reset(); return -1; }
  if (!ShiftOnMassShells())        { Reset(); return -1; }
  if (!FormOriginalClusters())     { Reset(); return -1; }
  if (!ApplyColourReconnections()) { Reset(); return 0; }
  if (!MergeClusterListsIntoOne()) { Reset(); return 0; }
  if (!ClustersToHadrons(blob))    { 
    Reset(); return -1; 
  }
  
  // Vec4D blobmom(0.,0.,0.,0.);
  // for (size_t i(0);i<blob->NOutP();i++) 
  //   blobmom+=blob->OutParticle(i)->Momentum();
  // msg_Out()<<"____________________________________________________\n"
  // 	   <<"____________________________________________________\n"
  // 	   <<"Cluster list after all merging etc.:\n"<<(*p_clulist)
  // 	   <<"____________________________________________________\n"
  // 	   <<"Blob momentum: "<<blobmom<<";\n"<<(*blob)<<"\n"
  // 	   <<"____________________________________________________\n"
  // 	   <<"____________________________________________________\n";
  
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
      while(!clist.empty()) clist.pop_back();
      delete &clist;
    }
    m_clulists.clear();
  }
}


bool Cluster_Formation_Handler::ExtractSinglets(Blob * blob)
{
  Particle   * part(NULL);
  // Incoming are partons et al from the shower --- no outgoing particles at that point
  // We look for closed colour structures
  //msg_Out() << (*blob);

  //// Filtering
  std::vector<int> temp;
  for (unsigned int i=0;i<blob->NInP();i++) {
    part = blob->InParticle(i); 
    if ((part->Status()!=part_status::active  && part->Status()!=part_status::fragmented)
            || (part->GetFlow(1)==0 && part->GetFlow(2)==0)) continue;
    temp.push_back(i);
  }
 

  /// One Particle_List for each closed colour structure
  Proto_Particle_List * pli = new Proto_Particle_List;
  
  /// Tagging colour
  unsigned int col = blob->InParticle(temp[0])->GetFlow(2);
  /// Iterate over particles, except the first, that's always good
  /// We look for "colour rings" here
  unsigned int beg(0);
  for (unsigned int i=1;i<temp.size();i++) {
    part = blob->InParticle(temp[i]);
    // Found a closed colour ring 
    if (col==part->GetFlow(1)) {
      // append copies of the ring's particles to a particle list
      for (unsigned int j=beg;j<=i;j++) {
        /// This bit catches weirdly ordered colour partners and triggers a retry event
        if (j<i && !(blob->InParticle(j+1)->GetFlow(2) == blob->InParticle(j)->GetFlow(1))) return false;
        /// Copy of particle for further usage
        Proto_Particle * copy = new Proto_Particle(
                blob->InParticle(j)->Flav(),
                blob->InParticle(j)->Momentum(),
                blob->InParticle(j)->Info()=='B' ?'B':'L'
                );
        SetInfoTagForPrimaryParticle(copy);
        pli->push_back(copy);
      }
      /// Append the closed colour ring to m_partlists
      m_partlists.push_back(pli);
      /// Prepare a new colour ring
      pli = new Proto_Particle_List;
      beg=i+1; // Make sure we skip the first i of all particles for the next ring
      if (i+1 <temp.size()) { /// Reset tag colour
          col=blob->InParticle(temp[i+1])->GetFlow(2);
      }
    }

  }
  //for (LPPL_Iterator pli=m_partlists.begin();pli!=m_partlists.end();pli++) {
    //msg_Out()<<(*pli)->size()<<"\n";
  //}
  //if (ana) { 
    //if (over) {
      //msg_Out()<<"\n\n"<<Nover<<" interesting particles: "<<id<<"\n"
                 //<<(*blob)<<"\n"<<"\n";
      //m_histograms[string("Forward_Number")]->Insert(Nover);
    //}
    //if (under) {
      //m_histograms[string("Central_Number")]->Insert(Nunder);
    //}
  //}
  return true;
}

void Cluster_Formation_Handler::
SetInfoTagForPrimaryParticle(Proto_Particle * proto) const {
  if (proto->m_info=='B') return;
  proto->m_info=(proto->m_flav.IsQuark() || proto->m_flav.IsDiQuark())?'L':'l';
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
      mass += hadpars->GetConstituents()->Mass((*pit)->m_flav);
    }
    
    if(mom.Abs2()>sqr(mass)) 
      shiftables.push_back(new Proto_Particle_List(**pplit));
    else 
      nonshiftables.push_back(new Proto_Particle_List(**pplit));
  }
  Proto_Particle_List * pplin;
  while(!nonshiftables.empty()) {
    bool takefromshift(false);
    if(nonshiftables.size()==1) {
      if(shiftables.empty()) {
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
      mass += hadpars->GetConstituents()->Mass((*pit)->m_flav);
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

  assert(nonshiftables.empty());

  while(!shiftables.empty()) {
    Proto_Particle_List * pplist=shiftables.front();
    if(!ShiftList(pplist)) {
      msg_Error()<<"Error in "<<METHOD<<" : \n"
		 <<"   Cannot shift particle list collectively.\n"
		 <<(*pplist)<<"\n"
		 <<"   Will trigger new event.\n";
      while(!shiftables.empty()) {
	delete shiftables.front();
	shiftables.pop_front();
      }
      return false;
    }
    shiftables.pop_front();
  }

  return true;
}

Proto_Particle_List * 
Cluster_Formation_Handler::SelectFromList(ListOfPPLs * lppl,
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
    masses[k]  = hadpars->GetConstituents()->Mass(flav);
  }
  if (!hadpars->AdjustMomenta(number,&momenta.front(),&masses.front()))  {
    msg_Error()<<"Warning in "<<METHOD<<".  Could not adjust momenta for:\n";
    for (pit=pl->begin();pit!=pl->end();++pit,++k) {
      msg_Error()<<"   "<<(*pit)->m_flav<<" "
		 <<(*pit)->m_mom<<" ("<<(*pit)->m_mom.Abs2()<<") vs. "
		 <<hadpars->GetConstituents()->Mass((*pit)->m_flav)<<".\n";
    }
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
    msg_Error()<<METHOD<<" yields momentum violation : \n"
	       <<checkbef<<" - "<<checkaft<<" --> "<<(checkbef-checkaft).Abs2()<<"\n";
  }
  else msg_Tracking()<<METHOD<<" conserves momentum.\n";
#endif

  return true;
}

bool Cluster_Formation_Handler::FormOriginalClusters()
{
  Cluster_List* clist=NULL;
  LPPL_Iterator pplit;

  while (!m_partlists.empty()) {
    pplit=m_partlists.begin();
   // msg_Tracking()<<"======= "<<METHOD<<" for :\n"<<(**pplit)<<"\n";
   // msg_Out()<<"========== before gluon splitting:\n"<<(**pplit)<<"\n";
    if(p_gludecayer->DecayList(*pplit)) {
      //msg_Out()<<"========== after gluon splitting:\n"<<(**pplit)<<"\n";
      clist = new Cluster_List;
      p_cformer->ConstructClusters(*pplit,clist);
     // msg_Out()<<"========== cluster list :\n"<<(*clist)<<"\n";
      m_clulists.push_back(clist);
      pplit=m_partlists.erase(pplit);
    }
    else {
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
    for (clit1=m_clulists.begin();clit1!=m_clulists.end();clit1++) 
      p_recons->Singlet_CR((*clit1));
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
    histomass = m_histograms[string("Cluster_Mass_Reconnections")];
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

  msg_Tracking()<<METHOD<<":\n"<<(*p_clulist)<<"\n";
  return true;
}


bool Cluster_Formation_Handler::ClustersToHadrons(Blob * blob)
{
  //msg_Out()<<"====== "<<METHOD<<" for: \n"<<(*p_clulist)<<"\n";
  if (!p_softclusters->TreatClusterList(p_clulist,blob)) return false;

  if (m_analyse) {
    Histogram * histomass, * histonumb;
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
