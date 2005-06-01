#include "Cluster_Formation_Handler.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Cluster_Formation_Handler::Cluster_Formation_Handler(bool ana) :
  m_single_cr(true), m_double_cr(false), 
  p_gludecayer(new Gluon_Decayer()), 
  p_cformer(new Cluster_Former()),
  p_recons(new Colour_Reconnections()), 
  p_ctransformer(new Cluster_Transformer()),
  p_clulist(new Cluster_List),
  m_analyse(ana)
{ 
  if (m_analyse) {
    m_histograms[string("Cluster_Mass_Formation")]     = new Histogram(0,0.,50.,100);
    m_histograms[string("Cluster_Mass_Reconnections")] = new Histogram(0,0.,50.,100);
    m_histograms[string("Cluster_Mass_Transformed")]   = new Histogram(0,0.,50.,100);
    m_histograms[string("Cluster_Number_Formation")]   = new Histogram(0,0.,20.,20);
    m_histograms[string("Cluster_Number_Transformed")] = new Histogram(0,0.,20.,20);
  }
}


Cluster_Formation_Handler::~Cluster_Formation_Handler()
{
  if (m_analyse) {
    Histogram * histo;
    string name;
    for (map<string,Histogram *>::iterator hit=m_histograms.begin();
	 hit!=m_histograms.end();hit++) {
      histo = hit->second;
      name  = hit->first+string(".dat");
      histo->Output(name);
      delete histo;
    }
    m_histograms.clear();
  }

  Reset();
  if (p_gludecayer)   { delete p_gludecayer;   p_gludecayer = NULL;      }
  if (p_cformer)      { delete p_cformer;      p_cformer    = NULL;      }
  if (p_recons)       { delete p_recons;       p_recons    = NULL;       }
  if (p_ctransformer) { delete p_ctransformer; p_ctransformer    = NULL; }
  if (p_clulist)      { delete p_clulist;      p_clulist    = NULL;      }
}

void Cluster_Formation_Handler::Reset()
{
  if (m_partlists.size()>0) {
    for (std::vector<Part_List *>::iterator plit=m_partlists.begin();
	 plit!=m_partlists.end();plit++) {
      if ((*plit)->size()>0) {
	do { delete (*plit)->back(); (*plit)->pop_back(); } while (!(*plit)->empty());
      }
      delete (*plit); 
    }
  }
  m_partlists.clear();

  if (m_clulists.size()>0) {
    for (std::vector<Cluster_List *>::iterator clit=m_clulists.begin();
	 clit!=m_clulists.end();clit++) {
      if ((*clit)->size()>0) {
	do { delete (*clit)->back(); (*clit)->pop_back(); } while (!(*clit)->empty());
      }
      delete (*clit); 
    }
  }
  m_clulists.clear();
}

Blob * Cluster_Formation_Handler::FormClusters(Blob_List * bl) 
{
  if (bl==NULL) return false;
  Reset();
  p_blob = new Blob();
  p_blob->SetType(btp::Cluster_Formation);
  p_blob->SetId();
  bl->push_back(p_blob);
  msg.Tracking()<<"Extract -----------------------------------------------------------"<<endl;
  ExtractSinglets(bl);
  msg.Tracking()<<"Form Original -----------------------------------------------------"<<endl;
  FormOriginalClusters();
  msg.Tracking()<<"Color Recons ------------------------------------------------------"<<endl;
  ApplyColourReconnections();
  msg.Tracking()<<"Clusters2Hadrons --------------------------------------------------"<<endl;
  ClustersToHadrons();
  msg.Tracking()<<"One List ----------------------------------------------------------"<<endl;
  MergeClusterListsIntoOne();
  msg.Tracking()<<"Leave Formation ---------------------------------------------------"<<endl;
  return p_blob;
}


void Cluster_Formation_Handler::ExtractSinglets(Blob_List * bl)
{
  Particle  * part1, * part2;
  Part_List * pl = new Part_List;
  for (Blob_List::iterator blit=bl->begin();blit!=bl->end();++blit) {
    if ((*blit)->Type()==btp::FS_Shower || (*blit)->Type()==btp::IS_Shower) {
      for (int i=0;i<(*blit)->NOutP();i++) {
	part2 = (*blit)->OutParticle(i); 
	if (part2->Status()==1 && 
	    (part2->GetFlow(1)!=0 || part2->GetFlow(2)!=0)) {
	  p_blob->AddToInParticles(part2);
	  part1 = new Particle(-1,part2->Flav(),part2->Momentum(),'L');
	  part1->SetNumber(0);
	  part1->SetFlow(1,part2->GetFlow(1));
	  part1->SetFlow(2,part2->GetFlow(2));
	  pl->push_back(part1);
	}
      }
    }
  }


  int  col1, col2;
  bool hit1, hit2;
  Part_List * pli=NULL;
  do {
    hit1 = false;
    for (Part_Iterator pit=pl->begin();pit!=pl->end();++pit) {
      col1 = (*pit)->GetFlow(1);
      col2 = (*pit)->GetFlow(2);
      if (col1!=0 && col2==0) {
	hit1 = true;
	pli  = new Part_List;
	pli->push_back((*pit));
	pit  = pl->erase(pit);
	m_partlists.push_back(pli);
	do {
	  hit2 = false;
	  for (Part_Iterator pit1=pl->begin();pit1!=pl->end();++pit1) {
	    if ((int)((*pit1)->GetFlow(2))==col1) {
	      col1 = (*pit1)->GetFlow(1);
	      pli->push_back((*pit1));
	      pit1 = pl->erase(pit1);
	      hit2 = true;
	      break;
	    }
	  }
	} while (hit2 && col1!=0);
      }
      if (hit1) break;
    }
  } while(pl->size()>0);

  pl->clear();
  delete pl;
}


void Cluster_Formation_Handler::FormOriginalClusters() 
{
  Cluster_List * clist=NULL;
  for (std::vector<Part_List *>::iterator plit=m_partlists.begin();
       plit!=m_partlists.end();plit++) {
    clist = new Cluster_List;
    p_gludecayer->DecayList(*plit);
    p_cformer->ConstructClusters(*plit,clist);
    m_clulists.push_back(clist);    
  }

  Histogram * histomass, * histonumb;
  if (m_analyse) {
    histomass = (m_histograms.find(string("Cluster_Mass_Formation")))->second;
    histonumb = (m_histograms.find(string("Cluster_Number_Formation")))->second;
    histonumb->Insert(clist->size());
    for (Cluster_Iterator cit=clist->begin();cit!=clist->end();cit++) {
      histomass->Insert((*cit)->Mass());
    }
  }
}


void Cluster_Formation_Handler::ApplyColourReconnections()
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
}

void Cluster_Formation_Handler::ClustersToHadrons()
{
  std::vector<Cluster_List *>::iterator clit,clit1;
  Cluster * clu;
  msg.Tracking()<<"   Start with "<<m_clulists.size()<<" lists."<<endl;
  for (clit=m_clulists.begin();clit!=m_clulists.end();) {
    if ((*clit)->size()==1) {
      msg.Tracking()<<"      List with 1 cluster only."<<endl;
      if (p_ctransformer->TreatSingleCluster((*clit),p_blob)) {
	if ((*clit)->size()!=0) {
	  Cluster_List * clist = NULL;
	  int maxsize = 10000;
	  for (clit1=m_clulists.begin();clit1!=m_clulists.end();clit1++) {
	    if ((*clit1)->size()<maxsize) clist = (*clit1);
	  }
	  clist->push_back(*(*clit)->begin());
	}
	clit=m_clulists.erase(clit);
      }
      else clit++;
    }
    else clit++;
  }
  msg.Tracking()<<"   Continue with "<<m_clulists.size()<<" lists."<<endl;
  for (clit=m_clulists.begin();clit!=m_clulists.end();clit++) p_ctransformer->TreatClusterList((*clit),p_blob);

  Histogram * histomass, * histonumb;
  if (m_analyse) {
    histomass = (m_histograms.find(string("Cluster_Mass_Transformed")))->second;
    histonumb = (m_histograms.find(string("Cluster_Number_Transformed")))->second;
    int numb  = 0;
    for (clit1=m_clulists.begin();clit1!=m_clulists.end();clit1++) {
      numb+=(*clit1)->size();
      for (Cluster_Iterator cit=(*clit1)->begin();cit!=(*clit1)->end();cit++) {
	histomass->Insert((*cit)->Mass());
      }
    }
    histonumb->Insert(numb);
  }
}

void Cluster_Formation_Handler::MergeClusterListsIntoOne()
{
  if (!p_clulist->empty()) {
    do {
      if (p_clulist->back()) delete p_clulist->back();
      p_clulist->pop_back();
    } while (!p_clulist->empty());
    p_clulist->clear();
  }

  if (m_clulists.size()>0) {
    for (std::vector<Cluster_List *>::iterator clit=m_clulists.begin();
	 clit!=m_clulists.end();clit++) {
      msg.Tracking()<<"New cluster list : "<<(*clit)->size()<<" : ";
      if (!(*clit)->empty()) {
	do {
	  p_clulist->push_front((*clit)->back());
	  (*clit)->pop_back();
	  msg.Tracking()<<(*clit)->size()<<" ";
	} while (!(*clit)->empty());
	msg.Tracking()<<endl;
      }
    }
  }
  Reset();
}
