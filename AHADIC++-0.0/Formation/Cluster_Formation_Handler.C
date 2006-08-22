#include "Cluster_Formation_Handler.H"
#include "Hadronisation_Parameters.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Cluster_Formation_Handler::Cluster_Formation_Handler(bool ana) :
  m_single_cr(true), m_double_cr(false),
  p_gludecayer(new Gluon_Decayer()), 
  p_cformer(new Cluster_Former()),
  p_recons(new Colour_Reconnections(0,0,1.)), 
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

Return_Value::code Cluster_Formation_Handler::FormClusters(Blob * blob,Blob_List * bl) 
{
  p_blob = blob;
  if (bl==NULL) {
    msg.Error()<<"ERROR in "<<METHOD<<":"<<std::endl
	       <<"   Continue with error and hope for the best."<<std::endl;
    return Return_Value::Error;
  }

  Return_Value::code success;
  Reset();
  msg.Tracking()<<"Extract -----------------------------------------------------------"<<endl;
  switch (int(ExtractSinglets(bl))) {
  case int(Return_Value::Retry_Method) : return Return_Value::Retry_Method;
  case int(Return_Value::Success) : 
  default:
    break;
  }
  msg.Tracking()<<"Form Original -----------------------------------------------------"<<endl;
  switch (int(FormOriginalClusters())) {
  case int(Return_Value::Retry_Method) : return Return_Value::Retry_Method;
  case int(Return_Value::Success) : 
  default:
    break;
  }
  msg.Tracking()<<"Color Recons ------------------------------------------------------"<<endl;
  switch (int(ApplyColourReconnections())) {
  case int(Return_Value::Retry_Method) : return Return_Value::Retry_Method;
  case int(Return_Value::Success) : 
  default:
    break;
  }
  msg.Tracking()<<"Clusters2Hadrons --------------------------------------------------"<<endl;
  switch (int(ClustersToHadrons())) {
  case int(Return_Value::Retry_Method) : return Return_Value::Retry_Method;
  case int(Return_Value::Success) : 
  default:
    break;
  }
  msg.Tracking()<<"One List ----------------------------------------------------------"<<endl;
  switch (int(MergeClusterListsIntoOne())) {
  case int(Return_Value::Retry_Method) : return Return_Value::Retry_Method;
  case int(Return_Value::Success) : 
  default:
    break;
  }
  msg.Tracking()<<"Leave Formation ---------------------------------------------------"<<endl;
  bl->push_back(p_blob);
  return Return_Value::Success;
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


Return_Value::code Cluster_Formation_Handler::ExtractSinglets(Blob_List * bl)
{
  Particle  * part1, * part2;
  Part_List * pl = new Part_List;
  for (Blob_List::iterator blit=bl->begin();blit!=bl->end();++blit) {
    if ((*blit)->Type()==btp::FS_Shower || 
	(*blit)->Type()==btp::IS_Shower ||
	(*blit)->Type()==btp::Shower) {
      for (int i=0;i<(*blit)->NOutP();i++) {
	part2 = (*blit)->OutParticle(i); 
	if (part2->Status()==part_status::active && 
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

  return Return_Value::Success;
}


Return_Value::code Cluster_Formation_Handler::FormOriginalClusters() 
{
  Cluster_List * clist=NULL;
  std::vector<Part_List *>::iterator help;
  Vec4D  totvec;
  double totmass;
  int k(0);
  Flavour flav;
  bool rearrange;
  do {
    rearrange = false;
    for (std::vector<Part_List *>::iterator plit=m_partlists.begin();
	 plit!=m_partlists.end();plit++,k++) {
      totmass = 0.;
      totvec  = Vec4D(0.,0.,0.,0.);
      for (Part_Iterator pit=(*plit)->begin();pit!=(*plit)->end();++pit) {
	flav     = (*pit)->Flav();
	totmass += hadpars.GetConstituents()->Mass(flav);
	totvec  += (*pit)->Momentum();
      }
      if (sqr(totmass)>totvec.Abs2()) {
	rearrange = true;
	help = plit;
	help++;
	if (help==m_partlists.end()) {
	  help = plit;
	  help--;
	  (*help)->merge((**plit));
	  (*plit)->clear();
	  delete (*plit);
	  m_partlists.pop_back();	  
	  break;
	}
	else {
	  (*plit)->merge((**help));
	  (*help)->clear();
	  delete (*help);
	  for (int j=k+1;j<m_partlists.size()-1;j++) m_partlists[j]=m_partlists[j+1];
	  m_partlists.pop_back();
	}
      }
    }
  } while(rearrange);
  
  
  for (std::vector<Part_List *>::iterator plit=m_partlists.begin();
       plit!=m_partlists.end();plit++,k++) {
    clist = new Cluster_List;
    //std::cout<<"Part list with "<<(*plit)->size()<<"."<<std::endl;
    if (!p_gludecayer->DecayList(*plit)) {
      msg.Info()<<"WARNING in "<<METHOD<<":"<<std::endl
		<<"   Not enough energy to move partons on their mass shell."<<std::endl
		<<"   Retry the formation procedure."<<std::endl;
      rvalue.IncRetryMethod(METHOD);
      return Return_Value::Retry_Method;
    }
    else {
      p_cformer->ConstructClusters(*plit,clist);
      m_clulists.push_back(clist);    
    }
  }

  Histogram * histomass, * histonumb;
  if (m_analyse) {
    histomass = (m_histograms.find(string("Cluster_Mass_Formation")))->second;
    histonumb = (m_histograms.find(string("Cluster_Number_Formation")))->second;
    //std::cout<<"Insert cluster number : "<<clist->size()<<"."<<std::endl;
    histonumb->Insert(clist->size());
    for (Cluster_Iterator cit=clist->begin();cit!=clist->end();cit++) {
      //std::cout<<"   Insert cluster mass : "<<(*cit)->Mass()<<"."<<std::endl;
      histomass->Insert((*cit)->Mass());
    }
  }
  return Return_Value::Success;
}


Return_Value::code Cluster_Formation_Handler::ApplyColourReconnections()
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
  return Return_Value::Success;
}

Return_Value::code Cluster_Formation_Handler::ClustersToHadrons()
{
  std::vector<Cluster_List *>::iterator clit,clit1;
  Cluster * clu;
  msg.Tracking()<<"   Start with "<<m_clulists.size()<<" lists."<<endl;
  for (clit=m_clulists.begin();clit!=m_clulists.end();) {
    if ((*clit)->size()==1) {
      msg.Tracking()<<"      List with 1 cluster only."<<endl;
      switch (int(p_ctransformer->TreatSingleCluster((*clit),p_blob))) {
      case int(Return_Value::Success):
      case int(Return_Value::Warning):
      case int(Return_Value::Error):
	if ((*clit)->size()!=0) {
	  Cluster_List * clist = NULL;
	  int maxsize = 10000;
	  for (clit1=m_clulists.begin();clit1!=m_clulists.end();clit1++) {
	    if ((*clit1)->size()<maxsize) clist = (*clit1);
	  }
	  clist->push_back(*(*clit)->begin());
	}
	clit=m_clulists.erase(clit);
	break;
      case int(Return_Value::Nothing): 
      default:
	clit++;
	break;
      }
    }
    else clit++;
  }
  msg.Tracking()<<"   Continue with "<<m_clulists.size()<<" lists ..."<<std::endl;
  for (clit=m_clulists.begin();clit!=m_clulists.end();clit++) {
    switch (int(p_ctransformer->TreatClusterList((*clit),p_blob))) {
    case int(Return_Value::Error):
      rvalue.IncRetryMethod(METHOD);
      return Return_Value::Retry_Method;
    case int(Return_Value::Success): 
    default:continue;
    }
  }
  msg.Tracking()<<"                                                   ... done."<<endl;

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
  return Return_Value::Success;
}

Return_Value::code Cluster_Formation_Handler::MergeClusterListsIntoOne()
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
  return Return_Value::Success;
}
