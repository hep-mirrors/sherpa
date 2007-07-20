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
  p_softclusters(new Soft_Cluster_Handler()),
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
  if (p_gludecayer)   { delete p_gludecayer;   p_gludecayer   = NULL; }
  if (p_cformer)      { delete p_cformer;      p_cformer      = NULL; }
  if (p_recons)       { delete p_recons;       p_recons       = NULL; }
  if (p_softclusters) { delete p_softclusters; p_softclusters = NULL; }
  if (p_clulist)      { delete p_clulist;      p_clulist      = NULL; }
}

Return_Value::code Cluster_Formation_Handler::FormClusters(Blob * blob,
							   Blob_List * bl)
{
  Reset();
  if (blob==NULL) Return_Value::Error;
  p_blob = blob;

  switch (int(ExtractSinglets())) {
  case int(Return_Value::Retry_Method) : return Return_Value::Retry_Method;
  case int(Return_Value::Success) : 
  default:
    break;
  }
  switch (int(FormOriginalClusters())) {
  case int(Return_Value::Retry_Method) : return Return_Value::Retry_Method;
  case int(Return_Value::Success) : 
  default:
    break;
  }

  switch (int(ApplyColourReconnections())) {
  case int(Return_Value::Retry_Method) : return Return_Value::Retry_Method;
  case int(Return_Value::Success) : 
  default:
    break;
  }
  switch (int(ClustersToHadrons(bl))) {
  case int(Return_Value::Retry_Method) : return Return_Value::Retry_Method;
  case int(Return_Value::Success) : 
  default:
    break;
  }
  switch (int(MergeClusterListsIntoOne())) {
  case int(Return_Value::Retry_Method) : return Return_Value::Retry_Method;
  case int(Return_Value::Success) : 
  default:
    break;
  }

  return Return_Value::Success;
}


void Cluster_Formation_Handler::Reset()
{
  if (m_partlists.size()>0) {
    for (VPPL_Iterator plit=m_partlists.begin();
	 plit!=m_partlists.end();plit++) delete (*plit); 
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


Return_Value::code Cluster_Formation_Handler::ExtractSinglets()
{
  Particle            * part;
  Proto_Particle_List * pli(NULL);
  bool        construct(false);
  int         col1, col2;

  for (int i=0;i<p_blob->NInP();i++) {
    part = p_blob->InParticle(i); 
    if ((part->Status()!=part_status::active && part->Status()!=part_status::fragmented) || 
	(part->GetFlow(1)==0 && part->GetFlow(2)==0)) continue;
    if (construct) {
      if (part->GetFlow(2)==col1) {
	Proto_Particle * copy = new Proto_Particle(part->Flav(),part->Momentum(),'L');
	control::s_AHAprotoparticles++;
	pli->push_back(copy);
	col1 = part->GetFlow(1);
	if (col1==col2) construct = false;
      }
      else {
	msg_Error()<<"ERROR in "<<METHOD<<":"<<std::endl
		   <<"   Assumed everything okay with blob."
		   <<std::endl<<(*p_blob)<<std::endl;
	abort();
	return Return_Value::Error;
      }
    }
    else {
      col1 = part->GetFlow(1);
      col2 = part->GetFlow(2);
      pli  = new Proto_Particle_List;
      Proto_Particle * copy = new Proto_Particle(part->Flav(),part->Momentum(),'L');
      control::s_AHAprotoparticles++;
      pli->push_back(copy);
      m_partlists.push_back(pli);
      construct = true;
    }
  }
  return Return_Value::Success;
}


Return_Value::code Cluster_Formation_Handler::FormOriginalClusters() 
{
  Cluster_List * clist=NULL;
  VPPL_Iterator pplit,help;
  Vec4D  totvec;
  double totmass;
  Flavour flav;
  bool rearrange;
  do {
    rearrange = false;
    for (pplit=m_partlists.begin();pplit!=m_partlists.end();pplit++) {
      totmass = 0.;
      totvec  = Vec4D(0.,0.,0.,0.);
      for (PPL_Iterator pit=(*pplit)->begin();pit!=(*pplit)->end();++pit) {
	flav     = (*pit)->m_flav;
	totmass += hadpars.GetConstituents()->Mass(flav);
	totvec  += (*pit)->m_mom;
      }
      if (sqr(totmass)>totvec.Abs2()) {
	rearrange = true;
	help = pplit;
	help++;
	if (help==m_partlists.end()) {
	  help = pplit;
	  help--;
	  while (!(*pplit)->empty()) { 
	    (*help)->push_back((*pplit)->front()); 
	    (*pplit)->pop_front(); 
	  }
	  m_partlists.pop_back();	  
	  break;
	}
	else {
	  while (!(*help)->empty())  { 
	    (*pplit)->push_back((*help)->front()); 
	    (*help)->pop_front(); 
	  }
	  pplit=m_partlists.erase(help);
	  break;
	}
      }
    }
  } while(rearrange);
  
  
  while (!m_partlists.empty()) {
    pplit=m_partlists.begin();
    clist = new Cluster_List;
    if (!p_gludecayer->DecayList(*pplit)) {
      msg_Info()<<"WARNING in "<<METHOD<<":"<<std::endl
		<<"   Not enough energy to move partons on their mass shell."<<std::endl
		<<"   Retry the formation procedure."<<std::endl;
      rvalue.IncRetryMethod(METHOD);
      return Return_Value::Retry_Method;
    }
    else {
      p_cformer->ConstructClusters(*pplit,clist);
      m_clulists.push_back(clist);    
    }
    delete (*pplit);
    pplit=m_partlists.erase(pplit);
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

Return_Value::code Cluster_Formation_Handler::ClustersToHadrons(ATOOLS::Blob_List * bl)
{
  std::vector<Cluster_List *>::iterator clit=m_clulists.begin();
  while(clit!=m_clulists.end()) {
    if (!p_softclusters->TreatClusterList((*clit),p_blob)) {
      Cluster_List * clist = NULL;
      int maxsize = 10000;
      for (std::vector<Cluster_List *>::iterator clit1=m_clulists.begin();
	   clit1!=m_clulists.end();clit1++) {
	if ((*clit1)->size()<maxsize && (*clit1)!=(*clit)) clist = (*clit1);
      }
      clist->push_back(*(*clit)->begin());
      clit=m_clulists.erase(clit);
    }
    if ((*clit)->empty()) clit = m_clulists.erase(clit);
    else {
      Cluster_Iterator cit=(*clit)->begin();
      Cluster_Iterator endit=(*clit)->end();
      while (cit!=endit) {
	if (!(*cit)->GetLeft()) {
	  p_blob->AddToOutParticles((*cit)->GetSelf());
	  cit++;
	}
	else {
	  if ((*cit)->GetLeft()->GetSelf()->Flav()==Flavour(kf::none)) {
	    (*clit)->push_back((*cit)->GetLeft());
	    (*clit)->push_back((*cit)->GetRight());
	  }
	  else {
	    p_blob->AddToOutParticles((*cit)->GetSelf());
	    bl->push_back((*cit)->CHHDecayBlob());
	  }
	  if (*cit) delete (*cit);
	  cit = (*clit)->erase(cit);
	}
      }
      clit++;
    }
  }

  

  Histogram * histomass, * histonumb;
  if (m_analyse) {
    histomass = (m_histograms.find(string("Cluster_Mass_Transformed")))->second;
    histonumb = (m_histograms.find(string("Cluster_Number_Transformed")))->second;
    int numb  = 0;
    for (clit=m_clulists.begin();clit!=m_clulists.end();clit++) {
      numb+=(*clit)->size();
      for (Cluster_Iterator cit=(*clit)->begin();cit!=(*clit)->end();cit++) {
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
      if (!(*clit)->empty()) {
	do {
	  p_clulist->push_front((*clit)->back());
	  (*clit)->pop_back();
	} while (!(*clit)->empty());
      }
    }
  }
  Reset();
  return Return_Value::Success;
}
