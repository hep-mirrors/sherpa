#include "Final_Selector.H"
#include "Primitive_Analysis.H"
#include "Message.H"
#include "MyStrStream.H"
#include "Durham_Algorithm.H"

#include <algorithm>

using namespace ANALYSIS;
using namespace ATOOLS;

Final_Selector::Final_Selector(const std::string & inlistname,
			       const std::string & outlistname,
			       int mode) :
  m_inlistname(inlistname),m_outlistname(outlistname),m_ownlist(false), m_extract(false),
  m_mode(mode), p_jetalg(NULL)
{
  m_splitt_flag = false;
  if (mode) {
    p_jetalg = new Durham_Algorithm();  
  }
  else {
    p_jetalg = new Kt_Algorithm();
  }
}



void Final_Selector::AddSelector(const Flavour & fl, const Final_Selector_Data & fs) 
{
  Final_Data_Map::iterator it = m_fmap.find(fl);
  if (it==m_fmap.end()) {
    m_fmap.insert(std::make_pair(fl,fs));
    if (m_extract) m_fmap[fl].keep = false;
  }
  else {
    it->second.eta_min = fs.eta_min; 
    it->second.eta_max = fs.eta_max;
    it->second.et_min  = fs.et_min;
    it->second.pt_min  = fs.pt_min;
    it->second.r_min   = fs.r_min;
  }
  if (fl==Flavour(kf::jet) && fs.r_min>0.) AddSelector(fl,fl,fs);
}

void Final_Selector::AddSelector(const Flavour & flav1, const Flavour & flav2, 
				 const Final_Selector_Data & fs) 
{
  std::pair<Flavour,Flavour> flavs(flav1,flav2);
  Final_Correlator_Map::iterator it = m_cmap.find(flavs);
  if (it==m_cmap.end()) {
    m_cmap.insert(std::make_pair(flavs,fs));
    if (m_extract) m_cmap[flavs].keep = false;
  }
  else {
    std::pair<Flavour,Flavour> flavs1(flav2,flav1);
    Final_Correlator_Map::iterator it1 = m_cmap.find(flavs1);
    if (it1==m_cmap.end()) {
      m_cmap.insert(std::make_pair(flavs,fs));
      if (m_extract) m_cmap[flavs].keep = false;
    }
    else {
      it->second.mass_min = fs.mass_min; 
      it->second.mass_max = fs.mass_max;
      it->second.r_min    = fs.r_min;
    }
  }
}

void Final_Selector::AddSelector(const Flavour & fl, int min, int max) 
{
  Final_Data_Map::iterator it = m_fmap.find(fl);
  if (it==m_fmap.end()) {
    Final_Selector_Data fs;
    fs.min_n = min;  
    fs.max_n = max;  
    if (m_extract) fs.keep = false;
    m_fmap.insert(std::make_pair(fl,fs));
  }
  else {
    it->second.min_n = min;  
    it->second.max_n = max;  
  }
}

void Final_Selector::AddKeepFlavour(const Flavour & fl) 
{
  if (fl==Flavour(kf::lepton)) {
    for (int i=0;i<fl.Size();++i) AddKeepFlavour(fl[i]);
  }

  if(!m_extract) {
    Final_Data_Map::iterator it;
    for (it=m_fmap.begin();it!=m_fmap.end();++it) it->second.keep = false;
    m_extract = true;
  }
  m_fmap[fl].keep = true;
}

void Final_Selector::Output()
{
  std::cout<<"Final_Selector : "<<m_fmap.size()<<"/"<<m_cmap.size()<<":"<<std::endl;
  for (Final_Data_Map::iterator it=m_fmap.begin();it!=m_fmap.end();++it) {
    if (it->first!=Flavour(kf::jet)) 
      std::cout<<" "<<it->first<<" : pt_min = "<<it->second.pt_min<<", eta = "
	       <<it->second.eta_min<<" ... "<<it->second.eta_max<<std::endl;
    else
      std::cout<<" "<<it->first<<" : pt_min = "<<it->second.pt_min<<", eta = "
	       <<it->second.eta_min<<" ... "<<it->second.eta_max
	       <<", jets with ktRunII, r_min = "<<it->second.r_min<<std::endl;
  }
  for (Final_Correlator_Map::iterator it=m_cmap.begin();it!=m_cmap.end();++it) {
    std::cout<<" "<<it->first.first<<" "<<it->first.second<<" : "<<it->second.r_min<<std::endl;
  }
  for (Final_Data_Map::iterator it=m_fmap.begin();it!=m_fmap.end();++it) {
    if ((it->second.min_n>-1) && (it->second.max_n>-1)) {
      std::cout<<" "<<it->first<<" : min = "<<it->second.min_n<<", max = "<<it->second.min_n<<std::endl;
    }
  }
}


bool Final_Selector::PtSelect(const Vec4D & mom, double ptmin) 
{
  if (mom.PPerp()<ptmin) return true;
  return false;
}

bool Final_Selector::EtSelect(const Vec4D & mom, double etmin) 
{
  if (mom.EPerp()<etmin) return true;
  return false;
}

bool Final_Selector::EtaSelect(const Vec4D & mom, double etamin,double etamax) 
{
  double eta = mom.Eta();
  if (eta<etamin || etamax<eta ) return true;
  return false;
}

bool Final_Selector::DeltaRSelect(const Vec4D & p1,const Vec4D & p2,double rmin) 
{
  double deta12 = p1.Eta()-p2.Eta();
  double dphi12 = acos( (p1[1]*p2[1]+p1[2]*p2[2])/(p1.PPerp()*p2.PPerp()) );
  if (sqrt(sqr(deta12)+sqr(dphi12))<rmin) return true;
  return false;
}

bool Final_Selector::MassSelect(const Vec4D & p1,const Vec4D & p2,
				double massmin,double massmax) 
{
  double mass = (p1+p2).Abs2();
  if (mass<massmin || mass>massmax) return true;
  return false;
}

double Final_Selector::DeltaR(const Vec4D & p1,const Vec4D & p2) 
{
  double deta12 = p1.Eta() - p2.Eta();

  double pt1=sqrt(p1[1]*p1[1]+p1[2]*p1[2]);
  double pt2=sqrt(p2[1]*p2[1]+p2[2]*p2[2]);
  double dphi12=acos((p1[1]*p2[1]+p1[2]*p2[2])/(pt1*pt2));
  
  return sqrt(sqr(deta12) + sqr(dphi12));
}




void Final_Selector::Select(Particle_List * pl,Final_Data_Map::iterator it) 
{
  bool hit;
  for (Particle_List::iterator pit=pl->begin();pit!=pl->end();) {
    if ((*pit)->Flav()==it->first) {
      hit = false;
      if (it->second.eta_min!=it->second.eta_max)  
	hit=EtaSelect((*pit)->Momentum(),it->second.eta_min,it->second.eta_max);
      if (it->second.et_min!=0. && !hit) hit=EtSelect((*pit)->Momentum(),it->second.et_min);
      if (it->second.pt_min!=0. && !hit) hit=PtSelect((*pit)->Momentum(),it->second.pt_min);
      if (!hit) ++pit;
      else {
	if (m_ownlist) delete *pit;
	pit = pl->erase(pit);
      }
    }
    else {
      ++pit;
    }
  }
}

void Final_Selector::Select2(Particle_List * pl,Final_Correlator_Map::iterator it) 
{
  if (it->second.r_min<=0.) return;

  Flavour flav1 = it->first.first;
  Flavour flav2 = it->first.second;

  bool hit = false;
  for (Particle_List::iterator pit=pl->begin();pit!=pl->end();++pit) {
    for (Particle_List::iterator pit2=pl->begin();pit2!=pl->end();++pit2) {
      if (flav1.Includes((*pit)->Flav()) && flav2.Includes((*pit2)->Flav()) && pit!=pit2) {
	hit = DeltaRSelect((*pit)->Momentum(),(*pit2)->Momentum(),it->second.r_min);
	break;
      }
    }
  } 
  if (hit) {
    for (Particle_List::iterator pit=pl->begin();pit!=pl->end();) {
      if (m_ownlist) delete *pit;
      pit=pl->erase(pit);
    }
  }
}

void Final_Selector::SelectN(Particle_List * pl,Final_Data_Map::iterator it) 
{
  if (pl->size()==0) return;
  if (it->second.min_n==-1 && it->second.max_n==-1) return;

  int counter=0;
  for (Particle_List::iterator pit=pl->begin();pit!=pl->end();++pit) {
    if ((*pit)->Flav()==it->first) ++counter;
  }
  if ((it->second.min_n>counter && it->second.min_n!=-1) ||
      (it->second.max_n<counter && it->second.max_n!=-1)) {
    // delete list
    for (Particle_List::iterator pit=pl->begin();pit!=pl->end();) {
      if (m_ownlist) delete *pit;
      pit=pl->erase(pit);
    }   
  }
}



void Final_Selector::Extract(Particle_List * pl) 
{
  if (!m_extract) return;
  if (pl->size()==0) return;
  for (Particle_List::iterator pit=pl->begin();pit!=pl->end();) {
    bool remove = true;
    for (Final_Data_Map::iterator it=m_fmap.begin();it!=m_fmap.end();++it) {
      if ((*pit)->Flav()==it->first && it->second.keep) {
	remove = false;
	break;
      }
    }
    if (remove) {
      if (m_ownlist) delete *pit;
      pit=pl->erase(pit);
    }
    else ++pit;
  }
}


void Final_Selector::Evaluate(const Blob_List &,double value, int ncount) {
  Particle_List * pl_in = p_ana->GetParticleList(m_inlistname);
  if (pl_in==NULL) {
    msg.Error()<<" WARNING: particle list "<<m_inlistname<<" not found "<<std::endl;
    return;
  }
//   if (m_fmap.empty() && m_cmap.empty()) {
//     p_ana->AddParticleList(m_outlistname,pl_in);
//     return;
//   }

  Particle_List * pl_out = new Particle_List;
  
  
  // look for kt and after for other selectors
  Final_Data_Map::iterator it =m_fmap.find(Flavour(kf::jet));
  if (it!=m_fmap.end()) {
    if (it->second.r_min>0.) {
      std::vector<double> * diffrates=new std::vector<double>();
      p_jetalg->ConstructJets(pl_in,pl_out,diffrates,it->second.r_min);
      // add leptons
      for (Particle_List::iterator pit=pl_in->begin();pit!=pl_in->end();++pit) {
	if ((*pit)->Flav().IsLepton()) pl_out->push_back(new Particle(*pit));
      }
      m_ownlist=true;
      MyStrStream str;
      str<<"KtJetrates("<<it->second.r_min<<")"<<m_listname;
      std::string key;
      str>>key;
      //      std::cout<<" creating : "<<key<<" with "<<diffrates->size()<<" elements "<<std::endl;
      p_ana->AddData(key,new Blob_Data<std::vector<double> *>(diffrates));
    }
    else {
      // else look only for other selectors
      std::copy(pl_in->begin(),pl_in->end(),back_inserter(*pl_out));
      m_ownlist=false;
    }
  }
  else {
    // else look only for other selectors
    std::copy(pl_in->begin(),pl_in->end(),back_inserter(*pl_out));
    m_ownlist=false;
  }

  // one particle select
  for (it=m_fmap.begin();it!=m_fmap.end();++it) Select(pl_out,it);  

  // two particle corr.
  Final_Correlator_Map::iterator ct;
  for (ct=m_cmap.begin();ct!=m_cmap.end();++ct) Select2(pl_out,ct);  

  // event conditions.
  for (it=m_fmap.begin();it!=m_fmap.end();++it) SelectN(pl_out,it);  

  // particle extraction
  Extract(pl_out);  

  if (!m_ownlist) {
    for (Particle_List::iterator itp=pl_out->begin(); itp!=pl_out->end();++itp) {
      *itp = new Particle(*itp);
    }
  }
  
  p_ana->AddParticleList(m_outlistname,pl_out);
}


Primitive_Observable_Base * Final_Selector::Copy() const 
{
  Final_Selector *fs = new Final_Selector(m_inlistname,m_outlistname,m_mode);
  for (Final_Data_Map::const_iterator it=m_fmap.begin();it!=m_fmap.end();++it) {
    fs->AddSelector(it->first,it->second);
  }

  for (Final_Data_Map::const_iterator it=m_fmap.begin();it!=m_fmap.end();++it) {
    if (m_extract && it->second.keep) fs->AddKeepFlavour(it->first);
  }

  for (Final_Correlator_Map::const_iterator ct=m_cmap.begin();ct!=m_cmap.end();++ct) {
    fs->AddSelector(ct->first.first,ct->first.second,ct->second);
  }
  return fs;
}    


Final_Selector::~Final_Selector() {
  if (p_jetalg) { delete p_jetalg; p_jetalg = NULL; }
}
