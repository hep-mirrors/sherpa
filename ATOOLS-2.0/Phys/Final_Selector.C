#include "Final_Selector.H"
#include "Primitive_Analysis.H"
#include "Kt_Algorithm.H"

#include <algorithm>

using namespace ATOOLS;

Final_Selector::Final_Selector(const std::string & ilist,const std::string & olist) :
  m_ilist(ilist),m_olist(olist),p_ktalg(NULL)
{
  m_splitt_flag=false;
  m_ownlist=false;
  m_extract=false;
  p_ktalg=new Kt_Algorithm();
}

void Final_Selector::AddSelector(const Flavour & fl, const Final_Selector_Data & fs) 
{
  Final_Data_Map::iterator it =m_fmap.find(fl);
  if (it==m_fmap.end()) {
    m_fmap[fl]=fs;
    if (m_extract) m_fmap[fl].keep=false;
  }
  else {
    it->second.eta_min = fs.eta_min; 
    it->second.eta_max = fs.eta_max;
    it->second.et_min  = fs.et_min;
    it->second.pt_min  = fs.pt_min;
    it->second.rmin    = fs.rmin;
  }
}

void Final_Selector::AddSelector(const Flavour & flav1, const Flavour & flav2, const Final_Selector_Data & fs) 
{
  std::pair<Flavour,Flavour> flavs(flav1,flav2);
  m_cmap[flavs]=fs;
}

void Final_Selector::AddSelector(const Flavour & fl, int min, int max) 
{
  Final_Data_Map::iterator it =m_fmap.find(fl);
  if (it!=m_fmap.end()) {
    it->second.min_n = min;  
    it->second.max_n = max;  
  }
  else {
    Final_Selector_Data fs;
    fs.min_n = min;  
    fs.max_n = max;  
    if (m_extract) fs.keep=false;
    m_fmap[fl]=fs;
  }
}

void Final_Selector::AddKeepFlavour(const Flavour & fl) 
{
  if (fl==Flavour(kf::lepton)) {
    for (int i=0;i<fl.Size();++i) AddKeepFlavour(fl[i]);
  }

  if(!m_extract) {
    Final_Data_Map::iterator it;
    for (it=m_fmap.begin();it!=m_fmap.end();++it) it->second.keep=false;
    m_extract=true;
  }
  m_fmap[fl].keep=true;
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
  double eta=mom.Eta();
  if (eta<etamin || etamax<eta ) return true;
  return false;
}

void Final_Selector::Select(Particle_List * jets,Final_Data_Map::iterator it) 
{
  for (Particle_List::iterator pit=jets->begin();pit!=jets->end();) {
    if ((*pit)->Flav()==it->first) {
      bool hit=false;
      if (it->second.eta_min!=it->second.eta_max) hit=EtaSelect((*pit)->Momentum(),it->second.eta_min,it->second.eta_max);
      if (it->second.et_min!=0. && !hit) hit=EtSelect((*pit)->Momentum(),it->second.et_min);
      if (it->second.pt_min!=0. && !hit) hit=PtSelect((*pit)->Momentum(),it->second.pt_min);
      if (!hit) ++pit;
      else {
	if (m_ownlist) delete *pit;
	pit=jets->erase(pit);
      }
    }
    else {
      ++pit;
    }
  }
}

void Final_Selector::Select2(Particle_List * jets,Final_Correlator_Map::iterator it) 
{
  if (it->second.rmin<=0.) return;

  Flavour flav1=it->first.first;
  Flavour flav2=it->first.second;

  bool hit=false;
  for (Particle_List::iterator pit=jets->begin();pit!=jets->end();++pit) {
    for (Particle_List::iterator pit2=jets->begin();pit2!=jets->end();++pit2) {
      if (flav1.Includes((*pit)->Flav()) && flav2.Includes((*pit2)->Flav()) && pit!=pit2) {
	if (DeltaR((*pit)->Momentum(),(*pit2)->Momentum()) < it->second.rmin && it->second.rmin>0.) {
	  hit=true;
	  break;
	}
      }
    }
  } 

  // delete list
  if (hit) {
    for (Particle_List::iterator pit=jets->begin();pit!=jets->end();) {
      if (m_ownlist) delete *pit;
      pit=jets->erase(pit);
    }
  }

}

void Final_Selector::SelectN(Particle_List * jets,Final_Data_Map::iterator it) 
{
  if (jets->size()==0) return;
  if (it->second.min_n==-1 && it->second.max_n==-1) return;

  int counter=0;
  for (Particle_List::iterator pit=jets->begin();pit!=jets->end();++pit) {
    if ((*pit)->Flav()==it->first) ++counter;
  }
  if ((it->second.min_n>counter && it->second.min_n!=-1) ||
      (it->second.max_n<counter && it->second.max_n!=-1)) {
    // delete list
    for (Particle_List::iterator pit=jets->begin();pit!=jets->end();) {
      if (m_ownlist) delete *pit;
      pit=jets->erase(pit);
    }   
  }
}

void Final_Selector::Extract(Particle_List * jets) 
{
  if (!m_extract) return;
  if (jets->size()==0) return;
  for (Particle_List::iterator pit=jets->begin();pit!=jets->end();) {
    bool remove=true;
    for (Final_Data_Map::iterator it=m_fmap.begin();it!=m_fmap.end();++it) {
      if ((*pit)->Flav()==it->first && it->second.keep) {
	remove=false;
	break;
      }
    }
    if (remove) {
      if (m_ownlist) delete *pit;
      pit=jets->erase(pit);
    }
    else ++pit;
  }

}


double Final_Selector::DeltaR(const Vec4D & p1,const Vec4D & p2) 
{
  /*
  double c1=p1[3]/Vec3D(p1).Abs();
  double c2=p2[3]/Vec3D(p2).Abs();
  double deta12 = 0.5 *log( (1 + c1)*(1 - c2)/((1-c1)*(1+c2)));
  */
  double deta12 = p1.Eta() - p2.Eta();

  double pt1=sqrt(p1[1]*p1[1]+p1[2]*p1[2]);
  double pt2=sqrt(p2[1]*p2[1]+p2[2]*p2[2]);
  double dphi12=acos((p1[1]*p2[1]+p1[2]*p2[2])/(pt1*pt2));
  
  return sqrt(sqr(deta12) + sqr(dphi12));
}

void Final_Selector::Evaluate(const Blob_List &,double value, int ncount) {
  Particle_List * pl=p_ana->GetParticleList(m_ilist);
  if (pl==NULL) {
    msg.Error()<<" WARNING: particle list "<<m_ilist<<" not found "<<std::endl;
    return;
  }

  Particle_List * jets= new Particle_List;
  
  
  // look for kt and after for other selectors
  Final_Data_Map::iterator it =m_fmap.find(Flavour(kf::jet));
  if (it!=m_fmap.end()) {
    if (it->second.rmin>0.) {
      p_ktalg->ConstructJets(pl,jets,0,it->second.rmin);
      // add leptons
      for (Particle_List::iterator pit=pl->begin();pit!=pl->end();++pit) {
	if ((*pit)->Flav().IsLepton()) jets->push_back(new Particle(*pit));
      }
      m_ownlist=true;
    }
    else {
      // else look only for other selectors
      std::copy(pl->begin(),pl->end(),back_inserter(*jets));
      m_ownlist=false;
    }
  }
  else {
    // else look only for other selectors
    std::copy(pl->begin(),pl->end(),back_inserter(*jets));
    m_ownlist=false;
  }

  // one particle select
  for (it=m_fmap.begin();it!=m_fmap.end();++it) Select(jets,it);  

  // two particle corr.
  Final_Correlator_Map::iterator ct;
  for (ct=m_cmap.begin();ct!=m_cmap.end();++ct) Select2(jets,ct);  

  // event conditions.
  for (it=m_fmap.begin();it!=m_fmap.end();++it) SelectN(jets,it);  

  // particle extraction
  Extract(jets);  

  if (!m_ownlist) {
    for (Particle_List::iterator itp=jets->begin(); itp!=jets->end();++itp) {
      *itp = new Particle(*itp);
    }
  }
  
  p_ana->AddParticleList(m_olist,jets);
}


Primitive_Observable_Base * Final_Selector::Copy() const 
{
  Final_Selector *fs = new Final_Selector(m_ilist,m_olist);
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
  if (p_ktalg) delete p_ktalg;

}
