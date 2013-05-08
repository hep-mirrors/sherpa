#include "AddOns/Analysis/Triggers/Final_Selector.H"

using namespace ANALYSIS;

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "AddOns/Analysis/Tools/Particle_Qualifier.H"
#include "AddOns/Analysis/Triggers/Midpoint_Cone.H"
#include "AddOns/Analysis/Triggers/DIS_Algorithm.H"
#include "AddOns/Analysis/Triggers/MySISCone.H"
#include "AddOns/Analysis/Triggers/MCFMCone.H"
#include <iomanip>

DECLARE_GETTER(Final_Selector,"Trigger",
 	       Analysis_Object,Argument_Matrix);

void ATOOLS::Getter
<Analysis_Object,Argument_Matrix,Final_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"{\n"
     <<std::setw(width+7)<<" "<<"InList  list\n"
     <<std::setw(width+7)<<" "<<"OutList list\n"
     <<std::setw(width+7)<<" "<<"JetMode mode\n"
     <<std::setw(width+7)<<" "<<"Finder  kf [ptmin etamin etamax rmin bjets]\n"
     <<std::setw(width+7)<<" "<<"DRMin   kf1 kf2 drmin\n"
     <<std::setw(width+7)<<" "<<"Counts  kf min max\n"
     <<std::setw(width+7)<<" "<<"Keep    kf\n"
     <<std::setw(width+7)<<" "<<"Qual    qualifier\n"
     <<std::setw(width+4)<<" "<<"}";
}

Analysis_Object * ATOOLS::Getter
<Analysis_Object,Argument_Matrix,Final_Selector>::
operator()(const Argument_Matrix &parameters) const
{
  Final_Selector_Data data;
  int jetmode=0;
  if (ATOOLS::rpa->gen.Beam1().Kfcode()==kf_e || 
      ATOOLS::rpa->gen.Beam2().Kfcode()==kf_e) jetmode=3;
  if (ATOOLS::rpa->gen.Beam1().Kfcode()==kf_e && 
      ATOOLS::rpa->gen.Beam2().Kfcode()==kf_e) jetmode=1;
  std::string inlist="FinalState", outlist="Analysed";
  ATOOLS::Particle_Qualifier_Base *qualifier=NULL;
  for (size_t i=0;i<parameters.size();++i) {
    const std::vector<std::string> &cur=parameters[i];
    if (cur[0]=="InList" && cur.size()>1) inlist=cur[1];
    else if (cur[0]=="OutList" && cur.size()>1) outlist=cur[1];
    else if (cur[0]=="JetMode" && cur.size()>1) jetmode=ATOOLS::ToType<int>(cur[1]);
    else if (cur[0]=="Qual" && cur.size()>1) {
      if (!qualifier) delete qualifier;
      qualifier = ATOOLS::Particle_Qualifier_Getter::GetObject(cur[1],cur[1]);
    }
  }
  if (!qualifier) qualifier = new ATOOLS::Is_Not_Lepton(); 
  Final_Selector *selector = new Final_Selector(inlist,outlist,jetmode,qualifier);
  selector->SetAnalysis(parameters());
  for (size_t i=0;i<parameters.size();++i) {
    const std::vector<std::string> &cur=parameters[i];
    if (cur[0]=="Finder" && cur.size()>1) {
      int kf=ATOOLS::ToType<int>(cur[1]);
      // if (kf!=kf_jet || kf==kf_bjet) continue;
      ATOOLS::Flavour flavour((kf_code)abs(kf));
      if (kf<0) flavour=flavour.Bar();
      data.pt_min=0.;
      data.eta_min=-20.;
      data.eta_max=+20.;
      data.f=0.5;
      if (cur.size()>2) data.pt_min=ATOOLS::ToType<double>(cur[2]);
      if (data.pt_min<0.) { data.et_min=-data.pt_min;data.pt_min=0.; }
      if (cur.size()>3) data.eta_min=ATOOLS::ToType<double>(cur[3]);
      if (cur.size()>4) data.eta_max=ATOOLS::ToType<double>(cur[4]);
      if (cur.size()>5 && (kf==93 || kf==97)) data.r_min=ATOOLS::ToType<double>(cur[5]);
      if (cur.size()>6 && (kf==93 || kf==97)) data.bf=ATOOLS::ToType<int>(cur[6]);
      if (cur.size()>7 && (kf==93 || kf==97)) data.f=ATOOLS::ToType<double>(cur[7]);
      selector->AddSelector(flavour,data);
    }
    else if (cur[0]=="DRMin" && cur.size()>3) {
      int kf=ATOOLS::ToType<int>(cur[1]);
      ATOOLS::Flavour f1((kf_code)abs(kf));
      if (kf<0) f1=f1.Bar();
      kf=ATOOLS::ToType<int>(cur[2]);
      ATOOLS::Flavour f2((kf_code)abs(kf));
      if (kf<0) f2=f2.Bar();
      data.r_min=0.;
      if (cur.size()>3) data.r_min=ATOOLS::ToType<double>(cur[3]);
      selector->AddSelector(f1,f2,data);
    }
    else if (cur[0]=="Counts" && cur.size()>3) {
      int kf=ATOOLS::ToType<int>(cur[1]);
      ATOOLS::Flavour flavour((kf_code)abs(kf));
      if (kf<0) flavour=flavour.Bar();
      selector->AddSelector(flavour,ATOOLS::ToType<int>(cur[2]),
			    ATOOLS::ToType<int>(cur[3]));
    }
    else if (cur[0]=="Keep" && cur.size()>1) {
      int kf=ATOOLS::ToType<int>(cur[1]);
      ATOOLS::Flavour flavour((kf_code)abs(kf));
      if (kf<0) flavour=flavour.Bar();
      selector->AddKeepFlavour(flavour);
    }
  }
  return selector;
}

DECLARE_GETTER(Leading_Particle,"LeadingParticle",
 	       Analysis_Object,Argument_Matrix);

void ATOOLS::Getter
<Analysis_Object,Argument_Matrix,Leading_Particle>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"{\n"
     <<std::setw(width+7)<<" "<<"InList  list\n"
     <<std::setw(width+7)<<" "<<"OutList list\n"
     <<std::setw(width+7)<<" "<<"Keep    kf\n"
     <<std::setw(width+7)<<" "<<"Qual    qualifier\n"
     <<std::setw(width+4)<<" "<<"}";
}

Analysis_Object * ATOOLS::Getter
<Analysis_Object,Argument_Matrix,Leading_Particle>::
operator()(const Argument_Matrix &parameters) const
{
  std::string inlist("FinalState"), outlist("Analysed");
  int mode(0);
  ATOOLS::Particle_Qualifier_Base *qualifier(NULL);
  for (size_t i=0;i<parameters.size();++i) {
    const std::vector<std::string> &cur=parameters[i];
    if      (cur[0]=="InList"  && cur.size()>1) inlist=cur[1];
    else if (cur[0]=="OutList" && cur.size()>1) outlist=cur[1];
    else if (cur[0]=="Mode"    && cur.size()>1) {
      if (cur[1]=="PT") mode = 1;
    }
    else if (cur[0]=="Qual"    && cur.size()>1) {
      if (!qualifier) delete qualifier;
      qualifier = ATOOLS::Particle_Qualifier_Getter::GetObject(cur[1],cur[1]);
    }
  }
  if (!qualifier) qualifier = new ATOOLS::Is_Hadron(); 
  Leading_Particle *selector = new Leading_Particle(inlist,outlist,mode,qualifier);
  selector->SetAnalysis(parameters());
  return selector;
}

//##############################################################################
//##############################################################################
//##############################################################################




#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Org/Message.H"
#include "AddOns/Analysis/Triggers/Durham_Algorithm.H"
#include "AddOns/Analysis/Triggers/Calorimeter_Cone.H"

#include <algorithm>

using namespace ATOOLS;

namespace ANALYSIS {
  std::ostream & operator<<(std::ostream & s,const Final_Selector_Data & fd )
  {
    s<<"[keep("<<fd.keep<<","<<fd.bf<<"), n("<<fd.min_n<<","<<fd.max_n<<"),"
     <<" eta("<<fd.eta_min<<","<<fd.eta_max<<"), "<<"etpt("<<fd.et_min<<","<<fd.pt_min<<"), "
     <<" r_min("<<fd.r_min<<"), "<<" mass("<<fd.mass_min<<","<<fd.mass_max<<")]";
    return s;
  }
}

Final_Selector::Final_Selector(const std::string & inlistname,
			       const std::string & outlistname) :
  m_inlistname(inlistname),m_outlistname(outlistname),m_ownlist(false), m_extract(false),
  m_mode(-1), p_jetalg(NULL)
{
  m_name="Trigger";
}

Final_Selector::Final_Selector(const std::string & inlistname,
			       const std::string & outlistname,
			       int mode, SP(ATOOLS::Particle_Qualifier_Base) qualifier) :
  p_qualifier(qualifier), m_inlistname(inlistname), m_outlistname(outlistname),
  m_ownlist(false), m_extract(false), m_mode(mode), p_jetalg(NULL)
{
  msg_Tracking()<<" init Final_Selector("<<inlistname<<","<<outlistname<<","
		<<mode<<","<<qualifier<<")"<<std::endl;
  m_name="Trigger";
  switch (mode) {
  case 1: p_jetalg = new Durham_Algorithm(--p_qualifier); break;
  case 0: p_jetalg = new Kt_Algorithm(--p_qualifier); break;
  case 3: p_jetalg = new DIS_Algorithm(--p_qualifier); break;
  default:
    break;
  }
}

void Final_Selector::AddSelector(const Flavour & fl, const Final_Selector_Data & fs) 
{
  msg_Tracking()<<" AddSelector("<<fl<<","<<fs<<")"<<std::endl;
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
    it->second.bf      = fs.bf;
    it->second.f       = fs.f;
  }
  
  if (fl==kf_jet || fl==kf_bjet) {
    switch(m_mode) {
    case 2: p_jetalg = new 
	      Calorimeter_Cone(fs.pt_min,fs.eta_min,fs.eta_max);break;
    case 10: p_jetalg = new Midpoint_Cone(--p_qualifier,0,fs.f); break;
    case 11: p_jetalg = new Midpoint_Cone(--p_qualifier,1,fs.f); break;
    case 20: p_jetalg = new SISCone(--p_qualifier,fs.f); break;
    case 30: p_jetalg = new MCFMCone(--p_qualifier,fs.f); break;
    case 40: p_jetalg = new Kt_Algorithm(--p_qualifier); break;
    }
    if (p_jetalg) p_jetalg->Setbflag(fs.bf);
  }
}

void Final_Selector::AddSelector(const Flavour & flav1, const Flavour & flav2, 
				 const Final_Selector_Data & fs) 
{
  msg_Tracking()<<" AddSelector("<<flav1<<","<<flav2<<","<<fs<<")"<<std::endl;
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
  msg_Tracking()<<" AddSelector("<<fl<<", n("<<min<<","<<max<<") )"<<std::endl;
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
    it->second.ko    = false;
  }
}

void Final_Selector::AddSelector(const Flavour & fl, const Final_Selector_Data & fs,
				 Calorimeter_Cone * const cone) {
  msg_Tracking()<<" AddSelector : Cone."<<std::endl;
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
    it->second.ko      = false;
  }
  if (p_jetalg!=NULL) {
    msg_Error()<<"Error in Final_Selector::AddSelector("<<cone<<") : "<<std::endl
	       <<"   Tried to add a cone finder based on Hcal,"
	       <<" jet finder already present."<<std::endl
	       <<"   Abort the run."<<std::endl;
    abort();
  }
  p_jetalg = cone;
}

void Final_Selector::AddKeepFlavour(const Flavour & fl) 
{
  msg_Tracking()<<" AddKeepFlavour("<<fl<<")"<<std::endl;
  if (fl==Flavour(kf_lepton)) {
    for (size_t i=0;i<fl.Size();++i) AddKeepFlavour(fl[i]);
  }

  if(!m_extract) {
    Final_Data_Map::iterator it;
    for (it=m_fmap.begin();it!=m_fmap.end();++it) it->second.keep = false;
    m_extract = true;
  }
  if (m_fmap.find(fl)==m_fmap.end()) m_fmap[fl].ko=true;
  m_fmap[fl].keep = true;
}

void Final_Selector::Output()
{
  if (!msg_LevelIsTracking()) return;
  msg_Out()<<"Final_Selector : "<<m_fmap.size()<<"/"<<m_cmap.size()<<":"<<std::endl;
  for (Final_Data_Map::iterator it=m_fmap.begin();it!=m_fmap.end();++it) {
    if (it->first!=Flavour(kf_jet) && it->first!=Flavour(kf_bjet)) 
      msg_Out()<<" "<<it->first<<" : pt_min = "<<it->second.pt_min<<", eta = "
	       <<it->second.eta_min<<" ... "<<it->second.eta_max<<std::endl;
    else
      msg_Out()<<" "<<it->first<<" : pt_min = "<<it->second.pt_min<<", eta = "
	       <<it->second.eta_min<<" ... "<<it->second.eta_max
	       <<", jets with ktRunII, r_min = "<<it->second.r_min<<std::endl;
  }
  for (Final_Correlator_Map::iterator it=m_cmap.begin();it!=m_cmap.end();++it) {
    msg_Out()<<" "<<it->first.first<<" "<<it->first.second<<" : "<<it->second.r_min<<std::endl;
  }
  for (Final_Data_Map::iterator it=m_fmap.begin();it!=m_fmap.end();++it) {
    if ((it->second.min_n>-1) || (it->second.max_n>-1)) {
      msg_Out()<<" "<<it->first<<" : min = "<<it->second.min_n<<", max = "<<it->second.max_n<<std::endl;
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

void Final_Selector::JetSelect(Particle_List * pl,const Flavour& jf) 
{
  for (Particle_List::iterator pit=pl->begin();pit!=pl->end();) {
    if ((*pit)->Flav()!=jf) {
      if (m_ownlist) delete *pit;
      pit = pl->erase(pit);
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
	if (hit) break;
      }
    }
    if (hit) break;
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

void Final_Selector::Evaluate(const Blob_List &bl,double value, double ncount) {
  Particle_List * pl_in = p_ana->GetParticleList(m_inlistname);
  if (pl_in==NULL) {
    msg_Out()<<"WARNING in Final_Selector::Evaluate : particle list "<<m_inlistname<<" not found "<<std::endl;
    return;
  }
  Particle_List * pl_out = new Particle_List;
  // look for kt and after for other selectors
  Final_Data_Map::iterator it = m_fmap.find(Flavour(kf_jet));
  if (it==m_fmap.end() || it->second.r_min==0.0) it = m_fmap.find(Flavour(kf_bjet));
  if (it!=m_fmap.end()) {
    if (it->second.r_min>0.) {
      std::vector<double> * diffrates=new std::vector<double>();
      p_jetalg->SetBlobList(&bl);
      p_jetalg->ConstructJets(pl_in,pl_out,diffrates,it->second.r_min);
      // JetSelect(pl_out,it->first);
      // add leptons
      for (Particle_List::iterator pit=pl_in->begin();pit!=pl_in->end();++pit) {
	if (!(*p_qualifier)(*pit))  pl_out->push_back(new Particle(**pit));
      }
      m_ownlist=true;
      std::string key;
      /*
      MyStrStream str;
      str<<"KtJetrates("<<it->second.r_min<<")"<<m_outlistname;
      str>>key;
      */
      key="KtJetrates(1)"+m_outlistname;
      p_ana->AddData(key,new Blob_Data<std::vector<double> *>(diffrates));

      Blob_Data_Base * ktdrs=(*p_ana)["KtDeltaRs"];
      if (ktdrs) {
	ktdrs->Get<std::vector<double> *>()->push_back(it->second.r_min);
      }
      else {
	std::vector<double> * drs = new std::vector<double>;
	drs->push_back(it->second.r_min);
	p_ana->AddData("KtDeltaRs",new Blob_Data<std::vector<double> *>(drs));
      }

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
      *itp = new Particle(**itp);
    }
  }
  //  std::sort(pl_out->begin(),pl_out->end(),ATOOLS::Order_PT());
  p_ana->AddParticleList(m_outlistname,pl_out);
}

void Final_Selector::SetAnalysis(Primitive_Analysis  * ana)
{
  p_ana=ana;
  if (p_jetalg!=NULL) {
    Calorimeter_Cone *cc(dynamic_cast<Calorimeter_Cone*>(p_jetalg));
    if (cc!=NULL) cc->SetAnalysis(p_ana);
  }
}


Analysis_Object * Final_Selector::GetCopy() const 
{
  Final_Selector *fs = new Final_Selector(m_inlistname,m_outlistname,m_mode,p_qualifier);
  fs->SetAnalysis(p_ana);
  for (Final_Data_Map::const_iterator it=m_fmap.begin();it!=m_fmap.end();++it) {
    if (!it->second.ko) fs->AddSelector(it->first,it->second);
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
  if (p_jetalg) delete p_jetalg;
}



//##############################################################################
//##############################################################################
//##############################################################################

Leading_Particle::Leading_Particle(const std::string & inlistname,
				   const std::string & outlistname) :
  p_qualifier(NULL), m_inlistname(inlistname), m_outlistname(outlistname), m_mode(0)
{
  m_name="Leading_Particle(E)";
}

Leading_Particle::Leading_Particle(const std::string & inlistname,
				   const std::string & outlistname,
				   int mode, ATOOLS::Particle_Qualifier_Base * const qualifier) :
  p_qualifier(qualifier), m_inlistname(inlistname), m_outlistname(outlistname), m_mode(mode)
{
  msg_Out()<<" Init Leading_Particle("<<inlistname<<","<<outlistname<<","
	   <<mode<<","<<qualifier<<")"<<std::endl;
  m_name="Leading_Particle";
  switch (mode) {
  case 1:  m_name+="(PT)"; break;
  default: m_name+="(E)"; break;
  }
}

Leading_Particle::~Leading_Particle() {}

Analysis_Object * Leading_Particle::GetCopy() const 
{
  Leading_Particle * lp = new Leading_Particle(m_inlistname,m_outlistname,m_mode,p_qualifier);
  lp->SetAnalysis(p_ana);
  return lp;
}    

void Leading_Particle::Evaluate(const Blob_List &,double value, double ncount) {
  Particle_List * pl_in = p_ana->GetParticleList(m_inlistname);
  if (pl_in==NULL) {
    msg_Out()<<"WARNING in Leading_Particle::Evaluate : particle list "
	     <<m_inlistname<<" not found "<<std::endl;
    return;
  }
  Particle * winner(NULL);
  double crit(0.),test(0.);
  //std::cout<<"-----------------------------------------------------------"<<std::endl;
  for (Particle_List::iterator pit=pl_in->begin();pit!=pl_in->end();++pit) {
    if ((*p_qualifier)(*pit)) {
      //std::cout<<"Test "<<(*pit)->Flav()<<" successful."<<std::endl;
      switch (m_mode) {
      case 1:  test =(*pit)->Momentum().PPerp2();
      default: test =(*pit)->Momentum()[0];
      }
      if (test>crit) { winner=(*pit); test=crit; }
    }
  }    

  Particle_List * pl_out = new Particle_List;
  if (winner!=NULL) pl_out->push_back(new Particle(*winner));
  p_ana->AddParticleList(m_outlistname,pl_out);
}

void Leading_Particle::Output()
{
  if (!msg_LevelIsTracking()) return;
  msg_Out()<<m_name<<"."<<std::endl;
}

