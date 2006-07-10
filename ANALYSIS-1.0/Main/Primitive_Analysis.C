#include "Primitive_Analysis.H"
#include "Primitive_Observable_Base.H"
#include "Particle_Selector.H"
#include "Universal_Selector.H"
#include "Message.H"
#include "MyStrStream.H"
#include "Shell_Tools.H"
#include "Particle_Qualifier.H"

#ifdef PROFILE__Analysis_Phase
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif

using namespace ANALYSIS;
using namespace ATOOLS;

Primitive_Analysis::Primitive_Analysis(const std::string _name, const int mode) :
  p_bfinder(0) , m_active(true)
{
  m_nevt = 0;
  p_partner = this;
  m_mode = mode;

  m_name = std::string("Analysis : ") + _name;
  msg_Tracking()<<" Initializing Primitive_Analysis : "<<m_name<<std::endl;
}

Primitive_Analysis::Primitive_Analysis(const int mode) :
  m_nevt(0), p_partner(this), p_bfinder(0), m_active(true)
{
  m_mode = mode;

  m_name = std::string("Analysis : noname");
  msg_Tracking()<<" Initializing Primitive_Analysis : "<<m_name<<std::endl;
}

Primitive_Analysis::~Primitive_Analysis()
{
  for(int i=m_observables.size();i>0;i--) {
    if (m_observables[i-1]) delete m_observables[i-1];
  }
  m_observables.clear();

  for (Analysis_List::iterator it=m_subanalyses.begin();it!=m_subanalyses.end();++it) 
    delete it->second;
  m_subanalyses.clear();
  if (p_bfinder) delete p_bfinder;
}

void Primitive_Analysis::AddObservable(Primitive_Observable_Base * obs) 
{
  obs->SetAnalysis(p_partner);
  std::string oname=obs->Name();
  std::string id="_A";
  size_t pos=oname.find(".dat");
  for(size_t i=0;i<m_observables.size();++i) {
    if (m_observables[i]->Name()==obs->Name()) {
      std::string pname;
      if (pos!=std::string::npos)
	pname=oname.substr(0,pos)+id+oname.substr(pos);
      else 
	pname=oname+id;
      obs->SetName(pname);
      ++id[1];
    }
  }
  m_observables.push_back(obs);
}

void Primitive_Analysis::AddSubAnalysis(const std::string & key,Primitive_Analysis * ana)
{
  Analysis_List::const_iterator cit=m_subanalyses.find(key);
  if (cit!=m_subanalyses.end()) {
    msg.Out()<<"WARNING in Primitive_Analysis::AddSubAnalysis :"<<std::endl
	     <<"   Analysis "<<key<<" already existent;"<<std::endl
	     <<" sub analysis not added, will be deleted."<<std::endl;
    if (ana) delete ana;
    return;
  }

  m_subanalyses[key]=ana;  
}

Primitive_Analysis * Primitive_Analysis::GetSubAnalysis(const std::string & key, int mode) 
{
  Analysis_List::const_iterator cit=m_subanalyses.find(key);
  if (cit!=m_subanalyses.end()) return cit->second;

  bool master=true;
  if (key=="ME" || key=="MI" || key=="Shower" || key=="Hadron") {
    master=false;
    if (key!="ME"     && mode&ANALYSIS::do_me) mode=mode^ANALYSIS::do_me;
    if (key!="MI"     && mode&ANALYSIS::do_mi) mode=mode^ANALYSIS::do_mi;
    if (key!="Shower" && mode&ANALYSIS::do_shower) mode=mode^ANALYSIS::do_shower;
    if (key!="Hadron" && mode&ANALYSIS::do_hadron) mode=mode^ANALYSIS::do_hadron;
  }

  Primitive_Analysis * ana = new Primitive_Analysis(m_name.substr(11)+key,mode);
  if (master) ana->SetPartner(p_partner);

  for (size_t i=0;i<m_observables.size();i++) {
    if (m_observables[i]->Splittable() || !master) 
      ana->AddObservable(m_observables[i]->Copy());
  }
  m_subanalyses[key]=ana;
  return ana;
}

void Primitive_Analysis::CallSubAnalysis(const Blob_List * const bl, double value) 
{
  int nout=-1;
  std::string name;
  for (Blob_List::const_iterator bit=bl->begin();bit!=bl->end();++bit) {
    if ((*bit)->Type()==btp::Signal_Process) {
      nout  = (*bit)->NOutP();      
      name  = (*bit)->TypeSpec();    //orig: (*bit)->Type();
      break;
    }
  }
  if (nout==-1) {
    msg.Out()<<"WARNING in Primitive_Analysis::CallSubAnalysis: no Signal process found "<<std::endl;
    return;
  }

  Blob_Data_Base * extra_info=operator[]("OrderEWeak");

  std::string key;
  int mode;
  if (extra_info && m_mode&ANALYSIS::splitt_extra) {
    mode=m_mode^ANALYSIS::splitt_extra;
    if (m_mode&ANALYSIS::splitt_jetseeds)
      mode=m_mode^ANALYSIS::splitt_jetseeds;
    mode=mode|ANALYSIS::output_this;

    double dkey=extra_info->Get<double>();
    key="extra_"+ATOOLS::ToString(dkey);

    GetSubAnalysis(key,mode);

    for (Analysis_List::iterator it=m_subanalyses.begin();it!=m_subanalyses.end();++it) {
      if (it->first.find("extra_")!=std::string::npos) {
	if (it->first==key) {
	  it->second->DoAnalysis(bl,value);
	}
	else {
	  m_active=false;
	  it->second->DoAnalysis(bl,value);
	  m_active=true;
	}
      }
    }
  }
  if (m_mode&ANALYSIS::splitt_jetseeds) {
    mode=m_mode^ANALYSIS::splitt_jetseeds;
    if (mode&ANALYSIS::splitt_extra)
      mode=m_mode^ANALYSIS::splitt_extra;
    mode=mode|ANALYSIS::output_this;
    key="j"+ToString(nout);
  }
  else {
    mode=m_mode^ANALYSIS::splitt_process;
//     if (m_mode&ANALYSIS::output_process) mode=mode|ANALYSIS::output_this;
//     else 
    if (m_mode&ANALYSIS::output_this) mode=mode^ANALYSIS::output_this;
      key=name;
  }
  
  Primitive_Analysis * ana=GetSubAnalysis(key,mode);
  ana->DoAnalysis(bl,value);
}


void Primitive_Analysis::DoAnalysis(const Blob_List * const bl, const double value) {
  ++m_nevt;

  if (m_mode&ANALYSIS::splitt_phase) {
    m_mode=m_mode|ANALYSIS::output_this;
    int mode=m_mode^ANALYSIS::splitt_phase;
    if (m_mode&ANALYSIS::do_me)     GetSubAnalysis("ME",mode)->DoAnalysis(bl,value);
    if (m_mode&ANALYSIS::do_mi)     GetSubAnalysis("MI",mode)->DoAnalysis(bl,value);
    if (m_mode&ANALYSIS::do_shower) GetSubAnalysis("Shower",mode)->DoAnalysis(bl,value);
    if (m_mode&ANALYSIS::do_hadron) GetSubAnalysis("Hadron",mode)->DoAnalysis(bl,value);
    return;
  }

  ClearAllData();
  p_blobs = bl;

  // if (value!=1.) m_mode=m_mode|ANALYSIS::weighted;

  if (p_partner==this) {
    m_mode=m_mode|ANALYSIS::fill_helper;
    m_mode=m_mode|ANALYSIS::output_this;
  }
  else if (m_mode&ANALYSIS::fill_helper) m_mode=m_mode^ANALYSIS::fill_helper;
  if (p_partner==this &&  m_mode&ANALYSIS::weighted) {
    m_mode= m_mode|ANALYSIS::splitt_process;
  }
  else if (m_mode&ANALYSIS::splitt_process) {
    m_mode=m_mode|ANALYSIS::output_process;
  }
  if ((m_mode&ANALYSIS::splitt_all)==0) m_mode=m_mode|ANALYSIS::fill_histos;
  if (m_mode&ANALYSIS::weighted && m_mode&ANALYSIS::splitt_all) {
    m_mode=m_mode^(m_mode&ANALYSIS::fill_histos);
  }
  Init();
  double weight=(*p_partner)["ME_Weight"]->Get<double>();
  double procweight=1.;
  if (m_mode&ANALYSIS::weighted_ns || 
      !(m_mode&ANALYSIS::weighted)) 
    procweight=(*p_partner)["Process_Weight"]->Get<double>();
  weight/=procweight;
  int    ncount=(*p_partner)["ME_NumberOfTrials"]->Get<int>();
  if (!IsEqual(value/procweight,weight)) {
    if (p_partner==this) {
      msg.Out()<<"WARNING in Primitive_Analysis::DoAnalysis :"<<std::endl
	       <<"   Weight in Primitive_Analysis ambiguous! ("<<value/procweight<<","<<weight<<")"<<std::endl;
    }
    else if (value/procweight==0.) {
      weight=0.;
    }
    else {
      msg.Out()<<"WARNING something is wrong in Primitive_Analysis::DoAnalysis :"<<std::endl
	       <<"   Weight in Primitive_Analysis ambiguous! ("<<value/procweight<<","<<weight<<")"<<std::endl;
    }
  }
  double weight_one=weight;
  int    ncount_one=ncount;
  Blob_Data_Base * info = (*p_partner)["ME_Weight_One"];
  if (info) {
    weight_one = info->Get<double>();
    weight_one/= procweight;
    ncount_one = (*p_partner)["ME_NumberOfTrials_One"]->Get<int>();
  }
  if (weight==0.) weight_one=0.;
  m_stats.sum_weight     += weight;
  m_stats.nevt           += ncount;
  m_stats.sum_weight_one += weight_one;
  m_stats.nevt_one       += ncount_one;
  
  // do nonsplittable (helper and legacy observables) first
  if (m_mode&ANALYSIS::fill_helper) {
    for (size_t i=0;i<m_observables.size();i++) {
      if (!m_observables[i]->Splittable()) {
	m_observables[i]->Evaluate(*bl,value/procweight,ncount);
      }
    }
  }

  if (m_mode&ANALYSIS::fill_histos) {
    for (size_t i=0;i<m_observables.size();i++) {
      if (m_observables[i]->Splittable()) {
	m_observables[i]->Evaluate(*bl,value/procweight,ncount);
      }
    }
  }


  if (m_mode&ANALYSIS::splitt_all) CallSubAnalysis(bl,value);
  if (p_partner==this && msg.LevelIsTracking()) PrintStatus();

  ClearAllData();
}

void Primitive_Analysis::FinishAnalysis(const std::string & resdir,long ntotal, double xs) 
{
  if (ntotal==0) ntotal=m_nevt;
  if (m_mode&ANALYSIS::output_this) 
    ATOOLS::MakeDir(resdir+OutputPath(),448); 

  for (Analysis_List::iterator it=m_subanalyses.begin();
       it!=m_subanalyses.end();++it) {
    std::string dir=resdir+OutputPath()+std::string("/")+it->first;
    //std::cout<<"Subanalysis: ";
    it->second->FinishAnalysis(dir,ntotal,xs);
  }

  if (!(m_mode&ANALYSIS::splitt_phase)) {
    for (size_t i=0;i<m_observables.size();i++) {
      if (m_mode&ANALYSIS::weighted  && m_mode&ANALYSIS::splitt_all && m_observables[i]->Splittable()) {
	std::string key = m_observables[i]->Name();
	m_observables[i]->Reset();  // just in case someone has filled something in already
      
	for (Analysis_List::iterator it=m_subanalyses.begin();
	     it!=m_subanalyses.end();++it) {
	  if (it->first.find("extra_")==std::string::npos) {
	    Primitive_Observable_Base * ob = it->second->GetObservable(key);
	    if (ob)   (*m_observables[i])+=(*ob);
	  }
	}
      }
      else {
	if (m_mode&ANALYSIS::weighted_ns) {
	  m_observables[i]->EndEvaluation(double(m_stats.nevt)/double(ntotal)*xs);
	}
	else {
	  if ((m_mode&ANALYSIS::weighted)==0 ) {
	    m_observables[i]->EndEvaluation(double(m_nevt)/double(ntotal)*xs);
	  }
	  else {
	    if (m_stats.nevt_one>0 && m_stats.sum_weight!=0.) {
	      double xshist=m_stats.sum_weight/double(m_stats.nevt);
	      double xsreal=m_stats.sum_weight_one/double(m_stats.nevt_one);
	      m_observables[i]->EndEvaluation(xsreal/xshist);
	    }
	    else {
	      m_observables[i]->EndEvaluation();
	    }
	  }
	}
      }
      if (m_mode&ANALYSIS::output_this)
	m_observables[i]->Output(resdir+OutputPath());
    }
  }
}

void Primitive_Analysis::Init()
{
  if (m_mode&ANALYSIS::fill_helper)
    CreateFinalStateParticleList();
}

bool Primitive_Analysis::SelectBlob(const ATOOLS::Blob *blob) 
{
  if (m_mode&ANALYSIS::do_hadron) return true;
  if (m_mode&ANALYSIS::do_shower && 
      (blob->Type()==btp::IS_Shower || blob->Type()==btp::FS_Shower ||
       blob->Type()==btp::Shower)) return true;
  if (m_mode&ANALYSIS::do_mi && 
      (blob->Type()==btp::Hard_Collision ||
       blob->Type()==btp::Signal_Process)) return true;
  if (m_mode&ANALYSIS::do_me && blob->Type()==btp::Signal_Process) return true;
  return false;
}

void Primitive_Analysis::CreateFinalStateParticleList(bool markb)
{
  std::string key="FinalState";
  if (markb) {
    key="FinalStateB";
    if (!p_bfinder) p_bfinder = ATOOLS::Particle_Qualifier_Getter::GetObject("DecayedBHadron","DecayedBHadron");
  }

  PL_Container::const_iterator cit=m_pls.find(key);
  if (cit!=m_pls.end()) return;

  Particle_List * pl = new Particle_List;

  for (Blob_List::const_iterator blit=p_blobs->begin();blit!=p_blobs->end();++blit) {
    if (!markb) {
      if ((*blit)->Type()==btp::Signal_Process) {
	Blob_Data_Base * info=(*(*blit))["ME_Weight"];
	if (info) {
	  m_datacontainer["ME_Weight"]=new Blob_Data<double>(info->Get<double>());
	  info=(*(*blit))["Process_Weight"];
	  if (info) {
	    m_datacontainer["Process_Weight"]=new Blob_Data<double>(info->Get<double>());
	  }
	  info=(*(*blit))["ME_NumberOfTrials"];
	  if (info) {
	    m_datacontainer["ME_NumberOfTrials"]=new Blob_Data<int>(info->Get<int>());
	  }
	}
	info=(*(*blit))["ME_Weight_One"];
	if (info) {
	  m_datacontainer["ME_Weight_One"]=new Blob_Data<double>(info->Get<double>());
	  info=(*(*blit))["ME_NumberOfTrials_One"];
	  if (info) {
	    m_datacontainer["ME_NumberOfTrials_One"]=new Blob_Data<int>(info->Get<int>());
	  }
	}
      }
      if ((*blit)->Type()==ATOOLS::btp::ME_PS_Interface_FS) {
	Blob_Data_Base * info=(*(*blit))["OrderStrong"];
	if (info &&m_datacontainer.find("OrderStrong")==m_datacontainer.end()) {
	  m_datacontainer["OrderStrong"]=new Blob_Data<double>(info->Get<double>());
	  info=(*(*blit))["OrderEWeak"];
	  if (info) {
	    m_datacontainer["OrderEWeak"]=new Blob_Data<double>(info->Get<double>());
	  }
	}
      }
    }
    for (String_BlobDataBase_Map::const_iterator it=(*blit)->GetData().begin();
	   it!=(*blit)->GetData().end(); ++it) {
      if (it->first.length()>2 && it->first[0]=='d' && it->first[1]=='#') {
	m_datacontainer[it->first.substr(2)]=new Blob_Data<double>(it->second->Get<double>());
	m_datacontainer["NULL"+it->first.substr(2)]=new Blob_Data<double>(0.);
      }
    }

    if (SelectBlob(*blit)) {
      for (int i=0;i<(*blit)->NOutP();++i) {
	Particle * p = (*blit)->OutParticle(i);
	if (p->DecayBlob()==NULL || 
	    (m_mode&ANALYSIS::do_hadron)==0 && p->Info()!='G') {
	  if ((p->Info()!='G' &&  p->Info()!='H')
	      || (*blit)->Type()!=btp::IS_Shower) {
	    pl->push_back(new Particle(*p));
	    if (markb && p_bfinder && (*p_bfinder)(p)) {
	      pl->back()->SetFlav(Flavour(kf::bjet));
	    }
	  }
	}
      }
    }
  }

  if (!markb) {
    bool found=false;
    if (m_datacontainer.find("ME_Weight")!=m_datacontainer.end()) found=true;
    if (!found) {
      m_datacontainer["ME_Weight"]=new Blob_Data<double>(1.);
      m_datacontainer["ME_NumberOfTrials"]=new Blob_Data<int>(1);
      m_datacontainer["ME_Weight_One"]=new Blob_Data<double>(0.);
      m_datacontainer["ME_NumberOfTrials_One"]=new Blob_Data<int>(0);
    }
    else if (!m_datacontainer["ME_NumberOfTrials"]) {
      m_datacontainer["ME_NumberOfTrials"]=new Blob_Data<int>(1);
      m_datacontainer["ME_Weight_One"]=new Blob_Data<double>(0.);
      m_datacontainer["ME_NumberOfTrials_One"]=new Blob_Data<int>(0);
    }
  }
  m_pls[key]=pl;
  AddParticleList("NULL",new Particle_List);
}

void Primitive_Analysis::CreatePrimordialHadronsList()
{
  PL_Container::const_iterator cit=m_pls.find("PrimordialHadrons");
  if (cit!=m_pls.end()) return;

  Particle_List * pl = new Particle_List;

  if (m_mode&ANALYSIS::do_hadron) {

    for (Blob_List::const_iterator blit=p_blobs->begin();blit!=p_blobs->end();++blit) {
      if ((*blit)->Type()==btp::Fragmentation) {
	for (int i=0;i<(*blit)->NOutP();++i) {
	  Particle * p = (*blit)->OutParticle(i);
	  if (p->Flav().IsHadron()) pl->push_back(new Particle(*p));
	}
      }
    }
  }

  m_pls["PrimordialHadrons"]=pl;
}

void Primitive_Analysis::CreateIntermediateHadronsList()
{
  PL_Container::const_iterator cit=m_pls.find("IntermediateHadrons");
  if (cit!=m_pls.end()) return;

  Particle_List * pl = new Particle_List;

  if (m_mode&ANALYSIS::do_hadron) {

    for (Blob_List::const_iterator blit=p_blobs->begin();blit!=p_blobs->end();++blit) {
      if ((*blit)->Type()==btp::Hadron_Decay || (*blit)->Type()==btp::Fragmentation) {
	for (int i=0;i<(*blit)->NOutP();++i) {
	  Particle * p = (*blit)->OutParticle(i);
	  if (p->Flav().IsHadron()) pl->push_back(new Particle(*p));
	}
      }
    }
  }

  m_pls["IntermediateHadrons"]=pl;
}

void Primitive_Analysis::CreateChargedParticleList()
{
  PL_Container::const_iterator cit=m_pls.find("ChargedParticle");
  if (cit!=m_pls.end()) return;

  CreateFinalStateParticleList();
  Particle_List * pl_fs=m_pls["FinalState"];
  
  Particle_List * pl = new Particle_List;
  copy_if(pl_fs->begin(),pl_fs->end(),
	  back_inserter(*pl),Is_Charged());

  m_pls["ChargedParticle"]=pl;
}


Particle_List * Primitive_Analysis::GetParticleList(const std::string & key) 
{
  if (!m_active) {
    PL_Container::const_iterator cit=m_pls.find("NULL");
    if (cit!=m_pls.end()) return cit->second;
  }
  PL_Container::const_iterator cit=m_pls.find(key);
  if (cit!=m_pls.end()) return cit->second;

  if (key=="FinalState")               CreateFinalStateParticleList();
  else if (key=="FinalStateB")         CreateFinalStateParticleList(true);
  else if (key=="IntermediateHadrons") CreateIntermediateHadronsList();
  else if (key=="PrimordialHadrons")   CreatePrimordialHadronsList();
  //  else if (key=="ChargedParticle") CreateChargedParticleList();
  if (key=="Analysed") return 0;

  cit=m_pls.find(key);
  if (cit!=m_pls.end()) return cit->second;

  Particle_Selector * ps=0;
  Universal_Selector * us=0;
  std::string testname1=std::string("ParticleSelector_")+key;
  std::string testname2=std::string("UniversalSelector_")+key;
  for(size_t i=0;i<m_observables.size();++i) {
    if (m_observables[i]->Name()==testname1) {
      msg.Error()<<"WARNING in Primitive_Analysis::GetParticleList:"<<std::endl
		 <<"   "<<testname1<<" already present, will continue."<<std::endl;
      ps = static_cast<Particle_Selector*>(m_observables[i]);
      break;
    } 
    if (m_observables[i]->Name()==testname2) {
      msg_Tracking()<<"found matching universal selector"<<std::endl;
      us = static_cast<Universal_Selector*>(m_observables[i]);
      break;
    } 
  }

  if (ps==0 && us==0) {
    ps = new Particle_Selector("FinalState","IntermediateHadrons",key,0);
    AddObservable(ps);
  }
  if (ps) ps->CreateParticleList();
  if (us) us->CreateParticleList();

  cit=m_pls.find(key);
  if (cit!=m_pls.end()) return cit->second;
  msg.Error()<<"WARNING in Primitive_Analysis::GetParticleList:"<<std::endl
	     <<"   "<<key<<" not found, return 0."<<std::endl;

  return 0;
}

void Primitive_Analysis::AddParticleList(const std::string & key,Particle_List * pl) 
{
  PL_Container::const_iterator cit=m_pls.find(key);
  if (cit!=m_pls.end()) {
    for (Particle_List::iterator pit=cit->second->begin(); 
	 pit!=cit->second->end();++pit) delete (*pit);
    delete cit->second;
  }

  m_pls[key]=pl;
}

ATOOLS::Blob_Data_Base * Primitive_Analysis::operator[](const std::string name) 
{
  ATOOLS::String_BlobDataBase_Map::const_iterator cit;
  if (!m_active) {
    cit=m_datacontainer.find("NULL"+name);
    if (cit!=m_datacontainer.end()) return cit->second;
  }
  cit=m_datacontainer.find(name);
  if (cit==m_datacontainer.end()) return 0;
  return cit->second;
} 

void Primitive_Analysis::AddData(const std::string name, Blob_Data_Base * data) 
{
  String_BlobDataBase_Map::iterator it=m_datacontainer.find(name);
  if (it==m_datacontainer.end()) {
    m_datacontainer[name]=data;
  }
  else {
    delete it->second;
    it->second=data;
  }
}

void Primitive_Analysis::ClearAllData() 
{
  for (PL_Container::iterator it=m_pls.begin();
       it!=m_pls.end(); ++it) {
    if (it->second->size()>0) {
      for (Particle_List::iterator pit=it->second->begin(); 
	   pit!=it->second->end();++pit) delete (*pit);
    }
    delete it->second;
  }
  m_pls.clear();

  for (String_BlobDataBase_Map::iterator it=m_datacontainer.begin();
       it!=m_datacontainer.end(); ++it) delete it->second;
  m_datacontainer.clear();
}

void Primitive_Analysis::PrintStatus() 
{

  msg.Out()<<"Particle_Lists:"<<std::endl;
  for (PL_Container::iterator it=m_pls.begin();
       it!=m_pls.end(); ++it) {
    msg.Out()<<"   * "<<it->first<<" ("<<it->second->size()<<")"<<std::endl;
  }
  for (PL_Container::iterator it=m_pls.begin();
       it!=m_pls.end(); ++it) {
    msg.Out()<<"   * "<<it->first<<std::endl<<*it->second<<std::endl;
  }

  msg.Out()<<"Data_Container:"<<std::endl;
  for (String_BlobDataBase_Map::iterator it=m_datacontainer.begin();
       it!=m_datacontainer.end(); ++it) {
    msg.Out()<<"   * "<<it->first<<" ("<<*(it->second)<<")"<<std::endl;
  }
}

void Primitive_Analysis::SetPartner(Primitive_Analysis * const ana)
{
  p_partner=ana;
}

Primitive_Observable_Base * Primitive_Analysis::GetObservable(const std::string & key)
{
  for (size_t i=0;i<m_observables.size();i++) {
    if (m_observables[i]->Name()==key) return m_observables[i];
  }
  return 0;
}



template <>
std::ostream & Blob_Data<std::vector<double> *>::operator>>(std::ostream & s) const 
{
  if (m_data->size()>0) 
    s<<(*m_data)[0];
  for (size_t i=1;i<m_data->size();++i) 
    s<<","<<(*m_data)[i];
  return s;
}

template <> Blob_Data<std::vector<double> *>::~Blob_Data() 
{
  delete m_data;
}

template class Blob_Data<std::vector<double> *>;
