#include "Primitive_Analysis.H"
#include "Primitive_Observable_Base.H"
#include "Message.H"
#include "List_Algorithms.H"
#include "MyStrStream.H"

#ifdef PROFILE__Analysis_Phase
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif

using namespace ANALYSIS;
using namespace ATOOLS;

Primitive_Analysis::Primitive_Analysis(const std::string _name, const int mode) :
  m_nevt(0), p_partner(this)
{
  m_mode = mode;

  m_name = std::string("Analysis : ") + _name;
  msg.Tracking()<<" Initializing Primitive_Analysis : "<<m_name<<std::endl;
}

Primitive_Analysis::Primitive_Analysis(const int mode) :
  m_nevt(0), p_partner(this)
{
  m_mode = mode;

  m_name = std::string("Analysis : noname");
  msg.Tracking()<<" Initializing Primitive_Analysis : "<<m_name<<std::endl;
}

Primitive_Analysis::~Primitive_Analysis()
{
  for(int i=m_observables.size();i>0;i--) {
    if (m_observables[i-1]) delete m_observables[i-1];
  }
  m_observables.clear();

  //
  for (Analysis_List::iterator it=m_subanalyses.begin();
       it!=m_subanalyses.end();++it) {
    delete it->second;
  }
  m_subanalyses.clear();
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
    msg.Out()<<" WARNING: Analysis "<<key<<" already existent;"<<std::endl
	     <<"          sub analysis not added!"<<std::endl;
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
  if (key=="ME" || key=="Shower" || key=="Hadron") {
    master=false;
    if (key!="ME"     && mode&ANALYSIS::do_me) mode=mode^ANALYSIS::do_me;
    if (key!="Shower" && mode&ANALYSIS::do_shower) mode=mode^ANALYSIS::do_shower;
    if (key!="Hadron" && mode&ANALYSIS::do_hadron) mode=mode^ANALYSIS::do_hadron;
  }

  Primitive_Analysis * ana=new Primitive_Analysis(m_name.substr(11)+key,mode);
  if (master) ana->SetPartner(p_partner);

  for (size_t i=0;i<m_observables.size();i++) {
    if (m_observables[i]->Splittable() || !master) 
      ana->AddObservable(m_observables[i]->Copy());
  }
  m_subanalyses[key]=ana;
  return ana;
}

void Primitive_Analysis::CallSubAnalysis(Blob_List * const bl, double value) 
{
  int nout=-1;
  std::string name;
  for (Blob_Const_Iterator bit=bl->begin();bit!=bl->end();++bit) {
    if ((*bit)->Type()==btp::Signal_Process) {
      nout  = (*bit)->NOutP();
      // orig: (*bit)->Type();
      name  = (*bit)->TypeSpec();
      break;
    }
  }
  if (nout==-1) {
    msg.Out()<<" no Signal process found "<<std::endl;
    return;
  }

  std::string key;
  int mode;
  if (m_mode&ANALYSIS::splitt_jetseeds) {
    mode=m_mode^ANALYSIS::splitt_jetseeds;
    mode=mode|ANALYSIS::output_this;
    MyStrStream str;
    str<<'j';
    str<<nout;
    str>>key;
  }
  else {
    mode=m_mode^ANALYSIS::splitt_process;
//     if (m_mode&ANALYSIS::output_process) mode=mode|ANALYSIS::output_this;
//     else 
      if (m_mode&ANALYSIS::output_this) mode=mode^ANALYSIS::output_this;

      // orig:    key=name.substr(17);
      key=name;
  }
  
  Primitive_Analysis * ana=GetSubAnalysis(key,mode);
  ana->DoAnalysis(bl,value);
}


void Primitive_Analysis::DoAnalysis(Blob_List * const bl, double value) {
  ++m_nevt;

  if (m_mode&ANALYSIS::splitt_phase) {
    m_mode=m_mode|ANALYSIS::output_this;
    int mode=m_mode^ANALYSIS::splitt_phase;
    if (m_mode&ANALYSIS::do_me)     GetSubAnalysis("ME",mode)->DoAnalysis(bl,value);
    if (m_mode&ANALYSIS::do_shower) GetSubAnalysis("Shower",mode)->DoAnalysis(bl,value);
    if (m_mode&ANALYSIS::do_hadron) GetSubAnalysis("Hadron",mode)->DoAnalysis(bl,value);
    return;
  }

  ClearAllData();
  p_blobs = bl;

  if (value!=1.) m_mode=m_mode|ANALYSIS::weighted;
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
  int    ncount=(*p_partner)["ME_NumberOfTrials"]->Get<int>();
  if (!IsEqual(value,weight)) std::cout<<" WARNING: weight in Primitive_Analysis ambiguous! "<<std::endl;

  // do nonsplittable (helper and legacy observables) first
  if (m_mode&ANALYSIS::fill_helper) {
    for (size_t i=0;i<m_observables.size();i++) {
      if (!m_observables[i]->Splittable()) {
	m_observables[i]->Evaluate(*bl,value,ncount);
      }
    }
  }

  if (m_mode&ANALYSIS::fill_histos) {
    for (size_t i=0;i<m_observables.size();i++) {
      if (m_observables[i]->Splittable()) {
	m_observables[i]->Evaluate(*bl,value,ncount);
      }
    }
  }


  if (m_mode&ANALYSIS::splitt_all)   CallSubAnalysis(bl,value);
  PrintStatus();

  ClearAllData();
}

void Primitive_Analysis::FinishAnalysis(const std::string & resdir,long ntotal, double xs) 
{
  if (ntotal==0) ntotal=m_nevt;
  int  mode_dir = 448;
  if (m_mode&ANALYSIS::output_this) mkdir(resdir.c_str(),mode_dir); 

  for (Analysis_List::iterator it=m_subanalyses.begin();
       it!=m_subanalyses.end();++it) {
    std::string dir=resdir+std::string("/")+it->first;
    it->second->FinishAnalysis(dir,ntotal,xs);
  }

  if (!(m_mode&ANALYSIS::splitt_phase)) {
    for (size_t i=0;i<m_observables.size();i++) {
      if (m_mode&ANALYSIS::weighted  && m_mode&ANALYSIS::splitt_all && m_observables[i]->Splittable()) {
	std::string key =m_observables[i]->Name();
	m_observables[i]->Reset();  // just in case someone has filled something in already
      
	for (Analysis_List::iterator it=m_subanalyses.begin();
	     it!=m_subanalyses.end();++it) {
	  Primitive_Observable_Base * ob = it->second->GetObservable(key);
	  if (ob) {
	    (*m_observables[i])+=(*ob);
	  }
	}
      }
      else {
	if ((m_mode&ANALYSIS::weighted)==0 ) {
	  m_observables[i]->EndEvaluation(double(m_nevt)/double(ntotal)*xs);
	}
	else {
	  m_observables[i]->EndEvaluation();
	}
      }
      if (m_mode&ANALYSIS::output_this) m_observables[i]->Output(resdir);
    }
  }
}

void Primitive_Analysis::Init()
{
  if (m_mode&ANALYSIS::fill_helper)
    CreateFinalStateParticleList();
}

void Primitive_Analysis::CreateFinalStateParticleList()
{
  PL_Container::const_iterator cit=m_pls.find("FinalState");
  if (cit!=m_pls.end()) return;

  //  msg.Out()<<" creating FinalState particle list "<<std::endl;
  Particle_List * pl = new Particle_List;

  for (Blob_Const_Iterator blit=p_blobs->begin();blit!=p_blobs->end();++blit) {
    if ((*blit)->Type()==btp::Signal_Process) {
      Blob_Data_Base * info=(*(*blit))["ME_Weight"];
      if (info) {
	m_datacontainer["ME_Weight"]=new Blob_Data<double>(info->Get<double>());
	info=(*(*blit))["ME_NumberOfTrials"];
	if (info) {
	  m_datacontainer["ME_NumberOfTrials"]=new Blob_Data<int>(info->Get<int>());
	}
      }
    }
    if (m_mode&ANALYSIS::do_hadron ||
        (m_mode&ANALYSIS::do_shower && ((*blit)->Type()==btp::IS_Shower || (*blit)->Type()==btp::FS_Shower)) ||
        (m_mode&ANALYSIS::do_me && (*blit)->Type()==btp::Signal_Process)) {
      for (int i=0;i<(*blit)->NOutP();++i) {
	Particle * p = (*blit)->OutParticle(i);
	if (p->DecayBlob()==NULL || 
	    (m_mode&ANALYSIS::do_hadron)==0 && p->Info()!='G') {
	  if ((p->Info()!='G' &&  p->Info()!='H')
	      || (*blit)->Type()!=btp::IS_Shower)
	    pl->push_back(new Particle(p));
	}
      }
    }
  }

  if (!m_datacontainer["ME_Weight"]) {
    std::cout<<" ME_Weight not in data container! "<<std::endl;
    m_datacontainer["ME_Weight"]=new Blob_Data<double>(1.);
    m_datacontainer["ME_NumberOfTrials"]=new Blob_Data<int>(1);
  }
  else if (!m_datacontainer["ME_NumberOfTrials"]) {
    m_datacontainer["ME_NumberOfTrials"]=new Blob_Data<int>(1);
  }

  m_pls["FinalState"]=pl;
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

  for (Particle_List::iterator it=pl->begin(); it!=pl->end();++it) {
    (*it)= new Particle(*it);
  }
  
  m_pls["ChargedParticle"]=pl;
}


Particle_List * Primitive_Analysis::GetParticleList(const std::string & key) 
{
  PL_Container::const_iterator cit=m_pls.find(key);
  if (cit!=m_pls.end()) return cit->second;

  if (key=="FinalState") CreateFinalStateParticleList();
  else if (key=="ChargedParticle") CreateChargedParticleList();

  cit=m_pls.find(key);
  if (cit!=m_pls.end()) return cit->second;

  //  msg.Out()<<key<<" not found in particle list container, Primitive_Analysis::GetParticleList failed "<<std::endl;
  return 0;
}

void Primitive_Analysis::AddParticleList(const std::string & key,Particle_List * pl) 
{
  PL_Container::const_iterator cit=m_pls.find(key);
  if (cit!=m_pls.end()) delete cit->second;

  m_pls[key]=pl;
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
    for (Particle_List::iterator pit=it->second->begin(); 
	 pit!=it->second->end();++pit) delete (*pit);
    delete it->second;
  }
  m_pls.clear();

  for (String_BlobDataBase_Map::iterator it=m_datacontainer.begin();
       it!=m_datacontainer.end(); ++it) delete it->second;
  m_datacontainer.clear();
}

void Primitive_Analysis::PrintStatus() 
{
  if (!msg.LevelIsTracking()) return; 
  //  std::cout<<" ana name: "<<m_name<<" "<<m_mode<<std::endl;

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

// ----------------------------------------------------------------------
// probably obsolete:


void Primitive_Analysis::DoAnalysis(const Particle_List & pl, double value) 
{
  ++m_nevt;
  for (size_t i=0;i<m_observables.size();i++) m_observables[i]->Evaluate(pl,value);
}
