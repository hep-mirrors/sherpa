#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "AddOns/Analysis/Main/Analysis_Object.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "PHASIC++/Process/NLO_Helpers.H"

#ifdef PROFILE__Analysis_Phase
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif

using namespace ANALYSIS;
using namespace ATOOLS;

Primitive_Analysis::Primitive_Analysis(const std::string _name, const int mode) :
  m_active(true), m_splitjetconts(true)
{
  m_nevt = 0;
  p_partner = this;
  m_mode = mode;

  m_name = std::string("Analysis : ") + _name;
  msg_Tracking()<<" Initializing Primitive_Analysis : "<<m_name<<std::endl;
}

Primitive_Analysis::Primitive_Analysis(const int mode) :
  m_nevt(0), p_partner(this), m_active(true), m_splitjetconts(true)
{
  m_mode = mode;

  m_name = std::string("Analysis : noname");
  msg_Tracking()<<" Initializing Primitive_Analysis : "<<m_name<<std::endl;
}

Primitive_Analysis::~Primitive_Analysis()
{
  for(int i=m_objects.size();i>0;i--) {
    if (m_objects[i-1]) delete m_objects[i-1];
  }
  m_objects.clear();

  for (Analysis_List::iterator it=m_subanalyses.begin();it!=m_subanalyses.end();++it) 
    delete it->second;
  m_subanalyses.clear();
}

void Primitive_Analysis::AddObject(Analysis_Object * obs) 
{
  obs->SetAnalysis(p_partner);
  std::string oname=obs->Name();
  std::string id="_A";
  size_t pos=oname.find(".dat");
  for(size_t i=0;i<m_objects.size();++i) {
    if (m_objects[i]->Name()==obs->Name()) {
      std::string pname;
      if (pos!=std::string::npos)
	pname=oname.substr(0,pos)+id+oname.substr(pos);
      else 
	pname=oname+id;
      obs->SetName(pname);
      ++id[1];
    }
  }
  m_objects.push_back(obs);
}

void Primitive_Analysis::AddSubAnalysis(const std::string & key,Primitive_Analysis * ana)
{
  Analysis_List::const_iterator cit=m_subanalyses.find(key);
  if (cit!=m_subanalyses.end()) {
    msg_Out()<<"WARNING in Primitive_Analysis::AddSubAnalysis :"<<std::endl
	     <<"   Analysis "<<key<<" already existent;"<<std::endl
	     <<" sub analysis not added, will be deleted."<<std::endl;
    if (ana) delete ana;
    return;
  }

  m_subanalyses[key]=ana;  
}

Primitive_Analysis * Primitive_Analysis::GetSubAnalysis
(const Blob_List *const bl,const std::string & key, int mode) 
{
  Analysis_List::const_iterator cit=m_subanalyses.find(key);
  if (cit!=m_subanalyses.end()) return cit->second;

  bool master=true;
  if (key=="ME" || key=="MENLO" || key=="MI" || key=="Shower" || key=="Hadron") {
    master=false;
    if (key!="ME"     && mode&ANALYSIS::do_me) mode=mode^ANALYSIS::do_me;
    if (key!="MENLO"  && mode&ANALYSIS::do_menlo) mode=mode^ANALYSIS::do_menlo;
    if (key!="MI"     && mode&ANALYSIS::do_mi) mode=mode^ANALYSIS::do_mi;
    if (key!="Shower" && mode&ANALYSIS::do_shower) mode=mode^ANALYSIS::do_shower;
    if (key!="Hadron" && mode&ANALYSIS::do_hadron) mode=mode^ANALYSIS::do_hadron;
  }

  Primitive_Analysis * ana = new Primitive_Analysis(m_name.substr(11)+key,mode);
  if (master) ana->SetPartner(p_partner);
  ana->SetMaxJetTag(m_maxjettag);

  for (size_t i=0;i<m_objects.size();i++) {
    if (m_objects[i]->IsObservable() || !master) 
      ana->AddObject(m_objects[i]->GetCopy());
  }
  for (long int i(0);i<m_nevt-1;++i) ana->DoAnalysis(bl,0.0);
  m_subanalyses[key]=ana;
  return ana;
}

std::string Primitive_Analysis::JetID
(std::string name,std::string max) const
{
  size_t jets(1), maxjets(100);
  std::string subprocs;
  if (max.length()>0) {
    size_t pos(max.find('['));
    if (pos!=std::string::npos) {
      maxjets=ToType<int>(max.substr(0,pos));
      max=max.substr(pos);
    }
    else {
      maxjets=ToType<int>(max);
      max="";
    }
  }
  for (size_t i(1);i<name.length();++i) {
    if (name[i]=='_' && name[i-1]=='_') ++jets;
    else if (name[i]=='[') {
      std::string cmax;
      for (size_t j(0);j<max.length();++j) {
	if (max[j]=='[') {
	  int open(1);
	  for (size_t k(j+1);k<max.length();++k) {
	    if (max[k]=='[') ++open;
	    if (max[k]==']') --open;
	    if (open==0) {
	      cmax=max.substr(j+1,k-j-1);
	      max=max.substr(k+1);
	      break;
	    }
	  }
	}
      }
      int open(1);
      for (size_t j(i+1);j<name.length();++j) {
	if (name[j]=='[') ++open;
	if (name[j]==']') --open;
	if (open==0) {
	  if (jets>1) subprocs+=ToString(jets-1);
	  subprocs+="["+JetID(name.substr(i+1,j-i-1),cmax)+"]";
	  jets=0;
	  i=j;
	  break;
	}
      }
    }
  }
  if (maxjets<jets) return "X";
  return subprocs+ToString(jets);
}

void Primitive_Analysis::CallSubAnalysis(const Blob_List * const bl, double value) 
{
  Blob *sp(bl->FindFirst(btp::Signal_Process));
  std::string name("no_signal_process");
  if (sp==NULL) {
    msg_Debugging()<<"WARNING in Primitive_Analysis::CallSubAnalysis: no Signal process found "<<std::endl;
    // if no signal process (i.e. hadrons execs etc.), proceed anyways
    // do not split jetseeds
    if (m_mode&ANALYSIS::splitt_jetseeds)
      m_mode=m_mode^ANALYSIS::splitt_jetseeds;
  }
  else name = sp->TypeSpec();

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

    GetSubAnalysis(bl,key,mode);

    for (Analysis_List::iterator it=m_subanalyses.begin();it!=m_subanalyses.end();++it) {
      if (it->first.find("extra_")!=std::string::npos) {
	if (it->first==key) {
	  it->second->DoAnalysis(bl,value);
	  m_called.insert(it->second);
	}
	else {
	  m_active=false;
	  it->second->DoAnalysis(bl,value);
	  m_called.insert(it->second);
	  m_active=true;
	}
      }
    }
  }
  if (m_mode&ANALYSIS::splitt_jetseeds) {
    mode=m_mode^ANALYSIS::splitt_jetseeds;
    if (mode&ANALYSIS::splitt_extra)
      mode=m_mode^ANALYSIS::splitt_extra;
    if (!m_splitjetconts)
      mode=mode-(mode&ANALYSIS::output_this);
    std::string fsname(name.substr(name.find("__")+3));
    fsname=fsname.substr(fsname.find("__")+3);
    fsname=fsname.substr(fsname.find("__")+2);
    key=JetID(fsname,m_maxjettag);
    key="j"+key;
  }
  else {
    mode=m_mode^ANALYSIS::splitt_process;
//     if (m_mode&ANALYSIS::output_process) mode=mode|ANALYSIS::output_this;
//     else 
    if (m_mode&ANALYSIS::output_this) mode=mode^ANALYSIS::output_this;
      key=name;
  }
  if (key.find('X')!=std::string::npos) {
    msg_Debugging()<<METHOD<<"(): Max jet number reached in '"<<key<<"'\n";
  }

  Primitive_Analysis * ana=GetSubAnalysis(bl,key,mode);
  ana->DoAnalysis(bl,value);
  m_called.insert(ana);
}


void Primitive_Analysis::DoAnalysis(const Blob_List * const bl, const double value)
{
  ++m_nevt;
  m_called.clear();
  if (IsNan(value)) {
    msg_Error()<<METHOD<<"(): Event weight is nan. Skip."<<std::endl;
    return;
  }

  if (m_mode&ANALYSIS::splitt_phase) {
    m_mode=m_mode|ANALYSIS::output_this;
    int mode=m_mode^ANALYSIS::splitt_phase;
    if (m_mode&ANALYSIS::do_me)     {
      Primitive_Analysis *ana(GetSubAnalysis(bl,"ME",mode));
      ana->DoAnalysis(bl,value);
      m_called.insert(ana);
    }
    if (m_mode&ANALYSIS::do_menlo)     {
      Primitive_Analysis *ana(GetSubAnalysis(bl,"MENLO",mode));
      ana->DoAnalysis(bl,value);
      m_called.insert(ana);
    }
    if (m_mode&ANALYSIS::do_mi)     {
      Primitive_Analysis *ana(GetSubAnalysis(bl,"MI",mode));
      ana->DoAnalysis(bl,value);
      m_called.insert(ana);
    }
    if (m_mode&ANALYSIS::do_shower)     {
      Primitive_Analysis *ana(GetSubAnalysis(bl,"Shower",mode));
      ana->DoAnalysis(bl,value);
      m_called.insert(ana);
    }
    if (m_mode&ANALYSIS::do_hadron)     {
      Primitive_Analysis *ana(GetSubAnalysis(bl,"Hadron",mode));
      ana->DoAnalysis(bl,value);
      m_called.insert(ana);
    }
    return;
  }
  if (m_mode&ANALYSIS::do_menlo) {
    if (DoAnalysisNLO(bl,value)) return;
  }

  ClearAllData();
  p_blobs = bl;

  if (p_partner==this) {
    m_mode=m_mode|ANALYSIS::fill_helper;
    m_mode=m_mode|ANALYSIS::output_this;
  }
  else if (m_mode&ANALYSIS::fill_helper) m_mode=m_mode^ANALYSIS::fill_helper;
  if (m_mode&ANALYSIS::splitt_process) {
    m_mode=m_mode|ANALYSIS::output_process;
  }
  if ((m_mode&ANALYSIS::splitt_all)==0) m_mode=m_mode|ANALYSIS::fill_histos;
  Init();
  Blob *sp(bl->FindFirst(btp::Signal_Process));
  // if no signal process present (i.e. hadrons execs etc.),
  // assume weight=1, ncount=1
  double weight(1.), ncount(1.);
  if (sp) {
    weight=(*sp)["Weight"]->Get<double>();
    ncount=(*sp)["Trials"]->Get<double>();
  }
  if (!IsEqual(value,weight)) {
    if (p_partner==this) {
      msg_Out()<<"WARNING in Primitive_Analysis::DoAnalysis :"<<std::endl
	       <<"   Weight in Primitive_Analysis ambiguous! ("<<value<<","<<weight<<")"<<std::endl;
    }
    else if (value==0.) {
      weight=0.;
    }
    else {
      msg_Out()<<"WARNING something is wrong in Primitive_Analysis::DoAnalysis :"<<std::endl
	       <<"   Weight in Primitive_Analysis ambiguous! ("<<value<<","<<weight<<")"<<std::endl;
    }
  }
  // do nonsplittable (helper and legacy objects) first
  if (m_mode&ANALYSIS::fill_helper) {
    for (size_t i=0;i<m_objects.size();i++) {
      if (!m_objects[i]->IsObservable()) {
	m_objects[i]->Evaluate(*bl,value,ncount);
      }
    }
  }

  if (m_mode&ANALYSIS::fill_histos) {
    for (size_t i=0;i<m_objects.size();i++) {
      if (m_objects[i]->IsObservable()) {
	m_objects[i]->Evaluate(*bl,value,ncount);
      }
    }
  }


  if (m_mode&ANALYSIS::splitt_all) CallSubAnalysis(bl,value);
  if (p_partner==this && msg_LevelIsTracking()) PrintStatus();

  for (Analysis_List::const_iterator ait(m_subanalyses.begin());
       ait!=m_subanalyses.end();++ait)
    if (m_called.find(ait->second)==m_called.end()) 
      ait->second->DoAnalysis(bl,0.0);

  ClearAllData();
}


bool Primitive_Analysis::DoAnalysisNLO(const Blob_List * const bl, const double value) {
  if (IsNan(value)) {
    msg_Error()<<METHOD<<"(): Event weight is nan. Skip."<<std::endl;
    return 0;
  }

  ClearAllData();
  p_blobs = bl;
  Blob_Data_Base* info=0;
  Blob* signal=0;
  for (Blob_List::const_iterator blit=p_blobs->begin();
       blit!=p_blobs->end();++blit) 
    if ((*blit)->Type()==btp::Signal_Process) {
      signal=(*blit);
      info=(*signal)["NLO_subeventlist"];
      break;
    }

  if (!info) return 0;

  PHASIC::NLO_subevtlist* nlos = info->Get<PHASIC::NLO_subevtlist*>();

  // (*nlos)*=(*signal)["XS_Weight"]->Get<double>();
  double ncount=(*signal)["Trials"]->Get<double>();

  for (size_t j=0;j<nlos->size();j++) {
    m_pls[finalstate_list]=(*nlos)[j]->CreateParticleList();
  
    // do nonsplittable (helper and legacy objects) first
    if (m_mode&ANALYSIS::fill_helper) {
      for (size_t i=0;i<m_objects.size();i++) {
	if (!m_objects[i]->IsObservable()) {
 	  m_objects[i]->Evaluate(*bl,(*nlos)[j]->m_result,ncount);
	}
      }
    }
    
    if (m_mode&ANALYSIS::fill_histos) {
      for (size_t i=0;i<m_objects.size();i++) {
	if (m_objects[i]->IsObservable()) {
	  m_objects[i]->EvaluateNLOcontrib((*nlos)[j]->m_result,ncount);
	}
      }
    }
    ClearAllData();

  }
  ++m_nevt;
  for (size_t i=0;i<m_objects.size();i++) {
    if (m_objects[i]->IsObservable())
      m_objects[i]->EvaluateNLOevt();
  }

  return 1;
}


void Primitive_Analysis::FinishAnalysis(const std::string & resdir) 
{
  if (m_mode&ANALYSIS::output_this) 
    ATOOLS::MakeDir(resdir+OutputPath()); 

  if (m_mode&ANALYSIS::do_menlo) {
    for (size_t i=0;i<m_objects.size();i++) {
      m_objects[i]->EndEvaluation(1.);
      m_objects[i]->Output(resdir+OutputPath());
    }
    return;
  }
  for (Analysis_List::iterator it=m_subanalyses.begin();
       it!=m_subanalyses.end();++it) {
    std::string dir=resdir+OutputPath()+std::string("/")+it->first;
    it->second->FinishAnalysis(dir);
  }

  if (!(m_mode&ANALYSIS::splitt_phase)) {
    for (size_t i=0;i<m_objects.size();i++) {
      m_objects[i]->EndEvaluation();
      if (m_mode&ANALYSIS::output_this) 
	m_objects[i]->Output(resdir+OutputPath());
    }
  }
}

void Primitive_Analysis::RestoreAnalysis() 
{
  for (Analysis_List::iterator it=m_subanalyses.begin();
       it!=m_subanalyses.end();++it)
    it->second->RestoreAnalysis();
  for (size_t i(0);i<m_objects.size();++i)
    m_objects[i]->Restore();
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
      blob->Type()==btp::Shower) return true;
  if (m_mode&ANALYSIS::do_mi && 
      (blob->Type()==btp::Hard_Collision ||
       blob->Type()==btp::Signal_Process)) return true;
  if ((m_mode&ANALYSIS::do_me||m_mode&ANALYSIS::do_menlo) && 
      blob->Type()==btp::Signal_Process) return true;
  return false;
}

void Primitive_Analysis::CreateFinalStateParticleList()
{
  PL_Container::const_iterator cit=m_pls.find(finalstate_list);
  if (cit!=m_pls.end()) return;

  Particle_List * pl = new Particle_List;

  for (Blob_List::const_iterator blit=p_blobs->begin();
       blit!=p_blobs->end();++blit) {
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
	if ((p->DecayBlob()==NULL || (m_mode&ANALYSIS::do_hadron)==0) &&
	    p->Info()!='G')
	  pl->push_back(p);
      }
    }
  }

  m_pls[finalstate_list]=pl;
  AddParticleList("NULL",new Particle_List);
}

Particle_List * Primitive_Analysis::GetParticleList(const std::string & key,
						    const bool nocreate) 
{
  if (!m_active) {
    PL_Container::const_iterator cit=m_pls.find("NULL");
    if (cit!=m_pls.end()) return cit->second;
  }
  PL_Container::const_iterator cit=m_pls.find(key);
  if (cit!=m_pls.end()) return cit->second;
  if (nocreate) return NULL;
  if (key==finalstate_list) CreateFinalStateParticleList();
  cit=m_pls.find(key);
  if (cit!=m_pls.end()) return cit->second;
  msg_Error()<<METHOD<<"(): List '"<<key<<"' not found."<<std::endl;
  return NULL;
}

void Primitive_Analysis::AddParticleList(const std::string & key,Particle_List * pl) 
{
  PL_Container::const_iterator cit=m_pls.find(key);
  if (cit!=m_pls.end()) {
    for (Particle_List::iterator pit=cit->second->begin(); 
	 pit!=cit->second->end();++pit) 
      if ((*pit)->ProductionBlob()==NULL && (*pit)->DecayBlob()==NULL) delete *pit;
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
    if (!it->second->empty()) {
      for (Particle_List::iterator pit=it->second->begin(); 
      	   pit!=it->second->end();++pit) 
      if ((*pit)->ProductionBlob()==NULL && (*pit)->DecayBlob()==NULL) 
	delete *pit;
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

  msg_Out()<<"Particle_Lists:"<<std::endl;
  for (PL_Container::iterator it=m_pls.begin();
       it!=m_pls.end(); ++it) {
    msg_Out()<<"   * "<<it->first<<" ("<<it->second->size()<<")"<<std::endl;
  }
  for (PL_Container::iterator it=m_pls.begin();
       it!=m_pls.end(); ++it) {
    msg_Out()<<"   * "<<it->first<<std::endl<<*it->second<<std::endl;
  }

  msg_Out()<<"Data_Container:"<<std::endl;
  for (String_BlobDataBase_Map::iterator it=m_datacontainer.begin();
       it!=m_datacontainer.end(); ++it) {
    msg_Out()<<"   * "<<it->first<<" ("<<*(it->second)<<")"<<std::endl;
  }
}

void Primitive_Analysis::SetPartner(Primitive_Analysis * const ana)
{
  p_partner=ana;
}

Analysis_Object * Primitive_Analysis::GetObject(const std::string & key)
{
  for (size_t i=0;i<m_objects.size();i++) {
    if (m_objects[i]->Name()==key) return m_objects[i];
  }
  return 0;
}

void Primitive_Analysis::Test(const int mode) 
{
  std::cout<<"Number of objects: "<<m_objects.size()<<std::endl;
  for (size_t i=0;i<m_objects.size();i++) {
    if (!m_objects[i]->IsObservable()) m_objects[i]->Test(mode);
  }
  for (size_t i=0;i<m_objects.size();i++) {
    if (m_objects[i]->IsObservable()) m_objects[i]->Test(mode);
  }
}

namespace ATOOLS {

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

}
