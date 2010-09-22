#include "ATOOLS/Org/CXXFLAGS.H"
#include "SHERPA/Tools/RootNtuple_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include <iostream>

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;



RootNtuple_Reader::RootNtuple_Reader(const std::string & path,const std::string & file,int mode) :
  Event_Reader_Base(path,file), 
  m_eventmode(mode), m_evtid(0), m_subevtid(0), m_evtcnt(0), m_entries(0), m_evtpos(0)
{
  std::string filename=m_path+m_file;
  msg_Out()<<" Reading from "<<filename<<"\n";

#ifdef USING__ROOT
  p_f=new TChain("t3");
  p_f->Add(filename.c_str());
  m_entries=p_f->GetEntries();
  if(m_entries==0) {
    msg_Error()<<"ERROR: Event file "<<filename<<" does not contain any event."<<std::endl
	       <<"   Will abort the run."<<std::endl;
    abort();
  }
  p_f->SetBranchAddress("id",&m_id);
  p_f->SetBranchAddress("nparticle",&m_nparticle);
  p_f->SetBranchAddress("px",p_px);
  p_f->SetBranchAddress("py",p_py);
  p_f->SetBranchAddress("pz",p_pz);
  p_f->SetBranchAddress("E",p_E);

  p_f->SetBranchAddress("kf",p_kf);
  p_f->SetBranchAddress("weight",&m_wgt);
  p_f->SetBranchAddress("weight2",&m_wgt2);

  msg_Out()<<"Event mode = "<<m_eventmode<<std::endl;
#else
  msg_Error()<<"Sherpa must be linked with root to read in root files!"<<endl;
#endif
}



void RootNtuple_Reader::CloseFile() {
#ifdef USING__ROOT
  delete p_f;
#endif
}





bool RootNtuple_Reader::FillBlobs(Blob_List * blobs) 
{
  bool result;
  switch (m_eventmode) {
  case 1: result=ReadInFullEvent(blobs);
    break;
  default: result=ReadInSubEvent(blobs);
  }
  if (result==0) rpa.gen.SetNumberOfEvents(rpa.gen.NumberOfDicedEvents());
  
  long nev=rpa.gen.NumberOfEvents();
  if(nev==rpa.gen.NumberOfDicedEvents()) CloseFile();
  return result;
}

bool RootNtuple_Reader::ReadInEntry() 
{
  if (m_evtpos>=m_entries) return 0;
#ifdef USING__ROOT
  p_f->GetEntry(m_evtpos);
  m_evtpos++;
  m_evtid=m_id;
#endif
  return 1;
}

bool RootNtuple_Reader::ReadInSubEvent(Blob_List * blobs) 
{
  if (!ReadInEntry()) return 0;
  Blob         *signalblob=new Blob();
  signalblob->SetType(btp::Signal_Process);
  signalblob->SetTypeSpec("NLOsubevt");
  signalblob->SetId();
  signalblob->SetPosition(Vec4D(0.,0.,0.,0.));
  signalblob->SetStatus(blob_status::code(30));
  signalblob->SetBeam(-1);
#ifdef USING__ROOT
  m_weight=m_wgt;
  signalblob->SetWeight(m_weight);
  for (int i=0;i<m_nparticle;i++) {
    Particle *part=new Particle(i,Flavour(abs(p_kf[i]),p_kf[i]<0),Vec4D(p_E[i],p_px[i],p_py[i],p_pz[i]));
    signalblob->AddToOutParticles(part);
  }
  signalblob->AddData("Weight",new Blob_Data<double>(m_weight));
  signalblob->AddData("Trials",new Blob_Data<double>(1.));
#endif


  blobs->push_back(signalblob);
  m_evtcnt++;
  return 1;
}


bool RootNtuple_Reader::ReadInFullEvent(Blob_List * blobs) 
{
  if (!m_nlos.empty()) {
    for (size_t i=0;i<m_nlos.size();i++) {
      delete[] m_nlos[i]->p_fl;
      delete[] m_nlos[i]->p_mom;
    }
    m_nlos.clear();
  }
  if (m_evtid==0) if (!ReadInEntry()) return 0;
  size_t currentid=m_evtid;
  Blob         *signalblob=new Blob();
  m_weight = 0.;
  signalblob->SetType(btp::Signal_Process);
  signalblob->SetTypeSpec("NLO");
  signalblob->SetId();
  signalblob->SetPosition(Vec4D(0.,0.,0.,0.));
  signalblob->SetStatus(blob_status::code(30));
  signalblob->SetBeam(-1);
  
#ifdef USING__ROOT
  while (currentid==m_evtid) {
    Vec4D *moms = new Vec4D[2+m_nparticle];
    Flavour *flav = new Flavour[2+m_nparticle];
    for (int i=0;i<m_nparticle;i++) {
      moms[i+2]=Vec4D(p_E[i],p_px[i],p_py[i],p_pz[i]);
      flav[i+2]=Flavour(abs(p_kf[i]),p_kf[i]<0);
    }
    m_nlos.push_back(new NLO_subevt(m_nparticle+2,NULL,flav,moms));
    m_nlos.back()->m_result=m_wgt2;
    m_nlos.back()->m_flip=0;
    m_weight+=m_wgt2;
    if (!ReadInEntry()) m_evtid=0;
  }  
#endif
  signalblob->SetWeight(m_weight);
  signalblob->AddData("Weight",new Blob_Data<double>(m_weight));
  signalblob->AddData("Trials",new Blob_Data<int>(1.));
  signalblob->AddData("NLO_subeventlist",new Blob_Data<NLO_subevtlist*>(&m_nlos));
  blobs->push_back(signalblob);
  m_evtcnt++;  
  return 1;
}
