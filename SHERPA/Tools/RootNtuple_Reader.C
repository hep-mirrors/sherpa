#include "ATOOLS/Org/CXXFLAGS.H"
#include "SHERPA/Tools/RootNtuple_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include <iostream>

#ifdef USING__ROOT
#include "TChain.h"
#endif

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

namespace SHERPA {
  struct RootNTupleReader_Variables {
#ifdef USING__ROOT
    static const Int_t s_kMaxParticle = 100;
    UInt_t m_id;
    Int_t m_nparticle;
    Float_t p_px[s_kMaxParticle];
    Float_t p_py[s_kMaxParticle];
    Float_t p_pz[s_kMaxParticle];
    Float_t p_E[s_kMaxParticle];
    Int_t p_kf[s_kMaxParticle];

    Double_t m_wgt,m_wgt2;
    TChain* p_f;
#endif
  };
}

RootNtuple_Reader::RootNtuple_Reader(const std::string & path,const std::string & file,int mode) :
  Event_Reader_Base(path,file), 
  m_eventmode(mode), m_evtid(0), m_subevtid(0), m_evtcnt(0), m_entries(0), m_evtpos(0)
{
  std::string filename=m_path+m_file;
  msg_Out()<<" Reading from "<<filename<<"\n";

#ifdef USING__ROOT
  p_vars = new RootNTupleReader_Variables();
  p_vars->p_f=new TChain("t3");
  p_vars->p_f->Add(filename.c_str());
  m_entries=p_vars->p_f->GetEntries();
  if(m_entries==0) {
    msg_Error()<<"ERROR: Event file "<<filename<<" does not contain any event."<<std::endl
	       <<"   Will abort the run."<<std::endl;
    abort();
  }
  p_vars->p_f->SetBranchAddress("id",&p_vars->m_id);
  p_vars->p_f->SetBranchAddress("nparticle",&p_vars->m_nparticle);
  p_vars->p_f->SetBranchAddress("px",p_vars->p_px);
  p_vars->p_f->SetBranchAddress("py",p_vars->p_py);
  p_vars->p_f->SetBranchAddress("pz",p_vars->p_pz);
  p_vars->p_f->SetBranchAddress("E",p_vars->p_E);

  p_vars->p_f->SetBranchAddress("kf",p_vars->p_kf);
  p_vars->p_f->SetBranchAddress("weight",&p_vars->m_wgt);
  p_vars->p_f->SetBranchAddress("weight2",&p_vars->m_wgt2);

  msg_Out()<<"Event mode = "<<m_eventmode<<std::endl;
#else
  msg_Error()<<"Sherpa must be linked with root to read in root files!"<<endl;
#endif
}



void RootNtuple_Reader::CloseFile() {
#ifdef USING__ROOT
  delete p_vars->p_f;
  delete p_vars;
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
  if (result==0) rpa->gen.SetNumberOfEvents(rpa->gen.NumberOfGeneratedEvents());
  
  long nev=rpa->gen.NumberOfEvents();
  if(nev==rpa->gen.NumberOfGeneratedEvents()) CloseFile();
  return result;
}

bool RootNtuple_Reader::ReadInEntry() 
{
  if (m_evtpos>=m_entries) return 0;
#ifdef USING__ROOT
  p_vars->p_f->GetEntry(m_evtpos);
  m_evtpos++;
  m_evtid=p_vars->m_id;
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
#ifdef USING__ROOT
  m_weight=p_vars->m_wgt;
  signalblob->SetWeight(m_weight);
  for (int i=0;i<p_vars->m_nparticle;i++) {
    Particle *part=new Particle(i,Flavour(abs(p_vars->p_kf[i]),
					  p_vars->p_kf[i]<0),
				Vec4D(p_vars->p_E[i],
				      p_vars->p_px[i],
				      p_vars->p_py[i],
				      p_vars->p_pz[i]));
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
  Blob         *signalblob=new Blob();
  m_weight = 0.;
  signalblob->SetType(btp::Signal_Process);
  signalblob->SetTypeSpec("NLO");
  signalblob->SetId();
  signalblob->SetPosition(Vec4D(0.,0.,0.,0.));
  signalblob->SetStatus(blob_status::code(30));
  
#ifdef USING__ROOT
  size_t currentid=m_evtid;
  while (currentid==m_evtid) {
    Vec4D *moms = new Vec4D[2+p_vars->m_nparticle];
    Flavour *flav = new Flavour[2+p_vars->m_nparticle];
    for (int i=0;i<p_vars->m_nparticle;i++) {
      moms[i+2]=Vec4D(p_vars->p_E[i],p_vars->p_px[i],
		      p_vars->p_py[i],p_vars->p_pz[i]);
      flav[i+2]=Flavour(abs(p_vars->p_kf[i]),p_vars->p_kf[i]<0);
    }
    m_nlos.push_back(new NLO_subevt(p_vars->m_nparticle+2,NULL,flav,moms));
    m_nlos.back()->m_result=p_vars->m_wgt2;
    m_nlos.back()->m_flip=0;
    m_weight+=p_vars->m_wgt2;
    if (!ReadInEntry()) m_evtid=0;
  }  
  for (size_t i=0;i<m_nlos.back()->m_n;++i) {
    Particle *part=new Particle
      (i,m_nlos.back()->p_fl[i],m_nlos.back()->p_mom[i]);
    signalblob->AddToOutParticles(part);
  }
#endif
  signalblob->SetWeight(m_weight);
  signalblob->AddData("Weight",new Blob_Data<double>(m_weight));
  signalblob->AddData("Trials",new Blob_Data<double>(1.));
  signalblob->AddData("NLO_subeventlist",new Blob_Data<NLO_subevtlist*>(&m_nlos));
  blobs->push_back(signalblob);
  m_evtcnt++;  
  return 1;
}
