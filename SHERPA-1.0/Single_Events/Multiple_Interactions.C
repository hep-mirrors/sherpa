#include "Multiple_Interactions.H"
#include "Matrix_Element_Handler.H"

#ifdef PROFILE__Multiple_Interactions
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif

using namespace SHERPA;
using namespace ATOOLS;

Multiple_Interactions::Multiple_Interactions(MI_Handler *_p_mihandler):
  m_diced(false),
  p_mihandler(_p_mihandler)
{
  m_name = std::string("Multiple_Interactions:")+p_mihandler->Name();
  m_type = eph::Perturbative;
}

Multiple_Interactions::~Multiple_Interactions() {}

bool Multiple_Interactions::Treat(ATOOLS::Blob_List *bloblist,double &weight)
{
  PROFILE_HERE;
  if (p_mihandler->Type()==MI_Handler::None || m_diced) return false;
  if (bloblist->empty()) {
    ATOOLS::msg.Error()<<"Multiple_Interactions::Treat(): "
		       <<"Incoming blob list is empty!"<<std::endl;
    return false;
  }
  bool check=false;
  for (Blob_Iterator bit=bloblist->begin();bit!=bloblist->end();++bit) {
    if ((*bit)->Type()==btp::Hard_Collision) check=true;
  }
  if (check) return false;
  ATOOLS::Blob *mepsblob=NULL, *signalblob=NULL;
  for (Blob_Iterator bit=bloblist->begin();bit!=bloblist->end();++bit) {
    if ((*bit)->Type()==btp::ME_PS_Interface_FS) mepsblob=*bit;
    else if ((*bit)->Type()==btp::Signal_Process) signalblob=*bit;
  }
  if (mepsblob==NULL) return false;
  m_diced=true;
  m_pperpmax=1.0e37; 
  m_ecmsmax=ATOOLS::rpa.gen.Ecms();
  double val=0.0;
  m_xmax[1]=m_xmax[0]=1.0;
  ATOOLS::Vec4D cur, ptot;
  for (int i=0;i<mepsblob->NInP();++i) {
    ptot+=cur=mepsblob->InParticle(i)->Momentum();
    val=ATOOLS::Max(val,cur.PPerp());
  }
  for (int i=0;i<signalblob->NInP();++i) {
    m_xmax[i]-=2.0*signalblob->InParticle(i)->Momentum()[0]/ATOOLS::rpa.gen.Ecms();
  }
  m_ecmsmax-=sqrt(ptot.Abs2());
  m_pperpmax=ATOOLS::Min(m_pperpmax,val);
  p_mihandler->SetScaleMax(m_pperpmax,0);
  p_mihandler->SetScaleMax(m_ecmsmax,1);
  p_mihandler->SetScaleMax(m_xmax[0],2);
  p_mihandler->SetScaleMax(m_xmax[1],3);
  size_t size=bloblist->size();
  p_mihandler->GenerateHardEvent(bloblist);
  return (size!=bloblist->size());
}

void Multiple_Interactions::Finish(const std::string &) 
{
}

void Multiple_Interactions::CleanUp() 
{ 
  m_diced=false;
}

