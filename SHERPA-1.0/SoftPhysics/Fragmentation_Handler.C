#include "Fragmentation_Handler.H"

#include "Data_Read.H"
#include "Run_Parameter.H"
#include "Exception.H"

#ifdef PROFILE__all
#define PROFILE__Fragmentation_Handler
#endif
#ifdef PROFILE__Fragmentation_Handler
#include "prof.hh" 
#else
#define PROFILE_HERE
#endif

using namespace SHERPA;

Fragmentation_Handler::Fragmentation_Handler(std::string _dir,std::string _file):
  m_dir(_dir), m_file(_file), m_mode(0), p_lund(NULL)
{
  ATOOLS::Data_Read dr(m_dir+m_file);
  m_fragmentationmodel=dr.GetValue<std::string>("FRAGMENTATION",std::string("Pythiav6.214"));
  if (m_fragmentationmodel==std::string("Lund")) {
    std::string lundfile=dr.GetValue<std::string>("LUND_FILE",std::string("Lund.dat"));
    p_lund = new Lund_Interface(m_dir,lundfile,true);
    m_mode=1;
    return;
  }
  else if (m_fragmentationmodel==std::string("Off")) {
    return;
  }
  throw(ATOOLS::Exception(ATOOLS::ex::critical_error,"Fragmentation model not implemented.",
			  "Fragmentation_Handler","Fragmentation_Handler"));
}
   
Fragmentation_Handler::~Fragmentation_Handler() 
{
  if (p_lund!=NULL)   { delete p_lund;   p_lund=NULL;   }
}

bool Fragmentation_Handler::PerformFragmentation(ATOOLS::Blob_List *bloblist,
						 ATOOLS::Particle_List *particlelist) 
{
  PROFILE_HERE;
  if (m_mode==0 || bloblist->size()==0) return true;
  switch (m_mode) {
  case 1  : return p_lund->Hadronize(bloblist,particlelist);
  default : return false;
  }
}

