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
  m_dir(_dir), m_file(_file), m_mode(0), 
#ifdef USING__Ahadic
  p_ahadic(NULL),
#endif
  p_lund(NULL)
{
  ATOOLS::Data_Read dr(m_dir+m_file);
  m_fragmentationmodel=dr.GetValue<std::string>("FRAGMENTATION",std::string("Pythiav6.214"));
  if (m_fragmentationmodel==std::string("Lund")) {
    std::string lundfile=dr.GetValue<std::string>("LUND_FILE",std::string("Lund.dat"));
    p_lund = new Lund_Interface(m_dir,lundfile,true);
    m_mode=1;
    return;
  }
#ifdef USING__Ahadic
  else if (m_fragmentationmodel==std::string("Ahadic")) {
    std::string clusterfile=dr.GetValue<std::string>("AHADIC_FILE",std::string("Cluster.dat"));
    p_ahadic = new AHADIC::Ahadic(m_dir,clusterfile,false);
    m_mode=2;
    return;
  }
#endif
  else if (m_fragmentationmodel==std::string("Off")) return;
  THROW(critical_error,"Fragmentation model not implemented.");
}
   
Fragmentation_Handler::~Fragmentation_Handler() 
{
  if (p_lund!=NULL)   { delete p_lund;   p_lund   = NULL;   }
#ifdef USING__Ahadic
  if (p_ahadic!=NULL) { delete p_ahadic; p_ahadic = NULL;   }
#endif
}

bool Fragmentation_Handler::PerformFragmentation(ATOOLS::Blob_List *bloblist,
						 ATOOLS::Particle_List *particlelist) 
{
  PROFILE_HERE;
  if (m_mode==0 || bloblist->size()==0) return true;
  switch (m_mode) {
  case 1  : return p_lund->Hadronize(bloblist,particlelist);
#ifdef USING__Ahadic
  case 2  : return p_ahadic->Hadronize(bloblist,particlelist);
#endif
  default : return false;
  }
}

