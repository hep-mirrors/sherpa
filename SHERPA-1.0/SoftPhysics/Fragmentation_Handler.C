#include "Fragmentation_Handler.H"

#include "Data_Read.H"
#include "Run_Parameter.H"
#include "Exception.H"
#include "Return_Value.H"

#ifdef PROFILE__all
#define PROFILE__Fragmentation_Handler
#endif
#ifdef PROFILE__Fragmentation_Handler
#include "prof.hh" 
#else
#define PROFILE_HERE
#endif

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Fragmentation_Handler::Fragmentation_Handler(string _dir,string _file):
  m_dir(_dir), m_file(_file), m_mode(0), 
#ifdef USING__Ahadic
  p_ahadic(NULL),
#endif
  p_lund(NULL)
{
  Data_Read dr(m_dir+m_file);
  m_fragmentationmodel=dr.GetValue<string>("FRAGMENTATION",string("Pythiav6.214"));
  if (m_fragmentationmodel==string("Lund")) {
    string lundfile=dr.GetValue<string>("LUND_FILE",string("Lund.dat"));
    p_lund = new Lund_Interface(m_dir,lundfile,true);
    m_mode=1;
    return;
  }
#ifdef USING__Ahadic
  else if (m_fragmentationmodel==string("Ahadic")) {
    string clusterfile=dr.GetValue<string>("AHADIC_FILE",string("Cluster.dat"));
    p_ahadic = new AHADIC::Ahadic(m_dir,clusterfile,true);
    m_mode=2;
    return;
  }
#endif
  else if (m_fragmentationmodel==string("Off")) return;
  THROW(critical_error,"Fragmentation model not implemented.");
}
   
Fragmentation_Handler::~Fragmentation_Handler() 
{
  if (p_lund!=NULL)   { delete p_lund;   p_lund   = NULL;   }
#ifdef USING__Ahadic
  if (p_ahadic!=NULL) { delete p_ahadic; p_ahadic = NULL;   }
#endif
}

Return_Value::code Fragmentation_Handler::PerformFragmentation(Blob_List *bloblist,
							       Particle_List *particlelist) 
{
  PROFILE_HERE;
  if (m_mode==0 || bloblist->size()==0) return Return_Value::Nothing;
  switch (m_mode) {
  case 1  : return p_lund->Hadronize(bloblist,particlelist);
#ifdef USING__Ahadic
  case 2  : return p_ahadic->Hadronize(bloblist,particlelist);
#endif
  default : return Return_Value::Nothing;
  }
}

