#include "Fragmentation_Handler.H"

#include "Message.H"
#include "Data_Read.H"
#include "Exception.H"

using namespace SHERPA;

Fragmentation_Handler::Fragmentation_Handler(const std::string &_dir,const std::string &_file):
  m_dir(_dir), 
  m_file(_file),
  m_mode(0),
  p_lund(NULL)
{
  ATOOLS::Data_Read dr(m_dir+m_file);
  m_fragmentationmodel=dr.GetValue<std::string>("FRAGMENTATION",std::string("Lund"));
  std::string lundfile;
  if (m_fragmentationmodel==std::string("Lund")) {
    lundfile=dr.GetValue<std::string>("LUND_FILE",std::string("Lund.dat"));
    ATOOLS::msg.Events()<<"Fragmentation_Handler::Fragmentation_Handler(..): "
			<<"Initialize Lund Fragmentation according to "<<lundfile<<std::endl;
    p_lund = new Lund_Interface(m_dir,lundfile);
    m_mode = 1;
    return;
  }
  else if (m_fragmentationmodel==std::string("Off")) {
    ATOOLS::msg.Out()<<"Fragmentation_Handler::Fragmentation_Handler(..): "
		     <<"WARNING: The Fragmentation is switched off "<<std::endl;
    return;
  }
  throw(ATOOLS::Exception(ATOOLS::ex::fatal_error,std::string("Fragmentation model '")+
			  m_fragmentationmodel+std::string("' not implemented."),
			  "Fragmentation_Handler","Fragmentation_Handler"));
}
   
Fragmentation_Handler::~Fragmentation_Handler() 
{
  if (p_lund!=NULL) delete p_lund;
}

bool Fragmentation_Handler::PerformFragmentation(ATOOLS::Blob_List *bloblist,ATOOLS::Particle_List *list) 
{
  if (m_mode==0) return true;
  ATOOLS::Blob *fragmentation = new ATOOLS::Blob();
  bloblist->push_back(fragmentation);
  fragmentation->SetId(bloblist->size());
  fragmentation->SetStatus(0);
  return p_lund->Hadronize(fragmentation,bloblist,list);
}


