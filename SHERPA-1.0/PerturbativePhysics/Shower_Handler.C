#include "Shower_Handler.H"
#include "Tree.H"
#include "ISR_Handler.H"

using namespace SHERPA;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;

Shower_Handler::Shower_Handler(std::string _dir,std::string _file,
			       MODEL::Model_Base * _model,
			       ISR::ISR_Handler * _isr,int _maxjet) :
  m_dir(_dir), m_file(_file), p_isr(_isr), m_maxjetnumber(_maxjet),
  p_apacic(NULL)
{
  msg.Debugging()<<"Open file "<<m_dir+m_file<<endl;
  p_dataread = new Data_Read(m_dir+m_file);
  Switch::code flag = Switch::On;
  
  m_showergenerator = p_dataread->GetValue<std::string>("SHOWER_GENERATOR",std::string("Apacic"));
  m_isrshowerswitch = 0;
  if (p_isr) {
    msg.Debugging()<<"In Shower_Handler : "<<p_isr<<" "<<p_isr->On()<<" "<<p_isr->Type()<<endl;
    if (p_isr->On()>0) m_isrshowerswitch = p_dataread->GetValue<int>("ISR_SHOWER",1);
  }
  m_fsrshowerswitch = p_dataread->GetValue<int>("FSR_SHOWER",1);
  
  if (m_showergenerator==std::string("Apacic")) {
    p_apacic = new APACIC::Hard_Interface(p_isr,m_maxjetnumber,m_isrshowerswitch,m_fsrshowerswitch,p_dataread);
  }
  else {
    msg.Error()<<"Error in Shower_Handler::ReadInFile()."<<endl
    <<"   Showers needed, but no valid shower generator found !"<<endl
    <<"   Don't know, how to deal with SHOWER_GENERATOR = "<<m_showergenerator<<endl
    <<"   Abort."<<endl;
    abort();
  }

  delete p_dataread;
}


Shower_Handler::~Shower_Handler() 
{
  cout<<"in  Shower_Handler::~Shower_Handler() "<<endl;
  if (p_apacic) { delete p_apacic; p_apacic = NULL; }
  cout<<"out Shower_Handler::~Shower_Handler() "<<endl;
}


int Shower_Handler::PerformShowers(bool jetveto) {
  if (p_apacic) return p_apacic->PerformShowers(m_isrshowerswitch,m_fsrshowerswitch,jetveto);
  return 0;
}

void Shower_Handler::FillBlobs(APHYTOOLS::Blob_List * _bloblist) 
{
  if (p_apacic) {
    if (!(p_apacic->ExtractPartons(m_isrshowerswitch,m_fsrshowerswitch,_bloblist))) {
      msg.Error()<<"Error in Shower_Handler::FillBlobs()."<<endl
		 <<"   Did not succeed to fill bloblist any further."<<endl;
    }
  }
}

void Shower_Handler::CleanUp() {
  if (p_apacic) p_apacic->PrepareTrees();
}

APACIC::Tree * Shower_Handler::GetFinTree() { 
  if (p_apacic) return p_apacic->FinTree();
  msg.Error()<<"Error in Shower_Handler::FinTree()."<<endl
	     <<"   Apacic is not the shower handler."<<endl
	     <<"   Initialized "<<m_showergenerator<<". Abort run."<<endl;
  abort();
}

APACIC::Tree ** Shower_Handler::GetIniTrees() { 
  if (p_apacic) return p_apacic->IniTrees();
  msg.Error()<<"Error in Shower_Handler::FinTree()."<<endl
	     <<"   Apacic is not the shower handler."<<endl
	     <<"   Initialized "<<m_showergenerator<<". Abort run."<<endl;
  abort();
}
