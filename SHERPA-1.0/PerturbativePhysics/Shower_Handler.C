#include "Shower_Handler.H"
#include "Tree.H"
#include "ISR_Handler.H"

using namespace SHERPA;
using namespace ATOOLS;

Shower_Handler::Shower_Handler(std::string _dir,std::string _file,
			       MODEL::Model_Base * _model,
			       PDF::ISR_Handler * _isr,int _maxjet) :
  m_dir(_dir), m_file(_file), m_maxjetnumber(_maxjet),
  p_apacic(NULL),p_isr_handler(_isr)
{
  p_dataread        = new Data_Read(m_dir+m_file);
  m_showergenerator = p_dataread->GetValue<std::string>("SHOWER_GENERATOR",std::string("Apacic"));
  m_isrshowerswitch = 0;
  if (_isr) {
     if (_isr->On()>0) m_isrshowerswitch = p_dataread->GetValue<int>("ISR_SHOWER",1);
  }
  m_fsrshowerswitch = p_dataread->GetValue<int>("FSR_SHOWER",1);
  if (m_isrshowerswitch && !m_fsrshowerswitch) {
    msg.Out()<<"WARNING: final state shower is switch on, since initial state shower is turned on as well."<<std::endl;
    m_fsrshowerswitch=true;
  }

  if (m_showergenerator==std::string("Apacic")) {
    p_apacic = new APACIC::Hard_Interface(_isr,_model,m_maxjetnumber,
					  m_isrshowerswitch,m_fsrshowerswitch,p_dataread);
  }
  else {
    msg.Error()<<"Error in Shower_Handler::ReadInFile()."<<std::endl
	       <<"   Showers needed, but no valid shower generator found !"<<std::endl
	       <<"   Don't know, how to deal with SHOWER_GENERATOR = "<<m_showergenerator<<std::endl
	       <<"   Abort."<<std::endl;
    abort();
  }

  delete p_dataread;
}


Shower_Handler::~Shower_Handler() 
{
  if (p_apacic) { delete p_apacic; p_apacic = NULL; }
}


int Shower_Handler::PerformShowers(int jetveto,double _x1,double _x2) {
  if (p_apacic) return p_apacic->PerformShowers(m_isrshowerswitch,m_fsrshowerswitch,jetveto,_x1,_x2);
  return 0;
}

int Shower_Handler::PerformDecayShowers(bool jetveto) {
  if (p_apacic) return p_apacic->PerformShowers(0,m_fsrshowerswitch,jetveto,1.,1.);
  return 0;
}

void Shower_Handler::FillBlobs(ATOOLS::Blob_List * _bloblist) 
{
  if (p_apacic) {
    if (!(p_apacic->ExtractPartons(m_isrshowerswitch,m_fsrshowerswitch,_bloblist))) {
      msg.Error()<<"Error in Shower_Handler::FillBlobs()."<<std::endl
		 <<"   Did not succeed to fill bloblist any further."<<std::endl;
    }
  }
}

void Shower_Handler::FillDecayBlobs(ATOOLS::Blob_List * _bloblist) 
{
  if (p_apacic) {
    if (!(p_apacic->ExtractPartons(0,m_fsrshowerswitch,_bloblist))) {
      msg.Error()<<"Error in Shower_Handler::FillBlobs()."<<std::endl
		 <<"   Did not succeed to fill bloblist any further."<<std::endl;
    }
  }
}

void Shower_Handler::CleanUp() {
  if (p_apacic) p_apacic->PrepareTrees();
}

APACIC::Tree * Shower_Handler::GetFinTree() { 
  if (p_apacic) return p_apacic->FinTree();
  msg.Error()<<"Error in Shower_Handler::FinTree()."<<std::endl
	     <<"   Apacic is not the shower handler."<<std::endl
	     <<"   Initialized "<<m_showergenerator<<". Abort run."<<std::endl;
  abort();
}

APACIC::Tree ** Shower_Handler::GetIniTrees() { 
  if (p_apacic) return p_apacic->IniTrees();
  msg.Error()<<"Error in Shower_Handler::FinTree()."<<std::endl
	     <<"   Apacic is not the shower handler."<<std::endl
	     <<"   Initialized "<<m_showergenerator<<". Abort run."<<std::endl;
  abort();
}
