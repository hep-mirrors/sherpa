#include "Shower_Initialization.H"

using namespace SHERPA;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;

Shower_Initialization::Shower_Initialization(std::string _dir,std::string _file,
					     MODEL::Model_Base * _model,
					     ISR::ISR_Handler * _isr,int _maxjet) :
  m_dir(_dir), m_file(_file), p_isr(_isr), m_maxjetnumber(_maxjet) 
{
  msg.Debugging()<<"Initialized Shower_Initialization for "<<m_dir<<m_file<<endl;
  p_apacic = NULL;

  ReadInFile();
}

void Shower_Initialization::ReadInFile()
{
  msg.Debugging()<<"Open file "<<m_dir+m_file<<endl;
  p_dataread = new Data_Read(m_dir+m_file);
  Switch::code flag = Switch::On;

  std::string showergenerator = p_dataread->GetValue<std::string>("SHOWER_GENERATOR",std::string("Apacic"));
  bool isrshowerswitch        = 0;
  if (p_isr) {
    if (p_isr->On()>0) isrshowerswitch = p_dataread->GetValue<int>("ISR_SHOWER",1);
  }
  bool fsrshowerswitch        = p_dataread->GetValue<int>("FSR_SHOWER",1);

  bool okay;
  if (showergenerator==std::string("Apacic")) okay = InitializeApacic(isrshowerswitch,fsrshowerswitch);
  delete p_dataread;
  if (okay) return;

  msg.Error()<<"Error in Shower_Initialization::ReadInFile()."<<endl
	     <<"   Showers needed, but no valid shower generator found !"<<endl
	     <<"   Don't know, how to deal with SHOWER_GENERATOR = "<<showergenerator<<endl
	     <<"   Abort."<<endl;
  abort();
}


Shower_Initialization::~Shower_Initialization() {};

bool Shower_Initialization::InitializeApacic(bool _fsrshower,bool _isrshower) 
{
  p_apacic = new APACIC::Hard_Interface(p_isr,m_maxjetnumber,_fsrshower,_isrshower,p_dataread);
}


void Shower_Initialization::Output() {};
