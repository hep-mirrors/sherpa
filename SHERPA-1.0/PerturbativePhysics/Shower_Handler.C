#include "Shower_Handler.H"
#include "Tree.H"
#include "ISR_Handler.H"
#include "Message.H"
#include "Data_Read.H"
#include "MyStrStream.H"

using namespace SHERPA;
using namespace ATOOLS;

Shower_Handler::Shower_Handler(std::string dir,std::string file,
			       MODEL::Model_Base * _model,
			       PDF::ISR_Handler * _isr,int _maxjet) :
  m_maxjetnumber(_maxjet), m_showermi(true), 
  p_apacic(NULL), 
  p_jf(NULL), p_isr_handler(_isr)
{
  Data_Read dataread(dir+file);
  m_showergenerator = dataread.GetValue<std::string>("SHOWER_GENERATOR",std::string("Apacic"));
  m_isrshowerswitch = 0;
  if (_isr) {
    if (_isr->On()>0) m_isrshowerswitch = dataread.GetValue<int>("ISR_SHOWER",1);
  }
  m_fsrshowerswitch = dataread.GetValue<int>("FSR_SHOWER",1);
  if (m_isrshowerswitch && !m_fsrshowerswitch) {
    msg.Out()<<"WARNING in Shower_Handler : "<<std::endl
	     <<"   final state shower is switched on, since initial state shower is turned on as well."<<std::endl;
    m_fsrshowerswitch=true;
  }
  m_showermi = dataread.GetValue<int>("SHOWER_MI",1);
  rpa.gen.SetVariable("SHOWER_MODE",dataread.GetValue<std::string>("SHOWER_MODE","0"));
  int type(4);
  if (rpa.gen.Beam1().IsLepton() && rpa.gen.Beam2().IsLepton())          type = 1;
  else if ((!rpa.gen.Beam1().IsLepton() && !rpa.gen.Beam2().IsLepton())) type = 4;
  else type = 4;

  p_jf = new Jet_Finder(rpa.gen.Variable("Y_CUT"),type,false);
  p_jf->SetDeltaR(ToType<double>(rpa.gen.Variable("DELTA_R")));

  if (m_showergenerator==std::string("Apacic")) {
    p_apacic = new APACIC::Apacic(_isr,_model,p_jf,&dataread);
  }
  else {
    msg.Error()<<"Error in Shower_Handler::ReadInFile()."<<std::endl
	       <<"   Showers needed, but no valid shower generator found !"<<std::endl
	       <<"   Don't know, how to deal with SHOWER_GENERATOR = "<<m_showergenerator<<std::endl
	       <<"   Abort."<<std::endl;
    abort();
  }
}


Shower_Handler::~Shower_Handler() 
{
  if (p_apacic) { delete p_apacic; p_apacic = NULL; }
  if (p_jf) delete p_jf;
}


int Shower_Handler::PerformShowers(int jetveto,int losejv,double _x1,double _x2) {
  if (p_apacic) return p_apacic->PerformShowers(jetveto,losejv,_x1,_x2);
  return 0;
}

int Shower_Handler::PerformDecayShowers(bool jetveto) {
  if (p_apacic) return p_apacic->PerformShowers(jetveto,1,1.,1.);
  return 0;
}

void Shower_Handler::FillBlobs(ATOOLS::Blob_List * _bloblist) 
{
  if (p_apacic) {
    if (!(p_apacic->ExtractPartons(m_isrshowerswitch,m_fsrshowerswitch,_bloblist))) {
      msg.Error()<<"Error in Shower_Handler::FillBlobs()."<<std::endl
		 <<"   Did not succeed to fill bloblist any further."<<std::endl;
    }
//     std::cout<<"------------------------------------------------------------------------------"<<std::endl
// 	     <<(*_bloblist)
// 	     <<"------------------------------------------------------------------------------"<<std::endl;
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

