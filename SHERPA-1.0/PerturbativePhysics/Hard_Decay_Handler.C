#include"Hard_Decay_Handler.H"

#include "Full_Decay_Table.H"
#include "Data_Read.H"
#include "Message.H"
#include "MyStrStream.H"

#include <iostream>
#include <iomanip>
#include <sstream>

using namespace SHERPA;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Hard_Decay_Handler::Hard_Decay_Handler(std::string _path,std::string _file,std::string _pfile,
				       MODEL::Model_Base * _model) :
  m_path(_path), m_file(_file), p_amegic(NULL), m_amegicflag(0)
{
  ReadInDecays();
  EvaluateWidths(_pfile,_model);

  SetWidths(0);
}

Hard_Decay_Handler::~Hard_Decay_Handler() { }

void Hard_Decay_Handler::ReadInDecays()
{
  ifstream from((m_path+m_file).c_str());
  if (!from) {
    msg.Error()<<"Error in Amegic::InitializeProcesses : "<<endl
	       <<"   File : "<<(m_path+m_file).c_str()<<" not found ! Abort program execution."<<endl;
    abort();
  }

  char buffer[100];
  string buf,number;
  int pos,kfc;
  Decay_Table * dt = NULL;
  Flavour flav;
  for(;from;) {
    from.getline(buffer,100);
    if (buffer[0] != '\%' && strlen(buffer)>0) {
      buf    = string(buffer);
      pos    = buf.find(string("Decays :")); 
      if (pos>-1 && pos<buf.length()) {
	buf  = buf.substr(pos+8);
	while(buf.length()>1) {
	  if (buf[0]==' ') buf = buf.substr(1);
	  else {
	    pos = buf.find(string(" "));
	    if (pos>0) buf = buf.substr(0,pos);
	    MyStrStream sstream;
	    sstream<<buf;
	    sstream>>kfc;
	    break;
	  }
	}
	flav = Flavour(kf::code(int(abs(double(kfc)))));
	dt   = new Decay_Table(Flavour(kf::code(int(abs(double(kfc))))));
	m_decaytables.insert(dt);
      }
      pos     = buf.find(string("overwrite"));  
      if (pos>-1 && pos<buf.length() && (dt)) dt->SetOverwrite(); 
    }
  }
}

void Hard_Decay_Handler::EvaluateWidths(std::string _pfile,MODEL::Model_Base * _model)
{
  for (DecIt dit=m_decaytables.begin();dit!=m_decaytables.end();++dit) {
    if ((*dit)->Overwrite()) {
      if (_model->FillDecay((*dit))) { (*dit)->Output(); }
      else {
	(*dit)->Flav().SetWidth(-1.);
	if (!p_amegic) p_amegic = new AMEGIC::Amegic(m_path,_pfile,_model);
	if (!m_amegicflag) {
	  if (!p_amegic->InitializeDecays(0)) {
	    msg.Error()<<"Error in Hard_Decay_Handler::EvaluateWidths("<<_pfile<<")"<<endl
		       <<"   InitializeDecays of Amegic went wrong. Abort run."<<endl;
	    abort();
	  }
	  m_amegicflag = 1;
	}
	if (!p_amegic->GetAllDecays()->AddToDecays((*dit)->Flav())) {
	  msg.Error()<<"Error in Hard_Decay_Handler::EvaluateWidths("<<_pfile<<")"<<endl
		     <<"   Could not add "<<(*dit)->Flav()
		     <<" to list of all_decays, no vertex found. Abort run."<<endl;
	  abort();
	}
      }
    }
  }
  if (!p_amegic->GetAllDecays()->InitializeDecayTables()) {
    msg.Error()<<"Some libraries were missing ! Type make install and rerun."<<endl;
    abort();
  }
  p_amegic->GetAllDecays()->CalculateWidths();
}


void Hard_Decay_Handler::SetWidths(bool flag)
{
  AMEGIC::Full_Decay_Table * fdt;
  for (DecIt dit=m_decaytables.begin();dit!=m_decaytables.end();++dit) {
    if (!flag) {
      if ((*dit)->Overwrite()) {
	if ((*dit)->Flav().Width()<0.) {
	  fdt = p_amegic->GetAllDecays()->GetFullDecayTable((*dit)->Flav());
	  for (int i=0;i<fdt->NumberOfChannels();i++) (*dit)->AddDecayChannel(fdt->GetChannel(i));
	  (*dit)->Flav().SetWidth(fdt->Width());
	}
      }
    }
    if (flag && !(*dit)->Overwrite()) {
      fdt = p_amegic->GetAllDecays()->GetFullDecayTable((*dit)->Flav());
      for (int i=0;i<fdt->NumberOfChannels();i++) (*dit)->AddDecayChannel(fdt->GetChannel(i));
    }
    (*dit)->Output();
  }
} 

bool Hard_Decay_Handler::InitializeAllHardDecays(std::string _pfile,MODEL::Model_Base * _model) 
{
  bool newones = 0;
  for (DecIt dit=m_decaytables.begin();dit!=m_decaytables.end();++dit) {
    if (!p_amegic) p_amegic = new AMEGIC::Amegic(m_path,_pfile,_model);
    if (!m_amegicflag) {
      if (!p_amegic->InitializeDecays(0)) {
	msg.Error()<<"Error in Hard_Decay_Handler::EvaluateWidths("<<_pfile<<")"<<endl
		   <<"   InitializeDecays of Amegic went wrong. Abort run."<<endl;
	abort();
      }
      m_amegicflag = 1;
    }
    if (p_amegic->GetAllDecays()->AddToDecays((*dit)->Flav())) newones = 1;
  }
  if (newones) {
    p_amegic->GetAllDecays()->InitializeDecayTables();
    p_amegic->GetAllDecays()->CalculateWidths();
  }
  SetWidths(1);

  return 1;
}



std::string Hard_Decay_Handler::Name() 
{
  return std::string("");
}


