#include"Hard_Decay_Handler.H"

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
  m_path(_path), m_file(_file), p_amegic(NULL)
{
  ReadInDecays();
  EvaluateWidths(_pfile,_model);
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
    cout<<buffer<<endl;
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
	cout<<"Found "<<flav<<" from "<<kfc<<"/ "<<buf<<" "<<buf.size()<<endl;
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
	if (!p_amegic) p_amegic = new AMEGIC::Amegic(m_path,_pfile,_model);
      }
    }
  }
}

std::string Hard_Decay_Handler::Name() 
{
  return std::string("");
}


