#include "Hadron_Decay_Handler.H"
#include "Message.H"
#include "Random.H"
#include "Vector.H"
#include "Data_Read.H"

using namespace SHERPA;
using namespace HADRONS;
using namespace ATOOLS;
using namespace std;


Hadron_Decay_Handler::Hadron_Decay_Handler(std::string _dir,std::string _file,
					   Hadron_Decays * _hadrons) :
  m_dir(_dir), m_file(_file), m_decmodel(string("")), m_mode(-1),
  p_lund(NULL), p_hadrons(_hadrons)
{
  msg.Error()<<"Error in Hadron_Decay_Handler::Hadron_Decay_Handler(string,string)."<<std::endl
	     <<"   This form of the Hadron_Decay_Handler is not yet available."<<std::endl
	     <<"   Abort program."<<std::endl;
  abort();
}

Hadron_Decay_Handler::Hadron_Decay_Handler(std::string _dir,std::string _file,
					   Lund_Interface * _lund) :
  m_dir(_dir), m_file(_file), m_decmodel(string("")), m_mode(-1),
  p_lund(NULL), p_hadrons(NULL)
{
  Data_Read dr(m_dir+m_file);
  m_decmodel = dr.GetValue<string>("DECAYMODEL",string("Lund"));
  if (m_decmodel==string("Lund") && (_lund!=NULL)) {
    p_lund = _lund;
    m_mode = 0;
    return;
  }
  else if (m_decmodel==std::string("Hadrons")) {
    string decayfile = dr.GetValue<string>("DECAYFILE",string("HadronDecays.dat"));
    cout<<"|"<<m_dir<<"|"<<decayfile<<"|"<<endl;
    p_hadrons = new Hadron_Decays(m_dir,decayfile);
    m_mode = 1;
    return;
  }
  THROW(critical_error,"Fragmentation model not implemented.");
}


Hadron_Decay_Handler::~Hadron_Decay_Handler() 
{
}

bool Hadron_Decay_Handler::FillHadronDecayBlobs(ATOOLS::Blob_List *,ATOOLS::Particle_List *)
{
  return 1;
}

bool Hadron_Decay_Handler::ReconstructLundHadronDecays() 
{
  return 1;
}
