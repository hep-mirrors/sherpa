#include "Hadron_Decay_Handler.H"
#include "Message.H"
#include "Random.H"
#include "Vector.H"
#include "Data_Read.H"

using namespace SHERPA;
using namespace ATOOLS;


Hadron_Decay_Handler::Hadron_Decay_Handler(std::string,std::string) :
  p_lund(NULL)
{
  msg.Error()<<"Error in Hadron_Decay_Handler::Hadron_Decay_Handler(string,string)."<<std::endl
	     <<"   This form of the Hadron_Decay_Handler is not yet available."<<std::endl
	     <<"   Abort program."<<std::endl;
  abort();
}

Hadron_Decay_Handler::Hadron_Decay_Handler(std::string _dir,std::string _file,
					   Lund_Interface * _lund) :
  m_dir(_dir), m_file(_file), p_lund(_lund)
{
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
