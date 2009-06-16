#include "SHERPA/SoftPhysics/Soft_Photon_Handler.H"
#include "ATOOLS/Org/Data_Reader.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace PHOTONS;
using namespace std;


Soft_Photon_Handler::Soft_Photon_Handler(string path,string datfile) :
  m_mode(softphotons::off), p_yfs(NULL)
{
  Data_Reader * dataread = new Data_Reader(" ",";","!","=");
  dataread->AddWordSeparator("\t");
  dataread->SetInputPath(path);
  dataread->SetInputFile(datfile);

  m_mode = softphotons::code(dataread->GetValue<int>("YFS_MODE",2));
  if ((m_mode == softphotons::exp) || (m_mode == softphotons::exp_order1))
    p_yfs  = new Photons(dataread,true);
  else
    p_yfs  = new Photons(dataread,false);
  
  delete dataread;
}

Soft_Photon_Handler::~Soft_Photon_Handler() 
{
  if (p_yfs) { delete p_yfs; p_yfs = NULL; }
}

bool Soft_Photon_Handler::AddRadiation(Blob * blob)
{
  if (m_mode==softphotons::off) {
    blob->UnsetStatus(blob_status::needs_extraQED);
    return true;
  }
  p_yfs->AddRadiation(blob);
  blob->UnsetStatus(blob_status::needs_extraQED);
  return true;
}
