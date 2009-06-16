#include "SHERPA/Single_Events/Hard_Decays.H"
#include "SHERPA/PerturbativePhysics/Hard_Decay_Handler.H"
#include "ATOOLS/Org/Message.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Hard_Decays::Hard_Decays(Hard_Decay_Handler * _dechandler) :
  p_dechandler(_dechandler)
{
  m_name      = std::string("Hard_Decays");
  m_type      = eph::Perturbative;
}

Hard_Decays::~Hard_Decays() 
{
}

Return_Value::code Hard_Decays::Treat(Blob_List * bloblist, double & weight) 
{
  return Return_Value::Nothing;
}

void Hard_Decays::CleanUp() { 
}

void Hard_Decays::Finish(const std::string &) {}
