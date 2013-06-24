#include "SHERPA/Single_Events/Minimum_Bias.H"
#include "ATOOLS/Phys/Blob.H"
#include <string>

using namespace SHERPA;


Minimum_Bias::Minimum_Bias(Soft_Collision_Handler * schandler) :
  p_schandler(schandler)
{
  m_name      = std::string("Minimum_Bias:")+p_schandler->Soft_CollisionModel();
  m_type      = eph::Perturbative;
}

Minimum_Bias::~Minimum_Bias() {}

ATOOLS::Return_Value::code 
Minimum_Bias::Treat(ATOOLS::Blob_List * blobs, double & weight) {
  //msg_Out()<<METHOD<<":\n"<<blobs<<" with "<<p_schandler<<"\n";
  return p_schandler->GenerateMinimumBiasEvent(blobs,weight);
}

void Minimum_Bias::CleanUp(const size_t & mode) {
  p_schandler->CleanUp(mode);
}

void Minimum_Bias::Finish(const std::string &) {}

