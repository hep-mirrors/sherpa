#include "SHERPA/Single_Events/Minimum_Bias.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Org/Message.H"
#include <string>

using namespace SHERPA;


Minimum_Bias::Minimum_Bias(Soft_Collision_Handler_Map * schandlers) :
  p_schandlers(schandlers)
{
  msg_Out()<<METHOD<<":\n";
  for (Soft_Collision_Handler_Map::iterator scit=schandlers->begin();scit!=schandlers->end();scit++) {
    msg_Out()<<"  * "<<scit->first<<": "<<scit->second<<"\n";
  }
  m_type = eph::Perturbative;
  m_name = std::string("Minimum_Bias: ");
  if (p_schandlers->find(PDF::isr::hard_subprocess)!=p_schandlers->end())
    m_name += ( (*p_schandlers)[PDF::isr::hard_subprocess]->Soft_CollisionModel()+
		std::string(" + ") );		
  if (p_schandlers->find(PDF::isr::bunch_rescatter)!=p_schandlers->end())
    m_name += ( (*p_schandlers)[PDF::isr::bunch_rescatter]->Soft_CollisionModel()+
		std::string(" (rescatter)") );
  for (Soft_Collision_Handler_Map::iterator scit=p_schandlers->begin();
       scit!=p_schandlers->end();scit++) {
    msg_Out()<<METHOD<<"["<<this<<"], model["<<scit->first<<"] --> "<<scit->second<<"\n";
  }
}

Minimum_Bias::~Minimum_Bias() {}

ATOOLS::Return_Value::code Minimum_Bias::Treat(ATOOLS::Blob_List* blobs)
{
  for (ATOOLS::Blob_List::iterator bit=blobs->begin();bit!=blobs->end();bit++) {
    if ((*bit)->Has(ATOOLS::blob_status::needs_minBias))
      return (*p_schandlers)[PDF::isr::hard_subprocess]->GenerateMinimumBiasEvent(blobs);
    if ((*bit)->Has(ATOOLS::blob_status::needs_beamRescatter))
      return (*p_schandlers)[PDF::isr::bunch_rescatter]->GenerateBunchRescatter(blobs);
  }
  return ATOOLS::Return_Value::Nothing;
}

void Minimum_Bias::CleanUp(const size_t & mode) {
  for (Soft_Collision_Handler_Map::iterator scit=p_schandlers->begin();
       scit!=p_schandlers->end();scit++) {
    //msg_Out()<<METHOD<<" [MinBias = "<<this<<", model = "<<scit->first<<"] "
    //	     <<"--> SC_Handler = "<<scit->second<<".\n";
    scit->second->CleanUp();
  }
  msg_Out()<<"Out of "<<METHOD<<"\n";
}

void Minimum_Bias::Finish(const std::string &) {}

