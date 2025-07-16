#include "SHERPA/Single_Events/Minimum_Bias.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Org/Message.H"
#include <string>

using namespace SHERPA;


Minimum_Bias::Minimum_Bias(Soft_Collision_Handler_Map * schandlers) : m_on(false)
{
  m_type   = eph::Perturbative;
  m_name   = std::string("Minimum_Bias: ");
  if (schandlers->find(PDF::isr::hard_subprocess)!=schandlers->end() &&
      (*schandlers)[PDF::isr::hard_subprocess]!=NULL) {
    p_schandler = (*schandlers)[PDF::isr::hard_subprocess];
    m_name     += p_schandler->Soft_CollisionModel();
    m_on        = true;
  }
  else m_name += "None";
}

Minimum_Bias::~Minimum_Bias() {}

ATOOLS::Return_Value::code Minimum_Bias::Treat(ATOOLS::Blob_List* blobs)
{
  if (m_on) {
    for (ATOOLS::Blob_List::iterator bit=blobs->begin();bit!=blobs->end();bit++) {
      if ((*bit)->Has(ATOOLS::blob_status::needs_minBias))
	return p_schandler->GenerateMinimumBiasEvent(blobs);
    }
  }
  return ATOOLS::Return_Value::Nothing;
}

void Minimum_Bias::CleanUp(const size_t & mode) { if (m_on) p_schandler->CleanUp(); }

void Minimum_Bias::Finish(const std::string &) {}

