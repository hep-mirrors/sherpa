#include "SHERPA/Single_Events/Baryon_Recombination.H"
#include "SHERPA/SoftPhysics/Baryon_Recombination_Handler.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Baryon_Recombination::
Baryon_Recombination(Baryon_Recombination_Handler * BRhandler) :
  p_BRhandler(BRhandler)
{
  m_name      = std::string("Baryon_Recombination");
  m_type      = eph::Hadronization;
}

Baryon_Recombination::~Baryon_Recombination() {}

Return_Value::code Baryon_Recombination::Treat(Blob_List* bloblist)
{
  DEBUG_FUNC("bloblist->size()="<<bloblist->size());
  if(bloblist->empty()) return Return_Value::Nothing;
  // TODO: Fix interaction with the BRhandler, and, through it, with the
  // physics model in RECOMBINATION
  Return_Value::Success;
}

void Baryon_Recombination::CleanUp(const size_t & mode) {}

void Baryon_Recombination::Finish(const std::string &)  {}
