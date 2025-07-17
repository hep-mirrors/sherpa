#include "HADRON_RESCATTERING/Main/Hadron_Rescatterings.H"

using namespace HADRON_RESCATTERING;
using namespace ATOOLS;

Hadron_Rescatterings::Hadron_Rescatterings(const bool & on) :
  m_on(on), m_nfails(0), p_BB(NULL)
{ }

void Hadron_Rescatterings::Initialize() {
  p_BB = new BaryonBaryon();
}

Hadron_Rescatterings::~Hadron_Rescatterings() {
  if (p_BB) { delete p_BB; p_BB = NULL; }
}
    
Return_Value::code Hadron_Rescatterings::operator()(Blob_List *const blobs,
						    Particle_List *const) {
  return Return_Value::Nothing;
}

void Hadron_Rescatterings::Reset() {}
