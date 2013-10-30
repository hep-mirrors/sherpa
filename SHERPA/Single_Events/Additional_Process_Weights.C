#include "SHERPA/Single_Events/Additional_Process_Weights.H"
#include "ATOOLS/Org/Message.H"

using namespace SHERPA;
using namespace ATOOLS;

Additional_Process_Weights::Additional_Process_Weights() : 
  m_treated(false), p_signal(NULL)
{
  m_name="Additional_Process_Weights";
  m_type=eph::Perturbative;
}


Additional_Process_Weights::~Additional_Process_Weights() {}

Return_Value::code 
Additional_Process_Weights::Treat(ATOOLS::Blob_List * blobs,double & weight) 
{
  p_signal=blobs->FindFirst(btp::Signal_Process);
  if (m_treated || !p_signal) return Return_Value::Nothing;
  m_treated = true;
  msg_Out()<<"Success in "<<METHOD<<"!\n";
  weight *= 1.;
  return Return_Value::Success;
}

void Additional_Process_Weights::CleanUp(const size_t & mode) {
  m_treated = false;
}

void Additional_Process_Weights::Finish(const std::string &) {}
