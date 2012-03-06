#include "SHERPA/Single_Events/Hard_Decays.H"
#include "SHERPA/Single_Events/Decay_Handler_Base.H"
#include "ATOOLS/Org/Message.H"
#include "METOOLS/SpinCorrelations/Amplitude2_Tensor.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Hard_Decays::Hard_Decays(Decay_Handler_Base* dechandler) :
  p_dechandler(dechandler)
{
  m_name      = std::string("Hard_Decays");
  m_type      = eph::Perturbative;
}

Hard_Decays::~Hard_Decays() 
{
}

Return_Value::code Hard_Decays::Treat(Blob_List * bloblist, double & weight)
{
  if(bloblist->empty()) return Return_Value::Nothing;

  bool didit(false);
  for (size_t blit(0);blit<bloblist->size();++blit) {
    Blob* blob=(*bloblist)[blit];
    if (p_dechandler && blob->Has(blob_status::needs_harddecays)) {
      DEBUG_FUNC("Treating blob "<<blob->Id());
      didit = true;
      p_dechandler->SetBlobList(bloblist);
      try {
        if (p_dechandler->SpinCorr()) {
          Blob* signal=bloblist->FindFirst(btp::Signal_Process);
          if (signal) {
            METOOLS::Amplitude2_Tensor* amps(NULL);
            Blob_Data_Base* data = (*signal)["ATensor"];
            if (data) amps=data->Get<METOOLS::Amplitude2_Tensor*>();
            Particle_Vector outparts=blob->GetOutParticles();
            p_dechandler->TreatInitialBlob(blob, amps, outparts);
          }
        }
        else p_dechandler->TreatInitialBlob(blob, NULL, Particle_Vector());
      } catch (Return_Value::code ret) {
        return ret;
      }
      blob->SetStatus(blob_status::needs_showers);
    }
  }
  return (didit ? Return_Value::Success : Return_Value::Nothing);
}

void Hard_Decays::CleanUp(const size_t & mode)
{
  if (p_dechandler) p_dechandler->CleanUp();
}

void Hard_Decays::Finish(const std::string &)
{
}
