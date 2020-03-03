#include "SHERPA/Single_Events/Hadron_Decays.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "SHERPA/Single_Events/Decay_Handler_Base.H"
#include "METOOLS/SpinCorrelations/Amplitude2_Tensor.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Hadron_Decays::Hadron_Decays(Decay_Handler_Base* dechandler) :
  p_dechandler(dechandler)
{
  m_name      = std::string("Hadron_Decays");
  m_type      = eph::Hadronization;
}

Hadron_Decays::~Hadron_Decays()
{
}

Return_Value::code Hadron_Decays::Treat(Blob_List* bloblist)
{
  DEBUG_FUNC("bloblist->size()="<<bloblist->size());
  if(bloblist->empty()) return Return_Value::Nothing;

  bool didit(false);
  for (size_t blit(0);blit<bloblist->size();++blit) {
    Blob* blob=(*bloblist)[blit];
    if (p_dechandler && blob->Has(blob_status::needs_hadrondecays)) {
      didit = true;
      p_dechandler->SetBlobList(bloblist);
      blob->UnsetStatus(blob_status::needs_hadrondecays);
      try {
        if (p_dechandler->SpinCorr()) {
          Blob* signal=bloblist->FindFirst(btp::Signal_Process);
          if (signal) {
            METOOLS::Amplitude2_Tensor* amps(NULL);
            Blob_Data_Base* data = (*signal)["ATensor"];
            if (data) amps=data->Get<METOOLS::Amplitude2_Tensor*>();
            DEBUG_VAR(*amps);
            p_dechandler->TreatInitialBlob(blob, amps);
          }
          else {
            p_dechandler->TreatInitialBlob(blob, NULL, Particle_Vector());
          }
        }
        else p_dechandler->TreatInitialBlob(blob, NULL);
      } catch (Return_Value::code ret) {
        msg_Tracking()<<METHOD<<" Something went wrong for event: "<<*bloblist
		      <<endl<<" Will retry event."<<endl;
        return ret;
      }
    }
  }
  return (didit ? Return_Value::Success : Return_Value::Nothing);
}

void Hadron_Decays::CleanUp(const size_t & mode)
{
  if (p_dechandler) p_dechandler->CleanUp();
}

void Hadron_Decays::Finish(const std::string &) 
{
}
