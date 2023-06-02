#include "SHERPA/Single_Events/Hard_Decays.H"
#include "SHERPA/Single_Events/Decay_Handler_Base.H"
#include "ATOOLS/Org/Message.H"
#include "METOOLS/SpinCorrelations/Amplitude2_Tensor.H"
#include "METOOLS/SpinCorrelations/Polarized_CrossSections_Handler.H"
#include "METOOLS/SpinCorrelations/PolWeights_Map.H"

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

Return_Value::code Hard_Decays::Treat(Blob_List * bloblist)
{
  if(bloblist->empty()) return Return_Value::Nothing;

  bool didit(false);
  for (size_t blit(0);blit<bloblist->size();++blit) {
    Blob* blob=(*bloblist)[blit];
    if (p_dechandler && blob->Has(blob_status::needs_harddecays)) {
      DEBUG_FUNC("Treating blob "<<blob->Id());
      didit = true;
      p_dechandler->SetBlobList(bloblist);
      if (blob->OutParticle(0)->Flav().Kfcode() == kf_instanton) {
        blob->UnsetStatus(blob_status::needs_harddecays);
        continue;
      }
      try {
        if (p_dechandler->SpinCorr()) {
          Blob* signal=bloblist->FindFirst(btp::Signal_Process);
          if (signal) {
            METOOLS::Amplitude2_Tensor* amps(NULL);
            METOOLS::Amplitude2_Tensor* prod_amps(NULL);
            bool Pol_CrossSec = p_dechandler->PolCrossSec();
            bool Spin_Coor = p_dechandler->SpinCorr();
            Blob_Data_Base* data = (*signal)["ATensor"];
            if (data) {
              amps=data->Get<METOOLS::Amplitude2_Tensor*>();
              // save production amplitude tensor before it is changed by the spin correlation algorithm
              if (Pol_CrossSec){
                prod_amps = new METOOLS::Amplitude2_Tensor(*amps);
              }
            }
            Particle_Vector outparts=blob->GetOutParticles();
            p_dechandler->TreatInitialBlob(blob, amps, outparts);
            // writing polarisation fractions to Weights_Map of signal blob

            if (Pol_CrossSec) {
              if (!Spin_Coor){
                THROW(fatal_error, "Calculation of polarized cross sections only possible together with "
                                   "spin correlations")
              }
              METOOLS::Polarized_CrossSections_Handler* pol_crosssection_handler = p_dechandler->GetPolarizationHandler();
              if (prod_amps){
                std::vector<METOOLS::PolWeights_Map*> polweights = pol_crosssection_handler->Treat(signal, prod_amps,
                                                                                                   p_dechandler->GetDecayMatrices());
                std::vector<std::string> refsystem = pol_crosssection_handler->GetRefSystems();
                Weights_Map wgtmap = bloblist->WeightsMap();
                for (size_t p(0); p<polweights.size(); ++p){
                  std::string name(refsystem[p]);
                  if (!(refsystem[p]=="Lab" || refsystem[p]=="RestFrames" || refsystem[p]=="COM" ||
                  refsystem[p]=="PPFr")){
                    name = "refsystem" + ToString(p);
                  }
                  for (auto  &e: *(polweights[p])) {
                    wgtmap["PolWeight_"+name][e.first] = e.second.real();
                  }
                  delete polweights[p];
                }
                signal->AddData("WeightsMap", new Blob_Data<Weights_Map>(wgtmap));
                delete prod_amps;
              }
              else{
                THROW(fatal_error, "No Amplitude2_Tensor for calculation of polarized cross sections found")
              }
            }
          }
        }
        else p_dechandler->TreatInitialBlob(blob, NULL, Particle_Vector());
      } catch (Return_Value::code ret) {
        return ret;
      }
      blob->UnsetStatus(blob_status::needs_harddecays);
      if (!bloblist->FourMomentumConservation()) {
	msg_Tracking()<<METHOD<<" found four momentum conservation error.\n";
	return Return_Value::New_Event;
      }
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
