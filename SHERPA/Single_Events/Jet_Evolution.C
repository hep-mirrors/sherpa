#include "SHERPA/Single_Events/Jet_Evolution.H"

#include "REMNANTS/Main/Remnant_Handler.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/NLO_Types.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_MPI.H"
#include "MODEL/Main/Running_AlphaS.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace PDF;
using namespace std;

Jet_Evolution::Jet_Evolution(Shower_Handler* showerhandler,
                             REMNANTS::Remnant_Handler* remnants) :
  p_showerhandler(showerhandler), p_remnants(remnants), p_on_signalblob(NULL)
{
  m_name      = string("Jet_Evolution:")+showerhandler->ShowerGenerator();
  m_type      = eph::Perturbative;
}

Jet_Evolution::~Jet_Evolution() 
{ 
}


Return_Value::code Jet_Evolution::Treat(Blob_List * bloblist)
{
  if (p_showerhandler==NULL || !p_showerhandler->On()) return Return_Value::Nothing;
  
  bool hit(false);
  for (auto blob: *bloblist) {
    if (blob->Has(blob_status::needs_showers)) {
      hit=true;
      try {
        DefineInitialConditions(bloblist, blob);
        PerformShowers(blob);
        Event_Weights weights = GetShower()->Weights();
        if (weights.Nominal()!=1.0 || weights.Size()>1) {
          // TODO proper shower weight storage, see #225
          // MPI shower weights should not go into SP blob(?)
          Blob* weightblob=blob;
          if (p_on_signalblob) weightblob = p_on_signalblob;
          Blob_Data_Base * bdb((*weightblob)["Shower_Weights"]);
          if (!bdb) weightblob->AddData("Shower_Weights",new Blob_Data<Event_Weights>(weights));
          else {
            bdb->Get<Event_Weights>()*=weights;
          }
        }
      } catch (Return_Value::code ret) {
        CleanUp();
        return ret;
      }
      CleanUp();
    }
  }

  if (hit) return Return_Value::Success;
  
  // Capture potential problem with empty remnants here.
  // This should only happen after retrying an event has been called.  In this case
  // we find the last (and hopefully only) shower blob and extract its initiators.
  Blob * showerblob = bloblist->FindLast(btp::Shower);
  if (showerblob!=NULL && !p_remnants->ExtractShowerInitiators(showerblob)) return Return_Value::New_Event;
  return Return_Value::Nothing;
}



void Jet_Evolution::PerformShowers(Blob * blob) 
{
  CleanUp();
  
  Blob_Data_Base * bdb((*blob)["ClusterAmplitude"]);
  if (!bdb) THROW(fatal_error, "Internal Error");
  Cluster_Amplitude* ampl=bdb->Get<Cluster_Amplitude*>();
  
  int stat=GetShower()->PerformShowers(ampl);
  Reset();
  if (stat==1) {
    // No Sudakov rejection, all good
    blob->UnsetStatus(blob_status::needs_showers);
    blob->SetStatus(blob_status::needs_beams |
                    blob_status::needs_hadronization);
    blob->AddStatus(blob_status::needs_reconnections);
    // enable shower generator independent FS QED correction to ME
    blob->AddStatus(blob_status::needs_extraQED);

    blob->SetTypeSpec(GetShower()->Name());
    for (int i=0;i<blob->NInP();++i)  blob->InParticle(i)->SetStatus(part_status::decayed);
    for (int i=0;i<blob->NOutP();++i) blob->OutParticle(i)->SetStatus(part_status::decayed);


    ExtractPartons(ampl, blob);

    if (!p_remnants->ExtractShowerInitiators(blob)) {
      blob->SetStatus(blob_status::inactive);
      CleanUp();
      throw Return_Value::New_Event;
    }
  }
  else if (stat==0) {
    // Sudakov rejection
    throw Return_Value::New_Event;
  }
  else {
    DEBUG_INFO("Shower failed, will retry event.");
    throw Return_Value::New_Event;
  }
}

void Jet_Evolution::ExtractPartons(const Cluster_Amplitude* ampl, Blob* blob)
{
  DEBUG_FUNC(*blob);
  while (ampl->Prev()) ampl=ampl->Prev();
  DEBUG_VAR(*ampl);
  int beam_order[2] = {0,1};
  if (p_on_signalblob && p_on_signalblob->InParticle(0)->Momentum()[3]<0.0) {
    beam_order[0]=1; // TODO hack for Differential2
    beam_order[1]=0;
  }
  for (size_t i=0;i<ampl->NIn();++i) {
    Particle* part=
      new Particle(-1, ampl->Leg(i)->Flav().Bar(), -ampl->Leg(i)->Mom(),'I');
    part->SetNumber();
    part->SetFinalMass(ampl->MS()->Mass(part->Flav()));
    blob->AddToInParticles(part);
    part->SetFlow(1,ampl->Leg(i)->Col().m_j);
    part->SetFlow(2,ampl->Leg(i)->Col().m_i);
    part->SetBeam(beam_order[i]);
  }
  for (size_t i=ampl->NIn();i<ampl->Legs().size();++i) {
    Particle* part=
      new Particle(-1, ampl->Leg(i)->Flav(), ampl->Leg(i)->Mom(), 'F');
    part->SetNumber();
    part->SetFinalMass(ampl->MS()->Mass(part->Flav()));
    part->SetFlow(1,ampl->Leg(i)->Col().m_i);
    part->SetFlow(2,ampl->Leg(i)->Col().m_j);
    blob->AddToOutParticles(part);
  }

  /// TODO Stefan: implement this method properly
  ///   - move stuff from CSS/Dire generically here
  ///   - Check special features in CSS like SetOriginalPart which is missing(?) in Dire
  ///   - There is one additional complication: for a unified decay cascade
  ///     we might have particles with a position offset (e.g. tau decay products)
  ///     as parts of the Cluster_Amplitude.
  ///     They will never be colour-connected, so maybe we can take that into account
  ///     by creating multiple separate shower blobs (with correct positions) from the
  ///     one Cluster_Amplitude.

  /// Note: in case you need to handle errors, you can use
  // if (error) {
  //   PRINT_INFO("extract partons problem...");
  //   throw Return_Value::New_Event;
  // }
  DEBUG_VAR(*blob);
}

void Jet_Evolution::CleanUp(const size_t & mode) 
{
  p_showerhandler->CleanUp();
}

void Jet_Evolution::Reset()
{
  p_showerhandler->GetISRHandler()->Reset(0);
  p_showerhandler->GetISRHandler()->Reset(1);
}

void Jet_Evolution::DefineInitialConditions(const Blob_List* bloblist, Blob* showerblob)
{
  DEBUG_FUNC(bloblist->size());
  Reset();

  // TODO hack to keep compatibility with master
  // Any p_on_signalblob special casing should be implemented through Cluster_Amplitude
  p_on_signalblob=NULL;
  for (auto part: showerblob->GetInParticles()) {
    Blob* upstreamblob=part->ProductionBlob();
    if (upstreamblob->Type()==btp::Signal_Process ||
        upstreamblob->Type()==btp::Hard_Decay) {
      p_on_signalblob = bloblist->FindFirst(btp::Signal_Process);
    }
  }

  // because of ISRHandler->Reset all shower blobs have to be extracted again
  for (auto blob: *bloblist) {
    if (blob->Type()==btp::Shower) {
      for (size_t beam=0; beam<2; ++beam) {
        size_t cbeam=0;
        for (int i=0;i<blob->NInP();++i) {
          Particle *cur=blob->InParticle(i);
          if (!cur->Flav().Strong() || cur->ProductionBlob()) continue;
          if (cbeam==beam) {
            // TODO disabled to mimick master
            //DEBUG_INFO(*cur<<", beam = "<<beam);
            //p_remnants->Extract(cur,beam);
            continue;
          }
          ++cbeam;
        }
      }
    }
  }
}

void Jet_Evolution::Finish(const string &) 
{
}
