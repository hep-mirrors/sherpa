#include "HADRONS++/Run/Main.H"

#include "ATOOLS/Org/MyStrStream.H"
#include "HADRONS++/Main/Hadron_Decay_Map.H"
#include "HADRONS++/Main/Hadron_Decay_Table.H"
#include "HADRONS++/Main/Hadron_Decay_Channel.H"

#include "SHERPA/Main/Sherpa.H"
#include "SHERPA/Single_Events/Event_Handler.H"
#include "ATOOLS/Math/Random.H"

static Flavour mother_flav;
static SHERPA::Sherpa* p_sherpa;

using namespace SHERPA;

void InitialiseGenerator(int argc, char *argv[])
{
  p_sherpa = new Sherpa();
  p_sherpa->InitializeTheRun(argc,argv);
  p_sherpa->InitializeTheEventHandler();

  Data_Reader read(" ",";","!","=");
  int mother_kf(0);
  if (!read.ReadFromFile(mother_kf,"DECAYER")) {
    cout<<"Usage: ./FullDecay DECAYER=<PDG_CODE> [...]"<<endl;
    THROW(normal_exit,"you didn't specify the decaying particle by PDG code.");
  }
  mother_flav.FromHepEvt(mother_kf);
  mother_flav.SetStable(false);
  rpa.gen.SetEcms(mother_flav.HadMass());
  msg_Info()<<"Welcome. I am decaying a "<<mother_flav<<endl;
}


Blob_List* GenerateEvent()
{
  Blob_List* blobs=p_sherpa->GetEventHandler()->GetBlobs();
  if (!blobs->empty()) THROW(fatal_error, "Bloblist not empty.");

  Vec4D mom(mother_flav.HadMass(), 0., 0., 0.);
  Particle* mother_in_part=new Particle(1, mother_flav, mom);
  Particle* mother_part=new Particle(1, mother_flav, mom);
  mother_part->SetTime();
  mother_part->SetFinalMass(mother_flav.HadMass());
  mother_in_part->SetStatus(part_status::decayed);
  
  Blob* blob = blobs->AddBlob(btp::Hadron_Decay);
  blob->SetStatus(blob_status::needs_hadrondecays);
  blob->AddToInParticles(mother_in_part);
  blob->AddToOutParticles(mother_part);

  p_sherpa->GenerateOneEvent(false);
  return blobs;
}


void CleanUpEvent(Blob_List* blobs)
{
  p_sherpa->GetEventHandler()->Reset();
  ATOOLS::ran.SaveStatus();
}


void FinishGenerator()
{
  p_sherpa->SummarizeRun();
  delete p_sherpa;
}
