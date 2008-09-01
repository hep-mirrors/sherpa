#include "Main.H"

#include "MyStrStream.H"
#include "Decay_Map.H"
#include "Hadron_Decay_Table.H"
#include "Hadron_Decay_Channel.H"

#include "Initialization_Handler.H"
#include "Event_Handler.H"
#include "Multiple_Interactions.H"
#include "Jet_Evolution.H"
#include "Hadronization.H"
#include "Hadron_Decays.H"
#include "Analysis_Phase.H"
#include "Analysis_Handler.H"
#include "Library_Loader.H"

static Flavour mother_flav;
static SHERPA::Event_Handler* p_eventhandler;
static SHERPA::Initialization_Handler* p_inithandler;

using namespace SHERPA;

void InitialiseGenerator(int argc, char *argv[])
{
  std::string statuspath;
  for (int i(1);i<argc;++i) {
    std::string cur(argv[i]);
    size_t pos(cur.find("STATUS_PATH"));
    if (pos==0 && cur.length()>11 && cur[11]=='=') {
      statuspath=cur.substr(12);
      if (statuspath=="") continue;
      if (statuspath[statuspath.length()-1]!='/') statuspath+=std::string("/");
      break;
    }
  }
  ATOOLS::s_loader = new Library_Loader();

  p_inithandler  = new Initialization_Handler(argc, argv);
  p_inithandler->InitializeTheFramework();
  if (statuspath!="") exh->ReadInStatus(statuspath);
  
  p_eventhandler  = new Event_Handler();
  p_eventhandler->AddEventPhase(new Jet_Evolution(p_inithandler->GetMatrixElementHandlers(),
                                                  p_inithandler->GetShowerHandler()));
  p_eventhandler->AddEventPhase(new Multiple_Interactions(p_inithandler->GetMIHandler()));
  p_eventhandler->AddEventPhase(new Hadronization(p_inithandler->GetFragmentationHandler()));
  p_eventhandler->AddEventPhase(new Hadron_Decays(p_inithandler->GetHadronDecayHandlers(),
                                                  p_inithandler->GetSoftPhotonHandler()));
  ANALYSIS::Analysis_Handler * ana = p_inithandler->GetSampleAnalysis();
  if (ana) p_eventhandler->AddEventPhase(new Analysis_Phase(ana));

  Data_Reader read(" ",";","!","=");
  int mother_kf(0);
  if (!read.ReadFromFile(mother_kf,"DECAYER")) {
    cout<<"Usage: ./FullDecay DECAYER=<PDG_CODE> [...]"<<endl;
    THROW(normal_exit,"you didn't specify the decaying particle by PDG code.");
  }
  mother_flav = Flavour((kf_code) abs(mother_kf));
  if(mother_kf<0) mother_flav=mother_flav.Bar();
  mother_flav.SetStable(false);
  rpa.gen.SetEcms(mother_flav.PSMass());
  msg_Info()<<"Welcome. I am decaying a "<<mother_flav<<endl;
}


Blob_List* GenerateEvent()
{
  Blob_List* blobs = p_eventhandler->GetBlobs();
  
  Particle* mother_in_part = new Particle( 1,mother_flav,Vec4D(mother_flav.PSMass(),0.,0.,0.) );
  Particle* mother_part = new Particle( 1,mother_flav,Vec4D(mother_flav.PSMass(),0.,0.,0.) );
  mother_part->SetTime();
  mother_part->SetFinalMass(mother_flav.PSMass());
  
  Blob* blob = blobs->AddBlob(btp::Signal_Process);
  blob->SetTypeSpec("1_1__"+mother_flav.IDName()+"__"+mother_flav.IDName()+"__"+
                    mother_flav.IDName()+"__"+mother_flav.IDName());
  blob->SetStatus(blob_status::needs_hadrondecays);
  blob->AddToInParticles(mother_in_part);
  mother_in_part->SetStatus(part_status::decayed);
  blob->AddToOutParticles(mother_part);

  p_eventhandler->GenerateEvent(0, false);
  if (blobs->size()<2) { // retried event
    CleanUpEvent(blobs);
    return GenerateEvent();
  }

  return blobs;
}


void CleanUpEvent(Blob_List* blobs)
{
  for(size_t i=0;i<p_eventhandler->NumberOfEventPhases();i++) {
    p_eventhandler->GetEventPhase(i)->CleanUp();
  }
  blobs->Clear();
  Blob::Reset();
  Particle::Reset();
  Flow::ResetCounter();
  ATOOLS::ran.SaveStatus();
}


void FinishGenerator()
{
  p_eventhandler->Finish();
  delete p_inithandler;
  delete ATOOLS::s_loader;
}
