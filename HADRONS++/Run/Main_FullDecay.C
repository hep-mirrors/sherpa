#include "Main.H"

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
  mother_flav=Flavour(mother_kf);
  mother_flav.SetStable(false);
  rpa->gen.SetEcms(mother_flav.HadMass());
  m_analysis = read.GetValue<int>("ANALYSIS",1);
  msg_Info()<<"Welcome. I am decaying a "<<mother_flav<<endl;
}


Blob_List* GenerateEvent()
{
  p_sherpa->GenerateOneEvent();
  return p_sherpa->GetEventHandler()->GetBlobs();
}


void CleanUpEvent(Blob_List* blobs)
{
  ATOOLS::ran->SaveStatus();
}


void FinishGenerator()
{
  p_sherpa->SummarizeRun();
  delete p_sherpa;
}
