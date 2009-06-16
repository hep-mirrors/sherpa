#include "HADRONS++/Run/Main.H"

#include "PHOTONS++/Main/Photons.H"

static Flavour mother_flav;
static Blob* ref_blob;
static Hadrons* hadrons;
static PHOTONS::Photons* photons;


void InitialiseGenerator(int argc, char *argv[])
{
  if(argc<2) {
    cout<<"Usage: ./SingleDecay <PDG_CODE>"<<endl;
    THROW(normal_exit,"you didn't specify the decaying particle by PDG code.");
  }

  small_sherpa_init(argc, argv);

  hadrons = new Hadrons(/*"../../SHERPA-1.0/Run/Decaydata/",*/ "./Decaydata/",
                        string("HadronDecays.dat"),
                        string("HadronConstants.dat") ) ;

  photons = new PHOTONS::Photons(true);
  PHOTONS::Photons::s_ircutoff=1e-3;

  mother_flav = Flavour( (kf_code) abs(ToType<int>(argv[1])) );
  mother_flav.SetStable(false);
  if(ToType<int>(argv[1])<0) mother_flav=mother_flav.Bar();

  rpa.gen.SetEcms(mother_flav.HadMass());
  msg_Info()<<"Welcome. I am decaying a "<<mother_flav<<endl;

  Particle* mother_part = new Particle( 1,mother_flav,
                                        Vec4D(mother_flav.HadMass(),0.,0.,0.) );
  mother_part->SetTime();
  mother_part->SetFinalMass(mother_flav.HadMass());
  
  ref_blob = new Blob();
  ref_blob->SetType(btp::Hadron_Decay);
  ref_blob->SetStatus(blob_status::needs_hadrondecays);
  ref_blob->AddToInParticles(mother_part);

  if(!hadrons->FillDecayBlob(ref_blob, mother_part->Momentum())) {
    THROW(normal_exit,"Hadrons::PerformDecay didn't succeed.");
  }
}


void InitialiseAnalysis()
{
}


Blob_List* GenerateEvent()
{
  Blob_List* blobs = new Blob_List();

  Blob* blob = new Blob(ref_blob);
  blob->SetStatus(blob_status::needs_extraQED);
  blobs->push_back(blob);

  photons->AddRadiation(blob);

  return blobs;
}


void AnalyseEvent(Blob_List* blobs)
{
}


void CleanUpEvent(Blob_List* blobs)
{
  blobs->Clear();
  Blob::Reset();
  Particle::Reset();
  Flow::ResetCounter();
  delete blobs;
}


void FinishGenerator()
{
  hadrons->CleanUp();
  delete hadrons;
}


void FinishAnalysis()
{
}
