/*
This is an example setup for an analysis of photon multiplicities and 
radiative energy loss as well as an angular radiation pattern analysis.
It has to include Main_FullDecay.C, and a decay with fixed configurations 
has to be started with the following specialities:

  - ./PhotonAnalysis <kfcode> EVENTS=<nevents> ANALYSIS=1
*/

#include "HADRONS++/Run/Main.H"
#include "Shell_Tools.H"

#include "ATOOLS/Org/Data_Reader.H"
#include "PHOTONS++/Main/Photons.H"

static Flavour mother_flav;
static Blob* ref_blob;
static SHERPA::Hadron_Decay_Handler* hadrons;
static PHOTONS::Photons* photons;

#ifdef USING__ROOT
static TFile* rootfile;
static TH1D * photonmultiplicity;
static TH1D * decayframeenergy;
static TH1D * multipoleframeangles;
#endif


void InitialiseGenerator(int argc, char *argv[])
{
  if(argc<2) {
    cout<<"Usage: ./SingleDecay <PDG_CODE>"<<endl;
    THROW(normal_exit,"you didn't specify the decaying particle by PDG code.");
  }

  small_sherpa_init(argc, argv);

  hadrons = new SHERPA::Hadron_Decay_Handler(".", "Fragmentation.dat");

  Data_Reader * reader = new Data_Reader(" ",";","!","=");
  reader->AddWordSeparator("\t");
  reader->SetInputPath("./");
  reader->SetInputFile("YFS.dat");
  photons = new PHOTONS::Photons(reader);

  mother_flav = Flavour( (kf_code) abs(ToType<int>(argv[1])) );
  mother_flav.SetStable(false);
  if(ToType<int>(argv[1])<0) mother_flav=mother_flav.Bar();

  rpa->gen.SetEcms(mother_flav.HadMass());
  msg_Info()<<"Welcome. I am decaying a "<<mother_flav<<endl;

  Particle* mother_part = new Particle( 1,mother_flav,
                                        Vec4D(mother_flav.HadMass(),0.,0.,0.) );
  mother_part->SetTime();
  mother_part->SetFinalMass(mother_flav.HadMass());
  
  ref_blob = new Blob();
  ref_blob->SetType(btp::Hadron_Decay);
  ref_blob->SetStatus(blob_status::needs_hadrondecays);
  ref_blob->AddToInParticles(mother_part);

  try {
    hadrons->FillOnshellDecay(ref_blob, NULL);
  } catch (Return_Value::code ret) {
    msg_Error()<<METHOD<<" Something went wrong for blob: "<<ref_blob<<endl;
    return;
  }
}


void InitialiseAnalysis()
{
#ifdef USING__ROOT
  msg_Out()<<"initialising ROOT analysis ..."<<std::endl;
  ATOOLS::MakeDir("PhotonAnalysisDirectory/SingleDecay/"+mother_flav.ShellName()+"_decays",true,493);
  rootfile = new TFile(string("PhotonAnalysisDirectory/SingleDecay/"+
                              mother_flav.ShellName()+"_decays/"+
                              mother_flav.ShellName()+"__"+
                              "DAUGHTERS"+".root").c_str(), "RECREATE");
  photonmultiplicity   = makeTH1D("photon_multiplicity","",
                                  10, -0.5, 10.5,
                                  Flavour(kf_photon).RootName()+" multiplicity","Events");
  decayframeenergy     = makeTH1D("decayframeenergy","",
                                  1000, 0., mother_flav.HadMass(),
                                  "total energy radiated in decay frame","Events");
  multipoleframeangles = makeTH1D("multipoleframeangles","",
                                  1000, 0., M_PI,
                                  "angular radiation pattern","Events");
#endif
}


Blob_List* GenerateEvent()
{
  Blob_List* blobs = new Blob_List();

  Blob* blob = new Blob(ref_blob, true);
  blob->SetStatus(blob_status::needs_extraQED);
  blobs->push_back(blob);

  photons->AddRadiation(blob);

  return blobs;
}


void AnalyseEvent(Blob_List* blobs)
{
#ifdef USING__ROOT
//   int outgoing = 1;
//   int incoming = -1;
//   Particle_List outparts = blobs->ExtractParticles(part_status::active, outgoing);

  ///////////////////////////////////////////////////////////////////////////////////
  // analyse primary decay blob, ignore subsequent decays                          //
  ///////////////////////////////////////////////////////////////////////////////////
  Blob * primarydecayblob = blobs->FindFirst(btp::Hadron_Decay);
//   msg_Out()<<"primary decay blob:"<<endl<<*primarydecayblob<<endl;
  // photon multiplicity and decay frame radiated energy (total)
  unsigned int photmult = 0;
  double       photener = 0.;
  for (int i=0; i<primarydecayblob->NOutP(); i++) {
    if ((primarydecayblob->OutParticle(i)->Flav().IsPhoton() == true) && 
        (primarydecayblob->OutParticle(i)->Info() == 'S')) {
      photmult++;
      photener = photener + primarydecayblob->OutParticle(i)->Momentum()[0];
    }
  }
  photonmultiplicity->Fill(photmult);
  if (photener != 0.)   decayframeenergy->Fill(photener);
  // multipole rest frame angles
  Vec4D multipolesum  = Vec4D(0.,0.,0.,0.);
  Vec4D axis          = Vec4D(0.,0.,0.,1.);
  std::list<Vec4D> multipole;
  std::list<Vec4D> newphot;
  for (int i=0; i<primarydecayblob->NOutP(); i++) {
    if (primarydecayblob->OutParticle(i)->Flav().Charge() != 0.) {
      multipolesum = multipolesum + primarydecayblob->OutParticle(i)->Momentum();
      multipole.push_back(primarydecayblob->OutParticle(i)->Momentum());
    }
  }
  if (primarydecayblob->InParticle(0)->Flav().Charge() != 0) {
    multipolesum = multipolesum + primarydecayblob->InParticle(0)->Momentum();
    multipole.push_front(primarydecayblob->InParticle(0)->Momentum());
  }
  Poincare boost(multipolesum);
  Poincare rotate;
  // charged initial state: rotate such that initial state at theta = 0
  if (mother_flav.Charge() != 0.) {
    Vec4D inmom = *multipole.begin();
    boost.Boost(inmom);
    rotate = Poincare(inmom,axis);
  }
  // neutral initial state: rotate such that heaviest charged final state at theta = 0
  else {
    std::list<Vec4D>::iterator heaviest = multipole.begin();
    for (std::list<Vec4D>::iterator iter=multipole.begin(); iter!=multipole.end(); iter++) {
      if (abs((iter->Abs2() - heaviest->Abs2())/(iter->Abs2() + heaviest->Abs2())) > 1E-6) {
        heaviest = iter;
      }
    }
    boost.Boost(*heaviest);
    rotate = Poincare(*heaviest,axis);
  }
  for (int i=0; i<primarydecayblob->NOutP(); i++) {
    if (primarydecayblob->OutParticle(i)->Flav().IsPhoton() == true) {
      Vec4D mom = primarydecayblob->OutParticle(i)->Momentum();
      boost.Boost(mom);
      rotate.Rotate(mom);
      double theta = acos((Vec3D(mom)*Vec3D(axis))/(Vec3D(mom).Abs()*Vec3D(axis).Abs()));
      multipoleframeangles->Fill(theta);
    }
  }  

  ///////////////////////////////////////////////////////////////////////////////////
  // inclusive analysis of whole decay chain                                       //
  ///////////////////////////////////////////////////////////////////////////////////

  // to be done ..
#endif
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
#ifdef USING__ROOT
  photonmultiplicity->Write();
  decayframeenergy->Write();
  multipoleframeangles->Write();
#endif
}
