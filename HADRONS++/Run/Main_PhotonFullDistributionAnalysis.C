/*
This is an example setup for an analysis of photon multiplicities and 
radiative energy loss as well as an angular radiation pattern analysis.
It has to include Main_FullDecay.C, and a FullDecay has to be
started with the following specialities:

  - ./PhotonAnalysis DECAYER=<kfcode> EVENTS=100000 ANALYSIS=1
*/

#include "Main_FullDecay.C"
#include "ATOOLS/Org/Shell_Tools.H"

#ifdef USING__ROOT
static TFile* rootfile;
static TH1D * photonmultiplicity;
static TH1D * decayframeenergy;
static TH1D * multipoleframeangles;
#endif

void InitialiseAnalysis()
{
#ifdef USING__ROOT
  std::string adir = "";
  ATOOLS::MakeDir("PhotonAnalysisDirectory/"+mother_flav.ShellName()+"_decays",true,493);
  rootfile = new TFile(string("PhotonAnalysisDirectory/"+
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
  std::list<Vec4D>::iterator heaviest = multipole.begin();
  if (mother_flav.Charge() != 0.) {
    Vec4D inmom = *multipole.begin();
    boost.Boost(inmom);
    rotate = Poincare(inmom,axis);
  }
  // neutral initial state: rotate such that heaviest charged final state at theta = 0
  else {
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


void FinishAnalysis()
{
#ifdef USING__ROOT
  photonmultiplicity->Write();
  decayframeenergy->Write();
  multipoleframeangles->Write();
#endif
}
