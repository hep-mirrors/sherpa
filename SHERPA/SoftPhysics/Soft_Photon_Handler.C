#include "SHERPA/SoftPhysics/Soft_Photon_Handler.H"
#include "ATOOLS/Org/Default_Reader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Tensor.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Momenta_Stretcher.H"
#include "ATOOLS/Phys/Particle.H"
#include "PHOTONS++/Main/Photons.H"
#include "SHERPA/PerturbativePhysics/Matrix_Element_Handler.H"
#include "SHERPA/SoftPhysics/Resonance_Finder.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace PHASIC;
using namespace std;


Soft_Photon_Handler::Soft_Photon_Handler(string path,string datfile,
                                         Matrix_Element_Handler * meh) :
  m_photonsadded(false), m_name(""),
  p_yfs(NULL), p_clusterer(NULL), p_mehandler(meh)
{
  Default_Reader reader;
  reader.SetInputPath(path);
  reader.SetInputFile(datfile);

  p_yfs       = new PHOTONS::Photons(&reader);
  p_clusterer = new Resonance_Finder(&reader,meh);
  m_name      = p_yfs->Name();
}

Soft_Photon_Handler::~Soft_Photon_Handler() 
{
  if (p_yfs)       { delete p_yfs;       p_yfs = NULL; }
  if (p_clusterer) { delete p_clusterer; p_clusterer = NULL; }
}

bool Soft_Photon_Handler::AddRadiation(Blob * blob)
{
  DEBUG_FUNC("");
  p_yfs->AddRadiation(blob);
  blob->UnsetStatus(blob_status::needs_extraQED);
  m_photonsadded = p_yfs->AddedAnything();
  if (m_photonsadded)
    for (size_t i(0); i<blob->NOutP(); ++i)
      if (blob->OutParticle(i)->DecayBlob())
        BoostDecayBlob(blob->OutParticle(i)->DecayBlob());
  return p_yfs->DoneSuccessfully();
}

bool Soft_Photon_Handler::AddRadiation(Particle_Vector& leps, Blob_Vector& blobs)
{
  // build effective verteces for resonant production
  // use subprocess infos if possible
  p_clusterer->BuildResonantBlobs(leps,blobs);
  bool photonsadded(false);
  // add radiation
  for (Blob_Vector::iterator it=blobs.begin();it!=blobs.end();++it) {
    // do nothing if no resonance determined
    if ((*it)->InParticle(0)->Flav().Kfcode()!=kf_none) {
      (*it)->SetStatus(blob_status::needs_extraQED);
      if (!AddRadiation(*it)) return false;
      photonsadded+=m_photonsadded;
    }
  }
  m_photonsadded=photonsadded;
  for (Blob_Vector::iterator it=blobs.begin();it!=blobs.end();++it) {
    msg_Debugging()<<**it<<endl;
    (*it)->DeleteInParticles();
  }
  return true;
}

void Soft_Photon_Handler::BoostDecayBlob(Blob * blob)
{
  DEBUG_FUNC("");
  //msg_Out()<<*blob<<std::endl;
  const Vec4D& P((*blob)["p_original"]->Get<Vec4D>());
  const Vec4D& Pt(blob->InParticle(0)->Momentum());
  const Vec4D e(P-Pt);
  msg_Debugging()<<"P-Pt="<<e<<" ["<<e.Mass()<<"]"<<std::endl;
  //msg_Out()<<"Before builidng tensors.\n";
  const Lorentz_Ten2D lambda(MetricTensor()-2.*BuildTensor(e,e)/e.Abs2());
  msg_Debugging()<<"\\Lambda="<<std::endl<<lambda<<std::endl;
  for (size_t i(0);i<blob->NOutP();++i) {
    Vec4D mom(blob->OutParticle(i)->Momentum());
    msg_Debugging()<<blob->OutParticle(i)->Flav().IDName()<<" "
                   <<mom<<" ["<<mom.Mass()<<"] -> ";
    mom=Contraction(lambda,2,mom);
    blob->OutParticle(i)->SetMomentum(mom);
    msg_Debugging()<<mom<<" ["<<mom.Mass()<<"]"<<std::endl;
  }
  //msg_Out()<<"Before checkong onshell-ness.\n";
  CheckOnshellness(blob);
  if (msg_LevelIsDebugging()) {
    for (size_t i(0);i<blob->NOutP();++i) {
      Vec4D mom(blob->OutParticle(i)->Momentum());
      msg_Debugging()<<blob->OutParticle(i)->Flav().IDName()<<" "
                     <<mom<<" ["<<mom.Mass()<<"]"<<std::endl;
    }
  }
  msg_Debugging()<<"Momentum conservation in decay blob of "
                 <<blob->InParticle(0)->Flav()<<": "
                 <<blob->CheckMomentumConservation()<<std::endl;
}

bool Soft_Photon_Handler::CheckOnshellness(Blob* blob)
{
  std::vector<double> masses;
  bool allonshell(true);
  double accu(sqrt(Accu()));
  for (size_t i(0);i<blob->NOutP();++i) {
    masses.push_back(blob->OutParticle(i)->FinalMass());
    if (allonshell &&
        !IsEqual(blob->OutParticle(i)->Momentum().Abs2(),
                 sqr(blob->OutParticle(i)->FinalMass()),accu)) allonshell=false;
  }
  msg_Debugging()<<"masses="<<masses<<std::endl;
  if (allonshell) return true;
  msg_Debugging()<<"need to put on-shell"<<std::endl;
  Momenta_Stretcher momstretch;
  momstretch.StretchMomenta(blob->GetOutParticles(),masses);
  return false;
}
