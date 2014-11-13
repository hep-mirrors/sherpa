#include "SHERPA/SoftPhysics/Soft_Photon_Handler.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Flavour.H"
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
  m_photonsadded(false),
  m_name(""), m_mode(softphotons::off),
  p_yfs(NULL), p_clusterer(NULL), p_mehandler(meh)
{
  Data_Reader * dataread = new Data_Reader(" ",";","!","=");
  dataread->AddComment("#");
  dataread->AddWordSeparator("\t");
  dataread->SetInputPath(path);
  dataread->SetInputFile(datfile);

  m_mode = softphotons::code(dataread->GetValue<int>("YFS_MODE",2));
  if (m_mode != softphotons::off) {
    p_yfs       = new PHOTONS::Photons(dataread);
    p_clusterer = new Resonance_Finder(dataread,meh);
    m_name      = p_yfs->Name();
  }

  delete dataread;
}

Soft_Photon_Handler::~Soft_Photon_Handler() 
{
  if (p_yfs)       { delete p_yfs;       p_yfs = NULL; }
  if (p_clusterer) { delete p_clusterer; p_clusterer = NULL; }
}

bool Soft_Photon_Handler::AddRadiation(Blob * blob)
{
  if (m_mode==softphotons::off) {
    blob->UnsetStatus(blob_status::needs_extraQED);
    return true;
  }
  p_yfs->AddRadiation(blob);
  blob->UnsetStatus(blob_status::needs_extraQED);
  m_photonsadded=p_yfs->AddedAnything();
  return p_yfs->DoneSuccessfully();
}

bool Soft_Photon_Handler::AddRadiation(Particle_Vector& leps, Blob_Vector& blobs)
{
  // build effective verteces for resonant production
  // use subprocess infos if possible
  p_clusterer->BuildResonantBlobs(leps,blobs);
  // add radiation
  for (Blob_Vector::iterator it=blobs.begin();it!=blobs.end();++it) {
    // do nothing if no resonance determined
    if ((*it)->InParticle(0)->Flav().Kfcode()!=kf_none) {
      (*it)->SetStatus(blob_status::needs_extraQED);
      if (!AddRadiation(*it)) return false;
    }
  }
  for (Blob_Vector::iterator it=blobs.begin();it!=blobs.end();++it) {
    msg_Debugging()<<**it<<endl;
    (*it)->DeleteInParticles();
  }
  return true;
}

