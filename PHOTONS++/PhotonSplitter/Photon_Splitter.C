#include "PHOTONS++/PhotonSplitter/Photon_Splitter.H"

#include "ATOOLS/Org/Scoped_Settings.H"

using namespace PHOTONS;

Photon_Splitter::Photon_Splitter() : m_on(0), m_sudakov() {}

Photon_Splitter::Photon_Splitter(int mode) : m_on(mode),
  m_sudakov(mode)
{
  msg_Debugging()<<METHOD<<"(){\n"
                 <<"  on = "<<m_on
                 <<" ,  "<<m_sudakov.Info()
                 <<"\n}"<<std::endl;
}

Photon_Splitter::~Photon_Splitter()
{}

bool Photon_Splitter::SplitPhotons(ATOOLS::Blob * blob)
{
  DEBUG_FUNC(blob->ShortProcessName());
  m_sudakov.SetNInParticles(blob->NInP());
  for (size_t i=0; i<(blob->NInP()+blob->NOutP()); ++i) {
    // access particles using blob->GetParticle(i), automatically handles indexing 
    // to test if InParticle: i < NInP()
    if (blob->GetParticle(i)->Flav().Charge() != 0) {
      m_sudakov.AddChargedParticle(blob->GetParticle(i),i);
    }
    else if (i > blob->NInP() && blob->GetParticle(i)->Info() == 'S' && blob->GetParticle(i)->Flav().Kfcode() == 22) {
      m_sudakov.AddSplitter(blob->GetParticle(i),i);
    }
  }

  m_sudakov.SetCutoff();
  bool success = m_sudakov.Run();

  if (m_sudakov.AddedAnything())
  {
    // update spectator momenta
    YFS_Particle_Vector specs = m_sudakov.GetSpectators();
    for (size_t i=0; i<specs.size(); ++i) {
      blob->GetParticle(specs[i]->Id())->SetMomentum(specs[i]->Momentum());
    }
    // remove all photons 
    size_t i = 0;
    while (i<blob->NOutP())
    {
      if (blob->OutParticle(i)->Info() == 'S' &&
          blob->OutParticle(i)->Flav().Kfcode() == kf_photon)
      {
        blob->DeleteOutParticle(blob->OutParticle(i));
      }
      else ++i;
    }
    // add photons back
    YFS_Particle_List remainingphotons = m_sudakov.GetRemainingSoftPhotons();
    for (YFS_Particle_List::iterator pvit=remainingphotons.begin();pvit!=remainingphotons.end();++pvit)
    {
      Particle *newphoton = new Particle(-1,(*pvit)->GetFlavour(),(*pvit)->Momentum(),'S');
      newphoton->SetNumber(0);
      blob->AddToOutParticles(newphoton);
    }
    // add fermions 
    // create new particles, use info 's' to identify them later on
    YFS_Particle_Vector addedparticles = m_sudakov.GetAddedParticles();
    for (YFS_Particle_Vector::iterator pvit=addedparticles.begin();pvit!=addedparticles.end();++pvit)
    {
      Particle *newparticle = new Particle(-1,(*pvit)->GetFlavour(),(*pvit)->Momentum(),'s');
      newparticle->SetNumber(0);
      blob->AddToOutParticles(newparticle);
    }
  }

  bool cleared = m_sudakov.ClearAll();
  return success && cleared;
}
