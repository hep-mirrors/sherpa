#include "REMNANTS/Main/Electron_Remnant.H"

#include "ATOOLS/Org/Exception.H"

using namespace REMNANTS;
using namespace ATOOLS;

Electron_Remnant::
Electron_Remnant(PDF::PDF_Base * pdf,const unsigned int & beam,const unsigned int & tag):
  Remnant_Base(pdf->Bunch(),beam,tag), p_pdfbase(pdf)
{
  // this is a *** very *** specific ordering - lepton at front, photon at back.
  // And we assume that we do not do anything with the photon, really.
  // This will be used explicitly in methods TestExtract and Fillblob
  m_constituents.push_back(p_pdfbase->Bunch());
  m_constituents.push_back(Flavour(kf_photon));
}

Electron_Remnant::
Electron_Remnant(YFS::YFS_Handler * yfs,const unsigned int & beam,const unsigned int & tag):
Remnant_Base(yfs->GetInFlav(beam),beam,tag),p_yfs(yfs)
{
  // this is a *** very *** specific ordering - lepton at front, photon at back.
  // And we assume that we do not do anything with the photon, really.
  // This will be used explicitly in methods TestExtract and Fillblob
  m_constituents.push_back(yfs->GetInFlav(beam));
  m_constituents.push_back(Flavour(kf_photon));
}

bool Electron_Remnant::FillBlob(ParticleMomMap *ktmap,const bool & copy)
{
  if (m_extracted.size()==0) {
    THROW(critical_error,"No particles extracted from intact beam.");
  }
  else if (m_extracted.size()>1) {
    THROW(critical_error,"Too many particles extracted from intact beam.");
  }  
  Particle * particle = (*m_extracted.begin()), * remnant(NULL);
  p_beamblob->AddToOutParticles(particle);
  Vec4D diff = p_beamblob->InParticle(0)->Momentum()-particle->Momentum();
  if (diff[0]>0. && diff[0]/p_beamblob->InParticle(0)->Momentum()[0]>1.e-8) {
    remnant = new Particle(-1,m_constituents.back(),diff);
    p_beamblob->AddToOutParticles(remnant);
  }
  return true;
}

bool Electron_Remnant::TestExtract(const Flavour &flav,const Vec4D &mom) {
  if (m_extracted.size()==1) {
    msg_Error()<<"Error in "<<METHOD<<" already extracted\n"
	       <<"   "<<(**m_extracted.begin())<<"\n"
	       <<"   will ignore it.\n";
  }
  if (flav!=m_constituents.front() || (p_beam->OutMomentum()[0]-mom[0])<0.) {
    if(abs(p_beam->OutMomentum()[0]-mom[0]) > 1e-6) {
      msg_Error()<<"Error in "<<METHOD<<": "<<mom<<" vs. "<<p_beam->OutMomentum()<<", "
  	       <<m_constituents.front()<<" vs. " << flav << ".\n";
      return false;
    }
  }
  return true;
}

