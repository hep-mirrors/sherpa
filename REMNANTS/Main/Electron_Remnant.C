#include "REMNANTS/Main/Electron_Remnant.H"

#include "ATOOLS/Org/Exception.H"

using namespace REMNANTS;
using namespace ATOOLS;

Electron_Remnant::Electron_Remnant(PDF::PDF_Base* pdf, const unsigned int& beam,
                                   const unsigned int& tag)
    : Remnant_Base(pdf->Bunch(), beam, tag), p_pdfbase(pdf), p_yfs(nullptr)
{
  // this is a *** very *** specific ordering - lepton at front, photon at back.
  // And we assume that we do not do anything with the photon, really.
  // This will be used explicitly in methods TestExtract and Fillblob
  m_constituents.push_back(p_pdfbase->Bunch());
  m_constituents.push_back(Flavour(kf_photon));
}

Electron_Remnant::
Electron_Remnant(YFS::YFS_Handler * yfs,const unsigned int & beam,const unsigned int & tag):
Remnant_Base(yfs->GetInFlav(beam),beam,tag),p_pdfbase(nullptr),p_yfs(yfs)
{
  // this is a *** very *** specific ordering - lepton at front, photon at back.
  // And we assume that we do not do anything with the photon, really.
  // This will be used explicitly in methods TestExtract and Fillblob
  m_constituents.push_back(yfs->GetInFlav(beam));
  m_constituents.push_back(Flavour(kf_photon));
}

bool Electron_Remnant::Extract(Particle* parton, Colour_Generator* colours,
                               const double& spair)
{
  const bool wholebeam(parton->Flav() == m_beamflav &&
                       IsEqual(parton->Momentum(), IncomingMomentum(), 1.e-8));
  const bool ok(Remnant_Base::Extract(parton, colours, spair));
  if (ok && wholebeam) p_wholebeam = parton;
  return ok;
}

void Electron_Remnant::Reset(const bool& resc, const bool& DIS)
{
  p_wholebeam = nullptr;
  Remnant_Base::Reset(resc, DIS);
}

bool Electron_Remnant::FillBlob(Colour_Generator* colours, ParticleMomMap* ktmap, const bool& copy)
{
  // The lepton took the whole beam (no ISR photon on this side): it leaves the
  // beam blob as it is, with nothing left over for a photon remnant.
  if (m_extracted.empty() && p_wholebeam != nullptr) {
    p_beamblob->AddToOutParticles(p_wholebeam);
    return true;
  }
  if (m_extracted.size() != 1) {
    THROW(critical_error, "None or too many particles extracted from intact beam.");
  }
  Particle* particle = (*m_extracted.begin());
  p_beamblob->AddToOutParticles(particle);
  Vec4D diff = p_beamblob->InParticle(0)->Momentum()-particle->Momentum();
  if (diff[0]>0. && diff[0]/p_beamblob->InParticle(0)->Momentum()[0]>1.e-8) {
    p_beamblob->AddToOutParticles(new Particle(-1, m_constituents.back(), diff));
  }
  return true;
}

bool Electron_Remnant::TestExtract(const Flavour &flav,const Vec4D &mom,
                                   const double &spair) {
  if (m_extracted.size()==1) {
    msg_Error() << "Error in " << METHOD << " already extracted\n"
                << "   " << (**m_extracted.begin()) << "\n"
                << "   will ignore it.\n";
  }
  if (flav != m_constituents.front() ||
      mom[0] > (1. + 1.e-6) * p_beam->OutMomentum()[0]) {
    msg_Error() << "Error in " << METHOD << ": parton " << mom << " vs. beam "
                << p_beam->OutMomentum() << ", " << m_constituents.front() << " vs. " << flav << ".\n";
    return false;
  }
  return true;
}

