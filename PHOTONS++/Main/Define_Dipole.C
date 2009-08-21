#include "PHOTONS++/Main/Define_Dipole.H"
#include "ATOOLS/Org/Message.H"

using namespace PHOTONS;
using namespace ATOOLS;
using namespace std;

Define_Dipole::Define_Dipole(Blob * blob) :
m_success(true), m_photonsadded(false) {
  p_blob = blob;
  for (unsigned int i=0; i<(blob->GetInParticles().size()); i++)
    if (blob->InParticle(i)->Flav().Charge() == 0)
      m_neutralinparticles.push_back(blob->InParticle(i));
    else 
      m_chargedinparticles.push_back(blob->InParticle(i));
  for (unsigned int i=0; i<(blob->GetOutParticles().size()); i++)
    if (blob->OutParticle(i)->Flav().Charge() == 0)
      m_neutraloutparticles.push_back(blob->OutParticle(i));
    else
      m_chargedoutparticles.push_back(blob->OutParticle(i));

#ifdef PHOTONS_DEBUG
  msg_Info()<<m_chargedinparticles.size()<<" "
            <<m_neutralinparticles.size()<<" "
            <<m_chargedoutparticles.size()<<" "
            <<m_neutraloutparticles.size()<<endl;
#endif

  if ((m_chargedinparticles.size() == 0) &&
      (m_chargedoutparticles.size() >= 2))
    m_dtype = Dipole_Type::ff;
  else if (  (m_chargedinparticles.size() == 1) &&
           ( (m_chargedoutparticles.size() >= 2) ||
            ((m_chargedoutparticles.size() >= 1) &&
             (m_neutraloutparticles.size() >= 1))))
    m_dtype = Dipole_Type::fi;
  else if ((m_chargedinparticles.size() >= 2) &&
           (m_chargedoutparticles.size() == 0))
    m_dtype = Dipole_Type::ii;
  else
    m_dtype = Dipole_Type::unknown;

  if (msg_LevelIsDebugging()) {
    msg_Out()<<METHOD<<"(){\n";
    if      (m_dtype==Dipole_Type::ff) msg_Out()<<"  Dipole_FF(";
    else if (m_dtype==Dipole_Type::fi) msg_Out()<<"  Dipole_FI(";
    else if (m_dtype==Dipole_Type::ii) msg_Out()<<"  Dipole_II(";
    else 			       msg_Out()<<"  Undefined(";
    msg_Out()<<m_chargedinparticles.size()<<","
	     <<m_neutralinparticles.size()<<","
	     <<m_chargedoutparticles.size()<<","
	     <<m_neutraloutparticles.size()<<")\n}\n";
  }

  m_pvv.push_back(m_chargedinparticles);
  m_pvv.push_back(m_neutralinparticles);
  m_pvv.push_back(m_chargedoutparticles);
  m_pvv.push_back(m_neutraloutparticles);
}

Define_Dipole::~Define_Dipole() {
}

void Define_Dipole::AddRadiation() {
  if (m_dtype == Dipole_Type::ff) {
    Dipole_FF dipole(m_pvv);
    dipole.AddRadiation();
    m_success = dipole.DoneSuccessfully();
    m_photonsadded = dipole.AddedAnything();
    if ((m_success == true) && (m_photonsadded == true)) {
      for (int i=0; i<dipole.GetPhotonNumber(); i++) {
        p_blob->AddToOutParticles(dipole.GetPhoton(i));
      }
    }
  }
  else if (m_dtype == Dipole_Type::fi) {
    Dipole_FI dipole(m_pvv);
    dipole.AddRadiation();
    m_success = dipole.DoneSuccessfully();
    m_photonsadded = dipole.AddedAnything();
    if ((m_success == true) && (m_photonsadded == true)) {
      for (int i=0; i<dipole.GetPhotonNumber(); i++) {
        p_blob->AddToOutParticles(dipole.GetPhoton(i));
      }
    }
  }
#ifdef PHOTONS_DEBUG
  msg_Info()<<"momentum conservation in lab: "
            <<p_blob->CheckMomentumConservation()<<endl;
#endif
}
