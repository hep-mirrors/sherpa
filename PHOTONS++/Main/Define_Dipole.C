#include "PHOTONS++/Main/Define_Dipole.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Particle.H"
#include "PHOTONS++/Tools/Dipole_FF.H"
#include "PHOTONS++/Tools/Dipole_FI.H"
#include "PHOTONS++/Tools/Dress_Blob_Base.H"

using namespace PHOTONS;
using namespace ATOOLS;
using namespace std;

Define_Dipole::Define_Dipole(Blob * blob) :
m_success(true), m_photonsadded(false) {
  DEBUG_FUNC("");
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
    if      (m_dtype==Dipole_Type::ff) msg_Out()<<"  Dipole_FF(";
    else if (m_dtype==Dipole_Type::fi) msg_Out()<<"  Dipole_FI(";
    else if (m_dtype==Dipole_Type::ii) msg_Out()<<"  Dipole_II(";
    else 			       msg_Out()<<"  Undefined(";
    msg_Out()<<m_chargedinparticles.size()<<","
	     <<m_neutralinparticles.size()<<","
	     <<m_chargedoutparticles.size()<<","
	     <<m_neutraloutparticles.size()<<")\n";
  }

  m_pvv.push_back(m_chargedinparticles);
  m_pvv.push_back(m_neutralinparticles);
  m_pvv.push_back(m_chargedoutparticles);
  m_pvv.push_back(m_neutraloutparticles);
}

Define_Dipole::~Define_Dipole() {
}

void Define_Dipole::AddRadiation() {
  Dress_Blob_Base * dipole(NULL);
  if (m_dtype == Dipole_Type::ff) {
    dipole = new Dipole_FF(m_pvv);
  }
  else if (m_dtype == Dipole_Type::fi) {
    dipole = new Dipole_FI(m_pvv);
  }
  else return;
  // treat and add to blob, reset photon's number
  if (dipole) {
    dipole->AddRadiation();
    m_success = dipole->DoneSuccessfully();
    m_photonsadded = dipole->AddedAnything();
    Particle_Vector photons;
    if ((m_success == true) && (m_photonsadded == true)) {
      for (int i=0; i<dipole->GetPhotonNumber(); i++) {
	photons.push_back(dipole->GetPhoton(i));
      }
    }
    delete dipole;
    if ((m_success == true) && (m_photonsadded == true)) {
      for (Particle_Vector::iterator pvit=photons.begin();
	   pvit!=photons.end();++pvit) {
	(*pvit)->SetNumber(0);
	p_blob->AddToOutParticles(*pvit);
      }
    }
  }
  else m_success = false;
}

bool Define_Dipole::CheckMasses()
{
  for (size_t i(0);i<m_chargedinparticles.size();++i)
    if (m_chargedinparticles[i]->FinalMass()==0.) return false;
  for (size_t i(0);i<m_chargedoutparticles.size();++i)
    if (m_chargedoutparticles[i]->FinalMass()==0.) return false;
  return true;
}
