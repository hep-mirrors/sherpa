#include "PHOTONS++/Tools/Weight_Dipole.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Message.H"

using namespace ATOOLS;
using namespace PHOTONS;
using namespace std;

// public members
Weight_Dipole::Weight_Dipole
(const Particle_Vector& olddip, const Particle_Vector& newdip,
 const Particle_Vector& phots, Dipole_Type::code dtype) {
  m_dtype       = dtype;
  m_olddipole   = olddip;
  m_newdipole   = newdip;
  m_softphotons = phots;

  CalculateWeight();
  CalculateMax();
}

Weight_Dipole::~Weight_Dipole() {
}

// private members
void Weight_Dipole::CalculateWeight() {
  double prod = 1.;
  for (unsigned int k=0; k<m_softphotons.size(); k++) {
    double sump = 0.;
    double sumq = 0.;
    for (unsigned int j=0; j<m_olddipole.size(); j++) {
      for (unsigned int i=0; i<j; i++) {
        double Zi   = m_olddipole[i]->Flav().Charge();
        double Zj   = m_olddipole[j]->Flav().Charge();
        double titj = 0.;
          if (m_newdipole[i]->ProductionBlob()
                == m_newdipole[j]->ProductionBlob())
            titj = +1.;
          else if (m_newdipole[i]->DecayBlob()
                     == m_newdipole[j]->ProductionBlob())
            titj = -1.;
          else if (m_newdipole[i]->ProductionBlob()
                     == m_newdipole[j]->DecayBlob())
            titj = -1.;
          else if (m_newdipole[i]->DecayBlob()
                     == m_newdipole[j]->DecayBlob())
            titj = +1.;
          else 
            titj = 0.;
        sump = sump + Zi*Zj*titj*SMod(m_newdipole[i]->Momentum(),
                                      m_newdipole[j]->Momentum(),
                                      m_softphotons[k]->Momentum());
        sumq = sumq + Zi*Zj*titj*SMod(m_olddipole[i]->Momentum(),
                                      m_olddipole[j]->Momentum(),
                                      m_softphotons[k]->Momentum());
      }
    }
    prod = prod * sump/sumq;
  }
  m_weight = prod;
}

void Weight_Dipole::CalculateMax() {
  m_maxweight = 1.;
}

double Weight_Dipole::SMod(const Vec4D& p1, const Vec4D& p2, const Vec4D& k) {
  return (p1/(p1*k)-p2/(p2*k)).Abs2();
}
