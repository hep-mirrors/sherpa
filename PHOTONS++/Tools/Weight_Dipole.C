#include "PHOTONS++/Tools/Weight_Dipole.H"

using namespace ATOOLS;
using namespace PHOTONS;
using namespace std;

// public members
Weight_Dipole::Weight_Dipole(const Particle_Vector olddip, const Particle_Vector newdip, const Particle_Vector k, Dipole_Type::code dtype) {
  m_weight    = 1;
  m_maxweight = 1;
  m_dtype     = dtype;

  m_olddipole   = olddip;
  m_newdipole   = newdip;
  m_softphotons = k;

  CalculateWeight();
  CalculateMax();
}

Weight_Dipole::~Weight_Dipole() {
}

// private members
void Weight_Dipole::CalculateWeight() {
  double prod = 1;
  for (unsigned int k=0; k<m_softphotons.size(); k++) {
    double sump = 0;
    double sumq = 0;
    for (unsigned int j=0; j<m_olddipole.size(); j++) {
      for (unsigned int i=0; i<j; i++) {
        double Zi   = m_olddipole.at(i)->Flav().Charge();
        double Zj   = m_olddipole.at(j)->Flav().Charge();
        double titj = 0;
          if (m_newdipole.at(i)->ProductionBlob() == m_newdipole.at(j)->ProductionBlob())
            titj = +1;
          else if (m_newdipole.at(i)->DecayBlob() == m_newdipole.at(j)->ProductionBlob())
            titj = -1;
          else if (m_newdipole.at(i)->ProductionBlob() == m_newdipole.at(j)->DecayBlob())
            titj = -1;
          else if (m_newdipole.at(i)->DecayBlob() == m_newdipole.at(j)->DecayBlob())
            titj = +1;
          else 
            titj = 0;
        sump = sump + Zi*Zj*titj*SMod(m_newdipole.at(i)->Momentum(),m_newdipole.at(j)->Momentum(),m_softphotons.at(k)->Momentum());
        sumq = sumq + Zi*Zj*titj*SMod(m_olddipole.at(i)->Momentum(),m_olddipole.at(j)->Momentum(),m_softphotons.at(k)->Momentum());
      }
    }
    prod = prod * sump/sumq;
  }
  if (m_softphotons.size() > 0)
    m_weight = prod;
}

void Weight_Dipole::CalculateMax() {
  m_maxweight = 1;
}

double Weight_Dipole::SMod(const Vec4D p1, const Vec4D p2, const Vec4D k) {
  return ((p1/(p1*k)-p2/(p2*k))*(p1/(p1*k)-p2/(p2*k)));
}
