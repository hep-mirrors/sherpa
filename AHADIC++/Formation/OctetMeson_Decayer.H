#ifndef AHADIC_Formations_OctetMeson_Decayer_H
#define AHADIC_Formations_OctetMeson_Decayer_H

#include "AHADIC++/Tools/Singlet_Tools.H"
#include <list>

namespace AHADIC {
  class OctetMeson_Decayer : public Singlet_Tools {
  private:
    std::list<Singlet *>        * p_singlets;
    std::list<Proto_Particle *> * p_hadrons;
    int    m_offset;
    double m_kappa, m_minE;
    
    Proto_Particle * p_part1, * p_part2;
    ATOOLS::Vec4D    m_mom[3];

    bool FindOctet(Singlet * singlet,std::list<Proto_Particle *>::iterator & pit1);
    bool FixSpectator(Singlet * singlet,std::list<Proto_Particle *>::iterator & pit1);
    bool FixKinematics();
    void UpdateColouredObjectsAndAddHadron();
  public:
    OctetMeson_Decayer(std::list<Singlet *> * singlets,
		       std::list<Proto_Particle *> * hadrons);
    ~OctetMeson_Decayer() {}

    void Init() {};
    bool operator()();
  };
}

#endif
