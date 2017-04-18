#ifndef AHADIC_Formation_Trivial_Splitter_H
#define AHADIC_Formation_Trivial_Splitter_H

#include "AHADIC++/Tools/Flavour_Selector.H"
#include "AHADIC++/Tools/KT_Selector.H"
#include "AHADIC++/Tools/Z_Selector.H"
#include "AHADIC++/Tools/Constituents.H"
#include "AHADIC++/Tools/Proto_Particle.H"
#include "ATOOLS/Math/Poincare.H"
#include <list>

namespace AHADIC {
  class Trivial_Splitter {
  private:
    Constituents         * p_constituents;
    Flavour_Selector       m_flavourselector;
    KT_Selector            m_ktselector;
    Z_Selector             m_zselector;
    
    double m_minmass;

    std::list<Proto_Particle *> * p_singlet;
    Proto_Particle * p_part1, * p_part2;
    ATOOLS::Vec4D    m_spectmom;
    double           m_Q2, m_E;
    ATOOLS::Poincare m_boost, m_rotat;

    ATOOLS::Flavour m_newflav;
    double          m_popped_mass, m_popped_mass2;
    double          m_ktmax, m_kt, m_kt2, m_phi, m_beta, m_z;
    ATOOLS::Vec4D   m_ktvec, m_q1mom, m_q2mom, m_glumom;
    

    bool InitKinematics(bool rescue=false);
    void SelectFlavour();
    bool FixTrialKinematics();
    void FixTransverseMomentum(bool rescue=false);
    void FixBetaAndZ();
    bool ConstructMomenta();
    bool CheckKinematics();

    bool Rescue();
    bool ConstructRescueMomenta();
 public:
    Trivial_Splitter();
    ~Trivial_Splitter() {}  
    void Init();
    
    bool operator()(std::list<Proto_Particle *> * singlet);
    bool operator()(Proto_Particle * part1, Proto_Particle * part2);
  };
}

#endif
