#ifndef AHADIC_Tools_Dipole_H
#define AHADIC_Tools_Dipole_H

#include "AHADIC++/Tools/Cluster.H"
#include "AHADIC++/Tools/Constituents.H"
#include <cassert>
#include <list>

namespace AHADIC {
  class Hadronisation_Parameters;
  class Dipole;

  class Dipole {
  private:
    bool   m_mustdecay,m_switched;
    double m_mass2,m_massbar2;
    Proto_Particle * p_triplet, * p_antitriplet; 

    static size_t s_cnt;
    static std::list<Dipole * > s_actives;
  public:
    Dipole(Proto_Particle * trip,Proto_Particle * anti);
    ~Dipole();

    void Update();
    void Output();

    // simple access methods
    bool             MustDecay();
    bool             IsSwitched();
    void             SetSwitched(const bool=true);
    ATOOLS::Vec4D    Momentum() const;
    double           Mass2() const;
    double           MassBar2() const;
    Proto_Particle * Triplet();
    Proto_Particle * AntiTriplet();
    void             SetTriplet(Proto_Particle * pp);
    void             SetAntiTriplet(Proto_Particle * pp);

    bool CheckConsistency(std::ostream & s,std::string method);

    static void PrintActives(std::ostream & s) {
      s<<"   ----- Active Dipoles -----"<<std::endl;
      for (std::list<Dipole * >::iterator sit=s_actives.begin();
	   sit!=s_actives.end(); sit++) {
	s<<"   (mass^2 = "<<(*sit)->Mass2()
	 <<" /"<<(*sit)<<")"<<std::endl;
      }
      s<<"   ----- Total: "<<s_actives.size()<<" -----"<<std::endl;
    };
    
    static long int RemainingActives() {
      return s_actives.size();
    };
    
    static void DeleteActives() {
      if (s_actives.empty()) return;
      while (!s_actives.empty()) delete s_actives.front();
      s_actives.clear();
    };
  };

  inline bool   Dipole::MustDecay()  { return m_mustdecay; }
  inline bool   Dipole::IsSwitched() { return m_switched; }
  inline void   Dipole::SetSwitched(const bool sw) { m_switched=sw; }
  inline ATOOLS::Vec4D Dipole::Momentum() const {
    return p_triplet->m_mom+p_antitriplet->m_mom;
  }
  inline double Dipole::Mass2() const           { return m_mass2; }
  inline double Dipole::MassBar2() const        { return m_massbar2; }
  inline Proto_Particle * Dipole::Triplet()     { return p_triplet; }
  inline Proto_Particle * Dipole::AntiTriplet() { return p_antitriplet; }
  inline void   Dipole::SetTriplet(Proto_Particle * pp)     { p_triplet = pp; }
  inline void   Dipole::SetAntiTriplet(Proto_Particle * pp) { p_antitriplet = pp; }

  class DipoleList;
  class DipoleList : public std::list<Dipole * > {
  private:
    static std::list<DipoleList * > s_actives;
  public:
    DipoleList() {
      s_actives.push_back(this);
    }

    ~DipoleList() {
      s_actives.remove(this);
    }

    static void PrintActives(std::ostream & s) {
      s<<"   ----- Active DipoleLists -----"<<std::endl;
      for (std::list<DipoleList * >::iterator sit=s_actives.begin();
	   sit!=s_actives.end(); sit++) {
	s<<"   ("<<(*sit)->size()<<"/"<<(*sit)<<")"<<std::endl;
      }
      s<<"   ----- Total: "<<s_actives.size()<<" -----"<<std::endl;
    };
    
    static long int RemainingActives() {
      return s_actives.size();
    };
    
    static void DeleteActives() {
      if (s_actives.empty()) return;
      while (!s_actives.empty()) delete s_actives.front();
      s_actives.clear();
    };
  };

  typedef DipoleList::iterator DipIter;



  /*!
    \file 
    \brief Contains the class AHADIC::Dipole.
  */
  /*!
    \class Dipole
    \brief A wrapper around two AHADIC::Proto_Particle's to organize the decay of 
    AHADIC::Cluster's or the splitting of gluons through the AHADIC::Gluon_Decayer.
  */
  /*!
    \var bool Dipole::m_mustdecay
    \brief A flag that this dipole must decay, i.e. it contains a gluon.
  */
  /*!
    \var bool Dipole::m_switched
    \brief A flag to indicate whether the triplet or the antitriplet particle
    are the splitter or spectator in a dipole decay, respectively.
  */
  /*!
    \var double Dipole::m_mass2
    \brief The mass squared of the dipole, as given by the momenta.
  */
  /*!
    \var double Dipole::m_massbar2
    \brief The mass squared of the dipole minus the square of the sum of its
    constituent masses.
  */
  /*!
    \var Proto_Particle * Dipole::p_triplet 
    \brief A pointer to the triplet charge constituent.
  */
  /*!
    \var Proto_Particle * Dipole::p_antitriplet 
    \brief A pointer to the antitriplet charge constituent.
  */
  /*!
    \fn Dipole::Dipole(Proto_Particle * trip,Proto_Particle * anti)
    \brief
  */
  /*!
    \fn void Dipole::Update()
    \brief 
  */
  /*!
    \fn void Dipole::Output()
    \brief Output method.
  */
  /*!
    \fn bool Dipole::MustDecay()
    \brief Returns the must-decay switch, Dipole::m_mustdecay.
  */
  /*!
    \fn bool Dipole::IsSwitched()
    \brief Returns the switched flag, Dipole::m_switched, to decide which of the two 
    constituents decays.
  */
  /*!
    \fn void Dipole::SetSwitched(const bool=true)
    \brief Sets the switched flag, Dipole::m_switched, to decide which of the two 
    constituents decays.
  */
  /*!
    \fn ATOOLS::Vec4D Dipole::Momentum() const
    \brief Calculates and returns the total momentum of the dipole.
  */
  /*!
    \fn double Dipole::Mass2() const
    \brief Returns the mass squared of the dipole.
  */
  /*!
    \fn double Dipole::MassBar2() const
    \brief Returns the mass squared of the dipole minus the square of its 
    summed constituent masses.
  */
  /*!
    \fn Proto_Particle * Dipole::Triplet()
    \brief Returns the pointer to the triplet constituent, Dipole::p_triplet.
  */
  /*!
    \fn Proto_Particle * Dipole::AntiTriplet()
    \brief Returns the pointer to the antitriplet constituent, Dipole::p_antitriplet.
  */
  /*!
    \fn void Dipole::SetTriplet(Proto_Particle * pp)
    \brief Sets the pointer to the triplet constituent, Dipole::p_triplet.
  */
  /*!
    \fn void Dipole::SetAntiTriplet(Proto_Particle * pp)
    \brief Sets the pointer to the antitriplet constituent, Dipole::p_antitriplet.
  */
}

#endif
