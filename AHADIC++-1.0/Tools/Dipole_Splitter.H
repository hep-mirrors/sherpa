#ifndef Dipole_Splitter_H
#define Dipole_Splitter_H

#include "Strong_Coupling.H"
#include "Cluster.H"
#include "Constituents.H"
#include "Poincare.H"
#include "Return_Value.H"

namespace AHADIC {
  class Hadronisation_Parameters;

  struct recoils {
    enum code {
      spect_remains = 0,
      energy_dep    = 1,
      mass_dep      = 2,
      unknown       = 99
    };
  };

  struct DecaySpecs { double weight,massmin,popweight; };
  typedef std::map<ATOOLS::Flavour,DecaySpecs *> FlavDecayMap;
  typedef FlavDecayMap::iterator                 FDIter;
  


  class Dipole { 
  public:
    bool   mustdecay,switched;
    double mass2,massbar2;
    Proto_Particle * triplet, * antitriplet; 
    Dipole(Proto_Particle * trip,Proto_Particle * anti);
    void Update();
  };
  inline Dipole::Dipole(Proto_Particle * trip,Proto_Particle * anti) {
    triplet   = trip; antitriplet = anti;
    mustdecay = (triplet->m_flav.IsGluon() || antitriplet->m_flav.IsGluon()); 
    mass2     = (triplet->m_mom+antitriplet->m_mom).Abs2();
    massbar2  = ATOOLS::sqr(sqrt(mass2)-(triplet->m_mass+antitriplet->m_mass));
  }
  inline void Dipole::Update() {
    mustdecay = (triplet->m_flav.IsGluon() || triplet->m_flav.IsGluon());
    mass2     = (triplet->m_mom+antitriplet->m_mom).Abs2();
    massbar2  = ATOOLS::sqr(sqrt(mass2)-(triplet->m_mass+antitriplet->m_mass));
  }
  typedef std::list<Dipole *>  DipoleList;
  typedef DipoleList::iterator DipIter;

  struct Kin_Box {
    double s,W,yint;
    double sm1,sm2,sm3,sm1_2,sm2_2,sm3_2;
    double smom1,smom3,xt1,xt3,ymin,ymax;
    double xt2kinmax,xt2max,xt2cut;
    double xt2tmp,ytmp,x1,x2,x3;
    ATOOLS::Flavour flav;
  };

  class Dipole_Splitter {
  private:
    recoils::code      m_recoils;
    Strong_Coupling  * p_as; 
    FlavDecayMap     * p_options;

    Proto_Particle   * p_spect, * p_split, * p_out1, * p_out2;
    Kin_Box            m_kin;
    double             m_evolval, m_pt2min;

    ATOOLS::Vec4D      m_mom1, m_mom2, m_mom3;

    ATOOLS::Poincare   m_booster, m_rotator;
    Constituents     * p_constituents;

    void PrepareDipole(Dipole * dip,const double);
    void DetermineSplittingMode();
    bool BuildKinematicBounds(FDIter &);
    bool SelectPT_Y(double &,double &);
    bool Veto(const double,const double);

    double FixPsi(const double x1,const double x3,const double theta13);
  public:
    Dipole_Splitter(Strong_Coupling *);
    bool FixKinematics(Dipole *);
    void UpdateSpectator();
    void GetNewParticles(Proto_Particle *& out1,Proto_Particle *& out2);

    bool Split(Dipole *,const double=-1);
    void SetOptions(FlavDecayMap *);
  };

  inline void Dipole_Splitter::GetNewParticles(Proto_Particle *& out1,
					       Proto_Particle *& out2) {
    out1 = p_out1; out2 = p_out2;
  }

  inline void Dipole_Splitter::UpdateSpectator() {
    p_spect->m_mom = m_mom3;
  }

  inline void Dipole_Splitter::SetOptions(FlavDecayMap * options) {
    p_options = options;
  }
}

#endif
