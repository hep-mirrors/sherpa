#ifndef Particle_Qualifier_H
#define Particle_Qualifier_H 

#include "Particle.H"
#include "Getter_Function.H"

namespace ATOOLS {

  class Particle_Qualifier_Base {
  public:
    virtual bool operator() (const Particle *) const = 0;
    static void ShowQualifiers(const int mode);
  };

  typedef Primitive_Getter_Function<Particle_Qualifier_Base,std::string> 
  Particle_Qualifier_Getter;
  
  class Is_Photon : public Particle_Qualifier_Base {
  public:
    bool operator()(const Particle *) const;
  };

  class Is_Gluon : public Particle_Qualifier_Base {
  public:
    bool operator()(const Particle *) const;
  };

  class Is_Parton : public Particle_Qualifier_Base {
  public:
    bool operator()(const Particle *) const;
  };

  class Is_Charged : public Particle_Qualifier_Base {
  public:
    bool operator()(const Particle *) const;
  };

  class Is_Hadron : public Particle_Qualifier_Base {
  public:
    bool operator()(const Particle *) const;
  };

  class Is_Charged_Hadron : public Particle_Qualifier_Base {
  public:
    bool operator()(const Particle *) const;
  };

  class Is_Neutral_Hadron : public Particle_Qualifier_Base {
  public:
    bool operator()(const Particle *) const;
  };

  class Is_Final_State : public Particle_Qualifier_Base {
  public:
    bool operator()(const Particle *) const;
  };

  class Is_Neutral : public Particle_Qualifier_Base {
  public:
    bool operator()(const Particle *) const;
  };

  class Is_Charged_Pion : public Particle_Qualifier_Base {
  public:
    bool operator()(const Particle *) const;
  };

  class Is_Neutral_Pion : public Particle_Qualifier_Base {
  public:
    bool operator()(const Particle *) const;
  };

  class Is_Charged_Kaon : public Particle_Qualifier_Base {
  public:
    bool operator()(const Particle *) const;
  };

  class Is_Neutral_Kaon : public Particle_Qualifier_Base {
  public:
    bool operator()(const Particle *) const;
  };

  class Is_Charged_KStar : public Particle_Qualifier_Base {
  public:
    bool operator()(const Particle *) const;
  };

  class Is_Neutral_KStar : public Particle_Qualifier_Base {
  public:
    bool operator()(const Particle *) const;
  };

  class Is_Rho0 : public Particle_Qualifier_Base {
  public:
    bool operator()(const Particle *) const;
  };

  class Is_Eta : public Particle_Qualifier_Base {
  public:
    bool operator()(const Particle *) const;
  };

  class Is_EtaPrime : public Particle_Qualifier_Base {
  public:
    bool operator()(const Particle *) const;
  };

  class Is_Phi : public Particle_Qualifier_Base {
  public:
    bool operator()(const Particle *) const;
  };

  class Is_Omega : public Particle_Qualifier_Base {
  public:
    bool operator()(const Particle *) const;
  };

  class Is_Lambda : public Particle_Qualifier_Base {
  public:
    bool operator()(const Particle *) const;
  };

  class Is_Charged_Sigma : public Particle_Qualifier_Base {
  public:
    bool operator()(const Particle *) const;
  };

  class Is_Charged_Xi : public Particle_Qualifier_Base {
  public:
    bool operator()(const Particle *) const;
  };

  class Is_Neutral_Xi : public Particle_Qualifier_Base {
  public:
    bool operator()(const Particle *) const;
  };

  class Is_Proton_Antiproton : public Particle_Qualifier_Base {
  public:
    bool operator()(const Particle *) const;
  };

  class Is_Not_Lepton : public Particle_Qualifier_Base {
  public:
    bool operator()(const Particle *) const;
  };

  class Is_Not_Neutrino : public Particle_Qualifier_Base {
  public:
    bool operator()(const Particle *) const;
  };

  class Is_There : public Particle_Qualifier_Base {
  public:
    bool operator()(const Particle *) const;
  };

}// end of namespace ATOOLS

#endif
