#ifndef ATOOLS_Phys_Particle_Qualifier_H
#define ATOOLS_Phys_Particle_Qualifier_H

#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Org/Getter_Function.H"

namespace ATOOLS {

  class Particle_Qualifier_Base {
  public:
    virtual ~Particle_Qualifier_Base();
    virtual bool operator() (const Particle *) const = 0;
    static void ShowQualifiers(const int mode);
    void Keep(Particle_List *const list);
    void Erase(Particle_List *const list);
  };

  typedef Getter_Function<Particle_Qualifier_Base,std::string> 
  Particle_Qualifier_Getter;

  
  class Or_Particle_Qualifier : public Particle_Qualifier_Base {
    Particle_Qualifier_Base * p_qual_a;
    Particle_Qualifier_Base * p_qual_b;
  public:
    Or_Particle_Qualifier(Particle_Qualifier_Base * a,
			   Particle_Qualifier_Base * b)
      : p_qual_a(a),p_qual_b(b) {}
    bool operator()(const Particle *) const;
  };

  class And_Particle_Qualifier : public Particle_Qualifier_Base {
    Particle_Qualifier_Base * p_qual_a;
    Particle_Qualifier_Base * p_qual_b;
  public:
    And_Particle_Qualifier(Particle_Qualifier_Base * a,
			   Particle_Qualifier_Base * b)
      : p_qual_a(a),p_qual_b(b) {}
    bool operator()(const Particle *) const;
  };

  class Not_Particle_Qualifier : public Particle_Qualifier_Base {
    Particle_Qualifier_Base * p_qual_a;
  public:
    Not_Particle_Qualifier(Particle_Qualifier_Base * qual)
      : p_qual_a(qual) {}
    bool operator()(const Particle *) const;
  };

  class Is_ME_Particle : public Particle_Qualifier_Base {
  public:
    bool operator()(const Particle *) const;
  };

  class Is_KF : public Particle_Qualifier_Base {
  private:
    kf_code m_kfcode;
  public:
    Is_KF(const std::string &kfcode);
    bool operator()(const Particle *) const;
  };

  class Is_Flav : public Particle_Qualifier_Base {
  private:
    Flavour m_flav;
  public:
    Is_Flav(const std::string &kfcode);
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

  class Is_Charged_Lepton : public Particle_Qualifier_Base {
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
