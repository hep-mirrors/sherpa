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

  typedef Getter_Function<Particle_Qualifier_Base,std::string> 
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





  inline bool Is_Photon::operator() (const Particle * p) const {
    if ( p && p->Flav().IsPhoton() ) return 1;
    return 0;
  };

  inline bool Is_Gluon::operator() (const Particle * p)  const {
    if ( p && p->Flav().IsGluon() ) return 1;
    return 0;
  };

  inline bool Is_Parton::operator() (const Particle * p)  const {
    if ( p && ( p->Flav().IsGluon() || p->Flav().IsQuark() ) ) return 1;
    return 0;
  };

  inline bool Is_Charged::operator() (const Particle * p)  const{
      if ( p && (p->Flav().IntCharge() !=0) ) return 1;
      return 0;
  };

  inline bool Is_Hadron::operator() (const Particle * p)  const{
    if ( p && p->Flav().IsHadron() &&  !p->Flav().IsDiQuark()) return 1;
    return 0;
  };

  inline bool Is_Charged_Hadron::operator() (const Particle * p)  const{
    if ( p && p->Flav().IntCharge() !=0 &&  
	 p->Flav().IsHadron() &&  !p->Flav().IsDiQuark()) return 1;
    return 0;
  };

  inline bool Is_Neutral_Hadron::operator() (const Particle * p)  const{
    if ( p && p->Flav().IntCharge() ==0 &&  
	 p->Flav().IsHadron() &&  !p->Flav().IsDiQuark()) return 1;
    return 0;
  };

  inline bool Is_Final_State::operator() (const Particle * p)  const{
    if ( p && (p->Status() == 1) ) return 1;
    return 0;
  };



  inline bool Is_Neutral::operator() (const Particle * p)  const{
    if ( p && p->Flav().IntCharge()==0) return 1;
    return 0;
  };

  inline bool Is_Charged_Pion::operator() (const Particle * p) const {
    if ( p && p->Flav().Kfcode()==kf::pi_plus) return 1;
    return 0;
  };
  inline bool Is_Neutral_Pion::operator() (const Particle * p) const {
    if ( p && p->Flav().Kfcode()==kf::pi) return 1;
    return 0;
  };

  inline bool Is_Charged_Kaon::operator() (const Particle * p) const {
    if ( p && p->Flav().Kfcode()==kf::K_plus) return 1;
    return 0;
  };
  inline bool Is_Neutral_Kaon::operator() (const Particle * p) const {
    if ( p && p->Flav().Kfcode()==kf::K) return 1;
    return 0;
  };


  inline bool Is_Charged_KStar::operator() (const Particle * p) const {
    if ( p && p->Flav().Kfcode()==kf::K_star_892_plus) return 1;
    return 0;
  };
  inline bool Is_Neutral_KStar::operator() (const Particle * p) const {
    if ( p && p->Flav().Kfcode()==kf::K_star_892) return 1;
    return 0;
  };

  inline bool Is_Rho0::operator() (const Particle * p) const {
    if ( p && p->Flav().Kfcode()==kf::rho_770) return 1;
    return 0;
  };
  inline bool Is_Eta::operator() (const Particle * p) const {
    if ( p && p->Flav().Kfcode()==kf::eta) return 1;
    return 0;
  };
  inline bool Is_EtaPrime::operator() (const Particle * p) const {
    if ( p && p->Flav().Kfcode()==kf::eta_prime_958) return 1;
    return 0;
  };

  inline bool Is_Phi::operator() (const Particle * p) const {
    if ( p && p->Flav().Kfcode()==kf::phi_1020) return 1;
    return 0;
  };
  inline bool Is_Omega::operator() (const Particle * p) const {
    if ( p && p->Flav().Kfcode()==kf::omega_782) return 1;
    return 0;
  };
  inline bool Is_Lambda::operator() (const Particle * p) const {
    if ( p && p->Flav().Kfcode()==kf::Lambda) return 1;
    return 0;
  };
  inline bool Is_Charged_Sigma::operator() (const Particle * p) const {
    if ( p && p->Flav().Kfcode()==kf::Sigma_minus) return 1;
    return 0;
  };
  inline bool Is_Charged_Xi::operator() (const Particle * p) const {
    if ( p && p->Flav().Kfcode()==kf::Xi_minus) return 1;
    return 0;
  };
  inline bool Is_Neutral_Xi::operator() (const Particle * p) const {
    if ( p && p->Flav().Kfcode()==kf::Xi) return 1;
    return 0;
  };

  inline bool Is_Proton_Antiproton::operator() (const Particle * p) const {
    if ( p && p->Flav().Kfcode()==2212) return 1;
    return 0;
  };

  inline bool Is_Not_Lepton::operator() (const Particle * p) const {
    if ( p && !p->Flav().IsLepton() ) return 1;
    return 0;
  };

  inline bool Is_Not_Neutrino::operator() (const Particle * p) const {
    if ( p && !(p->Flav().IsLepton() && p->Flav().IntCharge()==0) ) return 1;
    return 0;
  };


  inline bool Is_There::operator() (const Particle * p) const {
    if ( p ) return 1;
    return 0;
  };

}// end of namespace ATOOLS

#endif
