//bof
//Version: 3 ADICIC++-0.0/2005/09/08

//Inline methods of Dipole_Parameter.H.





#include <cassert>
#include <cstdlib>
#include "Run_Parameter.H"





namespace ADICIC {



  //===========================================================================



  inline const bool Dipole_Parameter::Sud::RunAlphaS() const {
    return s_runalphas;
  }
  inline const double Dipole_Parameter::Sud::AlphaSFix() const {
    return s_alphasfix;
  }
  inline const unsigned short Dipole_Parameter::Sud::NfFix() const {
    return s_nffix;
  }
  inline const Radiation::Type Dipole_Parameter::Sud::RadiationType() const {
    return s_radiatype;
  }
  inline const double Dipole_Parameter::Sud::MinK2t() const {
    return s_k2tmin;
  }
  inline const double Dipole_Parameter::Sud::MaxK2t() const {
    return s_k2tmax;
  }
  inline const double Dipole_Parameter::Sud::MinIIK2t() const {
    return s_k2tiimin;
  }
  inline const double Dipole_Parameter::Sud::MaxIIK2t() const {
    return s_k2tiimax;
  }
  inline const double Dipole_Parameter::Sud::IIEffExp() const {
    return s_iieffexp;
  }



  //---------------------------------------------------------------------------



  inline void Dipole_Parameter::Sud::SetMaxIIScale(const double d) {
    s_k2tiivarscale=d;
    if(s_k2tiivarscale>s_k2tiimin) s_k2tiimax=s_k2tiifac*s_k2tiivarscale;
    else s_k2tiimax=s_k2tiifac*s_k2tiifixscale;
    assert(s_k2tiimax>s_k2tiimin);
  }



  //===========================================================================



  inline const int Dipole_Parameter::Kin::ShowerMode() const {
    return s_dsmode;
  }
  inline const std::vector<Recoil_Strategy::Type>&
  Dipole_Parameter::Kin::RecoilStrategy() const {
    return v_recostrat;
  }



  //===========================================================================



  inline const std::vector<Chain_Evolution_Strategy::Type>&
  Dipole_Parameter::Evo::ChainEvolutionStrategy() const {
    return v_chevostrat;
  }



  //---------------------------------------------------------------------------



  inline void
  Dipole_Parameter::Evo::SetChainParticleLimit(std::size_t n) {
    assert(n>=2);
    s_chpartlim=n;
  }
  inline void
  Dipole_Parameter::Evo::SetChainCorrelationLimit(std::size_t n) {
    assert(n>=2);
    s_chcorrlim=n;
  }
  inline const std::size_t
  Dipole_Parameter::Evo::ChainParticleLimit() const {
    return s_chpartlim;
  }
  inline const std::size_t
  Dipole_Parameter::Evo::ChainCorrelationLimit() const {
    return s_chcorrlim;
  }



  //===========================================================================



  inline const double Dipole_Parameter::MaxIISHat(double s) const {
    assert(s>0.0);
    double val=sud.MaxIIK2t();
    return s+2*(val+sqrt(val*(s+val)));
  }
  inline const double Dipole_Parameter::MaxIIInvScale(double s) const {
    assert(s>0.0);
    double val=sud.MaxIIK2t()/s;
    val=2*(val+sqrt(val*(1.0+val)));
    return ATOOLS::sqr(val)*s/4.0;
  }
  inline const double Dipole_Parameter::HighestIIInvScale(double s) const {
    assert(s>0.0);
    double val=ATOOLS::sqr(ATOOLS::rpa.gen.Ecms())-s;
    assert(val>0.0);
    return ATOOLS::sqr(val)/4.0/s;
  }



  //===========================================================================



}    //eo namespace ADICIC





//eof
