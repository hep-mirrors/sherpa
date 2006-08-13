//bof
//Version: 4 ADICIC++-0.0/2006/07/27

//Inline methods of Dipole_Parameter.H.





#include <cassert>
#include <cstdlib>
#include "MathTools.H"
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
  inline const xbool Dipole_Parameter::Sud::GsplitRule() const {
    return s_gsplit;
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



  inline const fascat::code Dipole_Parameter::Evo::FactScaleType() const {
    return s_fascatype;
  }
  inline const double
  Dipole_Parameter::Evo::GetFactScaleFrom(const Multidouble& mudo) const {
    if(mudo.size()>=fascat::stop) return mudo[s_fascatype]+s_scaoffset;
    else return s_scaoffset;
  }
  inline const double
  Dipole_Parameter::Evo::GetFactScaleFrom(const Sudakov_Result& sure) const {
    switch(s_fascatype) {
    case fascat::p2t : return sure.P2t+s_scaoffset;
    case fascat::k2t : if(sure.Isr.empty()) return s_scaoffset;
      return ATOOLS::sqr(sure.Isr[sr::kt])+s_scaoffset;
    case fascat::m2t : if(sure.Isr.empty()) return s_scaoffset;
      return ATOOLS::sqr(sure.Isr[sr::mt])+s_scaoffset;
    case fascat::shat: if(sure.Isr.empty()) return s_scaoffset;
      return sure.Isr[sr::shat]+s_scaoffset;
    default: assert(s_fascatype==fascat::p2t ||
		    s_fascatype==fascat::k2t ||
		    s_fascatype==fascat::m2t ||
		    s_fascatype==fascat::shat);
    }
    return 0.0;
  }



  inline const double Dipole_Parameter::Evo::FactScaleFactor() const {
    return s_fmuf;
  }
  inline const double Dipole_Parameter::Evo::RenoScaleFactor() const {
    return s_fmur;
  }



  inline const std::vector<Chain_Evolution_Strategy::Type>&
  Dipole_Parameter::Evo::ChainEvolutionStrategy() const {
    return v_chevostrat;
  }



  inline const std::size_t
  Dipole_Parameter::Evo::ChainParticleLimit() const {
    return s_chpartlim;
  }
  inline const std::size_t
  Dipole_Parameter::Evo::ChainCorrelationLimit() const {
    return s_chcorrlim;
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



  inline void Dipole_Parameter::Evo::SetFactScaleOffset(double offset) {
    //Once explicitly set to zero, no chance to turn it on again.
    static bool active=true;
    if(!active) return;
    if(offset<=0.0) active=false;
    s_scaoffset=offset;
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
