//bof
//Version: 4 ADICIC++-0.0/2006/07/02

//Inline methods of Chain_Handler.H.





#include <cassert>
#include <cstdlib>
#include "Evolution_Strategy.hpp"
#include "Sudakov_Calculator.H"





namespace ADICIC {



  //===========================================================================



  inline const bool Chain_Handler::HasNewChain() const {
    return bool(p_cix);
  }


  inline const bool Chain_Handler::IsDocked() const {
    return bool(p_cha);
  }


  inline const bool Chain_Handler::IsDockedAt(const Chain& cha) const {
    return (p_cha==&cha);
  }


  inline const double Chain_Handler::CompScale() const {
    return m_k2tcomp;
  }



  //===========================================================================



  inline const bool Chain_Handler::AttachChain(Chain* pC) {
    if(p_cix) return false;
    if(p_cha==NULL && pC->IsHandledBy(*this)) {
      p_cha=pC;
      m_box.Mup.Clear();
      m_code=ATOOLS::Flavour();
      f_below=false;
      this->PresetCompScale();
      this->CleanUpDHs();
      return true;
    }
    return false;
  }


  inline const bool Chain_Handler::DetachChain(const Chain* pC) {
    if(p_cha==pC && pC->IsHandled()==false) {
      this->FreeChain();
      //Decoupling of new results should still be guaranteed!
      p_cha=NULL;
      return true;
    }
    return false;
  }



  //===========================================================================



  inline const bool Chain_Handler::ProductionStrategyInfo() {    //Static.
    std::cout<<"\nFor the purpose of confirmation: "
	     <<"Chain_Evolution_Strategy::Production has been chosen!"
	     <<std::endl;
    return false;
  }
  inline const bool Chain_Handler::EmissionStrategyInfo() {    //Static.
    std::cout<<"\nFor the purpose of confirmation: "
	     <<"Chain_Evolution_Strategy::Emission has been chosen!"
	     <<std::endl;
    return false;
  }



  //===========================================================================



  inline void Chain_Handler::PresetCompScale() {
    m_k2tcomp=ATOOLS::Min(dpa.sud.MinK2t(),
			  dpa.sud.MinIIK2t());
  }





  template<> const bool
  Chain_Handler::FindTheDipole<Chain_Evolution_Strategy::Production>();
  template<> const bool
  Chain_Handler::FindTheDipole<Chain_Evolution_Strategy::Emission>();
  //template<> const bool
  //Chain_Handler::FindTheDipole<Chain_Evolution_Strategy::Mass>();



  inline const bool Chain_Handler::FindDipole() {
    //return FindTheDipole<Chain_Evolution_Strategy::Ret>();
#ifdef __GNUC__
#if __GNUC__ >2
    return (this->*v_finddip[dpa.evo.ChainEvolutionStrategy()[cel::def]])();
#else
#error
    switch(dpa.evo.ChainEvolutionStrategy()[cel::def]) {
    case Chain_Evolution_Strategy::Production:
      return (this->FindTheDipole<Chain_Evolution_Strategy::Production>());
    case Chain_Evolution_Strategy::Emission:
      return (this->FindTheDipole<Chain_Evolution_Strategy::Emission>());
    case Chain_Evolution_Strategy::Mass:
      return (this->FindTheDipole<Chain_Evolution_Strategy::Mass>());
    default:
      return (this->FindTheDipole<Chain_Evolution_Strategy::Unknown>());
    }
#endif
#endif
  }



  //===========================================================================



  inline Chain_Handler::Carrier::Carrier()
    : ChaOrder(0), EmitFlav(),
      pCha(NULL), Mup() {}


  inline void Chain_Handler::Carrier::Reset() {
    ChaOrder=0; EmitFlav=ATOOLS::Flavour();
    pCha=NULL; Mup.Clear();
  }



  //===========================================================================



}    //eo namespace ADICIC





//eof
