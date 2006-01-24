//bof
//Version: 3 ADICIC++-0.0/2005/09/21

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
      m_code=ATOOLS::Flavour();
      f_below=false;
      this->PresetCompScale();
      this->CleanUp();
      return true;
    }
    return false;
  }


  inline const bool Chain_Handler::DetachChain(const Chain* pC) {
    if(p_cha==pC && pC->IsHandled()==false) {
      this->FreeChain();
      p_cha=NULL;
      return true;
    }
    return false;
  }



  //===========================================================================



  inline void Chain_Handler::DecoupleNew(Recoil_Tool*& pR, Chain*& pC,
					 ATOOLS::Flavour& code, bool& below) {
    if(pR || pC) return;
    pR=p_rec; p_rec=NULL;
    pC=p_cix; p_cix=NULL;
    code=m_code; m_code=ATOOLS::Flavour();
    below=f_below; f_below=false;
  }


  inline void Chain_Handler::RemoveNewProducts() {    //Resets the news.
    f_below=false;
    m_code=ATOOLS::Flavour();
    if(p_cix) { delete p_cix; p_cix=NULL;}
    if(p_rec) { delete p_rec; p_rec=NULL;}
  }



  //===========================================================================



  inline const bool Chain_Handler::ProductionStrategyInfo() {    //Static.
    std::cout<<"\nFor the purpose of confirmation: "
	     <<"Chain_Evolution_Strategy::Production has been chosen!"
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
  //template<> const bool
  //Chain_Handler::FindTheDipole<Chain_Evolution_Strategy::Emission>();
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



}    //eo namespace ADICIC





//eof
