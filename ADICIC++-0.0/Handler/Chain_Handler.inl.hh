//bof
//Version: 2 ADICIC++-0.0/2005/01/31

//Inline methods of Chain_Handler.H.





#include <cassert>
#include <cstdlib>
#include "Evolution_Strategy.hpp"
#include "Sudakov_Calculator.H"





namespace ADICIC {



  //===========================================================================



  inline const int Chain_Handler::ChainEvolutionStrategy() {    //Static.
    return s_param;
  }



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
      m_code=ATOOLS::kf::none;
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



  inline void Chain_Handler::DecoupleNew(Chain*& pC, ATOOLS::kf::code& code,
					 bool& below) {
    if(pC) return;
    pC=p_cix; p_cix=NULL;
    code=m_code; m_code=ATOOLS::kf::none;
    below=f_below; f_below=false;
  }


  inline void Chain_Handler::RemoveNewProducts() {
    //Resets the news.
    f_below=false;
    m_code=ATOOLS::kf::none;
    if(p_cix) { delete p_cix; p_cix=NULL;}
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
    m_k2tcomp=Sudakov_Calculator::MinOfK2t();
  }


  inline void Chain_Handler::CleanUp() {
    assert(!m_dh1.IsDocked() && !m_dh2.IsDocked());
    m_dh1.RemoveNewProducts();
    m_dh2.RemoveNewProducts();
    p_dhwait=&m_dh1;
    p_dhaciv=&m_dh2;
    p_dhtemp=NULL;
    i_fix=NULL;
    i_run=NULL;
  }





  template<> const bool
  Chain_Handler::FindTheDipole<Chain_Evolution_Strategy::Production>();
  //template<> const bool
  //Chain_Handler::FindTheDipole<Chain_Evolution_Strategy::Emission>();
  //template<> const bool
  //Chain_Handler::FindTheDipole<Chain_Evolution_Strategy::Mass>();



  inline const bool Chain_Handler::FindDipole() {
    //return FindDecDipole<Chain_Evolution_Strategy::Ret>();
#ifdef __GNUC__
#if __GNUC__ >2
    return (this->*fp_finddip[s_param])();
#else
    switch(s_param) {
    case 1:
      return (this->FindTheDipole<Chain_Evolution_Strategy::Production>());
    case 2:
      return (this->FindTheDipole<Chain_Evolution_Strategy::Emission>());
    case 3:
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
