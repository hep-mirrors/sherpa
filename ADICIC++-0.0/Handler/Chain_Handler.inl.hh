//bof
//Version: 1 ADICIC++-0.0/2004/06/03

//Inline methods of Chain_Handler.H.





#include <cassert>
#include <cstdlib>
#include "Sudakov_Calculator.H"





namespace ADICIC {



  //===========================================================================



  inline const bool Chain_Handler::IsNewChain() const {
    return bool(p_cix);
  }


  inline const bool Chain_Handler::IsDocked() const {
    return bool(p_cha);
  }


  inline const bool Chain_Handler::IsDockedAt(const Chain& cha) const {
    if(p_cha==&cha) return true; return false;
  }


  inline const double Chain_Handler::CompScale() const {
    return m_k2tcomp;
  }



  //===========================================================================



  inline const bool Chain_Handler::AttachChain(Chain* pC) {
    if( p_cha==NULL && pC->IsHandledBy(*this) ) {
      p_cha=pC;
      return true;
    }
    return false;
  }


  inline const bool Chain_Handler::DetachChain(const Chain* pC) {
    if(p_cha==pC && pC->IsHandled()==false) {
      this->Reset();
      p_cha=NULL;
      return true;
    }
    return false;
  }



  //===========================================================================



  inline void Chain_Handler::DecoupleNewChain(Chain*& pC, bool& below) {
    if(pC || !p_cix) return;
    pC=p_cix; p_cix=NULL;
    below=f_below;
  }



  //===========================================================================



  inline void Chain_Handler::PresetCompScale() {
    m_k2tcomp=Sudakov_Calculator::MinOfK2t();
  }



  template<>
  const bool Chain_Handler::FindDipole<Chain_Evolution_Strategy::Production>();
  //template<>
  //const bool Chain_Handler::FindDipole<Chain_Evolution_Strategy::Emission>();
  //template<>
  //const bool Chain_Handler::FindDipole<Chain_Evolution_Strategy::Mass>();



  inline const bool Chain_Handler::FindDipole() {
    return FindDipole<Chain_Evolution_Strategy::Ret>();
  }



  //===========================================================================



}    //eo namespace ADICIC





//eof
