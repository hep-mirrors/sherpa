//bof
//Version: 3 ADICIC++-0.0/2005/09/21

//Inline methods of Recoil_Strategy.hpp.





#include <cassert>





namespace ADICIC {



  //===========================================================================



  inline void Recoil_Result::Reset() {
    Poc=both;
    Vec.clear();
  }



  //===========================================================================



  inline const rdt::code Recoil_Tool::Mode() const {
    return m_mode;
  }


  inline bool Recoil_Tool::KeepThisBoost(ATOOLS::Poincare& B) {
    if(p_fly) {
      std::cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": "
	       <<"Warning: Already kept a boost!\n"<<std::endl;
      return false;
    }
    p_fly=&B;
    return true;
  }


  inline bool Recoil_Tool::KeepThisBackBoost(ATOOLS::Poincare& BB) {
    if(p_flyprime) {
      std::cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": "
	       <<"Warning: Already kept a backboost!\n"<<std::endl;
      return false;
    }
    p_flyprime=&BB;
    return true;
  }


  inline ATOOLS::Poincare& Recoil_Tool::GetBoost() const {
    assert(p_fly); return *p_fly;
  }


  inline ATOOLS::Poincare& Recoil_Tool::GetBackBoost() const {
    assert(p_flyprime); return *p_flyprime;
  }



  //===========================================================================



}    //eo namespace ADICIC





//eof
