//bof
//Version: 2 ADICIC++-0.0/2004/09/09

//Inline methods of Cascade_Handler.H.





#include <cassert>
#include <cstdlib>





namespace ADICIC {



  //===========================================================================



  inline const bool Cascade_Handler::IsEvolving() const {
    //List of chain iterators is not empty.
    return bool(!l_itt.empty());
  }
  inline const bool Cascade_Handler::HasEvolved() const {
    if(!p_cas) return false;
    return (p_cas->IsEvolved() && l_itt.empty());
  }





  inline const bool Cascade_Handler::IsDocked() const {
    return bool(p_cas);
  }
  inline const bool Cascade_Handler::IsDockedAt(const Cascade& cas) const {
    return (p_cas==&cas);
  }



  //===========================================================================



  inline const bool Cascade_Handler::AttachCascade(Cascade* pC) {
    //It is p_cas==NULL, if Cascade_Handler() or DetachCascade(...) have run.
    if(p_cas==NULL && pC->IsHandledBy(*this)) {
      p_cas=pC;
      //this->ResetCounter();//????????????????????????????????????????????????
      return true;
    }
    return false;
  }


  inline const bool Cascade_Handler::DetachCascade(const Cascade* pC) {
    if(p_cas==pC && pC->IsHandled()==false) {
      this->FreeCascade();
      p_cas=NULL;
      return true;
    }
    return false;
  }



  //===========================================================================



  inline
  std::size_t Cascade_Handler::ReadoutCounter(const Counter::code c) const {
    return v_count[c];
  }


  inline void Cascade_Handler::ResetCounter() {
    for(std::size_t c=Counter::start; c<Counter::stop; ++c) v_count[c]=0;
  }



  //===========================================================================



}    //eo namespace ADICIC





//eof
