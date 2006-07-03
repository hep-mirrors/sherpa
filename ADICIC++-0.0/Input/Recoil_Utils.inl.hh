//bof
//Version: 4 ADICIC++-0.0/2006/05/21

//Inline methods of Recoil_Utils.hpp.





#include <cassert>





namespace ADICIC {



  //===========================================================================



  inline const rdt::code Multipoincare::Mode() const {
    return m_mode;
  }


  inline
  const std::vector<Multipoincare::Set>& Multipoincare::GetItsVec() const {
    return v_trafs;
  }


  inline void Multipoincare::Swap(Multipoincare& mup) {
    rdt::code temp=m_mode;
    m_mode=mup.m_mode;
    mup.m_mode=temp;
    v_trafs.swap(mup.v_trafs);
  }


  inline void Multipoincare::Keep(ATOOLS::Poincare& traf, bool rot, bool bak) {
    //The last one in the vector will be applied first.
    std::vector<Set>::iterator it=v_trafs.begin();
    v_trafs.insert(it,Set(0,NULL));
    Set& pair=v_trafs.front();
    pair.first=2*rot+bak;
    pair.second=&traf;
  }


  inline void Multipoincare::Destruct() {
    for(size_t i=0; i<v_trafs.size(); ++i) {
      if(v_trafs[i].second) delete v_trafs[i].second;
    }
    v_trafs.clear();
  }



  //===========================================================================



  inline void Recoil_Result::Reset(rdt::code c) {
    Poc=both;
    Vec.clear();
    Mup.Clear(c);
  }



  //===========================================================================



}    //eo namespace ADICIC





//eof
