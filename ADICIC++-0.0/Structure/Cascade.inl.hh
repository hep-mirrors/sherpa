//bof
//Version: 2 ADICIC++-0.0/2004/09/09

//Inline methods of Cascade.H.





#include <cassert>
#include <cstdlib>
#include "Cascade_Handler.H"





namespace ADICIC {



  //===========================================================================



  inline const Bool Cascade::Status() const {
    return caset.f_active;
  }
  inline const bool Cascade::IsClear() const {
    return caset.IsInit();
  }
  inline const bool Cascade::IsEmpty() const {
    return (caset.l_ifo.empty() && caset.l_cha.empty() && caset.m_nroot==0);
  }





  inline const double Cascade::Mass() const {
    return caset.m_mass;
  }
  inline const double Cascade::InvMass() const {
    return caset.m_invmass;
  }
  inline const ATOOLS::Vec4D& Cascade::Momentum() const {
    return caset.m_momentum;
  }





  inline const bool Cascade::IsHandled() const {
    return bool(caset.p_hdl);
  }
  inline const bool Cascade::IsHandledBy(const Cascade_Handler& CH) const {
    return (caset.p_hdl==&CH);
  }





  inline const bool Cascade::IsLine() const {
    if(caset.l_cha.size()!=1) return false;
    return (caset.l_cha.front()->ChainType()==Chain::line);
  }


  inline const bool Cascade::IsRing() const {
    if(caset.l_cha.size()!=1) return false;
    return (caset.l_cha.front()->ChainType()==Chain::ring);
  }


  inline const bool Cascade::IsLines() const {
    if(caset.l_cha.size()<2) return false;
    for(std::list<Chain*>::const_iterator cat=caset.l_cha.begin();
	cat!=caset.l_cha.end(); ++cat)
      if((*cat)->ChainType()!=Chain::line) return false;
    return true;
  }





  inline const bool Cascade::IsEvolved() const {
    if(caset.l_cha.empty()) return false;
    for(std::list<Chain*>::const_iterator cat=caset.l_cha.begin();
	cat!=caset.l_cha.end(); ++cat)
      if((*cat)->Status()!=Off) return false;
    return true;
  }


  inline const bool Cascade::HasBlockedChain() const {
    if(caset.l_cha.empty()) return false;
    for(std::list<Chain*>::const_iterator cat=caset.l_cha.begin();
	cat!=caset.l_cha.end(); ++cat)
      if((*cat)->Status()==Blocked) return true;
    return false;
  }



  //===========================================================================



  inline std::size_t Cascade::RootChainNumber() const {
    return caset.m_nroot;
  }
  inline std::size_t Cascade::ChainNumber() const {
    return caset.l_cha.size();
  }
  inline std::size_t Cascade::MaxChainNumber() const {
    return caset.l_cha.max_size();
  }
  inline std::size_t Cascade::DipoleNumber() const {
    std::size_t num=0;
    for(std::list<Chain*>::const_iterator cat=caset.l_cha.begin();
	cat!=caset.l_cha.end(); ++cat)
      num+=(*cat)->DipoleNumber();
    return num;
  }
  inline std::size_t Cascade::ParticleNumber() const {
    std::size_t num=0;
    for(std::list<Chain*>::const_iterator cat=caset.l_cha.begin();
	cat!=caset.l_cha.end(); ++cat)
      num+=(*cat)->ParticleNumber();
    return num;
  }



  //===========================================================================



  inline Bool& Cascade::SetStatus() {
    //if(caset.p_hdl) caset.p_hdl->"Reset()";??????????????????????????????????
    return caset.f_active;
  }



  //===========================================================================



  inline const bool Cascade::operator|(Cascade_Handler& CH) {
    if(caset.p_hdl || CH.IsDocked()) return false;
    caset.p_hdl=&CH;
    if(caset.p_hdl->AttachCascade(this)) return true;
    caset.p_hdl=NULL; return false;
  }


  inline void Cascade::operator|(bool fake) {
    if(!caset.p_hdl) return;
    Cascade_Handler* phand=caset.p_hdl;
    if(caset.p_hdl->IsDockedAt(*this)) {
      caset.p_hdl=NULL;
      phand->DetachCascade(this);
      return;
    }
    std::cerr<<"\nBug: Wrong Cascade-Cascade_Handler connection emerged!\n";
    assert(caset.p_hdl->IsDockedAt(*this));
  }



  //===========================================================================



  inline void Cascade::Print() const {
    s_print=true;
    std::cout<<(*this)<<std::endl;
    s_print=false;
  }



  //===========================================================================



  inline const Chain* const Cascade::ChainPreparationPointer() const {
    return caset.p_add;
  }
  inline std::list<Cascade::Mirror>& Cascade::MirrorList() {
    return caset.l_ifo;
  }
  inline const std::list<Cascade::Mirror>& Cascade::MirrorList() const {
    return caset.l_ifo;
  }
  inline std::list<Chain*>& Cascade::ChainPointerList() {
    return caset.l_cha;
  }
  inline const std::list<Chain*>& Cascade::ChainPointerList() const {
    return caset.l_cha;
  }



  //===========================================================================



}    //eo namespace ADICIC





//eof
