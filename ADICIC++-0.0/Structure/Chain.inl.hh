//bof
//Version: 2 ADICIC++-0.0/2004/09/10

//Inline methods of Chain.H.





#include <cassert>
#include <cstdlib>
#include "Chain_Handler.H"





namespace ADICIC {



  //===========================================================================



  inline const Bool Chain::Status() const {
    return varset.f_active;
  }
  inline const bool Chain::IsClear() const {
    if(!varset.f_clear) return false;
    assert(varset.l_dip.empty() && !varset.p_quab && !varset.p_1glu &&
	   !varset.p_hdl);    //Checking.
    return true;
  }
  inline const bool Chain::IsEmpty() const {
    return (varset.l_dip.empty() && !varset.p_quab && !varset.p_1glu);
  }





  inline const double Chain::LastScale() const {
    return varset.m_k2tlast;
  }


  inline const double Chain::Mass() const {
    return varset.m_mass;
  }
  inline const double Chain::InvMass() const {
    return varset.m_invmass;
  }
  inline const ATOOLS::Vec4D& Chain::Momentum() const {
    return varset.m_momentum;
  }





  inline const bool Chain::IsHandled() const {
    return bool(varset.p_hdl);
  }
  inline const bool Chain::IsHandledBy(const Chain_Handler& CH) const {
    return (varset.p_hdl==&CH);
  }





  inline const bool Chain::IsLine() const {
    if(varset.l_dip.empty()) return false;
    return bool(varset.l_dip.front()->GetTopBranchPointer().operator->()
		==varset.p_quab	&&
		varset.l_dip.back()->GetBotBranchPointer().operator->()
		==varset.p_atib);
  }


  inline const bool Chain::IsRing() const {
    if(varset.l_dip.empty()) return false;
    return bool(varset.l_dip.front()->GetTopBranchPointer().operator->()
		==varset.p_1glu	&&
		varset.l_dip.back()->GetBotBranchPointer().operator->()
		==varset.p_1glu);
  }





  inline const Dipole& Chain::ChainRoot() const {
    assert(varset.p_root);
    return *varset.p_root;
  }



  //===========================================================================



  inline std::size_t Chain::ParticleNumber() const {
    return (varset.l_glub.size()+bool(varset.p_quab)+bool(varset.p_atib));
  }
  inline std::size_t Chain::MaxParticleNumber() const {
    return (varset.l_glub.max_size()+2);
  }
  inline std::size_t Chain::DipoleNumber() const {
    return varset.l_dip.size();
  }
  inline std::size_t Chain::MaxDipoleNumber() const {
    return varset.l_dip.max_size();
  }



  //===========================================================================



  inline Bool& Chain::SetStatus() {
    //if(varset.p_hdl) varset.p_hdl->"Reset()";????????????????????????????????
    //else varset.f_clear=false;???????????????????????????????????????????????
    varset.f_clear=false;
    return varset.f_active;
  }


  inline double& Chain::SetLastScale() {
    //if(varset.p_hdl) varset.p_hdl->"Reset()";????????????????????????????????
    //else varset.f_clear=false;???????????????????????????????????????????????
    varset.f_clear=false;
    return varset.m_k2tlast;
  }



  //===========================================================================



  inline const bool Chain::operator|(Chain_Handler& CH) {
    if(varset.p_hdl || CH.IsDocked()) return false;
    varset.p_hdl=&CH;
    if(varset.p_hdl->AttachChain(this)) { varset.f_clear=false; return true;}
    varset.p_hdl=NULL; return false;
  }


  inline void Chain::operator|(bool s) {
    if(!varset.p_hdl) return;
    Chain_Handler* phand=varset.p_hdl;
    if(varset.p_hdl->IsDockedAt(*this)) {
      varset.p_hdl=NULL;    //Could be a clear chain right now.
      phand->DetachChain(this);
      return;
    }
    std::cerr<<"\nBug: Wrong Chain-Chain_Handler connection emerged!\n";
    assert(varset.p_hdl->IsDockedAt(*this));
  }



  //===========================================================================



  inline void Chain::Print() const {
    if(!s_print) { s_print=true; std::cout<<(*this)<<std::endl; s_print=false;}
    else std::cout<<(*this)<<std::endl;
  }



  //===========================================================================



  inline const Dipole*& Chain::RootPointer() {
    varset.f_clear=false;
    return varset.p_root;
  }
  inline Dipole::Branch*& Chain::BranchPointer() {
    varset.f_clear=false;
    return varset.p_quab;
  }
  inline const Dipole::Branch* Chain::BranchPointer() const {
    return varset.p_quab;
  }
  inline Dipole::Antibranch*& Chain::AntibranchPointer() {
    varset.f_clear=false;
    return varset.p_atib;
  }
  inline const Dipole::Antibranch* Chain::AntibranchPointer() const {
    return varset.p_atib;
  }
  inline const Dipole::Glubranch*& Chain::FirstGlubranchPointer() {
    varset.f_clear=false;
    return varset.p_1glu;
  }
  inline const Dipole::Glubranch* Chain::FirstGlubranchPointer() const {
    return varset.p_1glu;
  }
  inline std::list<Dipole::Glubranch*>& Chain::GlubranchPointerList() {
    varset.f_clear=false;
    return varset.l_glub;
  }
  inline
  const std::list<Dipole::Glubranch*>& Chain::GlubranchPointerList() const {
    return varset.l_glub;
  }
  inline std::list<Dipole*>& Chain::DipolePointerList() {
    varset.f_clear=false;
    return varset.l_dip;
  }
  inline const std::list<Dipole*>& Chain::DipolePointerList() const {
    return varset.l_dip;
  }



  //===========================================================================



}    //eo namespace ADICIC





//eof
