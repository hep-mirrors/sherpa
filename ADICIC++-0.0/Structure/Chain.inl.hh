//bof
//Version: 1 ADICIC++-0.0/2004/05/21

//Inline methods of Chain.H.





#include <cassert>
#include <cstdlib>
//#include "Chain_Handler.H"





namespace ADICIC {



  //===========================================================================



  inline const Bool Chain::Status() const {
    return varset.f_active;
  }
  inline const bool Chain::IsEmpty() const {
    return varset.l_dip.empty();
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
  //inline const bool Chain::IsHandledBy(const Chain_Handler& CH) const {
  //  if(varset.p_hdl==&CH) return true; return false;
  //}





  inline const Chain::Type Chain::ChainType() const {
    if( varset.p_quab &&  varset.p_atib) return Chain::line;
    if(!varset.p_quab && !varset.p_atib &&
       !varset.l_dip.empty()) return Chain::ring;
    return Chain::incorrect;
  }
  inline const bool Chain::IsRing() const {
    return bool(!varset.p_quab && !varset.p_atib &&
		!varset.l_dip.empty());
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
    return varset.f_active;
  }


  inline double& Chain::SetLastScale() {
    return varset.m_k2tlast;
  }



  //===========================================================================



  //inline const bool Dipole::operator|(Dipole_Handler& DH) {
  //  if(varset.p_hdl || DH.IsDocked()) return false;
  //  varset.p_hdl=&DH;
  //  return varset.p_hdl->AttachDipole(this);
  //}


  //inline void Dipole::operator|(bool s) {
  //  if(!varset.p_hdl) return;
  //  Dipole_Handler* phand=varset.p_hdl;
  //  if(varset.p_hdl->IsDockedAt(*this)) {
  //    varset.p_hdl=NULL; phand->DetachDipole(this); return;
  //  }
  //  std::cerr<<"\nBug: Wrong Dipole-Dipole_Handler connection emerged!\n";
  //  assert(varset.p_hdl->IsDockedAt(*this));
  //}



  //===========================================================================



  inline Dipole::Branch*& Chain::ChainBranchPointer() {
    return varset.p_quab;
  }
  inline Dipole::Antibranch*& Chain::ChainAntibranchPointer() {
    return varset.p_atib;
  }
  inline std::list<Dipole::Glubranch*>& Chain::GlubranchPointerList() {
    return varset.l_glub;
  }
  inline std::list<Dipole*>& Chain::DipolePointerList() {
    return varset.l_dip;
  }
  inline const std::list<Dipole*>& Chain::DipolePointerList() const {
    return varset.l_dip;
  }



  //===========================================================================



}    //eo namespace ADICIC





//eof
