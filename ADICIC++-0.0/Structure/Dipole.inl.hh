//bof
//Version: 1 ADICIC++-0.0/2004/07/05

//Inline methods of Dipole.H.





#include <cassert>
#include <cstdlib>
#include "Dipole_Handler.H"





namespace ADICIC {



  //===========================================================================



  inline const short Dipole::PointerHandling() const {
    return (2*f_top+1*f_bot);
  }


  inline const int Dipole::Source() const {
    return m_memory;
  }


  inline const Bool Dipole::Status() const {
    return f_active;
  }


  inline const double Dipole::ProdScale() const {
    return m_p2t;
  }
  inline const double Dipole::BootScale() const {
    return m_k2t;
  }
  inline const double Dipole::EmitScale() const {
    return m_l2t;
  }


  inline const Dipole::Type Dipole::IsType() const {
    return m_type;
  }


  inline const double Dipole::Mass() const {
    return m_mass;
  }
  inline const double Dipole::InvMass() const {
    return m_invmass;
  }
  inline const ATOOLS::Vec4D& Dipole::TotP() const {
    return m_momentum;
  }


  inline const bool Dipole::IsHandled() const {
    return bool(p_hdl);
  }


  inline const bool Dipole::IsHandledBy(const Dipole_Handler& DH) const {
    if(p_hdl==&DH) return true; return false;
  }





  inline void Dipole::PrintTowers() const {
    p_top->ShowDipoles();
    std::cout<<"  >Dipole:"<<m_name<<"<";
    p_bot->ShowDipoles();
  }





  inline int& Dipole::SetSource() {
    ++m_nchg;
    return m_memory;
  }


  inline Bool& Dipole::SetStatus() {
    ++m_nchg;
    return f_active;
  }


  inline double& Dipole::SetProdScale() {
    ++m_nchg;
    return m_p2t;
  }
  inline double& Dipole::SetBootScale() {
    ++m_nchg;
    return m_k2t;
  }
  inline double& Dipole::SetEmitScale() {
    ++m_nchg;
    return m_l2t;
  }


  inline const bool Dipole::operator|(Dipole_Handler& DH) {
    if(p_hdl || DH.IsDocked()) return false;
    p_hdl=&DH;
    if(p_hdl->AttachDipole(this)) return true;
    p_hdl=NULL; return false;
  }


  inline void Dipole::operator|(bool s) {
    if(!p_hdl) return;
    Dipole_Handler* phand=p_hdl;
    if(p_hdl->IsDockedAt(*this)) {
      p_hdl=NULL; phand->DetachDipole(this); return;
    }
    std::cerr<<"\nBug: Wrong Dipole-Dipole_Handler connection emerged!\n";
    assert(p_hdl->IsDockedAt(*this));
  }





  inline Dipole::Particle_Pointer Dipole::GetTopBranchPointer() const {
    Particle_Pointer pB;
    if(p_top->OrgType()==Positive) pB.p_dipa=static_cast<Branch*>(p_top);
    else pB.p_dipa=static_cast<Glubranch*>(p_top);
    return pB;
  }
  inline Dipole::Particle_Pointer Dipole::GetBotBranchPointer() const {
    Particle_Pointer pA;
    if(p_bot->OrgType()==Negative) pA.p_dipa=static_cast<Antibranch*>(p_bot);
    else pA.p_dipa=static_cast<Glubranch*>(p_bot);
    return pA;
  }



  //===========================================================================



  inline std::list<Dipole*>& Dipole::AccessTopBranch() const {
    return (p_top->m_tow);
  }
  inline std::list<Dipole*>& Dipole::AccessBotBranch() const {
    return (p_bot->m_tow);
  }
  inline void Dipole::AddDipoleToTowers() {
    std::list<Dipole*>& toptower=AccessTopBranch();
    std::list<Dipole*>& bottower=AccessBotBranch();
    toptower.remove(this);    //Not sure if this is really necessary.
    bottower.remove(this);
    toptower.push_front(this); bottower.push_front(this);
  }
  inline void Dipole::RemoveDipoleFromTowers() {
    std::list<Dipole*>& toptower=AccessTopBranch();
    std::list<Dipole*>& bottower=AccessBotBranch();
    toptower.remove(this); bottower.remove(this);
  }



  //===========================================================================



  inline bool Dipole::Gate::operator()(const Dipole* dip,
				       const Dipole_Particle* P) const {
    if(dip->p_top==P) return true;
    if(dip->p_bot==P) return false;
    std::cerr<<"\nBug: Dipole_Particle does not belong to this Dipole!\n";
    assert(0);
  }



  //===========================================================================



  inline void Dipole_Particle::SetPacNum() {
    m_pac.SetNumber(m_num);
  }



  //===========================================================================



  inline void Dipole_Particle::sc_info() const {
#ifdef DIPOLE_OUTPUT
    std::cout<<"    construct standard Dipole_Particle ["<<m_num<<"]\n";
#endif
  }
  inline void Dipole_Particle::cc_info() const {
#ifdef DIPOLE_OUTPUT
    std::cout<<"    construct copied Dipole_Particle ["<<m_num<<"]\n";
#endif
  }
  inline void Dipole_Particle::bc_info() const {
#ifdef DIPOLE_OUTPUT
    std::cout<<"    construct Dipole_Particle ["<<m_num<<"]\n";
#endif
  }
  inline void Dipole_Particle::sd_info() const {
#ifdef DIPOLE_OUTPUT
    std::cout<<"    destruct Dipole_Particle ["<<m_num<<"]\n";
#endif
  }
  inline void Dipole_Particle::cp_info() const {
#ifdef DIPOLE_OUTPUT
    std::cout<<"    renew Dipole_Particle ["<<m_num<<"]\n";
#endif
  }
  inline void Dipole_Particle::nm_info() const {
    std::cout<<"[dipa"<<m_num<<"]";
  }



  //===========================================================================



  inline Dipole_Particle::Dipole_Particle(Trio i, const ATOOLS::Particle& par)
    : m_num(++s_maxcount),
      m_typ(i), m_tag(i), m_pac(par),
      m_tow(std::list<Dipole*>()),
      Name(m_num), Parton(m_pac) {

    ++s_count; this->SetPacNum(); this->if_info();
  }


  inline short& Dipole_Particle::SetTag() {
    return m_tag;
  }


  inline void Dipole_Particle::if_info() const {
#ifdef DIPOLE_OUTPUT
    std::cout<<"    construct Dipole_Particle ["<<m_num<<"] from interface\n";
#endif
  }



  //===========================================================================



  inline void Dipole_Particle::WhatIsIt() const {
    std::cout<<"Dipole_Particle."<<std::endl;
  }


  inline void Dipole_Particle::ShowParticle() const {
    std::cout<<std::endl<<&m_pac<<std::endl;
  }





  inline const Trio Dipole_Particle::OrgType() const {
    return m_typ;
  }


  inline const short Dipole_Particle::Tag() const {
    return m_tag;
  }


  inline const ATOOLS::Particle& Dipole_Particle::operator*() const {
    return m_pac;
  }


  inline const ATOOLS::Flavour& Dipole_Particle::Flav() const {
    return m_pac.RefFlav();
  }


  inline const ATOOLS::Vec4D& Dipole_Particle::Momentum() const {
    return m_pac.Momentum();
  }





  inline void Dipole_Particle::SetMomentum(const ATOOLS::Vec4D& vec) {
    m_pac.SetMomentum(vec);
    //Towering!
    for(std::list<Dipole*>::iterator it=m_tow.begin(); it!=m_tow.end(); ++it)
      (*it)->UpdateMass();
    return;
  }



  //===========================================================================



  inline void Dipole_Branch::cb_info() const {
#ifdef DIPOLE_OUTPUT
    std::cout<<"      Branch ["<<Name<<"]"<<std::endl;
#endif
  }
  inline void Dipole_Branch::db_info() const {
#ifdef DIPOLE_OUTPUT
    std::cout<<"      ~Branch ["<<Name<<"]";
#endif
  }
  inline void Dipole_Branch::nm_info() const {
    std::cout<<"[dipa"<<Name<<".branch]";
  }





  inline void Dipole_Antibranch::ca_info() const {
#ifdef DIPOLE_OUTPUT
    std::cout<<"      Antibranch ["<<Name<<"]"<<std::endl;
#endif
  }
  inline void Dipole_Antibranch::da_info() const {
#ifdef DIPOLE_OUTPUT
    std::cout<<"      ~Antibranch ["<<Name<<"]";
#endif
  }
  inline void Dipole_Antibranch::nm_info() const {
    std::cout<<"[dipa"<<Name<<".antibranch]";
  }





  inline void Dipole_Glubranch::cg_info() const {
#ifdef DIPOLE_OUTPUT
    std::cout<<"      Glubranch ["<<Name<<"]"<<std::endl;
#endif
  }
  inline void Dipole_Glubranch::dg_info() const {
#ifdef DIPOLE_OUTPUT
    std::cout<<"      ~Glubranch ["<<Name<<"]";
#endif
  }
  inline void Dipole_Glubranch::nm_info() const {
    std::cout<<"[dipa"<<Name<<".glubranch]";
  }



  //===========================================================================



  inline Dipole::Particle_Pointer::Particle_Pointer() : p_dipa(NULL) {}


  inline Dipole::Particle_Pointer::~Particle_Pointer() { p_dipa=NULL;}


  inline Dipole::Particle_Pointer::operator bool() const {
    return bool(p_dipa);
  }


  inline const bool
  Dipole::Particle_Pointer::operator==(const Particle_Pointer& B) const {
    return (p_dipa==B.p_dipa);
  }


  inline Dipole_Particle* Dipole::Particle_Pointer::operator->() const {
    assert(p_dipa); return p_dipa;
  }


  inline Dipole::Particle_Pointer::Particle_Pointer(const Particle_Pointer& B)
    : p_dipa(B.p_dipa) {}    //Private!



  //===========================================================================



}    //eo namespace ADICIC





//eof
