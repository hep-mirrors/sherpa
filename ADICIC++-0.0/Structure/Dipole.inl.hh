//bof
//Version: 1 ADICIC++-0.0/2004/03/05

//Inline methods of Dipole.H.





#include <cassert>
#include <cstdlib>





namespace ADICIC {



  inline const short Dipole::PointerHandling() const {
    return (2*f_top+1*f_bot);
  }


  inline const int Dipole::Source() const {
    return m_memory;
  }


  inline const Bool Dipole::Status() const {
    return f_active;
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





  inline void Dipole::PrintTowers() const {
    p_top->ShowDipoles();
    std::cout<<"  >Dipole:"<<m_name<<"<";
    p_bot->ShowDipoles();
  }





  inline int& Dipole::SetSource() {
    return m_memory;
  }


  inline Bool& Dipole::SetStatus() {
    return f_active;
  }





  inline const Dipole::Type Dipole::UpdateType() {
    if(p_top->OrgType()==-1 || p_bot->OrgType()==1) {    //Cast protector!
      std::cerr<<"\nError: Demand will produce invalid Dipole type!\n";
      assert( (p_top->OrgType()==-1 || p_bot->OrgType()==1) == false );
    }
#ifdef STRICT_VERSION
    m_type=Type( (10*p_top->Tag()+p_bot->Tag()) );
    return m_type;
    //Does that sufficiently work? No.
    //Only if invalid Dipoles are really excluded.
#else
    short t = 10*p_top->Tag() + p_bot->Tag();
    if(t==9 || t==10 || t==-1 || t==0); else {
      std::cerr<<"\nError: Algorithm used produces invalid Dipole type!\n";
      assert(t==9 || t==10 || t==-1 || t==0);
    }
    m_type=Type(t);
    return m_type;
#endif
  }


  inline const double Dipole::UpdateMass() {
    m_momentum=p_top->Momentum()+p_bot->Momentum();
    m_invmass=m_momentum.Abs2();
    if(m_invmass<0.0) {
      std::cerr<<"\nMethod: const double ADICIC::Dipole::UpdateMass(): "
	       <<"Warning: Negative invariant mass!\n"<<std::endl;
      m_mass=-1*sqrt(-1*m_invmass);
    }
    else m_mass=sqrt(m_invmass);
    return m_mass;
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





  inline void Dipole::Test(const std::list<Dipole*>& L,
			   const Dipole_Particle* P,
			   const std::string& s) const {
    for(std::list<Dipole*>::const_iterator it=L.begin(); it!=L.end(); ++it) {
      if( (*it)->p_top==P ) {
	if( (*it)->f_top ) {
	  std::cerr<<"\nBug: "<<s
		   <<"Branch physically belongs to several Dipoles!\n";
	  assert((*it)->f_top==false);}
      } else {
	if( (*it)->p_bot==P ) {
	  if( (*it)->f_bot ) {
	    std::cerr<<"\nBug: "<<s
		     <<"Branch physically belongs to several Dipoles!\n";
	    assert((*it)->f_bot==false);}
	} else {
	  std::cerr<<"\nBug: Dipole_Particle's tower carries external Dipole!";
	  std::cerr<<std::endl; assert(0);
	}
      }
    }
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



  inline void Dipole_Particle::Test() const {
    Dipole::Gate gate;
    for(std::list<Dipole*>::const_iterator it=m_tow.begin();
	it!=m_tow.end(); ++it) {
      short ph=(*it)->PointerHandling();
      if(gate(*it,this)) {    //top branch
	if(ph==2 || ph==3) {
	  std::cerr<<"\nBug: TopBranch physically belongs to a Dipole!\n";
	  assert((ph==2 || ph==3)==false);
	}
      } else {    //bot branch
	if(ph==1 || ph==3) {
	  std::cerr<<"\nBug: BotBranch physically belongs to a Dipole!\n";
	  assert((ph==1 || ph==3)==false);
	}
      }
    }
  }





  inline void Dipole_Particle::SetPacNum() {
    m_pac.SetNumber(m_num);
  }



  //===========================================================================



  inline void Dipole_Particle::sc_info() const {
    std::cout<<"    construct standard Dipole_Particle ["<<m_num<<"]\n";
  }
  inline void Dipole_Particle::cc_info() const {
    std::cout<<"    construct copied Dipole_Particle ["<<m_num<<"]\n";
  }
  inline void Dipole_Particle::bc_info() const {
    std::cout<<"    construct Dipole_Particle ["<<m_num<<"]\n";
  }
  inline void Dipole_Particle::sd_info() const {
    std::cout<<"    destruct Dipole_Particle ["<<m_num<<"]\n";
  }
  inline void Dipole_Particle::cp_info() const {
    std::cout<<"    renew Dipole_Particle ["<<m_num<<"]\n";
  }
  void Dipole_Particle::nm_info() const {
    std::cout<<"[dipa"<<m_num<<"]";
  }



  //===========================================================================



  inline Dipole_Particle::Dipole_Particle()
    : m_num(++s_maxcount), Name(m_num),
      m_typ(Nil), m_tag(info.gluon.g.Tag()),
      m_pac( ATOOLS::Particle(m_num,info.gluon.g()) ), Parton(m_pac),
      m_tow(std::list<Dipole*>()) {

    ++s_count; this->sc_info();
  }





  void Dipole_Particle::WhatIsIt() const {
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
    std::cout<<"      Branch ["<<Name<<"]"<<std::endl;
  }
  inline void Dipole_Branch::db_info() const {
    std::cout<<"      ~Branch ["<<Name<<"]";
  }
  inline void Dipole_Branch::nm_info() const {
    std::cout<<"[dipa"<<Name<<".branch]";
  }





  inline void Dipole_Antibranch::ca_info() const {
    std::cout<<"      Antibranch ["<<Name<<"]"<<std::endl;
  }
  inline void Dipole_Antibranch::da_info() const {
    std::cout<<"      ~Antibranch ["<<Name<<"]";
  }
  inline void Dipole_Antibranch::nm_info() const {
    std::cout<<"[dipa"<<Name<<".antibranch]";
  }





  inline void Dipole_Glubranch::cg_info() const {
    std::cout<<"      Glubranch ["<<Name<<"]"<<std::endl;
  }
  inline void Dipole_Glubranch::dg_info() const {
    std::cout<<"      ~Glubranch ["<<Name<<"]";
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



}    //eo namespace ADICIC





//eof
