//bof
//Version: 1 ADICIC++-0.0/2004/03/04

//The Dipole structure is an carrier of information.



#ifndef _Dipole_H_
#define _Dipole_H_ _Dipole_H_


#include <iostream>
#include <list>
#include "Particle.H"
#include "Vector.H"
#include "Information.hpp"


#define STRICT_VERSION STRICT_VERSION
#undef STRICT_VERSION





namespace ADICIC {



  class Dipole;

  class Dipole_Particle;

  class Dipole_Branch;
  class Dipole_Antibranch;
  class Dipole_Glubranch;





  std::ostream& operator<<(std::ostream&, const Dipole&);





  class Dipole {

  public:

    typedef Dipole_Branch     Branch;
    typedef Dipole_Antibranch Antibranch;
    typedef Dipole_Glubranch  Glubranch;

  private:

    enum Type {    //Blocked: qq=11, qbarq=-9, qbarqbar=-11, qbarg=-10, gq=1.
      qqbar    = 9,
      qg       = 10,
      gqbar    =-1,
      gg       = 0,
      incorrect=-9999
    };

  public:

    class Gate {
      friend class Dipole_Particle;
      Gate() {}
      ~Gate() {}
      inline bool operator()(const Dipole*, const Dipole_Particle*) const;
      void operator()(Dipole*, Dipole_Particle*, Dipole_Particle*, bool) const;
    };

    class Particle_Pointer;

  private:

    friend std::ostream& operator<<(std::ostream&, const Dipole&);

    friend bool Gate::operator()(const Dipole*, const Dipole_Particle*) const;
    friend void Gate::operator()(Dipole*, Dipole_Particle*, Dipole_Particle*,
				 bool) const;

  private:

    static int s_count;
    static int s_maxcount;

    bool f_top;
    bool f_bot;
    int  m_name;
    int  m_copy;

    int  m_memory;
    Bool f_active;

    Dipole_Particle* p_top;
    Dipole_Particle* p_bot;

    Type          m_type;
    double        m_mass;
    double        m_invmass;
    ATOOLS::Vec4D m_momentum;

    //Private inline methods:
    //=======================
    inline std::list<Dipole*>& AccessTopBranch() const;
    inline std::list<Dipole*>& AccessBotBranch() const;
    inline void AddDipoleToTowers();
    inline void RemoveDipoleFromTowers();
    inline void Test(const std::list<Dipole*>&, const Dipole_Particle*,
		     const std::string&) const;

  public:

    Dipole();
    Dipole(const Dipole&, bool phdl=false);
    Dipole(Branch&, Antibranch&, int source=0, bool phdl=false);
    Dipole(Branch&, Glubranch&, int source=0, bool phdl=false);
    Dipole(Glubranch&, Antibranch&, int source=0, bool phdl=false);
    Dipole(Glubranch&, Glubranch&, int source=0, bool phdl=false);
    ~Dipole();

    Dipole& operator=(const Dipole&);
    const bool operator==(const Dipole&) const;//-

    //References: just return values
    //==============================
    static const int& InStore;
    const int&        Name;
    const int&        CopyOf;

    //Public inline methods:
    //======================
    inline const short          PointerHandling() const;
    inline const int            Source() const;
    inline const Bool           Status() const;
    inline const Type           IsType() const;
    inline const double         Mass() const;
    inline const double         InvMass() const;
    inline const ATOOLS::Vec4D& TotP() const;

    inline void PrintTowers() const;

    inline int&  SetSource();
    inline Bool& SetStatus();

    inline const Type   UpdateType();
    inline const double UpdateMass();

    inline Particle_Pointer GetTopBranchPointer() const;
    inline Particle_Pointer GetBotBranchPointer() const;

    //Public methods:
    //===============
    void Print() const;
    void RenewBranch(Branch&);
    void RenewBranch(Antibranch&);
    void RenewBranch(bool, Glubranch&);
    const bool RenewBranch(bool, const Dipole&, bool);//-

  };





  class Dipole_Particle {

    friend std::list<Dipole*>& Dipole::AccessTopBranch() const;
    friend std::list<Dipole*>& Dipole::AccessBotBranch() const;

  private:

    static int s_count;
    static int s_maxcount;

    int m_num;

    Trio               m_typ;
    short              m_tag;
    ATOOLS::Particle   m_pac;
    std::list<Dipole*> m_tow;

    //Private inline methods:
    //=======================
    inline void Test() const;
    inline void SetPacNum();

    //Private methods:
    //================
    Dipole_Particle* Copy() const;

  protected:

    //Monitoring via output functions:
    //================================
    inline  void sc_info() const;
    inline  void cc_info() const;
    inline  void bc_info() const;
    inline  void sd_info() const;
    inline  void cp_info() const;
    virtual void nm_info() const;

  public:

    inline Dipole_Particle();    //Standard particle is the SM gluon!
    Dipole_Particle(const Dipole_Particle&);

    Dipole_Particle(const Dipole_Quark_Base&, const ATOOLS::Vec4D&);
    Dipole_Particle(const Dipole_Antiquark_Base&, const ATOOLS::Vec4D&);
    Dipole_Particle(const Dipole_Gluon_Base&, const ATOOLS::Vec4D&);

    virtual ~Dipole_Particle();

    Dipole_Particle& operator=(const Dipole_Particle&);
    const bool operator==(const Dipole_Particle&) const;//-

    //References:
    //===========
    static const int&       InStore;
    const int&              Name;
    const ATOOLS::Particle& Parton;

    //Public inline methods:
    //======================
    virtual void WhatIsIt() const;
    inline  void ShowParticle() const;

    inline const Trio              OrgType() const;
    inline const short             Tag() const;
    inline const ATOOLS::Particle& operator*() const;
    inline const ATOOLS::Flavour&  Flav() const;
    inline const ATOOLS::Vec4D&    Momentum() const;

    inline void SetMomentum(const ATOOLS::Vec4D&);    //Form is constrained :o(

    //Public methods:
    //===============
    void ShowDipoles() const;

    virtual const ATOOLS::Particle& Quarkize(const Dipole_Quark_Base&);
    virtual const ATOOLS::Particle& Antiquarkize(const Dipole_Antiquark_Base&);
    virtual const ATOOLS::Particle& Gluonize();
    //Later:
    //virtual const ATOOLS::Particle& Gluonize(const Dipole_Gluon_Base&);

  };





  struct Tag_Real {};    //Maybe these are not necessary.
  struct Tag_Anti {};
  struct Tag_Glue {};





  class Dipole_Branch : public Dipole_Particle, public Tag_Real {

  private:

    inline void cb_info() const;
    inline void db_info() const;
    inline void nm_info() const;

    Dipole_Branch(const Dipole_Antiquark_Base&, const ATOOLS::Vec4D&);
    Dipole_Branch(const Dipole_Gluon_Base&, const ATOOLS::Vec4D&);
    Dipole_Branch& operator=(const Dipole_Particle&);
    Dipole_Branch& operator=(const Dipole_Antibranch&);
    Dipole_Branch& operator=(const Dipole_Glubranch&);

    //Blocking of group-changing methods:
    const ATOOLS::Particle& Quarkize(const Dipole_Quark_Base&) {}
    const ATOOLS::Particle& Antiquarkize(const Dipole_Antiquark_Base&) {}
#ifdef STRICT_VERSION
    const ATOOLS::Particle& Gluonize() {}
#endif

  public:

    Dipole_Branch()
      : Dipole_Particle(info.quark.d,ATOOLS::Vec4D()) { cb_info();}
    Dipole_Branch(const Dipole_Branch& DB)
      : Dipole_Particle(DB) { cb_info();}

    Dipole_Branch(const ATOOLS::Particle&);//-
    Dipole_Branch(const Dipole_Quark_Base& Q, const ATOOLS::Vec4D& P)
      : Dipole_Particle(Q,P) { cb_info();}

    ~Dipole_Branch() { db_info();}

    Dipole_Branch& operator=(const Dipole_Branch& DB) {
      this->Dipole_Particle::operator=(DB); return *this;
    }
    void WhatIsIt() const { std::cout<<"Dipole_Branch."<<std::endl;}

  };





  class Dipole_Antibranch : public Dipole_Particle, public Tag_Anti {

  private:

    inline void ca_info() const;
    inline void da_info() const;
    inline void nm_info() const;

    Dipole_Antibranch(const Dipole_Quark_Base&, const ATOOLS::Vec4D&);
    Dipole_Antibranch(const Dipole_Gluon_Base&, const ATOOLS::Vec4D&);
    Dipole_Antibranch& operator=(const Dipole_Particle&);
    Dipole_Antibranch& operator=(const Dipole_Branch&);
    Dipole_Antibranch& operator=(const Dipole_Glubranch&);

    //Blocking of group-changing methods:
    const ATOOLS::Particle& Quarkize(const Dipole_Quark_Base&) {}
    const ATOOLS::Particle& Antiquarkize(const Dipole_Antiquark_Base&) {}
#ifdef STRICT_VERSION
    const ATOOLS::Particle& Gluonize() {}
#endif

  public:

    Dipole_Antibranch()
      : Dipole_Particle(info.antiquark.d,ATOOLS::Vec4D()) { ca_info();}
    Dipole_Antibranch(const Dipole_Antibranch& DA)
      : Dipole_Particle(DA) { ca_info();}

    Dipole_Antibranch(const ATOOLS::Particle&);//-
    Dipole_Antibranch(const Dipole_Antiquark_Base& A, const ATOOLS::Vec4D& P)
      : Dipole_Particle(A,P) { ca_info();}

    ~Dipole_Antibranch() { da_info();}

    Dipole_Antibranch& operator=(const Dipole_Antibranch& DA) {
      this->Dipole_Particle::operator=(DA); return *this;
    }
    void WhatIsIt() const { std::cout<<"Dipole_Antibranch."<<std::endl;}

  };





  class Dipole_Glubranch : public Dipole_Particle, public Tag_Glue {

  private:

    inline void cg_info() const;
    inline void dg_info() const;
    inline void nm_info() const;

    Dipole_Glubranch(const Dipole_Quark_Base&, const ATOOLS::Vec4D&);
    Dipole_Glubranch(const Dipole_Antiquark_Base&, const ATOOLS::Vec4D&);
    Dipole_Glubranch& operator=(const Dipole_Particle&);
    Dipole_Glubranch& operator=(const Dipole_Branch&);
    Dipole_Glubranch& operator=(const Dipole_Antibranch&);

    //Blocking of group-changing methods:
#ifdef STRICT_VERSION
    const ATOOLS::Particle& Quarkize(const Dipole_Quark_Base&) {}
    const ATOOLS::Particle& Antiquarkize(const Dipole_Antiquark_Base&) {}
#endif
    const ATOOLS::Particle& Gluonize() {}

  public:

    Dipole_Glubranch()
      : Dipole_Particle() { cg_info();}
    //Dipole_Glubranch()    //Later?
    //  : Dipole_Particle(info.gluon.g,ATOOLS::Vec4D()) { cg_info();}
    Dipole_Glubranch(const Dipole_Glubranch& DG)
      : Dipole_Particle(DG) { cg_info();}

    Dipole_Glubranch(const ATOOLS::Particle&);//-
    Dipole_Glubranch(const ATOOLS::Vec4D& P)
      : Dipole_Particle(info.gluon.g,P) { cg_info();}
    //Dipole_Glubranch(const Dipole_Gluon_Base& G, const ATOOLS::Vec4D& P)
    //  : Dipole_Particle(G,P) { cg_info();}    //Later.

    ~Dipole_Glubranch() { dg_info();}

    Dipole_Glubranch& operator=(const Dipole_Glubranch& DG) {
      this->Dipole_Particle::operator=(DG); return *this;
    }
    void WhatIsIt() const { std::cout<<"Dipole_Glubranch."<<std::endl;}

  };





  class Dipole::Particle_Pointer {

    friend class Dipole;

  private:

    Dipole_Particle* p_dipa;

    inline Particle_Pointer(const Particle_Pointer&);

    Particle_Pointer& operator=(const Particle_Pointer&);

  public:

    inline Particle_Pointer();
    inline ~Particle_Pointer();

    inline operator bool() const;
    inline const bool operator==(const Particle_Pointer&) const;
    inline Dipole_Particle* operator->() const;

  };



}    //eo namespace ADICIC





#include "Dipole.inl.hh"


#endif    //eo _Dipole_H_



//eof
