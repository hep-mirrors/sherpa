//bof
//Version: 1 ADICIC++-0.0/2004/06/06

//This structure supports the evolution of a Chain of Dipoles.



#ifndef _Chain_H_
#define _Chain_H_ _Chain_H_


#include "Particle_List.H"
#include "Dipole.H"


#define CHAIN_OUTPUT CHAIN_OUTPUT
#undef  CHAIN_OUTPUT





namespace ADICIC {



  class Chain;

  class Chain_Handler;





  std::ostream& operator<<(std::ostream&, const Chain&);





  class Chain {

  public:

    //typedef

    enum Type {
      ring     = 1,
      line     = 0,
      incorrect=-9999
    };

    struct Initiator {
      static const struct Simple_Epem {} simple_epem;
    };

  private:

    friend class Chain_Handler;
    friend std::ostream& operator<<(std::ostream&, const Chain&);

    //typedef,enum,class

    static int s_count;

    //In-struct:
    //==========
    struct {
    private:
      friend class Chain;

      Bool          f_active;
      double        m_k2tlast;
      double        m_mass;
      double        m_invmass;
      ATOOLS::Vec4D m_momentum;

      const Dipole* p_root;

      Chain_Handler* p_hdl;

      Dipole::Branch*               p_quab;
      Dipole::Antibranch*           p_atib;
      const Dipole::Glubranch*      p_1glu;
      std::list<Dipole::Glubranch*> l_glub;
      std::list<Dipole*>            l_dip;

      void Init() {
	f_active=Blocked;
	m_k2tlast=m_mass=m_invmass=0.0;	m_momentum=ATOOLS::Vec4D();
	p_root=NULL;
	p_hdl=NULL;
	p_quab=NULL; p_atib=NULL; p_1glu=NULL;
	l_glub.clear();	l_dip.clear();
      }

    } varset;

    //Blocked methods:
    //================
    const bool operator==(const Chain&) const;

    //Private inline methods:
    //=======================
    inline Dipole::Branch*&               ChainBranchPointer();
    inline Dipole::Antibranch*&           ChainAntibranchPointer();
    inline const Dipole::Glubranch*       ChainFirstGlubranchPointer() const;
    inline std::list<Dipole::Glubranch*>& GlubranchPointerList();
    inline std::list<Dipole*>&            DipolePointerList();
    inline const std::list<Dipole*>&      DipolePointerList() const;

    //Private methods:
    //================
    const ATOOLS::Vec4D& UpdateMomentum(double, const ATOOLS::Vec4D&);

  public:

    Chain();
    Chain(const Chain&);
    Chain(const Dipole::Branch&, const Dipole::Antibranch&,
	  const Initiator::Simple_Epem);
    Chain(const Dipole::Glubranch&, const Dipole::Glubranch&,
	  const Initiator::Simple_Epem);
    ~Chain();

    Chain& operator=(const Chain&);

    //References: just return values
    //==============================
    static const int& InStore;

    //Public inline methods:
    //======================
    inline const Bool           Status() const;
    inline const bool           IsEmpty() const;
    inline const double         LastScale() const;
    inline const double         Mass() const;
    inline const double         InvMass() const;
    inline const ATOOLS::Vec4D& Momentum() const;
    inline const bool           IsHandled() const;
    inline const bool           IsHandledBy(const Chain_Handler&) const;
    inline const Type           ChainType() const;
    inline const bool           IsLine() const;
    inline const bool           IsRing() const;
    inline const Dipole&        ChainRoot() const;

    inline std::size_t ParticleNumber() const;
    inline std::size_t MaxParticleNumber() const;
    inline std::size_t DipoleNumber() const;
    inline std::size_t MaxDipoleNumber() const;

    inline Bool&      SetStatus();
    inline double&    SetLastScale();

    inline const bool operator|(Chain_Handler&);
    inline void       operator|(bool);

    //Public methods:
    //===============
    void Print() const;
    const bool CheckMomentumConservation(ATOOLS::Vec4D&) const;
    const bool ExtractPartons(ATOOLS::Particle_List&) const;

    const bool Clear();

    const bool Initialize(const Dipole::Branch&, const Dipole::Antibranch&);

    const bool StartInitChain(const Dipole::Branch&);//////////////////////////
    const bool StartInitChain(const Dipole::Glubranch&);///////////////////////
    const bool AddToInitChain(const Dipole::Glubranch&, const double&);////////
    const bool FinishInitChain(const Dipole::Antibranch&, const double&);//////
    const bool FinishInitChain(const Dipole::Glubranch&, const double&);///////

  };



}    //eo namespace ADICIC





#include "Chain.inl.hh"


#endif    //eo _Chain_H_



//eof
