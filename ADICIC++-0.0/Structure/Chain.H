//bof
//Version: 4 ADICIC++-0.0/2006/06/11

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
      static const struct EpEm {} epem;    //Simple & ordinary e+e- collision.
    };

  private:

    friend class Chain_Handler;
    friend std::ostream& operator<<(std::ostream&, const Chain&);

    //typedef,enum

    struct Kernel {
      //An intra-struct!
    private:
      friend class Chain;

      Bool          f_active;
      double        m_k2tlast;
      double        m_mass;
      double        m_invmass;
      ATOOLS::Vec4D m_momentum;
      Multidouble   m_mufac;

      const Dipole* p_root;

      Chain_Handler* p_hdl;

      Dipole::Branch*               p_quab;
      Dipole::Antibranch*           p_atib;
      const Dipole::Glubranch*      p_1glu;
      std::list<Dipole::Glubranch*> l_glub;
      std::list<Dipole*>            l_dip;
      std::list<ATOOLS::Particle*>  l_corr;

      Kernel();
      ~Kernel() {}
      const bool IsInit() const;
      //Of course, the proper construction and destruction has to be
      //controlled by the chain itself.
      void Init();
      void Destruct();
    };

    static bool        s_print;
    static std::size_t s_count;
    static std::size_t s_maxcount;

    std::size_t m_name;
    std::size_t m_memo;

    //Responsibility for all particles related to chain is with the chain.
    Kernel                        varset;
    std::list<ATOOLS::Particle*>  l_fill;

    //Blocked methods:
    //================
    const bool operator==(const Chain&) const;

    //Private inline methods:
    //=======================
    inline       Multidouble&                   FactScaleBox();
    inline const Dipole*&                       RootPointer();
    inline       Dipole::Branch*&               BranchPointer();
    inline const Dipole::Branch*                BranchPointer()const;
    inline       Dipole::Antibranch*&           AntibranchPointer();
    inline const Dipole::Antibranch*            AntibranchPointer()const;
    inline const Dipole::Glubranch*&            FirstGlubranchPointer();
    inline const Dipole::Glubranch*             FirstGlubranchPointer()const;
    inline       std::list<Dipole::Glubranch*>& GlubranchPointerList();
    inline const std::list<Dipole::Glubranch*>& GlubranchPointerList()const;
    inline       std::list<Dipole*>&            DipolePointerList();
    inline const std::list<Dipole*>&            DipolePointerList()const;
    inline       std::list<ATOOLS::Particle*>&  CorrParticlePointerList();
    inline const std::list<ATOOLS::Particle*>&  CorrParticlePointerList()const;

    //Private methods:
    //================
    void                 Copy(const Chain&, const Type);
    const ATOOLS::Vec4D& UpdateMomentum(double, const ATOOLS::Vec4D&);

  public:

    Chain();
    Chain(const Chain&);
    Chain(const Dipole::Branch&, const Dipole::Antibranch&,
	  const Initiator::EpEm);
    Chain(const Dipole::Glubranch&, const Dipole::Glubranch&,
	  const Initiator::EpEm);
    ~Chain();

    Chain& operator=(const Chain&);

    //References: just return values
    //==============================
    static bool&              MoreCout;
    static const std::size_t& InStore;
    static const std::size_t& MaxCount;
    const std::size_t&        Name;
    const std::size_t&        Source;

    //Public inline methods:
    //======================
    inline const Bool           Status() const;
    inline const bool           IsClear() const;
    inline const bool           IsEmpty() const;
    inline const double         LastScale() const;
    inline const double         FactScale() const;
    inline const double         Mass() const;
    inline const double         InvMass() const;
    inline const ATOOLS::Vec4D& Momentum() const;
    inline const bool           IsHandled() const;
    inline const bool           IsHandledBy(const Chain_Handler&) const;
    inline const bool           IsLine() const;
    inline const bool           IsRing() const;
    inline const Dipole&        ChainRoot() const;

    inline std::size_t ParticleNumber() const;
    inline std::size_t MaxParticleNumber() const;
    inline std::size_t DipoleNumber() const;
    inline std::size_t MaxDipoleNumber() const;
    inline std::size_t CorrParticleNumber() const;

    inline std::size_t& SetSource();
    inline Bool&        SetStatus();
    inline double&      SetLastScale();

    inline const bool operator|(Chain_Handler&);
    inline void       operator|(bool);

    inline void Print() const;

    //Public methods:
    //===============
    const Type     ChainType() const;
    const unsigned INumber() const;
    const bool     CheckMomentumConservation(ATOOLS::Vec4D&) const;
    const boolint  ExtractPartons(ATOOLS::Particle_List&) const;

    const bool Clear();

    const bool AddCorrParticleToInit(const ATOOLS::Particle&);
    const bool Initialize(const Dipole::Branch&, const Dipole::Antibranch&,
			  double scale=false);
    //const bool Initialize(const Dipole::Glubranch&, const Dipole::Glubranch&,
    //                      double scale=false);

    const bool StartInit(const Dipole::Branch&);///////////////////////////////
    const bool StartInit(const Dipole::Glubranch&);////////////////////////////
    const bool AddToInit(const Dipole::Glubranch&, const double&);/////////////
    const bool FinishInit(const Dipole::Antibranch&, const double&);///////////
    const bool FinishInit(const Dipole::Glubranch&, const double&);////////////

  };



}    //eo namespace ADICIC





#include "Chain.inl.hh"


#endif    //eo _Chain_H_



//eof
