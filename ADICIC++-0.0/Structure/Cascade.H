//bof
//Version: 4 ADICIC++-0.0/2006/05/23

//This structure supports the cascading - the evolution of a list of Chains.



#ifndef _Cascade_H_
#define _Cascade_H_ _Cascade_H_


#include "Chain.H"


#define CASCADE_OUTPUT CASCADE_OUTPUT
#undef  CASCADE_OUTPUT





namespace ADICIC {



  class Cascade;

  class Cascade_Handler;





  std::ostream& operator<<(std::ostream&, const Cascade&);





  class Cascade {

  public:

    //typedef

    enum Type {
      ring     = 1,
      line     = 0,
      //rings    = 11,
      lines    = 10,
      incorrect=-9999
    };

    //struct Initiator { static const struct EpEm {} epem;};

  private:

    friend class Cascade_Handler;
    friend std::ostream& operator<<(std::ostream&, const Cascade&);

    //1st: name of the initiating root chain.
    //2nd: number of splittings of this chain.
    typedef std::pair<std::size_t,std::size_t> Mirror;

    //enum,class

    static bool s_print;
    static int  s_count;

    //In-struct:
    //==========
    struct {
    private:
      friend class Cascade;

      Bool          f_active;
      std::size_t   m_nroot;
      double        m_mass;
      double        m_invmass;
      ATOOLS::Vec4D m_momentum;

      Cascade_Handler* p_hdl;

      Chain* p_add;

      std::list<Mirror> l_ifo;    //Root chains are indicated by 1st=chainname.
      std::list<Chain*> l_cha;

      const bool IsInit() const {
	return (f_active==Blocked && m_nroot==0 &&
		m_momentum==ATOOLS::Vec4D() && p_hdl==NULL && p_add==NULL &&
		l_ifo.empty() && l_cha.empty());
      }
      void Init() {
	f_active=Blocked; m_nroot=0;
	m_mass=m_invmass=0.0; m_momentum=ATOOLS::Vec4D();
	p_hdl=NULL; p_add=NULL;	l_ifo.clear(); l_cha.clear();
      }
    } caset;

    //Blocked methods:
    //================
    const bool operator==(const Cascade&) const;

    //Private inline methods:
    //=======================
    inline const Chain* const       ChainPreparationPointer() const;
    inline std::list<Mirror>&       MirrorList();
    inline const std::list<Mirror>& MirrorList() const;
    inline std::list<Chain*>&       ChainPointerList();
    inline const std::list<Chain*>& ChainPointerList() const;

    //Private methods:
    //================
    void                 Destruct();
    void                 Copy(const Cascade&);
    const ATOOLS::Vec4D& UpdateMomentum(double, const ATOOLS::Vec4D&);

  public:

    Cascade();
    Cascade(const Cascade&);
    ~Cascade();

    Cascade& operator=(const Cascade&);

    //References: just return values
    //==============================
    static const int& InStore;

    //Public inline methods:
    //======================
    inline const Bool           Status() const;
    inline const bool           IsClear() const;
    inline const bool           IsEmpty() const;
    inline const double         Mass() const;
    inline const double         InvMass() const;
    inline const ATOOLS::Vec4D& Momentum() const;
    inline const bool           IsHandled() const;
    inline const bool           IsHandledBy(const Cascade_Handler&) const;
    inline const bool           IsLine() const;
    inline const bool           IsRing() const;
    inline const bool           IsLines() const;
    inline const bool           IsEvolved() const;
    inline const bool           HasBlockedChain() const;

    inline std::size_t RootChainNumber() const;
    inline std::size_t ChainNumber() const;
    inline std::size_t MaxChainNumber() const;
    inline std::size_t DipoleNumber() const;
    inline std::size_t ParticleNumber() const;
    inline std::size_t RelatedParticleNumber() const;

    inline Bool& SetStatus();

    inline const bool operator|(Cascade_Handler&);
    inline void       operator|(bool);

    inline void Print() const;

    //Public methods:
    //===============
    const Type     CascadeType() const;
    const unsigned INumber() const;
    const bool     CheckMomentumConservation(ATOOLS::Vec4D&) const;
    const boolint  ExtractPartons(std::list<ATOOLS::Particle_List>&) const;

    const bool Clear();

    const bool AddChain(const Dipole::Branch&, const Dipole::Antibranch&,
			const ATOOLS::Particle_List* pcorrs=NULL,
			double scale=false, bool ontop=false);
    //const bool AddChain(const Dipole::Glubranch&, const Dipole::Glubranch&,
    //                    const ATOOLS::Particle_List*,
    //                    double scale, bool ontop);
    //const bool StartAddingChain(const Dipole::Branch&);
    //const bool ContinueAddingChain(const Dipole::Glubranch&);
    //const bool StopAddingChain(const Dipole::Antibranch&);

  };



}    //eo namespace ADICIC





#include "Cascade.inl.hh"


#endif    //eo _Cascade_H_



//eof
