//bof
//Version: 2 ADICIC++-0.0/2004/08/10    (surely preliminarily)

//Only the Chain_Handler structure can reasonably modify the state of a Chain.



#include "Chain.H"


#ifndef _Chain_Handler_H_
#define _Chain_Handler_H_ _Chain_Handler_H_


//#include <...>
//#include "..."


#define CHAIN_HANDLER_OUTPUT CHAIN_HANDLER_OUTPUT
#undef  CHAIN_HANDLER_OUTPUT





namespace ADICIC {



  //class;





  //std::ostream& operator<<(std::ostream&, const ...&);





  class Chain_Handler {

  public:

    //typedef,enum,class

  private:

    //friend std::ostream& operator<<(std::ostream&, const ...&);

    //enum,class

    typedef const bool (ADICIC::Chain_Handler::*FindDecDip_Func)();

    //Blocked methods:
    //================
    Chain_Handler(const Chain_Handler&);
    Chain_Handler& operator=(const Chain_Handler&);
    const bool operator==(const Chain_Handler&) const;    //What to compare?

  private:

    static int s_count;
    static int s_param;

    bool f_below;
    Chain* p_cix;

    double m_k2tcomp;
    Chain* p_cha;

    Dipole_Handler m_dh1;
    Dipole_Handler m_dh2;

    Dipole_Handler* p_dhwait;
    Dipole_Handler* p_dhaciv;
    Dipole_Handler* p_dhtemp;

    std::list<Dipole*>::iterator i_fix;
    std::list<Dipole*>::iterator i_run;

    FindDecDip_Func fp_finddecdip[4];

    //Private static methods:
    //=======================
    static inline const bool ProductionStrategyInfo();

    //Private inline methods:
    //=======================
    inline void PresetCompScale();
    inline const bool FindDipole();

    //Private methods:
    //================
    template<class _Strategy> const bool FindDecDipole();
    const bool SplitDipole();

  public:

    Chain_Handler();
    Chain_Handler(Chain&);
    ~Chain_Handler();

    //References: just return values
    //==============================
    static const int& InStore;

    //Public static methods:
    //======================
    static inline const int ChainEvolutionStrategy();

    static       void ShowParameters();
    static const int  AdjustParameters();

    //Public inline methods:
    //======================
    inline const bool   IsNewChain() const;
    inline const bool   IsDocked() const;
    inline const bool   IsDockedAt(const Chain&) const;
    inline const double CompScale() const;

    //Only useful in interplay with the Chain::operator|'s.
    inline const bool AttachChain(Chain*);
    inline const bool DetachChain(const Chain*);

    inline void DecoupleNewChain(Chain*&, bool&);

    //Public methods:
    //===============
    void Reset();
    const bool EvolveChainByOneStep();
    const bool EvolveChain();

  };



}    //eo namespace ADICIC





#include "Chain_Handler.inl.hh"


#endif    //eo _Chain_Handler_H_



//eof
