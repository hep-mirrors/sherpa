//bof
//Version: 2 ADICIC++-0.0/2004/09/09    (surely preliminarily)

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

    typedef const bool (ADICIC::Chain_Handler::*FindDipole_Func)();

    //Blocked methods:
    //================
    Chain_Handler(const Chain_Handler&);
    Chain_Handler& operator=(const Chain_Handler&);
    const bool operator==(const Chain_Handler&) const;    //What to compare?

  private:

    static int s_count;
    static int s_param;

    bool             f_below;
    ATOOLS::kf::code m_code;
    Chain*           p_cix;

    double m_k2tcomp;
    Chain* p_cha;

    Dipole_Handler m_dh1;
    Dipole_Handler m_dh2;

    Dipole_Handler* p_dhwait;
    Dipole_Handler* p_dhaciv;
    Dipole_Handler* p_dhtemp;

    std::list<Dipole*>::iterator i_fix;
    std::list<Dipole*>::iterator i_run;

    FindDipole_Func fp_finddip[4];

    //========================
    bool                f_bot;
    Trio                m_rec;
    Chain::Type         m_typ;
    Dipole*             p_dw;
    Dipole::Glubranch*  p_gw;
    Dipole::Antibranch* p_aw;
    Dipole::Branch*     p_bw;
    Dipole*             p_win;
    ATOOLS::Vec4D       m_vec;
    ATOOLS::Vec4D       m_old;
    //========================

    //Private static methods:
    //=======================
    static inline const bool ProductionStrategyInfo();

    //Private inline methods:
    //=======================
    inline void       PresetCompScale();
    inline void       CleanUp();
    inline const bool FindDipole();

    //Private methods:
    //================
    void FreeChain();    //Needed since here we may have handlers acting.

    template<class _Strategy> const bool FindTheDipole();

    const bool ModifyChain();

    void EmitGluon();
    void TopSplitLine();
    void BotSplitLine();
    void TopSplitRing();
    void BotSplitRing();

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

    static void      ShowParameters();
    static const int AdjustParameters();

    //Public inline methods:
    //======================
    inline const bool   HasNewChain() const;
    inline const bool   IsDocked() const;
    inline const bool   IsDockedAt(const Chain&) const;
    inline const double CompScale() const;

    //Only useful in interplay with the Chain::operator|'s.
    inline const bool AttachChain(Chain*);
    inline const bool DetachChain(const Chain*);

    inline void DecoupleNew(Chain*&, ATOOLS::kf::code&, bool&);
    inline void RemoveNewProducts();

    //Public methods:
    //===============
    const bool EvolveChainByOneStep();
    const bool EvolveChain();//////////////////////////////////////////////////

  };



}    //eo namespace ADICIC





#include "Chain_Handler.inl.hh"


#endif    //eo _Chain_Handler_H_



//eof
