//bof
//Version: 4 ADICIC++-0.0/2006/07/02    (surely preliminarily)

//Only the Chain_Handler structure can reasonably modify the state of a Chain.



#include <vector>
#include "Chain.H"
#include "Evolution_Strategy.hpp"


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

    //typedef,enum

    struct Carrier {
      bool            ChaOrder;
      ATOOLS::Flavour EmitFlav;
      Chain*          pCha;
      Multipoincare   Mup;
      inline Carrier();
      ~Carrier() {}
      inline void Reset();
    };

  private:

    //friend std::ostream& operator<<(std::ostream&, const ...&);

    //enum,class

#ifdef __GNUC__
#if __GNUC__ >2
    typedef const bool (ADICIC::Chain_Handler::*FindDipole_Func)();
#endif
#endif

    //Blocked methods:
    //================
    Chain_Handler(const Chain_Handler&);
    Chain_Handler& operator=(const Chain_Handler&);
    const bool operator==(const Chain_Handler&) const;    //What to compare?

  private:

    static int s_count;

    bool                    f_below;
    ATOOLS::Flavour         m_code;
    Chain*                  p_cix;

    Dipole_Handler::Carrier m_box;

    int    m_nii;
    double m_k2tcomp;
    Chain* p_cha;

    Dipole_Handler m_dh1;
    Dipole_Handler m_dh2;

    Dipole_Handler* p_dhwait;
    Dipole_Handler* p_dhaciv;
    Dipole_Handler* p_dhtemp;

    std::list<Dipole*>::iterator i_fix;
    std::list<Dipole*>::iterator i_run;

#ifdef __GNUC__
#if __GNUC__ >2
    std::vector<FindDipole_Func> v_finddip;
#endif
#endif

    //============================
    Chain::Type             m_typ;
    Dipole*                 p_win;
    double                  m_dm2;
    ATOOLS::Vec4D           m_emi;
    ATOOLS::Vec4D           m_vec;
    ATOOLS::Vec4D           m_old;
    //============================

    //Private static methods:
    //=======================
    static inline const bool ProductionStrategyInfo();
    static inline const bool EmissionStrategyInfo();

    //Private inline methods:
    //=======================
    inline void       PresetCompScale();
    inline const bool FindDipole();

    //Private methods:
    //================
    void InitStrategies();
    void CleanUpDHs();
    void FreeChain();    //Needed since here we may have handlers acting.

    template<Chain_Evolution_Strategy::Type _Strategy>
    const bool FindTheDipole();

    const bool ModifyChain();

    void DistributeIIRecoil();
    void EmitGluon();
    void EmitQuark();
    void TopSplitLine();    //Still as in version 2!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    void BotSplitLine();    //Still as in version 2!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    void TopSplitRing();    //Still as in version 2!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    void BotSplitRing();    //Still as in version 2!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  public:

    Chain_Handler();
    Chain_Handler(Chain&);
    ~Chain_Handler();

    //References: just return values
    //==============================
    static const int& InStore;

    //Public inline methods:
    //======================
    inline const bool   HasNewChain() const;
    inline const bool   IsDocked() const;
    inline const bool   IsDockedAt(const Chain&) const;
    inline const double CompScale() const;

    //Only useful in interplay with the Chain::operator|'s.
    inline const bool AttachChain(Chain*);
    inline const bool DetachChain(const Chain*);

    //Public methods:
    //===============
    void DecoupleNew(Carrier&);
    void RemoveNewProducts();
    const bool EvolveChainByOneStep();
    const bool EvolveChain();//////////////////////////////////////////////////
    const bool CorrectChain(const Multipoincare&);

  };



}    //eo namespace ADICIC





#include "Chain_Handler.inl.hh"


#endif    //eo _Chain_Handler_H_



//eof
