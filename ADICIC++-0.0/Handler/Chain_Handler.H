//bof
//Version: 1 ADICIC++-0.0/2004/05/27    (surely preliminarily)

//Only the Chain_Handler structure can reasonably modify the state of a Chain.



#include "Chain.H"


#ifndef _Chain_Handler_H_
#define _Chain_Handler_H_ _Chain_Handler_H_


//#include <map>


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

    //enum,class,typedef

    //Blocked methods:
    //================
    Chain_Handler(const Chain_Handler&);
    Chain_Handler& operator=(const Chain_Handler&);
    const bool operator==(const Chain_Handler&) const;    //What to compare?

  private:

    static int s_count;

    bool f_below;
    Chain* p_cix;

    Chain* p_cha;

    Dipole_Handler* p_dhwait;
    Dipole_Handler* p_dhaciv;

    Dipole_Handler m_dh1;
    Dipole_Handler m_dh2;

    //Private static methods:
    //=======================

    //Private inline methods:
    //=======================

    //Private methods:
    //================

  public:

    Chain_Handler();
    Chain_Handler(Chain&);
    ~Chain_Handler();

    //References: just return values
    //==============================
    static const int& InStore;

    //Public static methods:
    //======================

    //Public inline methods:
    //======================
    inline const bool IsNewChain() const;
    inline const bool IsDocked() const;
    inline const bool IsDockedAt(const Chain&) const;

    //Only useful in interplay with the Chain::operator|'s.
    inline const bool AttachChain(Chain*);
    inline const bool DetachChain(const Chain*);

    inline void DecoupleNewChain(Chain*&, bool&);

    //Public methods:
    //===============

  };



}    //eo namespace ADICIC





#include "Chain_Handler.inl.hh"


#endif    //eo _Chain_Handler_H_



//eof
