//bof
//Version: 2 ADICIC++-0.0/2004/09/09    (surely preliminarily)

//Only the Cascade_Handler structure can reasonably modify
//the state of a Cascade.



#include "Cascade.H"


#ifndef _Cascade_Handler_H_
#define _Cascade_Handler_H_ _Cascade_Handler_H_


//#include <...>
//#include "..."


#define CASCADE_HANDLER_OUTPUT CASCADE_HANDLER_OUTPUT
#undef  CASCADE_HANDLER_OUTPUT





namespace ADICIC {



  //class;





  class Cascade_Handler {

  public:

    //typedef,enum,class

    struct Counter {
      enum code {
	start    =0,
	total    =0,
	emission =1,
	splitting=2,
	gluon    =3,
	dquark   =4,
	uquark   =5,
	squark   =6,
	cquark   =7,
	bquark   =8,
	stop     =9
      };
    };

  private:

    //friend,typedef,enum,class

    //Blocked methods:
    //================
    Cascade_Handler(const Cascade_Handler&);
    Cascade_Handler& operator=(const Cascade_Handler&);
    const bool operator==(const Cascade_Handler&) const;    //What to compare?

  private:

    static int s_count;
 
    Cascade* p_cas;

    Chain_Handler m_chh;

    //Each Cascade_Handler counts on its own.
    std::vector<std::size_t> v_count;

    std::list<std::list<Cascade::Mirror>::iterator> l_mit;
    std::list<std::list<Chain*>::iterator>          l_itt;

    //========================
    bool                f_bot;
    ATOOLS::kf::code    m_kfc;
    Cascade::Type       m_typ;
    Chain*              p_nec;    //New and current chain.
    Chain*              p_cuc;    //Independent chains -> no winner.
    ATOOLS::Vec4D       m_vec;
    ATOOLS::Vec4D       m_old;
    //========================

    //Private static methods:
    //=======================

    //Private inline methods:
    //=======================

    //Private methods:
    //================
    void       FreeCascade();    //Needed since we may have handlers acting.
    const bool EvolveCurrentChain();

  public:

    Cascade_Handler();
    Cascade_Handler(Cascade&);
    ~Cascade_Handler();

    //References: just return values
    //==============================
    static const int& InStore;

    //Public static methods:
    //======================

    //Public inline methods:
    //======================
    inline const bool IsEvolving() const;
    inline const bool HasEvolved() const;
    inline const bool IsDocked() const;
    inline const bool IsDockedAt(const Cascade&) const;

    //Only useful in interplay with the Cascade::operator|'s.
    inline const bool AttachCascade(Cascade*);
    inline const bool DetachCascade(const Cascade*);

    inline std::size_t ReadoutCounter(const Counter::code) const;
    inline void        ResetCounter();

    //Public methods:
    //===============
    void       PrintCounter() const;
    const bool EvolveCascadeByOneStep();///////////////////////////////////////
    const bool EvolveCascade();

  };



}    //eo namespace ADICIC





#include "Cascade_Handler.inl.hh"


#endif    //eo _Cascade_Handler_H_



//eof
