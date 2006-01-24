//bof
//Version: 3 ADICIC++-0.0/2005/09/22    (surely preliminarily)

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
	emission =10,
	gluon    =7,
	dquark   =1,
	uquark   =2,
	squark   =3,
	cquark   =4,
	bquark   =5,
	antid    =11,
	antiu    =12,
	antis    =13,
	antic    =14,
	antib    =15,
	splitting=20,
	dad      =21,
	uau      =22,
	sas      =23,
	cac      =24,
	bab      =25,
	stop     =26
	//Meaningful code. Do not change numbers!
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
    unsigned            m_nin;
    Cascade::Type       m_typ;
    bool                f_bot;
    ATOOLS::Flavour     m_fla;
    Recoil_Tool*        p_ret;
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
