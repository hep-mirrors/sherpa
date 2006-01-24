//bof
//Version: 3 ADICIC++-0.0/2005/08/22    (surely preliminarily)

//The Dipole_Handler structure has the task of a Dipole modifier.



#include "Dipole.H"


#ifndef _Dipole_Handler_H_
#define _Dipole_Handler_H_ _Dipole_Handler_H_


#include <map>
#include <utility>
#include "Sudakov_Utils.hpp"
#include "Recoil_Strategy.hpp"


#define DIPOLE_HANDLER_OUTPUT DIPOLE_HANDLER_OUTPUT
#undef  DIPOLE_HANDLER_OUTPUT
//Dipole::PointerHandling() restriction:
#define STRICT_DIPOLE_HANDLER STRICT_DIPOLE_HANDLER
//#undef  STRICT_DIPOLE_HANDLER





namespace ADICIC {



  class Sudakov_Calculator;
  class Recoil_Calculator;





  //std::ostream& operator<<(std::ostream&, const ...&);





  class Dipole_Handler {

  public:

    typedef Dipole::Type         Dkey;
    typedef Radiation::Group     Rkey;
    typedef std::pair<Dkey,Rkey> Key;

  private:

    typedef std::map<Dkey,Sudakov_Calculator*> Sudakovbox;
    typedef std::map<Key,Recoil_Calculator*>   Recoilbox;

    class Calcbox {
      friend class Dipole_Handler;
      std::vector<Sudakov_Calculator*> v_psud;
      std::vector<Recoil_Calculator*>  v_prec;
      Calcbox();
      ~Calcbox();
    };

    //Blocked methods:
    //================
    Dipole_Handler(const Dipole_Handler&);
    Dipole_Handler& operator=(const Dipole_Handler&);
    const bool operator==(const Dipole_Handler&) const;    //What to compare?

  private:

    static bool sf_1stadjust;
    static int  s_count;

    static Sudakovbox s_sumap;
    static Recoilbox  s_remap;

    Key m_key;

    Sudakov_Calculator* p_sudakov;
    Recoil_Calculator*  p_recoil;

    Dipole* p_dip;

    Dipole*             p_dix;
    Dipole::Branch*     p_ban;
    Dipole::Antibranch* p_ati;
    Dipole::Glubranch*  p_glu;

    bool f_below;    //g emit: new dip's place, g split: place of splitting
    int  f_gate;

    Sudakov_Result m_sur;
    Recoil_Result  m_rer;

    //Private static methods:
    //=======================

    //Private inline methods:
    //=======================

    //Private methods:
    //================
    const bool GenerateSplitting();

  public:

    Dipole_Handler();
    Dipole_Handler(Dipole&);
    ~Dipole_Handler();

    //References: just return values
    //==============================
    static const int& InStore;

    //Public static methods:
    //======================
    static       void ListCalcBox();
    static       void ShowCalcBox();
    static const bool AdjustCalcBox();

    //Public inline methods:
    //======================
    inline const char Status() const;
    inline const int  IsDocked() const;
    inline const bool IsDockedAt(const Dipole&) const;
    inline const bool IsWaiting() const;

    //Only useful in interplay with the Dipole::operator|'s.
    inline const bool AttachDipole(Dipole*);
    inline const bool DetachDipole(const Dipole*);

    inline void DecoupleNew(Dipole*&, Dipole::Glubranch*&,
			    Dipole::Antibranch*&, Dipole::Branch*&,
			    bool&, xbool&);

    //Public methods:
    //===============
    void RemoveNewProducts();
    void ShowSudakov() const;
    void ShowRecoil() const;
    const bool InduceDipoleRadiation(bool t=false, bool o=false, bool f=false);
    const bool FinishDipoleRadiation();

  };



}    //eo namespace ADICIC





#include "Dipole_Handler.inl.hh"


#endif    //eo _Dipole_Handler_H_



//eof
