//bof
//Version: 1 ADICIC++-0.0/2004/06/01    (surely preliminarily)

//The Dipole_Handler structure has the task of a Dipole modifier.



#include "Dipole.H"


#ifndef _Dipole_Handler_H_
#define _Dipole_Handler_H_ _Dipole_Handler_H_


#include <map>


#define DIPOLE_HANDLER_OUTPUT DIPOLE_HANDLER_OUTPUT
#undef  DIPOLE_HANDLER_OUTPUT





namespace ADICIC {



  class Sudakov_Calculator;
  class Recoil_Calculator;





  //std::ostream& operator<<(std::ostream&, const ...&);





  class Dipole_Handler {

  public:

    //typedef,enum,class

  private:

    //friend std::ostream& operator<<(std::ostream&, const ...&);

    //enum

    class Calcpair {
      friend class Dipole_Handler;
      Sudakov_Calculator* p_sud;
      Recoil_Calculator*  p_rec;
      Calcpair() : p_sud(NULL), p_rec(NULL) {}
      ~Calcpair();
    };

    typedef std::map<Dipole::Type,Calcpair*> Calcbox;

    //Blocked methods:
    //================
    Dipole_Handler(const Dipole_Handler&);
    Dipole_Handler& operator=(const Dipole_Handler&);
    const bool operator==(const Dipole_Handler&) const;    //What to compare?

  private:

    static int     s_count;
    static Calcbox s_map;

    static const bool sf_init;

    Calcpair* p_tempa;

    Sudakov_Calculator* p_sudakov;
    Recoil_Calculator*  p_recoil;

    Dipole* p_dip;

    Dipole*             p_dix;
    Dipole::Branch*     p_ban;
    Dipole::Antibranch* p_ati;
    Dipole::Glubranch*  p_glu;

    bool f_below;
    Trio f_recoil;    //Positive = top branch, Negative = bottom branch, Nil
    int  f_gate;

    double m_p2t;
    double m_x1;
    double m_x3;

    ATOOLS::Vec4D m_p1;
    ATOOLS::Vec4D m_p2;
    ATOOLS::Vec4D m_p3;

    //Private static methods:
    //=======================
    static const bool InitCalcBox();

    //Private inline methods:
    //=======================

    //Private methods:
    //================
    const bool GenerateMomenta();
    const bool GenerateSplitting();
    void CrossProductTest(const ATOOLS::Vec4D&) const;

  public:

    Dipole_Handler();
    Dipole_Handler(Dipole&);
    ~Dipole_Handler();

    //References: just return values
    //==============================
    static const int& InStore;

    //Public static methods:
    //======================
    static void ShowCalcBox();

    //Public inline methods:
    //======================
    inline const bool4::level Status() const;
    inline const int          IsDocked() const;
    inline const bool         IsDockedAt(const Dipole&) const;
    inline const bool         IsWaiting() const;

    //Only useful in interplay with the Dipole::operator|'s.
    inline const bool AttachDipole(Dipole*);
    inline const bool DetachDipole(const Dipole*);

    inline void ResetStatus();
    inline void DecoupleNewDipole(Dipole*&, bool&);
    inline void DecoupleGlubranch(Dipole::Glubranch*&);

    //Public methods:
    //===============
    void ShowSudakov() const;
    void ShowRecoilStrategy() const;
    const bool InduceGluonEmission();
    const bool FinishGluonEmission();
    const bool ManageGluonEmission();

  };



}    //eo namespace ADICIC





#include "Dipole_Handler.inl.hh"


#endif    //eo _Dipole_Handler_H_



//eof
