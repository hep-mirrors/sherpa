//bof
//Version: 1 ADICIC++-0.0/2004/03/12    (surely preliminarily)

//The Dipole_Handler structure has the task of a Dipole modifier.



#include "Dipole.H"


#ifndef _Dipole_Handler_H_
#define _Dipole_Handler_H_ _Dipole_Handler_H_


//#include <...>


#define DIPOLE_HANDLER_OUTPUT DIPOLE_HANDLER_OUTPUT
#undef  DIPOLE_HANDLER_OUTPUT





namespace ADICIC {



  //class;





  //std::ostream& operator<<(std::ostream&, const ...&);





  class Dipole_Handler {

  public:

    //typedef

    //class

  private:

    //enum

    //friend std::ostream& operator<<(std::ostream&, const ...&);

    //Blocked methods:
    //================
    Dipole_Handler(const Dipole_Handler&);
    Dipole_Handler& operator=(const Dipole_Handler&);
    const bool operator==(const Dipole_Handler&) const;    //What to compare?

  private:

    static int s_count;

    static const double s_k2tmin;
    static const double s_k2tmax;
    static const double s_alphasfix;

    Dipole* p_dip;

    Dipole*             p_dix;
    Dipole::Branch*     p_ban;
    Dipole::Antibranch* p_ati;
    Dipole::Glubranch*  p_glu;

    bool f_quarkrecoil;

    double m_p2t;
    double m_x1;
    double m_x3;

    ATOOLS::Vec4D m_p1;
    ATOOLS::Vec4D m_p2;
    ATOOLS::Vec4D m_p3;

    //Private inline methods:
    //=======================

    //Private methods:
    //================
    const bool GenerateEfracs();
    const bool TestEfracs() const;
    const bool GenerateMomenta();
    const bool GenerateSplitting();

  public:

    Dipole_Handler();
    Dipole_Handler(Dipole&);
    ~Dipole_Handler();

    //References: just return values
    //==============================
    static const int& InStore;

    //Public inline methods:
    //======================
    inline const bool4::level Status() const;
    inline const int          IsDocked() const;
    inline const bool         IsDockedAt(const Dipole&) const;

    inline const bool AttachDipole(Dipole*);
    inline const bool DetachDipole(const Dipole*);

    inline void DecoupleNewDipole(Dipole*&);
    inline void DecoupleGlubranch(Dipole::Glubranch*&);

    //Public methods:
    //===============
    const bool InduceGluonEmission();

  };



}    //eo namespace ADICIC





#include "Dipole_Handler.inl.hh"


#endif    //eo _Dipole_Handler_H_



//eof
