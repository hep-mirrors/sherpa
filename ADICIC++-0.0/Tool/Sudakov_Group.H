//bof
//Version: 4 ADICIC++-0.0/2006/02/03    (surely preliminarily)

//The Sudakov classes are provided, which manage the Sudakov calculation needed
//for the gluon emission or splitting of a specific FF dipole.
//                                                  =========



#ifndef _Sudakov_Group_H_
#define _Sudakov_Group_H_ _Sudakov_Group_H_


//#include <...>
#include "Sudakov_Calculator.H"





namespace ADICIC {



  template<Dipole::Type _Diptype>
  class Sudakov_Group : public Sudakov_Calculator {    //i.e. FFSudakov_Group

  public:

    //typedef, enum, class;

  private:

    //enum, friend;

    //Blocked methods:

    Sudakov_Group(const Sudakov_Group<_Diptype>&);
    Sudakov_Group<_Diptype>& operator=(const Sudakov_Group<_Diptype>&);
    const bool operator==(const Sudakov_Group<_Diptype>&) const;

  private:

    const Radiation::Type m_radtype;

    double m_s;
    double m_x2tmin;
    double m_x2tmax;
    double m_x2t;
    double m_ymax;
    double m_rap;
    double m_corr;

    std::list<Sudakov_Base*> l_sud;    //Put the Sudakov's together.

    //Private methods:

    inline const bool TestEfracs(const double, const double) const;
    inline       void Reset();

    const bool InitRadiation() const;
    const bool InitWithCurrentDipole();

  public:

    Sudakov_Group(const Radiation::Type=Radiation::g);
    ~Sudakov_Group();

    //Public (inline) methods:

    inline void Which() const;
    inline void ShowSpecification() const;

    inline const double Sdip() const;
    inline const double X2tmin() const;
    inline const double X2t() const;
    inline const double Ymax() const;

    const bool GenerateVariablesFor(const Dipole&, Sudakov_Result&);

  };



  //###########################################################################



  //Blocked classes:

  template<Dipole::Type _Diptype, Radiation::Group _Radgroup>
  class Sudakov : public Sudakov_Base {
    Sudakov();
    void ShowSpecification() const;
    const double CalculateRapLimit();
    const double GenerateRap();
    const double GenerateCorr(const Multidouble&);
    const double GenerateX2t();
  };
  template<>
  class Sudakov<Dipole::qqbar,Radiation::quark> : public Sudakov_Base {
    Sudakov();
    void ShowSpecification() const;
    const double CalculateRapLimit();
    const double GenerateRap();
    const double GenerateCorr(const Multidouble&);
    const double GenerateX2t();
  };





  //Specified class: gluon emission

  template<Dipole::Type _Diptype>
  class Sudakov<_Diptype,Radiation::gluon> : public Sudakov_Base {

    //i.e. FFSudakov

  private:

    Sudakov();
    Sudakov(const Sudakov<_Diptype,Radiation::gluon>&);
    Sudakov<_Diptype,Radiation::gluon>&
    operator=(const Sudakov<_Diptype,Radiation::gluon>&);
    const bool operator==(const Sudakov<_Diptype,Radiation::gluon>&) const;

  private:

    //Compile time information.
    static const short  s_x1pow;
    static const short  s_x3pow;
    static const double s_colfac;

    const Sudakov_Group<_Diptype>& m_sgroup;

  public:

    Sudakov(const Sudakov_Group<_Diptype>&, const Sudakov_Flavour&);
    ~Sudakov();

    inline void Which() const;
    inline void ShowSpecification() const;

    inline const double CalculateRapLimit();
    inline const double GenerateRap();
    inline const double GenerateCorr(const Multidouble&);
    const double GenerateX2t();

  };





  //Specified class: quark-antiquark splitting

  template<Dipole::Type _Diptype>
  class Sudakov<_Diptype,Radiation::quark> : public Sudakov_Base {

    //i.e. FFSudakov

  private:

    Sudakov();
    Sudakov(const Sudakov<_Diptype,Radiation::quark>&);
    Sudakov<_Diptype,Radiation::quark>&
    operator=(const Sudakov<_Diptype,Radiation::quark>&);
    const bool operator==(const Sudakov<_Diptype,Radiation::quark>&) const;

  private:

    //Compile time information.
    static const double s_colfac;

    const Sudakov_Group<_Diptype>& m_sgroup;

    inline void SetSplitBranch();

  public:

    Sudakov(const Sudakov_Group<_Diptype>&, const Sudakov_Flavour&);
    ~Sudakov();

    inline void Which() const;
    inline void ShowSpecification() const;

    inline bool NfVeto(const double) const;

    inline const double CalculateRapLimit();
    inline const double GenerateRap();
    inline const double GenerateCorr(const Multidouble&);
    const double GenerateX2t();

  };



}    //eo namespace ADICIC





#include "Sudakov_Group.inl.hh"


#endif    //eo _Sudakov_Group_H_



//eof
