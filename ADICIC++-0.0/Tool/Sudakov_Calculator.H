//bof
//Version: 2 ADICIC++-0.0/2004/10/28    (surely preliminarily)

//Classes are provided, which manage the Sudakov calculation needed for the
//gluon emission or splitting of a specific dipole.



#ifndef _Sudakov_Calculator_H_
#define _Sudakov_Calculator_H_ _Sudakov_Calculator_H_


//#include <...>
#include "Model_Base.H"
#include "Function_Base.H"
#include "Running_AlphaS.H"
#include "Sudakov_Strategy.hpp"
#include "Sudakov_Result.hpp"
#include "Dipole.H"


#define SUDAKOV_CALCULATOR_OUTPUT SUDAKOV_CALCULATOR_OUTPUT
#undef  SUDAKOV_CALCULATOR_OUTPUT
#define TEMP_OUTPUT TEMP_OUTPUT////////////////////////////////////////////////
#undef  TEMP_OUTPUT////////////////////////////////////////////////////////////





namespace ADICIC {



  //Possibility of having information at compile time.
  template<Dipole::Type _Diptype, Radiation::Group _Radgroup>
  struct Sudakov_Info {};
  template<Dipole::Type _Diptype>
  struct Sudakov_Info<_Diptype,Radiation::gluon> {
    static const Dipole::Type     Dipoletype;
    static const Radiation::Group Radiationgroup;
    static const short            X1power;
    static const short            X3power;
    static const double           Colourfactor;
  };
  template<Dipole::Type _Diptype>
  struct Sudakov_Info<_Diptype,Radiation::quark> {
    static const Dipole::Type     Dipoletype;
    static const Radiation::Group Radiationgroup;
    static const double           Colourfactor;
  };



  //###########################################################################



  class Sudakov_Calculator {    //i.e. Sudakov_Group_Base

  private:

    typedef const double (*Double_Double_Func)(const double);
    typedef const int    (*Int_Double_Func)(const double);

  private:

    Sudakov_Calculator(const Sudakov_Calculator&);
    Sudakov_Calculator& operator=(const Sudakov_Calculator&);
    const bool operator==(const Sudakov_Calculator&) const;

  private:

    static const bool sf_start;
    static const bool sf_ariadne;

    static int s_count;

    static bool   s_isalphasrun;
    static double s_alphasfix;
    static double s_k2tmin;
    static double s_k2tmax;

    static double s_approx;

    static const int s_nffix;

    static ATOOLS::Function_Base* s_pas;

    static Double_Double_Func  GetAlphaSCorr;
    static inline const double FixAlphaSCorr(const double);
    static inline const double RunAlphaSCorr(const double);

    static Int_Double_Func  GetNf;
    static inline const int FixNf(const double);
    static inline const int RunNf(const double);

  protected:

    xbool m_spl;

    ATOOLS::kf::code m_rad;

    double m_p2t;
    double m_x1;
    double m_x3;

  public:

    inline  Sudakov_Calculator();
    virtual ~Sudakov_Calculator();

    static const bool& Ariadne;
    static const int&  InStore;

    static inline const bool   IsAlphaSRunning();
    static inline const int    NfFix();
    static inline const double MinOfK2t();
    static inline const double MaxOfK2t();
    static inline const double AlphaSFix();
    static inline const double AlphaSApprox();
    static inline const double AlphaSCorr(const double);
    static inline const int    Nf(const double);

    static       void ShowParameters();
    static const bool AdjustParameters();
    static const bool Init(MODEL::Model_Base*);

    inline void GetResult(Sudakov_Result&) const;

    virtual void Which() const;

    virtual       void ShowSpecification() const=0;
    virtual const bool GenerateEfracsFor(const Dipole&,
					 Sudakov_Strategy::Factorization)=0;
    virtual const bool GenerateEfracsFor(const Dipole&,
					 Sudakov_Strategy::Distribution)=0;

  };





  class Sudakov_Base {

  private:

    Sudakov_Base(const Sudakov_Base&);
    Sudakov_Base& operator=(const Sudakov_Base&);
    const bool operator==(const Sudakov_Base&) const;

  private:

    static int s_count;

  protected:

    const ATOOLS::kf::code m_code;

    xbool  m_split;
    double m_mass;

  public:

    inline  Sudakov_Base(ATOOLS::kf::code);
    virtual ~Sudakov_Base();

    static const int& InStore;

    inline const ATOOLS::kf::code RadCode() const;
    inline const xbool            RadPart() const;
    inline const double           RadMass() const;

    inline void InitRadParticle();

    virtual void Which() const;

    virtual       void   ShowSpecification() const=0;
    virtual const double CalculateRapLimit()=0;
    virtual const double GenerateRap()=0;
    virtual const double GenerateCorr(const double, const double)=0;
    virtual const double GenerateX2t()=0;

  };



  //###########################################################################



  template<Dipole::Type _Diptype>
  class Sudakov_Group : public Sudakov_Calculator {

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

    const ATOOLS::kf::code InitRadiation() const;
    void                   InitWith(const Dipole&);

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

    const bool GenerateEfracsFor(const Dipole&,
				 Sudakov_Strategy::Factorization);
    const bool GenerateEfracsFor(const Dipole&,
				 Sudakov_Strategy::Distribution);//////////////

  };



  //###########################################################################



  //Blocked classes:

  template<Dipole::Type _Diptype, Radiation::Group _Radgroup>
  class Sudakov : public Sudakov_Base {
    Sudakov();
    void ShowSpecification() const;
    const double CalculateRapLimit();
    const double GenerateRap();
    const double GenerateCorr(const double, const double);
    const double GenerateX2t();
  };
  template<>
  class Sudakov<Dipole::qqbar,Radiation::quark> : public Sudakov_Base {
    Sudakov();
    void ShowSpecification() const;
    const double CalculateRapLimit();
    const double GenerateRap();
    const double GenerateCorr(const double, const double);
    const double GenerateX2t();
  };





  //Specified class: gluon emission

  template<Dipole::Type _Diptype>
  class Sudakov<_Diptype,Radiation::gluon> : public Sudakov_Base {

  private:

    Sudakov();
    Sudakov(const Sudakov<_Diptype,Radiation::gluon>&);
    Sudakov<_Diptype,Radiation::gluon>&
    operator=(const Sudakov<_Diptype,Radiation::gluon>&);
    const bool operator==(const Sudakov<_Diptype,Radiation::gluon>&) const;

  private:

    short m_step;

    const Sudakov_Group<_Diptype>& m_sgroup;

  public:

    Sudakov(const Sudakov_Group<_Diptype>&, ATOOLS::kf::code);
    ~Sudakov() {
#ifdef TEMP_OUTPUT
      std::cout<<"~Sudakov(gluon)"<<std::endl;/////////////////////////////////
#endif
    }

    inline void Which() const;
    inline void ShowSpecification() const;

    inline const double CalculateRapLimit();
    inline const double GenerateRap();
    inline const double GenerateCorr(const double, const double);
    const double GenerateX2t();

  };





  //Specified class: quark-antiquark splitting

  template<Dipole::Type _Diptype>
  class Sudakov<_Diptype,Radiation::quark> : public Sudakov_Base {

  private:

    Sudakov();
    Sudakov(const Sudakov<_Diptype,Radiation::quark>&);
    Sudakov<_Diptype,Radiation::quark>&
    operator=(const Sudakov<_Diptype,Radiation::quark>&);
    const bool operator==(const Sudakov<_Diptype,Radiation::quark>&) const;

  private:

    short m_step;

    const Sudakov_Group<_Diptype>& m_sgroup;

    inline void SetSplitBranch();

  public:

    Sudakov(const Sudakov_Group<_Diptype>&, ATOOLS::kf::code);
    ~Sudakov() {
#ifdef TEMP_OUTPUT
      std::cout<<"~Sudakov(quark)"<<std::endl;/////////////////////////////////
#endif
    }

    inline void Which() const;
    inline void ShowSpecification() const;

    inline const double CalculateRapLimit();
    inline const double GenerateRap();
    inline const double GenerateCorr(const double, const double);
    const double GenerateX2t();

  };



}    //eo namespace ADICIC





#include "Sudakov_Calculator.inl.hh"


#endif    //eo _Sudakov_Calculator_H_



//eof
