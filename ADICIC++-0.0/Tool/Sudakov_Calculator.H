//bof
//Version: 2 ADICIC++-0.0/2004/08/09    (surely preliminarily)

//Classes are provided, which manage the Sudakov calculation needed for the
//gluon emission or splitting of a specific dipole.



#ifndef _Sudakov_Calculator_H_
#define _Sudakov_Calculator_H_ _Sudakov_Calculator_H_


//#include <...>
#include "Model_Base.H"
#include "Function_Base.H"
#include "Dipole.H"


#define SUDAKOV_CALCULATOR_OUTPUT SUDAKOV_CALCULATOR_OUTPUT
#undef  SUDAKOV_CALCULATOR_OUTPUT





namespace ADICIC {



  //Possibility of having information at compile time.
  //Some globally defined parameter sets.
  template<Dipole::Type _DipType> struct Sudakov_Info {
    static const Dipole::Type  Dipoletype;
    static const short         X1power;
    static const short         X3power;
    static const double        Colourfactor;
  };





  class Sudakov_Calculator {

  private:

    typedef const double (*AlphaSCorr_Func)(const double);

  private:

    Sudakov_Calculator(const Sudakov_Calculator&);
    Sudakov_Calculator& operator=(const Sudakov_Calculator&);
    const bool operator==(const Sudakov_Calculator&) const;

  private:

    static const bool sf_start;

    static int s_count;

    static bool   s_isalphasrun;
    static double s_alphasfix;
    static double s_k2tmin;
    static double s_k2tmax;

    static double s_approx;

    static ATOOLS::Function_Base* s_pas;

    static AlphaSCorr_Func     GetAlphaSCorr;
    static inline const double FixAlphaSCorr(const double);
    static inline const double RunAlphaSCorr(const double);

  protected:

    bool f_gsplit;

    double m_p2t;
    double m_x1;
    double m_x3;

  public:

    inline Sudakov_Calculator();
    virtual ~Sudakov_Calculator();

    static const int& InStore;

    static inline const bool   IsAlphaSRunning();
    static inline const double MinOfK2t();
    static inline const double MaxOfK2t();
    static inline const double AlphaSFix();
    static inline const double AlphaSApprox();
    static inline const double AlphaSCorr(const double);

    static       void ShowParameters();
    static const bool AdjustParameters();
    static const bool Init(MODEL::Model_Base*);

    inline void GetResult(bool&, double&, double&, double&) const;

    virtual void Which() const;
    virtual void ShowSpecification() const=0;
    virtual const bool GenerateEfracsFor(const Dipole&)=0;

  };





  template<Dipole::Type _DipType> class Sudakov : public Sudakov_Calculator {

  public:

    //typedef, enum, class;

  private:

    //enum, friend;

    //Blocked methods:

    Sudakov(const Sudakov<_DipType>&);
    Sudakov<_DipType>& operator=(const Sudakov<_DipType>&);
    const bool operator==(const Sudakov<_DipType>&) const;

  private:

    double m_s;
    double m_x2tmin;
    double m_x2t;
    double m_rap;
    double m_corr;

    //Private methods:

    inline void InitWith(const Dipole&);
    inline void GenerateRap();
    inline void GenerateCorr();
    inline const bool TestEfracs() const;

    const bool GenerateX2t();

  public:

    Sudakov();
    ~Sudakov() {}

    //Public (inline) methods:

    inline void Which() const;
    inline void ShowSpecification() const;

    const bool GenerateEfracsFor(const Dipole&);

  };



}    //eo namespace ADICIC





#include "Sudakov_Calculator.inl.hh"


#endif    //eo _Sudakov_Calculator_H_



//eof
