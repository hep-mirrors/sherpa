//bof
//Version: 1 ADICIC++-0.0/2004/05/07    (surely preliminarily)

//Classes are provided, which manage the Sudakov calculation needed for the
//gluon emission or splitting of a specific dipole.



#ifndef _Sudakov_Calculator_H_
#define _Sudakov_Calculator_H_ _Sudakov_Calculator_H_


//#include <...>
#include "Dipole.H"
#include "Sudakov_Calculator.dat.hh"


#define SUDAKOV_CALCULATOR_OUTPUT SUDAKOV_CALCULATOR_OUTPUT
#undef  SUDAKOV_CALCULATOR_OUTPUT





namespace ADICIC {



  struct Alpha_S_Fix {};
  struct Alpha_S_Run {};




 
  class Sudakov_Calculator {

  private:

    Sudakov_Calculator(const Sudakov_Calculator&);
    Sudakov_Calculator& operator=(const Sudakov_Calculator&);
    const bool operator==(const Sudakov_Calculator&) const;

  private:

    static int s_count;

    static const bool sf_alphasrun;

    static const double s_k2tmin;
    static const double s_k2tmax;
    static const double s_alphasfix;

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
    static inline const double FixedAlphaS();

    inline void GetResult(bool&, double&, double&, double&) const;

    virtual void Which() const;
    virtual void ShowSpecification() const=0;
    virtual const bool GenerateEfracsFor(const Dipole&)=0;

  };





  template<Dipole::Type _DipType, class _AlphaSTreatment> class Sudakov
  : public Sudakov_Calculator {

  public:

    //typedef, enum, class;

  private:

    //enum, friend;

    //Blocked methods:

    Sudakov(const Sudakov<_DipType,_AlphaSTreatment>&);
    Sudakov<_DipType,_AlphaSTreatment>&
    operator=(const Sudakov<_DipType,_AlphaSTreatment>&);
    const bool operator==(const Sudakov<_DipType,_AlphaSTreatment>&) const;

  private:

    double m_alphas;

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

    inline Sudakov();
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
