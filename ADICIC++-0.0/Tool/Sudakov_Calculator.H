//bof
//Version: 3 ADICIC++-0.0/2005/09/30    (surely preliminarily)

//Base classes are provided, which manage the Sudakov calculations.



#ifndef _Sudakov_Calculator_H_
#define _Sudakov_Calculator_H_ _Sudakov_Calculator_H_


#include <string>
#include "Model_Base.H"
#include "Function_Base.H"
#include "Running_AlphaS.H"
#include "PDF_Base.H"
#include "Sudakov_Utils.hpp"
#include "Dipole.H"


#define SUDAKOV_CALCULATOR_OUTPUT SUDAKOV_CALCULATOR_OUTPUT
#undef  SUDAKOV_CALCULATOR_OUTPUT
#define TEMP_OUTPUT TEMP_OUTPUT////////////////////////////////////////////////
//#undef  TEMP_OUTPUT//////////////////////////////////////////////////////////





namespace ADICIC {



  class Sudakov_Calculator {    //i.e. Sudakov_Group_Base

  private:

    typedef std::vector<ATOOLS::Function_Base*> Funcs;
    typedef std::vector<PDF::PDF_Base*>         Pdfs;

    typedef const double (*Double_Double_Func)(const double);
    typedef const int    (*Int_Double_Func)(const double);
    typedef const double (*PDF_Corr_Func)(bool, const Multiflavour&,
					  const Multidouble&);

    class Toolbox {
      friend class Sudakov_Calculator;
      Funcs m_ras;
      Pdfs  m_pdf;
      Toolbox();
      ~Toolbox();
    };

  private:

    Sudakov_Calculator(const Sudakov_Calculator&);
    Sudakov_Calculator& operator=(const Sudakov_Calculator&);
    const bool operator==(const Sudakov_Calculator&) const;

  private:

    static const bool sf_ariadne;
    static const bool sf_pdf;    //Temporary entry.

    static int s_count;

    static double s_asapprox;    //AlphaS approximation.

    static Toolbox s_box;

    static Double_Double_Func  GetAlphaSCorr;
    static inline const double FixAlphaSCorr(const double);
    static inline const double RunAlphaSCorr(const double);

    static Int_Double_Func  GetNf;
    static inline const int FixNf(const double);
    static inline const int RunNf(const double);

    static PDF_Corr_Func       GetPDFCorr;
    static inline const double NoPDFCorr(bool, const Multiflavour&,
					 const Multidouble&);
    static        const double IsPDFCorr(bool, const Multiflavour&,
					 const Multidouble&);

  protected:

    const Dipole*   p_dip;
    Sudakov_Result* p_sur;

  public:

    inline  Sudakov_Calculator();
    virtual ~Sudakov_Calculator();

    static const bool& Ariadne;
    static const int&  InStore;

    static inline const bool   IsAlphaSRunning();
    static inline const double AlphaSApprox();
    static inline const double AlphaSCorr(const double);
    static inline const int    Nf(const double);
    static inline const bool   ArePDFsInitialized();
    static inline const double PlusPDFCorr(const Multiflavour&,
					   const Multidouble&);
    static inline const double MinusPDFCorr(const Multiflavour&,
					    const Multidouble&);

    static       void ShowEnvironment();
    static const Trio AdjustEnvironment(const std::string& path="default",
					MODEL::Model_Base* pmb=NULL);

    inline const Dipole& CurrentDipole() const;
    virtual void         Which() const;

    virtual       void ShowSpecification() const=0;
    virtual const bool GenerateVariablesFor(const Dipole&, Sudakov_Result&)=0;

  };



  //###########################################################################



  class Sudakov_Base {

  private:

    Sudakov_Base(const Sudakov_Base&);
    Sudakov_Base& operator=(const Sudakov_Base&);
    const bool operator==(const Sudakov_Base&) const;

  private:

    static int s_count;

  protected:

    static const double s_average;

    Sudakov_Flavour  m_code;
    xbool            m_split;
    Radiation::Group m_rg;
    short            m_step;
    double           m_mass;

  public:

    inline  Sudakov_Base(const Sudakov_Flavour&, const Radiation::Group);
    virtual ~Sudakov_Base();

    static const int& InStore;

    inline const Sudakov_Flavour& RadCode() const;
    inline const xbool            RadPart() const;
    inline const Radiation::Group RadGroup() const;
    inline const double           RadMass() const;
    inline void                   InitRadParticle();

    virtual void Which() const;
    virtual void ShowSpecification() const=0;

    virtual bool   NfVeto(const double) const { return false;}
    virtual void   SetPDFFlavs(bool) {}
    virtual void   SetGenX2tFac() {}
    virtual double GeneratePDFCorr(const Multidouble&) { return 1.0;}

    virtual const double CalculateRapLimit()=0;
    virtual const double GenerateRap()=0;
    virtual const double GenerateCorr(const Multidouble&)=0;
    virtual const double GenerateX2t()=0;

  };



}    //eo namespace ADICIC





#include "Sudakov_Calculator.inl.hh"


#endif    //eo _Sudakov_Calculator_H_



//eof
