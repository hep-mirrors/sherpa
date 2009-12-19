#ifndef PDF_Main_PDF_Base_H
#define PDF_Main_PDF_Base_H

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Getter_Function.H"

namespace ATOOLS {

  class Integration_Info;
  class Data_Reader;

}

namespace PDF {

  struct PDF_Arguments {
    ATOOLS::Flavour m_bunch;
    std::string m_path;
    ATOOLS::Data_Reader *p_read;
    PDF_Arguments(const ATOOLS::Flavour &bunch,
		  const std::string &path,
		  ATOOLS::Data_Reader *const read):
      m_bunch(bunch), m_path(path), p_read(read) {}
  };// end of struct PDF_Arguments

  class PDF_Base {
  public:

    typedef ATOOLS::Getter_Function<PDF_Base,PDF_Arguments> 
    PDF_Getter_Function;

  protected:

    std::vector<PDF_Base*> m_copies;

    ATOOLS::Flavour_Set m_partons;    
    ATOOLS::Flavour     m_bunch;

    std::string m_type;    

    int m_orderas, m_member;

    double m_exponent, m_rescale, m_asmz;
    double m_xmin, m_xmax, m_q2min, m_q2max;

    std::vector<std::pair<ATOOLS::Flavour,double> > m_extracted;

    void SingleExtract(const ATOOLS::Flavour flavour,const double x);
    void SingleReset();
    virtual void CalculateSpec(double x,double Q2) = 0;

  private:

    double m_fac_scale_factor, m_rq2min, m_rq2max;

  public:

    // constructor
    PDF_Base();

    // destructor
    virtual ~PDF_Base();

    // member functions
    virtual void AssignKeys(ATOOLS::Integration_Info *const info);

    virtual bool   Collinear(const double kp2) const;
    virtual double Cut(const std::string &type);
    
    void   Calculate(double x,double Q2);
    virtual double GetXPDF(const ATOOLS::Flavour flavour) = 0;

    virtual PDF_Base *GetCopy() = 0;
    virtual PDF_Base *GetBasicPDF();

    void SetBounds();

    virtual void SetPDFMember();

    static void ShowSyntax(const size_t i);

    // inline functions
    inline double XMin() const  { return m_xmin;  }
    inline double XMax() const  { return m_xmax;  }
    inline double Q2Min() const { return m_rq2min; }
    inline double Q2Max() const { return m_rq2max; }

    inline const ATOOLS::Flavour_Set &Partons() const { return m_partons; }
    inline const ATOOLS::Flavour      Bunch() const   { return m_bunch;   }

    inline const std::string &Type() const { return m_type;     }  
    inline double         Exponent() const { return m_exponent; }

    inline PDF_Base *operator->() { return GetBasicPDF(); }

    inline void Extract(const ATOOLS::Flavour flavour,const double x)
    { GetBasicPDF()->SingleExtract(flavour,x); }
    inline void Reset()
    { GetBasicPDF()->SingleReset(); }

    inline bool Contains(const ATOOLS::Flavour &a) const
    { return m_partons.find(a)!=m_partons.end(); }

    inline double RescaleFactor() const { return m_rescale; }

    inline int    OrderAS() const { return m_orderas; }
    inline double ASMZ() const    { return m_asmz;    }

    inline int Member() const { return m_member; }

  };
  
  /*!
    \class PDF_Base
    \brief is the base class for parton distributions and structure functions.

    This class is the base class for both parton distributions and structure functions
    (e.g. of leptons). It hosts some information about both the incoming bunch particle 
    and the partons that constitute it. Furthermore, the corresponding weights for 
    given energy fractions and scales are calculated and made available through this class.
  */
  /*!
    \var std::vector<ATOOLS::Flavour> PDF_Base::m_partons    
    The vector listing the constituents. for protons this list has usually eleven members,
    one gluon and the five lightest quarks and anti-quarks. For leptons, due to the very nature
    of the structure function, only the lepton in question is included as legitimate parton.
  */
  /*!
    \var ATOOLS::Flavour PDF_Base::m_bunch
    The incoming particle. 
  */
  /*!
    \var std::string PDF_Base::m_type
    The type of the PDF_Base. It is given by PDF_flavourname.
  */
  /*!
    \fn virtual PDF_Base * PDF_Base::GetCopy()
    A method to initialize another PDF as exact copy of the current one.
    This is needed for the initial state shower of APACIC.
  */
  /*!
    \fn void PDF_Base::Calculate(const double, const double)
    This calculates all parton densities inside a hadron or the structure function
    depending on the energy fraction \f$x\f$ and the scale squared \f$Q^2\f$. It should
    be noted that the result will be given - as usual - as \f$x\f$ times the function. 
  */
  /*!
    \fn double PDF_Base::GetXPDF(const ATOOLS::Flavour &)
    This returns x times one specific pdf out of the set of all calculated ones - it checks 
    whether the flavour passed over as argument is included in the list of legitimate constituents, 
    PDF_Base::m_partons and selects the appropriate one.
  */
} 

#define DECLARE_PDF_GETTER(NAME)					\
  									\
  class NAME: public ATOOLS::Getter_Function				\
  <PDF_Base,PDF_Arguments,std::less<std::string> > {			\
  protected:								\
    std::string m_key;							\
    void PrintInfo(std::ostream &str,const size_t width) const;		\
    Object_Type *							\
      operator()(const Parameter_Type &parameters) const;		\
  public:								\
    NAME(const std::string &name):					\
      ATOOLS::Getter_Function<PDF_Base,PDF_Arguments,			\
      std::less<std::string> >(name), m_key(name) {}			\
  }

#endif
