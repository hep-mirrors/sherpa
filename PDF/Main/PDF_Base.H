#ifndef PDF_Main_PDF_Base_H
#define PDF_Main_PDF_Base_H

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Org/Settings.H"

namespace PDF {

  struct PDF_Defaults {
    std::map<kf_code,std::string> m_deflib, m_defset;

    PDF_Defaults();

    inline std::string DefaultPDFLibrary(kf_code kf) { return m_deflib[kf]; }
    inline std::string DefaultPDFSet(kf_code kf)     { return m_defset[kf]; }
  };

  extern PDF_Defaults * pdfdefs;

  struct PDF_Arguments {
    ATOOLS::Flavour       m_bunch;
    int m_ibeam;
    std::string m_set;
    int m_member;
    int m_order, m_scheme;
    PDF_Arguments(const ATOOLS::Flavour &bunch,
                  int ibeam,
                  std::string set,
                  int member,
                  int order=-1,
                  int scheme=-1):
      m_bunch(bunch),
      m_ibeam(ibeam),
      m_set(set),
      m_member(member),
      m_order(order),
      m_scheme(scheme)
    {}
  };

  struct PDF_Flavour: public ATOOLS::Flavour {
    double m_mass, m_thres;
    PDF_Flavour(const kf_code &kf=kf_none,const bool anti=false):
      Flavour(kf,anti), m_mass(0.0), m_thres(0.0) {}
  };

  struct PDF_AS_Info {
    int m_order,m_nf;
    double m_asmz,m_mz2;
    std::vector<PDF_Flavour> m_flavs;
    inline PDF_AS_Info(): m_order(-1), m_nf(-1), m_asmz(0.0), m_mz2(0.0) {}
  };

  std::ostream &operator<<(std::ostream &ostr,const PDF_AS_Info &asi);

  class PDF_Base {
  public:
    typedef ATOOLS::Getter_Function<PDF_Base, PDF_Arguments>
            PDF_Getter_Function;

  protected:

    std::vector<PDF_Base*> m_copies;

    ATOOLS::Flavour_Set m_partons;
    ATOOLS::Flavour     m_bunch;

    std::string m_type, m_set;
    int m_member, m_lhef_number, m_nf;

    double m_exponent, m_rescale;
    double m_xmin, m_xmax, m_q2min, m_q2max;
    bool   m_rescX;

    PDF_AS_Info m_asinfo;

    virtual void CalculateSpec(const double& x,const double& Q2) = 0;
  private:
    double m_rq2min, m_rq2max;

    static void RegisterDefaults() ;

  public:
    PDF_Base();
    virtual ~PDF_Base();

    void   Calculate(double x,double Q2);
    virtual double GetXPDF(const ATOOLS::Flavour& flavour) = 0;
    virtual double GetXPDF(const kf_code& kf, bool anti) = 0;
    virtual double AlphaSPDF(const double &q2);
    inline  double GetXPDF(const int& id) { return GetXPDF(abs(id),id<0); }

    virtual PDF_Base *GetCopy() = 0;
    virtual PDF_Base *GetBasicPDF();

    void SetBounds();

    virtual void SetAlphaSInfo();
    virtual void SetPDFMember();
    virtual double GetDefaultAlpha();
    virtual double GetDefaultScale();
    virtual int GetFlavourScheme();

    static void ShowSyntax(size_t i);

    inline const std::string& Type() const { return m_type; }
    inline std::string Set() const                    { return m_set; }
    inline ATOOLS::Flavour      Bunch() const   { return m_bunch;   }
    inline const ATOOLS::Flavour_Set &Partons() const { return m_partons; }
    inline bool Contains(const ATOOLS::Flavour & flav) const {
      if (m_partons.find(flav) != m_partons.end()) return true;
      for (size_t i(0); i < flav.Size(); ++i)
        if (m_partons.find(flav[i]) == m_partons.end()) return false;
      return true;
    }
    inline int Member() const         { return m_member; }
    inline int LHEFNumber() const     { return m_lhef_number; }
    inline PDF_Base *operator->()     { return GetBasicPDF(); }
    inline PDF_AS_Info ASInfo() const { return m_asinfo; }
    inline double Exponent() const    { return m_exponent; }
    inline double XMin() const        { return m_xmin;  }
    inline double XMax() const        { return m_xmax;  }
    inline double Q2Min() const       { return m_rq2min; }
    inline double Q2Max() const       { return m_rq2max; }
    inline void SetRescaleFactor(const double & fac) { m_rescale = fac; }
    inline void ResetRescaleFactor()                 { m_rescale = 1.; }
    inline void SetRescaling(const bool & resc=true) { m_rescX   = resc; }
    inline double RescaleFactor() const              { return m_rescale; }
  };

  /*!
    \class PDF_Base
    \brief is the base class for parton distributions and structure functions.

    This class is the base class for both parton distributions and structure
    functions (e.g. of leptons). It hosts some information about both the
    incoming bunch particle and the partons that constitute it. Furthermore, the
    corresponding weights for given energy fractions and scales are calculated
    and made available through this class.
  */
  /*!
    \var std::vector<ATOOLS::Flavour> PDF_Base::m_partons
    The vector listing the constituents. for protons this list has usually
    eleven members, one gluon and the five lightest quarks and anti-quarks. For
    leptons, due to the very nature of the structure function, only the lepton
    in question is included as legitimate parton.
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
    This calculates all parton densities inside a hadron or the structure
    function depending on the energy fraction \f$x\f$ and the scale squared
    \f$Q^2\f$. It should be noted that the result will be given - as usual - as
    \f$x\f$ times the function.
  */
  /*!
    \fn double PDF_Base::GetXPDF(const ATOOLS::Flavour &)
    This returns x times one specific pdf out of the set of all calculated ones
    - it checks whether the flavour passed over as argument is included in the
    list of legitimate constituents, PDF_Base::m_partons and selects the
    appropriate one.
  */
}

#define DECLARE_PDF_GETTER(NAME)					\
  									\
  class NAME: public ATOOLS::Getter_Function				\
  <PDF::PDF_Base,PDF::PDF_Arguments,std::less<std::string> > {			\
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
