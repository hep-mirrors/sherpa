#ifndef Doubly_Unintegrated_PDF_H
#define Doubly_Unintegrated_PDF_H

#include "PDF_Base.H"
#include "LL_Sudakov.H"

namespace PDF {

  struct kps {
    
    enum type {
      constant   = 0,
      function   = 1,
      derivative = 2
    };
    
  };// end of struct kps
  
  class Doubly_Unintegrated_PDF: public PDF_Base {
  private:

    PDF_Base   *p_pdf;
    LL_Sudakov *p_sudakov;

    MODEL::Running_AlphaS *p_alphas;

    std::map<ATOOLS::Flavour,Splitting_Kernel*> m_splitting;
    std::map<ATOOLS::Flavour,LL_Branching*>     m_branching;

    double m_x, m_z, m_kperp2, m_mu2, m_mu02, m_epsilon;
    double m_integrated, m_unintegrated;

    std::vector<std::pair<double,Splitting_Kernel*> > m_partsums;

    kps::type m_kperpscheme;
    double    m_fixedktexponent, m_cweight;
    bool      m_calculate;
    int       m_sudmode, m_splitmode, m_ordermode, m_weightmode;

    Splitting_Kernel *p_jkernel;

    double SmoothIntegrated(const ATOOLS::Flavour flavour);
    double ConstantIntegrated(const ATOOLS::Flavour flavour);

    bool Unintegrate(const ATOOLS::Flavour flavour);

  public:

    // constructor
    Doubly_Unintegrated_PDF(PDF_Base *pdf,MODEL::Running_AlphaS *alphas,
			    const double mu02);

    // destructor
    ~Doubly_Unintegrated_PDF();

    // member functions
    bool Initialize();
    void AssignKeys(ATOOLS::Integration_Info *const info);
    bool Collinear(const double kp2=0.0) const;

    void   Calculate(double x,double z,double kperp2,double mu2);
    double GetXPDF(const ATOOLS::Flavour flavour);

    PDF_Base *GetCopy();
    PDF_Base *GetBasicPDF();

    void    Output();
    double  Cut(const std::string &type);

    bool SelectJetFlavour(ATOOLS::Flavour &a,ATOOLS::Flavour &c,
			  const double &rn);

    void Test();

    // inline functions
    inline kps::type KPerpScheme() const     { return m_kperpscheme;     }
    inline double    FixedKtExponent() const { return m_fixedktexponent; }
    inline int       SudakovMode() const     { return m_sudmode;         }
    inline int       SplitMode() const       { return m_splitmode;       }
    inline int       OrderingMode() const    { return m_ordermode;       }
    inline int       WeightMode() const      { return m_weightmode;      }

    inline LL_Sudakov *Sudakov() const { return p_sudakov; }

    inline void SetKPerpScheme(const kps::type scheme) 
    { m_kperpscheme=scheme; }
    inline void SetFixedKtExponent(const double exponent) 
    { m_fixedktexponent=exponent; }
    inline void SetSudakovMode(const int mode) 
    { m_sudmode=mode; }
    inline void SetSplitMode(const int mode) 
    { m_splitmode=mode; }
    inline void SetOrderingMode(const int mode) 
    { m_ordermode=mode; }
    inline void SetWeightMode(const int mode) 
    { m_weightmode=mode; }

    inline Splitting_Kernel *JetKernel() const { return p_jkernel; }

    inline double CorrectionWeight() const { return m_cweight; }

  };// end of class Doubly_Unintegrated_PDF

}// end of namespace PDF

#endif
