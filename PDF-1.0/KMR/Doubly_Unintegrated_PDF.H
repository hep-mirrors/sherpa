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

    std::map<ATOOLS::Flavour,APACIC::Splitting_Function*> m_splitting;
    std::map<ATOOLS::Flavour,LL_Branching*>               m_branching;

    double m_x, m_z, m_kperp2, m_mu2, m_mu02, m_epsilon;
    double m_integrated, m_unintegrated;

    kps::type m_kperpscheme;
    double    m_cutrescale;
    bool      m_calculate;

    double SmoothIntegrated(const ATOOLS::Flavour flavour);
    double ConstantIntegrated(const ATOOLS::Flavour flavour);

    bool Unintegrate(const ATOOLS::Flavour flavour);

  public:

    // constructor
    Doubly_Unintegrated_PDF(PDF_Base *pdf,MODEL::Running_AlphaS *alphas,
			    const double mu02);

    // destructor
    ~Doubly_Unintegrated_PDF();

    // getters 
    void    Output();
    double  Cut(const std::string &type);

    // member functions
    void AssignKeys(ATOOLS::Integration_Info *const info);
    bool Collinear(const double kp2=0.0) const;

    void   Calculate(double x,double z,double kperp2,double mu2);
    double GetXPDF(const ATOOLS::Flavour flavour);

    PDF_Base *GetCopy();

    PDF_Base *GetBasicPDF();

    // inline functions
    inline kps::type KPerpScheme() const { return m_kperpscheme; }
    inline double    CutRescale() const  { return m_cutrescale;  }

    inline void SetKPerpScheme(const kps::type scheme) 
    { m_kperpscheme=scheme; }
    inline void SetCutRescale(const double rescale) 
    { m_cutrescale=rescale; }

  };// end of class Doubly_Unintegrated_PDF

}// end of namespace PDF

#endif
