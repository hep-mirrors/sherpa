#ifndef Continued_PDF_H
#define Continued_PDF_H

#include "PDF_Base.H"
#include <map>

namespace PDF {

  struct pcs {
    
    enum type {
      linear = 1,
      cubic  = 2
    };

  };// end of struct pcs

  class Continued_PDF: public PDF_Base {
  private:

    PDF_Base *p_pdf;

    pcs::type m_scheme;
    
    bool m_continue;

    double m_x, m_z, m_kperp2, m_mu2, m_epsilon;

    double ContinueLinear(const ATOOLS::Flavour flavour);
    double ContinueSquarish(const ATOOLS::Flavour flavour);

  public:

    // constructor
    Continued_PDF(PDF_Base *const pdf,const pcs::type scheme);

    // destructor
    ~Continued_PDF();
    
    // member functions
    void   Output();
    double Cut(const std::string &type);

    void AssignKeys(ATOOLS::Integration_Info *const info);
    bool Collinear(const double kp2=0.0) const;

    void   Calculate(double x,double z,double kperp2,double mu2);
    double GetXPDF(const ATOOLS::Flavour flavour);

    PDF_Base *GetCopy();

    PDF_Base *GetBasicPDF();

  };// end of class Contibued_PDF

}// end of namespace PDF

#endif
