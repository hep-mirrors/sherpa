#ifndef Continued_PDF_H
#define Continued_PDF_H

#include "PDF_Base.H"

namespace PDF {

  struct pcs {
    
    enum type {
      none  = 0
    };

  };// end of struct pcs

  class Continued_PDF: public PDF_Base {
  private:

    PDF_Base *p_pdf;

    pcs::type m_scheme;

    double m_q2cut, m_xcut;

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


    // inline functions
    inline void SetQ2Cut(const double q2cut) { m_q2cut=q2cut; }
    inline void SetXCut(const double xcut)   { m_xcut=xcut;   }

    inline double Q2Cut() const { return m_q2cut; }
    inline double XCut() const  { return m_xcut;  }
    
  };// end of class Contibued_PDF

}// end of namespace PDF

#endif
