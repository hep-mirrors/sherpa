#ifndef PDF_Main_ISR_Base_H
#define PDF_Main_ISR_Base_H

#include "PDF/Main/PDF_Base.H"

namespace PDF {

  struct isrtype {
    enum code {
      intact  = 0,
      lepton  = 1,
      hadron  = 2,
      yfs     = 3,
      unknown = 99
    };
  };
  std::ostream& operator<<(std::ostream&,isrtype::code type);
  
  class ISR_Base { 
  protected:
    PDF_Base *p_pdf;
    //std::string     m_type;
    ATOOLS::Flavour m_bunch;
    isrtype::code   m_type;
    
    double m_weight, m_exponent, m_xmax;
    bool   m_on;

  public:
    ISR_Base(PDF_Base *pdf);
    virtual ~ISR_Base();

    bool On() const   { return m_on;   }

    inline isrtype::code   Type()    const { return m_type;  }
    inline ATOOLS::Flavour Flavour() const { return m_bunch; }
    inline double Exponent()         const { return m_exponent; }
    inline double XMax()             const { return m_xmax;     } 
    inline double XMin()             const { return (p_pdf!=NULL ? p_pdf->XMin() : 0.); }
    inline PDF_Base * PDF()          const { return p_pdf; }
    inline void SetPDF(PDF_Base *pdf)      { p_pdf = pdf;  }
    inline virtual void SetRescaleFactor(const double & rescale) {
      p_pdf->SetRescaleFactor(rescale);
    }
    inline virtual void SetRescaling(const bool & resc=true) {
      p_pdf->SetRescaling(resc);
    }      

    virtual bool   CalculateWeight(double x,double z,double kp2,
				   double Q2,int warn) = 0;
    virtual double Weight(ATOOLS::Flavour) = 0;

  };
}

#endif
