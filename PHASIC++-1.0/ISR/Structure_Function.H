#ifndef Structure_Function_H
#define Structure_Function_H

#include "ISR_Base.H"
#include "PDF_Base.H"

namespace ISR {  
  class Structure_Function : public ISR_Base {
  private:
    /*!
      The pdf for the structure function. It is initialized in the constructor,
      depending on the flavour of the incoming beam. The PDF_Handler uses the beam 
      as arguments to select the corresponding PDF's. I'm sure we could do this
      better ... .
    */
    PDF::PDF_Base * pdf;
  public:
    Structure_Function(APHYTOOLS::Flavour); 
    ~Structure_Function();


    /*!
      Calculates the PDF's at the corresponding x-values and the scale.
    */
    bool   CalculateWeight(double,double);
    /*!
      The parton-flavour dependent weight;
    */
    double Weight(APHYTOOLS::Flavour);
    bool            On()       { return 1;  }
    PDF::PDF_Base * PDF()      { return pdf; }
    double          Exponent() { return pdf->Exponent(); }
    double          Xmax()     { return pdf->GetXMax(); }
  };
}


#endif
