#ifndef Structure_Function_H
#define Structure_Function_H

#include "ISR_Base.H"
#include "PDF_Base.H"

namespace ISR {  
  class Structure_Function : public ISR_Base {
  private:
    PDF::PDF_Base * p_pdf;
  public:
    Structure_Function(PDF::PDF_Base *,APHYTOOLS::Flavour); 
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
    PDF::PDF_Base * PDF()      { return p_pdf; }
    double          Exponent() { return p_pdf->Exponent(); }
    double          Xmax()     { return p_pdf->GetXMax(); }
  };
  /*!
    \file 
    \brief contains the class ISR::Structure_Function.
  */
  /*!
    \class ISR::Structure_Function
    \brief This class is used for initial state radiation off particles characterized
           by energy fractions and scales of the hard interaction, i.e. structure functions
	   and PDFs.

    This class is used for treatment of initial state radiation through structure functions 
    or PDFs that are consequently fully incorporated via a pointer to a suitable PDF_Base.
    It should be stressed here that for instance electromagnetic radiation off an electron
    is described through a structure function that is cast into the form of a PDF_Base.
    In any case, this treatment is clearly constrained to collinear initial state radiation.

    The type related to this ISR_Base is "(SF)".
  */
  /*!
    \var PDF::PDF_Base * Structure_Function::pdf
    The pdf for the structure function. It is handed over to the constructor,
    after having been instantiated by the PDF::PDF_Handler. Note that also 
    structure functions like the one for the electron are handed over as PDF_Base.
  */
  /*! 
    \var bool Structure_Function::CalculateWeight(double,double);
    After checking the ranges of the arguments, the energy fraction of the parton and
    the scale, this calculates the parton densities/weights for the structure function
    through a call to p_pdf->Calculate(x,q2). However the result of the pdf is not
    the parton density but \f$x\f$-times the parton density, therefore the actual internal
    weight is set to \f$1/x\f$ to be multiplied later on with the - eventuall flavour
    dependent - pdf piece.
  */
}


#endif
