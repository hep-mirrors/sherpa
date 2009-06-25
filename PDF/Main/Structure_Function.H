#ifndef PDF_Main_Structure_Function_H
#define PDF_Main_Structure_Function_H

#include "PDF/Main/ISR_Base.H"
#include "PDF/Main/PDF_Base.H"

namespace PDF {  

  class Structure_Function : public ISR_Base {
  public:

    // constructor
    Structure_Function(PDF::PDF_Base *,ATOOLS::Flavour); 

    // meber functions
    bool   CalculateWeight(double,double,double,double);
    double Weight(ATOOLS::Flavour);

  };// end of class Structure_Function

  /*!
    \file 
    \brief contains the class PDF::Structure_Function.
  */
  /*!
    \class PDF::Structure_Function
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
    \var PDF_Base * Structure_Function::pdf
    The pointer to the pdf for the structure function. It is handed over to the 
    constructor, after having been instantiated by the PDF::PDF_Handler. Note that also 
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

}// end of namespace PDF

#endif
