#ifndef APACIC_Showers_Backward_Splitting_Group_H
#define APACIC_Showers_Backward_Splitting_Group_H

#include "APACIC++/Showers/Splitting_Group.H"
#include "PDF/Main/PDF_Base.H"

namespace APACIC {

  class Backward_Splitting_Group : public Splitting_Group {
  protected:

    PDF::PDF_Base *p_pdf;

  public:

    // constructor
    Backward_Splitting_Group(ATOOLS::Mass_Selector *&ms,
			     Splitting_Function *const spl, 
			     PDF::PDF_Base *const pdf);

    // member functions
    double CrudeInt(double zmin, double zmax);

  };

  /*!
    \file
    \brief contains the class APACIC::Backward_Splitting_Group
  */

  /*! 
     \class Backward_Splitting_Group
     \brief A group of splitting functions specialised for the 
     backward evolution of the inital state shower
  */

} 

#endif 


