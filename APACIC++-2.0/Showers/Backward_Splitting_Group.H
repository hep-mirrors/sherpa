#ifndef Backward_Splitting_Group_H
#define Backward_Splitting_Group_H


#include "PDF_Base.H"
#include "Splitting_Group.H"

//----------------------------------------------------------------------
//     Backward Initial State Shower - Splitting Function Group 
//      ( includes PDFs)
//----------------------------------------------------------------------

namespace APACIC {
  class Backward_Splitting_Group : public Splitting_Group {
  protected:
    PDF::PDF_Base * p_pdf;
  public:
    Backward_Splitting_Group(Splitting_Function * spl, PDF::PDF_Base * pdf); 
    ~Backward_Splitting_Group() {}
    // summing over all splittings
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


