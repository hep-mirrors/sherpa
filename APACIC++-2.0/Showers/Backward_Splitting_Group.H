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
    Backward_Splitting_Group(Splitting_Function * spl, PDF::PDF_Base * _pdf); 
    ~Backward_Splitting_Group() {}
    // summing over all splittings
    double CrudeInt(double _zmin, double _zmax);
  };

} 

#endif 


