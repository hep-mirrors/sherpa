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
    PDF::PDF_Base * pdf;
  public:
    Backward_Splitting_Group(Splitting_Function * spl, PDF::PDF_Base * _pdf): 
      pdf(_pdf),Splitting_Group(spl) {
      //    pdf = new PDF::PDF_MRST();  // initialising proton;
    }

    // summing over all splittings
    double CrudeInt(double _zmin, double _zmax) {
      if (!partsums) partsums = new double[group.GetLength()];
      lastint = 0;
      int i   = 0;
      for (SplFunIter iter(group);iter();++iter,++i) {
	if (pdf->GetXPDF(iter()->GetFlB())==0.) {
	  partsums[i]=0;
	} 
	else {
	  partsums[i] = lastint += 
	    iter()->CrudeInt(_zmin,_zmax) *  
	    pdf->GetXPDF(iter()->GetFlA()) /
	    pdf->GetXPDF(iter()->GetFlB());
	}
	/*
	std::cout<<"In CrudeInt("<<_zmin<<","<<_zmax<<") : "<<std::endl
	    <<"   "<<iter()->CrudeInt(_zmin,_zmax)<<" * "  
	    <<pdf->GetXPDF(iter()->GetFlA())<<" / "
	    <<pdf->GetXPDF(iter()->GetFlB())<<" for "
	    <<iter()->GetFlA()<<" "<<iter()->GetFlB()<<std::endl;
	*/

      }
      return lastint;
    };        
  };

} // namespace APACIC

#endif // Backward_Splitting_Group_H


