#include "Backward_Splitting_Group.H"

using namespace APACIC;

Backward_Splitting_Group::Backward_Splitting_Group(Splitting_Function * spl, PDF::PDF_Base * _pdf): 
  p_pdf(_pdf),Splitting_Group(spl) {
  //    pdf = new PDF::PDF_MRST();  // initialising proton;
}

double Backward_Splitting_Group::CrudeInt(double _zmin, double _zmax) {
  if (!p_partsums) p_partsums = new double[m_group.GetLength()];
  m_lastint = 0;
  int i     = 0;
  for (SplFunIter iter(m_group);iter();++iter,++i) {
    if (p_pdf->GetXPDF(iter()->GetFlB())==0.) {
      p_partsums[i]=0;
    } 
    else {
      p_partsums[i] = m_lastint += 
	iter()->CrudeInt(_zmin,_zmax) *  
	p_pdf->GetXPDF(iter()->GetFlA()) /
	p_pdf->GetXPDF(iter()->GetFlB());
    }
    
//     std::cout<<"In CrudeInt("<<_zmin<<","<<_zmax<<") : "<<iter()->GetFlA()<<" ->"<<iter()->GetFlB()<<std::endl
// 	     <<"   "<<iter()->CrudeInt(_zmin,_zmax)<<" * "  
// 	     <<p_pdf->GetXPDF(iter()->GetFlA())<<" / "
// 	     <<p_pdf->GetXPDF(iter()->GetFlB())<<" for "
// 	     <<iter()->GetFlA()<<" "<<iter()->GetFlB()<<std::endl;
    
    
  }
  return m_lastint;
}        
