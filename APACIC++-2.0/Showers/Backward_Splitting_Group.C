#include "Backward_Splitting_Group.H"

using namespace APACIC;

Backward_Splitting_Group::Backward_Splitting_Group(Splitting_Function * spl, PDF::PDF_Base * pdf): 
  Splitting_Group(spl), p_pdf(pdf) {
}

double Backward_Splitting_Group::CrudeInt(double zmin, double zmax) {
  if (!p_partsums) p_partsums = new double[m_group.GetLength()];
  m_lastint = 0;
  int i     = 0;
  for (SplFunIter iter(m_group);iter();++iter,++i) {
    if (p_pdf->GetXPDF(iter()->GetFlB())==0.) {
      p_partsums[i]=0.;
    } 
    else {
      p_partsums[i] = m_lastint += 
	iter()->CrudeInt(zmin,zmax) *  
	p_pdf->GetXPDF(iter()->GetFlA()) /
	p_pdf->GetXPDF(iter()->GetFlB());
    }
  }
  return m_lastint;
}        
