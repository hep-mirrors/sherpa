#include "AddOns/Apacic++/Showers/Backward_Splitting_Group.H"

using namespace APACIC;

Backward_Splitting_Group::
Backward_Splitting_Group(ATOOLS::Mass_Selector *&ms,
			 Splitting_Function *const spl,
			 PDF::PDF_Base *const pdf): 
  Splitting_Group(ms,spl), p_pdf(pdf) {
  if (spl) m_flavs[1]=spl->GetFlB();
}

double Backward_Splitting_Group::CrudeInt(double zmin, double zmax) 
{
  if (m_partsums.empty()) m_partsums.resize(m_splittings.size());
  m_lastint=0.0;
  for (size_t size(m_splittings.size()), i(0);i<size;++i) {
    Splitting_Function *split(m_splittings[i]);
    double xpdfb(p_pdf->GetXPDF(split->GetFlB()));
    if (xpdfb==0.0) {
      m_partsums[i]=0.0;
    } 
    else {
      m_partsums[i]=m_lastint+=
	split->CrudeInt(zmin,zmax)*p_pdf->GetXPDF(split->GetFlA())/xpdfb;
    }
  }
  return m_lastint;
}        
