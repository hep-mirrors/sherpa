#include "PDF/Main/Structure_Function.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace PDF;
using namespace ATOOLS;

Structure_Function::Structure_Function(PDF::PDF_Base * _p_pdf,const ATOOLS::Flavour& _m_bunch):
  ISR_Base(_p_pdf)
{
  m_bunch = _m_bunch;
  if (m_bunch.IsChargedLepton()) m_type = isrtype::lepton;
  else if (m_bunch.IsPhoton() ||
	   m_bunch.IsHadron())   m_type = isrtype::hadron;
  //m_type  = std::string("(SF)");
}

bool Structure_Function::CalculateWeight(double x,double z,double kp2,double q2,int warn) 
{
  if ( (x  > p_pdf->XMax()) || (x<= p_pdf->XMin()) ) {
    if (warn)
      msg_Error() << METHOD << ": x out of bounds: " << x << " at Q2 = " << q2
                  << ", "
                  << "xrange = " << p_pdf->XMin() << " ... " << p_pdf->XMax()
                  << std::endl;
    return false;
  }
  if ( (q2 >= p_pdf->Q2Max()) || (q2<= p_pdf->Q2Min()) ) {
    if (warn)
      msg_Error() << METHOD << ": q2 out of bounds " << x << " at " << q2
                  << ", "
                  << "q2range = " << p_pdf->Q2Min() << " ... " << p_pdf->Q2Max()
                  << std::endl;
    return false;
  }
  p_pdf->Calculate(x,q2);
  m_weight=1.0/x;
  return true;
}

double Structure_Function::Weight(ATOOLS::Flavour flin)
{
  if (!p_pdf->Contains(flin)) return 0.0;
  return m_weight * p_pdf->GetXPDF(flin); 
}
