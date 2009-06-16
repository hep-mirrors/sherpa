#include "PDF/Main/Structure_Function.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace PDF;
using namespace ATOOLS;

Structure_Function::Structure_Function(PDF::PDF_Base * _p_pdf,ATOOLS::Flavour _m_bunch):
  ISR_Base(_p_pdf), m_lastx(-1.0), m_lastq2(-1.0)
{
  m_bunch = _m_bunch;
  m_type  = std::string("(SF)");
}

bool Structure_Function::CalculateWeight(double x,double z,double kp2,double q2) 
{
  if (IsEqual(x,m_lastx) && IsEqual(q2,m_lastq2)) return m_lastres;
  m_lastx=x;
  m_lastq2=q2;  
  if ( (x  > p_pdf->XMax()) || (x<= p_pdf->XMin()) ) {
    msg_Error()<<"SF::CalculateWeight : x out of bounds "<<x<<" at "<<q2<<", "
	       <<"xrange = "<<p_pdf->XMin()<<" ... "<<p_pdf->XMax()<<std::endl;
    return m_lastres=0; 
  }
  if ( (q2 >= p_pdf->Q2Max()) || (q2<= p_pdf->Q2Min()) ) { 
    msg_Error()<<"SF::CalculateWeight : q2 out of bounds "<<x<<" at "<<q2<<", "
	       <<"q2range = "<<p_pdf->Q2Min()<<" ... "<<p_pdf->Q2Max()<<std::endl;
    return m_lastres=0; 
  }
  p_pdf->Calculate(x,q2);
  m_weight = 1./x;
  return m_lastres=1;
}

double Structure_Function::Weight(ATOOLS::Flavour flin)
{
  return m_weight * p_pdf->GetXPDF(flin); 
}
