#include "PDF/Main/Structure_Function.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace PDF;
using namespace ATOOLS;

Structure_Function::Structure_Function(PDF::PDF_Base * _p_pdf,ATOOLS::Flavour _m_bunch):
  ISR_Base(_p_pdf), m_x(0.), m_q2(0.)
{
  m_bunch = _m_bunch;
  m_type  = std::string("(SF)");
}

bool Structure_Function::CalculateWeight(double x,double z,double kp2,double q2) 
{
  if ( (x  > p_pdf->XMax()) || (x<= p_pdf->XMin()) ) {
    msg_Error()<<"SF::CalculateWeight : x out of bounds "<<x<<" at "<<q2<<", "
	       <<"xrange = "<<p_pdf->XMin()<<" ... "<<p_pdf->XMax()<<std::endl;
    return 0; 
  }
  if ( (q2 >= p_pdf->Q2Max()) || (q2<= p_pdf->Q2Min()) ) { 
    msg_Error()<<"SF::CalculateWeight : q2 out of bounds "<<x<<" at "<<q2<<", "
	       <<"q2range = "<<p_pdf->Q2Min()<<" ... "<<p_pdf->Q2Max()<<std::endl;
    return 0; 
  }
  if (x!=m_x||q2!=m_q2) {
    p_pdf->Calculate(x,q2);
    m_x=x;
    m_q2=q2;
    m_weight = 1./m_x;
  }
  return 1;
}

double Structure_Function::Weight(ATOOLS::Flavour flin)
{
  return m_weight * p_pdf->GetXPDF(flin); 
}
