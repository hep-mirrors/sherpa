#include "Structure_Function.H"
#include "PDF_Handler.H"
#include "Message.H"

using namespace ATOOLS;
using namespace PDF;


Structure_Function::Structure_Function(PDF::PDF_Base * _pdf,Flavour _bunch) :
  p_pdf(_pdf)
{
  m_bunch = _bunch;
  m_type  = std::string("(SF)");
  p_pdf->Output();
}

Structure_Function::~Structure_Function() {
  if (p_pdf) { delete p_pdf; p_pdf = NULL; }   
}

bool Structure_Function::CalculateWeight(double x,double q2) 
{

  //cout<<x<<" : "<<p_pdf->GetXMax()<<endl;

  if ( (x  > p_pdf->GetXMax()) || (x<= p_pdf->GetXMin()) ) {
    msg.Error()<<"SF : x out of bounds "<<x<<" at "<<q2<<", "
	       <<"xrange = "<<p_pdf->GetXMin()<<" ... "<<p_pdf->GetXMax()<<std::endl;
    return 0; 
  }
  if ( (q2 >= p_pdf->GetQ2Max()) || (q2<= p_pdf->GetQ2Min()) ) { 
    msg.Error()<<"SF : q2 out of bounds "<<x<<" at "<<q2<<", "
	       <<"q2range = "<<p_pdf->GetQ2Min()<<" ... "<<p_pdf->GetQ2Max()<<std::endl;
    return 0; 
  }
  
  p_pdf->Calculate(x,q2);
  m_weight = 1./x;
  return 1;
};

double Structure_Function::Weight(Flavour flin)
{
  return m_weight * p_pdf->GetXPDF(flin); 
}
