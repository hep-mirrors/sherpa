#include "Structure_Function.H"
#include "PDF_Handler.H"
#include "Message.H"

using namespace ISR;
using namespace PDF;
using namespace APHYTOOLS;
using namespace AORGTOOLS;


Structure_Function::Structure_Function(Flavour beam)
{
  type = std::string("(SF)");

  PDF_Handler pdf_hdl;
  pdf  = pdf_hdl.GetPDFLib(beam);

  msg.Tracking()<<"Initialised structure function for beam "<<beam<<std::endl;
}

Structure_Function::~Structure_Function() {
  if (pdf) { delete pdf; pdf = NULL; }   
}

bool Structure_Function::CalculateWeight(double x,double q2) 
{
  if ( (x  > pdf->GetXMax()) || (x<= pdf->GetXMin()) ) {
    msg.Error()<<"SF : x out of bounds "<<x<<" at "<<q2<<", "
	       <<"xrange = "<<pdf->GetXMin()<<" ... "<<pdf->GetXMax()<<std::endl;
    return 0; 
  }
  if ( (q2 >= pdf->GetQ2Max()) || (q2<= pdf->GetQ2Min()) ) { 
    msg.Error()<<"SF : q2 out of bounds "<<x<<" at "<<q2<<", "
	       <<"q2range = "<<pdf->GetQ2Min()<<" ... "<<pdf->GetQ2Max()<<std::endl;
    return 0; 
  }
  
  pdf->Calculate(x,q2);
  weight = 1./x;
  return 1;
};

double Structure_Function::Weight(Flavour flin)
{
  return weight * pdf->GetXPDF(flin); 
}
