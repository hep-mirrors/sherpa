#include "SHRiMPS/Tools/Kernels.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MathTools.H"

using namespace SHRIMPS;
using namespace ATOOLS;

std::vector<double> & DEQ_Kernel_NoKT::
operator()(const std::vector<double> & input,const double param) {
  if (input.size()!=2) {
    msg_Error()<<"Error in "<<METHOD<<" : "<<std::endl
	       <<"   Wrong input size, expected 2, received "
	       <<input.size()<<"."<<std::endl
	       <<"   Will abort the run."<<std::endl;
    abort();
  }
  double absorpterm(1.);
  if (dabs(m_lambda)>1.e-12) {
    double term0(ATOOLS::Max(1.e-12,m_lambda*m_expfactor*input[0]));
    double term1(ATOOLS::Max(1.e-12,m_lambda*m_expfactor*input[1]));
    switch (m_absorp) {
    case absorption::exponential:
      absorpterm = exp(-(term0+term1));
      break;
    case absorption::factorial:
    default:
      absorpterm = (1.-exp(-term0))/term0 * (1.-exp(-term1))/term1;
      break;
    }
  }
  m_output[0] = +input[0]*m_Delta*absorpterm; 
  m_output[1] = -input[1]*m_Delta*absorpterm; 
  return m_output;
}

Integration_Kernel_B2::
Integration_Kernel_B2(Eikonal_Contributor * Omega1, 
		      Eikonal_Contributor * Omega2, 
		      const int & test) : 
  m_kernel(Integration_Kernel_Theta(Omega1,Omega2,m_test)), 
  m_integrator(ATOOLS::Gauss_Integrator((&m_kernel))),
  m_accu(0.01), m_test(test)
{ }

double Integration_Kernel_B2::operator()(double b1) {
  m_kernel.Setb1(b1);
  double integral(m_integrator.Integrate(0.,M_PI,m_accu,1));
  return 2.*b1*integral;
}

Integration_Kernel_Theta::    
Integration_Kernel_Theta(Eikonal_Contributor * Omega1, 
			 Eikonal_Contributor * Omega2, 
			 const int & test) : 
  p_Omega1(Omega1), p_Omega2(Omega2), m_test(test),
  m_errmax12(0.), m_errmax21(0.), m_maxvalue(0.) {}

double Integration_Kernel_Theta::operator()(double theta1) {
  double b2((m_b==0.)?m_b1:sqrt(m_b*m_b+m_b1*m_b1-2.*m_b*m_b1*cos(theta1)));
  double eik12((*p_Omega1)(m_b1,b2,m_y));
  double eik21((*p_Omega2)(m_b1,b2,m_y));
  double value(eik12*eik21);
  if (m_b1*value>m_maxvalue) m_maxvalue = m_b1*value;
  if (value<0.) {
    msg_Error()<<"Warning in "<<METHOD<<"(B="<<m_b<<", b1="<<m_b1<<", "
	       <<"b2="<<b2<<", theta="<<theta1<<") = "
	       <<value<<" (y="<<m_y<<")"<<std::endl
	       <<"   (eikonals = "
	       <<eik12<<" vs. "<<(p_Omega1->FF1()->FourierTransform(m_b1)*
				  exp(0.3*(p_Omega1->Y()+m_y)))<<", "
	       <<eik21<<" vs. "<<(p_Omega2->FF2()->FourierTransform(b2)*
				  exp(0.3*(p_Omega2->Y()-m_y)))<<")."
	       <<std::endl;
  }
  return value;
}

void Integration_Kernel_Theta::PrintErrors() { 
  msg_Info()<<"Maximal errors in evaluating product of single terms: "
	    <<std::endl<<"    "
	    <<"delta_max{Omega_12} = "<<m_errmax12<<", "
	    <<"delta_max{Omega_21} = "<<m_errmax12<<"."<<std::endl;
}

