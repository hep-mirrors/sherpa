#include "AMISIC++/Perturbative/MI_Processes.H"
#include "AMISIC++/Tools/Impact_Parameter.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Message.H"

using namespace AMISIC;
using namespace ATOOLS;


// All equations in this file refer to 
// Sjostrand-van der Zijl, PRD 36 (1987) 2019.

Impact_Parameter::Impact_Parameter() : p_pint(new Interaction_Probability()),
				       p_mo(p_pint->GetOverlap()),
				       m_test(false) {}

void Impact_Parameter::Initialize(const double & xsecratio) {
  p_pint->Initialize(xsecratio);
  m_oexp = p_pint->OverlapExpectation();
  m_fc   = m_oexp/p_mo->Integral()*p_pint->Integral();
  m_bmax = p_mo->Bmax();
  if (m_test) Test();
}

double Impact_Parameter::operator()(const double & b) {
  // This is f(b), the enhancement factor
  return m_enhancement = (b<m_bmax? (*p_mo)(b)/m_oexp : 0.);
}

double Impact_Parameter::SelectB(const double & pt2) {
  // Select b according to f(b) and accept or reject b with probability given by
  // "factorized Sudakov form factor", Eq. (37)
  double hardpart = p_procs->SudakovArgument(pt2);
  double sudakov, softpart;
  int    trials   = 1000;
  do {
    softpart = m_fc * (*this)(m_b = p_mo->SelectB());
    sudakov  = exp(-softpart * hardpart);
  } while (sudakov<ran->Get() && (trials--)>0);
  if (trials<=0)
    msg_Error()<<METHOD<<" throws warning:\n"
	       <<"   no impact parameter in accordance with Sudakov.\n"
	       <<"   Return b = "<<m_b<<" for pt^2 = "<<pt2
	       <<" without Sudakov argument.\n";
  return m_b;
}

//##########################################################################################
//##########################################################################################
//##########################################################################################
//##########################################################################################
//##########################################################################################

void Impact_Parameter::Test() {
  msg_Out()<<METHOD<<" starts testing enhancement factor.\n";
  Histogram histoOverlap(0,0.,m_bmax,100);
  double b(0.), btot(0.), bstep(m_bmax/100.);
  while (b<m_bmax) {
    histoOverlap.Insert(b+bstep/2.,(*p_mo)(b+bstep/2.));
    b+= bstep;
  }
  histoOverlap.Output("Overlap.dat");

  Histogram histoPInt(0,0.,m_bmax,100);
  b = 0.;
  while (b<m_bmax) {
    histoPInt.Insert(b+bstep/2.,(*p_pint)(b+bstep/2.));
    b+= bstep;
  }
  histoPInt.Output("PInt.dat");

  Histogram histoBWeight(0,0.,m_bmax,100);
  b = 0.;
  while (b<m_bmax) {
    histoBWeight.Insert(b+bstep/2.,(*this)(b+bstep/2.));
    b+= bstep;
  }
  histoBWeight.Output("Enhancement_Factor.dat");
  
  msg_Out()<<METHOD<<" starts testing b selection.\n";
  double ntrials = 2.5e7;
  Histogram histoB(0,0.,1.1*m_bmax,100);
  for (long int i=0;double(i)<ntrials;i++) {
    b = p_mo->SelectB();
    histoB.Insert(b);
  }
  histoB.Finalize();
  histoB.Output("B_Distribution.dat");

  Histogram histoB10(0,0.,1.1*m_bmax,100);
  for (long int i=0;double(i)<ntrials;i++) {
    b = SelectB(100.);
    histoB10.Insert(b,1.);
  }
  histoB10.Finalize();
  histoB10.Output("B_Distribution_10.dat");

  Histogram histoB100(0,0.,1.1*m_bmax,100);
  for (long int i=0;double(i)<ntrials;i++) {
    b = SelectB(10000.);
    histoB100.Insert(b,1.);
  }
  histoB100.Finalize();
  histoB100.Output("B_Distribution_100.dat");
}

