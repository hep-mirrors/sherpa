#include "Strong_Coupling.H"
#include "Hadronisation_Parameters.H"
#include "Running_AlphaS.H"
#include "Run_Parameter.H"
#include "Message.H"

using namespace AHADIC;
using namespace ATOOLS;

Strong_Coupling::Strong_Coupling(const asform::code form) :
  m_form(form)
{
  m_pt02    = hadpars.Get(std::string("pt02"));
  m_beta0   = 12.*M_PI/(33.-2.*3.);
  MODEL::Running_AlphaS * as = static_cast<MODEL::Running_AlphaS *>
    (rpa.gen.GetScalarFunction(std::string("alpha_S")));
  m_Lambda2 = m_pt02*exp(-1./(m_beta0*(*as)(m_pt02)));
  m_kappa2  = m_pt02/exp(1.)*m_beta0*log(m_pt02/m_Lambda2);
  if (msg_LevelIsDebugging()) {
    std::ofstream was;
    was.open("as_in_ahadic_test.dat");
    for (double Q(0.01);Q<10.;Q*=1.001) {
      PRINT_INFO(Q<<" "<<(*this)(sqr(Q))<<" : "<<(*this)(m_pt02)<<" ("<<sqrt(m_pt02)<<").");
      was<<Q<<" "<<(*this)(sqr(Q))<<" : "<<(*this)(m_pt02)<<" ("<<sqrt(m_pt02)<<").\n";
    }
    was.close();
    //abort();
  }
  m_asmax = (*this)(m_pt02);
}

const double Strong_Coupling::operator()(double q2) const {
  if (q2<m_pt02) {
    switch (m_form) {
    case (asform::fall_off) : return q2/m_kappa2*exp(-q2/m_pt02);
    case (asform::constant) :
    default:                  return m_asmax;
    }
  }
  return 1./(m_beta0*log(q2/m_Lambda2));
}
