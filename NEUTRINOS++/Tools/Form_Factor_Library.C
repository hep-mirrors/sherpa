#include "NEUTRINOS++/Tools/Form_Factor_Library.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include <iomanip>

using namespace NEUTRINOS;
using namespace ATOOLS;
using namespace std;

std::ostream & NEUTRINOS::operator<<(std::ostream & s,const ff_type::code & type) {
  if (type==ff_type::none)             s<<setw(18)<<"none";
  if (type==ff_type::dipole)           s<<setw(18)<<"dipole";
  if (type==ff_type::neutron_electric) s<<setw(18)<<"neutron_electric";
  if (type==ff_type::exponential)      s<<setw(18)<<"exponential";
  if (type==ff_type::Gaussian)         s<<setw(18)<<"Gaussian";
  if (type==ff_type::Kelly)            s<<setw(18)<<"Kelly";
  if (type==ff_type::BBBA)             s<<setw(18)<<"BBBA";
  if (type==ff_type::ArringtonHill)    s<<setw(18)<<"Arrington-Hill";
  if (type==ff_type::Helm)             s<<setw(18)<<"Helm";
  if (type==ff_type::Lovato)           s<<setw(18)<<"Lovato";
  if (type==ff_type::unknown)          s<<setw(18)<<"unknown";
  return s;
}

std::ostream & NEUTRINOS::operator<<(std::ostream & s,const cpl_info::code & cpl) {
  if (cpl==cpl_info::unknown)      s<<setw(12)<<"unknown";
  if (cpl==cpl_info::scalar)       s<<setw(12)<<"scalar";
  if (cpl==cpl_info::pseudoscalar) s<<setw(12)<<"pseudoscalar";
  if (cpl==cpl_info::vector)       s<<setw(12)<<"vector";
  if (cpl==cpl_info::axialvector)  s<<setw(12)<<"axialvector";
  if (cpl==cpl_info::tensor)       s<<setw(12)<<"tensor";
  return s;
}

std::ostream & NEUTRINOS::operator<<(std::ostream & s,const ff_info & info) {
  s<<"     "<<info.m_cpl<<" ["<<info.m_type<<", "<<info.m_params.size()<<" parameters]: \n";
  s<<"                  ";
  if (info.m_type!=ff_type::none && info.m_type!=ff_type::unknown) {
    for (size_t i=0;i<info.m_params.size();i++) s<<info.m_params[i]<<" "; s<<"\n";
  }
  return s;
}
  

/////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////////////
Dipole_Form_Factor::Dipole_Form_Factor(const ff_info & info) :
  Form_Factor_Base("Dipole", info),
  m_norm(info.m_params[0]), m_invlambda2(sqr(1./info.m_params[1])) {}

double Dipole_Form_Factor::Calc(const double & q2) { return m_norm / sqr(1.+m_invlambda2*q2); }




/////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////////////
Neutron_Electric_Form_Factor::Neutron_Electric_Form_Factor(const ff_info & info) :
  Form_Factor_Base("Neutron_Electric", info),
  m_norm(info.m_params[0]), m_invlambda2(sqr(1./info.m_params[1])),
  m_mass2(sqr(info.m_params[2])), m_pref(info.m_params[3]) {}

double Neutron_Electric_Form_Factor::Calc(const double & q2) {
  return m_norm /sqr(1.+m_invlambda2*q2) * q2/(1.+m_pref*q2/(4.*sqr(m_mass2)));
}





/////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////////////
Exponential_Form_Factor::Exponential_Form_Factor(const ff_info & info) :
  Form_Factor_Base("Exponential", info),
  m_norm(info.m_params[0]), m_invlambda2(sqr(1./info.m_params[1])) {}

double Exponential_Form_Factor::Calc(const double & q2) { return m_norm * exp(-m_invlambda2*q2); }





/////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////////////
Gaussian_Form_Factor::Gaussian_Form_Factor(const ff_info & info) :
  Form_Factor_Base("Gaussian", info),
  m_norm(info.m_params[0]), m_invlambda2(sqr(1./info.m_params[1])) {}


double Gaussian_Form_Factor::Calc(const double & q2) { return m_norm * exp(-m_invlambda2*q2); }






/////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////////////
Kelly_Form_Factor::Kelly_Form_Factor(const ff_info & info) :
  Form_Factor_Base("Kelly", info),
  m_norm(info.m_params[0]), m_invlambda2(sqr(1./info.m_params[1])),
  m_mass2(sqr(info.m_params[2])) {
  for (size_t i=0;i<4;i++) m_A[i]  = info.m_params[i+3];
  for (size_t i=0;i<2;i++) m_En[i] = info.m_params[i+7];
}

double Kelly_Form_Factor::Calc(const double & q2) { return ModifiedDipole(q2) + Polynomial(q2); }

const double Kelly_Form_Factor::ModifiedDipole(const double & q2) const {
  double tau = q2/(4.*m_mass2);
  return (m_En[0]*tau/(1.+m_En[1]*tau)) * 1./sqr(1.+m_invlambda2*q2);
}

const double Kelly_Form_Factor::Polynomial(const double & q2) const {
  double tau = q2/(4.*m_mass2);
  return m_norm*(1.+m_A[0]*tau)/(1.+m_A[1]*tau+m_A[2]*tau*tau+m_A[3]*tau*tau*tau);
}






/////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////////////
BBBA_Form_Factor::BBBA_Form_Factor(const ff_info & info) :
  Form_Factor_Base("BBBA", info),
  m_norm(info.m_params[0]), m_mass2(sqr(info.m_params[1])) {
  for (size_t i=0;i<4;i++) { m_A[i]  = info.m_params[i+2]; m_A[i]  = info.m_params[i+6]; }
}

double BBBA_Form_Factor::Calc(const double & q2) {
  double tau = q2/(4.*m_mass2);
  return m_norm * Numerator(tau)/Denominator(tau);
}

const double BBBA_Form_Factor::Numerator(const double & tau) const {
  return m_A[0]+tau*(m_A[1]+tau*(m_A[2]+tau*m_A[3]));
}

const double BBBA_Form_Factor::Denominator(const double & tau) const {
  return 1+tau*(m_B[0]+tau*(m_B[1]+tau*(m_B[2]+tau*m_B[3])));
}





/////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////////////
ArringtonHill_Form_Factor::ArringtonHill_Form_Factor(const ff_info & info) :
  Form_Factor_Base("ArringtonHill", info),
  m_mass2(sqr(info.m_params[0])), m_tcut(info.m_params[1]), m_t0(info.m_params[2]) {
  for (size_t i=0;i<13;i++) { m_A[i]  = info.m_params[i+3]; }
}

double ArringtonHill_Form_Factor::Calc(const double & q2) {
  double z = sqr(sqrt(m_tcut+q2)-sqrt(m_tcut-m_t0))/(q2+m_t0);
  return ZExpand(z);
}

const double ArringtonHill_Form_Factor::ZExpand(const double & z) const {
  double result = 0.;
  for (size_t i=12;i>0;i--) { result += m_A[i]; result *= z; }
  return result + m_A[0];
}




/////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////////////
Helm_Form_Factor::Helm_Form_Factor(const ff_info & info) :
  Form_Factor_Base("Helm", info),
  m_s(info.m_params[0]),
  m_r(sqrt(sqr(1.2*std::cbrt(info.m_params[1]))-5.*sqr(m_s))),
  m_kappa(info.m_params[2]),
  m_hbarc(rpa->hBar()*rpa->c()) {}

double Helm_Form_Factor::Calc(const double & q2) {
  double kappa = sqrt(q2)/m_hbarc, kappaR = kappa*m_r;
  return 3.*exp(-sqr(kappa*m_s)/2.)*(sin(kappaR)-kappaR*cos(kappaR))/pow(kappaR,3);
}




/////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////////////
Lovato_Form_Factor::Lovato_Form_Factor(const ff_info & info) :
  Form_Factor_Base("Lovato", info),
  m_b(info.m_params[0]), m_hbarc(rpa->hBar()*rpa->c()) {
  for (size_t i=0;i<5;i++) m_c[i] = info.m_params[i+1];
}

double Lovato_Form_Factor::Calc(const double & q2) {
  double x = sqrt(q2)/m_hbarc, bx = m_b*x;
  return exp(-sqr(bx)/2.)*Polynomial(bx);
}

const double Lovato_Form_Factor::Polynomial(const double & bx) const {
  double result = 0.;
  for (size_t i=4;i>0;i--) {result += m_c[i]; result *= bx; }
  return (m_c[0]+result)/6.;
}

