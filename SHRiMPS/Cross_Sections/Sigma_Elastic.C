#include "SHRiMPS/Cross_Sections/Sigma_Elastic.H"
#include "SHRiMPS/Tools/Special_Functions.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

using namespace SHRIMPS;
using namespace ATOOLS;

double Sigma_Elastic::dSigma_dt::operator()(double B) {
  // factors:
  // * B         from B integration in polar coordinates
  // * 2 pi J0   from integral of angle between Q and B
  // * (arg)     usual argument of elastic scattering
  return B * 2.*M_PI*SF.Jn(0,B*m_Q) * p_sigma_el->GetDiffArgument(B);
}

Sigma_Elastic::Sigma_Elastic() :
  m_tmin(0), m_tmax(1.), m_summed(0.), m_steps(1000),
  m_delta((m_tmax-m_tmin)/m_steps)
{ }


double Sigma_Elastic::GetValue(const double & B) { 
  return ATOOLS::sqr(p_eikonal->Prefactor()*(1.-exp(-(*p_eikonal)(B)/2.))); 
}

double Sigma_Elastic::GetCombinedValue(const double & B) { 
  return ATOOLS::sqr(GetDiffArgument(B));
}

double Sigma_Elastic::GetDiffArgument(const double & B) { 
  double value(0.);
  for (size_t i=0;i<p_eikonals->size();i++) {
    for (size_t j=0;j<(*p_eikonals)[i].size();j++) {
      Omega_ik * eikonal = (*p_eikonals)[i][j];
      value += eikonal->Prefactor()*(1.-exp(-(*eikonal)(B)/2.)); 
    }
  }
  return value;
}

void Sigma_Elastic::FillGrids() {
  m_diffgrid.clear();
  m_intgrid.clear();
  FillDiffQGrid();
  FillIntQGridAndNormalize();
}

void Sigma_Elastic::FillDiffQGrid() {
  msg_Out()<<METHOD<<" for ["<<m_tmin<<", "<<m_tmax<<"] in "<<m_steps<<" steps of "
	   <<"size = "<<m_delta<<"\n";
  dSigma_dt differential(this);
  Gauss_Integrator integrator(&differential);

  double value, t;
  for (size_t i=0;i<m_steps;i++) {
    t = m_tmin + m_delta*i;
    differential.SetQ(sqrt(t));
    value = rpa->Picobarn()/(4.*M_PI) *
      ATOOLS::sqr(integrator.Integrate(0.,MBpars.GetEikonalParameters().bmax,
				       MBpars.GetEikonalParameters().accu,1.));
    if (dabs(value<0.)) value = 0.;
    m_diffgrid.push_back(value);
  }
}

double Sigma_Elastic::GetXSvsT(double t) {
  size_t i = 0;
  size_t l_ind(i), r_ind;
  double t_current(m_tmin + m_delta*i);
  double l_dist;
  double r_dist;
  if (t >= m_tmin && t <= m_tmax) {
    while (t_current <= t) {
      l_dist = dabs(t_current - t);
      l_ind = i;
      i++;
      t_current = m_tmin + m_delta*i;
    }
    r_dist = dabs(t_current - t);
    r_ind = i;
    double value;
    if (r_dist < l_dist) value = m_diffgrid[l_ind];
    else if (l_dist < r_dist) value = m_diffgrid[r_ind];
    else value = (m_diffgrid[l_ind] + m_diffgrid[r_ind])/2.;
    return value;
  }
  else {
    msg_Error() << "Error in " << METHOD << " t value out of range";
    return 0.;
  }
}

void Sigma_Elastic::FillIntQGridAndNormalize() {
  m_intgrid.push_back(0.);
  m_summed = 0.;
  double average, binvalue;
  for (size_t i=1;i<m_steps;i++) {
    average   = (m_diffgrid[i-1]+m_diffgrid[i])/2.;
    binvalue  = average * m_delta;
    if (binvalue<0.) binvalue = 0.;
    m_summed += binvalue;
    m_intgrid.push_back(m_summed);
  }
  for (size_t i=0;i<m_steps;i++) m_intgrid[i] /= m_summed;
}

double Sigma_Elastic::SelectT() const {
  double random(ran->Get());
  unsigned int i(0);
  while (random-m_intgrid[i]>=0) i++;
  return m_tmin+(i-1)*m_delta + m_delta *(random-m_intgrid[i-1])/(m_intgrid[i]-m_intgrid[i-1]);
}

double Sigma_Elastic::Test() {
  const Eikonal_Parameters & eikparams(MBpars.GetEikonalParameters());
  const FormFactor_Parameters & ffparams(MBpars.GetFFParameters());
  const double EulerGamma= 0.577215664901532860606512090082 ;
  double a(ffparams.Lambda2/(8.*(1.+ffparams.kappa)));
  double c(eikparams.beta02*ffparams.Lambda2*(1.+ffparams.kappa)*
	   exp(2.*eikparams.Delta*eikparams.Ymax)/(8.*M_PI));
  double alpha(2.*M_PI*ffparams.norm);
  ExpInt expint;
  double ei(expint.GetExpInt(-c)), ei2(expint.GetExpInt(-c/2.));
  return alpha*(EulerGamma+ei-ei2+log(c/4.))/(2.*a)*rpa->Picobarn();
}


