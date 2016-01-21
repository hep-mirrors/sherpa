#include "SHRiMPS/Cross_Sections/Sigma_Elastic.H"
#include "SHRiMPS/Tools/Special_Functions.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace SHRIMPS;
using namespace ATOOLS;

double Sigma_Elastic::dSigma_dt::operator()(double B) {
  return 2.*M_PI*B*SF.Jn(0,B*m_Q)*p_sigma_el->GetDiffArgument(B);
}

double Sigma_Elastic::GetValue(const double & B) { 
  return ATOOLS::sqr(p_eikonal->Prefactor()*(1.-exp(-(*p_eikonal)(B)/2.))); 
}

double Sigma_Elastic::GetCombinedValue(const double & B) { 
  return ATOOLS::sqr(GetDiffArgument(B));
}

double Sigma_Elastic::GetDiffArgument(const double & B) { 
  double value(0.);
  for (std::list<Omega_ik *>::iterator eikonal=p_eikonals->begin();
       eikonal!=p_eikonals->end(); eikonal++) {
    value += (*eikonal)->Prefactor()*(1.-exp(-(**eikonal)(B)/2.)); 
  }
  return value;
}

void Sigma_Elastic::FillDifferentialGrids() {
  m_diffgrid.clear();
  double Qmax(5.);
  FillDiffQGrid(Qmax);
  m_intgrid.clear();
  double cumul = FillIntQGridAndNormalize();
}

void Sigma_Elastic::FillDiffQGrid(const double & Qmax) {
  dSigma_dt differential(this);
  Gauss_Integrator integrator(&differential);

  double value(1.), pref(0.), Q(Qmax);
  while (dabs((value-pref)/(value+pref))>1.e-12 ||
	 m_diffgrid.size()<10) {
    differential.SetQ(Q);
    pref  = value;
    value =
      sqr(integrator.Integrate(0.,MBpars.GetEikonalParameters().bmax,
			       MBpars.GetEikonalParameters().accu,1.)) *
      rpa->Picobarn()/(4.*M_PI);
    if (dabs(value<1.e-10)) value = 0.;
    m_diffgrid[Q] = value;
    Q *= exp(-1./m_logdelta);
  }
  differential.SetQ(0.);
  m_diffgrid[0.] =
    sqr(integrator.Integrate(0.,MBpars.GetEikonalParameters().bmax,
			     MBpars.GetEikonalParameters().accu,1.)) *
    rpa->Picobarn()/(4.*M_PI);
}

double Sigma_Elastic::FillIntQGridAndNormalize() {
  m_intgrid[0] = 0.;
  double prefQ(0.),Q(0.),prefval(0.),val(0.),cumul(0.);
  for (std::map<double,double>::iterator diffiter=m_diffgrid.begin();
       diffiter!=m_diffgrid.end();diffiter++) {
    if (diffiter->first==0.) {
      prefval = diffiter->second;
      prefQ   = 0.;
      continue;
    }
    Q   = diffiter->first;
    val = diffiter->second;
    cumul += (val+prefval)/2. * (Q-prefQ)*(Q+prefQ);
    m_intgrid[Q] = cumul;
    prefval = val;
    prefQ   = Q;
  }
  for (std::map<double,double>::iterator intiter=m_intgrid.begin();
       intiter!=m_intgrid.end();intiter++) {
    intiter->second /= cumul;
  }
  return cumul;
}


double Sigma_Elastic::SelectPT2() const {
  double random(ran->Get());
  // unsigned int i(0);
  // while (random-m_intgrid[i]>=0) i++;

  // double Q1(sqr(m_Qmax*exp(-double(i-1)/m_logdelta)));
  // double Q2(sqr(i==m_intgrid.size()-1?0.:m_Qmax*exp(-double(i)/m_logdelta)));
  // return ((Q2*(m_intgrid[i-1]-random)+Q1*(random-m_intgrid[i]))/
  // 	  (m_intgrid[i-1]-m_intgrid[i]));
}

void Sigma_Elastic::Test() {
  const double EulerGamma= 0.577215664901532860606512090082 ;
  double m_a,m_c,m_alpha,m_res,m_ei,m_ei2;
  double m_Delta,m_prefactor,m_Lambda2,m_beta0,m_kappa;
  ExpInt m_expint;
  // m_Delta     = (*p_eikonals).front()->Delta();
  // m_prefactor = (*p_eikonals).front()->Prefactor();
  // m_kappa     = (*p_eikonals).front()->Kappa_i();
  // m_Lambda2   = (*p_eikonals).front()->Lambda2();
  // m_beta0     = (*p_eikonals).front()->FF1()->Beta0();
  // m_a         = m_Lambda2/(8.*(1.+m_kappa));
  // m_c         = ATOOLS::sqr(m_beta0)*m_Lambda2*(1.+m_kappa)*
  //   exp(2.*m_Delta*m_Y)/(8.*M_PI);
  // m_alpha     = 2.*M_PI*m_prefactor;
  // m_ei        = m_expint.GetExpInt(-m_c);
  // m_ei2       = m_expint.GetExpInt(-m_c/2.);
  // m_res       = m_alpha*(EulerGamma+m_ei-m_ei2+log(m_c/4.))/(2.*m_a);
  msg_Out() << "In " << METHOD << " sigma_elas = "<< m_res <<" 1/GeV^2 = "
	    <<m_res*rpa->Picobarn()/1.e9<<" mb ."<<std::endl;
}


