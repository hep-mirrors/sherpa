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

Sigma_Elastic::
Sigma_Elastic(std::list<Omega_ik *> * eikonals,const double & energy,
	      const int & test) : 
  Sigma_Base(eikonals),
  m_Bmin(MBpars("bmin")), m_Bmax(MBpars("bmax")),
  m_Qmax(energy/2.), m_logQsteps(300), m_logdelta(20.),m_test(test) 
{ 
  FillGrid(); 
}



void Sigma_Elastic::FillGrid() {
  if (m_test==10) PrintDifferentialelasticXsec(true);
  msg_Tracking()<<"In "<<METHOD<<": Integrate from "
		<<m_Bmin<<" to "<<m_Bmax<<"."<<std::endl
		<<"   Maximal sqrt{|t|} = "<<m_Qmax<<"."<<std::endl;
  m_intgrid.clear();
  m_intgrid.push_back(0.);
  m_diffgrid.clear();


  dSigma_dt differential(this);
  ATOOLS::Gauss_Integrator integrator(&differential);
  double value(1.), pref(0.), prefQ(0.), cumul(0.), Q;
  size_t step(0);
  while (step<m_logQsteps || (step>1 && (value-pref)/(value+pref)>1.e-12)) {
    Q     = m_Qmax*exp(-double(step)/m_logdelta);
    differential.SetQ(Q);
    value = 
      ATOOLS::sqr(integrator.Integrate(m_Bmin,m_Bmax,m_accu,1.))*
      ATOOLS::rpa->Picobarn()/(4.*M_PI);
    msg_Debugging()<<"   Q = "<<Q<<" --> dsigma/dt = "
		   <<(value/1.e9)<<" mbarn"<<std::endl;
    m_diffgrid.push_back(value);
    if (step>0) {
      cumul += (value+pref)/2. * (prefQ-Q);
//     msg_Out()<<"   Q = "<<Q<<" --> cumul = "
// 		   <<(cumul/1.e9)<<std::endl;
      m_intgrid.push_back(cumul);
    }
    prefQ = Q;
    pref  = value;
    step++;
  }
  Q = 0.;
  value = 
    ATOOLS::sqr(integrator.Integrate(m_Bmin,m_Bmax,m_accu,1.))*
    ATOOLS::rpa->Picobarn()/(4.*M_PI);
  m_diffgrid.push_back(value);  
  cumul += (value+pref)/2. * (prefQ-Q);
  m_intgrid.push_back(cumul);
  m_sigma = cumul;
  msg_Debugging()<<"   Q = "<<Q<<" --> dsigma/dt = "
		 <<(value/1.e9)<<" mbarn"<<std::endl;

  for (size_t i=0;i<m_intgrid.size();i++) {
    m_intgrid[i]/=cumul;
    msg_Debugging()<<i<<"  Q = "<<m_Qmax*exp(-double(i)/m_logdelta)
		   <<" --> "<<m_intgrid[i]<<"."<<std::endl;
  }
}

void Sigma_Elastic::PrintDifferentialelasticXsec(const bool & onscreen,
			const bool & tuning, std::string dirname) {
  if(!tuning){
    std::ofstream was;
    std::string Estring(ATOOLS::ToString(2.*m_Qmax));
    std::string filename(dirname+std::string("/Dsigma_el_by_dt_"+Estring+".dat"));
    was.open(filename.c_str());

    dSigma_dt differential(this);
    ATOOLS::Gauss_Integrator integrator(&differential);
    double value(1.), Q(m_Qmax);
    size_t step(0);
    if (onscreen) msg_Out()<<"---------------------------------------------\n";
    while (Q>1.e-3) {
      Q     = m_Qmax*exp(-double(step)/20.);
      differential.SetQ(Q);
      value = 
        ATOOLS::sqr(integrator.Integrate(m_Bmin,m_Bmax,m_accu,1.))/(4.*M_PI) *
        ATOOLS::rpa->Picobarn()/1.e9;
      was<<" "<<(Q*Q)<<"   "<<value/(2.*Q)<<std::endl;
      if (onscreen) msg_Out()<<" "<<(Q*Q)<<"   "<<value/(2.*Q)<<" mbarn/GeV^2\n";
      step++;
    }
    Q = 0.;
    differential.SetQ(Q);
    value = 
      ATOOLS::sqr(integrator.Integrate(m_Bmin,m_Bmax,m_accu,1.))/(4.*M_PI) *
      ATOOLS::rpa->Picobarn()/1.e9;
    was<<" "<<(Q*Q)<<"   "<<value<<std::endl;
    was.close();
    if (onscreen) 
      msg_Out()<<" "<<(Q*Q)<<"   "<<value<<" mbarn/GeV^2\n"
	     <<"---------------------------------------------\n";
  }
  else{
    std::string Estring(ATOOLS::ToString(2.*m_Qmax));
    std::vector<double> tvals;
    std::string infile(std::string("tvals_dsigma_el_dt_"+Estring+".dat"));
    std::ifstream input;
    input.open(infile.c_str());
    std::string test;
    while (!input.eof()) {
      input>>test;
      tvals.push_back(std::atof(test.c_str()));
    }
    input.close();
    
    std::ofstream was;
    std::string filename(dirname+std::string("/Dsigma_el_by_dt_tuning_"+Estring+".dat"));
    was.open(filename.c_str());
    was<<"# BEGIN HISTOGRAM tmp"<<std::endl;
    double value(1.), Q(m_Qmax),Qlow,Qhigh,vallow,valhigh,a,b;
    unsigned int ilow,ihigh;
//     msg_Out()<<"Calculating differential elastic cross sections for tuning."<<std::endl;
    for (int i=0; i<tvals.size(); i++) {
//       msg_Out()<<"calculating for t = "<<tvals[i]<<" GeV^2"<<std::endl;
      Q=sqrt(tvals[i]);
      ilow=int(m_logdelta*log(m_Qmax/Q))+1;
      if(ilow>m_logQsteps) ilow=m_diffgrid.size();
      ihigh=int(m_logdelta*log(m_Qmax/Q));
      Qlow=(ilow==m_diffgrid.size()?0.:m_Qmax*exp(-double(ilow)/m_logdelta));
      Qhigh=m_Qmax*exp(-double(ihigh)/m_logdelta);
      if(ilow>m_logQsteps) ilow=m_logQsteps;
      vallow=m_diffgrid[ilow];
      valhigh=m_diffgrid[ihigh];
      a=(valhigh-vallow)/(Qhigh-Qlow);
      b=vallow-a*Qlow;
      value=a*Q+b;
      was<<tvals[i]<<"   "<<tvals[i]<<"   "<<value/(2.*Q*1.e9)<<"   0.0\n";
    }
    was<<"# END HISTOGRAM"<<std::endl;
    was.close();
  }
}


double Sigma_Elastic::GetValue(const double & B) { 
  return ATOOLS::sqr(p_eikonal->Prefactor()*(1.-exp(-(*p_eikonal)(B)/2.))); 
}

double Sigma_Elastic::GetCombinedValue(const double & B) { 
  double value(0.);
  for (std::list<Omega_ik *>::iterator eikonal=p_eikonals->begin();
       eikonal!=p_eikonals->end(); eikonal++) {
    value += (*eikonal)->Prefactor()*(1.-exp(-(**eikonal)(B)/2.)); 
  }
  return ATOOLS::sqr(value);
}

double Sigma_Elastic::GetDiffArgument(const double & B) { 
  double value(0.);
  for (std::list<Omega_ik *>::iterator eikonal=p_eikonals->begin();
       eikonal!=p_eikonals->end(); eikonal++) {
    value += (*eikonal)->Prefactor()*(1.-exp(-(**eikonal)(B)/2.)); 
  }
  return value;
}

double Sigma_Elastic::PT() const {
  double random(ran->Get());
  unsigned int i(0);
  while (random-m_intgrid[i]>=0) i++;

  double Q1(m_Qmax*exp(-double(i-1)/m_logdelta));
  double Q2(i==m_intgrid.size()-1?0.:m_Qmax*exp(-double(i)/m_logdelta));
  return ((Q2*(m_intgrid[i-1]-random)+Q1*(random-m_intgrid[i]))/
	  (m_intgrid[i-1]-m_intgrid[i]));
}

void Sigma_Elastic::TestElasticCrossSection(){
  const double EulerGamma= 0.577215664901532860606512090082 ;
  double m_a,m_c,m_alpha,m_res,m_ei,m_ei2;
  double m_Delta,m_prefactor,m_Lambda2,m_beta0,m_kappa;
  ExpInt m_expint;
  m_Delta     = (*p_eikonals).front()->Delta();
  m_prefactor = (*p_eikonals).front()->Prefactor();
  m_kappa     = (*p_eikonals).front()->Kappa_i();
  m_Lambda2   = (*p_eikonals).front()->Lambda2();
  m_beta0     = (*p_eikonals).front()->FF1()->Beta0();
  m_a         = m_Lambda2/(8.*(1.+m_kappa));
  m_c         = ATOOLS::sqr(m_beta0)*m_Lambda2*(1.+m_kappa)*
    exp(2.*m_Delta*m_Y)/(8.*M_PI);
  m_alpha     = 2.*M_PI*m_prefactor;
  m_ei        = m_expint.GetExpInt(-m_c);
  m_ei2       = m_expint.GetExpInt(-m_c/2.);
  m_res       = m_alpha*(EulerGamma+m_ei-m_ei2+log(m_c/4.))/(2.*m_a);
  msg_Out() << "In " << METHOD << " sigma_elas = "<< m_res <<" 1/GeV^2 = "
	    <<m_res*rpa->Picobarn()/1.e9<<" mb ."<<std::endl;
}


