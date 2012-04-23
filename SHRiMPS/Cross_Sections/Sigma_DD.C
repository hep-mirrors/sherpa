#include "SHRiMPS/Cross_Sections/Sigma_DD.H"
#include "SHRiMPS/Tools/Special_Functions.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace SHRIMPS;
using namespace ATOOLS;

double Sigma_DD::dSigma_dt_Kernel::operator()(double B) {
//   msg_Out()<<METHOD<<"("<<B<<"), eikonal = "<<(*p_eikonal)(B)<<".\n";
  return 2.*M_PI*B*SF.Jn(0,B*m_Q)*(1.-exp(-(*p_eikonal)(B)/2.));
}

Sigma_DD::Sigma_DD(Sigma_Elastic * sigma_el,Sigma_SD * sigma_sd,const int & test) :
  Sigma_Base(sigma_el->Eikonals()),
  p_sigma_el(sigma_el),p_sigma_sd(sigma_sd),
  m_Bmin(p_sigma_el->Bmin()), m_Bmax(p_sigma_el->Bmax()),
  m_Qmax(p_sigma_el->Qmax()), 
  m_logQsteps(p_sigma_el->Steps()), m_logdelta(p_sigma_el->Delta()), 
  m_test(test) 
{ 
  FillGrids();
}

void Sigma_DD::FillGrids() {
  if (m_test==10) PrintDifferentialElasticAndDiffXsec(true);
  msg_Tracking()<<"In "<<METHOD<<": Integrate from "
		<<m_Bmin<<" to "<<m_Bmax<<"."<<std::endl
		<<"   Maximal sqrt{|t|} = "<<m_Qmax<<"."<<std::endl;
  m_intgrid_DD.clear();
  m_intgrid_DD.push_back(0.);
  m_diffgrid_DD.clear();

  int noFF(0);
  for (std::list<Omega_ik *>::iterator eikiter=p_eikonals->begin();
       eikiter!=p_eikonals->end();eikiter++) {
    if ((*eikiter)->FF1()->Number()>noFF) noFF = (*eikiter)->FF1()->Number();
  }
  noFF++;
  std::vector<std::vector<Omega_ik *> > eikonals;
  std::vector<double> prefs;
  eikonals.resize(noFF);
  prefs.resize(noFF);
  for (int i=0;i<noFF;i++) {
    eikonals[i].resize(noFF);
    for (int j=0;j<noFF;j++) {
      for (std::list<Omega_ik *>::iterator eikiter=p_eikonals->begin();
	   eikiter!=p_eikonals->end();eikiter++) {
	if ((*eikiter)->FF1()->Number()==i && (*eikiter)->FF2()->Number()==j)
	  eikonals[i][j] = (*eikiter);
	if ((*eikiter)->FF1()->Number()==i) 
	  prefs[i] = (*eikiter)->FF1()->Prefactor();
      }
    }
  }

  dSigma_dt_Kernel differential;
  ATOOLS::Gauss_Integrator integrator(&differential);
  double pref(0.), cumul(0.), value(0.);
  double prefQ(0.), Q(0.), sigmael(p_sigma_el->Sigma()),sigmasd(p_sigma_sd->Sigma());
  std::vector<std::vector<double> > values;
  values.resize(noFF);
  for (int i=0;i<noFF;i++) values[i].resize(noFF);

  size_t step(0);
  while (step<m_logQsteps || 
	 (step>1 && (value-pref)/(value+pref)>1.e-12)) {
    Q     = m_Qmax*exp(-double(step)/m_logdelta);
    differential.SetQ(Q);
    for (int i=0;i<noFF;i++) {
      for (int j=0;j<noFF;j++) {
// 	msg_Out()<<METHOD<<"("<<i<<", "<<j<<") : "<<eikonals[i][j]<<".\n";
	differential.SetEikonal(eikonals[i][j]);
	values[i][j] = integrator.Integrate(m_Bmin,m_Bmax,m_accu,1.);
      }
    }
    value = 0.;
    for (int i=0;i<noFF;i++) {
      for (int j=0;j<noFF;j++) {
	value += 
	values[i][j]*values[i][j]*sqr(prefs[i]*prefs[j]);
      }
    }
    value *= ATOOLS::rpa->Picobarn()/(4.*M_PI);
    m_diffgrid_DD.push_back(value);

    msg_Tracking()<<"   Q = "<<Q<<" --> dsigma/dt = "
		  <<(value/1.e9)<<" mbarn."<<std::endl;

    if (step>0) {
      cumul += (value+pref)/2. * (prefQ-Q)*(prefQ+Q);
      m_intgrid_DD.push_back(cumul);
    }
    prefQ = Q;
    pref = value;
    step++;
  }
  Q = 0.;

  differential.SetQ(Q);
  for (int i=0;i<noFF;i++) {
    for (int j=0;j<noFF;j++) {
      differential.SetEikonal(eikonals[i][j]);
      values[i][j] = integrator.Integrate(m_Bmin,m_Bmax,m_accu,1.);
    }
  }
  value = 0.;
  for (int i=0;i<noFF;i++) {
    for (int j=0;j<noFF;j++) {
      value += values[i][j]*values[i][j]*
	  sqr(prefs[i])*sqr(prefs[j])/(4*M_PI);
    }
  }
  value *= ATOOLS::rpa->Picobarn()/(4.*M_PI);
  m_diffgrid_DD.push_back(value);

  cumul += (value+pref)/2. * (prefQ-Q)*(prefQ+Q);
  m_intgrid_DD.push_back(cumul);

  m_sigma_DD = cumul-sigmasd-sigmael;
  m_sigma     = m_sigma_DD;

  msg_Tracking()<<"Sigma_{el} = "<<(sigmael/1.e9)<<", "<<"Sigma_{SD} = "<<(sigmasd/1.e9)
		<<", Sigma_{DD} = "<<(m_sigma_DD/1.e9)<<", sum = "<<(cumul/1.e9)<<" mbarn."<<std::endl;

  const std::vector<double> & gridel = *p_sigma_el->Grid();
  const std::vector<double> & gridsd = *p_sigma_sd->Grid1();
  for (size_t i=0;i<m_intgrid_DD.size();i++) {
    m_intgrid_DD[i] -= gridel[i]*sigmael + gridsd[i]*sigmasd;
    m_intgrid_DD[i] /= (cumul-sigmasd-sigmael);
    msg_Tracking()<<i<<"  Q = "<<m_Qmax*exp(-double(i)/m_logdelta)<<" --> "
		  <<m_intgrid_DD[i]<<"."<<std::endl;
  }
}


double Sigma_DD::GetValue(const double & B) { 
  return 0.;
}

double Sigma_DD::GetCombinedValue(const double & B) { 
  double sdvalue(0.),elvalue(0.),ddvalue(0.),fac(1.);
  for (std::list<Omega_ik *>::iterator eikonal1=p_eikonals->begin();
       eikonal1!=p_eikonals->end(); eikonal1++) {
    for (std::list<Omega_ik *>::iterator eikonal2=p_eikonals->begin();
	 eikonal2!=p_eikonals->end(); eikonal2++) {
      fac = 1.;
      if ((*eikonal1)->GetSingleTerm(0)->FF1()->Number()==
	  (*eikonal2)->GetSingleTerm(0)->FF1()->Number()
	  && (*eikonal1)->GetSingleTerm(0)->FF2()->Number()==
	  (*eikonal2)->GetSingleTerm(0)->FF2()->Number()) {
	fac += 1./(sqr((*eikonal1)->GetSingleTerm(0)->FF1()->Prefactor())*
	   sqr((*eikonal1)->GetSingleTerm(0)->FF2()->Prefactor()));
      }
      if ((*eikonal1)->GetSingleTerm(0)->FF1()->Number()==
	  (*eikonal2)->GetSingleTerm(0)->FF1()->Number()) {
	fac -= 1./sqr((*eikonal1)->GetSingleTerm(0)->FF1()->Prefactor());
      }
      if ((*eikonal1)->GetSingleTerm(0)->FF2()->Number()==
	  (*eikonal2)->GetSingleTerm(0)->FF2()->Number()) {
	fac -= 1./sqr((*eikonal1)->GetSingleTerm(0)->FF2()->Prefactor());
      }
      ddvalue += fac * (*eikonal1)->Prefactor()*(1.-exp(-(**eikonal1)(B)/2.))
	            * (*eikonal2)->Prefactor()*(1.-exp(-(**eikonal2)(B)/2.));
    }  
  }      
  return ddvalue;
}


void Sigma_DD::PrintDifferentialElasticAndDiffXsec(const bool & onscreen,std::string dirname) {
  std::ofstream was;
  std::string Estring(ATOOLS::ToString(2.*m_Qmax));
  std::string filename(dirname+std::string("/Dsigma_DD_by_dt_"+Estring+".dat"));
  was.open(filename.c_str());

  double Q(m_Qmax);
  if (onscreen) msg_Out()<<"---------------------------------------------\n";

  for (size_t i=0;i<m_diffgrid_DD.size();i++) {
    Q     = m_Qmax*exp(-double(i)/m_logdelta);
    was<<" "<<(Q*Q)<<"   "<<m_diffgrid_DD[i]/1.e9<<std::endl;
    if (onscreen) msg_Out()<<" "<<(Q*Q)<<"   "<<m_diffgrid_DD[i]/1.9<<" mbarn/GeV^2\n";
  }
  was.close();
  if (onscreen) msg_Out()<<"---------------------------------------------\n";
  }


double Sigma_DD::PT() {  
  const std::vector<double> & grid = m_intgrid_DD;

  double random(ran->Get());
  size_t i(0);
  while (random-grid[i]>=0) i++;

  double Q1(sqr(m_Qmax*exp(-double(i-1)/m_logdelta)));
  double Q2(sqr(i==grid.size()-1?0.:m_Qmax*exp(-double(i)/m_logdelta)));
  return ((Q2*(grid[i-1]-random)+Q1*(random-grid[i]))/(grid[i-1]-grid[i]));
}
