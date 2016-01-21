#include "SHRiMPS/Cross_Sections/Sigma_SD.H"
#include "SHRiMPS/Tools/Special_Functions.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace SHRIMPS;
using namespace ATOOLS;

double Sigma_SD::dSigma_dt_Kernel::operator()(double B) {
//   msg_Out()<<METHOD<<"("<<B<<"), eikonal = "<<(*p_eikonal)(B)<<".\n";
  return 2.*M_PI*B*SF.Jn(0,B*m_Q)*(1.-exp(-(*p_eikonal)(B)/2.));
}

Sigma_SD::Sigma_SD(Sigma_Elastic * sigma_el) :
  Sigma_Base(sigma_el->Eikonals()),
  p_sigma_el(sigma_el),
  m_Bmin(p_sigma_el->Bmin()), m_Bmax(p_sigma_el->Bmax()),
  m_Qmax(p_sigma_el->Qmax()), 
  m_logQsteps(p_sigma_el->Steps()), m_logdelta(p_sigma_el->Delta())
{ 
  FillGrids();
}

void Sigma_SD::FillGrids() {
  msg_Tracking()<<"In "<<METHOD<<": Integrate from "
		<<m_Bmin<<" to "<<m_Bmax<<"."<<std::endl
		<<"   Maximal sqrt{|t|} = "<<m_Qmax<<"."<<std::endl;
  m_intgrid_SD1.clear();
  m_intgrid_SD1.push_back(0.);
  m_intgrid_SD2.clear();
  m_intgrid_SD2.push_back(0.);
  m_diffgrid_SD1.clear();
  m_diffgrid_SD2.clear();

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
  double pref1(0.), cumul1(0.), value1(0.);
  double pref2(0.), cumul2(0.), value2(0.);
  double prefQ(0.), Q(0.), sigmael(p_sigma_el->Sigma());
  std::vector<std::vector<double> > values;
  values.resize(noFF);
  for (int i=0;i<noFF;i++) values[i].resize(noFF);

  size_t step(0);
  while (step<m_logQsteps || 
	 (step>1 && (value1-pref1)/(value1+pref1)>1.e-12 && 
	  (value2-pref2)/(value2+pref2)>1.e-12)) {
    Q     = m_Qmax*exp(-double(step)/m_logdelta);
    differential.SetQ(Q);
    for (int i=0;i<noFF;i++) {
      for (int j=0;j<noFF;j++) {
	//msg_Out()<<METHOD<<"("<<i<<", "<<j<<") : "<<eikonals[i][j]<<".\n";
	differential.SetEikonal(eikonals[i][j]);
	values[i][j] = integrator.Integrate(m_Bmin,m_Bmax,m_accu,1.);
      }
    }
    value1 = value2 = 0.;
    for (int i=0;i<noFF;i++) {
      for (int j=0;j<noFF;j++) {
	for (int k=0;k<noFF;k++) {
	  value1 += 
	    values[i][j]*values[i][k]*sqr(prefs[i]*prefs[j]*prefs[k]);
	  value2 += 
	    values[j][i]*values[k][i]*sqr(prefs[i]*prefs[j]*prefs[k]);
	}
      }
    }
    value1 *= ATOOLS::rpa->Picobarn()/(4.*M_PI);
    value2 *= ATOOLS::rpa->Picobarn()/(4.*M_PI);
    m_diffgrid_SD1.push_back(value1);
    m_diffgrid_SD2.push_back(value2);

    msg_Tracking()<<"   Q = "<<Q<<" --> dsigma/dt = "
		  <<(value1/1.e9)<<"/"<<(value2/1.e9)<<" mbarn."<<std::endl;

    if (step>0) {
      cumul1 += (value1+pref1)/2. * (prefQ-Q)*(prefQ+Q);
      cumul2 += (value2+pref2)/2. * (prefQ-Q)*(prefQ+Q);
      m_intgrid_SD1.push_back(cumul1);
      m_intgrid_SD2.push_back(cumul2);
    }
    prefQ = Q;
    pref1 = value1;
    pref2 = value2;
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
  value1 = value2 = 0.;
  for (int i=0;i<noFF;i++) {
    for (int j=0;j<noFF;j++) {
      for (int k=0;k<noFF;k++) {
	value1 += 
	  values[i][j]*values[i][k]*
	  sqr(prefs[i])*sqr(prefs[j])*sqr(prefs[k])/(4*M_PI);
	value2 += 
	  values[j][i]*values[k][i]*
	  sqr(prefs[i])*sqr(prefs[j])*sqr(prefs[k])/(4*M_PI);
      }
    }
  }
  value1 *= ATOOLS::rpa->Picobarn()/(4.*M_PI);
  value2 *= ATOOLS::rpa->Picobarn()/(4.*M_PI);
  m_diffgrid_SD1.push_back(value1);
  m_diffgrid_SD2.push_back(value2);

  cumul1 += (value1+pref1)/2. * (prefQ-Q)*(prefQ+Q);
  cumul2 += (value2+pref2)/2. * (prefQ-Q)*(prefQ+Q);
  m_intgrid_SD1.push_back(cumul1);
  m_intgrid_SD2.push_back(cumul2);

  m_sigma_SD1 = cumul1-sigmael;
  m_sigma_SD2 = cumul2-sigmael;
  m_sigma     = m_sigma_SD1+m_sigma_SD2;

  msg_Tracking()<<"Sigma_{el} = "<<(sigmael/1.e9)<<", "
		<<"Sigma_{SD} = "<<(m_sigma_SD1/1.e9)<<" + "
		<<(m_sigma_SD2/1.e9)<<", total = "<<(cumul1+cumul2)/1.e9<<" mbarn."<<std::endl;

  const std::vector<double> & grid = *p_sigma_el->Grid();
  for (size_t i=0;i<m_intgrid_SD1.size();i++) {
    m_intgrid_SD1[i] -= grid[i]*sigmael;
    m_intgrid_SD2[i] -= grid[i]*sigmael;
    m_intgrid_SD1[i] /= (cumul1-sigmael);
    m_intgrid_SD2[i] /= (cumul2-sigmael);
    msg_Debugging()<<i<<"  Q = "<<m_Qmax*exp(-double(i)/m_logdelta)<<" --> "
		  <<m_intgrid_SD1[i]<<", "<<m_intgrid_SD2[i]<<"."<<std::endl;
  }
}


double Sigma_SD::GetValue(const double & B) { 
  return 0.;
}

double Sigma_SD::GetCombinedValue(const double & B) { 
  double sdvalue(0.),elvalue(0.);
  for (std::list<Omega_ik *>::iterator eikonal1=p_eikonals->begin();
       eikonal1!=p_eikonals->end(); eikonal1++) {
    for (std::list<Omega_ik *>::iterator eikonal2=p_eikonals->begin();
	 eikonal2!=p_eikonals->end(); eikonal2++) {
      if ((*eikonal1)->GetSingleTerm(0)->FF1()->Number()==
	  (*eikonal2)->GetSingleTerm(0)->FF1()->Number()) {
	sdvalue += 
	  (*eikonal1)->Prefactor()*(1.-exp(-(**eikonal1)(B)/2.)) *
	  (*eikonal2)->Prefactor()*(1.-exp(-(**eikonal2)(B)/2.)) /
	  ATOOLS::sqr((*eikonal1)->GetSingleTerm(0)->FF1()->Prefactor());
      }
      if ((*eikonal1)->GetSingleTerm(0)->FF2()->Number()==
	  (*eikonal2)->GetSingleTerm(0)->FF2()->Number()) {
	sdvalue += 
	  (*eikonal1)->Prefactor()*(1.-exp(-(**eikonal1)(B)/2.)) *
	  (*eikonal2)->Prefactor()*(1.-exp(-(**eikonal2)(B)/2.)) /
	  ATOOLS::sqr((*eikonal1)->GetSingleTerm(0)->FF2()->Prefactor());
      }
    }
  }
  double value(0.);
  for (std::list<Omega_ik *>::iterator eikonal=p_eikonals->begin();
       eikonal!=p_eikonals->end(); eikonal++) {
    value += (*eikonal)->Prefactor()*(1.-exp(-(**eikonal)(B)/2.)); 
  }
  elvalue = 2.*ATOOLS::sqr(value);
  return sdvalue-elvalue;
}


void Sigma_SD::PrintDifferentialElasticAndSDXsec(const bool & onscreen,std::string dirname) {
  std::ofstream was;
  std::string Estring(ATOOLS::ToString(2.*m_Qmax));
  std::string filename(dirname+std::string("/Dsigma_SD_by_dt_"+Estring+".dat"));
 was.open(filename.c_str());

  double Q(m_Qmax);
  if (onscreen) msg_Out()<<"---------------------------------------------\n";

  for (size_t i=0;i<m_diffgrid_SD1.size();i++) {
    Q     = m_Qmax*exp(-double(i)/m_logdelta);
    was<<" "<<(Q*Q)<<"   "<<(m_diffgrid_SD1[i]+m_diffgrid_SD2[i])/1.e9<<std::endl;
    if (onscreen) msg_Out()<<" "<<(Q*Q)<<"   "<<(m_diffgrid_SD1[i]+m_diffgrid_SD2[i])/1.e9<<" mbarn/GeV^2\n";
  }
  was.close();
  if (onscreen) msg_Out()<<"---------------------------------------------\n";
    
}


double Sigma_SD::PT2(bool & mode) {  
  const std::vector<double> & grid = m_intgrid_SD1;
  if (m_sigma_SD1/m_sigma>ATOOLS::ran->Get()) mode = false;
  else mode = true;

  double random(ran->Get());
  size_t i(0);
  while (random-grid[i]>=0) i++;
  
  double Q1(sqr(m_Qmax*exp(-double(i-1)/m_logdelta)));
  double Q2(sqr(i==grid.size()-1?0.:m_Qmax*exp(-double(i)/m_logdelta)));
  return ((Q2*(grid[i-1]-random)+Q1*(random-grid[i]))/(grid[i-1]-grid[i]));
}
