#include "SHRiMPS/Cross_Sections/Sigma_SD.H"
#include "SHRiMPS/Tools/Special_Functions.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace SHRIMPS;
using namespace ATOOLS;

double Sigma_SD::SD_Term::operator()(double B) {
  return B * 2.*M_PI*SF.Jn(0,B*m_Q) * (1.-exp(-(*p_eikonal)(B)/2.));
}

Sigma_SD::Sigma_SD() :
  m_tmin(0.), m_tmax(1.), m_steps(1), m_delta(1.) {
  for (size_t i=0;i<2;i++) m_summed[i] = 0.;
}

double Sigma_SD::GetValue(const double & B)         { return 0.; }
double Sigma_SD::GetCombinedValue(const double & B) { return 0.; }

void Sigma_SD::FillGrids(Sigma_Elastic * sigma_el) {
  m_tmin  = sigma_el->Tmin();
  m_tmax  = sigma_el->Tmax();
  m_steps = sigma_el->Steps();
  m_delta = (m_tmax-m_tmin)/double(m_steps);
  msg_Out()<<METHOD<<" for ["<<m_tmin<<", "<<m_tmax<<"] in "<<m_steps<<" steps of "
       <<"size = "<<m_delta<<"\n";
  m_tgrids.resize(p_eikonals->size());
  for (size_t i=0;i<p_eikonals->size();i++) m_tgrids[i].resize(p_eikonals->size());

  FillTGrids();
  for (size_t diff=0;diff<3;diff++) {
    CombineTGrids(diff);
    CreateIntGrids(diff,sigma_el);
  }
}

void Sigma_SD::FillTGrids() {
  SD_Term term;
  Gauss_Integrator integrator(&term);
  double t,value;
  for (size_t k=0;k<m_steps;k++) {
    t = m_tmin + m_delta*k;
    term.SetQ(sqrt(t));
    for (size_t i=0;i<p_eikonals->size();i++) {
      for (size_t j=0;j<(*p_eikonals)[i].size();j++) {
	term.SetEikonal((*p_eikonals)[i][j]);
	value = integrator.Integrate(0.,MBpars.GetEikonalParameters().bmax,
				     MBpars.GetEikonalParameters().accu,1.);
    if (dabs(value<0.)) value = 0.;
	m_tgrids[i][j].push_back(value);
      }
    }
  }
}

void Sigma_SD::CombineTGrids(const size_t diff) {
  double pref, value, t;
  for (size_t q=0;q<m_steps;q++) {
    t     = m_tmin + m_delta*q;
    value = 0.;
    for (size_t i=0;i<p_eikonals->size();i++) {
      for (size_t j=0;j<(*p_eikonals)[i].size();j++) {
	for (size_t k=0;k<(*p_eikonals)[i].size();k++) {
	  if (diff==0) {
	    pref  = ((*p_eikonals)[i][j]->Prefactor()*(*p_eikonals)[i][k]->Prefactor()/
		     sqrt((*p_eikonals)[i][i]->Prefactor())/
		     (4.*M_PI));
	    value += pref * m_tgrids[i][j][q] * m_tgrids[i][k][q] * rpa->Picobarn();
	  }
	  else if (diff==1) {
	    pref  = ((*p_eikonals)[j][i]->Prefactor()*(*p_eikonals)[k][i]->Prefactor()/
		     sqrt((*p_eikonals)[i][i]->Prefactor())/
		     (4.*M_PI));
	    value += pref * m_tgrids[j][i][q] * m_tgrids[k][i][q] * rpa->Picobarn();
	  }
	  else if (diff==2) {
	    for (size_t l=0;l<p_eikonals->size();l++) {
	      pref  = ((*p_eikonals)[i][j]->Prefactor()*(*p_eikonals)[l][k]->Prefactor()/
		       (4.*M_PI));
	      value += pref * m_tgrids[i][j][q] * m_tgrids[l][k][q] * rpa->Picobarn();
	    }
	  }
	}
      }
    }
    m_diffgrids[diff].push_back(value);
  }
}

void Sigma_SD::CreateIntGrids(const size_t diff,Sigma_Elastic * sigma_el) {
  m_summed[diff] = 0.;
  m_intgrids[diff].push_back(0.);
  std::vector<double> el_grid = sigma_el->GetDiffGrid();
  for (size_t i=0;i<m_diffgrids[diff].size();i++) m_diffgrids[diff][i] -= el_grid[i];
  for (size_t i=1;i<m_diffgrids[diff].size();i++) {
    m_summed[diff] += (m_diffgrids[diff][i]+m_diffgrids[diff][i-1])/2. * m_delta;
    m_intgrids[diff].push_back(m_summed[diff]);
  }
  for (size_t i=0;i<m_intgrids[diff].size();i++) m_intgrids[diff][i] /= m_summed[diff];
}

double Sigma_SD::SelectT(const size_t & mode) const {
  double random(ran->Get());
  unsigned int i(0);
  while (random-m_intgrids[mode][i]>=0) i++;
  return m_tmin+(i-1)*m_delta +
    m_delta *(random-m_intgrids[mode][i-1])/(m_intgrids[mode][i]-m_intgrids[mode][i-1]);
}

