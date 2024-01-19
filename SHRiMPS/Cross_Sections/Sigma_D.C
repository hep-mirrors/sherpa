#include "SHRiMPS/Cross_Sections/Sigma_D.H"
#include "SHRiMPS/Tools/Special_Functions.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace SHRIMPS;
using namespace ATOOLS;

double Sigma_D::D_Term::operator()(double B) {
  return B * 2.*M_PI*SF.Jn(0,B*m_Q) * (1.-exp(-(*p_eikonal)(B)/2.));
}

Sigma_D::Sigma_D() :
  m_tmin(0.), m_tmax(1.), m_steps(1), m_delta(1.) {
  for (size_t i=0;i<2;i++) m_summed[i] = 0.;
}

double Sigma_D::GetValue(const double & B) {
    return p_eikonal->Prefactor()*sqr(1.-exp(-(*p_eikonal)(B)/2.));
}

double Sigma_D::GetCombinedValue(const double & B) {
  double value(0.);
  for (size_t i=0;i<p_eikonals->size();i++) {
    for (size_t j=0;j<(*p_eikonals)[i].size();j++) {
      Omega_ik * eikonal = (*p_eikonals)[i][j];
      value += eikonal->Prefactor()*sqr(1.-exp(-(*eikonal)(B)/2.));
    }
  }
  return value;
}

double Sigma_D::GetValuePerChannel(const int k, const int l, const double &B) {
    if (k >= int(p_eikonals->size()) || l >= int(p_eikonals->size())) {
        msg_Error()<<"Error in "<<METHOD<<": requested GW state does not exist, will return zero for cross section.\n";
        return 0.;
    }
    double value(0.);
    if (k >= 0 && l >= 0) {
        Omega_ik * eikonal = (*p_eikonals)[k][l];
        value = eikonal->Prefactor()*sqr(1.-exp(-(*eikonal)(B)/2.));
        value *= sqr(p_eikonals->size());
    }
    if (k >= 0 && l < 0) {
        for (size_t j=0;j<(*p_eikonals)[k].size();j++) {
            Omega_ik * eikonal = (*p_eikonals)[k][j];
            value += eikonal->Prefactor()*sqr(1.-exp(-(*eikonal)(B)/2.));
        }
        value *= p_eikonals->size();
    }
    if (k < 0 && l >= 0) {
        for (size_t i=0;i<(*p_eikonals)[l].size();i++) {
            Omega_ik * eikonal = (*p_eikonals)[i][l];
            value += eikonal->Prefactor()*sqr(1.-exp(-(*eikonal)(B)/2.));
        }
        value *= p_eikonals->size();
    }
    if (k < 0 && l < 0) value = GetCombinedValue(B);
    return value;
}


double Sigma_D::GetCombinedValueEL(const double & B) {
    double value(0.);
    for (size_t i=0;i<p_eikonals->size();i++) {
      for (size_t j=0;j<(*p_eikonals)[i].size();j++) {
        Omega_ik * eikonal = (*p_eikonals)[i][j];
        value += eikonal->Prefactor()*(1.-exp(-(*eikonal)(B)/2.));
      }
    }
    return sqr(value);
}

double Sigma_D::GetCombinedValueSD0(const double & B) {
  double value(0.),kaverage(0.);
  for (size_t i=0;i<p_eikonals->size();i++) {
      kaverage = 0.;
      for (size_t j=0;j<(*p_eikonals)[i].size();j++) {
          Omega_ik * eikonal = (*p_eikonals)[i][j];
          kaverage += sqr(eikonal->FF2()->Prefactor())*(1.-exp(-(*eikonal)(B)/2.));
      }
      value += sqr((*p_eikonals)[i][i]->FF1()->Prefactor())*sqr(kaverage);
  }
  return value - GetCombinedValueEL(B);
}

double Sigma_D::GetCombinedValueSD1(const double & B) {
    double value(0.),kaverage(0.);
    for (size_t i=0;i<p_eikonals->size();i++) {
        kaverage = 0.;
        for (size_t j=0;j<(*p_eikonals)[i].size();j++) {
            Omega_ik * eikonal = (*p_eikonals)[i][j];
            kaverage += sqr(eikonal->FF1()->Prefactor())*(1.-exp(-(*eikonal)(B)/2.));
        }
        value += sqr((*p_eikonals)[i][i]->FF2()->Prefactor())*sqr(kaverage);
    }
    return value - GetCombinedValueEL(B);
}

double Sigma_D::GetCombinedValueDD(const double & B) {
  return GetCombinedValue(B) - GetCombinedValueSD0(B) - GetCombinedValueSD1(B) - GetCombinedValueEL(B);
}


void Sigma_D::FillGrids(Sigma_Elastic * sigma_el) {
  m_tgrids.clear();
  for (size_t i = 0; i < 3; ++i) {
    m_intgrids[i].clear();
    m_diffgrids[i].clear();
  }
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

void Sigma_D::FillTGrids() {
  D_Term term;
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
        //if (dabs(value<0.)) value = 0.;
      	m_tgrids[i][j].push_back(value);
      }
    }
  }
}

double Sigma_D::GetXSvsT(const size_t diff, double t) {
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
    if (r_dist < l_dist) value = m_diffgrids[diff][l_ind];
    else if (l_dist < r_dist) value = m_diffgrids[diff][r_ind];
    else value = (m_diffgrids[diff][l_ind] + m_diffgrids[diff][r_ind])/2.;
    return value;
  }
  else {
    msg_Error() << "Error in " << METHOD << " t value out of range";
    return 0.;
  }
}

void Sigma_D::CombineTGrids(const size_t diff) {
  double pref, value, t;
  for (size_t q=0;q<m_steps;q++) {
    t     = m_tmin + m_delta*q;
    value = 0.;
    for (size_t i=0;i<p_eikonals->size();i++) {        
        for (size_t j=0;j<(*p_eikonals)[i].size();j++) {
            if (diff == 2) {
                pref = (*p_eikonals)[i][j]->Prefactor()/(4.*M_PI);
                value += pref * sqr(m_tgrids[i][j][q]) * rpa->Picobarn();
            }
            else {
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
                }
            }
        }
    }
    m_diffgrids[diff].push_back(value);
  }
}

void Sigma_D::CreateIntGrids(const size_t diff,Sigma_Elastic * sigma_el) {
  m_summed[diff] = 0.;
  m_intgrids[diff].push_back(0.);
  std::vector<double> el_grid = sigma_el->GetDiffGrid();
  if (diff < 2) {
    for (size_t i=0;i<m_diffgrids[diff].size();i++) {
      m_diffgrids[diff][i] -= el_grid[i];
    }
  }
  else if (diff == 2) {
    for (size_t i=0;i<m_diffgrids[diff].size();i++) {
      m_diffgrids[diff][i] -= (m_diffgrids[0][i] + m_diffgrids[1][i] + el_grid[i]);
    }
  }
  for (size_t i=1;i<m_diffgrids[diff].size();i++) {
    m_summed[diff] += (m_diffgrids[diff][i]+m_diffgrids[diff][i-1])/2. * m_delta;
    m_intgrids[diff].push_back(m_summed[diff]);
  }
  for (size_t i=0;i<m_intgrids[diff].size();i++) m_intgrids[diff][i] /= m_summed[diff];
}

double Sigma_D::SelectT(const size_t & mode) const {
  double random(ran->Get());
  unsigned int i(0);
  while (random-m_intgrids[mode][i]>=0) i++;
  return m_tmin+(i-1)*m_delta +
    m_delta *(random-m_intgrids[mode][i-1])/(m_intgrids[mode][i]-m_intgrids[mode][i-1]);
}

