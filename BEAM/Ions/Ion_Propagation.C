#include "BEAM/Ions/Ion_Propagation.H"
#include "BEAM/Ions/Ion_Parameters.H"
#include "BEAM/Ions/Ion_Base.H"
#include "ATOOLS/Math/Special_Functions.H"

using namespace BEAM;
using namespace ATOOLS;
using namespace std;

Ion_Propagation::Ion_Propagation(Ion_Base * ion) :
  p_ion(ion),
  m_A(p_ion->A()), m_Z(ion->Z()), m_R(p_ion->R()),
  m_V(4./3.*M_PI*pow(m_R,3.)),
  m_potentials(Ion_Potentials(m_A)),
  m_analyse(true)
{
  msg_Out()<<METHOD<<"\n"; ran->Get();
  InitParameters();
  InitLists();
}

void Ion_Propagation::InitParameters() {
  m_sigma    = ionpars->Get("Sigma_Psi");
  m_sigma2   = sqr(m_sigma);
  m_Nsteps   = ionpars->GetInt("N_steps");
  m_method   = ionpars->GetInt("Evolver");
  m_deltat   = ionpars->Get("delta_t");
}

void Ion_Propagation::InitLists() {
  p_constituents = p_ion->GetConstituents();
  size_t levels = (m_method==0 ? 2 : 5);
  m_mom.resize(levels);
  m_pos.resize(levels);
  m_deltarvec.resize(levels);
  m_deltar.resize(levels);
  m_deltar2norm.resize(levels);
  for (size_t l=0;l<levels;l++) {
    m_mom[l].resize(m_A);
    m_pos[l].resize(m_A);
    m_deltarvec[l].resize(m_A);
    m_deltar[l].resize(m_A);
    m_deltar2norm[l].resize(m_A);
    for (size_t i=0;i<m_A;i++) {
      m_deltarvec[l][i].resize(m_A);
      m_deltar[l][i].resize(m_A);
      m_deltar2norm[l][i].resize(m_A);
    }
  }
  m_potentials.SetDistances(&m_deltarvec,&m_deltar,&m_deltar2norm);
  
  if (m_analyse) InitAnalysis();
}


void Ion_Propagation::Initialise() {
  HarvestPositionsAndMomenta();
  UpdateDistances();
  m_potentials.UpdateDensities();
  m_potentials.UpdatePotentials();
  CalculateRadii();
}

void Ion_Propagation::HarvestPositionsAndMomenta() {
  for (size_t i=0;i<m_A;i++) {
    m_mom[0][i] = (*p_constituents)[i]->m_p;
    m_pos[0][i] = (*p_constituents)[i]->m_r;
  }
}

void Ion_Propagation::UpdateDistances(size_t level) {
  double pref = 1./(4.*m_sigma2);
  for (size_t i=0;i<m_A-1;i++) {
    for (size_t j=i+1;j<m_A;j++) {
      m_deltarvec[level][i][j]   =  m_pos[level][i]-m_pos[level][j];
      m_deltarvec[level][j][i]   = -m_deltarvec[level][i][j];
      m_deltar2norm[level][i][j] = m_deltar2norm[level][j][i] =
	pref*m_deltarvec[level][i][j].Sqr();
      m_deltar[level][i][j]      = m_deltar[level][j][i] =
	m_deltarvec[level][i][j].Abs();
    }
  }
}

const double Ion_Propagation::Pmax(size_t i) const {
  return sqrt(2.*(*p_constituents)[i]->m_flav.Mass()*dabs(m_potentials.U(i))); 
  return ATOOLS::rpa->hBar_c()*(3.*M_PI*M_PI*m_tilderho[i]);
}


void Ion_Propagation::Evolve(int Nsteps,double dt) {
  if (Nsteps<0) Nsteps = m_Nsteps;
  if (dt<0.)    dt     = m_deltat;
  msg_Out()<<METHOD<<" for "<<Nsteps<<" time steps of length "<<dt<<" fm/c.\n";
  for (int i=0;i<Nsteps;i++) {
    HarvestPositionsAndMomenta();
    UpdateDistances(0);
    if (m_method==1) RungeKuttaStep(dt);
    else             NewtonEulerStep(dt);
    CalculateRadii(i);
  }
  if (m_analyse) FinishAnalysis();
}

void Ion_Propagation::RungeKuttaStep(const double & dt) {
  double inct=0.;
  for (size_t l=1;l<5;l++) {
    if (l==1 || l==2) inct = dt/2.;
    if (l==3 || l==4) inct = dt;
    for (size_t i=0;i<m_A;i++) {
      m_mom[l][i] = m_mom[0][i]-m_potentials.Nabla(l-1,i)*inct;
      Vec3D dr    = m_mom[l][i]/(*p_constituents)[i]->m_flav.HadMass();
      m_pos[l][i] = m_pos[0][i]+dr*inct;
    }
    UpdateDistances(l);
  }
  for (size_t i=0;i<m_A;i++) {
    (*p_constituents)[i]->m_r =
      Vec4D(0.,1./6.*(m_pos[1][i]+2.*m_pos[2][i]+2.*m_pos[3][i]+m_pos[4][i])); 
    (*p_constituents)[i]->m_p =
      Vec4D(0.,1./6.*(m_mom[1][i]+2.*m_mom[2][i]+2.*m_mom[3][i]+m_pos[4][i])); 
  }
}

void Ion_Propagation::NewtonEulerStep(const double & dt) {
  for (int i=0;i<m_A;i++) {
    m_mom[1][i] = m_mom[0][i]-m_potentials.Nabla(0,i)*dt;
    m_pos[1][i] = m_pos[0][i]+m_mom[1][i]/(*p_constituents)[i]->m_flav.Mass()*dt;
    (*p_constituents)[i]->m_p = Vec4D(0.,m_mom[1][i]);
    (*p_constituents)[i]->m_r = Vec4D(0.,m_pos[1][i]);
  }
}

void Ion_Propagation::CalculateRadii(int step) {
  double RA2 = 0., RZ2 = 0., rr;
  Vec4D R    = Vec4D(0.,0.,0.,0.);
  for (size_t i=0;i<m_A;i++) {
    Vec4D r = (*p_constituents)[i]->m_r;
    R      += r*(*p_constituents)[i]->m_flav.Mass();
    RA2    += rr = Vec3D(r).Sqr();
    if ((*p_constituents)[i]->m_flav.Kfcode()==kf_p_plus)
      RZ2  += rr;
  }
  if (step!=-1 && m_analyse) {
    m_histos["R"]->Insert(step,sqrt(-R.Abs2()));
    m_histos["R_A"]->Insert(step,sqrt(RA2/double(m_A)));
    m_histos["R_Z"]->Insert(step,sqrt(RZ2/double(m_Z)));
    m_histos["R_A_norm"]->Insert(step,sqrt(RA2/double(m_A))/m_R);
    m_histos["R_Z_norm"]->Insert(step,sqrt(RZ2/double(m_Z))/m_R);
  }
  msg_Out()<<"Check[t = "<<(step*m_deltat)<<"]: R = "<<R<<" (nom: "<<m_R<<"), "
	   <<"r_A/R = "<<(sqrt(RA2/double(m_A))/m_R)<<", "
	   <<"r_Z/R = "<<(sqrt(RZ2/double(m_Z))/m_R)<<".\n";
  if (sqrt(RA2/double(m_A))/m_R>2.) {
    for (size_t i=0;i<m_A;i++) {
      msg_Out()<<"- "
	       <<Vec3D((*p_constituents)[i]->m_r).Abs()<<" / "
	       <<Vec3D((*p_constituents)[i]->m_p).Abs()<<" from "
	       <<(*p_constituents)[i]->m_r<<" / "
	       <<(*p_constituents)[i]->m_p<<".\n";
    }
    FinishAnalysis();
    exit(1);
  }
}
  

void Ion_Propagation::Step(double dt) {
}

void Ion_Propagation::InitAnalysis() {
  m_histos["R"]        = new Histogram(0,0,m_Nsteps*m_deltat,m_Nsteps);
  m_histos["R_A"]      = new Histogram(0,0,m_Nsteps*m_deltat,m_Nsteps);
  m_histos["R_Z"]      = new Histogram(0,0,m_Nsteps*m_deltat,m_Nsteps);
  m_histos["R_A_norm"] = new Histogram(0,0,m_Nsteps*m_deltat,m_Nsteps);
  m_histos["R_Z_norm"] = new Histogram(0,0,m_Nsteps*m_deltat,m_Nsteps);
}

void Ion_Propagation::FinishAnalysis() {
  Histogram * histo;
  string name;
  for (map<string,Histogram*>::iterator hit=m_histos.begin();
       hit!=m_histos.end();hit++) {
      histo = hit->second;
      name  = string("Ion_Analysis/")+hit->first+string(".dat");
      //histo->Finalize();
      histo->Output(name);
      delete histo;
  }
  m_histos.clear();
}
