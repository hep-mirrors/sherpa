#include "PHOTONS++/PhaseSpace/Avarage_Photon_Number.H"
#include "MODEL/Main/Model_Base.H"

using namespace PHOTONS;
using namespace ATOOLS;
using namespace std;

Avarage_Photon_Number::Avarage_Photon_Number(Particle_Vector dip, double wmax, double wmin) {
  m_nbar     = 0;
  m_omegaMax = wmax;
  m_omegaMin = wmin;
  m_dipole   = dip;
#ifdef PHOTONS_DEBUG
  Vec4D sum(0.,0.,0.,0.);
  for (unsigned int i=0; i<m_dipole.size(); i++) {
    msg_Info()<<*m_dipole[i]<<endl;
    sum+=m_dipole[i]->Momentum();
  }
  msg_Info()<<sum<<endl;
#endif

  CalculateAvaragePhotonNumber();

#ifdef PHOTONS_DEBUG
  msg_Events()<<endl<<"nbar: "<<m_nbar<<" = ";
  for (unsigned int i=0; i<m_nbars.size(); i++) {
    msg_Events()<<m_nbars.at(i)<<" + ";
  }
  msg_Events()<<endl;
#endif
}

Avarage_Photon_Number::~Avarage_Photon_Number() {
}

void Avarage_Photon_Number::CalculateAvaragePhotonNumber() {
  double alphaQED = MODEL::s_model->ScalarConstant("alpha_QED(0)");
  double sum      = 0;
  for(unsigned int j=0; j<m_dipole.size(); j++) {
    for(unsigned int i=0; i<j; i++) {
      double Zi       = m_dipole.at(i)->Flav().Charge();
      double Zj       = m_dipole.at(j)->Flav().Charge();
      double titj     = 0;
      if ((m_dipole.at(i)->ProductionBlob() == m_dipole.at(j)->ProductionBlob()) && (m_dipole.at(i)->ProductionBlob() != NULL))       titj = +1;
      else if ((m_dipole.at(i)->DecayBlob() == m_dipole.at(j)->ProductionBlob()) && (m_dipole.at(i)->DecayBlob() != NULL))            titj = -1;
      else if ((m_dipole.at(i)->ProductionBlob() == m_dipole.at(j)->DecayBlob()) && (m_dipole.at(i)->ProductionBlob() != NULL))       titj = -1;
      else if ((m_dipole.at(i)->DecayBlob() == m_dipole.at(j)->DecayBlob()) && (m_dipole.at(i)->DecayBlob() != NULL))            titj = +1;
      double betai    = CalculateBeta(m_dipole.at(i)->Momentum());
      double betaj    = CalculateBeta(m_dipole.at(j)->Momentum());
      Vec3D  pi       = m_dipole.at(i)->Momentum();
      Vec3D  pj       = m_dipole.at(j)->Momentum();
      double alpha    = 0.;
      if (!IsZero(pi.Abs()) && !IsZero(pj.Abs())) {
        double arg = pi*pj/(pi.Abs()*pj.Abs());
        if (IsEqual(arg,-1.))     alpha = 0.;
        else if (IsEqual(arg,1.)) alpha = M_PI/2.;
        else                      alpha = (M_PI - acos(arg))/2.;
      }
      double ai       = betai*sin(alpha);
      double aj       = betaj*sin(alpha);
      double bi       = betai*cos(alpha);
      double bj       = betaj*cos(alpha);
#ifdef PHOTONS_DEBUG
      msg_Info()<<"DEBUG"<<endl;
      msg_Info()<<"i,j: "<<i<<","<<j<<endl;
      msg_Info()<<"betai: "<<betai<<" ,  betaj: "<<betaj<<endl;
      msg_Info()<<"pi:    "<<pi<<" ,  pj:    "<<pj<<endl;
      msg_Info()<<"pi*pj: "<<pi*pj<<" ,  |pi|: "<<pi.Abs()<<" ,  |pj|: "<<pj.Abs()<<endl;
      msg_Info()<<"arg:   "<<pi*pj/(pi.Abs()*pj.Abs())<<endl;
      msg_Info()<<acos(-1.)<<endl;
      msg_Info()<<"alpha: "<<alpha<<endl;
      msg_Info()<<"ai:    "<<ai<<" ,  aj:    "<<aj<<endl;
      msg_Info()<<"bi:    "<<bi<<" ,  bj:    "<<bj<<endl;
#endif
      double dipoleij = alphaQED/M_PI*Zi*Zj*titj*log(m_omegaMax/m_omegaMin)
                          *(2 - (1-ai*aj+bi*bj)*InterferenceTerm(ai,aj,bi,bj));
#ifdef PHOTONS_DEBUG
      msg_Info()<<"ana: "<<dipoleij<<endl;
      msg_Info()<<"int: "<<(1-ai*aj+bi*bj)*InterferenceTerm(ai,aj,bi,bj)<<endl;
#endif
      m_nbars.push_back(dipoleij);
      sum = sum + dipoleij;
    }
  }
  m_nbar = sum;
}

double Avarage_Photon_Number::CalculateBeta(Vec4D p) {
  return Vec3D(p).Abs()/p[0];
}

double Avarage_Photon_Number::InterferenceTerm(double ai, double aj, double bi, double bj) {
  if (IsZero(ai) && IsZero(aj) && IsZero(bi) && !IsZero(bj))
    return 1./bj*log((1.+bj)/(1.-bj));
  if (IsZero(ai) && IsZero(aj) && !IsZero(bi) && IsZero(bj))
    return 1./bi*log((1.+bi)/(1.-bi));
  double A  = bi-bj;
  double B  = 2*bi*bj;
  double Ci = 1-ai*ai;
  double Cj = 1-aj*aj;
  double Di = -2*bi;
  double Dj = 2*bj;
  double Ei = ai*ai+bi*bi;
  double Ej = aj*aj+bj*bj;
  double rooti = sqrt(B*B*Ci-A*B*Di+A*A*Ei);
  double rootj = sqrt(B*B*Cj-A*B*Dj+A*A*Ej);
#ifdef PHOTONS_DEBUG
      msg_Info()<<"A:  "<<A<<endl;
      msg_Info()<<"B:  "<<B<<endl;
      msg_Info()<<"Ci: "<<Ci<<" , Cj: "<<Cj<<endl;
      msg_Info()<<"Di: "<<Di<<" , Dj: "<<Dj<<endl;
      msg_Info()<<"Ei: "<<Ei<<" , Ej: "<<Ej<<endl;
      msg_Info()<<"rooti: "<<rooti<<" , rootj: "<<rootj<<endl;
#endif
  // avoid divergences in A+-B when alpha = 0 and beta1=1 and beta2=0.3333 
  if (ai/bi < 1E-4) // ai/bi = tan(alpha)
    return 1./(bi+bj)*log(((1+bi)*(1+bj))/((1-bi)*(1-bj)));
  else
    return bi*bj*(1./(bj*rooti)
                    *log(abs((A+B)/(A-B))
                         *(2*sqrt(Ci-Di+Ei)+(B*(2*Ci-Di)-A*(Di-2*Ei))/rooti)
                          /(2*sqrt(Ci+Di+Ei)+(B*(2*Ci+Di)-A*(Di+2*Ei))/rooti))
                 -1./(bi*rootj)
                    *log(abs((A+B)/(A-B))
                         *(2*sqrt(Cj-Dj+Ej)+(B*(2*Cj-Dj)-A*(Dj-2*Ej))/rootj)
                          /(2*sqrt(Cj+Dj+Ej)+(B*(2*Cj+Dj)-A*(Dj+2*Ej))/rootj)));
}

