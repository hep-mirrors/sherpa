#include "PHOTONS++/PhaseSpace/Avarage_Photon_Number.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Particle.H"
#include "PHOTONS++/Main/Photons.H"
#include "PHOTONS++/Main/Dipole_Type.H"

using namespace PHOTONS;
using namespace ATOOLS;
using namespace std;

Avarage_Photon_Number::Avarage_Photon_Number
(const Particle_Vector& dip, const double& wmax, const double& wmin) :
  m_omegaMax(wmax), m_omegaMin(wmin), m_dipole(dip), m_nbar(0.)
{
  DEBUG_FUNC("");
  CalculateAvaragePhotonNumber();
  for (unsigned int i=0; i<m_nbars.size(); i++) {
    msg_Debugging()<<m_nbars[i]<<std::endl;
  }
  msg_Debugging()<<"-> "<<m_nbar<<endl;
  return;
  // set those which cancel to zero if individually larger than sum
  // and they cancel to within 10% of the resultant nbar
  // or they cancel to within 5% of each other
  double nbardiff(0.);
  size_t nonzeronbars(m_nbars.size());
  for (size_t i(0);i<m_nbars.size();++i) {
    for (size_t j(0);j<m_nbars.size();++j) {
      if (m_dipole[m_nbars[i].ij.i]->Flav()!=
            m_dipole[m_nbars[i].ij.j]->Flav().Bar() &&
          m_dipole[m_nbars[i].ij.i]->Flav()!=
            m_dipole[m_nbars[i].ij.j]->Flav().Bar() &&
          m_nbars[i].nbar>0.1*m_nbar && m_nbars[j].nbar<-0.1*m_nbar &&
          ((m_nbars[i].nbar+m_nbars[j].nbar)<0.1*m_nbar ||
           (m_nbars[i].nbar+m_nbars[j].nbar)<0.025*(m_nbars[i].nbar-m_nbars[j].nbar))) {
        msg_Debugging()<<"removing "
                       <<m_nbars[i]<<", "<<m_nbars[j]<<" -> "
                       <<m_nbars[i].nbar+m_nbars[j].nbar<<std::endl;
        nbardiff+=m_nbars[i].nbar+m_nbars[j].nbar;
        m_nbars[i].nbar=m_nbars[j].nbar=0.;
        nonzeronbars-=2;
      }
    }
  }
  // adjust other nbars
  if (nbardiff!=0.) {
    for (size_t i(0); i<m_nbars.size(); ++i) {
      if (m_nbars[i].nbar!=0.) m_nbars[i].nbar+=(nbardiff)/nonzeronbars;
    }
    for (unsigned int i=0; i<m_nbars.size(); i++) {
      msg_Debugging()<<m_nbars[i]<<std::endl;
    }
    msg_Debugging()<<"-> "<<m_nbar<<endl;
  }
}

Avarage_Photon_Number::~Avarage_Photon_Number() {
}

void Avarage_Photon_Number::CalculateAvaragePhotonNumber() {
  double sum      = 0.;
  for(unsigned int j=0; j<m_dipole.size(); j++) {
    for(unsigned int i=0; i<j; i++) {
      if ((j==1 &&  i==0) ||
	  (j==2 && (i==0 || i==1)) ||
	  (j==3 && (i==0 || i==2))) continue;
      double Zi       = m_dipole[i]->Flav().Charge();
      double Zj       = m_dipole[j]->Flav().Charge();
      double titj     = 0;
      if ((m_dipole[i]->ProductionBlob() == m_dipole[j]->ProductionBlob()) &&
          (m_dipole[i]->ProductionBlob() != NULL))       titj = +1.;
      else if ((m_dipole[i]->DecayBlob() == m_dipole[j]->ProductionBlob()) &&
               (m_dipole[i]->DecayBlob() != NULL))       titj = -1.;
      else if ((m_dipole[i]->ProductionBlob() == m_dipole[j]->DecayBlob()) &&
               (m_dipole[i]->ProductionBlob() != NULL))  titj = -1.;
      else if ((m_dipole[i]->DecayBlob() == m_dipole[j]->DecayBlob()) &&
               (m_dipole[i]->DecayBlob() != NULL))       titj = +1.;
      double betai    = CalculateBeta(m_dipole[i]->Momentum());
      double betaj    = CalculateBeta(m_dipole[j]->Momentum());
      Vec3D  pi       = m_dipole[i]->Momentum();
      Vec3D  pj       = m_dipole[j]->Momentum();
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
      double dipoleij = Photons::s_alpha/M_PI*Zi*Zj*titj
			  *log(m_omegaMax/m_omegaMin)
                          *(2.-(1.-ai*aj+bi*bj)*InterferenceTerm(ai,aj,bi,bj));
#ifdef PHOTONS_DEBUG
      msg_Debugging()<<"ana: "<<dipoleij<<endl;
      msg_Debugging()<<"alpha/pi*ZiZjtitj: "
                     <<Photons::s_alpha/M_PI*Zi*Zj*titj<<endl;
      msg_Debugging()<<"log(Emax/Emin): "<<log(m_omegaMax/m_omegaMin)<<endl;
      msg_Debugging()<<"int: "
                     <<(1-ai*aj+bi*bj)*InterferenceTerm(ai,aj,bi,bj)<<endl;
      msg_Debugging()<<"Emax,Emin: "<<m_omegaMax<<" ,  "<<m_omegaMin<<endl;
#endif
      m_nbars.push_back(IdPairNbar(IdPair(i,j),dipoleij));
      sum += dipoleij;
    }
  }
  m_nbar = sum;
}

double Avarage_Photon_Number::CalculateBeta(const Vec4D& p) {
  return Vec3D(p).Abs()/p[0];
}

double Avarage_Photon_Number::InterferenceTerm
(const double& ai, const double& aj, const double& bi, const double& bj) {
  if (IsZero(ai) && IsZero(aj) && IsZero(bi) && !IsZero(bj))
    return 1./bj*log((1.+bj)/(1.-bj));
  if (IsZero(ai) && IsZero(aj) && !IsZero(bi) && IsZero(bj))
    return 1./bi*log((1.+bi)/(1.-bi));
  double A  = bi-bj;
  double B  = 2.*bi*bj;
  double Ci = 1.-ai*ai;
  double Cj = 1.-aj*aj;
  double Di = -2.*bi;
  double Dj = 2.*bj;
  double Ei = ai*ai+bi*bi;
  double Ej = aj*aj+bj*bj;
  double rooti = sqrt(B*B*Ci-A*B*Di+A*A*Ei);
  double rootj = sqrt(B*B*Cj-A*B*Dj+A*A*Ej);
  // avoid divergences in A+-B when alpha = 0 and beta1=1 and beta2=0.3333 
  if (ai/bi < 1E-4) // ai/bi = tan(alpha)
    return 1./(bi+bj)*log(((1.+bi)*(1.+bj))/((1.-bi)*(1.-bj)));
  else
    return bi*bj*(1./(bj*rooti)
                  *log(abs((A+B)/(A-B))
                       *(2.*sqrt(Ci-Di+Ei)+(B*(2.*Ci-Di)-A*(Di-2.*Ei))/rooti)
                       /(2.*sqrt(Ci+Di+Ei)+(B*(2.*Ci+Di)-A*(Di+2.*Ei))/rooti))
                 -1./(bi*rootj)
                  *log(abs((A+B)/(A-B))
                       *(2.*sqrt(Cj-Dj+Ej)+(B*(2.*Cj-Dj)-A*(Dj-2.*Ej))/rootj)
                       /(2.*sqrt(Cj+Dj+Ej)+(B*(2.*Cj+Dj)-A*(Dj+2.*Ej))/rootj)));
}

