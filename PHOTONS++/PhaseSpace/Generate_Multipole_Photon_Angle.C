#include "PHOTONS++/PhaseSpace/Generate_Multipole_Photon_Angle.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Particle.H"
#include "PHOTONS++/Main/Photons.H"
#include "PHOTONS++/Main/Dipole_Type.H"
#include "PHOTONS++/PhaseSpace/Generate_Dipole_Photon_Angle.H"


using namespace PHOTONS;
using namespace ATOOLS;
using namespace std;

Generate_Multipole_Photon_Angle::Generate_Multipole_Photon_Angle
(const Particle_Vector& dip, const IdPairNbarVector& nbars) :
  m_dipole(dip), m_nbars(nbars), m_posnbar(0.), m_theta(0.), m_phi(0.)
{
  DEBUG_FUNC("");
  msg_Debugging()<<"positive dipoles:\n";
  for (unsigned int i=0; i<nbars.size(); i++) {
    if (nbars[i].nbar>0) {
      m_posnbar += nbars[i].nbar;
      m_posnbars.push_back(nbars[i]);
      msg_Debugging()<<i<<" "<<nbars[i]<<": ["
                     <<m_dipole[m_posnbars.back().ij.i]->Flav()<<","
                     <<m_dipole[m_posnbars.back().ij.j]->Flav()<<"]"<<std::endl;
    }
  }
  msg_Debugging()<<"-> "<<m_posnbar<<std::endl;
  GenerateMultipoleAngle();
}

Generate_Multipole_Photon_Angle::~Generate_Multipole_Photon_Angle() {
}

double Generate_Multipole_Photon_Angle::CalculateBeta(const Vec4D& p) {
  return p.PSpat()/p[0];
}

void Generate_Multipole_Photon_Angle::GenerateMultipoleAngle() {
  // choose according to which eikonal factor angle should be generated
  // p_i and p_j

  // only use pos. dipoles to generate, then reweight to full distribution
  DEBUG_FUNC("");
  bool accept(false);
  while (!accept) {
    double rdm = ran->Get();
    double s   = 0.;
    size_t k   = 0;
    for (size_t i(0); i<m_nbars.size(); i++) {
      s += m_posnbars[i].nbar/m_posnbar;
      if (rdm < s) { k = i; break; }
    }
    msg_Debugging()<<"select "<<m_nbars[k]<<", ["
                              <<m_dipole[m_posnbars[k].ij.i]->Flav()<<","
                              <<m_dipole[m_posnbars[k].ij.j]->Flav()<<"]\n";

    Vec4D p1 = m_dipole[m_nbars[k].ij.i]->Momentum();
    Vec4D p2 = m_dipole[m_nbars[k].ij.j]->Momentum();
    Generate_Dipole_Photon_Angle gdpa(p1,p2);
    double theta = gdpa.GetTheta();
    double phi   = gdpa.GetPhi();
    msg_Debugging()<<"theta = "<<theta<<", phi = "<<phi<<std::endl;

    // reweight pos-dipoles to all dipoles
    // double wgt(CalculateWeightByThetaPhi(theta,phi)); // untested
    double wgt(CalculateWeightByVector(gdpa.GetVector()));
    msg_Debugging()<<"weight = "<<wgt;
    if (wgt>1. || wgt<0.) {
      msg_Debugging()<<" -> reject"<<std::endl;
      msg_Error()<<METHOD<<"(): Internal error. Retry point.\n";
      continue;
    }
    if (wgt > ran->Get()) {
      accept  = true;
      m_theta = theta;
      m_phi   = phi;
      msg_Debugging()<<" -> accept"<<std::endl;
    }
    else msg_Debugging()<<" -> reject"<<std::endl;
  }
}

double Generate_Multipole_Photon_Angle::CalculateWeightByThetaPhi
(const double& theta, const double& phi) {
  // calculate \sum_{i,j} QiQj*(  (1-bi^2)/(1-bi*ei*ek)^2
  //                            + (1-bj^2)/(1-bj*ej*ek)^2
  //                            -2(1-bi*bj*ei*ej)/((1-bi*ei*ek)*(1-bj*ej*ek)))
  // untested
  DEBUG_FUNC("");
  size_t precision(msg->Out().precision());
  msg->SetPrecision(16);
  double all(0.),pos(0.);
  double stk(sin(theta)), spk(sin(phi)), ctk(cos(theta)), cpk(cos(phi));
  for (size_t i(0);i<m_nbars.size();++i) {
    const Vec4D& pi(m_dipole[m_nbars[i].ij.i]->Momentum());
    const Vec4D& pj(m_dipole[m_nbars[i].ij.j]->Momentum());
    double QiQj(m_dipole[m_nbars[i].ij.i]->Flav().Charge()
                *m_dipole[m_nbars[i].ij.j]->Flav().Charge());
    double bi(CalculateBeta(pi)), bj(CalculateBeta(pj));
    msg_Debugging()<<bi<<" "<<bj<<std::endl;
    double thetai(acos(pi[3]/pi.PSpat())), thetaj(acos(pj[3]/pj.PSpat()));
    msg_Debugging()<<thetai<<" "<<thetaj<<std::endl;
    double sti(sin(thetai)), stj(sin(thetaj));
    double cti(cos(thetai)), ctj(cos(thetaj));
    double phii(IsZero(sti)?0.:acos(pi[1]/sqrt(sqr(pi[1])+sqr(pi[2])))),
           phij(IsZero(stj)?0.:acos(pj[1]/sqrt(sqr(pj[1])+sqr(pj[2]))));
    if (pi[2]<0.) phii=-phii;
    if (pj[2]<0.) phij=-phij;
    msg_Debugging()<<phii<<" "<<phij<<std::endl;
    double spi(sin(phii)), spj(sin(phij));
    double cpi(cos(phii)), cpj(cos(phij));
    double eiej(sti*stj*(cpi*cpj+spi*spj)+cti*ctj);
    double eiek(sti*stk*(cpi*cpk+spi*spk)+cti*ctk);
    double ejek(stj*stk*(cpj*cpk+spj*spk)+ctj*ctk);
    double fii((1.-sqr(bi))/sqr(1.-bi*eiek));
    double fjj((1.-sqr(bj))/sqr(1.-bj*ejek));
    double fij((1.-bi*bj*eiej)/((1.-bi*eiek)*(1.-bj*ejek)));
    msg_Debugging()<<fii<<" "<<fjj<<" "<<fij<<std::endl;
    double val(QiQj*(fii+fjj-2.*fij));
    if (IsBad(val)) {
      msg_Error()<<METHOD<<"(): encountered "<<val<<", set to 0"<<std::endl;
      val=0.;
    }
    msg_Debugging()<<i<<" ("<<m_nbars[i].ij.i<<","
                            <<m_nbars[i].ij.j<<"), ["
                            <<m_dipole[m_nbars[i].ij.i]->Flav()<<","
                            <<m_dipole[m_nbars[i].ij.j]->Flav()<<"]: "
                   <<val<<std::endl;
    all+=val;
    if (QiQj<0.) pos+=val;
  }
  msg_Debugging()<<all<<" / "<<pos<<" = "<<all/pos<<std::endl;
  msg->SetPrecision(precision);
  return all/pos;
}

double Generate_Multipole_Photon_Angle::CalculateWeightByVector
(const Vec4D& ek) {
  // calculate \sum_{i,j} QiQj * (pi/pi*ek - pj/pj*ek)^2
  DEBUG_FUNC("");
  size_t precision(msg->Out().precision());
  msg->SetPrecision(16);
  double all(0.),pos(0.);
  bool close(false);
  for (size_t i(0);i<m_nbars.size();++i) {
    if (m_nbars[i].nbar==0.) continue;
    const Vec4D& pi(m_dipole[m_nbars[i].ij.i]->Momentum());
    const Vec4D& pj(m_dipole[m_nbars[i].ij.j]->Momentum());
    double QiQj(m_dipole[m_nbars[i].ij.i]->Flav().Charge()
                *m_dipole[m_nbars[i].ij.j]->Flav().Charge());
    Vec4D eikvec(pi/(pi*ek)-pj/(pj*ek));
    double titj(TiTj(m_nbars[i].ij.i,m_nbars[i].ij.j));
    double val(QiQj*titj*eikvec.Abs2());
    if (IsBad(val)) {
      msg_Error()<<METHOD<<"(): encountered "<<val<<", set to 0"<<std::endl;
      val=0.;
    }
    msg_Debugging()<<i<<" ("<<m_nbars[i].ij.i<<","
                            <<m_nbars[i].ij.j<<"), ["
                            <<m_dipole[m_nbars[i].ij.i]->Flav()<<","
                            <<m_dipole[m_nbars[i].ij.j]->Flav()<<"]: "
                   <<val<<std::endl;
    if (pi.DR(ek)<Photons::s_drcut || pj.DR(ek)<Photons::s_drcut) close=true;
    all+=val;
    if (val>0.) pos+=val;
  }
  msg_Debugging()<<all<<" / "<<pos<<" = "<<all/pos<<std::endl;
  msg->SetPrecision(precision);
  if (!close) return 0.;
  if (all==0. && pos==0.) pos=1.;
  return all/pos;
}

double Generate_Multipole_Photon_Angle::TiTj(const size_t& i, const size_t& j)
{
  if      ((m_dipole[i]->ProductionBlob() == m_dipole[j]->ProductionBlob()) &&
           (m_dipole[i]->ProductionBlob() != NULL))  return  1.;
  else if ((m_dipole[i]->DecayBlob() == m_dipole[j]->ProductionBlob()) &&
           (m_dipole[i]->DecayBlob() != NULL))       return -1.;
  else if ((m_dipole[i]->ProductionBlob() == m_dipole[j]->DecayBlob()) &&
           (m_dipole[i]->ProductionBlob() != NULL))  return -1.;
  else if ((m_dipole[i]->DecayBlob() == m_dipole[j]->DecayBlob()) &&
           (m_dipole[i]->DecayBlob() != NULL))       return  1.;
  else                                               return  0.;
}
