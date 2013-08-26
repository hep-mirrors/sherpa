# include "PHOTONS++/Tools/Generate_One_Photon.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Particle.H"
#include "PHOTONS++/PhaseSpace/Generate_Dipole_Photon_Angle.H"
#include "PHOTONS++/PhaseSpace/Generate_Multipole_Photon_Angle.H"

#include "ATOOLS/Math/Histogram_2D.H"

using namespace ATOOLS;
using namespace PHOTONS;
using namespace std;

Histogram_2D *Generate_One_Photon::s_histo = new Histogram_2D(1,0.,M_PI,100,0.,2.*M_PI,100);

// public members
Generate_One_Photon::Generate_One_Photon
(const double& b1, const double& b2, const double& E, Dipole_Type::code dtype) :
  m_w(E), m_dtype(dtype)
{
  // dipole radiation
  if (m_dtype == Dipole_Type::ff) {
    m_beta1 = b1;
    m_beta2 = b2;
  }
  else if (m_dtype == Dipole_Type::fi) {
    // define in different order, because b1 in -z direction in Dipole_FI
    m_beta1 = b2;
    m_beta2 = b1;
  }
  else if (m_dtype == Dipole_Type::ii) {
    m_beta1 = b1;
    m_beta2 = b2;
  }
//   if (s_histo) { 
//     size_t precision(msg->Out().precision());
//     msg->SetPrecision(16);
//     msg_Info()<<"\\beta_1 = "<<m_beta1<<" ,  \\beta_2 = "<<m_beta2<<std::endl;
//     msg->SetPrecision(precision);
//     for (size_t i(0); i<100000000; ++i) {
//       Generate_Dipole_Photon_Angle gdpat(m_beta1,m_beta2);
//       s_histo->Insert(gdpat.GetTheta(),gdpat.GetPhi(),1.,1);
//       if ((i<10000 && (i+1)%1000==0) ||
//           (i<100000 && (i+1)%10000==0) ||
//           (i<1000000 && (i+1)%100000==0) ||
//           (i<10000000 && (i+1)%1000000==0) ||
//           ((i+1)%10000000==0)) msg_Info()<<"Generated "<<(i+1)<<" angles\n";
//     }
//     s_histo->Finalize();
//     s_histo->Output("theta-phi.dat");
//     delete s_histo;
//     s_histo=NULL;
//   }
  Generate_Dipole_Photon_Angle gdpa(m_beta1,m_beta2);
  m_theta = gdpa.GetTheta();
  m_phi   = gdpa.GetPhi();
  GeneratePhoton();
}

Generate_One_Photon::Generate_One_Photon
(const Particle_Vector& dip, const IdPairNbarVector& nbars,
 double E, Dipole_Type::code dtype) :
  m_w(E), m_dtype(dtype)
{
  // multipole radiation
//   if (s_histo) {
//     for (size_t j(0); j<dip.size(); ++j) {
//       msg_Info()<<dip[j]->Momentum()<<std::endl;
//     }
//     for (size_t j(0); j<dip.size(); ++j) {
//       double pspat(dip[j]->Momentum().PSpat());
//       double theta(acos(dip[j]->Momentum()[3]/pspat));
//       double phi(IsZero(theta)?0.:acos(dip[j]->Momentum()[1]/(pspat*sin(theta))));
//       double phi1(IsZero(theta)?0.:acos(dip[j]->Momentum()[1]/sqrt(sqr(dip[j]->Momentum()[1])+sqr(dip[j]->Momentum()[2]))));
//       size_t precision(msg->Out().precision());
//       msg->SetPrecision(16);
//       msg_Info()<<j<<": \\beta = "<<pspat/dip[j]->Momentum()[0]<<" ,  \\theta = "<<theta<<" ,  \\phi = "<<phi<<std::endl;
//       msg->SetPrecision(precision);
//     }
//     for (size_t i(0); i<100000000; ++i) {
//       Generate_Multipole_Photon_Angle gmpat(dip,nbars);
//       s_histo->Insert(gmpat.GetTheta(),gmpat.GetPhi(),1.,1);
//       if ((i<10000 && (i+1)%1000==0) ||
//           (i<100000 && (i+1)%10000==0) ||
//           (i<1000000 && (i+1)%100000==0) ||
//           (i<10000000 && (i+1)%1000000==0) ||
//           ((i+1)%10000000==0)) msg_Info()<<"Generated "<<(i+1)<<" angles\n";
//     }
//     s_histo->Finalize();
//     s_histo->Output("theta-phi.dat");
//     delete s_histo;
//     s_histo=NULL;
//   }
  Generate_Multipole_Photon_Angle gmpa(dip,nbars);
  m_phi   = gmpa.GetPhi();
  m_theta = gmpa.GetTheta();
  GeneratePhoton();
}

Generate_One_Photon::~Generate_One_Photon() {
}

// private members
void Generate_One_Photon::GeneratePhotonAngleMassless() {
// Generation of theta for two massless particles
  double num = ran->Get();
  m_theta = acos(sqrt(1.-(sqr(sin(m_delta))/((1.-num)*sqr(sin(m_delta))+num))));
  if (ran->Get()>=0.5) m_theta = M_PI - m_theta;
}

void Generate_One_Photon::GeneratePhoton() {
  p_photon = new Particle(-1,kf_photon,Vec4D(0.,0.,0.,0.),'S');
  p_photon->SetMomentum(Vec4D(  m_w ,
                                m_w * sin(m_theta) * cos(m_phi) ,
                                m_w * sin(m_theta) * sin(m_phi) ,
                                m_w * cos(m_theta)  ));
}

void Generate_One_Photon::SetMinimalPhotonAngle(const double& del) {
  m_delta = del;
  delete p_photon;
  GeneratePhoton();
}
