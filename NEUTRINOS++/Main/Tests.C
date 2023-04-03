#include "SHERPA/Main/Sherpa.H"
#include "NEUTRINOS++/Main/Tests.H"
#include "NEUTRINOS++/Current_Library/Lepton_Lepton.H"
#include "NEUTRINOS++/Current_Library/Nucleon_Nucleon.H"
#include "NEUTRINOS++/Tools/Form_Factor_Library.H"
#include "HADRONS++/ME_Library/Current_ME.H"
#include "ATOOLS/Org/Message.H"

using namespace NEUTRINOS;
using namespace ATOOLS;
using namespace std;

Tests::Tests(int argc,char* argv[]) {
  small_sherpa_init(argc,argv);
  msg_Info()<<METHOD<<" is starting:\n";
  ep_ep();
}

void Tests::InitVectors(const size_t n) {
  m_flavs.resize(n);
  m_indices.resize(n);
  for (size_t i=0;i<n;i++) m_indices[i] = i;
  m_momenta.resize(n);
}

void Tests::MakeSimpleScatterKinematics(const double & E,const double & costheta,const double & phi) {
  double massNin  = m_flavs[m_indices[1]].HadMass();       // TODO: different masses
  double massNout = m_flavs[m_indices[3]].HadMass();       // TODO: different masses
  double omega    = E*massNout/(E*(1.-costheta)+massNout); // TODO: Fix it!
  double sinphi   = sin(phi), cosphi = cos(phi), sintheta = sqrt(1.-sqr(costheta));
  m_momenta[0]    = E * Vec4D(1.,0.,0.,1.);
  m_momenta[1]    = massNin * Vec4D(1.,0.,0.,0.);
  m_momenta[2]    = omega * Vec4D(1., sintheta*cosphi, sintheta*sinphi, costheta);
  m_momenta[3]    = Vec4D(E+massNout-omega, -sintheta*cosphi, -sintheta*sinphi, E-omega*costheta);
}

void Tests::ep_ep() {
  /////////////////////////////////////////////////////////////////////////////
  // Tests for the electromagnetic current-current interaction, using ep -> ep
  // Sequence is: outgoing then incoming in currents LL and NN
  /////////////////////////////////////////////////////////////////////////////
  InitVectors(4);
  m_flavs.push_back(Flavour(kf_e));      
  m_flavs.push_back(Flavour(kf_e));      
  m_flavs.push_back(Flavour(kf_p_plus)); 
  m_flavs.push_back(Flavour(kf_p_plus)); 

  p_JLL = new Lepton_Lepton(m_flavs,m_indices,"ee");
  p_JNN = new Nucleon_Nucleon(m_flavs,m_indices,"pp");
  p_ME  = new Current_ME(m_flavs, m_indices, "Current_ME");
  p_ME->SetCurrent1(p_JLL);
  p_ME->SetCurrent2(p_JNN);
  
  p_ME->Calculate(m_momenta);
  
  delete p_JLL;
  delete p_JNN;
}

int main(int argc,char* argv[]) { Tests(argc,argv); }

