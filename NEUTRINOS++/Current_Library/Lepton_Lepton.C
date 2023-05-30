#include "NEUTRINOS++/Current_Library/Lepton_Lepton.H"
#include "ATOOLS/Org/Exception.H"

using namespace NEUTRINOS;
using namespace ATOOLS;
using namespace std;

///////////////////////////////////////////////////////////////////////////////
// The first index must denote the "barred spinor", the second the non-barred one
///////////////////////////////////////////////////////////////////////////////
Lepton_Lepton::Lepton_Lepton(const ATOOLS::Flavour_Vector& flavs,
			     const std::vector<int>& indices,
			     const std::string& name) :
  Scatter_Current_Base(flavs, indices, name),
  m_anti(m_flavs[m_indices[0]].IsAnti()) {
  /////////////////////////////////////////////////////////////////////////////
  // Relevant parameters here.
  /////////////////////////////////////////////////////////////////////////////

  double alphaQED   = 1./137.;
  double sin2thetaW = 0.22290, cos2thetaW = 1. - sin2thetaW;

  double e_coupling = sqrt(4.*M_PI*alphaQED);
  double gz_coupling = e_coupling/(sqrt(sin2thetaW*cos2thetaW));
  double gw_coupling = e_coupling/(sqrt(sin2thetaW));

  kf_code N1        = m_flavs[m_indices[0]].Kfcode(), N2 = m_flavs[m_indices[1]].Kfcode();

  /////////////////////////////////////////////////////////////////////////////
  // fixing the coupling constants:
  // - if both flavours identical, check if neutral weak/electromagnetic current
  //   and fix the couplings accordingly as a mixture of left- and right-handed
  //   couplings
  // - if lepton-neutrino, it is only left-handed.
  /////////////////////////////////////////////////////////////////////////////
  if (m_flavs[m_indices[0]]==m_flavs[m_indices[1]]) {
    if (m_flavs[m_indices[0]].IsNeutrino() ) {
      /////////////////////////////////////////////////////////////////////////
      // Neutrino => Neutrino
      /////////////////////////////////////////////////////////////////////////

      /////////////////////////////////////////////////////////////////////////
      // QED coupling:  0
      /////////////////////////////////////////////////////////////////////////
      QED_coupling = QED_cR = QED_cL = Complex( 0., 0.);

      /////////////////////////////////////////////////////////////////////////
      // Weak Neutral coupling: -i g_Z gamma^{mu} (cL P^{L} + cR P^{R})
      /////////////////////////////////////////////////////////////////////////
      Weak_NC_coupling = (-Complex( 0., 1.) * gz_coupling );
      Weak_NC_cR = Complex(0.,0.);
      Weak_NC_cL = Complex(1.,0.) * (1./2.);

      /////////////////////////////////////////////////////////////////////////
      // Weak Charged (left-handed) coupling: 0
      /////////////////////////////////////////////////////////////////////////
      Weak_CC_coupling = Weak_CC_cR = Weak_CC_cL = Complex(0.,0.);
    }
    else if (m_flavs[m_indices[0]].IsChargedLepton() || m_flavs[m_indices[0]].IsBaryon()) {
      /////////////////////////////////////////////////////////////////////////
      // Lepton => Lepton
      /////////////////////////////////////////////////////////////////////////

      /////////////////////////////////////////////////////////////////////////
      // QED coupling: -i e e_f gamma^{mu}
      /////////////////////////////////////////////////////////////////////////
      QED_coupling = ( -Complex( 0., 1.) * m_flavs[m_indices[0]].Charge() * e_coupling );
      QED_cR = QED_cL = Complex(1.,0.);

      /////////////////////////////////////////////////////////////////////////
      // Weak Neutral coupling: -i g_Z gamma^{mu} (cL P^{L} + cR P^{R})
      /////////////////////////////////////////////////////////////////////////
      Weak_NC_coupling = (-Complex( 0., 1.) * gz_coupling );
      Weak_NC_cR = Complex(1.,0.) * -(m_flavs[m_indices[0]].Charge())*sin2thetaW;
      Weak_NC_cL = Complex(1.,0.) * ((-1./2.) - (m_flavs[m_indices[0]].Charge())*sin2thetaW);

      /////////////////////////////////////////////////////////////////////////
      // Weak Charged (left-handed) coupling: 0
      /////////////////////////////////////////////////////////////////////////
      Weak_CC_coupling = Weak_CC_cR = Weak_CC_cL = Complex(0.,0.);
    }
  }
  else if (m_flavs[m_indices[0]].LeptonFamily()==m_flavs[m_indices[1]].LeptonFamily()) {
    /////////////////////////////////////////////////////////////////////////
    // Lepton_i => Lepton_f 
    /////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////
    // QED coupling: 0
    /////////////////////////////////////////////////////////////////////////    
    QED_coupling = QED_cR = QED_cL = Complex( 0., 0.);

    /////////////////////////////////////////////////////////////////////////
    // Weak Neutral coupling: 0
    /////////////////////////////////////////////////////////////////////////
    Weak_NC_coupling = Weak_NC_cR = Weak_NC_cL = Complex(0.,0.);

    ///////////////////////////////////////////////////////////////////////////
    // Weak Charged (left-handed) coupling: -i g_W gamma^{mu} / (sqrt(2)) * (cL P^{L}) 
    ///////////////////////////////////////////////////////////////////////////
    Weak_CC_coupling = (-Complex( 0., 1.) * gw_coupling) / (sqrt(2.));
    Weak_CC_cR = Complex(0.,0.);
    Weak_CC_cL = Complex(1.,0.);

  }
  else THROW(fatal_error,"family non-diagonal lepton interaction not yet implemented.")
};

void Lepton_Lepton::Calc(const ATOOLS::Vec4D_Vector& moms, METOOLS::XYZFunc * F)
{

  // , std::string Diagram_Type
  std::string Diagram_Type = "QED";

  /////////////////////////////////////////////////////////////////////////
  // This assumes the momentum transfer from the other (lepton) current
  // taken as incoming, i.e. p_0 = p_1 + q.
  // Propagator term in Nucleon side.
  /////////////////////////////////////////////////////////////////////////

  const int N  = m_flavs.size();
  const int pf = 0; 
  const int pi = 1;

  Complex Zero = Complex(0.,0.);
  Complex One = Complex(1.,0.);

  /////////////////////////////////////////////////////////////////////////
  // J^mu = ubar(0) gamma^mu (P_L c_L + P_R c_R) u(1) 
  // We separate the left and right handed terms 
  /////////////////////////////////////////////////////////////////////////

  Vec4C amp;
  for(int hf=0; hf<2; hf++) {
    for(int hi=0; hi<2; hi++) {
      amp *= 0.0;

      if ( Diagram_Type == "QED" && fabs(QED_coupling) > 0.0 ) {
        amp += F->L(pi,hi, pf,hf, QED_cR, QED_cL);
        amp = amp * QED_coupling;
      }
      else if ( Diagram_Type == "NC" && fabs(Weak_NC_coupling) > 0.0 ) {
        amp += F->L(pi,hi, pf,hf, Weak_NC_cR, Weak_NC_cL);
        amp = amp * Weak_NC_coupling;
      }
      else if ( Diagram_Type == "CC" && fabs(Weak_CC_coupling) > 0.0 ) {
        amp += F->L(pi,hi, pf,hf, Weak_CC_cR, Weak_CC_cL);
        amp = amp * Weak_CC_coupling;
      }

      // Factor of two to undo spin averaging.
      amp = amp / 2.0;

      vector<pair<int,int> > spins;
      spins.push_back(make_pair(pf,hf));
      spins.push_back(make_pair(pi,hi));
      Insert( amp ,spins );
    }
  }
}
