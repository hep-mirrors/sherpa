#include "NEUTRINOS++/Current_Library/Lepton_Lepton.H"
#include "ATOOLS/Org/Exception.H"

#include "NEUTRINOS++/Tools/Form_Factor_Parameter_Maps.H"

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

  double alphaQED   = 1.0/137.;
  double sin2thetaW = 0.22290, cos2thetaW = 1.0 - sin2thetaW;

  double e_coupling = sqrt(4.*M_PI*alphaQED);
  double gz_coupling = e_coupling/(sqrt(sin2thetaW*cos2thetaW));
  double gw_coupling = e_coupling/(sqrt(sin2thetaW));
  
  //kf_code of nucleons IN vs OUT
  kf_code IN = m_flavs[m_indices[1]].Kfcode();
  kf_code OUT = m_flavs[m_indices[0]].Kfcode();

  //Turn contributions from currents on or off...
  double QED_ON = ffs->GetModelParms("Bosons", "gamma");
  double Weak_NC_ON = ffs->GetModelParms("Bosons", "Z");
  double Weak_CC_ON = ffs->GetModelParms("Bosons", "W");

  /////////////////////////////////////////////////////////////////////////////
  // fixing the coupling constants:
  // - if both flavours identical, check if neutral weak/electromagnetic current
  //   and fix the couplings accordingly as a mixture of left- and right-handed
  //   couplings
  // - if lepton-neutrino, it is only left-handed.
  /////////////////////////////////////////////////////////////////////////////
  if (IN==OUT) {
    if (Flavour(OUT).IsNeutrino() ) {
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
    else if (Flavour(OUT).IsChargedLepton() || Flavour(OUT).IsBaryon()) {
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
  else if (Flavour(IN).LeptonFamily()==Flavour(OUT).LeptonFamily()) {
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
    Weak_CC_coupling = (-Complex( 0., 1.) * gw_coupling ) / ( sqrt(2.));
    Weak_CC_cR = Complex(0.,0.);
    Weak_CC_cL = Complex(1.,0.);
  }
  else THROW(fatal_error,"family non-diagonal lepton interaction not yet implemented.")

  QED_coupling *= QED_ON;
  Weak_NC_coupling *= Weak_NC_ON;
  Weak_CC_coupling *= Weak_CC_ON;
};

void Lepton_Lepton::Calc(const ATOOLS::Vec4D_Vector& moms, METOOLS::XYZFunc * F)
{
  /////////////////////////////////////////////////////////////////////////
  // This assumes the momentum transfer from the other (lepton) current
  // taken as incoming, i.e. p_0 = p_1 + q.
  // Propagator term in Nucleon side.
  /////////////////////////////////////////////////////////////////////////

  const int N  = m_flavs.size();
  int pf = 0;
  int pi = 1;

  if (m_anti) F->Set_m_Anti(true);
  else        F->Set_m_Anti(false);

  Complex Zero = Complex(0.0,0.0);
  Complex One = Complex(1.0,0.0);

  /////////////////////////////////////////////////////////////////////////
  // J^mu = ubar(0) gamma^mu (P_L c_L + P_R c_R) u(1) 
  // We separate the left and right handed terms 
  /////////////////////////////////////////////////////////////////////////


  Vec4C QED_amp, Weak_NC_amp, Weak_CC_amp;
  for(int hf=0; hf<2; hf++) {
    for(int hi=0; hi<2; hi++) {
      QED_amp *= 0.0;
      Weak_NC_amp *= 0.0;
      Weak_CC_amp *= 0.0;

      QED_amp = F->L(pf,hf, pi,hi, QED_cR, QED_cL);
      if (!m_anti) QED_amp = QED_amp * QED_coupling;
      else         QED_amp = QED_amp * conj(QED_coupling);
      
      Weak_NC_amp = F->L(pf,hf, pi,hi, Weak_NC_cR, Weak_NC_cL);
      if (!m_anti) Weak_NC_amp = Weak_NC_amp * Weak_NC_coupling;
      else         Weak_NC_amp = Weak_NC_amp * conj(Weak_NC_coupling);
    

      Weak_CC_amp = F->L(pf,hf, pi,hi, Weak_CC_cR, Weak_CC_cL);
      if (!m_anti) Weak_CC_amp = Weak_CC_amp * Weak_CC_coupling;
      else         Weak_CC_amp = Weak_CC_amp * conj(Weak_CC_coupling);
    

      // Factor of two to undo spin averaging.
      QED_amp = QED_amp / 2.0;
      Weak_NC_amp = Weak_NC_amp / 2.0;
      Weak_CC_amp = Weak_CC_amp / 2.0;

      // msg_Out() << "Lepton_Lepton\n";
      // msg_Out() 
      //   << "Anti?: " << m_anti << "\n"
      //   << "QED: \n     " 
      //   << "Coupling: " << QED_coupling << " \n     "
      //   << "Left: " << QED_cL << " \n     "
      //   << "Right: " << QED_cR << " \n"
      //   << "Weak_NC: \n     " 
      //   << "Coupling: " << Weak_NC_coupling << " \n     "
      //   << "Left: " << Weak_NC_cL << " \n     "
      //   << "Right: " << Weak_NC_cR << " \n"
      //   << "Weak_CC: \n     " 
      //   << "Coupling: " << Weak_CC_coupling << " \n     "
      //   << "Left: " << Weak_CC_cL << " \n     "
      //   << "Right: " << Weak_CC_cR << " \n\n\n";

      vector<pair<int,int> > spins;
      spins.push_back(make_pair(pf,hf));
      spins.push_back(make_pair(pi,hi));

      //Old Method QED ONLY
      Insert(QED_amp, spins);

      //New Method
      Insert_ProcessType("QED", QED_amp, spins);
      Insert_ProcessType("Weak_NC", Weak_NC_amp, spins);
      Insert_ProcessType("Weak_CC", Weak_CC_amp, spins);
    }
  }
}