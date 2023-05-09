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
  // As a quick fix, add relevant parameters here.
  // TODO: Will have to make them part of an overall "reduced" model or input 
  //       structure at a later stage.
  /////////////////////////////////////////////////////////////////////////////

  double alphaQED   = 1./137.;
  double sin2thetaW = 0.22290, cos2thetaW = 1.-sin2thetaW;
  double I_f        = -1.0/2.0;
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
      // Weak Neutral coupling: -i g_Z/(2) (gamma^{mu L} + gamma^{mu R}) To neutinos
      /////////////////////////////////////////////////////////////////////////
      
      //TODO Think how to add these if cL and cR are asymmetric...

      // m_cL = ( -Complex( 0., 1.) * 
      //     ((I_f) - (m_flavs[m_indices[0]].Charge())*sin2thetaW) *
		  //     sqrt(4.*M_PI*alphaQED/(2*sin2thetaW*cos2thetaW))
      // );
      // m_cR = Complex( 0., 0.);   

    }
    else if (m_flavs[m_indices[0]].IsChargedLepton()) {
      /////////////////////////////////////////////////////////////////////////
      // QED coupling only for the time being:  -i e e_f gamma^mu To leptons
      /////////////////////////////////////////////////////////////////////////
      m_cL = m_cR = ( -Complex( 0., 1.) *
		      m_flavs[m_indices[0]].Charge() *
		      sqrt(4.*M_PI*alphaQED) );

      /////////////////////////////////////////////////////////////////////////
      // Weak Neutral coupling:  -i g_Z/(2) (gamma^{mu L} + gamma^{mu R}) To leptons
      // TODO Add diagrams together...
      /////////////////////////////////////////////////////////////////////////

      // m_cL = ( -Complex( 0., 1.) * 
      //     ((I_f) - (m_flavs[m_indices[0]].Charge())*sin2thetaW) *
		  //     sqrt(4.*M_PI*alphaQED/(2*sin2thetaW*cos2thetaW))
      // );

      // m_cR = ( -Complex( 0., 1.) * 
      //     (-(m_flavs[m_indices[0]].Charge())*sin2thetaW) *
		  //     sqrt(4.*M_PI*alphaQED/(2*sin2thetaW*cos2thetaW))
      // );    
    }
  }
  else if (m_flavs[m_indices[0]].LeptonFamily()==m_flavs[m_indices[1]].LeptonFamily()) {
    ///////////////////////////////////////////////////////////////////////////
    // Weak left-handed coupling: -i g_W/(2 sqrt(2)) gamma^{mu L}
    ///////////////////////////////////////////////////////////////////////////
    // m_cL = Complex( 0., 1. ) * sqrt(4.*M_PI*alphaQED/(sqrt(8.)*sin2thetaW)) ;
    // m_cR = Complex( 0., 0. );
  }
  else THROW(fatal_error,"family non-diagonal lepton interaction not yet implemented.")
};

void Lepton_Lepton::Calc(const ATOOLS::Vec4D_Vector& moms, METOOLS::XYZFunc * F)
{

  /////////////////////////////////////////////////////////////////////////
  // This assumes the momentum transfer from the other (lepton) current
  // taken as incoming, i.e. p_0 = p_1 + q.
  // Propagator term in Nucleon side.
  /////////////////////////////////////////////////////////////////////////

  const int N  = m_flavs.size();
  const int pf = 0; 
  const int pi = 1;

  /////////////////////////////////////////////////////////////////////////
  // J^mu = ubar(0) gamma^mu (P_L c_L + P_R c_R) u(1) 
  /////////////////////////////////////////////////////////////////////////
  Vec4C amp;
  for(int hf=0; hf<2; hf++) {
    for(int hi=0; hi<2; hi++) {
      /////////////////////////////////////////////////////////////////////////
      // L() = ubar(0) gamma^mu (P_L c_L + P_R c_R) u(1) 
      /////////////////////////////////////////////////////////////////////////
      amp = F->L(pi,hf, pf,hi, m_cR,m_cL);
      vector<pair<int,int> > spins;
      spins.push_back(make_pair(0,hf));
      spins.push_back(make_pair(1,hi));
      Insert( amp ,spins );
    }
  }
}
