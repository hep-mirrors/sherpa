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
  /////////////////////////////////////////////////////////////////////////////
  // fixing the coupling constants:
  // - if both flavours identical, check if neutral weak/electromagnetic current
  //   and fix the couplings accordingly as a mixture of left- and right-handed
  //   couplings
  // - if lepton-neutrino, it is only left-handed.
  /////////////////////////////////////////////////////////////////////////////
  if (m_flavs[m_indices[0]]==m_flavs[m_indices[1]]) {
    if (m_flavs[m_indices[0]].IsNeutrino()) {
      /////////////////////////////////////////////////////////////////////////
      // TODO: Weak coupling only, to be implemented later  
      /////////////////////////////////////////////////////////////////////////
      THROW(fatal_error,"weak neutral current not yet implemented.")
    }
    else if (m_flavs[m_indices[0]].IsChargedLepton()) {
      /////////////////////////////////////////////////////////////////////////
      // assume QED coupling only for the time being:  -i e e_f gamma^mu
      // TODO: add weak neutral coupling
      /////////////////////////////////////////////////////////////////////////
      m_cL = m_cR = ( -Complex( 0., 1.) *
		      m_flavs[m_indices[0]].Charge() *
		      sqrt(4.*M_PI*alphaQED) );
    }
  }
  else if (m_flavs[m_indices[0]].LeptonFamily()==m_flavs[m_indices[1]].LeptonFamily()) {
    ///////////////////////////////////////////////////////////////////////////
    // Weak left-handed coupling: -i g_W/(2 sqrt(2)) gamma^{mu L}
    ///////////////////////////////////////////////////////////////////////////
    m_cL = Complex( 0., 1. ) * sqrt(4.*M_PI*alphaQED/(sqrt(8.)*sin2thetaW)) ;
    m_cR = Complex( 0., 0. );
  }
  else THROW(fatal_error,"family non-diagonal lepton interaction not yet implemented.")
};

void Lepton_Lepton::Calc(const ATOOLS::Vec4D_Vector& moms, METOOLS::XYZFunc * F)
{
  /////////////////////////////////////////////////////////////////////////
  // J^mu = ubar(0) gamma^mu (P_L c_L + P_R c_R) u(1) 
  /////////////////////////////////////////////////////////////////////////
  Vec4C amp;
  for(int h0=0; h0<2; h0++) {
    for(int h1=0; h1<2; h1++) {
      /////////////////////////////////////////////////////////////////////////
      // L() = ubar(0) gamma^mu (P_L c_L + P_R c_R) u(1) 
      /////////////////////////////////////////////////////////////////////////
      amp = F->L(0,h0, 1,h1, m_cR,m_cL);
      vector<pair<int,int> > spins;
      spins.push_back(make_pair(0,h0));
      spins.push_back(make_pair(1,h1));
      Insert( amp ,spins );
    }
  }
}
