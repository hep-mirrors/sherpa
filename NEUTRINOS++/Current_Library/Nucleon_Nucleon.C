#include "NEUTRINOS++/Current_Library/Nucleon_Nucleon.H"
#include "METOOLS/Main/XYZFuncs.H"
#include "ATOOLS/Org/Exception.H"

using namespace NEUTRINOS;
using namespace ATOOLS;
using namespace std;

///////////////////////////////////////////////////////////////////////////////
// The first index must denote the "barred spinor", the second the non-barred one
///////////////////////////////////////////////////////////////////////////////
Nucleon_Nucleon::Nucleon_Nucleon(const ATOOLS::Flavour_Vector& flavs,
				 const std::vector<int>& indices,
				 const std::string& name) :
  Scatter_Current_Base(flavs, indices, name) {
  /////////////////////////////////////////////////////////////////////////////
  // As a quick fix, add relevant parameters here.
  // TODO: Will have to make them part of an overall "reduced" model or input 
  //       structure at a later stage.
  /////////////////////////////////////////////////////////////////////////////
  double alphaQED   = 1./137.;
  double sin2thetaW = 0.22290, cos2thetaW = 1.-sin2thetaW;
  if (m_flavs[m_indices[2]]==m_flavs[m_indices[3]]) {
    ///////////////////////////////////////////////////////////////////////////
    // Electromagnetic interaction (ignoring neutral weak interaction for the
    // time being).
    // TODO: Add weak neutral interaction & form factors - we will have to
    //       find a way to make this "switchable" with an input/model file. 
    ///////////////////////////////////////////////////////////////////////////
    m_cL = m_cR = ( -Complex( 0., 1.) *
		    m_flavs[m_indices[0]].Charge() *
		    sqrt(4.*M_PI*alphaQED) );
    if (dabs(m_flavs[m_indices[0]].Charge())>0.) {
      /////////////////////////////////////////////////////////////////////////
      // TODO: Charged baryon (most likely proton)
      /////////////////////////////////////////////////////////////////////////
      //THROW(fatal_error,"weak neutral current not yet implemented.")
    }
    else {
      /////////////////////////////////////////////////////////////////////////
      // TODO: Neutral baryon (most likely neutron)
      /////////////////////////////////////////////////////////////////////////
      //THROW(fatal_error,"weak neutral current not yet implemented.")
    }
  }
  else THROW(fatal_error,"current not yet implemented.")
};

void Nucleon_Nucleon::Calc(const ATOOLS::Vec4D_Vector& moms,METOOLS::XYZFunc * F)
{
  /////////////////////////////////////////////////////////////////////////////
  // J^mu = ubar(0) [ gamma^mu F_1(q^2) + i/2 sigma^{mu nu} q_nu F_2(q^2)] u(1) 
  /////////////////////////////////////////////////////////////////////////////
  Vec4C amp;
  for(int h0=0; h0<2; h0++) {
    for(int h1=0; h1<2; h1++) {
      /////////////////////////////////////////////////////////////////////////
      // L() = ubar(0) gamma^mu (c_L+c_R) u(1) 
      /////////////////////////////////////////////////////////////////////////
      amp = F->L(2,h0, 3,h1, m_cR,m_cL);
      vector<pair<int,int> > spins;
      spins.push_back(make_pair(0,h0));
      spins.push_back(make_pair(1,h1));
      Insert( amp ,spins );
    }
  }
}
