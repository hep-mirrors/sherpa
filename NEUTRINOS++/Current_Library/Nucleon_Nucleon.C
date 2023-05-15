#include "NEUTRINOS++/Current_Library/Nucleon_Nucleon.H"
#include "METOOLS/Main/XYZFuncs.H"
#include "ATOOLS/Org/Exception.H"

#include "NEUTRINOS++/Tools/Form_Factor_Parameter_Maps.H"
#include "NEUTRINOS++/Tools/Propagator_Maps.H"

using namespace NEUTRINOS;
using namespace ATOOLS;
using namespace std;

///////////////////////////////////////////////////////////////////////////////
// The first index must denote the "barred spinor", the second the non-barred one
///////////////////////////////////////////////////////////////////////////////
Nucleon_Nucleon::Nucleon_Nucleon(const ATOOLS::Flavour_Vector& flavs,
				 const std::vector<int>& indices,
				 const std::string& name) :
  Scatter_Current_Base(flavs, indices, name),
  m_massin(m_flavs[m_indices[0]].HadMass()),
  m_massout(m_flavs[m_indices[1]].HadMass()) {
  /////////////////////////////////////////////////////////////////////////////
  // As a quick fix, add relevant parameters here.
  // TODO: Will have to make them part of an overall "reduced" model or input 
  //       structure at a later stage.
  /////////////////////////////////////////////////////////////////////////////

  double alphaQED   = 1./137.;
  double sin2thetaW = 0.22290, cos2thetaW = 1.-sin2thetaW;
  double I_f        = -1.0/2.0;
  kf_code N1        = m_flavs[m_indices[0]].Kfcode(), N2 = m_flavs[m_indices[1]].Kfcode();

  //Read in model parameters
  double Vud = ffs->GetModelParms("Vud");
  double Vuc = ffs->GetModelParms("Vuc");
  double Vub = ffs->GetModelParms("Vub");
  double Vcs = ffs->GetModelParms("Vcs");
  double Vcb = ffs->GetModelParms("Vcb");
  double Vtb = ffs->GetModelParms("Vtb");

  if (m_flavs[m_indices[0]]==m_flavs[m_indices[1]]) {
    ///////////////////////////////////////////////////////////////////////////
    // Electromagnetic interaction (ignoring neutral weak interaction for the
    // time being).
    // TODO: Add weak neutral interaction & form factors - we will have to
    //       find a way to make this "switchable" with an input/model file. 
    ///////////////////////////////////////////////////////////////////////////
    m_cL = m_cR = ( -Complex( 0., 1.) *
		    m_flavs[m_indices[0]].Charge() *
		    sqrt(4.*M_PI*alphaQED) );

    kf_code prop_kf   = kf_photon;
    cpl_info::code GE = cpl_info::GE, GM = cpl_info::GM;
    cpl_info::code A  = cpl_info::axialvector, P = cpl_info::pseudoscalar;
    m_ffs["GE"] = ffs->GetFF(N1,N2,prop_kf,GE); 
    m_ffs["GM"] = ffs->GetFF(N1,N2,prop_kf,GM); 
    m_ffs["A"]  = ffs->GetFF(N1,N2,prop_kf,A);
    m_ffs["P"]  = ffs->GetFF(N1,N2,prop_kf,P);

    prop_type::code prop_type = prop_type::unstable;
    m_ffprops["prop"] = ffprops->GetProp(prop_kf, prop_type);

    /////////////////////////////////////////////////////////////////////////
    // Weak Neutral coupling:  -i g_Z/(2) (gamma^{mu L} + gamma^{mu R})
    /////////////////////////////////////////////////////////////////////////
  
    //TODO Think how to add these if cL and cR are asymmetric...

    // wn_cL = ( -Complex( 0., 1.) * 
    //     ((I_f) - (m_flavs[m_indices[0]].Charge())*sin2thetaW) *
    //     sqrt(4.*M_PI*alphaQED/(2*sin2thetaW*cos2thetaW))
    // );

    // wn_cR = ( -Complex( 0., 1.) * 
    //     (-(m_flavs[m_indices[0]].Charge())*sin2thetaW) *
    //     sqrt(4.*M_PI*alphaQED/(2*sin2thetaW*cos2thetaW))
    // );  
  }
  else {

    /////////////////////////////////////////////////////////////////////////
    // Weak Charged coupling
    /////////////////////////////////////////////////////////////////////////

    m_cL = Complex( 0., 1. ) * sqrt(4.*M_PI*alphaQED/(8.0*sin2thetaW)) ;
    m_cR = Complex( 0., 0. );
    // TODO Vckm factor
    THROW(fatal_error,"current not yet implemented.")
  }
};

void Nucleon_Nucleon::Calc(const ATOOLS::Vec4D_Vector& moms,METOOLS::XYZFunc * F)
{
  /////////////////////////////////////////////////////////////////////////////
  // J^mu = ubar(0) [ gamma^mu F_1(q^2) + i/2 sigma^{mu nu} q_nu F_2(q^2) + q^{mu} F_3(q^2)] u(1) 
  /////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////
  // This assumes the momentum transfer from the other (lepton) current
  // taken as incoming, i.e. p_0 = p_1 + q.
  /////////////////////////////////////////////////////////////////////////

  const int N  = m_flavs.size();
  const int pf = 2; const int pf_bar = 2+N; 
  const int pi = 3; const int pi_bar = 3+N; 

  const ATOOLS::Vec4<Complex> qmom = (F->P(pf)-F->P(pi));
  const double q2  = qmom.Abs2().real();
  const complex prop_factor = m_ffprops["prop"]->Calc(q2);

  const double G_E = m_ffs["GE"]->Calc(-q2), G_M = m_ffs["GM"]->Calc(-q2);
  const double F_A = m_ffs["A"]->Calc(-q2),  F_P = m_ffs["P"]->Calc(-q2);
  const double tau = -q2/(4.*m_massin*m_massout); 

  // msg_Out() << q2 << "\n";
  // msg_Out() << "GE: " << m_ffs["GE"]->Name() << " "<< m_ffs["GE"]->Type() << " " << m_ffs["GE"]->Cpl() << "\n\n";

  // msg_Out() << "qsq= " << 0 << " " << m_ffs["GE"]->Calc(0) << "\n";
  // msg_Out() << "qsq= " << -1 << " " << m_ffs["GE"]->Calc(-1) << "\n";
  // msg_Out() << "qsq= " << -2 << " " << m_ffs["GE"]->Calc(-2) << "\n";
  // msg_Out() << "qsq= " << -3 << " " << m_ffs["GE"]->Calc(-3) << "\n";
  // msg_Out() << "qsq= " << -4 << " " << m_ffs["GE"]->Calc(-4) << "\n";
  // msg_Out() << "qsq= " << -5 << " " << m_ffs["GE"]->Calc(-5) << "\n";
  // exit(1);

  
  //double ff1  = (G_E+tau*G_M)/(1.+tau), ff2 = (G_E-G_M)/(1.+tau);
  double ff1  = (G_E+tau*G_M)/(1.+tau), ff2 = (G_M-G_E)/(1.+tau);
  double ff3  = F_A, ff4 = F_P;
  Complex Zero = Complex(0.,0.);
  Complex One = Complex(1.,0.);

  ff1 = 1.0;
  ff2 = 0.0;
  ff3 = 0.0;
  ff4 = 0.0;
  
  Vec4C amp;
  for(int hf=0; hf<2; hf++) {
    for(int hi=0; hi<2; hi++) {
      if ( ff1 != 0.0 ) {
        amp = ff1 * F->L(pf,hf, pi,hi, m_cR,m_cL);
      } 
      
      if ( ff2 != 0.0 ) {
        for (int hq=0;hq<2;hq++) {
          amp += ff2/(4.*m_massin) * 0.5 *
        (
              -F->Y(pf,hf, pi,hq, m_cR,m_cL) * F->L(pi,hq, pi,hi, One,One) +
              F->Y(pi,hq, pi,hi, m_cR,m_cL) * F->L(pf,hf, pi,hq, One,One) +
              -F->Y(pf,hq, pi,hi, m_cR,m_cL) * F->L(pf,hf, pf,hq, One,One) +
              F->Y(pf,hf, pf,hq, m_cR,m_cL) * F->L(pf,hq, pi,hi, One,One) +

              -F->Y(pf,hf, pi_bar,hq, m_cR,m_cL) * F->L(pi_bar,hq, pi,hi, One,One) +
              F->Y(pi_bar,hq, pi,hi, m_cR,m_cL) * F->L(pf,hf, pi_bar,hq, One,One) +
              -F->Y(pf_bar,hq, pi,hi, m_cR,m_cL) * F->L(pf,hf, pf_bar,hq, One,One) +
              F->Y(pf,hf, pf_bar,hq, m_cR,m_cL) * F->L(pf_bar,hq, pi,hi, One,One) 
        );
        }
      }

      if ( ff3 != 0.0 ) {
        amp += -ff3 * qmom * F->Y(pf,hf, pi,hi, m_cR,m_cL) / m_massin;
      }
      //TODO Add ff4?

      //Divide by propagator on Nucleon Current side
      amp = amp * prop_factor;

      vector<pair<int,int> > spins;
      spins.push_back(make_pair(0,hf));
      spins.push_back(make_pair(1,hi));
      Insert( amp, spins );
    }
  }
}
