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

  cpl_info::code GE = cpl_info::GE, GM = cpl_info::GM;
  cpl_info::code A  = cpl_info::axialvector, P = cpl_info::pseudoscalar;

  prop_type::code prop_type = prop_type::unstable;

  if (m_flavs[m_indices[0]]==m_flavs[m_indices[1]]) {
    ///////////////////////////////////////////////////////////////////////////
    // Electromagnetic interaction (ignoring neutral weak interaction for the
    // time being).
    // TODO: Add weak neutral interaction & form factors - we will have to
    //       find a way to make this "switchable" with an input/model file. 
    ///////////////////////////////////////////////////////////////////////////
    QED_coupling = ( -Complex( 0., 1.) * m_flavs[m_indices[0]].Charge() * sqrt(4.*M_PI*alphaQED) );
    QED_cR = QED_cL = Complex(1.,0.);

    kf_code prop_kf_P   = kf_photon;
    m_ffs["GE_P"] = ffs->GetFF(N1,N2,prop_kf_P,GE); 
    m_ffs["GM_P"] = ffs->GetFF(N1,N2,prop_kf_P,GM); 
    m_ffs["A_P"]  = ffs->GetFF(N1,N2,prop_kf_P,A);
    m_ffs["P_P"]  = ffs->GetFF(N1,N2,prop_kf_P,P);

    m_ffprops["prop_P"] = ffprops->GetProp(prop_kf_P, prop_type);

    /////////////////////////////////////////////////////////////////////////
    // Weak Neutral coupling:  -i g_Z/(2) (gamma^{mu L} + gamma^{mu R})
    /////////////////////////////////////////////////////////////////////////

    Weak_coupling = (-Complex( 0., 1.) * sqrt(4.*M_PI*alphaQED/(2*sin2thetaW*cos2thetaW)));
    Weak_cR = Complex(1.,0.) * -(m_flavs[m_indices[0]].Charge())*sin2thetaW;
    Weak_cL = Complex(1.,0.) * ((I_f) - (m_flavs[m_indices[0]].Charge())*sin2thetaW);

    kf_code prop_kf_Z   = kf_Z;
    m_ffs["GE_Z"] = ffs->GetFF(N1,N2,prop_kf_Z,GE); 
    m_ffs["GM_Z"] = ffs->GetFF(N1,N2,prop_kf_Z,GM); 
    m_ffs["A_Z"]  = ffs->GetFF(N1,N2,prop_kf_Z,A);
    m_ffs["P_Z"]  = ffs->GetFF(N1,N2,prop_kf_Z,P);

    m_ffprops["prop_Z"] = ffprops->GetProp(prop_kf_Z, prop_type);
  }
  else {

    /////////////////////////////////////////////////////////////////////////
    // Weak Charged coupling
    /////////////////////////////////////////////////////////////////////////
    QED_coupling = Complex( 0., 0.);
    QED_cR = QED_cL = Complex(0.,0.);

    //TODO How to pick this...? Check quark content between in and out particles?
    double Vckm = 1;

    Weak_coupling = sqrt(4.*M_PI*alphaQED/(8.0*sin2thetaW))*Vckm;
    Weak_cR = 0.;
    Weak_cL = 1.;

    kf_code prop_kf_W   = kf_Wplus;
    m_ffs["GE_W"] = ffs->GetFF(N1,N2,prop_kf_W,GE); 
    m_ffs["GM_W"] = ffs->GetFF(N1,N2,prop_kf_W,GM); 
    m_ffs["A_W"]  = ffs->GetFF(N1,N2,prop_kf_W,A);
    m_ffs["P_W"]  = ffs->GetFF(N1,N2,prop_kf_W,P);

    m_ffprops["prop_W"] = ffprops->GetProp(prop_kf_W, prop_type);
  }
};

void Nucleon_Nucleon::Calc(const ATOOLS::Vec4D_Vector& moms,METOOLS::XYZFunc * F)
{
  //JW: TODO. I think we need to specify the diagrams carefully here so that we're only multiplying the currents between the same diagrams?
  
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
  const complex prop_factor = m_ffprops["prop_P"]->Calc(q2);

  const double G_E = m_ffs["GE_P"]->Calc(-q2), G_M = m_ffs["GM_P"]->Calc(-q2);
  const double F_A = m_ffs["A_P"]->Calc(-q2),  F_P = m_ffs["P_P"]->Calc(-q2);
  const double tau = -q2/(4.*m_massin*m_massout); 
 
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
      amp *= 0.0;

      if ( ff1 != 0.0 ) {
        amp += ff1 * F->L(pf,hf, pi,hi, One,One);
      } 
      
      if ( ff2 != 0.0 ) {
        for (int hq=0;hq<2;hq++) {
          amp += (ff2/(4.*m_massin)) * 0.5 *
            (
                F->L(pf,hq,pi,hi,One,One) * F->Y(pf,hf,pf,hq,One,One) -
                F->L(pi,hq,pi,hi,One,One) * F->Y(pf,hf,pi,hq,One,One) - 
                F->L(pf,hf,pf,hq,One,One) * F->Y(pf,hq,pi,hi,One,One) +
                F->L(pf,hf,pi,hq,One,One) * F->Y(pi,hq,pi,hi,One,One) +

                F->L(pf_bar,hq,pi,hi,One,One) * F->Y(pf,hf,pf_bar,hq,One,One) -
                F->L(pi_bar,hq,pi,hi,One,One) * F->Y(pf,hf,pi_bar,hq,One,One) - 
                F->L(pf,hf,pf_bar,hq,One,One) * F->Y(pf_bar,hq,pi,hi,One,One) +
                F->L(pf,hf,pi_bar,hq,One,One) * F->Y(pi_bar,hq,pi,hi,One,One)
          );
        }
      }

      if ( ff3 != 0.0 ) {
        amp += ff3 * qmom * F->Y(pf,hf, pi,hi, One,One) / m_massin;
      }
      //TODO Add ff4?

      //Divide by propagator on Nucleon Current side
      amp = amp * QED_coupling * prop_factor / 2.0; //JW: Factor of two required...

      vector<pair<int,int> > spins;
      spins.push_back(make_pair(pf,hf));
      spins.push_back(make_pair(pi,hi));
      Insert( amp, spins );
    }
  }
}
