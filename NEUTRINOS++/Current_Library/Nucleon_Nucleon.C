#include "NEUTRINOS++/Current_Library/Nucleon_Nucleon.H"
#include "METOOLS/Main/XYZFuncs.H"
#include "ATOOLS/Org/Exception.H"

#include "NEUTRINOS++/Tools/Form_Factor_Parameter_Maps.H"
#include "NEUTRINOS++/Tools/Propagator_Maps.H"

using namespace NEUTRINOS;
using namespace ATOOLS;
using namespace std;

vector<int> getQuarkContent(int x)
{
    //Pull out quark content from kf_code
    vector<int> digits = {};
    for(int di=0; di<7; di++) {
      digits.push_back(x % 10);
      x /= 10;
      if (x == 0) break;
    }
    int length = digits.size();
    return {digits[length-4],digits[length-3],digits[length-2]};
}

std::string compareQuarkContent(vector<int> x1, vector<int> x2)
{
  std::sort(x1.begin(), x1.end());
  std::sort(x2.begin(), x2.end());

  //Get Intersection of quark content
  vector<int> inter;
  set_intersection(
    x1.begin(), x1.end(),
    x2.begin(), x2.end(),
    back_inserter(inter)
  );

  //Remove mutual quark content
  x1.erase(
    set_difference(
      x1.begin(), x1.end(),
      inter.begin(), inter.end(),
      x1.begin()
    ),
    x1.end()
  );

  x2.erase(
    set_difference(
      x2.begin(), x2.end(),
      inter.begin(), inter.end(),
      x2.begin()
    ),
    x2.end()
  );
  if ( x1.size() == 0 && x2.size() == 0 ) return "Vqq";
  if ( x1.size() != 1 | x2.size() != 1 ) THROW(fatal_error,"Flavour changing of more or less than one quarks not implemented");
  if ( x1[0] == 2 | x2[0] == 2 ){ // u
    if ( x1[0] == 1 | x2[0] == 1 ) return "Vud";
    if ( x1[0] == 3 | x2[0] == 3 ) return "Vus";
    if ( x1[0] == 5 | x2[0] == 5 ) return "Vub";
  } else if ( x1[0] == 4 | x2[0] == 4 ){ // c
    if ( x1[0] == 1 | x2[0] == 1 ) return "Vcd";
    if ( x1[0] == 3 | x2[0] == 3 ) return "Vcs";
    if ( x1[0] == 5 | x2[0] == 5 ) return "Vcb";
  } else if ( x1[0] == 6 | x2[0] == 6 ){ // c
    if ( x1[0] == 1 | x2[0] == 1 ) return "Vtd";
    if ( x1[0] == 3 | x2[0] == 3 ) return "Vts";
    if ( x1[0] == 5 | x2[0] == 5 ) return "Vtb";
  }
  return "Vuu_Vdd";
}

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
  // Relevant parameters here.
  /////////////////////////////////////////////////////////////////////////////
  
  double alphaQED   = 1./137.;
  double sin2thetaW = 0.22290, cos2thetaW = 1.-sin2thetaW;

  double e_coupling = sqrt(4.*M_PI*alphaQED);
  double gz_coupling = e_coupling/(sqrt(sin2thetaW*cos2thetaW));
  double gw_coupling = e_coupling/(sqrt(sin2thetaW));

  kf_code N1        = m_flavs[m_indices[0]].Kfcode(), N2 = m_flavs[m_indices[1]].Kfcode();

  //kf_code of nucleons IN vs OUT
  kf_code IN = m_flavs[m_indices[0]].Kfcode();
  kf_code OUT = m_flavs[m_indices[1]].Kfcode();

  //Read in model parameters
  std::string Vckm_string = compareQuarkContent(getQuarkContent(IN), getQuarkContent(OUT)); 
  double Vckm = ffs->GetModelParms(Vckm_string);
  
  cpl_info::code GE = cpl_info::GE, GM = cpl_info::GM;
  cpl_info::code A  = cpl_info::axialvector, P = cpl_info::pseudoscalar;

  /////////////////////////////////////////////////////////////////////////////
  // Propagator info
  // We define all propagators here by QED (Photon), NC (Z boson), CC (W boson)
  /////////////////////////////////////////////////////////////////////////////
  prop_type::code prop_type = prop_type::unstable;
  kf_code prop_kf_P   = kf_photon;
  kf_code prop_kf_Z   = kf_Z;
  kf_code prop_kf_W   = kf_Wplus;

  m_ffs["QED_GE"] = ffs->GetFF(N1,N2,prop_kf_P,GE); 
  m_ffs["QED_GM"] = ffs->GetFF(N1,N2,prop_kf_P,GM); 
  m_ffs["QED_A"]  = ffs->GetFF(N1,N2,prop_kf_P,A);
  m_ffs["QED_P"]  = ffs->GetFF(N1,N2,prop_kf_P,P);
  m_ffprops["QED"] = ffprops->GetProp(prop_kf_P, prop_type);

  m_ffs["NC_GE"] = ffs->GetFF(N1,N2,prop_kf_Z,GE); 
  m_ffs["NC_GM"] = ffs->GetFF(N1,N2,prop_kf_Z,GM); 
  m_ffs["NC_A"]  = ffs->GetFF(N1,N2,prop_kf_Z,A);
  m_ffs["NC_P"]  = ffs->GetFF(N1,N2,prop_kf_Z,P);
  m_ffprops["NC"] = ffprops->GetProp(prop_kf_Z, prop_type);

  m_ffs["CC_GE"] = ffs->GetFF(N1,N2,prop_kf_W,GE); 
  m_ffs["CC_GM"] = ffs->GetFF(N1,N2,prop_kf_W,GM); 
  m_ffs["CC_A"]  = ffs->GetFF(N1,N2,prop_kf_W,A);
  m_ffs["CC_P"]  = ffs->GetFF(N1,N2,prop_kf_W,P);
  m_ffprops["CC"] = ffprops->GetProp(prop_kf_W, prop_type);

  if (m_flavs[m_indices[0]]==m_flavs[m_indices[1]]) {
    /////////////////////////////////////////////////////////////////////////
    // Hadron => Hadron
    /////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////
    // QED coupling: -i e e_f gamma^{mu}
    /////////////////////////////////////////////////////////////////////////
    QED_coupling = ( -Complex( 0., 1.) * m_flavs[m_indices[0]].Charge() *  e_coupling);
    QED_cR = QED_cL = Complex(1.,0.);

    /////////////////////////////////////////////////////////////////////////
    // Weak Neutral coupling: -i g_Z gamma^{mu} (cL P^{L} + cR P^{R})
    /////////////////////////////////////////////////////////////////////////
    double I_f = 1./2.; //TODO: Check the sign of different nucleons
    Weak_NC_coupling = (-Complex( 0., 1.) * gz_coupling );
    Weak_NC_cR = Complex(1.,0.) * -(m_flavs[m_indices[0]].Charge())*sin2thetaW;
    Weak_NC_cL = Complex(1.,0.) * ((I_f) - (m_flavs[m_indices[0]].Charge())*sin2thetaW);

    /////////////////////////////////////////////////////////////////////////
    // Weak Charged (left-handed) coupling: 0
    /////////////////////////////////////////////////////////////////////////
    Weak_CC_coupling = Weak_CC_cR = Weak_CC_cL = Complex(0.,0.);
  }
  else {
    /////////////////////////////////////////////////////////////////////////
    // Hadron_A => Hadron_B
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
    Weak_CC_coupling = (-Complex( 0., 1.) * gw_coupling * Vckm ) / (sqrt(2.));
    Weak_CC_cR = Complex( 0., 0.);
    Weak_CC_cL = Complex( 1., 0.);
  }
};

void Nucleon_Nucleon::Calc(const ATOOLS::Vec4D_Vector& moms,METOOLS::XYZFunc * F)
{
  /////////////////////////////////////////////////////////////////////////////
  // J^mu = ubar(0) [ gamma^mu F_1(q^2) + i/2 sigma^{mu nu} q_nu F_2(q^2) / m + q^{mu} F_3(q^2) / m] u(1) 
  /////////////////////////////////////////////////////////////////////////////
  // , std::string Diagram_Type
  std::string Diagram_Type = "QED";

  /////////////////////////////////////////////////////////////////////////
  // This assumes the momentum transfer from the other (lepton) current
  // taken as incoming, i.e. p_0 = p_1 + q.
  // Propagator term in Nucleon side. (HERE)
  /////////////////////////////////////////////////////////////////////////

  const int N  = m_flavs.size();
  const int pf = 2; const int pf_bar = 2+N; 
  const int pi = 3; const int pi_bar = 3+N;

  Complex Zero = Complex(0.,0.);
  Complex One = Complex(1.,0.); 

  const ATOOLS::Vec4<Complex> qmom = (F->P(pf)-F->P(pi));
  const double q2  = qmom.Abs2().real();

  const complex prop_factor = m_ffprops[Diagram_Type]->Calc(q2);

  const double G_E = m_ffs[Diagram_Type+"_GE"]->Calc(-q2), G_M = m_ffs[Diagram_Type+"_GM"]->Calc(-q2);
  const double F_A = m_ffs[Diagram_Type+"_A"]->Calc(-q2),  F_P = m_ffs[Diagram_Type+"_P"]->Calc(-q2);
  const double tau = -q2/(4.*m_massin*m_massout); 
 
  double ff1  = (G_E+tau*G_M)/(1.+tau), ff2 = (G_M-G_E)/(1.+tau);
  double ff3  = F_A, ff4 = F_P;
  

  // double mass_gordon = sqrt((m_massin*m_massin + m_massout*m_massout) / 2.0);

  // ff1 = 0.0;
  // ff2 = 1.0;
  // ff3 = 0.0;
  // ff4 = 0.0;
  
  Vec4C amp;
  for(int hf=0; hf<2; hf++) {
    for(int hi=0; hi<2; hi++) {
      amp *= 0.0;
      
      if ( ff1 != 0.0 ) {
        amp += ff1 * F->L(pf,hf, pi,hi, One,One);
      } 
      
      if ( ff2 != 0.0 ) {
        for (int hq=0;hq<2;hq++) {
          amp += (ff2/(4.*m_massin)) *
            (
                F->L(pf,hq,pi,hi,One,One) * F->Y(pf,hf,pf,hq,One,One) -
                F->L(pi,hq,pi,hi,One,One) * F->Y(pf,hf,pi,hq,One,One) - 
                F->L(pf,hf,pf,hq,One,One) * F->Y(pf,hq,pi,hi,One,One) +
                F->L(pf,hf,pi,hq,One,One) * F->Y(pi,hq,pi,hi,One,One)
          );
        }
      }

      if ( ff3 != 0.0 ) {
        amp += ff3 * qmom * F->Y(pf,hf, pi,hi, One,One) / m_massin;
      }
      //TODO: Add ff4?

      if ( Diagram_Type == "QED" && fabs(QED_coupling) > 0.0 ) {
        amp = amp * QED_coupling;
      } 
      else if ( Diagram_Type == "NC" && fabs(Weak_NC_coupling) > 0.0 ) {
        amp = amp * Weak_NC_coupling;
      }
      else if ( Diagram_Type == "CC" && fabs(Weak_CC_coupling) > 0.0 ) {
        amp = amp * Weak_CC_coupling;
      }
      
      //Propagator on Nucleon Current side
      amp = amp * prop_factor;

      vector<pair<int,int> > spins;
      spins.push_back(make_pair(pf,hf));
      spins.push_back(make_pair(pi,hi));
      Insert( amp, spins );
    }
  }
}