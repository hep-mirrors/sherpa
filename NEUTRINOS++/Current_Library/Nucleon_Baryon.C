#include "NEUTRINOS++/Current_Library/Nucleon_Baryon.H"
#include "METOOLS/Main/XYZFuncs.H"
#include "ATOOLS/Org/Exception.H"

#include "NEUTRINOS++/Tools/Form_Factor_Parameter_Maps.H"
#include "NEUTRINOS++/Tools/Propagator_Maps.H"
#include "NEUTRINOS++/Tools/Transition_Maps.H"

using namespace NEUTRINOS;
using namespace ATOOLS;
using namespace std;

///////////////////////////////////////////////////////////////////////////////
// The first index must denote the "barred spinor", the second the non-barred one
///////////////////////////////////////////////////////////////////////////////
Nucleon_Baryon::Nucleon_Baryon(const ATOOLS::Flavour_Vector& flavs,
				 const std::vector<int>& indices,
				 const std::string& name) :
  Scatter_Current_Base(flavs, indices, name),
  m_anti(m_flavs[m_indices[0]].IsAnti()),
  m_massin(m_flavs[m_indices[1]].HadMass()),
  m_massout(m_flavs[m_indices[0]].HadMass()) {
  /////////////////////////////////////////////////////////////////////////////
  // Relevant parameters here.
  /////////////////////////////////////////////////////////////////////////////
  double alphaQED   = 1./137.;
  double sin2thetaW = 0.22290, cos2thetaW = 1.-sin2thetaW;

  double e_coupling = sqrt(4.*M_PI*alphaQED);
  double gz_coupling = e_coupling/(sqrt(sin2thetaW*cos2thetaW));
  double gw_coupling = e_coupling/(sqrt(sin2thetaW));

  //kf_code of nucleons IN vs OUT
  kf_code IN = m_flavs[m_indices[1]].Kfcode();
  kf_code OUT = m_flavs[m_indices[0]].Kfcode();

  //Read in model parameters
  std::string Vckm_string = compareQuarkContent(getQuarkContent(IN), getQuarkContent(OUT)); 
  double Vckm = ffs->GetModelParms("CKM", Vckm_string);

  //Turn contributions from currents on or off...
  double QED_ON = ffs->GetModelParms("Bosons", "gamma");
  double Weak_NC_ON = ffs->GetModelParms("Bosons", "Z");
  double Weak_CC_ON = ffs->GetModelParms("Bosons", "W");

  //Get the Clebsch–Gordan coefficients
  a_CG = ffs->GetModelParms("CG", "a");
  b_CG = ffs->GetModelParms("CG", "b");

  /////////////////////////////////////////////////////////////////////////////
  // Form factor info
  /////////////////////////////////////////////////////////////////////////////
  
  cpl_info::code GE = cpl_info::GE, GM = cpl_info::GM, f3 = cpl_info::f3;
  cpl_info::code g1 = cpl_info::g1, g2 = cpl_info::g2, g3 = cpl_info::g3;
  cpl_info::code ffsnull = cpl_info::unknown;

  /////////////////////////////////////////////////////////////////////////////
  // Propagator info
  // We define all propagators here by QED (Photon), NC (Z boson), CC (W boson)
  /////////////////////////////////////////////////////////////////////////////
  prop_type::code prop_type = prop_type::massive;
  kf_code prop_kf_P   = kf_photon;
  kf_code prop_kf_Z   = kf_Z;
  kf_code prop_kf_W   = kf_Wplus;

  kf_code proton_pid  = kf_p_plus;
  kf_code neutron_pid = kf_n;

  kf_code null_pid = kf_none; 

  //Using the EM form factors...
  m_ffs["GE_proton"]  = ffs->GetFF(proton_pid,proton_pid,prop_kf_P,GE); 
  m_ffs["GM_proton"]  = ffs->GetFF(proton_pid,proton_pid,prop_kf_P,GM); 
  m_ffs["f3_proton"]  = ffs->GetFF(proton_pid,proton_pid,prop_kf_P,f3); 

  m_ffs["GE_neutron"]  = ffs->GetFF(neutron_pid,neutron_pid,prop_kf_P,GE); 
  m_ffs["GM_neutron"]  = ffs->GetFF(neutron_pid,neutron_pid,prop_kf_P,GM); 
  m_ffs["f3_neutron"]  = ffs->GetFF(neutron_pid,neutron_pid,prop_kf_P,f3); 

  //Only valid for X->Y & X != Y
  if ( IN != OUT ) {
    m_ffs["g1_protonneutron"]  = ffs->GetFF(proton_pid,neutron_pid,prop_kf_P,g1); 
    m_ffs["g2_protonneutron"]  = ffs->GetFF(proton_pid,neutron_pid,prop_kf_P,g2); 
    m_ffs["g3_protonneutron"]  = ffs->GetFF(proton_pid,neutron_pid,prop_kf_P,g3); 
  } else {
    m_ffs["g1_protonneutron"]  = ffs->GetFF(proton_pid,neutron_pid,prop_kf_P,ffsnull); 
    m_ffs["g2_protonneutron"]  = ffs->GetFF(proton_pid,neutron_pid,prop_kf_P,ffsnull); 
    m_ffs["g3_protonneutron"]  = ffs->GetFF(proton_pid,neutron_pid,prop_kf_P,ffsnull); 
  }

  m_ffprops["QED"] = ffprops->GetProp(prop_kf_P, prop_type);
  m_ffprops["NC"] = ffprops->GetProp(prop_kf_Z, prop_type);
  m_ffprops["CC"] = ffprops->GetProp(prop_kf_W, prop_type);

  if (IN==OUT) {
    /////////////////////////////////////////////////////////////////////////
    // Hadron => Hadron
    /////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////
    // QED coupling: -i e e_f gamma^{mu}
    /////////////////////////////////////////////////////////////////////////
    QED_coupling = ( -Complex( 0., 1.) * Flavour(IN).Charge() *  e_coupling);
    QED_cR = QED_cL = Complex(1.,0.);

    /////////////////////////////////////////////////////////////////////////
    // Weak Neutral coupling: -i g_Z gamma^{mu} (cL P^{L} + cR P^{R})
    /////////////////////////////////////////////////////////////////////////
    double I_f = 1./2.; //TODO: Check the sign of different nucleons
    Weak_NC_coupling = (-Complex( 0., 1.) * gz_coupling );
    Weak_NC_cR = Complex(1.,0.) * -(Flavour(IN).Charge())*sin2thetaW;
    Weak_NC_cL = Complex(1.,0.) * ((I_f) - (Flavour(IN).Charge())*sin2thetaW);

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
    Weak_CC_coupling = (-Complex( 0., 1.) * gw_coupling * fabs(Vckm) ) / (2.0 * sqrt(2.));
    Weak_CC_cR = Complex( 0., 0.);
    Weak_CC_cL = Complex( 1., 0.);
  }
  QED_coupling *= QED_ON;
  Weak_NC_coupling *= Weak_NC_ON;
  Weak_CC_coupling *= Weak_CC_ON;
};

void Nucleon_Baryon::Calc(const ATOOLS::Vec4D_Vector& moms,METOOLS::XYZFunc * F)
{
  /////////////////////////////////////////////////////////////////////////////
  // J^mu = ubar(0) [ gamma^mu F_1(q^2) + i/2 sigma^{mu nu} q_nu F_2(q^2) / m + q^{mu} F_3(q^2) / m] u(1) 
  /////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////
  // This assumes the momentum transfer from the other (lepton) current
  // taken as incoming, i.e. p_0 = p_1 + q.
  // Propagator term in Nucleon side. (HERE)
  /////////////////////////////////////////////////////////////////////////

  const int N  = m_flavs.size();
  int pf = 2; 
  int pi = 3; 

  if (m_anti) F->Set_m_Anti(true);
  else        F->Set_m_Anti(false);
  
  const int pf_bar = pf+N; 
  const int pi_bar = pi+N;

  Complex Zero = Complex(0.,0.);
  Complex One = Complex(1.,0.); 

  const ATOOLS::Vec4<Complex> qmom = (F->P(pf)-F->P(pi));
  const double q2  = qmom.Abs2().real();
  const double tau = -q2/(4.*m_massin*m_massout); 

  const complex QED_prop_factor = m_ffprops["QED"]->Calc(q2);
  const complex Weak_NC_prop_factor = m_ffprops["NC"]->Calc(q2);
  const complex Weak_CC_prop_factor = m_ffprops["CC"]->Calc(q2);

  const double GE_p = m_ffs["GE_proton"]->Calc(-q2), GM_p = m_ffs["GM_proton"]->Calc(-q2);
  const double GE_n = m_ffs["GE_neutron"]->Calc(-q2), GM_n = m_ffs["GM_neutron"]->Calc(-q2);
  
  //f1, f2, f3 for Proton and Neutron
  double f1_p = (GE_p+tau*GM_p)/(1.+tau), f2_p = (GM_p-GE_p)/(1.+tau);
  double f3_p = m_ffs["f3_proton"]->Calc(-q2);
  double f1_n = (GE_n+tau*GM_n)/(1.+tau), f2_n = (GM_n-GE_n)/(1.+tau);
  double f3_n = m_ffs["f3_neutron"]->Calc(-q2);

  //g1, g2, g3 for Proton and Neutron
  double g1_pn = m_ffs["g1_protonneutron"]->Calc(-q2);
  double g2_pn = m_ffs["g2_protonneutron"]->Calc(-q2);
  double g3_pn = m_ffs["g3_protonneutron"]->Calc(-q2);

  /////////////////////////////////////////////////////////////////////////
  // Derive Fi^V and Di^V, Fi^A and Di^A
  // Using definition from neutron and proton EM processes e.g. p -> p and n -> n
  // Fi^V = fi^p + 0.5*fi^n, Di^V = (-3/2)*fi^n
  // Fi^A = gi^p + 0.5*gi^n, Di^A = (-3/2)*gi^n
  /////////////////////////////////////////////////////////////////////////

  //Vector using f1, f2, f3
  double F1V = f1_p + 0.5*f1_n, F2V = f2_p + 0.5*f2_n, F3V = f3_p + 0.5*f3_n;
  double D1V = (-3.0/2.0)*f1_n, D2V = (-3.0/2.0)*f2_n, D3V = (-3.0/2.0)*f3_n;
  //Axial using g1, g2
  double x1_np = 0.364, x2_np = x1_np;
  double F1A = g1_pn * x1_np, F2A = g2_pn * x2_np;
  double D1A = g1_pn * (1-x1_np), D2A = g2_pn * (1-x2_np);

  /////////////////////////////////////////////////////////////////////////
  // Now using Clebsch-Gordan coefficients get relevant form factors for N -> Y.
  /////////////////////////////////////////////////////////////////////////  
  //Vector 
  double f1_NY = a_CG*F1V  + b_CG*D1V, f2_NY = a_CG*F2V  + b_CG*D2V, f3_NY = a_CG*F3V  + b_CG*D3V;
  //Axial
  double g1_NY = a_CG*F1A  + b_CG*D1A, g2_NY = a_CG*F2A  + b_CG*D2A;

  double masskaon = .497648; //TODO Check this
  double g3_NY = g1_NY*sqr(m_massin+m_massout)/(2*(sqr(masskaon)-q2)); //Nambu...

  /////////////////////////////////////////////////////////////////////////////
  // J^mu =  
  //  ubar(0) [ 
  //    (f_1(q^2) gamma^mu  + f_2(q^2) (i/2) sigma^{mu nu} q_nu  / m + f_3(q^2) q^{mu} / m) -
  //    (g_1(q^2) gamma^mu  + g_2(q^2) (i/2) sigma^{mu nu} q_nu  / m + g_3(q^2) q^{mu} / m) gamma^5
  //  ] u(1) 
  /////////////////////////////////////////////////////////////////////////////

  //Remove Z boson...
  Weak_NC_coupling = 0.0;


  Vec4C Gen_amp, QED_amp, Weak_NC_amp, Weak_CC_amp;
  Vec4C f1_term, f2_term, f3_term, g1_term, g2_term, g3_term;
  for(int hf=0; hf<2; hf++) {
    for(int hi=0; hi<2; hi++) {
      QED_amp *= 0.0;
      Weak_NC_amp *= 0.0;
      Weak_CC_amp *= 0.0;

      f1_term *= 0.0;
      f2_term *= 0.0;
      f3_term *= 0.0;
      g1_term *= 0.0;
      g2_term *= 0.0;
      g3_term *= 0.0;

      f1_term = F->L(pf,hf, pi,hi, One,One);
      g1_term = -F->L(pf,hf, pi,hi, One,-One);

      for (int hq=0;hq<2;hq++) {
          f2_term += (
                F->L(pf,hq,pi,hi,One,One) * F->Y(pf,hf,pf,hq,One,One) -
                F->L(pi,hq,pi,hi,One,One) * F->Y(pf,hf,pi,hq,One,One) - 
                F->L(pf,hf,pf,hq,One,One) * F->Y(pf,hq,pi,hi,One,One) +
                F->L(pf,hf,pi,hq,One,One) * F->Y(pi,hq,pi,hi,One,One) 
          );
          g2_term += (
                F->L(pf,hq,pi,hi,One,-One) * F->Y(pf,hf,pf,hq,One,One) -
                F->L(pi,hq,pi,hi,One,-One) * F->Y(pf,hf,pi,hq,One,One) - 
                F->L(pf,hf,pf,hq,One,One) * F->Y(pf,hq,pi,hi,One,-One) +
                F->L(pf,hf,pi,hq,One,One) * F->Y(pi,hq,pi,hi,One,-One)
          );
        }
      f2_term *=  (1/(4.*m_massout));
      g2_term *= -(1/(4.*m_massout));

      f3_term = qmom * F->Y(pf,hf, pi,hi, One,One) / m_massout;
      g3_term = -qmom * F->Y(pf,hf, pi,hi, One,-One) / m_massout;

      //Propagator on Nucleon Current side
      Gen_amp = f1_term*f1_NY + f2_term*f2_NY + f3_term*f3_NY - g1_term*g1_NY - g2_term*g2_NY - g3_term*g3_NY;
      //QED TERM
      if ( fabs(QED_coupling) > 0.0 ) {
        if (!m_anti) QED_amp = Gen_amp * QED_coupling * QED_prop_factor;
        else         QED_amp = Gen_amp * conj(QED_coupling) * QED_prop_factor;
      } 
      //NC TERM
      else if ( fabs(Weak_NC_coupling) > 0.0 ) {
        if (!m_anti) Weak_NC_amp = Gen_amp * Weak_NC_coupling * Weak_NC_prop_factor;
        else         Weak_NC_amp = Gen_amp * conj(Weak_NC_coupling) * Weak_NC_prop_factor;
      }
      //CC TERM
      else if ( fabs(Weak_CC_coupling) > 0.0 ) {
        if (!m_anti) Weak_CC_amp = Gen_amp * Weak_CC_coupling * Weak_CC_prop_factor;
        else         Weak_CC_amp = Gen_amp * conj(Weak_CC_coupling) * Weak_CC_prop_factor;
      }

      // msg_Out() << "Nucleon_Baryon\n";
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
      //   << "Right: " << Weak_CC_cR << " \n"
      //   << "Form factors: \n     " 
      //   << "f1_p: " << f1_p << " \n     "
      //   << "f2_p: " << f2_p << " \n     "
      //   << "f3_p: " << f3_p << " \n     "
      //   << "g1_p: " << g1_p << " \n     "
      //   << "g2_p: " << g2_p << " \n     "
      //   << "g3_p: " << g3_p << " \n\n     "
      //   << "f1_n: " << f1_n << " \n     "
      //   << "f2_n: " << f2_n << " \n     "
      //   << "f3_n: " << f3_n << " \n     "
      //   << "g1_n: " << g1_n << " \n     "
      //   << "g2_n: " << g2_n << " \n     "
      //   << "g3_n: " << g3_n << " \n\n     "
      //   << "f1_NY: " << f1_NY << " \n     "
      //   << "f2_NY: " << f2_NY << " \n     "
      //   << "f3_NY: " << f3_NY << " \n     "
      //   << "g1_NY: " << g1_NY << " \n     "
      //   << "g2_NY: " << g2_NY << " \n     "
      //   << "g3_NY: " << g3_NY << " \n\n\n";
      // exit(1);
      
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