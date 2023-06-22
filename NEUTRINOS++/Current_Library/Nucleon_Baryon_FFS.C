#include "NEUTRINOS++/Current_Library/Nucleon_Baryon_FFS.H"
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
Nucleon_Baryon_FFS::Nucleon_Baryon_FFS(const ATOOLS::Flavour_Vector& flavs,
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

  /////////////////////////////////////////////////////////////////////////////
  // Form factor info
  /////////////////////////////////////////////////////////////////////////////
  
  cpl_info::code f1 = cpl_info::f1, f2 = cpl_info::f2, f3 = cpl_info::f3;
  cpl_info::code g1 = cpl_info::g1, g2 = cpl_info::g2, g3 = cpl_info::g3;

  /////////////////////////////////////////////////////////////////////////////
  // Propagator info
  // We define all propagators here by QED (Photon), NC (Z boson), CC (W boson)
  /////////////////////////////////////////////////////////////////////////////
  //prop_type::code prop_type = prop_type::unstable;
  prop_type::code prop_type = prop_type::massive;
  kf_code prop_kf_P   = kf_photon;
  kf_code prop_kf_Z   = kf_Z;
  kf_code prop_kf_W   = kf_Wplus;

  m_ffs["QED_f1"]  = ffs->GetFF(IN,OUT,prop_kf_P,f1); 
  m_ffs["QED_f2"]  = ffs->GetFF(IN,OUT,prop_kf_P,f2); 
  m_ffs["QED_f3"]  = ffs->GetFF(IN,OUT,prop_kf_P,f3);
  m_ffs["QED_g1"]  = ffs->GetFF(IN,OUT,prop_kf_P,g1); 
  m_ffs["QED_g2"]  = ffs->GetFF(IN,OUT,prop_kf_P,g2); 
  m_ffs["QED_g3"]  = ffs->GetFF(IN,OUT,prop_kf_P,g3);
  m_ffprops["QED"] = ffprops->GetProp(prop_kf_P, prop_type);

  m_ffs["NC_f1"]  = ffs->GetFF(IN,OUT,prop_kf_Z,f1); 
  m_ffs["NC_f2"]  = ffs->GetFF(IN,OUT,prop_kf_Z,f2); 
  m_ffs["NC_f3"]  = ffs->GetFF(IN,OUT,prop_kf_Z,f3);
  m_ffs["NC_g1"]  = ffs->GetFF(IN,OUT,prop_kf_Z,g1); 
  m_ffs["NC_g2"]  = ffs->GetFF(IN,OUT,prop_kf_Z,g2); 
  m_ffs["NC_g3"]  = ffs->GetFF(IN,OUT,prop_kf_Z,g3);
  m_ffprops["NC"] = ffprops->GetProp(prop_kf_Z, prop_type);

  m_ffs["CC_f1"]  = ffs->GetFF(IN,OUT,prop_kf_W,f1); 
  m_ffs["CC_f2"]  = ffs->GetFF(IN,OUT,prop_kf_W,f2); 
  m_ffs["CC_f3"]  = ffs->GetFF(IN,OUT,prop_kf_W,f3);
  m_ffs["CC_g1"]  = ffs->GetFF(IN,OUT,prop_kf_W,g1); 
  m_ffs["CC_g2"]  = ffs->GetFF(IN,OUT,prop_kf_W,g2); 
  m_ffs["CC_g3"]  = ffs->GetFF(IN,OUT,prop_kf_W,g3);
  m_ffprops["CC"] = ffprops->GetProp(prop_kf_W, prop_type);

  // double GF = (sqrt(2)/8) * sqr(gw_coupling)/sqr(Flavour(prop_kf_W).Mass());
  // msg_Out() << Flavour(prop_kf_W).Mass() << " " << GF << " " << gw_coupling << "\n";

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
    double Fudge = 1;//sqrt(sqrt(sqrt(sin2thetaW)));
    Weak_CC_coupling = (-Complex( 0., 1.) * gw_coupling * fabs(Vckm) ) / (2.0 * sqrt(2.) * Fudge);
    Weak_CC_cR = Complex( 0., 0.);
    Weak_CC_cL = Complex( 1., 0.);
  }
  QED_coupling *= QED_ON;
  Weak_NC_coupling *= Weak_NC_ON;
  Weak_CC_coupling *= Weak_CC_ON;
};

void Nucleon_Baryon_FFS::Calc(const ATOOLS::Vec4D_Vector& moms,METOOLS::XYZFunc * F)
{
  /////////////////////////////////////////////////////////////////////////////
  // J^mu =  
  //  ubar(0) [ 
  //    (f_1(q^2) gamma^mu  + f_2(q^2) (i/2) sigma^{mu nu} q_nu  / m + f_3(q^2) q^{mu} / m) -
  //    (g_1(q^2) gamma^mu  + g_2(q^2) (i/2) sigma^{mu nu} q_nu  / m + g_3(q^2) q^{mu} / m) gamma^5
  //  ] u(1) 
  /////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////
  // This assumes the momentum transfer from the other (lepton) current
  // taken as incoming, i.e. p_0 = p_1 + q.
  // Propagator term in Nucleon side. (HERE)
  /////////////////////////////////////////////////////////////////////////

  const int N  = m_flavs.size();
  int pf = 2; 
  int pi = 3; 
  const int pf_bar = pf+N; 
  const int pi_bar = pi+N;

  if (m_anti) F->Set_m_Anti(true);
  else        F->Set_m_Anti(false);

  Complex Zero = Complex(0.,0.);
  Complex One = Complex(1.,0.); 

  const ATOOLS::Vec4<Complex> qmom = (F->P(pf)-F->P(pi));
  const double q2  = qmom.Abs2().real();

  const complex QED_prop_factor = m_ffprops["QED"]->Calc(q2);
  const complex Weak_NC_prop_factor = m_ffprops["NC"]->Calc(q2);
  const complex Weak_CC_prop_factor = m_ffprops["CC"]->Calc(q2);

  const double f1_QED = m_ffs["QED_f1"]->Calc(-q2);
  const double f2_QED = m_ffs["QED_f2"]->Calc(-q2);
  const double f3_QED = m_ffs["QED_f3"]->Calc(-q2);
  const double g1_QED = m_ffs["QED_g1"]->Calc(-q2);
  const double g2_QED = m_ffs["QED_g2"]->Calc(-q2);
  const double g3_QED = m_ffs["QED_g3"]->Calc(-q2);

  const double f1_NC = m_ffs["NC_f1"]->Calc(-q2);
  const double f2_NC = m_ffs["NC_f2"]->Calc(-q2);
  const double f3_NC = m_ffs["NC_f3"]->Calc(-q2);
  const double g1_NC = m_ffs["NC_g1"]->Calc(-q2);
  const double g2_NC = m_ffs["NC_g2"]->Calc(-q2);
  const double g3_NC = m_ffs["NC_g3"]->Calc(-q2);

  const double f1_CC = m_ffs["CC_f1"]->Calc(-q2);
  const double f2_CC = m_ffs["CC_f2"]->Calc(-q2);
  const double f3_CC = m_ffs["CC_f3"]->Calc(-q2);
  const double g1_CC = m_ffs["CC_g1"]->Calc(-q2);
  const double g2_CC = m_ffs["CC_g2"]->Calc(-q2);
  const double g3_CC = m_ffs["CC_g3"]->Calc(-q2);
  
  // double mass_gordon = sqrt((m_massin*m_massin + m_massout*m_massout) / 2.0);

  Vec4C QED_amp, Weak_NC_amp, Weak_CC_amp;
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
      //QED TERM
      if ( fabs(QED_coupling) > 0.0 ) {
        QED_amp = f1_term*f1_QED + f2_term*f2_QED + f3_term*f3_QED - g1_term*g1_QED - g2_term*g2_QED - g3_term*g3_QED;

        if (!m_anti) QED_amp = QED_amp * QED_coupling * QED_prop_factor;
        else         QED_amp = QED_amp * conj(QED_coupling) * QED_prop_factor;
      } 
      //NC TERM
      else if ( fabs(Weak_NC_coupling) > 0.0 ) {
        Weak_NC_amp = f1_term*f1_NC + f2_term*f2_NC + f3_term*f3_NC - g1_term*g1_NC - g2_term*g2_NC - g3_term*g3_NC;

        if (!m_anti) Weak_NC_amp = Weak_NC_amp * Weak_NC_coupling * Weak_NC_prop_factor;
        else         Weak_NC_amp = Weak_NC_amp * conj(Weak_NC_coupling) * Weak_NC_prop_factor;
      }
      //CC TERM
      else if ( fabs(Weak_CC_coupling) > 0.0 ) {
        Weak_CC_amp = f1_term*f1_CC + f2_term*f2_CC + f3_term*f3_CC - g1_term*g1_CC - g2_term*g2_CC - g3_term*g3_CC;

        if (!m_anti) Weak_CC_amp = Weak_CC_amp * Weak_CC_coupling * Weak_CC_prop_factor;
        else         Weak_CC_amp = Weak_CC_amp * conj(Weak_CC_coupling) * Weak_CC_prop_factor;
      }

      // msg_Out() << "Nucleon_Baryon_FFS\n";
      // msg_Out() 
      //   << "Anti?: " << m_anti << "\n"
      //   << "QED: \n     " 
      //   << "Coupling: " << QED_coupling << " \n     "
      //   << "Left: " << QED_cL << " \n     "
      //   << "Right: " << QED_cR << " \n     "
      //   << "f1: " << f1_QED << " \n     "
      //   << "f2: " << f2_QED << " \n     "
      //   << "f3: " << f3_QED << " \n     "
      //   << "g1: " << g1_QED << " \n     "
      //   << "g2: " << g2_QED << " \n     "
      //   << "g3: " << g3_QED << " \n"
      //   << "Weak_NC: \n     " 
      //   << "Coupling: " << Weak_NC_coupling << " \n     "
      //   << "Left: " << Weak_NC_cL << " \n     "
      //   << "Right: " << Weak_NC_cR << " \n     "
      //   << "f1: " << f1_NC << " \n     "
      //   << "f2: " << f2_NC << " \n     "
      //   << "f3: " << f3_NC << " \n     "
      //   << "g1: " << g1_NC << " \n     "
      //   << "g2: " << g2_NC << " \n     "
      //   << "g3: " << g3_NC << " \n"
      //   << "Weak_CC: \n     " 
      //   << "Coupling: " << Weak_CC_coupling << " \n     "
      //   << "Left: " << Weak_CC_cL << " \n     "
      //   << "Right: " << Weak_CC_cR << " \n     "
      //   << "f1: " << f1_CC << " \n     "
      //   << "f2: " << f2_CC << " \n     "
      //   << "f3: " << f3_CC << " \n     "
      //   << "g1: " << g1_CC << " \n     "
      //   << "g2: " << g2_CC << " \n     "
      //   << "g3: " << g3_CC << " \n\n\n";
      // exit(1);

      vector<pair<int,int> > spins;
      spins.push_back(make_pair(pf,hf));
      spins.push_back(make_pair(pi,hi));

      //Old Method QED ONLY
      Insert(QED_amp, spins);

      Insert_ProcessType("QED", QED_amp, spins);
      Insert_ProcessType("Weak_NC", Weak_NC_amp, spins);
      Insert_ProcessType("Weak_CC", Weak_CC_amp, spins);
    }
  }
}