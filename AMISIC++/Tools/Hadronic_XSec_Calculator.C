#include "AMISIC++/Tools/Hadronic_XSec_Calculator.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace AMISIC;
using namespace ATOOLS;
using namespace std;

// will have to make sure that pions are initialised below.
//
// All cross sections here are returned in GeV^{-2}, with s in GeV^2
//
// sigma(s) = X_pomeron s^{0.0808} + X_reggeon s^{-0.4525}
//
// The X_{pomeron/reggeon} (m_xsecpom and m_xsecregge) are given in units of mb,
// to make comparison with literature simpler.
// We use units of mb^{1/2} for the (constant) triple pomeron vertex.
Hadronic_XSec_Calculator::
Hadronic_XSec_Calculator(MODEL::Model_Base * model,
			 const Flavour & fl1,const Flavour & fl2) :
  m_mass_proton(Flavour(kf_p_plus).Mass()),m_mass_proton2(sqr(m_mass_proton)),
  m_mass_pi(Flavour(kf_pi).Mass()), m_mres(2.), m_cres(2.), m_s1(sqr(20.)),
  m_Ypp(-1.), m_c0(2.24), m_c1(2.1), m_testmode(0)
{
  m_flavs[0] = fl1; m_flavs[1] = fl2;
  for (size_t i=0;i<2;i++) {
    m_masses[i] = m_flavs[i].HadMass(); m_masses2[i] = sqr(m_masses[i]);
  }
  m_alphaQED       = (dynamic_cast<MODEL::Running_AlphaQED *>
		      (model->GetScalarFunction("alpha_QED")))->AqedThomson();
  m_eps_pomeron    = (*mipars)("PomeronIntercept");
  m_alphaP_pomeron = (*mipars)("PomeronSlope");
  m_triple_pomeron = (*mipars)("TriplePomeronCoupling");
  m_eta_reggeon    = (*mipars)("ReggeonIntercept");
  m_xsnd_norm      = (*mipars)("SigmaND_Norm");
  m_s0             = 1./m_alphaP_pomeron;

  //////////////////////////////////////////////////////////////////////////////////////////
  // Prefactors, converted to mb, 1/{mb^1/2 GeV^2}, 1/mb for elastic, SD, and DD
  // elastic:            1. / (16 pi)
  // single diffractive:  g_3P    * s_1^{3 eps_P/2}/ (16 pi)
  // double diffractive: (g_3P)^2 * s_1^{2 eps_P/2}/ (16 pi)
  // for single & double diffractive add a factor to compensate the scaling of the
  // pomeron-hadron couplings beta0 with s, from a safe E_cm = 20 GeV down to zero.
  //////////////////////////////////////////////////////////////////////////////////////////
  for (size_t i=0;i<4;i++) m_beta0[i] = sqrt(s_X[i][i]);
  m_prefElastic    = 1.e9/(rpa->Picobarn()*16.*M_PI);
  m_prefSD         = ( m_triple_pomeron * pow(m_s1,3.*m_eps_pomeron/2.)/(16.*M_PI) /
		       (rpa->Picobarn()/1e9) );
  m_prefDD         = ( sqr(m_triple_pomeron) * pow(m_s1,m_eps_pomeron)/(16.*M_PI) /
		       (rpa->Picobarn()/1e9) );
  FixType();
  FixTables();
  if (m_testmode>0) TestXSecs();
}

void Hadronic_XSec_Calculator::FixType() {
  m_type = xsec_type::none;
  if      (m_flavs[0].IsPhoton()  && m_flavs[1].IsPhoton())  m_type = xsec_type::photon_photon;
  else if (m_flavs[0].IsPhoton()  && m_flavs[1].IsNucleon()) m_type = xsec_type::photon_nucleon;
  else if (m_flavs[0].IsNucleon() && m_flavs[1].IsPhoton())  m_type = xsec_type::nucleon_photon;
  else if (m_flavs[0].IsNucleon() && m_flavs[1].IsNucleon()) {
    m_type = xsec_type::nucleon_nucleon;
    if (m_flavs[0].IsAnti() != m_flavs[1].IsAnti()) m_Ypp = 98.39;
  }
  if (m_type==xsec_type::none) THROW(fatal_error,"Unknown type of hadronic cross section.");
}

void Hadronic_XSec_Calculator::TestXSecs() {
  list<double> Es = { 23.5, 62.5, 546., 1800., 16000., 40000. };
  for (size_t i=0;i<2;i++) {
    switch (m_testmode) {
    case 3:
      m_flavs[i] = Flavour(kf_photon);
      m_type     = xsec_type::nucleon_photon;
      break;
    case 2:
      m_flavs[i] = (i==0) ? Flavour(kf_p_plus) : Flavour(kf_photon);
      m_type     = xsec_type::nucleon_photon;
      break;
    case 1:
      m_flavs[i] = Flavour(kf_p_plus);
      m_type     = xsec_type::nucleon_nucleon;
      break;
    default:
      return;
    }
    m_masses[i] = m_flavs[i].HadMass(); m_masses2[i] = sqr(m_masses[i]);
  }
  for (auto E : Es) {
    (*this)(sqr(E));
    Output();
  }
  THROW(normal_exit,"testing complete");
}

void Hadronic_XSec_Calculator::operator()(double s)
{
  m_s     = s;
  m_xstot = m_xsel = m_xssdA = m_xssdB = m_xsdd = 0.;
  ////////////////////////////////////////////////////////////////////////////////////////////
  // All cross sections in mb so far.
  ////////////////////////////////////////////////////////////////////////////////////////////
  switch (m_type) {
  case xsec_type::nucleon_nucleon: CalculateHHXSecs();           break;
  case xsec_type::photon_nucleon:  CalculateHGammaXSecs(0);      break;
  case xsec_type::nucleon_photon:  CalculateHGammaXSecs(1);      break;
  case xsec_type::photon_photon:   CalculatePhotonPhotonXSecs(); break;
  default:
    THROW(fatal_error, "Not yet implemented for unknown type");
  }
  m_xsnd    = m_xstot - m_xsel - m_xssdA - m_xssdB - m_xsdd;
  // convert non-diffractive cross section from millibarn to 1/GeV^2
  m_xsnd   *= 1.e9/rpa->Picobarn();
}

void Hadronic_XSec_Calculator::CalculateHHXSecs() {
  size_t hadtags[2];
  hadtags[0] = hadtags[1] = 0;
  double masses[2];
  masses[0] = m_masses[0];masses[1] = m_masses[1];
  m_xstot   = TotalXSec(hadtags);
  m_xsel    = m_sigmaEl[0][0]  = IntElXSec(hadtags,m_xstot);
  m_xssdA   = m_sigmaSDA[0][0] = IntSDXSec(hadtags,0,masses);
  m_xssdB   = m_sigmaSDB[0][0] = IntSDXSec(hadtags,1,masses);
  m_xsdd    = m_sigmaDD[0][0]  = IntDDXSec(hadtags,masses);
}

void Hadronic_XSec_Calculator::CalculateHGammaXSecs(const size_t photon) {
  size_t hadtags[2];
  hadtags[1-photon] = 0;
  double xstot, prefV, masses[2];
  masses[1-photon] = m_masses[1-photon];
  // Iterate over VMD hadrons and add cross sections
  for (map<Flavour, double>::const_iterator flit=m_fVs.begin();
       flit!=m_fVs.end();flit++) {
    hadtags[photon] = m_indexmap[flit->first];
    masses[photon]  = flit->first.Mass();
    prefV           = m_alphaQED/m_fVs[flit->first];
    m_xstot        += prefV * (xstot = TotalXSec(hadtags));
    m_xsel         += prefV * IntElXSec(hadtags,xstot);
    m_xssdA        += prefV * IntSDXSec(hadtags,0,masses);
    m_xssdB        += prefV * IntSDXSec(hadtags,1,masses);
    m_xsdd         += prefV * IntDDXSec(hadtags,masses);
  }
}

void Hadronic_XSec_Calculator::CalculatePhotonPhotonXSecs() {
  size_t hadtags[2];
  double xstot, prefVV, masses[2];
  ////////////////////////////////////////////////////////////////////////////////////////////
  // Iterate over VMD hadrons and add cross sections
  ////////////////////////////////////////////////////////////////////////////////////////////
  for (map<Flavour, double>::const_iterator flit0=m_fVs.begin();
       flit0!=m_fVs.end();flit0++) {
    hadtags[0] = m_indexmap[flit0->first];
    masses[0]  = flit0->first.Mass();
    for (map<Flavour, double>::const_iterator flit1=m_fVs.begin();
	 flit1!=m_fVs.end();flit1++) {
      hadtags[1] = m_indexmap[flit1->first];
      masses[1]  = flit0->first.Mass();
      prefVV     = sqr(m_alphaQED)/(m_fVs[flit0->first] * m_fVs[flit1->first]);
      m_xstot   += prefVV * (xstot = TotalXSec(hadtags));
      m_xsel    += prefVV * IntElXSec(hadtags,xstot);
      m_xssdA   += prefVV * IntSDXSec(hadtags,0,masses);
      m_xssdB   += prefVV * IntSDXSec(hadtags,1,masses);
      m_xsdd    += prefVV * IntDDXSec(hadtags,masses);
    }
  }
}

double Hadronic_XSec_Calculator::TotalXSec(const size_t hadtags[2]) const {
  ////////////////////////////////////////////////////////////////////////////////////////////
  // Eq.(4) in Schuler and Sjostrand, PRD 49 (Donnachie-Landshof fit)
  ////////////////////////////////////////////////////////////////////////////////////////////
  return ( s_X[hadtags[0]][hadtags[1]]                     * pow(m_s,m_eps_pomeron) +
	   (m_Ypp>0 ? m_Ypp : s_Y[hadtags[0]][hadtags[1]]) * pow(m_s,m_eta_reggeon) );
}

double Hadronic_XSec_Calculator::IntElXSec(const size_t hadtags[2],const double & xstot) const {
  ////////////////////////////////////////////////////////////////////////////////////////////
  // Eq.(7) in Schuler and Sjostrand, PRD 49 (Donnachie-Landshof fit) with elastic slope taken
  // from Eq.(11) ibidem.
  ////////////////////////////////////////////////////////////////////////////////////////////
  double b_elastic =  2.*(s_slopes[hadtags[0]] + s_slopes[hadtags[1]] +
			  m_c0*pow(m_s/m_s0,m_eps_pomeron)-m_c1);
  return m_prefElastic * sqr(xstot) / b_elastic;
}

double Hadronic_XSec_Calculator::IntSDXSec(const size_t hadtags[2],const size_t & diff,
					   const double masses[2]) const {
  ////////////////////////////////////////////////////////////////////////////////////////////
  // Eq.(24) in Schuler and Sjostrand, PRD 49 (Donnachie-Landshof fit) and
  // Eq. (19) in Schuler and Sjostrand Z fuer Physik C 73, with parameters
  // for pp taken from PRD 49 and for VMD states from Table 1 in ZfP C 73.
  // Note: To arrive at the correct prefactor, need to rescale the beta according to Eq.(27).
  ////////////////////////////////////////////////////////////////////////////////////////////
  size_t nodiff = 1-diff;
  double mmin2  = sqr(masses[diff]+2.*m_mass_pi),      mmin  = sqrt(mmin2);
  double mres   = masses[nodiff]-m_mass_proton+m_mres, mres2 = sqr(mres);
  double mmax2  = ( s_c[hadtags[0]][hadtags[1]][0+4*diff]*m_s +
		    s_c[hadtags[0]][hadtags[1]][1+4*diff] );
  if (m_s<=mmin2 || m_s<=mmin*mres) return 0.;
  double bA     = s_slopes[hadtags[nodiff]];
  double B_AX   = ( s_c[hadtags[0]][hadtags[1]][2+4*diff]     +
		    s_c[hadtags[0]][hadtags[1]][3+4*diff]/m_s );
  double J_AX   = Max(0., ( 1./(2.*m_alphaP_pomeron) *
			    log((bA+m_alphaP_pomeron*log(m_s/mmin2))/
				(bA+m_alphaP_pomeron*log(m_s/mmax2))) +
			    m_cres/(2.*(bA+m_alphaP_pomeron*log(m_s/(mmin*mres)))+B_AX) *
			    log(1.+mres2/mmin2) ) );
  return m_prefSD * m_beta0[hadtags[nodiff]] * s_X[hadtags[0]][hadtags[1]] * J_AX;
}

double Hadronic_XSec_Calculator::IntDDXSec(const size_t hadtags[2],
					   const double masses[2]) const {
  ////////////////////////////////////////////////////////////////////////////////////////////
  // Eq.(24) in Schuler and Sjostrand, PRD 49 (Donnachie-Landshof fit), with parameters
  // for pp taken from PRD 49 and for VMD states from Table 1 in Schuler and Sjostrand
  // Z fuer Physik C 73
  ////////////////////////////////////////////////////////////////////////////////////////////
  double logs   = log(m_s),                       log2s  = sqr(logs);
  double s0     = 1./m_alphaP_pomeron,            ss0    = m_s*s0;
  double m1min2 = sqr(masses[0]+2.*m_mass_pi),    m1min  = sqrt(m1min2);
  double m2min2 = sqr(masses[1]+2.*m_mass_pi),    m2min  = sqrt(m2min2);
  double m1res  = masses[0]-m_mass_proton+m_mres, m1res2 = sqr(m1res);
  double m2res  = masses[1]-m_mass_proton+m_mres, m2res2 = sqr(m2res);
  double mmax2  = ( s_d[hadtags[0]][hadtags[1]][3] +
		    s_d[hadtags[0]][hadtags[1]][4]/logs +
		    s_d[hadtags[0]][hadtags[1]][5]/log2s ) * m_s;
  if (m_s    <= sqr(m1min+m2min)   ||
      m1min2 > mmax2               || m2min2 > mmax2||
      ss0    <= m1min2*m2res*m2min || ss0 <= m2min2*m1res*m1min ||
      ss0    <= mmax2*m1res*m1min  || ss0 <= mmax2*m2res*m2min  ||
      ss0    <= m1res*m1min*m2res*m2min ||
      m_s    <= (m1min2*m2min2)/m_mass_proton2) return 0.;
  double arg11  = Max(1.001, ss0/(m1min2*m2res*m2min));
  double arg12  = Max(1.001, ss0/(mmax2*m2res*m2min));
  double arg21  = Max(1.001, ss0/(m2min2*m1res*m1min));
  double arg22  = Max(1.001, ss0/(mmax2*m1res*m1min));
  double Delta0 = ( s_d[hadtags[0]][hadtags[1]][0] +
		    s_d[hadtags[0]][hadtags[1]][1]/logs +
		    s_d[hadtags[0]][hadtags[1]][2]/log2s );
  if (Delta0 <= 0.) return 0.;
  double Bxx    = ( s_d[hadtags[0]][hadtags[1]][6] +
		    s_d[hadtags[0]][hadtags[1]][7]/sqrt(m_s) +
		    s_d[hadtags[0]][hadtags[1]][8]/m_s );
  double Deltay = log((m_s * m_mass_proton2)/(m1min2 * m2min2));
  double J_XX   = Max(0.,
                    ( 1./(2.*m_alphaP_pomeron) * ( Deltay*( log(Deltay/Delta0) - 1.) + Delta0 ) +
                     m_cres/(2.*m_alphaP_pomeron) * (log(1.+m2res2/m2min2) *
                                                             log(log(arg11)/log(arg12)) +
                                                         log(1.+m1res2/m1min2) *
                                                             log(log(arg21)/log(arg22))) +
                     sqr(m_cres)/(2.*m_alphaP_pomeron*log(ss0/(m1res*m1min*m2res*m2min))+Bxx) *
                         log(1.+m1res2/m1min2) * log(1.+m2res2/m2min2) ) );
  return m_prefDD * s_X[hadtags[0]][hadtags[1]] * J_XX;
}

void Hadronic_XSec_Calculator::Output() const {
  msg_Out()<<METHOD<<": Results for "<<m_flavs[0]<<" "<<m_flavs[1]<<" collisions "
	   <<"at E_cm = "<<sqrt(m_s)<<" GeV are {\n"
	   <<"   \\sigma_{tot}   = "<<m_xstot<<" mb\n"
	   <<"   \\sigma_{el}    = "<<m_xsel<<" mb\n"
	   <<"   \\sigma_{sd}(A) = "<<m_xssdA<<" mb\n"
	   <<"   \\sigma_{sd}(B) = "<<m_xssdB<<" mb\n"
	   <<"   \\sigma_{dd}    = "<<m_xsdd<<" mb\n"
	   <<"   \\sigma_{nd}    = "<<m_xsnd/1.e9*rpa->Picobarn()<<" mb = "
	   <<m_xsnd<<" GeV^-2\n}"<<endl;
}


void Hadronic_XSec_Calculator::FixTables() {
  m_hadroncomponents = {
    { Flavour(kf_p_plus),       { Flavour(kf_p_plus) }                          },
    { Flavour(kf_n),            { Flavour(kf_n) }                               },
    { Flavour(kf_photon),       { Flavour(kf_rho_770), Flavour(kf_omega_782),
				  Flavour(kf_phi_1020), Flavour(kf_J_psi_1S) }  }
  };
  ////////////////////////////////////////////////////////////////////////////////////////////
  // Cross section parametrisation taken from Donnachie and Landshoff,
  // and from Schuler and Sjöstrand, Z Phys C 73 677-688 (1997).
  // Following their papers we assume that for pomeron/regeeon fits, there is no difference
  // between nucleons (i.e. we treat protons and neutrons as if they were the same), and between
  // rho(770) and omega(782).
  ////////////////////////////////////////////////////////////////////////////////////////////
  m_indexmap = {
    { Flavour(kf_p_plus),    0 },
    { Flavour(kf_n),         0 },
    { Flavour(kf_rho_770),   1 },
    { Flavour(kf_omega_782), 1 },
    { Flavour(kf_phi_1020),  2 },
    { Flavour(kf_J_psi_1S),  3 }
  };
  ////////////////////////////////////////////////////////////////////////////////////////////
  // Critical values for the total and elastic cross section fit: VMD
  // factors f_V^2/(4 pi) have no units.
  ////////////////////////////////////////////////////////////////////////////////////////////
  m_fVs = {
    { Flavour(kf_rho_770),    2.20 },
    { Flavour(kf_omega_782), 23.60 },
    { Flavour(kf_phi_1020),  18.40 },
    { Flavour(kf_J_psi_1S),  11.50 }
  };
}

////////////////////////////////////////////////////////////////////////////////////////////
// Slopes for the elastic cross section fit in units of GeV^{-2}.
// b slope parameters below eq 19 in Schuler and Sjostrand, Z fuer Physik C 73
////////////////////////////////////////////////////////////////////////////////////////////
double Hadronic_XSec_Calculator::s_slopes[4] = {
  2.30, 1.40, 1.40, 0.23
};
////////////////////////////////////////////////////////////////////////////////////////////
// Critical values for the total/elastic cross section fit in mb.
// For nucleon-nucleon: the Y-values depend on particle/anti-particle, will use the
// particle-particle one here. If necessary, m_Ypp is set for p pbar collisions,
// c.f. the intialization.
// Table 1 in Schuler and Sjostrand, Z fuer Physik C 73
////////////////////////////////////////////////////////////////////////////////////////////
double Hadronic_XSec_Calculator::s_X[4][4] = {
  { 21.700,  13.630,  10.010,  0.970  },
  { 13.630,   8.560,   6.290,  0.609  },
  { 10.010,   6.290,   4.620,  0.447  },
  {  0.970,   0.609,   0.447,  0.0434 }
};

double Hadronic_XSec_Calculator::s_Y[4][4] = {
  { 56.080,  31.790,  -1.510,  -0.146  },
  { 31.790,  13.080,   -0.62,  -0.060  },
  { -1.510,   -0.62,    0.03,  -0.0028 },
  { -0.146,  -0.060, -0.0028,  0.00028 }
};

////////////////////////////////////////////////////////////////////////////////////////////
// The parametrisation of single diffractive events in hadron-hadron collisions.
// First four entries for diffraction of A, second four entries for diffraction of B.
// p p values taken from eq. (26) in Schuler and Sjostrand PRD 49, otherwise
// Table 1 in Schuler and Sjostrand, Z fuer Physik C 73
////////////////////////////////////////////////////////////////////////////////////////////
double Hadronic_XSec_Calculator::s_c[4][4][8] = {             // A    B
  { { 0.213, 0.0, -0.47, 150., 0.213, 0.0, -0.47, 150.},      // p    p
    { 0.213, 0.0, -0.47, 150., 0.267, 0.0, -0.47, 100.},      // p    rho
    { 0.213, 0.0, -0.47, 150., 0.232, 0.0, -0.47, 110.},      // p    phi
    { 0.213, 7.0, -0.55, 800., 0.115, 0.0, -0.47, 110.} },    // p    Jpsi
  { { 0.213, 0.0, -0.47, 150., 0.267, 0.0, -0.47, 100.},      // rho  p
    { 0.267, 0.0, -0.46,  75., 0.267, 0.0, -0.46,  75.},      // rho  rho
    { 0.267, 0.0, -0.48, 100., 0.232, 0.0, -0.46,  85.},      // rho  phi
    { 0.267, 6.0, -0.56, 420., 0.115, 0.0, -0.50,  90.} },    // rho  Jpsi
  { { 0.213, 0.0, -0.47, 150., 0.232, 0.0, -0.47, 110.},      // phi  p
    { 0.267, 0.0, -0.48, 100., 0.232, 0.0, -0.46,  85.},      // phi  rho
    { 0.232, 0.0, -0.48, 110., 0.232, 0.0, -0.48, 110.},      // phi  phi
    { 0.232, 6.0, -0.56, 470., 0.115, 0.0, -0.52, 120.} },    // phi  Jpsi
  { { 0.213, 7.0, -0.55, 800., 0.115, 0.0, -0.47, 110.},      // Jpsi p
    { 0.267, 6.0, -0.56, 420., 0.115, 0.0, -0.50,  90.},      // Jpsi rho
    { 0.232, 6.0, -0.56, 470., 0.115, 0.0, -0.52, 120.},      // Jpsi phi
    { 0.115, 5.5, -0.58, 570., 0.115, 5.5, -0.58, 570.} }     // Jpsi Jpsi
};

double Hadronic_XSec_Calculator::s_d[4][4][9] = {                       // A    B
  { { 3.20,  -9.00, 17.40, 0.070, -0.44, 1.36, -1.05,  40.0, 8000. },   // p    p
    { 3.11,  -7.10, 10.60, 0.073, -0.41, 1.17, -1.41,  31.6,   95. },   // p    rho
    { 3.12,  -7.43,  9.21, 0.067, -0.44, 1.41, -1.35,  36.5,  132. },   // p    phi
    { 3.13,  -8.18, -4.20, 0.056, -0.71, 3.12, -1.12,  55.2, 1298. } }, // p    Jpsi
  { { 3.11,  -7.10, 10.60, 0.073, -0.41, 1.17, -1.41,  31.6,   95. },   // rho  p
    { 3.11,  -6.90, 11.40, 0.078, -0.40, 1.05, -1.40,  28.4,   78. },   // rho  rho
    { 3.11,  -7.13, 10.00, 0.071, -0.41, 1.23, -1.34,  33.1,  105. },   // rho  phi
    { 3.12,  -7.90, -1.49, 0.054, -0.64, 2.72, -1.13,  53.1,  995. } }, // rho  Jpsi
  { { 3.12,  -7.43,  9.21, 0.067, -0.44, 1.41, -1.35,  36.5,  132. },   // phi  p
    { 3.11,  -7.13, 10.00, 0.071, -0.41, 1.23, -1.34,  33.1,  105. },   // phi  rho
    { 3.11,  -7.39,  8.22, 0.065, -0.44, 1.45, -1.36,  38.1,  148. },   // phi  phi
    { 3.18,  -8.95, -3.37, 0.057, -0.76, 3.32, -1.12,  55.6, 1472. } }, // phi  Jpsi
  { { 3.13,  -8.18, -4.20, 0.056, -0.71, 3.12, -1.12,  55.2, 1298. },   // Jpsi p
    { 3.12,  -7.90, -1.49, 0.054, -0.64, 2.72, -1.13,  53.1,  995. },   // Jpsi rho
    { 3.18,  -8.95, -3.37, 0.057, -0.76, 3.32, -1.12,  55.6, 1472. },   // Jpsi phi
    { 4.18, -29.20, 56.20, 0.074, -1.36, 6.67, -1.14, 116.2, 6532. } }  // Jpsi Jpsi
};

