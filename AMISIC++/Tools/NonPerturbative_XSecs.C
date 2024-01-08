#include "AMISIC++/Tools/NonPerturbative_XSecs.H"
#include "BEAM/Main/Beam_Base.H"
#include "REMNANTS/Main/Remnant_Handler.H"
#include "ATOOLS/Phys/Weights.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace AMISIC;
using namespace ATOOLS;
using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////
// Equations below refer mainly to Schuler and Sjostrand, PRD 49 (1994) 2257
///////////////////////////////////////////////////////////////////////////////////////////

NonPerturbative_XSecs::
NonPerturbative_XSecs(REMNANTS::Remnant_Handler * remnants,Hadronic_XSec_Calculator * xsecs) :
  p_remnants(remnants), p_xsecs(xsecs), m_variableS(false),
  m_evttype(mipars->GetEvtType()),
  m_inflav(p_xsecs->GetFlavs()), 
  m_smin(p_xsecs->Smin()), 
  m_eps_pomeron(p_xsecs->EpsPomeron()), m_alphaP_pomeron(p_xsecs->AlphaPPomeron()),
  m_triple_pomeron(p_xsecs->TriplePomeron()), m_alphaQED(p_xsecs->AlphaQED()),
  m_s0(1./m_alphaP_pomeron), 
  m_mpi(Flavour(kf_pi).Mass()), m_mpi2(m_mpi*m_mpi),
  m_deltaMres(p_xsecs->Diffractive_Mres()), m_cres(p_xsecs->Diffractive_cres()),
  m_mrho(Flavour(kf_rho_770).Mass()), m_mrho2(sqr(m_mrho)),
  m_mrho_min(0.3), m_q_rho(sqrt(m_mrho2-4.*m_mpi2)), 
  m_momega(Flavour(kf_omega_782).Mass()), m_momega2(sqr(m_momega)),
  m_q_omega(sqrt(m_momega2-4.*m_mpi2)), 
  m_Grho(Flavour(kf_rho_770).Width()), m_Grho2(sqr(m_Grho)),
  m_Gomega(Flavour(kf_omega_782).Width()), m_Gomega2(sqr(m_Gomega)), m_A2max(0.), 
  m_twopions(two_pions::none),
  m_calls(0), m_fails(0),
  m_ana(false)
{
  m_twopions   = mipars->GetTwoPionTreatment();
  m_f_omega    = (*mipars)("f_omega"); 
  m_phi_omega  = (*mipars)("phi_omega");
  m_f_nr       = (*mipars)("f_nr");
  m_Lambda2_nr = sqr((*mipars)("Lambda_nr"));
  m_delta_nr   = (*mipars)("delta_nr");
  for (size_t i=0;i<1000;i++) {
    double A2test = RhoMassModifier(sqr(m_mrho_min+double(i)/100)); 
    if (A2test>m_A2max) m_A2max = A2test;
  }
  for (size_t i=0;i<2;i++) p_beams[i] = NULL;
  for (size_t i=0;i<2;i++) {
    m_masses[i] = m_inflav[i].HadMass(); m_masses2[i] = sqr(m_masses[i]);
  }
  if (m_ana) { Tests(); exit(0); }
  //msg_Out()<<METHOD<<": f = "<<m_f_nr<<", Lambda^2 = "<<m_Lambda2_nr<<", delta = "<<m_delta_nr<<", "
  //	   <<"m(pi)^2 = "<<m_mpi2<<", Gamma(omega) = "<<m_Gomega<<".\n"; exit(1);
}

NonPerturbative_XSecs::~NonPerturbative_XSecs() {
  if (m_fails>0)
    msg_Out()<<METHOD<<" with "<<m_fails<<" fails in "<<m_calls<<" calls.\n";
}

void NonPerturbative_XSecs::SetBeams(BEAM::Beam_Base * beam1,BEAM::Beam_Base * beam2) {
  /////////////////////////////////////////////////////////////////////////////////////////
  // This essentially fixes whether we have to select a c.m. frame of the collision (in
  // the case of EPA, i.e. photon-photon or photon-nucleon collisions) or not (for pure
  // hadronic initial states, i.e. nucleon-nucleon collisions).
  /////////////////////////////////////////////////////////////////////////////////////////
  p_beams[0] = beam1; p_beams[1] = beam2;
  for (size_t i=0;i<2;i++) {
    if (p_beams[i]!=NULL && p_beams[i]->Type()==BEAM::beamspectrum::EPA) {
      m_variableS = true;
      break;
    }
  }
}

void NonPerturbative_XSecs::CalculateSDependentCrossSections() {
  /////////////////////////////////////////////////////////////////////////////////////////
  // Pre-calculating the s-dependent cross sections by folding the EPA spectrum with the
  // non-perturbative cross sections enoced in the p_xsecs, the Hadronic_XSec_Calculator
  // class.  The initial-state integration is performed within the Beam_Integrator class.
  /////////////////////////////////////////////////////////////////////////////////////////
  if (!m_variableS) return;	
  m_integrator.Init(p_xsecs,p_beams[0],p_beams[1],m_evttype);
}

Blob * NonPerturbative_XSecs::MakeScatter() {
  /////////////////////////////////////////////////////////////////////////////////////////
  // Creating an initial state kinematics, i.e. its c.m. energy squared s and its
  // rapidity y - the latter also fixing the boost from the c.m. to the lab frame.  In
  // the second step the mode (elastic/diffractive) is fixed and the final state flavours
  // are selected.  For initial photons these can be either of the four lowest vector
  // mesons: rho (770), omega (782), phi (1020), or J/Psi.  Then the final state kinematics
  // is created according to the differential cross sections of the Schuler-Sjostrand
  // model (relevant equations are documented in the respective methods).
  /////////////////////////////////////////////////////////////////////////////////////////
  m_calls++;
  if (m_variableS) {
    if (!m_integrator()) {
      msg_Error()<<"Error in "<<METHOD<<": "
		 <<"could not create initial-state phase-space.\n";
      return NULL;
    }
  }
  for (size_t beam=0;beam<2;beam++)
    m_inmom[beam] = p_remnants->GetRemnant(beam)->InMomentum();
  m_s     = (m_inmom[0]+m_inmom[1]).Abs2();
  m_boost = Poincare(m_inmom[0]+m_inmom[1]);
  for (size_t beam=0;beam<2;beam++) m_boost.Boost(m_inmom[beam]);
  if (m_s<m_smin) THROW(fatal_error,"Insufficient s = "+ToString(m_s)+
			" vs "+ToString(m_smin)+
			" in Non-Perturbative Events in AMISIC++");
  (*p_xsecs)(m_s);
  array<Flavour, 2> flavs;
  switch(SelectMode()) {
  case event_mode::elastic:
    if (!p_xsecs->SelectEl(flavs) || !FixFS(flavs))
      THROW(fatal_error,"Could not define elastic FS");
    return ElasticScatter();       
  case event_mode::SDA:
    if (!p_xsecs->SelectSDA(flavs) || !FixFS(flavs))
      THROW(fatal_error,"Could not define SD(A) FS");
    return SingleDiffractiveScatter(0);
  case event_mode::SDB:
    if (!p_xsecs->SelectSDB(flavs) || !FixFS(flavs))
      THROW(fatal_error,"Could not define SD(B) FS");
    return SingleDiffractiveScatter(1);
  case event_mode::DD:
    if (!p_xsecs->SelectDD(flavs) || !FixFS(flavs))
      THROW(fatal_error,"Could not define DD FS");
    return DoubleDiffractiveScatter();
  case event_mode::unknown:
  default:                  break;
  }
  THROW(fatal_error,"Unknown Event Mode in Non-Perturbative Events in AMISIC++");
  return NULL;
}

const event_mode::code NonPerturbative_XSecs::SelectMode() {
  /////////////////////////////////////////////////////////////////////////////////////////
  // Trivial initialisation of outgoing flavours as being either "none" or the
  // initial flavour.  In principle we can have up to four outgoing particles in a
  // two-body scatter, when each of the initial particles decomposes diffractively into
  // two constituents.  This means we map the initial i = {0,1} onto i -> {2i, 2i+1}. 
  /////////////////////////////////////////////////////////////////////////////////////////
  m_weight = 0.;
  for (size_t i=0;i<2;i++) {
    m_outflav[2*i]     = m_inflav[i];
    m_outflav[2*i+1]   = Flavour(kf_none);
    m_outmasses[2*i]   = m_masses[i];
    m_outmasses2[2*i]  = m_masses2[i];
    m_outmasses[2*i+1] = m_outmasses2[2*i+1] = 0.;
  }
  double xsecs[4]; for (size_t i=0;i<3;i++) xsecs[i] = 0.;

  if (m_evttype==evt_type::Elastic || m_evttype==evt_type::QuasiElastic)
    m_weight += xsecs[0] = p_xsecs->SigmaEl();
  if (m_evttype==evt_type::DiffractiveA || m_evttype==evt_type::QuasiElastic)
    m_weight += xsecs[1] = p_xsecs->SigmaSDA();
  if (m_evttype==evt_type::DiffractiveB || m_evttype==evt_type::QuasiElastic)
    m_weight += xsecs[2] = p_xsecs->SigmaSDB();
  if (m_evttype==evt_type::DiffractiveAB || m_evttype==evt_type::QuasiElastic)
    m_weight += xsecs[3] = p_xsecs->SigmaDD();
  double total = m_weight * ran->Get() * 0.999999999;
  size_t i = 0;
  for (i=3;i>0;i--) {
    total -= xsecs[i];
    if (total<=0.) break;
  }
  switch (i) {
  case 3: return event_mode::DD;
  case 2: return event_mode::SDB;
  case 1: return event_mode::SDA;
  case 0: return event_mode::elastic;
  default:
    break;
  }
  THROW(fatal_error,"Could not define event_mode and/or outgoing flavours.")
}

bool NonPerturbative_XSecs::FixFS(array<ATOOLS::Flavour, 2> & flavs) {
  /////////////////////////////////////////////////////////////////////////////////////////
  // This method is mainly to meant to initialise the outgoing flavours and masses for 
  // initial photons, which will turn into vector mesons through the VMD model, and to
  // define the correct associated hadron tags for the various pomeron couplings.
  /////////////////////////////////////////////////////////////////////////////////////////
  for (size_t i=0;i<2;i++) {
    m_outflav[2*i]    = flavs[i];
    m_outmasses[2*i]  = m_outflav[2*i].HadMass();
    m_outmasses2[2*i] = sqr(m_outmasses[2*i]);
    m_hadtags[i]      = p_xsecs->Index(flavs[i],i);
    if (m_hadtags[i]>99) return false;
  }
  return true;
}

double NonPerturbative_XSecs::DiffElXSec(const double & s,const double & t) const
{
  /////////////////////////////////////////////////////////////////////////////////////////
  // Differential elastic scatter dsigma_el/dt as given in (7).  We consciously ignore the
  // effect of Reggeon exchange that would be relevant for lower c.m. energies, and operate
  // solely in the high-energy limit.  The couplings of the pomeron to the hadrons, beta, are
  // encoded through intercept and slope, and realised as inline funtions in the .H file.
  // We only use this method for testing purposes.
  /////////////////////////////////////////////////////////////////////////////////////////
  return ( 1./(16.*M_PI) * sqr(beta(m_hadtags[0],t) * beta(m_hadtags[1],t)) *
	   exp((2.*m_eps_pomeron + 2.*m_alphaP_pomeron*t) * log(s/m_s0)) ); 
}

Blob * NonPerturbative_XSecs::ElasticScatter() {
  /////////////////////////////////////////////////////////////////////////////////////////
  // MC realisation of selecting the 4-momentum transfer t in (quasi-)elastic scatters,
  // again in the high-energy limit also used in DiffElXSec.  Once t is chosen, we create outgoing
  // momenta in FixOutMomenta() and fill the scatter blob.
  // In the Tests method we explicitly verify that the MC reproduces the differential
  // cross section above in DiffElXSec.
  /////////////////////////////////////////////////////////////////////////////////////////
  double arg  = 2.*(m_eps_pomeron + m_alphaP_pomeron * log(m_s/m_s0) +
		    p_xsecs->s_slopes[m_hadtags[0]] + p_xsecs->s_slopes[m_hadtags[1]]);
  double t;
  size_t trials=0;
  do { t = ExponentialDist(-m_s,0.,arg); } while (!FixOutMomenta(t) && (trials++)<1000);
  if (trials>=1000) {
    if (m_fails<5)
      msg_Error()<<METHOD<<" fails: t = "<<t<<" yields no momenta, s = "<<m_s<<".\n";
    m_fails++;
    return NULL;
  }
  Blob * blob   = InitBlob();
  return blob;
}

double NonPerturbative_XSecs::DiffSDXSec(const double & s, const double & t,
					 const double & M2,const size_t & pos) const
{
  /////////////////////////////////////////////////////////////////////////////////////////
  // Differential single diffractive scatter dsigma_el/dt dM^2 as given in (12), again
  // in the high-energy limit.   In the spirit of Schuler-Sjostrand we modify this however with
  // the F_SD function of (22), to arrive at the differential form underpinning (24) and (25),
  // but without the slope modifications of (20).
  // F_SD has been realised as inline function in the .H file.
  // Again, we only use this method for testing purposes.
  /////////////////////////////////////////////////////////////////////////////////////////
  double m     = GetFlavour(pos,m_hadtags[pos]).Mass(), m2 = sqr(m);
  double M2res = sqr(m+m_deltaMres), M2min = sqr(m+2.*m_mpi);
  if (M2<M2min) return 0.;
  return ( F_SD(M2,M2res) / (16.*M_PI*M2) *
	   m_triple_pomeron * sqr(beta(m_hadtags[pos],t)) * beta(m_hadtags[1-pos],0.) *
	   exp((2.*m_eps_pomeron + 2.*m_alphaP_pomeron*t) * log(s/M2)) ); 
}

ATOOLS::Blob * NonPerturbative_XSecs::SingleDiffractiveScatter(const size_t & pos) {
  /////////////////////////////////////////////////////////////////////////////////////////
  // MC realisation of selecting the 4-momentum transfer t and the diffractive mass M in single
  // diffractive scattering of particle at pos = {0,1}, again in the high-energy limit.  This
  // is the MC realisation of the DiffSDXSec above.
  // After t and M are fixed and the outgoing momenta are constructed in FixOutMomenta, the
  // diffracted hadron is split into constitutents with SplitDiffractiveState(pos) and the
  // scatter blob is filled.
  // In the Tests method we explicitly verify that the MC reproduces the differential
  // cross section above in DiffSDXSec.
  /////////////////////////////////////////////////////////////////////////////////////////
  double M2min  = sqr(m_outmasses2[2*pos]+2.*m_mpi), M2res = sqr(m_outmasses2[2*pos]+m_deltaMres);
  double argt   = 2.*(p_xsecs->s_slopes[m_hadtags[1-pos]] + m_alphaP_pomeron*log(m_s/M2min));
  double argM   = -(1. + 2.*m_eps_pomeron);
  double maxval = F_SD(M2min,M2res);
  double t, M2, value;
  size_t attempts = 0;
   do {
    if (attempts++>10000) return NULL;
    do {
      t      = ExponentialDist(-m_s,0,argt);
      M2     = PowerDist(M2min,m_s/4.,argM);
      value  = ( F_SD(M2,M2res) *
		 exp(2.*(p_xsecs->s_slopes[m_hadtags[1-pos]]+m_alphaP_pomeron*log(m_s/M2))*t) /
		 exp(argt * t) );
    } while (value/maxval<ran->Get());
    m_outmasses2[2*pos]     = M2;
  } while (!FixOutMomenta(t));
  SplitDiffractiveState(pos);
  Blob * blob = InitBlob(M2,M2);
  return blob;
}

double NonPerturbative_XSecs::DiffDDXSec(const double & s,const double & t,
					 const std::array<double, 2> & M2) const
{
  /////////////////////////////////////////////////////////////////////////////////////////
  // Differential double diffractive scatter dsigma_el/dt dM_1^2 dM_2^2 as given in (13), again
  // in the high-energy limit.   In the spirit of Schuler-Sjostrand we modify this however with
  // the F_DD function of (22), to arrive at the differential form underpinning (24) and (25),
  // but again without the slope modifications of (23).
  // F_DD has been realised as inline function in the .H file.
  // Again, we only use this method for testing purposes.
  /////////////////////////////////////////////////////////////////////////////////////////
  array<double ,2> m, m2, M2res, M2min;
  for (size_t i=0;i<2;i++) {
    m[i]     = GetFlavour(i,m_hadtags[i]).Mass();
    m2[i]    = sqr(m[i]);
    M2res[i] = sqr(m[i]+m_deltaMres);
    M2min[i] = sqr(m[i]+2.*m_mpi);
    if (M2[i]<M2min[i]) return 0.;
  }
  return ( F_DD(M2,M2res) / (16.*M_PI*M2[0]*M2[1]) *
	   sqr(m_triple_pomeron) * beta(m_hadtags[0],0.) * beta(m_hadtags[1],0.) *
	   exp((2.*m_eps_pomeron + 2.*m_alphaP_pomeron*t) * log((s*m_s0)/(M2[0]*M2[1]))) ); 
}

ATOOLS::Blob * NonPerturbative_XSecs::DoubleDiffractiveScatter() {
  /////////////////////////////////////////////////////////////////////////////////////////
  // MC realisation of selecting the 4-momentum transfer t and the diffractive masses M_1 and
  // M_2 in double diffractive scattering of particle, again in the high-energy limit.  This
  // is the MC realisation of the DiffDDXSec above.
  // After t, M_1 and M_2 are fixed and the outgoing momenta are constructed in FixOutMomenta, the
  // diffracted hadron is split into constitutents with SplitDiffractiveState() and the
  // scatter blob is filled.
  // In the Tests method we explicitly verify that the MC reproduces the differential
  // cross section above in DiffDDXSec.
  /////////////////////////////////////////////////////////////////////////////////////////
  array<double, 2> M2min, M2res, M2;
  for (size_t i=0;i<2;i++) {
    M2min[i] = sqr(m_outmasses2[2*i]+2.*m_mpi);
    M2res[i] = sqr(m_outmasses2[2*i]+m_deltaMres);
  }
  double argt   = 2.*m_alphaP_pomeron * log((m_s*m_s0)/(M2min[0]*M2min[1]));
  double argM   = -(1. + 2.*m_eps_pomeron);
  double maxval = F_DD(M2min,M2res);
  double t, value;
  size_t attempts = 0;
  do {
    if (attempts++>10000) return NULL;
    do {
      t       = ExponentialDist(-m_s,0,argt);
      for (size_t i=0;i<2;i++) M2[i] = PowerDist(M2min[i],m_s,argM);
      value   = ( F_DD(M2, M2res) *
		  exp(2.*m_alphaP_pomeron * log((m_s*m_s0)/(M2[0]*M2[1])) * t) /
		  exp(argt*t) );
    } while (value/maxval<ran->Get());
    for (size_t i=0;i<2;i++) m_outmasses2[2*i] = M2[i];
  } while (!FixOutMomenta(t));
  for (size_t i=0;i<2;i++) SplitDiffractiveState(i);
  Blob * blob = InitBlob(max(M2[0],M2[1]),max(M2[0],M2[1]));
  return blob;
}

bool NonPerturbative_XSecs::FixOutMomenta(const double & t) {
  /////////////////////////////////////////////////////////////////////////////////////////
  // Fixing the outgoing momenta of the 2->2 scatter, given the masses of the outgoing,
  // possibly diffracted, systems.  As we map the incoming particles i = {0,1} like
  // i -> {2i, 2i+1}, we only have to create outmomenta m_outmom[0] and m_outmom[2].
  // This method also includes a boost back into the lab-frame.
  /////////////////////////////////////////////////////////////////////////////////////////
  if (m_s<sqr(sqrt(m_outmasses2[0])+sqrt(m_outmasses2[2]))) return false;
  if ((m_outflav[0].Kfcode()==kf_rho_770 || m_outflav[2].Kfcode()==kf_rho_770) &&
      !SetRhoMasses2()) return false;
  double p22  = ( (sqr(m_s-m_outmasses2[0]-m_outmasses2[2]) - 4.*m_outmasses2[0]*m_outmasses2[2] )/
		  (4.*m_s) );
  double p2   = sqrt(p22), E[4];
  for (size_t i=0;i<2;i++) E[2*i] = sqrt(p22+m_outmasses2[2*i]);
  double cost = (t-m_masses2[0]-m_outmasses2[0] + 2.*m_inmom[0][0]*E[0])/(2.*m_inmom[0][3]*p2);
  if (dabs(cost)>1.) return false;
  double sint = sqrt(1.-cost*cost);
  double phi  = 2.*M_PI*ran->Get();
  m_outmom[0] = Vec4D(E[0],  p2*sint*cos(phi),  p2*sint*sin(phi),  p2*cost);
  m_outmom[2] = Vec4D(E[2], -p2*sint*cos(phi), -p2*sint*sin(phi), -p2*cost);
  for (size_t beam=0;beam<2;beam++) {
    m_boost.BoostBack(m_inmom[beam]);
    m_boost.BoostBack(m_outmom[2*beam]);
  }
  for (size_t i=0;i<2;i++)
    if (m_outflav[2*i].Kfcode()==kf_rho_770) SplitRhoIntoPions(2*i);
  return true;
}

void NonPerturbative_XSecs::SplitRhoIntoPions(const size_t & pos) {
  /////////////////////////////////////////////////////////////////////////////////////////
  // Decaying a potentially off-shell rho isotropically into two pions.
  /////////////////////////////////////////////////////////////////////////////////////////
  if (m_twopions==two_pions::none || m_outflav[pos].Kfcode()!=kf_rho_770) return;
  Vec4D  rhomom    = m_outmom[2*pos];
  double rhoM2     = m_outmasses2[2*pos];
  double p2        = rhoM2/4.-m_mpi2,  E = sqrt(p2+m_mpi2), p = sqrt(p2);
  double costh     = 1.-2.*ran->Get(),   sinth = sqrt(1.-costh*costh);
  double phi       = 2.*M_PI*ran->Get(), cosph = cos(phi), sinph = sin(phi); 
  m_outmom[pos]    = Vec4D(E, p*sinth*cosph, p*sinth*sinph, p*costh);
  m_outmom[pos+1]  = Vec4D(E,-p*sinth*cosph,-p*sinth*sinph,-p*costh);
  m_outflav[pos]   = Flavour(kf_pi_plus);
  m_outflav[pos+1] = Flavour(kf_pi_plus).Bar();
  for (size_t i=0;i<2;i++) m_outmasses[pos+i] = m_outflav[pos+i].HadMass();
  Poincare boost(rhomom);
  boost.BoostBack(m_outmom[pos]);   
  boost.BoostBack(m_outmom[pos+1]);  
}

bool NonPerturbative_XSecs::SetRhoMasses2() {
  /////////////////////////////////////////////////////////////////////////////////////////
  // If there are rho mesons in the final state and we apply any form of in-class distribution
  // of virtual rho masses in their decay into two pions, we fix the rho mass here
  /////////////////////////////////////////////////////////////////////////////////////////
  if (m_twopions==two_pions::none) return true;
  size_t trials=0;
  do {
    for (size_t i=0;i<2;i++) {
      if (m_outflav[2*i]!=Flavour(kf_rho_770)) continue;
      do {
	m_outmasses2[2*i] = sqr(m_mrho_min + (5.-m_mrho_min)*ran->Get());
      } while(RhoMassModifier(m_outmasses2[2*i]) < m_A2max*ran->Get());
    }
  } while (trials++ < int(m_A2max)*10000 &&
	   m_s<sqr(sqrt(m_outmasses2[0])+sqrt(m_outmasses2[2])));
  return (trials<int(m_A2max)*10000);
}

double NonPerturbative_XSecs::RhoMassModifier(const double & M2) {
  /////////////////////////////////////////////////////////////////////////////////////////
  // Reweighting of the original Breit-Wigner Mass distribution with fixed rho-width to
  // the actual mass distribution/
  /////////////////////////////////////////////////////////////////////////////////////////
  double A_nr = ( (m_twopions==two_pions::cont_only || m_twopions==two_pions::rho_omega_cont) ?
		  m_f_nr/pow(M2-4.*m_mpi2+m_Lambda2_nr, m_delta_nr) : 0. );
  /////////////////////////////////////////////////////////////////////////////////////////
  // correct for the initial fixed-width Breit-Wigner when only using the continuum
  /////////////////////////////////////////////////////////////////////////////////////////
  if (m_twopions==two_pions::cont_only)
    return sqr(A_nr); // * (sqr(M2-m_mrho2)+m_mrho2*m_Grho2)/(m_mrho2*m_Grho2);
  /////////////////////////////////////////////////////////////////////////////////////////
  // Including a more complicated structure
  /////////////////////////////////////////////////////////////////////////////////////////
  double  qrhoratio = pow(sqrt(M2-4.*m_mpi2)/m_q_rho,3.);
  double  GrhoM2    = m_Grho*qrhoratio*m_mrho/sqrt(M2), Grho2M2 = GrhoM2*GrhoM2;
  Complex BWrho     = (m_mrho*m_Grho/(sqr(M2-m_mrho2)+m_mrho2*Grho2M2) *
		       Complex(M2-m_mrho2, -m_mrho*GrhoM2) );
  Complex A         = Complex(1.,0.);
  /////////////////////////////////////////////////////////////////////////////////////////
  // adding the omega-rho-interference, and, potentially, the continuum
  /////////////////////////////////////////////////////////////////////////////////////////
  if (m_twopions==two_pions::rho_omega || m_twopions==two_pions::rho_omega_cont) {
    double  qomegaratio = pow(sqrt(M2-4.*m_mpi2)/m_q_omega,3.);
    double  GomegaM2    = m_Gomega*qomegaratio*m_momega/sqrt(M2), Gomega2M2 = GomegaM2*GomegaM2;
    Complex BWomega     = ( m_momega*m_Gomega/(sqr(M2-m_momega2)+m_momega2*Gomega2M2) *
			    Complex(M2-m_momega2, -m_momega*GomegaM2) );
    A += m_f_omega * M2/m_momega2 * Complex(cos(m_phi_omega),sin(m_phi_omega)) * BWomega;
  }
  /////////////////////////////////////////////////////////////////////////////////////////
  // Correct for the initial fixed-width Breit-Wigner when only using the rho:
  // this is mainly the effect of running vs. fixed width - the fixed-width is part of
  // the selection and therefore the weight must be divided out.
  /////////////////////////////////////////////////////////////////////////////////////////
  double rel = (sqr(std::abs(BWrho * A + Complex(A_nr,0.)))/sqr(std::abs(BWrho * A)));
  if (rel<1.) 
    msg_Out()<<"Gotcha! "<<METHOD<<"(M = "<<sqrt(M2)<<", A = "<<(BWrho*A)<<", "<<A_nr
	     <<" --> wt = "<<sqr(std::abs(BWrho * A + Complex(A_nr, 0.)))<<" ["<<rel<<"].\n";
  // *
  //			      (sqr(M2-m_mrho2) + m_mrho2*m_Grho2)/(m_mrho2*m_Grho2))<<"\n";
  return sqr(std::abs(BWrho * A + Complex(A_nr, 0.) ));
  // * (sqr(M2-m_mrho2) + m_mrho2*m_Grho2)/(m_mrho2*m_Grho2);
}

bool NonPerturbative_XSecs::SplitDiffractiveState(const size_t & pos) {
  /////////////////////////////////////////////////////////////////////////////////////////
  // Given the momentum and mass of the diffractive system, we select its (constituent)
  // flavours in SelectFlavoursOfDiffraction(pos,M2).  As an additional safety measure we
  // check if the diffractive mass is larger than the sum of the constituent masses - if this
  // is not the case (which may happen for nucleons) we explicitly split the excited system
  // into a nucleon and a pion, in SplitDiffractiveStateIntoHadrons(pos,M2).
  // Once we have create two outgoing particles we let the diffractive system decay
  // isotropically in its own rest frame.
  // This method also includes a boost back of the constituent momenta into the lab-frame.
  /////////////////////////////////////////////////////////////////////////////////////////
  Vec4D  diffmom  = m_outmom[2*pos];
  double M2       = diffmom.Abs2();
  if (!SelectFlavoursOfDiffraction(pos,M2)) {
    if (!SplitDiffractiveStateIntoHadrons(pos,M2)) exit(1);
  }
  for (size_t i=0;i<2;i++) {
    m_outmasses[2*pos+i]  = m_outflav[2*pos+i].HadMass();
    m_outmasses2[2*pos+i] = sqr(m_outmasses[2*pos+i]); 
  }
  double p2   = ( (sqr(M2-m_outmasses2[2*pos]-m_outmasses2[2*pos+1]) -
		   4.*m_outmasses2[2*pos]*m_outmasses2[2*pos+1]) / (4.*M2) ), p = sqrt(p2);
  double E[2];
  for (size_t i=0;i<2;i++) E[i] = sqrt(p2+m_outmasses2[2*pos+i]);
  double cost = -1.+2.*ran->Get(), sint = sqrt(1.-cost*cost), phi = 2.*M_PI*ran->Get();
  m_outmom[2*pos]   = Vec4D(E[0],p*sint*cos(phi),p*sint*sin(phi),p*cost);
  m_outmom[2*pos+1] = Vec4D(E[1],-Vec3D(m_outmom[2*pos]));
  Poincare boost(diffmom);
  boost.BoostBack(m_outmom[2*pos]);   
  boost.BoostBack(m_outmom[2*pos+1]);
  return true;
}

bool NonPerturbative_XSecs::
SelectFlavoursOfDiffraction(const size_t & pos,const double & M2) {
  /////////////////////////////////////////////////////////////////////////////////////////
  // This method splits a diffractive sysmte of mass M into its (di-)quark constituents:
  // - rho / omega mesons are split with equal probability either into u-ubar or d-dbar pairs
  // - phi mesons are split into s-sbar pairs
  // - J/psi's are split into c-cbar pairs
  // - protons are split into
  //   d-uu_1 with probability P = 1/3, u-ud_1 with P = 1/6, and into u-ud_0 with P = 1/2
  //   provided that the diffractive mass is large enough.  If this is not the case, the
  //   rescue system will be invoked.
  // TODO: Also add neutrons
  /////////////////////////////////////////////////////////////////////////////////////////
  bool    anti = m_outflav[2*pos].IsAnti();
  kf_code diff = m_outflav[2*pos].Kfcode();
  if (diff==2212) {
    double random = ran->Get();
    if (random<1./3. && M2>sqr(Flavour(kf_d).HadMass()+Flavour(kf_uu_1).HadMass())) {
      m_outflav[2*pos]=Flavour(kf_d); m_outflav[2*pos+1]=Flavour(kf_uu_1);
    }
    else if (random<1./2. && M2>sqr(Flavour(kf_u).HadMass()+Flavour(kf_ud_1).HadMass())) {
      m_outflav[2*pos]=Flavour(kf_u); m_outflav[2*pos+1]=Flavour(kf_ud_1);
    }
    else if (M2>sqr(Flavour(kf_u).HadMass()+Flavour(kf_ud_0).HadMass())) {
      m_outflav[2*pos]=Flavour(kf_u); m_outflav[2*pos+1]=Flavour(kf_ud_0);
    }
    else return false;
  }
  else if ((diff==113 || diff==223) &&
	   M2>sqr(2.*Max(Flavour(kf_u).HadMass(),Flavour(kf_u).HadMass())) ) {
    m_outflav[2*pos]   = ran->Get()>0.5 ? Flavour(kf_u) : Flavour(kf_d);
    m_outflav[2*pos+1] = m_outflav[2*pos].Bar();
  }
  else if (diff==333 && M2>sqr(2.*Flavour(kf_s).HadMass())) {
    m_outflav[2*pos]   = Flavour(kf_s);
    m_outflav[2*pos+1] = m_outflav[2*pos].Bar();
  }
  else if (diff==433 && M2>sqr(2.*Flavour(kf_c).HadMass())) {
    m_outflav[2*pos]   = Flavour(kf_c);
    m_outflav[2*pos+1] = m_outflav[2*pos].Bar();
  }
  else {
    if (m_fails<5)
      msg_Error()<<METHOD<<" couldn't split diffractive state ["<<diff<<"] "
		 <<"with mass = "<<sqrt(M2)<<".\n";
    m_fails++;
    return false;
  }
  if (anti) for (size_t i=0;i<2;i++) m_outflav[2*pos+i]=m_outflav[2*pos+i].Bar(); 
  return true;
}

bool NonPerturbative_XSecs::
SplitDiffractiveStateIntoHadrons(const size_t & pos,const double & M2) {
  /////////////////////////////////////////////////////////////////////////////////////////
  // This is a safety measure, currently only encoded for protons, which is only
  // invoked if a diffractive sysmte cannot be split into its (di-)quark constituents.
  // We tentatively split the diffracted proton p^* with equal probability either
  // p* -> p + pi^0   or   p^* -> n + pi^+ 
  // and return whether the combined decay masses are larger or smaller than the mass
  // of the diffractive system.
  // TODO: We need to extend this method also to neutrons and make sure we are safe for
  // diffracted vector mesons (although this should work with default parameters).
  /////////////////////////////////////////////////////////////////////////////////////////
  Flavour split = m_inflav[pos];
  kf_code diff = m_outflav[2*pos].Kfcode();
  if (diff==2212) {
    if (ran->Get()<0.5) {
      m_outflav[2*pos]=Flavour(kf_p_plus); m_outflav[2*pos+1]=Flavour(kf_pi);
    }
    else {
      m_outflav[2*pos]=Flavour(kf_n); m_outflav[2*pos+1]=Flavour(kf_pi_plus);
    }
  }
  if (split.IsAnti()) for (size_t i=0;i<2;i++) m_outflav[2*pos+i]=m_outflav[2*pos+i].Bar(); 
  return (m_outflav[2*pos].HadMass()+m_outflav[2*pos+1].HadMass()<sqrt(M2));
}

Blob * NonPerturbative_XSecs::InitBlob(const double & muR2,const double & muQ2) {
  /////////////////////////////////////////////////////////////////////////////////////////
  // Initialising and filling the scatter blob with the outgoing particles; their flavours
  // and momenta have been fixed before.  Incoming flavours and momenta are taken directly
  // from the remnants, accessed through the Remnant_Handler p_remnants.
  // Coloured outgoing particles are assigned colour indices, and in their presence the
  // blob will also be given the "needs_showers" status, which would normally be absent and
  // replaced with needs_beams | needs_hadrondecays.
  /////////////////////////////////////////////////////////////////////////////////////////
  Blob * blob = new Blob();
  blob->SetId();
  blob->AddData("WeightsMap",new Blob_Data<Weights_Map>( m_integrator.TotalXSec()*1e9 ));
  for (size_t i=0;i<2;i++) {
    Particle * part = new Particle(-1,m_inflav[i],m_inmom[i],'I');
    part->SetNumber();
    part->SetBeam(p_remnants->GetRemnant(i)->Beam());
    blob->AddToInParticles(part);
  }
  bool needs_showers = false;
  for (size_t i=0;i<4;i++) {
    if (m_outflav[i]==Flavour(kf_none)) continue;
    Particle * part = new Particle(-1,m_outflav[i],m_outmom[i],'F');
    if (part->Flav().IsQuark())        part->SetFlow(part->Flav().IsAnti()?2:1,500+i/2);
    else if (part->Flav().IsDiQuark()) part->SetFlow(part->Flav().IsAnti()?1:2,500+i/2);
    if (part->Flav().IsQuark() || part->Flav().IsDiQuark()) needs_showers = true;
    part->SetFinalMass(m_outmasses[i]);
    part->SetNumber();
    blob->AddToOutParticles(part);
  }
  if (needs_showers) {
    blob->AddData("Renormalization_Scale",new Blob_Data<double>(sqrt(muR2)));
    blob->AddData("Factorization_Scale",new Blob_Data<double>(sqrt(muQ2)));
    blob->AddData("Resummation_Scale",new Blob_Data<double>(sqrt(muQ2)));
    blob->SetType(btp::Hard_Collision);
    blob->SetStatus(blob_status::needs_showers);
    m_muf2 = muQ2;
    m_mur2 = muR2;
  }
  else {
    blob->AddData("Renormalization_Scale",new Blob_Data<double>(sqrt(0.)));
    blob->AddData("Factorization_Scale",new Blob_Data<double>(sqrt(0.)));
    blob->AddData("Resummation_Scale",new Blob_Data<double>(sqrt(0.)));
    blob->SetType(btp::Elastic_Collision);
    blob->SetStatus(blob_status::needs_beams | blob_status::needs_hadrondecays);
    m_muf2 = 0.;
    m_mur2 = 0.;
  }
  return blob;
}

void NonPerturbative_XSecs::Tests() {
  InitHistos();
  double E = 7000;
  double p = sqrt(sqr(sqr(E)-m_masses2[0]-m_masses2[1])-4.*m_masses2[0]*m_masses2[1])/(2.*E);
  m_s = sqr(E);
  (*p_xsecs)(m_s);
  m_hadtags[0] = m_hadtags[1] = 0;
  m_inmom[0]   = Vec4D(E/2., 0., 0.,  p);
  m_inmom[1]   = Vec4D(E/2., 0., 0., -p);
  for (size_t i=0;i<2;i++) {
    m_outflav[2*i]    = m_inflav[i];
    m_outflav[2*i+1]  = Flavour(kf_none);
    m_outmasses[2*i]  = m_masses[i];
    m_outmasses2[2*i] = m_masses2[i];
  }
  TestElastic();
  TestSingleD();
  TestDoubleD();
  WriteHistos();
}

void NonPerturbative_XSecs::TestElastic() {
  double t, dsigmaEl;
  double conv     = 1.e9/rpa->Picobarn(), sigmaEl = 0.;
  m_weight        = p_xsecs->SigmaEl();
  double binwidth = m_maxT/double(m_Nsteps);
  for (size_t i=0;i<m_Nsteps;i++) {
    t        = (double(i)+0.5) * m_maxT/double(m_Nsteps);
    sigmaEl += dsigmaEl = DiffElXSec(m_s,-t);
    m_histos[string("dsigmaEl_by_dt(Ex)")]->Insert(t,conv*dsigmaEl/binwidth);
  }
  for (size_t i=0;i<m_Nevents;i++) {
    Blob * blob = ElasticScatter();
    t = dabs((m_inmom[0] - m_outmom[0]).Abs2());
    m_histos[string("dsigmaEl_by_dt(MC)")]->Insert(t,m_weight);
  }
}

void NonPerturbative_XSecs::TestSingleD() {
  double t, M, dsigmaSD, prev, act, conv = 1.e9/rpa->Picobarn();
  unsigned long int intsteps = 1000.;
  m_weight       = p_xsecs->SigmaSDA();
  double Mmin    = sqr(m_masses2[0]+2.*m_mpi), Mmax = 20.+Mmin;

  double sigmaSD = 0.;
  double bint    = m_maxT/double(m_Nsteps);
  double binM    = (Mmax-Mmin)/double(intsteps);
  for (size_t i=0;i<m_Nsteps;i++) {
    t        = (double(i)+0.5) * bint;
    prev     = 0.;
    dsigmaSD = 0.;
    for (long int j=0;j<intsteps;j++) {
      M         = Mmin + double(j)*binM;
      act       = DiffSDXSec(m_s,-t, M*M, 0);
      dsigmaSD += M*(prev+act)/2.*binM;
      prev      = act;
    }
    sigmaSD    += dsigmaSD;
    m_histos[string("dsigmaSD1_by_dt(Ex)")]->Insert(t,dsigmaSD/bint);
  }

  sigmaSD = 0.;
  binM    = (Mmax-Mmin)/double(m_Nsteps);
  bint    = m_maxT/double(intsteps);
  for (size_t i=0;i<m_Nsteps;i++) {
    M        = Mmin + double(i)*binM;
    prev     = 0.;
    dsigmaSD = 0.;
    for (long int j=0;j<intsteps;j++) {
      t         = (double(j)+0.5) * bint;
      act       = DiffSDXSec(m_s,-t, M*M, 0);
      dsigmaSD += (prev+act)/2.*bint;
      prev      = act;
    }
    sigmaSD    += dsigmaSD;
    m_histos[string("dsigmaSD1_by_dM(Ex)")]->Insert(M,M*dsigmaSD/binM);
  }
  
  Vec4D  mom;
  for (size_t i=0;i<m_Nevents/100;i++) {
    Blob * blob = SingleDiffractiveScatter(0);
    mom = m_outmom[0]+m_outmom[1];
    t   = dabs((m_inmom[0] - mom).Abs2());
    M   = sqrt(mom.Abs2());
    m_histos[string("dsigmaSD1_by_dt(MC)")]->Insert(t,m_weight);
    m_histos[string("dsigmaSD1_by_dM(MC)")]->Insert(M,m_weight);
  }
}

void NonPerturbative_XSecs::TestDoubleD() {
  double t, m1, m2, dsigmaDD, prev, act, conv = 1.e9/rpa->Picobarn();
  unsigned long int intsteps = 1000.;
  m_weight       = p_xsecs->SigmaDD();
  double Mmin    = sqr(m_masses2[0]+2.*m_mpi), Mmax = 20.+Mmin;

  double sigmaDD = 0.;
  double bint    = m_maxT/double(m_Nsteps);
  double binM1   = (Mmax-Mmin)/double(intsteps);
  double binM2   = (Mmax-Mmin)/double(intsteps);
  std::array<double, 2> M2;
  for (size_t i=0;i<m_Nsteps;i++) {
    t        = (double(i)+0.5) * bint;
    prev     = 0.;
    dsigmaDD = 0.;
    for (long int j=0;j<intsteps;j++) {
      m1          = Mmin + double(j)*binM1;
      M2[0]       = m1*m1;
      for (long int k=0;k<intsteps;k++) {
	m2        = Mmin + double(k)*binM2;
	M2[1]     = m2*m2;
	act       = DiffDDXSec(m_s,-t, M2);
	dsigmaDD += m1*m2*(prev+act)/2.*binM1*binM2;
	prev      = act;
      }
    }
    sigmaDD    += dsigmaDD;
    m_histos[string("dsigmaDD_by_dt(Ex)")]->Insert(t,dsigmaDD/bint);
  }

  sigmaDD = 0.;
  binM1   = (Mmax-Mmin)/double(m_Nsteps);
  binM2   = (Mmax-Mmin)/double(intsteps);
  bint    = m_maxT/double(intsteps);
  for (size_t i=0;i<m_Nsteps;i++) {
    m1       = Mmin + double(i)*binM1;
    M2[0]    = m1*m1;
    prev     = 0.;
    dsigmaDD = 0.;
    for (long int k=0;k<intsteps;k++) {
      m2     = Mmin + double(k)*binM2;
      M2[1]  = m2*m2;
      for (long int j=0;j<intsteps;j++) {
	t         = (double(j)+0.5) * bint;
	act       = DiffDDXSec(m_s,-t, M2);
	dsigmaDD += m1*m2*(prev+act)/2.*bint*binM2;
	prev      = act;
      }
    }
    sigmaDD    += dsigmaDD;
    m_histos[string("dsigmaDD_by_dM1(Ex)")]->Insert(m1,dsigmaDD/binM1);
  }

  Vec4D  mom1, mom2;
  for (size_t i=0;i<m_Nevents/100;i++) {
    Blob * blob = DoubleDiffractiveScatter();
    mom1 = m_outmom[0]+m_outmom[1];
    mom2 = m_outmom[2]+m_outmom[3];
    t    = dabs((m_inmom[0] - mom1).Abs2());
    m1   = sqrt(mom1.Abs2());
    m2   = sqrt(mom2.Abs2());
    m_histos[string("dsigmaDD_by_dt(MC)")]->Insert(t,m_weight);
    m_histos[string("dsigmaDD_by_dM1(MC)")]->Insert(m1,m_weight);
  }
}

void NonPerturbative_XSecs::InitHistos() {
  if (m_ana) {
    m_Nevents = 10000000;
    m_Nsteps  = 200;
    m_maxT    =  2. ;
    m_minM    =  1.2;
    m_maxM    = 21.2;
    m_histos[string("dsigmaEl_by_dt(MC)")]   = new Histogram(0,    0.,m_maxT,m_Nsteps);
    m_histos[string("dsigmaEl_by_dt(Ex)")]   = new Histogram(0,    0.,m_maxT,m_Nsteps);
    m_histos[string("dsigmaSD1_by_dt(MC)")]  = new Histogram(0,    0.,m_maxT,m_Nsteps);
    m_histos[string("dsigmaSD1_by_dt(Ex)")]  = new Histogram(0,    0.,m_maxT,m_Nsteps);
    m_histos[string("dsigmaSD1_by_dM(MC)")]  = new Histogram(0,m_minM,m_maxM,m_Nsteps);
    m_histos[string("dsigmaSD1_by_dM(Ex)")]  = new Histogram(0,m_minM,m_maxM,m_Nsteps);
    m_histos[string("dsigmaDD_by_dt(MC)")]   = new Histogram(0,    0.,m_maxT,m_Nsteps);
    m_histos[string("dsigmaDD_by_dt(Ex)")]   = new Histogram(0,    0.,m_maxT,m_Nsteps);
    m_histos[string("dsigmaDD_by_dM1(MC)")]  = new Histogram(0,m_minM,m_maxM,m_Nsteps);
    m_histos[string("dsigmaDD_by_dM1(Ex)")]  = new Histogram(0,m_minM,m_maxM,m_Nsteps);
  }
}

void NonPerturbative_XSecs::WriteHistos() {
  if (m_ana) {
    std::string name  = std::string("NonPert_Analysis/");
    std::string dat   = std::string(".dat");
    for (std::map<std::string, ATOOLS::Histogram * >::iterator hit=m_histos.begin();
	 hit!=m_histos.end();hit++) {
      hit->second->Finalize();
      hit->second->Output(name+hit->first+dat);
      delete hit->second;
    }
  }
}

