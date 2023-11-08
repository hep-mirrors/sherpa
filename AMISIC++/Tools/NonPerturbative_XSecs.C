#include "AMISIC++/Tools/NonPerturbative_XSecs.H"
#include "BEAM/Main/Beam_Base.H"
#include "REMNANTS/Main/Remnant_Handler.H"
#include "ATOOLS/Phys/Weights.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace AMISIC;
using namespace ATOOLS;
using namespace std;

NonPerturbative_XSecs::
NonPerturbative_XSecs(REMNANTS::Remnant_Handler * remnants,Hadronic_XSec_Calculator * xsecs) :
  p_remnants(remnants), p_xsecs(xsecs), m_variableS(false),
  m_evttype(mipars->GetEvtType()),
  m_smin(p_xsecs->Smin()), 
  m_eps_pomeron(p_xsecs->EpsPomeron()), m_alphaP_pomeron(p_xsecs->AlphaPPomeron()),
  m_triple_pomeron(p_xsecs->TriplePomeron()), m_alphaQED(p_xsecs->AlphaQED()),
  m_s0(1./m_alphaP_pomeron),
  m_mpi(Flavour(kf_pi).HadMass()), m_deltaMres(1.), m_cres(2.),
  m_inflav(p_xsecs->GetFlavs()),
  m_ana(false)
{
  for (size_t i=0;i<2;i++) p_beams[i] = NULL;
  for (size_t i=0;i<2;i++) {
    m_masses[i] = m_inflav[i].HadMass(); m_masses2[i] = sqr(m_masses[i]);
  }
  if (m_ana) { Tests(); exit(0); }
}

void NonPerturbative_XSecs::SetBeams(BEAM::Beam_Base * beam1,BEAM::Beam_Base * beam2) {
  p_beams[0] = beam1; p_beams[1] = beam2;
  for (size_t i=0;i<2;i++) {
    if (p_beams[i]!=NULL && p_beams[i]->Type()==BEAM::beamspectrum::EPA) {
      m_variableS = true;
      break;
    }
  }
}

void NonPerturbative_XSecs::CalculateSDependentCrossSections() {
  if (!m_variableS) return;	
  m_integrator.Init(p_xsecs,p_beams[0],p_beams[1],m_evttype);
}

Blob * NonPerturbative_XSecs::MakeScatter() {
  if (m_variableS) {
    if (!m_integrator()) {
      msg_Error()<<"Error in "<<METHOD<<": could not create initial-state phase-space.\n";
      return NULL;
    }
  }
  for (size_t beam=0;beam<2;beam++) m_inmom[beam] = p_remnants->GetRemnant(beam)->InMomentum();
  m_s     = (m_inmom[0]+m_inmom[1]).Abs2();
  m_boost = Poincare(m_inmom[0]+m_inmom[1]);
  for (size_t beam=0;beam<2;beam++) m_boost.Boost(m_inmom[beam]);
  if (m_s<m_smin) THROW(fatal_error,"Insufficient s = "+ToString(m_s)+" vs "+ToString(m_smin)+
			" in Non-Perturbative Events in AMISIC++");
  (*p_xsecs)(m_s);
  array<Flavour, 2> flavs;
  switch(SelectMode()) {
  case event_mode::elastic:
    if (!p_xsecs->SelectEl(flavs) || !FixFS(flavs))  THROW(fatal_error,"Could not define elastic FS");
    return ElasticScatter();       
  case event_mode::SDA:
    if (!p_xsecs->SelectSDA(flavs) || !FixFS(flavs)) THROW(fatal_error,"Could not define SD(A) FS");
    return SingleDiffractiveScatter(0);
  case event_mode::SDB:
    if (!p_xsecs->SelectSDB(flavs) || !FixFS(flavs)) THROW(fatal_error,"Could not define SD(B) FS");
    return SingleDiffractiveScatter(1);
  case event_mode::DD:
    if (!p_xsecs->SelectDD(flavs) || !FixFS(flavs))  THROW(fatal_error,"Could not define DD FS");
    return DoubleDiffractiveScatter();
  case event_mode::unknown:
  default:                  break;
  }
  THROW(fatal_error,"Unknown Event Mode in Non-Perturbative Events in AMISIC++");
  return NULL;
}

const event_mode::code NonPerturbative_XSecs::SelectMode() {
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
  for (size_t i=0;i<2;i++) {
    m_outflav[2*i]    = flavs[i];
    m_outmasses[2*i]  = m_outflav[2*i].Mass();
    m_outmasses2[2*i] = ATOOLS::sqr(m_outmasses[2*i]);
    m_hadtags[i]      = p_xsecs->Index(flavs[i],i);
    if (m_hadtags[i]>99) return false;
  }
  return true;
}

double NonPerturbative_XSecs::DiffElXSec(const double & s,const double & t) const
{
  return ( 1./(16.*M_PI) * sqr(beta(m_hadtags[0],t) * beta(m_hadtags[1],t)) *
	   exp((2.*m_eps_pomeron + 2.*m_alphaP_pomeron*t) * log(s/m_s0)) ); 
}

Blob * NonPerturbative_XSecs::ElasticScatter() {
  double arg  = 2.*(m_eps_pomeron + m_alphaP_pomeron * log(m_s/m_s0) +
		    p_xsecs->s_slopes[m_hadtags[0]] + p_xsecs->s_slopes[m_hadtags[1]]);
  double t;
  size_t trials=0;
  do { t = ExponentialDist(-m_s,0.,arg); trials++; } while (!FixOutMomenta(t) && trials<100);
  if (trials>=100) {
    msg_Error()<<METHOD<<" fails: t = "<<t<<" yields no momenta.\n";
    return NULL;
  }
  Blob * blob   = InitBlob();
  blob->SetType(btp::Soft_Collision);
  blob->SetStatus(blob_status::needs_beams);
  return blob;
}

double NonPerturbative_XSecs::DiffSDXSec(const double & s, const double & t,
					 const double & M2,const size_t & pos) const
{
  double m     = GetFlavour(pos,m_hadtags[pos]).Mass(), m2 = sqr(m);
  double M2res = sqr(m+m_deltaMres), M2min = sqr(m+2.*m_mpi);
  if (M2<M2min) return 0.;
  return ( F_SD(M2,M2res) / (16.*M_PI*M2) *
	   m_triple_pomeron * sqr(beta(m_hadtags[pos],t)) * beta(m_hadtags[1-pos],0.) *
	   exp((2.*m_eps_pomeron + 2.*m_alphaP_pomeron*t) * log(s/M2)) ); 
}

ATOOLS::Blob * NonPerturbative_XSecs::SingleDiffractiveScatter(const size_t & pos) {
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
  blob->SetType(btp::Hard_Collision);
  blob->SetStatus(blob_status::needs_showers);
  return blob;
}

double NonPerturbative_XSecs::DiffDDXSec(const double & s,const double & t,
					 const std::array<double, 2> & M2) const
{
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
  blob->SetType(btp::Hard_Collision);
  blob->SetStatus(blob_status::needs_showers);
  return blob;
}

bool NonPerturbative_XSecs::FixOutMomenta(const double & t) {
  if (m_s<sqr(sqrt(m_outmasses2[0])+sqrt(m_outmasses2[2]))) return false;
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
  return true;
}

bool NonPerturbative_XSecs::SplitDiffractiveState(const size_t & pos) {
  Vec4D  diffmom  = m_outmom[2*pos];
  double M2       = diffmom.Abs2();
  if (!SelectFlavoursOfDiffraction(pos,M2)) return false;
  if (M2<sqr(m_outmasses[2*pos]+m_outmasses[2*pos+1])) return false;
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

bool NonPerturbative_XSecs::SelectFlavoursOfDiffraction(const size_t & pos,const double & M2) {
  kf_code diff = m_outflav[2*pos].Kfcode();
  bool    anti = m_outflav[2*pos].IsAnti();
  if (diff==2212) {
    double random = ran->Get();
    if (random<1./3. && M2>sqr(Flavour(kf_d).Mass()+Flavour(kf_uu_1).Mass())) {
      m_outflav[2*pos]=Flavour(kf_d); m_outflav[2*pos+1]=Flavour(kf_uu_1);
    }
    else if (random<1./2. && M2>sqr(Flavour(kf_u).Mass()+Flavour(kf_ud_1).Mass())) {
      m_outflav[2*pos]=Flavour(kf_u); m_outflav[2*pos+1]=Flavour(kf_ud_1);
    }
    else if (M2>sqr(Flavour(kf_d).Mass()+Flavour(kf_uu_1).Mass())) {
      m_outflav[2*pos]=Flavour(kf_d); m_outflav[2*pos+1]=Flavour(kf_uu_1);
    }
    else {
      msg_Error()<<METHOD<<" couldn't split diffractive state of mass = "<<sqrt(M2)<<".\n";
      return false;
    }
  }
  else if (diff==113 || diff==223) {
    m_outflav[2*pos]   = ran->Get()>0.5 ? Flavour(kf_u) : Flavour(kf_d);
    m_outflav[2*pos+1] = m_outflav[2*pos].Bar();
  }
  else if (diff==333) {
    m_outflav[2*pos]   = Flavour(kf_s);
    m_outflav[2*pos+1] = m_outflav[2*pos].Bar();
  }
  else if (diff==433) {
    m_outflav[2*pos]   = Flavour(kf_c);
    m_outflav[2*pos+1] = m_outflav[2*pos].Bar();
  }
  else THROW(fatal_error,"No meaningful hadronic state for "+m_outflav[2*pos].IDName());
  if (anti) { for (size_t i=0;i<2;i++) m_outflav[2*pos+i]=m_outflav[2*pos+i].Bar(); }
  for (size_t i=0;i<2;i++) {
    m_outmasses[2*pos+i]  = m_outflav[2*pos+i].Mass();
    m_outmasses2[2*pos+i] = sqr(m_outmasses[2*pos+i]); 
  }
  return true;
}

Blob * NonPerturbative_XSecs::InitBlob(const double & muR2,const double & muQ2) {
  Blob * blob = new Blob();
  blob->SetId();
  blob->AddData("WeightsMap",new Blob_Data<Weights_Map>({}));
  blob->AddData("Renormalization_Scale",new Blob_Data<double>(sqrt(muR2)));
  blob->AddData("Factorization_Scale",new Blob_Data<double>(sqrt(muQ2)));
  blob->AddData("Resummation_Scale",new Blob_Data<double>(sqrt(muQ2)));
  for (size_t i=0;i<2;i++) {
    Particle * part = new Particle(-1,m_inflav[i],m_inmom[i],'I');
    part->SetNumber();
    part->SetBeam(p_remnants->GetRemnant(i)->Beam());
    blob->AddToInParticles(part);
  }
  for (size_t i=0;i<4;i++) {
    if (m_outflav[i]==Flavour(kf_none)) continue;
    Particle * part = new Particle(-1,m_outflav[i],m_outmom[i],'F');
    if (part->Flav().IsQuark())        part->SetFlow(part->Flav().IsAnti()?2:1,500+i/2);
    else if (part->Flav().IsDiQuark()) part->SetFlow(part->Flav().IsAnti()?1:2,500+i/2);
    part->SetNumber();
    blob->AddToOutParticles(part);
  }
  m_muf2 = muQ2;
  m_mur2 = muR2;
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

