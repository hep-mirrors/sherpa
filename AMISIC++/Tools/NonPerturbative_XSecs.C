#include "AMISIC++/Tools/NonPerturbative_XSecs.H"
#include "REMNANTS/Main/Remnant_Handler.H"
#include "ATOOLS/Phys/Weights.H"
#include "ATOOLS/Math/Random.H"

using namespace AMISIC;
using namespace ATOOLS;
using namespace std;

NonPerturbative_XSecs::NonPerturbative_XSecs(Hadronic_XSec_Calculator * xsecs) :
  p_xsecs(xsecs), m_evttype(mipars->GetEvtType()),
  m_eps_pomeron(p_xsecs->EpsPomeron()), m_alphaP_pomeron(p_xsecs->AlphaPPomeron()),
  m_triple_pomeron(p_xsecs->TriplePomeron()), m_alphaQED(p_xsecs->AlphaQED()),
  m_s0(1./m_alphaP_pomeron),
  m_mpi(Flavour(kf_pi).HadMass()), m_deltaMres(1.), m_cres(2.),
  m_inflav(p_xsecs->GetFlavs())
{
  for (size_t i=0;i<2;i++) {
    m_masses[i] = m_inflav[i].HadMass(); m_masses2[i] = sqr(m_masses[i]);
  }
}

Blob * NonPerturbative_XSecs::MakeScatter(REMNANTS::Remnant_Handler * remnants) {
  msg_Out()<<"=============== "<<METHOD<<":\n";
  m_inmom[0] = remnants->GetRemnant(0)->InMomentum();
  m_inmom[1] = remnants->GetRemnant(1)->InMomentum();
  m_s        = (m_inmom[0]+m_inmom[1]).Abs2();
  (*p_xsecs)(m_s);
  array<size_t,2> hadpars; hadpars[0] = hadpars[1] = 0;
  switch(SelectMode(hadpars)) {
  case event_mode::elastic:
    return ElasticScatter(remnants,hadpars);
  case event_mode::SDA:
    return SingleDiffractiveScatter(remnants,hadpars,0);
  case event_mode::SDB:
    return SingleDiffractiveScatter(remnants,hadpars,1);
  case event_mode::DD:
    return DoubleDiffractiveScatter(remnants,hadpars);
  default:
  case event_mode::unknown:
    exit(1);
  }
  return NULL;
}

const event_mode::code NonPerturbative_XSecs::SelectMode(array<size_t,2> & hadpars) {
  double total = 0.;
  if (m_evttype==evt_type::Elastic ||
      m_evttype==evt_type::QuasiElastic) total += p_xsecs->SigmaEl(hadpars);
  if (m_evttype==evt_type::DiffractiveA ||
      m_evttype==evt_type::QuasiElastic) total += p_xsecs->SigmaSDA(hadpars);
  if (m_evttype==evt_type::DiffractiveB ||
      m_evttype==evt_type::QuasiElastic) total += p_xsecs->SigmaSDB(hadpars);
  if (m_evttype==evt_type::DiffractiveAB ||
      m_evttype==evt_type::QuasiElastic) total += p_xsecs->SigmaDD(hadpars);
  total *= ran->Get();
  msg_Out()<<"*  "<<METHOD<<": total = "<<total<<"\n";
  for (size_t i=0;i<2;i++) {
    m_outflav[2*i]    = m_inflav[i];
    m_outflav[2*i+1]  = Flavour(kf_none);
    m_outmasses[2*i]  = m_masses[i];
    m_outmasses2[2*i] = m_masses2[i];
  }
  if (m_evttype==evt_type::DiffractiveAB ||
      m_evttype==evt_type::QuasiElastic) total -= p_xsecs->SigmaDD(hadpars);
  if (total<0) return event_mode::DD;
  if (m_evttype==evt_type::DiffractiveB ||
      m_evttype==evt_type::QuasiElastic) total -= p_xsecs->SigmaSDB(hadpars);
  if (total<0) return event_mode::SDB;
  if (m_evttype==evt_type::DiffractiveA ||
      m_evttype==evt_type::QuasiElastic) total -= p_xsecs->SigmaSDA(hadpars);
  if (total<0) return event_mode::SDA;
  return event_mode::elastic;
}

Blob * NonPerturbative_XSecs::ElasticScatter(REMNANTS::Remnant_Handler * remnants,
					     const array<size_t,2> hadtags) {
  double arg    = 2.*(m_alphaP_pomeron * log(m_s/m_s0) +
		      p_xsecs->s_slopes[hadtags[0]] + p_xsecs->s_slopes[hadtags[1]]);
  double t      = ExponentialDist(-m_s,0.,arg); 
  double p22    = ( sqr(m_s-m_masses2[2]-m_masses2[3]) - 4.*m_masses2[2] *m_masses2[3] )/(4.*m_s);
  double E2     = sqrt(p22+m_masses2[2]), E3 = sqrt(p22+m_masses2[3]);
  double p2     = sqrt(p22);
  double cost   = (t - m_masses2[0] - m_masses2[2] + 2.*m_inmom[0][0]*E2)/(2.*m_inmom[0][3]*p2);
  double sint   = sqrt(1.-cost*cost);
  double phi    = 2.*M_PI*ran->Get();
  //msg_Out()<<METHOD<<" yields t = "<<t<<" for s = "<<m_s<<" "
  //	   <<"("<<m_alphaP_pomeron<<", "<<p_xsecs->s_slopes[hadtags[0]]<<", "
  //	   <<p_xsecs->s_slopes[hadtags[1]]<<")\n"
  //	   <<"*   cos(theta) = "<<cost<<" and sin(theta) = "<<sint<<"\n";
  m_outmom[0]   = Vec4D(E2,p2*sint*cos(phi),p2*sint*sin(phi),p2*cost);
  m_outmom[2]   = Vec4D(E3,-Vec3D(m_outmom[0]));
  Blob * blob   = InitBlob(remnants);
  //msg->SetPrecision(12);
  //msg_Out()<<"*   p2 = "<<p2<<" from s = "<<m_s<<", t = "<<t<<" "
  //	   <<"and masses = "<<m_masses2[2]<<" & "<<m_masses2[3]<<"\n"
  //	   <<"*   moms = "<<m_outmom[0]<<" + "<<m_outmom[1]<<"\n"
  //	   <<(*blob)<<"\n";
  blob->SetType(btp::Soft_Collision);
  blob->SetStatus(blob_status::needs_beams);
  return blob;
}

ATOOLS::Blob * NonPerturbative_XSecs::
SingleDiffractiveScatter(REMNANTS::Remnant_Handler * remnants,
			 const std::array<size_t,2> hadtags,const size_t & diffbeam) {
  double M2min  = sqr(m_masses[2+diffbeam]+2.*m_mpi), M2res = sqr(m_masses[2+diffbeam]+m_deltaMres);
  double argt   = 2.*(m_alphaP_pomeron * log(m_s/M2min) + p_xsecs->s_slopes[hadtags[1-diffbeam]]);
  double argM   = -(1. + m_eps_pomeron);
  double maxval = ExponentialWeight(-m_s,0,argt)*PowerWeight(M2min,m_s,argM)*F_SD(M2min,M2res,m_cres);
  double t, M2, value;
  do {
    t     = ExponentialDist(-m_s,0,argt);
    M2    = PowerDist(M2min,m_s,argM);
    value = ExponentialWeight(-m_s,0,argt)*PowerWeight(M2min,m_s,argM)*F_SD(M2,M2res,m_cres);
    msg_Out()<<METHOD<<": t = "<<t<<", M = "<<sqrt(M2)<<" > "<<sqrt(M2min)<<", "
	     <<"Mres = "<<sqrt(M2res)<<", "
	     <<"value = "<<value<<"/"<<maxval<<"\n";
  } while (value/maxval<ran->Get());
  double p22    = ( sqr(m_s-m_masses2[3-diffbeam]-M2) - 4.*m_masses2[3-diffbeam]*M2 )/(4.*m_s);
  double E2     = sqrt(p22+(diffbeam==0?M2:m_masses2[2])), E3 = sqrt(p22+(diffbeam==1?M2:m_masses2[3]));
  double p2     = sqrt(p22);
  double cost   = ((t - m_masses2[0] - (diffbeam==0?M2:m_masses2[2]) + 2.*m_inmom[0][0]*E2)/
		   (2.*m_inmom[0][3]*p2) );
  double sint   = sqrt(1.-cost*cost);
  double phi    = 2.*M_PI*ran->Get();
  //msg_Out()<<METHOD<<" yields t = "<<t<<" for s = "<<m_s<<" "
  //	   <<"("<<m_alphaP_pomeron<<", "<<p_xsecs->s_slopes[hadtags[0]]<<", "
  //	   <<p_xsecs->s_slopes[hadtags[1]]<<")\n"
  //	   <<"*   cos(theta) = "<<cost<<" and sin(theta) = "<<sint<<"\n";
  m_outmom[0]   = Vec4D(E2,p2*sint*cos(phi),p2*sint*sin(phi),p2*cost);
  m_outmom[2]   = Vec4D(E3,-Vec3D(m_outmom[0]));
  msg_Out()<<"*   p2 = "<<p2<<" from s = "<<m_s<<" and masses = "
	   <<(diffbeam==0?sqrt(M2):m_masses2[2])<<" & "<<(diffbeam==1?sqrt(M2):m_masses2[3])<<"\n"
  	   <<"*   moms = "<<m_outmom[0]<<"("<<sqrt(m_outmom[0].Abs2())<<") +\n"
	   <<"           "<<m_outmom[2]<<"("<<sqrt(m_outmom[2].Abs2())<<")\n";
  SplitDiffractiveState(diffbeam);
  Blob * blob   = InitBlob(remnants,M2,M2);
  blob->SetType(btp::Hard_Collision);
  blob->SetStatus(blob_status::needs_showers);
  msg_Out()<<(*blob)<<"\n";
  return blob;
}

ATOOLS::Blob * NonPerturbative_XSecs::DoubleDiffractiveScatter(REMNANTS::Remnant_Handler * remnants,
							       const std::array<size_t,2> hadtags) {
  return NULL;
}

bool NonPerturbative_XSecs::SplitDiffractiveState(const size_t & pos) {
  Flavour diff     = m_outflav[2*pos];
  Vec4D   diffmom  = m_outmom[2*pos];
  if (diff==Flavour(kf_p_plus) || diff==Flavour(kf_p_plus).Bar()) {
    msg_Out()<<METHOD<<" for a diffracted proton! Mom = "<<diffmom<<"\n";
    double random=ran->Get();
    if (random<1./2.)            { m_outflav[2*pos]=Flavour(kf_u); m_outflav[2*pos+1]=Flavour(kf_ud_0);}
    else if (random<1./2.+1./6.) { m_outflav[2*pos]=Flavour(kf_u); m_outflav[2*pos+1]=Flavour(kf_ud_1);}
    else                         { m_outflav[2*pos]=Flavour(kf_d); m_outflav[2*pos+1]=Flavour(kf_uu_1);}
  }
  if (diff.IsAnti()) { for (size_t i=0;i<2;i++) m_outflav[2*pos+i]=m_outflav[2*pos+i].Bar(); }
  for (size_t i=0;i<2;i++) {
    m_outmasses[2*pos+i]  = m_outflav[2*pos+i].Mass();
    m_outmasses2[2*pos+i] = sqr(m_outmasses[2*pos+i]); 
    msg_Out()<<"Fix masses for "<<m_outflav[2*pos+i]<<" --> "<<m_outmasses[2*pos+i]<<"\n";
  }
  double M2   = diffmom.Abs2();
  if (M2<sqr(m_outmasses[2*pos]+m_outmasses[2*pos+1])) return false;
  double p2   = ( (sqr(M2-m_outmasses2[2*pos]-m_outmasses2[2*pos+1]) -
		   4.*m_outmasses2[2*pos]*m_outmasses2[2*pos+1]) / (4.*M2) ), p = sqrt(p2);
  double E[2];
  for (size_t i=0;i<2;i++) E[i] = sqrt(p2+m_outmasses2[2*pos+i]);
  double cost = -1.+2.*ran->Get(), sint = sqrt(1.-cost*cost), phi = 2.*M_PI*ran->Get();
  m_outmom[2*pos]   = Vec4D(E[0],p*sint*cos(phi),p*sint*sin(phi),p*cost);
  m_outmom[2*pos+1] = Vec4D(E[1],-Vec3D(m_outmom[2*pos]));
  msg_Out()<<METHOD<<"("<<pos<<"): "<<sqrt(M2)<<" --> "
	   <<m_outmom[2*pos]<<" ("<<sqrt(m_outmom[2*pos].Abs2())<<"), "
	   <<m_outmom[2*pos+1]<<" ("<<sqrt(m_outmom[2*pos+1].Abs2())<<")\n";
  Poincare boost(diffmom);
  boost.BoostBack(m_outmom[2*pos]); //m_outmom[2*pos+1]=diffmom-m_outmom[2*pos];
  boost.BoostBack(m_outmom[2*pos+1]); //m_outmom[2*pos+1]=diffmom-m_outmom[2*pos];
  msg_Out()<<diffmom<<" --> "
	   <<m_outmom[2*pos]<<" ("<<sqrt(m_outmom[2*pos].Abs2())<<") + "
	   <<m_outmom[2*pos+1]<<" ("<<sqrt(m_outmom[2*pos+1].Abs2())<<")\n";
  return true;
}

double NonPerturbative_XSecs::DiffElXSec(const array<size_t,2> hadtags,
					 const double & s,const double & t) const
{
  return ( 1./(16.*M_PI) * sqr(beta(hadtags[0],t) * beta(hadtags[1],t)) *
	   exp((2.*m_eps_pomeron + 2.*m_alphaP_pomeron*t) * log(s/m_s0)) ); 
}

Blob * NonPerturbative_XSecs::InitBlob(REMNANTS::Remnant_Handler * remnants,
				       const double & muR,const double & muQ) {
  Blob * blob = new Blob();
  blob->SetId();
  blob->AddData("WeightsMap",new Blob_Data<Weights_Map>({}));
  blob->AddData("Renormalization_Scale",new Blob_Data<double>(muR));
  blob->AddData("Factorization_Scale",new Blob_Data<double>(0.));
  blob->AddData("Resummation_Scale",new Blob_Data<double>(muQ));
  for (size_t i=0;i<2;i++) {
    Particle * part = new Particle(-1,m_inflav[i],m_inmom[i],'I');
    part->SetNumber();
    part->SetBeam(remnants->GetRemnant(i)->Beam());
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
  return blob;
}

