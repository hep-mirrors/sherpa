#include "AHADIC++/Tools/Splitter_Base.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Splitter_Base::Splitter_Base(list<Cluster *> * cluster_list,
			     Soft_Cluster_Handler * softclusters) :
  p_cluster_list(cluster_list), p_softclusters(softclusters)
{}

void Splitter_Base::Init() {
  p_singletransitions = hadpars->GetSingleTransitions();
  p_doubletransitions = hadpars->GetDoubleTransitions();
  p_constituents      = hadpars->GetConstituents();
  m_flavourselector.InitWeights();
  m_ktselector.Init();
  m_zselector.Init();
  m_minmass           = m_flavourselector.MinimalMass();
}

bool Splitter_Base::
operator()(Proto_Particle * part1,Proto_Particle * part2,
	   Proto_Particle * part3) {
  if (!InitSplitting(part1,part2,part3)) return false;
  size_t attempts=100;
  do { attempts--; } while(attempts>0 && !MakeSplitting());
  //msg_Out()<<"************* "<<METHOD<<" with "<<attempts<<" attempts "
  //	   <<"(z1,2 = "<<m_z1<<", "<<m_z2<<" and kt = "
  //	   <<sqrt(m_kt2)<<" < "<<m_ktmax<<").\n";
  return (attempts>0);
}

bool Splitter_Base::
InitSplitting(Proto_Particle * part1,Proto_Particle * part2,
	      Proto_Particle * part3)
{
  p_part1 = part1;
  p_part2 = part2;
  p_part3 = part3;
  FillMasses();
  ConstructLightCone();
  ConstructPoincare();
  return (m_Emax>2.*m_minmass);
}

void Splitter_Base::FillMasses() {
  m_barrd = ((p_part1->Flavour().IsQuark() && p_part1->Flavour().IsAnti()) ||
	     (p_part1->Flavour().IsDiQuark() && !p_part1->Flavour().IsAnti()));
  m_flavs1.first  = m_barrd?p_part1->Flavour().Bar():p_part1->Flavour();
  m_flavs2.second = m_barrd?p_part2->Flavour().Bar():p_part2->Flavour();
  m_Q2    = (p_part1->Momentum()+p_part2->Momentum()).Abs2();
  m_E     = sqrt(m_Q2)/2.;
  m_mass1 = p_constituents->Mass(p_part1->Flavour());
  m_mass2 = p_constituents->Mass(p_part2->Flavour());
  m_m12   = sqr(m_mass1);
  m_m22   = sqr(m_mass2);
  m_Emax  = 2.*m_E-m_mass1-m_mass2;
  if (p_part3!=0) {
    m_mass3 = p_constituents->Mass(p_part3->Flavour());
    m_m32   = sqr(m_mass1);
    m_flavs2.second = m_barrd?p_part3->Flavour().Bar():p_part3->Flavour();
  }
}

void Splitter_Base::ConstructLightCone(const double & kt2) {
  m_alpha = m_beta = 0;
  if (m_m12>1.e-6 && m_m22>1.e-6) {
    double lambda = Lambda(m_Q2,m_m12,m_m22,kt2);
    m_alpha = (m_Q2+m_m12-m_m22)/(2.*m_Q2)+lambda;
    m_beta  = (m_Q2+m_m22-m_m12)/(2.*m_Q2)+lambda;
  }
  else if (m_m12>1.e-6) {
    m_alpha = 1.;
    m_beta  = 1.-m_m12/m_Q2;
  }
  else if (m_m22>1.e-6) {
    m_alpha = 1.-m_m22/m_Q2;
    m_beta  = 1.;
  }
}

void Splitter_Base::ConstructPoincare() {
  Vec4D mom1(p_part1->Momentum()), mom2(p_part2->Momentum());
  m_boost = Poincare(mom1+mom2);
  m_boost.Boost(mom1);
  m_rotat = Poincare(mom1,m_E*s_AxisP); 
}

bool Splitter_Base::MakeSplitting() {
  PopFlavours();
  DetermineMinimalMasses();
  if (MakeKinematics()) {
    FillParticlesInLists();
    return true;
  }
  return false;
}

void Splitter_Base::PopFlavours() {
  // Here we should set vetodi = false -- but no heavy baryons (yet)
  Flavour flav    = m_flavourselector(m_Emax/2.,false);
  m_newflav1      = m_barrd?flav:flav.Bar();
  m_newflav2      = m_newflav1.Bar();
  m_popped_mass   = p_constituents->Mass(flav);
  m_popped_mass2  = sqr(m_popped_mass);
  // m_barrd = true  if part1 = AntiQuark or DiQuark
  // m_barrd = false if part1 = Quark or AntiDiQuark
  m_flavs1.second = m_barrd?m_newflav1.Bar():m_newflav1;
  m_flavs2.first  = m_barrd?m_newflav2.Bar():m_newflav2;
}

void Splitter_Base::DetermineMinimalMasses() {
  if (!m_flavs1.first.IsGluon() && !m_flavs1.second.IsGluon()) {
    m_minQ_1 = p_singletransitions->GetLightestMass(m_flavs1)+0.1;
  }
  else {
    m_minQ_1 = (p_constituents->Mass(m_flavs1.first)+
		p_constituents->Mass(m_flavs1.second));
  }
  if (!m_flavs2.first.IsGluon() && !m_flavs2.second.IsGluon()) {
    m_minQ_2 = p_singletransitions->GetLightestMass(m_flavs2)+0.1;
  }
  else {
    m_minQ_2 = (p_constituents->Mass(m_flavs2.first)+
		p_constituents->Mass(m_flavs2.second));
  }
}

bool Splitter_Base::MakeKinematics() {
  MakeTransverseMomentum();
  return (MakeLongitudinalMomenta() && CheckKinematics());
}

void Splitter_Base::MakeTransverseMomentum() {
  m_ktmax        = (m_Emax-2.*m_popped_mass)/2.;
  if (m_ktmax<0.) {
    msg_Error()<<METHOD<<" yields error ktmax = "<<m_ktmax
	       <<" from "<<m_Emax<<", "<<m_popped_mass<<" vs. "
	       <<" min = "<<m_minmass<<".\n";
    abort();
  }
  m_kt    = m_ktselector(m_ktmax);
  m_kt2   = m_kt*m_kt;
  m_phi   = 2.*M_PI*ran->Get();
  m_ktvec = m_kt * Vec4D(0.,cos(m_phi),sin(m_phi),0.);
}
