#include "AHADIC++/Tools/Splitter_Base.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Splitter_Base::Splitter_Base(list<Cluster *> * cluster_list,
			     Soft_Cluster_Handler * softclusters) :
  p_cluster_list(cluster_list), p_softclusters(softclusters),
  m_ktorder(false), m_ktfac(1.),
  m_attempts(100)
{
  m_analyse = false; //hadpars->Switch("Analysis");
}

Splitter_Base::~Splitter_Base() {
  Histogram * histo;
  string name;
  for (map<string,Histogram *>::iterator hit=m_histograms.begin();
       hit!=m_histograms.end();hit++) {
    histo = hit->second;
    name  = string("Fragmentation_Analysis/")+hit->first+string(".dat");
    histo->Output(name);
    delete histo;
  }
  m_histograms.clear();
}

void Splitter_Base::Init(const bool & isgluon) {
  p_singletransitions = hadpars->GetSingleTransitions();
  p_doubletransitions = hadpars->GetDoubleTransitions();
  p_constituents      = hadpars->GetConstituents();
  m_flavourselector.InitWeights();
  m_ktorder  = (hadpars->Switch("KT_Ordering")>0);
  m_ktmax    = hadpars->Get("kT_max");
  m_ktselector.Init(isgluon);
  m_zselector.Init(this);
  m_minmass  = m_flavourselector.MinimalMass();
}

bool Splitter_Base::
operator()(Proto_Particle * part1,Proto_Particle * part2,
	   Proto_Particle * part3) {
  if (!InitSplitting(part1,part2,part3)) {
    return false;
  }
  size_t attempts(m_attempts);
  do { attempts--; } while(attempts>0 && !MakeSplitting());
  return (attempts>0);
}

bool Splitter_Base::
InitSplitting(Proto_Particle * part1,Proto_Particle * part2,
	      Proto_Particle * part3)
{
  p_part[0] = part1; p_part[1] = part2; p_part[2] = part3;
  FillMasses();
  ConstructLightCone();
  ConstructPoincare();
  return (m_Emax>2.*m_minmass);
}

void Splitter_Base::FillMasses() {
  m_barrd = ((p_part[0]->Flavour().IsQuark() && p_part[0]->Flavour().IsAnti()) ||
	     (p_part[0]->Flavour().IsDiQuark() && !p_part[0]->Flavour().IsAnti()));
  m_flavs1.first  = m_barrd?p_part[0]->Flavour().Bar():p_part[0]->Flavour();
  m_flavs2.second = m_barrd?p_part[1]->Flavour().Bar():p_part[1]->Flavour();
  if (p_part[2]!=0) {
    m_flavs2.second = m_barrd?p_part[2]->Flavour().Bar():p_part[2]->Flavour();
  }
  m_Q2    = (p_part[0]->Momentum()+p_part[1]->Momentum()).Abs2();
  m_Q     = sqrt(m_Q2);
  m_E     = m_Q/2.;
  m_Emax  = m_Q;
  for (size_t i=0;i<3;i++) {
    m_mass[i] = (p_part[i]==0)?0.:p_constituents->Mass(p_part[i]->Flavour());
    m_m2[i]   = sqr(m_mass[i]);
    if (i!=2) m_Emax -= m_mass[i];
  }
}

void Splitter_Base::ConstructLightCone(const double & kt2) {
  m_lc[0] = m_lc[1] = 0.;
  if (m_m2[0]>1.e-6 && m_m2[1]>1.e-6) {
    double lambda = Lambda(m_Q2,m_m2[0],m_m2[1],kt2);
    for (size_t i=0;i<2;i++)
      m_lc[i] = (m_Q2+m_m2[i]-m_m2[1-i])/(2.*m_Q2)+lambda;
  }
  else {
    for (size_t i=0;i<2;i++) {
      if (m_m2[i]>1.e-6) {
	m_lc[i]   = 1.;
	m_lc[1-i] = 1.-m_m2[i]/m_Q2;
      }
    }
  }
}

void Splitter_Base::ConstructPoincare() {
  Vec4D mom1(p_part[0]->Momentum()), mom2(p_part[1]->Momentum());
  m_boost = Poincare(mom1+mom2);
  m_boost.Boost(mom1);
  m_rotat = Poincare(mom1,m_E*s_AxisP); 
}

bool Splitter_Base::MakeSplitting() {
  PopFlavours();
  DetermineMinimalMasses();
  return (MakeKinematics() && FillParticlesInLists());
}

void Splitter_Base::PopFlavours() {
  // Here we should set vetodi = false -- but no heavy baryons (yet)
  Flavour flav    = m_flavourselector(m_Emax/2.,false);
  // m_barrd = true  if part1 = AntiQuark or DiQuark
  // m_barrd = false if part1 = Quark or AntiDiQuark
  m_newflav[0]    = m_barrd?flav:flav.Bar();
  m_newflav[1]    = m_newflav[0].Bar();
  m_popped_mass   = p_constituents->Mass(flav);
  m_popped_mass2  = sqr(m_popped_mass);
  m_flavs1.second = m_barrd?m_newflav[0].Bar():m_newflav[0];
  m_flavs2.first  = m_barrd?m_newflav[1].Bar():m_newflav[1];
}

void Splitter_Base::DetermineMinimalMasses() {
  m_msum[0]  = (p_constituents->Mass(m_flavs1.first)+
		p_constituents->Mass(m_flavs1.second));
  m_msum[1]  = (p_constituents->Mass(m_flavs2.first)+
		p_constituents->Mass(m_flavs2.second));
  m_mdec[0] = p_doubletransitions->GetLightestMass(m_flavs1);
  m_mdec[1] = p_doubletransitions->GetLightestMass(m_flavs2);
  if (!m_flavs1.first.IsGluon() && !m_flavs1.second.IsGluon()) {
    if (!(m_flavs1.first.IsDiQuark() && m_flavs1.second.IsDiQuark())) {
      m_minQ[0] = Min(m_mdec[0],
		      Max(0.,p_singletransitions->GetLightestMass(m_flavs1)));
    }
    else {
      m_minQ[0] = m_mdec[0];
    }
  }
  if (!m_flavs2.first.IsGluon() && !m_flavs2.second.IsGluon()) {
    if (!(m_flavs2.first.IsDiQuark() && m_flavs2.second.IsDiQuark())) {
      m_minQ[1] = Min(m_mdec[1],
		      Max(0.,p_singletransitions->GetLightestMass(m_flavs2)));
    }
    else {
      m_minQ[1] = m_mdec[1];
    }
  }
  else {
    for (size_t i=0;i<2;i++) m_minQ[i] = m_msum[i];
  }
  for (size_t i=0;i<2;i++) {
    m_minQ2[i] = sqr(m_minQ[i]);
    m_msum2[i] = sqr(m_msum[i]);
    m_mdec2[i] = sqr(m_mdec[i]);
  }
}

bool Splitter_Base::MakeKinematics() {
  MakeTransverseMomentum();
  return (MakeLongitudinalMomenta() && CheckKinematics());
}

void Splitter_Base::MakeTransverseMomentum() {
  m_ktfac  = Max(1.,m_Q2/(4.*m_minQ[0]*m_minQ[1]));
  m_kt2max = Min(p_part[0]->KT2_Max(),p_part[1]->KT2_Max());
  double ktmax  = Min(m_ktmax,
		      (m_ktorder?
		       Min(sqrt(m_kt2max),(m_Emax-2.*m_popped_mass)/2.):
		       (m_Emax-2.*m_popped_mass)/2.));
  // have to make this a parameter for the beam breakup?
  //if (p_part[0]->IsBeam() || p_part[1]->IsBeam()) ktmax = Min(5.0,m_ktmax);
  if (ktmax<0.) {
    msg_Error()<<METHOD<<" yields error ktmax = "<<ktmax
	       <<" from "<<m_Emax<<", "<<m_popped_mass<<" vs. "
	       <<" min = "<<m_minmass<<".\n";
    abort();
  }
  m_ktfac = 1.;
  bool islead = p_part[0]->IsLeading() || p_part[1]->IsLeading();
  m_kt    = m_ktselector(ktmax,m_ktfac);
  m_kt2   = m_kt*m_kt;
  m_phi   = 2.*M_PI*ran->Get();
  m_ktvec = m_kt * Vec4D(0.,cos(m_phi),sin(m_phi),0.);
}
