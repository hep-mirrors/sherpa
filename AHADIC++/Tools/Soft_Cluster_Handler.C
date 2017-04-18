#include "AHADIC++/Tools/Soft_Cluster_Handler.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Soft_Cluster_Handler::Soft_Cluster_Handler(list<Proto_Particle *> * hadrons) :
  p_hadrons(hadrons)
{ }

Soft_Cluster_Handler::~Soft_Cluster_Handler() 
{ }

void Soft_Cluster_Handler::Init() {
  p_constituents       = hadpars->GetConstituents();
  p_singletransitions  = hadpars->GetSingleTransitions(); 
  p_doubletransitions  = hadpars->GetDoubleTransitions();
  m_light              = hadpars->Get("decay_threshold");
  m_piphoton_threshold = hadpars->Get("piphoton_threshold");
  m_dipion_threshold   = hadpars->Get("dipion_threshold");
  m_open_threshold     = (2.*p_constituents->MinMass()+
			  hadpars->Get("open_threshold"));
  m_ktselector.Init();
}

bool Soft_Cluster_Handler::MustPromptDecay(Cluster * cluster) {
  FillFlavours(cluster);
  // will assume clusters have to decay, if they are lighter than heaviest
  // single (one-hadron) transition or lighter than heaviest decay into
  // two hadrons
  return (m_mass < p_doubletransitions->GetHeaviestMass(m_flavs) ||
	  m_mass < p_singletransitions->GetHeaviestMass(m_flavs));
}

bool Soft_Cluster_Handler::Treat(Cluster * cluster,bool force)
{
  FillFlavours(cluster);
  if (!force && CheckOutsideRange()) return false;
  return Decay();
}

bool Soft_Cluster_Handler::CheckOutsideRange() {
  // we may want to check if we want to take the full range of possible
  // cluster decays into two hadrons
  double mass_l   = p_singletransitions->GetLightestMass(m_flavs);
  if (mass_l<0. && m_mass<p_doubletransitions->GetLightestMass(m_flavs))
    return false;
  if (m_mass<=mass_l) {
    Flavour had = p_singletransitions->GetLightestTransition(m_flavs);
    msg_Error()<<"Problem spotted in "<<METHOD<<":\n"
	       <<"   test mass = "<<m_mass<<" too light for lightest hadron, "
	       <<had<<" with mass = "<<mass_l<<"\n"
	       <<"   made from <"<<m_flavs.first<<", "<<m_flavs.second<<"> "
	       <<" with masses "<<m_flavs.first.HadMass()<<" + "
	       <<m_flavs.second.HadMass()<<".\n";
    exit(1);
    return true;
  }
  double mass_dec =
    p_doubletransitions->GetLightestMass(m_flavs) * m_light       + 
    p_doubletransitions->GetLightestMass(m_flavs) * (1.-m_light); 
  if (m_mass>mass_dec) return true;
  return false;
}

bool Soft_Cluster_Handler::RadiativeDecay(Cluster * cluster) {
  FillFlavours(cluster);
  if (m_mass>p_singletransitions->GetLightestMass(m_flavs) &&
      RadiationWeight()>0.) {
    m_hads[0] = p_singletransitions->GetLightestTransition(m_flavs);
    m_hads[1] = Flavour(kf_photon);
    return FixKinematics();
  }
  return false;
}

bool Soft_Cluster_Handler::TreatTwoGluons(Cluster * cluster) {
  FillFlavours(cluster);
  // below pi0 + gamma threshold
  if (m_mass < m_piphoton_threshold) {
    m_hads[0] = m_hads[1] = Flavour(kf_photon);
  }
  // below two-pion threshold
  else if (m_mass <  m_dipion_threshold) {
    size_t i(2.*ran->Get());
    m_hads[i]   = Flavour(kf_photon);
    m_hads[1-i] = Flavour(kf_pi);
  }
  // below two-pion threshold
  else {
    if (ran->Get()>0.5) {
      m_hads[0] = m_hads[1] = Flavour(kf_pi);
    }
    else {
      size_t i(2.*ran->Get());
      m_hads[i]   = Flavour(kf_pi);
      m_hads[1-i] = Flavour(kf_pi).Bar();
    }
  }
  return FixKinematics();
}
  
bool Soft_Cluster_Handler::FillFlavours(Cluster * cluster) {
  p_cluster      = cluster;
  m_mass2        = cluster->Momentum().Abs2();
  m_mass         = sqrt(m_mass2);
  m_flavs.first  = (*cluster)[0]->Flavour();
  m_flavs.second = (*cluster)[1]->Flavour();
}

bool Soft_Cluster_Handler::Decay() {
  m_hads[0] = m_hads[1] = Flavour(kf_none);
  double decweight(DecayWeight());
  if (decweight>0.) return FixKinematics();
  m_hads[0] = Flavour(kf_none); m_hads[1] = Flavour(kf_photon);
  double radweight = RadiationWeight();
  if (radweight>0.) return FixKinematics();
  return false;
}

bool Soft_Cluster_Handler::FixKinematics() {
  Vec4D mom1((*p_cluster)[0]->Momentum()), mom2((*p_cluster)[1]->Momentum());
  Poincare boost = Poincare(mom1+mom2);
  boost.Boost(mom1);
  Poincare rotat = Poincare(mom1,s_AxisP); 

  double M2(m_mass*m_mass);
  double m12(sqr(m_hads[0].Mass())),m22(sqr(m_hads[1].Mass()));
  double E1((M2+m12-m22)/(2.*m_mass)), p12(sqr(E1)-m12);
  double ktmax = sqrt(p12);
  double kt    = m_ktselector(ktmax), kt2 = kt*kt;
  double pl    = sqrt(E1*E1-m12-kt2);
  double phi   = 2.*M_PI*ran->Get();
  m_moms[0] = Vec4D(E1,kt*cos(phi),kt*sin(phi),pl);
  m_moms[1] = Vec4D(m_mass-E1,-kt*cos(phi),-kt*sin(phi),-pl);
  for (size_t i=0;i<2;i++) {
    rotat.RotateBack(m_moms[i]);
    boost.BoostBack(m_moms[i]);
    Proto_Particle * part = new Proto_Particle(m_hads[i],m_moms[i],false);
    p_hadrons->push_back(part);
  }
  return true;
}

double Soft_Cluster_Handler::RadiationWeight() {
  // no radiation for diquark-diquark clusters -- must annihilate
  if (m_flavs.first.IsDiQuark() && m_flavs.second.IsDiQuark())
    return Annihilation();
  Single_Transition_List * radiations = (*p_singletransitions)[m_flavs];
  // this should ** NEVER ** happen ..... unless stable BSM coloured particles
  if (radiations==NULL) return 0.;
  m_hads[0] = (--radiations->end())->first;
  // everything is fine - get on with your life and just decay.
  map<Flavour,double> weights;
  double totweight(0.), weight;
  for (Single_Transition_List::reverse_iterator sit=radiations->rbegin();
       sit!=radiations->rend();sit++) {
    double m2(sit->first.Mass());
    if (m2>m_mass) break;
    // wave-function overlap * phase-space (units of 1 in total)
    weight     = sit->second * PhaseSpace(m2,0.);
    totweight += weights[sit->first] = weight;
  }
  double disc = totweight * ran->Get();
  map<Flavour,double>::iterator wit=weights.begin();
  do {
    disc -= wit->second;
    if (disc<=1.e-12) break;
    wit++;
  } while (wit!=weights.end());
  if (wit!=weights.end()) m_hads[0] = wit->first;
  return totweight;
}

double Soft_Cluster_Handler::DecayWeight() {
  Double_Transition_List * decays = (*p_doubletransitions)[m_flavs];
  // this should ** NEVER ** happen ..... unless stable BSM coloured particles
  if (decays==NULL) {
    msg_Error()<<"No decays found for "
	       <<m_flavs.first<<"/"<<m_flavs.second<<".\n";
    return 0.;
  }
  // lightest possible pair of hadrons.
  m_hads[0] = (--decays->end())->first.first;
  m_hads[1] = (--decays->end())->first.second;

  // last resort: if cluster is light, but consists of two diquarks -
  // may have to "cut open" the diquarks and form two mesons out of
  // two quarks and two anti-quarks.
  if (m_hads[0].Mass()+m_hads[1].Mass()>m_mass) return Annihilation();

  // everything is fine - get on with your life and just decay.
  map<Flavour_Pair,double> weights;
  double totweight(0.), weight;
  for (Double_Transition_List::reverse_iterator dit=decays->rbegin();
       dit!=decays->rend();dit++) {
    double m2(dit->first.first.Mass()), m3(dit->first.second.Mass());
    if (m2+m3>m_mass) break;
    // wave-function overlap * phase-space (units of 1 in total)
    weight     = dit->second * PhaseSpace(m2,m3);
    totweight += weights[dit->first] = weight;
  }

  double disc = totweight * ran->Get();
  map<Flavour_Pair,double>::iterator wit=weights.begin();
  do {
    disc -= wit->second;
    if (disc<=1.e-12) break;
    wit++;
  } while (wit!=weights.end());
  if (wit!=weights.end()) {
    m_hads[0] = wit->first.first;
    m_hads[1] = wit->first.second;
  }
  return totweight;
}

double Soft_Cluster_Handler::Annihilation() {
  Flavour_Pair one, two;
  Flavour one1, one2, two1, two2;
  if (!(DiQuarkToQuarks(m_flavs.first,one1,one2) &&
	DiQuarkToQuarks(m_flavs.second,two1,two2))) return 0.;
  bool disc(ran->Get()>0.5);
  one.first  = disc?two1:two2;
  one.second = one1;
  two.first  = disc?two2:two1;
  two.second = one2;
  if (DefineHadronsInAnnihilation(one,two)) return true;
  msg_Error()<<METHOD<<" yields error - no annihilation defined.\n"
	     <<"   Will return false and hope for the best.\n";
  return false;
}

bool Soft_Cluster_Handler::
DiQuarkToQuarks(const Flavour & di,Flavour & q1,Flavour & q2) {
  if (!di.IsDiQuark()) return false;
  int kfdi(int(di.Kfcode())), kf1(kfdi/1000), kf2((kfdi-kf1*1000)/100);
  q1 = Flavour(kf1);
  q2 = Flavour(kf2);
  if (di.IsAnti()) { q1 = q1.Bar(); q2 = q2.Bar(); }
  return true;
}

double Soft_Cluster_Handler::
DefineHadronsInAnnihilation(const Flavour_Pair & one,const Flavour_Pair & two) {
  Single_Transition_List * ones = (*p_singletransitions)[one];
  Single_Transition_List * twos = (*p_singletransitions)[two];
  map<Flavour_Pair,double> weights;
  double m2, m3, totweight(0.), weight;
  for (Single_Transition_List::reverse_iterator oit=ones->rbegin();
       oit!=ones->rend();oit++) {
    m2 = oit->first.Mass();
    if (m2>m_mass) break;
    for (Single_Transition_List::reverse_iterator tit=twos->rbegin();
       tit!=twos->rend();tit++) {
      m3 = tit->first.Mass();
      if (m2+m3>m_mass) break;
      // wave-function overlap * phase-space (units of 1 in total)
      weight     = oit->second * tit->second * PhaseSpace(m2,m3);
      Flavour_Pair flpair;
      flpair.first = oit->first; flpair.second = tit->first;
      totweight += weights[flpair] = weight;
    }
  }
  double disc = totweight*ran->Get()*0.9999999;
  map<Flavour_Pair,double>::iterator wit=weights.begin();
  while (disc>0. && wit!=weights.end()) {
    disc-=wit->second;
    wit++;
  }
  // extra safety net
  if (wit==weights.end()) wit = weights.begin();
  m_hads[0] = wit->first.first;
  m_hads[1] = wit->first.second;
  return totweight;
}
  
double Soft_Cluster_Handler::PhaseSpace(const double & m2,const double & m3) {
  double m22(m2*m2),m32(m3*m3);
  return sqrt(sqr(m_mass2-m22-m32)-4.*m22*m32)/(8.*M_PI*m_mass2);
  // extra weight to possible steer away from phase space only ... may give
  // preference to higher or lower mass pairs
  // * pow((m2+m3)/m_mass,m_chi);
}

double Soft_Cluster_Handler::
MinSingleMass(const Flavour & fl1,const Flavour & fl2) {
  m_flavs.first  = fl1;
  m_flavs.second = fl2;
  return p_singletransitions->GetLightestMass(m_flavs);
}

ATOOLS::Flavour Soft_Cluster_Handler::
LowestTransition(const ATOOLS::Flavour & fl1,const ATOOLS::Flavour & fl2) {
  m_flavs.first  = fl1;
  m_flavs.second = fl2;
  return p_singletransitions->GetLightestTransition(m_flavs);
}

double Soft_Cluster_Handler::
MinDoubleMass(const Flavour & fl1,const Flavour & fl2) {
  m_flavs.first  = fl1;
  m_flavs.second = fl2;
  return p_doubletransitions->GetLightestMass(m_flavs);
}

void Soft_Cluster_Handler::PrintHadrons() {
  for (list<Proto_Particle *>::iterator hit=p_hadrons->begin();
       hit!=p_hadrons->end();hit++) {
    Proto_Particle * hadron = (*hit);
  }
}
