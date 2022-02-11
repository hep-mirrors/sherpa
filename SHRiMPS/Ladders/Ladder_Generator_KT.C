#include "SHRiMPS/Ladders/Ladder_Generator_KT.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <list>

using namespace SHRIMPS;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Ladder_Generator_KT::Ladder_Generator_KT() : Ladder_Generator_Base() {}

Ladder_Generator_KT::~Ladder_Generator_KT() {}

Ladder * Ladder_Generator_KT::operator()(const Vec4D & pos) {
  InitLadder(pos);
  if (MakeTrialInitialPartons() && MakeTrialLadder()) {
    ConstructISKinematics();
    SelectPropagatorColours();
  }
  else { delete p_ladder; p_ladder = NULL; }
  return p_ladder;
}

bool Ladder_Generator_KT::MakeTrialInitialPartons() {
  //msg_Out()<<METHOD<<": E1, 2 = "<<m_E[0]<<", "<<m_E[1]
  //	   <<" --> smax = "<<(4.*m_E[0]*m_E[1])<<"\n";
  do {
    m_shat     = m_partonic.MakeEvent();
  } while (m_shat>(4.*m_E[0]*m_E[1]));
  if (m_shat<0.) return false;
  m_sumQ     = Vec4D(0.,0.,0.,0.);
  m_sigmahat = m_partonic.SigmaHat();
  for (size_t beam=0;beam<2;beam++) {
    m_ylimits[beam] = (beam==0 ? 1.: -1.) * (m_Ymax + ran->Get()*m_deltaY);
    m_sumQ         += m_qini[beam] = 2.*m_partonic.X(beam) * rpa->gen.PBeam(beam);
    m_flavs[beam]   = m_partonic.Flav(beam);
  }
  return true;
}

bool Ladder_Generator_KT::MakeTrialLadder() {
  double y = 0.5*log(m_shat/m_qt2min), yhat = m_partonic.YHat();
  size_t trials = 0;
  size_t N = m_density.NGluons(-y+yhat, y+yhat);
  do {
    p_ladder->ResetFS();
    if ((trials++)>1000) return false;
    MakeTrialRapidities(N);
    if (SelectKTs()) {
      FillRapidities();
      ConstructISKinematics();
      ConstructPropagators();
      SelectPropagatorColours();
      CalculateWeight();
    }
    else m_weight = 0.;
  } while(m_weight<ran->Get());
  return true;
}

void Ladder_Generator_KT::MakeTrialRapidities(const size_t & N) {
  m_rapidities.clear();
  for (size_t beam=0;beam<2;beam++) {
    m_rapidities[m_ylimits[beam]] = Vec4D(0.,0.,0.,0.);
  }
  for (size_t i=0;i<N;i++) {
    double y = m_density.SelectRapidity(m_ylimits[1], m_ylimits[0]);
    m_rapidities[y] = Vec4D(0.,0.,0.,0.);
  }
}

bool Ladder_Generator_KT::SelectKTs() {
  size_t trials = 0;
  double kt2, kt, weight;
  do {
    if (trials++>1000) return false;
    m_sumK = Vec4D(0.,0.,0.,0.);
    for (map<double,Vec4D>::iterator rit=m_rapidities.begin();
	 rit!=m_rapidities.end();rit++) {
      if (rit->first!=m_rapidities.rbegin()->first) {
	m_sumK += rit->second = SelectKT(rit->first,dabs((m_sumQ-m_sumK).Abs2()));
      }
    }
    kt2    = m_sumK.PPerp(); kt = sqrt(kt2);
    weight = (kt2>m_kt2min ?
	      AlphaSWeight(kt2)/kt2 :
	      m_rapidities.size()<4 ? 1. : 0.);
  } while (weight<ran->Get());
  double y = m_rapidities.rbegin()->first;
  m_sumK += m_rapidities.rbegin()->second =
    kt * Vec4D(cosh(y),0.,0.,sinh(y)) - m_sumK.Perp();
  return true;
}

Vec4D Ladder_Generator_KT::SelectKT(const double & y,const double & seff) {
  double kt2max = seff/sqr(2.*cosh(y));
  double kt2min = (dabs(y)>m_Ymax?m_qt2minFF:m_kt2min);
  if (kt2max<kt2min) return Vec4D(0.,0.,0.,0.);
  double kt2 = 0.;
  MakeTransverseUnitVector();
  if (y>=m_Ymax)       kt2 = p_eikonal->FF(0)->SelectQT2(kt2max,m_qt2minFF);
  else if (y<=-m_Ymax) kt2 = p_eikonal->FF(1)->SelectQT2(kt2max,m_qt2minFF);
  else  {
    do {
      kt2 = m_kt2min * pow(kt2max/m_kt2min, ATOOLS::ran->Get());
    } while (AlphaSWeight(kt2)<ran->Get());
  }
  return sqrt(kt2) * (Vec4D(cosh(y),0.,0.,sinh(y)) + m_eqt);
}

void Ladder_Generator_KT::FillRapidities() {
  for (map<double,Vec4D>::iterator rit=m_rapidities.begin();
       rit!=m_rapidities.end();rit++) {
    p_ladder->AddRapidity(rit->first,Flavour(kf_gluon),rit->second);
  }
}

void Ladder_Generator_KT::ConstructPropagators() {
  if (p_emissions->empty()) return;
  for (size_t i=0;i<p_emissions->size()-1;i++)
    p_props->push_back(T_Prop(colour_type::octet,Vec4D(0.,0.,0.,0.),m_qt2min));
  
  TPropList::iterator pit = p_props->begin();
  Vec4D qin = p_ladder->InPart(0)->Momentum();
  for (LadderMap::iterator lit=p_emissions->begin();
       lit!=(--p_emissions->end());lit++) {
    qin -= lit->second.Momentum();
    pit->SetQ(qin[0]<0. ? -qin: qin);
    pit->SetQ2(qin.Abs2());
    pit->SetQT2(qin.PPerp2());
    pit->SetQ02(m_qt2min);
    pit++; 
  }
  //msg_Out()<<"check : "<<(p_ladder->InPart(1)->Momentum()+qin)<<" = "
  //	   <<(--p_emissions->end())->second.Momentum()<<"\n";
}

void Ladder_Generator_KT::SelectPropagatorColours() {
  LadderMap::iterator lit1=p_emissions->begin(),  lit2=p_emissions->end();  lit2--;
  TPropList::iterator pit1=p_props->begin(), pit2=p_props->end(); pit2--;
  double y1, y2, wt1, wt8;
  size_t dir;
  while (lit1->first>lit2->first) {
    dir = (dabs(lit1->first) > dabs(lit2->first))?1:0;
    if (dir) { y1 = lit1->first; lit1++; y2 = lit1->first; }
    else     { y2 = lit2->first; lit2--; y1 = lit2->first; }
    wt1 = m_density.SingletWeight(y1,y2);
    wt8 = m_density.OctetWeight(y1,y2);
    if (wt1/(wt1+wt8)>ran->Get()) {
      if (dir) { pit1->SetCol(colour_type::singlet); pit1++; }
      else     { pit2->SetCol(colour_type::singlet); pit2--; }
    }
  }
  pit1 = p_props->begin(); pit2 = pit1; pit2++;
  while (pit2!=p_props->end()) {
    if (pit1->Col()==colour_type::singlet && pit2->Col()==colour_type::singlet) {
      if (ran->Get()>0.5) pit1->SetCol(colour_type::octet); 
      else pit2->SetCol(colour_type::octet);
    }
    pit1++; pit2++;
  }
}

void Ladder_Generator_KT::CalculateWeight() {
  Vec4D Pcms       = (p_ladder->InPart(0)->Momentum() +
		      p_ladder->InPart(1)->Momentum()); 
  double Y         = Pcms.Y(), SHat = Pcms.Abs2();
  double sigma_act = m_partonic.dSigma(SHat,Y);
  double sigma_ratio, me, regge;
  m_weight = ( (sigma_ratio = sigma_act/m_sigmahat) *
	       (me          = m_me(p_ladder,m_qt2min)) *
	       (regge       = TotalReggeWeight(p_ladder)) );
}

