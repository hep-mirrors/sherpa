#include "SHRiMPS/Ladders/Ladder_Generator_QT.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <list>

using namespace SHRIMPS;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Ladder_Generator_QT::Ladder_Generator_QT() :
  Ladder_Generator_Base(),
  m_partonic(Sigma_Partonic(xs_mode::Regge)),
  m_S(sqr(rpa->gen.Ecms())),
  m_qt2min(m_Q02), m_qt2minFF(0.), m_shatmin(m_qt2min),
  m_weight(0.),
  m_fixflavour(true)
{
  for (size_t i=0;i<2;i++) m_Ebeam[i] = rpa->gen.PBeam(i)[0];
}

Ladder * Ladder_Generator_QT::operator()(const Vec4D & pos) {
  m_weight = 0.;
  p_ladder = new Ladder(pos);
  if (FixInitialPartons() &&
      MakeTrialLadder()) {
    ConstructISKinematics();
    SelectPropagatorColours();
    CalculateWeight();
  }
  else { delete p_ladder; p_ladder = NULL; }
  return p_ladder;
}

bool Ladder_Generator_QT::FixInitialPartons() {
  m_shat = m_partonic.MakeEvent(m_fixflavour);
  if (m_shat<0.) return false;
  for (size_t beam=0;beam<2;beam++) {
    m_ylimits[beam] = (beam==0 ? 1.: -1.) * (m_Ymax + ran->Get()*m_deltaY);
    m_q[beam]       = m_partonic.X(beam) * rpa->gen.PBeam(beam);
    m_flavs[beam]   = m_partonic.Flav(beam);
  }
  return true;
}

bool Ladder_Generator_QT::MakeTrialLadder() {
  for (size_t beam=0;beam<2;beam++) {
    m_ynew[beam]   = m_yold[beam] = m_ylimits[beam];
    m_qt2old[beam] = m_q[beam].PPerp2();  
  }
  TPropList::iterator pit = p_ladder->GetProps()->begin(); pit++;
  size_t dir;
  do {
    m_seff = (m_q[0]+m_q[1]).Abs2();
    if (m_seff<4.*m_qt2min) break;
    dir    = dabs(m_yold[0])>dabs(m_yold[1]) ? 0 : 1;
    if (TrialEmission(dir)) AddEmission(dir,pit);
  } while (m_ynew[0]>m_ynew[1]);
  if (LastEmissions()) {
    for (size_t dir=0;dir<2;dir++) p_ladder->AddRapidity(m_yold[dir],m_flavs[dir],m_k[dir]);
    p_ladder->GetProps()->insert(pit, T_Prop(colour_type::octet, m_qT, m_Q02));
    return true;
  }
  return false;
}

void Ladder_Generator_QT::AddEmission(size_t dir, TPropList::iterator & pit) {
  p_ladder->AddRapidity(m_yold[dir],m_flavs[dir],m_k[dir]);
  p_ladder->GetProps()->insert(pit, T_Prop(colour_type::octet, m_qT, m_Q02));
  if (dir==1) pit--;
  m_q[dir]     -= m_k[dir];	
  m_yold[dir]   = m_ynew[dir];
  m_qt2old[dir] = m_qt2; 
  if (dabs(m_ynew[dir])<m_Ymax) { m_ynew[dir] == (dir==0 ? m_Ymax : -m_Ymax); }
}

bool Ladder_Generator_QT::TrialEmission(size_t dir) {
  double qt2min    = QT2Min(dir), weight;
  double arg       = M_PI/(3.*AlphaS(0.)*log(m_seff/m_qt2min));
  Vec4D  eqT       = MakeQTNorm();
  Form_Factor * ff = dabs(m_yold[dir])>m_Ymax ? p_eikonal->FF(dir) : NULL;
  do {
    m_ynew[dir] += (dir==1? -1. : 1.) * arg * log(ran->Get());
    m_qT     = MakePropMomentum(qt2min, m_seff, eqT, ff);
    m_k[dir] = MakeFSMomentum(dir);
    weight   = ( ReggeWeight(m_qt2, dabs(m_ynew[dir]-m_yold[dir])) *
		 LDCWeight(m_qt2, m_qt2old[dir]) *
		 EmissionWeight(m_k[dir].PPerp2(), ff) *
		 AbsorptionWeight(m_k[dir], m_ynew[dir]) );
    if ((dir==0 && m_ynew[0]<m_ynew[1]) || (dir==1 && m_ynew[1]>m_ynew[0])) {
      return false;
    }
  } while (m_q[dir][0]-m_k[dir][0]<0. ||
	   weight<ran->Get());
  return true;
}

bool Ladder_Generator_QT::LastEmissions() {
  if (p_ladder->GetEmissions()->size()==0) return FixSimpleKinematics();
  double qt2min  = QT2Min(), qt2max = QT2Max(), weight;
  if (qt2max<qt2min) return false;
  Vec4D eqT      = MakeQTNorm();
  long int trials = 10000;
  Form_Factor * ff[2];
  for (size_t beam=0;beam<2;beam++) {
    ff[beam] = dabs(m_yold[beam])>m_Ymax ? p_eikonal->FF(beam) : NULL;
  }
  do {
    m_qT = MakePropMomentum(qt2min, qt2max, eqT,
			    dabs(m_yold[0])>m_Ymax ? ff[0] :
			    dabs(m_yold[1])>m_Ymax ? ff[1] :
			    NULL);
    double weight = ReggeWeight(m_qt2, dabs(m_yold[0]-m_yold[1]));
    for (size_t beam=0;beam<2;beam++) {
      m_k[beam] = MakeFSMomentum(beam);
      weight   *= LDCWeight(m_qt2, m_qt2old[beam], dabs(m_yold[1-beam])<m_Ymax);
      weight   *= EmissionWeight(m_k[beam].PPerp2(), ff[beam]);
    }
    if ((m_k[0]+m_k[1]).Abs2() < m_seff && weight > ran->Get()) return true;
  } while ((trials--)>0); 
  return false;
}

bool Ladder_Generator_QT::FixSimpleKinematics() {
  double qt2max = m_shat/cosh(Max(m_ylimits[0],m_ylimits[1]));
  double factor = (sqr(cosh(m_ylimits[0])+cosh(m_ylimits[1]))-
		   sqr(sinh(m_ylimits[0])+sinh(m_ylimits[1])));
  size_t trials = 0;
  do {
    m_qt2 = sqrt(p_eikonal->FF(0)->SelectQT2(qt2max)*p_eikonal->FF(0)->SelectQT2(qt2max));
    if ((trials++)>1000000) return false;
  } while (m_qt2*factor>m_shat);
  double qt = sqrt(m_qt2), phi = 2.*M_PI*ran->Get();
  m_qT      = Vec4D(0.,cos(phi),sin(phi),0.);
  for (size_t i=0;i<2;i++) {
    m_k[i]  = qt * (Vec4D(cosh(m_ylimits[i]),0,0,sinh(m_ylimits[i])) +
		    (i==0?1.:-1.) * m_qT);
  }
  return true;
}

double Ladder_Generator_QT::LDCWeight(const double & qt2,const double & qt2prev,const bool & apply) {
  return  apply ? qt2 / Max(qt2,qt2prev) : 1.;
}

double Ladder_Generator_QT::EmissionWeight(const double & kt2, Form_Factor * ff) {
  return (ff ? 1. : AlphaS(kt2)/AlphaS(0.) );	   
}

double Ladder_Generator_QT::AbsorptionWeight(const Vec4D & k,const double & y) {
  //double kt2 = k.PPerp2();
  // double wt = kt2/(kt2+m_qt2min);
  return m_density.AbsorptionWeight(y);
}

double Ladder_Generator_QT::ReggeWeight(const double & q2,const double & deltay) {
  return exp(-3.*AlphaS(dabs(q2))/M_PI * log(q2/m_qt2min) * dabs(deltay));
}


void Ladder_Generator_QT::SelectPropagatorColours() {
  // Iterate over propagators and assign colours different than octet
  LadderMap::iterator lit1=p_ladder->GetEmissions()->begin(), lit2 = lit1; lit2++;
  TPropList::iterator pit1=p_ladder->GetProps()->begin(), pit2=pit1;
  double y1,y2,wt1,wt8,ratio1,ratio2=0.;
  while (lit2!=p_ladder->GetEmissions()->end() && pit1!=p_ladder->GetProps()->end()) {
    y1     = lit1->first;
    y2     = lit2->first;
    wt1    = m_density.SingletWeight(y2,y1);
    wt8    = m_density.OctetWeight(y2,y1);
    ratio1 = wt1/(wt1+wt8); 
    if (ratio1>ran->Get()) pit1->SetCol(colour_type::singlet);
    if (pit1!=p_ladder->GetProps()->begin() &&
	pit1->Col()==colour_type::singlet &&
	pit2->Col()==colour_type::singlet) {
      if (ratio1>ratio2) { pit2->SetCol(colour_type::octet); }
      else               { pit1->SetCol(colour_type::octet); }
    }
    ratio2 = ratio1;
    pit2   = pit1;
    lit1++;lit2++;pit1++;
  }
}

void Ladder_Generator_QT::CalculateWeight() {
  double qt2max(m_Q02), tmax(m_Q02);
  for (TPropList::iterator pit=p_ladder->GetProps()->begin();
       pit!=p_ladder->GetProps()->end();pit++) {
    double qt2 = pit->Q().PPerp2();
    if (qt2>qt2max) {
      qt2max = qt2;  tmax = Max(tmax,dabs(pit->Q().Abs2()));
    }
  }
  m_weight = qt2max*m_Q02/sqr(tmax);
}


