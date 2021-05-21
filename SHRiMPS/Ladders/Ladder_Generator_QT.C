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
  m_S(sqr(rpa->gen.Ecms())), m_qt2minFF(0.), m_shatmin(m_qt2min), m_seff(0.), 
  m_fixflavour(true)
{
  for (size_t i=0;i<2;i++) m_Ebeam[i] = rpa->gen.PBeam(i)[0];
}

Ladder * Ladder_Generator_QT::operator()(const Vec4D & pos) {
  InitLadder(pos);
  if (FixInitialPartons() && MakeTrialLadder()) {
    ConstructISKinematics();
    SelectPropagatorColours();
    CalculateWeight();
  }
  else { delete p_ladder; p_ladder = NULL; }
  return p_ladder;
}

bool Ladder_Generator_QT::FixInitialPartons() {
  m_shat     = m_partonic.MakeEvent(m_fixflavour);
  if (m_shat<0.) return false;
  m_sigmahat = m_partonic.SigmaHat();
  for (size_t beam=0;beam<2;beam++) {
    m_ylimits[beam] = (beam==0 ? 1.: -1.) * (m_Ymax + ran->Get()*m_deltaY);
    m_qini[beam]    = m_q[beam] = m_partonic.X(beam) * rpa->gen.PBeam(beam);
    m_flavs[beam]   = m_partonic.Flav(beam);
  }
  return true;
}

bool Ladder_Generator_QT::MakeTrialLadder() {
  for (size_t beam=0;beam<2;beam++) {
    m_y[beam][0]    = m_y[beam][1] = m_ylimits[beam];
    m_qt2prev[beam] = m_q[beam].PPerp2();  
  }
  TPropList::iterator pit = p_ladder->GetProps()->begin(); pit++;
  size_t dir;
  do {
    m_seff = (m_q[0]+m_q[1]).Abs2();
    if (m_seff<4.*m_qt2min) break;
    dir    = dabs(m_y[0][1])>dabs(m_y[1][1]) ? 0 : 1;
    if (TrialEmission(dir)) AddEmission(dir,pit);
  } while (m_y[0][0]>m_y[1][0]);
  if (LastEmissions()) {
    for (size_t dir=0;dir<2;dir++)
      p_ladder->AddRapidity(m_y[dir][1],m_flavs[dir],m_k[dir]);
    p_ladder->GetProps()->insert(pit,
				 T_Prop(colour_type::octet, m_qT, m_qt2min));
    return true;
  }
  return false;
}

void Ladder_Generator_QT::AddEmission(size_t dir, TPropList::iterator & pit) {
  p_ladder->AddRapidity(m_y[dir][1],m_flavs[dir],m_k[dir]);
  p_ladder->GetProps()->insert(pit,
			       T_Prop(colour_type::octet, m_qT, m_qt2min));
  if (dir==1) pit--;
  m_q[dir]      -= m_k[dir];	
  m_y[dir][1]    = m_y[dir][0];
  m_qt2prev[dir] = m_qt2; 
  if (dabs(m_y[dir][0])<m_Ymax) { m_y[dir][0] == (dir==0 ? m_Ymax : -m_Ymax); }
}

bool Ladder_Generator_QT::TrialEmission(size_t dir) {
  MakeTransverseUnitVector();
  double qt2min    = QT2Min(dir), weight, dy;
  double arg       = M_PI/(3.*AlphaS(0.)*log(m_seff/m_qt2min));
  Form_Factor * ff = dabs(m_y[dir][1])>m_Ymax ? p_eikonal->FF(dir) : NULL;
  do {
    dy           = arg * log(ran->Get());
    m_y[dir][0] += (dir==1? -1. : 1.) * dy;
    m_qT         = MakePropMomentum(qt2min, m_seff, ff);
    m_k[dir]     = MakeFSMomentum(dir);
    weight       = ( ReggeWeight(m_qt2, m_y[dir][0], m_y[dir][1]) *
		     LDCWeight(m_qt2, m_qt2prev[dir]) *
		     AlphaSWeight(m_k[dir].PPerp2()) *
		     AbsorptionWeight(m_k[dir], m_y[dir][0]) );
    if ((dir==0 && m_y[0][0]<m_y[1][0]) || (dir==1 && m_y[1][0]>m_y[0][0])) {
      return false;
    }
  } while (m_q[dir][0]-m_k[dir][0]<0. || weight<ran->Get());
  return true;
}

bool Ladder_Generator_QT::LastEmissions() {
  if (p_ladder->GetEmissions()->size()==0) return FixSimpleKinematics();
  double qt2min  = QT2Min(), qt2max = QT2Max();
  if (qt2max<qt2min) return false;
  MakeTransverseUnitVector();
  //long int trials = 1000;
  Form_Factor * ff = NULL;
  if (dabs(m_y[0][1])>dabs(m_y[1][1]) && dabs(m_y[0][1])>m_Ymax)
    ff = p_eikonal->FF(0);
  if (dabs(m_y[1][1])>dabs(m_y[0][1]) && dabs(m_y[1][1])>m_Ymax)
    ff = p_eikonal->FF(1);
  size_t trials = 1000;
  do {
    m_qT = MakePropMomentum(qt2min, qt2max,ff); 
    for (size_t beam=0;beam<2;beam++) m_k[beam] = MakeFSMomentum(beam);
    if ((m_k[0]+m_k[1]).Abs2()<m_seff) {
      double weight = (ReggeWeight(m_qt2, m_y[0][1], m_y[1][1]) *
		       LDCWeight(m_qt2, m_qt2prev[0]) *
		       LDCWeight(m_qt2, m_qt2prev[1]) *
		       AlphaSWeight(m_k[0].PPerp2())  *
		       AlphaSWeight(m_k[1].PPerp2()) *
		       AbsorptionWeight(m_k[0], m_y[0][1]) *
		       AbsorptionWeight(m_k[1], m_y[1][1]));
      if (weight>ran->Get()) return true;
    }
  } while ((trials--)>0); 
  return false;
}

bool Ladder_Generator_QT::FixSimpleKinematics() {
  double qt2max = m_shat/cosh(Max(m_ylimits[0],m_ylimits[1]));
  double factor = (sqr(cosh(m_ylimits[0])+cosh(m_ylimits[1]))-
		   sqr(sinh(m_ylimits[0])+sinh(m_ylimits[1])));
  size_t trials = 0;
  do {
    m_qt2 = sqrt(p_eikonal->FF(0)->SelectQT2(qt2max,0.)*
		 p_eikonal->FF(0)->SelectQT2(qt2max,0.));
    if ((trials++)>1000) return false;
  } while (m_qt2*factor>m_shat);
  MakeTransverseUnitVector();
  double qt = sqrt(m_qt2);
  for (size_t i=0;i<2;i++)
    m_k[i]  = qt * (Vec4D(cosh(m_ylimits[i]),0,0,sinh(m_ylimits[i])) +
		    (i==0?1.:-1.) * m_eqt);
  return true;
}

double Ladder_Generator_QT::
AbsorptionWeight(const Vec4D & k,const double & y) {
  return m_density.AbsorptionWeight(y);
  //return (1.-exp(-m_kt2min/(k.PPerp2()+m_kt2min)))/(1.-exp(-1.)) *
  //  m_density.AbsorptionWeight(y);
}

void Ladder_Generator_QT::SelectPropagatorColours() {
  // Iterate over propagators and assign colours different than octet
  LadderMap::iterator lit1=p_ladder->GetEmissions()->begin(),
    lit2 = lit1; lit2++;
  TPropList::iterator pit1=p_ladder->GetProps()->begin(), pit2=pit1;
  double y1,y2,wt1,wt8,ratio1,ratio2=0.;
  while (lit2!=p_ladder->GetEmissions()->end() &&
	 pit1!=p_ladder->GetProps()->end()) {
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
  Vec4D Pcms       = (p_ladder->InPart(0)->Momentum() +
  		      p_ladder->InPart(1)->Momentum()); 
  double Y         = Pcms.Y(), SHat = Pcms.Abs2();
  double sigma_act = m_partonic.dSigma(SHat,Y);
  m_weight  = (sigma_act/m_sigmahat);
  m_weight *= TWeight();
}


