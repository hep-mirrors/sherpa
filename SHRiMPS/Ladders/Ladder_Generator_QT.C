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
  Ladder_Generator_Base(), m_seff(4.)
{
  for (size_t beam=0;beam<2;beam++)
    m_ylimits[beam] = (beam==0 ? 1.:-1.) * (m_Ymax-m_deltaY);
}

Ladder * Ladder_Generator_QT::operator()(const Vec4D & pos) {
  InitLadder(pos);
  do {
    p_props->clear();
    p_emissions->clear();
    CreateInitialMomenta();
    if (!CreateLadder()) { delete p_ladder; return NULL; }
    CalculateWeight();
  } while (m_weight<ran->Get());
  m_colourgenerator(p_ladder);
  msg_Out()<<METHOD<<" constructs ladder (weight = "<<m_weight<<")\n"<<(*p_ladder);
  return p_ladder;
}

void Ladder_Generator_QT::CreateInitialMomenta() {
  double x_ini[2], s_ini;
  do {
    for (size_t beam=0;beam<2;beam++) {
      x_ini[beam]  = m_partonic.GetPDF(beam)->SelectX(Flavour(kf_gluon),m_seff);
      m_qini[beam] = x_ini[beam] * m_E[beam] * Vec4D(1.,0,0,(beam==0)?1.:-1.);
    }
    m_cms  = m_qini[0]+m_qini[1];
    m_sact = m_cms.Abs2();
  } while (m_sact<m_seff);
  for (size_t beam=0;beam<2;beam++) {
    m_y[beam] = m_ylimits[beam];
    m_q[beam] = m_qini[beam];
    msg_Out()<<"Selected x for "<<beam<<": "
    	     <<x_ini[beam]<<" --> "<<m_qini[beam]<<", y = "<<m_y[beam]<<"\n";
  }
}

bool Ladder_Generator_QT::CreateLadder() {
  do { } while (YStep() && m_sact>m_seff);
  CreatePropagators();
  MergeSinglets();
  if (p_emissions->size()<2) return false;
  m_weight = RescaleLadder(p_ladder, m_cms);
  if (!(p_ladder->FSMomentum()-(m_qini[0]+m_qini[1])).IsZero()) return false;
  for (size_t beam=0;beam<2;beam++) p_ladder->InPart(beam)->SetMomentum(m_qini[beam]);
  return true;
}

bool Ladder_Generator_QT::YStep() {
  size_t step  = dabs(m_y[0])<dabs(m_y[1]);
  double arg   = 3./M_PI*AlphaSMax()*log(1.+m_sact/(4.*m_kt2min));
  double sud   = exp(-arg*dabs(m_y[step]-m_y[1-step]));
  msg_Out()<<METHOD<<": step = "<<step<<", arg = "<<arg<<", sud = "<<sud<<"\n";
  if (sud>ran->Get()) {
    msg_Out()<<"No Sudakov left for splitting.\n";
    m_y[step] = m_y[1-step];
    return false;
  }
  double y, kt2, wt, phi, s;
  Vec4D  q, k;
  bool   success = true;
  size_t trials  = 0;
  do {
    y   = m_y[step] + (step==0 ? 1.:-1.) * log(ran->Get())/arg;
    kt2 = m_kt2min * ( pow(1.+m_sact/(4.*m_kt2min),ran->Get()) - 1.);
    msg_Out()<<"- y = "<<y<<", kt2 = "<<kt2<<" vs "<<(m_kt2min/cosh(y))<<"s = "<<s<<"\n";
    if (kt2<m_kt2min/cosh(y)) { trials++; continue; }
    phi = 2.*M_PI*ran->Get();
    k   = sqrt(kt2)*Vec4D(cosh(y),cos(phi),sin(phi),sinh(y));
    q   = m_q[step]-k;
    s   = (q+m_q[1-step]).Abs2();
    if (y<m_y[1] || y>m_y[0] || trials++>1000) { success = false; break; }
    wt  = AlphaS(kt2+m_kt2min)/AlphaSMax();
    wt *= ReggeWeight(dabs(m_q[step].Abs2()),y,m_y[step]);
    wt *= m_density.AbsorptionWeight(y);
  } while (s<m_seff || wt<ran->Get());
  if (success) {
    m_q[step] = q;
    m_sact    = s;
    p_ladder->AddRapidity(y,Flavour(kf_gluon),k);
  }
  //msg_Out()<<"   "<<m_y[step]<<" --> m_y["<<step<<"] = "<<y<<", kt^2 = "<<kt2<<" "
  //	   <<" --> "<<k<<"\n"
  //	   <<"   new s = "<<m_sact<<" ("<<m_q[0]<<"+"<<m_q[1]<<")\n";
  m_y[step] = y;
  return success;				  
}

void Ladder_Generator_QT::CreatePropagators() {
  Vec4D q = m_qini[0];
  LadderMap::iterator eit1=p_ladder->GetEmissions()->begin(), eit2=eit1; eit2++;
  do {
    q          -= eit1->second.Momentum();
    double y1   = eit1->first, y2 = eit2->first;
    double wt1  = m_density.SingletWeight(y1,y2), wt8 = m_density.OctetWeight(y1,y2);
    T_Prop prop = T_Prop((wt1/(wt1+wt8)>ran->Get()?colour_type::singlet:colour_type::octet),q,m_qt2min);
    p_ladder->AddPropagator(prop);
    eit1++; eit2++;
  } while (eit2!=p_ladder->GetEmissions()->end());
}

void Ladder_Generator_QT::MergeSinglets() {
  TPropList::iterator pit1=p_ladder->GetProps()->begin(), pit2=pit1, pithelp; pit2++;
  LadderMap::iterator lit=p_ladder->GetEmissions()->begin(), lithelp; lit++;
  bool hit = false;
  do {
    if (pit1->Col()==colour_type::singlet && pit2->Col()==colour_type::singlet) {
      lithelp = lit;  lit++;
      pithelp = pit2; pit2++;
      p_ladder->DeleteRapidity(lithelp);
      p_ladder->DeletePropagator(pithelp);
      hit = true;
    }
    else { pit1++; pit2++; lit++; }
  } while (pit2!=p_ladder->GetProps()->end());
}

void Ladder_Generator_QT::CalculateWeight() {}
