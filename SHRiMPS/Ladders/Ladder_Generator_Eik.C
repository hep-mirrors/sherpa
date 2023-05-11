//emission k_t's can become small, maybe such emissions need to be kicked out
#include "SHRiMPS/Ladders/Ladder_Generator_Eik.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <list>

using namespace SHRIMPS;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Ladder_Generator_Eik::Ladder_Generator_Eik() : Ladder_Generator_Base() {}

Ladder_Generator_Eik::~Ladder_Generator_Eik() {}

Ladder * Ladder_Generator_Eik::operator()(const Vec4D & pos,Sigma_Elastic * sigma_el,Sigma_D * sigma_sd) {
  InitLadder(pos);
  FillGluons();
  SelectPropagatorColours();
  if (p_ladder->GetEmissions()->size()==2) {
    ConstructSimpleLadder();
    ConstructISKinematics();
  }
  else {
    size_t trials = 0;
    do {
      if ((trials++)>1000) { delete p_ladder; return NULL; }
      if (SelectPropagatorQTs()) {
	ConstructISKinematics();
	CalculateWeight();
      }
      else m_weight = 0.;
    } while (m_weight<ran->Get());
  }
  //msg_Out() << *p_ladder << endl;
  return p_ladder;
}

void Ladder_Generator_Eik::FillGluons() {
  for (size_t beam=0;beam<2;beam++) {
    m_ylimits[beam] = (beam==0? 1. : -1.) * (m_Ymax + ran->Get()*m_deltaY);
    p_ladder->AddRapidity(m_ylimits[beam]);
  }
  size_t N = 0; //m_density.NGluons(m_ylimits[1], m_ylimits[0]);
  for (size_t i=0;i<N;i++) {
    p_ladder->AddRapidity(m_density.SelectRapidity(m_ylimits[1], m_ylimits[0]));
  }
}

void Ladder_Generator_Eik::SelectPropagatorColours() {
  for (size_t i=0;i<p_emissions->size()-1;i++)
    //    p_props->push_back(T_Prop(colour_type::octet,Vec4D(0.,0.,0.,0.),m_qt2min));
    p_props->push_back(T_Prop(colour_type::singlet,Vec4D(0.,0.,0.,0.),m_qt2min));
  return;
  
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

bool Ladder_Generator_Eik::SelectPropagatorQTs() {
  LadderMap::iterator flit[2];
  LadderMap::reverse_iterator rlit[2];
  TPropList::iterator         fpit=p_ladder->GetProps()->begin();
  TPropList::reverse_iterator rpit=p_ladder->GetProps()->rbegin();
  for (size_t i=0;i<2;i++) {
    flit[i] = p_ladder->GetEmissions()->begin(); 
    rlit[i] = p_ladder->GetEmissions()->rbegin(); 
    if (i==1) { flit[i]++; rlit[i]++; }
    m_y[0][i]    = flit[i]->first;
    m_y[1][i]    = rlit[i]->first;
    m_qt2prev[i] = 0.;
    m_qtprev[i]  = Vec4D(0.,0.,0.,0.);
  }
  size_t lastdir;
  do {
    size_t dir=dabs(m_y[0][0])>dabs(m_y[1][0])?0:1;
    if (!SelectPropagatorQT(dir,dir==0?(*fpit):(*rpit))) return false;
    Vec4D kt = m_qtprev[dir]-m_qt;
    Vec4D k  = kt.PPerp()*Vec4D(cosh(m_y[dir][0]),0.,0.,sinh(m_y[dir][0]))-kt;
    if (dir==0) {
      flit[0]->second.SetMomentum(k);
      flit[0] = flit[1]; flit[1]++; fpit++;
      for (size_t i=0;i<2;i++) m_y[dir][i] = flit[i]->first;
    }
    else {
      rlit[0]->second.SetMomentum(k);
      rlit[0] = rlit[1]; rlit[1]++; rpit++;
      for (size_t i=0;i<2;i++) m_y[dir][i] = rlit[i]->first;
    }
    m_qtprev[dir] = m_qt;
    lastdir       = dir;
  } while (m_y[0][0]>m_y[1][0]);
  Vec4D kt   = m_qtprev[lastdir]+m_qtprev[1-lastdir];
  m_lastk    = kt.PPerp()*Vec4D(cosh(m_y[1-lastdir][0]),0.,0.,sinh(m_y[1-lastdir][0]))-kt;
  if (lastdir==0) rlit[0]->second.SetMomentum(m_lastk);
  else            flit[0]->second.SetMomentum(m_lastk);
  return true;
}

bool Ladder_Generator_Eik::SelectPropagatorQT(const size_t dir,T_Prop & prop) {
  double qt2max = m_shat/(4.*sqr(cosh(m_y[dir][0]))), weight;
  double qt2min = m_qt2min/sqr(cosh(m_y[dir][0]));
  size_t trials = 0;
  m_qt2 = 0.;
  MakeTransverseUnitVector();
  do {
    if (trials++ > 1000) return false;
    if (dabs(m_y[dir][0])>m_Ymax) {
      m_qt2  = p_eikonal->FF(m_y[dir][0]>0.?0:1)->SelectQT2(qt2max,0.);
      weight = 1.;
    }
    else {
      m_qt2  = (prop.Col()==colour_type::octet ?
        qt2min * pow(qt2max/qt2min,ran->Get()) : //p(qt2) = 1/qt2
        //qt2min * pow((qt2max+qt2min)/qt2min,ran->Get()) - qt2min : //p(qt2) = 1/(qt2+qt2min)
        qt2min*qt2max / (qt2max - ran->Get()*(qt2max-qt2min)) ); //p(qt2) = ln(qt2)
      weight = ( AlphaSWeight((m_qtprev[dir]-sqrt(m_qt2)*m_eqt).PPerp2()) *
         LDCWeight(m_qt2,m_qtprev[dir].PPerp2()) );
    }
    weight *= (prop.Col()==colour_type::octet ?
	       ReggeWeight(m_qt2,m_y[dir][0],m_y[1-dir][1]) :
	       sqr(ReggeWeight(m_qt2,m_y[dir][0],m_y[1-dir][1])) );
  } while (weight < ran->Get());
  prop.SetQ(m_qt = sqrt(m_qt2)*m_eqt);
  prop.SetQT2(m_qt2);
  return true;
}

void Ladder_Generator_Eik::CalculateWeight() {
  m_weight  = AlphaSWeight(m_lastk.PPerp2());
  m_weight *= m_me(p_ladder,m_qt2min); 
}

