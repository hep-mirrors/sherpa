//emission k_t's can become small, maybe such emissions need to be kicked out
#include "SHRiMPS/Ladders/Ladder_Generator_Seeded.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <list>

using namespace SHRIMPS;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Ladder_Generator_Seeded::Ladder_Generator_Seeded() : Ladder_Generator_Base() {
  m_E[0] = m_E[1] = rpa->gen.Ecms()/2.;
  for (size_t beam=0;beam<2;beam++) {
    m_yseed[beam==0?0:3] = m_ylimits[beam] = (beam==0?1.:-1.)*(m_Ymax-m_deltaY);
  }
}

Ladder_Generator_Seeded::~Ladder_Generator_Seeded() {}

Ladder * Ladder_Generator_Seeded::operator()(const Vec4D & pos) {
  size_t trials = 0;
  do {
    SeedLadder(pos);
    FillIntervals();
    CompensateKTs();
    ConstructFSMomenta();
    FillPropagators();
    ConstructISKinematics();
    CalculateWeight();
    if (ran->Get()<m_weight) {
      m_colourgenerator(p_ladder);
      //msg_Out()<<"   "<<METHOD<<"["<<m_emissions[0]<<", "<<m_emissions[1]<<", "<<m_emissions[2]<<"], "
      //       <<"ktsum = "<<m_ktsum<<" for incoming E's "<<p_ladder->InPart(0)->Momentum()[0]<<" and "
      //       <<p_ladder->InPart(1)->Momentum()[0]<<" \n"<<(*p_ladder)<<"\n";
      break;
    }
    delete p_ladder; p_ladder = NULL;
  } while ((trials++)<=100);
  return p_ladder;
}

void Ladder_Generator_Seeded::SeedLadder(const Vec4D & pos) {
  p_ladder     = new Ladder(pos);
  p_emissions  = p_ladder->GetEmissions();
  p_props      = p_ladder->GetProps();
  do { m_shat  = m_partonic.MakeEvent();
  } while (m_partonic.Y(0)>m_yseed[0] || m_partonic.Y(1)<m_yseed[3]);
  m_kt2max     = m_partonic.PT2();
  m_phi        = m_partonic.Phi();
  m_ktsum      = Vec4D(0.,0.,0.,0.);
  double kt    = sqrt(m_kt2max);
  Vec4D  ktvec = kt*Vec4D(0.,cos(m_phi),sin(m_phi),0.);
  for (size_t i=0;i<2;i++) {
    m_yseed[i+1] = m_partonic.Y(i);
    m_xbeam[i]   = pow(m_partonic.X(i),ran->Get());
    m_oldxpdf[i] = m_partonic.XPDF(i);
    p_ladder->AddRapidity(m_yseed[i+1],m_partonic.Flav(i),(i==0?1.:-1.)*ktvec);
  }
  p_props->push_back(T_Prop(colour_type::octet,-ktvec,m_kt2max));
}

void Ladder_Generator_Seeded::FillIntervals() {
  for (size_t i=0;i<3;i++) {
    if (m_yseed[i]<=m_yseed[i+1]) continue;
    double wt1 = m_density.SingletWeight(m_yseed[i],m_yseed[i+1]);
    double wt8 = m_density.OctetWeight(m_yseed[i],m_yseed[i+1]);
    if (wt1/(wt1+wt8)>ran->Get()) {
      m_cols[i]      = colour_type::singlet;
      m_emissions[i] = 0;
    }
    else {
      m_cols[i]      = colour_type::octet;
      m_emissions[i] = m_density.NGluons(m_yseed[i], m_yseed[i+1]);
      Vec4D  ktvec;
      if (m_emissions[i]>0) {
	double kt2min = m_kt2min;
	for (size_t j=0;j<m_emissions[i];j++) {
	  double y      = m_density.SelectRapidity(m_yseed[i], m_yseed[i+1]);
	  double kt2max = (i==1)? m_E[0]*m_E[1]/sqr(2.*cosh(y)) : m_kt2max;
	  m_ktsum += ktvec = SelectKT(y,kt2min,kt2max);
	  p_ladder->AddRapidity(y,Flavour(kf_gluon),ktvec);
	}
      }
    }
  }
  if (m_cols[1]==colour_type::singlet) {
    if (m_cols[0]==colour_type::singlet && m_cols[2]==colour_type::singlet) m_cols[1]=colour_type::octet;
    else if (m_cols[0]==colour_type::singlet) {
      double wt01 = m_density.SingletWeight(m_yseed[0],m_yseed[1]);
      double wt12 = m_density.SingletWeight(m_yseed[1],m_yseed[2]);
      m_cols[wt12/(wt01+wt12)>ran->Get() ? 0 : 1] = colour_type::octet;
    }
    else if (m_cols[2]==colour_type::singlet) {
      double wt12 = m_density.SingletWeight(m_yseed[1],m_yseed[2]);
      double wt23 = m_density.SingletWeight(m_yseed[2],m_yseed[3]);
      m_cols[wt12/(wt12+wt23)>ran->Get() ? 2 : 1] = colour_type::octet;
    }
  }
}

void Ladder_Generator_Seeded::CompensateKTs() {
  for (size_t i=1;i<3;i++) {
    Vec4D before = (*p_emissions)[m_yseed[i]].Momentum();
    (*p_emissions)[m_yseed[i]].SetMomentum(before - m_ktsum/2.);
    Vec4D after  = (*p_emissions)[m_yseed[i]].Momentum();
  }
}

void Ladder_Generator_Seeded::ConstructFSMomenta() {
  m_ktsum = Vec4D(0.,0.,0.,0.);
  for (LadderMap::iterator pit=p_emissions->begin();pit!=p_emissions->end();pit++) {
    double y     = pit->first;
    Vec4D  ktvec = pit->second.Momentum();
    double kt    = sqrt(-ktvec.Abs2());
    pit->second.SetMomentum(ktvec+kt*Vec4D(cosh(y),0.,0.,sinh(y)));
    m_ktsum     += pit->second.Momentum();
  }  
}

void Ladder_Generator_Seeded::FillPropagators() {
  Vec4D qt(0.,0.,0.,0.);
  LadderMap::iterator pit1=p_emissions->begin(), pit2=pit1; pit2++;
  for (size_t i=0;i<p_emissions->size();i++) {
    if (m_cols[i]==colour_type::singlet) {
      qt -= pit1->second.Momentum();
      p_props->push_back(T_Prop(colour_type::singlet,qt,m_qt2min));
      pit2++; pit1++;
    }
    else if (m_cols[i]==colour_type::octet) {
      double y1, y2, wt1, wt8, ract;
      double rprev = (i>0 && m_cols[i]==colour_type::singlet)?1.:0.;
      for (size_t j=0;j<1+m_emissions[i];j++) {
	qt  -= pit1->second.Momentum();
	y1   = pit1->first, y2 = pit2->first;
	wt1  = m_density.SingletWeight(y1,y2);
	wt8  = m_density.OctetWeight(y1,y2);
	ract = wt1/(wt1+wt8);
	if (ract<ran->Get()) p_props->push_back(T_Prop(colour_type::octet,qt,m_qt2min));
	else {
	  if (ract>rprev && !(i<2 && m_cols[i+1]==colour_type::singlet)) {
	    if (rprev!=1. && p_props->back().Col()==colour_type::singlet) {
	      p_props->back().SetCol(colour_type::octet);
	    }
	    p_props->push_back(T_Prop(colour_type::singlet,qt,m_qt2min));
	  }
	  else p_props->push_back(T_Prop(colour_type::octet,qt,m_qt2min));
	}
	pit2++; pit1++;
      }
    }
  }
}

ATOOLS::Vec4D Ladder_Generator_Seeded::
SelectKT(const double & y, const double & kt2min, const double & kt2max) {
  double kt2    = 0., rand, weight; 
  MakeTransverseUnitVector();
  if (y>=m_Ymax)       kt2 = p_eikonal->FF(0)->SelectQT2(kt2max,m_qt2minFF);
  else if (y<=-m_Ymax) kt2 = p_eikonal->FF(1)->SelectQT2(kt2max,m_qt2minFF);
  else  {
    do {
      rand   = ATOOLS::ran->Get();
      kt2    = kt2min * ( pow(kt2max/kt2min+1., rand) - 1.);
      weight = AlphaSWeight(kt2+m_kt2min);
      if (kt2>m_kt2max) weight *= (*p_alphaS)(kt2+m_kt2min)/(*p_alphaS)(m_kt2max+m_kt2min);
      //kt2  = 1./((1.-rand)/kt2max+rand/kt2min);
    } while (weight<ran->Get());
  }
  return sqrt(kt2) * m_eqt;
}

void Ladder_Generator_Seeded::CalculateWeight() {
  m_weight  = 1.; //TotalReggeWeight(p_ladder);
  // adding extra factors of alphaS for singlets
  LadderMap::iterator lit1=p_ladder->GetEmissions()->begin(), lit2=lit1; lit2++;
  TPropList::iterator pit=p_ladder->GetProps()->begin(), hardest;
  while (lit2!=p_ladder->GetEmissions()->end() && pit!=p_ladder->GetProps()->end()) {
    if (pit->Col()==colour_type::singlet)
      m_weight *= ( pit->Q02()/(pit->QT2()+pit->Q02()) *
		    (*p_alphaS)(lit1->second.Momentum().PPerp2()) *
		    (*p_alphaS)(lit2->second.Momentum().PPerp2()) );
    pit++; lit1++; lit2++;
  }
  double x[2], ratio=1., qt2=0.;
  if (p_ladder->ExtractHardest(hardest,qt2)) {
    if (qt2>m_kt2max) m_weight *= sqr((m_kt2max+m_kt2min)/(qt2+m_kt2min));
  }
  for (size_t beam=0;beam<2;beam++) {
    x[beam]         = 2.*p_ladder->InPart(beam)->Momentum()[0]/rpa->gen.Ecms();
    if (x[beam]>0.99 || x[beam]<m_partonic.MinX(beam)) continue;
    m_newxpdf[beam] = Max(0.,m_partonic.PDF(beam,x[beam],m_kt2min));
    ratio           = m_newxpdf[beam]/m_oldxpdf[beam];
    m_weight       *= ratio;
  }
}

