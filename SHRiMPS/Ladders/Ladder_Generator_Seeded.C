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
  for (size_t i=0;i<2;i++) m_pbeam[i] = rpa->gen.Ecms()/2. * Vec4D(1.,0.,0.,i==0?1:-1);
}

Ladder_Generator_Seeded::~Ladder_Generator_Seeded() {}

Ladder * Ladder_Generator_Seeded::operator()(const Vec4D & pos) {
  SeedLadder(pos);
  AddBeamPartons();
  FillIntervals();
  CompensateKTs();
  ConstructFSMomenta();
  FillPropagators();
  ConstructISKinematics();
  //msg_Out()<<METHOD<<"["<<m_emissions[0]<<", "<<m_emissions[1]<<", "<<m_emissions[2]<<"], "
  //	   <<"ktsum = "<<m_ktsum<<" for incoming E's "<<p_ladder->InPart(0)->Momentum()[0]<<" and "
  //	   <<p_ladder->InPart(1)->Momentum()[0]<<" \n";//<<(*p_ladder)<<"\n";
  return p_ladder;
}

void Ladder_Generator_Seeded::SeedLadder(const Vec4D & pos) {
  p_ladder     = new Ladder(pos);
  p_emissions  = p_ladder->GetEmissions();
  p_props      = p_ladder->GetProps();
  m_shat       = m_partonic.MakeEvent();
  m_kt2max     = m_partonic.PT2();
  m_phi        = m_partonic.Phi();
  m_ktsum      = Vec4D(0.,0.,0.,0.);
  double kt    = sqrt(m_kt2max);
  Vec4D  ktvec = kt*Vec4D(0.,cos(m_phi),sin(m_phi),0.);
  for (size_t i=0;i<2;i++)
    p_ladder->AddRapidity(m_yseed[i+1] = m_partonic.Y(i),m_partonic.Flav(i),
			  (i==0?1.:-1.)*ktvec);
}

void Ladder_Generator_Seeded::AddBeamPartons() {
  Vec4D ktvec;
  for (size_t beam=0;beam<2;beam++) {
    if ((beam==0 && m_yseed[1]< m_Ymax) || (beam==1 && m_yseed[2]>-m_Ymax)) {
      m_yseed[beam==0?0:3] = m_ylimits[beam] = (beam==0?1.:-1.) * (m_Ymax + ran->Get()*m_deltaY);
      m_ktsum += ktvec = SelectKT(m_ylimits[beam]);
      //msg_Out()<<METHOD<<" adds beam["<<beam<<"] at y = "
      //       <<m_ylimits[beam]<<" = "<<m_yseed[beam==0?0:3]<<"\n";
      p_ladder->AddRapidity(m_ylimits[beam],Flavour(kf_gluon),ktvec);
    }
    else {
      //msg_Out()<<METHOD<<" identifies beam["<<beam<<"] at y = "
      //       <<m_ylimits[beam]<<" = "<<m_yseed[beam==0?0:3]<<"\n";
      m_ylimits[beam] = m_yseed[beam==0?0:3] = beam==0?m_yseed[1]:m_yseed[2];
    }
  }
}

void Ladder_Generator_Seeded::FillIntervals() {
  double wt1, wt8;
  for (size_t i=0;i<3;i++) {
    if (m_yseed[i]>m_yseed[i+1]+1.e-6) {
      wt1       = m_density.SingletWeight(m_yseed[i],m_yseed[i+1]);
      wt8       = m_density.OctetWeight(m_yseed[i],m_yseed[i+1]);
      if (wt1/(wt1+wt8)>ran->Get()) {
	m_cols[i]      = colour_type::singlet;
	m_emissions[i] = 0;
      }
      else {
	m_cols[i]      = colour_type::octet;
	m_emissions[i] = m_density.NGluons(m_yseed[i], m_yseed[i+1]);
	Vec4D  ktvec;
	if (m_emissions[i]>0) {
	  //msg_Out()<<METHOD<<" needs to add "<<m_emissions[i]<<" emissions in "<<i<<"th interval:\n";
	  for (size_t j=0;j<m_emissions[i];j++) {
	    double y = m_density.SelectRapidity(m_yseed[i], m_yseed[i+1]);
	    m_ktsum += ktvec = SelectKT(y);
	    p_ladder->AddRapidity(y,Flavour(kf_gluon),ktvec);
	    //msg_Out()<<"   * y = "<<y<<"\n";
	  }
	}
      }
      //msg_Out()<<METHOD<<"("<<m_yseed[i]<<", "<<m_yseed[i+1]<<"): "
      //       <<"wt1 = "<<wt1<<", wt8 = "<<wt8<<" --> "<<m_cols[i]<<", N = "<<m_emissions[i]<<"\n";
    }
    else {
      m_cols[i]      = colour_type::none;
      m_emissions[i] = 0;
      //msg_Out()<<METHOD<<"("<<m_yseed[i]<<", "<<m_yseed[i+1]<<") --> "<<m_cols[i]<<"\n";
    }
  }
}

void Ladder_Generator_Seeded::CompensateKTs() {
  //msg_Out()<<METHOD<<" must compensate for "<<m_ktsum<<", distribute over hard process.\n";
  for (size_t i=1;i<3;i++) {
    Vec4D before = (*p_emissions)[m_yseed[i]].Momentum();
    (*p_emissions)[m_yseed[i]].SetMomentum(before - m_ktsum/2.);
    Vec4D after  = (*p_emissions)[m_yseed[i]].Momentum();
    //msg_Out()<<METHOD<<"("<<i<<"): "<<before<<" --> "<<after<<"\n";
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
  //msg_Out()<<METHOD<<" yields total mom = "<<m_ktsum<<"\n";
}

void Ladder_Generator_Seeded::FillPropagators() {
  Vec4D qt(0.,0.,0.,0.);
  LadderMap::iterator pit1=p_emissions->begin(), pit2=pit1; pit2++;
  for (size_t i=0;i<3;i++) {
    if (m_cols[i]==colour_type::singlet) {
      qt -= pit1->second.Momentum();
      //msg_Out()<<METHOD<<"["<<i<<"] adds singlet propagator: qt = "<<qt<<"\n";
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
	//msg_Out()<<METHOD<<"["<<i<<j<<"("<<m_emissions[i]<<")]: "
	//	 <<"wt1 = "<<wt1<<", wt8 = "<<wt8<<", ract = "<<ract<<"("<<rprev<<")\n";
	if (ract<ran->Get()) {
	  //msg_Out()<<"   adds octet propagator: qt = "<<qt<<"\n";
	  p_props->push_back(T_Prop(colour_type::octet,qt,m_qt2min));
	}
	else {
	  if (ract>rprev && !(i<2 && m_cols[i+1]==colour_type::singlet)) {
	    if (rprev!=1. && p_props->back().Col()==colour_type::singlet) {
	      //msg_Out()<<"   overwrites previous propagator colour.\n";
	      p_props->back().SetCol(colour_type::octet);
	    }
	    //msg_Out()<<"   adds singlet propagator: qt = "<<qt<<".\n";
	    p_props->push_back(T_Prop(colour_type::singlet,qt,m_qt2min));
	  }
	  else {
	    //msg_Out()<<"   overwrites weight and adds octet propagator: qt = "<<qt<<".\n";
	    p_props->push_back(T_Prop(colour_type::octet,qt,m_qt2min));
	  }
	}
	pit2++; pit1++;
      }
    }
  }
}

ATOOLS::Vec4D Ladder_Generator_Seeded::SelectKT(const double & y) {
  double kt2 = 0., kt2max = Min(m_kt2max, m_E[0]*m_E[1]/sqr(2.*cosh(y)));
  //msg_Out()<<METHOD<<"(y = "<<y<<", cosh(y) = "<<cosh(y)<<", "
  //	   <<"kt2 < "<<kt2max<<" from "<<m_kt2max<<" and "
  //	   <<m_E[0]*m_E[1]/(4.*cosh(y))<<")\n";
  MakeTransverseUnitVector();
  if (y>=m_Ymax)       kt2 = p_eikonal->FF(0)->SelectQT2(kt2max,m_qt2minFF);
  else if (y<=-m_Ymax) kt2 = p_eikonal->FF(1)->SelectQT2(kt2max,m_qt2minFF);
  else  {
    do {
      kt2 = m_kt2min * ( pow(kt2max/m_kt2min+1., ATOOLS::ran->Get()) - 1.);
    } while (AlphaSWeight(kt2+m_kt2min)<ran->Get());
  }
  return sqrt(kt2) * m_eqt;
}

void Ladder_Generator_Seeded::CalculateWeight() {
  m_weight  = 1.;
  // we may want to add some Regge/Angular Odering weights etc. here.  
}

