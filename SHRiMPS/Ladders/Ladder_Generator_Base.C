#include "SHRiMPS/Ladders/Ladder_Generator_Base.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <list>

using namespace SHRIMPS;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Ladder_Generator_Base::Ladder_Generator_Base() :
  m_Ymax(MBpars.GetEikonalParameters().Ymax),
  m_deltaY(MBpars.GetEikonalParameters().cutoffY),
  m_qt2min(MBpars.GetLadderParameters().Q02),
  m_kt2min(MBpars.GetLadderParameters().Q02),
  m_kt2minShower(MBpars.GetShowerLinkParameters().KT2min),
  p_alphaS(new Strong_Coupling(static_cast<Running_AlphaS *>
			       (s_model->GetScalarFunction(string("alpha_S"))),
			       asform::smooth,
			       MBpars.GetLadderParameters().Qas2)),
  m_density(MBpars.GetEikonalParameters().Delta,
	    MBpars.GetEikonalParameters().lambda,m_Ymax,
	    MBpars.GetEikonalParameters().absorp),
  m_mecorrection(ME_Correction(m_kt2min)),
 p_ladder(0)
{}

Ladder_Generator_Base::~Ladder_Generator_Base() {
  delete p_alphaS;
}


void Ladder_Generator_Base::InitLadder(const Vec4D & pos) {
  m_shat      = 4.*m_E[0]*m_E[1];
  m_weight    = 1.;
  p_ladder    = new Ladder(pos);
  p_emissions = p_ladder->GetEmissions();
  p_props     = p_ladder->GetProps();
}

void Ladder_Generator_Base::ConstructSimpleLadder() {
  size_t dir    = (dabs(p_emissions->begin()->first) >
		   dabs(p_emissions->rbegin()->first) ? 0 : 1);
  double qt2max = sqr(m_E[dir]/(cosh(dir==0 ?
				     p_ladder->GetEmissions()->begin()->first :
				     p_ladder->GetEmissions()->rbegin()->first)));
  do {
    m_qt2 = p_eikonal->FF(dir)->SelectQT2(qt2max);
  } while (ReggeWeight(m_qt2,m_ylimits[0],m_ylimits[1]) < ran->Get());
  MakeTransverseUnitVector();
  Vec4D k[2];
  k[0] = sqrt(m_qt2)*(Vec4D(cosh(m_ylimits[0]),0.,0.,sinh(m_ylimits[0])) + m_eqt);
  k[1] = sqrt(m_qt2)*(Vec4D(cosh(m_ylimits[1]),0.,0.,sinh(m_ylimits[1])) - m_eqt);
  p_emissions->begin()->second.SetMomentum(k[0]);
  p_emissions->rbegin()->second.SetMomentum(k[1]);
  T_Prop & prop = *p_props->begin();
  prop.SetQT2(m_qt2);
  prop.SetQ02(m_qt2min);
  prop.SetQ(sqrt(m_qt2)*m_eqt);
}

void Ladder_Generator_Base::ConstructISKinematics() {
  p_ladder->UpdatePropagatorKinematics();
  Vec4D Ksum = p_ladder->FSMomentum();
  for (size_t i=0;i<2;i++) {
    p_ladder->InPart(i)->SetMomentum(i==0 ?
				     Ksum.PPlus()/2.  * Vec4D(1.,0.,0., 1.) :
				     Ksum.PMinus()/2. * Vec4D(1.,0.,0.,-1.) );
    p_ladder->InPart(i)->SetBeam(i);
  }
  p_emissions->begin()->second.SetBeam(p_ladder->InPart(0)->Beam());
  p_emissions->rbegin()->second.SetBeam(p_ladder->InPart(1)->Beam());
}

void Ladder_Generator_Base::MakeTransverseUnitVector() {
  double phi = 2.*M_PI*ran->Get();
  m_eqt = Vec4D(0.,cos(phi),sin(phi),0.);
}

void Ladder_Generator_Base::ResetFSFlavours() {
  for (LadderMap::iterator lit=p_emissions->begin();lit!=p_emissions->end();lit++) {
    lit->second.SetFlavour(Flavour(kf_gluon));
  }
  for (TPropList::iterator pit=p_props->begin();pit!=p_props->end();pit++) {
    if (pit->Col()==colour_type::triplet) pit->SetCol(colour_type::octet);
  }
}

void Ladder_Generator_Base::QuarkReplace() {
  LadderMap::iterator lit1=p_emissions->begin(), lit2=lit1; lit2++;
  TPropList::iterator pit=p_props->begin();
  double pp = p_ladder->InPart(0)->Momentum().PPlus(), ppold=pp;
  bool last = false;
  do {
    double pp = pit->Q().PPlus();
    if (last) last=false;
    else {
      if (pit->Col()!=colour_type::singlet) {
	if (pp > ran->Get()*ppold) {
	  pit->SetCol(colour_type::triplet);
	  Flavour flav = Flavour(int(1+ran->Get()*3.)); 
	  lit1->second.SetFlavour(flav);
	  lit1->second.SetFlow(2,0);
	  lit2->second.SetFlavour(flav.Bar());
	  lit2->second.SetFlow(1,0);
	  last = true;
	}
      }
    }
    ppold = pp;
    lit1++; lit2++; pit++; 
  } while (pit!=p_props->end());
}

double Ladder_Generator_Base::AlphaSWeight(const double & kt2) {
  return AlphaS(kt2)/AlphaS(0.);
}

double Ladder_Generator_Base::ReggeWeight(const double & qt2, const double & y1,
					  const double y2) {
  return (qt2>m_qt2min ? 
	  exp(-3.*AlphaS(qt2)/(2.*M_PI) * dabs(y1-y2) * log(qt2/m_qt2min))  : 1.);
}

double Ladder_Generator_Base::LDCWeight(const double & qt2, const double & qt2prev) {
  return qt2/Max(qt2,qt2prev);
}

double Ladder_Generator_Base::TWeight() {
  if (p_ladder->GetProps()->size()==1) return 1.;
  double qt2max = m_qt2min, qt2, Q2;
  colour_type::code ctype=colour_type::octet;
  for (TPropList::iterator pit=p_ladder->GetProps()->begin();
       pit!=p_ladder->GetProps()->end();pit++) {
    qt2      = pit->QT2();
    if (qt2>qt2max) {
      qt2max = qt2;
      ctype  = pit->Col();
    }
  }
  return (ctype==colour_type::triplet ? 1:
	  ctype==colour_type::singlet ? sqr(m_qt2min/qt2max) : m_qt2min/qt2max);
}

		 
void Ladder_Generator_Base::Test() {
  vector<vector<Omega_ik *> > * eikonals(MBpars.GetEikonals());
  double b1, b2, y, asym12, asym34, d1, d2, d3, d4;
  for (size_t i=0;i<eikonals->size();i++) {
    for (size_t j=i;j<(*eikonals)[i].size();j++) {
      msg_Out()<<"=================================\n"
	       <<"Testing eikonals["<<i<<"]["<<j<<"]\n";
      if (i==j) {
	InitCollision((*eikonals)[i][j],0.);
	for (int k=0;k<3;k++) {
	  for (int l=k;l<3;l++) {
	    b1 = double(k)*2.;
	    b2 = double(l)*2.;
	    msg_Out()<<"   for b1 = "<<b1<<", b2 = "<<b2<<"\n";
	    for (int m=0;m<8;m++) {
	      y    = double(m);
	      m_density.SetImpactParameters(b1,b2);
	      d1     = m_density(y);
	      d2     = m_density(-y);
	      asym12 = (d1-d2)/(d1+d2);
	      if (l!=k) {
		m_density.SetImpactParameters(b2,b1);
		d3     = m_density(y);
		d4     = m_density(-y);
		asym34 = (d3-d4)/(d3+d4);
		msg_Out()<<"  y = "<<y<<", asym = "<<(asym12+asym34)<<" ["<<asym12<<" and "<<asym34<<"] "
			 <<"from d's = "<<d1<<", "<<d2<<", "<<d3<<", and "<<d4<<"\n";
	      }
	      else {
		msg_Out()<<"  y = "<<y<<", asym = "<<asym12<<" from d's = "<<d1<<", and "<<d2<<"\n";
	      }
	    }
	  }
	}
      }
      else {
	for (int m=0;m<8;m++) {
	  y      = double(m);
	  for (int k=0;k<2;k++) {
	    for (int l=0;l<2;l++) {
	      b1 = double(k)*2.;
	      b2 = double(l)*2.;
	      m_density.SetImpactParameters(b1,b2);
	      InitCollision((*eikonals)[i][j],0.);
	      d1     = m_density(y);
	      InitCollision((*eikonals)[j][i],0.);
	      d2     = m_density(-y);
	      asym12 = (d1-d2)/(d1+d2);
	      if (k!=l) {
		m_density.SetImpactParameters(b2,b1);
		InitCollision((*eikonals)[i][j],0.);
		d3     = m_density(y);
		InitCollision((*eikonals)[j][i],0.);
		d4     = m_density(-y);
		asym34 = (d3-d4)/(d3+d4);
	      }
	      else {
		msg_Out()<<"  y = "<<y<<", asym = "<<asym12<<" from d's = "<<d1<<", and "<<d2<<"\n";
	      }
	    }
	  }
	}
      }
    }
  }
  exit(1);
}
