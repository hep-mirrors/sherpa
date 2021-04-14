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

Ladder_Generator_KT::Ladder_Generator_KT() :
  Ladder_Generator_Base(),
  m_S(sqr(rpa->gen.Ecms())), m_shatmin(4.*m_kt2min)
{
  for (size_t i=0;i<2;i++) m_Ebeam[i] = rpa->gen.PBeam(i)[0];
}

Ladder * Ladder_Generator_KT::operator()(const Vec4D & pos) {
  InitLadder(pos);
  if (FixInitialState() &&
      FillForwardPartons() &&
      FillEmissions()) { 
    SelectPropagatorColours();
    CalculateWeight();
  }
  else { delete p_ladder; p_ladder = NULL; }
  return p_ladder; 
}

bool Ladder_Generator_KT::FixInitialState() {
  double xmax[2];
  for (size_t i=0;i<2;i++) {
    xmax[i] = m_E[i]/m_Ebeam[i];
    if (xmax[i] < Max(2./m_Ebeam[i], 1.e-3)) return false;
  }
  m_shat          = (p_pdf[0]->XMaxPDF(kf_gluon)*p_pdf[1]->XMaxPDF(kf_gluon)) * m_S;
  double eta      = m_density.EffectiveIntercept(m_b1,m_b2,0.);
  double maxvalue = ( (p_pdf[0]->XPDFMax(kf_gluon)*p_pdf[1]->XPDFMax(kf_gluon)) *
		      pow(m_shat/m_shatmin,eta) / (2.*m_shat)  );
  double x[2], value;
  do {
    value  = 1.;
    m_shat = m_S;
    for (size_t i=0;i<2;i++) {
      x[i]    = ran->Get() * 0.99*xmax[i];
      if (x[i]<0.) {
	msg_Error()<<METHOD<<" is funny: "<<x[i]<<" from "<<xmax[i]<<" and "<<m_E[i]<<"\n";
	return false;
      }
      p_pdf[i]->Calculate(x[i],0.);
      m_shat *= x[i];   
      value  *= p_pdf[i]->XPDF(Flavour(kf_gluon));
    }
    value    *= pow(m_shat,eta)/(2.*m_shat);
  } while(m_shat<m_shatmin && value/maxvalue<ran->Get());
  for (size_t i=0;i<2;i++) {
    m_q[i] = x[i]*rpa->gen.PBeam(i);
    p_ladder->InPart(i)->SetMomentum(m_q[i]);
    p_ladder->InPart(i)->SetBeam(m_dir==1?1-i:i);
    if (m_q[i][0]>m_E[i]) {
      msg_Error()<<METHOD<<" gives funny momentum for x = "<<x[i]<<": "<<m_q[i]<<"\n";
      return false;
    }
  }
  return (m_shat>m_shatmin); 
}

bool Ladder_Generator_KT::FillForwardPartons() {
  for (size_t i=0;i<2;i++) {
    m_ylimits[i]  = (i==0?1.:-1.) * (m_Ymax + ran->Get()*m_deltaY);
    double qt2max = m_shat/(4.*sqr(cosh(m_ylimits[i])));
    double qt2    = p_eikonal->FF(i)->SelectQT2(qt2max);
    double phi    = 2.*M_PI*ran->Get();
    Vec4D  mom    = sqrt(qt2) * Vec4D(cosh(m_ylimits[i]),cos(phi),sin(phi),sinh(m_ylimits[i]));
    p_ladder->AddRapidity(m_ylimits[i],p_ladder->InPart(i)->Flavour(),mom);
    m_q[i]       -= mom;
    if (m_q[i][0]<0.) return false;
  }
  return true;
}

bool Ladder_Generator_KT::FillEmissions() {
  return (MakeTrialLadder() && p_ladder->CheckFourMomentum());
}

bool Ladder_Generator_KT::MakeTrialLadder() {
  size_t dlast;
  double ylast, wt1, wt8, ratio1;
  for (size_t i=0;i<2;i++) m_y[i][0] = m_y[i][1] = m_ylimits[i];
  m_pit[0] = p_ladder->GetProps()->begin(); m_pit[0]++;
  m_pit[1] = p_ladder->GetProps()->end();
  do {
    m_seff     = (m_q[0]+m_q[1]).Abs2();
    if (m_seff<m_kt2min) break;
    double arg = M_PI/(3.*AlphaS(0.)*log(m_seff/m_kt2min));
    size_t dir = (dabs(m_y[0][0]) > dabs(m_y[1][0])) ? 0 : 1;
    if (TrialEmission(m_y[dir][0],arg,dir)                       &&
	AbsorptionWeight(m_k,m_y[dir][0]) > ran->Get()           &&
	ReggeWeight(m_q[dir].PPerp2(),m_y[dir][0],m_y[dir][1]) > ran->Get() &&
	Ordering(dir)) {
      AddEmission(dir,m_y[dir][0]);
      m_y[dir][1] = ylast = m_y[dir][0];
      dlast       = dir;
    }
  } while (m_y[0][0]>m_y[1][0]);
  if (p_ladder->GetEmissions()->size()==2) return FixSimpleKinematics();
  return AdjustLastEmission(dlast,ylast,m_y[1-dlast][1]);
}

bool Ladder_Generator_KT::
TrialEmission(double & y,const double & arg,const size_t dir) {
  m_k = Vec4D(0.,0.,0.,0.);
  size_t attempts = 0;
  do {
    y    -= ((dir) ? arg : -arg) * log(ran->Get());
    m_kt2 = m_kt2min * pow(m_seff/m_kt2min,ran->Get());
    if (attempts > 1000) return false;
    attempts++;
  } while (m_kt2>m_seff/(4.*sqr(cosh(y))) &&
	   AlphaS(m_kt2)/AlphaS(0)<ran->Get());
    
  double phi = 2.*M_PI*ran->Get();
  m_k        = sqrt(m_kt2)*Vec4D(cosh(y),cos(phi),sin(phi),sinh(y));
  return true;
}

double Ladder_Generator_KT::AbsorptionWeight(const Vec4D & k,const double & y) {
  return m_density.AbsorptionWeight(y);
}

bool Ladder_Generator_KT::Ordering(const size_t dir) {
  Vec4D  qtest = m_q[dir]-m_k;
  return (qtest[0]>0 && qtest.PPerp2()>m_q[dir].PPerp2());
}

void Ladder_Generator_KT::AddEmission(const size_t dir,const double & y) {
  p_ladder->AddRapidity(y,Flavour(kf_gluon),m_k);
  p_ladder->GetProps()->insert(m_pit[dir],T_Prop(colour_type::octet,m_q[dir],m_qt2min));
  if (dir==1) { for (size_t i=0;i<2;i++) m_pit[i]--; }
  m_q[dir] -= m_k;
}

bool Ladder_Generator_KT::FinishKinematics() {
  return p_ladder->CheckFourMomentum();
}

void Ladder_Generator_KT::CalculateWeight() {
  m_weight *= KTMaxWeight();
}

double Ladder_Generator_KT::KTMaxWeight() {
  if (!p_ladder) return 0.;
  double kt2max = m_kt2min, kt2test;
  for (LadderMap::iterator lit=p_ladder->GetEmissions()->begin();
       lit!=p_ladder->GetEmissions()->end();lit++) {
    kt2test = lit->second.Momentum().PPerp2();
    if (kt2test>kt2max) kt2max = kt2test;
  }
  return m_kt2min/kt2max;
}

void Ladder_Generator_KT::SelectPropagatorColours() {
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

bool Ladder_Generator_KT::FixSimpleKinematics() {
  p_ladder->GetEmissions()->clear();
  Vec4D  in     = p_ladder->InPart(0)->Momentum()+p_ladder->InPart(1)->Momentum();
  m_seff        = in.Abs2();
  if (m_seff<m_kt2min) return false;
  double E2     = sqr(in[0]), P2 = sqr(in[3]), s = E2-P2, y[2];
  double qt2    = 0.;
  size_t trials = 0;
  do {
    for (size_t i=0;i<2;i++) qt2 += p_eikonal->FF(i)->SelectQT2(m_shatmin);
    trials++;
    if (trials>1000) return false;
  } while (4.*qt2>s);
  double phi = 2.*M_PI*ran->Get();
  Vec4D  qt  = sqrt(qt2/2.) * Vec4D(0.,cos(phi),sin(phi),0.);
  for (size_t i=0;i<2;i++) {
    double k_long = 1./2. * ( in[3] + (i==0 ?1.:-1.) * in[0] * sqrt(1.-4.*qt2/s) );
    double omega  = sqrt(k_long*k_long + qt2);
    Vec4D k       = Vec4D(omega,0.,0.,k_long) + (i==0 ?1:-1.) * qt;
    y[i]          = k.Y();
    p_ladder->AddRapidity(y[i],p_ladder->InPart(0)->Flavour(),k);
  }
  p_ladder->AddPropagator(T_Prop(colour_type::octet,qt,m_qt2min));
  m_weight = ReggeWeight(qt.PPerp2(),y[0],y[1]);
  return true;
}

bool Ladder_Generator_KT::
AdjustLastEmission(const bool dlast, const double & ylast,const double & y0) {
  if (dlast==0) m_pit[dlast]--;
  Vec4D oldk = (*p_ladder->GetEmissions())[ylast].Momentum();
  Vec4D QT   = (m_q[1-dlast]+m_pit[dlast]->Q()).Perp();
  Vec4D newk = QT.PPerp()*Vec4D(cosh(ylast),0.,0.,sinh(ylast))+QT;  
  p_ladder->GetProps()->insert(m_pit[1-dlast],T_Prop(colour_type::octet,m_q[1-dlast],m_qt2min));
  (*p_ladder->GetEmissions())[ylast].SetMomentum(newk);
  m_weight = ReggeWeight(m_q[1-dlast].PPerp2(),ylast,y0);
  return UpdateInitialState();
}

bool Ladder_Generator_KT::UpdateInitialState() {
  ConstructISKinematics();
  return (p_ladder->InPart(0)->Momentum()[0] >
	  p_ladder->GetEmissions()->begin()->second.Momentum()[0] &&
	  p_ladder->InPart(0)->Momentum()[0] < m_E[0] &&
	  p_ladder->InPart(1)->Momentum()[0] >
	  p_ladder->GetEmissions()->rbegin()->second.Momentum()[0] &&
	  p_ladder->InPart(1)->Momentum()[0] < m_E[1]);
}
