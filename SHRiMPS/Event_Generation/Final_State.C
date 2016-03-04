#include "SHRiMPS/Event_Generation/Final_State.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/Random.H"

using namespace SHRIMPS;
using namespace ATOOLS;
using namespace MODEL;
using namespace std;

Final_State::Final_State():
  p_alphaS(new Strong_Coupling(static_cast<Running_AlphaS *>
			       (s_model->GetScalarFunction(string("alpha_S"))),
			       asform::smooth,0.25)),
  m_weights(Eikonal_Weights()),
  m_origY(MBpars.GetEikonalParameters().originalY),
  m_Ymax(MBpars.GetEikonalParameters().Ymax),
  m_DeltaY(dabs(m_origY-m_Ymax))
{}

Final_State::~Final_State() {
  if (p_alphaS) { delete p_alphaS; p_alphaS = NULL; }
}

void Final_State::SetEikonal(Omega_ik * eikonal) {
  m_weights.SetEikonal(eikonal);
  p_omegaik = eikonal->GetSingleTerm(0);
  p_omegaki = eikonal->GetSingleTerm(1);
}

void Final_State::SetLadder(Ladder * ladder) {
  p_ladder = ladder;
}

void Final_State::
SetImpactParameters(const double & b1, const double & b2) {
  m_b1 = b1; m_b2 = b2;
  m_weights.SetImpactParameters(m_b1,m_b2);
}

void Final_State::SetAvailableEnergy(const double & E) {
  m_availableE = E;
}

void Final_State::FillPrimaryLadder() {
  m_E = m_pz = 0.;
  AddInitialRapiditiesToLadder();
  m_weights.AddRapidities(p_ladder,-m_Ymax,m_Ymax);
  AddPropagators();
  SanitizePropagators();
  ConstructKinematics();
}

void Final_State::AddInitialRapiditiesToLadder() {
  p_ladder->AddRapidity(-m_origY+ran->Get()*m_DeltaY);
  p_ladder->AddRapidity(m_origY-ran->Get()*m_DeltaY);
}

void Final_State::AddPropagators() {
  LadderMap::iterator lit1(p_ladder->GetEmissions()->begin()),lit2(lit1);
  lit2++;
  while(lit2!=p_ladder->GetEmissions()->end()) {
    T_Prop * prop(new T_Prop(m_weights.PropColour(lit1->first,lit2->first)));
    p_ladder->GetProps()->push_back(*prop);
    if (prop->m_col==colour_type::singlet) p_ladder->SetDiffractive(true);
    lit1 = lit2;
    lit2++;
  }
}

void Final_State::SanitizePropagators() {
  LadderMap::iterator lit1(p_ladder->GetEmissions()->begin()), lit2(lit1);
  lit2++;
  TPropList::iterator prop(p_ladder->GetProps()->begin());
  while (lit2!=p_ladder->GetEmissions()->end() &&
	 prop!=p_ladder->GetProps()->end()) {
    if (prop->m_col==colour_type::singlet) {
      CheckNextPropagator(prop,lit2);
    }
    else prop++;
    lit1++; lit2++;
  }
}

void Final_State::
CheckNextPropagator(TPropList::iterator & prop,LadderMap::iterator & lit2) {
  prop++;
  if (prop==p_ladder->GetProps()->end() ||
      prop->m_col!=colour_type::singlet) return;
  LadderMap::iterator & lit1(lit2); lit1--;
  double wt12(m_weights.WeightSingletOverOctet(lit1->first,lit2->first));
  LadderMap::iterator & lit3(lit2); lit3++;
  double wt23(m_weights.WeightSingletOverOctet(lit2->first,lit3->first));
  msg_Out()<<"Found double singlet in around "<<lit2->first<<": ";
  if (wt12>wt23) {
    prop->m_col = colour_type::octet;
    msg_Out()<<"made second interval octet.\n";
  }
  else {
    TPropList::iterator & prev(prop); prev--;
    prev->m_col = colour_type::octet;
    msg_Out()<<"made first interval octet.\n";
  }
}

void Final_State::ConstructKinematics() {
  LadderMap::iterator forward(p_ladder->GetEmissions()->begin());
  LadderMap::iterator backward(p_ladder->GetEmissions()->end());backward--;
  LadderMap::iterator next;
  TPropList::iterator down(p_ladder->GetProps()->begin());
  TPropList::iterator up(p_ladder->GetProps()->end());up--;
  Vec4D upqt(0.,0.,0.,0.),downqt(0.,0.,0.,0.);
  double qt,y,E,pz;
  while (forward!=backward) {
    if (dabs(forward->first)>dabs(backward->first)) {
      next   = forward; next++;
      y      = forward->first;
      ConstructPropagator(y,next->first,&(*down),1);
      downqt = -downqt+down->m_q;
      qt     = downqt.PPerp();
      m_E   += E  = qt*cosh(y);
      m_pz  += pz = qt*sinh(y);
      forward->second.m_mom = Vec4D(E,0.,0.,pz)-downqt;
      forward->second.m_flav = Flavour(kf_gluon);
      forward++; down++;
    }
    else {
      next   = backward; next--;
      y      = backward->first;
      ConstructPropagator(y,next->first,&(*up),-1);
      upqt   = -upqt+up->m_q;
      qt     = upqt.PPerp();
      m_E   += E  = qt*cosh(y);
      m_pz  += pz = qt*sinh(y);
      backward->second.m_mom = Vec4D(E,0.,0.,pz)+upqt;
      backward->second.m_flav = Flavour(kf_gluon);
      backward--; up--;
    }
  }
  Vec4D qtvec(down->m_q-up->m_q);
  qt    = qtvec.PPerp();
  y     = next->first;
  m_E  += E  = qt*cosh(y);
  m_pz += pz = qt*sinh(y);
  next->second.m_mom = qt*Vec4D(E,0.,0.,pz)-qtvec;
  next->second.m_flav = Flavour(kf_gluon);
}

void Final_State::ConstructPropagator(const double & y1,const double & y2,
					   T_Prop* prop,const int & sign) {
  double qt2max(sqr(m_availableE/(2.*cosh(y1))));
  double qt2(QT2(qt2max,y1,y2,prop->m_col)), qt(sqrt(qt2));
  double phi(ran->Get()*2.*M_PI);
  Vec4D qtvec(sign*qt*Vec4D(0.,cos(phi),sin(phi),0.));
  prop->m_qt2=qt2;
  prop->m_q=qtvec;
}

double Final_State::QT2(const double & qt2max,
			const double & y1,const double & y2,
			const colour_type::code & col) {
  double m_Ymax(MBpars.GetEikonalParameters().Ymax);
  double qt2(-1.);
  if (dabs(y1)>m_Ymax) {
    if (y1<y2) qt2 = p_omegaik->FF1()->SelectQT2(1.);
    else qt2 = p_omegaki->FF1()->SelectQT2(1.);
  }
  else {
    double reggefac(0.);
    if (col==colour_type::octet)
      reggefac = 3./M_PI*dabs(y2-y1);
    qt2 = SelectQT2(qt2max,-(1+reggefac*p_alphaS->MaxValue()));
  } 
  return qt2;
}

double Final_State::SelectQT2(const double & qt2max,const double & expo) {
  double mu2(1.),rand(ran->Get());
  return pow(rand*pow(qt2max+mu2,expo)+
	     (1.-rand)*pow(mu2,expo),1./expo)-mu2;
}


void Final_State::Test(const std::string & dirname) {
  m_weights.Test(dirname);
}
