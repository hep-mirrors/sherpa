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

void Final_State::SetNLadders(const size_t & N) { m_Nladders = N; }

void Final_State::SetLadder(Ladder * ladder) {
  p_ladder = ladder;
}

void Final_State::
SetImpactParameters(const double & b1, const double & b2) {
  m_b1 = b1; m_b2 = b2;
  m_weights.SetImpactParameters(m_b1,m_b2);
}

void Final_State::SetAvailableEnergies(const double E[2]) {
  m_E[0] = E[0]; m_E[1] = E[1];
}

void Final_State::FillPrimaryLadder() {
  m_pp = m_pm = 0.;
  m_totkt = Vec4D(0.,0.,0.,0.);
  m_weights.AddRapidities(p_ladder,-m_Ymax,m_Ymax);
  AddInitialRapiditiesToLadder();
  AddPropagators();
  SanitizePropagators();
  ConstructKinematics();
}

void Final_State::AddInitialRapiditiesToLadder() {
  double miny(0.), maxy(0.);
  if (p_ladder->Size()>0) {
    miny = Min(miny,p_ladder->GetEmissions()->begin()->first);
    maxy = Max(maxy,p_ladder->GetEmissions()->rbegin()->first);
  }
  double y, m_Delta(0.2);
  do {
    y = miny+ran->Get()*(-m_origY-miny);
  } while (dabs(y)>m_Ymax && exp(-m_Delta*dabs(y-miny))>ran->Get());
  p_ladder->AddRapidity(y);
  do {
    y = maxy+ran->Get()*(m_origY-maxy);
  } while (dabs(y)>m_Ymax && exp(-m_Delta*dabs(y-maxy))>ran->Get());
  p_ladder->AddRapidity(y);
  //msg_Out()<<"-------------------------------------------------------\n"
  //	   <<(*p_ladder->GetEmissions())
  //	   <<"-------------------------------------------------------\n";
  //p_ladder->AddRapidity(-m_origY+ran->Get()*m_DeltaY);
  //p_ladder->AddRapidity(m_origY-ran->Get()*m_DeltaY);
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
    prop++; lit1++; lit2++;
  }
}

void Final_State::
CheckNextPropagator(TPropList::iterator & prop,LadderMap::iterator & lit2) {
  TPropList::iterator prop1(prop); prop1++;
  if (prop1==p_ladder->GetProps()->end() ||
      prop1->m_col!=colour_type::singlet) return;
  LadderMap::iterator & lit1(lit2); lit1--;
  double wt12(m_weights.WeightSingletOverOctet(lit1->first,lit2->first));
  LadderMap::iterator & lit3(lit2); lit3++;
  double wt23(m_weights.WeightSingletOverOctet(lit2->first,lit3->first));
  if (wt12>wt23) {
    prop1->m_col = colour_type::octet;
    //msg_Out()<<"made second interval octet.\n";
  }
  else {
    prop->m_col = colour_type::octet;
    //msg_Out()<<"made first interval octet.\n";
  }
}

void Final_State::ConstructKinematics() {
  LadderMap::iterator forward(p_ladder->GetEmissions()->begin());
  LadderMap::iterator backward(p_ladder->GetEmissions()->end());backward--;
  LadderMap::iterator next;
  TPropList::iterator down(p_ladder->GetProps()->begin());
  TPropList::iterator up(p_ladder->GetProps()->end());up--;
  Vec4D downqt(0.,0.,0.,0.),upqt(0.,0.,0.,0.),qt;
  bool FWBW(false);
  while (forward!=backward) {
    if (dabs(forward->first)>dabs(backward->first)) {
      next   = forward; next++;
      ConstructPropagator(forward->first,next->first,down);
      ConstructEmission(forward,downqt-down->m_q);
      downqt = down->m_q;
      forward++; down++;
      FWBW = true;
    }
    else {
      next   = backward; next--;
      ConstructPropagator(backward->first,next->first,up);
      ConstructEmission(backward,-up->m_q-upqt);
      upqt   = -up->m_q;
      backward--; up--;
      FWBW = false;
    }
  }
  ConstructEmission(next,downqt-upqt);
  //msg_Out()<<METHOD<<"("<<FWBW<<"): "<<(downqt-upqt)<<".\n";
  if (dabs(m_totkt.Abs2())>1.e-6) {
    msg_Out()<<METHOD<<" violates transverse momentum:"<<m_totkt<<"\n";
    msg_Out()<<(*p_ladder)<<"\n";
    exit(1);
  }
}

void Final_State::ConstructPropagator(const double & y1,const double & y2,
				      TPropList::iterator & prop) {
  double qt2max(sqr((y1<0?m_E[0]:m_E[1])/cosh(y1))/m_Nladders);
  double qt2(QT2(qt2max,y1,y2,prop)), qt(sqrt(qt2));
  double phi(ran->Get()*2.*M_PI);
  Vec4D  qtvec(qt*Vec4D(0.,cos(phi),sin(phi),0.));
  prop->m_qt2 = prop->m_q2 = qt2;
  prop->m_q   = qtvec;
  prop->m_q02 = 1.;
  // we need m_q02 and m_q2 when we extract the kinematics of the hardest
  // propagator to reweight the ladder with an extra factor of 1/t or so.
  // This happens in the Ladder_Generator.
}

void Final_State::ConstructEmission(const LadderMap::iterator & emit,
				    const Vec4D & ktvec) {
  double kt(ktvec.PPerp()), y(emit->first), E(kt*cosh(y)), pz(kt*sinh(y));
  m_totkt += ktvec;
  m_pp += E+pz;
  m_pm += E-pz;
  emit->second.m_mom = Vec4D(E,0.,0.,pz)-ktvec;
  emit->second.m_flav = Flavour(kf_gluon);
}


double Final_State::
QT2(const double & qt2max,const double & y1,const double & y2,
    TPropList::iterator & prop) {
  double qt2(-1.);
  /*
    if (dabs(y1-p_ladder->GetEmissions()->begin()->first)<1.e-12) {
    qt2 = p_omegaik->FF1()->SelectQT2(Min(qt2max,sqr(1./cosh(y1))));
    }
    else if (dabs(y1-p_ladder->GetEmissions()->begin()->first)<1.e-12) {
    qt2 = p_omegaki->FF1()->SelectQT2(Min(qt2max,sqr(1./cosh(y1))));
    }
    else {
  */
  double reggefac(0.);
  if (prop->m_col==colour_type::octet) reggefac = 3./(4.*M_PI)*dabs(y2-y1);
  qt2 = SelectQT2(qt2max,reggefac*p_alphaS->MaxValue()); 
  return qt2;
}

double Final_State::SelectQT2(const double & qt2max,const double & expo) {
  // assume a probablity density for qt^2 of the form
  //  1/(qt^2+mu^2) * exp(-expo*alphaS(qt2)),
  // where expo is the typical regge exponent - it may need fixing -
  // which is set in QT2 and is either 0 for singlets or Nc/pi*|Delta y|
  double qt2(0.),mu2(1.);
  double arg((qt2max+mu2)/mu2),expmax(1.);
  do {
    qt2 = mu2*(pow(arg,ran->Get())-1.);
  } while (dabs(expo)>1.e-12 && qt2>mu2 && 
	   exp(-expo*(*p_alphaS)(qt2)*log(qt2/mu2))/expmax<ran->Get());
  return qt2;
}


void Final_State::Test(const std::string & dirname) {
  m_weights.Test(dirname);
}
