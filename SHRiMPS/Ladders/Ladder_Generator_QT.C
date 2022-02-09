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
  Ladder_Generator_Base(), m_seff(0.)
{ }

Ladder * Ladder_Generator_QT::operator()(const Vec4D & pos) {
  InitLadder(pos);
  size_t trials = 0;
  do {  
    if ((trials++)>1000) { delete p_ladder; return NULL; }
    if (FixInitialPartons() && MakeTrialLadder()) {
      ConstructISKinematics();
      SelectPropagatorColours();
      CalculateWeight();
      //msg_Out()<<METHOD<<"(weight = "<<m_weight<<"):\n"<<(*p_ladder)<<"\n";
    }
    else m_weight = 0.;
  } while (m_weight<ran->Get());    
  AddRescatters();
  return p_ladder;
}

bool Ladder_Generator_QT::FixInitialPartons() {
  p_emissions->clear(); p_props->clear();
  m_shat     = m_partonic.MakeEvent();
  if (m_shat<0.) return false;
  m_sigmahat = m_partonic.SigmaHat();
  for (size_t beam=0;beam<2;beam++) {
    m_ylimits[beam] = (beam==0 ? 1.: -1.) * (m_Ymax + ran->Get()*m_deltaY);
    m_qini[beam]    = m_q[beam] = 2.*m_partonic.X(beam) * rpa->gen.PBeam(beam);
    m_flavs[beam]   = m_partonic.Flav(beam);
  }
  return true;
}

bool Ladder_Generator_QT::MakeTrialLadder() {
  m_weight = 1.;
  for (size_t beam=0;beam<2;beam++) {
    m_y[beam][0]    = m_y[beam][1] = m_ylimits[beam];
    m_qt2prev[beam] = m_q[beam].PPerp2();  
  }
  size_t dir;
  TPropList::iterator pit = p_props->begin(); pit++;
  do {
    m_seff = (m_q[0]+m_q[1]).Abs2();
    if (m_seff<4.*m_qt2min) break;
    dir    = dabs(m_y[0][1])>dabs(m_y[1][1]) ? 0 : 1;
    if (TrialEmission(dir)) AddEmission(p_ladder,dir,pit);
  } while (m_y[0][0]>m_y[1][0]);
  if (LastEmissions(p_ladder)) {
    for (dir=0;dir<2;dir++) {
      p_ladder->AddRapidity(m_y[dir][1],m_flavs[dir],m_k[dir]);
    }
    p_props->insert(pit,T_Prop(colour_type::octet,m_qt, m_qt2min));
    m_weight *= RescaleLadder(p_ladder,m_qini[0]+m_qini[1]);
    return true;
  }
  return false;
}

bool Ladder_Generator_QT::TrialEmission(size_t dir,const double & yhat,const bool isrescatter) {
  double qt2min    = QT2Min(dir), qt2max = isrescatter ? m_seff/(4.*sqr(cosh(m_y[dir][1]))) : m_seff;
  double arg       = M_PI/(3.*AlphaS(0.)*dabs(log(qt2max/m_qt2min)));
  Form_Factor * ff = dabs(m_y[dir][1])>m_Ymax ? p_eikonal->FF(dir) : NULL;
  double weight, dy;
  do {
    MakeTransverseUnitVector();
    dy           = arg * log(ran->Get());
    m_y[dir][0] += (dir==1? -1. : 1.) * dy;
    m_qt         = MakePropMomentum(qt2min, qt2max, ff);
    m_k[dir]     = MakeFSMomentum(dir);
    weight       = ( LDCWeight(m_qt2, m_qt2prev[dir]) *
		     AlphaSWeight(m_k[dir].PPerp2()) *
		     (isrescatter ? 1. : AbsorptionWeight(m_k[dir], m_y[dir][0]+yhat)) );
    if ((m_q[0]+m_q[1]-m_k[dir]).Abs2()<4.*m_qt2min ||
	(dir==0 && m_y[0][0]<m_y[1][0]) ||
	(dir==1 && m_y[1][0]>m_y[0][0])) return false;
  } while (m_q[dir][0]-m_k[dir][0]<0. ||
	   weight<ran->Get());    
  return true;
}

void Ladder_Generator_QT::AddEmission(Ladder * ladder,size_t dir, TPropList::iterator & pit) {
  ladder->AddRapidity(m_y[dir][1],m_flavs[dir],m_k[dir]);
  ladder->GetProps()->insert(pit,T_Prop(colour_type::octet, m_qt, m_qt2min));
  if (dir==1) pit--;
  Vec4D old      = m_q[dir];
  m_q[dir]      -= m_k[dir];	
  m_y[dir][1]    = m_y[dir][0];
  m_qt2prev[dir] = m_qt2;
  if (dabs(m_y[dir][0])<m_Ymax) { m_y[dir][0] == (dir==0 ? m_Ymax : -m_Ymax); }
}

bool Ladder_Generator_QT::LastEmissions(Ladder * ladder) {
  if (ladder->GetEmissions()->size()==0) return FixSimpleKinematics();
  double qt2min    = QT2Min();
  Form_Factor * ff = NULL;
  if (dabs(m_y[0][1])>dabs(m_y[1][1]) && dabs(m_y[0][1])>m_Ymax)
    ff = p_eikonal->FF(0);
  if (dabs(m_y[1][1])>dabs(m_y[0][1]) && dabs(m_y[1][1])>m_Ymax)
    ff = p_eikonal->FF(1);
  size_t trials = 1000; 
  do {
    MakeTransverseUnitVector();
    m_qt = MakePropMomentum(qt2min, m_seff, ff); 
    for (size_t beam=0;beam<2;beam++) m_k[beam] = MakeFSMomentum(beam);
    if ((m_k[0]+m_k[1]).Abs2()<m_seff) {
      double weight = (AlphaSWeight(m_k[0].PPerp2()) *
		       AlphaSWeight(m_k[1].PPerp2())  *
		       LDCWeight(m_qt2, m_qt2prev[0]) *
		       LDCWeight(m_qt2, m_qt2prev[1]));
      if (weight>ran->Get()) return true; 
    }
  } while ((trials--)>0); 
  return false;
}

bool Ladder_Generator_QT::FixLastEmissions() {
  if (p_emissions->size()==0) return FixSimpleKinematics();
  Vec4D Q    = m_q[0]+m_q[1];
  double Q0  = Q[0], Q1 = Q[1], Q2 = Q[2], Q3 = Q[3], kt[2];
  double QT2 = Q1*Q1+Q2*Q2, QT = sqrt(QT2);
  for (size_t i=0;i<2;i++) {
    kt[i] = ( (Q0/cosh(m_y[1-i][1]) - Q3/sinh(m_y[1-i][1]))/
	      (cosh(m_y[i][1])/cosh(m_y[1-i][1]) -
	       sinh(m_y[i][1])/sinh(m_y[1-i][1])) );
  }
  if (QT<dabs(kt[0]-kt[1]) || QT>dabs(kt[0]+kt[1])) return false;
  double QTtilde  = (QT2+kt[0]*kt[0]-kt[1]*kt[1])/(2.*kt[0]);
  double QT2tilde = sqr(QTtilde), disc = QT2-QT2tilde;
  if (disc<0.) return false;
  double cosphi[2], sinphi[2];
  cosphi[0] = 1./QT2 * (QTtilde*Q1 + sqrt(disc)*Q2);
  cosphi[1] = (Q1-kt[0]*cosphi[0])/kt[1];
  for (size_t i=0;i<2;i++) sinphi[i] = sqrt(1.-sqr(cosphi[i]));
  for (size_t i=0;i<2;i++) {
    sinphi[0] *= -1.;
    for (size_t j=0;j<2;j++) {
      sinphi[1] *= -1.;
      if (dabs(Q2-(kt[0]*sinphi[0]+kt[1]*sinphi[1]))<1.e-6) break;
    }
  }
  m_qt     = (m_q[0]-m_k[0]).Perp(); 
  m_weight = (AlphaSWeight(m_k[0].PPerp2()) *
	      AlphaSWeight(m_k[1].PPerp2()) *
	      LDCWeight(m_qt2, m_qt2prev[0]) *
  	      LDCWeight(m_qt2, m_qt2prev[1]));
  return true;
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
}

void Ladder_Generator_QT::SelectPropagatorColours() {
  //msg_Out()<<"       - "<<METHOD<<"\n";
  // Iterate over propagators and assign colours different than octet
  LadderMap::iterator lit1=p_emissions->begin(),
    lit2 = lit1; lit2++;
  TPropList::iterator pit1=p_props->begin(), pit2=pit1;
  double y1,y2,wt1,wt8,ratio1,ratio2=0.;
  while (lit2!=p_emissions->end() && pit1!=p_props->end()) {
    y1     = lit1->first;
    y2     = lit2->first;
    wt1    = m_density.SingletWeight(y2,y1);
    wt8    = m_density.OctetWeight(y2,y1);
    ratio1 = wt1/(wt1+wt8); 
    if (ratio1>ran->Get()) pit1->SetCol(colour_type::singlet);
    if (pit1!=p_props->begin() &&
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
  m_weight = 1.; return;
  Vec4D Pcms       = (p_ladder->InPart(0)->Momentum() +
  		      p_ladder->InPart(1)->Momentum()); 
  double Y         = Pcms.Y(), SHat = Pcms.Abs2();
  double sigma_act = m_partonic.dSigma(SHat,Y);
  double start     = m_weight, me, regge;
  m_weight *= sigma_act/m_sigmahat;
  //m_weight *= me    = m_me(p_ladder,m_qt2min);
  m_weight *= regge = TotalReggeWeight();
}


void Ladder_Generator_QT::AddRescatters() {
  msg_Out()<<"====================================================================\n"
  	   <<(*p_ladder)<<"\n";
  LadderMap::iterator lit1=p_emissions->begin(), litN=lit1; litN++;
  TPropList::iterator pit=p_props->begin();
  while (lit1!=p_emissions->end() && litN!=p_emissions->end() &&
	 pit!=p_props->end()) {
    double prob = (pit->Col()!=colour_type::singlet ?
		   m_density.RescatterProbability(lit1->first,litN->first) : 0.);
    if (ran->Get()<prob && MakeRescatterLadder(lit1,pit)) {
      LadderMap::iterator olit1=lit1, olitN=litN;
      lit1 = litN; lit1++; 
      p_ladder->DeleteRapidity(olit1);p_ladder->DeleteRapidity(olitN);
      pit = p_ladder->GetProps()->erase(pit);
      p_ladder->UpdatePropagatorKinematics();
      msg_Out()<<"Inserted another ladder:\n"<<(*p_ladder)<<"\n";
    }
    else {
      lit1++; pit++;
    }
    litN=lit1; litN++;
  }
  msg_Out()<<"====================================================================\n\n";
}

bool Ladder_Generator_QT::MakeRescatterLadder(LadderMap::iterator & lit1,
					      TPropList::iterator & pit1) {
  LadderMap::iterator litN=lit1; litN++;
  TPropList::iterator pitN=pit1; pitN++;
  m_qini[0] = m_q[0] = lit1->second.Momentum();
  m_qini[1] = m_q[1] = litN->second.Momentum();
  Vec4D  Pcms  = m_q[0]+m_q[1];
  double shat  = Pcms.Abs2(), yhat = Pcms.Y();
  Poincare cms = Poincare(Pcms);
  for (size_t i=0;i<2;i++) cms.Boost(m_q[i]);
  Poincare zax = Poincare(m_q[0],Vec4D(1.,0.,0.,1.)); 
  Ladder rescatter(lit1->second.Position());
  TPropList::iterator pit=rescatter.GetProps()->begin(); pit++;
  for (size_t i=0;i<2;i++) {
    m_y[i][0]  = m_y[i][1] = m_q[i].Y();
    zax.Rotate(m_q[i]);
    rescatter.InPart(i)->SetFlavour(i==0?lit1->second.Flavour():litN->second.Flavour());
    rescatter.InPart(i)->SetMomentum(m_q[i]);
    m_flavs[i] = rescatter.InPart(i)->Flavour();
  }
  //msg_Out()<<"Add rescatter emissions from\n"
  //	   <<lit1->second<<litN->second
  //	   <<"to the ladder:\n"<<rescatter<<"\n";
  size_t dir;
  do {
    m_seff = (m_q[0]+m_q[1]).Abs2();
    if (m_seff<4.*m_qt2min) break;
    dir = dabs(m_y[0][1])>dabs(m_y[1][1]) ? 0 : 1;
    //msg_Out()<<" * Trial(dir = "<<dir<<") in ["<<m_y[0][1]<<", "<<m_y[1][1]<<"], "
    //	     <<"actual y = "<<m_y[dir][0]<<"\n"; 
    if (TrialEmission(dir,yhat,true)) {
      //msg_Out()<<"   --> new y = "<<m_y[dir][0]<<", qT = "<<m_qt<<"\n";
      AddEmission(&rescatter,dir,pit);
    }
  } while (m_y[0][0]>m_y[1][0]);
  if (rescatter.GetEmissions()->size()==0 || !LastEmissions(&rescatter)) return false;
  for (dir=0;dir<2;dir++) rescatter.AddRapidity(m_y[dir][1],m_flavs[dir],m_k[dir]);
  rescatter.GetProps()->insert(pit,T_Prop(colour_type::octet,m_qt, m_qt2min));
  
  m_weight *= RescaleLadder(&rescatter,m_qini[0]+m_qini[1]);
  for (LadderMap::iterator lit=rescatter.GetEmissions()->begin();
       lit!=rescatter.GetEmissions()->end();lit++) {
    Vec4D mom = lit->second.Momentum();
    zax.RotateBack(mom);
    cms.BoostBack(mom);
    lit->second.SetMomentum(mom);
  }
  for (size_t i=0;i<2;i++) rescatter.InPart(i)->SetMomentum(m_qini[i]);
  
  for (LadderMap::iterator lit=rescatter.GetEmissions()->begin();
       lit!=rescatter.GetEmissions()->end();lit++) {
    p_ladder->AddRapidity(lit->second.Momentum().Y(),
			  lit->second.Flavour(),lit->second.Momentum());
  }
  for (TPropList::iterator pit=rescatter.GetProps()->begin();
       pit!=rescatter.GetProps()->end();pit++) {
    p_ladder->GetProps()->insert(pit1,T_Prop(colour_type::octet,pit->Q(),pit->Q02()));
  }
  msg_Out()<<METHOD<<" adds "<<rescatter.GetEmissions()->size()<<" partons.\n";
  return true;
}

double Ladder_Generator_QT::QT2Min(size_t dir) {
  double qt2min = m_qt2min;
  for (size_t beam=0;beam<2;beam++) {
    if ((dir==2 || dir==beam) &&
	(ATOOLS::dabs(m_y[beam][1])>m_Ymax)) qt2min = m_qt2minFF;
  }
  return qt2min;
}

double Ladder_Generator_QT::QT2Max() {
  return m_seff/ATOOLS::sqr(cosh(ATOOLS::Max(ATOOLS::dabs(m_y[0][1]),
					     ATOOLS::dabs(m_y[1][1]))));
}

ATOOLS::Vec4D Ladder_Generator_QT::MakeFSMomentum(size_t dir) {
  double kt = sqrt((m_q[dir]+(dir==0 ? -m_qt : m_qt)).PPerp2());
  return ( kt * ATOOLS::Vec4D(cosh(m_y[dir][1]),0.,0.,sinh(m_y[dir][1])) +
	   (m_q[dir] + (dir==0 ? -m_qt : m_qt)).Perp());
}

ATOOLS::Vec4D Ladder_Generator_QT::MakePropMomentum(const double & qt2min,
						    const double & qt2max,
						    Form_Factor * ff) {
  if (qt2max<0.) exit(1);
  m_qt2 = ATOOLS::Max(0.,(ff ?
			  ff->SelectQT2(qt2max,qt2min) :
			  qt2min * pow(qt2max/qt2min, ATOOLS::ran->Get())));
  return sqrt(m_qt2) * m_eqt;
}

