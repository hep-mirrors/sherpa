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
    if (FixInitialPartons() && MakeTrialLadder(p_ladder)) {
      ConstructISKinematics();
      SelectPropagatorColours(p_ladder);
      CalculateWeight();
    }
    else m_weight = 0.;
  } while (m_weight<ran->Get());    
  AddRescatters(p_ladder);
  return p_ladder;
}

void Ladder_Generator_QT::AddRescatters(Ladder * ladder) {
  //msg_Out()<<"====================================================================\n";
  //	   <<(*ladder)<<"\n";
  LadderMap::iterator lit[2]; 
  TPropList::iterator pit;
  bool success = true;
  while (success) {
    success = false;
    lit[0]=ladder->GetEmissions()->begin();
    while (lit[0]!=ladder->GetEmissions()->end()) {
      lit[1] = ladder->GetEmissions()->end(); lit[1]--;
      pit    = ladder->GetProps()->end(); pit--;
      while (lit[1]!=lit[0] && pit!=ladder->GetProps()->begin()) {
	double prob = 0.;
	if (!(lit[0]->second.IsMarked() || lit[1]->second.IsMarked()) &&
	    !(dabs(lit[0]->first)>m_Ymax && dabs(lit[1]->first)>m_Ymax) ) {
	  prob = m_density.RescatterProbability(lit[0]->first,lit[1]->first);
	}
	//msg_Out()<<"==== P_resc("<<lit[0]->first<<", "<<lit[1]->first<<") = "<<prob<<" --> ";
	if (ran->Get()<prob) {
	  //msg_Out()<<"try rescatter.\n";
	  Ladder * rescatter = MakeRescatterLadder(lit,pit);
	  if (rescatter) {
	    //msg_Out()<<"     Insert another ladder with "<<rescatter->GetEmissions()->size()<<" partons.\n";
	    MergeLadders(ladder,rescatter,lit,pit);
	    delete rescatter;
	    //if (ladder->CheckFourMomentum()) msg_Out()<<"Merge successful: "<<(*ladder)<<"\n";
	    //else exit(1);
	    success = true;
	    break;
	  }
	  //else msg_Out()<<"     Couldn't construct a rescatter ladder.\n";
	}
	//else msg_Out()<<"no rescatter.\n";
	lit[1]--;
	pit--;
      }
      if (success) break;
      lit[0]++;
    }
  }
  //msg_Out()<<"====================================================================\n\n";
}

void Ladder_Generator_QT::MergeLadders(Ladder * ladder,Ladder * rescatter,
				       LadderMap::iterator lit[2],
				       TPropList::iterator pit) {
  if (ladder==NULL || rescatter==NULL) return;
  for (LadderMap::iterator rlit=rescatter->GetEmissions()->begin();
       rlit!=rescatter->GetEmissions()->end();rlit++) {
    p_ladder->AddRapidity(rlit->second.Momentum().Y(),
			  rlit->second.Flavour(),rlit->second.Momentum())->SetMark(true);
  }
  for (TPropList::iterator rpit=rescatter->GetProps()->begin();
       rpit!=rescatter->GetProps()->end();rpit++) {
    p_ladder->GetProps()->insert(pit,T_Prop(colour_type::octet,rpit->Q(),rpit->Q02()));
  }
  //msg_Out()<<METHOD<<" adds "<<rescatter->GetEmissions()->size()<<" partons.\n";
  LadderMap::iterator olit[2];
  for (size_t i=0;i<2;i++) olit[i] = lit[i];
  lit[0] = lit[1]; lit[0]++; 
  for (size_t i=0;i<2;i++) ladder->DeleteRapidity(olit[i]);
  pit = ladder->GetProps()->erase(pit);
  ladder->UpdatePropagatorKinematics();
}


bool Ladder_Generator_QT::FixInitialPartons() {
  p_emissions->clear(); p_props->clear();
  m_shat     = m_partonic.MakeEvent();
  m_yhat     = 0;
  if (m_shat<0.) return false;
  m_sigmahat = m_partonic.SigmaHat();
  for (size_t beam=0;beam<2;beam++) {
    m_ylimits[beam] = (beam==0 ? 1.: -1.) * (m_Ymax + ran->Get()*m_deltaY);
    m_y[beam][0]    = m_y[beam][1] = m_ylimits[beam];
    m_qini[beam]    = m_q[beam] = 2.*m_partonic.X(beam) * rpa->gen.PBeam(beam);
    m_qt2prev[beam] = m_q[beam].PPerp2();  
    m_flavs[beam]   = m_partonic.Flav(beam);
  }
  return true;
}

Ladder * Ladder_Generator_QT::InitializeRescatterLadder(Ladder_Particle * parts[2]) {
  Ladder * rescatter = new Ladder(parts[0]->Position(),true);
  Vec4D Pcms(0.,0.,0.,0.);
  for (size_t i=0;i<2;i++) Pcms += m_qini[i] = m_q[i] = parts[i]->Momentum();
  m_shat = Pcms.Abs2();
  m_yhat = 0.; 
  for (size_t i=0;i<2;i++) {
    m_y[i][0]  = m_y[i][1] = m_q[i].Y();
    rescatter->InPart(i)->SetFlavour(m_flavs[i] = parts[i]->Flavour());
    rescatter->InPart(i)->SetMomentum(m_q[i]);
  }
  return rescatter;
}

Ladder * Ladder_Generator_QT::MakeRescatterLadder(LadderMap::iterator lit[2],
						  TPropList::iterator & pit) {
  Ladder_Particle * parts[2];
  for (size_t i=0;i<2;i++) parts[i] = &lit[i]->second;
  size_t trials = 100;
  Ladder * rescatter = InitializeRescatterLadder(parts);
  while (true) {
    for (size_t i=0;i<2;i++) {
      m_q[i]    = parts[i]->Momentum();
      m_y[i][0] = m_y[i][1] = m_q[i].Y();
    }
    m_weight = 1.;
    rescatter->GetEmissions()->clear();
    rescatter->GetProps()->clear();
    if (MakeTrialLadder(rescatter)) {
      m_weight *= m_me(rescatter);
      if (ran->Get()<m_weight) break;
    }
    if ((trials--)<=0) { delete rescatter; return NULL; }
  }
  for (size_t i=0;i<2;i++) rescatter->InPart(i)->SetMomentum(m_qini[i]);
  if (!rescatter->CheckFourMomentum()) exit(1);
  SelectPropagatorColours(rescatter);
  return rescatter;
}

bool Ladder_Generator_QT::MakeTrialLadder(Ladder * ladder) {
  m_weight = 1.;
  size_t dir;
  TPropList::iterator pit = ladder->GetProps()->begin(); pit++;
  //msg_Out()<<"     trial ladder for "<<m_q[0]<<"+ "<<m_q[1]<<".\n";
  do {
    m_seff = (m_q[0]+m_q[1]).Abs2();
    if (m_seff<4.*m_qt2min) break;
    dir    = dabs(m_y[0][1])>dabs(m_y[1][1]) ? 0 : 1;
    if (TrialEmission(ladder,dir,m_yhat)) AddEmission(ladder,dir,pit);
  } while (m_y[0][0]>m_y[1][0]);
  if ((ladder->IsRescatter() && ladder->GetEmissions()->size()==0) ||
      !LastEmissions(ladder)) return false;
  for (dir=0;dir<2;dir++) ladder->AddRapidity(m_y[dir][1],m_flavs[dir],m_k[dir]);
  ladder->GetProps()->insert(pit,T_Prop(colour_type::octet,m_qt, m_qt2min));
  m_weight *= RescaleLadder(ladder,m_qini[0]+m_qini[1]);
  return true;
}

bool Ladder_Generator_QT::TrialEmission(Ladder * ladder,size_t dir,const double & yhat) {
  bool   isrescatter = ladder->IsRescatter();
  double qt2min      = isrescatter ? QT2Min(dir) : QT2Min(dir);
  double qt2max      = (isrescatter ? m_seff : m_seff)/cosh(m_y[dir][1]);
  double arg         = M_PI/(3.*AlphaSMax()*dabs(log(qt2max/m_qt2min)));
  Form_Factor * ff   = dabs(m_y[dir][1])>m_Ymax ? p_eikonal->FF(dir) : NULL;
  double weight, dy, pdfratio, random;
  do {
    MakeTransverseUnitVector();
    dy           = arg * log(ran->Get());
    m_y[dir][0] += (dir==1? -1. : 1.) * dy;
    m_qt         = MakePropMomentum(qt2min, qt2max, ff, isrescatter);
    m_k[dir]     = MakeFSMomentum(dir);
    pdfratio     = (isrescatter ? 1. :
		    m_me.PDFratio(m_q[dir],Flavour(kf_gluon),m_q[dir]-m_k[dir],Flavour(kf_gluon),dir));
    weight       = ( LDCWeight(m_qt2, m_qt2prev[dir]) *
		     AlphaSWeight(m_k[dir].PPerp2()) *
		     pdfratio *
		     AbsorptionWeight(m_k[dir], m_y[dir][0]+yhat) );
    /*
      if (ladder->GetEmissions()->size()>=3 && m_k[dir].PPerp2()>4.)
      msg_Out()<<"      - try("<<ladder->GetEmissions()->size()<<") "<<m_y[dir][0]<<": check = "
    	       <<((m_q[0]+m_q[1]-m_k[dir]).Abs2()<4.*m_qt2min ||
    		  (dir==0 && m_y[0][0]<m_y[1][0]) ||
    		  (dir==1 && m_y[1][0]>m_y[0][0]))<<", "
    	       <<"weight = "<<weight<<" for "<<m_k[dir]<<"\n"
    	       <<"                weight = "<<AlphaSWeight(m_k[dir].PPerp2())<<" * "
    	       <<pdfratio<<" * "<<AbsorptionWeight(m_k[dir], m_y[dir][0]+yhat)<<"\n";
    */
    if ((m_q[0]+m_q[1]-m_k[dir]).Abs2()<4.*m_qt2min ||
	(dir==0 && m_y[0][0]<m_y[1][0]) ||
	(dir==1 && m_y[1][0]>m_y[0][0])) return false;
  } while (m_q[dir][0]-m_k[dir][0]<0. || weight<ran->Get());
  //if (ladder->GetEmissions()->size()>=3 && m_k[dir].PPerp2()>4.)
  //msg_Out()<<"       --> accept emission: "<<m_q[dir]<<" -> "<<m_k[dir]<<".\n";
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
  if (dabs(m_y[dir][0])<m_Ymax) { m_y[dir][0] = (dir==0 ? m_Ymax : -m_Ymax); }
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

void Ladder_Generator_QT::SelectPropagatorColours(Ladder * ladder) {
  //msg_Out()<<"       - "<<METHOD<<"\n";
  // Iterate over propagators and assign colours different than octet
  LadderMap::iterator lit1=ladder->GetEmissions()->begin(),
    lit2 = lit1; lit2++;
  TPropList::iterator pit1=ladder->GetProps()->begin(), pit2=pit1;
  double y1,y2,wt1,wt8,ratio1,ratio2=0.;
  while (lit2!=ladder->GetEmissions()->end() && pit1!=p_props->end()) {
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
  //Vec4D Pcms       = (p_ladder->InPart(0)->Momentum() +
  //		      p_ladder->InPart(1)->Momentum()); 
  //double Y         = Pcms.Y(), SHat = Pcms.Abs2();
  //double sigma_act = m_partonic.dSigma(SHat,Y);
  //double start     = m_weight, me, regge;
  //m_weight *= sigma_act/m_sigmahat;
  m_weight *= m_me(p_ladder,0.)/m_sigmahat;
  //m_weight *= TotalReggeWeight(p_ladder);
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

ATOOLS::Vec4D Ladder_Generator_QT::
MakePropMomentum(const double & qt2min,const double & qt2max,Form_Factor * ff,const bool rescatter) {
  if (qt2max<0.) exit(1);
  if (ff) m_qt2 = ff->SelectQT2(qt2max,qt2min);
  else    m_qt2 = qt2min * (pow(1.+qt2max/qt2min, ATOOLS::ran->Get())-1.);
  return sqrt(m_qt2) * m_eqt;
}

