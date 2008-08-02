#include "Dipole_Splitter.H"
#include "Hadronisation_Parameters.H"
#include "Random.H"

using namespace AHADIC;
using namespace ATOOLS;

size_t Dipole::s_cnt=0;

Dipole_Splitter::Dipole_Splitter(Strong_Coupling * as,const double ptmax) :
  p_as(as),p_constituents(hadpars.GetConstituents()), p_options(NULL),
  p_spect(0), p_split(0), p_out1(0), p_out2(0), 
  m_pt2veto(sqr(ptmax)), m_asmax(p_as->MaxValue()), 
  m_2p3min(2.*p_constituents->MinMass())
{ }


bool Dipole_Splitter::SplitCluster(SP(Cluster) cluster,const double pt2max) {
  SP(Dipole) dip1 = new Dipole(new Proto_Particle((*cluster->GetTrip())),
			       new Proto_Particle((*cluster->GetAnti())));
#ifdef memchecker
  msg_Out()<<"### Two new Proto_Particles ("
  	   <<dip1->Triplet()<<"/"<<dip1->AntiTriplet()<<") from "<<METHOD<<"."<<std::endl;
#endif

  p_dip = dip1;
  if (msg->LevelIsDebugging()) {
    msg_Out()<<METHOD<<" : "<<cluster->Momentum()<<" ("<<cluster->Mass2()<<")."<<std::endl;
    p_dip->Output();
  }
  
  SP(Dipole) dip2 = EmitGluon(pt2max);
  if (dip2==NULL) {
    msg_Tracking()<<"ERROR in "<<METHOD<<":"<<std::endl
		  <<"   Could not emit gluon, hence no dipole splitting."<<std::endl
		  <<"   Return false, may lead to new event."<<std::endl;
#ifdef memchecker
    msg_Out()<<"### Delete Proto_Particle ("<<dip1->Triplet()
    	     <<"/"<<dip1->Triplet()->m_flav<<") in "<<METHOD<<"."<<std::endl;
#endif
#ifdef memchecker
    msg_Out()<<"### Delete Proto_Particle ("<<dip1->AntiTriplet()
    	     <<"/"<<dip1->AntiTriplet()->m_flav<<") in "<<METHOD<<"."<<std::endl;
#endif
    return false;
  }
  dip1->SetAntiTriplet(dip2->Triplet());
  dip1->Update();
  dip2->Update();

  if (msg->LevelIsDebugging()) {
    msg_Out()<<"       after gluon emission: "<<std::endl;
    dip1->Output();
    dip2->Output();
  }


  SP(Dipole) first(dip1), second(dip2);
  if (ran.Get()<(dip1->Mass2()-dip1->Triplet()->m_mom.Abs2())/
      (dip1->Mass2()-dip1->Triplet()->m_mom.Abs2()+
       dip2->Mass2()-dip2->AntiTriplet()->m_mom.Abs2())) {
    first  = dip2;
    second = dip1;
  }
  if (!SplitDipole(first,pt2max) && !SplitDipole(second,pt2max)) {
    msg_Tracking()<<"Error in "<<METHOD<<" :"<<std::endl
		  <<"   Two unsplittable dipoles emerging from :"<<std::endl
		  <<(*cluster)<<std::endl;
    return false;
  }


  dip1->SetAntiTriplet(p_out1);
  dip2->SetTriplet(p_out2);
  
  dip1->Update();
  dip2->Update();

  if (msg->LevelIsDebugging()) {
    msg_Out()<<"       after gluon splitting: "<<std::endl;
    dip1->Output();
    dip2->Output();
    msg_Out()<<"       Total momentum : "
	     <<(dip1->Momentum()+dip2->Momentum())<<" ("
	     <<(dip1->Momentum()+dip2->Momentum()).Abs2()<<")."<<std::endl;
  }

  SP(Cluster) left  = new Cluster(dip1->Triplet(),dip1->AntiTriplet());
  SP(Cluster) right = new Cluster(dip2->Triplet(),dip2->AntiTriplet());
#ifdef memchecker
  msg_Out()<<"@@@ Two new clusters "<<left<<"/"<<right<<" from "<<METHOD<<"."<<std::endl;
#endif

  cluster->SetLeft(left);
  cluster->SetRight(right);
  left->SetPrev(cluster);
  right->SetPrev(cluster);

  return true;
}

Dipole *Dipole_Splitter::EmitGluon(const double pt2max) {
  SetSpectatorAndSplitter();

  m_M2      = p_dip->Mass2();
  m_flav    = Flavour(kf_gluon);
  m_m1      = p_constituents->Mass(p_spect->m_flav);
  m_m12     = sqr(m_m1);
  m_m2      = p_constituents->Mass(m_flav);
  m_m22     = sqr(m_m2);
  m_m3      = p_constituents->Mass(p_split->m_flav);
  m_m32     = sqr(m_m3);
  m_2p3     = m_m3;
  m_xt2min  = sqr(m_m1+m_2p3min/2.)*sqr(m_m3+m_2p3min/2.)/m_M2;
  if (pt2max>0.) m_pt2veto = pt2max;
  m_xt2     = m_xt2max = 1./4.;
  m_y       = m_ybound = -1./2. * log(m_xt2min);  
  m_pref    = 4./3.;

  if (!DetermineSplitting(false)) {
    msg_Tracking()<<"Error in "<<METHOD<<":"<<std::endl
		  <<"   Could not split dipole, must retry event."<<std::endl;
    return NULL;
  }
  if (FixKinematics(false)) {
    return new Dipole(p_out1,p_dip->AntiTriplet());
  }
  return NULL;
}

bool Dipole_Splitter::SplitDipole(SP(Dipole) dip,const double pt2max) {
  p_dip = dip;

  SetSpectatorAndSplitter();
  m_M2      = p_dip->Mass2();
  m_m1      = p_constituents->Mass(p_spect->m_flav);
  if (sqrt(m_M2)-m_m1-m_2p3min<=1.e-8) return false;
  m_m12     = sqr(m_m1);
  m_xt2min  = sqr((m_m1+m_2p3min/2.)*(m_2p3min)/m_M2);
  if (pt2max>0.) m_pt2veto = pt2max;
  m_xt2     = m_xt2max = 1./4.;
  m_y       = m_ybound = -1./2. * log(m_xt2min);  
  m_pref    = 0.;
  m_flav    = Flavour(kf_none);
  double maxwt(0.), wt;
  for (FDIter fdit=p_options->begin();fdit!=p_options->end();fdit++) {
    if (fdit->second->popweight<=0. || 
	m_m1+fdit->second->massmin>sqrt(m_M2)) continue;
    m_pref += wt = fdit->second->popweight;
    if (wt>maxwt) maxwt = wt;
  }
  m_pref /= maxwt;
  if (maxwt==0.) {
    if (msg->LevelIsDebugging()) {
      msg_Out()<<METHOD<<" : cannot decay this dipole: "<<std::endl;
      p_dip->Output();
    }
    return false;
  }

  if (!DetermineSplitting()) {
    msg_Tracking()<<"Error in "<<METHOD<<":"<<std::endl
		  <<"   Could not split dipole, must retry event."<<std::endl;
    return false;
  }
  if (FixKinematics()) {
#ifdef memchecker
    msg_Out()<<"### Delete Proto_Particle ("<<p_split
    	     <<"/"<<p_split->m_flav<<") in "<<METHOD<<"."<<std::endl;
#endif
    return true;
  }
  msg_Tracking()<<"ERROR in "<<METHOD<<" FixKinematics did not work out!"<<std::endl;
  return false;
}

void Dipole_Splitter::SetSpectatorAndSplitter() {
  if (p_dip->Triplet()->m_flav.IsGluon() && 
      (p_dip->AntiTriplet()->m_flav.IsQuark() ||
       p_dip->AntiTriplet()->m_flav.IsDiQuark())) {
    p_spect = p_dip->AntiTriplet(); 
    p_split = p_dip->Triplet();
    p_dip->SetSwitched(true);
  }
  else if (p_dip->AntiTriplet()->m_flav.IsGluon() && 
	   (p_dip->Triplet()->m_flav.IsQuark() ||
	    p_dip->Triplet()->m_flav.IsDiQuark())) {
    p_split = p_dip->AntiTriplet(); 
    p_spect = p_dip->Triplet();
    p_dip->SetSwitched(false);
  }
  else if ((p_dip->Triplet()->m_flav.IsGluon() && 
	    p_dip->AntiTriplet()->m_flav.IsGluon()) ||
	   (!p_dip->Triplet()->m_flav.IsGluon() && 
	    !p_dip->AntiTriplet()->m_flav.IsGluon())) {
    if (ran.Get()>0.5) {
      p_spect = p_dip->AntiTriplet(); 
      p_split = p_dip->Triplet();
      p_dip->SetSwitched(true);
    }
    else {
      p_split = p_dip->AntiTriplet(); 
      p_spect = p_dip->Triplet();
      p_dip->SetSwitched(false);
    }
  }
}

bool Dipole_Splitter::DetermineSplitting(bool glusplit) {
  do { 
    if (!SelectPT_Y(glusplit)) {
      if (!MinimalDecay(glusplit)) return false;
      break; 
    } 
  } while (Veto());
  return true;
}

bool Dipole_Splitter::SelectPT_Y(bool glusplit) {
  double expo = sqr(log(m_xt2))-log(ran.Get())/(m_asmax*m_pref);
  if (expo<0.) {
    return false;
  }
  m_xt2 = exp(-sqrt(expo));
  if (m_xt2<m_xt2min) {
    return false;
  }
  m_y   = m_ybound*(2.*ran.Get()-1.);
  if (glusplit) {
    FDIter fdit;
    double disc = 0.9999999*m_pref*ran.Get();
    for (fdit=p_options->begin();fdit!=p_options->end();fdit++) {
      if (fdit->second->popweight<=0. || 
	  m_m1+fdit->second->massmin>sqrt(m_M2)) continue;
      disc -= fdit->second->popweight;
      if (disc<0) { 
	m_flav = fdit->first; 
	m_2p3 = fdit->second->massmin; 
	m_m2  = m_m3 = m_2p3/2.;
	m_m22 = m_m32 = sqr(m_m2);
	break; 
      }
    }
    if (m_flav.IsDiQuark()) m_flav = m_flav.Bar();
  }
  return true;
}

bool Dipole_Splitter::Veto() {
  if (m_flav==Flavour(kf_none))        return true;
  double ytruebound(log(1./(2.*sqrt(m_xt2))*(1.+sqrt(1.-4.*m_xt2))));
  if (dabs(m_y)>ytruebound)             return true;
  CalculateInvariants();
  if (!ConstructKinematics())           return true;
  if (m_asmax*ran.Get()>(*p_as)(m_pt2)) return true;
  if (m_pt2>m_pt2veto)                  return true;
  return false;
}

bool Dipole_Splitter::MinimalDecay(bool glusplit) {
  if (glusplit) {
    double disc = 0.9999999*m_pref*ran.Get();
    for (FDIter fdit=p_options->begin();fdit!=p_options->end();fdit++) {
      if (fdit->second->popweight<=0. || 
	  m_m1+fdit->second->massmin>sqrt(m_M2)) continue;
      disc -= fdit->second->popweight;
      if (disc<0) { m_flav = fdit->first; m_2p3 = fdit->second->massmin; break; }
    }
    if (m_flav.IsDiQuark()) m_flav = m_flav.Bar();
    m_m2  = m_m3 = m_2p3/2.;
    m_m22 = m_m32 = sqr(m_m2);
  }
  else {
    m_2p3 = m_2p3min+m_m3;
  }

  m_s23 = sqr(m_2p3);
  m_s12 = m_m12+m_m22+Max(m_m2,m_2p3min/2.)/m_2p3*(m_M2-m_m12-m_s23);

  if (!ConstructKinematics()) {
    if (msg->LevelIsTracking()) {
      msg_Out()<<"ERROR in "<<METHOD<<" : "<<std::endl
	       <<"   Could not construct any decay for dipole below."<<std::endl
	       <<"   This will lead to a new event."<<std::endl;
      p_dip->Output();
    }
    return false;
  }
  return true;
}

void Dipole_Splitter::CalculateInvariants() {
  double Mp(sqrt(m_xt2)*m_M2);
  m_s12 = Mp*exp(-m_y);
  m_s23 = Mp*exp(+m_y);
}

bool Dipole_Splitter::ConstructKinematics() {
  if (m_s12<sqr(m_m1+m_m2) || m_s23<sqr(m_m2+m_m3)) {
    return false;
  }
  double Mred2 = m_M2-(m_s23+m_m12);     
  if (sqr(Mred2)<4.*m_s23*m_m12 || Mred2<0.) {
    return false;
  }
  double gamma = (Mred2+sqrt(sqr(Mred2)-4.*m_s23*m_m12))/2.;
  double a1    = Max(m_m12,0.)/gamma;
  double a23   = m_s23/gamma;
  double lpm   = sqrt(gamma)/2.;
  Vec4D lplus  = lpm*Vec4D(1.,0.,0.,1.);
  Vec4D lminus = lpm*Vec4D(1.,0.,0.,-1.);
  Vec4D mom23  = lplus+a23*lminus;
  double z2    = (m_s12-(m_m12+m_m22)-a1*(m_s23-(m_m32-m_m22)))/(gamma-a1*m_s23);
  double z3    = 1.-z2;
  m_pt2        = z2*(1.-z2)*m_s23-(1.-z2)*m_m22-z2*m_m32;
  if (ATOOLS::IsZero(m_pt2)) m_pt2 = 0.;
  if (m_pt2<0. || m_pt2>m_pt2veto) {
    return false;
  }                                   
  double phi   = 2.*M_PI*ran.Get();
  Vec4D  kperp = sqrt(m_pt2)*Vec4D(0.,cos(phi),sin(phi),0.);
  double y2    = (m_m22+m_pt2)/(z2*m_s23);
  double y3    = (m_m32+m_pt2)/(z3*m_s23);

  m_mom1       = lminus+a1*lplus;
  m_mom2       = z2*lplus+y2*a23*lminus+kperp;
  m_mom3       = z3*lplus+y3*a23*lminus-kperp;
  return true;
}

bool Dipole_Splitter::FixKinematics(bool glusplit) {
  Vec4D q1(p_split->m_mom),q2(p_spect->m_mom),q(q1+q2);
  Vec4D k(m_mom1+m_mom2+m_mom3);

#ifdef AHAmomcheck
  Vec4D checkbef(q);
#endif

  Poincare boostk(k);
  boostk.Boost(m_mom1);
  boostk.Boost(m_mom2);
  boostk.Boost(m_mom3);
  
  Poincare boost(q);
  boost.Boost(q1);
  boost.Boost(q2);
  Poincare rotate(q1,Vec4D(1.,Vec3D::ZVEC));
  rotate.Rotate(q1);
  rotate.Rotate(q2);

  rotate.RotateBack(m_mom1);
  rotate.RotateBack(m_mom2);
  rotate.RotateBack(m_mom3);
  boost.BoostBack(m_mom1);
  boost.BoostBack(m_mom2);
  boost.BoostBack(m_mom3);

#ifdef AHAmomcheck
  Vec4D checkaft = m_mom1+m_mom2+m_mom3;

  if (dabs((checkbef-checkaft).Abs2())>1.e-12) {
    msg_Out()<<METHOD<<" yields momentum violation : "<<std::endl
 	     <<checkbef<<" - "<<checkaft<<" --> "<<(checkbef-checkaft).Abs2()<<std::endl;
  }
  else msg_Debugging()<<METHOD<<" conserves momentum."<<std::endl;
#endif

  if (glusplit) {
    p_out1 = new Proto_Particle(m_flav.Bar(),m_mom2,'l');
    p_out2 = new Proto_Particle(m_flav,m_mom3,'l');
#ifdef memchecker
    msg_Out()<<"### Two new Proto_Particles ("<<p_out1<<"/"<<p_out2
    	     <<" -- "<<m_flav.Bar()<<"/"<<m_flav<<") from "<<METHOD<<"."<<std::endl;
#endif
  }
  else {
    p_out1 = new Proto_Particle(m_flav,m_mom2,'l');
#ifdef memchecker
    msg_Out()<<"### New Proto_Particle ("<<p_out1<<"/"<<m_flav<<") from "<<METHOD<<"."<<std::endl;
#endif
    p_out2 = NULL;
    p_split->m_mom = m_mom3;
  }
  p_spect->m_mom = m_mom1;
  return true;
}
