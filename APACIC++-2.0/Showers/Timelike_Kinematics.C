#include "Timelike_Kinematics.H"
#include "Run_Parameter.H"
#include "Poincare.H"
#include "Tree.H"
#include "Data_Read.H"
#include "Exception.H"
#include <iomanip>

using namespace APACIC;
using namespace ATOOLS;


Timelike_Kinematics::Timelike_Kinematics(double pt2min, Data_Read * const dataread) : 
  m_pt2min(pt2min), m_pt_scheme(1), m_zrange_scheme(1), m_mass_scheme(1)
{
  m_type = 1;
  if (rpa.gen.Beam1().IsLepton() && rpa.gen.Beam2().IsLepton()) {
    m_type = 1;
  }
  else if ((!rpa.gen.Beam1().IsLepton() && !rpa.gen.Beam2().IsLepton())) {
    m_type = 4;
  }
  else {
    msg.Error()<<"ERROR in Timelike_Kinematics : "<<std::endl
	       <<"   DIS is not yet implemented in the jetfinder, continue with hh-mode."<<std::endl;
    m_type = 4;
    throw(ATOOLS::Exception(ATOOLS::ex::not_implemented,"DIS is not implemented yet",
			    "Timelike_Kinematics","Timelike_Kinematics"));
  }
  p_jf = new ATOOLS::Jet_Finder(rpa.gen.Ycut(),m_type); 
  p_jf->SetDeltaR(rpa.gen.DeltaR()); 


    /*   0=constrained z,  1=unconstrained z)  */
  m_zrange_scheme = dataread->GetValue<int>("FS_Z_RANGES",1);      

  m_losejet_veto = dataread->GetValue<int>("FS_LOSEJETVETO",0);
}


//-----------------------------------------------------------------------
//------------------- Checks for kinematics : The shuffles --------------
//----------------------------------------------------------------------- 

bool Timelike_Kinematics::CheckZRange(Knot const * const  mo,
				      Flavour const * const d1_flav,
				      Flavour const * const d2_flav) const
{
  Knot * d1=mo->left;
  Knot * d2=mo->right;
  
  if ((d1->stat == 3) || (d2->stat == 3)) return 1; 

  double t   = mo->t;  
  double t1  = d1->t, t2  = d2->t;
  double t01 = 0.,    t02 = 0.;

  if (m_zrange_scheme==0) {
    t01 = d1->tout;
    t02 = d2->tout;
  }

  if (t < t1+t2+2.*sqrt(t1*t2)) {
    if (d1->stat==0 && d2->stat==0) return 0;
    if (d1->stat==0 && d2->stat!=0) {
      d2->stat = 3;
      return 0;
    }
    if (d1->stat!=0 && d2->stat==0) {
      d1->stat = 3;
      return 0;
    }
    if (d1->t > d2->t) d1->stat=3;
    else d2->stat = 3;
    return 0;
  }

  double z  = CalcZShift(mo->z,t,t1,t2,t01,t02);
  double e1 = z*sqrt(mo->E2);
  double e2 = (1-z)*sqrt(mo->E2);

  bool do1 = 0, do2 = 0;
  if (d1->stat) {
    if (m_zrange_scheme==1  && 
	!CheckZRange(d1->z,sqr(e1),t1,0.,0.)) do1=1;
    if (m_zrange_scheme==0  &&
	!CheckZRange(d1->z,sqr(e1),t1,d1_flav[0].PSMass(),d1_flav[1].PSMass())) do1=1;
  }

  // check z-range of daugther two (unconstrained)
  if (d2->stat) {
    if (m_zrange_scheme==1 && 
	!CheckZRange(d2->z,sqr(e2),t2,0.,0.)) do2=1;
    if (m_zrange_scheme==0  &&
	!CheckZRange(d2->z,sqr(e2),t2,d2_flav[0].PSMass(),d2_flav[1].PSMass())) do2=1;
  }

  // if necessary select a daughter to be diced again.
  if (!do1 && !do2) return 1; // all fine

  if (!do1 && do2) {
    d2->stat=3;
    return 0; // dice d2 again!
  }

  if (do1 && !do2) {
    d1->stat=3;
    return 0; // dice d1 again!
  }

  if (d1->t > d2->t) d1->stat=3;
  else d2->stat=3;

  return 0;
}

int Timelike_Kinematics::Shuffle(Knot * const mo, const int first) const
{
  if (first) return ShuffleMoms(mo);
  return ShuffleZ(mo);
} 

int Timelike_Kinematics::ShuffleZ(Knot * const mo) const
{
  double t   = mo->t;
  double t1  = mo->left->t;
  double t2  = mo->right->t;
  double t01 = 0.,    t02 = 0.;

  if (m_zrange_scheme==0) {
    t01 = mo->left->tout;
    t02 = mo->right->tout;
  }

  if (t - (t1+t2+2.*sqrt(t1*t2)) < rpa.gen.Accu()) {
    msg_Debugging()<<"Timelike_Kinematics::ShuffleZ() : not enough virtuality: "<<sqrt(t)<<" < "<<sqrt(t1)<<" + "<<sqrt(t2)<<std::endl;
    return 0; 
  }

  double z      = mo->z;
  mo->z         = CalcZShift(z,t,t1,t2,t01,t02);

  //  double  pt2 = mo->z*(1.-mo->z)*t - (1.-mo->z)*t1 - mo->z*t2;
  double  pt2 = 0.25*CalcKt2(mo->z,mo->E2,t,t1,t2);
  /*
  double pt2 = z*(1.-z)*t;
  if (m_pt_scheme==1) pt2 = mo->z*(1.-mo->z)*t - (1.-mo->z)*t1 - mo->z*t2;
 //       pt2 -= (1.-z)*t1 + z*t2;
  else if (m_pt_scheme==2)  {
    pt2 = 0.25*Min((1.-z)/z,z/(1.-z))*t;
    double kt2 = 0.25*CalcKt2(mo->z,mo->E2,t,t1,t2);
    double pt2a = z*(1.-z)*t;
    double pt2b = mo->z*(1.-mo->z)*t - (1.-mo->z)*t1 - mo->z*t2;
    std::cout<<" pt2="<<pt2<<"   kt2="<<kt2<<std::endl;
    std::cout<<" at2="<<pt2a<<"   bt2="<<pt2b<<std::endl;
  }
  */
  if (pt2<m_pt2min) {
    msg_Debugging()<<" WARNING Timelike_Kinematics::ShuffleZ() pt2 limit"<<std::endl;
    mo->z = z;
    return 0;
  }
  int stat = KinCheck(0,mo);
  if (stat) {
    msg_Debugging()<<" WARNING Timelike_Kinematics::ShuffleZ() : kinCheck failed "<<mo->z<<" "<<z<<std::endl;
    mo->z = z;
    if (stat==3) return 3;
    return 0;
  }
  mo->left->E2  = mo->z*mo->z*mo->E2;
  mo->right->E2 = (1.-mo->z)*(1.-mo->z)*mo->E2;
  return 1;
}

int Timelike_Kinematics::ShuffleMoms(Knot * const mo) const
{ 
  // modifies already existing momenta
  Knot * d1       = mo->left;
  Knot * d2       = mo->right;

  double t        = mo->t;
  double t1       = d1->t;
  double t2       = d2->t;

  if (dabs(t - mo->part->Momentum().Abs2())>rpa.gen.Accu()) {
    msg.Out()<<" WARNING Timelike_Kinematics::ShuffleMoms : \n"
	     <<"    strong mass deviation "<<t<<" vs. "<<mo->part->Momentum().Abs2()<<std::endl;
  }

  if (t - (t1+t2+2.*sqrt(t1*t2)) < rpa.gen.Accu()) {
    msg_Debugging()<<" WARNING Timelike_Kinematics::ShuffleMoms() : not enough virtuality: "<<sqrt(t)<<" < "<<sqrt(t1)<<" + "<<sqrt(t2)<<std::endl;
    return 0; 
  }

  double r1 = 0., r2 = 0.;
  double z  = mo->z;

  Vec4D p1 = d1->part->Momentum();
  Vec4D p2 = d2->part->Momentum();

  if ((p1.Abs2()/d1->E2 > rpa.gen.Accu()) ||
      (p2.Abs2()/d2->E2 > rpa.gen.Accu())) {
    double t1n = p1.Abs2();
    double t2n = p2.Abs2();
    double A   =  ((t2 - t2n) - (t1 -t1n)) / (t + t1n - t2n);
    double B   =  (t + t2n-t1n)/ (t + t1n -t2n);
    double C   =  t - t1n - t2n;
    double D   =  (2.*t2n - 2. * A*B*t1n + (A - B)*C)/(2.*(t2n + B*B*t1n - B*C));
    double E   =  (t2n - t2 + A*A*t1n + A*C)/(t2n + B*B*t1n - B*C);

    r2 = D - sqrt(D*D- E);
    r1 = A + r2*B;
    Vec4D p1a( (1.-r1)*p1 + r2*p2 );
    Vec4D p =mo->part->Momentum();
    mo->z = p1a[0]/p[0];

    if (KinCheck(1,mo)) {
      msg_Debugging()<<"Timelike_Kinematics::ShuffleMoms() : kinCheck failed "<<mo->z<<" "<<z<<std::endl;
      mo->z = z;
      return 0;
    }
  }
  else {
    double lambda = sqrt(sqr(t-t1-t2)-4.*t1*t2); 
    r1     = (t+t2-t1-lambda)/(2.*t);
    r2     = (t-t2+t1-lambda)/(2.*t);
    mo->z         = z - r1*z + r2*(1.-z);
    /*
      mo->z         = CalcZShift(z,t,t1,t2);
    */

    if (dabs(mo->z-z) < rpa.gen.Accu()) return 1;

    if (KinCheck(1,mo)) {
      //    msg_Debugging()<<" WARNING Timelike_Kinematics::ShuffleMoms() : kinCheck failed "<<mo->z<<" "<<z<<std::endl;
      mo->z = z;
      return 0;
    }
  }

  msg_Debugging()<<"Timelike_Kinematics::ShuffleMoms(["<<mo->kn_no<<"])"<<std::endl;

  d1->part->SetMomentum( (1.-r1)*p1 + r2*p2 );
  d2->part->SetMomentum( (1.-r2)*p2 + r1*p1 );
  d1->E2   = mo->z*mo->z*mo->E2;
  d2->E2   = (1.-mo->z)*(1.-mo->z)*mo->E2;
  msg_Debugging()<<" d1 : "<<sqr(d1->part->Momentum()[0])<<" == "<<d1->E2<<std::endl;
  msg_Debugging()<<" d2 : "<<sqr(d2->part->Momentum()[0])<<" == "<<d2->E2<<std::endl;

  // boost daughters if existent
  BoostDaughters(mo);

  // update daughter E2,z 
  Tree::UpdateDaughters(mo);

  return 1;
}

int Timelike_Kinematics::KinCheck(const int first,Knot * const mo) const
{
  // KinCheck returns 1 in case the kinematics does not work out,
  //                  0 in case everything is fine.
  Knot * d1      = mo->left;
  Knot * d2      = mo->right;
  // no daughters no checks -> OK!
  if ((d1==0) || (d2==0)) return 0;
  
  double w1      = mo->z*mo->z*mo->E2;
  double w2      = (1.-mo->z)*(1.-mo->z)*mo->E2;
  double t1      = d1->t;
  double t2      = d2->t;
  if ((t1>w1) || (t2>w2)) return 1;

  double p1p2    = sqrt((w1-t1)*(w2-t2));

  // triangular three momementum relation     
  if (mo->E2-mo->t - (w1-t1 + w2-t2 + 2.*p1p2) > first*rpa.gen.Accu() ) return 1;

  double cosreal = (2.*mo->z*(1.-mo->z)*mo->E2-mo->t+t1+t2)/(2.*p1p2); 
  // physical opening angle
  if ((dabs(cosreal) > 1.) && !(first)) return 1;
  if (cosreal > 1.)  cosreal = 1.; 
  if (cosreal < -1.) cosreal = -1.; 
  mo->costh = cosreal;

  // physical deflection angle
  if (cosreal >  1.) mo->costh = 1.;
  if (cosreal < -1.) mo->costh = -1.;
  double coth1;
  if (dabs(mo->E2-mo->t)<rpa.gen.Accu()) coth1 = -1.;
  else coth1 = (mo->costh*p1p2+w1-t1)/(sqrt((mo->E2-mo->t)*(w1-t1)));

  if (dabs(coth1) > 1.+rpa.gen.Accu()) return 1;


  if (!first) {
    // Test for extra Jet
    if (m_jetveto) {      
      double pt2 = mo->z*(1.-mo->z)*mo->t;
      double tb  = d1->tout;         
      double tc  = d2->tout;         
      if (m_pt_scheme == 1) 
	pt2       -= (1.-mo->z)*tb + mo->z*tc;
      else if (m_pt_scheme == 2)
	pt2 = 0.25*Min((1.-mo->z)/mo->z,mo->z/(1.-mo->z))*mo->t;
      double pt2th    = sqrt(pt2/mo->E2)/(mo->z*(1.- mo->z));
      // double crudeth  = sqrt( mo->t/(mo->z*(1.- mo->z)*mo->E2) );      
      double coscrude = cos(pt2th);

      if (p_jf->TwoJets(mo->E2,mo->z,coscrude,0)) {
	// if (p_jf->TwoJets(mo->E2,mo->z,cosreal,0)) {
	msg_Debugging()<<" JetVeto in Timelike_Kinematics::KinCheck "<<std::endl;
	return 3;
      }
    }
    return 0;
  }
  return 0;
}


bool Timelike_Kinematics::ExtraJetCheck(Knot const * const mo, Knot const * const d1, Knot const * const d2) const 
{
  if (!m_losejet_veto) return 1;

  if (m_type==4) {
    if (d1==0) {
      if (! (p_jf->TwoJets(d2->part->Momentum()))) return 0;
      return 1;
    }
    
    if (d2==0) {
      if (! (p_jf->TwoJets(d1->part->Momentum()))) return 0;
      return 1;
    }
    
    if (! (p_jf->TwoJets(d1->part->Momentum(),d2->part->Momentum()))) return 0;
    return 1;
  }
  
  double z, E2, t;
  if (mo) {
    z  = mo->z;
    E2 = mo->E2;
    t  = mo->t;
  } 
  else {
    Vec4D mom = d1->part->Momentum() + d2->part->Momentum();
    E2   = sqr(mom[0]);
    z    = d1->part->Momentum()[0]/mom[0];
    t    = mom.Abs2();
  }

  double pt2 = z*(1.-z)*t;
  double tb  = d1->tout;         
  double tc  = d2->tout;         
  if (m_pt_scheme == 1) 
    pt2 -= (1.-z)*tb + z*tc;
  else if (m_pt_scheme == 2) {
    pt2 = 0.25*Min((1.-z)/z,z/(1.-z))*t;
  }
  double pt2th  = sqrt(pt2/E2)/(z*(1.- z));
  //double crudeth  = sqrt( t/(z*(1.- z)*E2) );
  double coscrude = cos(pt2th);


  /*
  double t1      = d1->t;
  double t2      = d2->t;
  double w1      = z*z*E2;
  double w2      = (1.-z)*(1.-z)*E2;
  double p1p2    = sqrt((w1-t1)*(w2-t2));
  double cosreal = (2.*z*(1.-z)*E2 - t+t1+t2)/(2.*p1p2); 
  */  

  //if (! (p_jf->TwoJets(E2,z,cosreal,0))) {
  if (! (p_jf->TwoJets(E2,z,coscrude,0))) {
    msg_Debugging()<<" ExtraJetVeto in Timelike_Kinematics::ExtraJetCheck "<<std::endl;

    return 0;
  }
  return 1;
}


bool Timelike_Kinematics::JetVeto(double t, double e2, double z, 
				  double t1, double t2) 
{
  bool flag=0;
  if (t<0.) {
    t=-t;
    flag=1;
  }
  if (m_jetveto) {      
    double pt2 = z*(1.-z)*t;
//     if (m_pt_scheme == 1) 
//       pt2       -= (1.-z)*t1 + z*t2;
//     else if (m_pt_scheme == 2)
//       pt2 = 0.25*Min((1.-z)/z,z/(1.-z))*t;

    double pt2th    = sqrt(pt2/e2)/(z*(1.- z));
    //double crudeth  = sqrt( t/(z*(1.- z)*e2) );
    double coscrude = cos(pt2th); 
    
    if (p_jf->TwoJets(e2,z,coscrude,0)) {
      return 1;
    }
  }
  return 0;
}

//-----------------------------------------------------------------------
//--------------------- Evaluation of the kinematics --------------------
//----------------------------------------------------------------------- 
 
bool Timelike_Kinematics::DoKinematics(Knot * const mo) const
{
  msg_Debugging()<<"Timelike_Kinematics::DoKinematics([<<"<<mo->kn_no<<"])"<<std::endl;
  if (!(mo)) return 1;
  if (!(mo->left)) {
    if (mo->part->Info()==' ') {
      mo->part->SetStatus(1);
      mo->part->SetInfo('F');
    }
    return 1;
  }
  
  double t      = mo->t, z = mo->z, E2 = mo->E2;
  double p      = sqrt(E2-t);
  double t1     = mo->left->t, w1 = mo->left->E2;
  double p1     = sqrt(w1-t1);
  double p1real = mo->left->part->Momentum()[1]; 
  double t2     = mo->right->t,w2 = mo->right->E2;
  double p2     = sqrt(w2-t2);
  
  if (p1real==0.) {
    mo->part->SetStatus(2);
    Knot * au = mo->prev->left;
    int sign  = 0;
    int mode  = 1;
    if (mo==au) {
      au      = mo->prev->right;
      sign    = 1;
      mode    = 3;
    }
    Vec3D na(au->part->Momentum()); // aunt
    Vec3D nm(mo->part->Momentum()); // mother
    na = na/na.Abs();
    nm = nm/nm.Abs();
    
    Vec3D n1     = cross(na,nm);
    double n1abs = n1.Abs();
    if (n1abs<1.e-5) {
      n1         = cross(Vec3D(0.,0.,1.),nm);
      n1abs      = n1.Abs();
      mode = mode|4;
    }
    if (n1abs<1.e-5) {
      n1         = cross(Vec3D(0.,1.,0.),nm);
      n1abs      = n1.Abs();
      mode = mode|8;
    }
    
    n1           = n1/n1abs;
    if (sign) n1 = -1.*n1; 
    Vec3D n2     = cross(nm,n1);
    
    double phi=mo->phi + mo->polinfo.Angle();
    //    double sph=sin(mo->phi),cph=cos(mo->phi);
    double bph=cos(phi),cph=-sin(phi);
    Vec3D es   = cph*n1 + bph*n2;
    bool test=false;
    if (test) {
      msg.Out()<<"==============="<<std::endl;
      msg.Out()<<" nm ="<<nm<<"\n"
	       <<" n1 ="<<n1<<"\n"
	       <<" n2 ="<<n2<<"\n"
	       <<" phi="<<mo->phi<<" ("<<180.*mo->phi/M_PI<<")\n"
	       <<" na="<<na<<"("<<mode<<")\n";
      msg.Out()<<" es="<<es<<"("<<acos(es*n2)<<")\n";
      msg.Out()<<"test = "<<cross(nm,n1)*n2<<"\n";
    }


    double cth1 = (p*p-p2*p2+p1*p1)/(2.*p*p1);
    double sth1 = sqrt(1.-sqr(cth1));
    double cth2 = (p*p+p2*p2-p1*p1)/(2.*p*p2);
    double sth2 = sqrt(1.-sqr(cth2));
    mo->costh   = cth1*cth2-sth1*sth2;
 
    Vec3D p1vec = p1*(cth1*nm - sth1*es);
    mo->left->part->SetMomentum(Vec4D(sqrt(w1),p1vec));
    
    Vec3D p2vec = p2*(cth2*nm + sth2*es);
    mo->right->part->SetMomentum(Vec4D(sqrt(w2),p2vec));
    
    bool error = (CheckVector(mo->part->Momentum()) || 
                  CheckVector(mo->left->part->Momentum()) || 
                  CheckVector(mo->right->part->Momentum()));
    if (error) msg.Error()<<"Error after CheckVector() !"<<std::endl;
    if ( dabs((mo->part->Momentum()-mo->left->part->Momentum()-
	       mo->right->part->Momentum()).Abs2()) > rpa.gen.Accu()) error = 1;

    double py=mo->left->part->Momentum()[2];
    if (!(py<0) && !(py>0)) error =1;


    if (error) {
      Vec3D pm(mo->part->Momentum());
      msg.Out()<<" p = "<<nm*pm<<" = "<<pm.Abs()<<" = "<<p<<std::endl;
      msg.Out()<<" p_l = "<<p1*cth1<<" + "<<p2*cth2<<" = "<<(p1*cth1 + p2*cth2)<<std::endl;
      msg.Out()<<" p_t = "<<p1*sth1<<" - "<<p2*sth2<<" = "<<(p1*sth1-p2*sth2)<<std::endl;

      msg.Error()<<"Error in Timelike_Kinematics."<<std::endl
		 <<"Timelike_Kinematics::DoKinematics : After moms set"<<std::endl
		 <<"    Mother & Daughters :"<<std::endl
		 <<"    mother : "<<t<<","<<z<<","<<E2<<std::endl
		 <<"    d1     : "<<t1<<","<<w1<<" ("<<(z*z*E2)<<")"<<std::endl
		 <<"    d2     : "<<t2<<","<<w2<<" ("<<((1.-z)*(1.-z)*E2)<<")"<<std::endl
		 <<"    no, p "<<mo->kn_no<<", "<<mo->part->Momentum()<<", "<<std::endl
		 <<"      "<<mo->part->Momentum().Abs2()<<"("<<mo->t<<")"<<std::endl
		 <<"    flav : "<<mo->part->Flav()<<std::endl
		 <<"    no, p "<<mo->left->kn_no<<", "<<mo->left->part->Momentum()<<", "<<std::endl
		 <<"      "<<mo->left->part->Momentum().Abs2()
		 <<"("<<mo->left->t<<", "<<mo->left->E2<<", "<<p1<<") "<<std::endl
		 <<"    flav : "<<mo->left->part->Flav()<<std::endl
		 <<"    no, p "<<mo->right->kn_no<<", "<<mo->right->part->Momentum()<<", "<<std::endl
		 <<"      "<<mo->right->part->Momentum().Abs2()
		 <<"("<<mo->right->t<<", "<<mo->right->E2<<", "<<p2<<") "<<std::endl
		 <<"    flav : "<<mo->right->part->Flav()<<std::endl
		 <<"    Three Vectors :"<<std::endl
		 <<"    "<<es<<std::endl<<"    "<<na<<std::endl<<"    "<<nm<<std::endl
		 <<"    Just for fun : prev :"<<std::endl
		 <<"    no, p "<<mo->prev->kn_no<<", "<<mo->prev->part->Momentum()<<", "<<std::endl
		 <<"      "<<mo->prev->part->Momentum().Abs2()<<"("<<mo->prev->t<<")"<<std::endl
		 <<"    flav : "<<mo->prev->part->Flav()<<std::endl;
      if (mo->prev->left == mo) {
	msg.Error()<<"    Just for fun : prev->right :"<<std::endl
		   <<"    no, p "<<mo->prev->right->kn_no<<", "
		   <<mo->prev->right->part->Momentum()<<", "<<std::endl
		   <<"      "<<mo->prev->right->part->Momentum().Abs2()
		   <<"("<<mo->prev->right->t<<")"<<std::endl
		   <<"    flav : "<<mo->prev->right->part->Flav()<<std::endl;
      }
      else {
	msg.Error()<<"    Just for fun : prev->left :"<<std::endl
		   <<"    no, p "<<mo->prev->left->kn_no<<", "
		   <<mo->prev->left->part->Momentum()<<", "<<std::endl
		   <<"      "<<mo->prev->left->part->Momentum().Abs2()
		   <<"("<<mo->prev->left->t<<")"<<std::endl
		   <<"    flav : "<<mo->prev->left->part->Flav()<<std::endl;
      }
      return 0;
    }
  }

  if (!DoKinematics(mo->left))  return 0;
  if (!DoKinematics(mo->right)) return 0;

  BoostDaughters(mo);
  return 1;
}


bool Timelike_Kinematics::ArrangeColourPartners(Particle const * const aup,Knot const * const d1,Knot const * const d2) const
{
  if (!aup) return 0;
  if (!d1) return 0;
  if (!d2) return 0;
  if (p_jf->MTij2(aup->Momentum(),d1->part->Momentum()) <
      p_jf->MTij2(aup->Momentum(),d2->part->Momentum()) ) {
    return 0;
  }
  return 1;
}




bool Timelike_Kinematics::CheckVector(const Vec4D mom) const 
{
  if ( (mom.Abs2() > 0) && (mom.Abs2() < 0) ) return 1;
  if (mom[0] < 0) return 1;
  return 0;
}
 
 
void Timelike_Kinematics::BoostDaughters(Vec4D pold, Vec4D pnew, 
					 const Vec4D & pmom, Knot * const mo) const
{
  msg_Debugging()<<"Timelike_Kinematics::BoostDaughters(...,["<<mo->kn_no<<"]"<<std::endl;

  int bigboost=0;
  Vec3D prot = cross(Vec3D(pnew),Vec3D(pold));

  Poincare bmom;  
  if (prot.Abs()>rpa.gen.Accu()) { 
    bigboost = 1;
    bmom = Poincare(pmom);

    Tree::BoRo(bmom,mo->left);
    Tree::BoRo(bmom,mo->right);
    bmom.Boost(pold);
    bmom.Boost(pnew);

    bmom=Poincare(Vec4D(pmom[0],-1.*Vec3D(pmom)));
  }
  Poincare bos1(pnew);
  Poincare bos(bos1*pold);
  Tree::BoRo(bos,mo->left);
  Tree::BoRo(bos,mo->right);
  
  if (bigboost) {
    Tree::BoRo(bmom,mo->left);
    Tree::BoRo(bmom,mo->right);
  }
}

void Timelike_Kinematics::BoostDaughters(Knot * const mo) const
{
  Knot * d1= mo->left;
  Knot * d2= mo->right;
  Vec4D p=d1->part->Momentum()+d2->part->Momentum();
  if (d1->left) {
    Vec4D p1old = d1->left->part->Momentum() + d1->right->part->Momentum();
    Vec4D p1new = d1->part->Momentum();
    if (p1new!=p1old && p1old!=Vec4D(0.,0.,0.,0.)) {
      BoostDaughters(p1old,p1new,p,d1);
    }
  }
  if (d2->left) {
    Vec4D p2old = d2->left->part->Momentum() + d2->right->part->Momentum();
    Vec4D p2new = d2->part->Momentum();
    if (p2new!=p2old && p2old!=Vec4D(0.,0.,0.,0.)) {
      BoostDaughters(p2old,p2new,p,d2);
    }
  }

}


bool  Timelike_Kinematics::CheckZRange(double z, double E2, double t, double t1, double t2)
{
  double mean_z  = 0.5 * ( 1. + (t1-t2)/t ); 
  double delta_z = 0.5 * sqrt((1.-t/E2) * ( sqr(t-t1-t2) - 4.*t1*t2 ))/t;
  
  if ((z<mean_z - delta_z) || (mean_z + delta_z<z)) {
    msg_Debugging()<<" Timelike_Kinematics::CheckZRange:"<<mean_z - delta_z<<" < "<<z<<" < "<<mean_z + delta_z<<std::endl;
    return false;
  }

  return true;
}

double Timelike_Kinematics::CalcZShift(double z, double t, double t1, double t2, double t01, double t02)
{
  double lambda1 = sqr(t-t1-t2)-4.*t1*t2; 
  if (lambda1<0) {
    msg_Tracking()<<" WARNING: Timelike_Kinematics::CalcZShift kinematics does not fit!"<<std::endl;
    return 0.5;
  }
  if (t01==0. && t02==0.) 
    return ((2.*z-1.)*sqrt(lambda1) + (t+ t1-t2))/(2.*t);

  double lambda0 = sqr(t-t01-t02)-4.*t01*t02; 
  return  (z-(t+t01-t02)/(2.*t))*sqrt(lambda1/lambda0) + (t+t1-t2)/(2.*t);
}

double Timelike_Kinematics::CalcKt2(double z, double E2, double t, double t1, double t2)
{
  double w1      = z*z*E2;
  double w2      = (1.-z)*(1.-z)*E2;
  double p1p2    = sqrt((w1-t1)*(w2-t2));
  double cosreal = (2.*z*(1.-z)*E2 - t+t1+t2)/(2.*p1p2); 
  
  double kt2 = 2.*ATOOLS::Min(w1,w2)*(1. - cosreal);

  return kt2;
}
