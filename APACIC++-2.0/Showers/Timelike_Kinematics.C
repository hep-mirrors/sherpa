#include "Timelike_Kinematics.H"
#include "Run_Parameter.H"
#include "Poincare.H"
#include "Tree.H"
#include "Data_Read.H"
#include <iomanip>


using std::endl;

using namespace APACIC;
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace AORGTOOLS;




Timelike_Kinematics::Timelike_Kinematics(double _pt2min, Data_Read * const dataread) : 
  pt2min(_pt2min), t0(4.*pt2min),pt_scheme(1),mass_scheme(1)
{
  double ycut   = AORGTOOLS::rpa.gen.Ycut();
  m_type = 1;
  if (rpa.gen.Beam1().IsLepton() && rpa.gen.Beam2().IsLepton()) {
    msg.Debugging()<<" Jet_Finder in Timelike_Kinematics set up  to deal with lepton-lepton collisions "<<endl;
    m_type = 1;
  }
  else if ((!rpa.gen.Beam1().IsLepton() && !rpa.gen.Beam2().IsLepton())) {
    msg.Debugging()<<" Jet_Finder in Timelike_Kinematics set up  to deal with hadron-hadron collisions "<<endl;
    m_type = 4;
  }
  else {
    msg.Error()<<"Error in Timelike_Kinematics DIS is not yet implemented in the Jetfinder "<<endl;
    m_type = 4;
  }
  jf            = new APHYTOOLS::Jet_Finder(ycut,m_type); 
  if (m_type==1) m_losejet_veto = dataread->GetValue<int>("FS LOSEJETVETO",1);
  else m_losejet_veto = dataread->GetValue<int>("FS LOSEJETVETO",0);
}


//-----------------------------------------------------------------------
//------------------- Checks for kinematics : The shuffles --------------
//----------------------------------------------------------------------- 

bool Timelike_Kinematics::CheckZRange(Knot * mo) {
  Knot * d1=mo->left;
  Knot * d2=mo->right;
  
  if ((d1->stat == 3) || (d2->stat == 3)) return 1; 

  double t  = mo->t;  double t1 =d1->t;  double t2 =d2->t;

  if (t  < t1+t2+2.*sqrt(t1*t2)) {
    if (d1->stat==0 && d2->stat==0) return 0;
    if (d1->stat==0 && d2->stat!=0) {
      d2->stat=3;
      return 0;
    }
    if (d1->stat!=0 && d2->stat==0) {
      d1->stat=3;
      return 0;
    }
    if (d1->t > d2->t) d1->stat=3;
    else d2->stat=3;
    return 0;
  }

  double lambda   = sqrt(sqr(t-t1-t2)-4.*t1*t2);   
  double z        = mo->z;
  double r1 = (t+t2-t1-lambda)/(2.*t);
  double r2 = (t-t2+t1-lambda)/(2.*t);
  z     = z - r1*z + r2*(1.-z);
  double e1= z*sqrt(mo->E2);
  double p1=0;
  if (sqr(e1)>t1) p1=sqrt(sqr(e1)-t1);
  double e2= (1-z)*sqrt(mo->E2);
  double p2=0;
  if (sqr(e2)>t2) p2=sqrt(sqr(e2)-t2);

  double zp,zm,zdelta;  
  bool do1=0, do2=0;
  if (d1->stat) {
    zdelta=p1/e1;
    zm=0.5*(1.-zdelta);
    zp=0.5*(1.+zdelta);
    if ((d1->z<zm) || (zp<d1->z))  do1=1;
  }

  // check z-range of daugther two (unconstrained)
  if (d2->stat) {
    zdelta=p2/e2;
    zm=0.5*(1.-zdelta);
    zp=0.5*(1.+zdelta);
    if ((d2->z<zm) || (zp<d2->z))   do2=1;
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

bool Timelike_Kinematics::Shuffle(Knot * mo, int first)
{
  if (first) return ShuffleMoms(mo);
  return ShuffleZ(mo);
} 

bool Timelike_Kinematics::ShuffleZ(Knot * mo) 
{
  double t      = mo->t;
  double t1     = mo->left->t;
  double t2     = mo->right->t;

  if (t - (t1+t2+2.*sqrt(t1*t2)) < rpa.gen.Accu()) return 0; 

  double lambda = sqrt(sqr(t-t1-t2)-4.*t1*t2); 
  double r1     = (t+t2-t1-lambda)/(2.*t);
  double r2     = (t-t2+t1-lambda)/(2.*t);
  double z      = mo->z;
  mo->z         = z - r1*z + r2*(1.-z);

  double pt2 = z*(1.-z)*t;
  if (pt_scheme==1)       pt2 -= (1.-z)*t1 + z*t2;
  else if (pt_scheme==2)  pt2 = 0.25*Min((1.-z)/z,z/(1.-z))*t;
  if (pt2<pt2min) {
    mo->z = z;
  }
  if (KinCheck(0,mo)) {
    mo->z = z;
    return 0;
  }
  mo->left->E2  = mo->z*mo->z*mo->E2;
  mo->right->E2 = (1.-mo->z)*(1.-mo->z)*mo->E2;
  return 1;
}

bool Timelike_Kinematics::ShuffleMoms(Knot * mo) 
{ 
  Knot * d1       = mo->left;
  Knot * d2       = mo->right;

  double t        = mo->t;
  double t1       = d1->t;
  double t2       = d2->t;

  if (t - (t1+t2+2.*sqrt(t1*t2)) < rpa.gen.Accu()) return 0; 
  double lambda   = sqrt(sqr(t-t1-t2)-4.*t1*t2);   
  double z        = mo->z;

  double r1 = (t+t2-t1-lambda)/(2.*t);
  double r2 = (t-t2+t1-lambda)/(2.*t);
  mo->z     = z - r1*z + r2*(1.-z);

  if (dabs(mo->z-z) < rpa.gen.Accu()) return 1;

  if (KinCheck(1,mo)) {
    mo->z = z;
    return 0;
  }

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
    Vec4D p2a( (1.-r2)*p2 + r1*p1 );
    Vec4D p =mo->part->Momentum();
    mo->z= p1a[0]/p[0];

    if (KinCheck(1,mo)) {
      mo->z = z;
      return 0;
    }
  }

  d1->part->SetMomentum( (1.-r1)*p1 + r2*p2 );
  d2->part->SetMomentum( (1.-r2)*p2 + r1*p1 );
  d1->E2   = mo->z*mo->z*mo->E2;
  d2->E2   = (1.-mo->z)*(1.-mo->z)*mo->E2;
  return 1;
}

bool Timelike_Kinematics::KinCheck(int first,Knot * mo) 
{
  // KinCheck returns 1 in case the kinematics does not work out,
  //                  0 in case everything is fine.
  Knot * d1      = mo->left;
  Knot * d2      = mo->right;
  // no daughters no checks
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
    if (jetveto) {      
      double pt2 = mo->z*(1.-mo->z)*mo->t;
      double tb  = d1->tout;         
      double tc  = d2->tout;         
      if (pt_scheme == 1) 
	pt2       -= (1.-mo->z)*tb + mo->z*tc;
      else if (pt_scheme == 2)
	pt2 = 0.25*Min((1.-mo->z)/mo->z,mo->z/(1.-mo->z))*mo->t;
      double pt2th    = sqrt(pt2/mo->E2)/(mo->z*(1.- mo->z));
      double crudeth  = sqrt( mo->t/(mo->z*(1.- mo->z)*mo->E2) );
      
      double coscrude = cos(pt2th);

      int hit=0;
      if (jf->TwoJets(mo->E2,mo->z,coscrude,0)) return 1;
    }
    return 0;
  }
  else return 0;


  // already known daughters to be checked
  if (first) { 
    // one might have grand daughters with already fixed momenta
    if (d1->stat!=0) { if (KinCheck(first,d1)) return 1; }
    if (d2->stat!=0) { if (KinCheck(first,d2)) return 1; }
    return 0;
  }
 
  if (KinCheck(first,d1)) return 1; 
  if (KinCheck(first,d2)) return 1;
  return 0;
}


bool Timelike_Kinematics::ExtraJetCheck(Knot * mo, Knot * d1, Knot * d2) 
{
  if (!m_losejet_veto) return 1;

  if (m_type==4) {
    if (d1==0) {
      if (! (jf->TwoJets(d2->part->Momentum()))) return 0;
      return 1;
    }
    
    if (d2==0) {
      if (! (jf->TwoJets(d1->part->Momentum()))) return 0;
      return 1;
    }
    
    if (! (jf->TwoJets(d1->part->Momentum(),d2->part->Momentum()))) return 0;
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
  if (pt_scheme == 1) 
    pt2 -= (1.-z)*tb + z*tc;
  else if (pt_scheme == 2)
    pt2 = 0.25*Min((1.-z)/z,z/(1.-z))*t;
  double pt2th  = sqrt(pt2/E2)/(z*(1.- z));
  double crudeth  = sqrt( t/(z*(1.- z)*E2) );
  
  double coscrude = cos(pt2th);
  
  int hit=0;
  if (! (jf->TwoJets(E2,z,coscrude,0))) return 0;
  return 1;
}


bool Timelike_Kinematics::JetVeto(double mo_t, double mo_e2, double mo_z, 
				  double tb, double tc) 
{
  if (jetveto) {      
    double pt2 = mo_z*(1.-mo_z)*mo_t;
    if (pt_scheme == 1) 
      pt2       -= (1.-mo_z)*tb + mo_z*tc;
    else if (pt_scheme == 2)
      pt2 = 0.25*Min((1.-mo_z)/mo_z,mo_z/(1.-mo_z))*mo_t;

    double pt2th    = sqrt(pt2/mo_e2)/(mo_z*(1.- mo_z));
    double crudeth  = sqrt( mo_t/(mo_z*(1.- mo_z)*mo_e2) );
   
    double coscrude = cos(pt2th); 
    
    if (jf->TwoJets(mo_e2,mo_z,coscrude,0)) {
      return 1;
    }
  }
  return 0;
}

//-----------------------------------------------------------------------
//--------------------- Evaluation of the kinematics --------------------
//----------------------------------------------------------------------- 
 
bool Timelike_Kinematics::DoKinematics(Knot * mo) 
{
  if (!(mo)) return 1;
  if (!(mo->left)) {
    if (mo->part->Info()==' ') {
      mo->part->SetStatus(1);
      mo->part->SetInfo('F');
    }
    return 1;
  }
  
  double t      = mo->t, z=mo->z, E2=mo->E2;
  double p      = sqrt(E2-t);
  double t1     = mo->left->t, w1 = mo->left->E2;
  double p1     = sqrt(w1-t1);
  double p1real = mo->left->part->Momentum()[1]; 
  double t2     = mo->right->t,w2 = mo->right->E2;
  double p2     = sqrt(w2-t2);
  double p2real = mo->right->part->Momentum()[1]; 
  
  if (p1real==0.) {
    mo->part->SetStatus(2);
    Knot * au = mo->prev->left;
    int sign  = 0;
    if (mo==au) {
      au      = mo->prev->right;
      sign    = 1;
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
    }
    if (n1abs<1.e-5) {
      n1         = cross(Vec3D(0.,1.,0.),nm);
      n1abs      = n1.Abs();
    }
    
    n1           = n1/n1abs;
    if (sign) n1 = -1.*n1; 
    Vec3D n2     = cross(nm,n1);
    
    double sph=sin(mo->phi),cph=cos(mo->phi);
    Vec3D es   = cph*n1 + sph*n2;
    
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
    if (error) msg.Error()<<"Error after CheckVector() !"<<endl;
    if ( (mo->part->Momentum()-mo->left->part->Momentum()-
	  mo->right->part->Momentum()).Abs2() > rpa.gen.Accu()) error = 1;

    double py=mo->left->part->Momentum()[2];
    if (!(py<0) && !(py>0)) error =1;


    if (error) {
      msg.Error()<<"Error in Timelike_Kinematics."<<endl
		 <<"Timelike_Kinematics::DoKinematics : After moms set"<<endl
		 <<"    Mother & Daughters :"<<endl
		 <<"    mother : "<<t<<","<<z<<","<<E2<<endl
		 <<"    d1     : "<<t1<<","<<w1<<" ("<<(z*z*E2)<<")"<<endl
		 <<"    d2     : "<<t2<<","<<w2<<" ("<<((1.-z)*(1.-z)*E2)<<")"<<endl
		 <<"    no, p "<<mo->kn_no<<", "<<mo->part->Momentum()<<", "<<endl
		 <<"      "<<mo->part->Momentum().Abs2()<<"("<<mo->t<<")"<<endl
		 <<"    flav : "<<mo->part->Flav()<<endl
		 <<"    no, p "<<mo->left->kn_no<<", "<<mo->left->part->Momentum()<<", "<<endl
		 <<"      "<<mo->left->part->Momentum().Abs2()
		 <<"("<<mo->left->t<<", "<<mo->left->E2<<", "<<p1<<") "<<endl
		 <<"    flav : "<<mo->left->part->Flav()<<endl
		 <<"    no, p "<<mo->right->kn_no<<", "<<mo->right->part->Momentum()<<", "<<endl
		 <<"      "<<mo->right->part->Momentum().Abs2()
		 <<"("<<mo->right->t<<", "<<mo->right->E2<<", "<<p2<<") "<<endl
		 <<"    flav : "<<mo->right->part->Flav()<<endl
		 <<"    Three Vectors :"<<endl
		 <<"    "<<es<<endl<<"    "<<na<<endl<<"    "<<nm<<endl
		 <<"    Just for fun : prev :"<<endl
		 <<"    no, p "<<mo->prev->kn_no<<", "<<mo->prev->part->Momentum()<<", "<<endl
		 <<"      "<<mo->prev->part->Momentum().Abs2()<<"("<<mo->prev->t<<")"<<endl
		 <<"    flav : "<<mo->prev->part->Flav()<<endl;
      if (mo->prev->left == mo) {
	msg.Error()<<"    Just for fun : prev->right :"<<endl
		   <<"    no, p "<<mo->prev->right->kn_no<<", "
		   <<mo->prev->right->part->Momentum()<<", "<<endl
		   <<"      "<<mo->prev->right->part->Momentum().Abs2()
		   <<"("<<mo->prev->right->t<<")"<<endl
		   <<"    flav : "<<mo->prev->right->part->Flav()<<endl;
      }
      else {
	msg.Error()<<"    Just for fun : prev->left :"<<endl
		   <<"    no, p "<<mo->prev->left->kn_no<<", "
		   <<mo->prev->left->part->Momentum()<<", "<<endl
		   <<"      "<<mo->prev->left->part->Momentum().Abs2()
		   <<"("<<mo->prev->left->t<<")"<<endl
		   <<"    flav : "<<mo->prev->left->part->Flav()<<endl;
      }
      return 0;
    }
  }

  if (!DoKinematics(mo->left))  return 0;
  if (!DoKinematics(mo->right)) return 0;

  BoostDaughters(mo);
  return 1;
}


bool Timelike_Kinematics::ArrangeColourPartners(Knot * au,Knot * d1,Knot * d2) {
  if (!au) return 0;
  if (!d1) return 0;
  if (!d2) return 0;
  if (jf->PTij(au->part->Momentum(),d1->part->Momentum()) <
      jf->PTij(au->part->Momentum(),d2->part->Momentum()) ) {
    return 0;
  }
  return 1;
}




bool Timelike_Kinematics::CheckVector(Vec4D vec) {
  if ( (vec.Abs2() > 0) && (vec.Abs2() < 0) ) return 1;
  if (vec[0] < 0) return 1;
  return 0;
}
 
 
void Timelike_Kinematics::BoostDaughters(Vec4D pold, Vec4D pnew, 
					 const Vec4D & pmom, Knot * mo) {

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

void Timelike_Kinematics::BoostDaughters(Knot * mo) {
  Knot * d1= mo->left;
  Knot * d2= mo->right;
  Vec4D p=d1->part->Momentum()+d2->part->Momentum();
  if (d1->left) {
    Vec4D p1old = d1->left->part->Momentum() + d1->right->part->Momentum();
    Vec4D p1new = d1->part->Momentum();
    if (p1new != p1old) {
      BoostDaughters(p1old,p1new,p,d1);
    }
  }
  if (d2->left) {
    Vec4D p2old = d2->left->part->Momentum() + d2->right->part->Momentum();
    Vec4D p2new = d2->part->Momentum();
    if (p2new != p2old) {
      BoostDaughters(p2old,p2new,p,d2);
    }
  }

}
