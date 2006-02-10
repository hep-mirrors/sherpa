#include "Kinematics_Base.H"
#include "Random.H"

using namespace CS_SHOWER;
using namespace ATOOLS;
using namespace std;

const Vec3D Kinematics_Base::s_ex(Vec3D(1.,0.,0.));
const Vec3D Kinematics_Base::s_ey(Vec3D(0.,1.,0.));
const Vec3D Kinematics_Base::s_ez(Vec3D(0.,0.,1.));


Parton * Kinematics_FF::MakeKinematics(Parton * split, ATOOLS::Flavour & newfl)
{
  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum();
  Vec3D n1 = Vec3D(p1), n2 = Vec3D(p2), n3 = cross(n1,n2);
  if (n3.Sqr()<1.e-12) {
    double phi = 2.*M_PI*ran.Get();
    double cost = 1.-2.*ran.Get(), sint = sqrt(1.-cost*cost);
    Vec3D axis = sint*(cos(phi)*s_ex+sin(phi)*s_ey)+cost*s_ez;
    n3 = cross(n1,axis);
    if (n3.Sqr()<1.e-12) {
      n3 = cross(n2,axis);
      if (n3.Sqr()<1.e-12) {
	Vec3D axis = sint*(cos(phi)*s_ex-sin(phi)*s_ey)+cost*s_ez;
	n3 = cross(n1,axis);
	if (n3.Sqr()<1.e-12) {
	  msg.Error()<<"Error in Kinematics_FF::MakeKinematics : "<<endl
		     <<"   Could not construct an axis for the splitting : "<<endl
		     <<"   "<<n1<<" "<<n2<<"."<<endl;
	  abort();
	}
      }
    }
  }
  Vec4D n_perp = Vec4D(0.,n3/n3.Abs());

  double kt = sqrt(split->KtTest()), z = split->ZTest(), y = split->YTest();

  Vec4D 
    q1 = z*p1      + y*(1.-z)*p2 + kt*n_perp, 
    q2 =               (1.-y)*p2,
    q3 = (1.-z)*p1 +      y*z*p2 - kt*n_perp;

  if ((IsZero(q1[0]) && IsZero(q1.Abs2())) ||
      (IsZero(q2[0]) && IsZero(q2.Abs2())) ||
      (IsZero(q3[0]) && IsZero(q3.Abs2()))) {
    cout<<"Error in  Kinematics_FF::MakeKinematics("<<kt<<","<<z<<","<<y<<"):"<<endl
	<<"Splitter : "<<split<<":"<<(*split)<<"Spectator :"<<(*spect);
  }

  split->SetMomentum(q1);
  spect->SetMomentum(q2);
  return new Parton(newfl,q3,pst::FS);
}


Parton * Kinematics_FI::MakeKinematics(Parton * split, ATOOLS::Flavour & newfl)
{
  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum();
  Vec3D n1 = Vec3D(p1), n2 = Vec3D(p2), n3 = cross(n1,n2);
  if (n3.Sqr()<1.e-12) {
    double phi = 2.*M_PI*ran.Get();
    double cost = 1.-2.*ran.Get(), sint = sqrt(1.-cost*cost);
    Vec3D axis = sint*(cos(phi)*s_ex+sin(phi)*s_ey)+cost*s_ez;
    n3 = cross(n1,axis);
    if (n3.Sqr()<1.e-12) {
      n3 = cross(n2,axis);
      if (n3.Sqr()<1.e-12) {
	Vec3D axis = sint*(cos(phi)*s_ex-sin(phi)*s_ey)+cost*s_ez;
	n3 = cross(n1,axis);
	if (n3.Sqr()<1.e-12) {
	  msg.Error()<<"Error in Kinematics_FI::MakeKinematics : "<<endl
		     <<"   Could not construct an axis for the splitting : "<<endl
		     <<"   "<<n1<<" "<<n2<<"."<<endl;
	  abort();
	}
      }
    }
  }
  Vec4D n_perp = Vec4D(0.,n3/n3.Abs());
  
  double kt = sqrt(split->KtTest()), z = split->ZTest(), y = split->YTest();
  
  std::cout<<" New momenta for FI splitting : "<<z<<" "<<y<<std::endl;
  
  //sign of p2 terms changed!!!!
  Vec4D 
    q1 = z*p1      - (1.-z)*y/(1.-y)*p2 + kt*n_perp, 
    q2 =           +       1./(1.-y)*p2,
    q3 = (1.-z)*p1 -      z*y/(1.-y)*p2 - kt*n_perp;

  std::cout<<q1<<std::endl;
  std::cout<<q2<<std::endl;
  std::cout<<q3<<std::endl;
  
  if ((IsZero(q1[0]) && IsZero(q1.Abs2())) ||
      (IsZero(q2[0]) && IsZero(q2.Abs2())) ||
      (IsZero(q3[0]) && IsZero(q3.Abs2()))) {
    cout<<"Error in  Kinematics_FI::MakeKinematics("<<kt<<","<<z<<","<<y<<"):"<<endl
	<<"Splitter : "<<split<<":"<<(*split)<<"Spectator :"<<(*spect);
  }
  
  split->SetMomentum(q1);
  spect->SetMomentum(q2);
  double xp = spect->Xbj()/(1.-y);
  spect->SetXbj(xp);
  
  return new Parton(newfl,q3,pst::FS);
}


Parton * Kinematics_IF::MakeKinematics(Parton * split, ATOOLS::Flavour & newfl)
{
  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum();
  Vec3D n1 = Vec3D(p1), n2 = Vec3D(p2), n3 = cross(n1,n2);
  if (n3.Sqr()<1.e-12) {
    double phi = 2.*M_PI*ran.Get();
    double cost = 1.-2.*ran.Get(), sint = sqrt(1.-cost*cost);
    Vec3D axis = sint*(cos(phi)*s_ex+sin(phi)*s_ey)+cost*s_ez;
    n3 = cross(n1,axis);
    if (n3.Sqr()<1.e-12) {
      n3 = cross(n2,axis);
      if (n3.Sqr()<1.e-12) {
	Vec3D axis = sint*(cos(phi)*s_ex-sin(phi)*s_ey)+cost*s_ez;
	n3 = cross(n1,axis);
	if (n3.Sqr()<1.e-12) {
	  msg.Error()<<"Error in Kinematics_IF::MakeKinematics : "<<endl
		     <<"   Could not construct an axis for the splitting : "<<endl
		     <<"   "<<n1<<" "<<n2<<"."<<endl;
	  abort();
	}
      }
    }
  }
  Vec4D n_perp = Vec4D(0.,n3/n3.Abs());

  double kt = sqrt(split->KtTest()), z = split->ZTest(), y = split->YTest();

  //signs changed !!!!
  Vec4D 
    q1 =             1./z*p1,
    q2 =      -y*(1.-z)/z*p1 + (1.-y)*p2 - kt*n_perp, 
    q3 = -(1.-y)*(1.-z)/z*p1 +      y*p2 + kt*n_perp;

  if ((IsZero(q1[0]) && IsZero(q1.Abs2())) ||
      (IsZero(q2[0]) && IsZero(q2.Abs2())) ||
      (IsZero(q3[0]) && IsZero(q3.Abs2()))) {
    cout<<"Error in  Kinematics_IF::MakeKinematics("<<kt<<","<<z<<","<<y<<"):"<<endl
	<<"Splitter : "<<split<<":"<<(*split)<<"Spectator :"<<(*spect);
  }

  split->SetMomentum(q1);
  spect->SetMomentum(q2);
  double xp = split->Xbj()/z;
  split->SetXbj(xp);
  
  return new Parton(newfl,q3,pst::FS);
}
