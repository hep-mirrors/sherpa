#include "Polarization_Tools.H"
#include "Message.H"

#ifndef SQRT_05
#define SQRT_05 0.70710678118654757
#endif

using namespace HELICITIES;
using namespace ATOOLS;
using namespace std;

Polarization_Vector::Polarization_Vector(Vec4D p, double m2, bool anti,
                                         bool out) : std::vector<Vec4C>()
{
  double pAbs2(p.Abs2()), pPSpat(p.PSpat());
  if (pAbs2 < -sqrt(Accu()) || pAbs2 < m2-sqrt(Accu())) {
    msg_Error()<<METHOD<<": p^2 = "<<pAbs2<<" < 0 or < m2 ="<<m2<<". "
               <<"Please check definition of complex boson polarization vector "
               <<"for correctness (signs) and remove this warning."<<endl;
  }

  if (pPSpat < Accu()) {      // decay from rest
    push_back(Vec4C(Complex(0.,0.), Complex(SQRT_05,0.), Complex(0.,-SQRT_05),
                    Complex(0.,0.)));
    push_back(Vec4C(Complex(0.,0.), Complex(SQRT_05,0.), Complex(0.,SQRT_05),
                    Complex(0.,0.)));
    push_back(Vec4C(Complex(0.,0.), Complex(0.,0.), Complex(0.,0.),
                    Complex(1.,0.)));
    push_back(Vec4C(Complex(0.,0.), Complex(0.,0.), Complex(0.,0.),
                    Complex(0.,0.)));
  }
  else {
    double ct(p.CosTheta()), st(p.SinTheta()), cp(p.CosPhi()), sp(p.SinPhi());
    
    push_back(Vec4C(Complex(0.,0.),SQRT_05*Complex(ct*cp,-sp),
                    SQRT_05*Complex(ct*sp,cp),SQRT_05*Complex(-st,0.0)));
    push_back(Vec4C(Complex(0.,0.),SQRT_05*Complex(ct*cp,sp),
                    SQRT_05*Complex(ct*sp,-cp),SQRT_05*Complex(-st,0.0)));
    push_back(1.0/p.Mass()*Vec4C(Complex(pPSpat,0.0),
                                 (p[0]/pPSpat)*Vec3C(p[1],p[2],p[3]))); // 0
    Vec4D real= (abs(pAbs2-m2)<Accu()) ? 
      Vec4D(0.,0.,0.,0.):sqrt((pAbs2-m2)/(pAbs2*m2))*p; // s
    push_back(Vec4C(Complex(real[0],0.0),Complex(real[1],0.0),
                    Complex(real[2],0.0),Complex(real[3],0.0)));
  }
}


Polarization_Vector::Polarization_Vector(Vec4D p, bool anti,
                                         bool out) : std::vector<Vec4C>()
{
  double pPSpat(p.PSpat());
  if (pPSpat < Accu()) {      // decay from rest
    push_back(Vec4C(Complex(0.,0.), Complex(SQRT_05,0.), Complex(0.,-SQRT_05),
                    Complex(0.,0.)));
    push_back(Vec4C(Complex(0.,0.), Complex(SQRT_05,0.), Complex(0.,SQRT_05),
                    Complex(0.,0.)));
    push_back(Vec4C(Complex(0.,0.), Complex(0.,0.), Complex(0.,0.),
                    Complex(1.,0.)));
    push_back(Vec4C(Complex(0.,0.), Complex(0.,0.), Complex(0.,0.),
                    Complex(0.,0.)));
  }
  else {
    double ct(p.CosTheta()), st(p.SinTheta()), cp(p.CosPhi()), sp(p.SinPhi());
    
    push_back(Vec4C(Complex(0.,0.),SQRT_05*Complex(ct*cp,-sp),
                    SQRT_05*Complex(ct*sp,cp),SQRT_05*Complex(-st,0.0)));
    push_back(Vec4C(Complex(0.,0.),SQRT_05*Complex(ct*cp,sp),
                    SQRT_05*Complex(ct*sp,-cp),SQRT_05*Complex(-st,0.0)));
    push_back(1.0/p.Mass()*Vec4C(Complex(pPSpat,0.0),
                                 (p[0]/pPSpat)*Vec3C(p[1],p[2],p[3]))); // 0
    push_back(Vec4C(Complex(0.,0.), Complex(0.,0.), Complex(0.,0.),
                    Complex(0.,0.)));
  }
}


Polarization_Tensor::Polarization_Tensor(Vec4D p, double m2, bool anti,
                                         bool out) : std::vector<CMatrix>()
{
  Polarization_Vector eps(p,m2);
  
  CMatrix tensor(4);
  for(int mu = 0; mu < 4; mu++){
    for(int nu = 0; nu < 4; nu++){
      Complex eps_munu = (eps[0][mu]*eps[0][nu]);
      tensor[mu][nu] = eps_munu;
    }
  }
  push_back(tensor);

  for(int mu = 0; mu < 4; mu++){
    for(int nu = 0; nu < 4; nu++){
      Complex eps_munu = (eps[0][mu]*eps[2][nu]+eps[2][mu]*
                          eps[0][nu])/sqrt(2.0);
      tensor[mu][nu] = eps_munu;
    }
  }
  push_back(tensor);

  for(int mu = 0; mu < 4; mu++){
    for(int nu = 0; nu < 4; nu++){
      Complex eps_munu = (eps[0][mu]*eps[1][nu]+eps[1][mu]*
                          eps[0][nu]
                          -2.0*eps[2][mu]*eps[2][nu])/sqrt(6.0);
      tensor[mu][nu] = eps_munu;
    }
  }
  push_back(tensor);

  for(int mu = 0; mu < 4; mu++){
    for(int nu = 0; nu < 4; nu++){
      Complex eps_munu = (eps[1][mu]*eps[2][nu]+eps[2][mu]*
                          eps[1][nu])/sqrt(2.0);
      tensor[mu][nu] = eps_munu;
    }
  }
  push_back(tensor);

  for(int mu = 0; mu < 4; mu++){
    for(int nu = 0; nu < 4; nu++){
      Complex eps_munu = (eps[1][mu]*eps[1][nu]);
      tensor[mu][nu] = eps_munu;
    }
  }
  push_back(tensor);
}


double g(int mu, int nu)
{
  if (mu<0 || mu>3 || nu<0 || nu>3)
    std::cout<<"wrong indices in g(mu, nu)."<<std::endl;
  if (mu!=nu) return 0.0;
  if (mu==0) return 1.0;
  if (mu>0 && mu<4) return -1.0;
  return 0.0;
}

double S(int mu, int nu, Vec4D p)
{
  return p[mu]*p[nu]/p.Abs2()-g(mu, nu);
}

void Polarization_Vector::Test(Vec4D p)
{
  std::cout<<METHOD<<": Testing transversality..."<<std::endl;
  bool success=true;
  for(int s=0;s<3;s++) {
    Complex d = (*this)[s]*p;
    if(!IsZero(d)) {
      success=false;
      msg_Out()<<"s="<<s<<std::endl;
      msg_Out()<<"  d="<<d<<" should be: 0.0"<<std::endl;
    }
  }
  std::cout<<METHOD<<"... "<<(success?"passed":"failed")<<std::endl<<std::endl;

  std::cout<<METHOD<<": Testing orthogonality..."<<std::endl;
  success=true;
  for(int s=0;s<3;s++) {
    for(int sp=0;sp<3;sp++) {
      Complex d = conj((*this)[s])*(*this)[sp];
      if( (s==sp && !IsEqual(d,Complex(-1.0,0.0))) ||
          (s!=sp && !IsZero(d))) {
        success=false;
        msg_Out()<<"s="<<s<<" sp="<<sp<<std::endl;
        msg_Out()<<"  d="<<d<<" should be: "<<(s==sp?-1.0:0.0)<<std::endl;
      }
    }
  }
  std::cout<<METHOD<<"... "<<(success?"passed":"failed")<<std::endl<<std::endl;

  std::cout<<METHOD<<": Testing completeness..."<<std::endl;
  success=true;
  for(int mu=0;mu<4;mu++) {
    for(int nu=0;nu<4;nu++) {
      Complex d(0.0,0.0);
      for(int s=0;s<3;s++) {
        d+=conj((*this)[s])[nu]*(*this)[s][mu];
      }
      Complex ref=p[mu]*p[nu]/p.Abs2()-g(mu, nu);
      if(!IsEqual(d,ref)) {
        success=false;
        msg_Out()<<"  mu="<<mu<<std::endl;
        msg_Out()<<"    nu="<<nu<<std::endl;
        msg_Out()<<"      d="<<d<<std::endl;
        msg_Out()<<"      r="<<ref<<std::endl;  
      }
    }
  }
  std::cout<<METHOD<<"... "<<(success?"passed":"failed")<<std::endl<<std::endl;
}

void Polarization_Tensor::Test(Vec4D p)
{
  bool success = true;
  std::cout<<METHOD<<": Testing symmetry..."<<std::endl;
  for(int s=0;s<5;s++) {
    for(int mu=0;mu<4;mu++) {
      for(int nu=mu+1;nu<4;nu++) {
        if(!IsEqual((*this)[s][mu][nu],(*this)[s][nu][mu])) {
        success=false;
        msg_Out()<<"s="<<s<<std::endl;
        msg_Out()<<"  eps["<<mu<<"]["<<nu<<"]="<<(*this)[s][mu][nu]<<std::endl
                 <<"  eps["<<nu<<"]["<<mu<<"]="<<(*this)[s][nu][mu]<<std::endl;
        }
      }
    }
  }
  std::cout<<"... "<<(success?"passed":"failed")<<std::endl<<std::endl;

  success=true;
  std::cout<<METHOD<<": Testing transversality..."<<std::endl;
  for(int s=0;s<5;s++) {
    Vec4C test1;
    for(int nu=0;nu<4;nu++) {
      for(int mu=0;mu<4;mu++) {
        double pmu = mu==0 ? p[mu] : -p[mu];
        test1[nu] += (*this)[s][mu][nu]*pmu;
      }
      if(!IsZero(test1[nu])) {
        success=false;
        msg_Out()<<"s="<<s<<std::endl;
        msg_Out()<<"  test1["<<nu<<"]="<<test1[nu]<<std::endl;
      }
    }
  }
  std::cout<<"... "<<(success?"passed":"failed")<<std::endl<<std::endl;

  success=true;
  std::cout<<METHOD<<": Testing tracelessness..."<<std::endl;
  for(int s=0;s<5;s++) {
    Complex d(0.0,0.0);
    for(int mu=0;mu<4;mu++) {
      double facmu = mu==0 ? 1.0 : -1.0;
      // should the metric go here, or is it really a trace?
      d += facmu*(*this)[s][mu][mu];
    }
    if(!IsZero(d)) {
      success=false;
      msg_Out()<<"s="<<s<<std::endl;
      msg_Out()<<"  d="<<d<<std::endl;
    }
  }
  std::cout<<"... "<<(success?"passed":"failed")<<std::endl<<std::endl;

  success=true;
  std::cout<<METHOD<<": Testing orthonormality..."<<std::endl;
  for(int s=0;s<5;s++) {
    for(int sp=0;sp<5;sp++) {
      Complex d(0.0,0.0);
      for(int nu=0;nu<4;nu++) {
        for(int mu=0;mu<4;mu++) {
          double facmu = mu==0 ? 1.0 : -1.0;
          double facnu = nu==0 ? 1.0 : -1.0;
          d += facmu*facnu*(*this)[s][mu][nu]*(*this)[sp].Conjugate()[mu][nu];
        }
      }
      if( (s==sp && !IsEqual(d,Complex(1.0,0.0))) ||
          (s!=sp && !IsZero(d))) {
        success=false;
        msg_Out()<<"s="<<s<<" sp="<<sp<<std::endl;
        msg_Out()<<"  d="<<d<<" should be: "<<(s==sp?1.0:0.0)<<std::endl;
      }
    }
  }
  std::cout<<"... "<<(success?"passed":"failed")<<std::endl<<std::endl;

  success=true;
  std::cout<<METHOD<<": Testing completeness..."<<std::endl;
  for(int mu=0;mu<4;mu++) {
    for(int nu=0;nu<4;nu++) {
      for(int a=0;a<4;a++) {
        for(int b=0;b<4;b++) {
          Complex d(0.0,0.0);
          for(int s=0;s<5;s++) {
            d+=(*this)[s][mu][nu]*(*this)[s].Conjugate()[a][b];
          }
          Complex ref=0.5*(S(mu,a,p)*S(nu,b,p)+S(mu,b,p)*S(nu,a,p))-
            1.0/3.0*S(mu,nu,p)*S(a,b,p);
          if(!IsEqual(d,ref)) {
            success=false;
            msg_Out()<<"  mu="<<mu<<std::endl;
            msg_Out()<<"    nu="<<nu<<std::endl;
            msg_Out()<<"      a="<<a<<std::endl;
            msg_Out()<<"        b="<<b<<std::endl;
            msg_Out()<<"          d="<<d<<std::endl;
            msg_Out()<<"          r="<<ref<<std::endl;  
          }
        }
      }
    }
  }
  std::cout<<"... "<<(success?"passed":"failed")<<std::endl<<std::endl;
}
