#include "METOOLS/Main/Polarization_Tools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "assert.h"

#ifndef SQRT_05
#define SQRT_05 0.70710678118654757
#endif

#ifndef sqrttwo
#define sqrttwo 1.4142135623730950488
#endif

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;


Polarization_Vector::Polarization_Vector(Vec4D p, bool anti,
                                         bool out) : std::vector<Vec4C>()
{
  Init(p);
}


Polarization_Vector::Polarization_Vector(Vec4D p, double m2, bool anti,
                                         bool out) : std::vector<Vec4C>()
{
  if(!IsZero((p.Abs2()-m2)/(fabs(p[0])+fabs(m2)), 1e-6)) {
    PRINT_INFO("Undefined for p.Abs2()="<<p.Abs2()<<" and m2="<<m2);
  }
  Init(p);
}

void Polarization_Vector::Init(Vec4D p)
{
  m_k=Vec4D(1.0,SQRT_05,SQRT_05,0.0);
  m_kp=SpinorType(1,m_k);
  m_km=SpinorType(-1,m_k);

  push_back(EMP(p));
  push_back(EMM(p));
  push_back(EML(p));
  push_back(Vec4C(Complex(0.,0.), Complex(0.,0.), Complex(0.,0.),
                  Complex(0.,0.)));
}

Vec4C Polarization_Vector::VT(const SpinorType &a,const SpinorType &b)
{
  Vec4C e;
  e[0]=a.U1()*b.U1()+a.U2()*b.U2();
  e[SpinorType::R3()]=a.U1()*b.U1()-a.U2()*b.U2();
  e[SpinorType::R1()]=a.U1()*b.U2()+a.U2()*b.U1();
  e[SpinorType::R2()]=Complex(0.0,1.0)*(a.U1()*b.U2()-a.U2()*b.U1());
  return e;
}

Vec4C Polarization_Vector::EM(const Vec4D &p)
{
  SpinorType pp(1,p);
  Vec4C e(VT(pp,m_km));
  return e/(sqrttwo*std::conj(m_kp*pp));
}

Vec4C Polarization_Vector::EP(const Vec4D &p)
{
  SpinorType pm(-1,p);
  Vec4C e(VT(m_kp,pm));
  return e/(sqrttwo*std::conj(m_km*pm));
}

Vec4C Polarization_Vector::EMM(const Vec4D &p)
{
  return EM(p-p.Abs2()/(2.0*m_k*p)*m_k);
}

Vec4C Polarization_Vector::EMP(const Vec4D &p)
{
  return EP(p-p.Abs2()/(2.0*m_k*p)*m_k);
}

Vec4C Polarization_Vector::EML(const Vec4D &p)
{
  double a(p.Abs2()/(2.0*m_k*p));
  Vec4D b(p-a*m_k);
  SpinorType bm(-1,b), bp(1,b), am(-1,m_k), ap(1,m_k);
  Vec4C e(VT(bp,bm)-a*VT(ap,am));
  return e/(2.0*p.Mass());
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
  bool massless(IsZero(p.Mass(),1e-6));
  
  std::cout<<METHOD<<": Testing transversality..."<<std::endl;
  bool success=true;
  for(int s=0;s<(massless?2:3);s++) {
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
  for(int s=0;s<(massless?2:3);s++) {
    for(int sp=0;sp<(massless?2:3);sp++) {
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
      Complex ref;
      if (massless) {
        for(int s=0;s<2;s++) {
          d+=conj((*this)[s])[nu]*(*this)[s][mu];
        }
        Vec4D q(p[0], -p[1], -p[2], -p[3]);
        ref=(p[mu]*q[nu]+p[nu]*q[mu])/(p*q)-g(mu, nu);
      }
      else {
        for(int s=0;s<3;s++) {
          d+=conj((*this)[s])[nu]*(*this)[s][mu];
        }
        ref=p[mu]*p[nu]/p.Abs2()-g(mu, nu);
      }
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
  bool massless(IsZero(p.Mass(),1e-6));
  if (massless) {
    THROW(not_implemented, "Massless polarisation tensors not implemented.");
  }
  
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
      for(int rho=0;rho<4;rho++) {
        for(int sigma=0;sigma<4;sigma++) {
          Complex d(0.0,0.0);
          for(int s=0;s<5;s++) {
            d+=(*this)[s][mu][nu]*(*this)[s].Conjugate()[rho][sigma];
          }
          Complex ref=S(mu,rho,p)*S(nu,sigma,p)/2.0+
                      S(mu,sigma,p)*S(nu,rho,p)/2.0-
                      S(mu,nu,p)*S(rho,sigma,p)/3.0;
          if(!IsEqual(d,ref)) {
            success=false;
            msg_Out()<<"  mu="<<mu<<std::endl;
            msg_Out()<<"    nu="<<nu<<std::endl;
            msg_Out()<<"      rho="<<rho<<std::endl;
            msg_Out()<<"        sigma="<<sigma<<std::endl;
            msg_Out()<<"          d="<<d<<std::endl;
            msg_Out()<<"          r="<<ref<<std::endl;  
          }
        }
      }
    }
  }
  std::cout<<"... "<<(success?"passed":"failed")<<std::endl<<std::endl;
}
