//bof
//Version: 1 ADICIC++-0.0/2004/03/12

//Implementation of Dipole_Handler.H.



#include "Random.H"
#include "Poincare.H"
#include "Dipole_Handler.H"
#include "Dipole_Handler.dat.cc"





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





//#include "....exa.cc"





//=============================================================================



//ostream& ADICIC::operator<<(ostream& ost, const ADICIC::...&) {
//}



//=============================================================================



int Dipole_Handler::s_count=0;    //So far there is no static Dipole_Handler.
const int& Dipole_Handler::InStore=Dipole_Handler::s_count;



//=============================================================================



Dipole_Handler::Dipole_Handler()
  : p_dip(NULL), p_dix(NULL), p_ban(NULL), p_ati(NULL), p_glu(NULL) {
  ++s_count;
}





Dipole_Handler::Dipole_Handler(Dipole& dip)
  : p_dip(NULL), p_dix(NULL), p_ban(NULL), p_ati(NULL), p_glu(NULL) {

  ++s_count;

  if(dip|*this) {
    if(dip.IsHandledBy(*this)); else {
      cerr<<"\nBug: Wrong Dipole-Dipole_Handler connection emerged!\n";
      assert(dip.IsHandledBy(*this));
    }
  } else {
    cerr<<"\nMethod: ADICIC::Dipole_Handler::Dipole_Handler(ADICIC::Dipole&): "
	<<"Warning: Attaching Dipole failed!\n"<<endl;
  }

}





Dipole_Handler::~Dipole_Handler() {

  --s_count;

  if(!p_dip) return;
  if(p_dip->IsHandledBy(*this)==false) {
    cerr<<"\nBug: Wrong Dipole-Dipole_Handler connection emerged!\n";
    assert(p_dip->IsHandledBy(*this));
  }

  *p_dip|0;

  if(p_dix) delete p_dix;
  if(p_ban) delete p_ban;
  if(p_ati) delete p_ati;
  if(p_glu) delete p_glu;

}





const bool Dipole_Handler::InduceGluonEmission() {

  if(!p_dip) return false;

  if(p_dip->Status() && p_dip->PointerHandling()==0 &&
     p_dip->IsType()==Dipole::qqbar);
  else return false;

  if(p_dix) { delete p_dix; p_dix=NULL;}
  if(p_ban) { delete p_ban; p_ban=NULL;}
  if(p_ati) { delete p_ati; p_ati=NULL;}
  if(p_glu) { delete p_glu; p_glu=NULL;}

  //No testing of global parameters.

  assert(p_dip->InvMass() > s_k2tmin);

  //assert(GenerateEfracs());
  //assert(TestEfracs());

  if( !GenerateEfracs() || !TestEfracs() ) {
    //cerr<<"\n"
    //<<"Method: const bool ADICIC::Dipole_Handler::InduceGluonEmission(): "
    //<<" Warning: Could not generate energy fractions!\n"<<endl;
    return false;
  }

#ifdef DIPOLE_HANDLER_OUTPUT
  cout<<"\ttransverse momentum and energy fractions:\n\t\t p2t=";
  cout<<m_p2t<<endl;
  cout<<"\t\t x1="<<m_x1<<endl;
  cout<<"\t\t x3="<<m_x3<<endl;
#endif

  assert((p_dip->TotP())[0] > 0.0);

  m_p1=p_dip->GetTopBranchPointer()->Momentum();
  m_p3=p_dip->GetBotBranchPointer()->Momentum();
  assert(m_p1.Abs2() >= 0.0);
  assert(m_p3.Abs2() >= 0.0);
  assert(m_p1[0] > 0.0);
  assert(m_p3[0] > 0.0);

  assert(GenerateMomenta());
  assert(GenerateSplitting());

  return true;

}



//=============================================================================



const bool Dipole_Handler::GenerateEfracs() {

  const double s=p_dip->InvMass();
  const double x2tmin=s_k2tmin/s;

  double x2t=Min(1.0,s_k2tmax/s);

  while(x2t>x2tmin) {
    double ran=ATOOLS::ran.Get();
#ifdef DIPOLE_HANDLER_OUTPUT
    cout<<"\t\t\tran="<<ran<<endl;
#endif
    double coeff=std::log(ran)*/*1.5*/0.75*M_PI/s_alphasfix;
    double A=sqr(std::log(x2t));
    if( coeff < A-sqr(std::log(x2tmin)) ) return false;
    x2t=std::exp(-sqrt(A-coeff));

    double ymax=-0.5*std::log(x2t);
    double rap=ymax*(-1.0+2.0*ATOOLS::ran.Get());

    m_x3=sqrt(x2t);
    m_x1=1.0-m_x3*std::exp(rap);
    m_x3=1.0-m_x3*std::exp(-rap);

    if( ATOOLS::ran.Get() < 0.5*(sqr(m_x1)+sqr(m_x3)) ) {
      m_p2t=s*x2t; return true;
    }
  }

  return false;

}





const bool Dipole_Handler::TestEfracs() const {
  double sum=m_x1+m_x3;
  if(sum>1.0 && sum<2.0) return true;
  return false;
}





const bool Dipole_Handler::GenerateMomenta() {

  const Vec4D& Plab=p_dip->TotP();
  //Preliminary approach - already done in EmitGluon():
  //m_p1=p_dip->GetTopBranchPointer()->Momentum();
  //m_p3=p_dip->GetBotBranchPointer()->Momentum();

#ifdef DIPOLE_HANDLER_OUTPUT
  cout<<"\tlab frame - before:\n\t\t P =";
  cout<<Plab<<"\t"<<p_dip->InvMass()<<"  "<<p_dip->Mass()<<endl;
  cout<<"\t\t q1="<<m_p1<<"\t "<<m_p1.Abs2()<<endl;
  cout<<"\t\t q3="<<m_p3<<"\t "<<m_p3.Abs2()<<endl;
#endif

  Vec4D& axis=m_p2;

  double E2=sqrt(p_dip->InvMass());
  double E1=0.5*E2*m_x1;
  double E3=0.5*E2*m_x3;

  E2=E2-E1-E3;

  double costheta=(sqr(E2)-sqr(E1)-sqr(E3))/(2.0*E1*E3);
  if(dabs(costheta)>1.0) {
    cerr<<"\nError: `costheta' value is out of range!\n";
    assert(0);
  }
  double sintheta=sqrt(1.0-sqr(costheta));
  double phi=2.0*M_PI*ATOOLS::ran.Get();

  if( ATOOLS::ran.Get() < sqr(m_x1)/(sqr(m_x1)+sqr(m_x3)) ) {
    f_quarkrecoil=false;
    axis=m_p1;    //lab frame
    m_p1=Vec4D(E1,0.0,0.0,E1);    //z-axis frame
    m_p3=Vec4D(E3,E3*sintheta*cos(phi),E3*sintheta*sin(phi),E3*costheta);
  } else {
    f_quarkrecoil=true;
    axis=m_p3;    //lab frame
    m_p1=Vec4D(E1,E1*sintheta*cos(phi),E1*sintheta*sin(phi),E1*costheta);
    m_p3=Vec4D(E3,0.0,0.0,E3);    //z-axis frame
  }

  Poincare fly(Plab);
  fly.Boost(axis);

#ifdef DIPOLE_HANDLER_OUTPUT
  cout<<"\tcm frame - before:\n\t\t ax=";
  cout<<axis<<"\t "<<axis.Abs2()<<endl;

  cout<<"\tz-axis frame - before (quarkrecoil="<<f_quarkrecoil<<"):\n";
  cout<<"\t\t p1="<<m_p1<<"\t "<<m_p1.Abs2()<<endl;
  cout<<"\t\t p3="<<m_p3<<"\t "<<m_p3.Abs2()<<endl;
#endif

  Poincare rot(Vec4D(1.0,0.0,0.0,1.0),axis);
  rot.Rotate(m_p1);
  rot.Rotate(m_p3);

#ifdef DIPOLE_HANDLER_OUTPUT
  cout<<"\tcm frame - after:\n";
  cout<<"\t\t p1="<<m_p1<<"\t "<<m_p1.Abs2()<<endl;
  cout<<"\t\t p3="<<m_p3<<"\t "<<m_p3.Abs2()<<endl;
#endif

  fly.BoostBack(m_p1);
  fly.BoostBack(m_p3);

  m_p2=Plab+(-1.0)*(m_p1+m_p3);

#ifdef DIPOLE_HANDLER_OUTPUT
  cout<<"\tlab frame - after:\n";
  cout<<"\t\t p1="<<m_p1<<"\t "<<m_p1.Abs2()<<endl;
  cout<<"\t\t p2="<<m_p2<<"\t "<<m_p2.Abs2()<<endl;
  cout<<"\t\t p3="<<m_p3<<"\t "<<m_p3.Abs2()<<endl;
#endif

  return true;

}





const bool Dipole_Handler::GenerateSplitting() {

  p_dip->GetTopBranchPointer()->SetMomentum(m_p1);
  p_dip->GetBotBranchPointer()->SetMomentum(m_p3);

  p_glu=new Dipole::Glubranch(m_p2); assert(p_glu);
  p_dix=new Dipole(*p_dip); assert(p_dix);

  if(f_quarkrecoil) {
    p_dix->RenewBranch(false,*p_glu);
    p_dip->RenewBranch(true,*p_glu);
  } else {
    p_dip->RenewBranch(false,*p_glu);
    p_dix->RenewBranch(true,*p_glu);
  }

  p_dix->SetSource()=p_dip->Name;
  p_dip->SetProdScale()=m_p2t;
  p_dix->SetProdScale()=m_p2t;

  return true;

}



//=============================================================================





//eof
