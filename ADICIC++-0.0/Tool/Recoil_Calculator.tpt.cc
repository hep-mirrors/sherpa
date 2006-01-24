//bof
//Version: 3 ADICIC++-0.0/2005/08/22

//Implementation of template structures of Recoil_Calculator.H.



#include "Run_Parameter.H"





//using;





//=============================================================================



template<Recoil_Strategy::Type ST> Recoil<ST>::Recoil()
  : Recoil_Calculator(),
    m_costheta(-1.0), m_sintheta(0.0), m_phi(0.0), m_e() {}





template<Recoil_Strategy::Type ST> Recoil<ST>::~Recoil() {
#ifdef RECOIL_CALCULATOR_OUTPUT
  cout<<"  ~Recoil\n";
#endif
}



//-----------------------------------------------------------------------------



template<Recoil_Strategy::Type ST>
void Recoil<ST>::CrossProductTest(const Vec4D& axis) const {
  Vec3D q1(p_rer->Vec[rr::p1]);
  Vec3D q3(p_rer->Vec[rr::p3]);
  Vec3D ax(axis);
  Vec3D B=cross(q1,q3);
  Vec3D A=cross(ax,q1);
  cout<<"  Cross product test."<<endl;
  cout<<"  ax="<<ax<<endl;
  cout<<"  q1="<<q1<<endl;
  cout<<"  q3="<<q3<<endl;
  cout<<"   A="<<A<<endl;
  cout<<"   B="<<B<<endl;
  cout<<"     ";
  for(char i=1; i<4; ++i) cout<<A[i]/B[i]<<" : ";
  cout<<"     : ";
  for(char i=1; i<4; ++i) cout<<B[i]/A[i]<<" : ";
  cout<<endl;
  cout<<"  +++++++++++++++++++"<<endl;
}





template<Recoil_Strategy::Type ST>
void Recoil<ST>::RotateOnto(const ATOOLS::Vec4D& axis) {
  ATOOLS::Poincare rot(Recoil_Calculator::ZAxis,axis);
  rot.Rotate(p_rer->Vec[rr::p1]);
  rot.Rotate(p_rer->Vec[rr::p3]);
}


template<Recoil_Strategy::Type ST>
void Recoil<ST>::Rotate(const ATOOLS::Vec4D& axis,
			const ATOOLS::Vec4D& newaxis) {
  ATOOLS::Poincare rot(axis,newaxis);
  rot.Rotate(p_rer->Vec[rr::p1]);
  rot.Rotate(p_rer->Vec[rr::p3]);
}





template<Recoil_Strategy::Type ST>
const bool Recoil<ST>::CmsMomenta() {

  assert(p_sur->X1>0.0 && p_sur->X3>0.0);

  const Vec4D& Plab=p_dip->TotP();
  assert(Plab[0]>0.0);

#ifdef RECOIL_CALCULATOR_OUTPUT
  cout<<"\tlab frame - before:\n\t\t P =";
  cout<<Plab<<"\t"<<p_dip->InvMass()<<"  "<<p_dip->Mass()<<"\n";
  cout<<"\t\t q1="<<p_rer->Vec[rr::p1]<<"\t "<<p_rer->Vec[rr::p1].Abs2()<<"\n";
  cout<<"\t\t q3="<<p_rer->Vec[rr::p3]<<"\t "<<p_rer->Vec[rr::p3].Abs2()<<"\n";
#endif

  m_e.resize(4,0.0);
  m_e[2]=sqrt(p_dip->InvMass());
  m_e[1]=0.5*m_e[2]*p_sur->X1;
  m_e[3]=0.5*m_e[2]*p_sur->X3;
  m_e[2]=m_e[2]-m_e[1]-m_e[3];

  p_rer->Vec[rr::axis]=p_rer->Vec[rr::p1];    //Lab frame.
  Poincare fly(Plab);
  fly.Boost(p_rer->Vec[rr::axis]);    //This is the initial cms frame axis.

  m_costheta=( sqr(m_e[2])-sqr(m_e[1])-sqr(m_e[3]) )/( 2.0*m_e[1]*m_e[3] );
  if(dabs(m_costheta)>1.0) {
    cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": "
	<<"Error: `m_costheta' value is out of range!\n"<<endl;
    assert(dabs(m_costheta)<=1.0);
    p_rer->Reset();
    return false;
  }
  m_sintheta=sqrt(1.0-sqr(m_costheta));
  m_phi=2.0*M_PI*ATOOLS::ran.Get();

  if(Calculate()); else { p_rer->Reset(); return false;}

#ifdef RECOIL_CALCULATOR_OUTPUT
  cout<<"\tcms frame - after:\n";
  cout<<"\t\t p1="<<p_rer->Vec[rr::p1]<<"\t "<<p_rer->Vec[rr::p1].Abs2()<<endl;
  cout<<"\t\t p3="<<p_rer->Vec[rr::p3]<<"\t "<<p_rer->Vec[rr::p3].Abs2()<<endl;
#endif

  if(TEMP::CPTEST) CrossProductTest(p_rer->Vec[rr::axis]);/////////////////////

  fly.BoostBack(p_rer->Vec[rr::p1]);
  fly.BoostBack(p_rer->Vec[rr::p3]);

  p_rer->Vec[rr::p2]=Plab+(-1.0)*(p_rer->Vec[rr::p1]+p_rer->Vec[rr::p3]);
  //And hence the axis information is dead.

#ifdef RECOIL_CALCULATOR_OUTPUT
  cout<<"\tlab frame - after:\n";
  cout<<"\t\t p1="<<p_rer->Vec[rr::p1]<<"\t "<<p_rer->Vec[rr::p1].Abs2()<<endl;
  cout<<"\t\t p2="<<p_rer->Vec[rr::p2]<<"\t "<<p_rer->Vec[rr::p2].Abs2()<<endl;
  cout<<"\t\t p3="<<p_rer->Vec[rr::p3]<<"\t "<<p_rer->Vec[rr::p3].Abs2()<<endl;
#endif

  return true;

}





template<Recoil_Strategy::Type ST>
const bool Recoil<ST>::LabMomenta() {

  assert(p_sur->Isr.size()==sr::stop);
  assert(p_sur->Isr[sr::kt]>0.0);
  assert(p_sur->Isr[sr::expy]>0.0);
  assert(p_sur->Isr[sr::xpfin]>0.0 && p_sur->Isr[sr::xmfin]>0.0);

  //cout<<p_sur->Dir<<endl;
  //cout<<p_rer->Vec[rr::p1]<<endl;
  //cout<<p_rer->Vec[rr::p3]<<endl;

  extern Run_Parameter ATOOLS::rpa;
  double sqrtS=rpa.gen.Ecms();

  //Check.
  //double u=
  //  -p_sur->Isr[sr::xmfin]*sqrtS*p_sur->Isr[sr::kt]*p_sur->Isr[sr::expy];
  //double t=
  //  -p_sur->Isr[sr::xpfin]*sqrtS*p_sur->Isr[sr::kt]/p_sur->Isr[sr::expy];
  //cout<<u<<endl<<t<<endl;

  //Annihilation: 3..quark,        2..gluon,        1..antiquark.
  //Compton     : 3..quark(gluon), 2..f(anti)quark, 1..gluon(antiquark).

  m_e.resize(4,0.0);
  if(p_sur->Dir) {
    m_e[3]=p_sur->Isr[sr::xpfin]*sqrtS/2.0;
    m_e[1]=p_sur->Isr[sr::xmfin]*sqrtS/2.0;
    p_rer->Vec[rr::p3]=Vec4D(m_e[3],0.0,0.0,m_e[3]);
    p_rer->Vec[rr::p1]=Vec4D(m_e[1],0.0,0.0,-m_e[1]);
  } else {
    m_e[3]=p_sur->Isr[sr::xmfin]*sqrtS/2.0;
    m_e[1]=p_sur->Isr[sr::xpfin]*sqrtS/2.0;
    p_rer->Vec[rr::p3]=Vec4D(m_e[3],0.0,0.0,-m_e[3]);
    p_rer->Vec[rr::p1]=Vec4D(m_e[1],0.0,0.0,m_e[1]);
  }

  double g=p_sur->Isr[sr::expy]/2.0;
  double h=0.25/g;
  m_phi=2.0*M_PI*ATOOLS::ran.Get();
  p_rer->Vec[rr::p2]=Vec4D(g+h,cos(m_phi),sin(m_phi),g-h);
  p_rer->Vec[rr::p2]*=p_sur->Isr[sr::kt];

  p_rer->Poc=both;

  return true;

}



//=============================================================================



template<Recoil_Strategy::Type ST>
const bool Recoil<ST>::Calculate() {
  cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": "
      <<"Warning: Recoil strategy has not been specified!\n"<<endl;
  return false;
}



//-----------------------------------------------------------------------------



template<>
const bool Recoil<Recoil_Strategy::Kleiss>::Calculate() {
  if(TEMP::CPTEST) cout<<"Kleiss strategy."<<endl;/////////////////////////////

  if( ATOOLS::ran.Get() < sqr(m_e[1])/(sqr(m_e[1])+sqr(m_e[3])) ) {
    p_rer->Poc=end;
    //p_rer->Vec[rr::axis]=p_rer->Vec[rr::axis];    //cms frame
    p_rer->Vec[rr::p1]=Vec4D(m_e[1],0.0,0.0,m_e[1]);    //z-axis frame
    p_rer->Vec[rr::p3]=Vec4D(m_e[3],
			     m_e[3]*m_sintheta*cos(m_phi),
			     m_e[3]*m_sintheta*sin(m_phi),
			     m_e[3]*m_costheta);
  } else {
    p_rer->Poc=front;
    p_rer->Vec[rr::axis][1]*=-1;
    p_rer->Vec[rr::axis][2]*=-1;
    p_rer->Vec[rr::axis][3]*=-1;    //cms frame
    p_rer->Vec[rr::p1]=Vec4D(m_e[1],
			     m_e[1]*m_sintheta*cos(m_phi),
			     m_e[1]*m_sintheta*sin(m_phi),
			     m_e[1]*m_costheta);
    p_rer->Vec[rr::p3]=Vec4D(m_e[3],0.0,0.0,m_e[3]);    //z-axis frame
  }

#ifdef RECOIL_CALCULATOR_OUTPUT
  cout<<"\tcms frame - before:\n\t\t ax=";
  cout<<p_rer->Vec[rr::axis]<<"\t "<<p_rer->Vec[rr::axis].Abs2()<<endl;
  cout<<"\tz-axis frame - before (recoiled direction="<<p_rer->Poc<<"):\n";
  cout<<"\t\t p1="<<p_rer->Vec[rr::p1]<<"\t "<<p_rer->Vec[rr::p1].Abs2()<<endl;
  cout<<"\t\t p3="<<p_rer->Vec[rr::p3]<<"\t "<<p_rer->Vec[rr::p3].Abs2()<<endl;
#endif

  RotateOnto(p_rer->Vec[rr::axis]);

  return true;

}





template<>
const bool Recoil<Recoil_Strategy::FixDir1>::Calculate() {
  if(TEMP::CPTEST) cout<<"1fix strategy."<<endl;///////////////////////////////

  p_rer->Poc=end;
  //p_rer->Vec[rr::axis]=p_rer->Vec[rr::axis];    //cms frame
  p_rer->Vec[rr::p1]=Vec4D(m_e[1],0.0,0.0,m_e[1]);    //z-axis frame
  p_rer->Vec[rr::p3]=Vec4D(m_e[3],
			   m_e[3]*m_sintheta*cos(m_phi),
			   m_e[3]*m_sintheta*sin(m_phi),
			   m_e[3]*m_costheta);

#ifdef RECOIL_CALCULATOR_OUTPUT
  cout<<"\tcms frame - before:\n\t\t ax=";
  cout<<p_rer->Vec[rr::axis]<<"\t "<<p_rer->Vec[rr::axis].Abs2()<<endl;
  cout<<"\tz-axis frame - before (recoiled direction="<<p_rer->Poc<<"):\n";
  cout<<"\t\t p1="<<p_rer->Vec[rr::p1]<<"\t "<<p_rer->Vec[rr::p1].Abs2()<<endl;
  cout<<"\t\t p3="<<p_rer->Vec[rr::p3]<<"\t "<<p_rer->Vec[rr::p3].Abs2()<<endl;
#endif

  RotateOnto(p_rer->Vec[rr::axis]);

  return true;

}





template<>
const bool Recoil<Recoil_Strategy::FixDir3>::Calculate() {
  if(TEMP::CPTEST) cout<<"3fix strategy."<<endl;///////////////////////////////

  p_rer->Poc=front;
  p_rer->Vec[rr::axis][1]*=-1;
  p_rer->Vec[rr::axis][2]*=-1;
  p_rer->Vec[rr::axis][3]*=-1;    //cms frame
  p_rer->Vec[rr::p1]=Vec4D(m_e[1],
			   m_e[1]*m_sintheta*cos(m_phi),
			   m_e[1]*m_sintheta*sin(m_phi),
			   m_e[1]*m_costheta);
  p_rer->Vec[rr::p3]=Vec4D(m_e[3],0.0,0.0,m_e[3]);    //z-axis frame

#ifdef RECOIL_CALCULATOR_OUTPUT
  cout<<"\tcms frame - before:\n\t\t ax=";
  cout<<p_rer->Vec[rr::axis]<<"\t "<<p_rer->Vec[rr::axis].Abs2()<<endl;
  cout<<"\tz-axis frame - before (recoiled direction="<<p_rer->Poc<<"):\n";
  cout<<"\t\t p1="<<p_rer->Vec[rr::p1]<<"\t "<<p_rer->Vec[rr::p1].Abs2()<<endl;
  cout<<"\t\t p3="<<p_rer->Vec[rr::p3]<<"\t "<<p_rer->Vec[rr::p3].Abs2()<<endl;
#endif

  RotateOnto(p_rer->Vec[rr::axis]);

  return true;

}





template<>
const bool Recoil<Recoil_Strategy::MinimizePt>::Calculate() {
  if(TEMP::CPTEST) cout<<"Minimize strategy."<<endl;///////////////////////////
  double psi;
  const double theta=acos(m_costheta);
  if(m_e[3]<m_e[1]) {
    p_rer->Poc=end;
    //p_rer->Vec[rr::axis]=p_rer->Vec[rr::axis];    //cms frame
    p_rer->Vec[rr::p1]=Vec4D(m_e[1],0.0,0.0,m_e[1]);    //z-axis frame
    p_rer->Vec[rr::p3]=Vec4D(m_e[3],
			     m_e[3]*m_sintheta,
			     0.0,
			     m_e[3]*m_costheta);
    psi=atan2(-sqr(m_e[3])*sin(2*theta),
	      sqr(m_e[1])+sqr(m_e[3])*cos(2*theta))/2.0;
  } else {
    p_rer->Poc=front;
    p_rer->Vec[rr::axis][1]*=-1;
    p_rer->Vec[rr::axis][2]*=-1;
    p_rer->Vec[rr::axis][3]*=-1;    //cms frame
    p_rer->Vec[rr::p1]=Vec4D(m_e[1],
			     m_e[1]*m_sintheta,
			     0.0,
			     m_e[1]*m_costheta);
    p_rer->Vec[rr::p3]=Vec4D(m_e[3],0.0,0.0,m_e[3]);    //z-axis frame
    psi=atan2(-sqr(m_e[1])*sin(2*theta),
	      sqr(m_e[3])+sqr(m_e[1])*cos(2*theta))/2.0;
  }

  //m_phi=0.0;/////////////////////////////////////////////////////////////////

  Vec4D newaxis(1.0,sin(psi)*cos(m_phi),sin(psi)*sin(m_phi),cos(psi));
  RotateOnto(newaxis);

#ifdef RECOIL_CALCULATOR_OUTPUT
  cout<<"\tcms frame - before:\n\t\t ax=";
  cout<<p_rer->Vec[rr::axis]<<"\t "<<p_rer->Vec[rr::axis].Abs2()<<endl;
  cout<<"\tz-axis frame - before (recoiled direction="<<p_rer->Poc<<"):\n";
  cout<<"\t\t p1="<<p_rer->Vec[rr::p1]<<"\t "<<p_rer->Vec[rr::p1].Abs2()<<endl;
  cout<<"\t\t p3="<<p_rer->Vec[rr::p3]<<"\t "<<p_rer->Vec[rr::p3].Abs2()<<endl;
#endif

  RotateOnto(p_rer->Vec[rr::axis]);

  return true;

}





template<>
const bool Recoil<Recoil_Strategy::Lonnblad>::Calculate() {
  if(TEMP::CPTEST) cout<<"Lonnblad strategy."<<endl;///////////////////////////

  p_rer->Poc=both;
  //p_rer->Vec[rr::axis]=p_rer->Vec[rr::axis];    //cms frame
  p_rer->Vec[rr::p1]=Vec4D(m_e[1],0.0,0.0,m_e[1]);    //z-axis frame
  p_rer->Vec[rr::p3]=Vec4D(m_e[3],
			   m_e[3]*m_sintheta,
			   0.0,
			   m_e[3]*m_costheta);
  double psi=(M_PI-acos(m_costheta))*sqr(m_e[3])/(sqr(m_e[1])+sqr(m_e[3]));

  //m_phi=0.0;/////////////////////////////////////////////////////////////////

  Vec4D newaxis(1.0,sin(psi)*cos(m_phi),sin(psi)*sin(m_phi),cos(psi));
  RotateOnto(newaxis);

#ifdef RECOIL_CALCULATOR_OUTPUT
  cout<<"\tcms frame - before:\n\t\t ax=";
  cout<<p_rer->Vec[rr::axis]<<"\t "<<p_rer->Vec[rr::axis].Abs2()<<endl;
  cout<<"\tz-axis frame - before (recoiled direction="<<p_rer->Poc<<"):\n";
  cout<<"\t\t p1="<<p_rer->Vec[rr::p1]<<"\t "<<p_rer->Vec[rr::p1].Abs2()<<endl;
  cout<<"\t\t p3="<<p_rer->Vec[rr::p3]<<"\t "<<p_rer->Vec[rr::p3].Abs2()<<endl;
#endif

  RotateOnto(p_rer->Vec[rr::axis]);

  return true;

}





template<>
const bool Recoil<Recoil_Strategy::OldAdicic>::Calculate() {
  if(TEMP::CPTEST) cout<<"OldAdicic strategy."<<endl;//////////////////////////
  double psi;
  const double theta=acos(m_costheta);
  if(m_e[3]<m_e[1]) {
    p_rer->Poc=end;
    //p_rer->Vec[rr::axis]=p_rer->Vec[rr::axis];    //cms frame
    p_rer->Vec[rr::p1]=Vec4D(m_e[1],0.0,0.0,m_e[1]);    //z-axis frame
    p_rer->Vec[rr::p3]=Vec4D(m_e[3],
			     m_e[3]*m_sintheta,
			     0.0,
			     m_e[3]*m_costheta);
    psi=atan2(-sqr(m_e[3])*sin(2*theta),
	      sqr(m_e[1])+sqr(m_e[3])*cos(2*theta))/2.0;
  } else {
    p_rer->Poc=front;
    p_rer->Vec[rr::axis][1]*=-1;
    p_rer->Vec[rr::axis][2]*=-1;
    p_rer->Vec[rr::axis][3]*=-1;    //cms frame
    p_rer->Vec[rr::p1]=Vec4D(m_e[1],
			     m_e[1]*m_sintheta,
			     0.0,
			     m_e[1]*m_costheta);
    p_rer->Vec[rr::p3]=Vec4D(m_e[3],0.0,0.0,m_e[3]);    //z-axis frame
    psi=atan2(-sqr(m_e[1])*sin(2*theta),
	      sqr(m_e[3])+sqr(m_e[1])*cos(2*theta))/2.0;
  }

  RotateOnto(p_rer->Vec[rr::axis]);

  //m_phi=0.0;/////////////////////////////////////////////////////////////////

  Vec4D newaxis(1.0,sin(psi)*cos(m_phi),sin(psi)*sin(m_phi),cos(psi));
  RotateOnto(newaxis);

  return true;

}





template<>
const bool Recoil<Recoil_Strategy::Test>::Calculate() {
  if(TEMP::CPTEST) cout<<"Test strategy."<<endl;///////////////////////////////
  double psi, eta, f, g;
  const double theta=acos(m_costheta);
  if(m_e[3]<m_e[1]) {
    p_rer->Poc=end;
    //p_rer->Vec[rr::axis]=p_rer->Vec[rr::axis];    //cms frame
    p_rer->Vec[rr::p1]=Vec4D(m_e[1],0.0,0.0,m_e[1]);    //z-axis frame
    p_rer->Vec[rr::p3]=Vec4D(m_e[3],
			     m_e[3]*m_sintheta,
			     0.0,
			     m_e[3]*m_costheta);
    psi=atan2(-sqr(m_e[3])*sin(2*theta),
	      sqr(m_e[1])+sqr(m_e[3])*cos(2*theta))/2.0;
  } else {
    p_rer->Poc=front;
    p_rer->Vec[rr::axis][1]*=-1;
    p_rer->Vec[rr::axis][2]*=-1;
    p_rer->Vec[rr::axis][3]*=-1;    //cms frame
    p_rer->Vec[rr::p1]=Vec4D(m_e[1],
			     m_e[1]*m_sintheta,
			     0.0,
			     m_e[1]*m_costheta);
    p_rer->Vec[rr::p3]=Vec4D(m_e[3],0.0,0.0,m_e[3]);    //z-axis frame
    psi=atan2(-sqr(m_e[1])*sin(2*theta),
	      sqr(m_e[3])+sqr(m_e[1])*cos(2*theta))/2.0;
  }

  Vec4D psiaxis(1.0,sin(psi),0.0,cos(psi));
  Rotate(Vec4D::ZVEC,psiaxis);

  Vec4D phiaxis(1.0,cos(m_phi),sin(m_phi),0.0);
  Rotate(Vec4D::XVEC,phiaxis);

  do {
    eta=asin(2.0*ATOOLS::ran.Get()-1.0);
    eta=0.0;//0.17;////////////////////////////////////////////////////////////
    f=0.098032909*sqr(eta-1.570796327)*sqr(eta+1.570796327);
    g=0.6*cos(eta);
  } while(ATOOLS::ran.Get()<f/g);
  Vec4D etaaxis(1.0,0.0,sin(eta),cos(eta));
  Rotate(Vec4D::ZVEC,etaaxis);

#ifdef RECOIL_CALCULATOR_OUTPUT
  cout<<"\tcms frame - before:\n\t\t ax=";
  cout<<p_rer->Vec[rr::axis]<<"\t "<<p_rer->Vec[rr::axis].Abs2()<<endl;
  cout<<"\tz-axis frame - before (recoiled direction="<<p_rer->Poc<<"):\n";
  cout<<"\t\t p1="<<p_rer->Vec[rr::p1]<<"\t "<<p_rer->Vec[rr::p1].Abs2()<<endl;
  cout<<"\t\t p3="<<p_rer->Vec[rr::p3]<<"\t "<<p_rer->Vec[rr::p3].Abs2()<<endl;
#endif

  RotateOnto(p_rer->Vec[rr::axis]);

  return true;

}



//=============================================================================





//eof
