//bof
//Version: 1 ADICIC++-0.0/2004/05/10

//Implementation of template structures of Recoil_Calculator.H.



//#include ""





//using;





//=============================================================================



template<class ST>
const bool Recoil<ST>::Calculate() {    //Kleiss prescription.
  if(TEMP::CPTEST) cout<<"Kleiss strategy."<<endl;/////////////////////////////
  const double& E1=p_ini->GetE1();
  const double& E3=p_ini->GetE3();

  if( ATOOLS::ran.Get() < sqr(E1)/(sqr(E1)+sqr(E3)) ) {
    f_recoil=Negative;
    //m_cmsaxis=m_cmsaxis;    //cms frame
    m_p1=Vec4D(E1,0.0,0.0,E1);    //z-axis frame
    m_p3=Vec4D(E3,E3*m_sintheta*cos(m_phi),E3*m_sintheta*sin(m_phi),
	       E3*m_costheta);
  } else {
    f_recoil=Positive;
    m_cmsaxis[1]*=-1; m_cmsaxis[2]*=-1; m_cmsaxis[3]*=-1;    //cms frame
    m_p1=Vec4D(E1,E1*m_sintheta*cos(m_phi),E1*m_sintheta*sin(m_phi),
	       E1*m_costheta);
    m_p3=Vec4D(E3,0.0,0.0,E3);    //z-axis frame
  }

#ifdef RECOIL_CALCULATOR_OUTPUT
  cout<<"\tcms frame - before:\n\t\t ax=";
  cout<<m_cmsaxis<<"\t "<<m_cmsaxis.Abs2()<<endl;
  cout<<"\tz-axis frame - before (recoiled direction="<<f_recoil<<"):\n";
  cout<<"\t\t p1="<<m_p1<<"\t "<<m_p1.Abs2()<<endl;
  cout<<"\t\t p3="<<m_p3<<"\t "<<m_p3.Abs2()<<endl;
#endif

  RotateOnto(m_cmsaxis);

  return true;

}



//=============================================================================



template<>
const bool Recoil<FixDir1_Strategy>::Calculate() {    //FixDir1 prescription.
  if(TEMP::CPTEST) cout<<"1fix strategy."<<endl;///////////////////////////////
  const double& E1=p_ini->GetE1();
  const double& E3=p_ini->GetE3();

  f_recoil=Negative;
  //m_cmsaxis=m_cmsaxis;    //cms frame
  m_p1=Vec4D(E1,0.0,0.0,E1);    //z-axis frame
  m_p3=Vec4D(E3,E3*m_sintheta*cos(m_phi),E3*m_sintheta*sin(m_phi),
	     E3*m_costheta);

#ifdef RECOIL_CALCULATOR_OUTPUT
  cout<<"\tcms frame - before:\n\t\t ax=";
  cout<<m_cmsaxis<<"\t "<<m_cmsaxis.Abs2()<<endl;
  cout<<"\tz-axis frame - before (recoiled direction="<<f_recoil<<"):\n";
  cout<<"\t\t p1="<<m_p1<<"\t "<<m_p1.Abs2()<<endl;
  cout<<"\t\t p3="<<m_p3<<"\t "<<m_p3.Abs2()<<endl;
#endif

  RotateOnto(m_cmsaxis);

  return true;

}





template<>
const bool Recoil<FixDir3_Strategy>::Calculate() {    //FixDir3 prescription.
  if(TEMP::CPTEST) cout<<"3fix strategy."<<endl;///////////////////////////////
  const double& E1=p_ini->GetE1();
  const double& E3=p_ini->GetE3();

  f_recoil=Positive;
  m_cmsaxis[1]*=-1; m_cmsaxis[2]*=-1; m_cmsaxis[3]*=-1;    //cms frame
  m_p1=Vec4D(E1,E1*m_sintheta*cos(m_phi),E1*m_sintheta*sin(m_phi),
	     E1*m_costheta);
  m_p3=Vec4D(E3,0.0,0.0,E3);    //z-axis frame

#ifdef RECOIL_CALCULATOR_OUTPUT
  cout<<"\tcms frame - before:\n\t\t ax=";
  cout<<m_cmsaxis<<"\t "<<m_cmsaxis.Abs2()<<endl;
  cout<<"\tz-axis frame - before (recoiled direction="<<f_recoil<<"):\n";
  cout<<"\t\t p1="<<m_p1<<"\t "<<m_p1.Abs2()<<endl;
  cout<<"\t\t p3="<<m_p3<<"\t "<<m_p3.Abs2()<<endl;
#endif

  RotateOnto(m_cmsaxis);

  return true;

}





template<>
const bool Recoil<MinimizePt_Strategy>::Calculate() {    //MinimizePt
                                                         //prescription.
  if(TEMP::CPTEST) cout<<"Minimize strategy."<<endl;///////////////////////////
  double psi;
  const double theta=acos(m_costheta);

  const double& E1=p_ini->GetE1();
  const double& E3=p_ini->GetE3();

  if(E3<E1) {
    f_recoil=Negative;
    //m_cmsaxis=m_cmsaxis;    //cms frame
    m_p1=Vec4D(E1,0.0,0.0,E1);    //z-axis frame
    m_p3=Vec4D(E3,E3*m_sintheta,0.0,E3*m_costheta);
    psi=atan2(-sqr(E3)*sin(2*theta),sqr(E1)+sqr(E3)*cos(2*theta))/2.0;
  } else {
    f_recoil=Positive;
    m_cmsaxis[1]*=-1; m_cmsaxis[2]*=-1; m_cmsaxis[3]*=-1;    //cms frame
    m_p1=Vec4D(E1,E1*m_sintheta,0.0,E1*m_costheta);
    m_p3=Vec4D(E3,0.0,0.0,E3);    //z-axis frame
    psi=atan2(-sqr(E1)*sin(2*theta),sqr(E3)+sqr(E1)*cos(2*theta))/2.0;
  }

  //m_phi=0.0;/////////////////////////////////////////////////////////////////

  Vec4D newaxis(1.0,sin(psi)*cos(m_phi),sin(psi)*sin(m_phi),cos(psi));
  RotateOnto(newaxis);

#ifdef RECOIL_CALCULATOR_OUTPUT
  cout<<"\tcms frame - before:\n\t\t ax=";
  cout<<m_cmsaxis<<"\t "<<m_cmsaxis.Abs2()<<endl;
  cout<<"\tz-axis frame - before (recoiled direction="<<f_recoil<<"):\n";
  cout<<"\t\t p1="<<m_p1<<"\t "<<m_p1.Abs2()<<endl;
  cout<<"\t\t p3="<<m_p3<<"\t "<<m_p3.Abs2()<<endl;
#endif

  RotateOnto(m_cmsaxis);

  return true;

}





template<>
const bool Recoil<Lonnblad_Strategy>::Calculate() {    //Lonnblad prescription.
  if(TEMP::CPTEST) cout<<"Lonnblad strategy."<<endl;///////////////////////////
  const double& E1=p_ini->GetE1();
  const double& E3=p_ini->GetE3();

  f_recoil=Nil;
  //m_cmsaxis=m_cmsaxis;    //cms frame
  m_p1=Vec4D(E1,0.0,0.0,E1);    //z-axis frame
  m_p3=Vec4D(E3,E3*m_sintheta,0.0,E3*m_costheta);
  double psi=(M_PI-acos(m_costheta))*sqr(E3)/(sqr(E1)+sqr(E3));

  //m_phi=0.0;/////////////////////////////////////////////////////////////////

  Vec4D newaxis(1.0,sin(psi)*cos(m_phi),sin(psi)*sin(m_phi),cos(psi));
  RotateOnto(newaxis);

#ifdef RECOIL_CALCULATOR_OUTPUT
  cout<<"\tcms frame - before:\n\t\t ax=";
  cout<<m_cmsaxis<<"\t "<<m_cmsaxis.Abs2()<<endl;
  cout<<"\tz-axis frame - before (recoiled direction="<<f_recoil<<"):\n";
  cout<<"\t\t p1="<<m_p1<<"\t "<<m_p1.Abs2()<<endl;
  cout<<"\t\t p3="<<m_p3<<"\t "<<m_p3.Abs2()<<endl;
#endif

  RotateOnto(m_cmsaxis);

  return true;

}





template<>
const bool Recoil<OldAdicic_Strategy>::Calculate() {    //OldAdicic
                                                        //prescription.
  if(TEMP::CPTEST) cout<<"OldAdicic strategy."<<endl;//////////////////////////
  double psi;
  const double theta=acos(m_costheta);

  const double& E1=p_ini->GetE1();
  const double& E3=p_ini->GetE3();

  if(E3<E1) {
    f_recoil=Negative;
    //m_cmsaxis=m_cmsaxis;    //cms frame
    m_p1=Vec4D(E1,0.0,0.0,E1);    //z-axis frame
    m_p3=Vec4D(E3,E3*m_sintheta,0.0,E3*m_costheta);
    psi=atan2(-sqr(E3)*sin(2*theta),sqr(E1)+sqr(E3)*cos(2*theta))/2.0;
  } else {
    f_recoil=Positive;
    m_cmsaxis[1]*=-1; m_cmsaxis[2]*=-1; m_cmsaxis[3]*=-1;    //cms frame
    m_p1=Vec4D(E1,E1*m_sintheta,0.0,E1*m_costheta);
    m_p3=Vec4D(E3,0.0,0.0,E3);    //z-axis frame
    psi=atan2(-sqr(E1)*sin(2*theta),sqr(E3)+sqr(E1)*cos(2*theta))/2.0;
  }

  RotateOnto(m_cmsaxis);

  //m_phi=0.0;/////////////////////////////////////////////////////////////////

  Vec4D newaxis(1.0,sin(psi)*cos(m_phi),sin(psi)*sin(m_phi),cos(psi));
  RotateOnto(newaxis);

  return true;

}





template<>
const bool Recoil<Test_Strategy>::Calculate() {    //Test prescription.
  if(TEMP::CPTEST) cout<<"Test strategy."<<endl;///////////////////////////////
  double psi, eta, f, g;
  const double theta=acos(m_costheta);

  const double& E1=p_ini->GetE1();
  const double& E3=p_ini->GetE3();

  if(E3<E1) {
    f_recoil=Negative;
    //m_cmsaxis=m_cmsaxis;    //cms frame
    m_p1=Vec4D(E1,0.0,0.0,E1);    //z-axis frame
    m_p3=Vec4D(E3,E3*m_sintheta,0.0,E3*m_costheta);
    psi=atan2(-sqr(E3)*sin(2*theta),sqr(E1)+sqr(E3)*cos(2*theta))/2.0;
  } else {
    f_recoil=Positive;
    m_cmsaxis[1]*=-1; m_cmsaxis[2]*=-1; m_cmsaxis[3]*=-1;    //cms frame
    m_p1=Vec4D(E1,E1*m_sintheta,0.0,E1*m_costheta);
    m_p3=Vec4D(E3,0.0,0.0,E3);    //z-axis frame
    psi=atan2(-sqr(E1)*sin(2*theta),sqr(E3)+sqr(E1)*cos(2*theta))/2.0;
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
  cout<<m_cmsaxis<<"\t "<<m_cmsaxis.Abs2()<<endl;
  cout<<"\tz-axis frame - before (recoiled direction="<<f_recoil<<"):\n";
  cout<<"\t\t p1="<<m_p1<<"\t "<<m_p1.Abs2()<<endl;
  cout<<"\t\t p3="<<m_p3<<"\t "<<m_p3.Abs2()<<endl;
#endif

  RotateOnto(m_cmsaxis);

  return true;

}



//=============================================================================





//eof
