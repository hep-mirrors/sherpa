//bof
//Version: 2 ADICIC++-0.0/2004/10/28

//Implementation of the template structures of Sudakov_Calculator.H.



//#include ""





//using;    //is already done in Sudakov_Calculator.C





//=============================================================================



template<Dipole::Type DT>
Sudakov_Group<DT>::Sudakov_Group(const Radiation::Type ratyp)
  : Sudakov_Calculator(),
    m_radtype(ratyp),
    m_s(Sudakov_Calculator::MaxOfK2t()),
    m_x2tmin(Sudakov_Calculator::MinOfK2t()/m_s),
    m_x2tmax(Sudakov_Calculator::MaxOfK2t()/m_s),
    m_x2t(1.0), m_ymax(0.0), m_rap(0.0), m_corr(1.0),
    l_sud() {

  if(ratyp>5) {
    Sudakov_Base* gsud=
      new Sudakov<DT,Radiation::gluon>(*this,kf::gluon);
    assert(gsud);
    l_sud.push_back(gsud);
  }

  kf::code kfc[6];
  kfc[1]=kf::d;
  kfc[2]=kf::u;
  kfc[3]=kf::s;
  kfc[4]=kf::c;
  kfc[5]=kf::b;
  size_t stop=ratyp%10;
  for(size_t i=1; i<=stop; ++i) {
    Sudakov_Base* qsud=
      new Sudakov<DT,Radiation::quark>(*this,kfc[i]);
    assert(qsud);
    l_sud.push_back(qsud);
  }

}


template<>
Sudakov_Group<Dipole::qqbar>::Sudakov_Group(const Radiation::Type ratyp)
  : Sudakov_Calculator(),
    m_radtype(ratyp),
    m_s(Sudakov_Calculator::MaxOfK2t()),
    m_x2tmin(Sudakov_Calculator::MinOfK2t()/m_s),
    m_x2tmax(Sudakov_Calculator::MaxOfK2t()/m_s),
    m_x2t(1.0), m_ymax(0.0), m_rap(0.0), m_corr(1.0),
    l_sud() {

  if(ratyp>5) {
    Sudakov_Base* gsud=
      new Sudakov<Dipole::qqbar,Radiation::gluon>(*this,kf::gluon);
    assert(gsud);
    l_sud.push_back(gsud);
  }
}





template<Dipole::Type DT>
Sudakov_Group<DT>::~Sudakov_Group() {
  for(list<Sudakov_Base*>::iterator bit=l_sud.begin(); bit!=l_sud.end(); ++bit)
    if(*bit) delete (*bit);
#ifdef TEMP_OUTPUT
  cout<<"~Sudakov_Group"<<endl;////////////////////////////////////////////////
#endif
}



//-----------------------------------------------------------------------------



template<Dipole::Type DT> const bool
Sudakov_Group<DT>::GenerateEfracsFor(const Dipole& dip,
				     Sudakov_Strategy::Factorization) {

  InitWith(dip);

  if(l_sud.empty()) return false;

  double x1, x3;

  for(list<Sudakov_Base*>::const_iterator bit=l_sud.begin();
      bit!=l_sud.end(); ++bit) {

#ifdef SUDAKOV_CALCULATOR_OUTPUT
    cout<<"----------\n";
#endif

    Reset();
    Sudakov_Base& suda=**bit;

    while(m_x2t > m_x2tmin) {
      m_x2t=suda.GenerateX2t();
#ifdef SUDAKOV_CALCULATOR_OUTPUT
      cout<<" -> "<<m_x2t<<"\n";
#endif
      if(m_x2t<0.0) break;
      m_ymax=suda.CalculateRapLimit();
      m_rap=suda.GenerateRap();
      x3=sqrt(m_x2t);
      x1=1.0-x3*std::exp(m_rap);
      x3=1.0-x3*std::exp(-m_rap);
      m_corr=suda.GenerateCorr(x1,x3);
      if(ATOOLS::ran.Get() < m_corr) {
	if(TestEfracs(x1,x3)) {
	  if(m_p2t > m_x2t) break;    //Keep the old values.
	  m_spl=suda.RadPart();
	  m_rad=suda.RadCode();
	  m_p2t=m_x2t;    //Temporarily keeping the new values.
	  m_x1=x1;
	  m_x3=x3;
	  break;
	}
	break;//???????????????????????????????????????????????????????????????
      }
    }

  }

#ifdef SUDAKOV_CALCULATOR_OUTPUT
  cout<<"==========\n";
  cout<<" = "<<m_p2t<<"\n";
  cout<<"=========="<<endl;
#endif

  m_p2t*=m_s;    //Making m_p2t a real p2t.

  return m_rad!=kf::none;

}





template<Dipole::Type DT> const bool
Sudakov_Group<DT>::GenerateEfracsFor(const Dipole& dip,
				     Sudakov_Strategy::Distribution) {

  InitWith(dip);

  return false;

}



//-----------------------------------------------------------------------------



template<Dipole::Type DT>
const kf::code Sudakov_Group<DT>::InitRadiation() const {
  std::cout<<"RADINIT"<<std::endl;/////////////////////////////////////////////
  for(list<Sudakov_Base*>::const_iterator cit=l_sud.begin();
      cit!=l_sud.end(); ++cit)
    if(*cit) (*cit)->InitRadParticle();
  return kf::none;
}





template<Dipole::Type DT> void Sudakov_Group<DT>::InitWith(const Dipole& dip) {

  static const kf::code code_none=InitRadiation();    //Warning removal trick.

  bool dipole_type_allowed=
    dip.IsType()==Sudakov_Info<DT,Radiation::gluon>::Dipoletype;
  assert(dipole_type_allowed);

  m_spl=between; m_rad=code_none; m_p2t=0.0; m_x1=1.0; m_x3=1.0;

  m_s=dip.InvMass();
  m_x2tmin=Sudakov_Calculator::MinOfK2t()/m_s;

  //m_x2tmax=Min(1.0,Sudakov_Calculator::MaxOfK2t()/m_s);
  //m_x2tmax=Min(dip.ProdScale()/m_s,Sudakov_Calculator::MaxOfK2t()/m_s);
  //m_x2tmax=Min(dip.BootScale()/m_s,Sudakov_Calculator::MaxOfK2t()/m_s);//Old.

  m_x2tmax=Min(0.25,Sudakov_Calculator::MaxOfK2t()/m_s);
  m_x2tmax=Min(m_x2tmax,dip.BootScale()/m_s);

}



//=============================================================================



//===============
//Gluon emission.
//===============


template<Dipole::Type DT>
Sudakov<DT,Radiation::gluon>::Sudakov(const Sudakov_Group<DT>& sg, kf::code go)
  : Sudakov_Base(go), m_step(-1), m_sgroup(sg) {
  assert(go==kf::gluon);
}





template<Dipole::Type DT>
const double Sudakov<DT,Radiation::gluon>::GenerateX2t() {
  m_step=0;
  double ran=ATOOLS::ran.Get();
#ifdef SUDAKOV_CALCULATOR_OUTPUT
  cout<<"\t\t\tg_random="<<ran<<"\n";
#endif
  double coeff=
    std::log(ran)*
    Sudakov_Info<DT,Radiation::gluon>::Colourfactor/
    Sudakov_Calculator::AlphaSApprox();
  if(Sudakov_Calculator::Ariadne) coeff/=2.0;
  double A=sqr(std::log(m_sgroup.X2t()));
  if( coeff < A-sqr(std::log(m_sgroup.X2tmin())) ) return -1.0;
  return std::exp(-sqrt(A-coeff));
}





//==========================
//Quark-antiquark splitting.
//==========================


template<Dipole::Type DT>
Sudakov<DT,Radiation::quark>::Sudakov(const Sudakov_Group<DT>& sg, kf::code qo)
  : Sudakov_Base(qo), m_step(-1), m_sgroup(sg) {
  assert(qo==kf::d || qo==kf::u || qo==kf::s || qo==kf::c || qo==kf::b);
  m_split=
    xbool(-1+2*(Sudakov_Info<DT,Radiation::quark>::Dipoletype==Dipole::gqbar));
  assert(m_split!=between);////////////////////////////////////////////////////
}





template<Dipole::Type DT>
const double Sudakov<DT,Radiation::quark>::GenerateX2t() {
  m_step=0;
  double ran=ATOOLS::ran.Get();
#ifdef SUDAKOV_CALCULATOR_OUTPUT
  cout<<"\t\t\tq_random="<<ran<<"\n";
#endif
  double coeff=
    Sudakov_Info<DT,Radiation::quark>::Colourfactor/
    Sudakov_Calculator::AlphaSApprox();
  if(Sudakov_Calculator::Ariadne) coeff/=sqrt(m_sgroup.Sdip());
  double x2tcal=m_sgroup.X2t();
#ifdef SUDAKOV_CALCULATOR_OUTPUT
  cout<<"    "<<x2tcal<<"; "<<coeff<<"; "<<coeff*std::log(ran)<<", "
      <<std::log(m_sgroup.X2tmin()/x2tcal)<<"\n";
#endif
  if( coeff*std::log(ran) < std::log(m_sgroup.X2tmin()/x2tcal) )
    return -1.0;
  x2tcal*=std::pow(ran,coeff);
  if( m_code > Sudakov_Calculator::Nf(m_sgroup.Sdip()*x2tcal) )
    return -1.0;
  SetSplitBranch();
  return x2tcal;
}



//=============================================================================





//eof
