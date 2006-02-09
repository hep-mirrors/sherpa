//bof
//Version: 4 ADICIC++-0.0/2006/02/03

//Implementation of the template structures of Sudakov_Group.H.



//#include ""





//using;    //is already done in Sudakov_Group.C





//=============================================================================



template<Dipole::Type DT>
Sudakov_Group<DT>::Sudakov_Group(const Radiation::Type ratyp)
  : Sudakov_Calculator(),
    m_radtype(ratyp),
    m_s(dpa.sud.MaxK2t()),
    m_x2tmin(dpa.sud.MinK2t()/m_s),
    m_x2tmax(dpa.sud.MaxK2t()/m_s),
    m_x2t(1.0), m_ymax(0.0), m_rap(0.0), m_corr(1.0),
    l_sud() {

  if(ratyp>Radiation::duscb) {
    Sudakov_Flavour sfc; sfc.Glu=&info.gluon.g;
    Sudakov_Base* gsud=new Sudakov<DT,Radiation::gluon>(*this, sfc);
    assert(gsud);
    l_sud.push_back(gsud);
  }

  Sudakov_Flavour sfc[6];
  sfc[1].Qua=&info.quark.d; sfc[1].Aqu=&info.antiq.d;
  sfc[2].Qua=&info.quark.u; sfc[2].Aqu=&info.antiq.u;
  sfc[3].Qua=&info.quark.s; sfc[3].Aqu=&info.antiq.s;
  sfc[4].Qua=&info.quark.c; sfc[4].Aqu=&info.antiq.c;
  sfc[5].Qua=&info.quark.b; sfc[5].Aqu=&info.antiq.b;
  size_t stop=ratyp%10;
  for(size_t i=1; i<=stop; ++i) {
    Sudakov_Base* qsud=new Sudakov<DT,Radiation::quark>(*this, sfc[i]);
    assert(qsud);
    l_sud.push_back(qsud);
  }

}


template<>
Sudakov_Group<Dipole::qqbar>::Sudakov_Group(const Radiation::Type ratyp)
  : Sudakov_Calculator(),
    m_radtype(ratyp),
    m_s(dpa.sud.MaxK2t()),
    m_x2tmin(dpa.sud.MinK2t()/m_s),
    m_x2tmax(dpa.sud.MaxK2t()/m_s),
    m_x2t(1.0), m_ymax(0.0), m_rap(0.0), m_corr(1.0),
    l_sud() {

  if(ratyp>Radiation::duscb) {
    Sudakov_Flavour sfc; sfc.Glu=&info.gluon.g;
    Sudakov_Base* gsud=new Sudakov<Dipole::qqbar,Radiation::gluon>(*this, sfc);
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



template<Dipole::Type DT>
const bool Sudakov_Group<DT>::GenerateVariablesFor(const Dipole& dip,
						   Sudakov_Result& sur) {

  sur.Reset();

  if(l_sud.empty()) { assert(p_dip==NULL && p_sur==NULL); return false;}

  assert(p_dip==NULL);
  p_dip=&dip;
  assert(p_sur==NULL);
  p_sur=&sur;

  if(InitWithCurrentDipole());
  else { p_dip=NULL; p_sur=NULL; return false;}

  xbool ksplt=between;
  Multidouble x(2);

  for(list<Sudakov_Base*>::const_iterator bit=l_sud.begin();
      bit!=l_sud.end(); ++bit) {

#ifdef SUDAKOV_CALCULATOR_OUTPUT
    cout<<"----------\n";
#endif

    Reset();
    Sudakov_Base& suda=**bit;

    while(true) {
      m_x2t=suda.GenerateX2t();
#ifdef SUDAKOV_CALCULATOR_OUTPUT
      cout<<" -> "<<m_x2t<<"\n";
#endif
      if(m_x2t<0.0) break;
      assert(m_s*m_x2t>dpa.sud.MinK2t());
      m_ymax=suda.CalculateRapLimit();
      m_rap=suda.GenerateRap();
      x[1]=sqrt(m_x2t);
      x[0]=1.0-x[1]*std::exp(m_rap);
      x[1]=1.0-x[1]*std::exp(-m_rap);
      m_corr=suda.GenerateCorr(x);
      if(ATOOLS::ran.Get() < m_corr) {
	if(TestEfracs(x[0],x[1])) {
	  if(p_sur->P2t > m_x2t) break;    //Keep the old values.
	  ksplt=suda.RadPart();
	  p_sur->Sfc=suda.RadCode();
	  p_sur->P2t=m_x2t;    //Temporarily keeping the new values.
	  p_sur->Y=m_rap;
	  p_sur->X1=x[0];
	  p_sur->X3=x[1];
	  break;
	}
	//break;//?????????????????????????????????????????????????????????????
      }
    }

  }

#ifdef SUDAKOV_CALCULATOR_OUTPUT
  cout<<"==========\n";
  cout<<" = "<<p_sur->P2t<<"\n";
  cout<<"=========="<<endl;
#endif

  p_sur->P2t*=m_s;    //Making p_sur->P2t a real p2t.

  if(p_sur->Sfc.Gluic() || p_sur->Sfc.Bqqic()) {
    switch(ksplt) {
    case between:
      p_sur->Rad=Radiation::gluon; assert(p_sur->Sfc.Gluic()); break;
    case front  :
      p_sur->Rad=Radiation::qtop; assert(p_sur->Sfc.Bqqic()); break;
    case end    :
      p_sur->Rad=Radiation::qbarbot; assert(p_sur->Sfc.Bqqic()); break;
    default     :
      assert(0);
    }
    p_dip=NULL;
    p_sur=NULL;
    return true;
  }
  else {
    p_dip=NULL;
    p_sur->Reset();
    p_sur=NULL;
    return false;
  }

}



//-----------------------------------------------------------------------------



template<Dipole::Type DT>
const bool Sudakov_Group<DT>::InitRadiation() const {
  cout<<__PRETTY_FUNCTION__<<endl;/////////////////////////////////////////////
  for(list<Sudakov_Base*>::const_iterator cit=l_sud.begin();
      cit!=l_sud.end(); ++cit)
    if(*cit) (*cit)->InitRadParticle();
  return true;
}





template<Dipole::Type DT>
const bool Sudakov_Group<DT>::InitWithCurrentDipole() {

  static bool check=InitRadiation();    //Warning removal trick.

  check=p_dip->IsType()==DT; assert(check);

  m_s=p_dip->InvMass();
  m_x2tmin=dpa.sud.MinK2t()/m_s;

  //m_x2tmax=Min(1.0,dpa.sud.MaxK2t()/m_s);
  //m_x2tmax=Min(p_dip->ProdScale()/m_s,dpa.sud.MaxK2t()/m_s);
  //m_x2tmax=Min(p_dip->BootScale()/m_s,dpa.sud.MaxK2t()/m_s);//Old.

  m_x2tmax=Min(0.25,dpa.sud.MaxK2t()/m_s);
  m_x2tmax=Min(m_x2tmax,p_dip->BootScale()/m_s);

  if(m_x2tmin<m_x2tmax) return true;
  //cout<<m_s<<" : "<<m_x2tmin<<" > "<<m_x2tmax<<" -> block evolution"<<endl;
  return false;

}



//=============================================================================



//===============
//Gluon emission.
//===============


template<Dipole::Type DT>
Sudakov<DT,Radiation::gluon>::Sudakov(const Sudakov_Group<DT>& sg,
				      const Sudakov_Flavour& go)
  : Sudakov_Base(go,Radiation::gluon), m_sgroup(sg) {

  assert(go.Gluic());
}





template<Dipole::Type DT>
Sudakov<DT,Radiation::gluon>::~Sudakov() {
#ifdef TEMP_OUTPUT
  std::cout<<"~Sudakov(gluon)"<<std::endl;/////////////////////////////////////
#endif
}





template<Dipole::Type DT>
const double Sudakov<DT,Radiation::gluon>::GenerateX2t() {
  m_step=0;
  double ran=ATOOLS::ran.Get();
#ifdef SUDAKOV_CALCULATOR_OUTPUT
  cout<<"\t\t\tg_random="<<ran<<"\n";
#endif
  double coeff=std::log(ran)*s_colfac/Sudakov_Calculator::AlphaSApprox();
  if(Sudakov_Calculator::Ariadne) coeff/=2.0;
  double A=sqr(std::log(m_sgroup.X2t()));
  if( coeff < A-sqr(std::log(m_sgroup.X2tmin())) ) return -1.0;
  return std::exp(-sqrt(A-coeff));
}





//==========================
//Quark-antiquark splitting.
//==========================


template<Dipole::Type DT>
Sudakov<DT,Radiation::quark>::Sudakov(const Sudakov_Group<DT>& sg,
				      const Sudakov_Flavour& qo)
  : Sudakov_Base(qo,Radiation::quark), m_sgroup(sg) {

  assert(qo.Bqqic());
  m_split=xbool(-1+2*(DT==Dipole::gqbar));
  assert(m_split!=between);////////////////////////////////////////////////////
}





template<Dipole::Type DT>
Sudakov<DT,Radiation::quark>::~Sudakov() {
#ifdef TEMP_OUTPUT
  std::cout<<"~Sudakov(quark)"<<std::endl;/////////////////////////////////////
#endif
}





template<Dipole::Type DT>
const double Sudakov<DT,Radiation::quark>::GenerateX2t() {
  m_step=0;
  double ran=ATOOLS::ran.Get();
#ifdef SUDAKOV_CALCULATOR_OUTPUT
  cout<<"\t\t\tq_random="<<ran<<"\n";
#endif
  double coeff=s_colfac/Sudakov_Calculator::AlphaSApprox();
  if(Sudakov_Calculator::Ariadne) coeff/=sqrt(m_sgroup.Sdip());
  double x2tcal=m_sgroup.X2t();
#ifdef SUDAKOV_CALCULATOR_OUTPUT
  cout<<"    "<<x2tcal<<"; "<<coeff<<"; "<<coeff*std::log(ran)<<", "
      <<std::log(m_sgroup.X2tmin()/x2tcal)<<"\n";
#endif
  if( coeff*std::log(ran) < std::log(m_sgroup.X2tmin()/x2tcal) )
    return -1.0;
  x2tcal*=std::pow(ran,coeff);
  if(NfVeto(m_sgroup.Sdip()*x2tcal)) return -1.0;
  //Subsequent x2tcal will be descending, and therefore not change the NfVeto()
  //result once it has been given true.
  SetSplitBranch();
  return x2tcal;
}



//=============================================================================





//eof
