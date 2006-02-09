//bof
//Version: 4 ADICIC++-0.0/2006/02/03

//Implementation of the template structures of IISudakov_Group.H.



#include "Run_Parameter.H"





//using;    //is already done in IISudakov_Group.C





//=============================================================================



template<Dipole::Type DT>
IISudakov_Group<DT>::IISudakov_Group(const Radiation::Type ratyp)
  : Sudakov_Calculator(),
    m_radtype(ratyp),
    m_sdip(dpa.sud.MaxIIK2t()),
    m_a(2*sqrt(2.0)+3.0),
    m_s(sqr(m_a-1.0)*dpa.sud.MaxIIK2t()),
    m_x2tmin(dpa.sud.MinIIK2t()/m_s),
    m_x2tmax(0.25),
    m_x2t(1.0), m_ymax(0.0), m_rap(0.0), m_corr(1.0),
    l_sud() {

  if(ratyp>Radiation::duscb) {
    Sudakov_Flavour sfc; sfc.Glu=&info.gluon.g;
    Sudakov_Base* gsud=new IISudakov<DT,Radiation::gluon>(*this,sfc);
    assert(gsud);
    l_sud.push_back(gsud);
  }

  if(ratyp!=Radiation::g) {
    Sudakov_Flavour sfc;
    Sudakov_Base* qsud;
    if(DT==Dipole::iiqbarg) {
      sfc.Qua=&info.quark.d;
      qsud=new IISudakov<DT,Radiation::igluon>(*this,sfc);
      assert(qsud);
      l_sud.push_back(qsud);
    }
    if(DT==Dipole::iigq) {
      sfc.Aqu=&info.antiq.d;
      qsud=new IISudakov<DT,Radiation::igluon>(*this,sfc);
      assert(qsud);
      l_sud.push_back(qsud);
    }
  }

  Sudakov_Flavour sfc[6];
  sfc[1].Qua=&info.quark.d; sfc[1].Aqu=&info.antiq.d;
  sfc[2].Qua=&info.quark.u; sfc[2].Aqu=&info.antiq.u;
  sfc[3].Qua=&info.quark.s; sfc[3].Aqu=&info.antiq.s;
  sfc[4].Qua=&info.quark.c; sfc[4].Aqu=&info.antiq.c;
  sfc[5].Qua=&info.quark.b; sfc[5].Aqu=&info.antiq.b;
  size_t stop=ratyp%10;
  for(size_t i=1; i<=stop; ++i) {
    Sudakov_Base* qsud=new IISudakov<DT,Radiation::quark>(*this,sfc[i]);
    assert(qsud);
    l_sud.push_back(qsud);
  }

}



template<>
IISudakov_Group<Dipole::iiqbarq>::IISudakov_Group(const Radiation::Type ratyp)
  : Sudakov_Calculator(),
    m_radtype(ratyp),
    m_sdip(dpa.sud.MaxIIK2t()),
    m_a(2*sqrt(2.0)+3.0),
    m_s(sqr(m_a-1.0)*dpa.sud.MaxIIK2t()),
    m_x2tmin(dpa.sud.MinIIK2t()/m_s),
    m_x2tmax(0.25),
    m_x2t(1.0), m_ymax(0.0), m_rap(0.0), m_corr(1.0),
    l_sud() {

  if(ratyp>Radiation::duscb) {
    Sudakov_Flavour sfc; sfc.Glu=&info.gluon.g;
    Sudakov_Base* gsud=
      new IISudakov<Dipole::iiqbarq,Radiation::gluon>(*this,sfc);
    assert(gsud);
    l_sud.push_back(gsud);
  }

  if(ratyp!=Radiation::g) {
    Sudakov_Flavour sfc;
    Sudakov_Base* qsud;
    sfc.Qua=&info.quark.d;
    qsud=new IISudakov<Dipole::iiqbarq,Radiation::igluon>(*this,sfc);
    assert(qsud);
    l_sud.push_back(qsud);
    sfc.Qua=NULL;
    sfc.Aqu=&info.antiq.d;
    qsud=new IISudakov<Dipole::iiqbarq,Radiation::igluon>(*this,sfc);
    assert(qsud);
    l_sud.push_back(qsud);
  }

}





template<>
IISudakov_Group<Dipole::iigg>::IISudakov_Group(const Radiation::Type ratyp)
  : Sudakov_Calculator(),
    m_radtype(ratyp),
    m_sdip(dpa.sud.MaxIIK2t()),
    m_a(2*sqrt(2.0)+3.0),
    m_s(sqr(m_a-1.0)*dpa.sud.MaxIIK2t()),
    m_x2tmin(dpa.sud.MinIIK2t()/m_s),
    m_x2tmax(0.25),
    m_x2t(1.0), m_ymax(0.0), m_rap(0.0), m_corr(1.0),
    l_sud() {

  if(ratyp>Radiation::duscb) {
    Sudakov_Flavour sfc; sfc.Glu=&info.gluon.g;
    Sudakov_Base* gsud=new IISudakov<Dipole::iigg,Radiation::gluon>(*this,sfc);
    assert(gsud);
    l_sud.push_back(gsud);
  }

}





template<Dipole::Type DT>
IISudakov_Group<DT>::~IISudakov_Group() {
  for(list<Sudakov_Base*>::iterator bit=l_sud.begin(); bit!=l_sud.end(); ++bit)
    if(*bit) delete (*bit);
#ifdef TEMP_OUTPUT
  cout<<"~IISudakov_Group"<<endl;//////////////////////////////////////////////
#endif
}



//-----------------------------------------------------------------------------



template<Dipole::Type DT>
const bool IISudakov_Group<DT>::GenerateVariablesFor(const Dipole& dip,
						     Sudakov_Result& sur) {

  static IISudakov_Stats efstat(DT,"EfracTest",false);
  static IISudakov_Stats nfstat(DT,"NfVeto",false);
  static IISudakov_Stats xmstat(DT,"XMinus",false);
  static IISudakov_Stats xpstat(DT,"XPlus",false);
  static IISudakov_Stats costat(DT,"CorrWeight");
  static IISudakov_Stats pfstat(DT,"PDFWeight",false);

  sur.Reset();

  if(l_sud.empty()) { assert(p_dip==NULL && p_sur==NULL); return false;}

  assert(p_dip==NULL);
  p_dip=&dip;
  assert(p_sur==NULL);
  p_sur=new Sudakov_Result; assert(p_sur);
  p_sur->Isr.resize(sr::stop,0.0);
  sur.Isr.resize(sr::stop,0.0);
  if(InitWithCurrentDipole()==false) {
    p_dip=NULL; delete p_sur; p_sur=NULL; sur.Reset(); return false;
  }

  xbool            kspl=between;
  Radiation::Group kgrp=Radiation::gluon;
  extern Run_Parameter ATOOLS::rpa;
  double Ecm=rpa.gen.Ecms();
  p_sur->Isr.push_back(0.0);    //Trick to have easy carriage.
  double& x1=p_sur->Isr.back();
  p_sur->Isr.push_back(0.0);
  double& x3=p_sur->Isr.back();

  for(list<Sudakov_Base*>::const_iterator bit=l_sud.begin();
      bit!=l_sud.end(); ++bit) {

    Sudakov_Base& suda=**bit;
    suda.SetPDFFlavs(p_sur->Dir);    //The radiation flavour is now fixed.
    xbool emi=suda.RadPart();
    if(emi==between) emi=xbool(p_sur->Dir);
    else {
      if(suda.RadGroup()==Radiation::quark) {
	if(emi==front) emi=xbool(p_sur->Dir); else emi=xbool(!p_sur->Dir);
      } else {
	size_t stop=m_radtype%10;
	if(emi==front) {
	  if(size_t((*suda.RadCode().Qua)().Kfcode()) > stop) continue;
	  emi=xbool(p_sur->Dir);
	} else {
	  if(size_t((*suda.RadCode().Aqu)().Kfcode()) > stop) continue;
	  emi=xbool(!p_sur->Dir);
	}
      }
    }

    suda.RadCode().Print(); cout<<endl;
    this->Reset();
    suda.SetGenX2tFac();

    while(true) {
      m_x2t=suda.GenerateX2t();
      if(m_x2t < 0.0) break;
      m_ymax=suda.CalculateRapLimit();
      m_rap=suda.GenerateRap();
      double h[2]; h[0]=std::exp(m_rap); h[1]=1.0/h[0];    //=std::exp(-m_rap);
      x3=m_sdip*m_x2t*m_s;    //This is m2D*p2t.
      p_sur->Isr[sr::shat]=m_sdip+sqrt(x3)*(h[0]+h[1]);
      if(p_sur->Isr[sr::shat]>sqr(Ecm)) continue;
      p_sur->Isr[sr::kt]=x3/p_sur->Isr[sr::shat];    //This is k2t!
      //If wished for, the explicit k2t cutoff must be commented in here:
      //if(p_sur->Isr[sr::kt] < dpa.sud.MinIIK2t()) break;
      p_sur->Isr[sr::mt]=sqrt(m_sdip+p_sur->Isr[sr::kt]);
      p_sur->Isr[sr::kt]=sqrt(p_sur->Isr[sr::kt]);    //Make it now linear.
      //Scales can be p2t, k2t or m2t, that determines the place of the veto!
      //Subsequent scales will be descending!
      if(suda.NfVeto(m_s*m_x2t)) { nfstat.Include(0.77); break;}
      //Dipole direction dependence if g emit type emission:
      //-+dir: h[1] & p_sur->Dir=1, +-dir: h[0] & p_sur->Dir=0
      p_sur->Isr[sr::expy]=
	p_sur->Isr[sr::expydip] *
	(sqrt(x3)*h[/*1*/emi]-sqr(p_sur->Isr[sr::kt])) /
	(p_sur->Isr[sr::kt]*p_sur->Isr[sr::mt]);
      assert(p_sur->Isr[sr::expy]>0.0);
      p_sur->Isr[sr::xpfin]=
	(p_sur->Isr[sr::kt]*p_sur->Isr[sr::expy]+
	 p_sur->Isr[sr::mt]*p_sur->Isr[sr::expydip])/Ecm;
      p_sur->Isr[sr::xmfin]=
	p_sur->Isr[sr::shat]/(p_sur->Isr[sr::xpfin]*sqr(Ecm));
      xpstat.Include(p_sur->Isr[sr::xpfin]);/////////////////////////////////
      xmstat.Include(p_sur->Isr[sr::xmfin]);/////////////////////////////////
      if(p_sur->Isr[sr::xpfin]>1.0 || p_sur->Isr[sr::xmfin]>1.0) continue;
      x3=sqrt(x3/sqr(m_sdip));
      //anni:xq(xg)  //comp:xq(b) (always remaining flav)  //split:x_i_keep
      x1=1.0+x3*h[1];
      //anni:xqb(xg)  //comp:xg (always new gluon)  //split:x_i_new
      x3=1.0+x3*h[0];
      //p_sur->Print();//////////////////////////////////////////////////////
      m_corr=suda.GenerateCorr(p_sur->Isr);
      double ran=ATOOLS::ran.Get();
      if(ran > m_corr) continue;
      double pfcorr=suda.GeneratePDFCorr(p_sur->Isr);
      m_corr*=pfcorr;
      //Without analysing: m_corr*=suda.GeneratePDFCorr(p_sur->Isr);
      pfstat.Include(pfcorr);////////////////////////////////////////////////
      costat.Include(m_corr);////////////////////////////////////////////////
      if(ran > m_corr) continue;
      if(TestEfracs(x1,x3)) {
	if(sur.P2t > m_x2t) break;    //Keep the old values.
	//Otherwise: Temporarily keeping the new values.
	kspl=suda.RadPart();
	kgrp=suda.RadGroup();
	for(size_t i=0; i<sr::stop; ++i) sur.Isr[i]=p_sur->Isr[i];
	sur.Dir=p_sur->Dir;
	sur.Sfc=suda.RadCode();
	sur.P2t=m_x2t;
	sur.Y=m_rap;
	sur.X1=x1;    //ann:xq  //comp:xq(b) (always remaining flav) //split:xk
	sur.X3=x3;    //ann:xqb //comp:xg (always new gluon)         //split:xn
	break;
      }
      efstat.Include(7.77);///////////////////////////////////////////////////
    }

  }

  p_dip=NULL;
  delete p_sur; p_sur=NULL;
  sur.P2t*=m_s;    //Making sur.P2t a real p2t.

  if(sur.Sfc.Gluic() || sur.Sfc.Quaic() || sur.Sfc.Aquic()) {
    switch(kgrp) {
    case Radiation::gluon:
      sur.Rad=Radiation::gluon;
      assert(kspl==between); assert(sur.Sfc.Gluic()); break;
    case Radiation::igluon:
      if(kspl==front) {
	sur.Rad=Radiation::qfront;
	assert(sur.Sfc.Quaic()); break;
      } else {
	sur.Rad=Radiation::qbarend;
	assert(kspl==end); assert(sur.Sfc.Aquic()); break;
      }
    case Radiation::quark:
      assert(!sur.Sfc.Gluic());
      if(kspl==front) { sur.Rad=Radiation::qitop; break;}
      else { sur.Rad=Radiation::qibot; assert(kspl==end); break;}
    default:
      assert(0);
    }
    return true;
  }
  else {
    sur.Reset();
    return false;
  }

}



//-----------------------------------------------------------------------------



template<Dipole::Type DT>
const bool IISudakov_Group<DT>::InitRadiation() const {
  cout<<__PRETTY_FUNCTION__<<endl;/////////////////////////////////////////////
  for(list<Sudakov_Base*>::const_iterator cit=l_sud.begin();
      cit!=l_sud.end(); ++cit)
    if(*cit) (*cit)->InitRadParticle();
  return true;
}





template<Dipole::Type DT>
const bool IISudakov_Group<DT>::InitWithCurrentDipole() {

  static bool chk=InitRadiation();    //Warning removal trick, use chk again.

  chk=(p_dip->IsType()==DT); assert(chk);

  extern Run_Parameter ATOOLS::rpa;
  double sqrtS=rpa.gen.Ecms();
  const Vec4D& quap=p_dip->GetBotBranchPointer()->Momentum();//Incoming q/g.
  const Vec4D& antp=p_dip->GetTopBranchPointer()->Momentum();//Incoming qbar/g.

  double pplus=quap.PPlus(), pminus;
  if(pplus<quap.PMinus()) {    //Dipole direction is then +(qbar) to -(q) beam.
    p_sur->Dir=false;
    pplus=antp.PPlus();
    pminus=quap.PMinus();
  } else {    //Dipole direction is then -(qbar) to +(q) beam.
    pminus=antp.PMinus();
  }
  p_sur->Isr[sr::xpini]=pplus/sqrtS;
  p_sur->Isr[sr::xmini]=pminus/sqrtS;
  assert(p_sur->Isr[sr::xpini]<=1.0);
  assert(p_sur->Isr[sr::xmini]<=1.0);

  static IISudakov_Stats dipstat(DT,"DipoleMass",false,40.0,270.0);
  static IISudakov_Stats saxstat(DT,"-/SHatMax",false,70.0,sqrtS);

  m_s=m_sdip=p_dip->InvMass();
  p_sur->Isr[sr::mdip]=sqrt(m_sdip);
  p_sur->Isr[sr::expydip]=pplus/p_sur->Isr[sr::mdip];
  dipstat.Include(p_sur->Isr[sr::mdip]);

  //dpv.sud.SetMaxIIScale(/*m_sdip*/);    //Static treatment!!!!!!!!!!!!!!!!!!!

  m_a=dpa.sud.MaxIIK2t()/m_sdip;
  m_a=1.0+2*m_a+2*sqrt(m_a*(1.0+m_a));
  assert(m_a>1.0);
  p_sur->Isr[sr::shatmax]=m_a*m_sdip;    //Only for control purposes.
  saxstat.Include(sqrt(p_sur->Isr[sr::shatmax]));//////////////////////////////
  if(p_sur->Isr[sr::shatmax]>sqr(sqrtS)) {
    p_sur->Isr[sr::shatmax]=sqr(sqrtS);
    m_a=p_sur->Isr[sr::shatmax]/m_sdip;
  }
  //saxstat.Include(sqrt(p_sur->Isr[sr::shatmax]));////////////////////////////

  m_s*=sqr(m_a-1.0);

  m_x2tmin=dpa.sud.MinIIK2t()/m_s;    //Use the direct treatment.
  m_x2tmax=Min(0.25,p_dip->BootScale()/m_s);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //m_x2tmax=0.25;//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  assert(m_x2tmin<m_x2tmax);

  //cout<<"-/shatmax="<<sqrt(p_sur->Isr[sr::shatmax])//////////////////////////
  //    <<"  m_a="<<m_a<<"  k2tmax="<<dpa.sud.MaxIIK2t()<<endl;////////////////
  //if(m_x2tmax<0.249999) cout<<"m_x2tmax="<<m_x2tmax<<"\n"<<*p_dip<<endl;/////
  //cout<<"p_dip->BootScale()/m_s="<<p_dip->BootScale()/m_s<<endl;/////////////

  if(m_x2tmin<m_x2tmax) return true;
  return false;

}



//=============================================================================



//===============
//Gluon emission.
//===============


template<Dipole::Type DT>
IISudakov<DT,Radiation::gluon>::IISudakov(const IISudakov_Group<DT>& sg,
					  const Sudakov_Flavour& go)
  : Sudakov_Base(go,Radiation::gluon),
    m_genx2tfac(0.0), m_mfl(), m_sgroup(sg) {

  assert(go.Gluic());
}





template<Dipole::Type DT>
IISudakov<DT,Radiation::gluon>::~IISudakov() {
#ifdef TEMP_OUTPUT
  std::cout<<"~IISudakov(gluon)"<<std::endl;///////////////////////////////////
#endif
}





template<Dipole::Type DT>
const double IISudakov<DT,Radiation::gluon>::GenerateX2t() {
  m_step=0;
  double ran=ATOOLS::ran.Get();
  double coeff=
    std::log(ran)*m_genx2tfac*s_colfac/Sudakov_Calculator::AlphaSApprox();
  double A=sqr(std::log(m_sgroup.X2t()));
  if( coeff < A-sqr(std::log(m_sgroup.X2tmin())) ) return -1.0;
  return std::exp(-sqrt(A-coeff));
}





//=============================================
//I gluon with associated (anti)quark emission.
//=============================================



template<Dipole::Type DT>
IISudakov<DT,Radiation::igluon>::IISudakov(const IISudakov_Group<DT>& sg,
					   const Sudakov_Flavour& qo)
  : Sudakov_Base(qo,Radiation::igluon),
    m_genx2tfac(0.0), m_mfl(), m_sgroup(sg) {

  if(qo.Quaic()) { m_split=front; return;}
  if(qo.Aquic()) { m_split=end; return;}
  assert(qo.Quaic() || qo.Aquic());
}





template<Dipole::Type DT>
IISudakov<DT,Radiation::igluon>::~IISudakov() {
#ifdef TEMP_OUTPUT
  std::cout<<"~IISudakov(igluon)"<<std::endl;//////////////////////////////////
#endif
}





template<Dipole::Type DT>
const double IISudakov<DT,Radiation::igluon>::GenerateX2t() {
  m_step=0;
  double ran=ATOOLS::ran.Get();
  double coeff=
    std::log(ran)*m_genx2tfac*s_colfac/Sudakov_Calculator::AlphaSApprox();
  double A=sqr(std::log(m_sgroup.X2t()));
  if( coeff < A-sqr(std::log(m_sgroup.X2tmin())) ) return -1.0;
  return std::exp(-sqrt(A-coeff));
}





//==========================
//Quark-antiquark splitting.
//==========================



template<Dipole::Type DT>
IISudakov<DT,Radiation::quark>::IISudakov(const IISudakov_Group<DT>& sg,
					  const Sudakov_Flavour& qo)
  : Sudakov_Base(qo,Radiation::quark),
    m_genx2tfac(0.0), m_inicode(qo), v_mfl(4), m_sgroup(sg) {

  assert(qo.Bqqic());
  m_split=xbool(-1+2*(DT==Dipole::iigq));
  assert(m_split!=between);////////////////////////////////////////////////////
}





template<Dipole::Type DT>
IISudakov<DT,Radiation::quark>::~IISudakov() {
  for(size_t i=0; i<v_mfl.size(); ++i) v_mfl[i].clear();
  v_mfl.clear();
#ifdef TEMP_OUTPUT
  std::cout<<"~IISudakov(quark)"<<std::endl;///////////////////////////////////
#endif
}





template<Dipole::Type DT>
const double IISudakov<DT,Radiation::quark>::GenerateX2t() {
  m_step=0;
  double ran=ATOOLS::ran.Get();
  double coeff=m_genx2tfac*s_colfac/Sudakov_Calculator::AlphaSApprox();
  double x2tcal=m_sgroup.X2t();
  if(coeff*std::log(ran) < std::log(m_sgroup.X2tmin()/x2tcal)) return -1.0;
  x2tcal*=std::pow(ran,coeff);
  return x2tcal;
}



//=============================================================================





//eof
