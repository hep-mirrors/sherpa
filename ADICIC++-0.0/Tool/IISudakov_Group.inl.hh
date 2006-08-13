//bof
//Version: 4 ADICIC++-0.0/2006/08/12

//Inline methods of IISudakov_Group.H.





#include <cassert>
#include <cstdlib>
#include <mathextra>
#include "Random.H"


#define IISUDAKOV_K2T_RSCALE IISUDAKOV_K2T_RSCALE
#undef  IISUDAKOV_K2T_RSCALE





using ATOOLS::sqr;



namespace ADICIC {



  //===========================================================================



  template<> inline void IISudakov_Group<Dipole::iiqbarq>::Which() const {
    std::cout<<"Sudakov_Calculator for an ii qbar-q dipole using ";
    if(IsAlphaSRunning()) std::cout<<"running alpha_s.\n";
    else std::cout<<"fixed alpha_s.\n";
  }
  template<> inline void IISudakov_Group<Dipole::iiqbarg>::Which() const {
    std::cout<<"Sudakov_Calculator for an ii qbar-g dipole using ";
    if(IsAlphaSRunning()) std::cout<<"running alpha_s.\n";
    else std::cout<<"fixed alpha_s.\n";
  }
  template<> inline void IISudakov_Group<Dipole::iigq>::Which() const {
    std::cout<<"Sudakov_Calculator for an ii g-q dipole using ";
    if(IsAlphaSRunning()) std::cout<<"running alpha_s.\n";
    else std::cout<<"fixed alpha_s.\n";
  }
  template<> inline void IISudakov_Group<Dipole::iigg>::Which() const {
    std::cout<<"Sudakov_Calculator for an ii g-g dipole using ";
    if(IsAlphaSRunning()) std::cout<<"running alpha_s.\n";
    else std::cout<<"fixed alpha_s.\n";
  }


  template<Dipole::Type DT>
  inline void IISudakov_Group<DT>::ShowSpecification() const {
    std::cout<<"IISudakov_Group for dipole and radiation type: "
	     <<DT<<" and "<<m_radtype<<".\n";
    if(l_sud.empty()) { std::cout<<" No IISudakov's assigned.\n"; return;}
    short j=1;
    for(std::list<Sudakov_Base*>::const_iterator bit=l_sud.begin();
	bit!=l_sud.end(); ++bit) {
      std::cout<<" "<<j<<"."; (*bit)->Which();
      std::cout<<" "<<j<<"."; (*bit)->ShowSpecification();
      ++j;
    }
  }



  //---------------------------------------------------------------------------



  template<Dipole::Type DT>
  inline const double IISudakov_Group<DT>::A() const {
    return m_a;
  }

  template<Dipole::Type DT>
  inline const double IISudakov_Group<DT>::Scale() const {
    return m_s;
  }

  template<Dipole::Type DT>
  inline const double IISudakov_Group<DT>::X2tmin() const {
    return m_x2tmin;
  }

  template<Dipole::Type DT>
  inline const double IISudakov_Group<DT>::X2t() const {
    return m_x2t;
  }

  template<Dipole::Type DT>
  inline const double IISudakov_Group<DT>::Ymax() const {
    return m_ymax;
  }



  //---------------------------------------------------------------------------



  template<Dipole::Type DT>
  inline const bool IISudakov_Group<DT>::TestEfracs(const double x1,
						    const double x3) const {
    double sum=x1+x3;
    if(x1>=1.0 && x3>=1.0 && sum>2.0 && sum<=1.0+m_a) return true;
    //assert(0);
    return false;
  }

  template<Dipole::Type DT>
  inline void IISudakov_Group<DT>::Reset() {
    m_x2t=m_x2tmax;
    m_ymax=0.0;
    m_rap=0.0;
    m_corr=1.0;
  }



  //===========================================================================



  template<>
  inline void IISudakov<Dipole::iiqbarq,Radiation::gluon>::Which() const {
    std::cout<<"Sudakov_Base for: ii qbar-q dipole producing g emission.\n";
  }
  template<>
  inline void IISudakov<Dipole::iiqbarq,Radiation::igluon>::Which() const {
    std::cout<<"Sudakov_Base for: ii qbar-q dipole producing (anti)q "
	     <<"emission.\n";
  }


  template<>
  inline void IISudakov<Dipole::iiqbarg,Radiation::gluon>::Which() const {
    std::cout<<"Sudakov_Base for: ii qbar-g dipole producing g emission.\n";
  }
  template<>
  inline void IISudakov<Dipole::iiqbarg,Radiation::igluon>::Which() const {
    std::cout<<"Sudakov_Base for: ii qbar-g dipole producing q emission.\n";
  }
  template<>
  inline void IISudakov<Dipole::iiqbarg,Radiation::quark>::Which() const {
    std::cout<<"Sudakov_Base for: ii qbar-g dipole producing q-qbar"
	     <<" splitting.\n";
  }


  template<>
  inline void IISudakov<Dipole::iigq,Radiation::gluon>::Which() const {
    std::cout<<"Sudakov_Base for: ii g-q dipole producing g emission.\n";
  }
  template<>
  inline void IISudakov<Dipole::iigq,Radiation::igluon>::Which() const {
    std::cout<<"Sudakov_Base for: ii g-q dipole producing antiq emission.\n";
  }
  template<>
  inline void IISudakov<Dipole::iigq,Radiation::quark>::Which() const {
    std::cout<<"Sudakov_Base for: ii g-q dipole producing q-qbar splitting.\n";
  }


  template<>
  inline void IISudakov<Dipole::iigg,Radiation::gluon>::Which() const {
    std::cout<<"Sudakov_Base for: ii g-g dipole producing g emission.\n";
  }



  template<Dipole::Type DT>
  inline void IISudakov<DT,Radiation::gluon>::ShowSpecification() const {
    std::cout<<"IISudakov for: dip.type="<<DT
	     <<",  rad.group="<<Radiation::gluon
	     <<" and flav.code="<<(*m_code.Glu)().Kfcode()
	     <<" ("<<m_mass<<" GeV)"
	     <<",  ME.powers="<<s_x1pow<<","<<s_x3pow
	     <<",  col.factor="<<s_colfac
	     <<",  eff.base="<<s_iieffbas
	     <<",  PDF.approx="<<s_pdfapprox<<".\n";
  }
  template<Dipole::Type DT>
  inline void IISudakov<DT,Radiation::igluon>::ShowSpecification() const {
    Radiation::Group rg=Radiation::qfront;
    if(m_code.Aquic()) rg=Radiation::qbarend;
    std::cout<<"IISudakov for: dip.type="<<DT
	     <<",  rad.group="<<Radiation::igluon<<"("<<rg<<")"
	     <<",  col.factor="<<s_colfac
	     <<",  eff.base="<<s_iieffbas
	     <<",  PDF.approx="<<s_pdfapprox<<".\n";
  }
  template<Dipole::Type DT>
  inline void IISudakov<DT,Radiation::quark>::ShowSpecification() const {
    std::cout<<"IISudakov for: dip.type="<<DT
	     <<",  rad.group="<<Radiation::quark
	     <<" and flav.code="<<(*m_inicode.Qua)().Kfcode()
	     <<" ("<<m_mass<<" GeV)"
	     <<",  col.factor="<<s_colfac
	     <<",  eff.base="<<s_iieffbas
	     <<",  PDF.approx="<<s_pdfapprox<<".\n";
  }



  //---------------------------------------------------------------------------



  //===============
  //Gluon emission.
  //===============



  template<Dipole::Type DT> inline void
  IISudakov<DT,Radiation::gluon>::SetPDFFlavs(bool dipdir) {
    //dipdir=true: +beam & iantibranch(q), -beam & ibranch(qbar).
    m_mfl.clear();
    ATOOLS::Flavour flav[2];
    flav[!dipdir]=m_sgroup.CurrentDipole().GetBotBranchPointer()->Flav();
    flav[dipdir]=m_sgroup.CurrentDipole().GetTopBranchPointer()->Flav();
    m_mfl.push_back(flav[0]);    //+ini
    m_mfl.push_back(flav[1]);    //-ini
    m_mfl.push_back(flav[0]);    //+fin
    m_mfl.push_back(flav[1]);    //-fin
    //std::cout<<m_mfl[sf::plusini]<<", "<<m_mfl[sf::miusini]<<", "////////////
    //	       <<m_mfl[sf::plusfin]<<", "<<m_mfl[sf::miusfin]<<"\n";///////////
  }



  template<> inline void IISudakov<Dipole::iiqbarq,Radiation::gluon>
  ::SetGenX2tFac() {
    static const double eff=std::pow(s_iieffbas,dpa.sud.IIEffExp());
    static Sudakov_Stats effstat(Dipole::iiqbarq,eff,"gEff");
    m_genx2tfac=eff; if(m_pdfr!=on) return;
    const Multidouble& is=m_sgroup.CurrentSudakovResult().Isr;
    m_pdfapprox[0]=s_pdf1stat.GiveEstimate(is[sr::xpini],s_pdfapprox);    //QQ
    m_pdfapprox[1]=s_pdf1stat.GiveEstimate(is[sr::xmini],s_pdfapprox);    //QQ
    //std::cout<<m_pdfapprox[0]<<" "<<m_pdfapprox[1]<<std::endl;
  }
  template<Dipole::Type DT> inline void IISudakov<DT,Radiation::gluon>
  ::SetGenX2tFac() {
    static const double eff=std::pow(s_iieffbas,dpa.sud.IIEffExp());
    static Sudakov_Stats effstat(DT,eff,"gEff");
    if(m_sgroup.CurrentDipole().SpinCorr()==false) m_genx2tfac=eff;
    else m_genx2tfac=2.0*eff/(m_sgroup.A()+1.0);
    if(m_pdfr!=on) return;    //resort fixed to 10^0.602..=4
    static const double resort=std::pow(s_pdfapprox,0.60205999);
    const Multidouble& is=m_sgroup.CurrentSudakovResult().Isr;
    if(m_mfl[sf::plusfin].Kfcode()==ATOOLS::kf::gluon) {    //+:GG=2 //-:QQ=1
      m_pdfapprox[0]=s_pdf2stat.GiveEstimate(is[sr::xpini],s_pdfapprox);
      m_pdfapprox[1]=s_pdf1stat.GiveEstimate(is[sr::xmini],resort);
    } else {    //+:QQ=1//-:GG=2
      m_pdfapprox[0]=s_pdf1stat.GiveEstimate(is[sr::xpini],resort);
      m_pdfapprox[1]=s_pdf2stat.GiveEstimate(is[sr::xmini],s_pdfapprox);
    }
    //std::cout<<m_pdfapprox[0]<<" "<<m_pdfapprox[1]<<std::endl;
  }
  template<> inline void IISudakov<Dipole::iigg,Radiation::gluon>
  ::SetGenX2tFac() {
    static const double eff=std::pow(s_iieffbas,dpa.sud.IIEffExp());
    static Sudakov_Stats effstat(Dipole::iigg,eff,"gEff");
    if(m_sgroup.CurrentDipole().SpinCorr()==false) m_genx2tfac=eff;
    else m_genx2tfac=2.0*eff/(m_sgroup.A()+1.0);
    if(m_pdfr!=on) return;    //resort fixed to 25^0.715..=10
    static const double resort=std::pow(s_pdfapprox,0.71533828);
    const Multidouble& is=m_sgroup.CurrentSudakovResult().Isr;
    m_pdfapprox[0]=s_pdf2stat.GiveEstimate(is[sr::xpini],resort);    //GG
    m_pdfapprox[1]=s_pdf2stat.GiveEstimate(is[sr::xmini],resort);    //GG
    //std::cout<<m_pdfapprox[0]<<" "<<m_pdfapprox[1]<<std::endl;
  }



  template<> inline double IISudakov<Dipole::iiqbarq,Radiation::gluon>
  ::GeneratePDFCorr() {
    assert(m_step==3); ++m_step;
    double pwgt=m_sgroup.PlusPDFCorr(m_mfl), mwgt=m_sgroup.MinusPDFCorr(m_mfl);
    if(m_pdfr==on) {
      const Multidouble& is=m_sgroup.CurrentSudakovResult().Isr;
      s_pdf1stat.Include(pwgt,m_pdfapprox[0],is[sr::xpini]);
      s_pdf1stat.Include(mwgt,m_pdfapprox[1],is[sr::xmini]);
      return pwgt*mwgt/(m_pdfapprox[0]*m_pdfapprox[1]);
    }
    pwgt/=sqrt(m_pdfapx); mwgt/=sqrt(m_pdfapx);
    s_pdf1stat.Include(pwgt); s_pdf1stat.Include(mwgt);
    return pwgt*mwgt;
  }
  template<Dipole::Type DT> inline double IISudakov<DT,Radiation::gluon>
  ::GeneratePDFCorr() {
    assert(m_step==3); ++m_step;
    double pwgt=m_sgroup.PlusPDFCorr(m_mfl), mwgt=m_sgroup.MinusPDFCorr(m_mfl);
    if(m_pdfr==on) {
      const Multidouble& is=m_sgroup.CurrentSudakovResult().Isr;
      if(m_mfl[sf::plusfin].Kfcode()==ATOOLS::kf::gluon) {
	s_pdf2stat.Include(pwgt,m_pdfapprox[0],is[sr::xpini]);    //GG
	s_pdf1stat.Include(mwgt,m_pdfapprox[1],is[sr::xmini]);    //QQ
	return pwgt*mwgt/(m_pdfapprox[0]*m_pdfapprox[1]);
      }
      s_pdf1stat.Include(pwgt,m_pdfapprox[0],is[sr::xpini]);    //QQ
      s_pdf2stat.Include(mwgt,m_pdfapprox[1],is[sr::xmini]);    //GG
      return pwgt*mwgt/(m_pdfapprox[0]*m_pdfapprox[1]);
    }
    if(m_mfl[sf::plusfin].Kfcode()==ATOOLS::kf::gluon) {
      pwgt/=std::pow(m_pdfapx,0.69897); mwgt/=std::pow(m_pdfapx,0.30103);
      s_pdf2stat.Include(pwgt); s_pdf1stat.Include(mwgt);
      return pwgt*mwgt;
    }
    pwgt/=std::pow(m_pdfapx,0.30103); mwgt/=std::pow(m_pdfapx,0.69897);
    s_pdf1stat.Include(pwgt); s_pdf2stat.Include(mwgt);
    return pwgt*mwgt;
  }
  template<> inline double IISudakov<Dipole::iigg,Radiation::gluon>
  ::GeneratePDFCorr() {
    assert(m_step==3); ++m_step;
    double pwgt=m_sgroup.PlusPDFCorr(m_mfl), mwgt=m_sgroup.MinusPDFCorr(m_mfl);
    if(m_pdfr==on) {
      const Multidouble& is=m_sgroup.CurrentSudakovResult().Isr;
      s_pdf2stat.Include(pwgt,m_pdfapprox[0],is[sr::xpini]);
      s_pdf2stat.Include(mwgt,m_pdfapprox[1],is[sr::xmini]);
      return pwgt*mwgt/(m_pdfapprox[0]*m_pdfapprox[1]);
    }
    pwgt/=sqrt(m_pdfapx); mwgt/=sqrt(m_pdfapx);
    s_pdf2stat.Include(pwgt); s_pdf2stat.Include(mwgt);
    return pwgt*mwgt;
  }



  template<Dipole::Type DT>
  inline const double IISudakov<DT,Radiation::gluon>::CalculateRapLimit() {
    assert(m_step==0); ++m_step;
    double ratio=0.25/m_sgroup.X2t();
    return std::log( sqrt(ratio)+sqrt(ratio-1.0) );
  }



  template<Dipole::Type DT>
  inline const double IISudakov<DT,Radiation::gluon>::GenerateRap() {
    assert(m_step==1); ++m_step;
    return m_sgroup.Ymax()*(-1.0+2.0*ATOOLS::ran.Get());
  }



  template<Dipole::Type DT>
  inline const double IISudakov<DT,Radiation::gluon>::GenerateCorr() {
    static Sudakov_Stats mestat(DT,"MEgWeight",false);
    assert(m_step==2); ++m_step;
    const double& x1=m_sgroup.CurrentSudakovResult().X1;    //xqb(xg).
    const double& x3=m_sgroup.CurrentSudakovResult().X3;    //xq(xg).
    double fac;
    double x2t=m_sgroup.X2t();
    const Multidouble& isrx=m_sgroup.CurrentSudakovResult().Isr;
    assert(isrx.size()==sr::stop);
    if(m_sgroup.CurrentDipole().SpinCorr()==false) fac=2/s_average;
    //else fac=(power<s_x1pow>(x1) + power<s_x3pow>(x3)) *    //gcc4.0 dislike.
    else fac=(power(x1,s_x1pow) + power(x3,s_x3pow)) *
	   ATOOLS::sqr(ATOOLS::sqr(isrx[sr::mdip])/isrx[sr::shat])
	   //This is the additional flux correction.
	   ;
    x2t=
#ifdef IISUDAKOV_K2T_RSCALE
      Sudakov_Calculator::AlphaSCorr(ATOOLS::sqr(isrx[sr::kt])) *
#else
      Sudakov_Calculator::AlphaSCorr(m_sgroup.Scale()*x2t) *
#endif
      (-m_sgroup.Ymax()/std::log(x2t)) *
      //The missing "2" is either compensated by the missing "1/2" of the
      //power correction or directly introduced.
      m_genx2tfac * fac;
    mestat.Include(x2t);/////////////////////////////////////////////////////
    //std::cout<<fac<<std::endl;///////////////////////////////////////////////
    return x2t;
  }





  //=============================================
  //I gluon with associated (anti)quark emission.
  //=============================================



  template<Dipole::Type DT>
  inline bool IISudakov<DT,Radiation::igluon>::NfVeto(const double sl) const {
    assert(m_step==2);
    switch(m_split) {
    case front:
      if((*m_code.Qua)().Kfcode() > Sudakov_Calculator::Nf(sl)) return true;
      return false;
    case   end:
      if((*m_code.Aqu)().Kfcode() > Sudakov_Calculator::Nf(sl)) return true;
      return false;
    default   : assert(m_split==front || m_split==end);
    }
    return false;
  }



  template<Dipole::Type DT>
  inline void IISudakov<DT,Radiation::igluon>::SetPDFFlavs(bool beamdipdir) {
    //beamdipdir=true, i.e. dipole goes from -(qbar) to +(q) beam.
    m_mfl.clear();
    ATOOLS::Flavour ifla[2];
    ifla[!beamdipdir]=m_sgroup.CurrentDipole().GetBotBranchPointer()->Flav();
    ifla[beamdipdir]=m_sgroup.CurrentDipole().GetTopBranchPointer()->Flav();
    m_mfl.push_back(ifla[0]);    //+ini
    m_mfl.push_back(ifla[1]);    //-ini
    switch(m_split) {
    case front: {    //TopBranch determines emission flavour (qg->Bq).
      m_code.Qua=info.quark.pkf[ifla[beamdipdir].Kfcode()];
      assert(m_code.Qua);
      m_mass=ifla[beamdipdir].Mass();
      ifla[beamdipdir]=ATOOLS::Flavour(ATOOLS::kf::gluon);
      m_mfl.push_back(ifla[0]);    //+fin
      m_mfl.push_back(ifla[1]);    //-fin
    } break;
    case   end: {    //BotBranch determines emission flavour (gqbar->qbarB).
      m_code.Aqu=info.antiq.pkf[ifla[!beamdipdir].Kfcode()];
      assert(m_code.Aqu);
      m_mass=ifla[!beamdipdir].Mass();
      ifla[!beamdipdir]=ATOOLS::Flavour(ATOOLS::kf::gluon);
      m_mfl.push_back(ifla[0]);    //+fin
      m_mfl.push_back(ifla[1]);    //-fin
    } break;
    default   : assert(m_split==front || m_split==end);
    }
    //std::cout<<m_mfl[sf::plusini]<<", "<<m_mfl[sf::miusini]<<", "////////////
    //	       <<m_mfl[sf::plusfin]<<", "<<m_mfl[sf::miusfin]<<"\n";///////////
  }



  template<Dipole::Type DT>
  inline void IISudakov<DT,Radiation::igluon>::SetGenX2tFac() {
    static const double eff=std::pow(s_iieffbas,dpa.sud.IIEffExp());
    static Sudakov_Stats effstat(DT,eff,"qEff");
    m_genx2tfac=eff;
    if(m_pdfr!=on) return;
    static const sf::code c=(DT==Dipole::iiqbarq ? sf::plusfin : sf::miusini);
    static const double pdfapprox=
      (DT==Dipole::iiqbarq ? s_pdfapprox : std::pow(s_pdfapprox,0.80103));
    static const double resort=
      (DT==Dipole::iiqbarq ?
       std::pow(s_pdfapprox,0.37580365) :    //fixed to 40^0.37..=4
       sqrt(s_pdfapprox));    //fixed to 100^0.5=10
    const Multidouble& is=m_sgroup.CurrentSudakovResult().Isr;
    if(m_mfl[c].Kfcode()==ATOOLS::kf::gluon) {    //+:GQ=2 //-:QQ/GG=1
      m_pdfapprox[0]=s_pdf2stat.GiveEstimate(is[sr::xpini],pdfapprox);
      m_pdfapprox[1]=s_pdf1stat.GiveEstimate(is[sr::xmini],resort);
    } else {    //+:QQ/GG=1 //-:GQ=2
      m_pdfapprox[0]=s_pdf1stat.GiveEstimate(is[sr::xpini],resort);
      m_pdfapprox[1]=s_pdf2stat.GiveEstimate(is[sr::xmini],pdfapprox);
    }
  }



  template<Dipole::Type DT>
  inline double IISudakov<DT,Radiation::igluon>::GeneratePDFCorr() {
    static const sf::code c=(DT==Dipole::iiqbarq ? sf::plusfin : sf::miusini);
    assert(m_step==3); ++m_step;
    double pwgt=m_sgroup.PlusPDFCorr(m_mfl);
    double mwgt=m_sgroup.MinusPDFCorr(m_mfl);
    if(m_pdfr==on) {
      const Multidouble& is=m_sgroup.CurrentSudakovResult().Isr;
      if(m_mfl[c].Kfcode()==ATOOLS::kf::gluon) {
	s_pdf2stat.Include(pwgt,m_pdfapprox[0],is[sr::xpini]);    //GQ
	s_pdf1stat.Include(mwgt,m_pdfapprox[1],is[sr::xmini]);    //QQ/GG
	return pwgt*mwgt/(m_pdfapprox[0]*m_pdfapprox[1]);
      }
      s_pdf1stat.Include(pwgt,m_pdfapprox[0],is[sr::xpini]);    //QQ/GG
      s_pdf2stat.Include(mwgt,m_pdfapprox[1],is[sr::xmini]);    //GQ
      return pwgt*mwgt/(m_pdfapprox[0]*m_pdfapprox[1]);
    }
    static const double expon=(DT==Dipole::iiqbarq ? 0.81209818 : 0.650515);
    if(m_mfl[c].Kfcode()==ATOOLS::kf::gluon) {
      pwgt/=std::pow(m_pdfapx,expon); mwgt/=std::pow(m_pdfapx,1.0-expon);
      s_pdf2stat.Include(pwgt); s_pdf1stat.Include(mwgt);
      return pwgt*mwgt;
    }
    pwgt/=std::pow(m_pdfapx,1.0-expon); mwgt/=std::pow(m_pdfapx,expon);
    s_pdf1stat.Include(pwgt); s_pdf2stat.Include(mwgt);
    return pwgt*mwgt;
  }



  template<Dipole::Type DT>
  inline const double IISudakov<DT,Radiation::igluon>::CalculateRapLimit() {
    assert(m_step==0); ++m_step;
    double ratio=0.25/m_sgroup.X2t();
    return std::log( sqrt(ratio)+sqrt(ratio-1.0) );
  }



  template<Dipole::Type DT>
  inline const double IISudakov<DT,Radiation::igluon>::GenerateRap() {
    assert(m_step==1); ++m_step;
    return m_sgroup.Ymax()*(-1.0+2.0*ATOOLS::ran.Get());
  }



  template<Dipole::Type DT>
  inline const double IISudakov<DT,Radiation::igluon>::GenerateCorr() {
    static Sudakov_Stats mestat(DT,"MEqWeight",false);
    assert(m_step==2); ++m_step;
    const double& x1=m_sgroup.CurrentSudakovResult().X1;    //xg.
    const double& x3=m_sgroup.CurrentSudakovResult().X3;    //xq(b).
    double fac;
    double x2t=m_sgroup.X2t();
    const Multidouble& isrx=m_sgroup.CurrentSudakovResult().Isr;
    assert(isrx.size()==sr::stop);
    if(m_sgroup.CurrentDipole().SpinCorr()==false)
      fac=2*(x1-1.)/(x1+x3-1.)/s_average;
    else
      fac=((power<2>(x3)+power<2>(x1+x3-2.))*(x1-1.)/(x1+x3-1.)) *
	ATOOLS::sqr(ATOOLS::sqr(isrx[sr::mdip])/isrx[sr::shat])
	//This is the additional flux correction.
	;
    x2t=
#ifdef IISUDAKOV_K2T_RSCALE
      Sudakov_Calculator::AlphaSCorr(ATOOLS::sqr(isrx[sr::kt])) *
      //#error
#else
      Sudakov_Calculator::AlphaSCorr(m_sgroup.Scale()*x2t) *
#endif
      (-m_sgroup.Ymax()/std::log(x2t)) *
      //The missing "2" is either compensated by the missing "1/2" of the
      //power correction or directly introduced.
      m_genx2tfac * fac;
    mestat.Include(x2t);/////////////////////////////////////////////////////
    //std::cout<<fac<<std::endl;///////////////////////////////////////////////
    return x2t;
  }





  //==========================
  //Quark-antiquark splitting.
  //==========================



  template<Dipole::Type DT>
  inline bool IISudakov<DT,Radiation::quark>::NfVeto(const double sl) const {
    assert(m_step==2);
    //The decision is not made yet whether we have the quark or antiquark being
    //the new initial leg. So, the Sudakov_Flavour is still a true Bqqic().
    if((*m_inicode.Qua)().Kfcode() > Sudakov_Calculator::Nf(sl)) return true;
    return false;
  }



  template<Dipole::Type DT>
  inline void IISudakov<DT,Radiation::quark>::SetPDFFlavs(bool beamdipdir) {
    //beamdipdir=true, i.e. dipole goes from -(qbar) to +(q) beam.
    for(size_t i=0; i<v_mfl.size(); ++i) v_mfl[i].clear();
    ATOOLS::Flavour ifla[2];
    ifla[!beamdipdir]=m_sgroup.CurrentDipole().GetBotBranchPointer()->Flav();
    ifla[beamdipdir]=m_sgroup.CurrentDipole().GetTopBranchPointer()->Flav();
    size_t i=0, j=v_mfl.size(); assert(j==4);
    if(DT==Dipole::iigq) j-=2; else if(DT==Dipole::iiqbarg) i=2;
    for(; i<j; ++i) {
      v_mfl[i].push_back(ifla[0]);    //+ini
      v_mfl[i].push_back(ifla[1]);    //-ini
      v_mfl[i].push_back(ifla[0]);    //pre +fin
      v_mfl[i].push_back(ifla[1]);    //pre -fin
    }
    if(!v_mfl[0].empty()) v_mfl[0][2+beamdipdir]=(*m_inicode.Qua)();
    if(!v_mfl[1].empty()) v_mfl[1][2+beamdipdir]=(*m_inicode.Aqu)();
    if(!v_mfl[2].empty()) v_mfl[2][2+!beamdipdir]=(*m_inicode.Qua)();
    if(!v_mfl[3].empty()) v_mfl[3][2+!beamdipdir]=(*m_inicode.Aqu)();
    for(size_t i=0; i<v_mfl.size(); ++i) {/////////////////////////////////////
      if(v_mfl[i].empty()) { std::cout<<i<<":\n"; continue;}///////////////////
      std::cout<<i<<": "///////////////////////////////////////////////////////
	       <<v_mfl[i][sf::plusini]<<", "<<v_mfl[i][sf::miusini]<<", "//////
    	       <<v_mfl[i][sf::plusfin]<<", "<<v_mfl[i][sf::miusfin]<<"\n";/////
    }
  }



  template<Dipole::Type DT>
  inline void IISudakov<DT,Radiation::quark>::SetGenX2tFac() {
    static double effpdfapprox=
      std::pow(s_iieffbas,dpa.sud.IIEffExp())/s_pdfapprox;
    static Sudakov_Stats effstat(DT,effpdfapprox,"gsplitEff");
    m_genx2tfac=0.5*effpdfapprox/(m_sgroup.A()-1.0);
  }



  template<Dipole::Type DT> inline double IISudakov<DT,Radiation::quark>
  ::GeneratePDFCorr() {
    static Sudakov_Stats pdfchstat(DT,"PDFQGSpl",true,0.0,s_pdfapprox/2);
    static Sudakov_Stats pdfkpstat(DT,"PDFQQSpl",true,0.0,s_pdfapprox/8);
    assert(m_step==3); ++m_step;
    size_t i; double pwgt[2], mwgt[2], wgt;
    if(m_split==front) i=0; else i=2;
    pwgt[0]=m_sgroup.PlusPDFCorr(v_mfl[i]);
    mwgt[0]=m_sgroup.MinusPDFCorr(v_mfl[i]);
    wgt=pwgt[0]*mwgt[0];
    pwgt[1]=m_sgroup.PlusPDFCorr(v_mfl[i+1]);
    mwgt[1]=m_sgroup.MinusPDFCorr(v_mfl[i+1]);
    bool win=wgt<(wgt+pwgt[1]*mwgt[1])*ATOOLS::ran.Get();
    m_code=m_inicode;
    if(win==1) m_code.Qua=NULL; else m_code.Aqu=NULL;
    if(v_mfl[i+win][sf::miusini].Kfcode()==ATOOLS::kf::gluon) {
      pdfkpstat.Include(pwgt[win]); pdfchstat.Include(mwgt[win]);
    } else {
      pdfkpstat.Include(mwgt[win]); pdfchstat.Include(pwgt[win]);
    }
    return (pwgt[win]*mwgt[win]/s_pdfapprox);
  }
  template<> inline double IISudakov<Dipole::iigg,Radiation::quark>
  ::GeneratePDFCorr() {
    static Sudakov_Stats
      pdfchstat(Dipole::iigg,"PDFQGSpl",true,0.0,s_pdfapprox/5);
    static Sudakov_Stats
      pdfkpstat(Dipole::iigg,"PDFGGSpl",true,0.0,s_pdfapprox/20);
    assert(m_step==3); ++m_step;
    //double pwgt=m_sgroup.PlusPDFCorr(m_mfl);
    //double mwgt=m_sgroup.MinusPDFCorr(m_mfl);
    //if(m_mfl[sf::plusfin].Kfcode()==ATOOLS::kf::gluon) {
    //  pdfglustat.Include(pwgt); pdfquastat.Include(mwgt);
    //} else {
    //  pdfglustat.Include(mwgt); pdfquastat.Include(pwgt);
    //}
    //return (pwgt*mwgt/s_pdfapprox);
    return 0.0;
  }



  template<Dipole::Type DT>
  inline const double IISudakov<DT,Radiation::quark>::CalculateRapLimit() {
    assert(m_step==0); ++m_step;
    double term=4*m_sgroup.X2t();
    return (1.0/sqrt(term))*(1.0+sqrt(1.0-term));    //This is Zmax=exp(Ymax).
  }



  template<Dipole::Type DT>
  inline const double IISudakov<DT,Radiation::quark>::GenerateRap() {
    assert(m_step==1); ++m_step;
    double zmaxr=1.0/m_sgroup.Ymax();
    return -std::log(zmaxr+(m_sgroup.Ymax()-zmaxr)*ATOOLS::ran.Get());
  }



  template<Dipole::Type DT>
  inline const double IISudakov<DT,Radiation::quark>::GenerateCorr() {
    static Sudakov_Stats mestat(DT,"MEgsplitWeight",false);
    assert(m_step==2); ++m_step;
    const double& x1=m_sgroup.CurrentSudakovResult().X1;    //x_i_new.
    const double& x3=m_sgroup.CurrentSudakovResult().X3;    //x_i_keep.
    double fac=1.0+sqr((1.0-x1)/(1.0-x1-x3));
    double x2t=m_sgroup.Scale()*m_sgroup.X2t();
    double zmax=m_sgroup.Ymax();    //
    const Multidouble& isrx=m_sgroup.CurrentSudakovResult().Isr;
    assert(isrx.size()==sr::stop);
    //if(m_sgroup.CurrentDipole().SpinCorr()==false); else;
    x2t=
#ifdef IISUDAKOV_K2T_RSCALE
      Sudakov_Calculator::AlphaSCorr(ATOOLS::sqr(isrx[sr::kt])) *
#else
      Sudakov_Calculator::AlphaSCorr(x2t) *
#endif
      (sqrt(x2t)/isrx[sr::mdip])*(zmax-1.0/zmax) *
      //This is the rapidity weight.
      (m_genx2tfac*s_pdfapprox) * fac;
    mestat.Include(x2t);/////////////////////////////////////////////////////
    return x2t;
  }



  //===========================================================================



}    //eo namespace ADICIC





//eof
