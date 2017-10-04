#include "PHASIC++/Process/Massive_Kernels.H"

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Default_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/NLO_Types.H"

#include "MODEL/Main/Model_Base.H"

using namespace ATOOLS;
using namespace PHASIC;

//for equations see hep-ph/0201036

Massive_Kernels::Massive_Kernels(ATOOLS::sbt::subtype st,
                                 const size_t &nf, const size_t nmf) :
  m_stype(st),
  m_nf(nf), m_nmf(nmf), m_NC(3.), m_CA(3.), m_CF(4./3.), m_TR(0.5),
  m_TRbyCA(m_TR/m_CA), m_CFbyCA(m_CF/m_CA), m_TRbyCF(m_TR/m_CF),
  m_g1t(0.), m_g2t(0.), m_g3t(0.), m_K1t(0.), m_K2t(0.), m_K3t(0.),
  m_beta0qcd(0.), m_beta0qed(0.),
  m_alpha_ff(1.), m_alpha_fi(1.), m_alpha_if(1.), m_alpha_ii(1.),
  m_kappa(2./3.),
  m_logaff(0.), m_logafi(0.), m_logaif(0.), m_logaii(0.),
  m_VNS(0.), m_gKterm(0.), m_aterm(0.), m_Vsubmode(1)
{
  m_subtype = 0;
  for (size_t i(1);i<=m_nf+m_nmf;i++) {
    Flavour flav((kf_code)(i));
    if (flav.IsMassive()) m_massflav.push_back(flav.Mass());
  }
  p_VS[0]=p_VS[1]=p_VS[2]=0.;
  p_Gammat[0]=p_Gammat[1]=0.;
}

void Massive_Kernels::SetNC(const double &nc)
{
  DEBUG_FUNC(m_stype<<((m_stype==sbt::qcd)?(" Nc="+ToString(nc)):""));
  // all massless charged particles in the model,
  //
  double sumQ2(0.),sumQ2quark(0.),sumQ2lepton(0.);
  Flavour_Vector flavs(MODEL::s_model->IncludedFlavours());
  for (size_t i(0);i<flavs.size();++i) {
    if (flavs[i].IsAnti()) continue;
    if (!flavs[i].IsMassive() && flavs[i].Charge()) {
      if (flavs[i].IsQuark()) {
        sumQ2quark += nc*sqr(flavs[i].Charge());
      }
      else if (flavs[i].IsLepton()) {
        sumQ2lepton += sqr(flavs[i].Charge());
      }
      else {
        THROW(fatal_error,"Massless charged flavour encountered. Cannot cope. Abort.");
      }
    }
  }
  sumQ2=sumQ2quark+sumQ2lepton;
  switch (m_stype) {
  case sbt::qcd:
    // T_j^2 for all quarks the same -> define m_g1, m_K1
    m_NC  = nc;
    m_TR  = 0.5;
    m_CA  = 2.*m_TR*m_NC;
    m_CF  = m_TR*(m_NC*m_NC-1.0)/m_NC;
    m_TRbyCA = m_TR/m_CA;
    m_CFbyCA = m_CF/m_CA;
    m_TRbyCF = m_TR/m_CF;
    m_g1t = 1.5;
    m_g2t = 11./6.-2./3.*m_TRbyCA*m_nf;
    m_g3t = 2.;
    m_K1t = 7./2.-sqr(M_PI)/6.;
    m_K2t = 67./18.-sqr(M_PI)/6.-10./9.*m_TRbyCA*m_nf;
    m_K3t = 4.-sqr(M_PI)/6.;
    break;
  case sbt::qed:
    // Q_j^2 different for all fermions -> do not define m_TR, m_CA, m_CF
    m_TR  = 1.;
    m_CA  = 1.;
    m_CF  = 1.;
    m_TRbyCA = 1.;
    m_CFbyCA = 1.;
    m_TRbyCF = nc; // this should be 1 for initial state leptons
    m_g1t = 1.5;
    m_g2t = -2./3.*sumQ2;
    m_g3t = 2.;
    m_K1t = 7./2.-sqr(M_PI)/6.;
    m_K2t = -10./9.*sumQ2;;
    m_K3t = 4.-sqr(M_PI)/6.;
  case sbt::none:
    break;
  }
  m_beta0qcd = 11./6.*m_CA-2./3.*m_TR*m_nf;
  m_beta0qed = -2./3.*sumQ2;
  if (msg_LevelIsDebugging()) {
    if (m_stype==sbt::qcd) msg_Out()<<"Nc="<<m_NC<<", CA="<<m_CA
                                    <<", CF="<<m_CF<<", TR="<<m_TR<<std::endl;
    if (m_stype==sbt::qed) msg_Out()<<"sum_q="<<sumQ2quark
                                    <<", sum_l="<<sumQ2lepton
                                    <<", sum_{q,l}="<<sumQ2<<std::endl;
    msg_Out()<<"g1t="<<m_g1t<<", g2t="<<m_g2t<<", g3t="<<m_g3t<<std::endl;
    msg_Out()<<"K1t="<<m_K1t<<", K2t="<<m_K2t<<", K3t="<<m_K3t<<std::endl;
  }
}

void Massive_Kernels::SetSubType(const int subtype)
{
  m_subtype=subtype;
}

void Massive_Kernels::SetAlpha(double aff, double afi, double aif, double aii)
{ 
  m_alpha_ff=aff; m_alpha_fi=afi; m_alpha_if=aif; m_alpha_ii=aii; 
  m_logaff=log(aff); m_logafi=log(afi); m_logaif=log(aif); m_logaii=log(aii);
}


// I operator

void Massive_Kernels::CalcVS(ist::itype type,double s,double mj,double mk)
// V^S(eps) as in (6.20)
{
  DEBUG_FUNC("type="<<type<<", s="<<s<<", mj="<<mj<<", mk="<<mk);
  p_VS[0]=p_VS[1]=p_VS[2]=0.;
  // does not exist for photons
  if (m_stype==sbt::qed && type==ist::g) return;
  // both massive
  if (mj>0.&&mk>0.) {
    double mj2=sqr(mj);
    double mk2=sqr(mk);
    double Q2=s+mj2+mk2;
    double vjk=sqrt(Lambda(Q2,mj2,mk2))/s;
    double lrhoj=log(sqrt(((1.-vjk)*s+2.*mj2)/((1.+vjk)*s+2.*mj2)));
    double lrhok=log(sqrt(((1.-vjk)*s+2.*mk2)/((1.+vjk)*s+2.*mk2)));
    p_VS[2]=0.;
    p_VS[1]=(lrhoj+lrhok)/vjk;
    p_VS[0]=(-sqr(lrhoj)-sqr(lrhok)-sqr(M_PI)/6.+(lrhoj+lrhok)*log(Q2/s))/vjk;
  }
  // one massive, one massless
  else if (mj>0.||mk>0.) {
    double m2=sqr(mj+mk);
    double Q2=s+m2;
    double lms=log(m2/s);
    p_VS[2]=0.5;
    p_VS[1]=0.5*lms;
    p_VS[0]=-0.25*sqr(lms)-sqr(M_PI)/12.-0.5*log(s/Q2)*(lms+log(m2/Q2));
  }
  // both massless
  else {
    p_VS[2]=1.;
    p_VS[1]=0.;
    p_VS[0]=0.;
  }
  msg_Debugging()<<"VS[0]="<<p_VS[0]
                 <<" ,  VS[1]="<<p_VS[1]
                 <<" ,  VS[2]="<<p_VS[2]<<std::endl;
}

void Massive_Kernels::CalcVNS(ist::itype type,double s,double mj,double mk,
                              bool ini)
// V^NS
{
  DEBUG_FUNC(type<<": s="<<s<<", mj="<<mj<<", mk="<<mk<<", ini="<<ini);
  m_VNS=0.;
  // does not exist for photons
  if (m_stype==sbt::qed && type==ist::g) return;
  switch (type) {
  case ist::q:
  case ist::Q:
  case ist::sG:
    CalcVNSq(s,mj,mk);
    break;
  case ist::g:
    CalcVNSg(s,mk,ini);
    break;
  case ist::sQ:
    CalcVNSs(s,mj,mk);
    break;
  case ist::V:
    if      (m_Vsubmode==0) CalcVNSs(s,mj,mk);
    else if (m_Vsubmode==1) CalcVNSq(s,mj,mk);
    break;
  default:
    THROW(fatal_error,"Unknown splitting type.");
    break;
  }
  msg_Debugging()<<"VNS="<<m_VNS<<std::endl;
}

void Massive_Kernels::CalcVNSq(double s,double mj,double mk)
// V^NS_q as defined in (6.21)-(6.23)
{
  if (mj==0.&&mk==0.) return;
  else if (mj==0.) {
    double Q2=s+sqr(mj)+sqr(mk);
    double Q=sqrt(Q2);
    m_VNS = m_g1t*(log(s/Q2)-2.*log(1.-mk/Q)-2.*mk/(Q+mk))
            +sqr(M_PI)/6.-DiLog(s/Q2);
    if (m_subtype==1) m_VNS += .25-(Q-mk)*(Q+3.*mk)/(4.*sqr(Q+mk));
  }
  else if (mk==0.) {
    double mj2=sqr(mj);
    double mk2=sqr(mk);
    double Q2=s+mj2+mk2;
    m_VNS = (m_g1t-2.)*log(s/Q2)+sqr(M_PI)/6.-DiLog(s/Q2)-mj2/s*log(mj2/Q2);
    if (m_subtype==1) m_VNS += .25*(1.-(Q2+mj2)/(Q2-mj2))
			-mj2/(Q2-mj2)*(1.+(2.*Q2+mj2)/(Q2-mj2)*log(mj2/Q2)/2.);
  }
  else {
    double mj2=sqr(mj);
    double mk2=sqr(mk);
    double Q2=s+mj2+mk2;
    double Q=sqrt(Q2);
    double vjk=sqrt(Lambda(Q2,mj2,mk2))/s;
    double rhoj2=((1.-vjk)*s+2.*mj2)/((1.+vjk)*s+2.*mj2);
    double rhok2=((1.-vjk)*s+2.*mk2)/((1.+vjk)*s+2.*mk2);
    double rho2=rhoj2*rhok2;
    m_VNS = m_g1t*log(s/Q2)
            +(log(rho2)*log(1.+rho2)
              +2.*DiLog(rho2)-DiLog(1.-rhoj2)-DiLog(1.-rhok2)
              -sqr(M_PI)/6.)/vjk
            +log(1.-mk/Q)-2.*log((sqr(Q-mk)-mj2)/Q2)-2.*mj2/s*log(mj/(Q-mk))
            -mk/(Q-mk)+2.*mk*(2.*mk-Q)/s+0.5*sqr(M_PI);
    if (m_subtype==1) {
      double muj=mj/Q, muk=mk/Q, muj2=mj2/Q2, muk2=mk2/Q2;
      m_VNS += .25+(muj2*((muj2*(.25+1./(1.-muk)-3.*log(muj/(1.-muk))))/(1.-muj2-muk2)
	-2.*log(muj/(1.-muk))))/(1.-muj2-muk2)-((1.+4.*muj2+2.*(1.-muk)*muk-muk2)*sqr(1.-muk))/(4.*sqr(1.-muj2-muk2));
    }
  }
}

void Massive_Kernels::CalcVNSg(double s,double mk,bool ini)
// V^NS_g as defined in Eqs.(6.24) and (6.26);
// Q_aux-terms canceled with Gamma_g
{
  size_t nfjk=0;
  if (!ini) 
    for (size_t i=0;i<m_nmf;i++)
      if (4.*m_massflav[i]*(m_massflav[i]+mk)<s) nfjk++;
  if (mk==0.) {
    for (size_t i=0;i<nfjk;i++) {
      double rho1=sqrt(1.-4.*sqr(m_massflav[i])/s);
      m_VNS+=log(0.5+0.5*rho1)-rho1*(1.+sqr(rho1)/3.)
             -0.5*log(sqr(m_massflav[i])/s);
      if (m_subtype==1) {
	double muj2=sqr(m_massflav[i])/s;
	m_VNS+=3./4.*(12.*muj2*log((1.-rho1)/(1.+rho1))*(-4.-17.*muj2+4.*sqr(muj2))+
	  (-1.-154.*muj2-152.*sqr(muj2)+64.*pow(muj2,3))*rho1)/(18.*sqr(1.-2.*muj2));
      }
    }
    m_VNS*=4./3.*m_TRbyCA;
  }
  else {
    bool simplev=ini||IsEqual(m_kappa,2./3.);
    double Q2=s+sqr(mk);
    double Q=sqrt(Q2);
    double muk=mk/Q, muk2=mk*mk/Q2;
    m_VNS=m_g2t*(log(s/Q2)-2.*log(1.-mk/Q)-2.*mk/(Q+mk))
          +sqr(M_PI)/6.-DiLog(s/Q2);
    if (!simplev)
      m_VNS+=(m_kappa-2./3.)*sqr(mk)/s
             *((2.*m_nf*m_TRbyCA-1.)*log(2.*mk/(Q+mk)));
    if (m_subtype==1) {
      m_VNS += (muk2*log((2.*muk)/(1.+muk)))/(3.*(1.-muk2))+muk2/(2.*pow(1.+muk,3))-pow(muk/(1.+muk),3)/18.;
      m_VNS += m_nf*m_TR/m_CA*(-muk2*(9.-muk)/(9.*pow(1.+muk,3))-2.*muk2/(3.*(1.-muk2))*log(2.*muk/(1.+muk)));
    }
    double nfc=0.;
    for (size_t i=0;i<nfjk;i++) {
      double rho1=sqrt(1.-4.*sqr(m_massflav[i])/sqr(Q-mk));
      double rho2=sqrt(1.-4.*sqr(m_massflav[i])/s);
      nfc+=4./3.*(log(1.-mk/Q)+mk*rho1*rho1*rho1/(Q+mk)+log(0.5+0.5*rho1)
                  -rho1*(1.+sqr(rho1)/3.)-0.5*log(sqr(m_massflav[i])/Q2));
      if (!simplev) {
        nfc+=(m_kappa-2./3.)*2.*sqr(mk)/s
             *(rho2*rho2*rho2*log((rho2-rho1)/(rho2+rho1))
               -log((1.-rho1)/(1.+rho1))-8.*rho1*sqr(m_massflav[i])/s);
      }
      if (m_subtype==1) {
	double muj2=sqr(m_massflav[i])/(Q2+sqr(m_massflav[i]));
	nfc += (-2.*muk2*((-8.*muj2*rho1)/(1-muk2)-log((1-rho1)/(1+rho1))+log((-rho1+rho2)/(rho1+rho2))*pow(rho2,3)))/(3.*(1.-muk2))
	  -log((1.-rho1)/(1.+rho1))*((3.-2.*muj2)/(3.*(1-muk2))-(muj2*(3.+2.*muj2-(2.*muj2*(9.+5.*muj2))/(1.-2.*muj2-muk2)))/(3.*(1.-2.*muj2-muk2)*muj2))
	  -(rho1*(65./6.-15*muj2+(8*muj2)/(3.*(1-muk))-4*muk+muk2/6.-(2.*muj2*(36.-148.*muj2-26.*muk+77.*muj2*muk+59.*sqr(muj2)))/(3.*(1.-2.*muj2-muk2)*muj2)
		  +((-49.*muj2+10.*(1.-muk)+39.*muj2*muk+127.*sqr(muj2)-77.*muk*sqr(muj2)-30.*pow(muj2,3))*sqr(2./(1.-2.*muj2-muk2)))/3.))/(3.*(1.-muk2));
      }
    }
    m_VNS+=m_TRbyCA*nfc;
  }
}

void Massive_Kernels::CalcVNSs(double s,double mj,double mk)
// V^NS_s as defined in (C.9)-(C.10)
{
  if (mk==0.) {
    double mj2=sqr(mj);
    double mk2=sqr(mk);
    double Q2=s+mj2+mk2;
    m_VNS = (m_g3t-2.)*log(s/Q2)+sqr(M_PI)/6.-DiLog(s/Q2);
  }
  else {
    double mj2=sqr(mj);
    double mk2=sqr(mk);
    double Q2=s+mj2+mk2;
    double Q=sqrt(Q2);
    double vjk=sqrt(Lambda(Q2,mj2,mk2))/s;
    double rhoj2=((1.-vjk)*s+2.*mj2)/((1.+vjk)*s+2.*mj2);
    double rhok2=((1.-vjk)*s+2.*mk2)/((1.+vjk)*s+2.*mk2);
    double rho2=rhoj2*rhok2;
    m_VNS = m_g3t*log(s/Q2)
            +(log(rho2)*log(1.+rho2)
              +2.*DiLog(rho2)-DiLog(1.-rhoj2)-DiLog(1.-rhok2)
              -sqr(M_PI)/6.)/vjk
            -2.*log((sqr(Q-mk)-mj2)/Q2)
            +4.*mk*(mk-Q)/s+0.5*sqr(M_PI);
  }
}

void Massive_Kernels::CalcGamma(ist::itype type, double mu2, double s, double m)
// Gamma(eps) as in Eq.(6.27)-(6.29), (C.13)
{
  DEBUG_FUNC(type<<", mu2="<<mu2<<", m="<<m);
  p_Gammat[0]=p_Gammat[1]=0.;
  switch (type) {
  case ist::q:
    p_Gammat[0]=0.;
    p_Gammat[1]=m_g1t;
    if (m_subtype==1) p_Gammat[0]+=-1./4.;
    break;
  case ist::g:
    p_Gammat[0]=0.;
    p_Gammat[1]=m_g2t;
    if (m_subtype==1) p_Gammat[0]+=1./36.-m_nf*m_TRbyCA/18.;
    break;
  case ist::Q:
  case ist::sG:
    p_Gammat[0]=0.5*log(sqr(m)/mu2)-2.;
    p_Gammat[1]=1.;
    if (m_subtype==1) p_Gammat[0]+=-1./4.;
    break;
  case ist::sQ:
    p_Gammat[0]=log(sqr(m)/mu2)-2.;
    p_Gammat[1]=1.;
    break;
  case ist::V:
    if      (m_Vsubmode==0) {
      p_Gammat[0]=log(sqr(m)/mu2)-2.;
      p_Gammat[1]=1.;
    }
    else if (m_Vsubmode==1) {
      p_Gammat[0]=0.5*log(sqr(m)/mu2)-2.;
      p_Gammat[1]=1.;
    }
    break;
  default:
    THROW(fatal_error,"Unknown splitting type.");
  }
  double lmus=log(mu2/s);
  p_Gammat[0]-=lmus*p_Gammat[1];

  msg_Debugging()<<"Gammat[0]="<<p_Gammat[0]
                 <<" ,  Gammat[1]="<<p_Gammat[1]<<std::endl;
}

void Massive_Kernels::CalcgKterm(ATOOLS::ist::itype type, double mu2, double s,
                                 double mj, bool mode)
// \gamma+K as in Eqs.(5.91), (6.17), (C.12)
{
  DEBUG_FUNC(type<<", mu2="<<mu2<<", s="<<s<<", mj="<<mj<<", mode"<<mode);
  m_gKterm=0.;
  double lmus=log(mu2/s);
  switch (type) {
  case ist::q:
  case ist::Q:
  case ist::sG:
    m_gKterm+=m_g1t*(1.+lmus)+m_K1t;
    if (IsZero(mj)) m_gKterm-=(mode?0.5:0.);
    break;
  case ist::g:
    m_gKterm+=m_g2t*(1.+lmus)+m_K2t-(mode?1./6.:0.);
    break;
  case ist::sQ:
    m_gKterm+=m_g3t*(1.+lmus)+m_K3t;
    break;
  case ist::V:
    if      (m_Vsubmode==0) m_gKterm+=m_g3t*(1.+lmus)+m_K3t;
    else if (m_Vsubmode==1) m_gKterm+=m_g1t*(1.+lmus)+m_K1t;
    break;
  default:
    THROW(fatal_error,"Unknown splitting type.");
  }
  msg_Debugging()<<"gKterm="<<m_gKterm<<std::endl;
}

/// alpha terms.
void Massive_Kernels::CalcAterms(ist::itype type, double mu2, double s,
                                 double mj, double mk, bool inij, bool inik)
{
  m_aterm=0.;
  if (inij || inik || m_alpha_ff==1.) return;
  DEBUG_FUNC(type<<": mu2="<<mu2<<", s="<<s<<", mj="<<mj<<", mk="<<mk);
  switch (type) {
  case ist::q:
  case ist::Q:
  case ist::sG:
    CalcAq(mu2,s,mj,mk);
    break;
  case ist::g:
    CalcAg(mu2,s,mk);
    break;
  case ist::sQ:
    CalcAs(mu2,s,mj,mk);
    break;
  case ist::V:
    if      (m_Vsubmode==0) CalcAs(mu2,s,mj,mk);
    else if (m_Vsubmode==1) CalcAq(mu2,s,mj,mk);
    break;
  default:
    THROW(fatal_error,"Unknown type.");
    break;
  }
  msg_Debugging()<<"aterm="<<m_aterm<<std::endl;
}

void Massive_Kernels::CalcAq(double mu2, double s,double mj,double mk)
{
  double Q2 = s+sqr(mj)+sqr(mk);
  double muj2 = mj*mj/Q2;
  double muk2 = mk*mk/Q2;
  double muk = sqrt(muk2);
  if (mj==0.) {
    if (mk==0.) {
      m_aterm+= -sqr(m_logaff)-m_g1t*(m_logaff+1.-m_alpha_ff);
    }
    else {
      double yp = (1.-muk)/(1.+muk);
      double xp = yp*(1.-m_alpha_ff)
                  + sqrt((1.-m_alpha_ff)*(1.-m_alpha_ff*yp*yp));
      m_aterm += sqr(log((1.-yp*yp+2.*xp*yp)/(1.+yp-xp)/(1.-yp+xp)))
                 -2.*sqr(log((1.+yp-xp)/(1.+yp)))
                 +4.*(log((1.+yp)/2.)*log((1.-yp+xp)/(1.-yp))
                      +log((1.+yp)/(2.*yp))
                         *log((1.-yp*yp+2.*xp*yp)/(1.-yp*yp))
                      +DiLog((1.-yp)/(1.+yp))
                      -DiLog((1.-yp*yp+2.*xp*yp)/sqr(1.+yp))
                      +DiLog((1.-yp+xp)/2.)
                      -DiLog((1.-yp)/2.));
      m_aterm += -1.5*(m_logaff+yp*(1.-m_alpha_ff));
    }
  }
  else {
    if (mk==0.) {
      double muq2 = muj2;
      double lmq2 = log(muj2);
      m_aterm+= 2.*(-m_logaff*(1.+lmq2)
                    -DiLog((muq2 -1.)/muq2)
                    +DiLog((m_alpha_ff*(muq2 -1.))/muq2));
      m_aterm+= 0.5*(-(1.-m_alpha_ff)*(3.*m_alpha_ff*(1.-muq2)+2.*muq2)
                                    /(m_alpha_ff+(1.-m_alpha_ff)*muq2)
                     +(1.+muq2)/(1.-muq2)*log(m_alpha_ff+(1.-m_alpha_ff)*muq2));
    }
    else {
      double d   = 0.5*(1.-muj2-muk2);
      double yp  = 0.5*(sqr(1.-muk)-muj2)/d;
      double ap  = m_alpha_ff*yp;
      double vjk = 0.5*sqrt(Lambda(1.,muj2,muk2))/d;
      double a   = muk/d;
      double b   = (1.-muk)/d;
      double c   = a*b*d;
      double x   = yp-ap+sqrt((1.-m_alpha_ff)*(1.-ap*yp-muj2*sqr(a)));
      double xp  = yp+vjk;
      double xm  = yp-vjk;
      m_aterm += 1.5*(1.+ap)
                 + 1./(1.-muk)
                 - (2.-2.*muj2-muk)/d
                 + 0.5*(1.-ap)*muj2/(muj2+2.*ap*d)
                 - 2.*m_logaff
                 + 0.5*(muj2+d)/d*log((muj2+2.*ap*d)/sqr(1.-muk));
                 ///eq A20
      m_aterm += 2.*(- DiLog((a+x)/(a+xp)) + DiLog(a/(a+xp))
                     + DiLog((xm-x)/(xm+a)) - DiLog(xm/(xm+a))
                     + DiLog((c+x)/(c+xp)) - DiLog(c/(c+xp))
                     - DiLog((xm-x)/(xm+c)) + DiLog(xm/(xm+c))
                     - DiLog((b-x)/(b-xm)) + DiLog(b/(b-xm))
                     + DiLog((xp-x)/(xp-b)) - DiLog(xp/(xp-b))
                     + DiLog((b-x)/(b+a)) - DiLog(b/(b+a))
                     - DiLog((c+x)/(c-a)) + DiLog(c/(c-a))
                     + log(c+x)*log((a-c)*(xp-x)/((a+x)*(c+xp)))
                     - log(c)*log((a-c)*xp/(a*(c+xp)))
                     + log(b-x)*log((a+x)*(xm-b)/((a+b)*(xm-x)))
                     - log(b)*log(a*(xm-b)/(xm*(a+b)))
                     - log((a+x)*(b-xp))*log(xp-x)
                     + log(a*(b-xp))*log(xp)
                     + log(d)*log((a+x)*xp*xm/(a*(xp-x)*(xm-x)))
                     + log((xm-x)/xm)*log((c+xm)/(a+xm))
                     + 0.5*log((a+x)/a)*log(a*(a+x)*(a+xp)*(a+xp)))/vjk;
    }
  }
}

void Massive_Kernels::CalcAg(double mu2, double s,double mk)
{
  double Q2 = s+sqr(mk);
  double muk2 = mk*mk/Q2;
  double muk = sqrt(muk2);
  if (IsZero(mk)) {
    for (size_t i=0;i<m_nmf;i++) if (4.*m_massflav[i]*(m_massflav[i]+mk)<s) {
      double muj2 = m_massflav[i]*m_massflav[i]/Q2;
      double muj4 = 4.*muj2*muj2;
      double a = sqrt(sqr(m_alpha_ff*(1.-2.*muj2)) -muj4);
      double b = sqrt(1.-4.*muj2);
      m_aterm -= m_TR*2./3.
                 *(2.*a/(2.*(m_alpha_ff-1.)*muj2-m_alpha_ff)+a
                   +(2.*muj2 -1.)*(-log(-2.*(a+m_alpha_ff*(2.*muj2-1.)))
                                   +2.*atan(2.*muj2/a)
                                   +log(-2.*(2.*muj2+b-1.))
                                   -2.*atan(2.*muj2/b))
                   +b); ///ref 0 eq A9
    }
    m_aterm+=-sqr(m_logaff)-m_g2t*(m_logaff+1.-m_alpha_ff);
  }
  else {
    double yp = 1. - 2.*muk*(1.-muk)/(1.-muk2);
    double yl = 1. - muk2 + m_alpha_ff*sqr(1. -muk)
                - (1.-muk)*sqrt(sqr(m_alpha_ff*(1.-muk)) + sqr(1.+muk)
                               - 2.*m_alpha_ff*(1.+muk2));
    m_aterm += -(11.*sqr(2.-2.*muk-yl)/((1.-muk2)*(2.-yl))
                 -22.*log(4.*sqr((1.-muk))*muk)
                 +24.*(sqr(log(2.))-sqr(log(1.+muk)))
                 -2.*((11.-15.*muk2)*log((2.-yl)/(2.*muk))
                      +4.*muk2*(log((sqr(2.-yl)-4.*muk2*(1.-yl))/
                                    (8.*(1.-muk)*muk2))))
                   /(1.-muk2)
                 +22.*log((2. - 2.*muk2 - yl)*yl)
                 -12.*(4.*log(1.-yl/2.)*log(-(yl/(-1.+muk2)))
                       -sqr(log(-(yl/(-1.+muk2))))
                       +sqr(log((-2.*(-2. + 2.*muk2 + yl))
                                /((-1. + muk2)*(-2. + yl))))
                       +2.*log(-(yl/(-1.+muk2)))
                         *(log((-2.*(-2.+2.*muk2+yl))/((-1.+muk2)*(-2.+yl)))
                           -2.*log(1.+yl/(-2.+2.*muk2))))
                 +48.*DiLog(1. - muk)
                 -48.*DiLog(1./(1. + muk))
                 -48.*DiLog(yl/2.)
                 +48.*DiLog(yl/(2.-2.*muk2)))/12.;
    m_aterm += m_TRbyCA*m_nf*(2./3.*((1.-muk-m_alpha_ff*yp*(1.+muk))/(1.+muk)
                                      +log(m_alpha_ff*yp*(1.+muk)/(1.-muk)))
                               +(m_kappa-2./3.)*2.*muk2/(1.-muk2)
                                 *log((1.-m_alpha_ff*yp)*(1.+muk)/(2.*muk)));
    /// ref 42 eq 21.
    for (size_t i=0;i<m_nmf;i++) if (4.*m_massflav[i]*(m_massflav[i]+mk)<s) {
      double muj2 = m_massflav[i]*m_massflav[i]/Q2;
      double muj4 = 4.*muj2*muj2;
      double a = sqrt(1.-muk2);
      double c = -1. + 2.*muj2 +muk2;
      double c2 = c*c;
      double b = sqrt(c2*yp*yp - muj4);
      double d = sqrt(sqr(m_alpha_ff*c*yp) - muj4);
      double e = sqrt(c2-muj4);
      yp = 1. - 2.*muk*(1.-muk)/(1.-2.*muj2 - muk2);
      m_aterm -= m_TR*(b*d*a
                       *((-8.*muj2*muj2 + 2.*c2 +2.*c+4.*muj2)
                           *log((m_alpha_ff*c2*yp-d*e -muj4)
                                /(-b*e+c2*yp-muj4))
                         +2.*a*(c2+c-muj4+2.*muj2)
                           *log((1.-yp)/(1.-m_alpha_ff*yp))
                         +(-3.*c2 +4.*c*muj2-2.*c)*sqrt(2.0*muj2-c)
                           *log((m_alpha_ff*c*yp+d)/(b+c*yp))
                         +2.*sqrt(2.0*muj2-c)*(c2-2.*(c+1.)*muj2+muj4)
                           *(atan(2.*muj2/d)-atan(2.*muj2/b)))
                       +(c*c2*yp*sqrt(2.0*muj2-c)
                           *(m_alpha_ff*m_alpha_ff*b*yp-2.*m_alpha_ff*b
                             -d*(yp-2.))
                         +4.*c*muj2*(b*(m_alpha_ff*yp-1.)-d*yp +d)
                         +muj4*(b-d)))
                 /(3.*c*pow(2.0*muj2-c,3.0/2.0)*b*d);
      ///ref 0 eq A6 reformulated.
    }
  }
}

void Massive_Kernels::CalcAs(double mu2, double s,double mj,double mk)
{
  DEBUG_FUNC(mu2<<" "<<s<<" "<<mj<<" "<<mk);
  double Q2 = s+sqr(mj)+sqr(mk);
  double muj2 = mj*mj/Q2;
  double muk2 = mk*mk/Q2;
  double muk = sqrt(muk2);
  if (mk==0.) {
    double muq2 = muj2;
    double lmq2 = log(muj2);
    m_aterm += 2.*(-m_logaff*(1.+lmq2)
                   -DiLog((muq2 -1.)/muq2)
                   +DiLog((m_alpha_ff*(muq2 -1.))/muq2)
                   -(1.-m_alpha_ff)*(m_alpha_ff*(1.-muq2)+muq2)
                                   /(m_alpha_ff+(1.-m_alpha_ff)*muq2));
  }
  else {
    double mjmk1 = 1.-muj2 - muk2;
    double yp = (sqr(1.-muk)-muj2)/(1.-muj2-muk2);
    double ap = m_alpha_ff*yp;
    double d = 0.5*(1.-muj2-muk2);
    double a = muk/d;
    double b = (1.-muk)/d;
    double c = a*b*d;
    double vjk = sqrt(Lambda(1.,muj2,muk2))/(1.-muj2-muk2);
    double x   = yp-ap+sqrt((1.-m_alpha_ff)*(1.-ap*yp-muj2*sqr(a)));
    double xp  = yp+vjk;
    double xm  = yp-vjk;
    m_aterm += 1.5*(1.+ap)
               +1./(1.-muk)
               -(2.-2.*muj2-muk)/d
               +0.5*(1.-ap)*muj2/(muj2+2.*ap*d)
               -2.*m_logaff
               +(-2.*(1.-ap)*d+2.*muk*(1.-muk))*(ap*(1.-muk)*d+muj2)
                /(2.*(1.-muk)*d*(2.*ap*d+muj2));
    msg_Debugging()<<m_aterm<<" "
                   <<1.5*(1.+ap)+1./(1.-muk)
                     -2.*(2.-2.*muj2-muk)/mjmk1
                     +(1.-ap)*muj2/(2.*(muj2+ap*mjmk1))
                     -2.*log(ap*mjmk1/(sqr(1.-muk)-muj2))
                     -(-ap*mjmk1-muj2+sqr(muk-1.))*(ap*(1.-muk)*mjmk1+2.*muj2)
                                          /(2.*(1.-muk)*mjmk1*(ap*mjmk1+muj2))
                   <<std::endl;
    ///eq A20
    m_aterm += 2.*(-DiLog((a+x)/(a+xp))+DiLog(a/(a+xp))
                   +DiLog((xp-x)/(xp-b))-DiLog(xp/(xp-b))
                   +DiLog((c+x)/(c+xp))-DiLog(c/(c+xp))
                   +DiLog((xm-x)/(xm+a))-DiLog(xm/(xm+a))
                   -DiLog((b-x)/(b-xm))+DiLog(b/(b-xm))
                   -DiLog((xm-x)/(xm+c))+DiLog(xm/(xm+c))
                   +DiLog((b-x)/(b+a))-DiLog(b/(b+a))
                   -DiLog((c+x)/(c-a))+DiLog(c/(c-a))
                   +log(c+x)*log((a-c)*(xp-x)/((a+x)*(c+xp)))
                   -log(c)*log((a-c)*xp/(a*(c+xp)))
                   +log(b-x)*log((a+x)*(xm-b)/((a+b)*(xm-x)))
                   -log(b)*log(a*(xm-b)/(xm*(a+b)))
                   -log((a+x)*(b-xp))*log(xp-x)
                   +log(a*(b-xp))*log(xp)
                   +log(d)*log((a+x)*xp*xm/(a*(xp-x)*(xm-x)))
                   +log((xm-x)/xm)*log((c+xm)/(a+xm))
                   +0.5*log((a+x)/a)*log(a*(a+x)*(a+xp)*(a+xp)))/vjk;
    msg_Debugging()<<m_aterm<<std::endl;
  }
}

void Massive_Kernels::Calculate(ist::itype type, double mu2, double s,
                                double mj, double mk,
                                bool inij, bool inik, bool mode)
{
  DEBUG_FUNC(type);
  CalcVS(type,s,mj,mk);
  CalcVNS(type,s,mj,mk,inij);
  CalcGamma(type,mu2,s,mj);
  CalcgKterm(type,mu2,s,mj,mode);
  CalcAterms(type,mu2,s,mj,mk,inij,inik);
}

double Massive_Kernels::I_Fin()
{
  return p_VS[0]+m_VNS-sqr(M_PI)/3.+p_Gammat[0]+m_gKterm+m_aterm;
}

double Massive_Kernels::I_E1()
{
  return p_VS[1]+p_Gammat[1];
}

double Massive_Kernels::I_E2()
{
  return p_VS[2];
}


// KP operator

double Massive_Kernels::Kb1(int type,double x)
{
  if (m_stype==sbt::qed && type==4) return 0.;
  switch(type) {
  case 1:
  case 4:
    return 2./(1.-x)*log((1.-x)/x);
  }
  return 0.;
}

double Massive_Kernels::Kb2(int type)
{
  if (m_stype==sbt::qed && type==4) return -8./3.*m_g2t;
  switch(type) {
  case 1:
    return -(m_g1t+m_K1t-5./6.*sqr(M_PI));
  case 4:
    return -(m_g2t+m_K2t-5./6.*sqr(M_PI));
  }
  return 0.;
}

double Massive_Kernels::Kb3(int type,double x)
{
  if (m_stype==sbt::qed && type==4) return 0.;
  double c=0.;
  if (m_subtype==1) c=-log(1./x+sqr(1./x-1.))/x;
  switch(type) {
  case 1:
    return (-(1.+x)*log((1.-x)/x)+(1.-x));
  case 2:
    return m_CFbyCA*((1.+sqr(1.-x))/x*log((1.-x)/x)+x+2.*c);
  case 3:
    return m_TRbyCF*((x*x+sqr(1.-x))*log((1.-x)/x)+2.*x*(1.-x));
  case 4:
    return 2.*(((1.-x)/x-1.+x*(1.-x))*log((1.-x)/x)+c);
  }
  return 0.;
}

double Massive_Kernels::Kb4(int type,double x)
{
  if (m_stype==sbt::qed && type==4) return 0.;
  switch(type) {
  case 1:
  case 4: {
    double l1x(log(1.-x));
    return 2.*(-0.5*sqr(l1x)+l1x*log(x)+DiLog(x));
  }
  }
  return 0.;
}

double Massive_Kernels::KFS1(int type,double x)
{
  if (m_stype==sbt::qed && type==4) return 0.;
  return 0.;
}

double Massive_Kernels::KFS2(int type)
{
  if (m_stype==sbt::qed && type==4) return 0.;
  return 0.;
}

double Massive_Kernels::KFS3(int type,double x)
{
  if (m_stype==sbt::qed && type==4) return 0.;
  return 0.;
}

double Massive_Kernels::KFS4(int type,double x)
{
  if (m_stype==sbt::qed && type==4) return 0.;
  return 0.;
}

double Massive_Kernels::Kc1(int type,int typej,double muja2,double x)
{
  if (m_stype==sbt::qed && type==4) return 0.;
  return 0.;
}

double Massive_Kernels::Kc2(int type,int typej,double muja2)
{
  if (m_stype==sbt::qed && type==4) return 0.;
  return 0.;
}

double Massive_Kernels::Kc3(int type,int typej,double muja2,double x)
{
  if (m_stype==sbt::qed && type==4) return 0.;
  return 0.;
}

double Massive_Kernels::Kc4(int type,int typej,double muja2,double x)
{
  if (m_stype==sbt::qed && type==4) return 0.;
  return 0.;
}

double Massive_Kernels::KbM1(int type,double muak2,double x)
{
  if (m_stype==sbt::qed && type==4) return 0.;
  return 0.;
}

double Massive_Kernels::KbM2(int type,double muak2)
{
  if (m_stype==sbt::qed && type==4) return 0.;
  return 0.;
}

double Massive_Kernels::KbM3(int type,double muak2,double x)
{
  if (m_stype==sbt::qed && type==4) return 0.;
  return 0.;
}

double Massive_Kernels::KbM4(int type,double muak2,double x)
{
  if (m_stype==sbt::qed && type==4) return 0.;
  return 0.;
}

double Massive_Kernels::Kt1(int type,double x)
{
  if (m_stype==sbt::qed && type==4) return 0.;
  switch(type) {
  case 1:
    return 2./(1.-x)*log(1.-x);
  case 4:
    return 2./(1.-x)*log(1.-x);
  }
  return 0.;
}

double Massive_Kernels::Kt2(int type)
{
  if (m_stype==sbt::qed && type==4) return 0.;
  switch(type) {
  case 1:
  case 4:
    return -sqr(M_PI)/3.;
  }
  return 0.;
}

double Massive_Kernels::Kt3(int type,double x)
{
  if (m_stype==sbt::qed && type==4) return 0.;
  double ax=0., c=0.;
  if (m_subtype==1) c=-log(1./x+sqr(1./x-1.))/x;
  if (m_alpha_ii<(1.-x)) ax=log(m_alpha_ii/(1.-x));
  switch(type) {
  case 1:
    ax*=(1.+x*x)/(1.-x);
    if (m_subtype==1) ax+=-(1.-x);
    return -(1.+x)*log(1.-x)+ax;
  case 2:
    ax*=(1.+sqr(1.-x))/x;
    if (m_subtype==1) ax+=(1.-x)+4.*log(x)/x/(1.+x)-2.*c;
    if (m_subtype==2) ax+=2.*log(x)/x;
    return m_CFbyCA*((1.+sqr(1.-x))/x*log(1.-x)+ax);
  case 3:
    ax*=(1.-2.*x*(1.-x));
    if (m_subtype==1) ax+=-(1.-x)*(1.-3.*x);
    return m_TRbyCF*((x*x+sqr(1.-x))*log(1.-x)+ax);
  case 4:
    ax*=x/(1.-x)+(1.-x)/x+x*(1.-x);
    if (m_subtype==1) ax+=0.5*(1.-x*(4.-3.*x)+4.*log(x)/x/(1.+x)-2.*c);
    if (m_subtype==2) ax+=log(x)/x;
    return 2.*((1.-x)/x-1.+x*(1.-x))*log(1.-x)+2.*ax;
  }
  return 0.;
}

double Massive_Kernels::Kt4(int type,double x)
{
  if (m_stype==sbt::qed && type==4) return 0.;
  switch(type) {
  case 1:
    return -sqr(log(1.-x));
  case 4:
    return -sqr(log(1.-x));
  }
  return 0.;
}

double Massive_Kernels::P1(int type,double x)
{
  if (m_stype==sbt::qed && type==4) return 0.;
  switch(type) {
  case 1:
    return (1.+x*x)/(1.-x);
  case 4:
    return 2./(1.-x);
  }
  return 0.;
}

double Massive_Kernels::P2(int type)
{
  if (type==4) return m_g2t;
  return 0.;
}

double Massive_Kernels::P3(int type,double x)
{
  if (m_stype==sbt::qed && type==4) return 0.;
  switch(type) {
  case 2:
    return m_CFbyCA*(1.+sqr(1.-x))/x;
  case 3:
    return m_TRbyCF*(x*x+sqr(1.-x));
  case 4:
    return 2.*((1.-x)/x-1.+x*(1.-x));
  }
  return 0.;
}

double Massive_Kernels::P4(int type,double x)
{
  if (m_stype==sbt::qed && type==4) return 0.;
  switch(type) {
  case 1:
    return -x-0.5*x*x-2.*log(1.-x);
  case 4:
    return -2.*log(1.-x);
  }
  return 0.;
}

//muq2=m_j^2/s_ja

double Massive_Kernels::t1(int type,int spin,double muq2,double x)
// g(x)
{
  if (m_stype==sbt::qed && type==4) return 0.;
  if (type==2||type==3) return 0.;
  double aterm(0.);
  if (m_alpha_fi<1.) aterm = -at1(type, spin, muq2, x);
  switch (spin) {
  case 0:
    return 2./(1.-x)*(1.+log((1.-x+muq2)/(1.-x)))+aterm;
  case 1:
    return 2./(1.-x)*(1.+log((1.-x+muq2)/(1.-x)))
           -0.5*(1.-x)/sqr(1.-x+muq2)+aterm;
  case 2:
    return m_g2t/(1.-x)+aterm;
  }
  return 0.;
}

double Massive_Kernels::t2(int type,int spin,double muq2)
// h; in case of gluon muq2 must be s_ja!!
{
  if (m_stype==sbt::qed && type==4) return 0.;
  if (type==2||type==3) return 0.;
  double aterm(0.);
  if (m_alpha_fi<1.) aterm = -at2(type, spin, muq2);
  switch (spin) {
  case 0: {
    double mx=muq2/(1.+muq2);
    if (IsZero(muq2)) return m_g1t+aterm;
    double res(0.);
    if (type==4)
      res=m_g1t
          -2.*m_g2t*(log(sqrt(1.+muq2)-sqrt(muq2))+1./(sqrt(1./mx)+1.))
          -muq2*log(mx)-0.5*mx + aterm;
    else
      res=m_g1t*(1.-2.*(log(sqrt(1.+muq2)-sqrt(muq2))+1./(sqrt(1./mx)+1.)))
          -muq2*log(mx)-0.5*mx + aterm;
    return res+(muq2*log(mx)+0.5*mx)-(m_g1t-m_g3t);
  }
  case 1: {
    double mx=muq2/(1.+muq2);
    if (IsZero(muq2)) return m_g1t+aterm;
    if (type==4)
      return m_g1t
             -2.*m_g2t*(log(sqrt(1.+muq2)-sqrt(muq2))+1./(sqrt(1./mx)+1.))
             -muq2*log(mx)-0.5*mx + aterm;
    return m_g1t*(1.-2.*(log(sqrt(1.+muq2)-sqrt(muq2))+1./(sqrt(1./mx)+1.)))
           -muq2*log(mx)-0.5*mx + aterm;
  }
  case 2: {
    double mgs=0.;
    for (size_t i=0;i<m_nmf;i++) {
      double xp=1.-sqr(2.*m_massflav[i])/muq2;
      if (xp>0.) mgs+=pow(xp,1.5); // TODO: fix for QED
    }
    return (m_g2t-m_TRbyCA*2./3.*mgs) + aterm;
  }
  }
  return aterm;
}

double Massive_Kernels::t3(int type,int spin,double muq2,double x)
// k(x)
{
  if (m_stype==sbt::qed && type==4) return 0.;
  double aterm(0.);
  if (m_alpha_fi<1. || m_alpha_if<1.) aterm = -at3(type, spin, muq2, x);
  if (IsZero(muq2)) return aterm;
  if (spin==2) return aterm;
  double mx=log((1.-x)/(1.-x+muq2));
  switch (type) {
  case 1:
    return (1.+x)*mx + aterm;
  case 2:
    return -m_CFbyCA*((1.+sqr(1.-x))*mx-2.*muq2*log((1.-x)/muq2+1.))/x + aterm;
  case 3:
    return -m_TRbyCF*(x*x+sqr(1.-x))*mx + aterm;
  case 4:
    return -2.*((1./x-2.+x*(1.-x))*mx-muq2/x*log((1.-x)/muq2+1.)) + aterm;
  }
  return aterm;
}

double Massive_Kernels::t4(int type,int spin,double muq2,double x)
// G(x)
{
  if (m_stype==sbt::qed && type==4) return 0.;
  if (type==2||type==3) return 0.;
  double aterm(0.);
  if (m_alpha_fi<1.) aterm = -at4( type, spin, muq2, x);
  if (IsZero(muq2)) {
    switch (spin) {
    case 1:
      return -m_g1t*log(1.-x) + aterm;
    case 2:
      return -m_g2t*log(1.-x) + aterm;
    }
  }
  double y=1.-x;
  double lny=log(y);
  switch (spin) {
  case 0:
    return sqr(lny)+2.*(DiLog(-y/muq2)-DiLog(-1./muq2)-log(muq2)*lny)-2.*lny;
  case 1:
    return sqr(lny)+2.*(DiLog(-y/muq2)-DiLog(-1./muq2)-log(muq2)*lny)
      +0.5*(muq2*x/((1.+muq2)*(y+muq2))-log((1.+muq2)/(y+muq2)))-2.*lny + aterm;
  case 2:
    return -m_g2t*lny + aterm;
  }
  return aterm;
}

double Massive_Kernels::t5(int type,double x,double xp)
// g^{xp}(x)
{
  if (m_stype==sbt::qed && type==4) return 0.;
  if (type==2||type==3) return 0.;
  if (x>xp) return 0.;
  x=1.-x;
  xp=1.-xp;
  return -2./3.*m_TRbyCA*(x+0.5*xp)/sqr(x)*sqrt(1.-xp/x);
}

double Massive_Kernels::t6(int type,double xp)
// h^{xp}
{
  if (m_stype==sbt::qed && type==4) return 0.;
  if (type==2||type==3) return 0.;
  double sxp=sqrt(xp);
  return -2./3.*m_TRbyCA*(log((1.-sxp)/(1.+sxp))+sxp/3.*(6.-xp));
}

double Massive_Kernels::t7(int type,double x,double xp)
// G^{xp}(x)
{
  if (m_stype==sbt::qed && type==4) return 0.;
  if (type==2||type==3) return 0.;
  if (x>xp) x=xp;
  return -2./3.*m_TRbyCA*((sqrt(1.-(1.-xp)/(1.-x))
                            *(5.+(1.-xp)/(1.-x))-sqrt(xp)*(6.-xp))/3.
                          -log((1.+xp)/2.-x+sqrt((1.-x)*(xp-x)))
                          +log((1.+xp)/2.+sqrt(xp)));
}

double Massive_Kernels::at1(int type,int spin,double muq2,double x)
// g(x)
{
  if (m_stype==sbt::qed && type==4) return 0.;
  if (type==2||type==3) return 0.;
  double res(0.);
  if (spin == 0) {
    if (x<1.-m_alpha_fi) res -= 2.*(log((1.+muq2)/muq2) - 1.)/(1.-x);
  }
  else if (spin == 1) {
    if (x<1.-m_alpha_fi) {
      if (IsZero(muq2)) res = 2.*log(1.-x)/(1.-x) + 1.5/(1.-x);
      else              res -= 2.*(log((1.+muq2)/muq2) - 1.)/(1.-x);
    }
  }
  else if (spin == 2) {
    /// final state is gluon - sum over all possible splittings
    // TODO: fix for QED
    if (x<1.-m_alpha_fi) res -= m_TRbyCA*m_nf*(2./3./(1.-x));
    if (x<1.-m_alpha_fi) res -= (-2./(1.-x)*log(1.-x)-11./6./(1.-x));
    size_t nfjk=0;
    for (size_t i=0;i<m_nmf;i++)
      if (4.*m_massflav[i]*(m_massflav[i])<muq2) nfjk++;
    for (size_t i=0; i<nfjk;i++) {
      double muQ2 = (m_massflav[i]*m_massflav[i])/muq2;
      if (x<1.-m_alpha_fi) res += 2./3.*((1.-x+2.*muQ2)/sqr(1.-x))
                                       *sqrt(1.-4.*muQ2/(1.-x));
    }
  }
  return res;
}

double Massive_Kernels::at2(int type,int spin,double muq2)
// h; in case of gluon muq2 must be s_ja!!
{
  if (m_stype==sbt::qed && type==4) return 0.;
  if (type==2||type==3) return 0.;
  double res(0.);
  if (spin == 0) { /// final state is scalar
    res += 2.*log(m_alpha_fi)*(log((1.+muq2)/muq2)-1.);
  }
  else if (spin == 1) {  /// final state is quark
    if (IsZero(muq2)) res += (-1.5*m_logafi - sqr(m_logafi));
    else              res += 2.*m_logafi*(log((1.+muq2)/muq2)-1.);
  }
  else if (spin == 2) {
    /// final state is gluon - sum over all possible splittings
    // TODO: fix for QED
    res += m_TRbyCA*m_nf*(m_logafi*2./3.);
    res -= /*2.**/(sqr(m_logafi) + 11./6.*m_logafi);
    double muQ2, a, b, c;
    c=sqrt(m_alpha_fi);
    size_t nfjk=0;
    for (size_t i=0;i<m_nmf;i++)
      if (4.*m_massflav[i]*(m_massflav[i])<muq2) nfjk++;
    for (size_t i=0; i<nfjk;i++) {
      muQ2 = (m_massflav[i]*m_massflav[i])/muq2;
      a = sqrt(1.-4.*muQ2);
      b = sqrt(m_alpha_fi-4.*muQ2);
      res += 2./9.*(-4.*muQ2*(b/m_alpha_fi/c+4./a)-5.*b/c
                    -sqr(4.*muQ2)/a +5./a +6.*log(b+c)-6.*log(a+1.));
    }
  }
  return res;
}


double Massive_Kernels::at3(int type,int spin,double muq2,double x)
// k(x)
{
  if (m_stype==sbt::qed && type==4) return 0.;
  if (spin==2) muq2 = muq2*x;
  else muq2 = muq2/x;
  double res(0.);
  if (type!=2 && type!=3) {
    if (spin == 0) { /// final state is scalar
      if (x<1.-m_alpha_fi) res -= /*(1.-x)/(2.*sqr(1.-x+muq2*x))+*/
                                  2./(1.-x)*log((2.-x+muq2*x)*muq2
                                                /(1.+muq2)/(1.-x+muq2*x));
    }
    else if (spin == 1) {  /// final state is quark
      if (x<1.-m_alpha_fi) {
        if (IsZero(muq2)) res += -2.*log(2.-x)/(1.-x);
        else              res -= (1.-x)/(2.*sqr(1.-x+muq2*x))
                                 +2./(1.-x)*log((2.-x+muq2*x)*muq2/(1.+muq2)
                                                  /(1.-x+muq2*x));
      }
    }
    else if (spin == 2) {
      /// final state is gluon - sum over all possible splittings
      /// (only one non-zero here is g->gg)
      if (x<1.-m_alpha_fi) res -= (2.*log(2.-x)/(1.-x));
    }
  }
  double zp=(1.-x)/(1.-x+muq2*x);
  if (spin==2) zp = 1.;
  switch (type) {
  case 1:
    if (zp>m_alpha_if)
      res-=2./(1.-x)*log(zp*(1.-x+m_alpha_if)/m_alpha_if/(1.-x+zp))
           -(1.+x)*log(zp/m_alpha_if);
    break;
  case 2:
    if (zp>m_alpha_if) {
      if (zp!=1.) res += -m_CFbyCA*((1.+sqr(1.-x))/x*(log(zp/m_alpha_if))
                                     +2.*muq2*log((1.-zp)/(1.-m_alpha_if)));
      else        res += -m_CFbyCA*(2.-2.*x+x*x)/x*log(zp/m_alpha_if);
    }
    break;
  case 3:
    if (zp>m_alpha_if) res += -m_TRbyCA*(x*x+sqr(1.-x))*log(zp/m_alpha_if);
    break;
  case 4:
    if (zp>m_alpha_if) {
      if (zp!=1.) res += -2.*((1./x-2.+x*(1.-x))*log(zp/m_alpha_if)
                              +muq2*log((1.-zp)/(1.-m_alpha_if))
                              -log(m_alpha_if*(1.-x+zp)/(zp*(1.-x+m_alpha_if)))
                                /(1.-x));
      else        res += -2.*((1./x-2.+x*(1.-x))*log(zp/m_alpha_if)
                              -log(m_alpha_if*(1.-x+zp)/(zp*(1.-x+m_alpha_if)))
                                /(1.-x));
    }
    break;
  }
  return res;
}

double Massive_Kernels::at4(int type,int spin,double muq2,double x)
// G(x)
{
  if (m_stype==sbt::qed && type==4) return 0.;
  if (type==2||type==3) return 0.;
  double res(0.);
  if (spin == 0) { /// final state is scalar
    if (x>1.-m_alpha_fi) res -= - 2.*(log((1.+muq2)/muq2) - 1.)*m_logafi;
    else                 res -= - 2.*(log((1.+muq2)/muq2) - 1.)*log(1.-x);
  }
  else if (spin == 1) {  ///final state is quark
    if (IsZero(muq2)) {
      if (x>1.-m_alpha_fi) res -= sqr(m_logafi) + 1.5*m_logafi;
      else                 res -= sqr(log(1.-x)) + 1.5*log(1.-x);
    }
    else {
      if (x>1.-m_alpha_fi) res -= - 2.*(log((1.+muq2)/muq2) - 1.)*m_logafi;
      else                 res -= - 2.*(log((1.+muq2)/muq2) - 1.)*log(1.-x);
    }
  }
  else if (spin == 2) {
    /// final state is gluon - sum over all possible splittings
    // TODO: fix for QED
    if (x>1.-m_alpha_fi) res -= (-m_TRbyCA*m_nf*2./3.+11./6.)*m_logafi
                                +sqr(m_logafi);
    else                 res -= (-m_TRbyCA*m_nf*2./3.+11./6.)*log(1.-x)
                                +sqr(log(1.-x));
    size_t nfjk=0;
    for (size_t i=0;i<m_nmf;i++)
      if (4.*m_massflav[i]*(m_massflav[i])<muq2) nfjk++;
    for (size_t i=0; i<nfjk;i++) {
      double muQ2 = (m_massflav[i]*m_massflav[i])/muq2;
      double rt = sqrt(1.-4.*muQ2);
      if (x>1.-m_alpha_fi) {
        double rta = sqrt(1.-4.*muQ2/m_alpha_fi);
        res += 2./3.*log(2.*m_alpha_fi*(rta +1.)-4.*muQ2)
               -2./9./m_alpha_fi*rta*(4.*muQ2 +5.*m_alpha_fi)
               +2./9.*(4.*rt*muQ2+5.*rt-3.*log(-2.*muQ2+rt+1.)-log(8.));
      }
      else {
        double rta = sqrt(1.-4.*muQ2/(1.-x));
        res += 2./3.*log(2.*(1.-x)*(rta +1.)-4.*muQ2)
               -2./9./(1.-x)*rta*(4.*muQ2 +5.*(1.-x))
               +2./9.*(4.*rt*muQ2+5.*rt-3.*log(-2.*muQ2+rt+1.)-log(8.));
      }
    }
  }
  return res;
}

