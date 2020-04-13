#include "Recola_Virtual.H"

#include "AddOns/Recola/Recola_Interface.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "METOOLS/Loops/Divergence_Array.H"

using namespace Recola;
using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

Ioperator::Ioperator(const Flavour_Vector &fls,const bool& on) :
  m_on(on), m_n(fls.size()), m_mu2(0.)
{
  if (!m_on) return;
  m_Q2i.resize(m_n);
  m_mi.resize(m_n);
  m_m2i.resize(m_n);
  m_gi.resize(m_n);
  m_Ki.resize(m_n);
  m_Gi.resize(m_n);
  for (size_t i(0);i<m_n;++i) m_Gi[i].resize(2);
  m_Qik.resize(m_n);
  m_VSik.resize(m_n);
  m_VNSik.resize(m_n);
  for (size_t i(0);i<m_n;++i) {
    m_Qik[i].resize(m_n);
    m_VSik[i].resize(m_n);
    m_VNSik[i].resize(m_n);
  }
}

Ioperator::~Ioperator() {}

void Ioperator::ComputeConstants(const Flavour_Vector &fls)
{
  for (size_t i(0);i<m_n;++i) {
    m_Q2i[i] = sqr(fls[i].Charge());
    // so far no photons
    if (m_Q2i[i]==0.) continue;
    m_mi[i] = fls[i].Mass();
    m_m2i[i] = sqr(m_mi[i]);
    if (fls[i].IsFermion()) {
      m_gi[i] = m_Q2i[i]*3./2.;
      m_Ki[i] = m_Q2i[i]*(7./2.-sqr(M_PI)/6.);
      if (m_mi[i]==0.) {
        m_Gi[i][0] = 0.;
        m_Gi[i][1] = m_gi[i];
      }
      else {
        msg_Debugging()<<m_m2i[i]<<" "<<m_mu2<<std::endl;
        m_Gi[i][0] = m_Q2i[i]*(0.5*log(m_m2i[i]/m_mu2)-2.);
        m_Gi[i][1] = m_Q2i[i];
      }
    }
    else THROW(not_implemented,"Piece missing.");
    for (size_t k(0);k<m_n;++k) if (i!=k) {
      m_Qik[i][k] = m_Qik[k][i] = fls[i].Charge()*fls[k].Charge();
    }
  }
  if (msg_LevelIsDebugging()) {
    for (size_t i(0);i<m_n;++i) {
      msg_Out()<<"Q2["<<i<<"]="<<m_Q2i[i]<<std::endl;
      msg_Out()<<" G["<<i<<"]="<<m_Gi[i]<<std::endl;
      msg_Out()<<" g["<<i<<"]="<<m_gi[i]<<std::endl;
      msg_Out()<<" K["<<i<<"]="<<m_Ki[i]<<std::endl;
    }
  }
}

double Ioperator::VS(const size_t& i,const size_t& k,
                     const double& sik,const double& q2ik,
                     const size_t& id)
{
  if (m_mi[i]==0. && m_mi[k]==0.) {
    if (id==0 || id==1) return 0.;
    else                return 1.;
  }
  else if (m_mi[i]==0. || m_mi[k]==0.) {
    double mui2(m_mi[i]/q2ik);
    double muk2(m_mi[k]/q2ik);
    double m2(m_m2i[i]+m_m2i[k]);
    double lms=log(m2/sik);
    if      (id==0) return -0.25*sqr(lms)-sqr(M_PI)/12.
                           -0.5*log(sik/q2ik)*(lms+log(m2/q2ik));
    else if (id==1) return 0.5*lms;
    else            return 0.5;
  }
  else {
    if (id==2) return 0.;
    double mui2(m_m2i[i]/q2ik);
    double muk2(m_m2i[k]/q2ik);
    double vik(q2ik/sik*sqrt(Lambda(1.,mui2,muk2)));
    double logrho(log(sqrt((1.-vik)/(1.+vik))));
    if (id==1) return logrho/vik;
    double logrhoi(log(sqrt((1.-vik+2.*mui2/(1.-mui2-muk2))/
                            (1.+vik+2.*mui2/(1.-mui2-muk2)))));
    double logrhok(log(sqrt((1.-vik+2.*muk2/(1.-mui2-muk2))/
                            (1.+vik+2.*muk2/(1.-mui2-muk2)))));
    if (id==0) return (-sqr(logrhoi)-sqr(logrhok)-sqr(M_PI)/6.
                       +logrho*log(q2ik/sik))/vik;
  }
  THROW(fatal_error,"You should not be here.");
  return 0.;
}

double Ioperator::VNS(const size_t& i,const size_t& k,
                      const double& sik,const double& q2ik)
{
  if (m_mi[i]==0. && m_mi[k]==0.) return 0.;
  else {
    if (m_mi[i]==0.) {
      double qik(sqrt(q2ik));
      return 3./2.*(log(sik/q2ik)
                     -2.*log(1.-m_mi[k]/qik)-2.*m_mi[k]/(qik+m_mi[k]))
              +sqr(M_PI)/6.-DiLog(sik/q2ik);
    }
    else if (m_mi[k]==0.) {
      return (3./2.-2.)*log(sik/q2ik)+sqr(M_PI)/6.-DiLog(sik/q2ik)
             -m_m2i[i]/sik*log(m_m2i[i]/q2ik);
    }
    else {
      double qik(sqrt(q2ik));
      double mui2(m_m2i[i]/q2ik);
      double muk2(m_m2i[k]/q2ik);
      double vik(q2ik/sik*sqrt(Lambda(1.,mui2,muk2)));
      double rhoj2(((1.-vik)*sik+2.*m_m2i[i])/((1.+vik)*sik+2.*m_m2i[i]));
      double rhok2(((1.-vik)*sik+2.*m_m2i[k])/((1.+vik)*sik+2.*m_m2i[k]));
      double rho2(rhoj2*rhok2);
      return 3./2.*log(sik/q2ik)
              +(log(rho2)*log(1.+rho2)
                +2.*DiLog(rho2)-DiLog(1.-rhoj2)-DiLog(1.-rhok2)
                -sqr(M_PI)/6.)/vik
              +log(1.-m_mi[k]/qik)
              -2.*log((sqr(qik-m_mi[k])-m_m2i[i])/q2ik)
              -2.*m_m2i[i]/sik*log(m_mi[i]/(qik-m_mi[k]))
              -m_mi[k]/(qik-m_mi[k])
              +2.*m_mi[k]*(2.*m_mi[k]-qik)/sik
              +0.5*sqr(M_PI);
    }
  }
  THROW(fatal_error,"You should not be here.");
  return 0.;
}

double Ioperator::I(const Vec4D_Vector &moms,const double &epsfac)
{
  if (!m_on) return 0.;
  DEBUG_FUNC("");
  double Ie2(0.),Ie1(0.),Ifin(0.);
  for (size_t i(0);i<m_n;++i) {
    for (size_t k(0);k<m_n;++k) if (m_Qik[i][k]) {
      double q2ik((moms[i]+moms[k]).Abs2());
      double sik(q2ik-m_m2i[i]-m_m2i[k]);
      double splfin(0.),sple1(0.),sple2(0.);
      sple2 +=m_Q2i[i]*VS(i,k,sik,q2ik,2);
      sple1 +=m_Q2i[i]*VS(i,k,sik,q2ik,1)+m_Gi[i][1];
      splfin+=m_Q2i[i]*(VS(i,k,sik,q2ik,0)+VNS(i,k,sik,q2ik)-sqr(M_PI)/3.)
              +m_Gi[i][0]+m_gi[i]*(1.+log(m_mu2/sik))+m_Ki[i];

      double fac(m_Qik[i][k]/m_Q2i[i]);
      sple2 *=fac;
      sple1 *=fac;
      splfin*=fac;

      double lsc(log(4.*M_PI*m_mu2/dabs(sik)/epsfac));
      Ie2 +=sple2;
      Ie1 +=sple1+sple2*lsc;
      Ifin+=splfin+sple1*lsc+0.5*sple2*sqr(lsc);
      msg_Debugging()<<"  i="<<i<<", k="<<k<<std::endl
                     <<"  I_e2 = "<<sple2<<std::endl
                     <<"  I_e1 = "<<sple1+sple2*lsc<<std::endl
                     <<"  Ifin = "<<splfin+sple1*lsc+0.5*sple2*sqr(lsc)<<std::endl;
    }
  }
  msg_Debugging()<<"I_e2 = "<<Ie2<<std::endl
                 <<"I_e1 = "<<Ie1<<std::endl
                 <<"Ifin = "<<Ifin<<std::endl;
  return -Ifin;
}

Recola_Virtual::Recola_Virtual(const Process_Info& pi,
			      const Flavour_Vector& flavs,
			      unsigned int recola_id) :
  Virtual_ME2_Base(pi, flavs), m_recola_id(recola_id),
  m_ewscheme(ToType<int>(rpa->gen.Variable("EW_SCHEME"))),
  m_voqcd(pi.m_maxcpl[0]),
  m_boqcd(pi.m_maxcpl[0]-(pi.m_fi.m_nloqcdtype==nlo_type::loop)),
  p_iop(Recola_Interface::UseIopInEWapprox()?new Ioperator(flavs,true):NULL)
{
  m_providespoles=false;
  m_fixedIRscale=true;

  m_IRscale=Recola_Interface::IRScale();
  m_UVscale=Recola_Interface::UVScale();
  if (p_iop) {
    p_iop->SetIRScale2(sqr(m_IRscale));
    p_iop->ComputeConstants(flavs);
  }

  // init associated contribs
  size_t n(0);
  if (pi.m_fi.m_asscontribs&asscontrib::EW) {
    ++n;
    if (pi.m_fi.m_asscontribs&asscontrib::LO1) {
      ++n;
      if (pi.m_fi.m_asscontribs&asscontrib::LO2) {
        ++n;
        if (pi.m_fi.m_asscontribs&asscontrib::LO3) {
          ++n;
        }
      }
    }
  }
  m_asscontribs.resize(n);

  if (m_asscontribs.size()>0 && m_voqcd!=m_boqcd+1)
    THROW(fatal_error,"Associated contribs only implemented for NLO QCD.");
}

Recola_Virtual::~Recola_Virtual()
{
  if (p_iop) { delete p_iop; p_iop=NULL; }
}

void Recola_Virtual::Calc(const Vec4D_Vector& momenta)
{
  DEBUG_FUNC("");
  if (!Recola_Interface::checkProcGeneration()) {
    // use different ew scheme in different interface functions
    switch (m_ewscheme) {
    case 3:
      use_gfermi_scheme_and_set_alpha_rcl(AlphaQED());
      break;
    case 2:
      //use_alphaz_scheme_and_set_alpha_rcl(AlphaQED());
      use_alphaz_scheme_rcl(AlphaQED());
      break;
    case 1:
      //use_alpha0_scheme_and_set_alpha_rcl(AlphaQED());
      use_alpha0_scheme_rcl(AlphaQED());
      break;
    default:
      msg_Out()<<"The EW scheme "<<m_ewscheme<<" is not available with the "
               <<"Sherpa+Recola interface. Valid options are:"<<std::endl
               <<"  1. alpha_QED(0)"<<std::endl
               <<"  2. alpha_QED(M_Z)"<<std::endl
               <<"  3. Gmu"<<std::endl;
      THROW(not_implemented,"Electroweak scheme not implemented in Sherpa+Recola.");
    }


    int nlight=0;
    set_mu_ir_rcl(m_IRscale);
    set_mu_uv_rcl(m_UVscale);
    int fixed(Recola_Interface::GetFixedFlav());
    if (Recola_Interface::GetDefaultFlav()==0) fixed=5;

    double alpha_mat;
    int default_flavscheme(fixed);
    if (default_flavscheme==16) default_flavscheme=-1;
    if (fixed>0 && fixed<10) {
      nlight=fixed;
    }
    else {
      if (default_flavscheme>10) {
        nlight=Recola_Interface::PDFnf(m_mur2,default_flavscheme-10);
      }
      if (default_flavscheme==-1)
        nlight=-1;
      if (default_flavscheme==-2 || default_flavscheme==0) {
        if (Flavour(kf_c).Mass()!=0)
          nlight=3;
        else if (Flavour(kf_b).Mass()!=0)
          nlight=4;
        else if (Flavour(kf_t).Mass()!=0)
          nlight=5;
        else {
          msg_Out()<<"WARNING: 6 light flavours detected.\n";
          nlight=6;
        }
      }
    }
    if (nlight==0){
      msg_Error()<<METHOD<<"(): Cannot determine number of flavours\n";
    }
    if (nlight>6){
      msg_Error()<<METHOD<<"(): Too many light flavours: "<<nlight<<"\n   Max is 6\n";
    }

    // Recola_Interface::SetDefaultFlav(nlight);

    double default_alphaQCD=Recola_Interface::GetDefaultAlphaQCD();
    double default_scale=Recola_Interface::GetDefaultScale();
    set_alphas_rcl(default_alphaQCD,sqrt(default_scale),nlight);
    msg_Debugging()<<"use AlphaQCD="<<AlphaQCD()
                   <<",  sqrt(m_mur2)="<<sqrt(m_mur2)<<std::endl;

    msg_Out()<<"processes in Recola are being generated..."<<std::endl;
    generate_processes_rcl();
    Recola_Interface::setProcGenerationTrue();
    msg_Out()<<"process generation in Recola completed..."<<std::endl;
  }

  // calculate
  m_res*=0.; m_born=0.;
  for (size_t i(0);i<m_asscontribs.size();++i) m_asscontribs[i]=0.;

  MyTiming* timing;
  if (msg_LevelIsDebugging()) {
    timing = new MyTiming();
    timing->Start();
    msg_Out()<<"virtual correction, ";
  }
  set_alphas_rcl(AlphaQCD(),sqrt(m_mur2),Recola_Interface::GetDefaultFlav());
  Recola_Interface::EvaluateProcess(m_recola_id, momenta, m_voqcd, m_boqcd,
                                    m_res,m_born,m_asscontribs);

  if (msg_LevelIsDebugging()) {
    timing->Stop();
    PRINT_INFO(momenta[2][0]<<" "<<m_flavs<<" = "<<m_res.Finite()<<" / "<<m_born
                            <<", user="<<timing->UserTime()
                            <<", real="<<timing->RealTime()
                            <<", sys="<<timing->SystemTime());
    for (size_t i(0);i<m_asscontribs.size();++i) {
      PRINT_INFO(momenta[2][0]<<" "<<m_flavs<<" = "<<m_asscontribs[i]
                              <<", user="<<timing->UserTime()
                              <<", real="<<timing->RealTime()
                              <<", sys="<<timing->SystemTime());
    }
  }

  // factor which by Sherpa convention has to be divided out at this stage
  if (m_born==0) m_born=1.;
  double factor=m_born*AlphaQCD()/2.0/M_PI;
  m_res.Finite()/=factor;
  m_res.IR()/=factor;
  m_res.IR2()/=factor;
  if (p_iop && m_asscontribs.size()>1) {
    double I(AlphaQED()/2.0/M_PI*m_born
             *p_iop->I(momenta,Eps_Scheme_Factor(momenta)));
    msg_Debugging()<<"V="<<m_asscontribs[0]<<", I="<<I<<std::endl;
    m_asscontribs[0]+=I;
  }
  for (size_t i(0);i<m_asscontribs.size();++i) m_asscontribs[i]/=factor;
  msg_Debugging()<<"V/B="<<m_res.Finite()<<std::endl;
  for (size_t i(0);i<m_asscontribs.size();++i)
    msg_Debugging()<<"ASS/B="<<m_asscontribs[i]<<std::endl;
}



DECLARE_VIRTUALME2_GETTER(Recola_Virtual,"Recola_Virtual")
Virtual_ME2_Base *ATOOLS::Getter<Virtual_ME2_Base,Process_Info,Recola_Virtual>::
operator()(const Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  if (pi.m_loopgenerator!="Recola") return NULL;

  if (pi.m_fi.m_nloqcdtype!=nlo_type::loop) return NULL;

  int procIndex=Recola_Interface::RegisterProcess(pi, amptype::treeloop);
  if (procIndex>0) {
    Flavour_Vector flavs = pi.ExtractFlavours();
    return new Recola_Virtual(pi, flavs, procIndex);
  }
  else {
    return NULL;
  }
}