#include "Recola_Virtual.H"

#include "AddOns/Recola/Recola_Interface.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Message.H"

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

namespace Recola {

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
    m_modebackup(m_mode), m_ismapped(false),
    p_iop(Recola_Interface::UseIopInEWapprox()!=0?new Ioperator(flavs,true):NULL)
  {
    m_procmap[m_recola_id]=pi;
    Settings& s = Settings::GetMainSettings();
    m_providespoles=false;
    m_fixedIRscale=true;

    m_IRscale=s["RECOLA_IR_SCALE"].Get<double>();
    m_UVscale=s["RECOLA_UV_SCALE"].Get<double>();
    m_modebackup=m_mode=Recola_Interface::s_vmode;
    m_voqcd = pi.m_maxcpl[0];
    m_boqcd = pi.m_maxcpl[0]-pi.m_fi.m_nlocpl[0];

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
    
    m_mode=m_modebackup;
    if (!Recola_Interface::checkProcGeneration()){
        std::cout<<"process generation started..."<<std::endl;
        Recola_Interface::GenerateProcesses(AlphaQED(),AlphaQCD(),
                                            m_IRscale,m_UVscale,m_mur2);
    }

    m_res*=0.; m_born=0.;
    for (size_t i(0);i<m_asscontribs.size();++i) m_asscontribs[i]=0.;

    MyTiming* timing;
    if (msg_LevelIsDebugging()) {
      timing = new MyTiming();
      timing->Start();
    }

    double aqcd=AlphaQCD(); 
    int flav=Recola_Interface::GetDefaultFlav();
    set_alphas_rcl(aqcd,sqrt(m_mur2),flav);
    Recola_Interface::EvaluateLoop(m_recola_id, momenta, m_born, m_res, m_asscontribs);
    if(Recola_Interface::UseIopInEWapprox()==2 && p_iop)
    {
      double I(AlphaQED()/2.0/M_PI*m_born
             *p_iop->I(momenta,Eps_Scheme_Factor(momenta)));
      msg_Debugging()<<"V="<<m_res.Finite()<<", I="<<I<<std::endl;
      m_res.Finite()+=I;
    }

    if (msg_LevelIsDebugging()) {
      timing->Stop();
      PRINT_INFO(momenta[2][0]<<" "<<m_flavs<<" user="<<timing->UserTime()
     <<" real="<<timing->RealTime()<<" sys="<<timing->SystemTime());
    }
 
    
    double coupling(1.);
    if (m_stype&sbt::qcd) coupling=AlphaQCD();
    else if (m_stype&sbt::qed) coupling=AlphaQED();
    else THROW(fatal_error,"Unknown coupling.");
    

    // if Born vanishes, do not divide by it, reset mode for this event
    if(!(m_mode&1) && m_born==0.) {
      m_mode|=1;
      msg_Tracking()<<METHOD<<"(): switch to mode 1, Born vanishes"<<std::endl;
    }
    double factor=((m_mode&1)?1.:m_born)*coupling/2.0/M_PI;
    msg_Debugging()<<"cpl="<<coupling/2.0/M_PI<<std::endl;
    // factor which by Sherpa convention has to be divided out at this stage
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


  bool Recola_Virtual::IsMappableTo(const PHASIC::Process_Info& pi){
    return false;
  }

}

using namespace Recola;

DECLARE_VIRTUALME2_GETTER(Recola_Virtual,"Recola_Virtual")
Virtual_ME2_Base *ATOOLS::Getter<Virtual_ME2_Base,Process_Info,Recola_Virtual>::
operator()(const Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  if (pi.m_loopgenerator!="Recola") return NULL;

  if (pi.m_fi.m_nlotype!=nlo_type::loop) return NULL;

  int procIndex=Recola_Interface::RegisterProcess(pi, 11);
  
  if (procIndex>0) {
    Flavour_Vector flavs = pi.ExtractFlavours();
    return new Recola_Virtual(pi, flavs, procIndex);
  }
  else {
    return NULL;
  }
}
