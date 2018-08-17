#include "PHASIC++/Process/KP_Terms.H"

#include "ATOOLS/Math/Random.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/Default_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"

using namespace PHASIC;
using namespace ATOOLS;

KP_Terms::KP_Terms(Process_Base *const proc,const sbt::subtype st,
                   const std::vector<size_t>& partonlist):
  m_stype(st), m_itype(cs_itype::K|cs_itype::P),
  m_kcontrib(cs_kcontrib::Kb|cs_kcontrib::KFS|cs_kcontrib::t|cs_kcontrib::Kt),
  p_proc(proc), p_nlomc(NULL), p_kernel(NULL),
  p_cpl(NULL), m_flavs(p_proc->Flavours()),
  m_massive(true), m_cemode(false), m_cpldef(0.), m_NC(3.),
  m_Vsubmode(1), m_collVFF(1), m_facscheme(0),
  m_sa(false), m_sb(false),
  m_typea(m_flavs[0].IntSpin()), m_typeb(m_flavs[1].IntSpin()),
  m_plist(partonlist)
{
  DEBUG_FUNC("");
  m_subtype=ToType<int>(rpa->gen.Variable("NLO_SUBTRACTION_SCHEME"));
  Default_Reader reader;
  reader.SetInputPath(rpa->GetPath());
  reader.SetInputFile(rpa->gen.Variable("ME_DATA_FILE"));

  const size_t nf(Flavour(kf_quark).Size()/2);
  int nfgs = ToType<int>(rpa->gen.Variable("DIPOLE_NF_GSPLIT", ToString(nf)));
  
  if (nfgs<nf) THROW(fatal_error,"Number of flavours in g->qq splitting ("
                                 +ToString(nfgs)
                                 +") smaller than number of light flavours ("
                                 +ToString(nf)+").");
  const size_t nmf(nfgs-nf);
  p_kernel=new Massive_Kernels(st,nf,nmf);
  p_kernel->SetSubType(m_subtype);

  m_kcontrib=ToType<cs_kcontrib::type>
             (reader.Get<std::string>("KP_KCONTRIB","BSGT"));
  msg_Tracking()<<"Set K-Term contribution to "<<m_kcontrib<<" .\n";

  m_cemode=reader.Get("KP_CHECK_ENERGY",0);
  msg_Tracking()<<"Set KP-term energy check mode "<<m_cemode<<" .\n";
  m_negativepdf=reader.Get("KP_ACCEPT_NEGATIVE_PDF",1);
  msg_Tracking()<<"Set KP-term accepts negative PDF "<<m_negativepdf<<" .\n";

  m_facscheme=reader.Get("KP_FACTORISATION_SCHEME",0);
  std::string fs("");
  switch (m_facscheme) {
  case 0:
    fs="MSbar";
    break;
  case 1:
    fs="DIS";
    THROW(not_implemented,"DIS factorisation scheme not implemented yet.");
    break;
  default:
    THROW(fatal_error,"Unknown factorisation scheme.");
    break;
  }
  msg_Tracking()<<"Set KP-term factorisation scheme to "<<m_facscheme<<".\n";

  m_NC = reader.Get("N_COLOR",3.);
  if (m_NC!=3.) msg_Out()<<"Set N_color="<<m_NC<<"."<<std::endl;
  p_kernel->SetNC(m_NC);
  SetColourFactors();

  double dalpha(reader.GetValue<double>("DIPOLE_ALPHA",1.));
  double dalpha_ff(reader.GetValue<double>("DIPOLE_ALPHA_FF",dalpha));
  double dalpha_fi(reader.GetValue<double>("DIPOLE_ALPHA_FI",dalpha));
  double dalpha_if(reader.GetValue<double>("DIPOLE_ALPHA_IF",dalpha));
  double dalpha_ii(reader.GetValue<double>("DIPOLE_ALPHA_II",dalpha));
  double kappa(reader.GetValue<double>("DIPOLE_KAPPA",2./3.));
  msg_Tracking()<<"Set FF dipole alpha="<<dalpha_ff<<" . "<<std::endl;
  msg_Tracking()<<"Set FI dipole alpha="<<dalpha_fi<<" . "<<std::endl;
  msg_Tracking()<<"Set IF dipole alpha="<<dalpha_if<<" . "<<std::endl;
  msg_Tracking()<<"Set II dipole alpha="<<dalpha_ii<<" . "<<std::endl;
  msg_Tracking()<<"Set dipole kappa="<<kappa<<" . "<<std::endl;
  SetAlpha(dalpha_ff,dalpha_fi,dalpha_if,dalpha_ii);
  SetKappa(kappa);

  m_Vsubmode=reader.Get("DIPOLE_V_SUBTRACTION_MODE",1);
  msg_Tracking()<<"Use "<<(m_Vsubmode?"fermionic":"scalar")
                <<" V->VP subtraction term."<<std::endl;
  if (!(m_Vsubmode==0 || m_Vsubmode==1)) THROW(fatal_error,"Unknown mode.");
  p_kernel->SetVSubtractionMode(m_Vsubmode);

  int cvffdef(1);
  if (m_stype==sbt::qed) cvffdef=0;
  m_collVFF=reader.Get("DIPOLE_COLLINEAR_VFF_SPLITTINGS",cvffdef);
  msg_Tracking()<<"Switch collinear VFF splittings "<<(m_collVFF?"on":"off")
                <<"."<<std::endl;
  p_kernel->SetCollinearVFFSplitting(m_collVFF);

  SetMassive();
  m_sa=((m_stype==sbt::qcd && m_flavs[0].Strong()) ||
        (m_stype==sbt::qed && (m_flavs[0].Charge() || m_flavs[0].IsPhoton())));
  m_sb=((m_stype==sbt::qcd && m_flavs[1].Strong()) ||
        (m_stype==sbt::qed && (m_flavs[1].Charge() || m_flavs[1].IsPhoton())));


  for (int i=0;i<8;i++) m_kpca[i]=0.;
  for (int i=0;i<8;i++) m_kpcb[i]=0.;
}

KP_Terms::~KP_Terms()
{
  if (p_kernel) { delete p_kernel; p_kernel=NULL; }
}

void KP_Terms::SetColourFactors()
{
  DEBUG_FUNC(m_stype);
  std::vector<double> cpls(m_flavs.size(),0.);
  for (size_t i(0);i<m_flavs.size();++i) {
    if (m_stype==sbt::qcd) {
      if      (abs(m_flavs[i].StrongCharge())==0) cpls[i]=0.;
      else if (abs(m_flavs[i].StrongCharge())==3) cpls[i]=p_kernel->CF();
      else if (abs(m_flavs[i].StrongCharge())==8) cpls[i]=p_kernel->CA();
      else THROW(fatal_error,"Unknown colour charge.");
    }
    else if (m_stype==sbt::qed) {
      if      (m_flavs[i].Charge())   cpls[i]=sqr(m_flavs[i].Charge());
      else if (m_flavs[i].IsPhoton()) cpls[i]=1.;
      else                            cpls[i]=0.;
    }
    msg_Debugging()<<m_flavs[i]<<": "<<cpls[i]<<std::endl;
  }
  p_kernel->SetCpls(cpls);
}

void KP_Terms::SetNLOMC(PDF::NLOMC_Base *const nlomc)
{
  p_nlomc=nlomc;
  m_subtype=p_nlomc->SubtractionType();
  p_kernel->SetSubType(m_subtype);
  if (m_subtype==1) p_kernel->SetKappa(1.0);
}

void KP_Terms::SetMassive()
{
  // determine whether massive,
  // count number of FS gluons/photons that can split into massive partons
  size_t fsgluons(0);
  for (size_t i=p_proc->NIn();i<m_flavs.size();i++) {
    if (m_flavs[i].IsMassive() &&
        ((m_stype==sbt::qcd && m_flavs[i].Strong()) ||
         (m_stype==sbt::qed && m_flavs[i].Charge()))) m_massive=true;
    if ((m_stype==sbt::qcd && m_flavs[i].IsGluon() ) ||
        (m_stype==sbt::qed && m_flavs[i].IsPhoton())) fsgluons++;
  }
  if (p_kernel->Nmf()>0) m_massive=true;

  m_xpa.resize(p_kernel->Nmf()*fsgluons);
  m_xpb.resize(p_kernel->Nmf()*fsgluons);
}

void KP_Terms::SetCoupling(const MODEL::Coupling_Map *cpls)
{
  std::string cplname("");
  if      (m_stype==sbt::qcd) cplname="Alpha_QCD";
  else if (m_stype==sbt::qed) cplname="Alpha_QED";
  else THROW(fatal_error,"Cannot set coupling for subtraction type"
                         +ToString(m_stype));

  if (cpls->find(cplname)!=cpls->end()) p_cpl=cpls->find(cplname)->second;
  else THROW(fatal_error,"Coupling not found");
  msg_Tracking()<<METHOD<<"(): "<<cplname<<" = "<<*p_cpl<<std::endl;
  m_cpldef = p_cpl->Default()/(2.*M_PI);
}

void KP_Terms::SetAlpha(const double &aff,const double &afi,
                        const double &aif,const double &aii)
{
  p_kernel->SetAlpha(aff,afi,aif,aii);
}

void KP_Terms::SetKappa(const double &kappa)
{
  p_kernel->SetKappa(kappa);
}

void KP_Terms::Calculate
(const Vec4D_Vector &mom,const std::vector<std::vector<double > > &dsij,
 const double &x0,const double &x1,const double &eta0,const double &eta1,
 const double &wgt, const bool &applyID)
{
  DEBUG_FUNC(m_itype<<": type(a)="<<m_typea<<", beam(a)="<<m_sa
                    <<", type(b)="<<m_typeb<<", beam(b)="<<m_sb);
  msg_Debugging()<<"x0="<<x0<<", x1="<<x1
                 <<", eta0="<<eta0<<", eta1="<<eta1<<std::endl;
  if (!m_sa && !m_sb) return;
  if ((m_sa && x0<eta0) || (m_sb && x1<eta1)) return;
  double cpl(Coupling());
  msg_Debugging()<<"cpl="<<Coupling()<<std::endl;
  size_t pls=1;
  if (m_sa&&m_sb) pls++;
  double muf2(p_proc->ScaleSetter()->Scale(stp::fac,1));
  for (int i=0;i<8;i++) m_kpca[i]=0.;
  for (int i=0;i<8;i++) m_kpcb[i]=0.;
  msg_Debugging()<<"parton list: "<<m_plist<<std::endl;

  // itype&K
  // K = as/2pi [ Ti^2*(Kbar+K_FS) - sum TiTa Kcal(i)
  //              - sum TaTk KbarM(k) + TaTb Ktilde ]
  // itype&P
  // P = as/2pi [ sum TaTk P(k)*log(muf2/xsak) + TaTb P*log(muf2/xsab) ]
  if (m_sa) {
    double w(1.-eta0);
    if (m_itype&cs_itype::K) {
      if (m_kcontrib&cs_kcontrib::Kb) {
        // Kbar terms
        msg_Debugging()<<"sa: Kbar"<<std::endl;
        m_kpca[0]=-w*p_kernel->Kb1(m_typea,x0)
                  +p_kernel->Kb2(m_typea)
                  -p_kernel->Kb4(m_typea,eta0);
        m_kpca[1]=w*(p_kernel->Kb1(m_typea,x0)
                     +p_kernel->Kb3(m_typea,x0));
        m_kpca[2]=-w*p_kernel->Kb1(m_typea+2,x0)
                  +p_kernel->Kb2(m_typea+2)
                  -p_kernel->Kb4(m_typea+2,eta0);
        m_kpca[3]=w*(p_kernel->Kb1(m_typea+2,x0)
                     +p_kernel->Kb3(m_typea+2,x0));
      }
      if (m_kcontrib&cs_kcontrib::KFS) {
        // KFS terms (only if not MSbar)
        if (m_facscheme>0) {
          msg_Debugging()<<"sa: KFS"<<std::endl;
          m_kpca[0]+=-w*p_kernel->KFS1(m_typea,x0)
                     +p_kernel->KFS2(m_typea)
                     -p_kernel->KFS4(m_typea,eta0);
          m_kpca[1]+=w*(p_kernel->KFS1(m_typea,x0)
                        +p_kernel->KFS3(m_typea,x0));
          m_kpca[2]+=-w*p_kernel->KFS1(m_typea+2,x0)
                     +p_kernel->KFS2(m_typea+2)
                     -p_kernel->KFS4(m_typea+2,eta0);
          m_kpca[3]+=w*(p_kernel->KFS1(m_typea+2,x0)
                        +p_kernel->KFS3(m_typea+2,x0));
        }
      }
      for (int i=0;i<4;i++) m_kpca[i]*=p_kernel->Cpl(0)*dsij[0][0];
      msg_Debugging()
          <<"    kpca[0]="<<m_kpca[0]<<" ,  kpca[1]="<<m_kpca[1]
          <<" ,  kpca[2]="<<m_kpca[2]<<" ,  kpca[3]="<<m_kpca[3]<<std::endl;

      if (m_kcontrib&cs_kcontrib::t) {
        size_t xpcnt(0);
        // (1/(1-x)_+ + delta(1-x)) terms
        msg_Debugging()<<"sa: t"<<std::endl;
        for (size_t i=pls;i<m_plist.size();i++) {
          int spin=m_flavs[m_plist[i]].IntSpin();
          if (m_stype==sbt::qed && m_flavs[m_plist[i]].IsMassive() &&
              m_flavs[m_plist[i]].IsVector()) {
            if (m_Vsubmode==0) spin=0;
            else               spin=1;
          }
          double saj=dabs(2.*mom[m_plist[0]]*mom[m_plist[i]]);
          double muq2=saj;
          double sajx = saj/x0;
          double muq2x=sajx;
          if (spin!=2) muq2=sqr(m_flavs[m_plist[i]].Mass())/saj;
          if (spin!=2) muq2x=sqr(m_flavs[m_plist[i]].Mass())/sajx;
          m_kpca[0]+=dsij[0][i]*(-w*p_kernel->t1(m_typea,spin,muq2,x0)
                                 +p_kernel->t2(m_typea,spin,muq2)
                                 -p_kernel->t4(m_typea,spin,muq2,eta0));
          m_kpca[1]+=dsij[0][i]*(w*(p_kernel->t1(m_typea,spin,muq2x,x0)
                                    +p_kernel->t3(m_typea,spin,muq2x,x0)));
          m_kpca[2]+=dsij[0][i]*(-w*p_kernel->t1(m_typea+2,spin,muq2,x0)
                                 +p_kernel->t2(m_typea+2,spin,muq2)
                                 -p_kernel->t4(m_typea+2,spin,muq2,eta0));
          m_kpca[3]+=dsij[0][i]*(w*(p_kernel->t1(m_typea+2,spin,muq2x,x0)
                                    +p_kernel->t3(m_typea+2,spin,muq2x,x0)));
          if (spin==2) {
            for (size_t j=0;j<p_kernel->Nmf();j++) {
              m_xpa[xpcnt].xp=1.-4.*sqr(p_kernel->FMass(j))/saj;
              if (m_xpa[xpcnt].xp>eta0) {
                m_kpca[0]+=dsij[0][i]*p_kernel->t6(m_typea,m_xpa[xpcnt].xp);
                m_kpca[1]+=dsij[0][i]*w*p_kernel->t5(m_typea,x0,m_xpa[xpcnt].xp);
                m_kpca[2]+=dsij[0][i]*p_kernel->t6(m_typea+2,m_xpa[xpcnt].xp);
                m_kpca[3]+=dsij[0][i]*w*p_kernel->t5(m_typea+2,x0,m_xpa[xpcnt].xp);

                m_xpa[xpcnt].kpc=dsij[0][i]*(-w*p_kernel->t5(m_typea,x0,m_xpa[xpcnt].xp)
                                             -p_kernel->t6(m_typea,m_xpa[xpcnt].xp)
                                             -p_kernel->t7(m_typea,eta0,m_xpa[xpcnt].xp));
              }
              xpcnt++;
            }
          }
        }
        msg_Debugging()
            <<"    kpca[0]="<<m_kpca[0]<<" ,  kpca[1]="<<m_kpca[1]
            <<" ,  kpca[2]="<<m_kpca[2]<<" ,  kpca[3]="<<m_kpca[3]<<std::endl;
      }

      if (m_kcontrib&cs_kcontrib::Kt) {
        // Ktilde terms
        if (m_sb) {
          msg_Debugging()<<"sa: Ktilde"<<std::endl;
          m_kpca[0]-=dsij[0][1]*(-w*p_kernel->Kt1(m_typea,x0)
                                 +p_kernel->Kt2(m_typea)
                                 -p_kernel->Kt4(m_typea,eta0));
          m_kpca[1]-=dsij[0][1]*w*(p_kernel->Kt1(m_typea,x0)
                                   +p_kernel->Kt3(m_typea,x0));
          m_kpca[2]-=dsij[0][1]*(-w*p_kernel->Kt1(m_typea+2,x0)
                                 +p_kernel->Kt2(m_typea+2)
                                 -p_kernel->Kt4(m_typea+2,eta0));
          m_kpca[3]-=dsij[0][1]*w*(p_kernel->Kt1(m_typea+2,x0)
                                   +p_kernel->Kt3(m_typea+2,x0));
        }
        msg_Debugging()
            <<"    kpca[0]="<<m_kpca[0]<<" ,  kpca[1]="<<m_kpca[1]
            <<" ,  kpca[2]="<<m_kpca[2]<<" ,  kpca[3]="<<m_kpca[3]<<std::endl;
      }
    }

    if (m_itype&cs_itype::P) {
      // P terms
      msg_Debugging()<<"sa: P"<<std::endl;
      double asum=0.,fsum=0.;
      for (size_t i=1;i<m_plist.size();i++) {
        fsum+=dsij[0][i];
        asum+=dsij[0][i]*log(muf2/dabs(2.*mom[m_plist[0]]*mom[m_plist[i]]));
      }
      double p4(-w*p_kernel->P1(m_typea,x0)
                +p_kernel->P2(m_typea)
                -p_kernel->P4(m_typea,eta0));
      double p5(w*(p_kernel->P1(m_typea,x0)
                   +p_kernel->P3(m_typea,x0)));
      double p6(-w*p_kernel->P1(m_typea+2,x0)
                +p_kernel->P2(m_typea+2)
                -p_kernel->P4(m_typea+2,eta0));
      double p7(w*(p_kernel->P1(m_typea+2,x0)
                   +p_kernel->P3(m_typea+2,x0)));
      // TODO: fix muf2 dependence
      m_kpca[0]+=asum*p4;
      m_kpca[1]+=asum*p5;
      m_kpca[2]+=asum*p6;
      m_kpca[3]+=asum*p7;
      m_kpca[4]=fsum*p4;
      m_kpca[5]=fsum*p5;
      m_kpca[6]=fsum*p6;
      m_kpca[7]=fsum*p7;
      msg_Debugging()
          <<"    kpca[0]="<<m_kpca[0]<<" ,  kpca[1]="<<m_kpca[1]
          <<" ,  kpca[2]="<<m_kpca[2]<<" ,  kpca[3]="<<m_kpca[3]<<std::endl;
    }

  if(applyID) // subtract only if pseudo-dipoles are used
    for (size_t i=0;i<m_plist.size();i++) {
      if(i==0 || i==1)  continue; // II-dipoles
      else {
      /* subtract contribution to K-terms, which comes from IF CS-dipoles
         and add contribution from IF pseudo-dipoles */
        ATOOLS::Vec4D n = mom[0]+mom[1];
        ATOOLS::Vec4D nx = mom[0]/x0+mom[1];
        /* convolved with quark-PDF */
//        if(m_typea==1){ // only for uu>WWbb and gg>WWbb test cases
          m_kpca[0]-=dsij[0][i]*(-w*p_kernel->K_if1(m_typea,x0)
                                 +p_kernel->K_if2(m_typea,mom[0],mom[m_plist[i]],n)
                                 -p_kernel->K_if4(m_typea,eta0));
          m_kpca[1]-=dsij[0][i]*(w*(p_kernel->K_if1(m_typea,x0)
                                    +p_kernel->K_if3(m_typea,x0,mom[0]/x0,mom[m_plist[i]],nx)));
//        }
        /* convolved with gluon-PDF */
//        else if(m_typea==2){ // only for uu>WWbb and gg>WWbb test cases
          m_kpca[2]-=dsij[0][i]*(-w*p_kernel->K_if1(m_typea+2,x0)
                                 +p_kernel->K_if2(m_typea+2,mom[0],mom[m_plist[i]],n)
                                 -p_kernel->K_if4(m_typea+2,eta0));
          m_kpca[3]-=dsij[0][i]*(w*(p_kernel->K_if1(m_typea+2,x0)
                                    +p_kernel->K_if3(m_typea+2,x0,mom[0],mom[m_plist[i]],n)));
//        }
      /* subtract contribution to K-terms, which comes from FI-CS-dipoles,
         i.e. eq.(8.29) in CS-paper */
        /* convolved with quark-PDF */
        if(m_typea==1){
          m_kpca[0]+=dsij[0][i]*(-w*p_kernel->K_fi1(1,x0)  // type (i.e. first argument) is always 1,
                                 +p_kernel->K_fi2(1)       // as the associated splitting is q > q
                                 -p_kernel->K_fi4(1,eta0));
          m_kpca[1]+=dsij[0][i]*(w*(p_kernel->K_fi1(1,x0)
                                   +p_kernel->K_fi3(1,x0)));
        }
        /* convolved with gluon-PDF */
        else if(m_typea==2){
          m_kpca[2]+=dsij[0][i]*(-w*p_kernel->K_fi1(1,x0)
                                 +p_kernel->K_fi2(1)
                                 -p_kernel->K_fi4(1,eta0));
          m_kpca[3]+=dsij[0][i]*(w*(p_kernel->K_fi1(1,x0)
                                   +p_kernel->K_fi3(1,x0)));
        }
        else THROW(fatal_error, "Invalid splitting type in KP-term subtraction.")
      }
    }
  }

  if (m_sb) {
    double w(1.-eta1);
    if (m_itype&cs_itype::K) {
      if (m_kcontrib&cs_kcontrib::Kb) {
        // Kbar terms
        msg_Debugging()<<"sb: Kbar"<<std::endl;
        m_kpcb[0]=-w*p_kernel->Kb1(m_typeb,x1)
                  +p_kernel->Kb2(m_typeb)
                  -p_kernel->Kb4(m_typeb,eta1);
        m_kpcb[1]=w*(p_kernel->Kb1(m_typeb,x1)
                     +p_kernel->Kb3(m_typeb,x1));
        m_kpcb[2]=-w*p_kernel->Kb1(m_typeb+2,x1)
                  +p_kernel->Kb2(m_typeb+2)
                  -p_kernel->Kb4(m_typeb+2,eta1);
        m_kpcb[3]=w*(p_kernel->Kb1(m_typeb+2,x1)
                     +p_kernel->Kb3(m_typeb+2,x1));
      }
      if (m_kcontrib&cs_kcontrib::KFS) {
        // KFS terms (only if not MSbar)
        if (m_facscheme>0) {
          msg_Debugging()<<"sb: KFS"<<std::endl;
          m_kpcb[0]+=-w*p_kernel->KFS1(m_typeb,x1)
                     +p_kernel->KFS2(m_typeb)
                     -p_kernel->KFS4(m_typeb,eta1);
          m_kpcb[1]+=w*(p_kernel->KFS1(m_typeb,x1)
                        +p_kernel->KFS3(m_typeb,x1));
          m_kpcb[2]+=-w*p_kernel->KFS1(m_typeb+2,x1)
                     +p_kernel->KFS2(m_typeb+2)
                     -p_kernel->KFS4(m_typeb+2,eta1);
          m_kpcb[3]+=w*(p_kernel->KFS1(m_typeb+2,x1)
                        +p_kernel->KFS3(m_typeb+2,x1));
        }
      }
      for (int i=0;i<4;i++) m_kpcb[i]*=p_kernel->Cpl(1)*dsij[0][0];
      msg_Debugging()
          <<"    kpcb[0]="<<m_kpcb[0]<<" ,  kpcb[1]="<<m_kpcb[1]
          <<" ,  kpcb[2]="<<m_kpcb[2]<<" ,  kpcb[3]="<<m_kpcb[3]<<std::endl;

      if (m_kcontrib&cs_kcontrib::t) {
        size_t xpcnt(0);
        // (1/(1-x)_+ + delta(1-x)) terms
        msg_Debugging()<<"sb: t"<<std::endl;
        for (size_t i=pls;i<m_plist.size();i++) {
          int spin=m_flavs[m_plist[i]].IntSpin();
          if (m_stype==sbt::qed && m_flavs[m_plist[i]].IsMassive() &&
              m_flavs[m_plist[i]].IsVector()) {
            if (m_Vsubmode==0) spin=0;
            else               spin=1;
          }
          double saj=dabs(2.*mom[m_plist[pls-1]]*mom[m_plist[i]]);
          double muq2=saj;
          double sajx = saj/x1;
          double muq2x=sajx;
          if (spin!=2) muq2=sqr(m_flavs[m_plist[i]].Mass())/saj;// mu-tilde
          if (spin!=2) muq2x=sqr(m_flavs[m_plist[i]].Mass())/sajx;// mu
          m_kpcb[0]+=dsij[pls-1][i]*(-w*p_kernel->t1(m_typeb,spin,muq2,x1)
                                     +p_kernel->t2(m_typeb,spin,muq2)
                                     -p_kernel->t4(m_typeb,spin,muq2,eta1));
          m_kpcb[1]+=dsij[pls-1][i]*(w*(p_kernel->t1(m_typeb,spin,muq2x,x1)
                                        +p_kernel->t3(m_typeb,spin,muq2x,x1)));
          m_kpcb[2]+=dsij[pls-1][i]*(-w*p_kernel->t1(m_typeb+2,spin,muq2,x1)
                                     +p_kernel->t2(m_typeb+2,spin,muq2)
                                     -p_kernel->t4(m_typeb+2,spin,muq2,eta1));
          m_kpcb[3]+=dsij[pls-1][i]*(w*(p_kernel->t1(m_typeb+2,spin,muq2x,x1)
                                        +p_kernel->t3(m_typeb+2,spin,muq2x,x1)));
          if (spin==2) {
            for (size_t j=0;j<p_kernel->Nmf();j++) {
              m_xpb[xpcnt].xp=1.-4.*sqr(p_kernel->FMass(j))/saj;
              if (m_xpb[xpcnt].xp>eta1) {
                m_kpcb[0]+=dsij[pls-1][i]*p_kernel->t6(m_typeb,m_xpb[xpcnt].xp);
                m_kpcb[1]+=dsij[pls-1][i]*w*p_kernel->t5(m_typeb,x1,m_xpb[xpcnt].xp);
                m_kpcb[2]+=dsij[pls-1][i]*p_kernel->t6(m_typeb+2,m_xpb[xpcnt].xp);
                m_kpcb[3]+=dsij[pls-1][i]*w*p_kernel->t5(m_typeb+2,x1,m_xpb[xpcnt].xp);

                m_xpb[xpcnt].kpc=dsij[pls-1][i]*(-w*p_kernel->t5(m_typeb,x1,m_xpb[xpcnt].xp)
                                                 -p_kernel->t6(m_typeb,m_xpb[xpcnt].xp)
                                                 -p_kernel->t7(m_typeb,eta1,m_xpb[xpcnt].xp));
              }
              xpcnt++;
            }
          }
        }
        msg_Debugging()
            <<"    kpcb[0]="<<m_kpcb[0]<<" ,  kpcb[1]="<<m_kpcb[1]
            <<" ,  kpcb[2]="<<m_kpcb[2]<<" ,  kpcb[3]="<<m_kpcb[3]<<std::endl;
      }

      if (m_kcontrib&cs_kcontrib::Kt) {
        // Ktilde terms
        if (m_sa) {
          msg_Debugging()<<"sb: Ktilde"<<std::endl;
          m_kpcb[0]-=dsij[0][1]*(-w*p_kernel->Kt1(m_typeb,x1)
                                 +p_kernel->Kt2(m_typeb)
                                 -p_kernel->Kt4(m_typeb,eta1));
          m_kpcb[1]-=dsij[0][1]*w*(p_kernel->Kt1(m_typeb,x1)
                                   +p_kernel->Kt3(m_typeb,x1));
          m_kpcb[2]-=dsij[0][1]*(-w*p_kernel->Kt1(m_typeb+2,x1)
                                 +p_kernel->Kt2(m_typeb+2)
                                 -p_kernel->Kt4(m_typeb+2,eta1));
          m_kpcb[3]-=dsij[0][1]*w*(p_kernel->Kt1(m_typeb+2,x1)
                                   +p_kernel->Kt3(m_typeb+2,x1));
        }
        msg_Debugging()
            <<"    kpcb[0]="<<m_kpcb[0]<<" ,  kpcb[1]="<<m_kpcb[1]
            <<" ,  kpcb[2]="<<m_kpcb[2]<<" ,  kpcb[3]="<<m_kpcb[3]<<std::endl;
      }
    }

    if (m_itype&cs_itype::P) {
      // P terms
      msg_Debugging()<<"sb: P"<<std::endl;
      double asum=0.,fsum=0.;
      for (size_t i=0;i<m_plist.size();i++) if (i!=pls-1) {
        fsum+=dsij[pls-1][i];
        asum+=dsij[pls-1][i]*log(muf2/dabs(2.*mom[m_plist[pls-1]]
                                             *mom[m_plist[i]]));
      }
      double p4(-w*p_kernel->P1(m_typeb,x1)
                +p_kernel->P2(m_typeb)
                -p_kernel->P4(m_typeb,eta1));
      double p5(w*(p_kernel->P1(m_typeb,x1)
                   +p_kernel->P3(m_typeb,x1)));
      double p6(-w*p_kernel->P1(m_typeb+2,x1)
                +p_kernel->P2(m_typeb+2)
                -p_kernel->P4(m_typeb+2,eta1));
      double p7(w*(p_kernel->P1(m_typeb+2,x1)
                   +p_kernel->P3(m_typeb+2,x1)));
      // TODO: fix muf2 dependence
      m_kpcb[0]+=asum*p4;
      m_kpcb[1]+=asum*p5;
      m_kpcb[2]+=asum*p6;
      m_kpcb[3]+=asum*p7;
      m_kpcb[4]=fsum*p4;
      m_kpcb[5]=fsum*p5;
      m_kpcb[6]=fsum*p6;
      m_kpcb[7]=fsum*p7;
      msg_Debugging()
          <<"    kpcb[0]="<<m_kpcb[0]<<" ,  kpcb[1]="<<m_kpcb[1]
          <<" ,  kpcb[2]="<<m_kpcb[2]<<" ,  kpcb[3]="<<m_kpcb[3]<<std::endl;
    }

  if(applyID) // subtract only if pseudo-dipoles are used
    for (size_t i=0;i<m_plist.size();i++) {
      if(i==0 || i==1)  continue;
      else {
        /* subtract contribution to K-terms, which comes from IF CS-dipoles
           and add contribution from IF pseudo-dipoles */
        ATOOLS::Vec4D n = mom[0]+mom[1];
        ATOOLS::Vec4D nx = mom[0]+mom[1]/x1;
        /* convolved with quark-PDF */
//        if(m_typeb==1){
          m_kpcb[0]-=dsij[1][i]*(-w*p_kernel->K_if1(m_typeb,x1)
                                 +p_kernel->K_if2(m_typeb,mom[1],mom[m_plist[i]],n)
                                 -p_kernel->K_if4(m_typeb,eta1));
          m_kpcb[1]-=dsij[1][i]*(w*(p_kernel->K_if1(m_typeb,x1)
                                   +p_kernel->K_if3(m_typeb,x1,mom[1]/x1,mom[m_plist[i]],nx)));
//        }
        /* convolved with gluon-PDF */
//        else if(m_typeb==2){
          m_kpcb[2]-=dsij[1][i]*(-w*p_kernel->K_if1(m_typeb+2,x1)
                                 +p_kernel->K_if2(m_typeb+2,mom[1],mom[m_plist[i]],n)
                                 -p_kernel->K_if4(m_typeb+2,eta1));
          m_kpcb[3]-=dsij[1][i]*(w*(p_kernel->K_if1(m_typeb+2,x1)
                                   +p_kernel->K_if3(m_typeb+2,x1,mom[1],mom[m_plist[i]],n)));
//        }
        /* subtract contribution to K-terms which comes from FI-dipoles,
           i.e. eq.(8.29) in CS-paper */
        /* convolved with quark-PDF */
        if(m_typeb==1){
          m_kpcb[0]+=dsij[1][i]*(-w*p_kernel->K_fi1(1,x1)
                                 +p_kernel->K_fi2(1)
                                 -p_kernel->K_fi4(1,eta1));
          m_kpcb[1]+=dsij[1][i]*(w*(p_kernel->K_fi1(1,x1)
                                   +p_kernel->K_fi3(1,x1)));
        }
        /* convolved with gluon-PDF */
        else if(m_typeb==2){
          m_kpcb[2]+=dsij[1][i]*(-w*p_kernel->K_fi1(1,x1)
                                 +p_kernel->K_fi2(1)
                                 -p_kernel->K_fi4(1,eta1));
          m_kpcb[3]+=dsij[1][i]*(w*(p_kernel->K_fi1(1,x1)
                                   +p_kernel->K_fi3(1,x1)));
        }
        else THROW(fatal_error, "Unvalid splitting type in KP-term subtraction.")
      }
    }
  }

  double cplfac(wgt*cpl);
  if (m_sa) cplfac/=(1.-eta0);
  if (m_sb) cplfac/=(1.-eta1);

  if (m_sa) {
    for (int i=0;i<8;i++) m_kpca[i]*=cplfac;
    for (size_t i=0;i<m_xpa.size();i++) m_xpa[i].kpc*=cplfac;
  }
  if (m_sb) {
    for (int i=0;i<8;i++) m_kpcb[i]*=cplfac;
    for (size_t i=0;i<m_xpb.size();i++) m_xpb[i].kpc*=cplfac;
  }

  msg_Debugging()<<"sa & sb: final"<<std::endl;
  msg_Debugging()
      <<"    kpca[0]="<<m_kpca[0]<<" ,  kpca[1]="<<m_kpca[1]
      <<" ,  kpca[2]="<<m_kpca[2]<<" ,  kpca[3]="<<m_kpca[3]<<std::endl;
  msg_Debugging()
      <<"    kpcb[0]="<<m_kpcb[0]<<" ,  kpcb[1]="<<m_kpcb[1]
      <<" ,  kpcb[2]="<<m_kpcb[2]<<" ,  kpcb[3]="<<m_kpcb[3]<<std::endl;
}

double KP_Terms::CalculateHP
(const Vec4D_Vector &mom, const double &cpl,
 const std::vector<size_t> &parton_indices, PHASIC::Color_Correlated_ME2* corr_me,
 const EXTAMP::DipoleCase &dipole_case, const ATOOLS::Flavour_Vector &flavs)
{ /* calculate H- and P- term for ID-dipole */

  double z = ATOOLS::ran->Get();
  m_H = 0;
  double H_temp; 
  ATOOLS::Vec4D n, pai, pb;  // momenta of BVI-kinematics
  ATOOLS::Vec4D nz, pa;      // z-dependent momenta
  ATOOLS::Vec4D pfs;         // momentum of other final state parton than ai (only used for IDin)
  size_t other_finalstate_parton; // other final state parton than ai (only used for IDin)

  switch(dipole_case){
  case EXTAMP::DipoleCase::IDa:
  case EXTAMP::DipoleCase::IDb:
    for(auto ai: parton_indices){    // sum over emitter-partons = ID-partons
      for(auto b: parton_indices){   // sum over all partons, except ai itself
        if(b==ai) continue;
  
        H_temp = 0;
        // labeling: 2=W+, 3=W-, 4=b, 5=bbar, is checked
        pa = z*mom[ai]; pb = mom[b]; pai = mom[ai];
        if(ai==4){                                    // b emitts
          switch(dipole_case){
          case EXTAMP::DipoleCase::IDa:
            n = mom[0]+mom[1]-pai-mom[3]-mom[5];
            nz = mom[0]+mom[1]-pa-mom[3]-mom[5];
            break;
          case EXTAMP::DipoleCase::IDb:
            n = mom[0]+mom[1]-pai-mom[2];
            nz = mom[0]+mom[1]-pa-mom[2];
          }
        }
        else if(ai==5){                               // bbar emitts
          switch(dipole_case){
          case EXTAMP::DipoleCase::IDa:
            n = mom[0]+mom[1]-pai-mom[2]-mom[4];
            nz = mom[0]+mom[1]-pa-mom[2]-mom[4];
            break;
          case EXTAMP::DipoleCase::IDb:
            n = mom[0]+mom[1]-pai-mom[3];
            nz = mom[0]+mom[1]-pa-mom[3]; 
          }
        }
        double Mab = corr_me->GetValue(ai,b);  // = <..|TbTai|..> = <..|TaiTb|..>
  
        double v = sqrt(1-(n.Abs2())*((pai+pb).Abs2())/pow((pai+pb)*n,2.));  // eq.(C.21) in CS
        H_temp += -( 0.25 + p_kernel->L2(1,pai,pb,n)
                          + p_kernel->L3(1,z,pa,pb,nz) );
  
        H_temp *= Mab;
        m_H    += H_temp;
      }
    }
    break;

  case EXTAMP::DipoleCase::IDin:
    /* parton_indices = (0,1,4,5) */
    for(auto ai: parton_indices){    // sum over emitter-partons = ID-partons
      if(ai<2) continue;             // emitter is in initial-state
      other_finalstate_parton = 9-ai;
      for(auto b: parton_indices){   // sum over all spectator-partons
        if(b==ai) continue;
  
        H_temp = 0;
        // labeling: 2=W+, 3=W-, 4=j, 5=j
        pa = z*mom[ai]; pai = mom[ai];
        pb = mom[b];    pfs = mom[other_finalstate_parton];
        if(flavs[ai].IsbQuark()) {
          n = mom[0]+mom[1]-pai-mom[3]-pfs;
          nz = mom[0]+mom[1]-pa-mom[3]-pfs;
          if(!flavs[other_finalstate_parton].IsbbarQuark()) 
            THROW(fatal_error,"Something went wrong in Assignment.");
        }
        else if (flavs[ai].IsbbarQuark()){
          n = mom[0]+mom[1]-pai-mom[2]-pfs;
          nz = mom[0]+mom[1]-pa-mom[2]-pfs;
          if(!flavs[other_finalstate_parton].IsbQuark()) 
            THROW(fatal_error,"Something went wrong in Assignment.");
        }
        else THROW(fatal_error,"Something went wrong in Assignment.");

        double Mab = corr_me->GetValue(ai,b);  // = <..|TbTai|..> = <..|TaiTb|..>
  
        double v = sqrt(1-(n.Abs2())*((pai+pb).Abs2())/pow((pai+pb)*n,2.));  // eq.(C.21) in CS
        H_temp += -( 0.25 + p_kernel->L2(1,pai,pb,n)
                          + p_kernel->L3(1,z,pa,pb,nz) );
  
        H_temp *= Mab;
        m_H    += H_temp;
      }
    }
  }
  return cpl*m_H;
}

  
double KP_Terms::Get(PDF::PDF_Base *pdfa, PDF::PDF_Base *pdfb,
                     const double &x0,const double &x1,
                     const double &eta0,const double &eta1,
                     const double &muf02,const double &muf12,
                     const double &muf02fac,const double &muf12fac,
                     const ATOOLS::Flavour &fl0,const ATOOLS::Flavour &fl1)
{
  DEBUG_FUNC(m_stype<<": "<<
             "fl0="<<fl0<<", x0="<<x0<<", eta0="<<eta0<<
                          ", muf02="<<muf02<<", muf02fac="<<muf02fac<<
             "\n                   "<<
             "fl1="<<fl1<<", x1="<<x1<<", eta1="<<eta1<<
                          ", muf12="<<muf12<<", muf12fac="<<muf12fac);
  if (m_sa && (pdfa==NULL || !pdfa->Contains(fl0))) return 0.0;
  if (m_sb && (pdfb==NULL || !pdfb->Contains(fl1))) return 0.0;
  if (!m_sa && !m_sb) return 0.;
  if ((m_sa && x0<eta0) || (m_sb && x1<eta1)) return 0.;
  size_t pls=1;
  if (m_sa&&m_sb) pls++;
  Flavour gluon(m_stype==sbt::qcd?kf_gluon:kf_photon);
  Flavour quark(kf_quark);
  double fa(0.),faq(0.),fag(0.),faqx(0.),fagx(0.);
  double fb(0.),fbq(0.),fbg(0.),fbqx(0.),fbgx(0.);
  double g2massq(0.);

  // assumption: a/a' = gluon/photon, quark only
  if (m_sa) {
    msg_Debugging()<<"sa"<<std::endl;
    if (m_cemode && eta0*rpa->gen.PBeam(0)[0]<fl0.Mass(true)) {
      msg_Tracking()<<METHOD<<"(): E < m ! ( "<<eta0*rpa->gen.PBeam(0)[0]
                    <<" vs. "<<fl0.Mass(true)<<" )"<<std::endl;
      return 0.0;
    }
    // 1/x f_a'(eta/x)
    pdfa->Calculate(eta0/x0,muf02*muf02fac);
    // a' = gluon/photon
    if (pdfa->Contains(gluon))                fagx = pdfa->GetXPDF(gluon)/eta0;
    // a' = quark (if a = quark)
    if (fl0.IsQuark() && pdfa->Contains(fl0)) faqx = pdfa->GetXPDF(fl0)/eta0;
    // a' = quark (if a = gluon,photon)
    // in QCD CF/CA already in kpca/b[1]
    // in QED -Q_i^2 needs to be multiplied in
    else {
      for (size_t i=0;i<quark.Size();i++)
        if (pdfa->Contains(quark[i])) {
          if (m_stype==sbt::qcd)              faqx += pdfa->GetXPDF(quark[i]);
          else                                faqx += sqr(quark[i].Charge())
                                                      *pdfa->GetXPDF(quark[i]);
        }
      faqx/=eta0;
    }
    // f_a(eta), f_a'(eta)
    pdfa->Calculate(eta0,muf02*muf02fac);
    if (pdfa->Contains(fl0)) fa = pdfa->GetXPDF(fl0)/eta0;
    // TODO: check whether this check is sensible (fa>0.)
    if (m_cemode && IsZero(fa,1.0e-16)) {
      msg_Tracking()<<METHOD<<"(): fa is zero, fa = "<<fa<<std::endl;
      return 0.;
    }
    if (!m_negativepdf && !(fa>0.)) {
      msg_Tracking()<<METHOD<<"(): fa is not pos. definite, fa = "<<fa<<std::endl;
      return 0.;
    }
    if (pdfa->Contains(gluon))                fag = pdfa->GetXPDF(gluon)/eta0;
    if (fl0.IsQuark() && pdfa->Contains(fl0)) faq = fa;
    // a' = quark (if a = gluon)
    // in QCD CF/CA already in kpca/b[1]
    // in QED -Q_i^2 needs to be multiplied in
    else {
      for (size_t i=0;i<quark.Size();i++)
        if (pdfa->Contains(quark[i])) {
          if (m_stype==sbt::qcd)              faq += pdfa->GetXPDF(quark[i]);
          else                                faq += sqr(quark[i].Charge())
                                                     *pdfa->GetXPDF(quark[i]);
        }
      faq/=eta0;
    }

    for (size_t i=0;i<m_xpa.size();i++) if (m_xpa[i].xp>eta0) {
      pdfa->Calculate(eta0/m_xpa[i].xp,muf02*muf02fac);
      g2massq+=m_xpa[i].kpc*pdfa->GetXPDF(fl0)/eta0/fa;
    }
  }

  if (m_sb) {
    msg_Debugging()<<"sb"<<std::endl;
    if (m_cemode && eta1*rpa->gen.PBeam(1)[0]<fl1.Mass(true)) {
      msg_Tracking()<<METHOD<<"(): E < m ! ( "<<eta1*rpa->gen.PBeam(1)[0]
                    <<" vs. "<<fl1.Mass(true)<<" )"<<std::endl;
      return 0.0;
    }
    // 1/x f_b'(eta/x)
    pdfb->Calculate(eta1/x1,muf12*muf12fac);
    // b' = gluon/photon
    if (pdfb->Contains(gluon))                fbgx = pdfb->GetXPDF(gluon)/eta1;
    // b' = quark (if b = quark)
    if (fl1.IsQuark() && pdfb->Contains(fl1)) fbqx = pdfb->GetXPDF(fl1)/eta1;
    // b' = quark (if b = gluon)
    // in QCD CF/CA already in kpca/b[1]
    // in QED -Q_i^2 needs to be multiplied in
    else {
      for (size_t i=0;i<quark.Size();i++)
        if (pdfb->Contains(quark[i])) {
          if (m_stype==sbt::qcd)              fbqx += pdfb->GetXPDF(quark[i]);
          else                                fbqx += sqr(quark[i].Charge())
                                                      *pdfb->GetXPDF(quark[i]);
        }
      fbqx/=eta1;
    }
    // f_a(eta), f_a'(eta)
    pdfb->Calculate(eta1,muf12*muf12fac);
    if (pdfb->Contains(fl1)) fb = pdfb->GetXPDF(fl1)/eta1;
    // TODO: check whether this check is sensible (fa>0.)
    if (m_cemode && IsZero(fb,1.0e-16)) {
      msg_Tracking()<<METHOD<<"(): fb is zero, fb = "<<fb<<std::endl;
      return 0.;
    }
    if (!m_negativepdf && !(fb>0.)) {
      msg_Tracking()<<METHOD<<"(): fb is not pos. definite, fb = "<<fb<<std::endl;
      return 0.;
    }
    if (pdfb->Contains(gluon))                fbg = pdfb->GetXPDF(gluon)/eta1;
    if (fl1.IsQuark() && pdfb->Contains(fl1)) fbq = fb;
    // b' = quark (if b = gluon)
    // in QCD CF/CA already in kpca/b[1]
    // in QED -Q_i^2 needs to be multiplied in
    else {
      for (size_t i=0;i<quark.Size();i++)
        if (pdfb->Contains(quark[i])) {
          if (m_stype==sbt::qcd)              fbq += pdfb->GetXPDF(quark[i]);
          else                                fbq += sqr(quark[i].Charge())
                                                     *pdfb->GetXPDF(quark[i]);
        }
      fbq/=eta0;
    }

    for (size_t i=0;i<m_xpb.size();i++) if (m_xpb[i].xp>eta1) {
      pdfb->Calculate(eta1/m_xpb[i].xp,muf12*muf12fac);
      g2massq+=m_xpb[i].kpc*pdfb->GetXPDF(fl1)/eta1/fb;
    }
  }

  msg_Debugging()<<"composing results"<<std::endl;
  double res(g2massq);
  // As this is intended to be a contribution to the *partonic* cross section,
  // multiplying with the two incoming parton PDFs should give the hadronic
  // cross section. Therefore, the following is used:
  //   (...)/fa * (fa * fb) = (...)*fb .
  // Note that this introduces an error, if fa or fb is 0.0. The only fix would
  // require this to be calculated elsewhere, e.g. from PHASIC's
  // Single_Process. But it would be semantically confusing to exclude the KP
  // from the ME generator's Partonic functions.
  if (m_sa && fa) {
    double a(m_kpca[0]*faq + m_kpca[1]*faqx + m_kpca[2]*fag + m_kpca[3]*fagx);
    if (muf02fac != 1.) {
      const double lF(log(muf02fac));
      a += (m_kpca[4]*faq + m_kpca[5]*faqx + m_kpca[6]*fag + m_kpca[7]*fagx) * lF;
    }
    a /= fa;
    res += a;
  }
  if (m_sb && fb) {
    double b(m_kpcb[0]*fbq + m_kpcb[1]*fbqx + m_kpcb[2]*fbg + m_kpcb[3]*fbgx);
    if (muf12fac != 1.) {
      const double lF(log(muf12fac));
      b += (m_kpcb[4]*fbq + m_kpcb[5]*fbqx + m_kpcb[6]*fbg + m_kpcb[7]*fbgx) * lF;
    }
    b /= fb;
    res += b;
  }
  if (msg_LevelIsDebugging()) {
    double precision(msg_Out().precision());
    msg->SetPrecision(16);
    msg_Out()<<"  "<<m_itype<<"-Terms Beam A:\n"
                   <<"    xa="<<x0<<" ,   etaa="<<eta0
                   <<"    fa="<<fa<<" ,  faq="<<faq<<" ,  faqx="<<faqx
                   <<" ,  fag="<<fag<<" ,  fagx="<<fagx<<"\n"
                   <<"    kpca[0]="<<m_kpca[0]<<" ,  kpca[1]="<<m_kpca[1]
                   <<" ,  kpca[2]="<<m_kpca[2]<<" ,  kpca[3]="<<m_kpca[3]
                   <<std::endl;
    msg_Out()<<"  => "<<(m_kpca[0]*faq+m_kpca[1]*faqx
                               +m_kpca[2]*fag+m_kpca[3]*fagx)*fb<<std::endl;
    msg_Out()<<"  "<<m_itype<<"-Terms Beam B:\n"
                   <<"    xb="<<x1<<" ,   etab="<<eta1
                   <<"    fb="<<fb<<" ,  fbq="<<fbq<<" ,  fbqx="<<fbqx
                   <<" ,  fbg="<<fbg<<" ,  fbgx="<<fbgx<<"\n"
                   <<"    kpcb[0]="<<m_kpcb[0]<<" ,  kpcb[1]="<<m_kpcb[1]
                   <<" ,  kpcb[2]="<<m_kpcb[2]<<" ,  kpcb[3]="<<m_kpcb[3]
                   <<std::endl;
    msg_Out()<<"  => "<<(m_kpcb[0]*fbq+m_kpcb[1]*fbqx
                               +m_kpcb[2]*fbg+m_kpcb[3]*fbgx)*fa<<std::endl;
    msg_Out()<<"res="<<res<<std::endl<<"  => "<<res*fa*fb<<std::endl;
    msg->SetPrecision(precision);
  }

  return res;
}

void KP_Terms::FillMEwgts(ATOOLS::ME_Weight_Info &wgtinfo)
{
  if (wgtinfo.m_type&mewgttype::KP) {
    size_t offset(m_stype==sbt::qed?16:0);
    for (int i=0;i<4;i++) {
      wgtinfo.m_wfac[i+offset]=wgtinfo.m_swap?m_kpcb[i]:m_kpca[i];
      wgtinfo.m_wfac[i+offset+4]=wgtinfo.m_swap?m_kpca[i]:m_kpcb[i];
      wgtinfo.m_wfac[i+offset+8]=wgtinfo.m_swap?m_kpcb[i+4]:m_kpca[i+4];
      wgtinfo.m_wfac[i+offset+12]=wgtinfo.m_swap?m_kpca[i+4]:m_kpcb[i+4];
    }
  }
}

