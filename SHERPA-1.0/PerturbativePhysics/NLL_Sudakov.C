#include "NLL_Sudakov.H"
#include "Message.H"
#include "MathTools.H"
#include "Run_Parameter.H"
#include "Running_AlphaS.H"
#include <iomanip> 
#include <stdio.h>  // sprintf

#include "NLL_Single_Sudakov.H"
#include "NLL_Combined_Sudakov.H"
#include "NLL_Branching_Probabilities.H"
#include "NLL_JetRate.H"
#include "MyStrStream.H"

namespace SHERPA {
  const double   Nc    = 3;
  const double   CA    = Nc;
  const double   CF    = (Nc*Nc-1.)/(2.*Nc);
  const double   Nf    = 5;
  const double   TR    =  1./2.;
  const double   BETA0 = (11.*CA-2.*Nf)/3.;
  const double   BETA1 = (17.*CA*CA- 3.*CF*Nf-5.*CA*Nf)/3.;
}

using namespace SHERPA;
using namespace ATOOLS;
using namespace MODEL;
using namespace std;

NLL_Sudakov::NLL_Sudakov(int mode, double _tmax,double _tmin,MODEL::Running_AlphaS * runas,int jetmode) 
  : m_mode(mode), p_runas(runas), p_jetrate(NULL)
{

  FixLambda2();

  msg.Debugging()<<"Init the NLL_Sudakov :"<<std::endl
		 <<"  lambda = "<<m_lambda<<"   --->  "<<std::endl
		 <<"  as(mue  ="<<sqrt(m_mu2)<<") = "<<m_asmu<<std::endl;
  // --------------------------------------------------


  // jetmode=0;  // 0, 3 ("1+2"), 4
  if (jetmode>=0)
    cout<<" jetmode="<<jetmode<<endl;
    p_jetrate = new NLL_JetRate(4,sqrt(_tmin),sqrt(_tmax),jetmode&7);
  //   p_jetrate = new NLL_JetRate(3,sqrt(0.001*_tmax),sqrt(_tmax),jetmode&7);

  if (jetmode==0 || jetmode==-2) PrepareMap(runas);
  else PrepareMassiveMap(runas,jetmode&7);

  CheckSudakovs(jetmode);
}


void NLL_Sudakov::PrepareMap(MODEL::Running_AlphaS * runas) 
{
  int  smode=Sudakov::analytic|(m_mode&3);
  BP::code    bpmode=(BP::code)(m_mode&124);
  if ((bpmode!=BP::gamma && bpmode!=BP::gamma_powercorr) || (runas!=0)) {
    smode=Sudakov::table;
    //     smode=Sudakov::numeric;
  }
  //  smode=Sudakov::numeric;

  //  BP::code    bpmode=BP::gamma_cut;
  // add dummy sudakov (const one)
  NLL_Dummy_Sudakov * dsud = new NLL_Dummy_Sudakov();
  m_sud_map[Flavour(kf::none)]=dsud;
  m_all_suds.push_back(dsud);  

  // add 5 massless quark sudakovs
  NLL_Single_Sudakov * ssud=0;
  ssud = new NLL_Single_Sudakov(new GammaQ_Lambda(bpmode,m_lambda,runas),smode);
  m_all_suds.push_back(ssud);
  for (int k=1;k<=5;++k) {
    Flavour fl = Flavour(kf::code(k));
    m_sud_map[fl]=ssud;
    m_sud_map[fl.Bar()]=ssud;
  }

  // add gluon sudakov (g->gg & g->qqb (5 fl) massless)
  NLL_Combined_Sudakov * csud=0;
  csud = new NLL_Combined_Sudakov(smode);
  csud->Add(new NLL_Single_Sudakov(new GammaG_Lambda(bpmode,m_lambda,runas),smode));
  csud->Add(new NLL_Single_Sudakov(new GammaF_Lambda(bpmode,m_lambda,runas),smode));
  m_all_suds.push_back(csud);
  m_sud_map[Flavour(kf::gluon)]=csud;

  // initialize nll jetrates
  if (p_jetrate) {
    p_jetrate->SetSudakovs(ssud,csud);
    p_jetrate->SetGammas(new GammaQ_Lambda(bpmode,m_lambda,runas),
			 new GammaG_Lambda(bpmode,m_lambda,runas),
			 new GammaF_Lambda(bpmode,m_lambda,runas));
    p_jetrate->Init();
    double r[6];
    p_jetrate->Rates(r[2],r[3],r[4],r[5]);
  }
}

void NLL_Sudakov::PrepareMassiveMap(MODEL::Running_AlphaS * runas, int mode) 
{
  NLL_Branching_Probability_Base * bp=NULL;
  //  std::cout<<" prepare massive map "<<endl;

  Sudakov::code smode=Sudakov::numeric;
  smode=Sudakov::table;
  //  BP::code    bpmode=BP::gamma_powercorr;
  //  BP::code    bpmode=BP::gamma_cut;
  //  BP::code    bpmode=BP::gamma;
  BP::code    bpmode=(BP::code)(m_mode&124);

  // ----------------------------------------
  // add dummy sudakov (const one)
  NLL_Dummy_Sudakov * dsud = new NLL_Dummy_Sudakov();
  m_sud_map[Flavour(kf::none)]=dsud;
  m_all_suds.push_back(dsud);  

  // ----------------------------------------
  // add 5 ( massive) quark sudakovs
  NLL_Single_Sudakov * ssud=0;
  for (int k=1;k<=5;++k) {
    Flavour fl = Flavour(kf::code(k));
    bp =    new GammaQ_Lambda_Massive(bpmode,m_lambda,runas,fl);
    //    cout<<" gamma("<<fl<<")="<<bp->Gamma(20.,2000.)<<endl;
    ssud = new NLL_Single_Sudakov(bp,smode);
    m_all_suds.push_back(ssud);
    m_sud_map[fl]=ssud;
    m_sud_map[fl.Bar()]=ssud;
  }
  NLL_Single_Sudakov * bsud=ssud;

  // ----------------------------------------
  // add gluon sudakov (g->gg & 5x g->qqb (5 fl))
  NLL_Combined_Sudakov * csud=0;
  csud = new NLL_Combined_Sudakov(smode);
  bp = new GammaG_Lambda_Massive(bpmode,m_lambda,runas);
  //  cout<<" gamma("<<Flavour(kf::gluon)<<")="<<bp->Gamma(20.,2000.)<<endl;
  csud->Add(new NLL_Single_Sudakov(bp,smode));
  for (int k=1;k<=5;++k) {
    Flavour fl = Flavour(kf::code(k));
    bp = new GammaF_Lambda_Massive(bpmode,m_lambda,runas,fl);
    //    cout<<" gamma("<<fl<<")="<<bp->Gamma(20.,2000.)<<endl;
    csud->Add(new NLL_Single_Sudakov(bp,smode));
  }
  m_all_suds.push_back(csud);
  m_sud_map[Flavour(kf::gluon)]=csud;

  //  cout<<" added "<<m_all_suds.size()<<" objects to sudakov list "<<endl;


  // ----------------------------------------
  // only 4 light flavours
  int nf4=4;
  // provide gluon sudakov (g->gg & 1x g->qq(4nf) + gqqb )
  NLL_Combined_Sudakov * gsudm=0;
  gsudm = new NLL_Combined_Sudakov(smode);
  gsudm->Add(new NLL_Single_Sudakov(new GammaG_Lambda(bpmode,m_lambda,runas),smode));
  gsudm->Add(new NLL_Single_Sudakov(new GammaF_Lambda(bpmode,m_lambda,runas,nf4),smode));
  gsudm->Add(new NLL_Single_Sudakov(new GammaF_Lambda_Massive(bpmode,m_lambda,runas,Flavour(kf::b)),smode));
  m_all_suds.push_back(gsudm);

  //  cout<<" added ("<<m_all_suds.size()<<") alternative gluon sudakov to sudakov list "<<endl;
  

  // ----------------------------------------
  // provide massless gluon sudakov
  NLL_Combined_Sudakov * gsud=0;
  gsud = new NLL_Combined_Sudakov(smode);
  gsud->Add(new NLL_Single_Sudakov(new GammaG_Lambda(bpmode,m_lambda,runas),smode));
  gsud->Add(new NLL_Single_Sudakov(new GammaF_Lambda(bpmode,m_lambda,runas),smode));
  m_all_suds.push_back(gsud);

  //  cout<<" added ("<<m_all_suds.size()<<") massless gluon sudakov to sudakov list "<<endl;


  // ----------------------------------------
  // provide massless quark sudakov
  NLL_Single_Sudakov * qsud=0;
  qsud = new NLL_Single_Sudakov(new GammaQ_Lambda(bpmode,m_lambda,runas),smode);
  m_all_suds.push_back(qsud);

  //  cout<<" added ("<<m_all_suds.size()<<") massless quark sudakov to sudakov list "<<endl;



  // initialize nll jetrates
  if (p_jetrate) {
    // we need four parts 
    // a) starting with light quark, and only light quark produced 
    //       (nf in Gamma_f replaced by four (explicitely and in DeltaG), all the rest like in the massless case)
    // b) starting with light quark, and only b quark produced 
    //       (Gamma_f only for a single massive particle)
    // c) starting with b quark, and only light quark produced
    // d) starting with b quark, and only b quark produced



    if (mode ==4) p_jetrate->SetSudakovs(qsud,gsud,bsud,gsudm); // dq, dg, dQ, dG
    else          p_jetrate->SetSudakovs(qsud,gsud,bsud,csud); // dq, dg, dQ, dtG
    p_jetrate->SetGammas(new GammaQ_Lambda(bpmode,m_lambda,runas),
			 new GammaG_Lambda(bpmode,m_lambda,runas),
			 new GammaF_Lambda(bpmode,m_lambda,runas,nf4),
			 new GammaQ_Lambda_Massive(bpmode,m_lambda,runas,Flavour(kf::b)),
			 new GammaF_Lambda_Massive(bpmode,m_lambda,runas,Flavour(kf::b)));
    p_jetrate->Init();
    double r[6];
    p_jetrate->Rates(r[2],r[3],r[4],r[5]);
  }

}


void NLL_Sudakov::FixLambda2() 
{
  m_mu2    = sqr(Flavour(kf::Z).Mass());
  m_asmu   = (*as)(m_mu2);
  m_lambda = sqrt( m_mu2 * exp(-4.*M_PI/(BETA0 * m_asmu)));
//   cout<<" mu2="<<m_mu2<<endl;
//   cout<<"asmu="<<m_asmu<<endl;
//   cout<<"lam ="<<m_lambda<<endl;
//   // test only:
//   double t=m_mu2;
//   cout<<" astest = "<<4.*M_PI/(BETA0*log(t/sqr(m_lambda)));
//   t=sqr(91.2);
//   cout<<" aastest = "<<4.*M_PI/(BETA0*log(t/sqr(m_lambda)));
}                 



NLL_Sudakov::~NLL_Sudakov()
{
  for (size_t i=0; i<m_all_suds.size();++i) 
    delete m_all_suds[i];

  m_all_suds.clear();
  m_sud_map.clear();
  
  if (p_jetrate) delete p_jetrate;
}




void NLL_Sudakov::CheckSudakovs(int jetmode)
{
  return;

  if (!rpa.gen.Debugging()) return;
  /*
  cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl;
  cout<<"kfactor for e+ e- -> Hadrons"<<endl;
  cout<<"============================"<<endl;
  double mz2 = sqr(Flavour(kf::Z).Mass());
  double alpi= (*as)(mz2)/M_PI;
  cout<<" alpi="<<alpi<<endl;
  double r0 =1.;
  double r1 =1.;
  // [1] W.~Celmaster and R.~J.~Gonsalves, Phys.\ Rev.\ Lett.\  {\bf 44} (1980) 560.
  // msbar R2 = (2 Zeta[3]/3 - 11/12) Nf + (365/24 - 11 Zeta[3]);
  double nf=5;
  double r2 =-0.1152953978936 * nf + 1.985707398578;

  double D = r0 + alpi * r1 + sqr(alpi) * r2;
  cout<<" "<<setw(10)<<r0<<" "<<setw(10)<<alpi * r1<<" "<<setw(10)<<sqr(alpi) * r2<<endl;
  cout<<" D = "<<D<<endl;
  */

  int steps=16;


  BP::code    bpmode=(BP::code)(m_mode&124);


  GammaQ_Lambda_Massive gammaq(BP::gamma,m_lambda,0,Flavour(kf::b));
  GammaQ_Lambda gammad(bpmode,m_lambda,p_runas);
  GammaG_Lambda gammag(bpmode,m_lambda,p_runas);
  GammaF_Lambda gammaf(bpmode,m_lambda,p_runas);
  GammaQ_AlphaS gammab(BP::gamma,m_asmu,m_mu2);
  cout<<endl;

  cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl;
  cout<<" m_lambda="<<m_lambda<<endl;
  cout<<gammaq(0.912,91.2)<<endl;
  cout<<gammaq(9.12,91.2)<<endl<<endl;
  cout<<gammaq.AlphaS(sqr(0.912))<<endl;
  cout<<gammaq.AlphaS(sqr(9.12))<<endl;
  cout<<gammaq.AlphaS(sqr(91.2))<<endl<<endl;
  cout<<gammab.AlphaS(sqr(0.912))<<endl;
  cout<<gammab.AlphaS(sqr(9.12))<<endl;
  cout<<gammab.AlphaS(sqr(91.2))<<endl;
  cout<<endl;
  cout<<endl;

  /*
  Sudakov::code  smode=Sudakov::table;
  NLL_Sudakov_Base * ssud = m_all_suds[8];   //  6,7,8 -- 5mass, 4ml+1mass, 5mless
    //  NLL_Single_Sudakov * ssud = new NLL_Single_Sudakov(new GammaQ_Lambda_Massive(BP::gamma,m_lambda,0,Flavour(kf::d)),smode); 
  std::cout<<"============================================================"<<std::endl;
  std::cout<<"  check table (massive) vs. internal (analytic/nummeric)"<<endl; 
  std::cout<<"============================================================"<<std::endl;
  //  std::cout<<"  y    q0     Q              test    sudq(q0,Q)  sudg(q0,Q) gammad(q0,Q) gamma(q0,Q) gammaf(q0,Q)"<<std::endl;
  //  std::cout<<"  y    q0     Q              test    sudq(q0,Q)  sudg(q0,Q) gammad(q0,Q) gamma(q0,Q) gammaf(q0,Q)"<<std::endl;
  std::cout<<"  y    q0     Q              sudg    sudg(4mls)  sudg(5ml) gammad(q0,Q) gamma(q0,Q) gammaf(q0,Q)"<<std::endl;
  double qmax=rpa.gen.Ecms(); 
  double qmin=sqrt(rpa.gen.Ycut())*rpa.gen.Ecms();
  for (int i=0; i<=steps;++i) {
    double qa=qmax*exp(-double(i)/steps*log(qmax/qmin));
    std::cout<<setiosflags(ios::scientific);
    std::cout.precision(4);
    std::cout<<" "<<setw(7)<<sqr(qa/qmax);
    std::cout<<resetiosflags(ios::scientific);
    std::cout<<" "<<setw(6)<<qmin<<" "<<setw(6)<<qa;
    std::cout.precision(6);
    std::cout<<" "<<setw(10)<<(*m_all_suds[6])(qa,qmin);   
    std::cout<<" "<<setw(10)<<(*m_all_suds[7])(qa,qmin);   
    std::cout<<" "<<setw(10)<<(*m_all_suds[8])(qa,qmin);   
//     std::cout<<" "<<setw(10)<<(*ssud)(qa,qmin);   
//     std::cout<<" "<<setw(10)<<Delta(Flavour(kf::d))(qa,qmin);
//     std::cout<<" "<<setw(10)<<Delta(Flavour(kf::gluon))(qa,qmin);
    std::cout<<" "<<setw(10)<<gammad(qmin,qa);
    std::cout<<" "<<setw(10)<<gammag(qmin,qa);
    std::cout<<" "<<setw(10)<<gammaf(qmin,qa);
    std::cout<<endl;
  }
  delete ssud;
  cout<<endl;
  cout<<endl;
  */


  double Q  = 2000.;
  double Qfix = Q;
  double q0 = 0.01*Q;


  MyStrStream s;
  s<<"sudakovs_"<<bpmode<<"_"<<Q<<".dat"<<endl;
  std::string fname;
  s>>fname;  

    //    std::ofstream  nllf("nll_rates_b0.dat");
  std::ofstream  sudf(fname.c_str());

  //  steps=64;
  steps=16;
  //  steps=40;
  std::cout<<"============================================================"<<std::endl
	   <<"   Selftest in NLL_Sudakov"<<std::endl
	   <<"------------------------------------------------------------"<<std::endl
	   <<"               q0     Q   sud_d(q0,Q) sud_g(q0,Q)   sud_b(q0,Q) gammad(q0,Q) gamma(q0,Q) gammaf(q0,Q)"<<std::endl;
  for (int i=0; i<=steps;++i) {
    Q=0.01*Qfix*exp(double(i)/steps*log(100.));
    //    q0=0.01*Q*exp(double(i)/steps*log(100.));
    std::cout<<setiosflags(ios::scientific);
    std::cout.precision(4);
    std::cout<<" "<<setw(7)<<sqr(q0/Q);
    std::cout<<resetiosflags(ios::scientific);
    std::cout<<" "<<setw(6)<<q0<<" "<<setw(6)<<Q;
    std::cout.precision(6);
    std::cout<<" "<<setw(10)<<Delta(Flavour(kf::d))(Q,q0);
    std::cout<<" "<<setw(10)<<Delta(Flavour(kf::gluon))(Q,q0);
    std::cout<<" "<<setw(10)<<Delta(Flavour(kf::b))(Q,q0);
    std::cout<<" "<<setw(10)<<gammad(q0,Q);
    std::cout<<" "<<setw(10)<<gammag(q0,Q);
    std::cout<<" "<<setw(10)<<gammaf(q0,Q);
    std::cout<<endl;

    sudf<<setiosflags(ios::scientific);
    sudf.precision(4);
    sudf<<" "<<setw(7)<<sqr(q0/Q);
    sudf<<resetiosflags(ios::scientific);
    sudf<<" "<<setw(6)<<q0<<" "<<setw(6)<<Q;
    sudf.precision(6);
    sudf<<" "<<setw(10)<<Delta(Flavour(kf::d))(Q,q0);
    sudf<<" "<<setw(10)<<Delta(Flavour(kf::gluon))(Q,q0);
    sudf<<" "<<setw(10)<<Delta(Flavour(kf::b))(Q,q0);
    sudf<<" "<<setw(10)<<gammad(q0,Q);
    sudf<<" "<<setw(10)<<gammag(q0,Q);
    sudf<<" "<<setw(10)<<gammaf(q0,Q);
    sudf<<endl;
  }
  sudf<<endl;
  sudf.close();
  cout<<endl;

  // partonlevel corrected ALEPH data from [hep-ph/9808364]
  const int ne =  21;
  const double ym = -1.;
  const double yp = -3.;
  double xx[ne] = {-1.0, -1.1, -1.2, -1.3, -1.4, -1.5, -1.6, -1.7, -1.8, -1.9, -2.0, -2.1, -2.2, -2.3, -2.4, -2.5, -2.6, -2.7, -2.8, -2.9, -3.0};
  
  double y2[ne] = { 0.952,0.929,0.906,0.878,0.849,0.814,0.780,0.734,0.695,0.654,0.608,0.562,0.516,0.4725,0.427,0.390,0.356,0.321,0.287,0.252,0.229};
  double y3[ne] = {0.0542,0.0745,0.0948,0.122,0.149,0.180,0.213,0.244,0.271,0.298, 0.325,0.348,.366,.379,.393,.402,.397,.397,.393,.379,.362};
  double y4[ne] = {0.0,0.0,0.0,.0017,.0042,.0083,.0125,.0208,.0308,.0442, .0608,.0808,.100,.125,.146,.167,.192,.222,.236,.254,.271};




  steps=40;
  int start=-10;
  int stop=steps;
  if (jetmode>7) {
    steps=8;
    start=2;
    stop=6;
  }
  if (p_jetrate) {
    BP::code    bpmode=(BP::code)(m_mode&124);
    //    int  smode=Sudakov::analytic|(m_mode&3);

    MyStrStream s;
    s<<"nll_rates_"<<bpmode<<"_"<<jetmode<<".dat"<<endl;
    std::string fname;
    s>>fname;
    

    //    std::ofstream  nllf("nll_rates_b0.dat");
    std::ofstream  nllf(fname.c_str());


    std::cout<<"===================================================================="<<std::endl
	     <<"   Selftest in NLL_Rates"<<std::endl
	     <<"--------------------------------------------------------------------"<<std::endl
	     <<"               q0     Q      R2         R3        R4        R5 "<<std::endl;
    nllf<<"#===================================================================="<<std::endl
	<<"#   Selftest in NLL_Rates"<<std::endl
	<<"#--------------------------------------------------------------------"<<std::endl
	<<"#               q0     Q      R2         R3        R4        R5 "<<std::endl;

    double chi2=0;

    for (int i=start; i<=stop;++i) {
      q0=0.01*Q*exp(double(i)/steps*log(100.));
      double r[6];
      p_jetrate->SetQ0(q0);
      p_jetrate->SetQ1(Q);
      p_jetrate->Rates(r[2],r[3],r[4],r[5]);
      if (jetmode>7) {
	double ytest=-4. + double(i)/steps*4.;
	int    xtest=int((ne-1)*(ytest - ym)/(yp - ym) + 0.5);
	if (xtest>=0 && xtest<ne) {
	  nllf<<"# "<<ytest<<" "<<xtest;
	  nllf<<setiosflags(ios::scientific);
	  nllf.precision(4);
	  nllf<<" "<<setw(7)<<xx[xtest]<<resetiosflags(ios::scientific)<<"              ";
	  nllf.precision(6);
	  nllf<<" "<<setw(10)<<y2[xtest]<<" "<<setw(10)<<y3[xtest]<<" "<<setw(10)<<y4[xtest]<<endl;

	  chi2+=sqr(r[2]-y2[xtest])+ sqr(r[3]-y3[xtest])+sqr(r[4]-y4[xtest]);
	}
      }

      std::cout<<setiosflags(ios::scientific);
      std::cout.precision(4);
      std::cout<<" "<<setw(7)<<sqr(q0/Q);
      std::cout<<resetiosflags(ios::scientific);
      std::cout<<" "<<setw(6)<<q0<<" "<<setw(6)<<Q;
      std::cout.precision(6);
      std::cout<<" "<<setw(10)<<r[2];
      std::cout<<" "<<setw(10)<<r[3];
      std::cout<<" "<<setw(10)<<r[4];
      std::cout<<" "<<setw(10)<<r[5];
      std::cout<<endl;      


      nllf<<setiosflags(ios::scientific);
      nllf.precision(4);
      nllf<<" "<<setw(7)<<sqr(q0/Q)<<resetiosflags(ios::scientific)<<" "<<setw(6)<<q0<<" "<<setw(6)<<Q;
      nllf.precision(6);
      nllf<<" "<<setw(10)<<r[2]<<" "<<setw(10)<<r[3]<<" "<<setw(10)<<r[4]<<" "<<setw(10)<<r[5]<<endl;
    }
    if (chi2>0.) {
      nllf<<" #  chi2 = "<<chi2<<endl;
    }
  }

  std::cout<<"============================================================"<<std::endl;
}





